import ast
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from scripts.modal_provenance import BOOLEAN_AUDIT_FILES, ODE_ARM_FILES, ProvenanceError, file_hashes, verify
from scripts.paths import PROJECT_ROOT

SHADOW = PROJECT_ROOT / "tests" / "fixtures" / "shadow"
# Importable in this environment without Modal: every shipped module except the wrappers themselves.
PLAIN = ("data_utils.py", "grn_coalition_sweep.py", "composition_scorer.py", "scripts/paths.py",
         "scripts/ode_coalition_sweep.py", "scripts/ode_engine.py", "scripts/run_ode_arm.py")


def first_party_imports(script: Path) -> set[str]:
    """Every project file reachable from `script` through absolute import statements, wherever
    they sit (inside functions, try blocks or `with image.imports()`), as paths from the root."""
    seen, stack = set(), [script.resolve()]
    while stack:
        path = stack.pop()
        relative = path.relative_to(PROJECT_ROOT).as_posix()
        if relative in seen:
            continue
        seen.add(relative)
        for node in ast.walk(ast.parse(path.read_text())):
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
                names = [node.module] + [f"{node.module}.{a.name}" for a in node.names]
            else:
                continue
            for name in names:
                candidate = PROJECT_ROOT / (name.replace(".", "/") + ".py")
                if candidate.exists():
                    stack.append(candidate.resolve())
    return seen


@pytest.mark.parametrize("wrapper,files", [("scripts/modal_ode_arm.py", ODE_ARM_FILES),
                                           ("scripts/modal_boolean_audit.py", BOOLEAN_AUDIT_FILES)])
def test_every_first_party_module_a_wrapper_imports_is_shipped_and_hashed(wrapper, files):
    imported = first_party_imports(PROJECT_ROOT / wrapper)
    assert "scripts/ode_engine.py" in imported or "scripts/audit_boolean_cycle_handling.py" in imported
    assert imported - set(files) == set()


def test_verify_accepts_the_hashed_tree():
    verify(PROJECT_ROOT, PLAIN, file_hashes(PROJECT_ROOT, PLAIN))


def test_verify_refuses_a_copy_whose_contents_changed(tmp_path):
    hashes = file_hashes(PROJECT_ROOT, PLAIN)
    for relative in PLAIN:
        (tmp_path / relative).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(PROJECT_ROOT / relative, tmp_path / relative)
    with (tmp_path / "data_utils.py").open("a") as f:
        f.write("\n# edited\n")
    with pytest.raises(ProvenanceError, match="data_utils.py"):
        verify(tmp_path, PLAIN, hashes)


def verify_in_a_fresh_interpreter(python_path: list[Path], cwd: Path, first: str = "") -> subprocess.CompletedProcess:
    """A neutral working directory, because Python searches it before PYTHONPATH. `first` is a
    module imported before the check runs."""
    env = {**os.environ, "PYTHONPATH": os.pathsep.join(str(p) for p in python_path), "SHADOWED_ROOT": str(PROJECT_ROOT)}
    code = ((f"import {first}\n" if first else "") +
            "from pathlib import Path\n"
            "from scripts.modal_provenance import file_hashes, verify\n"
            f"files = {PLAIN!r}\n"
            f"verify(Path({str(PROJECT_ROOT)!r}), files, file_hashes(Path({str(PROJECT_ROOT)!r}), files))\n")
    return subprocess.run([shutil.which("python"), "-c", code], env=env, cwd=cwd,
                          capture_output=True, text=True, timeout=300)


def test_a_shadowed_module_is_refused_although_every_file_hashes_correctly(tmp_path):
    clean = verify_in_a_fresh_interpreter([PROJECT_ROOT], tmp_path)
    assert clean.returncode == 0, clean.stderr
    shadowed = verify_in_a_fresh_interpreter([SHADOW, PROJECT_ROOT], tmp_path)
    assert shadowed.returncode != 0
    assert "ProvenanceError" in shadowed.stderr
    assert f"data_utils was imported from {SHADOW / 'data_utils.py'}" in shadowed.stderr


def first_party_import_order(script: Path) -> list[str]:
    names = []
    for node in ast.parse(script.read_text()).body:
        if isinstance(node, ast.Import):
            names += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom) and node.module:
            names.append(node.module)
    return [n for n in names if (PROJECT_ROOT / (n.replace(".", "/") + ".py")).exists()]


def test_a_module_that_puts_the_root_first_on_the_path_hides_the_shadow_if_imported_first(tmp_path):
    # scripts/ode_coalition_sweep.py inserts the repository root at the front of the import path,
    # so imported before data_utils it loads the real file, and the check has nothing to refuse.
    hidden = verify_in_a_fresh_interpreter([SHADOW, PROJECT_ROOT], tmp_path, first="scripts.ode_coalition_sweep")
    assert hidden.returncode == 0, hidden.stderr


def test_the_modal_wrapper_imports_data_utils_before_any_other_project_module():
    # What keeps the Modal negative control meaningful, given the test above.
    assert first_party_import_order(PROJECT_ROOT / "scripts/modal_ode_arm.py")[0] == "data_utils"
