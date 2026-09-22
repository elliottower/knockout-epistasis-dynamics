"""Checks that the code a Modal container runs is the code that was hashed at launch.

The Modal wrappers copy their first-party files into the image under one root and hash them at
launch. Two things must then hold in every container, and the first does not imply the second:

  1. each copy's SHA-256 equals the launch's;
  2. each first-party module Python actually imports was loaded from that copy.

A same-named module earlier on the import path, or a package Modal mounts on its own, is
imported instead of the copy while the copy still hashes correctly. Both checks cover every
first-party Python file the image carries; `tests/test_modal_provenance.py` checks that these
files include every first-party module the wrappers import.
"""

import hashlib
import importlib
from pathlib import Path

# Every first-party file each Modal image carries, relative to the repository root.
ODE_ARM_FILES = (
    "data_utils.py",
    "grn_coalition_sweep.py",
    "composition_scorer.py",
    "scripts/paths.py",
    "scripts/ode_coalition_sweep.py",
    "scripts/run_batch2b_extra_models.py",
    "scripts/run_batch2_blind_sweep.py",
    "scripts/structural_analysis.py",
    "scripts/ode_engine.py",
    "scripts/run_ode_arm.py",
    "scripts/walsh_estimators.py",
    "scripts/attractors.py",
    "scripts/modal_ode_arm.py",
    "scripts/modal_provenance.py",
    "scripts/registered_runs.py",
    "scripts/verify_freeze.py",
    "scripts/record_registered_identities.py",
    "results/paper_number_reconciliation.json",
    "experiments/ode_v2/solver_primary.json",
    "experiments/ode_v2/solver_legacy.json",
    "experiments/ode_v2/classifier_primary.json",
)
BOOLEAN_AUDIT_FILES = (
    "data_utils.py",
    "grn_coalition_sweep.py",
    "composition_scorer.py",
    "scripts/paths.py",
    "scripts/ode_coalition_sweep.py",
    "scripts/run_batch2b_extra_models.py",
    "scripts/run_batch2_blind_sweep.py",
    "scripts/structural_analysis.py",
    "scripts/audit_boolean_cycle_handling.py",
    "scripts/modal_boolean_audit.py",
    "scripts/modal_provenance.py",
    "results/paper_number_reconciliation.json",
)


class ProvenanceError(Exception):
    """The container's code is not the code hashed at launch."""


def file_hashes(root: Path, files: tuple[str, ...]) -> dict[str, str]:
    return {f: hashlib.sha256((root / f).read_bytes()).hexdigest() for f in files}


def module_name(relative: str) -> str:
    """scripts/ode_engine.py -> scripts.ode_engine"""
    return relative.removesuffix(".py").replace("/", ".")


def verify(root: Path, files: tuple[str, ...], expected_hashes: dict[str, str]) -> None:
    found = file_hashes(root, files)
    differing = [f for f in files if found[f] != expected_hashes.get(f)]
    if differing:
        raise ProvenanceError(f"the copies under {root} differ from the launch's files: {differing}")
    for relative in files:
        if not relative.endswith(".py"):
            continue
        module = importlib.import_module(module_name(relative))
        loaded = Path(module.__file__).resolve()
        if loaded != (root / relative).resolve():
            raise ProvenanceError(f"{module.__name__} was imported from {loaded}, not from the hashed copy {root / relative}")
