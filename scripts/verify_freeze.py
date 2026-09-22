"""Check that a launch runs from the frozen commit of the ODE rerun registrations.

    uv run python -m scripts.verify_freeze

The freeze (experiments/ode_v2/CONTEXT.md, Freeze and launch) is one sequence:

  1. commit the code, the shared context, the configuration files and the identities;
  2. commit the three registrations alone, the plan commit;
  3. run `prereg freeze` in each registration's directory, which names the plan commit;
  4. commit the three freeze headers alone;
  5. tag that commit with the signed annotated tag FREEZE_TAG, and push the branch and the tag.

A launch then runs from a checkout of the tagged commit. This check passes only if, in the
checkout it runs in:

  - the three registrations are frozen at one plan commit, and each passes `prereg check`;
  - HEAD is the tagged commit, the tag is annotated, its signature verifies, and origin holds
    the same tag object;
  - the plan commit is an ancestor of HEAD, and nothing but the three registrations changed
    between them;
  - every file the registrations rely on is tracked at the plan commit and unchanged since,
    in HEAD and in the working tree;
  - no relevant path holds a modified or untracked file;
  - experiments/ode_v2/registered_identities.json equals the identities computed from the tree.

`scripts/modal_ode_arm.py` runs the same check before every registered launch and writes the
tag and commits into every record's launch manifest; `scripts/analyze_ode_rerun.py` refuses a
record that names another commit. Writes results/audit/freeze_verification.json and raises
FreezeError if any check fails.
"""

import json
import re
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

from scripts.modal_provenance import ODE_ARM_FILES
from scripts.paths import AUDIT, PROJECT_ROOT
from scripts.record_registered_identities import OUT as IDENTITIES
from scripts.record_registered_identities import build

FREEZE_TAG = "ode-rerun-freeze"
REGISTRATIONS = (
    "experiments/2026-09-21_ode-primary-rerun/PREREG.md",
    "experiments/2026-09-21_graded-perturbation-rerun/PREREG.md",
    "experiments/2026-09-21_ode-sensitivity/PREREG.md",
)
# Everything the registrations name, run or read, beyond the files shipped to Modal: the
# environment the local steps run in, and the Boolean audit, whose per-trajectory shards are the
# Boolean arm's input. A directory entry covers every file under it.
BOOLEAN_AUDIT = "results/audit/boolean_cycle_handling"
FROZEN_FILES = tuple(sorted(set(ODE_ARM_FILES) | {
    "pyproject.toml",
    "uv.lock",
    BOOLEAN_AUDIT,
    "results/audit/attractor_grouping_validation.json",
    "results/audit/attractor_grouping_validation_cycle_grid_0.01-0.05.json",
    "experiments/ode_v2/CONTEXT.md",
    "experiments/ode_v2/registered_identities.json",
    "scripts/analyze_ode_rerun.py",
    "scripts/boolean_estimators.py",
    "scripts/audit_boolean_cycle_handling.py",
    "scripts/record_registered_identities.py",
    "scripts/validate_attractor_grouping.py",
    "scripts/verify_freeze.py",
    "tests/test_ode_engine.py",
    "tests/test_walsh_estimators.py",
    "tests/test_attractors.py",
    "tests/test_analyze_ode_rerun.py",
    "tests/test_modal_provenance.py",
    "tests/test_verify_freeze.py",
    "tests/fixtures/shadow/data_utils.py",
}))
# Paths that must hold no modified or untracked file, except the registrations' results.
RELEVANT_PATHS = ("scripts", "tests", "experiments/ode_v2", "data_utils.py", "grn_coalition_sweep.py",
                  "composition_scorer.py", "pyproject.toml", "uv.lock", BOOLEAN_AUDIT,
                  *(str(Path(r).parent) for r in REGISTRATIONS))
EXCLUDED_PATHS = tuple(f":(exclude){Path(r).parent}/results" for r in REGISTRATIONS)
FROZEN_AT = re.compile(r"^\*\*Status:\*\* FROZEN at `([0-9a-f]{12})`", re.M)
OUT = AUDIT / "freeze_verification.json"


class FreezeError(Exception):
    pass


@dataclass(frozen=True)
class RepoFacts:
    plan_commits: dict[str, str | None]
    prereg_check_exit: dict[str, int]
    head: str
    tag_commit: str | None
    tag_is_annotated: bool
    tag_signature_verifies: bool
    tag_object: str | None
    origin_tag_object: str | None
    plan_commit: str | None
    plan_is_ancestor_of_head: bool
    changed_since_plan: tuple[str, ...]
    frozen_files: dict[str, dict[str, bool]]
    dirty: tuple[str, ...]
    identities_match: bool


def git(root: Path, *args: str) -> subprocess.CompletedProcess:
    return subprocess.run(["git", *args], cwd=root, capture_output=True, text=True)


def tag_commit(root: Path, tag: str) -> str:
    found = git(root, "rev-parse", "--verify", "--quiet", f"refs/tags/{tag}^{{commit}}")
    if found.returncode:
        raise FreezeError(f"no tag {tag} in {root}")
    return found.stdout.strip()


def gather(root: Path = PROJECT_ROOT, tag: str = FREEZE_TAG) -> RepoFacts:
    plans = {r: (m.group(1) if (m := FROZEN_AT.search((root / r).read_text())) else None) for r in REGISTRATIONS}
    checks = {r: subprocess.run(["prereg", "check"], cwd=(root / r).parent, capture_output=True, text=True).returncode
              for r in REGISTRATIONS}
    head = git(root, "rev-parse", "HEAD").stdout.strip()
    found = git(root, "rev-parse", "--verify", "--quiet", f"refs/tags/{tag}^{{commit}}")
    tagged = found.stdout.strip() if found.returncode == 0 else None
    tag_object = git(root, "rev-parse", "--verify", "--quiet", f"refs/tags/{tag}").stdout.strip() or None
    remote = git(root, "ls-remote", "origin", f"refs/tags/{tag}").stdout.split()
    distinct = {c for c in plans.values() if c}
    plan = None
    if len(distinct) == 1 and None not in plans.values():
        resolved = git(root, "rev-parse", "--verify", "--quiet", f"{distinct.pop()}^{{commit}}")
        plan = resolved.stdout.strip() if resolved.returncode == 0 else None
    frozen = {}
    for path in FROZEN_FILES:
        tracked = plan is not None and git(root, "cat-file", "-e", f"{plan}:{path}").returncode == 0
        frozen[path] = {"tracked_at_plan": tracked,
                        "unchanged_since": tracked and git(root, "diff", "--quiet", plan, "--", path).returncode == 0}
    status = git(root, "status", "--porcelain", "--untracked-files=all", "--", *RELEVANT_PATHS, *EXCLUDED_PATHS)
    return RepoFacts(
        plan_commits=plans, prereg_check_exit=checks, head=head, tag_commit=tagged,
        tag_is_annotated=tag_object is not None and git(root, "cat-file", "-t", tag_object).stdout.strip() == "tag",
        tag_signature_verifies=tag_object is not None and git(root, "verify-tag", tag).returncode == 0,
        tag_object=tag_object, origin_tag_object=remote[0] if remote else None, plan_commit=plan,
        plan_is_ancestor_of_head=plan is not None and git(root, "merge-base", "--is-ancestor", plan, head).returncode == 0,
        changed_since_plan=tuple(git(root, "diff", "--name-only", plan, head).stdout.split()) if plan else (),
        frozen_files=frozen, dirty=tuple(line for line in status.stdout.splitlines() if line.strip()),
        identities_match=build() == json.loads((root / IDENTITIES.relative_to(PROJECT_ROOT)).read_text()),
    )


def problems(facts: RepoFacts, tag: str = FREEZE_TAG) -> list[str]:
    """Every way the checkout fails to be a launch from the frozen commit."""
    found = []
    for reg, commit in facts.plan_commits.items():
        if commit is None:
            found.append(f"{reg} is not frozen")
        elif facts.prereg_check_exit[reg] != 0:
            found.append(f"{reg} fails prereg check (exit {facts.prereg_check_exit[reg]})")
    named = {c for c in facts.plan_commits.values() if c}
    if len(named) > 1:
        found.append(f"the registrations name different plan commits: {sorted(named)}")
    elif named and facts.plan_commit is None and None not in facts.plan_commits.values():
        found.append("the plan commit the registrations name does not resolve")
    if facts.tag_commit is None:
        found.append(f"there is no tag {tag}")
    else:
        if facts.head != facts.tag_commit:
            found.append(f"HEAD {facts.head[:12]} is not the commit {tag} points to, {facts.tag_commit[:12]}")
        if not facts.tag_is_annotated:
            found.append(f"{tag} is not an annotated tag")
        if not facts.tag_signature_verifies:
            found.append(f"the signature on {tag} does not verify")
        if facts.origin_tag_object != facts.tag_object:
            found.append(f"origin does not hold {tag} as the same tag object")
    if facts.plan_commit is not None:
        if not facts.plan_is_ancestor_of_head:
            found.append("the plan commit is not an ancestor of HEAD")
        others = sorted(set(facts.changed_since_plan) - set(facts.plan_commits))
        if others:
            found.append(f"files other than the registrations changed since the plan commit: {others}")
        for path, state in facts.frozen_files.items():
            if not state["tracked_at_plan"]:
                found.append(f"{path} is not tracked at the plan commit")
            elif not state["unchanged_since"]:
                found.append(f"{path} has changed since the plan commit")
    if facts.dirty:
        found.append(f"relevant paths hold modified or untracked files: {list(facts.dirty)}")
    if not facts.identities_match:
        found.append(f"{IDENTITIES.relative_to(PROJECT_ROOT)} differs from the identities computed from the tree")
    return found


def require(root: Path = PROJECT_ROOT) -> dict:
    """Run the check, write its report, and return the launch's freeze record; raise if it fails."""
    facts = gather(root)
    found = problems(facts)
    report = {**asdict(facts), "problems": found, "passed": not found}
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2))
    if found:
        raise FreezeError("; ".join(found))
    return {"tag": FREEZE_TAG, "tag_object": facts.tag_object, "commit": facts.head, "plan_commit": facts.plan_commit}


def main() -> None:
    freeze = require()
    print(f"freeze verified: {FREEZE_TAG} at {freeze['commit'][:12]}, plan commit {freeze['plan_commit'][:12]}, "
          f"{len(FROZEN_FILES)} files, identities match")


if __name__ == "__main__":
    main()
