import dataclasses

import pytest

from scripts.boolean_estimators import SHARDS
from scripts.paths import PROJECT_ROOT
from scripts.verify_freeze import FROZEN_FILES, REGISTRATIONS, RELEVANT_PATHS, RepoFacts, plan_text, problems

PLAN, HEAD, TAG_OBJECT = "a" * 40, "b" * 40, "c" * 40


def clean(**changes) -> RepoFacts:
    facts = RepoFacts(
        plan_commits={r: PLAN[:12] for r in REGISTRATIONS}, prereg_check_exit={r: 0 for r in REGISTRATIONS},
        head=HEAD, tag_commit=HEAD, tag_is_annotated=True, tag_signature_verifies=True, tag_object=TAG_OBJECT,
        origin_tag_object=TAG_OBJECT, plan_commit=PLAN, plan_is_ancestor_of_head=True,
        changed_since_plan=REGISTRATIONS, registration_text_matches_plan={r: True for r in REGISTRATIONS},
        frozen_files={f: {"tracked_at_plan": True, "unchanged_since": True} for f in FROZEN_FILES},
        dirty=(), identities_match=True)
    return dataclasses.replace(facts, **changes)


def test_a_launch_from_the_tagged_freeze_commit_passes():
    assert problems(clean()) == []


@pytest.mark.parametrize("changes,expected", [
    ({"plan_commits": {**{r: PLAN[:12] for r in REGISTRATIONS}, REGISTRATIONS[1]: None}}, "is not frozen"),
    ({"plan_commits": {**{r: PLAN[:12] for r in REGISTRATIONS}, REGISTRATIONS[2]: "d" * 12}}, "different plan commits"),
    ({"prereg_check_exit": {**{r: 0 for r in REGISTRATIONS}, REGISTRATIONS[0]: 1}}, "fails prereg check"),
    ({"tag_commit": None}, "there is no tag"),
    ({"head": "e" * 40}, "is not the commit"),
    ({"tag_is_annotated": False}, "not an annotated tag"),
    ({"tag_signature_verifies": False}, "signature"),
    ({"origin_tag_object": None}, "origin does not hold"),
    ({"plan_is_ancestor_of_head": False}, "not an ancestor"),
    ({"changed_since_plan": (*REGISTRATIONS, "scripts/ode_engine.py")}, "other than the registrations"),
    ({"frozen_files": {"scripts/ode_engine.py": {"tracked_at_plan": False, "unchanged_since": False}}}, "not tracked"),
    ({"frozen_files": {"uv.lock": {"tracked_at_plan": True, "unchanged_since": False}}}, "has changed"),
    ({"dirty": ("?? scripts/data_utils.py",)}, "modified or untracked"),
    ({"identities_match": False}, "registered_identities.json differs"),
    ({"registration_text_matches_plan": {**{r: True for r in REGISTRATIONS}, REGISTRATIONS[1]: False}},
     "text above its log differs"),
])
def test_each_way_of_not_being_the_frozen_commit_is_refused(changes, expected):
    found = problems(clean(**changes))
    assert len(found) == 1
    assert expected in found[0]


def test_the_environment_the_local_steps_run_in_and_the_boolean_input_are_frozen():
    assert {"pyproject.toml", "uv.lock", "scripts/verify_freeze.py", "scripts/modal_ode_arm.py",
            "scripts/analyze_ode_rerun.py", "scripts/boolean_estimators.py"} <= set(FROZEN_FILES)
    assert str(SHARDS.parent.relative_to(PROJECT_ROOT)) in FROZEN_FILES
    assert str(SHARDS.parent.relative_to(PROJECT_ROOT)) in RELEVANT_PATHS
    assert {"results/ode_full", "results/grn_v2/ode_full"} <= set(FROZEN_FILES) & set(RELEVANT_PATHS)


def test_the_plan_text_ignores_the_freeze_header_and_the_log_and_nothing_else():
    draft = "# Q?\n\n**Status:** DRAFT — not frozen.\n\nH1 holds.\n\n---\n\n## Log\n\ncreated\n"
    frozen = ("# Q?\n\n**Status:** FROZEN at `851b7a96154f`\n**Plan sha256:** `ab`\n**Log:** 2 entries, head `09`\n"
              "**Frozen:** 2026-09-22\n\n"
              "H1 holds.\n\n---\n\n## Log\n\ncreated\nfrozen\n")
    assert plan_text(draft) == plan_text(frozen)
    assert plan_text(frozen.replace("H1 holds.", "H1 holds loosely.")) != plan_text(draft)
