# Provenance of pre-registrations and results

Fifteen commits dated 2026-08-02 add files written between 2026-07-28 and
2026-07-31. This document states plainly what that does and does not establish,
because a reader who notices the gap should not have to guess.

**What happened: the work was done on the dates the files carry, and it was not
committed to git until 2026-08-02.** There is no better explanation than that,
and none is offered.

---

## What each artifact's freeze actually rests on

The four groups differ, and conflating them would overstate three of them.

### 1. Cryptographically frozen, before any result existed

| artifact | first commit | date |
|---|---|---|
| `scripts/blind_prereg_brief.md` | `9e070ab` | 2026-07-28 08:13 |
| `prereg_shapiq_budget_sweep.md`, `_v2.md`, `_freeze_record.md` | `9e070ab` | 2026-07-28 08:13 |
| `prereg_grn_coalition_sweep.md`, `prereg_grn_v2/v3/v4.md` | `9e070ab` / `3966ce0` | 2026-07-28 |
| `prereg_auroc_addendum.md` | `3966ce0` | 2026-07-28 |
| `paper/prereg_extensions_v2.md` | `55cba8a` | 2026-07-28 |

These are ordinary pre-registrations. Git proves they existed before the work
they govern, and one of them is additionally hash-frozen.

**`prereg_shapiq_freeze_record.md` is a cryptographic freeze, and it verifies
today.** It was committed in `9e070ab` (the first commit) and records the
SHA-256 of `prereg_shapiq_budget_sweep_v2.md` together with version pins
(shapiq 1.6.0, numpy 2.2.6, scipy 1.15.3, scikit-learn 1.6.1) and verification
instructions:

```
recorded  2caaffff9590e906305eb4fc2c6ee1e8a797fcfec31d4bc86f6b7782b985a976
actual    2caaffff9590e906305eb4fc2c6ee1e8a797fcfec31d4bc86f6b7782b985a976
```

Anyone can re-run `sha256sum prereg_shapiq_budget_sweep_v2.md` and check it
against the hash in a commit that predates the analysis.

The `prereg_grn_*` documents carry no self-hash; their freeze is the git commit
itself, which is equivalent in force -- a commit hash covers its contents -- and
their commits also predate the analyses they govern.

**The blinding protocol is in this group.** `scripts/blind_prereg_brief.md` --
the brief instructing the predicting agent that it "has NOT seen any
experimental results" and must predict from Boolean rules "and NOTHING ELSE" --
was committed in the repository's **first commit**, at 08:13 on 2026-07-28.

### 2. Documentary only: the blind predictions

`paper/prereg_predictions.md` was committed on 2026-08-02 (`8235586`). Its
filesystem mtime is **2026-07-28 20:07**.

Supporting evidence for the ordering, all of which a reader can check:

| evidence | value |
|---|---|
| blinding protocol committed to git | 2026-07-28 **08:13** |
| predictions file written (mtime) | 2026-07-28 **20:07** |
| first result file written (mtime) | 2026-07-28 **20:19** |
| remaining results written (mtime) | 2026-07-28 20:43, 20:43, 21:00 |

So the protocol that defines the blinding was cryptographically frozen roughly
twelve hours before the predictions were made, and the predictions precede the
first result by twelve minutes on the filesystem.

**This is weaker than a git freeze and is not claimed to be equivalent.** File
mtimes can be altered. A reader who does not trust them should treat the four
blind predictions as unregistered, and the honest consequence is stated in the
next section.

### 3. Results, committed retrospectively

Everything under `results/` was committed on 2026-08-02 with mtimes from
2026-07-29 to 07-31. Results carry no freeze requirement -- what matters is that
the prediction preceded them, which is §2's question, not theirs.

### 4. Properly frozen, nothing run

`NEW PREREG 1`--`4` (`3095900`, `62d7a59`, `0537ff6`, `26c2b2a`) were committed
on 2026-08-02 **before any script for them exists**. Each commit message says
so. These have the freeze property the others in §2 lack, and they are the arms
whose results have not been produced.

---

## The consequence, stated rather than hidden

If a reader discounts §2 entirely, the affected claim is the paper's report that
pre-registered structural predictions achieved **48% direction accuracy,
indistinguishable from chance**.

That number survives the discount for a reason worth stating: **it is a negative
result about our own predictions.** Prediction 1 predicted a positive
correlation and measured Spearman rho = -0.509 -- wrong, and in the opposite
direction. Retrospective commitment creates an incentive to overstate
pre-registered success; here it would have been an incentive to *delete* the
file. The predictions were kept and reported as failures.

The reader is still entitled to treat that 48% as an unregistered analysis, and
the paper should not lean on it as though §1-grade evidence supported it.

---

## What was changed on 2026-08-02, beyond committing

Two corrections were made the same day, both recorded in
`docs/NOVELTY_AUDIT_v1.md`:

- The v12 title claimed "recoverable from 1% of knockout combinations". No
  recovery experiment exists in the paper, and Faure et al. (2024) had already
  published that result. Reverted (`3305ba3`).
- `PLAN.md` claimed head-to-head benchmarking against exact ground truth is
  rare. Refuted four ways; narrowed to the claim that survives (`3305ba3`).

Both are corrections against our own interest, made before submission rather
than after review.

---

## Practice going forward

Pre-registrations are committed before the script that implements them exists,
one commit per arm, and the commit message states that nothing has been run.
`NEW PREREG 1`--`4` follow this. Results are committed against the SHA of the
pre-registration they answer.
