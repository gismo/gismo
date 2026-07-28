# Paper Experiment Plan
## "An Unrefinement Algorithm for Planar THB-spline Parameterizations"

Last updated: 2026-07-28

---

## Overview

14 planned runs across 4 classes. Class 4 is complete.

| Class | Topic | Geometry | # runs | Status |
|---|---|---|---|---|
| 1 | ε tolerance variation | `tv_approximation_fine_L3.xml` | 3 | **READY** |
| 2 | Projection stage: LS / LO / NLO | `mask_L3_PatchOnly689_noR.xml` | 3 | **ready — rebuild required** |
| 3 | λ locality parameter | mask (existing) + **yeti L3** (new) | 6 | **global baseline done (unarchived); λ=0,1,2 pending** |
| 4 | Cell selection: Grenda vs lex | `joystick_approximation_fine_L3_NLO.xml` | 2 | **DONE** |

**Binary state:** Commit `2c926d01` (2026-07-09) changed default algorithm behaviour (Case 3 removed). Binary must be rebuilt before running Classes 2 and 3.

---

## Class 1 — Tolerance variation (ε_g)

**Purpose:** Show how tightening ε_g controls the tradeoff between DoF reduction and approximation quality.

**Geometry:** `filedata/generatedMPs/tv_approximation_fine_L3.xml`

Global baseline archived in `Dokumentierung/2026-07-25_tv_L3_global/` (34 accepted, unconstrained, 2026-07-25).

**Error profile:**
- lev-2 (13 steps): globalError 0.066–0.075, featureError ~1e-5 (negligible)
- lev-1 (5 steps): globalError 0.376 — **hard cliff** (5× jump from lev-2)
- lev-0 patches 0–10 (11 steps): globalError 0.376
- lev-0 patches 11–15 (5 steps): globalError 0.667–0.701

ε_f is the only active constraint axis — featureError is negligible throughout (max 0.049).
ε_f set to 0.1 in all runs to ensure it never binds.

**Runs:**

| Run | ε_g | ε_f | Expected effect |
|---|---|---|---|
| 1a tight | 0.10 | 0.10 | lev-2 only (13 accepted) — hard cliff blocks lev-1 |
| 1b medium | 0.40 | 0.10 | lev-2 + lev-1 + 10 lev-0 (28 accepted) |
| 1c loose | 0.75 | 0.10 | all 34 accepted |

**Flags:**
```
--epsilon-g <value> --epsilon-f 0.10
```

**Status: READY** (global baseline log available; constrained runs pending)

---

## Class 2 — Projection stage comparison (LS / LO / NLO)

**Purpose:** Show that each successive stage (FIT → LO → NLO) expands the set of successfully coarsened cases.

**Geometry:** `filedata/generatedMPs/mask_L3_PatchOnly689_noR.xml`

**Why:** Confirmed on an archived binary (2026-05-31) to produce 8 LO events + 3 NLO events at natural weights. With Case 3 removed (commit `2c926d01`), the current binary re-enables this behaviour. **Needs a fresh run to re-confirm counts on current code.**

**Runs:**

| Run | Command | Expected |
|---|---|---|
| 2a | `poissonTHB_example.exe --ls-only mask_L3_PatchOnly689_noR.xml` | Fewest accepted coarsenings |
| 2b | `poissonTHB_example.exe --lo-only mask_L3_PatchOnly689_noR.xml` | More accepted than 2a |
| 2c | `poissonTHB_example.exe mask_L3_PatchOnly689_noR.xml` | Most accepted; should reproduce 8 LO + 3 NLO |

**Note on --ls-only behaviour:** When FIT is irregular (minusnumber>0), `--ls-only` withdraws immediately (same as old Case 3 behaviour). Without the flag, the default full algorithm passes irregular FIT results to LO, and LO failures with minusnumber>0 to NLO.

---

## Class 3 — Locality parameter λ

**Purpose:** Show the effect of the locality extension parameter on which cells are accepted and on timing.

**Geometries (decision 2026-07-28):**
- `mask_L3_PatchOnly689_noR.xml` — compact, 10-patch geometry. Local region frequently covers all patches (fallback to global). Results exist in `Dokumentierung/2026-05-31_session4_PatchOnly689/`.
- `yeti_mp2_thb_approximation_fine_L3.xml` — **18-patch, spatially spread. Accepted 2026-07-28.**
  Local region is genuinely smaller than global at λ=0. Replaces joystick (freed for Class 4 only).

**Why yeti L3:** 18 patches give real spatial separation. Gradual error profile (no hard cliff).
Both ε_g and ε_f are meaningful constraints. Global run: 37 accepted, 18.4 s (fast, suitable for λ sweeps).
See "Geometry evaluation notes" section for full analysis and epsilon threshold table.

**Cell method:** Grenda (`g`) — default, no flag needed.

**All runs: unconstrained (no epsilon flags).** Epsilon would pre-filter cells and mask the locality effect.

**Yeti L3 runs:**

| Run | Command | Status | Archive |
|---|---|---|---|
| global (baseline) | `poissonTHB_example.exe filedata/generatedMPs/yeti_mp2_thb_approximation_fine_L3.xml` | **done 2026-07-28 (not archived — re-run to archive)** | — |
| 3a λ=0 | `poissonTHB_example.exe --local-fitting --lambda 0 filedata/generatedMPs/yeti_mp2_thb_approximation_fine_L3.xml` | pending | — |
| 3b λ=1 | `poissonTHB_example.exe --local-fitting --lambda 1 filedata/generatedMPs/yeti_mp2_thb_approximation_fine_L3.xml` | pending | — |
| 3c λ=2 | `poissonTHB_example.exe --local-fitting --lambda 2 filedata/generatedMPs/yeti_mp2_thb_approximation_fine_L3.xml` | pending | — |

Note: `--lambda` takes its value as the next argument (space-separated), not with `=`.

**Global baseline result (2026-07-28, unconstrained, unarchived):**
37 accepted, final `mpbes.size() = 117`, time ≈ **18.4 s**, all FIT (no LO/NLO), minusnumber = 0.
Files overwritten by constrained run — re-run to produce a clean archived baseline.

**Bonus: Grenda vs lex on yeti L3 (ε_g=0.1, ε_f=0.02) — 2026-07-28:**
Archived in `Dokumentierung/2026-07-28_yeti_L3_eg0.1_ef0.02/`. Not a Class 3 experiment —
run to validate geometry choice and check cell-selection behaviour before λ experiments.
- Grenda: **20.6 s**, mpbes.size()=264, 19 accepted
- Lex: **253.6 s**, mpbes.size()=340, 400 FIT attempts — **12.3× slower, worse result**
  Unlike joystick (same final DoF for both methods), here Grenda achieves better coarsening.
  Reason: tighter ε_f (0.02 vs 0.10) makes cell ordering matter — Grenda's bulk selection finds
  a more favourable greedy path than lex's sequential scan.

---

## Class 4 — Cell selection: Grenda vs lexicographic

**DONE.** Results archived in `Dokumentierung/2026-07-07_grenda_lex_comparison/`.

| Method | Time | FIT attempts | Final DoF |
|---|---|---|---|
| Grenda (`--cell-method g`) | 159.2 s | 43 | 585 |
| Lex (`--cell-method l`) | 3174.7 s | 630 | 585 |
| **Speedup** | **20×** | — | identical |

---

## Geometry evaluation notes

### Geometry assignments (resolved 2026-07-28)

| Class | Geometry | Decision |
|---|---|---|
| 1 | `tv_approximation_fine_L3.xml` | READY |
| 2 | `mask_L3_PatchOnly689_noR.xml` | READY (rebuild required) |
| 3 | `yeti_mp2_thb_approximation_fine_L3.xml` | **ACCEPTED 2026-07-28** |
| 4 | `joystick_approximation_fine_L3_NLO.xml` | DONE |

Classes 3 and 4 no longer share a geometry. Joystick remains in Class 4 (results archived).
Yeti L3 is the new Class 3 geometry.

---

### ACCEPTED: yeti_mp2_thb_approximation_fine_L3.xml (Class 3)

**Global baseline (unconstrained, 2026-07-28):**
37 accepted, mpbes.size()=117, 18.4 s, all FIT (no LO/NLO), minusnumber=0 throughout.

**Why accepted:**

1. **Both ε_g and ε_f are meaningful constraints** — unlike TV (featureError negligible).
   - lev-1: globalErrors (0.018–0.057) > featureErrors (0.013–0.023) → ε_g is the binding constraint
   - lev-0: both in the same decade (globalError 0.057–0.086, featureError 0.056–0.062)

2. **Gradual error profile** — no hard cliff between levels:
   - lev-1 steps 1–2: errors 0.013–0.022
   - lev-1 steps 3–5: errors 0.019–0.057
   - lev-0: errors 0.056–0.086

3. **18 patches** → local region at λ=0 is genuinely smaller than global (like joystick at 30 patches).

4. **Fast** — 18.4 s baseline means λ=0,1,2 sweep completes in under 2 minutes total.

**Epsilon decision thresholds (unconstrained baseline log):**

| What gets accepted | ε_g must be | ε_f must be |
|---|---|---|
| lev-2 only (machine precision) | any | any |
| + lev-1 step 1 | ≥ 0.019 | ≥ 0.014 |
| + lev-1 step 2 | ≥ 0.023 | ≥ 0.018 |
| + lev-1 step 3 | ≥ 0.052 | ≥ 0.019 |
| + all lev-1 | ≥ 0.057 | ≥ 0.023 |
| + lev-0 step 0 | ≥ 0.057 | ≥ 0.057 |
| + all lev-0 | ≥ 0.086 | ≥ 0.063 |

---

### Rejected candidate: yeti_mp2_thb_approximation_fine_L4.xml

Run 2026-07-28 (global, unconstrained):
- **~1787 log entries**, runtime unknown (very long)
- First ~1745 entries: lev-3 cells at machine precision (~1e-14). All trivially accepted.
- Last 41 entries: lev-1 (20 steps, errors 0.036–0.066) and lev-0 (21 steps, errors 0.058–0.108).
- **12-order-of-magnitude gap** between trivial and meaningful steps.

**Why L4 is unsuitable:**
- The 1700+ machine-precision lev-3 coarsenings carry no geometric content — they are meaningless
  as paper results and inflate runtime dramatically.
- ε cannot meaningfully distinguish between the trivial and meaningful regimes: any ε above 1e-13
  accepts all trivial steps; any ε below 0.036 rejects all meaningful steps. There is no useful
  middle ground.
- Grenda vs lex comparison is meaningless when almost all steps are trivial (Grenda offers no
  advantage over lex for machine-precision steps where all cells are geometrically equivalent).

**Verdict: do not use L4 for any paper experiment.**

---

## What is missing before running

| Item | Needed for | Status |
|---|---|---|
| Rebuild binary | Classes 2, 3 | Required (commit 2c926d01 changed behaviour) |
| Hexagon with wiggly boundary (XML) | Class 1 | Not created |
| Feature set F toggling | Class 1 | May need code change or geometry variant |
| Re-confirm 8 LO + 3 NLO on current binary | Class 2 | Run 2c first |
| Yeti L3 global baseline (archive) | Class 3 | Done but unarchived — re-run unconstrained |
| Yeti L3 λ=0, 1, 2 runs | Class 3 | Pending |

---

## CLI flags (full list as of 2026-07-09)

| Flag | Type | Default | Description |
|---|---|---|---|
| `--epsilon-g <v>` | real | 1e6 (effectively off) | Global approximation tolerance |
| `--epsilon-f <v>` | real | 1.0 | Feature boundary tolerance |
| `--cell-method <c>` | char | `g` | Cell selection: g=Grenda, l=lex, r=random, s=smallest |
| `--verbose` | bool | false | Verbose logging |
| `--local-fitting` | bool | false | Restrict LS/LO to local region around removed cell |
| `--lambda <v>` | int | 0 | Locality extension in cell widths (space-separated: `--lambda 1`) |
| `--ls-only` | bool | false | FIT only — withdraw if irregular or over tolerance; no LO/NLO |
| `--lo-only` | bool | false | LO enabled after FIT failure, NLO disabled |
| `--skip-lo-fallback` | bool | false | Skip LO; jump straight to NLO when FIT irregular |
| `--lo-weight-scale <v>` | real | 1.0 | Scale LO uniformity+length weights |

---

## Runner script

```powershell
# Dokumentierung/experiment_plan/run_experiments.ps1
$binary   = "build\bin\Release\poissonTHB_example.exe"
$joystick = "filedata/generatedMPs/joystick_approximation_fine_L3_NLO.xml"
$mask     = "filedata/generatedMPs/mask_L3_PatchOnly689_noR.xml"
$tv       = "filedata/generatedMPs/tv_approximation_fine_L3.xml"
$yeti     = "filedata/generatedMPs/yeti_mp2_thb_approximation_fine_L3.xml"

$experiments = @(
    # Class 1 — tolerance variation (tv_approximation_fine_L3.xml)
    @{ name="class1_tight";       args=@("--epsilon-g","0.10","--epsilon-f","0.10",$tv) },
    @{ name="class1_medium";      args=@("--epsilon-g","0.40","--epsilon-f","0.10",$tv) },
    @{ name="class1_loose";       args=@("--epsilon-g","0.75","--epsilon-f","0.10",$tv) },

    # Class 2 — projection stage (rebuild required)
    @{ name="class2_ls_only";     args=@("--ls-only",$mask) },
    @{ name="class2_lo_only";     args=@("--lo-only",$mask) },
    @{ name="class2_full";        args=@($mask) },

    # Class 3 — locality, yeti L3 (rebuild required; mask results already in Dokumentierung)
    @{ name="class3_yeti_global";  args=@($yeti) },
    @{ name="class3_yeti_lambda0"; args=@("--local-fitting","--lambda","0",$yeti) },
    @{ name="class3_yeti_lambda1"; args=@("--local-fitting","--lambda","1",$yeti) },
    @{ name="class3_yeti_lambda2"; args=@("--local-fitting","--lambda","2",$yeti) },
)

foreach ($exp in $experiments) {
    $outDir = "Dokumentierung\runs\$($exp.name)"
    if (Test-Path "$outDir\cmd_output.txt") { Write-Host "SKIP: $($exp.name)"; continue }
    New-Item -ItemType Directory -Path $outDir -Force | Out-Null
    $t0 = Get-Date
    & $binary @($exp.args) *> "$outDir\cmd_output.txt"
    $elapsed = [math]::Round(((Get-Date) - $t0).TotalSeconds, 1)
    Add-Content "$outDir\timing.txt" "Total: $elapsed s"
    Write-Host "DONE: $($exp.name) — $elapsed s"
}
```

After each run, copy `*_summary_log.txt` and `*_logFile_poissonTHB_example.txt` from the working directory into the run's output directory.

---

## Reproducibility note

Every run produces at minimum:
- `cmd_output.txt` — stdout+stderr (`*>` redirect)
- `*_summary_log.txt` — generated by binary in working directory
- `*_logFile_poissonTHB_example.txt` — generated by binary in working directory

Archive all three per run. The summary log is the primary result (FIT/LO/NLO counts, errors, timings). The logFile is the full diagnostic trace.
