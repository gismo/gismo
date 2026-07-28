# Paper Experiment Plan
## "An Unrefinement Algorithm for Planar THB-spline Parameterizations"

Last updated: 2026-07-28

---

## Overview

14 planned runs across 4 classes. Class 4 is complete.

| Class | Topic | Geometry | # runs | Status |
|---|---|---|---|---|
| 1 | ε tolerance variation | Hexagon with wiggly boundary | 3 | **geometry missing** |
| 2 | Projection stage: LS / LO / NLO | `mask_L3_PatchOnly689_noR.xml` | 3 | **ready — rebuild required** |
| 3 | λ locality parameter | mask (existing) + joystick (new) | 6 | **global + λ=0 done; λ=1,2 pending** |
| 4 | Cell selection: Grenda vs lex | `joystick_approximation_fine_L3_NLO.xml` | 2 | **DONE** |

**Geometry conflict:** Classes 3 and 4 both use joystick. Candidate replacement: `yeti_mp2_thb_approximation_fine_L3.xml` — see "Geometry evaluation notes" below.

**Binary state:** Commit `2c926d01` (2026-07-09) changed default algorithm behaviour (Case 3 removed). Binary must be rebuilt before running Classes 2 and 3.

---

## Class 1 — Tolerance variation (ε_g, ε_f, F)

**Purpose:** Show how tightening ε_g, ε_f, and the feature set F controls the tradeoff between DoF reduction and approximation quality.

**Geometry:** `filedata/generatedMPs/hexagon_wiggly.xml` — **RECREATED 2026-07-22**.
- Base: `hexagon_3p_4l.xml` (3-patch, level-4 THB-spline, degree 1, circumradius 1)
- Distortion: outer boundary wave, eps=0.03, k=9 half-sine half-periods (= 4.5 full oscillations)
- No interior radial distortion (-a 0)
- Min det J after distortion: 0.064 (all 6 outer sides pushing outward)
- Created by `fitting_mspline.cpp` with `-a 0 --boundary-wave-eps 0.03 --boundary-wave-k 9`
- Wave formula: `sin(k * π * s)` (half-sine). k=9 gives 4.5 full oscillations > level-3 Nyquist (4).
  The half-sine guarantees both CPs adjacent to each corner have positive (outward) amplitude —
  no inward dents near vertices. A centroid-based outward-normal check ensures all 6 outer sides
  are pushed consistently outward (the previous full-sine formula had 2 sides with inverted normals).

**Runs (ε_f calibrated to the half-sine wave, eps=0.03, k=9):**

The feature error when coarsening a level-3 boundary cell is the L∞ distance between the
original wave and its piecewise-linear approximation at level 3 (8 cells per side, knots at
s = 0, 1/8, …, 1). Measured empirically from the summary log.

| Run | ε_g | ε_f | Expected effect |
|---|---|---|---|
| 1a tight | 0.10 | 0.010 | ε_f below feature error → boundary cells rejected |
| 1b medium | 0.10 | 0.020 | ε_f at or near feature error → boundary at threshold |
| 1c loose_f | 0.10 | 0.050 | ε_f >> feature error → boundary cells accepted |

**Note:** ε_g = 0.10 is intentionally loose (interior has no radial distortion; interior coarsenings
produce small errors). The contrast across runs comes from ε_f only.

**Flags:**
```
--epsilon-g <value> --epsilon-f <value>
```
Feature set F = full outer boundary (default). Toggling F to interior-only requires either a flag
(check poissonTHB_example) or a second geometry variant — deferred.

**Status: READY**

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

**Decision (2026-07-09):** Two geometries.
- `mask_L3_PatchOnly689_noR.xml` — compact, 10-patch geometry. Local region frequently covers all patches (fallback to global). Results already exist in `Dokumentierung/2026-05-31_session4_PatchOnly689/` (λ=0 only; other λ values can be extracted if that run used `--local-fitting`).
- `joystick_approximation_fine_L3_NLO.xml` — 30-patch, spatially spread. Local region is genuinely local at λ=0. New runs needed for λ=0,1,2.

**Why two geometries:** Mask shows the fallback regime (local = global on compact geometry). Joystick shows genuine locality effect (local region actually smaller than global).

**Cell method:** Grenda (`g`) — default, no flag needed.

**Joystick runs:**

| Run | Command | Status | Archive |
|---|---|---|---|
| global (baseline) | `poissonTHB_example.exe filedata/generatedMPs/joystick_approximation_fine_L3_NLO.xml` | **DONE** | `Dokumentierung/2026-07-09_joystick_global/` |
| 3a λ=0 | `poissonTHB_example.exe --local-fitting --lambda 0 filedata/generatedMPs/joystick_approximation_fine_L3_NLO.xml` | **DONE** | `Dokumentierung/2026-07-09_joystick_local_lambda0/` |
| 3b λ=1 | `poissonTHB_example.exe --local-fitting --lambda 1 filedata/generatedMPs/joystick_approximation_fine_L3_NLO.xml` | pending | — |
| 3c λ=2 | `poissonTHB_example.exe --local-fitting --lambda 2 filedata/generatedMPs/joystick_approximation_fine_L3_NLO.xml` | pending | — |

Note: `--lambda` takes its value as the next argument (space-separated), not with `=`.

**Global baseline result (2026-07-09 12:33, commit 2c926d01):**
42/42 accepted, final DoF = **348**, time = **194.23 s**, minusnumber = 0 throughout.
FIT stage sufficient — no LO or NLO events. Archived in `Dokumentierung/2026-07-09_joystick_global/`.

**Local λ=0 result (2026-07-09 16:49, commit 2c926d01):**
42/42 accepted, final DoF = **348**, time = **193.11 s**. Mean basisSelected = **48.9 %** (range 6.3–97.8 %).
Genuinely local — joystick's 30 patches give real spatial separation at λ=0. Identical acceptance
decisions and final DoF as global; slight error differences at level 2 (local residual ≠ global, but
within tolerance). Archived in `Dokumentierung/2026-07-09_joystick_local_lambda0/`.

**Prior Grenda baselines (epsilon-constrained — for reference only):**
- `Dokumentierung/2026-07-06_grenda_joystick_baseline/`: 25/25 accepted, 585 DoF, 265.8 s (epsilon_g=-1.0 old default blocked level-0)
- `Dokumentierung/2026-07-07_grenda_lex_comparison/grenda/`: 25/25 accepted, 585 DoF, 159.2 s (`--epsilon-g 0.30 --epsilon-f 0.10` blocks level-0)

**Important:** The global reference for Class 3 is `2026-07-09_joystick_global` (unconstrained, new binary). The older baselines are NOT comparable because they used epsilon constraints that blocked level-0 cells.

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

### Open problem: Classes 3 and 4 both use joystick

Currently `joystick_approximation_fine_L3_NLO.xml` is assigned to both Class 3 (λ locality) and
Class 4 (Grenda vs lex). Using the same geometry for two different paper sections weakens the
presentation. A second geometry is needed for one of the two classes.

### Candidate: yeti_mp2_thb_approximation_fine_L3.xml

Run 2026-07-28 (global, unconstrained, no flags):
- **37 accepted coarsenings**, 18 patches, final `mpbes.size() = 117`, runtime **18.4 s**
- Log: `yeti_mp2_thb_approximation_fine_L3_summary_log.txt`

```
Command used:
poissonTHB_example.exe filedata/generatedMPs/yeti_mp2_thb_approximation_fine_L3.xml
```

**Why yeti L3 is promising:**

1. **Both ε_g and ε_f are meaningful constraints.** Unlike TV (featureError negligible) and hexagon
   (featureError dominant), yeti L3 shows genuine two-dimensional sensitivity:
   - At lev-1: globalErrors (0.018–0.057) > featureErrors (0.013–0.023) → ε_g binds
   - At lev-0: globalErrors (0.057–0.086) > featureErrors (0.056–0.062) → ε_g still binds, but both
     are in the same decade — either can be the limiting constraint depending on chosen values

2. **Gradual error profile.** No hard cliff between levels: errors rise monotonically from machine
   precision (~1e-14) through three distinct non-trivial regimes:
   - lev-1 steps 1–2: errors 0.013–0.022 (fine discrimination zone)
   - lev-1 steps 3–5: errors 0.019–0.057 (mid zone)
   - lev-0 steps 0–15: errors 0.056–0.086 (bulk zone)

3. **ε_g and ε_f act as independent controls.** Seven distinct (ε_g, ε_f) configurations span
   qualitatively different coarsening outcomes (see epsilon analysis 2026-07-28 in conversation log).

**Epsilon decision thresholds (from log):**

| Threshold | ε_g boundary | ε_f boundary |
|---|---|---|
| Below first non-trivial lev-1 | < 0.018 | < 0.013 |
| Adds lev-1 step 1 | 0.018 | 0.013 |
| Adds lev-1 step 2 | 0.023 | 0.018 |
| Adds lev-1 step 3 | 0.052 | 0.019 |
| Adds all lev-1 (steps 4–5) | 0.057 | 0.023 |
| Adds lev-0 step 0 (patch 5) | 0.057 | 0.057 |
| Adds all lev-0 | 0.086 | 0.063 |

**Candidate role:** Class 1 (tolerance variation, replacing hexagon or TV), or Class 3 (λ locality,
replacing the second joystick use). The 37-step run at 18.4 s is fast enough for multiple λ values.

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
| Joystick global baseline | Class 3 | **DONE** — `Dokumentierung/2026-07-09_joystick_global/` |
| Joystick λ=0, 1, 2 runs | Class 3 | Pending (global baseline done) |

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
$hexagon  = "filedata/generatedMPs/hexagon_wiggly.xml"   # not yet created

$experiments = @(
    # Class 1 — tolerance variation (geometry missing)
    @{ name="class1_tight";       args=@("--epsilon-g","0.10","--epsilon-f","0.10",$hexagon) },
    @{ name="class1_medium";      args=@("--epsilon-g","0.30","--epsilon-f","0.10",$hexagon) },
    @{ name="class1_loose_f";     args=@("--epsilon-g","0.30","--epsilon-f","0.30",$hexagon) },

    # Class 2 — projection stage (rebuild required)
    @{ name="class2_ls_only";     args=@("--ls-only",$mask) },
    @{ name="class2_lo_only";     args=@("--lo-only",$mask) },
    @{ name="class2_full";        args=@($mask) },

    # Class 3 — locality, joystick (rebuild required; mask results already in Dokumentierung)
    @{ name="class3_joy_lambda0"; args=@("--local-fitting","--lambda","0",$joystick) },
    @{ name="class3_joy_lambda1"; args=@("--local-fitting","--lambda","1",$joystick) },
    @{ name="class3_joy_lambda2"; args=@("--local-fitting","--lambda","2",$joystick) },
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
