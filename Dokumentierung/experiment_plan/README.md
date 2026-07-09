# Paper Experiment Plan
## "An Unrefinement Algorithm for Planar THB-spline Parameterizations"

Last updated: 2026-07-09

---

## Overview

14 planned runs across 4 classes. Class 4 is complete.

| Class | Topic | Geometry | # runs | Status |
|---|---|---|---|---|
| 1 | ε tolerance variation | Hexagon with wiggly boundary | 3 | **geometry missing** |
| 2 | Projection stage: LS / LO / NLO | `mask_L3_PatchOnly689_noR.xml` | 3 | **ready — rebuild required** |
| 3 | λ locality parameter | mask (existing) + joystick (new) | 6 | **ready — rebuild required** |
| 4 | Cell selection: Grenda vs lex | `joystick_approximation_fine_L3_NLO.xml` | 2 | **DONE** |

**Binary state:** Commit `2c926d01` (2026-07-09) changed default algorithm behaviour (Case 3 removed). Binary must be rebuilt before running Classes 2 and 3.

---

## Class 1 — Tolerance variation (ε_g, ε_f, F)

**Purpose:** Show how tightening ε_g, ε_f, and the feature set F controls the tradeoff between DoF reduction and approximation quality.

**Geometry needed:** Hexagon with wiggly boundary — does not exist yet. Must be created in `fitting_mspline.cpp` or as a hand-crafted XML. The wiggly boundary is essential to stress ε_f and F simultaneously; a smooth boundary makes the two indistinguishable.

**Proposed runs:**

| Run | ε_g | ε_f | F | Expected effect |
|---|---|---|---|---|
| 1a | 0.10 | 0.10 | full boundary | Tight: few coarsenings accepted |
| 1b | 0.30 | 0.10 | full boundary | Medium: baseline |
| 1c | 0.30 | 0.30 | interior only | Loose on feature: more accepted at boundary |

**Flags:**
```
--epsilon-g <value> --epsilon-f <value>
```
Feature set F is encoded in the geometry (which interfaces are marked as feature). Varying F may require a geometry variant with different feature markings.

**Missing before running:**
- Hexagon geometry XML
- Clarify whether F can be toggled by flag or requires geometry variants

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

**Joystick runs (new — to be executed):**

| Run | Command |
|---|---|
| 3a | `poissonTHB_example.exe --local-fitting --lambda 0 filedata/generatedMPs/joystick_approximation_fine_L3_NLO.xml` |
| 3b | `poissonTHB_example.exe --local-fitting --lambda 1 filedata/generatedMPs/joystick_approximation_fine_L3_NLO.xml` |
| 3c | `poissonTHB_example.exe --local-fitting --lambda 2 filedata/generatedMPs/joystick_approximation_fine_L3_NLO.xml` |

Note: `--lambda` takes its value as the next argument (space-separated), not with `=`.

**Baseline for comparison:** The Grenda baseline run without `--local-fitting` is archived in `Dokumentierung/2026-07-06_grenda_joystick_baseline/`. Use that as the global-fitting reference for the joystick.

---

## Class 4 — Cell selection: Grenda vs lexicographic

**DONE.** Results archived in `Dokumentierung/2026-07-07_grenda_lex_comparison/`.

| Method | Time | FIT attempts | Final DoF |
|---|---|---|---|
| Grenda (`--cell-method g`) | 159.2 s | 43 | 585 |
| Lex (`--cell-method l`) | 3174.7 s | 630 | 585 |
| **Speedup** | **20×** | — | identical |

---

## What is missing before running

| Item | Needed for | Status |
|---|---|---|
| Rebuild binary | Classes 2, 3 | Required (commit 2c926d01 changed behaviour) |
| Hexagon with wiggly boundary (XML) | Class 1 | Not created |
| Feature set F toggling | Class 1 | May need code change or geometry variant |
| Re-confirm 8 LO + 3 NLO on current binary | Class 2 | Run 2c first |
| Joystick λ=0, 1, 2 runs | Class 3 | Not yet run |

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
