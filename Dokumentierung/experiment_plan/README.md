# Paper Experiment Plan
## "An Unrefinement Algorithm for Planar THB-spline Parameterizations"

Last updated: 2026-07-08

---

## Overview

11 planned experiments across 4 classes + 2 already completed (Class 4).

| Class | Topic | Geometry | # runs | Status |
|---|---|---|---|---|
| 1 | ε tolerance variation | Hexagon with wiggly boundary | 3 | **geometry missing** |
| 2 | Projection stage: LS / LO / NLO | `mask_L3_PatchOnly689_noR.xml` | 3 | **flags missing** |
| 3 | λ locality parameter | TBD (mask or joystick) | 3 | ready |
| 4 | Cell selection: Grenda vs lex | `joystick_approximation_fine_L3_NLO.xml` | 2 | **DONE** |

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

**Flags (all existing):**
```
--epsilon-g <value> --epsilon-f <value>
```
Feature set F is currently encoded in the geometry (which interfaces are marked as feature). Varying F may require a code change or a geometry variant with a different feature marking.

**Missing before running:**
- Hexagon geometry XML
- Clarify whether F can be toggled by flag or requires geometry variants

---

## Class 2 — Projection stage comparison (LS / LO / NLO)

**Purpose:** Show that each successive stage (local LS → LO → NLO) expands the set of successfully coarsened cases on the same input geometry.

**Geometry:** `mask_L3_PatchOnly689_noR.xml`  
**Why:** Confirmed to produce 8 LO events + 3 NLO events at natural weights (base lowered to 2e-4). The NLO events are coarsenings where LS alone was insufficient and LO alone was insufficient but NLO restored regularity.

**Proposed runs:**

| Run | Flags | Expected |
|---|---|---|
| 2a | `--ls-only` | Fewest accepted coarsenings |
| 2b | `--lo-only` (LO but no NLO) | More accepted than 2a |
| 2c | (default — full LS+LO+NLO) | Most accepted; matches paper result |

**Flags — MISSING, need to add:**
- `--ls-only` — skip LO and NLO entirely; accept only when FIT is already regular and within tolerance
- `--lo-only` — allow LO fallback but not NLO fallback

Note: `--skip-lo-fallback` (already exists) goes the other way: skips LO and jumps straight to NLO — not what we need here.

---

## Class 3 — Locality parameter λ

**Purpose:** Show the effect of the locality extension parameter on fitting quality and/or timing.

**Geometry:** TBD. Options:
- `mask_L3_PatchOnly689_noR.xml` (main paper geometry)
- `joystick_approximation_fine_L3_NLO.xml` (30-patch, tests scaling)

**Proposed runs:**

| Run | Flags | Meaning |
|---|---|---|
| 3a | `--local-fitting --lambda=0` | Tight bounding box (default local) |
| 3b | `--local-fitting --lambda=1` | 1 cell-width extension |
| 3c | `--local-fitting --lambda=2` | 2 cell-width extension |

**Flags (all existing):**
```
--local-fitting --lambda=<value>
```

**Known issue:** `--local-fitting` currently has a performance regression: defect correction fires even when local region equals global region (fallback case), making it slower than global fitting (1.74 s vs 1.38 s mean). This must be fixed before timing comparisons in this class are meaningful. (Fix target: skip defect correction on fallback; see `feedback_local_fitting_perf.md`.)

**Decision needed:** Which geometry to use for Class 3?

---

## Class 4 — Cell selection: Grenda vs lexicographic

**DONE.** Results archived in `Dokumentierung/2026-07-07_grenda_lex_comparison/`.

| Method | Time | FIT attempts | Final DoF |
|---|---|---|---|
| Grenda (`--cell-method g`) | 159.2 s | 43 | 585 |
| Lex (`--cell-method l`) | 3174.7 s | 630 | 585 |
| **Speedup** | **20×** | — | identical |

---

## What is missing before running Classes 1–3

| Item | Needed for | Status |
|---|---|---|
| Hexagon with wiggly boundary (XML) | Class 1 | Not created |
| Feature set F toggling | Class 1 | May need code change |
| `--ls-only` flag | Class 2 | Not implemented |
| `--lo-only` flag | Class 2 | Not implemented |
| `--local-fitting` perf regression fix | Class 3 timing | Known issue, not yet fixed |
| Class 3 geometry choice (mask vs joystick) | Class 3 | Decision pending |

---

## Existing CLI flags (full list as of 2026-07-08)

| Flag | Type | Default | Description |
|---|---|---|---|
| `--epsilon-g <v>` | real | -1 (disabled) | Global approximation tolerance |
| `--epsilon-f <v>` | real | -1 (disabled) | Feature boundary tolerance |
| `--cell-method <c>` | char | g | Cell selection: g=Grenda, l=lex, r=random, s=smallest |
| `--verbose` | bool | false | Verbose logging (indexInTHB, functionDescription, per-patch/interface breakdown, etc.) |
| `--local-fitting` | bool | false | Restrict LS/LO to local region around removed cell |
| `--lambda=<v>` | int | 0 | Locality extension in cell widths (used with --local-fitting) |
| `--skip-lo-fallback` | bool | false | Skip LO; go straight to NLO when minusnumber > 0 |
| `--lo-weight-scale <v>` | real | 1.0 | Scale LO uniformity+length weights |

---

## Runner script skeleton

Once geometry and flags are ready, run all experiments via:

```powershell
# Dokumentierung/run_experiments.ps1
$binary  = "build\bin\Release\poissonTHB_example.exe"
$joystick = "filedata/generatedMPs/joystick_approximation_fine_L3_NLO.xml"
$mask     = "filedata/generatedMPs/mask_L3_PatchOnly689_noR.xml"
$hexagon  = "filedata/generatedMPs/hexagon_wiggly.xml"   # not yet created

$experiments = @(
    # Class 1 — tolerance variation
    @{ name="class1_tight";    args="--epsilon-g 0.10 --epsilon-f 0.10 $hexagon" },
    @{ name="class1_medium";   args="--epsilon-g 0.30 --epsilon-f 0.10 $hexagon" },
    @{ name="class1_loose_f";  args="--epsilon-g 0.30 --epsilon-f 0.30 $hexagon" },

    # Class 2 — projection stage
    @{ name="class2_ls_only";  args="--ls-only $mask" },        # flag not yet implemented
    @{ name="class2_lo_only";  args="--lo-only $mask" },        # flag not yet implemented
    @{ name="class2_full";     args="$mask" },

    # Class 3 — locality
    @{ name="class3_lambda0";  args="--local-fitting --lambda=0 $mask" },
    @{ name="class3_lambda1";  args="--local-fitting --lambda=1 $mask" },
    @{ name="class3_lambda2";  args="--local-fitting --lambda=2 $mask" },
)

foreach ($exp in $experiments) {
    $outDir = "Dokumentierung\runs\$($exp.name)"
    if (Test-Path "$outDir\cmd_output.txt") { Write-Host "SKIP: $($exp.name)"; continue }
    New-Item -ItemType Directory -Path $outDir -Force | Out-Null
    $t0 = Get-Date
    $argList = $exp.args -split '\s+'
    & $binary @argList *> "$outDir\cmd_output.txt"
    $elapsed = [math]::Round(((Get-Date) - $t0).TotalSeconds, 1)
    Add-Content "$outDir\timing.txt" "Total: $elapsed s"
    Write-Host "DONE: $($exp.name) — $elapsed s"
}
```

Skip-if-done (`Test-Path`) ensures you can re-run the script after adding a missing experiment without re-running completed ones.

---

## Reproducibility note

Every run should produce at minimum:
- `cmd_output.txt` — stdout+stderr (redirect both: `*>`)
- `summary_log.txt` — generated by binary in working directory
- `logFile.txt` — generated by binary in working directory

After each run, archive these three files to the experiment's output directory.
