# Section 6.4 Experiment — Grenda vs Lexicographic Cell Selection
## Joystick geometry, ε_g = 0.30, ε_f = 0.10

---

### Setup

| Parameter | Value |
|---|---|
| Input geometry | `filedata/generatedMPs/joystick_approximation_fine_L3_NLO.xml` |
| Patches | 30 |
| Input level | 3 |
| ε_g | 0.30 |
| ε_f | 0.10 |
| Binary | `build\bin\Release\poissonTHB_example.exe` |

---

### Results summary

| Metric | Grenda (`--cell-method g`) | Lexicographic (`--cell-method l`) |
|---|---|---|
| **Total time** | **159.2 s** | **3174.7 s** |
| **Speedup (Grenda/lex)** | **20×** | — |
| Total FIT attempts | 43 | 630 |
| Accepted coarsenings | 25 (lev-2 + lev-1) | 25 (lev-2 + lev-1) |
| Final DoF | **585** | **585** |
| lev-2 FIT entries | 17 | 480 |
| lev-1 FIT entries | 8 | 120 |
| lev-0 FIT entries (all rejected) | 18 | 30 |
| Max globalError (accepted) | 0.2946 | 0.2946 |
| Max featureError (accepted) | 0.0847 | 0.0847 |

---

### Analysis

**Same result, 20× cost difference.**

Both methods reach the same final parameterization: 585 DoF, identical max errors. The difference is entirely in the number of rejected attempts before each accepted coarsening.

#### Why lex is slow

Lex selects cells in lexicographic order of cell index, cycling through all 30 patches before repeating. At lev-2 it made **480 FIT attempts** to achieve the same 17 accepted coarsenings that Grenda achieved in 17 attempts:

- For patches 2–5: 16 attempts each producing *unchanged* globalError (0.0296) — these patches have no useful lev-2 cells to coarsen but lex cycles through all of them anyway.
- For patch 6: 16 attempts before finding the first lev-2 cell whose removal stays within tolerance.

Grenda picks cells geometry-guided (delta ∈ {0.2, 0.4, 0.6, 0.8, 1.0, 1.2}) and only tries patches where a geometry criterion is satisfied, skipping barren patches entirely.

#### Lev-0 — interesting asymmetry

Lex at lev-0 tries all 30 patches once (attempt 0 each). Patches 0–5 actually pass ε_g (globalError ≈ 0.294 < 0.30) but are rejected by ε_f (featureError ≈ 0.245 > 0.10). Patches 6–29 are rejected by both.

Grenda at lev-0 begins with patch 1 (globalError = 0.576), which already fails ε_g, and cycles through other patches in a different order — none pass.

This means lex discovers that patches 0–5 at lev-0 are *geometrically close* (small globalError) but have poor boundary continuity. Grenda never probes these patches first.

#### Per-level FIT attempt counts

| Level | Grenda attempts | Lex attempts | Ratio |
|---|---|---|---|
| 2 | 17 | 480 | 28× |
| 1 | 8 | 120 | 15× |
| 0 | 18 | 30 | 1.7× |

At lev-0, lex is actually more efficient in terms of attempt count (30 vs 18) because it tries each patch exactly once. Grenda spends more attempts at lev-0 because it re-evaluates the same patches multiple times across different delta values before giving up.

---

### Conclusion for the paper (Section 6.4)

Grenda's geometry-based cell selection is **20× faster** than lexicographic for this geometry under ε_g = 0.30, ε_f = 0.10, while producing the **identical final result**. The speedup comes from Grenda pruning geometrically unsuitable patches before attempting a fit+jack cycle, whereas lex exhaustively cycles through all patches at each level regardless of geometry.

---

### Files

| Directory | Contents |
|---|---|
| `grenda/` | 159.2 s run: summary_log, logFile, cmd_output |
| `lex/` | 3174.7 s run: summary_log, logFile, cmd_output |
