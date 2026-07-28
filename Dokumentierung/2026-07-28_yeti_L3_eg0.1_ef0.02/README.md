# Yeti L3 — Grenda vs Lexicographic (ε_g=0.1, ε_f=0.02)
## Candidate geometry for Class 3 (global vs local λ)

Date: 2026-07-28

---

### Setup

| Parameter | Value |
|---|---|
| Input geometry | `filedata/generatedMPs/yeti_mp2_thb_approximation_fine_L3.xml` |
| Patches | 18 |
| Input level | 3 |
| ε_g | **0.1** |
| ε_f | **0.02** |
| Fitting | Global (no `--local-fitting`) |
| Binary | `build\bin\Release\poissonTHB_example.exe` |

```
# Grenda (default)
poissonTHB_example.exe <geometry> --epsilon-g 0.1 --epsilon-f 0.02

# Lexicographic
poissonTHB_example.exe <geometry> --epsilon-g 0.1 --epsilon-f 0.02 --cell-method l
```

---

### Results summary

| Metric | Grenda (`--cell-method g`) | Lexicographic (`--cell-method l`) |
|---|---|---|
| FIT log entries | 39 | — (pending) |
| Runtime | — | — |
| Final DoF | — | — |

Full Grenda log in `grenda/yeti_mp2_thb_approximation_fine_L3_summary_log.txt`.

---

### Grenda log structure (from summary_log)

| Level | Entries | Notes |
|---|---|---|
| lev-2 | 15 | Machine precision (~3e-14), attempts 0–14 |
| lev-1 | 6 | Attempts 0,1,2,3,4,7 — gap at 5,6 indicates rejected attempts not logged |
| lev-0 | 18 | Attempts 0,1,1,3,… — repeated attempt numbers indicate bulk Grenda selection |

Note: The unconstrained baseline (no epsilon flags) produced 37 entries with sequential
attempts 0–36. The constrained run (ε_g=0.1, ε_f=0.02) produced 39 entries with
non-sequential attempts at lev-0, reflecting different Grenda bulk-selection behaviour
under epsilon constraints.

---

### Candidate role

**Primary candidate: Class 3 (global vs local λ)**
- 18 patches → local region at λ=0 will be a genuine subset of the full 18-patch basis
- Resolves the Class 3/4 geometry conflict (both currently use joystick)
- Next step: run `--local-fitting --lambda 0 / 1 / 2` to confirm locality is real

See `experiment_plan/README.md` → "Geometry evaluation notes" for full epsilon threshold analysis.

---

### Subdirectories

| Directory | Contents |
|---|---|
| `grenda/` | Grenda run output files |
| `lex/` | Lex run output files (pending) |
