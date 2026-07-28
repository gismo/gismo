# TV L3 — Global Baseline (Unconstrained)
## Class 1 geometry: tolerance variation (ε_g, ε_f, F)

Date: 2026-07-25

---

### Setup

| Parameter | Value |
|---|---|
| Input geometry | `filedata/generatedMPs/tv_approximation_fine_L3.xml` |
| ε_g | — (unconstrained, default 1e6) |
| ε_f | — (unconstrained, default 1.0) |
| Cell method | Grenda (default) |
| Binary | `build\bin\Release\poissonTHB_example.exe` |

```
poissonTHB_example.exe filedata/generatedMPs/tv_approximation_fine_L3.xml
```

---

### Results

| Metric | Value |
|---|---|
| Accepted coarsenings | **34** |
| Projection stages | FIT only (no LO/NLO) |
| minusnumber | 0 throughout |

| Level | Entries | globalError | featureError |
|---|---|---|---|
| lev-2 | 13 | 0.0659–0.0753 | 1.9e-05–3.3e-04 (negligible) |
| lev-1 | 5 | 0.376 | 0.0094 |
| lev-0 (patches 0–10) | 11 | 0.376 | 0.014–0.049 |
| lev-0 (patches 11–13) | 3 | 0.667 | 0.049 |
| lev-0 (patches 14–15) | 2 | 0.701 | 0.049 |

---

### Epsilon thresholds for Class 1 constrained runs

ε_g is the **only** meaningful axis — featureError is negligible throughout (max 0.049).
ε_f can be set to 0.1 to ensure it never binds.

| ε_g | Accepted | Expected result |
|---|---|---|
| 0.10 | 13 (lev-2 only) | Tight — hard cliff blocks lev-1 |
| 0.40 | 28 (lev-2 + lev-1 + 10 lev-0) | Medium — lev-0 patches 11–15 (globalError 0.667–0.701) blocked |
| 0.75 | 34 (all) | Loose — everything accepted |
