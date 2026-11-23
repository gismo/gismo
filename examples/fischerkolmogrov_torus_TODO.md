# FK torus (CMAME): to do next (written 2026-11-23)

Background and definitions: `fischerkolmogrov_torus_CMAME_MEASUREMENTS.md`.
Server runs: `run_fk_torus_cmame.sh`.

## 1. DSC and RMSE between two solutions (main gap)

The paper's local metrics (Lorenzo et al. 2024) are the Dice similarity
coefficient and the RMSE between N̂ maps. We need them between:
- body-fitted levels (convergence; reference = finest body-fitted run);
- immersed and body-fitted (method comparison).

### Definitions (continuous form of the paper's voxel metrics)
With H_x = H(N_x − 0.15) the tumour indicator (same smoothed H as for V_T,
k = 200, or the sharp step):

    DSC(a,b)  = 2 ∫_Ω H_a H_b dΩ / ( ∫_Ω H_a dΩ + ∫_Ω H_b dΩ )
    RMSE(a,b) = sqrt( (1/V_P) ∫_Ω (N_a − N_b)² dΩ )
    e_L2(a,b) = ‖N_a − N_b‖_{L²(Ω)} / ‖N_b‖_{L²(Ω)}

DSC = 1 means identical tumour regions. N̂ is unit-free, so the RMSE is too.
For reference, the paper reports DSC ≈ 0.8 and RMSE ≈ 0.02 (model vs MRI).

### Implementation plan
All integrals in ONE quadrature loop over the **body-fitted** mesh (exact
geometry, exact V_P).

1. **`--saveAt t1,t2,...`** in both codes: at the listed times (days), write
   the solution with `gsWrite` to `<out>/N_t<time>.xml`.
   - Body-fitted: `mp_nold.patch(0)`, the solution spline in parameter space
     (u,v,w). Note: with the periodic seam its w-knots are unclamped; write it
     as it is.
   - Immersed: `dbasis.basis(0).makeGeometry(Nold)`, a spline over [0,1]³ in
     **physical box coordinates**, because the background map is the identity.
2. **`--compareWith <file>`** and **`--compareType bodyfitted|immersed`** in the
   body-fitted code, plus **`--compareAt <time>`**. At t = compareAt, loop over
   the body-fitted quadrature points u with x = G(u):
   - N_a = own solution at u;
   - N_b = the file's spline:
     - body-fitted file: evaluate at the same u (all body-fitted meshes share
       the geometry map G);
     - immersed file: evaluate at x, i.e. `A.getCoeff(N_im, G)`.
   - Accumulate ∫(N_a − N_b)², ∫N_b², ∫H_aH_b, ∫H_a, ∫H_b and V_P.
   - Write `comparison.csv`: `time,DSC,RMSE,e_L2,V_T_a,V_T_b`.
   - Simplest workflow: rerun the coarse body-fitted case (cheap) with
     `--compareWith` pointing to the saved fine or immersed solution, instead
     of comparing two files offline.
3. Check:
   - DSC(a,a) = 1 and RMSE(a,a) = 0;
   - body-fitted r3 vs r4: DSC → 1 and RMSE → 0 under refinement;
   - immersed r5 vs body-fitted r4 lies in the same range as body-fitted r3
     vs r4.
4. Add the save/compare options to `run_fk_torus_cmame.sh`: `--saveAt 0,365,730`
   for all runs, and a final comparison lane.

Effort: about half a day.

## 2. Consistent initial rate dN/dt(0) (both codes)
Now dNold = 0 at t = 0, which gives a 20% mass-balance error in step 1
(`dintN_dt` against `growth`; from step 2 on they agree). Fix: solve once

    M dN0 = −D K N0 + ∫ ρ N0 (1 − N0) φ dΩ

with the same linear solver, before the time loop. Check that `dintN_dt` ≈
`growth` from step 1.

## 3. Body-fitted biomarkers in one pass
`writeBio` calls `computeInt` 6× per step (each builds its own assembler), and
V_P never changes. That is about 30% of the loop time at r3. Do one loop over
the elements/quadrature points (like the immersed `biomarkers` lambda), and
compute V_P once.

## 4. Same IC projection in both codes
- Body-fitted: L² projection onto the space over Ω.
- Immersed: `gsL2Projection` over the whole background box [0,1]³.

For an identical setup, assemble the immersed IC projection on the cut-cell
quadrature (`blocks`): mass matrix + regularization, rhs ∫ N0 φ over Ω.

## 5. `--plotEvery k`
With Δt = 1, `--plot` writes 731 frames (~2 GB). Write only every k-th frame,
and always the last one.

## 6. After the first server job
- Check immersed r6: peak memory (`sacct -j <id> --format=JobID,MaxRSS,Elapsed`)
  and time per step.
- Check body-fitted r5 with `cgdiag`: time per step, LinIts.
- κ against h from the `cond_*` runs. Table: κ(K) and κ(Jacobi) for both
  methods, plus the immersed active block.
- Update the TeX plots for the new file names (`biomarkers.csv`, t = 0 row,
  `TimeStep` from 0).

## 7. Open questions for the supervisor
- A_p: the paper reports 1/day. Is A_p divided by V_T?
- Which comparison times for DSC/RMSE (e.g. 1 and 2 years)?
- Is Δt = 1 day needed for all levels, or only for the reference?
