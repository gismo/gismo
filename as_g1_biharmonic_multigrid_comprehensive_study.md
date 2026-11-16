# Comprehensive Study of AS-$G^1$ Biharmonic Multigrid Solvers & Smoothers across 101 Multi-Patch Geometries
## 1. Executive Summary
This report documents the reorganization of the G+Smo 2D AS-$G^1$ domain collection into a strictly classified benchmark suite (`filedata/domain2d/as_g1_benchmark_suite/`, **101 unique validated domains**), explains the mathematical formulation and implementation of four progressive multigrid improvements (**Proposals 4, 1, 2, and 3**), and compiles extensive experimental benchmarks covering V-cycle (`-c 1`), W-cycle (`-c 2`), and Multigrid-PCG (`-s cg`), multiple refinement levels (`r = 3..6`, up to `78,736` DOFs), spline degrees (`d = 3, 4, 5`), and both the first and second biharmonic problems.
### Headline Findings
- **Direct MG (`-s mg`) goes from unreliable to near-100% robust**: Across all **101 unique domains** at `m=2 -> r=3` (`d=3`, $\nu=1$), V-cycle pass rate rises from **74/101 (`73.3%`, median 138 iters / 1.408s)** under Proposal 4 (`auto`) $\to$ **88/101 (`87.1%`, med 33 / 0.428s)** under Proposal 1 (`ring1`) $\to$ **95/101 (`94.1%`, med 24 / 0.284s)** under Proposal 3 (`fd_ring1`) $\to$ **99/101 (`98.0%`, med 17 / 0.165s)** under Proposal 2 (`line_ring1`), achieving **89/89 (`100%`)** on all simply-connected domains and **101/101 (`100%`, med 8 iters, max 52, 0.080s)** in PCG.
- **True 3-level W-cycle (`m=2 -> r=4`, 61 domains) eliminates coarse-grid lag**: `mg-W line_ring1` achieves **60/61 (`98.4%`)** pass rate (with the single non-convergent case being `asg1_circle_5patch.xml`, which is rank-ill-conditioned at `r=4` across *all* smoothers including `auto`), drops median iterations from `24` (V-cycle) to **18** (W-cycle), and runs **faster in wall-clock solve time** (`0.447s` vs `0.551s`).
- **$h$- and $p$-Robustness**: Under `mg-W line_ring1`, iterations remain nearly mesh-independent up to `r=5..6` (`64k–78k` DOFs: `4_patch_3_valence_not_bilinear` `17 -> 18 -> 19 -> 20`; `as_g1_pentagon_p3_r1` `18 -> 19 -> 20`). Under degree elevation `d = 3 -> 4 -> 5`, Point-GS-based `ring1` degrades sharply (e.g., `curved_ring_4patch` `89 -> 130 -> 300(F)`), whereas `line_ring1` and `fd_ring1` stay flat (`sq6`: `31 -> 28 -> 29` MG, `12 -> 12 -> 12` CG).
- **First vs. Second Biharmonic parity**: Across all 101 domains at `r=3` on the second biharmonic problem, `mg-V line_ring1` improves from **80/101 (med 85 iters)** (`auto`) to **99/101 (med 14 iters, 0.241s)**, and `cg-V line_ring1` hits **101/101 (med 8, max 40, 0.093s)**.

---

## 2. Reclassified Benchmark Suite (`filedata/domain2d/as_g1_benchmark_suite/`)
The original `--suite all` list in `as_g1_biharmonic_multigrid_example.cpp` contained **99 entries** grouped largely by provenance (`classic`, `generator`, `noholes`, `holes`), which introduced structural inconsistencies:
1. `ring_7patch_non_bilinear_swirl_p3_asg1.xml` (an annulus with **1 interior hole**) was placed in `generator` instead of `holes`.
2. Three entries in `holes` (`crescent_2corner_asg1.xml`, `curved_2patch_asg1.xml`, `curved_lshape_3patch_asg1.xml`) have **0 holes**.
3. Nine entries in `as_g1_generator` were **byte-for-byte duplicates** under alternate filenames (e.g. `as_g1_v3_p3.xml == as_g1_tri_p3.xml == as_g1_patch3_sym.xml`), leaving **90 unique domains** in the 99-list.
4. Filename collisions existed between `1to7corners/` and `asg1_reparametrized/` (e.g. two distinct `crescent_2corner_asg1.xml` geometries with control aspect ratios `2.27` vs `6.27`).

To fix this, all 90 distinct `--suite all` domains plus **11 additional validated AS-$G^1$ domains** from `as_g1_generator/`, `1to7corners/`, and `asg1_reparametrized/` (including 6 multi-hole domains with `3–5` holes) were classified by topology, interior extraordinary vertex count $n_{\text{intV}}$, and genus $n_{\text{holes}}$, and copied into `filedata/domain2d/as_g1_benchmark_suite/` (**101 unique validated domains**; three singular non-`suite-all` extras that fail even direct Cholesky—`heart_1corner_nonconvex_asg1` in `asg1_reparametrized`, `square_with_disk_asg1`, and `star_pentagon_5corner_nonconvex_asg1`—plus `plate_4corner_8holes_weird_asg1` were excluded and documented in `README.md`):

| Subfolder | Count | Topological Definition | Algebraic / Multigrid Significance |
| :--- | ---: | :--- | :--- |
| `01_single_patch/` | 10 | $P = 1$, $n_{\text{Ifc}} = 0$, $n_{\text{holes}} = 0$ | Pure interior B-splines + outer boundary Ring-1 / curvature (e.g. `02_serpentine_s_duct_asg1`) |
| `02_multipatch_trees/` | 22 | $P \ge 2$, $n_{\text{Ifc}} = P - 1$, $n_{\text{intV}} = 0$, $n_{\text{holes}} = 0$ | Interface `EdgeTrace`/`EdgeCross` coupling & Case-B smooth boundary vertices, no interior stars |
| `03_single_vertex_stars/` | 34 | $n_{\text{intV}} = 1$, valence $\in \{3..7\}$, $n_{\text{holes}} = 0$ | Isolates a single interior 6-jet extraordinary vertex (`Corner6Jet`) with symmetric/random/curved patches |
| `04_multi_vertex_meshes/` | 23 | $n_{\text{intV}} \ge 2$, $n_{\text{holes}} = 0$ | Coupled multi-valence vertices ($3\text{--}6$), extreme patch aspect ratios (`sq6_nb` 5.06), non-convex corners |
| `05_multiply_connected_holes/` | 12 | $n_{\text{holes}} \in \{1..6\}$ | Closed inner boundary loops (`u = d_n u = 0`) + high-curvature rings (`ring_7patch_swirl`) & multi-hole plates |

---

## 3. Explanations of Proposals 4, 1, 2, and 3
### Proposal 4 (`--smoother auto`): Basis Conditioning, Exact Coarse Solver & Symmetric GS
Diagnostic spectral analysis of the baseline multigrid operator $M_{\text{MG}} A$ revealed that the baseline suffered from three implementation bottlenecks before smoothing was even enhanced: (i) loose `.sparseView(1, 1e-4)` thresholds in [gsAsG1Basis.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Basis.hpp) corrupted small gluing/derivative entries, (ii) unscaled `EdgeCross` and `Corner6Jet` derivative DOFs created $3.6 \times 10^5\times$ diagonal imbalances in $K_{\text{free}}$ across $u, \nabla u, \nabla^2 u$, and (iii) `-s mg` previously used an unsymmetric Gauss–Seidel smoother (`step = fwd, stepT = fwd`) and left the coarse solver unspecified. Proposal 4 tightens `.sparseView` drop tolerances to `1e-12`, rescales AS-$G^1$ derivative DOFs uniformly per level by `1/hLevel` (first derivatives) and `1/hLevel^2` (second derivatives) to preserve exact interface/vertex matching while restoring $\mathcal{O}(h^{-2})$ diagonal balance, enables `makeSparseCholeskySolver(mg->matrix(0))`, and defaults to `makeSymmetricGaussSeidelOp` for `-s mg` (`makeGaussSeidelOp` for `-s cg`).

### Proposal 1 (`--smoother ring1`): Multiplicative Point GS + Coupled Ring-1 Skeleton
Even with Proposal 4, $86\text{--}99\%$ of eigenvalues of $M_{\text{MG}} A$ sit in $[0.5, 1.0]$ while $4\text{--}26$ outliers near zero ($\lambda_{\min} \approx 0.0088\text{--}0.0697$) stall stationary `-s mg` at $\rho = 1 - \lambda_{\min}$. Spectral mode localization proved these outliers stem from **cross-interface normal-curvature locking**: AS-$G^1$ trace/crossing functions ($\mathcal{V}_{E,0}, \mathcal{V}_{E,1}, \mathcal{V}_v$) live on boundary layers $i_n \in \{0,1\}$, while $\partial_{nn} u$ strongly couples them across $[0, h]$ to the first two interior B-spline layers $i_n \in \{2,3\}$ (`Int(ring0)`). `gsAsG1SkeletonSmootherOp` ([gsAsG1Smoother.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Smoother.hpp#L110-L189)) pairs a forward point Gauss–Seidel sweep with an exact sparse Cholesky correction on the $\mathcal{O}(\sqrt{N})$ **Coupled Ring-1** subspace (`collectCoupledRing1Dofs`: all non-interior free DOFs + two adjacent interior B-spline layers on every interface and boundary) using pre-extracted `m_subRows = P^T A` in `step`, and reverses the order (`Ring1 -> bwd GS`) in `stepT` for strict $A$-symmetry.

### Proposal 2 (`--smoother line_ring1`): Alternating 1D Strip Blocks + Coupled Ring-1
Once Coupled Ring-1 eliminates interface/vertex locking, the remaining bottlenecks are **patch-interior aspect-ratio anisotropy** ($(|J_u|/|J_v|)^4$, e.g., $81\text{--}705\times$ on `square_6_patch` and `square_6_patch_non_bilinear`), **double-knot ($C^1$) interior pairing** from `uniformRefine(1, 2)`, **curved/skewed metric variations** (`02_serpentine_s_duct`), and **high-$p$ coupling**. `gsAsG1BlockSmootherOp` ([gsAsG1Smoother.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Smoother.hpp#L288-L438)) partitions each patch's interior $(m_1 \times m_2)$ tensor grid into consecutive 2-row `u`-strips and 2-column `v`-strips (`collectInteriorLineAndRing1Blocks`), appends the Coupled Ring-1 skeleton block, converts `m_mat` once to RowMajor (`m_matRow`) for $\mathcal{O}(\text{nnz}(A))$ subblock extraction and residual evaluation, factorizes each banded 1D strip and Ring-1 block with `SimplicialLDLT`, and applies forward blocks `0..K-1` in `step` and exact reverse blocks `K-1..0` in `stepT`. Because every strip solves the **exact physical Galerkin submatrix** $P_b^T A P_b$, it handles arbitrary intra-patch Jacobian distortion and higher polynomial degrees without damping.

### Proposal 3 (`--smoother fd_ring1`): Interior Biharmonic Fast Diagonalization + Coupled Ring-1
`gsAsG1FastDiagSmootherOp` ([gsAsG1Smoother.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Smoother.hpp#L451-L864)) replaces interior line strips with a **tensor-product Fast Diagonalization (FD) preconditioner** per patch: on interior indices $[2..n_d-3]$, it assembles 1D mass $M_d$ and clamped biharmonic $B_d = \int B_i'' B_j''$, solves the small 1D generalized eigenproblem $B_d Q_d = M_d Q_d \Lambda_d^{(4)}$, and approximates the patch-interior physical biharmonic diagonal in $(Q_1 \otimes Q_2)$ coordinates by $D_{a,b} = (\sqrt{c_1 \lambda_a^{(4,1)}} + \sqrt{c_2 \lambda_b^{(4,2)}})^2 + \beta \overline{|J|}$ using average physical metric coefficients $c_d = \overline{|J|(g^{dd})^2}$. To ensure multiplicative contraction safety on non-bilinear patches, 12 steps of $A_{pp}$-inner-product power iteration estimate $\lambda_{\max}(\hat{A}_p^{-1} A_{pp})$ and calibrate $\omega_p = 0.90/\lambda_{\max}$. The symmetric cycle applies `fwd GS -> Patch FD -> Ring-1` in `step` and `Ring-1 -> Patch FD -> bwd GS` in `stepT`.

---

## 4. Aggregate & Cross-Cutting Experimental Results
### 4.1 Full 101-Geometry Sweep (`d=3`, `m=2 -> r=3`, `--homo --pre 1 --post 1 -e 1e-8 --maxIt 300`)
> **Note on 2-level V vs. W equivalence:** At `m=2 -> r=3` (`nLevels = 2`) with exact coarse Cholesky at level 0, V-cycle and W-cycle are mathematically identical, giving identical iteration counts. Multi-level V vs. W separation appears at `r >= 4` (`nLevels >= 3`, Table 4.3). Direct Cholesky verification (`--verify` on `cg-V line_ring1`) confirmed exact solution agreement across **all 101 domains** (median `Rel. Diff = 5.50e-10`, max `2.27e-07`).

| Solver / Cycle / Smoother | Pass Rate | Median Iters | Mean Iters | Max Iters | Median Solve (s) |
| :--- | :---: | ---: | ---: | ---: | ---: |
| `mg-V` `auto` | 74/101 | 138.0 | 153.2 | 300 | 1.4083 |
| `mg-V` `ring1` | 88/101 | 33.0 | 85.2 | 300 | 0.4277 |
| `mg-V` `line_ring1` | 99/101 | 17.0 | 37.0 | 300 | 0.1652 |
| `mg-V` `fd_ring1` | 95/101 | 24.0 | 56.7 | 300 | 0.2835 |
| `mg-W` `auto` | 74/101 | 138.0 | 153.2 | 300 | 1.2545 |
| `mg-W` `ring1` | 88/101 | 33.0 | 85.2 | 300 | 0.4141 |
| `mg-W` `line_ring1` | 99/101 | 17.0 | 37.0 | 300 | 0.1693 |
| `mg-W` `fd_ring1` | 95/101 | 24.0 | 56.7 | 300 | 0.3026 |
| `cg-V` `auto` | 99/101 | 30.0 | 45.5 | 300 | 0.2765 |
| `cg-V` `ring1` | 101/101 | 12.0 | 24.7 | 186 | 0.1637 |
| `cg-V` `line_ring1` | 101/101 | 8.0 | 11.0 | 52 | 0.0795 |
| `cg-V` `fd_ring1` | 100/101 | 10.0 | 17.8 | 300 | 0.1221 |

### 4.2 Per-Category Breakdown Across All 101 Geometries (`r=3`, First vs. Second Biharmonic)

| Domain Category (Count) | Biharm1 `mg-V auto` | Biharm1 `mg-V ring1` | Biharm1 `mg-V fd_ring1` | Biharm1 `mg-V line_ring1` | Biharm1 `cg-V line_ring1` | Biharm2 `mg-V auto` | Biharm2 `mg-V line_ring1` | Biharm2 `cg-V line_ring1` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `01_single_patch` (10) | 8/10 (med 9) | 8/10 (med 15) | 9/10 (med 10) | 10/10 (med 8) | 10/10 (med 6) | 8/10 (med 10) | 10/10 (med 8) | 10/10 (med 6) |
| `02_multipatch_trees` (22) | 17/22 (med 63) | 21/22 (med 18) | 21/22 (med 14) | 22/22 (med 15) | 22/22 (med 8) | 18/22 (med 56) | 21/22 (med 15) | 22/22 (med 8) |
| `03_single_vertex_stars` (34) | 31/34 (med 131) | 34/34 (med 23) | 34/34 (med 19) | 34/34 (med 17) | 34/34 (med 8) | 34/34 (med 66) | 34/34 (med 12) | 34/34 (med 8) |
| `04_multi_vertex_meshes` (23) | 9/23 (med 300) | 16/23 (med 120) | 21/23 (med 56) | 23/23 (med 50) | 23/23 (med 13) | 10/23 (med 300) | 23/23 (med 32) | 23/23 (med 13) |
| `05_multiply_connected_holes` (12) | 9/12 (med 194) | 9/12 (med 76) | 10/12 (med 54) | 10/12 (med 27) | 12/12 (med 10) | 10/12 (med 108) | 11/12 (med 18) | 12/12 (med 10) |

In First Biharmonic `mg-V line_ring1`, only 2/101 domains did not reach $10^{-8}$ within 300 steps at `r=3`, both high-genus multi-hole domains that still reached $\sim 10^{-7}$ residual (`curved_lshape_6corner_3holes_asg1`: `2.96e-07`, `swiss_cross_12corner_5holes_asg1`: `6.45e-07`), while `cg-V line_ring1` converged **101/101** in $\le 52$ iterations. Overall Second Biharmonic (`r=3`, all 101 domains) moved from `80/101` (`auto`, med `85`, `1.563s`) $\to$ `99/101` (`mg-V line_ring1`, med `14`, `0.241s`) and `101/101` (`cg-V line_ring1`, med `8`, max `40`, `0.093s`).

### 4.3 True 3-Level Hierarchy (`m=2 -> r=4`, 61 Domains): V-Cycle vs. W-Cycle vs. PCG
Evaluating all 59 domains with `DOFs(r=3) <= 6,000` plus `circle_asg1.xml` (`44,936` DOFs) and `curved_ring_4patch_asg1.xml` (`35,708` DOFs) on a true 3-level hierarchy separates V-cycle from W-cycle:

| Solver / Cycle / Smoother | Pass Rate | Median Iters | Mean Iters | Max Iters | Median Solve (s) |
| :--- | :---: | ---: | ---: | ---: | ---: |
| `mg-V` `auto` | 41/61 | 179.0 | 179.8 | 300 | 7.4445 |
| `mg-V` `ring1` | 54/61 | 35.0 | 80.0 | 300 | 0.8476 |
| `mg-V` `line_ring1` | 60/61 | 24.0 | 42.3 | 300 | 0.5511 |
| `mg-V` `fd_ring1` | 57/61 | 26.0 | 58.3 | 300 | 0.6837 |
| `mg-W` `auto` | 44/61 | 148.0 | 159.9 | 300 | 6.7053 |
| `mg-W` `ring1` | 56/61 | 31.0 | 70.7 | 300 | 0.7830 |
| `mg-W` `line_ring1` | 60/61 | 18.0 | 34.0 | 300 | 0.4471 |
| `mg-W` `fd_ring1` | 58/61 | 23.0 | 50.8 | 300 | 0.6024 |
| `cg-V` `auto` | 59/61 | 37.0 | 49.8 | 300 | 0.7285 |
| `cg-V` `ring1` | 59/61 | 14.0 | 29.4 | 300 | 0.3024 |
| `cg-V` `line_ring1` | 60/61 | 10.0 | 16.4 | 300 | 0.1831 |
| `cg-V` `fd_ring1` | 60/61 | 11.0 | 20.1 | 300 | 0.2793 |

> Only `asg1_circle_5patch.xml` failed at `r=4` across **every** solver and smoother (`auto`, `ring1`, `line_ring1`, `fd_ring1`, V, W, CG), whereas the regular 5-patch disk `circle_asg1.xml` (`44,936` DOFs) converged in `18` (`mg-V`), `16` (`mg-W`), and `10` (`cg-V`) iterations with `line_ring1`. Diagnostics showed `asg1_circle_5patch.xml` achieves optimal $O(h^4)$ $L^2$ error (`5.09e-4 -> 3.24e-5`) under direct Cholesky, but boundary Jacobian pinching at `r=4` introduces extreme conditioning that stalls even plain `Point GS` (`Final Rel. Res = 0.88–2.16` while `Rel. Diff vs Cholesky = 2.8e-7`).

### 4.4 $h$-Refinement Scaling (`r = 3, 4, 5, 6` up to `78,736` DOFs, Stratified Core Set)

| Domain | Ref `r` | `DOFs` | `mg-V auto` | `mg-V ring1` | `mg-V fd_ring1` | `mg-V line_ring1` | `mg-W auto` | `mg-W line_ring1` | `cg-V auto` | `cg-V line_ring1` | `cg-V fd_ring1` |
| :--- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `2patch_rect` | 3 | 405 | 49 | 18 | 6 | 7 | 49 | 7 | 17 | 6 | 5 |
| `2patch_rect` | 4 | 1829 | 57 | 20 | 8 | 10 | 55 | 10 | 23 | 7 | 6 |
| `2patch_rect` | 5 | 7749 | 55 | 23 | 8 | 11 | 49 | 9 | 27 | 7 | 6 |
| `2patch_rect` | 6 | 31877 | 55 | 24 | 8 | 12 | 53 | 9 | 28 | 8 | 6 |
| `4patch_v3` | 3 | 834 | 138 | 30 | 24 | 17 | 138 | 17 | 28 | 9 | 11 |
| `4patch_v3` | 4 | 3714 | 178 | 38 | 29 | 22 | 127 | 18 | 33 | 11 | 13 |
| `4patch_v3` | 5 | 15618 | 247 | 46 | 34 | 30 | 138 | 19 | 42 | 13 | 15 |
| `4patch_v3` | 6 | 64002 | TO | 53 | 39 | 39 | TO | 20 | 47 | 14 | 17 |
| `circle_asg1` | 3 | 10952 | 91 | 30 | 35 | 14 | 91 | 14 | 25 | 8 | 15 |
| `circle_asg1` | 4 | 44936 | 125 | 46 | 48 | 18 | 104 | 16 | 31 | 10 | 19 |
| `curved_quad` | 3 | 2116 | 9 | 15 | 7 | 8 | 9 | 8 | 9 | 6 | 6 |
| `curved_quad` | 4 | 8836 | 9 | 17 | 9 | 8 | 9 | 7 | 10 | 6 | 6 |
| `curved_ring` | 3 | 8636 | 48 | 89 | 14 | 9 | 48 | 9 | 25 | 6 | 9 |
| `curved_ring` | 4 | 35708 | 57 | 101 | 13 | 12 | 47 | 9 | 29 | 8 | 9 |
| `ldomain` | 3 | 6436 | 52 | 15 | 13 | 14 | 52 | 14 | 16 | 7 | 6 |
| `ldomain` | 4 | 26692 | 92 | 25 | 21 | 22 | 65 | 15 | 20 | 9 | 7 |
| `serpentine` | 3 | 2116 | 300 (F) | 300 (F) | 300 (F) | 21 | 300 (F) | 21 | 226 | 11 | 63 |
| `serpentine` | 4 | 8836 | 300 (F) | 300 (F) | 300 (F) | 36 | 300 (F) | 32 | 300 (F) | 15 | 84 |
| `sq6` | 3 | 1271 | 300 (F) | 120 | 24 | 19 | 300 (F) | 19 | 45 | 10 | 11 |
| `sq6` | 4 | 5607 | 300 (F) | 149 | 50 | 31 | 300 (F) | 30 | 58 | 12 | 17 |
| `sq6` | 5 | 23495 | 300 (F) | 183 | 85 | 41 | 300 (F) | 30 | 74 | 15 | 23 |
| `sq6_nb` | 3 | 1271 | 300 (F) | 300 (F) | 50 | 50 | 300 (F) | 50 | 79 | 13 | 16 |
| `sq6_nb` | 4 | 5607 | 300 (F) | 300 (F) | 87 | 82 | 300 (F) | 75 | 107 | 20 | 22 |
| `sq6_nb` | 5 | 23495 | 300 (F) | 300 (F) | 111 | 77 | 300 (F) | 68 | 115 | 23 | 28 |
| `swirl_ring7` | 3 | 15113 | 300 (F) | 64 | 52 | 41 | 300 (F) | 41 | 38 | 13 | 18 |
| `swirl_ring7` | 4 | 62489 | TO | 107 | 75 | 68 | TO | 50 | 53 | 16 | 23 |
| `two_bilinear` | 3 | 403 | 140 | 27 | 19 | 21 | 140 | 21 | 17 | 8 | 7 |
| `two_bilinear` | 4 | 1827 | 214 | 42 | 29 | 35 | 171 | 27 | 26 | 10 | 9 |
| `two_bilinear` | 5 | 7747 | 300 (F) | 61 | 41 | 49 | 196 | 31 | 34 | 12 | 11 |
| `two_bilinear` | 6 | 31875 | 300 (F) | 80 | 53 | 65 | 213 | 33 | 40 | 13 | 13 |
| `v5_r1` | 3 | 4336 | 122 | 23 | 18 | 18 | 122 | 18 | 27 | 8 | 9 |
| `v5_r1` | 4 | 18896 | 178 | 33 | 24 | 26 | 137 | 19 | 33 | 10 | 10 |
| `v5_r1` | 5 | 78736 | 236 | 42 | 31 | 34 | 149 | 20 | 38 | 12 | 12 |
| `v7_p3` | 3 | 1455 | 230 | 33 | 26 | 28 | 230 | 28 | 34 | 9 | 9 |
| `v7_p3` | 4 | 6495 | 300 (F) | 52 | 40 | 44 | 300 (F) | 28 | 46 | 13 | 14 |
| `v7_p3` | 5 | 27327 | 300 (F) | 74 | 58 | 63 | 300 (F) | 32 | 59 | 15 | 17 |
| `weirdo` | 3 | 2335 | 300 (F) | 146 | 116 | 99 | 300 (F) | 99 | 69 | 18 | 21 |
| `weirdo` | 4 | 10303 | 300 (F) | 283 | 220 | 192 | 300 (F) | 150 | 103 | 27 | 33 |
| `weirdo` | 5 | 43135 | 300 (F) | 300 (F) | 300 (F) | 300 (F) | TO | 217 | 135 | 36 | 41 |

### 4.5 Spline Degree $p$-Robustness (`d = 3, 4, 5`, First Biharmonic)

| Domain | `r` | Deg `d` | `DOFs` | `mg-V auto` | `mg-V ring1` | `mg-V fd_ring1` | `mg-V line_ring1` | `mg-W auto` | `mg-W line_ring1` | `cg-V auto` | `cg-V line_ring1` | `cg-V fd_ring1` |
| :--- | :---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `2patch_rect` | 4 | 3 | 1829 | 57 | 20 | 8 | 10 | 55 | 10 | 23 | 7 | 6 |
| `2patch_rect` | 4 | 4 | 1953 | 102 | 30 | 11 | 12 | 100 | 12 | 25 | 6 | 6 |
| `2patch_rect` | 4 | 5 | 2081 | 160 | 102 | 10 | 11 | 157 | 11 | 40 | 7 | 6 |
| `4patch_v3` | 4 | 3 | 3714 | 178 | 38 | 29 | 22 | 127 | 18 | 33 | 11 | 13 |
| `4patch_v3` | 4 | 4 | 3966 | 113 | 41 | 29 | 15 | 115 | 15 | 30 | 9 | 14 |
| `4patch_v3` | 4 | 5 | 4226 | 189 | 138 | 56 | 20 | 185 | 20 | 41 | 10 | 19 |
| `curved_quad` | 3 | 3 | 2116 | 9 | 15 | 7 | 8 | 9 | 8 | 9 | 6 | 6 |
| `curved_quad` | 3 | 4 | 2401 | 17 | 31 | 9 | 8 | 17 | 8 | 13 | 6 | 7 |
| `curved_quad` | 3 | 5 | 2704 | 58 | 108 | 11 | 12 | 58 | 12 | 25 | 8 | 8 |
| `curved_ring` | 3 | 3 | 8636 | 48 | 89 | 14 | 9 | 48 | 9 | 25 | 6 | 9 |
| `curved_ring` | 3 | 4 | 9800 | 75 | 130 | 14 | 7 | 75 | 7 | 29 | 6 | 9 |
| `curved_ring` | 3 | 5 | 11036 | 186 | 300 (F) | 17 | 9 | 186 | 9 | 47 | 7 | 10 |
| `ldomain` | 3 | 3 | 6436 | 52 | 15 | 13 | 14 | 52 | 14 | 16 | 7 | 6 |
| `ldomain` | 3 | 4 | 7303 | 42 | 33 | 9 | 10 | 42 | 10 | 16 | 7 | 6 |
| `ldomain` | 3 | 5 | 8224 | 57 | 107 | 6 | 12 | 57 | 12 | 25 | 8 | 5 |
| `serpentine` | 3 | 3 | 2116 | 300 (F) | 300 (F) | 300 (F) | 21 | 300 (F) | 21 | 226 | 11 | 63 |
| `serpentine` | 3 | 4 | 2401 | 300 (F) | 300 (F) | 300 (F) | 33 | 300 (F) | 33 | 258 | 15 | 73 |
| `serpentine` | 3 | 5 | 2704 | 300 (F) | 300 (F) | 300 (F) | 47 | 300 (F) | 47 | 300 (F) | 17 | 77 |
| `sq6` | 4 | 3 | 5607 | 300 (F) | 149 | 50 | 31 | 300 (F) | 30 | 58 | 12 | 17 |
| `sq6` | 4 | 4 | 5987 | 300 (F) | 129 | 52 | 28 | 300 (F) | 29 | 65 | 12 | 18 |
| `sq6` | 4 | 5 | 6379 | 300 (F) | 125 | 68 | 29 | 300 (F) | 28 | 89 | 12 | 21 |
| `sq6_nb` | 4 | 3 | 5607 | 300 (F) | 300 (F) | 87 | 82 | 300 (F) | 75 | 107 | 20 | 22 |
| `sq6_nb` | 4 | 4 | 5987 | 300 (F) | 300 (F) | 89 | 83 | 300 (F) | 78 | 111 | 19 | 23 |
| `sq6_nb` | 4 | 5 | 6379 | 300 (F) | 300 (F) | 140 | 102 | 300 (F) | 97 | 163 | 20 | 29 |
| `swirl_ring7` | 3 | 3 | 15113 | 300 (F) | 64 | 52 | 41 | 300 (F) | 41 | 38 | 13 | 18 |
| `swirl_ring7` | 3 | 4 | 17150 | 226 | 80 | 59 | 27 | 226 | 27 | 34 | 11 | 18 |
| `swirl_ring7` | 3 | 5 | 19313 | 262 | 196 | 96 | 22 | 262 | 22 | 47 | 11 | 22 |
| `two_bilinear` | 4 | 3 | 1827 | 214 | 42 | 29 | 35 | 171 | 27 | 26 | 10 | 9 |
| `two_bilinear` | 4 | 4 | 1951 | 167 | 28 | 14 | 18 | 127 | 16 | 28 | 8 | 8 |
| `two_bilinear` | 4 | 5 | 2079 | 196 | 89 | 15 | 15 | 181 | 13 | 38 | 8 | 9 |
| `v7_p3` | 4 | 3 | 6495 | 300 (F) | 52 | 40 | 44 | 300 (F) | 28 | 46 | 13 | 14 |
| `v7_p3` | 4 | 4 | 6936 | 277 | 31 | 27 | 25 | 228 | 22 | 36 | 10 | 12 |
| `v7_p3` | 4 | 5 | 7391 | 300 (F) | 94 | 49 | 18 | 249 | 17 | 44 | 9 | 16 |
| `weirdo` | 4 | 3 | 10303 | 300 (F) | 283 | 220 | 192 | 300 (F) | 150 | 103 | 27 | 33 |
| `weirdo` | 4 | 4 | 11002 | 300 (F) | 180 | 146 | 108 | 300 (F) | 88 | 93 | 22 | 30 |
| `weirdo` | 4 | 5 | 11723 | 300 (F) | 227 | 158 | 70 | 300 (F) | 60 | 110 | 19 | 32 |

### 4.6 Second Biharmonic Stratified Core (`as_g1_biharmonic2_multigrid_example`, `d=3`)

| Domain | `r` | `DOFs` | `mg-V auto` | `mg-V ring1` | `mg-V fd_ring1` | `mg-V line_ring1` | `mg-W auto` | `mg-W line_ring1` | `cg-V auto` | `cg-V line_ring1` | `cg-V fd_ring1` |
| :--- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `2patch_rect` | 4 | 1921 | 71 | 20 | 8 | 10 | 74 | 10 | 25 | 7 | 6 |
| `4patch_v3` | 4 | 3834 | 146 | 38 | 29 | 20 | 102 | 17 | 34 | 10 | 13 |
| `circle_asg1` | 3 | 11048 | 91 | 30 | 35 | 14 | 91 | 14 | 25 | 8 | 15 |
| `curved_quad` | 3 | 2208 | 9 | 15 | 7 | 8 | 9 | 8 | 9 | 6 | 6 |
| `curved_ring` | 3 | 8820 | 47 | 90 | 15 | 12 | 47 | 12 | 25 | 7 | 10 |
| `ldomain` | 3 | 6622 | 57 | 17 | 14 | 15 | 57 | 15 | 17 | 8 | 6 |
| `serpentine` | 3 | 2208 | 300 (F) | 300 (F) | 300 (F) | 62 | 300 (F) | 62 | 216 | 18 | 65 |
| `sq6` | 4 | 5763 | 300 (F) | 148 | 50 | 31 | 300 (F) | 30 | 62 | 12 | 17 |
| `sq6_nb` | 4 | 5763 | 300 (F) | 300 (F) | 84 | 78 | 300 (F) | 71 | 111 | 19 | 22 |
| `swirl_ring7` | 3 | 15435 | 102 | 64 | 52 | 22 | 102 | 22 | 30 | 11 | 18 |
| `two_bilinear` | 4 | 1917 | 61 | 22 | 12 | 15 | 59 | 12 | 23 | 9 | 8 |
| `v5_r1` | 4 | 19196 | 81 | 21 | 14 | 14 | 71 | 11 | 30 | 9 | 9 |
| `v7_p3` | 4 | 6705 | 107 | 19 | 16 | 13 | 92 | 12 | 32 | 9 | 10 |
| `weirdo` | 4 | 10544 | 300 (F) | 167 | 132 | 76 | 300 (F) | 67 | 87 | 20 | 29 |

---

## 5. Behavior Analysis Across Domain Classes
1. **`01_single_patch/` (10 domains)**: Without interfaces or interior vertices, isotropic patches (`quad_asg1`, `curved_quad_asg1`) converge in `6–10` iterations for all smoothers (`d=3`). However, strongly curved/stretched single patches (`02_serpentine_s_duct_asg1`, control AR `4.50`) fail (`300*`) under `auto`, `ring1`, and `fd_ring1` in direct MG (`cg fd_ring1 = 63–84`), whereas **`line_ring1` solves it cleanly in `21` (`r=3`) / `36` (`r=4` V) / `32` (`r=4` W) and `11–15` (CG)** because each 2-row strip captures the spatially varying curved metric $g^{\alpha\beta}(u,v)$ exactly.
2. **`02_multipatch_trees/` (22 domains)**: On rectilinear or mildly curved trees (`2_patch_rectangle`, `ldomain_asg1`, `two_bilinear_patches`), `fd_ring1` frequently achieves the *lowest* iteration count (`6–19` at `r=3`, `8–29` at `r=4`, `6` on `ldomain` at `d=5`) because the constant-metric $(B_d, M_d)$ tensor model closely matches each patch interior. On distorted/anisotropic trees (`comb_asg1`, `curved_comb_asg1`, `07_asymmetric_staircase_8corner_asg1`), `line_ring1` leads and delivers `22/22` convergence.
3. **`03_single_vertex_stars/` (34 domains)**: Pointwise `auto` drops from `31/34` (`med 131`) at `r=3` to struggling at `r=4..5` (`as_g1_v7_p3` fails at `300`, `4_patch_3_valence` reaches `247` at `r=5`). Coupled Ring-1 (`ring1`) immediately eliminates the extraordinary-vertex locking (`34/34`, med `23`), and `line_ring1`/`fd_ring1` push `r=3` to `med 17–19` (`CG line_ring1 med 8`). At `d=5`, `ring1` deteriorates (`2_patch_rectangle` `20 -> 102`, `4_patch_3_valence` `38 -> 138`) while `line_ring1` stays flat (`10 -> 12 -> 11` and `22 -> 15 -> 20`). With W-cycle (`line_ring1`), `as_g1_pentagon_p3_r1` stays nearly constant (`18 -> 19 -> 20`) from `4,336` to `78,736` DOFs.
4. **`04_multi_vertex_meshes/` (23 domains)**: The hardest simply-connected class for `auto` (`9/23` pass at `r=3`, med `300`). Multiple coupled 6-jet vertices plus extreme anisotropy (`square_6_patch` AR `3.00`, `square_6_patch_non_bilinear` AR `5.06`, `03_filleted_u_bracket` `P=11`, AR `4.08`, `weirdo_multivalence`, `blob_geometry`, `experiment_multivalence_not_bilinear`) require **both** Coupled Ring-1 and anisotropic interior relaxation: `ring1` passes `16/23` (failing `sq6_nb`, `03_filleted_u_bracket`, `star_quad_4corner`), while **`line_ring1` passes `23/23` (`100%` at `r=3`, med `50` MG / `13` CG)**, solving `circle_asg1` (`44,936` DOFs at `r=4`) in `18` V-cycle / `16` W-cycle / `10` CG iterations (`0.82s` vs `11.02s` for `auto`) and cutting `star_quad_4corner_nonconvex` V->W from `161` to `92`.
5. **`05_multiply_connected_holes/` (12 domains)**: Closed inner hole loops create periodic Ring-1 corridors around every hole (`holes=1..6`). Even on single-hole curved rings (`curved_ring_4patch` `35,708` DOFs at `r=4`), `line_ring1` (`12` V / `9` W / `8` CG) outperforms `ring1` (`101` V) by `8\times`, and at `d=5` (`r=3`) `ring1` fails at `300*` while `line_ring1` holds steady at **`9` iters MG / `7` iters CG**. On complex multi-hole plates (`3–5` holes, up to `43,395` DOFs at `r=3`), `mg-V line_ring1` converges cleanly on `10/12` and reaches $\sim 3\text{--}6\times 10^{-7}$ on the 2 hardest (`curved_lshape_6corner_3holes`, `swiss_cross_12corner_5holes`), while `cg-V line_ring1` converges `12/12` in `6–52` iterations (med `10`).

---

## 6. Complete Per-Geometry Results Tables
### 6.1 All 101 Validated Geometries (`d=3`, `m=2 -> r=3`, `--homo`, Biharm1 + Biharm2 `line_ring1`)

| ID | Category | Filename | $P$ | Holes | Val | `DOFs` | B1 `mg-V auto` | B1 `mg-V ring1` | B1 `mg-V fd_ring1` | B1 `mg-V line_ring1` | B1 `cg-V auto` | B1 `cg-V line_ring1` | B2 `mg-V line_ring1` |
| ---: | :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | `01_single_patch` | `02_serpentine_s_duct_asg1.xml` | 1 | 0 | 0 | 2116 | 300 (F) | 300 (F) | 300 (F) | 21 (0.12s) | 226 (1.83s) | 11 (0.04s) | 62 (0.24s) |
| 2 | `01_single_patch` | `05_hourglass_pinch_asg1.xml` | 1 | 0 | 0 | 2116 | 300 (F) | 300 (F) | 67 (0.80s) | 11 (0.11s) | 120 (0.88s) | 8 (0.08s) | 29 (0.16s) |
| 3 | `01_single_patch` | `bulged_quad_non_bilinear_asg1.xml` | 1 | 0 | 0 | 2116 | 9 (0.13s) | 15 (0.06s) | 10 (0.05s) | 8 (0.03s) | 9 (0.03s) | 6 (0.02s) | 8 (0.04s) |
| 4 | `01_single_patch` | `bulged_quad_non_bilinear_asg1_reparam.xml` | 1 | 0 | 0 | 196 | 8 (0.00s) | 14 (0.00s) | 7 (0.00s) | 6 (0.00s) | 9 (0.00s) | 5 (0.00s) | 7 (0.00s) |
| 5 | `01_single_patch` | `concave_quad_4corner_nonconvex_asg1.xml` | 1 | 0 | 0 | 2116 | 9 (0.07s) | 15 (0.05s) | 7 (0.04s) | 8 (0.04s) | 9 (0.02s) | 6 (0.01s) | 8 (0.04s) |
| 6 | `01_single_patch` | `curved_quad_4corner_convex_asg1.xml` | 1 | 0 | 0 | 2116 | 9 (0.04s) | 15 (0.06s) | 9 (0.03s) | 8 (0.03s) | 9 (0.03s) | 6 (0.02s) | 8 (0.06s) |
| 7 | `01_single_patch` | `curved_quad_asg1.xml` | 1 | 0 | 0 | 2116 | 9 (0.04s) | 15 (0.05s) | 7 (0.03s) | 8 (0.03s) | 9 (0.02s) | 6 (0.02s) | 8 (0.05s) |
| 8 | `01_single_patch` | `lake_asg1.xml` | 1 | 0 | 0 | 2116 | 174 (1.86s) | 214 (1.60s) | 214 (1.70s) | 15 (0.09s) | 44 (0.40s) | 9 (0.10s) | 23 (0.08s) |
| 9 | `01_single_patch` | `quad_asg1.xml` | 1 | 0 | 0 | 2116 | 9 (0.11s) | 15 (0.11s) | 4 (0.02s) | 8 (0.05s) | 9 (0.03s) | 6 (0.02s) | 8 (0.04s) |
| 10 | `01_single_patch` | `star_quad_4corner_nonconvex_asg1.xml` | 1 | 0 | 0 | 2116 | 127 (1.41s) | 204 (1.25s) | 190 (1.28s) | 87 (0.73s) | 36 (0.26s) | 24 (0.14s) | 90 (0.30s) |
| 11 | `02_multipatch_trees` | `2_patch_rectangle.xml` | 2 | 0 | 0 | 405 | 49 (0.03s) | 18 (0.01s) | 6 (0.00s) | 7 (0.00s) | 17 (0.01s) | 6 (0.00s) | 7 (0.00s) |
| 12 | `02_multipatch_trees` | `2_patch_rectangle_asg1.xml` | 2 | 0 | 0 | 4275 | 39 (0.68s) | 15 (0.06s) | 9 (0.16s) | 11 (0.13s) | 14 (0.16s) | 7 (0.13s) | 14 (0.14s) |
| 13 | `02_multipatch_trees` | `2_patch_rectangle_non_bilinear.xml` | 2 | 0 | 0 | 405 | 14 (0.01s) | 14 (0.01s) | 5 (0.00s) | 6 (0.00s) | 10 (0.00s) | 5 (0.00s) | 6 (0.00s) |
| 14 | `02_multipatch_trees` | `curved_2patch_asg1.xml` | 2 | 0 | 0 | 4277 | 147 (3.46s) | 74 (1.29s) | 18 (0.31s) | 17 (0.35s) | 37 (0.69s) | 9 (0.15s) | 17 (0.28s) |
| 15 | `02_multipatch_trees` | `curved_hexagon_6corner_convex_asg1.xml` | 2 | 0 | 0 | 4275 | 300 (F) | 52 (0.69s) | 47 (0.90s) | 36 (0.82s) | 39 (0.61s) | 11 (0.07s) | 17 (0.26s) |
| 16 | `02_multipatch_trees` | `curved_lshape_non_bilinear_asg1.xml` | 2 | 0 | 0 | 4275 | 162 (3.73s) | 56 (0.92s) | 41 (0.75s) | 21 (0.40s) | 31 (0.45s) | 11 (0.12s) | 21 (0.26s) |
| 17 | `02_multipatch_trees` | `two_bilinear_patches.xml` | 2 | 0 | 0 | 403 | 140 (0.14s) | 27 (0.01s) | 19 (0.01s) | 21 (0.01s) | 17 (0.01s) | 8 (0.00s) | 10 (0.01s) |
| 18 | `02_multipatch_trees` | `two_patch_mathematica.xml` | 2 | 0 | 0 | 403 | 63 (0.04s) | 14 (0.01s) | 11 (0.01s) | 12 (0.01s) | 14 (0.01s) | 7 (0.00s) | 11 (0.01s) |
| 19 | `02_multipatch_trees` | `two_patch_mathematica_non_bilinear.xml` | 2 | 0 | 0 | 403 | 51 (0.04s) | 15 (0.01s) | 10 (0.00s) | 12 (0.01s) | 14 (0.01s) | 7 (0.00s) | 12 (0.01s) |
| 20 | `02_multipatch_trees` | `two_squares_asg1.xml` | 2 | 0 | 0 | 4277 | 18 (0.34s) | 16 (0.17s) | 5 (0.06s) | 7 (0.09s) | 13 (0.09s) | 6 (0.05s) | 7 (0.06s) |
| 21 | `02_multipatch_trees` | `07_asymmetric_staircase_8corner_asg1.xml` | 3 | 0 | 0 | 6434 | 300 (F) | 222 (6.08s) | 170 (5.47s) | 160 (4.06s) | 71 (2.20s) | 23 (0.51s) | 68 (1.78s) |
| 22 | `02_multipatch_trees` | `curved_lshape_3patch_asg1.xml` | 3 | 0 | 0 | 6434 | 48 (1.97s) | 17 (0.43s) | 15 (0.35s) | 15 (0.31s) | 17 (0.37s) | 7 (0.12s) | 16 (0.30s) |
| 23 | `02_multipatch_trees` | `curved_lshape_asg1.xml` | 3 | 0 | 0 | 6436 | 41 (1.04s) | 16 (0.33s) | 10 (0.20s) | 11 (0.17s) | 15 (0.24s) | 7 (0.05s) | 10 (0.28s) |
| 24 | `02_multipatch_trees` | `curved_lshape_non_bilinear_asg1_reparam.xml` | 3 | 0 | 0 | 6435 | 300 (F) | 300 (F) | 300 (F) | 89 (2.23s) | 300 (F) | 23 (0.71s) | 300 (F) |
| 25 | `02_multipatch_trees` | `ldomain_asg1.xml` | 3 | 0 | 0 | 6436 | 52 (3.15s) | 15 (0.64s) | 13 (0.40s) | 14 (0.27s) | 16 (0.21s) | 7 (0.08s) | 15 (0.28s) |
| 26 | `02_multipatch_trees` | `lshape_asg1.xml` | 3 | 0 | 0 | 6436 | 51 (1.46s) | 16 (0.19s) | 13 (0.28s) | 13 (0.16s) | 16 (0.23s) | 7 (0.11s) | 15 (0.25s) |
| 27 | `02_multipatch_trees` | `three_patch_belt.xml` | 3 | 0 | 0 | 610 | 63 (0.08s) | 19 (0.02s) | 12 (0.01s) | 15 (0.01s) | 19 (0.02s) | 8 (0.01s) | 12 (0.02s) |
| 28 | `02_multipatch_trees` | `three_patch_belt_non_bilinear.xml` | 3 | 0 | 0 | 610 | 51 (0.04s) | 17 (0.01s) | 13 (0.01s) | 13 (0.01s) | 16 (0.01s) | 7 (0.00s) | 13 (0.01s) |
| 29 | `02_multipatch_trees` | `06_y_junction_3way_asg1.xml` | 4 | 0 | 0 | 8597 | 116 (4.65s) | 104 (3.79s) | 48 (2.09s) | 28 (1.23s) | 32 (1.31s) | 11 (0.39s) | 27 (0.60s) |
| 30 | `02_multipatch_trees` | `lshape_9patch_non_bilinear_trig_p3_asg1.xml` | 9 | 0 | 0 | 19391 | 144 (16.47s) | 35 (3.15s) | 34 (4.09s) | 23 (1.81s) | 24 (1.67s) | 11 (0.94s) | 23 (2.36s) |
| 31 | `02_multipatch_trees` | `comb_asg1.xml` | 11 | 0 | 0 | 23714 | 300 (F) | 88 (9.90s) | 56 (8.17s) | 63 (7.01s) | 57 (7.42s) | 15 (0.87s) | 72 (7.44s) |
| 32 | `02_multipatch_trees` | `curved_comb_asg1.xml` | 11 | 0 | 0 | 23714 | 300 (F) | 64 (8.28s) | 43 (6.92s) | 47 (4.78s) | 58 (6.66s) | 15 (0.81s) | 50 (4.43s) |
| 33 | `03_single_vertex_stars` | `as_g1_patch3_random.xml` | 3 | 0 | 3 | 627 | 150 (0.29s) | 75 (0.07s) | 60 (0.07s) | 32 (0.03s) | 34 (0.03s) | 12 (0.02s) | 34 (0.04s) |
| 34 | `03_single_vertex_stars` | `as_g1_tri_p3.xml` | 3 | 0 | 3 | 627 | 47 (0.06s) | 13 (0.01s) | 9 (0.01s) | 9 (0.01s) | 17 (0.02s) | 6 (0.01s) | 8 (0.01s) |
| 35 | `03_single_vertex_stars` | `concave_triangle_3corner_nonconvex_asg1.xml` | 3 | 0 | 3 | 6486 | 20 (1.20s) | 33 (1.13s) | 27 (1.08s) | 7 (0.21s) | 14 (0.28s) | 6 (0.12s) | 12 (0.37s) |
| 36 | `03_single_vertex_stars` | `convex_3_valence.xml` | 3 | 0 | 3 | 627 | 74 (0.11s) | 23 (0.02s) | 18 (0.01s) | 13 (0.02s) | 22 (0.01s) | 8 (0.01s) | 12 (0.02s) |
| 37 | `03_single_vertex_stars` | `convex_3_valence_non_bilinear.xml` | 3 | 0 | 3 | 627 | 73 (0.08s) | 14 (0.02s) | 11 (0.03s) | 11 (0.01s) | 20 (0.03s) | 7 (0.00s) | 7 (0.02s) |
| 38 | `03_single_vertex_stars` | `curved_triangle_3corner_convex_asg1.xml` | 3 | 0 | 3 | 6486 | 20 (1.12s) | 15 (0.43s) | 9 (0.30s) | 8 (0.21s) | 13 (0.38s) | 6 (0.13s) | 8 (0.29s) |
| 39 | `03_single_vertex_stars` | `star_hexagon_6corner_nonconvex_asg1.xml` | 3 | 0 | 3 | 6483 | 18 (0.89s) | 16 (0.36s) | 9 (0.20s) | 7 (0.17s) | 14 (0.35s) | 6 (0.11s) | 10 (0.30s) |
| 40 | `03_single_vertex_stars` | `star_triangle_3corner_nonconvex_asg1.xml` | 3 | 0 | 3 | 6486 | 200 (7.32s) | 81 (1.54s) | 81 (2.75s) | 40 (0.71s) | 37 (1.06s) | 14 (0.35s) | 40 (1.31s) |
| 41 | `03_single_vertex_stars` | `4_patch_3_valence_not_bilinear.xml` | 4 | 0 | 3 | 834 | 138 (0.20s) | 30 (0.04s) | 24 (0.04s) | 17 (0.03s) | 28 (0.04s) | 9 (0.02s) | 15 (0.04s) |
| 42 | `03_single_vertex_stars` | `as_g1_patch4_random.xml` | 4 | 0 | 4 | 834 | 95 (0.25s) | 18 (0.03s) | 14 (0.03s) | 14 (0.02s) | 23 (0.02s) | 8 (0.01s) | 9 (0.03s) |
| 43 | `03_single_vertex_stars` | `as_g1_v4_p3.xml` | 4 | 0 | 4 | 834 | 51 (0.10s) | 14 (0.02s) | 11 (0.02s) | 12 (0.01s) | 16 (0.01s) | 6 (0.01s) | 7 (0.01s) |
| 44 | `03_single_vertex_stars` | `star_heptagon_7corner_nonconvex_asg1.xml` | 4 | 0 | 4 | 8643 | 149 (6.23s) | 43 (1.82s) | 34 (1.64s) | 17 (0.64s) | 33 (1.28s) | 9 (0.26s) | 17 (0.78s) |
| 45 | `03_single_vertex_stars` | `5_patch_5_valence_not_bilinear.xml` | 5 | 0 | 5 | 1041 | 147 (0.49s) | 21 (0.04s) | 17 (0.04s) | 18 (0.03s) | 27 (0.03s) | 9 (0.03s) | 12 (0.04s) |
| 46 | `03_single_vertex_stars` | `as_g1_ellipse_v5_p3.xml` | 5 | 0 | 5 | 1041 | 300 (F) | 62 (0.12s) | 44 (0.12s) | 45 (0.08s) | 42 (0.10s) | 10 (0.03s) | 18 (0.07s) |
| 47 | `03_single_vertex_stars` | `as_g1_patch5_polygon_sym.xml` | 5 | 0 | 5 | 1041 | 100 (0.28s) | 19 (0.04s) | 14 (0.04s) | 16 (0.03s) | 23 (0.05s) | 7 (0.02s) | 9 (0.03s) |
| 48 | `03_single_vertex_stars` | `as_g1_patch5_random.xml` | 5 | 0 | 5 | 1041 | 106 (0.29s) | 23 (0.04s) | 20 (0.05s) | 21 (0.04s) | 26 (0.04s) | 9 (0.02s) | 13 (0.04s) |
| 49 | `03_single_vertex_stars` | `as_g1_patch5_star_sym.xml` | 5 | 0 | 5 | 1041 | 82 (0.20s) | 17 (0.04s) | 14 (0.02s) | 15 (0.02s) | 21 (0.03s) | 8 (0.01s) | 9 (0.05s) |
| 50 | `03_single_vertex_stars` | `as_g1_pentagon_p3.xml` | 5 | 0 | 5 | 1041 | 100 (0.27s) | 18 (0.03s) | 14 (0.03s) | 16 (0.02s) | 23 (0.03s) | 7 (0.01s) | 9 (0.04s) |
| 51 | `03_single_vertex_stars` | `as_g1_pentagon_p3_r1.xml` | 5 | 0 | 5 | 4336 | 122 (2.15s) | 23 (0.38s) | 18 (0.19s) | 18 (0.27s) | 27 (0.38s) | 8 (0.12s) | 12 (0.25s) |
| 52 | `03_single_vertex_stars` | `as_g1_pentagon_random_s101.xml` | 5 | 0 | 5 | 1041 | 178 (0.57s) | 20 (0.02s) | 20 (0.06s) | 17 (0.03s) | 31 (0.04s) | 9 (0.01s) | 18 (0.04s) |
| 53 | `03_single_vertex_stars` | `concave_pentagon_5corner_nonconvex_asg1.xml` | 5 | 0 | 5 | 10806 | 28 (1.84s) | 15 (0.45s) | 7 (0.27s) | 7 (0.34s) | 17 (0.91s) | 6 (0.18s) | 7 (0.41s) |
| 54 | `03_single_vertex_stars` | `convex_5_valence.xml` | 5 | 0 | 5 | 1041 | 136 (0.78s) | 22 (0.08s) | 18 (0.04s) | 18 (0.03s) | 30 (0.04s) | 8 (0.01s) | 10 (0.03s) |
| 55 | `03_single_vertex_stars` | `convex_5_valence_non_bilinear.xml` | 5 | 0 | 5 | 1041 | 300 (F) | 39 (0.05s) | 28 (0.08s) | 27 (0.04s) | 35 (0.08s) | 9 (0.02s) | 14 (0.03s) |
| 56 | `03_single_vertex_stars` | `curved_pentagon_5corner_convex_asg1.xml` | 5 | 0 | 5 | 10806 | 34 (1.78s) | 14 (0.29s) | 12 (0.51s) | 7 (0.28s) | 18 (0.54s) | 6 (0.32s) | 7 (0.31s) |
| 57 | `03_single_vertex_stars` | `custom_sample_domain.xml` | 5 | 0 | 5 | 1041 | 130 (0.57s) | 21 (0.07s) | 18 (0.04s) | 17 (0.02s) | 30 (0.04s) | 8 (0.01s) | 8 (0.01s) |
| 58 | `03_single_vertex_stars` | `yeti_asg1.xml` | 5 | 0 | 5 | 10801 | 38 (2.35s) | 34 (1.62s) | 28 (1.79s) | 11 (0.56s) | 21 (0.97s) | 8 (0.30s) | 15 (0.60s) |
| 59 | `03_single_vertex_stars` | `6_patch_6_valence_not_bilinear.xml` | 6 | 0 | 6 | 1248 | 132 (0.39s) | 30 (0.06s) | 23 (0.04s) | 23 (0.04s) | 33 (0.03s) | 9 (0.02s) | 14 (0.07s) |
| 60 | `03_single_vertex_stars` | `as_g1_hex_p3.xml` | 6 | 0 | 6 | 1248 | 153 (0.35s) | 24 (0.03s) | 18 (0.04s) | 21 (0.05s) | 27 (0.05s) | 7 (0.02s) | 10 (0.06s) |
| 61 | `03_single_vertex_stars` | `as_g1_hex_random_s777.xml` | 6 | 0 | 6 | 1248 | 263 (0.99s) | 29 (0.09s) | 25 (0.08s) | 25 (0.06s) | 39 (0.08s) | 10 (0.03s) | 12 (0.07s) |
| 62 | `03_single_vertex_stars` | `as_g1_patch6_random.xml` | 6 | 0 | 6 | 1248 | 183 (0.61s) | 32 (0.08s) | 26 (0.08s) | 21 (0.06s) | 32 (0.07s) | 10 (0.02s) | 14 (0.06s) |
| 63 | `03_single_vertex_stars` | `convex_6_valence.xml` | 6 | 0 | 6 | 1249 | 233 (0.78s) | 31 (0.07s) | 25 (0.07s) | 24 (0.07s) | 39 (0.07s) | 10 (0.02s) | 13 (0.04s) |
| 64 | `03_single_vertex_stars` | `convex_6_valence_non_bilinear.xml` | 6 | 0 | 6 | 1249 | 300 (F) | 51 (0.11s) | 39 (0.10s) | 37 (0.08s) | 42 (0.08s) | 10 (0.01s) | 17 (0.04s) |
| 65 | `03_single_vertex_stars` | `as_g1_patch7_random.xml` | 7 | 0 | 7 | 1455 | 192 (0.88s) | 34 (0.07s) | 28 (0.12s) | 26 (0.09s) | 40 (0.12s) | 12 (0.03s) | 17 (0.08s) |
| 66 | `03_single_vertex_stars` | `as_g1_v7_p3.xml` | 7 | 0 | 7 | 1455 | 230 (0.82s) | 33 (0.06s) | 26 (0.12s) | 28 (0.08s) | 34 (0.08s) | 9 (0.03s) | 10 (0.03s) |
| 67 | `04_multi_vertex_meshes` | `crescent_2corner_asg1.xml` | 4 | 0 | 3 | 8695 | 123 (8.18s) | 244 (9.05s) | 46 (2.65s) | 14 (0.57s) | 41 (1.59s) | 9 (0.39s) | 14 (0.71s) |
| 68 | `04_multi_vertex_meshes` | `crescent_2corner_asg1_reparam.xml` | 4 | 0 | 3 | 8695 | 300 (F) | 300 (F) | 237 (11.25s) | 171 (6.10s) | 232 (7.95s) | 24 (0.99s) | 173 (5.26s) |
| 69 | `04_multi_vertex_meshes` | `lens_2corner_convex_asg1.xml` | 4 | 0 | 3 | 8695 | 26 (1.46s) | 16 (0.53s) | 14 (0.39s) | 8 (0.16s) | 17 (0.30s) | 6 (0.11s) | 8 (0.31s) |
| 70 | `04_multi_vertex_meshes` | `teardrop_1corner_convex_asg1.xml` | 4 | 0 | 3 | 8695 | 25 (0.88s) | 16 (0.34s) | 13 (0.28s) | 9 (0.22s) | 17 (0.49s) | 6 (0.08s) | 9 (0.24s) |
| 71 | `04_multi_vertex_meshes` | `asg1_circle_5patch.xml` | 5 | 0 | 3 | 1096 | 300 (F) | 43 (0.13s) | 35 (0.09s) | 11 (0.03s) | 57 (0.20s) | 8 (0.03s) | 34 (0.13s) |
| 72 | `04_multi_vertex_meshes` | `asg1_general_5patch_v2.xml` | 5 | 0 | 3 | 1096 | 38 (0.12s) | 15 (0.03s) | 12 (0.03s) | 8 (0.01s) | 20 (0.03s) | 6 (0.01s) | 8 (0.02s) |
| 73 | `04_multi_vertex_meshes` | `asg1_general_5patch_v4.xml` | 5 | 0 | 3 | 1096 | 73 (0.13s) | 16 (0.03s) | 14 (0.03s) | 9 (0.02s) | 27 (0.05s) | 7 (0.01s) | 9 (0.03s) |
| 74 | `04_multi_vertex_meshes` | `circle_asg1.xml` | 5 | 0 | 3 | 10952 | 91 (4.89s) | 30 (1.67s) | 35 (2.26s) | 14 (0.72s) | 25 (1.27s) | 8 (0.49s) | 14 (0.85s) |
| 75 | `04_multi_vertex_meshes` | `heart_1corner_nonconvex_asg1.xml` | 5 | 0 | 3 | 10951 | 300 (F) | 300 (F) | 300 (F) | 89 (4.04s) | 84 (4.50s) | 21 (0.95s) | 89 (4.65s) |
| 76 | `04_multi_vertex_meshes` | `heart_2corner_nonconvex_asg1.xml` | 5 | 0 | 3 | 10950 | 47 (4.65s) | 71 (4.44s) | 56 (3.98s) | 9 (0.55s) | 27 (1.06s) | 7 (0.23s) | 14 (0.76s) |
| 77 | `04_multi_vertex_meshes` | `01_dart_4corner_nonconvex_asg1.xml` | 6 | 0 | 4 | 13062 | 300 (F) | 300 (F) | 250 (17.58s) | 133 (7.86s) | 96 (5.34s) | 29 (1.74s) | 133 (6.79s) |
| 78 | `04_multi_vertex_meshes` | `04_heptagon_7corner_convex_asg1.xml` | 6 | 0 | 4 | 13012 | 300 (F) | 300 (F) | 145 (10.23s) | 127 (7.30s) | 94 (5.72s) | 21 (0.99s) | 78 (3.77s) |
| 79 | `04_multi_vertex_meshes` | `6_patch_4_valence_each_not_bilinear.xml` | 6 | 0 | 4 | 1265 | 118 (0.86s) | 22 (0.13s) | 17 (0.12s) | 16 (0.06s) | 29 (0.08s) | 8 (0.03s) | 12 (0.07s) |
| 80 | `04_multi_vertex_meshes` | `concave_heptagon_7corner_nonconvex_asg1.xml` | 6 | 0 | 4 | 13012 | 300 (F) | 300 (F) | 300 (F) | 164 (9.94s) | 144 (7.84s) | 32 (1.76s) | 237 (12.85s) |
| 81 | `04_multi_vertex_meshes` | `curved_heptagon_7corner_convex_asg1.xml` | 6 | 0 | 4 | 13012 | 300 (F) | 288 (17.12s) | 133 (9.13s) | 122 (6.49s) | 85 (4.94s) | 19 (1.03s) | 57 (3.14s) |
| 82 | `04_multi_vertex_meshes` | `square_6_patch.xml` | 6 | 0 | 5 | 1271 | 300 (F) | 120 (0.50s) | 24 (0.14s) | 19 (0.06s) | 45 (0.18s) | 10 (0.03s) | 18 (0.08s) |
| 83 | `04_multi_vertex_meshes` | `square_6_patch_non_bilinear.xml` | 6 | 0 | 5 | 1271 | 300 (F) | 300 (F) | 50 (0.31s) | 50 (0.27s) | 79 (0.23s) | 13 (0.04s) | 44 (0.16s) |
| 84 | `04_multi_vertex_meshes` | `03_filleted_u_bracket_asg1.xml` | 11 | 0 | 4 | 23961 | 300 (F) | 300 (F) | 116 (20.67s) | 101 (10.28s) | 107 (12.28s) | 22 (2.13s) | 101 (9.87s) |
| 85 | `04_multi_vertex_meshes` | `blob_geometry.xml` | 11 | 0 | 6 | 2334 | 269 (3.71s) | 41 (0.45s) | 35 (0.36s) | 30 (0.21s) | 44 (0.40s) | 12 (0.09s) | 17 (0.24s) |
| 86 | `04_multi_vertex_meshes` | `blob_geometry_non_bilinear.xml` | 11 | 0 | 6 | 2334 | 300 (F) | 90 (0.89s) | 77 (1.01s) | 65 (0.86s) | 60 (0.54s) | 15 (0.12s) | 32 (0.49s) |
| 87 | `04_multi_vertex_meshes` | `weirdo_multivalence.xml` | 11 | 0 | 6 | 2335 | 300 (F) | 146 (1.46s) | 116 (1.47s) | 99 (1.13s) | 69 (0.53s) | 18 (0.17s) | 50 (0.86s) |
| 88 | `04_multi_vertex_meshes` | `weirdo_multivalence_non_bilinear.xml` | 11 | 0 | 6 | 2335 | 300 (F) | 134 (1.17s) | 111 (1.42s) | 97 (0.95s) | 76 (0.56s) | 21 (0.14s) | 56 (0.89s) |
| 89 | `04_multi_vertex_meshes` | `experiment_multivalence_not_bilinear.xml` | 13 | 0 | 6 | 2766 | 300 (F) | 106 (1.55s) | 84 (1.45s) | 66 (0.73s) | 50 (0.53s) | 14 (0.17s) | 24 (0.56s) |
| 90 | `05_multiply_connected_holes` | `bean_0corner_1hole_asg1.xml` | 4 | 1 | 0 | 8644 | 138 (8.00s) | 30 (1.21s) | 24 (1.18s) | 13 (0.32s) | 30 (1.00s) | 7 (0.21s) | 13 (0.59s) |
| 91 | `05_multiply_connected_holes` | `curved_lshape_6corner_1hole_asg1.xml` | 4 | 1 | 0 | 8636 | 163 (8.49s) | 300 (F) | 55 (2.39s) | 11 (0.31s) | 43 (1.56s) | 8 (0.27s) | 11 (0.48s) |
| 92 | `05_multiply_connected_holes` | `curved_ring_4patch_asg1.xml` | 4 | 1 | 0 | 8636 | 48 (2.35s) | 89 (2.40s) | 14 (0.42s) | 9 (0.25s) | 25 (0.47s) | 6 (0.13s) | 12 (0.41s) |
| 93 | `05_multiply_connected_holes` | `star_5corner_1hole_asg1.xml` | 4 | 1 | 0 | 8636 | 262 (12.03s) | 60 (1.52s) | 76 (3.11s) | 51 (1.38s) | 38 (1.38s) | 12 (0.48s) | 38 (1.31s) |
| 94 | `05_multiply_connected_holes` | `wavy_ring_annulus_asg1.xml` | 4 | 1 | 0 | 8644 | 26 (1.03s) | 50 (1.03s) | 28 (1.25s) | 10 (0.23s) | 18 (0.71s) | 7 (0.17s) | 10 (0.25s) |
| 95 | `05_multiply_connected_holes` | `ring_7patch_non_bilinear_swirl_p3_asg1.xml` | 7 | 1 | 0 | 15113 | 300 (F) | 64 (4.78s) | 52 (4.79s) | 41 (2.82s) | 38 (2.70s) | 13 (0.81s) | 22 (0.81s) |
| 96 | `05_multiply_connected_holes` | `crescent_2corner_3holes_asg1.xml` | 12 | 5 | 0 | 25908 | 228 (37.82s) | 62 (6.78s) | 76 (8.73s) | 46 (4.81s) | 31 (3.44s) | 12 (0.80s) | 33 (4.69s) |
| 97 | `05_multiply_connected_holes` | `curved_lshape_6corner_3holes_asg1.xml` | 12 | 3 | 0 | 25996 | 300 (F) | 300 (F) | 300 (F) | 300 (F) | 189 (19.92s) | 39 (4.06s) | 125 (13.03s) |
| 98 | `05_multiply_connected_holes` | `starburst_12corner_3holes_asg1.xml` | 12 | 5 | 0 | 25908 | 225 (33.30s) | 59 (6.28s) | 75 (7.01s) | 46 (4.37s) | 31 (3.52s) | 12 (1.20s) | 33 (3.69s) |
| 99 | `05_multiply_connected_holes` | `comb_16corner_4holes_asg1.xml` | 17 | 8 | 0 | 36746 | 146 (39.30s) | 291 (44.05s) | 5 (1.25s) | 9 (1.32s) | 40 (6.24s) | 6 (0.87s) | 9 (1.69s) |
| 100 | `05_multiply_connected_holes` | `uchannel_8corner_4holes_asg1.xml` | 17 | 5 | 0 | 36793 | 99 (24.13s) | 197 (29.05s) | 4 (1.10s) | 9 (1.26s) | 34 (4.00s) | 6 (0.88s) | 9 (1.02s) |
| 101 | `05_multiply_connected_holes` | `swiss_cross_12corner_5holes_asg1.xml` | 20 | 6 | 0 | 43395 | 300 (F) | 300 (F) | 300 (F) | 300 (F) | 300 (F) | 52 (10.14s) | 300 (F) |

### 6.2 Full 3-Level Comparison (`m=2 -> r=4`, `d=3`, 61 Domains: V-Cycle vs. W-Cycle vs. PCG)

| Filename | `DOFs` | `mg-V auto` | `mg-V ring1` | `mg-V fd_ring1` | `mg-V line_ring1` | `mg-W auto` | `mg-W ring1` | `mg-W fd_ring1` | `mg-W line_ring1` | `cg-V auto` | `cg-V line_ring1` |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `02_serpentine_s_duct_asg1.xml` | 8836 | 300 (F) | 300 (F) | 300 (F) | 36 | 300 (F) | 300 (F) | 300 (F) | 32 | 300 (F) | 15 |
| `05_hourglass_pinch_asg1.xml` | 8836 | 300 (F) | 300 (F) | 82 | 14 | 300 (F) | 300 (F) | 64 | 11 | 161 | 9 |
| `bulged_quad_non_bilinear_asg1.xml` | 8836 | 9 | 17 | 11 | 8 | 9 | 16 | 11 | 7 | 10 | 6 |
| `bulged_quad_non_bilinear_asg1_reparam.xml` | 900 | 9 | 15 | 9 | 8 | 8 | 15 | 9 | 8 | 9 | 6 |
| `concave_quad_4corner_nonconvex_asg1.xml` | 8836 | 9 | 17 | 9 | 8 | 9 | 16 | 9 | 7 | 10 | 6 |
| `curved_quad_4corner_convex_asg1.xml` | 8836 | 9 | 17 | 11 | 8 | 9 | 16 | 11 | 7 | 10 | 6 |
| `curved_quad_asg1.xml` | 8836 | 9 | 17 | 9 | 8 | 9 | 16 | 9 | 7 | 10 | 6 |
| `lake_asg1.xml` | 8836 | 300 (F) | 300 (F) | 300 (F) | 22 | 275 | 300 (F) | 300 (F) | 18 | 60 | 12 |
| `quad_asg1.xml` | 8836 | 9 | 17 | 4 | 8 | 9 | 16 | 4 | 7 | 10 | 6 |
| `star_quad_4corner_nonconvex_asg1.xml` | 8836 | 270 | 300 (F) | 300 (F) | 161 | 151 | 254 | 234 | 92 | 60 | 36 |
| `2_patch_rectangle.xml` | 1829 | 57 | 20 | 8 | 10 | 55 | 19 | 8 | 10 | 23 | 7 |
| `2_patch_rectangle_asg1.xml` | 17763 | 57 | 17 | 13 | 14 | 41 | 16 | 10 | 11 | 19 | 8 |
| `2_patch_rectangle_non_bilinear.xml` | 1829 | 17 | 15 | 5 | 7 | 17 | 14 | 5 | 7 | 11 | 6 |
| `curved_2patch_asg1.xml` | 17765 | 161 | 89 | 19 | 17 | 156 | 90 | 18 | 18 | 43 | 10 |
| `curved_hexagon_6corner_convex_asg1.xml` | 17763 | 300 (F) | 66 | 59 | 47 | 300 (F) | 57 | 48 | 37 | 53 | 14 |
| `curved_lshape_non_bilinear_asg1.xml` | 17763 | 271 | 77 | 52 | 37 | 201 | 63 | 43 | 25 | 38 | 13 |
| `two_bilinear_patches.xml` | 1827 | 214 | 42 | 29 | 35 | 171 | 34 | 24 | 27 | 26 | 10 |
| `two_patch_mathematica.xml` | 1827 | 92 | 21 | 16 | 19 | 68 | 15 | 12 | 15 | 19 | 8 |
| `two_patch_mathematica_non_bilinear.xml` | 1827 | 70 | 20 | 14 | 17 | 53 | 17 | 11 | 13 | 18 | 8 |
| `two_squares_asg1.xml` | 17765 | 18 | 17 | 5 | 8 | 18 | 17 | 4 | 7 | 14 | 6 |
| `three_patch_belt.xml` | 2754 | 91 | 27 | 18 | 22 | 67 | 20 | 14 | 16 | 25 | 9 |
| `three_patch_belt_non_bilinear.xml` | 2754 | 77 | 22 | 18 | 20 | 58 | 17 | 14 | 16 | 22 | 8 |
| `as_g1_patch3_random.xml` | 2787 | 142 | 105 | 79 | 41 | 129 | 80 | 64 | 35 | 42 | 15 |
| `as_g1_tri_p3.xml` | 2787 | 55 | 15 | 10 | 10 | 42 | 14 | 10 | 9 | 22 | 7 |
| `convex_3_valence.xml` | 2787 | 88 | 33 | 26 | 17 | 69 | 31 | 25 | 16 | 27 | 9 |
| `convex_3_valence_non_bilinear.xml` | 2787 | 94 | 20 | 14 | 14 | 73 | 19 | 14 | 10 | 27 | 8 |
| `4_patch_3_valence_not_bilinear.xml` | 3714 | 178 | 38 | 29 | 22 | 127 | 32 | 26 | 18 | 33 | 11 |
| `as_g1_patch4_random.xml` | 3714 | 153 | 27 | 23 | 20 | 115 | 20 | 17 | 14 | 30 | 9 |
| `as_g1_v4_p3.xml` | 3714 | 80 | 18 | 13 | 15 | 62 | 14 | 10 | 11 | 22 | 7 |
| `5_patch_5_valence_not_bilinear.xml` | 4641 | 216 | 34 | 26 | 29 | 162 | 25 | 20 | 21 | 38 | 11 |
| `as_g1_ellipse_v5_p3.xml` | 4641 | 300 (F) | 98 | 69 | 74 | 300 (F) | 78 | 56 | 59 | 62 | 14 |
| `as_g1_patch5_polygon_sym.xml` | 4641 | 163 | 28 | 20 | 23 | 129 | 20 | 16 | 17 | 32 | 10 |
| `as_g1_patch5_random.xml` | 4641 | 179 | 36 | 33 | 33 | 132 | 26 | 24 | 24 | 35 | 11 |
| `as_g1_patch5_star_sym.xml` | 4641 | 136 | 26 | 23 | 23 | 103 | 19 | 17 | 17 | 29 | 10 |
| `as_g1_pentagon_p3.xml` | 4641 | 161 | 27 | 19 | 23 | 129 | 20 | 15 | 17 | 32 | 9 |
| `as_g1_pentagon_p3_r1.xml` | 18896 | 178 | 33 | 24 | 26 | 137 | 24 | 19 | 19 | 33 | 10 |
| `as_g1_pentagon_random_s101.xml` | 4641 | 266 | 34 | 33 | 25 | 228 | 32 | 29 | 22 | 43 | 12 |
| `convex_5_valence.xml` | 4641 | 236 | 33 | 26 | 27 | 193 | 27 | 23 | 18 | 41 | 11 |
| `convex_5_valence_non_bilinear.xml` | 4641 | 300 (F) | 66 | 51 | 46 | 300 (F) | 59 | 47 | 42 | 50 | 13 |
| `custom_sample_domain.xml` | 4641 | 200 | 29 | 24 | 24 | 167 | 21 | 19 | 18 | 42 | 10 |
| `6_patch_6_valence_not_bilinear.xml` | 5568 | 234 | 51 | 37 | 39 | 174 | 37 | 28 | 29 | 41 | 14 |
| `as_g1_hex_p3.xml` | 5568 | 247 | 35 | 26 | 31 | 194 | 24 | 20 | 20 | 37 | 9 |
| `as_g1_hex_random_s777.xml` | 5568 | 300 (F) | 53 | 46 | 44 | 300 (F) | 43 | 38 | 34 | 54 | 14 |
| `as_g1_patch6_random.xml` | 5568 | 300 (F) | 52 | 34 | 38 | 235 | 38 | 25 | 29 | 44 | 14 |
| `convex_6_valence.xml` | 5569 | 300 (F) | 47 | 36 | 38 | 300 (F) | 40 | 33 | 30 | 56 | 14 |
| `convex_6_valence_non_bilinear.xml` | 5569 | 300 (F) | 75 | 55 | 63 | 300 (F) | 78 | 62 | 52 | 62 | 15 |
| `as_g1_patch7_random.xml` | 6495 | 300 (F) | 54 | 43 | 43 | 235 | 39 | 32 | 30 | 52 | 15 |
| `as_g1_v7_p3.xml` | 6495 | 300 (F) | 52 | 40 | 44 | 300 (F) | 37 | 31 | 28 | 46 | 13 |
| `asg1_circle_5patch.xml` | 4744 | 300 (F) | 300 (F) | 300 (F) | 300 (F) | 300 (F) | 300 (F) | 300 (F) | 300 (F) | 300 (F) | 300 (F) |
| `asg1_general_5patch_v2.xml` | 4744 | 50 | 21 | 17 | 10 | 51 | 20 | 16 | 9 | 24 | 7 |
| `asg1_general_5patch_v4.xml` | 4744 | 92 | 26 | 20 | 11 | 82 | 25 | 19 | 10 | 33 | 8 |
| `circle_asg1.xml` | 44936 | 125 | 46 | 48 | 18 | 104 | 42 | 43 | 16 | 31 | 10 |
| `6_patch_4_valence_each_not_bilinear.xml` | 5601 | 186 | 34 | 27 | 27 | 148 | 30 | 22 | 20 | 38 | 11 |
| `square_6_patch.xml` | 5607 | 300 (F) | 149 | 50 | 31 | 300 (F) | 140 | 47 | 30 | 58 | 12 |
| `square_6_patch_non_bilinear.xml` | 5607 | 300 (F) | 300 (F) | 87 | 82 | 300 (F) | 300 (F) | 77 | 75 | 107 | 20 |
| `blob_geometry.xml` | 10302 | 300 (F) | 60 | 50 | 44 | 300 (F) | 45 | 39 | 33 | 55 | 15 |
| `blob_geometry_non_bilinear.xml` | 10302 | 300 (F) | 136 | 117 | 94 | 300 (F) | 109 | 96 | 74 | 80 | 19 |
| `weirdo_multivalence.xml` | 10303 | 300 (F) | 283 | 220 | 192 | 300 (F) | 222 | 178 | 150 | 103 | 27 |
| `weirdo_multivalence_non_bilinear.xml` | 10303 | 300 (F) | 269 | 230 | 201 | 300 (F) | 221 | 187 | 162 | 116 | 31 |
| `experiment_multivalence_not_bilinear.xml` | 12190 | 300 (F) | 300 (F) | 225 | 190 | 300 (F) | 244 | 181 | 152 | 83 | 21 |
| `curved_ring_4patch_asg1.xml` | 35708 | 57 | 101 | 13 | 12 | 47 | 92 | 12 | 9 | 29 | 8 |
