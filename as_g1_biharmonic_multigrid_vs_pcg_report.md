# Root-Cause Analysis & Improvement Proposals: Direct Multigrid (`-s mg`) vs. Multigrid-PCG (`-s cg`) for AS-$G^1$ Biharmonic Problems

## 1. Executive Summary

Across the AS-$G^1$ multi-patch biharmonic benchmark suites in [as_g1_biharmonic_multigrid_example.cpp](file:///home/fhasanova/gismo/examples/as_g1_biharmonic_multigrid_example.cpp) and [as_g1_biharmonic2_multigrid_example.cpp](file:///home/fhasanova/gismo/examples/as_g1_biharmonic2_multigrid_example.cpp), **Multigrid-Preconditioned Conjugate Gradient (`-s cg`)** converges reliably across all benchmark geometries ($18\text{--}120$ iterations, $100\%$ pass rate at $\nu=1$ V-cycle), whereas **Direct Stationary Multigrid (`-s mg`)** with $\nu=1$ Gauss–Seidel smoothing requires $89\text{--}452$ iterations on simple geometries and exceeds $500$ iterations (failing to reach $\text{tol}=10^{-8}$) on stretched, non-bilinear, or multi-valence domains.

Spectral and mode-localization diagnostics of the V-cycle operator $M_{\text{MG}} A$ explain both *why* stationary MG stalls and *why* PCG dramatically outperforms it:

1. **Outlier-Dominated Spectrum of $M_{\text{MG}} A$**:
   - Over $96\text{--}99\%$ of eigenvalues on regular 2–4 patch geometries (and $86\text{--}87\%$ on 6-patch geometries) already lie in the well-conditioned bulk interval $[0.50, 1.00]$ ($\kappa_{\text{bulk}} \le 2$).
   - Only a **small cluster of $4\text{--}26$ outlier eigenvalues** lies near zero ($\lambda_{\min}(M_{\text{MG}} A) \in [0.0088, 0.0697]$).
   - Stationary MG (`-s mg`) is a memoryless fixed-point iteration $e^{(k)} = (I - M_{\text{MG}} A)^k e^{(0)}$ permanently governed by the single worst eigenvalue: $\rho(I - M_{\text{MG}} A) = 1 - \lambda_{\min} \approx 0.93\text{--}0.991$, requiring $\approx 18.42 / \lambda_{\min}$ iterations to reach $\text{tol}=10^{-8}$.
   - PCG (`-s cg`) provides both **$\mathcal{O}(\kappa^{-1/2})$ Chebyshev acceleration** and **Kaniel–Paige–Saad Krylov outlier deflation**: after placing polynomial roots on the few isolated small eigenvalues, PCG converges at the fast bulk rate.
2. **Primary Geometric/Algebraic Bottleneck — Cross-Interface Bending Locking (`EdgeTrace`/`Corner6Jet` $\leftrightarrow$ `Int(ring0)`)**:
   - In [gsAsG1Basis.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Basis.hpp), interface trace functions ($\mathcal{V}_{E,0}$) and crossing-derivative functions ($\mathcal{V}_{E,1}$) have support confined to the **first two boundary B-spline layers** ($i_n \in \{0,1\}$) and vanish on layer $i_n = 2$.
   - The biharmonic operator $\int_\Omega (\Delta u)^2$ involves the second normal derivative $\partial_{nn} u$, which across the first knot span $[0, h]$ strongly couples layers $i_n \in \{0,1\}$ (`EdgeTrace`, `EdgeCross`, `Corner6Jet`) with the **first two interior B-spline layers** $i_n \in \{2,3\}$ (`Int(ring0)`).
   - Pointwise Gauss–Seidel relaxes `EdgeTrace`/`Corner6Jet` and `Int(ring0)` one scalar DOF at a time while holding the adjacent ring fixed, causing $\mathcal{O}(h^{-4})$ normal-curvature locking across the interface.
3. **Secondary Bottlenecks — Unscaled 6-Jet Vertex Coupling, Patch Aspect-Ratio Anisotropy on $C^1$ Multiplicity-2 Interiors, and Boundary/Corner Regularity**:
   - At vertices (`Corner6Jet`), unscaled physical derivatives $(u, \nabla u, \nabla^2 u)$ create diagonal entries in $K_{\text{free}}$ spanning up to $3.6 \times 10^5\times$ on `square_6_patch.xml`, combined with strong $u_{xx}$–$u_{yy}$ coupling and coarse/fine vertex support mismatch.
   - Stretched or distorted patches (`square_6_patch.xml` with $|J_u|/|J_v| = 3.0$, `square_6_patch_non_bilinear.xml` with $|J_u|/|J_v| = 5.154$) amplify directional biharmonic stiffness as $(|J_u|/|J_v|)^4 \approx 81\text{--}705$, stalling pointwise Gauss–Seidel on interior $C^1$ (knot-multiplicity 2) modes.
   - Non-convex re-entrant corners and non-$G^2$ $180^\circ$ boundary junctions reduce dual $H^3(\Omega)$ regularity, inflating the two-grid approximation constant $C_A$.

Targeted $\mathcal{O}(N)$ smoother improvements — specifically adding a **1D interface-corridor (`Coupled Ring1`) subspace solve** and **alternating 1D line-block Gauss–Seidel** on patch interiors — reduce Direct MG iterations on `two_bilinear_patches.xml` ($r=2\to 4$) from **$338 \to 35$** ($9.6\times$ fewer iterations) and on `square_6_patch.xml` ($r=2\to 4$) from **$>500$ (failure) $\to 30$ iterations** ($0.179\text{ s}$) with $\nu = 1$ V-cycles.

---

## 2. Exact Spectral Analysis of $M_{\text{MG}} A$ vs. $I - M_{\text{MG}} A$

In [gsMultiGrid.hpp](file:///home/fhasanova/gismo/src/gsMultiGrid/gsMultiGrid.hpp#L148-L201), `gsMultiGridOp` uses symmetric pre/post Gauss–Seidel smoothing (`m_symmSmooth = true`, forward sweep in `step`, reverse sweep in `stepT`), Galerkin coarse-grid operators $K_{\ell-1} = (P_{\ell-1}^\ell)^T K_\ell P_{\ell-1}^\ell$, and exact coarse solves. Therefore, the V-cycle preconditioner $M_{\text{MG}}$ is symmetric positive definite, and $M_{\text{MG}} A$ has real eigenvalues in $(0, 1]$ (verified explicitly by assembling $M_{\text{MG}} A$ on the first biharmonic system $K_{\text{free}}$ with $\beta = 0, p = 3, \nu = 1$ V-cycle):

| Geometry (`filedata/domain2d/2patch/...`) | Levels | $N_{\text{free}}$ | $\lambda_{\min}(M_{\text{MG}} A)$ | $\lambda_{\max}(M_{\text{MG}} A)$ | $\kappa(M_{\text{MG}} A)$ | $\rho(I - M_{\text{MG}} A)$ | Eigs $< 0.10$ | Eigs $< 0.20$ | Eigs $\in [0.5, 1.0]$ | Direct MG Iters ($\nu=1$) | PCG Iters ($\nu=1$) |
| :--- | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `2_patch_rectangle.xml` | $2 \to 4$ | $1,829$ | $0.13989$ | $1.00000$ | $7.15$ | $0.86011$ | $0$ ($0.00\%$) | $3$ ($0.16\%$) | $1,802$ ($98.5\%$) | $103$ | $23$ |
| `two_bilinear_patches.xml` | $2 \to 3$ | $403$ | $0.06966$ | $1.00000$ | $14.35$ | $0.93034$ | $1$ ($0.25\%$) | $2$ ($0.50\%$) | $398$ ($98.8\%$) | — | $18$* |
| `two_bilinear_patches.xml` | $2 \to 4$ | $1,827$ | $0.04735$ | $1.00000$ | $21.12$ | $0.95265$ | $2$ ($0.11\%$) | $4$ ($0.22\%$) | $1,812$ ($99.2\%$) | $338$ | $26$ |
| `4_patch_3_valence_not_bilinear.xml` | $2 \to 3$ | $834$ | $0.06069$ | $1.00000$ | $16.48$ | $0.93931$ | $1$ ($0.12\%$) | $5$ ($0.60\%$) | $808$ ($96.9\%$) | — | $27$* |
| `square_6_patch.xml` | $2 \to 3$ | $1,271$ | $0.02387$ | $1.00000$ | $41.90$ | $0.97613$ | $12$ ($0.94\%$) | $37$ ($2.91\%$) | $1,106$ ($87.0\%$) | $>500$ ($r=2\to 4$) | $49$* |
| `square_6_patch_non_bilinear.xml` | $2 \to 3$ | $1,271$ | $0.00880$ | $1.00000$ | $113.58$ | $0.99120$ | $26$ ($2.05\%$) | $49$ ($3.86\%$) | $1,099$ ($86.5\%$) | $>500$ ($r=2\to 4$) | $101$ ($r=2\to 4$) |

*\*PCG counts marked with an asterisk for $r=2\to 3$ are from [as_g1_biharmonic_multigrid_results.tex](file:///home/fhasanova/gismo/as_g1_biharmonic_multigrid_results.tex) and [numerical_experiments_as_g1_multigrid.tex](file:///home/fhasanova/gismo/numerical_experiments_as_g1_multigrid.tex); all spectral quantities and $r=2\to 4$ counts were measured directly.*

### Why This Spectrum Crushes `-s mg` While `-s cg` Thrives

1. **Linear ($\kappa$) vs. Square-Root ($\sqrt{\kappa}$) Scaling**:
   - For stationary MG ($u^{(k+1)} = u^{(k)} + M_{\text{MG}}(f - A u^{(k)})$), the error reduction factor is bounded below by $\rho(I - M_{\text{MG}} A) = 1 - \lambda_{\min}$. Reaching $\text{tol} = 10^{-8}$ requires:
     $$k_{\text{MG}} \approx \frac{\ln(10^8)}{-\ln(1 - \lambda_{\min})} \approx \frac{18.42}{\lambda_{\min}} = 18.42 \, \kappa(M_{\text{MG}} A).$$
     For `two_bilinear_patches.xml` ($r=2\to 4$, $\lambda_{\min} = 0.04735$), $18.42 / 0.04735 = 389$ (matching the observed $338$ iterations). For `square_6_patch_non_bilinear.xml` ($\lambda_{\min} \le 0.00880$), $18.42 / 0.00880 \approx 2,093$ iterations — far past `--maxIt 300` or `500`.
   - For PCG, Chebyshev optimality improves the asymptotic worst-case bound from $\mathcal{O}(\kappa)$ to $\mathcal{O}(\sqrt{\kappa})$:
     $$k_{\text{PCG}} \le \frac{\ln(2 \cdot 10^8)}{2} \sqrt{\kappa} \approx 9.56 \sqrt{\kappa}.$$
2. **Krylov Outlier Deflation (Superlinear Convergence)**:
   - Because $98\text{--}99\%$ of eigenvalues on 2-patch geometries lie in $[0.50, 1.00]$ ($\kappa_{\text{bulk}} \le 2.0$), CG needs only a few iterations to place zeros of the residual polynomial $p_k(\lambda)$ near the $2\text{--}5$ smallest eigenvalues ($\{0.047, 0.075, 0.111, 0.169\}$ on `two_bilinear_patches.xml`), after which it converges at the bulk rate $\frac{\sqrt{2}-1}{\sqrt{2}+1} \approx 0.172$ per step.
   - Stationary MG has **no polynomial memory** and cannot deflate isolated bad modes: every iteration attenuates the slowest mode by the fixed factor $1 - \lambda_{\min}$.

---

## 3. Localization of the Slowest Eigenmodes & Stiffness Scaling

### 3.1 Diagonal Stiffness Imbalance Across DOF Types

Extracting the diagonal of $K_{\text{free}}$ across the four Argyris DOF classes at $r=3$ ($p=3$, first biharmonic) reveals extreme scale separation:

| Geometry ($r=3$) | `Interior` $\text{diag}(K)$ | `EdgeTrace` ($\mathcal{V}_{E,0}$) $\text{diag}(K)$ | `EdgeCross` ($\mathcal{V}_{E,1}$) $\text{diag}(K)$ | `Corner6Jet` ($\mathcal{V}_v$) $\text{diag}(K)$ | `Corner6Jet` Max/Min |
| :--- | :---: | :---: | :---: | :---: | ---: |
| `two_bilinear_patches.xml` ($N=403$) | $[9.258\text{e+}1,\, 1.946\text{e+}2]$ | $[2.361\text{e+}2,\, 2.812\text{e+}2]$ | $[3.064\text{e-}1,\, 3.357\text{e-}1]$ | (all clamped) | — |
| `two_bilinear_patches.xml` ($r=4, N=1827$) | $[3.620\text{e+}2,\, 8.368\text{e+}2]$ | $[9.204\text{e+}2,\, 1.210\text{e+}3]$ | $[3.042\text{e-}1,\, 3.435\text{e-}1]$ | (all clamped) | — |
| `2_patch_rectangle.xml` ($N=405$) | $1.884\text{e+}2$ | $3.196\text{e+}2$ | $4.368\text{e-}1$ | $5.243\text{e-}1$ (Case B $\partial_{nn} u$) | $1.0\times$ |
| `4_patch_3_valence_not_bilinear.xml` ($N=834$) | $[3.191\text{e+}1,\, 1.315\text{e+}2]$ | $[5.711\text{e+}1,\, 1.193\text{e+}2]$ | $[6.028\text{e-}2,\, 1.493\text{e-}1]$ | $[5.551\text{e-}1,\, 4.686\text{e+}2]$ | $844.1\times$ |
| `square_6_patch.xml` ($N=1271$) | $[1.899\text{e+}2,\, 1.709\text{e+}3]$ | $[4.272\text{e+}2,\, 1.648\text{e+}3]$ | $[3.500\text{e-}1,\, 3.468\text{e+}0]$ | $[1.983\text{e-}2,\, 7.170\text{e+}3]$ | $361,552.9\times$ |
| `square_6_patch_non_bilinear.xml` ($N=1271$) | $[4.969\text{e+}1,\, 7.818\text{e+}2]$ | $[1.458\text{e+}2,\, 9.283\text{e+}2]$ | $[9.783\text{e-}2,\, 1.303\text{e+}0]$ | $[3.447\text{e-}2,\, 2.817\text{e+}3]$ | $81,718.7\times$ |

Notice that `Interior` and `EdgeTrace` stiffness scales as $\mathcal{O}(h^{-2})$ (from $h^2 \cdot h^{-4}$ in 2D), whereas `EdgeCross` scales as $\mathcal{O}(1)$ (since crossing-derivative basis functions carry an implicit factor of $h$), and `Corner6Jet` spans **five orders of magnitude** because $u$, $\nabla u$, and $\nabla^2 u$ scale as $1$, $h$, and $h^2$ in amplitude and their $K$-diagonals scale as $h^{-2}$, $1$, and $h^2$ respectively.

### 3.2 Spatial & DOF-Class Decomposition of the Slowest Eigenmodes of $M_{\text{MG}} A$

Measuring the diagonal energy share $v_i^2 K_{ii} / \sum_j v_j^2 K_{jj}$ for the eigenvectors of $M_{\text{MG}} A$ corresponding to the smallest eigenvalues $\lambda_k$ exposes three distinct physical failure modes:

1. **Interface Bending Modes (`EdgeTrace` $\leftrightarrow$ `Int(ring0)`):**
   - On `two_bilinear_patches.xml` ($r=2\to 4$, $N=1827$, no free corner DOFs), **every single one of the 9 slowest modes** ($\lambda_1 = 0.0474$ through $\lambda_9 = 0.3778$) has its largest DOF amplitude on `EdgeTrace` (`1803..1811`), with $12.1\text{--}24.4\%$ energy on the adjacent interior ring `Int(ring0)` ($i_n \in \{2,3\}$) and $54.0\text{--}82.9\%$ extending into the interior `Int(ring>=1)`.
   - On `4_patch_3_valence_not_bilinear.xml` ($r=2\to 3$, $N=834$), 7 of the 10 slowest modes ($\lambda_1 = 0.0607$ to $\lambda_5 = 0.1673$, plus $\lambda_7, \lambda_{10}$) are dominated by `EdgeTrace` coupled to `Int(ring0)` ($25.8\text{--}34.8\%$).
2. **Vertex 6-Jet + Ring-0 Modes (`Corner6Jet` $\leftrightarrow$ `Int(ring0)` + `EdgeTrace`):**
   - On `square_6_patch.xml` ($r=2\to 3$, $N=1271$), the 4 slowest modes ($\lambda_1 = 0.0239$, $\lambda_2 = 0.0352$, $\lambda_4 = 0.0504$, $\lambda_5 = 0.0573$) peak directly on **`Corner6Jet`** DOFs (`1246, 1203, 1217, 1270`), carrying $17.4\text{--}28.5\%$ energy in `Corner6Jet` and $30.9\text{--}46.8\%$ in `Int(ring0)`.
   - Compare `2_patch_rectangle.xml` ($r=2\to 3$, $\lambda_1 = 0.1710$, $\lambda_2 = 0.2499$): the only difference in topology relative to `two_bilinear_patches.xml` ($\lambda_1 = 0.0697$) is that `2_patch_rectangle.xml` has two smooth boundary vertices (Case B free $\partial_{nn} u$) and symmetric bilinear patches without distortion — yet its mode 2 (`0.2499`) still concentrates $33.4\%$ energy in `Corner6Jet` + $36.6\%$ in `Int(ring0)`.
3. **Anisotropic Patch-Interior Modes (`Int(ring>=1)` on stretched patches):**
   - On `square_6_patch_non_bilinear.xml` ($r=2\to 3$), Modes 1–5 ($\lambda = 0.0088\text{--}0.0176$) are `EdgeTrace`/`Corner6Jet` modes, while Modes 6–10 ($\lambda = 0.0264\text{--}0.0317$) have **$98.6\text{--}99.9\%$ of their energy inside the patch interiors** (notably Patch 5, where the Jacobian aspect ratio $|J_u|/|J_v| = 5.154$).

---

## 4. Detailed Root-Cause Breakdown

### Root Cause 1: Cross-Interface Normal-Curvature Locking Between Argyris Edge/Corner Spaces and `Int(ring0)`
In `deriveArgyrisBasisEmbedding` ([gsAsG1Basis.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Basis.hpp#L540-L582)), each patch's tensor-product B-spline DOFs are partitioned into:
- `deriveInnerEmbedding(tensorBasis)` ([gsAsG1Basis.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Basis.hpp#L407-L441)), which keeps tensor indices $i_1, i_2 \in [2, n-3]$ (`Interior`), and
- `deriveEdgeEmbedding` / `deriveCornerEmbedding`, which act **exclusively** on the outer two B-spline layers ($i_n \in \{0, 1\}$ and $\{n-2, n-1\}$).

For a $C^1$ cubic spline ($p=3$, interior knot multiplicity $2$), the second normal derivative $\partial_{nn} u$ on the interface element $[0, h]$ depends on B-splines $B_0, B_1, B_2, B_3$. While `EdgeTrace` ($\mathcal{V}_{E,0}$) and `EdgeCross` ($\mathcal{V}_{E,1}$) reconstruct combinations of $B_0$ and $B_1$, the interior splines $B_2, B_3$ belong to `Int(ring0)`. Because `EdgeTrace` has $C^2$ continuity along the edge and steep normal decay across $[0, h]$, updating an `EdgeTrace` DOF while holding `Int(ring0)` fixed generates a large artificial $\partial_{nn} u$ jump inside $[0,h]$, causing strong algebraic coupling between `EdgeTrace` and `Int(ring0)`.

This is **quantitatively proven** by our smoother ablation (`task-81` / `task-85`):
- Adding an exact subspace solve on `(Ifc + Crn)` **alone** (without `Int(ring0)`) barely helps `two_bilinear_patches.xml`: Direct MG drops only from $338 \to 318$ iterations ($\kappa: 21.11 \to 19.79$).
- Adding the exact solve on `(Ifc + Crn + Int(ring0))` (`Coupled Ring1`) drops Direct MG from **$338 \to 43$ iterations** ($\kappa: 21.11 \to 3.27$)!

### Root Cause 2: Patch Geometric Anisotropy and Multiplicity-2 ($C^1$) Interior Basis Functions
Unlike Poisson ($H^1$) discretizations where stiffness scales with $(|J_u|/|J_v|)^2$, the biharmonic operator ($H^2$) scales with the **fourth power** of the Jacobian aspect ratio: $(|J_u|/|J_v|)^4$.
- On `square_6_patch.xml`, Patches 1 and 5 have aspect ratio $3.0$, giving $3^4 = 81\times$ directional stiffness contrast.
- On `square_6_patch_non_bilinear.xml`, Patch 5 has $|J_u|/|J_v| = 5.154$, giving $5.154^4 \approx 705\times$ directional stiffness contrast.
- On `weirdo_multivalence.xml`, Patch 4 has $|J_u|/|J_v| = 4.425$, giving $4.425^4 \approx 383\times$ contrast.

Furthermore, `dbasis.uniformRefine(1, 2)` inserts knots with **multiplicity 2**, so $C^1$ cubic splines have pairs of closely overlapping functions at every knot line ($2 \times 2$ blocks per cell). Pointwise Gauss–Seidel cannot smooth oscillatory modes along the weak direction of anisotropic patches or across $C^1$ double-knot pairs. Consequently, even after fixing the interface (`Coupled Ring1`), `square_6_patch_non_bilinear.xml` still has $\kappa = 64.21$ unless the interior smoother handles anisotropy (e.g., alternating 1D line-block Gauss–Seidel drops $\kappa$ to $6.81$ and Direct MG from $>500 \to 80$ iterations).

### Root Cause 3: Multi-Scale Vertex Coupling & Unscaled 6-Jet Basis
At each vertex (`Corner6Jet`), the 6 global DOFs represent physical values $(u, \partial_x u, \partial_y u, \partial_{xx} u, \partial_{yy} u, \partial_{xy} u)$ (or local $(\bn, \bt)$ derivatives).
1. **Unscaled diagonal entries**: As seen in Section 3.1, `diag(K)` varies by up to $3.6\times 10^5$ across the 6 DOFs at a vertex.
2. **Strong off-diagonal Hessian coupling**: Because the bilinear form integrates $(\Delta u)^2 = (u_{xx} + u_{yy})^2$, $u_{xx}$ and $u_{yy}$ at a vertex are strongly negatively correlated in the energy norm, and `deriveCornerEmbedding` ([gsAsG1Basis.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Basis.hpp#L448-L486)) couples all 6 functionals across all patches incident on a high-valence extraordinary vertex.
3. **Coarse-to-fine corner support mismatch**: Each level's corner basis functions are supported only on the two boundary layers $[0, 2h]$ of that level, so a coarse corner basis function (supported on $[0, 2H] = [0, 4h]$) cannot be represented solely by fine corner functions — its prolongation spills strongly into fine `Int(ring0)` and near-corner `EdgeTrace` DOFs.

### Root Cause 4: Implementation Details & First vs. Second Biharmonic Differences
1. **Half-strength Gauss–Seidel in the drivers**:
   - Line 495 of [as_g1_biharmonic_multigrid_example.cpp](file:///home/fhasanova/gismo/examples/as_g1_biharmonic_multigrid_example.cpp#L495) and line 523 of [as_g1_biharmonic2_multigrid_example.cpp](file:///home/fhasanova/gismo/examples/as_g1_biharmonic2_multigrid_example.cpp#L523) instantiate `makeGaussSeidelOp(mg->matrix(l))` (`gsGaussSeidel::forward`), not `makeSymmetricGaussSeidelOp`. Because `gsMultiGridOp` sets `m_symmSmooth = true`, $\nu=1$ executes **one forward sweep** on the way down and **one backward sweep** on the way up (preserving symmetry of $M_{\text{MG}}$, which is why PCG works cleanly), whereas `makeSymmetricGaussSeidelOp` doubles the smoothing work per call (reducing `two_bilinear_patches.xml` direct MG from $338 \to 214$ iterations).
   - Also, the smoother loop runs `for (index_t l = 0; l < nLevels; ++l)` even though level `0` uses the exact coarse solver, and `mg->setCoarseSolver(makeSparseCholeskySolver(mg->matrix(0)))` is omitted (so `gsMultiGridOp::initCoarseSolver()` in [gsMultiGrid.hpp](file:///home/fhasanova/gismo/src/gsMultiGrid/gsMultiGrid.hpp#L102-L122) falls back to `SparseLU`).
2. **Aggressive sparse drop tolerance (`1e-4`) in `gsAsG1Basis.hpp`**:
   - In `collocateCorners` ([gsAsG1Basis.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Basis.hpp#L287)), the first-biharmonic branch (`stride == 2`) uses `.sparseView(1, T(1e-4))`, whereas the second-biharmonic branch uses `T(1e-12)`. Similarly, `deriveEdgeEmbedding` ([gsAsG1Basis.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Basis.hpp#L304)) uses default `eps = 1e-4`.
   - Because derivative collocation rows scale as $\mathcal{O}(h)$ and $\mathcal{O}(h^2)$, thresholding entries below $10^{-4}$ can drop genuine entries on distorted patches (`4_patch_3_valence_not_bilinear.xml` exhibits a $7.83\times 10^{-11}$ prolongation nesting error $\|T_f P - R T_c\| / \|R T_c\|$ vs. $\sim 10^{-16}$ on `2_patch_rectangle.xml`).
3. **Why Biharmonic 1 (`338` MG iters) is harder than Biharmonic 2 (`89` MG iters) on `two_bilinear_patches.xml`**:
   - In the first biharmonic problem ($u = \partial_n u = 0$ on $\partial\Omega$), clamping $\partial_n u = 0$ eliminates all outer-boundary `EdgeCross` DOFs and forces boundary-adjacent interior bending onto the remaining `EdgeTrace` + `Int(ring0)` interface coupling, raising $\kappa(M_{\text{MG}} A)$ and shifting the lowest mode deeper into the interior-interface interaction.
4. **Non-convex / re-entrant corner regularity loss**:
   - On domains with re-entrant corners ($\omega > 180^\circ$, e.g., `weirdo_multivalence.xml`, `experiment_multivalence_not_bilinear.xml`) or $180^\circ$ boundary junctions that are only $G^1$ (not $G^2$), the biharmonic dual problem loses full $H^3(\Omega)$ regularity, enlarging the two-grid approximation constant $C_A$ and requiring stronger smoothing or Krylov acceleration.

---

## 5. Benchmark Comparison of Candidate Improvements

To verify which improvements actually fix standalone MG (`-s mg`) without breaking $\mathcal{O}(N)$ complexity, we benchmarked candidate smoothers on the first biharmonic problem ($p=3, r=2\to 4, \nu=1$, V-cycle, $\text{tol}=10^{-8}$, B-spline bubble RHS):

| Geometry | Smoother Variant ($\nu=1$ V-cycle unless noted) | Smoother Setup | $\kappa_{\text{Lanczos}}$ | Direct MG Iters (Time) | PCG Iters (Time) |
| :--- | :--- | ---: | ---: | ---: | ---: |
| `two_bilinear_patches.xml` ($r=2\to 4$) | Baseline Point GS (`makeGaussSeidelOp`) | $0.000\text{ s}$ | $21.11$ | $338$ ($0.402\text{ s}$) | $26$ ($0.028\text{ s}$) |
| | Point GS + `(Ifc+Crn)` only (no `Ring0`) | — | $19.79$ | $318$ ($0.528\text{ s}$) | $24$ ($0.037\text{ s}$) |
| | **Point GS + `Coupled Ring1` Block** | $0.001\text{ s}$ | **$3.27$** | **$43$ ($0.049\text{ s}$)** | **$13$ ($0.015\text{ s}$)** |
| | **Alternating 1D Line GS + `Coupled Ring1`** | $0.014\text{ s}$ | **$2.78$** | **$35$ ($0.043\text{ s}$)** | **$10$ ($0.013\text{ s}$)** |
| `two_bilinear_patches.xml` ($r=2\to 5$) | Baseline Point GS | $0.000\text{ s}$ | $29.04$ | $477$ ($3.636\text{ s}$) | $34$ ($0.270\text{ s}$) |
| | **Point GS + `Coupled Ring1` Block** | $0.004\text{ s}$ | **$4.39$** | **$62$ ($0.463\text{ s}$)** | **$15$ ($0.102\text{ s}$)** |
| | **Alternating 1D Line GS + `Coupled Ring1`** | $0.157\text{ s}$ | **$3.68$** | **$50$ ($0.359\text{ s}$)** | **$12$ ($0.092\text{ s}$)** |
| `4_patch_3_valence_not_bilinear.xml` | Baseline Point GS | $0.000\text{ s}$ | $23.67$ | $329$ ($1.001\text{ s}$) | $34$ ($0.119\text{ s}$) |
| | **Point GS + `Coupled Ring1` Block** | $0.004\text{ s}$ | **$3.07$** | **$41$ ($0.155\text{ s}$)** | **$16$ ($0.061\text{ s}$)** |
| | **Alternating 1D Line GS + `Coupled Ring1`** | $0.071\text{ s}$ | **$2.46$** | **$25$ ($0.087\text{ s}$)** | **$12$ ($0.042\text{ s}$)** |
| `square_6_patch.xml` | Baseline Point GS | $0.000\text{ s}$ | $87.15$ | $>500$ (FAIL, `1.1e-4`) | $58$ ($0.330\text{ s}$) |
| | **Point GS + `Coupled Ring1` Block** | $0.006\text{ s}$ | $21.59$ | $187$ ($1.184\text{ s}$) | $34$ ($0.188\text{ s}$) |
| | **Alternating 1D Line GS + `Coupled Ring1`** | $0.186\text{ s}$ | **$2.27$** | **$30$ ($0.179\text{ s}$)** | **$12$ ($0.066\text{ s}$)** |
| `square_6_patch_non_bilinear.xml` | Baseline Point GS | $0.000\text{ s}$ | $155.54$ | $>500$ (FAIL, `1.2e-4`) | $101$ ($0.555\text{ s}$) |
| | **Point GS + `Coupled Ring1` Block** | $0.007\text{ s}$ | $64.21$ | $500$ (`2.0e-8`) | $57$ ($0.307\text{ s}$) |
| | **Alternating 1D Line GS + `Coupled Ring1`** | $0.175\text{ s}$ | **$6.81$** | **$80$ ($0.415\text{ s}$)** | **$20$ ($0.104\text{ s}$)** |
| `weirdo_multivalence.xml` | Baseline Point GS | $0.000\text{ s}$ | $270.54$ | $>500$ (FAIL, `2.1e-4`) | $102$ ($1.092\text{ s}$) |
| | **Point GS + `Coupled Ring1` Block** | $0.011\text{ s}$ | **$29.05$** | **$308$ ($3.548\text{ s}$)** | **$43$ ($0.459\text{ s}$)** |
| | **Alternating 1D Line GS + `Coupled Ring1`** | $0.638\text{ s}$ | **$16.31$** | **$205$ ($1.971\text{ s}$)** | **$28$ ($0.282\text{ s}$)** |
| `experiment_multivalence_not_bilinear.xml` | Baseline Point GS | $0.000\text{ s}$ | $382.51$ | $>500$ (FAIL, `2.3e-5`) | $83$ ($1.126\text{ s}$) |
| | **Point GS + `Coupled Ring1` Block** | $0.015\text{ s}$ | **$48.66$** | **$372$ ($5.156\text{ s}$)** | **$42$ ($0.646\text{ s}$)** |
| | **Alternating 1D Line GS + `Coupled Ring1`** | $0.874\text{ s}$ | **$25.97$** | **$198$ ($2.440\text{ s}$)** | **$21$ ($0.272\text{ s}$)** |

*(Note: `makeIncompleteLUOp` with `fillfactor=1` diverged with `-nan` across all tested biharmonic matrices because $K_{\text{free}}$ is strongly off-diagonally coupled and contains orders-of-magnitude scale jumps at interfaces/vertices, causing unstable incomplete pivots.)*

---

## 6. Concrete Improvement Proposals

### Proposal 1 (Highest Return-on-Effort, $\mathcal{O}(N)$ Complexity): Multiplicative `Point GS + Coupled Ring1` Skeleton Smoother
- **Idea**: Keep the standard pointwise Gauss–Seidel sweep across all DOFs, and follow it (in reverse order during `stepT` so the smoother remains symmetric) with an **exact sparse Cholesky subspace correction** on the coupled skeleton index set $I_{\text{ring1}}$:
  $$I_{\text{ring1}} = \{\text{free `EdgeTrace`}\} \cup \{\text{free `EdgeCross`}\} \cup \{\text{free `Corner6Jet`}\} \cup \{\text{first 2 interior B-spline layers } i_n \in \{2, 3\} \text{ adjacent to every interface}\}.$$
- **Why it works & stays $\mathcal{O}(N)$**:
  - Including `Int(ring0)` ($i_n \in \{2,3\}$) directly resolves the normal-curvature ($\partial_{nn} u$) locking across every interface and vertex simultaneously.
  - Because $I_{\text{ring1}}$ lives on a 1D patch-interface skeleton, $|I_{\text{ring1}}| = \mathcal{O}(\sqrt{N})$, and its stiffness submatrix $A_{I_{\text{ring1}}, I_{\text{ring1}}}$ has 1D-graph sparsity. Factoring it takes $1\text{--}15\text{ ms}$ ($<2\%$ of setup time), and applying it takes $\mathcal{O}(\sqrt{N})$ per sweep.
- **Impact**: Cuts Direct MG iterations by $7\text{--}8\times$ on well-shaped domains ($338 \to 43$ on `two_bilinear_patches`, $329 \to 41$ on `4_patch_3_valence_not_bilinear`) and enables convergence on multi-valence domains (`weirdo_multivalence`, `experiment_multivalence_not_bilinear`) even with V-cycles and $\nu=1$.

### Proposal 2 (Best Overall for Direct MG on Stretched/Curved Domains, $\mathcal{O}(N)$ Complexity): Alternating 1D Line-Block GS on Patch Interiors + `Coupled Ring1`
- **Idea**: Replace pointwise Gauss–Seidel on each patch's interior DOFs with **alternating 1D line-strip solves** along $u$-lines and $v$-lines (each strip containing the 2 adjacent B-spline rows/columns corresponding to a multiplicity-2 knot line, solved via 1D banded Cholesky), followed by the `Coupled Ring1` skeleton solve.
- **Why it works & stays $\mathcal{O}(N)$**:
  - Each 1D strip has size $\mathcal{O}(h^{-1})$ and bandwidth $\mathcal{O}(p)$, so factoring all $\mathcal{O}(h^{-1})$ strips costs $\mathcal{O}(h^{-2} p^2) = \mathcal{O}(N)$.
  - Eliminates both the $(|J_u|/|J_v|)^4$ anisotropy bottleneck and the $C^1$ double-knot pairing inside patches.
- **Impact**: On `square_6_patch.xml`, $\kappa(M_{\text{MG}} A)$ drops from $87.15 \to 2.27$ and Direct MG converges in **30 iterations ($0.179\text{ s}$)**; on `square_6_patch_non_bilinear.xml`, $\kappa$ drops from $155.54 \to 6.81$ and Direct MG converges in **80 iterations ($0.415\text{ s}$)**; on `experiment_multivalence_not_bilinear.xml`, $\kappa$ drops from $382.51 \to 25.97$ ($198$ Direct MG iters / $21$ PCG iters).

### Proposal 3 (For Higher Spline Degrees $p \ge 4$): Biharmonic Subspace-Corrected Mass Smoother / Fast Diagonalization on Patch Interiors + `Coupled Ring1`
- **Idea**: For higher spline degrees $p \ge 4$, point and line Gauss–Seidel deteriorate exponentially in $p$ because B-spline support overlaps grow. Because interior basis functions map via the identity ($T_{\text{int}} = I$) and decouple from interfaces, one can apply a **biharmonic tensor-product Fast Diagonalization or Subspace-Corrected Mass Smoother** ([gsPatchPreconditionersCreator.h](file:///home/fhasanova/gismo/src/gsSolver/gsPatchPreconditionersCreator.h#L116-L149)) on each patch interior block, combined multiplicatively with the `Coupled Ring1` skeleton solver.

### Proposal 4 (Low-Effort Algebraic & Code Cleanups — No Smoother Redesign Required)
1. **Use `makeSymmetricGaussSeidelOp` (or $\nu \ge 2\text{--}4$ with W-cycles `-c 2`) when running `-s mg`**:
   - In [as_g1_biharmonic_multigrid_example.cpp](file:///home/fhasanova/gismo/examples/as_g1_biharmonic_multigrid_example.cpp#L495) and [as_g1_biharmonic2_multigrid_example.cpp](file:///home/fhasanova/gismo/examples/as_g1_biharmonic2_multigrid_example.cpp#L523), `makeGaussSeidelOp` executes 1 half-sweep per pre/post step. Switching to `makeSymmetricGaussSeidelOp` (or setting `--pre 2 --post 2` / `-c 2`) cuts Direct MG iterations dramatically (as already seen in [as_g1_biharmonic_multigrid_results.tex](file:///home/fhasanova/gismo/as_g1_biharmonic_multigrid_results.tex#L252-L267), where W-cycles with $\nu \ge 6\text{--}10$ converge on all tested non-convex multi-valence domains).
2. **Tighten sparse drop tolerances in [gsAsG1Basis.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Basis.hpp)**:
   - Change `stride == 2 ? T(1e-4) : T(1e-12)` in `collocateCorners` ([gsAsG1Basis.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Basis.hpp#L287)) and `eps = 1e-4` in `deriveEdgeEmbedding` ([gsAsG1Basis.hpp](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1Basis.hpp#L304)) to `1e-12` so small $\mathcal{O}(h)$ / $\mathcal{O}(h^2)$ derivative coefficients are not truncated on fine grids or curved geometries.
3. **Rescale vertex and crossing functionals by $h$**:
   - Scaling first-derivative functionals by $h$ and second-derivative functionals by $h^2$ inside `collocateCorners` / `collocateBoundaryCrossingDerivative` brings `EdgeCross` and `Corner6Jet` diagonals to $\mathcal{O}(h^{-2})$ matching `Interior` and `EdgeTrace`, eliminating the $10^5\times$ diagonal scale disparity.
4. **Explicit Cholesky coarse solver & smoother loop bounds**:
   - Call `mg->setCoarseSolver(makeSparseCholeskySolver(mg->matrix(0)))` instead of falling back to `SparseLU`, and start the smoother setup loop at `l = 1` instead of `l = 0`.
