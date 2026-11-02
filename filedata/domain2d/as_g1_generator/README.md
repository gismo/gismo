# AS-G1 Multi-Patch Geometries Catalog

This folder contains planar multi-patch B-spline geometries generated using the **Analysis-Suitable $G^1$ (AS-$G^1$) Multi-Patch Geometry Generator** ([`gsAsG1GeometryGenerator.h`](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1GeometryGenerator.h)), following the Bilinear Skeleton + Linearized Constrained Optimization pipeline (Kapl, Sangalli, Takacs (2018); Hasanova, Takacs, Takacs (2026)).

All geometries satisfy the AS-$G^1$ interface condition to machine precision ($\sim 10^{-15}$ residual) and have strictly regular Jacobian determinants ($\det \nabla G > 0$) across the entire domain.

---

## 1. Geometry Inventory

### A. Symmetric Geometries
Every patch is identical under cyclic rotation around the central extraordinary vertex.

| XML Geometry | Paraview PVD | Valence $V$ | Degree $p$ | Refinements | Shape / Description |
| :--- | :--- | :---: | :---: | :---: | :--- |
| `as_g1_patch3_sym.xml` | `as_g1_patch3_sym.pvd` | 3 | 3 | 0 | Curved equilateral triangle domain |
| `as_g1_patch4_sym.xml` | `as_g1_patch4_sym.pvd` | 4 | 3 | 0 | Curved symmetric 4-patch cross / square domain |
| `as_g1_patch5_sym.xml` | `as_g1_patch5_sym.pvd` | 5 | 3 | 0 | Curved regular pentagonal 5-patch domain |
| `as_g1_patch6_sym.xml` | `as_g1_patch6_sym.pvd` | 6 | 3 | 0 | Curved regular hexagonal 6-patch domain |
| `as_g1_patch7_sym.xml` | `as_g1_patch7_sym.pvd` | 7 | 3 | 0 | Curved regular heptagonal 7-patch domain |
| `as_g1_patch5_star_sym.xml` | `as_g1_patch5_star_sym.pvd` | 5 | 3 | 0 | 5-point star domain ($R_{\text{mid}} = 1.45$) |
| `as_g1_patch5_polygon_sym.xml` | `as_g1_patch5_polygon_sym.pvd` | 5 | 3 | 0 | Polygonal pentagon with straight chord edges |
| `as_g1_ellipse_v5_p3.xml` | `as_g1_ellipse_v5_p3.pvd` | 5 | 3 | 0 | Anisotropic elliptical 5-patch domain ($1.5x \times 0.8y$) |

### B. Randomized / Asymmetric Geometries
Constructed with non-uniform sector angles $\theta_i$, perturbed spoke radii $r_i$, and stochastic boundary curvature perturbations (as in Mathematica prototype `experiments/As-g1-mathematica attempt.txt`).

| XML Geometry | Paraview PVD | Valence $V$ | Degree $p$ | Seed | Max AS-$G^1$ Residual | Min $\det(J)$ |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| `as_g1_patch3_random.xml` | `as_g1_patch3_random.pvd` | 3 | 3 | 303 | $4.87 \times 10^{-15}$ | $0.570$ |
| `as_g1_patch4_random.xml` | `as_g1_patch4_random.pvd` | 4 | 3 | 404 | $8.23 \times 10^{-15}$ | $0.830$ |
| `as_g1_patch5_random.xml` | `as_g1_patch5_random.pvd` | 5 | 3 | 505 | $6.86 \times 10^{-15}$ | $0.626$ |
| `as_g1_patch6_random.xml` | `as_g1_patch6_random.pvd` | 6 | 3 | 606 | $1.84 \times 10^{-14}$ | $0.708$ |
| `as_g1_patch7_random.xml` | `as_g1_patch7_random.pvd` | 7 | 3 | 707 | $6.13 \times 10^{-15}$ | $0.419$ |

---

## 2. Generation CLI Commands

Using the compiled binary [`build/bin/as_g1_geometry_generation_example`](file:///home/fhasanova/gismo/build/bin/as_g1_geometry_generation_example):

### Symmetric Geometries:
```bash
./build/bin/as_g1_geometry_generation_example -v 5 -o filedata/domain2d/as_g1_generator/as_g1_patch5_sym
./build/bin/as_g1_geometry_generation_example -v 5 -s star -o filedata/domain2d/as_g1_generator/as_g1_patch5_star_sym
./build/bin/as_g1_geometry_generation_example -v 6 -o filedata/domain2d/as_g1_generator/as_g1_patch6_sym
```

### Randomized / Asymmetric Geometries:
```bash
./build/bin/as_g1_geometry_generation_example -v 5 -x --seed 123 -o filedata/domain2d/as_g1_generator/as_g1_patch5_random
./build/bin/as_g1_geometry_generation_example -v 6 -x --seed 606 -o filedata/domain2d/as_g1_generator/as_g1_patch6_random
./build/bin/as_g1_geometry_generation_example -v 7 -x --seed 707 --angPerturb 0.15 --radPerturb 0.2 -o filedata/domain2d/as_g1_generator/as_g1_patch7_random
```

---

## 3. PDE Verification (Biharmonic Problem)

Both symmetric and randomized geometries can be immediately used in the existing G+Smo AS-$G^1$ biharmonic solver:

```bash
# Direct solver with h-refinement convergence:
./build/bin/as_g1_biharmonic_example -f filedata/domain2d/as_g1_generator/as_g1_patch5_random.xml -m 2 -r 4

# Multigrid solver (PCG with V-cycles):
./build/bin/as_g1_biharmonic_multigrid_example -f filedata/domain2d/as_g1_generator/as_g1_patch5_random.xml -m 2 -r 4
```

Convergence rates on randomized geometries achieve the theoretical optimum:
- $L^2$ error rate: $\sim 3.8$ to $4.0$ ($O(h^{p+1})$)
- $H^1$ error rate: $\sim 3.5$
- $H^2$ error rate: $\sim 2.1$ ($O(h^p)$)
