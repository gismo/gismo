# AS-G1 Multi-Patch Geometries Catalog

This folder contains planar multi-patch B-spline geometries generated using the **Analysis-Suitable $G^1$ (AS-$G^1$) Multi-Patch Geometry Generator** ([`gsAsG1GeometryGenerator.h`](file:///home/fhasanova/gismo/src/gsModeling/gsAsG1GeometryGenerator.h)), following the Bilinear Skeleton + Linearized Constrained Optimization pipeline (Kapl, Sangalli, Takacs (2018); Hasanova, Takacs, Takacs (2026)).

All geometries satisfy the AS-$G^1$ interface condition to machine precision ($\sim 10^{-15}$ residual) and have strictly regular Jacobian determinants ($\det \nabla G > 0$) across the entire domain.

> [!NOTE]
> For a full mathematical derivation, literature background, step-by-step algorithms, and C++ software architecture, see the detailed reference document:
> **[`ALGORITHM.md`](file:///home/fhasanova/gismo/filedata/domain2d/as_g1_generator/ALGORITHM.md)**.

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
- $L^2$ error rate: $\sim 3.8$ to $4.7$ ($O(h^{p+1}) = O(h^4)$)
- $H^1$ error rate: $\sim 3.1$ to $4.0$
- $H^2$ error rate: $\sim 2.0$ to $3.0$ ($O(h^p) = O(h^2)$)

### Benchmark Summary Across All Geometries

| Geometry | Type | $L^2$ Rate ($r=3 \to 4$) | $H^1$ Rate | $H^2$ Rate | Multigrid V-cycles | Multigrid Solve Time |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `as_g1_patch3_sym.xml` | Triangle (V=3) | **4.68** | 3.44 | 2.16 | 20 | 0.065 s |
| `as_g1_patch3_random.xml` | Triangle (V=3) | **4.28** | 3.16 | 2.06 | 39 | 0.150 s |
| `as_g1_patch4_sym.xml` | Quad (V=4) | **4.41** | 3.89 | 2.55 | 20 | 0.085 s |
| `as_g1_patch4_random.xml` | Quad (V=4) | **4.50** | 3.49 | 2.14 | 25 | 0.115 s |
| `as_g1_patch5_sym.xml` | Pentagon (V=5) | **3.94** | 3.99 | 2.94 | 24 | 0.108 s |
| `as_g1_patch5_random.xml` | Pentagon (V=5) | **3.78** | 3.23 | 2.04 | 28 | 0.109 s |
| `as_g1_patch5_star_sym.xml` | Star (V=5) | **3.94** | 3.76 | 2.19 | 24 | 0.139 s |
| `as_g1_patch5_polygon_sym.xml` | Polygon (V=5) | **4.09** | 3.92 | 2.80 | 26 | 0.099 s |
| `as_g1_patch6_sym.xml` | Hexagon (V=6) | **4.18** | 4.04 | 2.99 | 34 | 0.181 s |
| `as_g1_patch6_random.xml` | Hexagon (V=6) | **3.75** | 3.20 | 2.06 | 37 | 0.155 s |
| `as_g1_patch7_sym.xml` | Heptagon (V=7) | **4.34** | 4.00 | 2.51 | 39 | 0.255 s |
| `as_g1_patch7_random.xml` | Heptagon (V=7) | **3.30** | 3.06 | 2.04 | 42 | 0.305 s |

