# Structured AS-G1 2D Multi-Patch Benchmark Suite (101 Unique Validated Domains)

Organizes all 90 unique domains from `--suite all` (99 entries including 9 byte-identical generator aliases) plus 11 additional validated AS-G1 domains into 5 non-overlapping topological classes.

## Categories

- `01_single_patch/` (10): Single-patch domains ($P=1$, no interfaces, $n_{\\text{holes}}=0$).
- `02_multipatch_trees/` (22): Simply-connected multi-patch trees ($P \\ge 2$, $n_{\\text{intV}}=0$, $n_{\\text{holes}}=0$).
- `03_single_vertex_stars/` (34): Single interior extraordinary vertex ($n_{\\text{intV}}=1$, valence $3\\text{--}7$, $n_{\\text{holes}}=0$).
- `04_multi_vertex_meshes/` (23): Multiple coupled interior vertices ($n_{\\text{intV}} \\ge 2$, $n_{\\text{holes}}=0$).
- `05_multiply_connected_holes/` (12): Multiply-connected domains with $1\\text{--}6$ topological boundary loops ($n_{\\text{holes}} \\ge 1$).

See `MANIFEST.csv` for full metadata per domain (`id`, `category`, `filename`, `orig_suite_tag`, `orig_rel_path`, `n_patches`, `n_interfaces`, `geom_degree`, `bilinear_param`, `n_holes`, `n_int_vertices`, `max_int_valence`, `int_valences`, `max_ctrl_aspect_ratio`).

## Byte-Identical Aliases Excluded from Copying (from original 99-suite)

- `domain2d/as_g1_generator/as_g1_v3_p3.xml` == `domain2d/as_g1_generator/as_g1_tri_p3.xml`
- `domain2d/as_g1_generator/as_g1_v5_p3_r1.xml` == `domain2d/as_g1_generator/as_g1_pentagon_p3_r1.xml`
- `domain2d/as_g1_generator/as_g1_v6_p3.xml` == `domain2d/as_g1_generator/as_g1_hex_p3.xml`
- `domain2d/as_g1_generator/as_g1_patch3_sym.xml` == `domain2d/as_g1_generator/as_g1_tri_p3.xml`
- `domain2d/as_g1_generator/as_g1_patch4_sym.xml` == `domain2d/as_g1_generator/as_g1_v4_p3.xml`
- `domain2d/as_g1_generator/as_g1_patch5_sym.xml` == `domain2d/as_g1_generator/as_g1_pentagon_p3.xml`
- `domain2d/as_g1_generator/as_g1_patch6_sym.xml` == `domain2d/as_g1_generator/as_g1_hex_p3.xml`
- `domain2d/as_g1_generator/as_g1_patch7_sym.xml` == `domain2d/as_g1_generator/as_g1_v7_p3.xml`
- `domain2d/as_g1_generator/as_g1_star_v5_p3.xml` == `domain2d/as_g1_generator/as_g1_patch5_star_sym.xml`

## Defective Non-Suite-All Extras Excluded (Direct Cholesky Ill-Posed)

- `domain2d/as_g1_domain_decomposition_results/asg1_reparametrized/heart_1corner_nonconvex_asg1.xml` (`NaN` in direct Cholesky)
- `domain2d/as_g1_generator/square_with_disk_asg1.xml` (direct Cholesky relative residual mismatch `1.0`)
- `domain2d/as_g1_no_holes_reparametrized/1to7corners/star_pentagon_5corner_nonconvex_asg1.xml` (divergent `2.79e20` / singular parameterization)
