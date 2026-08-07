# gsAutoDiff

`gsAutoDiff` integrates [autodiff](https://github.com/autodiff/autodiff) with G+Smo.

## Build

```bash
cmake -S . -B build -DGISMO_OPTIONAL=gsAutoDiff -Dautodiff_DIR=/path/to/autodiff/lib/cmake/autodiff -DCMAKE_CXX_STANDARD=17
make -C build -j1 gsAutoDiff-examples
```

## AD modes

- `GISMO_AUTODIFF_FORWARD=ON` enables `dual_t` and forward-mode examples.
- `GISMO_AUTODIFF_BACKWARD=ON` enables `var_t` and reverse/backward-mode examples.

## Examples

- **Forward AD**: `bSplineBasis_knot_derivative.cpp`, `poisson_rhs_sensitivity.cpp`
- **Reverse AD**: `poisson_geometry_sensitivity.cpp`, `knotVectorVar_example.cpp`, `test_var_basis.cpp`, `test_var_bspline.cpp`
- **Both**: `bSplineSurface_controlpoint_derivative.cpp`, `assemblyAutodiff.cpp`

The example sources use `#ifdef GISMO_AUTODIFF_FORWARD` and `#ifdef GISMO_AUTODIFF_BACKWARD` to keep the mode split explicit.
