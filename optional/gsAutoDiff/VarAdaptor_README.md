# VarAdaptor: Reverse-Mode Automatic Differentiation for GISMO

## Overview

`VarAdaptor` is a value-semantic wrapper around `autodiff::var` that provides reverse-mode automatic differentiation compatible with GISMO's template system.

## Why VarAdaptor?

The `autodiff` library's reverse-mode type (`autodiff::var`) uses expression templates where operations return `ExprPtr` (expression pointers) rather than `var` types. This breaks GISMO's template system which expects consistent value semantics.

`VarAdaptor` solves this by:
- Wrapping `autodiff::var` with value semantics
- Ensuring all operations return `VarAdaptor` types
- Providing transparent conversion to/from underlying `var` objects
- Implementing all necessary arithmetic and mathematical operations

## When to Use Reverse vs Forward Mode

### Forward-Mode (dual)
- **Efficient when**: #inputs << #outputs
- **Use for**: Computing Jacobians of vector functions
- **Example**: B-spline evaluation (few parameters, many spatial points)
- **Cost**: One forward pass per input variable

### Reverse-Mode (var/VarAdaptor)
- **Efficient when**: #inputs >> #outputs  
- **Use for**: Computing gradients of scalar functions
- **Example**: Optimization objectives, loss functions
- **Cost**: One reverse pass per output variable
- **Note**: This is the mode used in machine learning (backpropagation)

## Usage Example

```cpp
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsAutoDiff/gsVarAdaptor.h>

using namespace gismo;

// Create independent variables
VarAdaptorD x(1.0);
VarAdaptorD y(2.0);
VarAdaptorD z(3.0);

// Compute function
VarAdaptorD f = x*x + 2.0*y*z + sin(x*y);

// Get the value
double f_val = autodiff::val(f);

// Compute all gradients in one reverse sweep
auto [fx, fy, fz] = autodiff::derivatives(f, autodiff::wrt(x.get(), y.get(), z.get()));
```

## Supported Operations

### Arithmetic
- Binary: `+`, `-`, `*`, `/`
- Unary: `+x`, `-x`
- Compound: `+=`, `-=`, `*=`, `/=`

### Comparison
- `==`, `!=`, `<`, `<=`, `>`, `>=`
- Works with both `VarAdaptor` and scalars

### Mathematical Functions
- Trigonometric: `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`
- Hyperbolic: `sinh`, `cosh`, `tanh`
- Exponential/Logarithmic: `exp`, `log`, `log10`
- Power/Root: `pow`, `sqrt`
- Other: `abs`, `erf`

### Non-Differentiable Operations
- `floor`, `ceil`, `round`, `trunc` (return constant values)
- `min`, `max` (return one of the inputs)

### Utility Functions
- `isinf`, `isnan`, `isfinite`
- Stream operators: `<<`, `>>`

## Implementation Notes

1. **Not Precompiled**: `VarAdaptor` is header-only and not precompiled into the GISMO library. This avoids template instantiation issues.

2. **Access to Underlying var**: Use `.get()` to access the wrapped `autodiff::var` when needed for gradient computations:
   ```cpp
   VarAdaptorD x(1.0);
   autodiff::var& x_var = x.get();
   ```

3. **Value Extraction**: Use `autodiff::val()` or `.val()`:
   ```cpp
   double x_value = autodiff::val(x);  // or x.val()
   ```

4. **Gradient Computation**: The `autodiff` library's gradient functions work with the underlying `var` objects:
   ```cpp
   auto grad = autodiff::gradient(f.get(), {x.get(), y.get()});
   ```

## Limitations

- VarAdaptor works at the C++ level only, not available in the precompiled library
- For optimal performance with many parameters, reverse-mode is preferred
- Expression trees can grow large for complex functions

## See Also

- `helloAutoDiff.cpp`: Forward-mode examples
- `linearSolveAutoDiff.cpp`: AD through linear solvers  
- `reverseAutoDiff.cpp`: Reverse-mode examples (if created)
- [autodiff documentation](https://autodiff.github.io/)
