// Newton-Raphson method of finding roots

#include <math.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunction.h>

namespace gismo {

/// @brief Newton-Raphson method find a solution of the equation f(x) = y with starting vector x
///
/// \ingroup Utils
// function result: number of iterations >= 0 if root found, -1 if max_loop exceeded  
template <class T>
int newtonRaphson( const gsFunction<T> & f, gsVector<T> &x, const gsVector<T> &y, 
                   bool withSupport, const T accuracy = 1e-6, int max_loop = 100 )
{
    gsMatrix<T> supp;
    if (withSupport)
        supp = f.support();

    gsMatrix<T> delta, jac;
    int iter = 0;
    do {
        // compute residual y - f(x)
        f.eval_into (x, delta);
        delta = y - delta;

        // compute Jacobian and solve for next x
        f.deriv_into(x, jac);
        delta = jac.partialPivLu().solve( delta );
        x += delta;

        // clamp x to the support of the function
        if (withSupport)
            x = x.cwiseMax( supp.col(0) ).cwiseMin( supp.col(1) );

        if (delta.norm() <= accuracy)
            return iter;
    } while (++iter < max_loop);

    // no solution found within max_loop iterations
    return -1;
}


}; // namespace gismo
