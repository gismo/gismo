#include <gsCore/gsTemplateTools.h>

#include <gsUtils/gsQuadrature.h>
#include <gsUtils/gsQuadrature.hpp>

#define T real_t
#define uZ unsigned
#define Z int

namespace gismo
{

///////////////////////////////////////////////////////////////////
////// Quadrature
///////////////////////////////////////////////////////////////////

TEMPLATE_INST
void iteratedGaussRule<T>(gsMatrix<T>& ngrid, gsVector<T>& wgrid, int n, const std::vector<T>& intervals);

TEMPLATE_INST
void tensorGaussRule<T>(gsMatrix<T>& ngrid, gsVector<T>& wgrid, const gsVector<int>& numNodes, const gsVector<T>& lower, const gsVector<T>& upper);

TEMPLATE_INST
void tensorIteratedGaussRule<T>(gsMatrix<T>& ngrid, gsVector<T>& wgrid, const gsVector<int>& numNodes, const std::vector< std::vector<T> >& intervals);

TEMPLATE_INST
void uniformGaussRule<T>(gsMatrix<T>& ngrid, gsVector<T>& wgrid, int numEvals, int degree, const gsVector<T>& lower, const gsVector<T>& upper);

TEMPLATE_INST
void gsRefGaussRule<T> (gsVector<T>& x, gsVector<T>& w, int n, T a, T b);

} // namespace gismo

#undef T
#undef uZ
#undef Z
