#include <gsSolver/gsSimplePreconditioners.h>
#include <gsSolver/gsSimplePreconditioners.hpp>

namespace gismo
{

namespace internal
{

TEMPLATE_INST void gaussSeidelSweep(const gsSparseMatrix<real_t> & A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);
TEMPLATE_INST void reverseGaussSeidelSweep(const gsSparseMatrix<real_t> & A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);

} // namespace internal

} // namespace gismo
