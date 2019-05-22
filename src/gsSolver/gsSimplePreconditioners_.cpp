#include <gsSolver/gsSimplePreconditioners.h>
#include <gsSolver/gsSimplePreconditioners.hpp>

namespace gismo
{

namespace internal
{

TEMPLATE_INST void gaussSeidelSweep(const gsSparseMatrix<real_t> & A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);
TEMPLATE_INST void reverseGaussSeidelSweep(const gsSparseMatrix<real_t> & A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f);


TEMPLATE_INST void macroGaussSeidelSweep(const gsSparseMatrix<real_t> & A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, const std::vector< gsSparseMatrix<real_t> >& transfers, const std::vector< gsLinearOperator<real_t>::Ptr >& local_solvers);
TEMPLATE_INST void reverseMacroGaussSeidelSweep(const gsSparseMatrix<real_t> & A, gsMatrix<real_t>& x, const gsMatrix<real_t>& f, const std::vector< gsSparseMatrix<real_t> >& transfers, const std::vector<  gsLinearOperator<real_t>::Ptr >& local_solvers);




} // namespace internal

} // namespace gismo
