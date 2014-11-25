#include <gsCore/gsTemplateTools.h>

#include <gsUtils/gsCollocationMatrix.h>
#include <gsUtils/gsCollocationMatrix.hpp>

namespace gismo
{

TEMPLATE_INST
void gsCollocationMatrix_into (gsBasis<real_t> const& basis, gsMatrix<real_t> const& u, 
		       gsSparseMatrix<real_t> & res);


TEMPLATE_INST
gsSparseMatrix<real_t> * gsCollocationMatrix (gsBasis<real_t> const& basis, gsMatrix<real_t> const& u);

} // namespace gismo

