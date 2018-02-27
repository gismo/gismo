#include <gsTensor/gsTensorTools.hpp>

namespace gismo
{

TEMPLATE_INST
gsSparseMatrix<real_t> getKroneckerProduct(const gsSparseMatrix<real_t>& A, const gsSparseMatrix<real_t>& B);

TEMPLATE_INST
gsSparseMatrix<real_t, RowMajor> getKroneckerProduct(const gsSparseMatrix<real_t, RowMajor>& A, const gsSparseMatrix<real_t, RowMajor>& B);

TEMPLATE_INST
gsSparseMatrix<std::complex<real_t> > getKroneckerProduct(const gsSparseMatrix<std::complex<real_t> >& A, const gsSparseMatrix<std::complex<real_t> >& B);

TEMPLATE_INST
gsSparseMatrix<std::complex<real_t>, RowMajor> getKroneckerProduct(const gsSparseMatrix<std::complex<real_t>, RowMajor>& A, const gsSparseMatrix<std::complex<real_t>, RowMajor>& B);


TEMPLATE_INST
gsMatrix<real_t> getKroneckerProduct(const gsMatrix<real_t>& A, const gsMatrix<real_t>& B);

TEMPLATE_INST
gsMatrix<std::complex<real_t> > getKroneckerProduct(const gsMatrix<std::complex<real_t> >& A, const gsMatrix<std::complex<real_t> >& B);


}
