/** @file gsGridHierarchy_.cpp

    @brief Coarsening algorithms for knot vectors and bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#include <gsMultiGrid/gsGridHierarchy.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsGridHierarchy<real_t>;

TEMPLATE_INST
void uniformRefine_withTransfer(
    const gsMultiBasis<real_t>& mBasis,
    const gsBoundaryConditions<real_t>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    index_t numberOfKnotsToBeInserted,
    index_t multiplicityOfKnotsToBeInserted,
    gsMultiBasis<real_t> & refinedMBasis,
    gsSparseMatrix<real_t, RowMajor>& transferMatrix,
    std::vector< gsSparseMatrix<real_t, RowMajor> >& localTransferMatrices
    );

TEMPLATE_INST
gsKnotVector<real_t> coarsenKnotVector(const gsKnotVector<real_t>& kv, std::vector<real_t>& removedKnots);

TEMPLATE_INST
gsTensorBSplineBasis<1, real_t>::uPtr coarsenTensorBasis(const gsTensorBSplineBasis<1, real_t>& b, std::vector< std::vector<real_t> > & removedKnots);

TEMPLATE_INST
gsTensorBSplineBasis<2, real_t>::uPtr coarsenTensorBasis(const gsTensorBSplineBasis<2, real_t>& b, std::vector< std::vector<real_t> > & removedKnots);

TEMPLATE_INST
gsTensorBSplineBasis<3, real_t>::uPtr coarsenTensorBasis(const gsTensorBSplineBasis<3, real_t>& b, std::vector< std::vector<real_t> > & removedKnots);

TEMPLATE_INST
gsTensorBSplineBasis<4, real_t>::uPtr coarsenTensorBasis(const gsTensorBSplineBasis<4, real_t>& b, std::vector< std::vector<real_t> > & removedKnots);

TEMPLATE_INST
gsBasis<real_t>::uPtr coarsenBasis(const gsBasis<real_t> & b, std::vector< std::vector<real_t> >& removedKnots);

TEMPLATE_INST
gsBasis<real_t>::uPtr coarsenBasis_withTransfer(const gsBasis<real_t>& b, gsSparseMatrix<real_t, RowMajor>& transferMatrix);

TEMPLATE_INST
void coarsenMultiBasis_withTransfer(
    const gsMultiBasis<real_t>& mBasis,
    const gsBoundaryConditions<real_t>& boundaryConditions,
    const gsOptionList& assemblerOptions,
    gsMultiBasis<real_t> & refinedMBasis,
    gsSparseMatrix<real_t, RowMajor>& transferMatrix,
    std::vector< gsSparseMatrix<real_t, RowMajor> >& localTransferMatrices
    );

TEMPLATE_INST
void combineTransferMatrices(
    const std::vector< gsSparseMatrix<real_t, RowMajor> >& localTransferMatrices,
    const gsDofMapper& coarseMapper,
    const gsDofMapper& fineMapper,
    gsSparseMatrix<real_t, RowMajor>& transferMatrix
    );

TEMPLATE_INST
void copyIfNotIndexed(const std::vector<real_t>& data, const std::vector<index_t>& indices, std::vector<real_t>& result);

} // namespace gismo
