/** @file gsIncompleteLUT.h

    @brief Incomplete LU factorization with access to its factors and permutations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo {

/** \brief Incomplete LUT factorization (gsEigen::IncompleteLUT) exposing
    its combined L/U factor matrix and its fill-reducing permutation.

    Use gsEigen::IncompleteLUT directly when only compute()/solve() are
    needed (e.g. as a preconditioner); use this class when the factors
    themselves are needed, e.g. for block smoothers in multigrid.

    \tparam Scalar       coefficient type
    \tparam StorageIndex index type of the factor and the permutations

    \ingroup Matrix
*/
template <class Scalar, class StorageIndex = index_t>
class gsIncompleteLUT : public gsEigen::IncompleteLUT<Scalar, StorageIndex>
{
    typedef gsEigen::IncompleteLUT<Scalar, StorageIndex> Base;
public:
    typedef typename Base::FactorType FactorType;
    typedef gsEigen::PermutationMatrix<gsEigen::Dynamic, gsEigen::Dynamic, StorageIndex> PermutationType;

    gsIncompleteLUT() : Base() {}

    /// Computes the factorization of \a mat with drop tolerance \a droptol and fill factor \a fillfactor
    template <typename MatrixType>
    explicit gsIncompleteLUT(const MatrixType& mat,
                             const typename Base::RealScalar& droptol
                                 = gsEigen::NumTraits<Scalar>::dummy_precision(),
                             int fillfactor = 10)
    : Base(mat, droptol, fillfactor) {}

    /// Returns the combined factor matrix: strictly lower part L (unit diagonal implied), upper part U
    const FactorType& factors() const { return this->m_lu; }

    /// Returns the fill-reducing permutation P
    const PermutationType& fillReducingPermutation() const { return this->m_P; }

    /// Returns the inverse permutation P^{-1}
    const PermutationType& inversePermutation() const { return this->m_Pinv; }
};

} // namespace gismo
