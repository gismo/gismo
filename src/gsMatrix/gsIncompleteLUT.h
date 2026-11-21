/** @file gsIncompleteLUT.h

    @brief Subclass of Eigen::IncompleteLUT exposing internal members.

    ## Why this wrapper exists

    In Eigen 3.x the factor matrix and fill-reducing permutation were public
    members (`m_lu`, `m_P`, `m_Pinv`) of `Eigen::IncompleteLUT`, and G+Smo
    previously patched Eigen directly to add named accessors
    (`factors()`, `fillReducingPermutation()`, `inversePermutation()`).

    Eigen 5 moved those members to `protected` and did **not** introduce
    equivalent public accessor methods.  Code that only uses
    `Eigen::IncompleteLUT` as a black-box preconditioner (calling `compute()`
    and `solve()`) can continue to use `Eigen::IncompleteLUT` directly.
    Code that needs the explicit L/U factor matrix or the permutation
    matrices — e.g. for custom smoothers or multigrid setups — must use
    this wrapper, which re-exposes them through the three inline accessors
    below without requiring any patch to Eigen.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo {

/// \brief Eigen::IncompleteLUT subclass that exposes the ILUT factors
/// and fill-reducing permutation. Drop-in replacement — add any
/// missing forwarding constructors here if needed.
template <class Scalar, class StorageIndex = index_t>
class gsIncompleteLUT : public Eigen::IncompleteLUT<Scalar, StorageIndex>
{
    typedef Eigen::IncompleteLUT<Scalar, StorageIndex> Base;

public:
    typedef typename Base::FactorType FactorType;
    typedef Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, StorageIndex>
            PermutationType;

    gsIncompleteLUT() : Base() {}

    template <typename MatrixType>
    explicit gsIncompleteLUT(const MatrixType& mat,
                             const typename Base::RealScalar& droptol
                                 = Eigen::NumTraits<Scalar>::dummy_precision(),
                             int fillfactor = 10)
    : Base(mat, droptol, fillfactor) {}

    /// Returns the combined L/U factor matrix (L in strict lower, U in upper).
    const FactorType& factors() const { return this->m_lu; }

    /// Returns the fill-reducing permutation P.
    const PermutationType& fillReducingPermutation() const { return this->m_P; }

    /// Returns the inverse permutation P^{-1}.
    const PermutationType& inversePermutation() const { return this->m_Pinv; }
};

} // namespace gismo
