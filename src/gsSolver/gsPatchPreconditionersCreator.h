/** @file gsPatchPreconditionersCreator.h

    @brief Provides robust preconditioners for single patch geometries.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs, C. Hofreither
*/
#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsSolver/gsLinearOperator.h>


namespace gismo
{

/// @brief Provides robust preconditioners for single patch geometries.
///
/// This class provides efficient preconditioners for single patch geometries,
/// assuming to have \a gsTensorBasis.
///
/// @ingroup Solver
template<typename T>
class gsPatchPreconditionersCreator
{
    typedef typename gsLinearOperator<T>::uPtr OpUPtr;
    typedef typename gsLinearOperator<T>::Ptr  OpPtr;
public:

    /// Provieds the mass matrix on the parameter domain
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    static gsSparseMatrix<T> massMatrix(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions()
    );

    /// Provides \a gsLinearOperator representing the mass matrix (in a matrix-free way)
    /// on the parameter domain
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    static OpUPtr            massMatrixOp(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions()
    );

    /// Provides \a gsLinearOperator representing the inverse of the mass matrix (in a matrix-free way)
    /// on the parameter domain
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    static OpUPtr            massMatrixInvOp(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions()
    );

    /// Provieds stiffness matrix on the parameter domain
    ///
    /// The stiffness matrix represents \f$ -\Delta u + \alpha u \f$
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    /// \param alpha  Scaling parameter (see above)
    static gsSparseMatrix<T> stiffnessMatrix(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions(),
        T alpha = 0
    );

    /// Provides \a gsLinearOperator representing the stiffness matrix (in a matrix-free way)
    /// on the parameter domain
    ///
    /// The stiffness matrix represents \f$ -\Delta u + \alpha u \f$
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    /// \param alpha  Scaling parameter (see above)
    static OpUPtr            stiffnessMatrixOp(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions(),
        T alpha = 0
    );

    /// Provides \a gsLinearOperator representing the inverse stiffness matrix
    /// on the parameter domain based on the fast diagonalization approach
    /// (SIAM J. Sci. Comput., 38 (6), p. A3644 - A3671, 2016)
    ///
    /// The stiffness matrix represents \f$ -\Delta u + \alpha u \f$
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    /// \param alpha  Scaling parameter (see above)
    static OpUPtr            fastDiagonalizationOp(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions(),
        T alpha = 0
    );

    /// Provides \a gsLinearOperator representing the subspace corrected mass smoother
    /// on the parameter domain (SIAM J. on Numerical Analysis. 55 (4). p. 2004 - 2024, 2017)
    ///
    /// This operator is spectrally equivalent to the inverse of
    /// \f$ - \Delta u + \sigma h^{-2} u \f$
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    /// \param sigma  Scaling parameter (see above)
    static OpUPtr            subspaceCorrectedMassSmootherOp(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions(),
        T sigma = T(12)/T(100)
    );

};

} // namespace gismo
