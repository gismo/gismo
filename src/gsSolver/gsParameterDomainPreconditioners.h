/** @file gsParameterDomainPreconditioners.h

    @brief Provides preconditioners that live on the parameter domain.

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

/// @brief Provides preconditioners that live on the parameter domain.
///
/// This class provides efficient preconditioners based on the parameter domain,
/// assuming to have \a gsTensorBasis.
///
/// @ingroup Solver
template<typename T>
class gsParameterDomainPreconditioners
{
    typedef typename gsLinearOperator<T>::uPtr OpUPtr;
    typedef typename gsLinearOperator<T>::Ptr  OpPtr;
public:

    /// Constructor taking \a gsBasis, \a gsBoundaryConditions and the
    /// \a dirichtet::strategy
    gsParameterDomainPreconditioners(
        const gsBasis<T>& _basis,
        const gsBoundaryConditions<T>& _bc,
        dirichlet::strategy _dirichlet = dirichlet::elimination
    )
    : m_basis(_basis), m_bc(_bc), m_dirichlet(_dirichlet) { init(); }

    /// Constructor taking \a gsBasis, \a gsBoundaryConditions and
    /// \a gsOptionList object providing the Dirichlet strategy.
    gsParameterDomainPreconditioners(
        const gsBasis<T>& _basis,
        const gsBoundaryConditions<T>& _bc,
        const gsOptionList& _opt
    )
    : m_basis(_basis),
      m_bc(_bc),
      m_dirichlet( (dirichlet::strategy)_opt.askInt("DirichletStrategy", dirichlet::elimination) )
    { init(); }

    /// Assembles mass matrix on the parameter domain
    gsSparseMatrix<T> getMassMatrix()                      const;

    /// Provides \a gsLinearOperator representing the mass matrix (in a matrix-free way)
    OpUPtr            getMassMatrixOp()                    const;

    /// Provides \a gsLinearOperator representing the inverse of the mass matrix (in a matrix-free way)
    OpUPtr            getMassMatrixInvOp()                 const;

    /// Assembles stiffness matrix on the parameter domain
    ///
    /// The stiffness matrix represents \f$ -\Delta u + a u \f$
    gsSparseMatrix<T> getStiffnessMatrix(T a=0)            const;

    /// Provides \a gsLinearOperator representing the stiffness matrix (in a matrix-free way)
    ///
    /// The stiffness matrix represents \f$ -\Delta u + a u \f$
    OpUPtr            getStiffnessMatrixOp(T a=0)          const;

    /// Provides \a gsLinearOperator representing the inverse stiffness matrix
    /// based on the fast diagonalization approach
    /// (SIAM J. Sci. Comput., 38 (6), p. A3644 - A3671, 2016)
    ///
    /// The stiffness matrix represents \f$ -\Delta u + a u \f$
    OpUPtr            getFastDiagonalizationOp(T a=0)      const;

    // Will be provided in a followup pull request:
    // Provides \a gsLinearOperator representing the subspace corrected mass smoother
    // (SIAM J. on Numerical Analysis. 55 (4). p. 2004 - 2024, 2017)
    //
    // This operator is spectrally equivalent to the inverse of
    // \f$ - \Delta u + h^{-2} u \f$
    //OpUPtr            getSubspaceCorrectedMassSmootherOp() const;

    // Helper functions for implementation, which might be of interest also for use in
    // other functions
    static gsSparseMatrix<T>                assembleMass            (const gsBasis<T>& basis);
    static gsSparseMatrix<T>                assembleStiffness       (const gsBasis<T>& basis);
    static std::vector< gsSparseMatrix<T> > assembleTensorMass      (const gsBasis<T>& basis);
    static std::vector< gsSparseMatrix<T> > assembleTensorStiffness (const gsBasis<T>& basis);
    static void handleDirichletConditions(gsSparseMatrix<T>& matrix, const gsBoundaryConditions<T>& bc,
                                          const boxSide& west, const boxSide& east);

private:

    void init();

    const gsBasis<T> &      m_basis;
    gsBoundaryConditions<T> m_bc;
    dirichlet::strategy     m_dirichlet;
};

} // namespace gismo
