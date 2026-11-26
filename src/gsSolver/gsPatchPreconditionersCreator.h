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
    /// The stiffness matrix represents \f$ - \beta \Delta u + \alpha u \f$
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    /// \param alpha  Scaling parameter (see above)
    static gsSparseMatrix<T> stiffnessMatrix(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions(),
        T alpha = 0,
        T beta = 1
    );

    /// Provides \a gsLinearOperator representing the stiffness matrix (in a matrix-free way)
    /// on the parameter domain
    ///
    /// The stiffness matrix represents \f$ - \beta \Delta u + \alpha u \f$
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    /// \param alpha  Scaling parameter (see above)
    static OpUPtr            stiffnessMatrixOp(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions(),
        T alpha = 0,
        T beta = 1
    );

    /// Provides \a gsLinearOperator representing the inverse stiffness matrix
    /// on the parameter domain based on the fast diagonalization approach
    /// (SIAM J. Sci. Comput., 38 (6), p. A3644 - A3671, 2016)
    ///
    /// The stiffness matrix represents
    /// \f$ \beta (\nabla u, \nabla v)_{L_2} + \alpha (u, v)_{L_2} + \gamma (u, 1)_{L_2} (v,1)_{L_2} \f$
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    /// \param alpha  Scaling parameter (see above)
    /// \param beta   Scaling parameter (see above)
    /// \param gamma  Scaling parameter (see above). Only allowed for pure Neumann case.
    static OpUPtr            fastDiagonalizationOp(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions(),
        T alpha = 0,
        T beta = 1,
        T gamma = 0
    );

    /// Provides \a gsLinearOperator representing the inverse of an anisotropic
    /// stiffness matrix on the parameter domain based on the fast diagonalization
    /// approach (SIAM J. Sci. Comput., 38 (6), p. A3644 - A3671, 2016)
    ///
    /// The stiffness matrix represents
    /// \f$ \sum_{k=1}^d \beta_k (\partial_k u, \partial_k v)_{L_2} + \alpha (u, v)_{L_2} \f$
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    /// \param alpha  Scaling parameter (see above)
    /// \param beta   Vector with one positive scaling parameter per parametric direction
    static OpUPtr            fastDiagonalizationOp(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc,
        const gsOptionList& opt,
        T alpha,
        const gsVector<T>& beta
    );

    /// Provides \a gsLinearOperator representing a preconditioner for the
    /// stiffness matrix \f$ (\nabla u, \nabla v)_{L_2} \f$ on the physical
    /// domain of a single patch, based on the fast diagonalization approach
    /// (SIAM J. Sci. Comput., 38 (6), p. A3644 - A3671, 2016)
    ///
    /// The operator is the exact inverse of
    /// \f$ \sum_{k=1}^d (\tilde c_k \partial_k u, \partial_k v)_{L_2} \f$
    /// on the parameter domain, where
    /// \f$ \tilde c_k(\xi) = a_k(\xi_k) \prod_{j \ne k} m_j(\xi_j) \f$
    /// is a separable approximation of
    /// \f$ c_k = (|\det J| J^{-1} J^{-\top})_{kk} \f$ and \f$ J \f$ is the
    /// Jacobian of the geometry map. The weights \f$ a_k \f$ and
    /// \f$ m_j \f$ are those of \a elasticityFastDiagonalizationOp. If
    /// \f$ |\det J| J^{-1} J^{-\top} \f$ is diagonal and its entries are of
    /// the given form, the operator is the exact inverse of the stiffness
    /// matrix that is assembled with the quadrature rule given by \a opt.
    ///
    /// \param basis  A tensor basis
    /// \param geo    Geometry map of the patch
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    static OpUPtr            fastDiagonalizationOp(
        const gsBasis<T>& basis,
        const gsGeometry<T>& geo,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions()
    );

    /// Provides \a gsLinearOperator representing a preconditioner for the
    /// displacement formulation of linear elasticity on a single patch, based
    /// on the fast diagonalization approach
    /// (SIAM J. Sci. Comput., 38 (6), p. A3644 - A3671, 2016)
    ///
    /// The operator is the exact inverse of the block diagonal of the
    /// elasticity bilinear form
    /// \f$ 2 \mu (\varepsilon(u), \varepsilon(v))_{L_2} + \lambda (\nabla\cdot u, \nabla\cdot v)_{L_2} \f$
    /// on the parameter domain: the block of the displacement component \f$ c \f$
    /// represents
    /// \f$ \sum_{k=1}^d (\mu + (\lambda+\mu) \delta_{ck}) (\partial_k u_c, \partial_k v_c)_{L_2} \f$.
    /// By Korn's inequality, it is spectrally equivalent to the inverse of the
    /// elasticity stiffness matrix, with constants independent of the grid
    /// size and the spline degree, as long as each component has a Dirichlet
    /// condition on at least one side.
    ///
    /// The unknowns are ordered component by component, as for a vector
    /// valued space of \a gsExprAssembler. The boundary conditions of the
    /// component \f$ c \f$ are those registered for the unknown 0 with the
    /// component \f$ c \f$ or with all components.
    ///
    /// If a geometry is given, the block of the component \f$ c \f$ represents
    /// \f$ \sum_{k=1}^d (\tilde c_{c,k} \partial_k u_c, \partial_k v_c)_{L_2} \f$
    /// instead, where \f$ \tilde c_{c,k} \f$ is a separable approximation of
    /// \f$ c_{c,k} = \mu (|\det J| J^{-1} J^{-\top})_{kk} + (\lambda+\mu) |\det J| (J^{-1})_{kc}^2 \f$
    /// and \f$ J \f$ is the Jacobian of the geometry map. These coefficients
    /// are the ones of the terms \f$ \partial_k u_c \partial_k v_c \f$ of the
    /// diagonal block of the elasticity bilinear form on the physical
    /// domain; its terms with mixed derivatives are dropped. The
    /// approximation reads
    /// \f$ \tilde c_{c,k}(\xi) = a_k(\xi_k) \prod_{j \ne k} m_j(\xi_j) \f$,
    /// so the block is the sum over \f$ k \f$ of Kronecker products of a
    /// univariate stiffness matrix weighted with \f$ a_k \f$ and univariate
    /// mass matrices weighted with \f$ m_j \f$, which the fast diagonalization
    /// inverts exactly. The weights are means over the quadrature nodes of
    /// the basis: \f$ a_k(\xi_k) \f$ is the mean of \f$ c_{c,k} \f$ over all other
    /// directions and \f$ m_j(\xi_j) \f$ is the geometric mean over
    /// \f$ k \ne j \f$ of the mean of \f$ c_{c,k} \f$ over all directions but
    /// \f$ j \f$, divided by the mean of \f$ c_{c,k} \f$ over the patch. The
    /// approximation is exact if the coefficients are of the given form,
    /// like for a box with a different grading in each direction or for an
    /// annulus in two dimensions. For a constant Jacobian, each direction
    /// is scaled with one number, which accounts for the size and the
    /// aspect ratio of the patch. The setup evaluates the Jacobian once per
    /// quadrature node of the patch.
    ///
    /// If the determinant of the Jacobian vanishes or changes its sign at a
    /// quadrature node, each direction is scaled with the average of
    /// \f$ (|\det J| J^{-1} J^{-\top})_{kk} \f$ over the anchors of the
    /// geometry instead.
    ///
    /// \param basis   A tensor basis (of one displacement component)
    /// \param bc      Boundary conditions
    /// \param opt     Assembler options
    /// \param lambda  First Lame parameter
    /// \param mu      Second Lame parameter (shear modulus)
    /// \param geo     Geometry map of the patch (optional)
    static OpUPtr            elasticityFastDiagonalizationOp(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc,
        const gsOptionList& opt,
        T lambda,
        T mu,
        const gsGeometry<T>* geo = NULL
    );

    /// Provides \a gsLinearOperator representing a preconditioner for the
    /// displacement formulation of linear elasticity on a single patch, based
    /// on the fast diagonalization approach
    ///
    /// As above, with one set of scalar boundary conditions per displacement
    /// component. Use this version with assemblers that register the
    /// components as separate unknowns.
    ///
    /// \param basis   A tensor basis (of one displacement component)
    /// \param bc      Boundary conditions for each component
    /// \param opt     Assembler options
    /// \param lambda  First Lame parameter
    /// \param mu      Second Lame parameter (shear modulus)
    /// \param geo     Geometry map of the patch (optional)
    static OpUPtr            elasticityFastDiagonalizationOp(
        const gsBasis<T>& basis,
        const std::vector< gsBoundaryConditions<T> >& bc,
        const gsOptionList& opt,
        T lambda,
        T mu,
        const gsGeometry<T>* geo = NULL
    );

    /// Provides \a gsLinearOperator representing the subspace corrected mass smoother
    /// on the parameter domain (SIAM J. on Numerical Analysis. 55 (4). p. 2004 - 2024, 2017)
    ///
    /// This operator is spectrally equivalent to the inverse of
    ///
    /// \f$ \beta ( - \Delta u + \sigma h^{-2} u) + \alpha u \f$
    ///
    /// assuming \f$ \sigma = \mathcal{O}(1) \f$ to be large enough; the exact meaning of
    /// \f$ \sigma \f$ is explained in the abovementioned paper (\f$ \sigma \f$ from the
    /// paper equals \f$  1/(\sigma*h*h) \f$ here.)
    ///
    /// \param basis  A tensor basis
    /// \param bc     Boundary conditions
    /// \param opt    Assembler options
    /// \param sigma  Scaling parameter (see above)
    /// \param alpha  Scaling parameter (see above)
    /// \param beta   Scaling parameter (see above)
    static OpUPtr            subspaceCorrectedMassSmootherOp(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions(),
        T sigma = (T)(12) / (T)(100),
        T alpha = 0,
        T beta = 1
    );

    /// Provides matrices that represent the basis of the space \f$ \widetilde{S}_{p,h} \f$,
    /// as introduced in M^3AS. 26 (7), p. 1411 - 1445, 2016, and its \f$\ell^2\f$ - orthogonal
    /// complement.
    ///
    /// \param  basis  A tensor basis
    /// \param  bc     Boundary conditions
    /// \param  opt    Assembler options
    /// \return \a std::pair containing \a std::vector of d sparse matrices (for each dimension) that represent
    ///         the basis transformation for (a) \f$ \widetilde{S}_{p,h} \f$, and (b) its complement.
    static std::pair< std::vector< gsSparseMatrix<T> >, std::vector< gsSparseMatrix<T> > > getTildeSpaceBasisTransformation(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc = gsBoundaryConditions<T>(),
        const gsOptionList& opt = gsAssembler<T>::defaultOptions()
    );

private:

    // Fast diagonalization with one coefficient per parametric direction
    static OpUPtr            fastDiagonalization(
        const gsBasis<T>& basis,
        const gsBoundaryConditions<T>& bc,
        const gsOptionList& opt,
        T alpha,
        const gsVector<T>& beta,
        T gamma
    );

    // Fast diagonalization for given univariate matrices, where the i-th
    // matrix belongs to the direction d-1-i
    static OpUPtr            fastDiagonalization(
        const std::vector< gsSparseMatrix<T> >& local_stiff,
        const std::vector< gsSparseMatrix<T> >& local_mass,
        T alpha,
        T gamma
    );

};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPatchPreconditionersCreator.hpp)
#endif
