/** @file gsPoissonAssembler.h

    @brief Provides assembler for the Poisson equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsAssembler/gsAssemblerBase.h>
#include <gsAssembler/gsAssemblerOptions.h>

#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals

#include <gsPde/gsPoissonPde.h>

//#include <gsAssembler/gsAssemblerUtils.h>

namespace gismo
{

/** @brief
    Implementation of an (multiple right-hand side) Poisson assembler.

    The Poisson equation: \f$-\Delta\mathbf{u}=\mathbf{f} \f$

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).
    
    \ingroup Assembler
*/
template <class T>
class gsPoissonAssembler : public gsAssemblerBase<T>
{
public:
    typedef gsAssemblerBase<T> Base;

public:

    /* @brief
    Main Constructor of the assembler object.

    \param[in] pde A boundary value Poisson problem
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces

    gsPoissonAssembler( gsPoissonPde<T> const         & pde,
                        gsMultiBasis<T> const         & bases,
                        gsDirichletStrategy           dirStrategy,
                        gsInterfaceStrategy           intStrategy)
    :  Base(patches),
       m_rhsFun(&rhs),
       m_bConditions(bconditions),
       m_dirStrategy(dirStrategy),
       m_intStrategy(intStrategy)
    {
        m_bases.push_back(bases);

        const bool conforming = ( m_intStrategy == glue );

        if ( m_dirStrategy == elimination || m_dirStrategy == homogeneous)
            m_dofMapper.push_back( bases.makeMapper(conforming, bconditions) );
        else
            m_dofMapper.push_back( bases.makeMapper(conforming) );

        m_dofs = m_dofMapper.front()->freeSize();

        // Resize system matrix and right hand side
        m_matrix.resize(m_dofs, m_dofs);
        m_rhs.resize(m_dofs, rhs.targetDim() );
    }
//*/


    /** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] basis a multi-basis that contains patch-wise bases
    \param[in] bconditions is a gsBoundaryConditions object that holds all boundary conditions.
    \param[in] rhs is the right-hand side of the Poisson equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
    
    \ingroup Assembler
*/
    gsPoissonAssembler( gsMultiPatch<T> const         & patches,
                        gsMultiBasis<T> const         & basis,
                        gsBoundaryConditions<T> const & bconditions,
                        const gsFunction<T>           & rhs,
                        dirichlet::strategy           dirStrategy = dirichlet::elimination,
                        iFace::strategy               intStrategy = iFace::glue)
        : m_rhsFun(&rhs)
    {
        //set the values, options,...
        m_options.dirStrategy = dirStrategy;
        m_options.intStrategy = intStrategy;
        //m_options.dirValues = dirichlet::l2Projection;
        m_options.dirValues = dirichlet::interpolation;

        //initialize the values, options,...
        this->initialize(patches, basis, bconditions);
    }

    gsPoissonAssembler(const gsFunction<T> & rhs)
    : m_rhsFun(&rhs) //, m_initialized(false)
    { }

    /// Sets the Poisson assembler options
    void setOptions(const gsAssemblerOptions& options);

    /// Look at gsAssemblerBases for documentation
    bool isSymmertric() const { return true; }

    /// Main assembly routine.
    void assemble();

    /// Reconstruct solution field from computed solution vector
    gsField<T> * constructSolution(const gsMatrix<T> & solVector) const;

    /// Reconstruct solution from computed solution vector
    void constructSolution(const gsMatrix<T>& solVector,
                           gsMultiPatch<T>& result) const;

    /// Penalty constant for patch \a k, used for weak bounday conditions
    T penalty(int k) const
    {
        //return gsAssemblerUtils<T>::getMu(m_bases[0][k]);
        const int deg = m_bases[0][k].maxDegree();
        return (deg + m_bases[0][k].dim()) * (deg + 1) * T(2.0);
    }

    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() returns a lower diagonal matrix,
    /// since we exploit symmetry during assembly.
    Eigen::SparseSelfAdjointView< typename gsSparseMatrix<T>::Base, Lower> fullMatrix()
    {
        return m_matrix.template selfadjointView<Lower>();
    }

    /// Computes the Dirichlet DoF values
    void computeDirichletDofs();

protected:

    // initializes the pde specfic stuff
    void initializePdeSpecific()
    {
        applyOptions();
    }

    // Nitsche Dirichlet contributions
    void assembleNitsche();
    
    // Neumann contributions
    void assembleNeumann();

    /// Force Dirichlet boundary conditions by diagonal penalization
    void penalizeDirichlet();

    // Computes the Dirichlet DoF values by interpolation
    void computeDirichletDofsIntpl();

    // Computes the Dirichlet DoF values by L2 projection
    // S.Kleiss
    /** \brief Computes Dirichlet-boundary conditions by L2-projection.
     *
     * ...if the dirichlet::strategy is chosen as dirichlet::elimination.\n
     * A global \f$L_2\f$-projection is applied to the given Dirichlet data
     * and the eliminated coefficients are set to the corresponding values.
     * The projection is global in the sense that all Dirichlet-DOFs are
     * computed at once.
     *
     * \ingroup Assembler
     */
    void computeDirichletDofsL2Proj();

    // Applies the Poisson assembler options
    void applyOptions();

protected:

    /// Right hand side function
    const gsFunction<T> * m_rhsFun;

protected:

    // Members from gsAssemblerBase
    using gsAssemblerBase<T>::m_patches;
    using gsAssemblerBase<T>::m_bases;
    using gsAssemblerBase<T>::m_bConditions;
    using gsAssemblerBase<T>::m_dofMappers;
    using gsAssemblerBase<T>::m_ddof;
    using gsAssemblerBase<T>::m_options;
    using gsAssemblerBase<T>::m_matrix;
    using gsAssemblerBase<T>::m_rhs;
    using gsAssemblerBase<T>::m_dofs;
};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPoissonAssembler.hpp)
#endif
