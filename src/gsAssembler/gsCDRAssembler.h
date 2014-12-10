/** @file gsCDRAssembler.h

    @brief Provides assembler and solver for the Poisson equation, incl. adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

// NOT REALLY CHECKED YET, WORK IN PROGRESS!!!

#pragma once

#include <gsAssembler/gsAssemblerBase.h>
#include <gsAssembler/gsPoissonAssembler.h>
#include <gsAssembler/gsAssemblerOptions.h>

#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals
#include <gsAssembler/gsVisitorCDR.h> //
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals
#include <gsAssembler/gsVisitorDg.h>      // Disc. Galerkin interface integrals
#include <gsAssembler/gsVisitorResidual.h>// Residual error estimator

#include <gsPde/gsPoissonPde.h>

//#include <gsAssembler/gsAssemblerUtils.h>

namespace gismo
{

/** @brief
    WRONG, NOT UPDATED YET!\n
    Implementation of an (multiple righ-hand side) Poisson solver.

    The Poisson equation: \f$-\Delta\mathbf{u}=\mathbf{f} \f$

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).
*/
template <class T>
class gsCDRAssembler : public gsPoissonAssembler<T>
{
public:
    typedef gsAssemblerBase<T> Base;
    //typedef gsPoissonAssembler<T> PoissonAssembler;

public:
/** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] bconditions is a gsBoundaryConditions object that holds all boundary conditions.
    \param[in] rhs is the right-hand side of the Poisson equation, \f$\mathbf{f}\f$.
    \param[in] coeff_A diffusion coefficient.
    \param[in] coeff_b convection velocity.
    \param[in] coeff_c reaction coefficient.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
*/
    gsCDRAssembler( gsMultiPatch<T> const         & patches,
                        gsMultiBasis<T> const         & bases,
                        gsBoundaryConditions<T> const & bconditions,
                        const gsFunction<T>           & rhs,
                        const gsFunction<T>           & coeff_A,
                        const gsFunction<T>           & coeff_b,
                        const gsFunction<T>           & coeff_c,
                        dirichlet::strategy           dirStrategy,
                        iFace::strategy               intStrategy = iFace::glue)
    : gsPoissonAssembler<T>( patches,
                       bases,
                       bconditions,
                       rhs,
                       dirStrategy,
                       intStrategy),
       m_coeff_A( &coeff_A ),
       m_coeff_b( &coeff_b ),
       m_coeff_c( &coeff_c )
    {    }

    //void setOptions(const gsAssemblerOptions  & options)

    /// Main assembly routine
    void assemble()
    {
        // If we have a homogeneous Dirichlet problem fill boundary
        // DoFs with zeros
        if ( m_dirStrategy == dirichlet::homogeneous)
            m_ddof.setZero( m_dofMappers[0].boundarySize(), m_rhsFun->targetDim() );

        // If the Dirichlet strategy is elimination then precompute
        // Dirichlet dofs (m_dofMapper excludes these from the system)
        if ( m_dirStrategy == dirichlet::elimination)
            gsPoissonAssembler<T>::computeDirichletDofsL2Proj();
        

        if (m_dofs == 0 ) // Are there any interior dofs ?
        {
            gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
            return;
        }

        // Pre-allocate non-zero elements for each column of the
        // sparse matrix
        int nonZerosPerCol = 1;
        for (int i = 0; i < m_bases.front().dim(); ++i) // to do: improve
            nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;

        m_matrix = gsSparseMatrix<T>(m_dofs, m_dofs); // Clean matrices
        m_matrix.reserve( gsVector<int>::Constant(m_dofs, nonZerosPerCol) );

        // Resize the load vector
        m_rhs.setZero(m_dofs, m_rhsFun->targetDim() );


        // Assemble volume stiffness and load vector integrals
        gsVisitorCDR<T> visitCDR(*m_rhsFun,*m_coeff_A,*m_coeff_b,*m_coeff_c);
        for (unsigned np=0; np < m_patches.nPatches(); ++np )
        {
            //Assemble stiffness matrix and rhs for the local patch
            // with index np and add to m_matrix and m_rhs
            this->apply(visitCDR, np);
        }

        // If requested, force Dirichlet boundary conditions by Nitsche's method
        if ( m_dirStrategy == dirichlet::nitsche )
            gsPoissonAssembler<T>::assembleNitsche();

        // Enforce Neumann boundary conditions
        gsPoissonAssembler<T>::assembleNeumann();

        // If we are in in dg (Discontinuous Galerkin) mode: add
        // interface contributions
        if ( m_intStrategy == iFace::dg )
            gsPoissonAssembler<T>::assembleDg();
        
        // Assembly is done, compress the matrix
        m_matrix.makeCompressed();   
    }

protected:

    // Members from gsPoissonAssembler
    using gsPoissonAssembler<T>::m_patches;
    using gsPoissonAssembler<T>::m_bases;
    using gsPoissonAssembler<T>::m_dofMappers;
    using gsPoissonAssembler<T>::m_ddof;
    using gsPoissonAssembler<T>::m_matrix;
    using gsPoissonAssembler<T>::m_rhs;
    using gsPoissonAssembler<T>::m_dofs;
    using gsPoissonAssembler<T>::m_rhsFun;
    using gsPoissonAssembler<T>::m_bConditions;
    using gsPoissonAssembler<T>::m_dirStrategy;
    using gsPoissonAssembler<T>::m_intStrategy;

    const gsFunction<T> * m_coeff_A;
    const gsFunction<T> * m_coeff_b;
    const gsFunction<T> * m_coeff_c;
};


//////////////////////////////////////////////////
//////////////////////////////////////////////////


} // namespace gismo



