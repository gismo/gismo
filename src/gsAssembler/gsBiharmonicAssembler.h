/** @file gsBiharmonicAssembler.h

    @brief Provides assembler for a homogenius Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/


#pragma once

#include <gsAssembler/gsAssemblerBase.h>
#include <gsAssembler/gsPoissonAssembler.h>
#include <gsAssembler/gsVisitorBiharmonic.h>
#include <gsAssembler/gsVisitorNeumannBiharmonic.h>
#include <gsPde/gsPoissonPde.h>

namespace gismo
{

/** @brief
    Implementation of a homogenius Biharmonic Assembler.

    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can only be enforced strongly (i.e Nitsche is
    not implemented).
*/
template <class T>
class gsBiharmonicAssembler : public gsPoissonAssembler<T>
{
public:
    typedef gsAssemblerBase<T> Base;

public:
/** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] bconditions  is a gsBoundaryConditions object that holds boundary conditions on the form:
    \f[ \text{Dirichlet: } u = g \text{ on } \Gamma, \text{ and Neumann: } \nabla \Delta u \cdot \mathbf{n} = h \text{ on } \Gamma\f]
    \param[in] bconditions2 is a gsBoundaryConditions object that holds Neumann boundary conditions on the form:
    \f[\text{Neumann: } \nabla \Delta u \cdot \mathbf{n} = g\, \rightarrow \,(g,\nabla v \cdot \mathbf{n})_\Gamma, \f] where \f$ g \f$ is the Neumann data,
    \f$ v \f$ is the test function and \f$ \Gamma\f$ is the boundary side.
    \param[in] rhs is the right-hand side of the Biharmonic equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary in the \em bconditions object.
    \param[in] intStrategy option for the treatment of patch interfaces
*/
    gsBiharmonicAssembler( gsMultiPatch<T> const         & patches,
                           gsMultiBasis<T> const         & bases,
                           gsBoundaryConditions<T> const & bconditions,
                           gsBoundaryConditions<T> const & bconditions2,
                           const gsFunction<T>           & rhs,
                           dirichlet::strategy           dirStrategy,
                           iFace::strategy               intStrategy = iFace::glue)
    : gsPoissonAssembler<T>( patches,
                       bases,
                       bconditions,
                       rhs,
                       dirStrategy,
                       intStrategy),
      m_bConditions2(bconditions2)
    {    }

    /// Main assembly routine
    void assemble()
    {
         if (m_dofs == 0 ) // Are there any interior dofs ?
        {
            gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
            return;
        }

        gsPoissonAssembler<T>::computeDirichletDofs();

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
        GISMO_ASSERT(m_patches.nPatches() == 1, "Only valid for Multipatch");

        gsVisitorBiharmonic<T> visitBiHar(*m_rhsFun);
        for (unsigned np=0; np < m_patches.nPatches(); ++np )
            this->apply(visitBiHar, np);


        for ( typename gsBoundaryConditions<T>::const_iterator
                  it = m_bConditions.neumannBegin();
              it != m_bConditions.neumannEnd(); ++it )
        {
            gsVisitorNeumann<T> neumann(*it->function(), it->side());
            this->apply(neumann, it->patch(), it->side() );
        }
        for ( typename gsBoundaryConditions<T>::const_iterator
                  it = m_bConditions2.neumannBegin();
              it != m_bConditions2.neumannEnd(); ++it )
        {
            gsVisitorNeumannBiharmonic<T> neumann(*it->function(), it->side());
            this->apply(neumann, it->patch(), it->side() );
        }

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
    using gsPoissonAssembler<T>::m_options;

    // BC for the second kind of Neumann condition (See documentation in constructor.
    gsBoundaryConditions<T> m_bConditions2;

};


//////////////////////////////////////////////////
//////////////////////////////////////////////////


} // namespace gismo



