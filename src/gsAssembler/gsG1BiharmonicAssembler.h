/** @file gsG1BiharmonicAssembler.h

    @brief Provides assembler for a homogenius Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/


#pragma once

#include <gsAssembler/gsAssembler.h>

#include <gsPde/gsBiharmonicPde.h>

#include <gsAssembler/gsVisitorBiharmonic.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorLaplaceBiharmonic.h>
//#include <gsAssembler/gsVisitorNitscheBiharmonic.h>

namespace gismo
{

/** @brief
    Implementation of a homogeneous Biharmonic Assembler.

    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can only be enforced strongly (i.e Nitsche is
    not implemented).
*/
template <class T, class bhVisitor = gsVisitorBiharmonic<T> >
class gsG1BiharmonicAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

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
    gsG1BiharmonicAssembler( gsMultiPatch<T> const         & patches,
                           gsMappedBasis<2,T> const         & bases,
                           gsBoundaryConditions<T> const & bconditions,
                           gsBoundaryConditions<T> const & bconditions2,
                           const gsFunction<T>           & rhs,
                           dirichlet::strategy           dirStrategy = dirichlet::none,
                           iFace::strategy               intStrategy = iFace::none)
    : m_ppde(patches,bconditions,bconditions2,rhs), m_bases(bases)
    {
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        refresh();
    }

    void refresh();
    
    void assemble();

protected:

    // fixme: add constructor and remove this
    gsBiharmonicPde<T> m_ppde;

    // G1 Basis
    gsMappedBasis<2,T> const m_bases;

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
};

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::refresh()
{
    // We use predefined helper which initializes the system matrix
    // rows and columns using the same test and trial space
    //Base::scalarProblemGalerkinRefresh();
    gsVector<index_t> sz(m_bases.nPatches());
    for (size_t np = 0; np < m_bases.nPatches(); np++)
        sz[np] = m_bases.size(np);

    gsDofMapper map(sz);

    gsMatrix<index_t> boundaryDofs;
    typedef std::vector< patchSide >::const_iterator b_const_iter;
    for(b_const_iter iter = m_ppde.domain().bBegin();iter!=m_ppde.domain().bEnd();++iter)
    {
        boundaryDofs = m_bases.getBase(iter->patch).boundaryOffset(*iter, 0);
        map.markBoundary(iter->patch, boundaryDofs);
    }

    // TODO Matching Interface BFs

    map.finalize();

    map.print();

    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(map);
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::assemble()
{
    GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    index_t nz = 1;
    for (short_t i = 0; i != 2; ++i) // 2 == dim
        nz *= static_cast<index_t>(2 * m_bases.degree(0,i) + 1 + 0.5); // Patch 0
    nz = static_cast<index_t>(nz*(1.333333));
    m_system.reserve(nz, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    //Base::computeDirichletDofs();
    m_ddof.resize(m_system.numUnknowns());
    m_ddof[0].setZero(m_system.colMapper(0).boundarySize(), m_system.unkSize(0) * m_system.rhs().cols());


    // Assemble volume integrals
    Base::template push<bhVisitor >();
    
    // Newman conditions of first kind
    Base::template push<gsVisitorNeumann<T> >(
        m_ppde.bcFirstKind().neumannSides() );
    
    // Newman conditions of second kind
    Base::template push<gsVisitorLaplaceBiharmonic<T> >(
        m_ppde.bcSecondKind().laplaceSides() );

    if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
        gsWarn <<"DG option ignored.\n";

    /*
    // If requested, force Dirichlet boundary conditions by Nitsche's method
    this->template push<gsVisitorNitscheBiharmonic<T> >(
    m_ppde.bcSecondKind().dirichletSides() );
    */
    
    // Assembly is done, compress the matrix
    Base::finalize();
}


} // namespace gismo



