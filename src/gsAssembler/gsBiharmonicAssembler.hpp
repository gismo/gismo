/** @file gsBiharmonicAssembler.hpp

    @brief Provides assembler for a homogenius Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#include <gsAssembler/gsVisitorBiharmonic.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorNeumannBiharmonic.h>
//#include <gsAssembler/gsVisitorNitscheBiharmonic.h>

namespace gismo
{

template <class T, class bhVisitor>
void gsBiharmonicAssembler<T,bhVisitor>::refresh()
{
    // We use predefined helper which initializes the system matrix
    // rows and columns using the same test and trial space
    Base::scalarProblemGalerkinRefresh();
}

template <class T, class bhVisitor>
void gsBiharmonicAssembler<T,bhVisitor>::assemble()
{
    GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
    m_system.reserve(nz, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    Base::computeDirichletDofs();
    
    // Assemble volume integrals
    Base::template push<bhVisitor >();
    
    // Newman conditions of first kind
    Base::template push<gsVisitorNeumann<T> >(
        m_ppde.bcFirstKind().neumannSides() );
    
    // Newman conditions of second kind
    Base::template push<gsVisitorNeumannBiharmonic<T> >(
        m_ppde.bcSecondKind().neumannSides() );

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



