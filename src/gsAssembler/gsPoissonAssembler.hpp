/** @file gsPoissonAssembler.hpp

    @brief Provides assembler implementation for the Poisson equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris, J. Sogn
*/


#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals

namespace gismo
{

template<class T>
void gsPoissonAssembler<T>::refresh()
{
    // We use predefined helper which initializes the system matrix
    // rows and columns using the same test and trial space
    Base::scalarProblemGalerkinRefresh();
}

template<class T>
void gsPoissonAssembler<T>::assemble()
{
    GISMO_ASSERT(m_system.initialized(), 
                 "Sparse system is not initialized, call initialize() or refresh()");

    if ( 0 == Base::numDofs() ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
        return;
    }

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    Base::computeDirichletDofs();

    // Clean the sparse system
    m_system.setZero();

    // Assemble volume integrals
    Base::template push<gsVisitorPoisson<T> >();

    // Enforce Neumann boundary conditions
    Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );
     
     // If requested, enforce Dirichlet boundary conditions by Nitsche's method
     if ( m_options.dirStrategy == dirichlet::nitsche )
         Base::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());
     // If requested, enforce Dirichlet boundary conditions by diagonal penalization
     else if ( m_options.dirStrategy == dirichlet::penalize )
         Base::penalizeDirichletDofs();
     
    if ( m_options.intStrategy == iFace::dg )
        gsWarn <<"DG option is ignored.\n";
    
    // Assembly is done, compress the matrix
    Base::finalize();
}

}// namespace gismo
