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
void gsPoissonAssembler<T>::assemble()
{
    // Compute the Dirichlet Degrees of freedom
    this->computeDirichletDofs();

    if (m_dofs == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
        return;
    }

    // Allocate memory for the sparse matrix and right-hand side
    this->reserveSparseSystem();

    // Assemble volume integrals
    this->template push<gsVisitorPoisson<T> >();

    // Enforce Neumann boundary conditions
    this->template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides() );
     
     // If requested, enforce Dirichlet boundary conditions by Nitsche's method
     if ( m_options.dirStrategy == dirichlet::nitsche )
         this->template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());
    
     // If requested, enforce Dirichlet boundary conditions by diagonal penalization
     if ( m_options.dirStrategy == dirichlet::penalize )
         this->penalizeDirichlet();
     
    if ( m_options.intStrategy == iFace::dg )
        gsWarn <<"DG option is not available. Results will be incorrect.\n";
    
    // Assembly is done, compress the matrix
    this->finalize();
}

}// namespace gismo
