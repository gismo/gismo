/** @file gsCDRAssembler.hpp

    @brief Provides assembler and solver for the Poisson equation, incl. adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris
*/


#include <gsAssembler/gsAssemblerOptions.h>

#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals
#include <gsAssembler/gsVisitorCDR.h> //
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals

//#include <gsAssembler/gsAssemblerUtils.h>

namespace gismo
{

template <class T>
gsOptionList gsCDRAssembler<T>::defaultOptions()
{
    gsOptionList options = gsAssembler<T>::defaultOptions();
    options.update( gsVisitorNitsche<T>::defaultOptions(), gsOptionList::addIfUnknown );
    return options;
}

template<class T>
void gsCDRAssembler<T>::assemble()
{
    GISMO_ASSERT(m_system.initialized(),
                 "Sparse system is not initialized, call initialize() or refresh()");

    // Reserve sparse system
    m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    Base::computeDirichletDofs();

    if (0 == this->numDofs()) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" << "\n";
        return;
    }

    // Assemble volume integrals
    Base::template push<gsVisitorCDR<T> >();

    // Enforce Neumann boundary conditions
    Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides());

    // If requested, enforce Dirichlet boundary conditions by Nitsche's method
    if (m_options.getInt("DirichletStrategy") == dirichlet::nitsche)
        Base::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());

    // If requested, enforce Dirichlet boundary conditions by diagonal penalization
    else if (m_options.getInt("DirichletStrategy") == dirichlet::penalize)
            Base::penalizeDirichletDofs();

    // If we are in in dg (Discontinuous Galerkin) mode: add
    // interface contributions
    if (m_options.getInt("InterfaceStrategy") == iFace::dg)
        gsWarn << "DG option is ignored.\n";

    // Assembly is done, compress the matrix
    Base::finalize();
}

} // namespace gismo
