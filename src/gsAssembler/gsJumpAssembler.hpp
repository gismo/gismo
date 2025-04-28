/** @file gsJumpAssembler.hpp

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
#include <gsAssembler/gsVisitorDg.h>      // DG interface integrals

namespace gismo
{

template <class T>
gsOptionList gsJumpAssembler<T>::defaultOptions()
{
    gsOptionList options = gsAssembler<T>::defaultOptions();
    options.update( gsVisitorDg<T>::defaultOptions(), gsOptionList::addIfUnknown );
    options.setReal("DG.Alpha", 0);
    options.setReal("DG.Beta", 0);
    options.setReal("DG.Penalty", 1);
    return options;
}

template<class T>
void gsJumpAssembler<T>::refresh()
{
    // We use predefined helper which initializes the system matrix
    // rows and columns using the same test and trial space
    Base::scalarProblemGalerkinRefresh();
}

template<class T>
void gsJumpAssembler<T>::assemble()
{
    GISMO_ASSERT(m_system.initialized(), 
                 "Sparse system is not initialized, call initialize() or refresh()");

    // Reserve sparse system
    m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    Base::computeDirichletDofs();

    // Clean the sparse system
    // m_system.setZero(); //<< this call leads to a quite significant performance degrade!

    //Base::template pushInterface<gsVisitorDg<T> >();
    gsVisitorDg<T> visitor(*m_pde_ptr);

    const boundaryInterface & iFace = //recover master element
        ( m_bases[0][m_it->first() .patch].numElements(m_it->first() .side() ) <
            m_bases[0][m_it->second().patch].numElements(m_it->second().side() ) ?
            m_it->getInverse() : *m_it );

    this->apply(visitor, iFace);
        
    // Assembly is done, compress the matrix
    Base::finalize();
}


}// namespace gismo
