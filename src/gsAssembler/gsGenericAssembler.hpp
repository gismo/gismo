/** @file gsGenericAssembler.h

    @brief Provides an assembler for common IGA matrices

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsVisitorMass.h>
#include <gsAssembler/gsVisitorGradGrad.h>
#include <gsAssembler/gsVisitorMoments.h>
#include <gsAssembler/gsVisitorDg.h>

namespace gismo
{

template <class T>
gsOptionList gsGenericAssembler<T>::defaultOptions()
{
    gsOptionList options = gsAssembler<T>::defaultOptions();
    options.update( gsVisitorDg<T>::defaultOptions(), gsOptionList::addIfUnknown );
    return options;
}

template <class T>
const gsSparseMatrix<T> & gsGenericAssembler<T>::assembleMass(const index_t patchIndex)
{
    // Clean the sparse system
    gsGenericAssembler::refresh();
    const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
    m_system.matrix().reservePerColumn(nz);

    // Assemble mass integrals
    if ( patchIndex < 0 )
    {
        this->template push<gsVisitorMass<T> >();
    }
    else
    {
        gsVisitorMass<T> visitor(*m_pde_ptr);
        this->apply(visitor, patchIndex);
    }

    // Assembly is done, compress the matrix
    this->finalize();

    m_refresh = true;

    return m_system.matrix();
}

template <class T>
const gsSparseMatrix<T> & gsGenericAssembler<T>::assembleStiffness(const index_t patchIndex)
{
    // Clean the sparse system
    gsGenericAssembler::refresh();
    m_system.reserve(m_bases[0], m_options, this->pde().numRhs());
    Base::computeDirichletDofs();

    // Assemble stiffness integrals
    if ( patchIndex < 0 )
    {
        this->template push<gsVisitorGradGrad<T> >();
    }
    else
    {
        gsVisitorGradGrad<T> visitor(*m_pde_ptr);
        this->apply(visitor, patchIndex);
    }

    // Assembly is done, compress the matrix
    this->finalize();

    m_refresh = true;

    return m_system.matrix();
}

template <class T>
const gsMatrix<T> & gsGenericAssembler<T>::assembleMoments(const gsFunction<T> & func, const index_t patchIndex)
{
    // Reset the right-hand side vector
    m_system.rhs().setZero(m_system.cols(), 1);

    // Assemble moment integrals
    gsVisitorMoments<T> mom(func);
    if (patchIndex<0)
        this->push(mom);
    else
        this->apply(mom,patchIndex);

    // Assembly is done, compress the matrix
    this->finalize();

    m_refresh = true;

    return m_system.rhs();
}

template <class T>
const gsSparseMatrix<T> & gsGenericAssembler<T>::assembleDG(const boundaryInterface & iFace)
{
    GISMO_ENSURE( m_options.getInt("InterfaceStrategy")==iFace::dg,
        "Assembling DG terms only makes sense in corresponding setting." );

    // Clean the sparse system
    gsGenericAssembler::refresh();
    const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
    m_system.reserve(nz,1);

    gsVisitorDg<T> visitor(*m_pde_ptr);
    this->apply(visitor, iFace);

    // Assembly is done, compress the matrix
    this->finalize();

    m_refresh = true;

    return m_system.matrix();
}


} // namespace gismo
