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
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorNitsche.h>

namespace gismo
{

template <class T>
gsOptionList gsGenericAssembler<T>::defaultOptions()
{
    gsOptionList options = gsAssembler<T>::defaultOptions();
    options.update( gsVisitorDg<T>::defaultOptions(), gsOptionList::addIfUnknown );
    options.update( gsVisitorNitsche<T>::defaultOptions(), gsOptionList::addIfUnknown );
    return options;
}

template <class T>
const gsSparseMatrix<T> & gsGenericAssembler<T>::assembleMass(const index_t patchIndex, const bool refresh)
{
    // Clean the sparse system
    if (refresh)
    {
        gsGenericAssembler::refresh();
        m_system.reserve(m_bases[0], m_options, 1);
        Base::computeDirichletDofs();
    }

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

    return m_system.matrix();
}

template <class T>
const gsSparseMatrix<T> & gsGenericAssembler<T>::assembleStiffness(const index_t patchIndex, const bool refresh)
{
    // Clean the sparse system
    if (refresh)
    {
        gsGenericAssembler::refresh();
        m_system.reserve(m_bases[0], m_options, 1);
        Base::computeDirichletDofs();
    }

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

    return m_system.matrix();
}

template <class T>
const gsMatrix<T> & gsGenericAssembler<T>::assembleMoments(const gsFunction<T> & func, const index_t patchIndex, const bool refresh)
{
    // Clean the sparse system
    if (refresh)
    {
        gsGenericAssembler::refresh();
        m_system.rhs().setZero(m_system.cols(), 1);
    }

    // Assemble moment integrals
    gsVisitorMoments<T> mom(func);
    if (patchIndex<0)
        this->push(mom);
    else
        this->apply(mom,patchIndex);

    // Assembly is done, compress the matrix
    this->finalize();

    return m_system.rhs();
}

template <class T>
const gsSparseMatrix<T> & gsGenericAssembler<T>::assembleDG(const boundaryInterface & iFace, const bool refresh)
{
    GISMO_ENSURE( m_options.getInt("InterfaceStrategy")==iFace::dg,
        "Assembling DG terms only makes sense in corresponding setting." );

    // Clean the sparse system
    if (refresh)
    {
        gsGenericAssembler::refresh();
        m_system.reserve(m_bases[0], m_options, 1);
        Base::computeDirichletDofs();
    }

    gsVisitorDg<T> visitor(*m_pde_ptr);
    this->apply(visitor, iFace);

    // Assembly is done, compress the matrix
    this->finalize();

    return m_system.matrix();
}

template <class T>
const gsSparseMatrix<T> & gsGenericAssembler<T>::assembleNeumann(const boundary_condition<T> & bc, const bool refresh)
{
    GISMO_ASSERT( bc.ctype() == "Neumann", "gsGenericAssembler::assembleNeumann: Got " << bc.ctype() << "bc.");

    // Clean the sparse system
    if (refresh)
    {
        gsGenericAssembler::refresh();
        m_system.rhs().setZero(m_system.cols(), 1);
    }

    gsVisitorNeumann<T> visitor(*m_pde_ptr,bc);
    this->apply(visitor, bc.patch(), bc.side());

    // Assembly is done, compress the matrix
    this->finalize();

    return m_system.matrix();
}

template <class T>
const gsSparseMatrix<T> & gsGenericAssembler<T>::assembleNitsche(const boundary_condition<T> & bc, const bool refresh)
{
    GISMO_ASSERT( bc.ctype() == "Dirichlet", "gsGenericAssembler::assembleNitsche: Got " << bc.ctype() << "bc.");

    // Clean the sparse system
    if (refresh)
    {
        gsGenericAssembler::refresh();
        m_system.reserve(m_bases[0], m_options, 1);
        Base::computeDirichletDofs();
    }


    gsVisitorNitsche<T> visitor(*m_pde_ptr,bc);
    this->apply(visitor, bc.patch(), bc.side());

    // Assembly is done, compress the matrix
    this->finalize();

    return m_system.matrix();
}


} // namespace gismo
