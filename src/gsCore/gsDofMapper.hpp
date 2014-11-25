/** @file gsDofMapper.hpp

    @brief implementation file for the gsDofMapper

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
**/

#include <gsCore/gsMultiBasis.h>

namespace gismo 
{

template<class T>
void gsDofMapper::init( const gsMultiBasis<T> & bases)
{
    m_curElimId     = -1;
    m_curCouplingId =  1;
    m_offset.clear();

    const size_t nPatches = bases.nBases();

    // Initialize offsets and dof holder
    m_offset.reserve( nPatches );
    m_offset.push_back(0);
    for (size_t k = 1; k < nPatches; ++k)
    {
        m_offset.push_back( m_offset.back() + bases[k-1].size() );
    }

    m_numFreeDofs = m_offset.back() + bases.back().size();

    m_dofs.resize( m_numFreeDofs, 0);
}

template<class T>
void gsDofMapper::init(
        const gsMultiBasis<T>         &basis,
        const gsBoundaryConditions<T> &bc
        )
{
    init(basis);

    for (typename gsBoundaryConditions<T>::const_iterator
         it = bc.dirichletBegin() ; it != bc.dirichletEnd(); ++it )
    {
        gsMatrix<unsigned> * bnd = basis[it->ps.patch].boundary( it->ps.side );
        markBoundary(it->ps.patch, *bnd);
        delete bnd;
    }
}

template<class T>
void gsDofMapper::initSingle( const gsBasis<T> & basis)
{
    m_curElimId     = -1;
    m_curCouplingId =  1;

    m_offset.resize(1,0);
    m_numFreeDofs = basis.size();
    m_dofs.resize( m_numFreeDofs, 0);
}

}

