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
    m_curElimId   = -1;
    m_numCpldDofs =  1;
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
void gsDofMapper::init( std::vector<const gsMultiBasis<T> *> const & bases)
{
    m_curElimId   = -1;
    m_numCpldDofs =  1;
    m_offset.clear();

    const size_t nPatches = bases[0]->nBases();
    const index_t numComp = bases.size();

    //Checking if bases are same size in for components.
    std::vector<index_t> offsets(nPatches);
    for (index_t comp = 0; comp < numComp; ++comp)
    {
        for (size_t k = 0; k < nPatches; ++k)
        {
            if (comp != 0)
            {
                GISMO_ASSERT(offsets[k] == bases[comp]->basis(k).size(),
                         "The sizes of the bases are not the same for every component. Dofmapper requries this!");
            }
            offsets[k] =  bases[comp]->basis(k).size();
        }
    }
    // Initialize offsets and dof holder
    m_offset.reserve( nPatches );
    m_offset.push_back(0);
    for (size_t k = 1; k < nPatches; ++k)
    {
        m_offset.push_back( m_offset.back() + bases[0]->basis(k-1).size() );
    }

    if (nPatches == 1)
    {
        index_t dofsPatches = 0;
        for (index_t comp = 0; comp < numComp; ++comp)
        {
            dofsPatches += bases[comp]->back().size();
        }
        m_numFreeDofs = m_offset.back() + dofsPatches;
    }
    //Assuming each component are of same size;
    //i.e. bases[comp]->back().size() are equal for all comp
    else
    {
        m_numFreeDofs = (m_offset.back() + bases[0]->back().size())*nPatches;
    }

    m_dofs.resize( m_numFreeDofs, 0);
}

template<class T>
void gsDofMapper::init(const gsMultiBasis<T>         &basis,
                       const gsBoundaryConditions<T> &bc, int unk)
{
    init(basis);

    /// \todo move this code to gsMultiBasis::getMapper
    for (typename gsBoundaryConditions<T>::const_iterator
         it = bc.dirichletBegin() ; it != bc.dirichletEnd(); ++it )
    {
        if (unk == -1 || it->unknown() == unk) // special value -1 eliminates all BCs found
        {
            GISMO_ASSERT(it->ps.patch < static_cast<index_t>(m_offset.size()),
                         "Problem: a boundary condition is set on a patch id which does not exist.");

            gsMatrix<index_t> bnd = basis[it->ps.patch].boundary(it->ps.side());
            markBoundary(it->ps.patch, bnd);
        }
    }

    // corners
    for (typename gsBoundaryConditions<T>::const_citerator
         it = bc.cornerBegin() ; it != bc.cornerEnd(); ++it )
    {
        if (unk == -1 || it->unknown == unk)
        {
            GISMO_ASSERT(it->patch < m_offset.size(),
                         "Problem: a corner boundary condition is set on a patch id which does not exist.");

            eliminateDof(basis[it->patch].functionAtCorner(it->corner), it->patch);
        }
    }
}

template<class T>
void gsDofMapper::initSingle( const gsBasis<T> & basis)
{
    m_curElimId   = -1;
    m_numCpldDofs =  1;

    m_offset.resize(1,0);
    m_numFreeDofs = basis.size();
    m_dofs.resize( m_numFreeDofs, 0);
}

}

