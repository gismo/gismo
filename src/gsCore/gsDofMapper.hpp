/** @file gsDofMapper.hpp

    @brief implementation file for the gsDofMapper

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
**/

#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsBasis.h>
#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsMultiBasis.h>

#include <gsMSplines/gsMappedBasis.h>

#include <gsPde/gsBoundaryConditions.h>

namespace gismo
{

template<class T>
void gsDofMapper::init( std::vector<const gsFunctionSet<T> *> const & bases)
{
    const index_t numComp = bases.size();
    m_shift = m_bshift = 0;
    m_curElimId   = -1;
    m_numCpldDofs.assign(numComp+1,1); m_numCpldDofs.front()=0;
    m_offset.clear();

    const size_t nPatches = bases[0]->nPieces();

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
        m_offset.push_back( m_offset.back() + bases[0]->basis(k-1).size() );

    if (nPatches == 1)
    {
        index_t dofsPatches = 0;
        for (index_t comp = 0; comp < numComp; ++comp)
            dofsPatches += bases[comp]->basis(bases[comp]->nPieces()-1).size();
        m_numFreeDofs.assign(numComp+1,m_offset.back() + dofsPatches);
	m_numFreeDofs.front()=0;
    }
    //Assuming each component are of same size;
    //i.e. bases[comp]->back().size() are equal for all comp
    else
    {
      m_numFreeDofs.assign(numComp+1, (m_offset.back() + bases[0]->basis(bases[0]->nPieces()-1).size())); m_numFreeDofs.front()=0;
    }

    m_numElimDofs.assign(numComp+1,0);
    m_dofs.resize(numComp, std::vector<index_t>(m_numFreeDofs.back(), 0));
}

template<class T>
void gsDofMapper::init(const gsBasis<T>         &basis,
                       const gsBoundaryConditions<T> &bc,
                       index_t nComp,
                       int unk)
{
    gsMultiBasis<T> mbasis(basis);
    init(mbasis, gsBoxTopology(), bc, nComp, unk);
}

template<class T>
/// Initialize from a (possibly multi-patch) function set.
/// \param topology Interface topology between patches.
///        It is valid to pass an empty gsBoxTopology() for single-patch
///        bases or when no inter-patch conformity is required; in that
///        case no interface matching is performed.
void gsDofMapper::init(const gsFunctionSet<T> & bases,
                       const gsBoxTopology & topology,
                       index_t nComp, int unk, bool conforming)
{
    GISMO_ASSERT(nComp>0,"Zero components");
    m_shift = m_bshift = 0;
    m_curElimId   = -1;
    m_numCpldDofs.assign(nComp+1, 1); m_numCpldDofs.front()=0;
    m_numElimDofs.assign(nComp+1,0);
    m_offset.clear();

    const size_t nPatches = bases.nPieces();

    // Initialize offsets and dof holder
    m_offset.reserve( nPatches );
    m_offset.push_back(0);
    for (size_t k = 1; k < nPatches; ++k)
        m_offset.push_back( m_offset.back() + bases.basis(k-1).size() );

    m_numFreeDofs.assign(1+nComp,m_offset.back() + bases.basis(nPatches-1).size()); m_numFreeDofs.front()=0;

    m_dofs.resize(nComp, std::vector<index_t>(m_numFreeDofs.back(), 0));
    if (!conforming)
        return;

    gsMatrix<index_t> b1, b2;
    for ( gsBoxTopology::const_iiterator it = topology.iBegin();
            it != topology.iEnd(); ++it )
    {
        const gsBasis<T> & basis1 = bases.basis(it->first().patch);
        const gsBasis<T> & basis2 = bases.basis(it->second().patch);
        basis1.matchWith(*it, basis2, b1, b2);
        for (size_t i=0; i!=this->componentsSize(); ++i)
            this->matchDofs(it->first().patch, b1,it->second().patch, b2,i);
    }
}

template<class T>
void gsDofMapper::_addBCs(  const gsFunctionSet<T>         &basis,
                            const gsBoundaryConditions<T>  &bc,
                            index_t nComp,
                            int unk)
{
// Strong Dirichlet conditions
    gsMatrix<index_t> bnd, bnd1;
    for (typename gsBoundaryConditions<T>::const_iterator
         it = bc.begin("Dirichlet") ; it != bc.end("Dirichlet"); ++it )
    {
        if (unk!=-1 && it->unknown() != unk) continue;

        GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->numPatches(),
                    "Problem: a boundary condition is set on a patch id which does not exist.");

        bnd = basis.basis(it->ps.patch).boundary(it->ps.side());
        this->markBoundary(it->ps.patch, bnd, it->unkComponent());
    }

    // Clamped boundary condition (per DoF)
    for (typename gsBoundaryConditions<T>::const_iterator
         it = bc.begin("Clamped") ; it != bc.end("Clamped"); ++it )
    {
        if (unk!=-1 && it->unknown() != unk) continue;

        GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->numPatches(),
                        "Problem: a boundary condition is set on a patch id which does not exist.");

        const index_t cc = it->unkComponent();
        bnd = basis.basis(it->ps.patch).boundary(it->ps.side());
        bnd1= basis.basis(it->ps.patch).boundaryOffset(it->ps.side(), 1);
        if (!it->ps.parameter())
            bnd.swap(bnd1);
        for (index_t c = 0; c!=nComp; c++) // for all components
        {
            if (c==cc || cc==-1 )
                for (index_t k = 0; k < bnd.size(); ++k)
                    this->matchDof(it->ps.patch, (bnd)(k, 0),
                                   it->ps.patch, (bnd1)(k, 0), c);
        }
    }

    // Collapsed
    for (typename gsBoundaryConditions<T>::const_iterator
         it = bc.begin("Collapsed") ; it != bc.end("Collapsed"); ++it )
    {
        if (unk!=-1 && it->unknown() != unk) continue;

        GISMO_ASSERT(static_cast<size_t>(it->ps.patch) < this->numPatches(),
                        "Problem: a boundary condition is set on a patch id which does not exist.");
        const index_t cc = it->unkComponent();
        bnd = basis.basis(it->ps.patch).boundary(it->ps.side());
        // match all DoFs to the first one of the side
        for (index_t c = 0; c!=nComp; c++) // for all components
        {
            if (c==cc || cc==-1)
                for (index_t k = 0; k < bnd.size() - 1; ++k)
                    this->matchDof(it->ps.patch, (bnd)(0, 0),
                                   it->ps.patch, (bnd)(k + 1, 0), c);
        }
    }

    // Coupled boundary condition (per DoF)
    for (typename gsBoundaryConditions<T>::const_cpliterator
         it = bc.coupledBegin(); it != bc.coupledEnd(); ++it )
    {
        if (unk!=-1 && it->unknown != unk) continue;

        GISMO_ASSERT(static_cast<size_t>(it->ifc.first().patch) < this->numPatches(),
                        "Problem: a boundary condition is set on a patch id which does not exist.");
        GISMO_ASSERT(static_cast<size_t>(it->ifc.second().patch) < this->numPatches(),
                        "Problem: a boundary condition is set on a patch id which does not exist.");

        const index_t cc = it->component;
        bnd = basis.basis(it->ifc.first().patch).boundary(it->ifc.first().side());
        bnd1= basis.basis(it->ifc.second().patch).boundary(it->ifc.second().side());
        GISMO_ASSERT(bnd.rows() == bnd1.rows(),
                        "Problem: trying to couple boundaries of different size.");

        // match all DoFs to the first one of the side
        for (index_t c = 0; c!=nComp; c++) // for all components
        {
            if (c==cc || cc==-1)
            {
                for (index_t k = 0; k < bnd.size() -1; ++k)
                    this->matchDof(it->ifc.first() .patch, (bnd)(0, 0),
                                   it->ifc.first() .patch, (bnd)(k + 1, 0), c);
                for (index_t k = 0; k < bnd1.size(); ++k)
                    this->matchDof(it->ifc.second().patch, (bnd1)(k, 0),
                                   it->ifc.first().patch,  (bnd)(k, 0), c);
            }
        }
    }

    // Corners
    for (typename gsBoundaryConditions<T>::const_citerator
            it = bc.cornerBegin() ; it != bc.cornerEnd(); ++it )
    {
        if (unk!=-1 && it->unknown != unk) continue;
        for (index_t r = 0; r!=nComp; ++r)
        {
            if (it->component!=-1 && it->component!=r) continue;
            GISMO_ASSERT(static_cast<size_t>(it->patch) < this->numPatches(),
                            "Problem: a corner boundary condition is set on a patch id which does not exist.");
            this->eliminateDof(basis.basis(it->patch).functionAtCorner(it->corner), it->patch, it->component);
        }
    }
}

template<class T>
void gsDofMapper::init(const gsFunctionSet<T>        &basis,
                       const gsBoxTopology           &topology,
                       const gsBoundaryConditions<T> &bc,
                       index_t nComp,
                       int unk,
                       bool conforming)
{
    if (dynamic_cast<const gsMappedBasis<2,T>*>(&basis) !=nullptr || 
        dynamic_cast<const gsMappedBasis<3,T>*>(&basis) !=nullptr)
            this->setIdentity(basis.nPieces(), basis.size(), nComp);
    else
        init(basis,topology,nComp,unk,conforming);
    _addBCs(basis,bc,nComp,unk);
}

template<class T>
void gsDofMapper::initSingle( const gsBasis<T> & basis, index_t nComp)
{
    GISMO_ASSERT(nComp>0,"Zero components");
    m_shift = m_bshift = 0;
    m_curElimId   = -1;
    m_numFreeDofs.assign(nComp+1,basis.size()); m_numFreeDofs.front()=0;
    m_numCpldDofs.assign(nComp+1,1); m_numCpldDofs.front()=0;
    m_numElimDofs.assign(nComp+1,0);
    m_offset.resize(1,0);
    m_dofs.resize(nComp, std::vector<index_t>(m_numFreeDofs.back(), 0));
}

}//namespace gismo

