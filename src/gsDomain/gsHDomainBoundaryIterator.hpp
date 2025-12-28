/** @file gsHDomainBoundaryIterator.hpp

    @brief Provides implementation of iterator on boundary of hierarchical basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsHDomainBoundaryIterator.h>

namespace gismo
{

template<typename T, short_t d, typename Z>
gsHDomainBoundaryIterator<T,d,Z>::gsHDomainBoundaryIterator(const gsHTree<d,Z> & tree,
                              const gsHTensorBasis<d,T> & basis,
                              const boxSide & s,
                              const index_t patchId)
    :
    gsDomainIterator<T>(0, patchId),
    m_tree(tree),
    m_basis(basis),
    m_patchId(patchId)
{
    init(m_tree,s);
}

template<typename T, short_t d, typename Z>
gsHDomainBoundaryIterator<T,d,Z>::gsHDomainBoundaryIterator(const gsHDomain<d,T,Z> & domain,
                              const gsHTensorBasis<d,T> & basis,
                              const boxSide & s)
    : gsHDomainBoundaryIterator(domain.tree(),basis,s, domain.patchId())
{
}

template<typename T, short_t d, typename Z>
void gsHDomainBoundaryIterator<T,d,Z>::init(const gsHTree<d,Z> & tree, const boxSide & s)
{
    // Initialize mesh data
    m_meshStart.resize(d);
    m_meshEnd  .resize(d);
    m_curElement.resize(d);

    // Allocate breaks
    m_breaks = std::vector<std::vector<T> >(d, std::vector<T>());

    // Get the side information
    par = s.parameter();
    dir = s.direction();

    this->initLeaf(tree);
}

template<typename T, short_t d, typename Z>
void gsHDomainBoundaryIterator<T,d,Z>::next()
{
    bool isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);

    if (!isGood) // went through all elements in m_leaf
        isGood = nextLeaf();
}

template<typename T, short_t d, typename Z>
void gsHDomainBoundaryIterator<T,d,Z>::next(index_t increment)
{
    bool isGood(true);
    for (index_t i = 0; i < increment; i++)
        isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);

    if (isGood) // went through all elements in m_leaf
        isGood = nextLeaf();
}

template<typename T, short_t d, typename Z>
void gsHDomainBoundaryIterator<T,d,Z>::reset()
{
    initLeaf(m_tree);
}

template<typename T, short_t d, typename Z>
gsVector<T> gsHDomainBoundaryIterator<T,d,Z>::lowerCorner() const
{
    gsVector<T> lower;
    lower.resize(d);
    for (short_t i = 0; i < dir ; ++i)
        lower[i] = *m_curElement[i]; // in gsTensorDomainBoundaryIterator, we have: lower[i]  = m_curElement[i].lowerCorner().value();
    lower[dir] = (par ? *(m_curElement[dir]+1) : *m_curElement[dir] ); // in gsTensorDomainBoundaryIterator, we have: lower[dir] = (par ? m_curElement[dir].upperCorner().value() : m_curElement[dir].lowerCorner().value() );
    for (short_t i = dir+1; i < d; ++i)
        lower[i] = *m_curElement[i]; // in gsTensorDomainBoundaryiterator, we have: lower[i]  = m_curElement[i].lowerCorner().value();
    return lower;
}

template<typename T, short_t d, typename Z>
gsVector<T> gsHDomainBoundaryIterator<T,d,Z>::upperCorner() const
{
    gsVector<T> upper;
    upper.resize(d);
    for (short_t i = 0; i < dir ; ++i)
        upper[i] = *(m_curElement[i]+1); // in gsTensorDomainBoundaryiterator, we have: upper[i]  = m_curElement[i].upperCorner().value();
    upper[dir] = (par ? *(m_curElement[dir]+1) : *m_curElement[dir] ); // in gsTensorDomainBoundaryIterator, we have: upper[dir] = (par ? m_curElement[dir].upperCorner().value() : m_curElement[dir].upperCorner().value() );
    for (short_t i = dir+1; i < d; ++i)
        upper[i] = *(m_curElement[i]+1); // in gsTensorDomainBoundaryiterator, we have: upper[i]  = m_curElement[i].upperCorner().value();
    return upper;
}

template<typename T, short_t d, typename Z>
const T gsHDomainBoundaryIterator<T,d,Z>::getPerpendicularCellSize() const
{
    return *(m_curElement[dir]+1) - *m_curElement[dir];
}

template<typename T, short_t d, typename Z>
int gsHDomainBoundaryIterator<T,d,Z>::getLevel() const
{
    return m_leaf.level();
}

template<typename T, short_t d, typename Z>
void gsHDomainBoundaryIterator<T,d,Z>::initLeaf(const hDomain & tree_domain)
{
    // Get the first leaf
    m_leaf = tree_domain.beginLeafIterator();

    for (; m_leaf.good(); m_leaf.next() )
    {
        // Check if this leaf is on our side
        if ( leafOnBoundary() )
        {
            updateLeaf();
            return;
        }
    }
    GISMO_ERROR("No leaves.\n");
}

template<typename T, short_t d, typename Z>
bool gsHDomainBoundaryIterator<T,d,Z>::nextLeaf()
{
    for (m_leaf.next(); m_leaf.good(); m_leaf.next() )
    {
        // Check if this leaf is on our side
        if ( leafOnBoundary() )
        {
            updateLeaf();
            return true;
        }
    }
    return false;
}

template<typename T, short_t d, typename Z>
bool gsHDomainBoundaryIterator<T,d,Z>::leafOnBoundary() const
{
    if ( par )
    {
        // AM: a little ugly for now, to be improved
        size_t diadicSize;
        const gsHTensorBasis<d,T> * hbasis = dynamic_cast<const gsHTensorBasis<d,T> * >(&m_basis);
        if (hbasis->manualLevels() )
        {
            gsKnotVector<T> kv = hbasis->tensorLevel(m_leaf.level()).knots(dir);
            index_t start = 0;
            index_t end  = kv.uSize()-1;
            hbasis->_knotIndexToDiadicIndex(m_leaf.level(),dir,start);
            hbasis->_knotIndexToDiadicIndex(m_leaf.level(),dir,end);
            diadicSize = end - start;
        }
        else
            diadicSize = hbasis->tensorLevel(m_leaf.level()).knots(dir).uSize() - 1;

        return static_cast<size_t>(m_leaf.upperCorner().at(dir) ) == diadicSize;// todo: more efficient
    }
    else
    {
        return m_leaf.lowerCorner().at(dir) == 0;
    }
}

template<typename T, short_t d, typename Z>
void gsHDomainBoundaryIterator<T,d,Z>::updateLeaf()
{
    const point & lower = m_leaf.lowerCorner();
    const point & upper = m_leaf.upperCorner();
    // gsDebug<<"leaf "<<  lower.transpose() <<", "
    //        << upper.transpose() <<"\n";

    const int level2 = m_leaf.level();

    // Update leaf box
    for (short_t dim = 0; dim < d; ++dim)
    {
        index_t start = lower(dim);
        index_t end  = upper(dim) ;

        const gsHTensorBasis<d,T> * hbasis = dynamic_cast<const gsHTensorBasis<d,T> * >(&m_basis);
        if (hbasis->manualLevels() )
        {
            static_cast<const gsHTensorBasis<d,T>*>(&m_basis)->
                _diadicIndexToKnotIndex(level2,dim,start);
            static_cast<const gsHTensorBasis<d,T>*>(&m_basis)->
                _diadicIndexToKnotIndex(level2,dim,end);
        }

        const gsKnotVector<T> & kv =
            static_cast<const gsHTensorBasis<d,T>*>(&m_basis)
            ->tensorLevel(level2).component(dim).knots();

        m_breaks[dim].clear();
        if ( dim == dir )
        {
            if ( par )
            {
                m_breaks[dim].push_back( kv(end-1) );
                m_breaks[dim].push_back( kv(end  ) );
            }
            else
            {
                m_breaks[dim].push_back( kv(start)   );
                m_breaks[dim].push_back( kv(start+1) );
            }
        }
        else
        {
            for (index_t index = start; index <= end; ++index)
                m_breaks[dim].push_back( kv(index) );// unique index
        }

        m_curElement(dim) =
        m_meshStart(dim)  = m_breaks[dim].begin();


        // for n breaks, we have n - 1 elements (spans)
        m_meshEnd(dim) =  m_breaks[dim].end() - 1;
    }
    this->m_pside.patch = m_patchId;
}

template<typename T, short_t d, typename Z>
const gsHTensorBasis<d,T> & gsHDomainBoundaryIterator<T,d,Z>::basis() const { return *static_cast<const gsHTensorBasis<d,T>*>(&m_basis); }

} // namespace gismo
