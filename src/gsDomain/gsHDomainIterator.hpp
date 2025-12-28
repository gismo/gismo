/** @file gsHDomainIterator.hpp

    @brief Provides implementation of iterator of hierarchical domain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#pragma once

#include <gsDomain/gsHDomainIterator.h>

namespace gismo
{

template<typename T, short_t d, typename Z>
gsHDomainIterator<T,d,Z>::gsHDomainIterator(const gsHTree<d,Z> & tree,
                      const gsHTensorBasis<d,T> & basis,
                      const index_t patchId)
    :
    gsDomainIterator<T>(0, patchId),
    m_tree(tree),
    m_basis(basis),
    m_patchId(patchId)
{
    m_leaf = this->init(m_tree);
    updateLeaf();
}

template<typename T, short_t d, typename Z>
gsHDomainIterator<T,d,Z>::gsHDomainIterator(const gsHDomain<d,T,Z> & domain,
                      const gsHTensorBasis<d,T> & basis)
    :    gsHDomainIterator(domain.tree(),basis, domain.patchId())
{
}

template<typename T, short_t d, typename Z>
gsHDomainIterator<T,d,Z>::gsHDomainIterator(const gsHDomainIterator & other) = default;

template<typename T, short_t d, typename Z>
typename gsHDomainIterator<T,d,Z>::domainIter gsHDomainIterator<T,d,Z>::clone() const { return domainIter(new gsHDomainIterator(*this)); }

template<typename T, short_t d, typename Z>
typename gsHDomainIterator<T,d,Z>::leafIterator gsHDomainIterator<T,d,Z>::init(const gsHTree<d,Z> & tree)
{
    // Initialize mesh data
    m_meshStart.resize(d);
    m_meshEnd  .resize(d);

    // Initialize cell data
    m_curElement.resize(d);
    m_lower     .resize(d);
    m_upper     .resize(d);

    // Allocate breaks
    m_breaks = std::vector<std::vector<T> >(d, std::vector<T>());

    return tree.beginLeafIterator();
}

template<typename T, short_t d, typename Z>
void gsHDomainIterator<T,d,Z>::next()
{
    bool isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);
    if (!isGood) // went through all elements in m_leaf
        isGood = nextLeaf();
}

template<typename T, short_t d, typename Z>
void gsHDomainIterator<T,d,Z>::next(index_t increment)
{
    //todo: better implementation
    // compute the number of elements between curElement and meshEnd
    // use m_leaf.numElements() to skip leaves
    // arrive at the element or end
    bool isGood(m_leaf.good());
    for (index_t i = 0; i != increment && isGood; ++i)
    {
        isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);
        if (!isGood)
            isGood = nextLeaf();
    }
}

template<typename T, short_t d, typename Z>
void gsHDomainIterator<T,d,Z>::reset()
{
    m_leaf = m_tree.beginLeafIterator();
    updateLeaf();
}

template<typename T, short_t d, typename Z>
gsVector<T> gsHDomainIterator<T,d,Z>::lowerCorner() const
{
    gsVector<T> lower;
    lower.resize(d);
    for (short_t i = 0; i < d ; ++i)
        lower[i] = *m_curElement[i];
    return lower;
}

template<typename T, short_t d, typename Z>
gsVector<T> gsHDomainIterator<T,d,Z>::upperCorner() const
{
    gsVector<T> upper;
    upper.resize(d);
    for (short_t i = 0; i < d ; ++i)
        upper[i]  = *(m_curElement[i]+1);
    return upper;
}

template<typename T, short_t d, typename Z>
int gsHDomainIterator<T,d,Z>::getLevel() const
{
    return m_leaf.level();
}

template<typename T, short_t d, typename Z>
gsVector<index_t> gsHDomainIterator<T,d,Z>::elementMultiIndex() const
{
    gsVector<index_t> res(d);
    for (index_t i = 0; i!=d; ++i)
    {
        res[i] =  std::distance(m_breaks[i].begin(), m_curElement[i]);
    }
    return res;
}

template<typename T, short_t d, typename Z>
bool gsHDomainIterator<T,d,Z>::nextLeaf()
{
    bool isGood = m_leaf.next();

    if ( m_leaf.good() )
        updateLeaf();

    return isGood;
}

template<typename T, short_t d, typename Z>
void gsHDomainIterator<T,d,Z>::updateLeaf()
{
    const point & lower = m_leaf.lowerCorner();
    const point & upper = m_leaf.upperCorner();
    // gsDebug<<"leaf "<<  lower.transpose() <<", "
    //        << upper.transpose() <<"\n";

    const int level2 = m_leaf.level();

    // Update leaf box
    for (size_t dim = 0; dim < d; ++dim)
    {
        index_t start = lower(dim);
        index_t end  = upper(dim) ;

        if (basis().manualLevels() )
        {
            static_cast<const gsHTensorBasis<d,T>*>(&m_basis)->
                _diadicIndexToKnotIndex(level2,dim,start);
            static_cast<const gsHTensorBasis<d,T>*>(&m_basis)->
                _diadicIndexToKnotIndex(level2,dim,end);
        }

        const gsKnotVector<T> & kv =
            static_cast<const gsHTensorBasis<d,T>*>(&m_basis)
            ->tensorLevel(level2).component(dim).knots();

        // knotVals = kv.unique() 

        m_breaks[dim].clear();
        for (index_t index = start; index <= end; ++index)
            m_breaks[dim].push_back( kv(index) );// unique index

        m_curElement(dim) = 
        m_meshStart(dim)  = m_breaks[dim].begin();

        // for n breaks, we have n - 1 elements (spans)
        m_meshEnd(dim) =  m_breaks[dim].end() - 1;
    }
    this->m_pside.patch = m_patchId;
}

template<typename T, short_t d, typename Z>
const gsHTensorBasis<d,T> & gsHDomainIterator<T,d,Z>::basis() const { return *static_cast<const gsHTensorBasis<d,T>*>(&m_basis); }

} // namespace gismo
