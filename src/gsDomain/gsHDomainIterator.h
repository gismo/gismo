/** @file gsHDomainIterator.h

    @brief Provides declaration of iterator of hierarchical domain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/

#pragma once

#include <gsDomain/gsHTree.h>
#include <gsDomain/gsHDomain.h>
#include <gsDomain/gsKdNode.h>

#include <gsDomain/gsDomainIterator.h>

namespace gismo
{

// Documentation in gsDomainIterator
/** @brief Re-implements gsDomainIterator for iteration over all boundary
  * elements of a <b>hierarchical</b> parameter domain.
  *
  * <em>See
  * gsDomainIterator for more detailed documentation and an example of
  * the typical use!!!</em>\n Used, e.g., for basis of classes
  * gsHTensorBasis or gsTHBSplineBasis.
  *
  * \ingroup HSplines
  */

template<typename T, short_t d, typename Z>
class gsHDomainIterator: public gsDomainIterator<T>
{
public:

    typedef gsKdNode<d,Z> node;

    typedef typename node::point point;

    typedef typename std::vector<T>::const_iterator  uiter;

    typedef gsHTree<d,Z> hDomain;

    typedef typename hDomain::const_literator leafIterator;

    typedef typename gsDomainIterator<T>::uPtr domainIter;
public:

    gsHDomainIterator(const gsHTree<d,Z> & tree,
                      const gsHTensorBasis<d,T> & basis)
    :
    gsDomainIterator<T>(),
    m_tree(tree),
    m_basis(basis)
    {
        m_leaf = this->init(m_tree);
        updateLeaf();
    }


    gsHDomainIterator(const gsHDomain<d,T,Z> & domain,
                      const gsHTensorBasis<d,T> & basis)
    :
    gsHDomainIterator(domain.tree(),basis)
    {
    }

    gsHDomainIterator(const gsHDomainIterator & other) = default;
    domainIter clone() const override { return domainIter(new gsHDomainIterator(*this)); }

    leafIterator init(const gsHTree<d,Z> & tree)
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

    // ---> Documentation in gsDomainIterator.h
    void next() override
    {
        bool isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);
        if (!isGood) // went through all elements in m_leaf
            isGood = nextLeaf();
    }

    // ---> Documentation in gsDomainIterator.h
    void next(index_t increment) override
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

    /// Resets the iterator so that it can be used for another
    /// iteration through all boundary elements.
    void reset() override
    {
        m_leaf = m_tree.beginLeafIterator();
        updateLeaf();
    }

    gsVector<T> lowerCorner() const override
    {
        gsVector<T> lower;
        lower.resize(d);
        for (short_t i = 0; i < d ; ++i)
            lower[i] = *m_curElement[i];
        return lower;
    }

    gsVector<T> upperCorner() const override
    {
        gsVector<T> upper;
        upper.resize(d);
        for (short_t i = 0; i < d ; ++i)
            upper[i]  = *(m_curElement[i]+1);
        return upper;
    }

    int getLevel() const
    {
        return m_leaf.level();
    }

    // Returns the element multi-index at the current level
    // If you need the element at the level above, divide this all indices by 2
    gsVector<index_t> elementMultiIndex() const
    {
        gsVector<index_t> res(d);
        for (index_t i = 0; i!=d; ++i)
        {
            res[i] =  std::distance(m_breaks[i].begin(), m_curElement[i]);
        }
        return res;
    }

private:

    gsHDomainIterator();

    /// returns true if there is a another leaf with a boundary element
    bool nextLeaf()
    {
        bool isGood = m_leaf.next();

        if ( m_leaf.good() )
            updateLeaf();

        return isGood;
    }

    /// Computes lower, upper and center point of the current element, maps the reference
    /// quadrature nodes and weights to the current element, and computes the
    /// active functions.
    void updateLeaf()
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
    }

// =============================================================================
// members
// =============================================================================

public:
    // GISMO_DEPRECATED
    const gsHTensorBasis<d,T> & basis() const { return *static_cast<const gsHTensorBasis<d,T>*>(&m_basis); }

public:

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

private:

    const gsHTree<d,Z> & m_tree;
    const gsHTensorBasis<d,T> & m_basis;

    // The current leaf node of the tree
    leafIterator m_leaf;

    // Coordinates of the grid cell boundaries
    // \todo remove this member
    std::vector< std::vector<T> > m_breaks;

    // Extent of the tensor grid
    gsVector<uiter, d> m_meshStart, m_meshEnd;

    // Current element as pointers to it's supporting mesh-lines
    gsVector<uiter, d> m_curElement;

    // parameter coordinates of current grid cell
    gsVector<T> m_lower, m_upper;
};

} // end namespace gismo
