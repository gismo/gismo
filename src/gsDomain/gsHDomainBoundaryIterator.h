/** @file gsHDomainBoundaryIterator.h

    @brief Provides declaration of iterator on boundary of hierarchical basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsHTree.h>
#include <gsDomain/gsHDomain.h>
#include <gsDomain/gsKdNode.h>
#include <gsHSplines/gsHTensorBasis.h>

#include <gsDomain/gsDomainIterator.h>

namespace gismo
{

// Documentation in gsDomainIterator
/**
  * @brief
  * Re-implements gsDomainIterator for iteration over all boundary
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
class gsHDomainBoundaryIterator: public gsDomainIterator<T>
{
public:

    typedef gsKdNode<d,Z> node;

    typedef typename node::point point;

    typedef typename std::vector<T>::const_iterator  uiter;

    typedef gsHTree<d,Z> hDomain;

    typedef typename hDomain::const_literator leafIterator;

public:

    gsHDomainBoundaryIterator(const gsHTree<d,Z> & tree,
                              const gsHTensorBasis<d,T> & basis,
                              const boxSide & s)
    :
    gsDomainIterator<T>(),
    m_tree(tree),
    m_basis(basis)
    {
        init(m_tree,s);
    }

    gsHDomainBoundaryIterator(const gsHDomain<d,T,Z> & domain,
                              const gsHTensorBasis<d,T> & basis,
                              const boxSide & s)
    :
    gsHDomainBoundaryIterator(domain.tree(),basis,s)
    {
    }

    void init(const gsHTree<d,Z> & tree, const boxSide & s)
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
        bool isGood(true);
        for (index_t i = 0; i < increment; i++)
            isGood = nextLexicographic(m_curElement, m_meshStart, m_meshEnd);

        if (isGood) // went through all elements in m_leaf
            isGood = nextLeaf();
    }

    /// Resets the iterator so that it can be used for another
    /// iteration through all boundary elements.
    void reset() override
    {
        initLeaf(m_tree);
    }

    gsVector<T> lowerCorner() const override
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

    gsVector<T> upperCorner() const override
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

    const T getPerpendicularCellSize() const override
    {
        return *(m_curElement[dir]+1) - *m_curElement[dir];
    }

    int getLevel() const
    {
        return m_leaf.level();
    }

private:

    gsHDomainBoundaryIterator();

    /// Navigates to the first leaf on our side
    void initLeaf(const hDomain & tree_domain)
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


    /// returns true if there is a another leaf with a boundary element
    bool nextLeaf()
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

    /// returns true if the current leaf is on our side
    bool leafOnBoundary() const
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

    // Boundary parameters
    short_t dir; // direction normal to the boundary
    bool par;     // parameter value

    // The current leaf node of the tree
    leafIterator m_leaf;

    // Coordinates of the grid cell boundaries
    // \todo remove this member
    std::vector< std::vector<T> > m_breaks;

    // Extent of the tensor grid
    gsVector<uiter, d> m_meshStart, m_meshEnd;

    // Current element as pointers to it's supporting mesh-lines
    gsVector<uiter, d> m_curElement;
};

} // end namespace gismo
