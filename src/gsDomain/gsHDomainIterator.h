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
                      const gsHTensorBasis<d,T> & basis,
                      const index_t patchId = -1);


    gsHDomainIterator(const gsHDomain<d,T,Z> & domain,
                      const gsHTensorBasis<d,T> & basis);

    gsHDomainIterator(const gsHDomainIterator & other);
    domainIter clone() const override;

    leafIterator init(const gsHTree<d,Z> & tree);

    // ---> Documentation in gsDomainIterator.h
    void next() override;

    // ---> Documentation in gsDomainIterator.h
    void next(index_t increment) override;

    /// Resets the iterator so that it can be used for another
    /// iteration through all boundary elements.
    void reset() override;

    gsVector<T> lowerCorner() const override;

    gsVector<T> upperCorner() const override;

    int getLevel() const;

    // Returns the element multi-index at the current level
    // If you need the element at the level above, divide this all indices by 2
    gsVector<index_t> elementMultiIndex() const;

private:

    gsHDomainIterator();

    /// returns true if there is a another leaf with a boundary element
    bool nextLeaf();

    /// Computes lower, upper and center point of the current element, maps the reference
    /// quadrature nodes and weights to the current element, and computes the
    /// active functions.
    void updateLeaf();

// =============================================================================
// members
// =============================================================================

public:
    // GISMO_DEPRECATED
    const gsHTensorBasis<d,T> & basis() const;

public:

#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen

private:

    const gsHTree<d,Z> & m_tree;
    const gsHTensorBasis<d,T> & m_basis;
    index_t m_patchId;

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

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHDomainIterator.hpp)
#endif
