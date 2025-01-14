/** @file gsHDomain.h

    @brief Provides declaration of the HDomain class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/

#pragma once

#include <gsHSplines/gsHTree.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsHSplines/gsHDomainLeafIter.h>
#include <gsHSplines/gsHDomainIterator.h>
#include <gsHSplines/gsHDomainBoundaryIterator.h>
#include <gsCore/gsBoundary.h>

namespace gismo
{

template<typename T, short_t d, typename Z>
class gsHDomainIterator;

template<class T>
class gsSegment;

// template <short_t d, class Z>
// class gsKdNode;

template <class T>
class gsVSegment;

/**
\brief
Class with a <em>hierarchical domain structure</em> represented by a box
k-d-tree


The hierarchical domain structure represets a sequence of domains \f$\Omega^\ell  \subset \Omega\f$ which are nested in the sense that
\f[
\Omega = \Omega^0 \supset \Omega^1 \supset \Omega^2 \supset \ldots \supset \Omega^N \supset \Omega^{N+1} = \emptyset \f]
Each subdomain \f$\Omega^\ell\f$ is a (not necessarily connected) collection of axis-aligned elements/cells.

\remark In the context of HB-splines and THB-splines, these elements/cells are the knot spans of the tensor-product mesh at the respective level.

The information on the hierarchical domains is stored in a k-d-tree, where each leaf represents an axis-aligned box
\f$\omega \subseteq \Omega\f$, such that
\f$\omega \subseteq \Omega^\ell \land \omega \cap \Omega^{\ell+1} = \emptyset\f$ (i.e., each leaf of the tree can be assiciated with exactly one level of the hierarchy).

The implementation is, up to some technical differences, based on the technique described in the following publication

- G. Kiss, C. Giannelli, and B. Juettler.
Algorithms and data structures for truncated hierarchical B-splines. In M. Floater et al., editors, Mathematical Methods for Curves and Surfaces, volume 8177, pages 304-323.
Lecture Notes in Computer Science, 2014.

also available as a technical report

- G. Kiss, C. Giannelli, and B. Juettler.
Algorithms and Data Structures for
Truncated Hierarchical B-splines
DK-Report No. 2012-14, Doctoral Program Computational Mathematics: Numerical Anaylsis and Symbolic Computation, 2012.

Regarding the mentioned technical differences: A binary tree is used instead of a quad-tree (which was discussed in the above-mentioned publications). Also, the domains are not necessarily split at their middle, but according to the position of the domain of the next level.

Template parameters
\param d is the dimension
\param Z is the box-index type

\ingroup HSplines
*/


template<short_t d, class T, class Z = index_t>
class gsHDomain : public gsDomain<T> // is template correct?
{
public:

    typedef gsDomainIteratorWrapper<T> domainIter;

    template <class _T, short_t _d, class _Z>
    friend class gsHDomainIterator;

public:

    gsHDomain(const gsHTensorBasis<d,T> & basis)
    :
    m_tree(basis.tree())
    {
    }

    gsHDomain(const gsHTree<d,Z>& tree)
    :
    m_tree(tree)
    {
    }

    domainIter beginAll() const
    {
        return domainIter(new gsHDomainIterator<T,d,Z>(m_tree));
    }

    domainIter beginBdr(const boxSide bs) const override
    {
        return domainIter(new gsHDomainBoundaryIterator<T,d,Z>(m_tree, bs));
    }

    const gsHTree<d,Z> & tree() const { return m_tree; }
          gsHTree<d,Z> & tree()       { return m_tree; }

protected:
    const gsHTree<d,Z> & m_tree;

};

}// end namespace gismo

