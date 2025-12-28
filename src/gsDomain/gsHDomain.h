/** @file gsHDomain.h

    @brief Provides declaration of the HDomain class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): G. Kiss, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h> 
#include <gsDomain/gsDomain.h> 
#include <gsDomain/gsHTree.h> 
#include <gsDomain/gsHDomainLeafIter.h> 
#include <gsDomain/gsHDomainIterator.h> 
#include <gsDomain/gsHDomainBoundaryIterator.h> 
#include <gsDomain/gsCompositeDomain.h>
#include <gsDomain/gsTensorSubDomain.h>

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

template <class T>
struct LeafInfo 
{
    gsVector<index_t> lower;
    gsVector<index_t> upper;
    index_t global_id; // Keep original global ID for mapping if needed
    short_t level;     // Keep level for info
    size_t num_elements; // Number of elements in this leaf
};


template<short_t d, class T, class Z>
class gsHDomain : public gsDomain<T>
{
public:

    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef typename gsHTree<d,Z>::const_literator leafIterator;

    template <class _T, short_t _d, class _Z>
    friend class gsHDomainIterator;
    template <class _T, short_t _d, class _Z>
    friend class gsHDomainBoundaryIterator;

private: 
    const gsHTree<d,Z> & m_tree;
    const gsHTensorBasis<d,T> & m_basis;

    bool leafOnBoundary(const boxSide & s, const leafIterator leaf) const;

public:

    explicit gsHDomain(const gsHTree<d,Z>& tree,
                       const gsHTensorBasis<d,T>& basis);

    domainIter beginAll() const override;

    domainIter beginBdr(const boxSide bs) const override;

    size_t numElements() const override;

    size_t numElementsBdr(boxSide const & s = boundary::none) const override;

    short_t degree(short_t i) const override;

    short_t dim() const override;

    gsMatrix<T> boundingBox() const override;

    const gsHTree<d,Z> & tree() const;

    virtual typename gsDomain<T>::Ptr decompose(size_t npieces) const override;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHDomain.hpp)
#endif