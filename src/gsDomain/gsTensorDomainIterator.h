/** @file gsTensorDomainIterator.h

    @brief Iterator over the elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomainIterator.h>
#include <gsDomain/gsTensorDomain.h>
#include <gsUtils/gsCombinatorics.h>

namespace gismo
{
// Documentation in gsDomainIterator.h
// Class which enables iteration over all elements of a tensor product parameter domain

/**
 * @brief Re-implements gsDomainIterator for iteration over all elements of a <b>tensor product</b> parameter domain.\n
 * <em>See gsDomainIterator for more detailed documentation and an example of the typical use!!!</em>
 *
 * \ingroup Tensor
 */

template<class T, int D>
class gsTensorDomainIterator : public gsDomainIterator<T>
{
private:
    typedef typename gsDomainIterator<T>::uPtr domainIter;
    typedef gsDomainIteratorWrapper<T> domainIterWrapper;

public:

    explicit gsTensorDomainIterator(const gsTensorDomain<T,D> & domain)
    : gsDomainIterator<T>()
    {
        this->setPatchIndex(domain.patchIndex());

        // compute breaks and mesh size
        // meshStart.resize(D);
        // meshEnd.resize(D);
        // curElement.resize(D);

        for (int i=0; i < D; ++i)
        {
            meshEnd[i]    = give(domain.component(i)->endAll()  );
            meshStart[i]  = give(domain.component(i)->beginAll());
            curElement[i] = give(domain.component(i)->beginAll());
        }
    }

    gsTensorDomainIterator(const gsTensorDomainIterator & other) = default;
    domainIter clone() const override { return domainIter(new gsTensorDomainIterator(*this)); }

    // Documentation in gsDomainIterator.h
    void next() override
    {
        nextLexicographicIter(curElement, meshEnd);
    }

    /**
     * @brief Advances the iterator by \a increment elements in O(D) arithmetic.
     *
     * Elements are ordered lexicographically with direction 0 fastest, so the
     * element with component ids c_i has linear index sum_i c_i * prod_{m<i} n_m,
     * where n_i is the number of elements in direction i. The iterator moves to
     * the element of linear index (current + increment).
     *
     * If that index is not smaller than the number of elements N = prod_i n_i,
     * the iterator is set to the past-the-end state (0,...,0,n_{D-1}), which is
     * the state reached by the first lexicographic step that leaves the domain.
     * A further next(index_t) from this state leaves it unchanged; a single-step
     * next() does not stop here and must not be called on it.
     *
     * An \a increment <= 0 does nothing.
     *
     * Precondition: the current state is a valid element or the past-the-end
     * state, as obtained by construction, reset(), next() or next(index_t).
     *
     * Complexity: O(D) integer arithmetic plus one reset() and one advance of the
     * component iterators per direction (O(1) each for knot-span components).
     */
    void next(index_t increment) override
    {
        if (increment <= 0) return;

        // 64-bit arithmetic: the number of elements can exceed the range of index_t
        uint64_t stride[D+1];
        stride[0] = 1;
        uint64_t lin = 0;
        for (int i = 0; i < D; ++i)
        {
            const uint64_t n = static_cast<uint64_t>(meshEnd[i].id() - meshStart[i].id());
            lin += static_cast<uint64_t>(curElement[i].id() - meshStart[i].id()) * stride[i];
            stride[i+1] = stride[i] * n;
        }
        const uint64_t target = lin + static_cast<uint64_t>(increment);

        if (target >= stride[D])
        {
            for (int i = 0; i < D; ++i)
                curElement[i].reset();
            curElement[D-1] += static_cast<index_t>(meshEnd[D-1].id() - meshStart[D-1].id());
            return;
        }

        for (int i = 0; i < D; ++i)
        {
            const uint64_t n = stride[i+1] / stride[i];
            const index_t q = static_cast<index_t>((target / stride[i]) % n);
            curElement[i].reset();
            curElement[i] += q;
        }
    }

    // Documentation in gsDomainIterator.h
    void reset() override
    {
        for (index_t i = 0; i < D; ++i)
            curElement[i].reset();
    }

    /// return the tensor index of the current element
    gsVector<unsigned, D> index() const
    {
        gsVector<unsigned, D> curr_index(D);
        for (int i = 0; i < D; ++i)
            curr_index[i]  = curElement[i]->index();
        return curr_index;
    }

    void getVertices(gsMatrix<T>& result)
    {
        result.resize( D, 1 << D);

        const gsVector<T> lower = lowerCorner();
        const gsVector<T> upper = upperCorner();
        gsVector<T,D> v, l, u;
        l.setZero();
        u.setOnes();
        v.setZero();
        int r = 0;
        do {
            for ( int i = 0; i< D; ++i)
                result(i,r) = ( v[i] ? upper[i] : lower[i] );
        }
        while ( nextCubeVertex(v, l, u) );
    }

    gsVector<T> lowerCorner() const override
    {
        gsVector<T> lower(D);
        for (short_t i = 0; i < D ; ++i)
            lower[i]  = curElement[i].lowerCorner().value();
        return lower;
    }

    gsVector<T> upperCorner() const override
    {
        gsVector<T> upper(D);
        for (short_t i = 0; i < D ; ++i)
            upper[i]  = curElement[i].upperCorner().value();
        return upper;
    }

    bool isBoundaryElement() const override
    {
        for (int i = 0; i< D; ++i)
            if ( curElement[i].isBoundaryElement() )
                return true;
        return false;
    }

    index_t domainDim() const {return D;}


//    size_t numElements() const
//    {
//
//    }

// Data members
private:
    // Extent of the tensor grid and current element as pointers to
    // it's supporting mesh-lines
    gsVector<domainIterWrapper, D> meshStart, meshEnd, curElement;

public:
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen
}; // class gsTensorDomainIterator

} // namespace gismo
