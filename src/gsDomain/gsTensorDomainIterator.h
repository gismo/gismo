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

    // Documentation in gsDomainIterator.h
    void next(index_t increment) override
    {
        bool isGood(true);
        for (index_t i = 0; i < increment; i++)
            isGood = isGood && nextLexicographicIter(curElement, meshEnd);
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
