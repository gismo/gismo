/** @file gsTensorDomain.h

    @brief Iterator over the elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomain.h>
#include <gsNurbs/gsKnotVector.h>

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
class gsTensorDomain : public gsDomain<T>
{
private:
    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef typename gsKnotVector<T>::const_uiterator knotIter;

public: // constructors

    gsTensorDomain(const std::vector<typename gsDomain<T>::Ptr> & KVs)
    :
    m_knotVectors(give(KVs))
    {
        GISMO_ASSERT(KVs.size() == D, "Number of domains must match the dimension of the domain.");
    }

public: // iterators

    virtual domainIter beginAll() const override
    {
        return domainIter(new gsTensorDomainIterator<T,D>(*this));
    }

    domainIter beginBdr(const boxSide bs) const override
    { return domainIter(new gsTensorDomainBoundaryIterator<T,D,knotIter>(*this, bs)); }

public: // more members

    // Look at gsBasis class for a description
    size_t numElements() const override
    {
        size_t nElem = 1;
        for (short_t dim = 0; dim < D; ++dim)
            nElem *= m_knotVectors[dim]->numElements();
        return nElem;
    }

        // Look at gsBasis class for a description
    size_t numElementsBdr(boxSide const & s = boundary::none) const override
    {
        if(s==boundary::none)
        {
            GISMO_NO_IMPLEMENTATION
        }

        const short_t dir =  s.direction();
        size_t nElem = 1;
        for (short_t dim = 0; dim < D; ++dim)
        {
            if(dim == dir)
                continue;
            nElem *= m_knotVectors[dim]->numElements();
        }
        return nElem;
    }

    short_t degree(short_t i) const override
    {
        return m_knotVectors[i]->degree();
    }

    short_t dim() const override { return D; }

    gsMatrix<T> boundingBox() const override
    {
        gsMatrix<T> result(D, 2);
        for (short_t i = 0; i < D; ++i)
            result.row(i) = m_knotVectors[i]->boundingBox();
        return result;
    }

    virtual gsMesh<T> mesh() const override
    {
        // gsMesh<T> mesh;
        // mesh.setDimension(d);
        // mesh.setBasis(m_basis);
        // return mesh;
        GISMO_NO_IMPLEMENTATION
    }

// Specific for gsTensorDomain
public:

    typename gsDomain<T>::Ptr component(index_t i) const
    {
        return m_knotVectors[i];
    }

protected:
    // NOTE: change vector to array?
    std::vector< typename gsDomain<T>::Ptr> m_knotVectors;

};

} // namespace gismo
