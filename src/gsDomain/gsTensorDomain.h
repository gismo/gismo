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
#include <gsDomain/gsTensorDomainFaceIterator.h>
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

    domainIter beginSkeleton() const override
    { return domainIter(new gsTensorDomainFaceIterator<T,D,AllFaces>(breakGrid())); }

    size_t numSkeletonFaces() const override
    { return gsTensorDomainFaceIterator<T,D,AllFaces>::numFaces(breakGrid()); }

    /// A tensor domain carries no trimming information, so its ghost set is empty.
    domainIter beginGhost() const override
    { return domainIter(new gsDomainIteratorEnd<T>(0)); }

    size_t numGhostFaces() const override { return 0; }

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

private:
    /// Per-direction element boundaries, breaks[j].size() == numElements in j + 1
    std::vector< std::vector<T> > breakGrid() const
    {
        std::vector< std::vector<T> > breaks(D);
        for (short_t i = 0; i < D; ++i)
        {
            const gsKnotVector<T> * kv =
                dynamic_cast<const gsKnotVector<T>*>(m_knotVectors[i].get());
            GISMO_ENSURE(nullptr!=kv,
                         "gsTensorDomain: face iteration requires gsKnotVector components.");
            breaks[i] = kv->breaks();
        }
        return breaks;
    }

protected:
    // NOTE: change vector to array?
    std::vector< typename gsDomain<T>::Ptr> m_knotVectors;

};

} // namespace gismo
