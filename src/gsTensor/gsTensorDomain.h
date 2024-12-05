/** @file gsTensorDomain.h

    @brief Iterator over the elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDomain.h>
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
    typedef typename gsDomainIterator<T>::uPtr domainIter;

public:

    // gsTensorDomain(const std::vector<gsKnotVector<T> *> & KVs)
    // :
    // {
    //     GISMO_ASSERT(KVs.size() == D, "Number of knot vectors must match the dimension of the domain.");
    //     m_knotVectors.reserve(D);
    //     for (index_t i = 0; i < D; ++i)
    //         m_knotVectors.push_back(memory::make_shared_not_owned(KVs[i]));
    // }

    gsTensorDomain(const std::vector<typename gsDomain<T>::Ptr> & KVs)
    :
    m_knotVectors(give(KVs))
    {
        GISMO_ASSERT(KVs.size() == D, "Number of domains must match the dimension of the domain.");
    }

    typename gsDomainIterator<T>::uPtr begin(index_t i, const boxSide s = boundary::none) const override
    {
        return ( s == boundary::none ?
                 domainIter(new gsTensorDomainIterator<T,D>(*this)) :
                 domainIter(new gsTensorDomainBoundaryIterator<T,D>(*this, s))
        );
    }

    typename gsDomainIterator<T>::uPtr end(index_t i, const boxSide s = boundary::none) const override
    {
        return domainIter(new gsDomainIteratorEnd<T>(this->numElements(),s));
    }

    // Look at gsBasis class for a description
    size_t numElements(boxSide const & s = boundary::none) const override
    {
        const short_t dir =  s.direction();
        size_t nElem = 1;
        for (short_t dim = 0; dim < D; ++dim)
        {
            if(dim == dir && s!=boundary::none)
                continue;
            nElem *= m_knotVectors[dim]->numElements();
        }
        return nElem;
    }

    short_t dim() const override { return D; }

    gsMatrix<T> boundingBox() const override
    {
        gsMatrix<T> result(2,D);
        for (short_t i = 0; i < D; ++i)
            result.col(i) = m_knotVectors[i]->boundingBox();
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



    auto component(index_t i)
    {
        return m_knotVectors[i];
    }

    // auto breaksEnd(index_t i)
    // {
    //     return m_knotVectors[i]->uEnd();
    // }

protected:
    // NOTE: change vector to array?
    std::vector< typename gsDomain<T>::Ptr> m_knotVectors;

};

} // namespace gismo
