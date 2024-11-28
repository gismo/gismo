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
    // typedef typename std::vector<T>::const_iterator  uiter;
    typename gsTensorBasis<D,T>::domainIter domainIter;

public:

    gsTensorDomain(const std::vector< typename gsKnotVector<T>::Ptr> & KVs)
    :
    m_knotVectors(KVs)
    {
        GISMO_ASSERT(KVs.size() == D, "Number of knot vectors must match the dimension of the domain.");
    }

    typename gsDomainIterator<T>::uPtr domainIterator(index_t i, const boxSide s = boundary::none) override
    {
        return ( s == boundary::none ?
                 domainIter(new gsTensorDomainIterator<T,D>(*this)) :
                 domainIter(new gsTensorDomainBoundaryIterator<T,D>(*this, s))
        );
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
        {
            result(0,i) = m_knotVectors[i]->first();
            result(1,i) = m_knotVectors[i]->last();
        }
        return result;
    }

    virtual gsMesh<T> mesh() const override
    {
        // gsMesh<T> mesh;
        // mesh.setDimension(d);
        // mesh.setBasis(m_basis);
        // return mesh;
    }

// Specific for gsTensorDomain
public:

    auto breaksBegin(index_t i)
    {
        return m_knotVectors[i]->uBegin();
    }

    auto breaksEnd(index_t i)
    {
        return m_knotVectors[i]->uEnd();
    }

protected:
    std::vector< typename gsKnotVector<T>::Ptr> m_knotVectors;

};

} // namespace gismo
