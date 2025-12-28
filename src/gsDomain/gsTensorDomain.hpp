/** @file gsTensorDomain.hpp

    @brief Provides implementation of gsTensorDomain methods.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsDomain/gsTensorDomain.h>
#include <gsDomain/gsTensorSubDomain.h> // Keep this for gsTensorSubDomain references
#include <gsDomain/gsCompositeDomain.h> // Needed for gsCompositeDomain

namespace gismo
{

template<short_t d, class T>
gsTensorDomain<d,T>::gsTensorDomain(const gsKnotVector<T> & kvU)
{
    GISMO_ASSERT(d==1, "gsTensorDomain(const gsKnotVector<T> &) can only be used for 1D domains.");
    m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvU)));
}

template<short_t d, class T>
gsTensorDomain<d,T>::gsTensorDomain(const gsKnotVector<T> & kvU, const gsKnotVector<T> & kvV)
{
    GISMO_ASSERT(d==2, "gsTensorDomain(const gsKnotVector<T> &, const gsKnotVector<T> &) can only be used for 2D domains.");
    m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvU)));
    m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvV)));
}

template<short_t d, class T>
gsTensorDomain<d,T>::gsTensorDomain(const gsKnotVector<T> & kvU, const gsKnotVector<T> & kvV, const gsKnotVector<T> & kvW)
{
    GISMO_ASSERT(d==3, "gsTensorDomain(const gsKnotVector<T> &, ..., const gsKnotVector<T> &) can only be used for 3D domains.");
    m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvU)));
    m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvV)));
    m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvW)));
}

template<short_t d, class T>
gsTensorDomain<d,T>::gsTensorDomain(const std::vector<const gsKnotVector<T>*> & kvs)
{
    GISMO_ASSERT(kvs.size()==d, "Number of knot vectors must equal the dimension of the domain.");
    m_knotVectors.reserve(d);
    for (typename std::vector<const gsKnotVector<T>*>::const_iterator
         it = kvs.begin() ; it != kvs.end(); ++it )
        m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(**it)));
}

template<short_t d, class T>
gsTensorDomain<d,T>::gsTensorDomain(const std::vector<typename gsKnotVector<T>::Ptr> & kvs) : m_knotVectors(kvs)
{
    GISMO_ASSERT(kvs.size()==d, "Number of knot vectors must equal the dimension of the domain.");
}

template<short_t d, class T>
gsTensorDomain<d,T>::gsTensorDomain() : m_knotVectors(d) {}

template<short_t d, class T>
typename gsTensorDomain<d,T>::domainIter gsTensorDomain<d,T>::beginAll() const
{
    return domainIter(new gsTensorDomainIterator<T,d>(*this, this->patchId()));
}

template<short_t d, class T>
typename gsTensorDomain<d,T>::domainIter gsTensorDomain<d,T>::beginBdr(const boxSide btype) const
{
    return domainIter(new gsTensorDomainBoundaryIterator<T,d>(*this, btype));
}

template<short_t d, class T>
short_t gsTensorDomain<d,T>::dim() const { return d; }

template<short_t d, class T>
typename gsKnotVector<T>::Ptr gsTensorDomain<d,T>::knotVector(index_t i) const
{
    return m_knotVectors[i];
}

template<short_t d, class T>
typename gsKnotVector<T>::Ptr gsTensorDomain<d,T>::component(index_t i) const
{
    return m_knotVectors[i];
}

template<short_t d, class T>
typename gsDomain<T>::Ptr
gsTensorDomain<d,T>::decompose(size_t npieces) const
{
    // Define the parameter ranges
    std::vector<typename gsTensorSubDomain<d,T>::Range> fullRanges;
    for (short_t dim = 0; dim < d; ++dim) 
    {   
        fullRanges.push_back(typename gsTensorSubDomain<d,T>::Range(0, m_knotVectors[dim]->numElements()));
    }
    return gsTensorSubDomain<d,T>::decompose(*this, fullRanges, npieces, this->patchId());
}

template<short_t d, class T>
size_t gsTensorDomain<d,T>::numElements() const
{
    size_t nElem = 1;
    for (short_t dim = 0; dim < d; ++dim)
        nElem *= m_knotVectors[dim]->numElements();
    return nElem;
}

} // namespace gismo
