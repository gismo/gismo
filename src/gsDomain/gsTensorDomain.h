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

public:
    // constructors

    gsTensorDomain(const gsKnotVector<T> & kvU)
    {
        GISMO_ASSERT(D==1, "gsTensorDomain(const gsKnotVector<T> &) can only be used for 1D domains.");
        m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvU)));
    }

    gsTensorDomain(const gsKnotVector<T> & kvU, const gsKnotVector<T> & kvV)
    {
        GISMO_ASSERT(D==2, "gsTensorDomain(const gsKnotVector<T> &, const gsKnotVector<T> &) can only be used for 2D domains.");
        m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvU)));
        m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvV)));
    }

    gsTensorDomain(const gsKnotVector<T> & kvU, const gsKnotVector<T> & kvV, const gsKnotVector<T> & kvW)
    {
        GISMO_ASSERT(D==3, "gsTensorDomain(const gsKnotVector<T> &, ..., const gsKnotVector<T> &) can only be used for 3D domains.");
        m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvU)));
        m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvV)));
        m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(kvW)));
    }

    gsTensorDomain(const std::vector<const gsKnotVector<T>*> & kvs)
    {
        GISMO_ASSERT(kvs.size()==D, "Number of knot vectors must equal the dimension of the domain.");
        m_knotVectors.reserve(D);
        for (typename std::vector<const gsKnotVector<T>*>::const_iterator
             it = kvs.begin() ; it != kvs.end(); ++it )
            m_knotVectors.push_back(memory::make_shared(new gsKnotVector<T>(**it)));
    }

    gsTensorDomain(const std::vector<typename gsKnotVector<T>::Ptr> & kvs) : m_knotVectors(kvs)
    {
        GISMO_ASSERT(kvs.size()==D, "Number of knot vectors must equal the dimension of the domain.");
    }

    /// Default constructor
    gsTensorDomain() : m_knotVectors(D) {}


    // iterators

    virtual domainIter beginAll() const override
    {
        return domainIter(new gsTensorDomainIterator<T,D>(*this, this->patchId()));
    }

    /// Provides a new iterator which loops over boundary elements only.
    ///
    /// The iterator will run over boundary elements of type \a btype. If the
    /// given side is a corner of the domain, the boundary elements are cells
    /// of (dim-1)-dimensional boundaries. For an edge, these cells are
    /// 1-dimensional entities and for a volume they are 2-dimensional.
    virtual domainIter beginBdr(const boxSide btype = boundary::all) const override
    {
        return domainIter(new gsTensorDomainBoundaryIterator<T,D>(*this, btype));
    }


public: // more members

    short_t dim() const override { return D; }

    // Look at gsBasis class for a description

    virtual size_t numElements() const override;

// Specific for gsTensorDomain
public:
    typename gsKnotVector<T>::Ptr knotVector(index_t i) const
    {
        return m_knotVectors[i];
    }

    typename gsKnotVector<T>::Ptr component(index_t i) const
    {
        return m_knotVectors[i];
    }


    /// Specialized domain decomposition for tensor domains.
    /// Creates a smooth partitioning using axis-aligned boxes (tensor subdomain ranges).
    /// \param npieces Number of pieces to decompose into
    /// \return A composite domain made of tensor subdomains
    virtual typename gsDomain<T>::Ptr decompose(index_t npieces) const override;

    /// \brief Decomposes the domain into a given number of pieces using a specific strategy.
    /// \param npieces Number of pieces to decompose into.
    /// \param strategy The decomposition strategy to use.
    /// \return A shared pointer to the decomposed domain.
    virtual typename gsDomain<T>::Ptr decompose(index_t npieces, decompositionStrategy strategy) const override;

protected:
    // NOTE: change vector to array?
    std::vector< typename gsKnotVector<T>::Ptr > m_knotVectors;
}; // class gsTensorDomain

} // namespace gismo


#include <gsDomain/gsTensorDomain.hpp>
