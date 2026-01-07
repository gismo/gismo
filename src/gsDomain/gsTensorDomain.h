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
#include <gsDomain/gsTensorDomainIterator.h>
#include <gsDomain/gsTensorDomainBoundaryIterator.h>

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

template<short_t d, class T>
class gsTensorDomain : public gsDomain<T>
{
private:
    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef typename gsKnotVector<T>::const_uiterator knotIter;

public:
    // constructors

    gsTensorDomain(const gsKnotVector<T> & kvU);

    gsTensorDomain(const gsKnotVector<T> & kvU, const gsKnotVector<T> & kvV);

    gsTensorDomain(const gsKnotVector<T> & kvU, const gsKnotVector<T> & kvV, const gsKnotVector<T> & kvW);

    gsTensorDomain(const std::vector<const gsKnotVector<T>*> & kvs);

    gsTensorDomain(const std::vector<typename gsKnotVector<T>::Ptr> & kvs);

    /// Default constructor
    gsTensorDomain();


    // iterators

    virtual domainIter beginAll() const override;

    /// Provides a new iterator which loops over boundary elements only.
    ///
    /// The iterator will run over boundary elements of type \a btype. If the
    /// given side is a corner of the domain, the boundary elements are cells
    /// of (dim-1)-dimensional boundaries. For an edge, these cells are
    /// 1-dimensional entities and for a volume they are 2-dimensional.
    virtual domainIter beginBdr(const boxSide btype = boundary::all) const override;


public: // more members

    short_t degree(short_t i = 0) const override;

    short_t dim() const override;

    // Look at gsBasis class for a description
    virtual size_t numElements() const override;

    // Look at gsBasis class for a description
    size_t numElementsBdr(boxSide const & s = boundary::none) const override;

// Specific for gsTensorDomain
public:
    typename gsKnotVector<T>::Ptr knotVector(index_t i) const;

    typename gsKnotVector<T>::Ptr component(index_t i) const;

    /// \brief Decomposes the domain into a given number of pieces using a specific strategy.
    /// \param npieces Number of pieces to decompose into.
    /// \return A shared pointer to the decomposed domain.
    virtual typename gsDomain<T>::Ptr decompose(size_t npieces) const override;

protected:
    // NOTE: change vector to array?
    std::vector< typename gsKnotVector<T>::Ptr > m_knotVectors;
}; // class gsTensorDomain

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorDomain.hpp)
#endif