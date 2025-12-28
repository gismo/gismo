/** @file gsTensorDomainBoundaryIterator.h

    @brief Iterator over the boundary elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomainIterator.h>
#include <gsDomain/gsTensorDomain.h>
#include <gsUtils/gsCombinatorics.h>
#include <gsDomain/gsBreaksIterator.h>

namespace gismo
{

/**
 * @brief Re-implements gsDomainIterator for iteration over all elements of the boundary of a <b>tensor product</b> parameter domain.\n
 * <em>See gsDomainIterator for more detailed documentation and an example of the typical use!!!</em>
 *
 * \ingroup Tensor
 */

// Class which enables iteration over all elements of a tensor product parameter domain
// Documentation in gsDomainIterator.h

template<class T, int D, typename uiter>
class gsTensorDomainBoundaryIterator : public gsDomainIterator<T>
{
    typedef gsDomainIteratorWrapper<T> domainIterWrapper;
    typedef typename gsDomain<T>::uPtr domainPtr;
    using typename  gsDomainIterator<T>::uPtr;
    using typename  gsDomainIterator<T>::Ptr;
public:

    explicit gsTensorDomainBoundaryIterator(const gsTensorDomain<D,T> & domain,
                                            const boxSide & s);


    gsTensorDomainBoundaryIterator( const gsBasis<T>& b, const boxSide & s );

    // ---> Documentation in gsDomainIterator.h
    // proceed to the next element; returns true if end not reached yet
    void next() override;

    // ---> Documentation in gsDomainIterator.h
    // proceed to the next element (skipping #increment elements);
    // returns true if end not reached yet
    void next(index_t increment) override;

    /// Return the tensor index of the current element
    gsVector<unsigned, D> index() const;

    gsVector<T> lowerCorner() const override;

    gsVector<T> upperCorner() const override;

    const T getPerpendicularCellSize() const override;

    GISMO_DEPRECATED
    void adjacent( const gsVector<bool> & /* orient */,
                   gsDomainIterator<T>  & /* other */ ) override;

    /// Function to set the breakpoints in direction \a i manually
    void setBreaks(const std::vector<T> & newBreaks, index_t i); // i: direction


// Data members
private:

    // the dimension of the parameter space
    short_t d;

    // Boundary parameters
    short_t  dir;
    bool par;
    unsigned tindex;


    // First mesh-line on the tensor grid
    gsVector<domainIterWrapper, D> meshStart;

    // Last mesh-line on the tensor grid
    gsVector<domainIterWrapper, D> meshEnd;

    // Current element as pointers to it's supporting mesh-lines
    gsVector<domainIterWrapper, D> curElement;

public:
#   define Eigen gsEigen
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
#   undef Eigen
}; // class gsTensorDomainBoundaryIterator


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorDomainBoundaryIterator.hpp)
#endif
