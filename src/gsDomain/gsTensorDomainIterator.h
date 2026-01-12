/** @file gsTensorDomainIterator.h

    @brief Iterator over the elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h> // Centralized forward declarations
#include <gsDomain/gsDomainIterator.h>
#include <gsDomain/gsTensorDomain.h> // Full definition of gsTensorDomain needed
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

    explicit gsTensorDomainIterator(const gsTensorDomain<D,T> & domain, index_t patchId = -1);

    gsTensorDomainIterator(const gsTensorDomainIterator & other);
    domainIter clone() const override;

    // Documentation in gsDomainIterator.h
    void next() override;

    // Documentation in gsDomainIterator.h
    void next(index_t increment) override;

    // Documentation in gsDomainIterator.h
    void reset() override;

    /// return the tensor index of the current element
    gsVector<unsigned, D> index() const;

    void getVertices(gsMatrix<T>& result);

    gsVector<T> lowerCorner() const override;

    gsVector<T> upperCorner() const override;

    bool isBoundaryElement() const override;

    index_t domainDim() const;

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

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorDomainIterator.hpp)
#endif
