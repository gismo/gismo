/** @file gsDomain.h

    @brief Abstracgt Base class representing a domain. i.e. a
    collection of elements (triangles, rectangles, cubes, simplices.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsBoundary.h>
#include <gsUtils/gsMesh/gsMesh.h>

namespace gismo
{

/**
    @brief Class representing a domain. i.e. a collection of
    elements (triangles, rectangles, cubes, simplices.

    \warning  This interface is under development and is not used yet...

    \ingroup Core
*/


/*

    TODO (later):
        /// Begin iterator (pointer)
    // virtual gsDomainIterator<T> begin() const
    // { gsWarn << "gsDomain: begin() not defined at "<< *this << "\n"; return gsDomainIterator<T>(); }

    /// End iterator (pointer)
    // virtual gsDomainIterator<T> end() const
    // { gsWarn << "gsDomain: end() not defined at "<< *this << "\n"; return gsDomainIterator<T>(); }


 */




template<class T>
class gsDomain
{
public:

    typedef typename memory::shared_ptr<gsDomain<T> > Ptr;
    typedef typename memory::unique_ptr<gsDomain<T> > uPtr;

    typedef gsDomainIteratorWrapper<T> domainIterWrapper;

    virtual ~gsDomain() { }

#if EIGEN_HAS_RVALUE_REFERENCES && EIGEN_GNUC_AT_MOST(4,7) && !EIGEN_COMP_PGI
    // defaulted declaration required at least in Gcc 4.7.2
    gsDomain() = default;
    gsDomain(const gsDomain&) = default;
    gsDomain(gsDomain&&) = default;
    gsDomain & operator=(const gsDomain&) = default;
    gsDomain & operator=(gsDomain&&) = default;
#endif

public:


    // iterator(index_t i)

    // numSubdomains (for pieces)

    // From Basis:
    // -[X] gsDomainIterator
    // -[X] numElements
    // -[ ] elementIndex

    // From gsDomainIterator
    // - side
    // - numElements


    /*
        // Taking domain from basis (uses gsHDomain WHICH DOES NOT INHERIT FROM GSDOMAIN)
        gsTHBSplineBasis<d,T> thb;
        gsExprAssembler<> A;
        A.setIntegrationElements(thb.domain());

        // Taking domain from multiBasis
        gsMultiBasis<T> mb;
        gsExprAssembler<> A;
        A.setIntegrationElements(mb.domain());

        // Taking domain from set of points
        gsMatrix<T> points;
        gsPointDomain<T> pd(points);
        A.setIntegrationElements(pd);

     */

public: // Domain element iterators

    virtual domainIterWrapper beginAll() const = 0;

    /** @brief Iterator over the elements of the domain.
     *
     * @param i The index of the domain.
     * @param s The side of the domain (optional).
     */
    // Default implementation for one patch
    virtual domainIterWrapper beginAt(index_t k, const boxSide bs = boundary::none) const
    {
        GISMO_ASSERT(0==k, "This is a one-piece domain.");
        GISMO_UNUSED(k);
        if (bs == boundary::none)
            return beginAll();
        else
            return beginBdr(bs);
    }

    virtual domainIterWrapper beginBdr (const boxSide   bs) const = 0;

    virtual domainIterWrapper endAll() const
    {
        return domainIterWrapper(new gsDomainIteratorEnd<T>(this->numElements()));
    }

    // Default implementation for one patch
    virtual domainIterWrapper endAt(index_t k, const boxSide bs = boundary::none) const
    {
        GISMO_ASSERT(0==k, "This is a one-piece domain.");
        GISMO_UNUSED(k);
        if (bs == boundary::none)
            return endAll();
        else
            return endBdr(bs);
    }

    // Default implementation for one patch
    virtual domainIterWrapper endBdr(const boxSide bs) const
    {
        return domainIterWrapper(new gsDomainIteratorEnd<T>(this->numElements(bs), bs));
    }

    // for multipatch
    virtual domainIterWrapper beginIfc(const boundaryInterface bi) const
    {GISMO_NO_IMPLEMENTATION}
    virtual domainIterWrapper endIfc  (const boundaryInterface bi) const
    {GISMO_NO_IMPLEMENTATION}

    /** @brief Dimension of the domain
    */
    virtual short_t dim() const
    {GISMO_NO_IMPLEMENTATION}

    /** @brief Number of elements in the domain
    */
    virtual size_t numElements(boxSide const & s = boundary::none) const
    {GISMO_NO_IMPLEMENTATION}

    /** @brief Bounding box of the domain
    */
    virtual gsMatrix<T> boundingBox() const
    {GISMO_NO_IMPLEMENTATION}

    /** @brief Mesh of the domain
    */
    virtual gsMesh<T> mesh() const
    {GISMO_NO_IMPLEMENTATION }

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
        os<<"Domain of dimennsion "<<dim()<<".";
        return os;
    }
}; // class gsDomain

/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsDomain<T>& b)
{ return b.print(os); }


} // namespace gismo
