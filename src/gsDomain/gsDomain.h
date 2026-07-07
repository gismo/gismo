/** @file gsDomain.h

    @brief Abstract Base class representing a domain. i.e. a
    collection of elements (triangles, rectangles, cubes, simplices.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsDomain/gsDomainIterator.h>
#include <gsCore/gsBoundary.h>
#include <gsUtils/gsMesh/gsMesh.h>

namespace gismo
{

/**
    @brief Class representing a domain. i.e. a collection of elements
    (triangles, rectangles, cubes, simplices.

    gsDomain<> dom;

    // Iterate over all elements
    for (auto it = dom.beginAll(); it!=dom.endAll(); ++it)

    // Iterate over all elements of patch k
    for (auto it = dom.subdomain(k).beginAll(); it!=dom.subdomain(k).endAll(); ++it)

    // Iterate over all elements of the boundary of the domain
    for (auto it = dom.beginBdr(); it!=dom.endBdr(); ++it)

    // Iterate over all elements of the boundary of the subdomain k
    for (auto it = dom.subdomain(k).beginBdr(); it!=dom.subdomain(k).endBdr(); ++it)

    // Iterate over all elements of the boundary side bs of the subdomain k
    for (auto it = dom.subdomain(k).beginBdr(bs); it!=dom.subdomain(k).endBdr(bs); ++it)

    // Iterate over all elements of all the interfaces of the domain
    for (auto it = dom.beginIfc(); it!=dom.endIfc(); ++it)

    // Iterate over all elements of the interface \a bi
    for (auto it = dom.beginIfc(bi); it!=dom.endIfc(bi); ++it)

    // Number of elements of

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

    typedef gsDomainIteratorWrapper<T> iterator;

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

    /// Helper class for range-based for loops
    class ElementRange
    {
        iterator begin_, end_;
    public:
        ElementRange(iterator _begin, iterator _end)
        : begin_(give(_begin)), end_(give(_end)) { }
        iterator & begin() { return begin_; }
        iterator & end()   { return end_;   }
    };

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

public:

    /// Return the k-th subdomain, in the case that there are more than one
    virtual Ptr subdomain(index_t k) const
    {
        GISMO_ASSERT(0==k, "This is a single-piece domain.");
        return memory::make_shared_not_owned(this);
    }

    virtual size_t nPieces() const { return 1; }


public: // Domain element iterators

    /// Returns iterator over the elements in this domain
    virtual iterator beginAll() const = 0;

    /// Returns iterator at the past-to-end element in this domain
    virtual iterator endAll() const
    {
        return iterator(new gsDomainIteratorEnd<T>(this->numElements()));
    }

    /// Returns an iterator over the boundary.
    /// special value \a all: iterate over all boundaries
    virtual iterator beginBdr (const boxSide bs = boundary::all) const
    {GISMO_UNUSED(bs); GISMO_NO_IMPLEMENTATION}

    /// Returns an iterator to the end of the boundary elements
    /// special value \a all: iterate over all boundaries
    virtual iterator endBdr(const boxSide bs = boundary::all) const
    {
        return iterator(new gsDomainIteratorEnd<T>(this->numElementsBdr(bs)));
    }

    /// Returns a pair of two iterators, that define a chunk (range)
    /// per thread of elements that can be used in a parallel for
    /// loop. If a single thread is available then it returns the pair
    /// ( beginAll(), endAll() )
    inline ElementRange allElements() const
    {
#ifdef _OPENMP
        const int num_threads = omp_get_num_threads();
        const int num_elem    = numElements();
        const int chunk_size  = (num_elem + num_threads - 1) / num_threads;
        int chunk_start = chunk_size * omp_get_thread_num();
        if (chunk_start < num_elem)
        {
            iterator domIt = beginAll();
            domIt += chunk_start;
            chunk_start += chunk_size;
            iterator domItEnd = iterator(new gsDomainIteratorEnd<T>(
                     (chunk_start < num_elem) ? chunk_start : num_elem));
            return ElementRange( give(domIt), give(domItEnd) );
        }
        else // no work for this thread
        {
            iterator domIt = iterator(new gsDomainIteratorEnd<T>(0));
            return ElementRange( domIt, domIt );
        }
#else
        return ElementRange( beginAll(), endAll() );
#endif
    }

    // TO DO:
    //GISMO_ELEMENT_RANGE_LOOP(allElements, numElements, beginAll, endAll)
    //GISMO_ELEMENT_RANGE_LOOP(bdrElements, elementsBdr, beginBdr, endBdr)
    //GISMO_ELEMENT_RANGE_LOOP(ifcElements, elementsIfc, beginIfc, endIfc)

    // for multipatch
    virtual iterator beginIfc(const boundaryInterface bi) const
    {GISMO_UNUSED(bi); GISMO_NO_IMPLEMENTATION}
    virtual iterator endIfc  (const boundaryInterface bi) const
    {GISMO_UNUSED(bi); GISMO_NO_IMPLEMENTATION}

    /** @brief Number of elements in the domain
    */
    virtual size_t numElements() const = 0;

    /** @brief Number of elements in the domain
     */
    virtual size_t numElementsBdr(boxSide const & s = boundary::all) const
    {GISMO_UNUSED(s); GISMO_NO_IMPLEMENTATION}

    // NOTE: for immersed
    //virtual size_t numBackgroundElements() const;

    /** @brief Degree of the domain
    */
    virtual short_t degree(short_t i = 0) const
    {GISMO_UNUSED(i); GISMO_NO_IMPLEMENTATION}

    /** @brief Dimension of the domain
    */
    virtual short_t dim() const
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
        os<<"Domain of dimension "<<dim()<<", "<< "number of elements: "<< numElements()<<"\n";
        return os;
    }
}; // class gsDomain

/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsDomain<T>& b)
{ return b.print(os); }


} // namespace gismo
