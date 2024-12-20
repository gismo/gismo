/** @file gsDomainIterator.h

    @brief Provides declaration of DomainIterator abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither
*/

#pragma once

//#include <gsCore/gsBasis.h> // todo: remove
#include <gsCore/gsBoundary.h>


template<class T> class gsDomain;

namespace gismo
{

/** The gsDomainIterator \n
 * \brief Class which enables iteration over all elements of a parameter domain.
 *
 *
 * It also includes some additional functionality which is typically used
 * when assembling the system matrix for numerically solving a PDE.
 *
 * - <b>Iteration through the elements:</b>\n
 * The function next() jumps to the "next" element and sets up the quadrature
 * nodes and weights on that element.
 * The specific implementation
 * of how to select the next element depends on the structure of the underlying mesh.\n
 * The function good() indicates whether there still is a "next" element to be found.
 *
 *
 * Note that the features of the gsDomainIterator strongly depend on the underlying basis.
 * Hence the gsBasis is given as an input argument to the constructor.
 *
 * An example of the typical use of gsDomainIterator (remark: replace
 * the constructor by the constructor of the actually used derived
 * class):
 *
 * \verbatim
     gsDomainIterator domIter( basis );         // constructor

     for (; domIter.good(); domIter.next() )    // loop over all elements
     {
         // Your source code using
         domIter.centerPoint();
         domIter.lowerCorner();
         domIter.upperCorner();

     }
     \endverbatim

     \ingroup Core
 *
 *
 *
 */

template <class T>
class gsDomainIteratorWrapper
{
    typedef memory::unique_ptr< gsDomainIterator<T> > uPtr;
    uPtr m_domainIter;

public:
    gsDomainIteratorWrapper(uPtr _iter) : m_domainIter(give(_iter))
    { }

    gsDomainIteratorWrapper(gsDomainIterator<T> * _itptr) : m_domainIter(_itptr)
    { }

    gsDomainIteratorWrapper(gsDomainIteratorWrapper && _other) : m_domainIter(give(_other.m_domainIter))
    { }

    gsDomainIteratorWrapper & operator=(gsDomainIteratorWrapper && _other)
    {
        m_domainIter = give(_other.m_domainIter);
        return *this;
    }

    /// Equality operator to compare two iterators
    bool operator==(const gsDomainIteratorWrapper& other) const
    {
        return m_domainIter->id() == other.id();
    }

    /// Inequality operator to compare two iterators
    bool operator!=(const gsDomainIteratorWrapper& other) const
    {
        return m_domainIter->id() != other.id();
    }

    /// Increment operator to proceed to the next element
    gsDomainIteratorWrapper& operator++()
    {
        m_domainIter->next();
        m_domainIter->nextId();
        return *this;
    }

    /// Decrement operator to proceed to the next element
    gsDomainIteratorWrapper& operator--()
    {
        m_domainIter->prev();
        m_domainIter->prevId();
        return *this;
    }

    /// Incremental operator to proceed to the next element
    gsDomainIteratorWrapper& operator+(index_t k)
    {
        m_domainIter->next(k);
        m_domainIter->nextId(k);
        return *this;
    }

    /// Decrement operator to proceed to the next element
    gsDomainIteratorWrapper& operator-(index_t k)
    {
        m_domainIter->prev(k);
        m_domainIter->prevId(k);
        return *this;
    }

    /// Equality operator to compare two iterators
    size_t operator-(const gsDomainIteratorWrapper& other) const
    { return id() - other.id(); }

    void reset()
    {
        m_domainIter->reset();
        m_domainIter->resetId();
    }
public:

    const gsVector<T>& lowerCorner() const
    { return m_domainIter->lowerCorner(); }

    const gsVector<T>& upperCorner() const
    { return m_domainIter->upperCorner(); }

public:
    /// Returns the element id
    inline size_t id() const { return m_domainIter->id(); }

    inline boxSide side() const {return m_domainIter->m_side;}
};


template <class T>
class gsDomainIterator
{
    friend class gsDomainIteratorWrapper<T>;

public:
    /// Shared pointer for gsDomainIterator
    typedef memory::shared_ptr< gsDomainIterator > Ptr;
    /// Unique pointer for gsDomainIterator
    typedef memory::unique_ptr< gsDomainIterator > uPtr;

public:

    gsDomainIterator(index_t _id = 0) : m_basis(NULL), m_isGood( true ), m_id(_id) { }

    /// \brief Constructor using a basis
    gsDomainIterator( const gsDomain<T>& _dom, const boxSide & _bs = boundary::none)
    : m_domain(&_dom), m_side(_bs),
      m_isGood( true ), m_id(0)
    { }

    virtual ~gsDomainIterator() { }
    
private:

    /** @brief Proceeds to the next element.
     *
     * The function returns true if there are still elements remaining that have not been treated.\n
     * For the typical usage of this function, see the example in the
     * documentation of gsDomainIterator.
     */
    virtual bool next() = 0;

    /// \brief Proceeds to the next element (skipping \p increment elements).
    virtual bool next(index_t increment) = 0;

    virtual bool prev() { GISMO_NO_IMPLEMENTATION };

    virtual bool prev(index_t decrement) { GISMO_NO_IMPLEMENTATION };

    /// Resets the iterator so that it points to the first element
    virtual void reset()
    {
        GISMO_NO_IMPLEMENTATION
    }

protected:
    inline void resetId  () { m_id = 0;}
    inline void nextId(index_t _k = 1) { m_id += _k; }
    inline void prevId(index_t _k = 1) { m_id -= _k; }

public:

    /// Returns the element id
    size_t id() const   { return m_id; }

    /// Is the iterator still pointing to a valid element?
    // \todo use !=end() instead
    GISMO_DEPRECATED
    bool good() const   { return m_isGood; }

    /// Return dimension of the elements
    short_t dim() const   { return center.size(); }

    /// Updates \a other with and adjacent element
    /// \todo upgrade to return adjacent range instead
    // GISMO_DEPRECATED??
    virtual void adjacent( const gsVector<bool> & ,
                           gsDomainIterator & )
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// \brief Returns the center of the current element.
    ///
    /// The current element is a <em>d</em>-dimensional hypercube.
    /// The coordinates of its upper corner is returned as a gsVector of length \a d.\n
    /// \n
    /// E.g., if the current two-dimensional element is defined by <em>[a,b]x[c,d]</em>, then <em>[b,d]</em> is returned (see also lowerCorner()).
    const gsVector<T>& centerPoint () const
    { return center; }

    /// \brief Returns the lower corner of the current element.
    ///
    /// The current element is a <em>d</em>-dimensional hypercube.
    /// The coordinates of its lower corner is returned as a gsVector of length \a d.\n
    /// \n
    /// E.g., if the current two-dimensional element is defined by <em>[a,b]x[c,d]</em>, then <em>[a,c]</em> is returned (see also upperCorner()).
    virtual const gsVector<T>& lowerCorner() const
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// \brief Returns the upper corner of the current element.
    ///
    /// The current element is a <em>d</em>-dimensional hypercube.
    /// The coordinates of its upper corner is returned as a gsVector of length \a d.\n
    /// \n
    /// E.g., if the current two-dimensional element is defined by <em>[a,b]x[c,d]</em>, then <em>[b,d]</em> is returned (see also lowerCorner()).
    virtual const gsVector<T>& upperCorner() const
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// \brief Returns the perpendicular cell size of boundary iterator.
    ///
    /// Only works for boundary iterators. Returns the length from
    /// the boundary side to the parallel side not on the boundary.
    virtual const T getPerpendicularCellSize() const
    {
        GISMO_NO_IMPLEMENTATION
    }

    virtual bool isBoundaryElement() const
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Return the diagonal of the element
    T getCellSize() const
    {
        return (upperCorner() - lowerCorner()).norm();
    }

    /// Return the length of the smallest edge of the element
    T getMinCellLength() const
    {
        return (upperCorner() - lowerCorner()).minCoeff();
    }

    /// Return the length of the largest edge of the element
    T getMaxCellLength() const
    {
        return (upperCorner() - lowerCorner()).maxCoeff();
    }

    /// Return the volume of the element
    T volume() const
    { return (upperCorner() - lowerCorner()).prod(); }

    /// Returns the number of elements. --REMOVE 
    virtual size_t numElements() const
    {
        //\todo Remove this implementation. Probably using a shallow
        //copy, "reset" and "next" would do this better.

        // Buggy, and probably a terrible implementation,
        // but needed and therefore can be useful
        // sometimes.
        return m_domain->numElements();
    }


    inline boxSide side() const {return m_side;}

protected:

    const gsDomain<T> * m_domain;

    // \todo patchSide
    boxSide m_side;

private:
    size_t m_id;

protected:
    
    //// REMOVE

    /// Coordinates of a central point in the element (in the parameter domain).
    gsVector<T> center;


    /// Flag indicating whether the domain iterator is "good". If it
    /// is "good", the iterator can continue to the next element.
    bool m_isGood;

    /// The basis on which the domain iterator is defined.
    const gsBasis<T> * m_basis;

private:
    // disable copying
    gsDomainIterator( const gsDomainIterator& );
    gsDomainIterator& operator= ( const gsDomainIterator& );
}; // class gsDomainIterator


template <class T>
class gsDomainIteratorEnd : public gsDomainIterator<T>
{
    using gsDomainIterator<T>::m_side;

public:

    explicit gsDomainIteratorEnd(size_t id, boxSide s = boundary::none)
    :
    gsDomainIterator<T>(id)
    {
        m_side = s;
    }

    virtual bool next() override { GISMO_ERROR("Cannot proceed to next element. End iterator reached."); }

    virtual bool next(index_t increment) override { GISMO_ERROR("Cannot proceed to next element. End iterator reached."); }

};

/// Print (as string) operator to be used by all derived classes
//template<class T>
//std::ostream &operator<<(std::ostream &os, const gsDomainIterator<T>& b)
//{return b.print(os); }


} // namespace gismo
