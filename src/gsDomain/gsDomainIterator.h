/** @file gsDomainIterator.h

    @brief Provides declaration of DomainIterator abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

//#include <gsCore/gsBasis.h> // todo: remove
#include <gsCore/gsBoundary.h>


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
    uPtr m_domainIter;//change as Ptr

public:
    explicit gsDomainIteratorWrapper(gsDomainIterator<T> * _itptr = nullptr) : m_domainIter(_itptr)
    { }

    explicit gsDomainIteratorWrapper(uPtr _iter) : m_domainIter(give(_iter))
    { }

    gsDomainIteratorWrapper(const gsDomainIteratorWrapper & _other)
    {
        this->operator=(_other);
    }

#if EIGEN_HAS_RVALUE_REFERENCES
    /// Move constructor
    gsDomainIteratorWrapper(gsDomainIteratorWrapper && _other)
    : m_domainIter(give(_other.m_domainIter))
    { }

    /// Assignment operator
    gsDomainIteratorWrapper& operator= ( const gsDomainIteratorWrapper& _other )
    {
        m_domainIter = _other.m_domainIter->clone();
        return *this;
    }

    /// Move assignment operator
    gsDomainIteratorWrapper & operator=(gsDomainIteratorWrapper && _other) EIGEN_NOEXCEPT
    {
        m_domainIter = give(_other.m_domainIter);
        return *this;
    }
#else
    /// Assignment operator (uses copy-and-swap idiom)
    gsDomainIteratorWrapper& operator= ( gsDomainIteratorWrapper _other )
    {
        std::swap(m_domainIter, _other.m_domainIter);
        return *this;
    }
#endif

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

    /// Inequality operator to compare two iterators
    bool operator<(const gsDomainIteratorWrapper& other) const
    {
        return m_domainIter->id() < other.id();
    }

    /// Increment operator to proceed to the next element
    gsDomainIteratorWrapper& operator++()
    {
        m_domainIter->next();
        m_domainIter->nextId();
        return *this;
    }

    /// Post-increment operator to proceed to the next element
    size_t operator++(int)
    {
        size_t oid = id();
        ++(*this);
        return oid;
    }

    /// Decrement operator to proceed to the next element
    gsDomainIteratorWrapper& operator--()
    {
        m_domainIter->prev();
        m_domainIter->prevId();
        return *this;
    }

    /// Increment inplace by a number of steps
    gsDomainIteratorWrapper& operator+=(index_t k)
    {
        m_domainIter->next(k);
        m_domainIter->nextId(k);
        return *this;
    }

    /// Increment by a number of steps
    gsDomainIteratorWrapper operator+(index_t k) const
    {
        gsDomainIteratorWrapper result(m_domainIter->clone());
        result += k;
        return result;
    }

    /// Decrement operator to proceed to the next element
    gsDomainIteratorWrapper& operator-=(index_t k)
    {
        m_domainIter->prev(k);
        m_domainIter->prevId(k);
        return *this;
    }

    /// Decrement by a number of steps
    // gsDomainIteratorWrapper operator-(index_t k) const
    // {
    //     gsDomainIteratorWrapper result(m_domainIter->clone());
    //     result -= k;
    //     return result;
    // }

    /// Difference operator
    size_t operator-(const gsDomainIteratorWrapper& other) const
    { return id() - other.id(); }

    void reset()
    {
        m_domainIter->reset();
        m_domainIter->resetId();
    }

    gsDomainIterator<T> * get() { return m_domainIter.get(); }

    gsDomainIterator<T> & operator*() { return *m_domainIter; }

public:

    short_t dim() const
    { return this->centerPoint().rows(); }

    gsVector<T> lowerCorner() const
    { return m_domainIter->lowerCorner(); }

    gsVector<T> upperCorner() const
    { return m_domainIter->upperCorner(); }

    gsVector<T> centerPoint() const
    { return m_domainIter->centerPoint(); }

    const T getPerpendicularCellSize() const
    { return m_domainIter->getPerpendicularCellSize(); }

    bool isBoundaryElement() const
    { return m_domainIter->isBoundaryElement(); }

    T getCellSize() const
    { return m_domainIter->getCellSize(); }

    T getMinCellLength() const
    { return m_domainIter->getMinCellLength(); }

    T getMaxCellLength() const
    { return m_domainIter->getMaxCellLength(); }

    T volume() const
    { return m_domainIter->volume(); }

public:
    /// Returns the element id
    inline size_t id() const { return m_domainIter->id(); }

    inline boxSide side() const {return m_domainIter->side();}

    inline index_t patch() const {return m_domainIter->patch();}

    inline index_t localId() const {return m_domainIter->localId();}
    
    /// Fetches data of integer type based on string label
    const index_t & label(const std::string & _label)
    {return m_domainIter->label(_label); }

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

    explicit gsDomainIterator(index_t _id = 0, const boxSide & _bs = boundary::none)
    : m_id(_id), m_pside(0,_bs)
    { }

    virtual ~gsDomainIterator() { }

    virtual uPtr clone() const { GISMO_NO_IMPLEMENTATION }

    //void setPatch(index_t k) { m_pside.patch = k; }

private:

    /** @brief Proceeds to the next element.
     *
     * The function returns true if there are still elements remaining that have not been treated.\n
     * For the typical usage of this function, see the example in the
     * documentation of gsDomainIterator.
     *
     * @todo: remove the return value
     */

    virtual void next() = 0;

    virtual void prev() { GISMO_NO_IMPLEMENTATION }

    /// \brief Proceeds to the next element (skipping \p increment elements).
    virtual void next(index_t increment)
    {
        for(index_t i = 0; i!=increment; ++i)
            this->next();
    }

    virtual void prev(index_t decrement)
    {
        for(index_t i = 0; i!=decrement; ++i)
            this->prev();
    }

    /// Resets the iterator so that it points to the first element
    virtual void reset()
    {
        GISMO_NO_IMPLEMENTATION
        //*this = give(*m_domain->beginAll().get());
    }

protected:
    inline void resetId  () { m_id = 0;}
    inline void nextId(index_t _k = 1) { m_id += _k; }
    inline void prevId(index_t _k = 1) { m_id -= _k; }

public:

    /// Returns the element id -- see also patch() for the patch index
    size_t id() const   { return m_id; }

    /// Returns the local element id -- e.g. the id inside the patch
    virtual size_t localId() const { return m_id; }

    /// Return dimension of the elements
    short_t dim() const   { return centerPoint().size(); }

    /// fetches data of integer type based on string label
    virtual const index_t & label(const std::string & _label)
    {GISMO_ERROR("Cannot find property "<< _label); }

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
    gsVector<T> centerPoint () const
    {
        return (this->lowerCorner()+this->upperCorner()).array()/(T)2;
    }

    /// \brief Returns the lower corner of the current element.
    ///
    /// The current element is a <em>d</em>-dimensional hypercube.
    /// The coordinates of its lower corner is returned as a gsVector of length \a d.\n
    /// \n
    /// E.g., if the current two-dimensional element is defined by <em>[a,b]x[c,d]</em>, then <em>[a,c]</em> is returned (see also upperCorner()).
    virtual gsVector<T> lowerCorner() const
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// \brief Returns the upper corner of the current element.
    ///
    /// The current element is a <em>d</em>-dimensional hypercube.
    /// The coordinates of its upper corner is returned as a gsVector of length \a d.\n
    /// \n
    /// E.g., if the current two-dimensional element is defined by <em>[a,b]x[c,d]</em>, then <em>[b,d]</em> is returned (see also lowerCorner()).
    virtual gsVector<T> upperCorner() const
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

    inline boxSide side() const {return m_pside.side();}
    inline index_t patch() const {return m_pside.patch;}
    inline index_t & patch() {return m_pside.patch;}
protected:


    size_t m_id; ///< The element ID

    patchSide m_pside; ///< The patch side, when applicable

protected:
    // disable copying
    gsDomainIterator( const gsDomainIterator& ) = default;
    gsDomainIterator& operator= ( const gsDomainIterator& ) = default;

}; // class gsDomainIterator


template <class T>
class gsDomainIteratorEnd : public gsDomainIterator<T>
{
    typedef memory::unique_ptr< gsDomainIterator<T> > uPtr;
public:

    explicit gsDomainIteratorEnd(size_t id) : gsDomainIterator<T>(id) { }

    uPtr clone() const override { return uPtr(new gsDomainIteratorEnd(this->m_id)); }

    virtual void next() override
    { GISMO_ERROR("Cannot proceed to next element. End iterator reached."); }

    virtual void next(index_t /* increment */) override
    { GISMO_ERROR("Cannot proceed to next element. End iterator reached."); }

    virtual void prev() override
    {
        GISMO_NO_IMPLEMENTATION
    }

    virtual void prev(index_t /* decrement */) override
    {
        GISMO_NO_IMPLEMENTATION
    }


};

/// Print (as string) operator to be used by all derived classes
//template<class T>
//std::ostream &operator<<(std::ostream &os, const gsDomainIterator<T>& b)
//{return b.print(os); }


} // namespace gismo
