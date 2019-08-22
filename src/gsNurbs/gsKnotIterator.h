/** @file gsKnotIterator.h

    @brief Provides implementation of knot vector iterators.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mokris, A. Bressan
*/

#pragma once

namespace gismo {

namespace internal {

    template<typename T> class gsKnotIterator;

/**
   \brief A bi-directional knot iterator which provides extended
   information for the iterated knot

   The iteration is done over the unique knots in a knot-sequence, that is,
   the knots are visited **without** repetitions.
*/ 
template <typename T>
class gsUKnotIterator
{
public:
    friend class gsKnotVector<T>;
    friend class gsKnotIterator<T>;

    typedef typename gsKnotVector<T>::mult_t mult_t;    //index_t, gsKnotVector.h:85
    typedef std::random_access_iterator_tag iterator_category; 
    typedef const gsKnotVector<T> knotVector;
    typedef T value_type;
    typedef std::ptrdiff_t difference_type;
    typedef const T&       reference;
    typedef const T*       pointer;
    typedef const mult_t*  mltpointer;

private:
    mltpointer m_mlt ; ///< pointer to the beginning of the m_multSum sequence
    pointer    m_raw ; ///< pointer to the beginning of the m_repKnots sequence
    mult_t     m_upos; ///< unique index (without repetitions) of current knot

//#if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0
    mult_t m_dbg;// iteration limit: extra member for iterator debugging mode
//#endif

protected:

    // Change the value of the knot - access restricted to friend classes
    // Warning: call GISMO_ASSERT(check()) after using.
    void setValue(const T val)
    {
        std::fill(const_cast<T*>(m_raw) + firstAppearance(), 
                  const_cast<T*>(m_raw) + multSum(), val);
    }

public:

    /**
       \brief Default constructor (the iterator is initialized to
       NULL)
     */ 
    gsUKnotIterator()
    : m_mlt(NULL), m_raw(NULL), m_upos(0), m_dbg(0)
    { }

    /**
       \brief Constructs an iterator for the knot-vector \a KV.

       Optionally the iteration starts from from the knot with unique
       index (i.e. without repetitions) equal to \a upos
     */
    explicit gsUKnotIterator(knotVector & KV, const mult_t upos = 0)
    : m_mlt ( KV.multSumData() ),
      m_raw ( KV.data()        ),
      m_upos( upos             )
    {
        m_dbg = KV.uSize()+1;
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0
        //m_dbg = KV.uSize()+1;
        GISMO_ENSURE(upos < m_dbg, "Invalid iterator position "<< upos 
                     <<" for knot vector with "<<KV.uSize()<<" unique knots");
#       endif
    }

    /**
       \brief Returns an iterator pointing to the past-to-last
       position for the knot-vector \a KV
     */
    static inline gsUKnotIterator End(knotVector & KV)
    {   // the past-the-end position occurs for upos=KV.uSize()
        return gsUKnotIterator(KV, KV.uSize());
    }

public:

    /// \brief Dereferences the knot-iterator
    reference operator*  () const 
    {
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0 
        GISMO_ENSURE(m_upos >= 0 && m_upos + 1 < m_dbg, "Access to invalid knot position.");
#       endif 
        return m_raw[m_mlt[m_upos]-1];
    }

    pointer   operator-> () const {return m_raw+m_mlt[m_upos]-1 ;}

    gsUKnotIterator& operator++() 
    {
        ++m_upos;
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0
        GISMO_ENSURE(m_upos < m_dbg, "Invalid knot-iterator increment: "
                     << m_upos <<" is past the end position ("<<m_dbg-1<<").");
#       endif
        return *this;
    }

    gsUKnotIterator& operator--() 
    {
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0
        GISMO_ENSURE(m_upos > 0 && m_upos < m_dbg, "Invalid knot-iterator decrement");
#       endif
        --m_upos;
        return *this;
    }

    gsUKnotIterator operator++(int) { gsUKnotIterator tmp(*this); ++(*this); return tmp; }
    gsUKnotIterator operator--(int) { gsUKnotIterator tmp(*this); --(*this); return tmp; }

    /**
       \brief Equality operator for the iterator

       \note comparison between iterators on different (knot-)vectors
       is undefined behaviour, as of $ 24.2.5 of the C++ standard.
       However, this operator behaves as:
       \code
       this->uIndex() == other.uIndex();
       \endcode
       i.e. checks equality of the knot index positions (counted
       without repetitions) in two different knot-vectors.
     */
    bool operator == (const gsUKnotIterator& other) const 
    { return m_upos == other.m_upos;}// && m_raw == other.m_raw;}

    /**
       \brief Inequality operator for the iterator

       \note comparison between iterators on different (knot-)vectors
       is undefined behaviour, as of $ 24.2.5 of the C++ standard.
       However, this operator behaves as:
       \code
       this->uIndex() != other.uIndex();
       \endcode
       i.e. checks inequality of the knot index positions (counted
       without repetitions) in two different knot-vectors.
     */
    bool operator != (const gsUKnotIterator& other) const {return m_upos != other.m_upos;}
    bool operator <  (const gsUKnotIterator& other) const {return m_upos <  other.m_upos;}
    bool operator >  (const gsUKnotIterator& other) const {return m_upos >  other.m_upos;}
    bool operator <= (const gsUKnotIterator& other) const {return m_upos <= other.m_upos;}
    bool operator >= (const gsUKnotIterator& other) const {return m_upos >= other.m_upos;}

    /**
       \brief Returns the value of the knot which is \a a positions
       after the current knot (poistion and knots counted without repetitions).
     */
    reference operator [] (ptrdiff_t a) const
    {
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0
        GISMO_ENSURE(m_upos+a>=0 && m_upos+a+1 < m_dbg, "Invalid access to non-existent knot.");
#       endif
        return m_raw[m_mlt[m_upos+a]-1];
    }

    /**
       \brief Increments the iterator by \a a knot positions
       (counted on knot without repetitions)
    */
    gsUKnotIterator& operator+=(const difference_type & a)
    {
        m_upos += a;
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0 
        // Note: we allow invalid position for iterators on empty knot-vectors
        GISMO_ENSURE(m_dbg<2 || (m_upos >= 0 && m_upos < m_dbg),
                     "Iterator jumped to invalid knot position.");
#       endif
        return *this;
    }

    /**
       \brief Returns an iterator forwarded by \a a knot positions
       (counted on knots without repetitions)
     */
    gsUKnotIterator operator+(const difference_type & a) const
    {
        gsUKnotIterator tmp(*this);
        return tmp+=a;
    }

   /**
       \brief Decrements the iterator by \a a knot positions
       (counted without repetitions)
     */
    gsUKnotIterator& operator-=(const difference_type & a)
    {
        return operator+=(-a);
    }

    /**
       \brief Returns an iterator decremented by \a a knot positions
       (counted without repetitions)
     */
    gsUKnotIterator operator-(const difference_type & a) const
    {
        gsUKnotIterator tmp(*this);
        return tmp-=a;
    }

    friend difference_type operator-(const gsUKnotIterator & l, const gsUKnotIterator & r)
    {return l.m_upos - r.m_upos; }

    /**
       \brief Returns the value of the current knot
     */
    reference value() const 
    {
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0 
        GISMO_ENSURE(m_upos >= 0 && m_upos + 1< m_dbg, "Access to invalid knot position.");
#       endif 
        return this->operator*();
    }

    /**
       \brief Returns the multiplicity (i.e. the number of repetitions
       in the knot-sequence) of the current knot
     */
    mult_t multiplicity() const
    { 
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0 
        GISMO_ENSURE(m_upos >= 0 && m_upos + 1< m_dbg, "Access to invalid knot position.");
#       endif 
        if ( 0 == m_upos )//is it the first unique knot?
            return *m_mlt;
        else
        {
            mltpointer mp = m_mlt + m_upos;
            return *mp - *(mp-1);
        }
    }

    /**
       \brief Returns the index counted without repetitions (i.e. the
       unique index) of the current knot
     */
    mult_t uIndex() const {return m_upos;}

    /**
       \brief Returns the knot index (with repetitions) of the first
       knot in the knot sequence which is equal to the current knot
       value

       Use this function to get the position of the first appearance
       of the current knot
     */
     mult_t firstAppearance() const
    { 
        return 0 == m_upos ? 0 : m_mlt[m_upos-1];
    }

    /**
       \brief Returns the knot index of the last knot (with
       repetitions) in the knot sequence which is equal to the current
       knot value

       Use this function to get the position of the last appearance
       of the current knot
     */
    mult_t lastAppearance() const
    { 
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0 
        GISMO_ENSURE(m_upos >= 0 && m_upos + 1< m_dbg, "Access to invalid knot position.");
#       endif 
        return m_mlt[m_upos] - 1;
    }

    /**
       \brief Returns the first knot index of the next knot (counted with
       repetitions) in the knot sequence after the current knot

       Equivalently, this function returns the number of knot values
       (counted with multiplicitities) which are less or equal
       to the current knot.  Hence the name "multiplicity sum"
     */
    mult_t multSum() const
    { return m_mlt[m_upos];}
    
private:

    /*
       \brief Sets the iterator to the first (without repetitions)
       knot in the knot sequence

    // needed ?
    void reset()
    { 
        m_upos = 0;
    }
     */
};

/**
   \brief A bi-directional knot iterator which provides extended
   information for the iterated knot

   The iteration is done over all knots in a knot-sequence, that is,
   the knots are visited with repetitions.
*/ 
template <typename T>
class gsKnotIterator
{
public:
    friend class gsKnotVector<T>;

    typedef typename gsKnotVector<T>::mult_t mult_t;
    typedef std::random_access_iterator_tag iterator_category; 
    typedef const gsKnotVector<T> knotVector;
    typedef T value_type;
    typedef std::ptrdiff_t difference_type;
    typedef const T&       reference;
    typedef const T*       pointer;
    typedef const mult_t*  mltpointer;

private:
    gsUKnotIterator<T> m_uit; ///< iterator over unique knot-values
    mult_t             m_pos; ///< knot-index (with repetitions) of current knot

#if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0
    mult_t m_dbg;// iteration limit: extra member for iterator debugging mode
#endif

public:

    /**
       \brief Default constructor (the iterator is initialized to NULL)
    */ 
    gsKnotIterator()
    : m_uit(), m_pos(0)
    { }

    /**
       \brief Constructs an iterator for the knot-vector \a KV.

       Optionally the iteration starts from from the first appearance
       of the knot with unique index (i.e. without repetitions) equal to \a upos
    */
    explicit gsKnotIterator(knotVector & KV, const mult_t upos = 0)
    : m_uit(KV,upos), m_pos(firstAppearance())
    { 
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0
        m_dbg = KV.size()+1;
        GISMO_ENSURE(static_cast<size_t>(upos) <= KV.uSize(), 
                     "Invalid iterator position "<< upos 
                     <<" for knot vector with "<<KV.uSize()<<" unique knots");
#       endif
    }

    /**
       \brief Returns an iterator pointing to the past-to-last
       position for the knot-vector \a KV
     */
    static inline gsKnotIterator End(const gsKnotVector<T> & KV)
    {   // the past-the-end position occurs for upos=KV.uSize()
        return gsKnotIterator(KV, KV.uSize());
    }

public:

    /// \brief Dereferences the knot-iterator
    reference operator*  () const 
    {
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0 
        GISMO_ENSURE(m_pos >= 0 && m_pos + 1< m_dbg, 
                     "Access to invalid knot position.");
#       endif 
        return  m_uit.m_raw[m_pos];
    }

    pointer   get()         const {return  m_uit.m_raw+m_pos;}
    pointer   operator-> () const {return  get();}

    gsKnotIterator& operator++() 
    {
        if (++m_pos == m_uit.multSum())//crossing interval?
            ++m_uit;
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0
        GISMO_ENSURE(m_pos < m_dbg, "Invalid knot-iterator increment: "
                     << m_pos <<" is past the end position ("<<m_dbg-1<<").");
#       endif
        return *this;
    }

    gsKnotIterator& operator--() 
    {
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0
        GISMO_ENSURE(m_pos > 0 && m_pos < m_dbg, "Invalid knot-iterator decrement");
#       endif
        if ( m_pos-- == firstAppearance() )//crossing interval?
            --m_uit;
        return *this;
    }

    gsKnotIterator operator++(int) { gsKnotIterator tmp(*this); ++(*this); return tmp; }
    gsKnotIterator operator--(int) { gsKnotIterator tmp(*this); --(*this); return tmp; }

    /**
       \brief Equality operator for the iterator

       \note comparison between iterators on different (knot-)vectors
       is undefined behaviour, as of $ 24.2.5 of the C++ standard.
       However, this operator behaves as:
       \code
       this->index() == other.index();
       \endcode
       i.e. checks equality of the knot index positions (counted
       with repetitions) in two different knot-vectors.
     */
    bool operator == (const gsKnotIterator& other) const 
    { return m_pos == other.m_pos;}// && m_raw == other.m_raw;}

    /**
       \brief Inequality operator for the iterator

       \note comparison between iterators on different (knot-)vectors
       is undefined behaviour, as of $ 24.2.5 of the C++ standard.
       However, this operator behaves as:
       \code
       this->index() != other.index();
       \endcode
       i.e. checks inequality of the knot index positions (counted
       with repetitions) in two different knot-vectors.
     */
    bool operator != (const gsKnotIterator& other) const {return m_pos != other.m_pos;}
    bool operator <  (const gsKnotIterator& other) const {return m_pos <  other.m_pos;}
    bool operator >  (const gsKnotIterator& other) const {return m_pos >  other.m_pos;}
    bool operator <= (const gsKnotIterator& other) const {return m_pos <= other.m_pos;}
    bool operator >= (const gsKnotIterator& other) const {return m_pos >= other.m_pos;}

    /**
       \brief Returns the value of the knot which is \a a positions
       after the current knot (counted with repetitions).
     */
    reference operator [] (ptrdiff_t a) const
    {
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0
        GISMO_ENSURE(m_pos+a>=0 && m_pos+a+1 < m_dbg, "Invalid access to non-existent knot.");
#       endif
        return m_uit.m_raw[m_pos+a];
    }

    /**
       \brief Returns an iterator over the unique (without repetions)
       knot values of the knot sequence, which is set to the position
       of the current knot
     */
    const gsUKnotIterator<T> & uIterator() const
    { return m_uit; }

    /**
       \brief Sets the iterator to the position of the first
       appearance of the \em next unique knot value (ie. counted
       without repetitions) in the knot-sequence
     */
    void uNext()
    {
        m_pos = (m_uit++).multSum();//advance unique position and update m_pos
    }

    /**
       \brief Sets the iterator to the position of the first
       appearance of the \em previous unique knot value (ie. counted
       without repetitions) in the knot-sequence
     */
    void uPrev()
    {
        --m_uit;
        m_pos = firstAppearance();
    }

    /**
       \brief Increments the iterator by \a a knot positions
       (counted with repetitions)
     */
    gsKnotIterator& operator+=(const difference_type & a)
    {
        m_pos += a;
#       if defined(_GLIBCXX_DEBUG) || _SECURE_SCL != 0 
        GISMO_ENSURE((m_pos >= 0 && m_pos < m_dbg), 
                     "Iterator jumped to invalid knot position.");
#       endif 

        if (a<0) //substracting ?
        {
            mltpointer end = m_uit.m_mlt + m_uit.m_upos;
            mltpointer beg = end + a;
            if (beg < m_uit.m_mlt) beg = m_uit.m_mlt;
            //note: [beg, end) is a valid sorted range, complexity:  O(log a)
            m_uit.m_upos = std::upper_bound(beg, end, m_pos) - m_uit.m_mlt;
        }
        else    //incrementing
        {
            mltpointer beg = m_uit.m_mlt + m_uit.m_upos;
            // O(log a) efficient version
            mltpointer end = std::min(m_uit.m_mlt + m_uit.m_dbg-1, beg + a);
            m_uit.m_upos = std::upper_bound(beg, end, m_pos) - m_uit.m_mlt;

            // unsafe version withoud m_dbg
            //mltpointer end = beg + a; // note: could be over the end of m_mlt
            // while (beg!=end && (*beg)<=m_pos) { ++beg; }
            // m_uit.m_upos = beg - m_uit.m_mlt;
        }

        return *this;
    }

    /**
       \brief Returns an iterator forwarded by \a a knot positions
       (counted with repetitions)
     */
    gsKnotIterator operator+(const difference_type & a) const
    {
        gsKnotIterator tmp(*this);
        return tmp+=a;
    }


   /**
       \brief Decrements the iterator by \a a knot positions
       (counted with repetitions)
     */
    gsKnotIterator& operator-=(const difference_type & a)
    {
        return operator+=(-a);
    }

    /**
       \brief Returns an iterator decremented by \a a knot positions
       (counted with repetitions)
     */
    gsKnotIterator operator-(const difference_type & a) const
    {
        gsKnotIterator tmp(*this);
        return tmp-=a;
    }

    friend difference_type operator-(const gsKnotIterator & l, const gsKnotIterator & r)
    {return l.m_pos - r.m_pos; }

    /**
       \brief Returns the value of the current knot
     */
    reference value() const {return this->operator*();}

    /**
       \brief Returns the multiplicity (i.e. the number of repetitions
       in the knot-sequence) of the current knot
     */
    mult_t multiplicity() const
    { 
        return m_uit.multiplicity();
    }

    /**
       \brief Returns the knot-index (counted with repetitions) of the
       current knot
     */
    mult_t index() const {return m_pos;}

    /**
       \brief Returns the index counted without repetitions (i.e. the
       unique index) of the current knot
     */
    mult_t uIndex() const {return m_uit.m_upos;}

    /**
       \brief Returns the knot index of the first knot in the knot
       sequence which is equal to the current knot value

       Use this function to get the position of the first appearance
       of the current knot
     */
     mult_t firstAppearance() const
    { return m_uit.firstAppearance();}

    /**
       \brief Returns the knot index of the last knot in the knot
       sequence which is equal to the current knot value

       Use this function to get the position of the last appearance
       of the current knot
     */
    mult_t lastAppearance() const
    { return m_uit.lastAppearance();}

    /**
       \brief Returns the first knot index of the next knot (counted with
       repetitions) in the knot sequence after the current knot

       Equivalently, this function returns the number of knot values
       (counted with multiplicitities) which are less or equal
       to the current knot.  Hence the name "multiplicity sum"
     */
    mult_t multSum() const
    { return m_uit.multSum();}

    /**
       \brief Sets the iterator to the position of the first
       appearance of the current knot in the knot-sequence
     */
    void backToFirstAppearance()
    { m_pos = m_uit.firstAppearance();}

    /*
       \brief Sets the iterator to the first knot in the knot sequence
    // needed ?
    void reset()
    { 
        m_pos  = 0; 
        m_upos = 0;
    }
     */

};


} // namespace internal



}// namespace gismo

