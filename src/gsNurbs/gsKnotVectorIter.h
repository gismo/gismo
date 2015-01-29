/** @file gsKnotVectorIter.h

    @brief Provides declaration of knot vector iterator.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <iterator>

#include <gsCore/gsTemplateTools.h>

namespace gismo
{

template<class T>
class gsKnotVector;

  /** 
      @brief Iterator over the unique knot-values of a gsKnotVector.
      
      \ingroup Nurbs
  */
template <class T, bool isconst = false> 
class gsKnotVectorIter {
// Definitions
public:
    typedef std::random_access_iterator_tag iterator_category; 
    typedef T value_type;
    typedef std::ptrdiff_t difference_type;
    typedef typename choose<isconst, const T&, T&>::type reference;
    typedef typename choose<isconst, const T*, T*>::type pointer;
    typedef typename choose<isconst, typename gsKnotVector<T>::const_iterator,    
			    typename gsKnotVector<T>::iterator >::type knotptr;
    typedef typename choose<isconst, const gsKnotVector<T>, 
			    gsKnotVector<T> >::type containertype;

    friend class gsKnotVectorIter<T,true>;

// Constructors
public:
    gsKnotVectorIter()
    { }

    gsKnotVectorIter( containertype & c , bool start = true) 
        : m_current ( (start ? c.begin() : c.end() ) ),
          m_begin(c.begin()),
          m_end  (c.end()  )
    { 
            ignoreDublicates();
    }
    
    gsKnotVectorIter(const gsKnotVectorIter<T, false>& it) 
        : m_current(it.m_current), 
          m_begin(it.m_begin ),
          m_end(it.m_end)
    { 

    }
    
// Accessors
public:
    reference operator*() const { return *m_current; }
    pointer  operator->() const { return  m_current; }
    reference operator[](const int& rhs) {return *(*this+rhs);}

// Arithmetic
public:
    gsKnotVectorIter& operator=(value_type* rhs) 
    {m_current = rhs; return *this;}

    gsKnotVectorIter& operator=(const gsKnotVectorIter &rhs) 
    {
        m_current  = rhs.m_current;
        m_begin    = rhs.m_begin;
        m_end      = rhs.m_end;
        return *this;
    }
    
    gsKnotVectorIter& operator+=(const int& rhs) 
    {
        // to do: better implementation
        for (int i=0; i!= rhs; ++i)
            ++(*this);
        return *this;
    }

    gsKnotVectorIter& operator-=(const int& rhs) 
    {	
        // to do: better implementation
        for (int i=0; i!= rhs; ++i)
            --(*this);
        return *this;
    }

    // Advance to the next unique knot value
    gsKnotVectorIter& operator++() 
    { 
        ++m_current;
        ignoreDublicates();
        return *this;
    }

    gsKnotVectorIter operator+(const difference_type& n) const
    {
        gsKnotVectorIter tmp(*this);
        tmp+= n;
        return tmp;
    }
    
    //gsKnotVectorIter operator+(const gsKnotVectorIter& rhs) const
    //friend inline gsKnotVectorIter operator+(const int& lhs, const gsKnotVectorIter& rhs)
    
    gsKnotVectorIter operator++(int) 
    {
        gsKnotVectorIter tmp(*this);
        ++(*this);
        return tmp;
    }
    
    // Retreat to the previous unique knot value
    gsKnotVectorIter & operator--()
    { 
        if ( m_current != m_end)
        {
            const T val = *m_current;
            for (m_current = m_current-1; m_current >= m_begin; --m_current)
                if ( val != *m_current )
                    break;
        }
        else
            --m_current;
        return *this;
    }
    
    gsKnotVectorIter operator--(int) 
    {
        gsKnotVectorIter tmp(*this);
        --(*this);
        return tmp;
    }
    
    gsKnotVectorIter operator-(const difference_type& n)  const
    {
        gsKnotVectorIter tmp(*this);
        tmp-=n;
        return tmp;
    }
    
    //gsKnotVectorIter operator-(const gsKnotVectorIter& rhs) const
    //friend inline gsKnotVectorIter operator-(const int& lhs, const gsKnotVectorIter& rhs)
    
// Comparison
public:
    friend bool operator==(const gsKnotVectorIter& x, 
                           const gsKnotVectorIter& y) 
    {
        return (x.m_current == y.m_current );
    }
    
   friend bool operator!=(const gsKnotVectorIter& x, 
                          const gsKnotVectorIter& y) 
    {
        return (x.m_current != y.m_current);
    }

   friend bool operator> (const gsKnotVectorIter& x, 
                          const gsKnotVectorIter& y) 
    {
        return (x.m_current>y.m_current);
    }

   friend bool operator< (const gsKnotVectorIter& x, 
                          const gsKnotVectorIter& y) 
    {
        return (x.m_current<y.m_current);
    }
    
    friend bool operator>=(const gsKnotVectorIter& x, 
                           const gsKnotVectorIter& y) 
    {
        return (x.m_current>=y.m_current);
    }

    friend bool operator<=(const gsKnotVectorIter& x, 
                           const gsKnotVectorIter& y) 
    {
        return (x.m_current<=y.m_current);
    }
    
private:
    
    // Ignore iterator positions with the same value
    inline void ignoreDublicates()
    {
        if ( m_current != m_end)
        {
            const T val = *m_current;
            for (m_current = m_current+1; m_current != m_end; ++m_current)
                if ( val != *m_current )
                    break;
            --m_current;
        }
    }

// Data
protected:
    
    knotptr m_current;
    knotptr m_begin;
    knotptr m_end;
};


}// namespace gismo

