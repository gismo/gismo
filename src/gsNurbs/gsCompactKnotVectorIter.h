/** @file gsCompactKnotVectorIter.h

    @brief Iterator over the knots of a compact knotvector

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
class gsCompactKnotVector;

/** 
    @brief Iterator for gsCompactKnotVector.
      
    Beware of non-const version: changing one knot value does not
    insert a new simple knot, but just changes the value of all
    instances of the multiple knot.
      
    \ingroup Nurbs
*/
template <class T, bool isconst = false> 
class gsCompactKnotVectorIter {
// Definitions
public:
    typedef std::random_access_iterator_tag iterator_category; 
    typedef T value_type;
    typedef std::ptrdiff_t difference_type;
    typedef typename choose<isconst, const T&, T&>::type reference;
    typedef typename choose<isconst, const T*, T*>::type pointer;
    typedef typename choose<isconst, typename gsCompactKnotVector<T>::const_uiterator,    
                            typename gsCompactKnotVector<T>::uiterator >::type
    knotptr;
    typedef typename choose<isconst, typename std::vector<unsigned>::const_iterator,    
                            typename std::vector<unsigned>::iterator >::type
    multptr;
    typedef typename choose<isconst, const gsCompactKnotVector<T>,    
                            gsCompactKnotVector<T> >::type
    containertype;

    friend class  gsCompactKnotVectorIter<T, !isconst> ;

// Constructors
public:
    gsCompactKnotVectorIter() : m_index(0), m_span(0) { }
    
    gsCompactKnotVectorIter( containertype & c , bool start = true ) 
	: val   ( (start ? c.ubegin()  : c.uend() ) ), 
	  mult  ( (start ? c.mbegin()  : c.mend() ) ),
      m_index ( (start ? 0           : c.size() ) ),
      m_span  ( (start ? 0           : c.uSize()-1) ),
      multEnd ( c.mend() ){ }

    gsCompactKnotVectorIter( containertype & c , size_t pos)
    : val   ( c.ubegin()+pos ),
      mult  ( c.mbegin()+pos ),
      m_index ( *mult-1 ),
      m_span  ( pos ),
      multEnd ( c.mend() ){ }

    gsCompactKnotVectorIter(const gsCompactKnotVectorIter<T, false>& it) 
    : val(it.val), mult(it.mult), m_index(it.index()), m_span(it.span()), multEnd(it.multEnd) { }
    
// Accessors
public:
    reference operator*() const { return *val; }
    pointer  operator->() const { return  val; }
    reference operator[](const int& rhs) {return *(*this+rhs);}

// Arithmetic
public:

    gsCompactKnotVectorIter& operator=(value_type* rhs) 
    {val = rhs; return *this;}

    gsCompactKnotVectorIter& operator=(const gsCompactKnotVectorIter &rhs) 
    { val=rhs.val; mult= rhs.mult; m_index=rhs.m_index; m_span=rhs.m_span; multEnd=rhs.multEnd; return *this;}

    gsCompactKnotVectorIter& operator+=(const int& rhs) 
    {
        m_index+= rhs;
        while ( mult != multEnd &&
                m_index >= * mult )
        {
            ++ mult;
            ++ val ;
            ++ m_span ;
        }
        return *this;
    }

    gsCompactKnotVectorIter& operator-=(const int& rhs) 
    {	
        m_index-= rhs;
        while ( m_span > 0 &&
                m_index < *(mult-1) )
        {
            --mult;
            --val ;
            --m_span;
        }
        return *this;
    }

    gsCompactKnotVectorIter& operator++() 
    { 
        ++m_index;
        if ( mult != multEnd && m_index >= *mult )
	    {
            ++mult;
            ++val ;
            ++m_span;
	    }
        return *this;
    }

    gsCompactKnotVectorIter operator+(const difference_type& n) const
    {
        gsCompactKnotVectorIter tmp(*this);
        tmp.m_index+= n;
        while ( tmp.mult != tmp.multEnd && 
                tmp.m_index >= * tmp.mult )
        {  
            ++ tmp.mult;
            ++ tmp.val ;
            ++ tmp.m_span ;
        }
        return tmp;
    }

    //gsCompactKnotVectorIter operator+(const gsCompactKnotVectorIter& rhs) const
    //friend inline gsCompactKnotVectorIter operator+(const int& lhs, const gsCompactKnotVectorIter& rhs)

    gsCompactKnotVectorIter operator++(int) 
    {
        gsCompactKnotVectorIter tmp(*this);
        ++(*this);
        return tmp;
    }
    
    gsCompactKnotVectorIter & operator--()
    { 
        --m_index;
        if ( m_span > 0 && m_index < *(mult-1) )
	    {
            --mult;
            --val ;
            --m_span;
	    }
        return *this;
    }
    
    gsCompactKnotVectorIter operator--(int) 
    {
        gsCompactKnotVectorIter tmp(*this);
        --(*this);
        return tmp;
    }

    gsCompactKnotVectorIter operator-(const difference_type& n)  const
    {
        gsCompactKnotVectorIter tmp(*this);
        tmp.m_index -= n;
        while ( tmp.m_span > 0 &&
                tmp.m_index < *(tmp.mult-1) )
        {
            --tmp.mult;
            --tmp.val ;
            --tmp.m_span;
        }
        return tmp;
    }
    
    difference_type operator-(const gsCompactKnotVectorIter& rhs) const
    {
        return m_index - rhs.m_index;
    }

    //friend inline gsCompactKnotVectorIter operator-(const int& lhs, const gsCompactKnotVectorIter& rhs)

// Comparison
public:
    friend bool operator==(const gsCompactKnotVectorIter& x, 
                           const gsCompactKnotVectorIter& y) 
    {
        return ( 
            //x.val == y.val &&
            //x.mult==y.mult && 
            x.m_index==y.m_index);
    }
    
    friend bool operator!=(const gsCompactKnotVectorIter& x, 
                           const gsCompactKnotVectorIter& y) 
    {
        return (x.m_index != y.m_index);
    }

    friend bool operator> (const gsCompactKnotVectorIter& x, 
                           const gsCompactKnotVectorIter& y) 
    {
        return (x.m_index>y.m_index);
    }

    friend bool operator< (const gsCompactKnotVectorIter& x, 
                           const gsCompactKnotVectorIter& y) 
    {
        return (x.m_index<y.m_index);
    }

    friend bool operator>=(const gsCompactKnotVectorIter& x, 
                           const gsCompactKnotVectorIter& y) 
    {
        return (x.m_index>=y.m_index);
    }

    friend bool operator>=(const gsCompactKnotVectorIter& x, 
                           const difference_type& n) 
    {
        return (x.m_index>=n);
    }

    friend bool operator<=(const gsCompactKnotVectorIter& x, 
                           const gsCompactKnotVectorIter& y) 
    {
        return (x.m_index<=y.m_index);
    }

public:
    unsigned index() const
    { return m_index; }

    int span() const
    { return m_span; }

// Data
private:
    knotptr val;
    multptr mult;
    unsigned m_index;
    int m_span;
    multptr multEnd;
};



}// namespace gismo

