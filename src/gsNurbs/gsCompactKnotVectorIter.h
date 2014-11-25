
#pragma once

#include <iterator>

#include <gsCore/gsTemplateTools.h>


namespace gismo
{

template<class T>
class gsCompactKnotVector;

  /** 
      Iterator for gsCompactKnotVector.
      
      Beware of non-const version: changing one knot value does not
      insert a new simple knot, but just changes the value of all
      instances of the multiple knot.
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

// Constructors
public:
    gsCompactKnotVectorIter() : index(0), span(0) { }
    
    gsCompactKnotVectorIter( containertype & c , bool start = true ) 
	: val   ( (start ? c.ubegin()  : c.uend() ) ), 
	  mult  ( (start ? c.mbegin()  : c.mend() ) ),
	  index ( (start ? 0           : c.size() ) ),
      span  ( (start ? 0           : c.uSize()-1) ),
      multEnd ( c.mend() ){ }
    
    gsCompactKnotVectorIter(const gsCompactKnotVectorIter<T, false>& it) 
        : val(it.val), mult(it.mult), index(it.index), span(it.span), multEnd(it.multEnd) { }
    
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
    { val=rhs.val; mult= rhs.mult; index=rhs.index; span=rhs.span; multEnd=rhs.multEnd; return *this;}

    gsCompactKnotVectorIter& operator+=(const int& rhs) 
    {
        index+= rhs;
        while ( mult != multEnd &&
                index >= * mult )
        {
            ++ mult;
            ++ val ;
            ++span ;
        }
        return *this;
    }

    gsCompactKnotVectorIter& operator-=(const int& rhs) 
    {	
        index-= rhs;
        while ( span > 0 &&
                index < *(mult-1) )
        {
            --mult;
            --val ;
            --span;
        }
        return *this;
    }

    gsCompactKnotVectorIter& operator++() 
    { 
        ++index;
        if ( mult != multEnd && index >= *mult )
	    {
            ++mult;
            ++val ;
            ++span;
	    }
        return *this;
    }

    gsCompactKnotVectorIter operator+(const difference_type& n) const
    {
        gsCompactKnotVectorIter tmp(*this);
        tmp.index+= n;
        while ( tmp.mult != tmp.multEnd && 
                tmp.index >= * tmp.mult )
        {  
            ++ tmp.mult;
            ++ tmp.val ;
            ++ tmp.span ;
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
        --index;
        if ( span > 0 && index < *(mult-1) )
	    {
            --mult;
            --val ;
            --span;
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
        tmp.index -= n;
        while ( tmp.span > 0 &&
                tmp.index < *(tmp.mult-1) )
        {
            --tmp.mult;
            --tmp.val ;
            --tmp.span;
        }
        return tmp;
    }
    
    gsCompactKnotVectorIter operator-(const gsCompactKnotVectorIter& rhs) const
    {
        gsWarn<< "gsCompactKnotVectorIter difference not implemented.\n";
        gsCompactKnotVectorIter tmp(*this);
        tmp -= index - rhs.index;
        return tmp;
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
            x.index==y.index);
    }
    
   friend bool operator!=(const gsCompactKnotVectorIter& x, 
                          const gsCompactKnotVectorIter& y) 
   {
       return (x.index != y.index);
   }

   friend bool operator> (const gsCompactKnotVectorIter& x, 
                          const gsCompactKnotVectorIter& y) 
   {
       return (x.index>y.index);
   }

   friend bool operator< (const gsCompactKnotVectorIter& x, 
                          const gsCompactKnotVectorIter& y) 
    {
        return (x.index<y.index);
    }

   friend bool operator>=(const gsCompactKnotVectorIter& x, 
                          const gsCompactKnotVectorIter& y) 
   {
       return (x.index>=y.index);
   }

   friend bool operator>=(const gsCompactKnotVectorIter& x, 
                          const difference_type& n) 
   {
       return (x.index>=n);
   }

   friend bool operator<=(const gsCompactKnotVectorIter& x, 
                          const gsCompactKnotVectorIter& y) 
   {
       return (x.index<=y.index);
   }

// Data
public:
    knotptr val;
    multptr mult;
    unsigned index;
    int span;
    multptr multEnd;
};



}// namespace gismo

