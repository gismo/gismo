
#pragma once

#include <iterator>

namespace gismo
{

template<class T>
class gsHKnotVector;

// template <bool flag, class IsTrue, class IsFalse>
// struct choose;

// template <class IsTrue, class IsFalse>
// struct choose<true, IsTrue, IsFalse> {
//    typedef IsTrue type;
// };

// template <class IsTrue, class IsFalse>
// struct choose<false, IsTrue, IsFalse> {
//    typedef IsFalse type;
// };


  /** 
      Iterator over unique knots for gsHKnotVector.
  */
template <class T, bool isconst = false> 
class gsHKnotVectorUIter {
// Definitions
public:
    typedef std::random_access_iterator_tag iterator_category; 
    typedef T value_type;
    typedef std::ptrdiff_t difference_type;
    typedef typename choose<isconst, const T&, T&>::type reference;
    typedef typename choose<isconst, const T*, T*>::type pointer;
    typedef typename choose<isconst, typename gsHKnotVector<T>::const_global_uiterator,    			            typename gsHKnotVector<T>::global_uiterator >::type
      knotptr;
    typedef typename choose<isconst, const gsHKnotVector<T>, gsHKnotVector<T> >::type
      containertype;

// Constructors
public:
    gsHKnotVectorUIter() { }
    
    gsHKnotVectorUIter( containertype & c , bool start = true ) 
	: val ( (start ? c.global_ubegin()   : c.global_uend()) ), 
	  str ( c.stride() ) { }
    
    gsHKnotVectorUIter(const gsHKnotVectorUIter<T, false>& it) 
	: val(it.val), str(it.str) { }
    
// Accessors
public:
    reference operator*() const { return *val; }
    pointer  operator->() const { return  val; }
    reference operator[](const int& rhs) {return *(*this+rhs);}

// Arithmetic
public:

    gsHKnotVectorUIter& operator=(value_type* rhs) 
    {val = rhs; return *this;}

    gsHKnotVectorUIter& operator=(const gsHKnotVectorUIter &rhs) 
    { val=rhs.val; str=rhs.str; return *this;}

    gsHKnotVectorUIter& operator+=(const int& rhs) 
    {
	val += str*rhs;
	return *this;
    }

    gsHKnotVectorUIter& operator-=(const int& rhs) 
    {	
	val -= str*rhs;
	return *this;
    }

    gsHKnotVectorUIter& operator++() 
    { 
	val += str;
	return *this;
    }

    gsHKnotVectorUIter operator+(const difference_type& n) const
    {
	gsHKnotVectorUIter tmp(*this);
	tmp.val += str*n;
	return tmp;
    }

    //gsHKnotVectorUIter operator+(const gsHKnotVectorUIter& rhs) const
    //friend inline gsHKnotVectorUIter operator+(const int& lhs, const gsHKnotVectorUIter& rhs)

    gsHKnotVectorUIter operator++(int) 
    {
	gsHKnotVectorUIter tmp(*this);
	++(*this);
	return tmp;
    }
    
    gsHKnotVectorUIter & operator--()
    { 
	val -= str;
	return *this;
    }
    
    gsHKnotVectorUIter operator--(int) {
	gsHKnotVectorUIter tmp(*this);
	--(*this);
	return tmp;
    }

    gsHKnotVectorUIter operator-(const difference_type& n)  const
    {
	gsHKnotVectorUIter tmp(*this);
	tmp.val -= n*str;
	return tmp;
    }

    //gsHKnotVectorUIter operator-(const gsHKnotVectorUIter& rhs) const
    //friend inline gsHKnotVectorUIter operator-(const int& lhs, const gsHKnotVectorUIter& rhs)

// Comparison
public:
   friend bool operator==(const gsHKnotVectorUIter& x, 
			  const gsHKnotVectorUIter& y) 
    {
	return ( x.val == y.val );
    }
    
   friend bool operator!=(const gsHKnotVectorUIter& x, 
			  const gsHKnotVectorUIter& y) 
    {
	return ( x.val != y.val );
    }

   friend bool operator> (const gsHKnotVectorUIter& x, 
			  const gsHKnotVectorUIter& y) 
    {
	return (x.val>y.val);
    }

   friend bool operator< (const gsHKnotVectorUIter& x, 
			  const gsHKnotVectorUIter& y) 
    {
	return (x.val<y.val);
    }

   friend bool operator>=(const gsHKnotVectorUIter& x, 
			  const gsHKnotVectorUIter& y) 
    {
	return (x.val>=y.val);
    }

   friend bool operator<=(const gsHKnotVectorUIter& x, 
			  const gsHKnotVectorUIter& y) 
    {
	return (x.val<=y.val);
    }

// Data
protected:
    knotptr  val;
    unsigned str;
};


  /** 
      Iterator for gsHKnotVector.
  */

template <class T, bool isconst = false> 
class gsHKnotVectorIter {
// Definitions
public:
    typedef std::random_access_iterator_tag iterator_category; 
    typedef T value_type;
    typedef std::ptrdiff_t difference_type;
    typedef typename choose<isconst, const T&, T&>::type reference;
    typedef typename choose<isconst, const T*, T*>::type pointer;
    typedef typename choose<isconst, typename gsHKnotVector<T>::const_global_uiterator,    			            typename gsHKnotVector<T>::global_uiterator >::type
      knotptr;
    typedef typename choose<isconst, typename gsHKnotVector<T>::const_global_miterator,    			            typename gsHKnotVector<T>::global_miterator >::type
      multptr;
    typedef typename choose<isconst, const gsHKnotVector<T>,    
			    gsHKnotVector<T> >::type
      containertype;

// Constructors
public:
    gsHKnotVectorIter() : i(0) { }
    
    gsHKnotVectorIter( containertype & c , bool start = true ) 
	: val ( (start ? c.global_ubegin() : c.global_uend()) ), 
	  mult( (start ? c.global_mbegin() : c.global_mend()) ),
	  str (  c.stride()                                   ), 
	  i   ( ( start ? 0                : c.globalSize())  ) { }
    
    gsHKnotVectorIter(const gsHKnotVectorIter<T, false>& it) 
		: val(it.val), mult(it.mult), str(it.str), i(it.i) { }
    
// Accessors
public:
    reference operator*() const { return *val; }
    pointer  operator->() const { return  val; }
    reference operator[](const int& rhs) {return *(*this+rhs);}

// Arithmetic
public:

    gsHKnotVectorIter& operator=(value_type* rhs) 
    {val = rhs; return *this;}

    gsHKnotVectorIter& operator=(const gsHKnotVectorIter &rhs) 
    { val=rhs.val; mult= rhs.mult; str=rhs.str; return *this;}

    gsHKnotVectorIter& operator+=(const int& rhs) 
    {
	val += str*rhs;
	while ( i >= *mult )
	    {
		++mult;
		++val ;
	    }
	return *this;
    }

    gsHKnotVectorIter& operator-=(const int& rhs) 
    {	
	val-= str*rhs;
	while ( i < *(mult-1) )
	    {
		--mult;
		--val ;
	    }
	return *this;
    }

    gsHKnotVectorIter& operator++() // TO DO: test
    { 
	if ( ++i >= *mult )
	    {
		val  += str ;
		const int tmp = *mult;
		mult += str-1 ;
		i    += *mult - tmp;
		++mult;
	    }
	return *this;
    }

    gsHKnotVectorIter operator+(const difference_type& n) const
    {
	gsHKnotVectorIter tmp(*this);
	tmp.i+= n;
	while ( tmp.i >= * tmp.mult )
	    {
		++ tmp.mult;
		++ tmp.val ;
	    }
	return tmp;
    }

    //gsHKnotVectorIter operator+(const gsHKnotVectorIter& rhs) const
    //friend inline gsHKnotVectorIter operator+(const int& lhs, const gsHKnotVectorIter& rhs)

    gsHKnotVectorIter operator++(int) 
    {
	gsHKnotVectorIter tmp(*this);
	++(*this);
	return tmp;
    }
    
    gsHKnotVectorIter & operator--()
    { 
	if ( --i < *(mult-1) )
	    {
		--mult;
		--val ;
	    }
	return *this;
    }
    
    gsHKnotVectorIter operator--(int) {
	gsHKnotVectorIter tmp(*this);
	--(*this);
	return tmp;
    }

    gsHKnotVectorIter operator-(const difference_type& n)  const
    {
	gsHKnotVectorIter tmp(*this);
	tmp.i -= n;
	while ( tmp.i < *(tmp.mult-1) )
	    {
		--tmp.mult;
		--tmp.val ;
	    }
	return tmp;
    }

    //gsHKnotVectorIter operator-(const gsHKnotVectorIter& rhs) const
    //friend inline gsHKnotVectorIter operator-(const int& lhs, const gsHKnotVectorIter& rhs)

// Comparison
public:
   friend bool operator==(const gsHKnotVectorIter& x, 
			  const gsHKnotVectorIter& y) 
    {
	return (x.val == y.val && x.mult==y.mult && x.i==y.i);
    }
    
   friend bool operator!=(const gsHKnotVectorIter& x, 
			  const gsHKnotVectorIter& y) 
    {
	return (x.val != y.val || x.mult!=y.mult || x.i!=y.i);
    }

   friend bool operator> (const gsHKnotVectorIter& x, 
			  const gsHKnotVectorIter& y) 
    {
	return (x.i>y.i);
    }

   friend bool operator< (const gsHKnotVectorIter& x, 
			  const gsHKnotVectorIter& y) 
    {
	return (x.i<y.i);
    }

   friend bool operator>=(const gsHKnotVectorIter& x, 
			  const gsHKnotVectorIter& y) 
    {
	return (x.i>=y.i);
    }

   friend bool operator<=(const gsHKnotVectorIter& x, 
			  const gsHKnotVectorIter& y) 
    {
	return (x.i<=y.i);
    }

// Data
protected:
    knotptr val ;
    multptr mult;
    unsigned i  ;
    unsigned str;
};



}// namespace gismo
