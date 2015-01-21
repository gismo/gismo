/** @file gsKnotVector.hpp

    @brief Provides implementation of the KnotVector class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsNurbs/gsCompactKnotVector.h>
#include <gsNurbs/gsCompactKnotVectorIter.h>

#include <gsIO/gsXmlUtils.h>

#include <numeric>

namespace gismo
{

/// Private data member for knot vector
template <class T>
class gsKnotVectorPrivate
{
public:
  /// The knot values
  std::vector<T> knots;
  /// The degree (offset) of the knot vector
  int p;
};

template <class T>
gsKnotVector<T>::gsKnotVector() : gsDomain<T>(), my(new gsKnotVectorPrivate<T>) 
{ 
    my->p = 0;
}
    
template <class T>
gsKnotVector<T>::gsKnotVector( gsKnotVector const & other)
{
    my = new gsKnotVectorPrivate<T>(*other.my);
    my->knots= other.my->knots;
    my->p= other.my->p;
}

template <class T>
gsKnotVector<T>::gsKnotVector(int p) 
    : gsDomain<T>(), my(new gsKnotVectorPrivate<T>)
{ 
    my->p = p;
}

template <class T>
gsKnotVector<T>::gsKnotVector(int p, unsigned sz ) 
    : gsDomain<T>(), my(new gsKnotVectorPrivate<T>)
{ 
    my->p = p;
    my->knots.resize(sz); 
}

template <class T>
gsKnotVector<T>::gsKnotVector(T u0, T u1, unsigned interior, 
                              unsigned mult_ends, unsigned mult_interior, int degree)
    : gsDomain<T>(), my(new gsKnotVectorPrivate<T>)
{
    initUniform( u0, u1, interior, mult_ends, mult_interior, degree );
}


template <class T>
gsKnotVector<T>::gsKnotVector(std::vector<T> const& knots, int degree, int regularity)
    : gsDomain<T>(), my(new gsKnotVectorPrivate<T>)
{
  
  typename std::vector<T>::const_iterator itr;
  int mult= degree - regularity ;
  
  for ( int j=0; j<degree-1; j++ )
    my->knots.push_back( knots.front() );
  
  for ( itr= knots.begin(); itr != knots.end(); ++itr )
    for ( int j=0; j< mult; j++ )
      my->knots.push_back( *itr );
  
  for ( int j=0; j<degree-1; j++ )
    my->knots.push_back( knots.back() );
  
  my->p=degree;
}


template <class T>
gsKnotVector<T>::gsKnotVector(int degree, std::vector<T> const& knots)
    : gsDomain<T>(), my(new gsKnotVectorPrivate<T>)
{
  my->p = degree;
  my->knots = knots;
}

template <class T>
gsKnotVector<T>::gsKnotVector(gsCompactKnotVector<T> ckv)
    : gsDomain<T>(), my(new gsKnotVectorPrivate<T>)
{
  my->p = ckv.degree();
  for(gsCompactKnotVectorIter<T, true > it = ckv.begin();it!=ckv.end();it++)
  {
      my->knots.push_back(*it);
  }
}

template <class T>
gsKnotVector<T>::gsKnotVector(int deg, const_iterator start, const_iterator end)
    : gsDomain<T>(), my(new gsKnotVectorPrivate<T>)
{
  my->p = deg;
  my->knots.assign(start, end) ;
}


template <class T>
gsKnotVector<T>::~gsKnotVector() 
{ 
  delete my; 
  my = NULL;
}

template <class T> void 
gsKnotVector<T>::initUniform( T u0, T u1, unsigned interior, unsigned mult_ends, 
                             unsigned mult_interior, int degree)
{    
    my->knots.clear();
    my->knots.reserve( 2 * mult_ends + interior*mult_interior);
 
    const T h = (u1-u0) / (interior+1);

    my->knots.insert(my->knots.begin(), mult_ends, u0);
    
    for ( unsigned i=1; i<=interior; i++ )
        my->knots.insert(my->knots.end(), mult_interior, 
                         u0 + i*h );

    my->knots.insert(my->knots.end(), mult_ends, u1);

    if (degree == -1)   
        my->p = mult_ends-1;
    else
        my->p = degree;
}

template <class T> void 
gsKnotVector<T>::initUniform(unsigned numKnots, unsigned mult_ends,
                             unsigned mult_interior, int degree)
{
    initUniform(0.0, 1.0, numKnots - 2, mult_ends, mult_interior, degree );
}

template <class T> void 
gsKnotVector<T>::initGraded(T u0, T u1, unsigned interior, int degree, 
                           T grading, unsigned mult_interior )
{
// \todo: add grading target
    my->p = degree;
    my->knots.clear();
    my->knots.reserve( 2 *  (degree+1) + interior*mult_interior );

    const T h = (u1-u0) / (interior+1);

    my->knots.insert(my->knots.begin(), degree+1, u0);
    
    for ( unsigned i=1; i<=interior; i++ )
        my->knots.insert(my->knots.end(), mult_interior, 
                         math::pow(i*h, 1.0/grading) );

    my->knots.insert(my->knots.end(), degree+1, u1);
}


template <class T> void 
gsKnotVector<T>::initGraded(unsigned numKnots, int degree, 
                           T grading, unsigned mult_interior )
{
    initGraded( 0.0, 1.0, numKnots - 2, degree, grading, mult_interior);
}

template <class T> void 
gsKnotVector<T>::initClamped(T u0, T u1, int degree, 
                             unsigned interior,
                             unsigned mult_interior)
{
    initUniform(u0, u1, interior, degree + 1, mult_interior, degree );
}

template <class T> void 
gsKnotVector<T>::initClamped(int degree, unsigned numKnots, 
                             unsigned mult_interior)

{
    GISMO_ASSERT( numKnots > 1 , "Not enough knots.");
    initUniform(0.0, 1.0, numKnots - 2, degree + 1, mult_interior, degree );
}

template <class T>
typename gsKnotVector<T>::iterator
gsKnotVector<T>::begin()
{
    return my->knots.begin();
}

template <class T>
typename gsKnotVector<T>::reverse_iterator
gsKnotVector<T>::rbegin()
{
    return my->knots.rbegin();
}

template <class T>
typename gsKnotVector<T>::const_iterator
gsKnotVector<T>::begin() const 
{
    return my->knots.begin();
}

template <class T>
typename gsKnotVector<T>::uiterator
gsKnotVector<T>::ubegin()
{
    return gsKnotVectorIter<T,false>(*this,true);
}

template <class T>
typename gsKnotVector<T>::reverse_uiterator
gsKnotVector<T>::urbegin()
{
    return reverse_uiterator(uend());
}

template <class T>
typename gsKnotVector<T>::const_uiterator
gsKnotVector<T>::ubegin() const 
{
    return gsKnotVectorIter<T,true>(*this,true);
}

template <class T>
typename gsKnotVector<T>::const_reverse_uiterator
gsKnotVector<T>::urbegin() const 
{
    return const_reverse_uiterator(uend());
}

template <class T>
typename gsKnotVector<T>::const_reverse_iterator
gsKnotVector<T>::rbegin() const 
{
    return my->knots.rbegin();
}

template <class T>
typename gsKnotVector<T>::iterator
gsKnotVector<T>::end()
{
    return my->knots.end();
}

template <class T>
typename gsKnotVector<T>::reverse_iterator
gsKnotVector<T>::rend()
{
    return my->knots.rend();
}

template <class T>
typename gsKnotVector<T>::const_iterator
gsKnotVector<T>::end() const 
{
    return my->knots.end();
}

template <class T>
typename gsKnotVector<T>::const_reverse_iterator
gsKnotVector<T>::rend() const 
{
    return my->knots.rend();
}

template <class T>
typename gsKnotVector<T>::uiterator
gsKnotVector<T>::uend()
{
    return gsKnotVectorIter<T,false>(*this,false);
}

template <class T>
typename gsKnotVector<T>::reverse_uiterator
gsKnotVector<T>::urend()
{
    return reverse_uiterator(ubegin());
}

template <class T>
typename gsKnotVector<T>::const_uiterator
gsKnotVector<T>::uend() const 
{
    return gsKnotVectorIter<T,true>(*this,false);
}

template <class T>
typename gsKnotVector<T>::const_reverse_uiterator
gsKnotVector<T>::urend() const 
{
    return const_reverse_uiterator(ubegin());
}

template <class T>
void gsKnotVector<T>::swap(gsKnotVector &other)
{
  std::swap(my->p, other.my->p);
  my->knots.swap(other.my->knots);
}

template <class T>
void gsKnotVector<T>::resize(unsigned sz)
{
    my->knots.resize(sz);
}

template <class T>
void gsKnotVector<T>::clear()
{
    my->knots.clear();
}

template <class T>
T  gsKnotVector<T>::operator [] (size_t i) const 
{ return my->knots[i]; }
  
template <class T>
T& gsKnotVector<T>::operator [] (size_t i) { return my->knots[i]; }

template <class T>
bool gsKnotVector<T>::operator==(const gsKnotVector<T> &other) const
{
  // TODO: use tolerance?
  return (this->degree() == other.degree())
      && (this->size() == other.size())
      && std::equal( this->begin(), this->end(), other.begin() );
}

template <class T>
bool gsKnotVector<T>::operator!=(const gsKnotVector<T> &other) const
{
  return !(*this==other);
}

template <class T>
  T  gsKnotVector<T>::at (size_t i) const { return my->knots.at(i); }
template <class T>
T& gsKnotVector<T>::at (size_t i) { return my->knots.at(i); }
  
template <class T>
T gsKnotVector<T>::first () const { return my->knots.front(); }
template <class T>
T gsKnotVector<T>::last  () const { return my->knots.back(); }
  
template <class T>
void gsKnotVector<T>::push_back( T knot)  
{ 
    GISMO_ASSERT(my->knots.empty() || knot>=my->knots.back(), "Knot out of range");
    my->knots.push_back(knot); 
}

template <class T>
void gsKnotVector<T>::push_back( T knot, int mult)  
{ 
    GISMO_ASSERT(my->knots.empty() || knot>=my->knots.back(), "Knot out of range");
    my->knots.insert(my->knots.end(), mult, knot); 
}
  
template <class T>
void gsKnotVector<T>::push_front( T knot)  
{   
    GISMO_ASSERT( my->knots.empty() || knot<=my->knots.front(), "Knot out of range");
    my->knots.insert(my->knots.begin(), knot); 
}

template <class T>
void gsKnotVector<T>::push_front( T knot, int mult)  
{   
    GISMO_ASSERT( my->knots.empty() || knot<=my->knots.front(), "Knot out of range");
    my->knots.insert(my->knots.begin(), mult, knot); 
}
  
template <class T>
int gsKnotVector<T>::size() const 
{ return my->knots.size(); }
    
template <class T>
const std::vector<T>& gsKnotVector<T>::get() const 
{ return my->knots; }

template <class T>
gsVector<T>* gsKnotVector<T>::getVector() const
{
    std::vector<T> u = unique();
    gsVector<T> * res = new gsVector<T>(u.size());

    for ( unsigned i = 0; i < u.size(); ++i)
        (*res)[i] = u[i];

    return res; 
}

template <class T> // AM: we can do this more efficiently..
int gsKnotVector<T>::numKnotSpans() const
{ return this->unique().size()-1; }

// template <class T> // TO DO
// int gsKnotVector<T>::numDomainSpans() const
// { return number of elements; }

template <class T> // AM: we can do this more efficiently..
int gsKnotVector<T>::findElementIndex(T u) const
{
    const std::vector<T> uKnots = this->unique();

    unsigned low  = 0;// span index of domainStart (!)
    unsigned high = static_cast<unsigned>(uKnots.size() - 1); // span index of domainEnd (!)

    GISMO_ASSERT( (u >= uKnots[low]) && ( u  <= uKnots[high] ), 
                  "The requested abscissae u="<<u<<" is not in the knot vector." );

    if (u == uKnots[high]) // remove ?
        return high-1;

    do
    {
        const unsigned mid = (low + high) >> 1;
        if ( u < uKnots[mid] )
            high = mid;
        else if (u >= uKnots[mid+1])
            low  = mid;
        else
            return mid;
    }
    while (true);
}

template <class T>
std::vector<T> gsKnotVector<T>::knotSpanLengths() const
{
    const std::vector<T> u = this->unique();

    std::vector<T> spans;
    spans.reserve(u.size()-1);


    for ( typename  std::vector<T>::const_iterator
          it= u.begin()+1; it != u.end(); ++it )
        spans.push_back( *it - *(it-1) );

    return spans;
}


template <class T>
T gsKnotVector<T>::maxKnotSpanLength() const
{
    T hmax = 0.0;
    for (int i = 0; i < size() - 1; ++i)
        hmax = std::max(hmax, at(i+1) - at(i));
    return hmax;
}

template <class T>
T gsKnotVector<T>::minKnotSpanLength() const
{
    T tmp, hmin = std::numeric_limits<T>::max();
    for (int i = 0; i < size() - 1; ++i)
    {
        tmp = at(i+1) - at(i);
        if ( (tmp > 0) && (tmp < hmin) ) 
            hmin = tmp;
    }
    return hmin;
}


template <class T>
void gsKnotVector<T>::transform(T c, T d)
{
    T a = my->knots.front();
    T b = my->knots.back();
    my->knots.front() = c;
    my->knots.back()  = d;       

    for (typename std::vector<T>::iterator it = my->knots.begin()+1; it != my->knots.end()-1; ++it)
        (*it) = c + ((*it) - a) * (d - c) / (b - a) ;
}


template <class T>
void gsKnotVector<T>::merge(gsKnotVector<T> other)
{ 
        // check for degree
    if ( this->degree() != other.degree() )
    {std::cout<<"gsKnotVector: Cannot merge KnotVectors of different degree"<<"\n"; return;}

    other.setFirst( last() ) ;
    // TO DO: check that this->degree() = K.degree()
    // TO DO: check that last p+1 this->my->knots = first p+1 knots of K
    this->pop_back();
    this->append( other.begin()+my->p+1, other.end() );
}


template <class T>
std::vector<T> gsKnotVector<T>::unique() const 
{
  std::vector<T> result = my->knots;
  result.erase( std::unique( result.begin(), result.end() ),
                result.end() );
  return result;
}

template <class T>
std::vector<T> gsKnotVector<T>::unique(size_t const i, size_t const & j) const 
{
  typename std::vector<T>::const_iterator first = my->knots.begin() + i;
  typename std::vector<T>::const_iterator last  = my->knots.begin() + j;
  return std::vector<T>(first,last);
}


template <class T>
std::vector<T> gsKnotVector<T>::breaks() const
{
    std::vector<T> result(my->knots.begin() + my->p, my->knots.end() - my->p);
    result.erase( std::unique( result.begin(), result.end() ),
                  result.end() );
    return result; 
} 


template <class T>
typename gsKnotVector<T>::const_iterator 
gsKnotVector<T>::findspanIter (T u) const
{
    GISMO_ASSERT( ( u >= my->knots[my->p]) && ( u  <= *(my->knots.end()-my->p) ), 
                  "The requested abscissae u="<<u<<" is not in the knot vector." );
    /// \todo: reduce calls to findspan

    if ( u == *(my->knots.end()-my->p-1) )
        return my->knots.end()-my->p-2;
    else
        return std::upper_bound(my->knots.begin()+my->p, my->knots.end()-my->p, u ) -1 ; 
}

template <class T>
inline unsigned gsKnotVector<T>::findspan (T u) const
{
    // equivalent: return findspanIter(u) - my->knots.begin();

    // NB: version according to NURBS Book, Algorithm A2.1 (p. 68)

    /// \todo: reduce calls to findspan
    //gsInfo<<"u is " << u<<", "<<((u >= my->knots[my->p]) && ( u  <= my->knots[m-my->p])) <<"\n";

    const unsigned m = my->knots.size() - 1;        // last knot index
    unsigned low = my->p, high = m - my->p;

    GISMO_ASSERT( (u >= my->knots[low]) && ( u  <= my->knots[high]), 
		  "The requested abscissae u="<<u<<" is not in the knot vector." );

    if (u == my->knots[high])
        return high - 1;

    do
    {
        const unsigned mid = (low + high) >> 1;
        if ( u < my->knots[mid] )
            high = mid;
        else if (u >= my->knots[mid+1])
            low  = mid;
        else
            return mid;
    }
    while (true);
}


template <class T>
inline gsMatrix<unsigned,1> * gsKnotVector<T>::findspan (const gsMatrix<T,1> & u) const
{
    gsMatrix<unsigned,1> * fs = new gsMatrix<unsigned,1>(1, u.cols() );

    for( index_t i = 0; i < u.cols(); i++ )
        (*fs)(0,i) = findspan( u(0,i) );

    return fs;
}


// TODO
//template <class T>
//void gsKnotVector<T>::scale (T u0, T u1)
//{
//
//}

template <typename T>
void gsKnotVector<T>::reverse()
{
    // knot vector starts with first() and finish with last()
    const T ab = this->first() + this->last();

    std::reverse(my->knots.begin(), my->knots.end());

    for (iterator it = my->knots.begin(); it != my->knots.end(); ++it)
    {
        *it = ab - *it;
    }
}

template <class T>
void gsKnotVector<T>::insert(T knot, int mult)
{
    typename std::vector<T>::iterator itr
        = std::lower_bound(my->knots.begin(), my->knots.end(), knot);

    my->knots.insert(itr, mult, knot);
}

template <class T>
void gsKnotVector<T>::insert(std::vector<T> const & knots, int mult)
{
    for(typename std::vector<T>::const_iterator it=knots.begin();it!=knots.end();++it)
        insert(*it,mult);
}

template <class T>
bool gsKnotVector<T>::has(T knot) const
{return std::binary_search(my->knots.begin(), my->knots.end(), knot); };

template <class T>
void gsKnotVector<T>::append(const_iterator const & v0, const_iterator const & v1 ) 
{ my->knots.insert( my->knots.end(), v0, v1 ); }

template <class T>
void gsKnotVector<T>::pop_back(int i) 
{ my->knots.erase( my->knots.end()-i,my->knots.end() ); }
  
template <class T>
void gsKnotVector<T>::pop_front(int i) 
{ my->knots.erase( my->knots.begin(), my->knots.begin()+i ); }
  
template <class T>
void gsKnotVector<T>::addConstant(T t)  
{ std::transform(my->knots.begin(), my->knots.end(), my->knots.begin(),
                 std::bind2nd(std::plus<T>(), t));
}

template <class T>
void gsKnotVector<T>::setFirst(T t) 
{ addConstant( t - first() ); }

template <class T>
bool gsKnotVector<T>::isUniform() const
{
    //spans: std::adjacent_difference
    const std::vector<T> u = this->unique();
    T df = u[1]-u[0];
    for ( typename  std::vector<T>::const_iterator 
          it= u.begin()+2; it != u.end(); ++it )
        if ( *it - *(it-1) != df )
            return false;
    return true;
}

template <class T>
bool gsKnotVector<T>::isOpen() const
{
    return ( multiplicity(my->knots.front()) == my->p+1 &&
             multiplicity(my->knots.back() ) == my->p+1 );
}

template <class T>
T gsKnotVector<T>::firstInterval() const
{
    typename std::vector<T>::const_iterator it = my->knots.begin();
    const T u = *(it++);
    while ( u == *it ) 
        ++it;
    return *it-u;
}

template <class T>
void gsKnotVector<T>::uniformRefine(gsMatrix<T> const & interval, int numKnots)
{
  std::vector<T> u = this->unique( this->findspan(interval(0,0)), 
				   findspan( interval(0,1)) + 1 );

  for ( typename  std::vector<T>::iterator it= u.begin(); it != u.end(); ++it )
    this->insert ( *it, numKnots) ;
}

template <class T>
void gsKnotVector<T>::uniformRefine(int numKnots, int mul)
{
    T k0 = my->knots[my->p],
      k1 = *(my->knots.end()-my->p-1);

    std::vector<T> u = this->unique();

    for (std::size_t i = 0; i < u.size() - 1; ++i)
        for (int k = 1; k <= numKnots; ++k)
            this->insert(((numKnots+1-k) * u[i] + k * u[i+1]) / (numKnots + 1),mul);

    // trim extra knots
    typename  std::vector<T>::iterator it =
        std::upper_bound( my->knots.begin(), my->knots.end(), k0 );
    my->knots.erase( my->knots.begin(), it-my->p-1 );

    it = std::lower_bound( my->knots.begin(), my->knots.end(), k1 );
    my->knots.erase(it+my->p+1, my->knots.end() );
}

template <class T>
void gsKnotVector<T>::refineSpans(const std::vector<unsigned> & spanIndices, int numKnots)
{
    const std::vector<T> uK = this->unique();
    
    for ( typename  std::vector<unsigned>::const_iterator it= spanIndices.begin(); 
          it != spanIndices.end(); ++it )
        for (int k = 1; k <= numKnots; ++k)
            this->insert(((numKnots+1-k) * uK[*it] + k * uK[*it+1]) / (numKnots + 1));
}


template <class T>
void gsKnotVector<T>::getUniformRefinementKnots(int knotsPerSpan, std::vector<T>& result, int mul) const
{
    std::vector<T> u = this->unique();
    result.clear();
    result.reserve((u.size() - 1) * knotsPerSpan*mul);

    for (std::size_t i = 0; i < u.size() - 1; ++i)
        for (int k = 1; k <= knotsPerSpan; ++k)
            result.insert(result.end(),mul,((knotsPerSpan+1-k) * u[i] + k * u[i+1]) / (knotsPerSpan + 1));
}


template <class T>
void gsKnotVector<T>::degreeElevate(int const & i)
{
    GISMO_ASSERT(i>=0, "Degree elevation is only possible for non-negative numbers.");

    if (i==0) 
        return;

    T k0 = my->knots[my->p],
      k1 = *(my->knots.end()-my->p-1);

    std::vector<T> u = this->unique() ; 

    for ( typename  std::vector<T>::iterator it= u.begin(); it != u.end(); ++it )
        this->insert ( *it, i) ;
    my->p += i;

    // trim extra knots
    typename  std::vector<T>::iterator it =
        std::upper_bound( my->knots.begin(), my->knots.end(), k0 );
    my->knots.erase( my->knots.begin(), it-my->p-1 );

    it = std::lower_bound( my->knots.begin(), my->knots.end(), k1 );
    my->knots.erase(it+my->p+1, my->knots.end() );
}


template <class T>
void gsKnotVector<T>::degreeReduce(int const & i)
{
    reduceMultiplicity(i);
    my->p -= i;
}

template <class T>
void gsKnotVector<T>::degreeDecrease(int const & i)
{
    trim(i);
    my->p -= i;
}

template <class T>
void gsKnotVector<T>::degreeIncrease(int const & i)
{
    my->p += i;
    increaseMultFirst(i);
    increaseMultLast (i);
}

template <class T>
void gsKnotVector<T>::increaseMultiplicity(int const & i)
{
  const std::vector<T> u = this->unique();
  
  for (typename std::vector<T>::const_iterator itr = u.begin()+1;
       itr != u.end()-1; ++itr)
    this->insert( *itr, i);
}

template <class T>
void gsKnotVector<T>::reduceMultiplicity(int const & i)
{
    const std::vector<T> u = this->unique();
    
    for (const_iterator it = u.begin(); it != u.end(); ++it)
    {
        std::pair<iterator,iterator> itrs =
            equal_range(my->knots.begin(), my->knots.end(), *it);
        
        if ( itrs.second - itrs.first > i )
            my->knots.erase(itrs.first, itrs.first + i );
        else
            my->knots.erase(itrs.first, itrs.second);
    }
}

template <class T>
void gsKnotVector<T>::remove(T knot, int m)
{
    typename std::vector<T>::iterator itr =
        std::lower_bound(my->knots.begin(), my->knots.end(), knot);
    my->knots.erase(itr, itr+m);
}

template <class T>
gsMatrix<T> * gsKnotVector<T>::greville() const
{
  gsMatrix<T> * gr; 
  gr = new gsMatrix<T>( 1,this->size() - my->p - 1 );
  this->greville_into(*gr);
  return gr;
}

template <class T>
void gsKnotVector<T>::greville_into(gsMatrix<T> & result) const
{
  typename std::vector<T>::const_iterator itr = my->knots.begin() + 1;
  const int p = my->p;
  result.resize(1, this->size() -p-1 ) ; 
  unsigned i(1);
  
  if ( my->p!=0)
  {
      result(0,0)=  std::accumulate(itr, itr+p, T(0) ) / my->p;
      for (++itr; itr != my->knots.end()-p; ++itr, ++i )
      {
          result(0,i)=  std::accumulate( itr, itr+p, T(0) ) / my->p ;
          if ( result(0,i) == result(0,i-1) )
              result(0,i-1) -= 1e-10;// perturbe point to remain inside the needed support
      }
  }
  else
      std::copy(my->knots.begin(), my->knots.end()-1, result.data() );
}

template <class T>
T gsKnotVector<T>::greville(int i) const
{
    // TO DO: check special case as in greville_into
    typename std::vector<T>::const_iterator itr = my->knots.begin();
    return ( my->p!=0 ? 
	     std::accumulate( itr+1+i, itr+my->p+i+1, T(0.0) ) / my->p :
	     std::accumulate( itr+1+i, itr+my->p+i+1, T(0.0) )      );
}

template <class T>
void gsKnotVector<T>::set_degree(int p) 
{
    GISMO_ASSERT( my->knots.size() > static_cast<std::size_t>(2*p+1), "Not enough knots.");
    my->p=p;
}

template <class T>
int gsKnotVector<T>::degree() const 
{return my->p; }


template <class T>
std::vector<int> gsKnotVector<T>::multiplicities() const
{
    std::vector<int> mult;

    typename std::vector<T>::const_iterator itr = my->knots.begin();
    T last = *itr - 1;

    for (; itr != my->knots.end(); itr++)
    {
        if (last != *itr)
        {
            mult.push_back(1);
        }
        else
        {
            mult.back()++;
        }
        last = *itr;
    }

    return mult;
}


template <class T>
int gsKnotVector<T>::multiplicity(T knot) const
{
    std::pair<const_iterator,const_iterator> itrs =
        equal_range(my->knots.begin(), my->knots.end(), knot);

    return itrs.second - itrs.first ;
}

template <class T>
inline unsigned gsKnotVector<T>::multiplicityIndex(size_t const& i) const
{
    const T knot =  my->knots[i];
    size_t l = 0, r = 0;
    while (  l < i                && my->knots[i - (++l)] == knot ) { }
    while ( r+1< my->knots.size() && my->knots[i + (++r)] == knot ) { }
    
    return static_cast<unsigned>(l+r-1);
    
    /* equivalent:
       std::pair<const_iterator,const_iterator> itrs =
       equal_range(my->knots.begin(), my->knots.end(), my->knots[i] );
       return itrs.second - itrs.first ;
    */
}

template <class T>
int gsKnotVector<T>::multFirst() const
{
    typename std::vector<T>::const_iterator itr 
        = my->knots.begin();
    while ( *itr == *(itr+1) ) 
        itr++; 
    return itr - my->knots.begin() + 1;
}


template <class T>
int gsKnotVector<T>::multLast() const
{
    typename std::vector<T>::const_iterator itr 
        = my->knots.end()-1;
    while ( *itr == *(itr-1) ) 
        itr--; 
    return my->knots.end() - itr;
}


template <class T>
void gsKnotVector<T>::increaseMultFirst(int i)
{
    my->knots.insert( my->knots.begin() , i, my->knots.front() );
}

template <class T>
void gsKnotVector<T>::trim(int i)
{
    // same as:
    // pop_front(i); 
    // pop_back (i);
    my->knots.erase(my->knots.begin(), my->knots.begin()+i);
    my->knots.erase(my->knots.end()-i, my->knots.end()    );
}

template <class T>
void gsKnotVector<T>::increaseMultLast(int i)
{   
    my->knots.insert( my->knots.end()-1 , i, my->knots.back() );
    // Equivalent implementation
    // for (int k = 0; k!=i; ++k)
    //     my->knots.push_back( my->knots.back() );
}


/// Prints the object as a string.
template <class T>
std::ostream & gsKnotVector<T>::print(std::ostream &os) const
{ 
    os << "[ " ; 
    if ( size() > 2*my->p+8 )
    {
        for (const_iterator itr = this->begin(); itr != this->begin()+my->p+3; ++itr)
            os << *itr << " ";
        os << "... ";
        for (const_iterator itr = this->end()-my->p-3; itr != this->end(); ++itr)
            os << *itr << " ";
    }
    else
    {
        for (const_iterator itr = this->begin(); itr != this->end(); ++itr)
            os << *itr << " ";
    }
    os << "] (deg=" << degree()
       << ", size=" << size()
       << ", minSpan=" << minKnotSpanLength() 
       << ", maxSpan=" << maxKnotSpanLength()
       << ")";

    return os;
}

template <class T>
std::string gsKnotVector<T>::detail() const
{ 
    std::stringstream os;
    os << "[ " ; 
    for (const_iterator itr = this->begin(); itr != this->end(); ++itr)
    {
        os << *itr << " ";
    } 
    os << "]" <<"("<<my->p <<")";
    os << ". Size="<<size()<<", minSpan="<< minKnotSpanLength() <<", maxSpan="<< maxKnotSpanLength() <<"\n";
    return os.str();
}
  

namespace internal
{

/// Get a KnotVector from XML data
template<class T>
class gsXml< gsKnotVector<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsKnotVector<T>);
    static std::string tag () { return "KnotVector"; }
    //static std::string type() { return "Plain"; } // To do: add compact variant/or better make another tag CompactKnotVector
    static std::string type() { return ""; } // type=multiple / "default type"
    //static std::string type() { return "uniform"; } // start,end, num_interior, mult..,degre
    //static std::string type() { return "graded"; } // start,end, num_interior, mult..,degree, grading factor

    static gsKnotVector<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( ! strcmp( node->name(), "KnotVector"), "Invalid tag");
        // && node->first_attribute("type")->value() is "Plain");
        
        int p = atoi(node->first_attribute("degree")->value() );
        gsKnotVector<T> * kv = new gsKnotVector<T>;
        //gsWarn<<"Reading knots "<< node->value()<<"\n";
        
        std::istringstream str;
        str.str( node->value() );
        for (T knot; str >> knot;) 
            kv->push_back(knot);

        kv->set_degree(p);

        return kv;
    }

    static gsXmlNode * put (const gsKnotVector<T> & obj, gsXmlTree & data)
    {
        // Write the knot values
        std::ostringstream str;
        str << std::setprecision(FILE_PRECISION);

        for ( typename gsKnotVector<T>::const_iterator it = obj.begin();
              it != obj.end(); ++it )
        {
            str << *it <<" ";
        }

        // Make a new XML KnotVector node 
        gsXmlNode * tmp = internal::makeNode("KnotVector", str.str(), data);
        // Append the degree attribure
        str.str(std::string());// clean the ostream
        str<< obj.degree();
        tmp->append_attribute( makeAttribute("degree", str.str(),data) );
        return tmp;
    }

};

}// namespace internal


} // namespace gismo
