/** @file gsCompactKnotVector.h

    @brief Provides declaration of the CompactKnotVector class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/
 
#pragma once

#include <numeric>

#include <gsNurbs/gsKnotVector.h>
#include <gsNurbs/gsCompactKnotVectorIter.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsDomain.h>
#include <gsUtils/gsSortedVector.h>

namespace gismo
{

/** 
    Class representing a knot vector of a B-Spline space.
    Representation is based on unique knots plus multiplicity information.

    Instead of saving eg. [0 0 0 0.5 1 1 1] we implement it as the value
    vector [0 0.5 1] plus another vector of multiplicity sums, in this
    case [3 4 7].

    Template parameter
    \param T is the coefficient type
    
    \ingroup Nurbs
*/
    
template<class T>
class gsCompactKnotVector : public gsDomain<T>
{
     
public:
    //Type definitions
    typedef typename std::vector<T>::size_type             size_t;
    // Iterators over unique knots
    typedef typename std::vector<T>::iterator               uiterator;
    typedef typename std::vector<T>::const_iterator         const_uiterator;
    typedef typename std::vector<T>::reverse_iterator       reverse_uiterator;
    typedef typename std::vector<T>::const_reverse_iterator const_reverse_uiterator;

    // Iterators over multiplicity sums
    typedef typename std::vector<unsigned>::iterator       miterator;
    typedef typename std::vector<unsigned>::const_iterator const_miterator;
    // Iterators over all knots
    typedef gsCompactKnotVectorIter<T, false>              iterator;
    typedef gsCompactKnotVectorIter<T, true >              const_iterator;
    
public:
    /// Default empty constructor
    gsCompactKnotVector() : gsDomain<T>(), m_p(0) { }

    /// Construct by an expanded knot vector
    gsCompactKnotVector(const gsKnotVector<T> & KV) ;
    
    explicit gsCompactKnotVector(int p) : gsDomain<T>(), m_p(p) { }

    gsCompactKnotVector(int p, unsigned sz ) : gsDomain<T>(), m_p(p) 
    { 
        if ( sz != 0 )
        {
            m_knots.resize(sz);
            m_mult_sum.resize(sz);
            m_mult_sum[0]= 1.0;
            std::vector<unsigned>::iterator i = m_mult_sum.begin();
            while ( ++i != m_mult_sum.end() ) *i = *(i-1) + 1;
        }
    }

    ~gsCompactKnotVector() { } // Constructor

    /// Construct a knot vector
    /// \param u0 starting parameter
    /// \param u1 end parameter parameter
    /// \param interior number of interior knots
    /// \param mult_ends multiplicity at the two end knots
    /// \param mult_interior multiplicity at the interior knots
    /// \param degree degree of the spline space
    gsCompactKnotVector(T const& u0, T const& u1, unsigned const & interior, unsigned const& mult_ends=1, unsigned const& mult_interior=1, int const& degree=-1) ;

    /// @brief Construct an open knot vector from the given unique knots.
    /// \param knots sequence of distinct knots
    /// \param degree degree of a spline space
    /// \param regularity of spline space across the knots
    gsCompactKnotVector(std::vector<T> const& knots, int degree, int regularity);

    /// Construct a open knot vector from two iterators of a gsKnotVector.
    /// \param deg degree
    /// \param start iterator pointing the first knot value
    /// \param end iterator pointing the last knot value
    gsCompactKnotVector(int const & deg, const_uiterator start, const_uiterator end);

    /// Construct a open knot vector from two iterators of a gsCompactKnotVector
    /// \param deg degree
    /// \param start iterator pointing the first knot
    /// \param end iterator pointing the last knot
    gsCompactKnotVector(int const & deg, const_iterator start, const_iterator end);

public:


    void initClamped(int degree, unsigned numKnots = 2, 
                     unsigned mult_interior = 1)
    {
        gsWarn<<"not finished.";
        GISMO_ASSERT( numKnots > 1 , "Not enough knots.");
        //initUniform(0.0, 1.0, numKnots - 2, degree + 1, mult_interior, degree );
        m_knots.resize(2);
        m_knots[0]=0.0;
        m_knots[1]=1.0;
        m_mult_sum.resize(2);
        m_mult_sum[0] = 
            m_mult_sum[1] = degree+1;
    }


public:

    void swap(gsCompactKnotVector& other)
    {
        std::swap(m_p, other.m_p);
        m_knots.swap(other.m_knots);
        m_mult_sum.swap(other.m_mult_sum);
    }

    void resize(unsigned sz)
    {
        if ( sz != 0 )
        {
            m_knots.resize(sz);
            m_mult_sum.resize(sz);
            m_mult_sum[0]= 1;	  
            std::vector<unsigned>::iterator i = m_mult_sum.begin();
            while ( ++i != m_mult_sum.end() ) *i = *(i-1) + 1;
        }
    }

    /// Compresses the knot-vector by making the knot-sequence strictly
    /// increasing
    void makeCompressed(const T & tol = 1e-7)
    {
        // TO DO: Check and improve
        typename std::vector<T>::iterator k = m_knots.begin();
        std::vector<unsigned>::iterator m = m_mult_sum.begin();
        while( k+1 != m_knots.end() )
        {
            if ( fabs( *k - *(k+1) ) <= tol )
            {
                *m = *(m+1) ;
                m_knots.erase(k+1);  
                m_mult_sum.erase(m+1);
            }
            else
            {
                k++;
                m++;
            }
        }
        //m_knots.resize( k - m_knots.begin() );
        //m_mult_sum.resize( m - m_mult_sum.begin() );
    }

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

    /// Print the knot vector to the given stream.
    std::ostream &print(std::ostream &os) const;

    /// Return a string with detailed information on the knot vector.
    std::string detail() const;

    /// Clone function. Used to make a copy of the (derived) geometry
    gsCompactKnotVector * clone() const
    { return new gsCompactKnotVector(*this); }

    /// Get a const-iterator to the beginning of the knotvector
    /// \return an iterator to the beginning of the knotvector
    const_iterator begin() const
    { return const_iterator(*this); }

    /// \todo implement reverse gsCompactKnotVector iterators
    //const_reverse_iterator rbegin() const { }

    const_uiterator ubegin() const
    { return m_knots.begin(); }

    const_reverse_uiterator urbegin() const
    { return m_knots.rbegin(); }

    const_miterator mbegin() const
    { return m_mult_sum.begin(); }
  
    /// Get a const-iterator to the end of the knotvector
    /// \return an iterator to the end of the knotvector
    const_iterator end() const
    { return const_iterator(*this,false); }

    const_uiterator uend() const
    { return m_knots.end(); }

    const_reverse_uiterator urend() const
    { return m_knots.rend(); }

    const_miterator mend() const
    { return m_mult_sum.end(); }
  
    /// Get an iterator to the beginning of the knotvector
    /// \return an iterator to the beginning of the knotvector
    iterator begin()
    { return iterator(*this); }

    uiterator ubegin()
    { return m_knots.begin(); }

    reverse_uiterator urbegin()
    { return m_knots.rbegin(); }

    miterator mbegin()
    { return m_mult_sum.begin(); }

    /// Get an iterator to the end of the knotvector
    /// \return an iterator to the end of the knotvector
    iterator end()
    { return iterator(*this,false); }

    uiterator uend()
    { return m_knots.end(); }

    reverse_uiterator urend()
    { return m_knots.rend(); }

    miterator mend()
    { return m_mult_sum.end(); }

    inline T    operator [] (const size_t & i) const 
    {
        const std::vector<unsigned>::const_iterator itr
            = std::upper_bound(m_mult_sum.begin(), m_mult_sum.end(), i);
        return m_knots[ itr - m_mult_sum.begin() ]; 
    }
    
    inline T  & operator [] (const size_t & i) 
    {
        const std::vector<unsigned>::const_iterator itr
            = std::upper_bound(m_mult_sum.begin(), m_mult_sum.end(), i);
        return m_knots[ itr - m_mult_sum.begin() ]; 
    }

    inline T    at (const size_t & i) const 
    {
        std::vector<unsigned>::const_iterator itr
            = std::upper_bound(m_mult_sum.begin(), m_mult_sum.end(), i);
        return m_knots[ itr - m_mult_sum.begin() ]; 
    }
    
    inline T  & at (const size_t & i) 
    { 
        std::vector<unsigned>::const_iterator itr
            = std::upper_bound(m_mult_sum.begin(), m_mult_sum.end(), i);
        return m_knots[ itr - m_mult_sum.begin() ]; 
    }

    /// Returns the multiplicity sum at span-index \a i
    inline unsigned knotsUntilSpan(const size_t & i) const 
    { return m_mult_sum[i]; }

    /// Returns the value of a unique index
    inline T  uValue(const size_t & i) const
    { return m_knots[i]; }

    /// Returns the left-most knot-index that is equal to the knot-index \a i
    inline unsigned  firstAppearance(const size_t & i) const
    { 
        if ( i < m_mult_sum[0])
            return 0;
        else
        {
            return m_mult_sum[
                std::upper_bound(m_mult_sum.begin(),m_mult_sum.end(),i)
                -m_mult_sum.begin() - 1];
        }
    }

    /// Get the cardinal index of the knot-index \a i
    /// This is the index of the knot without counting multiplicities
    /// \param i index of the knot
    inline unsigned cardinalIndex(size_t i) const;

    // returns the span-index of a parametric value
    // to do: rename to findCardinalIndex ?? // similarities with findElementIndex ..
    inline unsigned Uniquefindspan (T u) const ; 

    /// Get the element index of the knot-index \a i
    /// Same as cardinalIndex but its result is always a valid element
    /// e.g. Ending knot is regarded in the last element
    /// \param i index of the knot
    inline unsigned elementIndex(size_t i) const;

    /// returns the span-index of a parameter value
    inline unsigned findElementIndex (T u) const ; 
    
    /// Returns the first knot-index of cardinal index i
    inline unsigned firstKnotIndex(const size_t & i) const
    { 
        return ( i == 0 ? 0 : m_mult_sum[i-1] ) ;
    }

    /// Returns the first knot-index of cardinal index i
    inline unsigned  findFirstKnotInstance(const T & knot) const 
    { return firstKnotIndex(Uniquefindspan(knot)); }
    
    /// Returns the last knot-index of cardinal index i
    inline unsigned lastKnotIndex(const size_t & i) const
    { return m_mult_sum[i]-1; }

    /// Returns the last knot-index of cardinal index i
    inline unsigned findLastKnotInstance(const T & knot) const
    { return lastKnotIndex(Uniquefindspan(knot)); }
     
    T    first () const { return m_knots.front(); }  
    T    last  () const { return m_knots.back(); }  
    
    void push_back( const T & knot)  
    {
        GISMO_ASSERT(m_knots.empty() || knot>=m_knots.back(), "Knot out of range");
        
        if ( m_knots.back() < knot )
        {
            m_knots.push_back(knot); 
            m_mult_sum.push_back(1); 
        }
        else if  ( m_knots.back() == knot )
            ++m_mult_sum.back(); 
        else
            gsWarn<<"Invalid knot value.\n";
    }

    void push_front( const T & knot)  
    {   assert( m_knots.empty() || knot<=m_knots.front());
        m_knots.insert(m_knots.begin(),knot); }
    
    unsigned uSize() const { return (unsigned)m_knots.size(); }
    int        size() const { return m_mult_sum.back() ; }
    
    unsigned spans() const { return this->m_knots.size()-1; }
       
    std::vector<T> unique() const ; 

    // Breaks (for integration)
    std::vector<T> breaks() const 
    {
        std::vector<T> result(m_knots.begin(),m_knots.end());
        return result; 
    } 

    inline unsigned findspan (T u) const ;
   
    const_iterator findspanIter (T u) const;

    inline gsMatrix<unsigned,1> * findspan (const gsMatrix<T,1> & u) const ;
    
    void scale (const T& u0, const T & u1);
    void mirror ();

    /// Reverse the knot vector.
    void reverse()
    {
        // knot vector starts with first() and finish with last()
        const T ab = first() + last();
        
        std::reverse(m_mult_sum.begin(), m_mult_sum.end());
        std::reverse(m_knots.begin()   , m_knots.end());
        
        for (uiterator it = m_knots.begin(); it != m_knots.end(); ++it)
        {
            *it = ab - *it;
        }
    }
    
    /// Insert a knot into the knot vector.
    /// \param knot parameter value of the new knot.
    /// \param mult multiplicity of the new knot
    void insert(const T & knot, const int & mult=1);

    /// @brief Insert knots into the knot vector.
    /// \param knot parameter values of the new knots, stored in a std::vector
    /// \param mult multiplicity of the new knots
    void insert(std::vector<T> const & knots, const int & mult=1);

    /// True iff the knot exists in the vector
    /// \param knot parameter value 
    inline bool has(T knot) const 
    {return std::binary_search(m_knots.begin(), m_knots.end(), knot); }

    /// Transforms the endpoints of the knot vector to [c , d]
    void transform(T const & c, T const & d) 
    { 
        T a = m_knots.front();
        T b = m_knots.back();
        m_knots.front() = c;
        m_knots.back()  = d;
        
        for ( typename gsSortedVector<T>::iterator it= 
                  m_knots.begin()+1; it!= m_knots.end()-1; ++it)
            (*it) = c + ((*it) - a )*(d -c) /( b - a) ;
    }
    
    /// Insert a knot range into the knot vector.
    void append(const_uiterator const & v0, const_uiterator const & v1 ) 
	{ 
	    // TO DO
        //m_knots.insert( m_knots.end(), v0, v1 ); 
    }

    /// Insert a knot range into the knot vector.
    /// \param other domain from which to insert
    void merge(gsCompactKnotVector other ) 
	{ 
	    // TO DO
	}
    
    /// Remove last knot from the knot vector.
    void pop_back(int const& i = 1) 
    { 
        m_knots.erase( m_knots.end()-i,m_knots.end() ); 
        m_mult_sum.erase( m_mult_sum.end()-i,m_mult_sum.end() ); 
    }
    
    /// Remove first knot from the knot vector.
    void pop_front(int const& i = 1) 
    { 
        m_knots.erase( m_knots.begin(), m_knots.begin()+i ); 
        m_mult_sum.erase( m_mult_sum.begin(), m_mult_sum.begin()+i );
    }
    
    /// Add a constant to all knots
    void addConstant(const T & t)  
	{ std::transform(m_knots.begin(), m_knots.end(), m_knots.begin(),
                     std::bind2nd(std::plus<T>(), t));
	}

    /// Shift the knot vector so that first knot is equal to t
    void setFirst(const T & t) { addConstant( t - first() ); }

    /// Returns the number of knot spans in the knot-vector
    int numKnotSpans() const { return m_knots.size() - 1; }

    /// Refine the knot vector by adding \a numKnots equally spaced knots of
    /// multiplicity \a mul in between every two distinct knots
    void uniformRefine(int numKnots = 1, int mul=1) ;

    /// Refine elements pointed by the indices in \a spanIndices
    void refineSpans(const std::vector<unsigned> & spanIndices, int numKnots = 1)
    {
        GISMO_NO_IMPLEMENTATION
    }


    /// Compute the new knots needed for uniform refinement with the
    /// given number of knots per span and return them in \a result.
    void getUniformRefinementKnots(int knotsPerSpan, std::vector<T>& result, int mul=1) const;

    /// Elevate the degree
    void degreeElevate(int const & i = 1) ;

    /// Elevate the degree
    void degreeReduce(int const & i = 1) 
    { 
        GISMO_NO_IMPLEMENTATION
    }

    /// Reduce the degree keeping interior knots intact
    void degreeDecrease(int const & i = 1);
    
    /// Increase the degree keeping interior knots intact (add clamped knots only)
    void degreeIncrease(int const & i = 1);
    
    /// Remove a knot from the knot vector.
    /// \param knot parameter value of the knot.
    void remove(T const& knot);

    /// Get the Greville abscissae
    gsMatrix<T> * greville() const;

    /// Get the Greville abscissae
    void greville_into(gsMatrix<T> & result) const;

    /// Get the i-th Greville abscissae
    T greville(int i) const;

    /// Set degree
    void set_degree(int const& p) {m_p=p;}

    /// Get degree
    int degree() const {return m_p; }

    /// Returns a knot vector with the knots, with multiplicities
    gsKnotVector<T> expand() const;

    /// Returns vector of multiplicities of the knots.
    std::vector<int> multiplicities() const;

    /// Get the multiplicity of the knot
    /// \param knot value of a knot
    inline int multiplicity(const T & knot) const;

    /// Get the multiplicity sum at knot value \a  knot
    /// \param knot value of a knot
    inline int multiplicitySum(const T & knot) const;

    /// Increase the multiplicity of all interior knots by i
    void increaseMultiplicity(int const & i = 1)
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Increase the multiplicity of the first knot by \a i.
    void increaseMultFirst(int i = 1)
    {
        GISMO_NO_IMPLEMENTATION
    }
    
    /// Increase the multiplicity of the last knot by \a i.
    void increaseMultLast(int i = 1)
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Get the multiplicity sum at knot index \a i
    /// \param i knot-index
    inline int multiplicitySumIndex(size_t const& i) const;

    /// Get the multiplicity of the unique knot indexed i
    /// \param i index of the knot
    inline unsigned u_multiplicityIndex(size_t const & i) const;

    /// Get the multiplicity of the knot indexed i
    /// \param i index of the knot
    inline unsigned multiplicityIndex(size_t const & i) const;

    /// Look at the supportIndex function.
    inline
    void supportIndex_into(const size_t& i, gsMatrix<unsigned>& result) const;

    /// Get the unique knot index of the beginning and end of support of the
    /// i-th basis function.
    /// \param i index of the basis function
    inline gsMatrix<unsigned> supportIndex(const size_t& i) const;


private:
    template <class It>
    void init(int deg, It start, It end);

// Data members
private:

    // Distinct knot values
    gsSortedVector<T> m_knots;

    // Contains the sums of multiplicities
    std::vector<unsigned> m_mult_sum;
    
    int m_p;// Should the knot vector have a degree?.. we need it in findspan

}; // class gsCompactKnotVector


//////////////////////////////////////////////////
//////////////////////////////////////////////////


template <class T>
gsCompactKnotVector<T>::gsCompactKnotVector(T const& u0, T const& u1, unsigned const& interior, unsigned const& mult_ends, unsigned const& mult_interior, int const& degree)
{
    T h= (u1-u0)/(interior+1);
    
    m_knots.push_back(u0);
    m_mult_sum.push_back(mult_ends);

    for ( unsigned i=1; i<=interior; i++ )
    {
        m_knots.push_back( u0 + i*h );
        m_mult_sum.push_back(mult_interior + m_mult_sum.back() );
    }
    m_knots.push_back(u1);
    m_mult_sum.push_back(mult_ends + m_mult_sum.back() );

    if (degree==-1)   
        m_p=mult_ends-1;
    else
        m_p=degree;
}


template <class T>
gsCompactKnotVector<T>::gsCompactKnotVector(std::vector<T> const& knots, int degree, int regularity)
{
    typename std::vector<T>::const_iterator itr;
    int mult= degree - regularity ;

    m_knots.push_back( knots.front() );
    m_mult_sum.push_back(degree+1);

    for ( itr= knots.begin()+1; itr != knots.end()-1; ++itr )
    {      
        m_knots.push_back( *itr );
        m_mult_sum.push_back(mult + m_mult_sum.back() );
    }
    m_knots.push_back( knots.back() );
    m_mult_sum.push_back(degree+1 + m_mult_sum.back() );

    m_p=degree;
}

template <class T>
template <class It>
void gsCompactKnotVector<T>::init(int deg, It start, It end)
{
    // ASSUMES that start..end is sorted
    m_p = deg;
    m_knots.push_back(*start);
    m_mult_sum.push_back(1);
    //size_t i   = 0;
    It itr;
    for ( itr= start+1; itr != end; ++itr )
        if ( *itr == m_knots.back() )
            m_mult_sum.back() += 1;
        else
        {
            m_knots.push_back(*itr);
            m_mult_sum.push_back(1 + m_mult_sum.back() );
        }
}



template <class T>
gsCompactKnotVector<T>::gsCompactKnotVector(int const & deg, const_uiterator start, const_uiterator end)
{ 
    init(deg,start,end); 
}

template <class T>
gsCompactKnotVector<T>::gsCompactKnotVector(int const & deg, const_iterator start, const_iterator end)
{
    init(deg,start,end);
}
    
template <class T>
gsCompactKnotVector<T>::gsCompactKnotVector(const gsKnotVector<T> & KV)
{ 
    init(KV.degree(), KV.begin(), KV.end() );  
}


template <class T>
std::vector<T> gsCompactKnotVector<T>::unique() const 
{
    return m_knots;
}

template <class T>
inline unsigned gsCompactKnotVector<T>::findspan (T u) const
{
    return m_mult_sum[ findElementIndex(u) ] - 1;
}

template <class T>
typename gsCompactKnotVector<T>::const_iterator 
gsCompactKnotVector<T>::findspanIter (T u) const
{
    GISMO_ERROR("not implemented");
}

// to do:
// rename to: findSpanIndex
template <class T>
inline unsigned gsCompactKnotVector<T>::Uniquefindspan (T u) const
{
    unsigned low  = 0;
    unsigned high = (unsigned)(m_knots.size() - 1);

    GISMO_ASSERT( (u >= m_knots[low]) && ( u  <= m_knots[high] ), 
                  "The requested abscissae u="<<u<<" is not in the knot vector." );

    if (u == m_knots[high]) // remove ?
        return high-1;

    do
    {
        const unsigned mid = (low + high) >> 1;
        if ( u < m_knots[mid] )
            high = mid;
        else if (u >= m_knots[mid+1])
            low  = mid;
        else
            return mid;
    }
    while (true);
}

template <class T>
inline unsigned gsCompactKnotVector<T>::findElementIndex(T u) const
{
    unsigned low  = cardinalIndex(m_p);
    unsigned high = cardinalIndex(size()-m_p-1);

    GISMO_ASSERT( (u >= m_knots[low]) && ( u  <= m_knots[high]), 
                  "The requested abscissae u="<<u<<" is not in the knot domain." );

    if (u == m_knots[high])
        return high-1;

    do
    {
        const unsigned mid = (low + high) >> 1;
        if ( u < m_knots[mid] )
            high = mid;
        else if (u >= m_knots[mid+1])
            low  = mid;
        else
            return mid;
    }
    while (true);
}

template <class T>
inline gsMatrix<unsigned,1> * gsCompactKnotVector<T>::findspan (const gsMatrix<T,1> & u) const
{
    gsMatrix<unsigned,1> * fs = new gsMatrix<unsigned,1>(1, u.cols() );
  
    for( index_t i=0; i<u.cols(); i++ )
        (*fs)(0,i)= findspan( u(0,i) );

    return fs;
}


template <class T>
void gsCompactKnotVector<T>::scale (const T& u0, const T & u1)
{
    // TODO

}
 
template <class T>
void gsCompactKnotVector<T>::mirror()
{
    // TODO

}


template <class T>
void gsCompactKnotVector<T>::insert(const T & knot, int const & mult)
{
    int i = Uniquefindspan(knot);
    if ( knot == m_knots[i])
        std::transform(m_mult_sum.begin()+i, m_mult_sum.end(), m_mult_sum.begin()+i,
                       std::bind2nd(std::plus<unsigned>(), mult) );
    // Equivalent implementation:
    // for(unsigned int j = i; j < m_mult_sum.size();j++)
    //     m_mult_sum[j]+=mult;
    else
    {
        m_knots.insert(m_knots.begin()+i+1, knot );
        m_mult_sum.insert(m_mult_sum.begin()+i+1,m_mult_sum[i]+mult );
        std::transform(m_mult_sum.begin()+i+2, m_mult_sum.end(), m_mult_sum.begin()+i+2,
                       std::bind2nd(std::plus<unsigned>(), mult) );
        // Equivalent implementation:
        //for(unsigned int j = i+2; j < m_knots.size();j++)
        //    m_mult_sum[j] +=mult;
    }
}

template <class T>
void gsCompactKnotVector<T>::insert(std::vector<T> const & knots, int const & mult)
{
    for(typename std::vector<T>::const_iterator it=knots.begin();it!=knots.end();++it)
        insert(*it,mult);
}


template <class T>
void gsCompactKnotVector<T>::uniformRefine(int numKnots, int mul)
{
    const unsigned s0 = elementIndex(m_p),
        s1 = elementIndex(size()-m_p-1);
    // assume s0 == 0 or s0 == m_p for now
    // otherwise m_mult_sum needs to be updated as well

    std::vector<T> newKnots;
    getUniformRefinementKnots(numKnots, newKnots);

    for ( unsigned i = s0*numKnots; i < (s1+1)*numKnots; ++i )
    {
        this->insert(newKnots[i],mul);// to do: more efficient
    }

/*
  std::vector<T> ghosts(m_knots.begin(), m_knots.begin()+s0);
  unsigned l = m_p-1;
  for ( unsigned i = s0; i != 0 && l>0 ; --i )
  {
  for (int k = numKnots-1; k >=0 && l>0 ; --k)
  m_knots[l--] = newKnots[i*numKnots+k];

  if (l>0)
  {
  l--;
  m_knots[l] = ghosts[i];
  }
  }

  ghosts = std::vector<T>(m_knots.end()-s1, m_knots.end());
  l = s1;
  for ( unsigned i = s1; i != m_knots.size()-2 && l<s1+m_p ; ++i )
  {
  for (int k = 0; k <numKnots && l<s1+m_p ; --k)
  m_knots[l++] = newKnots[i*numKnots+k];

  if (l<s1+m_p)
  {
  l++;
  m_knots[l] = ghosts[i];
  }
  }
*/
}


template <class T>
void gsCompactKnotVector<T>::getUniformRefinementKnots(int knotsPerSpan, std::vector<T>& result, int mul) const
{
    const std::vector<T> & u = m_knots;
    result.clear();
    result.reserve((u.size() - 1) * knotsPerSpan*mul);

    for (std::size_t i = 0; i < u.size() - 1; ++i)
        for (int k = 1; k <= knotsPerSpan; ++k)
            result.insert(result.end(),mul,((knotsPerSpan+1-k) * u[i] + k * u[i+1]) / (knotsPerSpan + 1));
}


template <class T>
void gsCompactKnotVector<T>::degreeElevate(int const & i )
{
    for ( typename  std::vector<unsigned>::iterator it= m_mult_sum.begin(); 
          it != m_mult_sum.end(); ++it )
        *it += i;
    m_p += i;
}

template <class T>
void gsCompactKnotVector<T>::degreeIncrease(int const & i)
{
    m_p += i;
    increaseMultFirst(i);
    increaseMultLast (i);
}

template <class T>
void gsCompactKnotVector<T>::remove(T const& knot)
{

    typename std::vector<T>::iterator itr =
        std::lower_bound(m_knots.begin(), m_knots.end(), knot);
    m_mult_sum.erase(m_mult_sum.begin() + itr-m_knots.begin() ) ;
    //to do: update mults tail
    m_knots.erase(itr);
}

template <class T>
gsMatrix<T> * gsCompactKnotVector<T>::greville() const
{
    gsMatrix<T> * gr; 
    gr = new gsMatrix<T>( 1,this->size() - m_p - 1 );
    this->greville_into(*gr);
    return gr;
}

template <class T>
void gsCompactKnotVector<T>::greville_into(gsMatrix<T> & result) const
{
    const_iterator itr;
    result.resize( 1,this->size() - m_p - 1 ) ; 
    unsigned i(0);
    
    if ( m_p!=0)
        for ( itr= begin(); itr != end()-m_p-1; ++itr )
            result(0,i++)=  std::accumulate( itr+1, itr+m_p+1, T(0) ) / m_p ;
    else
        for ( itr= begin(); itr != end()-1; ++itr )
            result(0,i++)=  std::accumulate( itr+1, itr+1, T(0) ) ;
}

template <class T>
T gsCompactKnotVector<T>::greville(int i) const
{
    int multiplicities = 0;
    T sum = 0;
    int j = 0;
    if (m_mult_sum[j] < static_cast<unsigned>(i))
    {
        while (m_mult_sum[j] < static_cast<unsigned>(i))
        {
            j++;
        }
    }

    multiplicities = m_mult_sum[j] - i;
    sum += multiplicities * m_knots[j];

    const int m_p_1 = m_p + 2;
    for (std::size_t jj = j + 1; jj != m_mult_sum.size(); ++jj)
    {
        int mult = static_cast<int>(m_mult_sum[jj] - m_mult_sum[jj - 1]);

        if (m_p_1 < multiplicities + mult)
        {
            mult = m_p_1 - multiplicities;
        }

        sum += mult * m_knots[jj];
        multiplicities += mult;

        if (multiplicities == m_p_1)
        {
            break;
        }
    }


    return (m_p != 0 ? sum / m_p_1 : sum);

//    OLD CODE does not work, because minus is not implemented in
//    gsCompactKnotVectorIter
//
//    const_iterator itr = begin();
//    return ( m_p!=0 ?
//             std::accumulate( itr+1+i, itr+m_p+i+1, T(0.0) ) / m_p :
//             std::accumulate( itr+1+i, itr+m_p+i+1, T(0.0) )      );
}

template <class T>
gsKnotVector<T> gsCompactKnotVector<T>::expand() const
{
    gsKnotVector<T> result(m_p);

    for (const_iterator it= begin(); it!=end(); ++it)
        result.push_back(*it);
    
    return result;
}


template <class T>
std::vector<int> gsCompactKnotVector<T>::multiplicities() const
{
    std::vector<int> mult;

    mult.push_back(m_mult_sum[0]);

    std::vector<unsigned>::size_type indx;
    for (indx = 1; indx != m_mult_sum.size(); indx++)
    {
        mult.push_back(m_mult_sum[indx] - m_mult_sum[indx - 1]);
    }

    return mult;
}


template <class T>
inline int gsCompactKnotVector<T>::multiplicity(T const& knot) const
{
    if( knot == m_knots[0] )
        return m_mult_sum[0];
    else
    {
        typedef typename gsSortedVector<T>::const_iterator iter_t;
        typedef std::pair<iter_t,iter_t> result_t;
        result_t result=std::equal_range(m_knots.begin(),m_knots.end(),knot);
        if ( result.first==result.second )
            return 0;
        size_t j=result.first-m_knots.begin();
        return m_mult_sum[j] - m_mult_sum[j-1];
    }
}

template <class T>
inline int gsCompactKnotVector<T>::multiplicitySum(T const& knot) const
{
    return m_mult_sum[ Uniquefindspan(knot)];
}

template <class T>
inline int gsCompactKnotVector<T>::multiplicitySumIndex(size_t const& i) const
{
    return m_mult_sum[std::upper_bound(m_mult_sum.begin(),m_mult_sum.end(),i)
                      - m_mult_sum.begin()];
}


template <class T>
inline unsigned gsCompactKnotVector<T>::u_multiplicityIndex(size_t const& i) const
{
    if ( i == 0)
        return m_mult_sum[0];
    else
        return m_mult_sum[i] - m_mult_sum[i-1] ;
}

template <class T>
inline unsigned gsCompactKnotVector<T>::multiplicityIndex(size_t const& i) const
{
    if ( i < m_mult_sum[0])
        return m_mult_sum[0];
    else
    {
        int k = std::upper_bound(m_mult_sum.begin(),m_mult_sum.end(),i)-m_mult_sum.begin();
        return m_mult_sum[k] - m_mult_sum[k-1];
    }
}

template <class T>
inline unsigned gsCompactKnotVector<T>::cardinalIndex(size_t i) const
{
    const std::vector<unsigned>::const_iterator it = 
        std::upper_bound(m_mult_sum.begin(), m_mult_sum.end(), i);
    return  it-m_mult_sum.begin();
}

template <class T>
inline unsigned gsCompactKnotVector<T>::elementIndex(size_t i) const
{
    const std::vector<unsigned>::const_iterator it = 
        std::upper_bound(m_mult_sum.begin(), m_mult_sum.end(), i);
    return (unsigned)( it == m_mult_sum.end()-1 ? 
             m_mult_sum.size()-2 : it-m_mult_sum.begin() );
}

/// Prints the object as a string.
template <class T>
std::ostream & gsCompactKnotVector<T>::print(std::ostream &os) const
{ 
    typename gsKnotVector<T>::const_iterator itr;
    unsigned i=0;

    os << "[ " ; 
    for ( itr= m_knots.begin(); itr != m_knots.end(); ++itr )
    {
        os << *itr << "("<< u_multiplicityIndex(i)<<")"<<" ";
        ++i;
    } 
    os << "]" <<"("<<m_p <<")";
  
    return os;
}
  
template <class T>
std::string gsCompactKnotVector<T>::detail() const
{ 
    std::stringstream os;
    os << "[ " ; 
    for (const_iterator itr = this->begin(); itr != this->end(); ++itr)
    {
        os << *itr << " ";
    } 
    os << "]" <<"("<<m_p <<")";
    //os << ". Size="<<size()<<", minSpan="<< minKnotSpanLength() <<", maxSpan="<< maxKnotSpanLength() <<"\n";
    return os.str();
}

template <typename T>
inline
void gsCompactKnotVector<T>::supportIndex_into(const size_t& i,
                                               gsMatrix<unsigned>& result) const
{
    result.resize(1, 2);

    GISMO_ASSERT(i < (*(m_mult_sum.end() - 1) - m_p - 1),
                 "Index i is out of range");

    std::vector<unsigned>::const_iterator begin
        = std::upper_bound(m_mult_sum.begin(), m_mult_sum.end(), i);

    std::vector<unsigned>::const_iterator tmp =
        m_mult_sum.end() - begin  > m_p + 1 ?
        begin + m_p + 1 : m_mult_sum.end();

    std::vector<unsigned>::const_iterator end
        = std::upper_bound(begin, tmp, i + m_p + 1);

    result(0, 0) = begin - m_mult_sum.begin();
    result(0, 1) = end - m_mult_sum.begin();
}


template <typename T>
inline
gsMatrix<unsigned> gsCompactKnotVector<T>::supportIndex(const size_t& i) const
{
    gsMatrix<unsigned> result;
    supportIndex_into(i, result);
    return result;
}


}// namespace gismo



#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCompactKnotVector.hpp)
#endif
