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
    \brief Class representing a knot vector of a B-Spline space.
    
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
            m_mult_sum[0]= 1;
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
    gsCompactKnotVector(T const u0, T const u1, unsigned const interior, 
                        unsigned const mult_ends=1, 
                        unsigned const mult_interior=1, int const degree=-1) ;

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

    void initUniform( T u0, T u1, unsigned interior, unsigned mult_ends, 
                      unsigned mult_interior = 1, int degree = -1);

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
    unsigned cardinalIndex(size_t i) const;

    // returns the span-index of a parametric value
    // to do: rename to findCardinalIndex ?? // similarities with findElementIndex ..
    unsigned Uniquefindspan (T u) const ; 

    /// Get the element index of the knot-index \a i
    /// Same as cardinalIndex but its result is always a valid element
    /// e.g. Ending knot is regarded in the last element
    /// \param i index of the knot
    unsigned elementIndex(size_t i) const;

    /// returns the span-index of a parameter value
    unsigned findElementIndex (T u) const ; 
    
    /// Returns the first knot-index of cardinal index i
    inline unsigned firstKnotIndex(const size_t & i) const
    { 
        return ( i == 0 ? 0 : m_mult_sum[i-1] ) ;
    }

    /// Returns the first knot-index of the knot value \a knot
    inline unsigned  findFirstKnotInstance(const T & knot) const 
    { return firstKnotIndex(Uniquefindspan(knot)); }
    
    /// Returns the last knot-index of cardinal index i
    inline unsigned lastKnotIndex(const size_t & i) const
    { return m_mult_sum[i]-1; }

    /// Returns the last knot-index of the knot value \a knot
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

    unsigned findspan (T u) const ;
   
    const_iterator findspanIter (T u) const;

    gsMatrix<unsigned,1> * findspan (const gsMatrix<T,1> & u) const ;
    
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
    /// \param knots parameter values of the new knots, stored in a std::vector
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

    /// Reduce the degree
    void degreeReduce(int const & i = 1) 
    { 
        reduceMultiplicity(i);
        m_p -= i;
    }

    /// Reduce the degree keeping interior knots intact
    void degreeDecrease(int const & i = 1);
    
    /// \brief Returns true if the knot vector contains all knots of
    /// \a other, taking into account multiplicities
    bool contains(gsCompactKnotVector<T> & other);


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
    int multiplicity(const T & knot) const;

    /// Get the multiplicity sum at knot value \a  knot
    /// \param knot value of a knot
    int multiplicitySum(const T & knot) const;

    /// Increase the multiplicity of all interior knots by \a i
    void increaseMultiplicity(int const & i = 1)
    {
        miterator m  = m_mult_sum.begin()+1;
        unsigned r = 0;
        for (; m != m_mult_sum.end()-1; ++m)
        {
            r  += i;
            *m += r;
        }
        
        *m += r; // update last mult. sum (no increase in mult.)
    }

    /// Reduce multiplicity of all interior knots by i
    void reduceMultiplicity(int i = 1)
    {
        miterator m  = m_mult_sum.begin() + 1;
        uiterator it = m_knots.begin()    + 1;

        unsigned r = math::min(static_cast<unsigned>(i), *(m-1) );
        *(m-1) -= r;

        for (; it != m_knots.end();)
        {
            const int mult = *m - *(m-1) - r;

            if ( mult > i )
            {
                r  += i;
                *m -= r;
                ++it   ;
                ++m    ;
            }
            else // delete the knot (new multiplicity == 0)
            {
                r  += mult;
                it  = m_knots.erase(it);
                m   = m_mult_sum.erase(m);
            }
        }

        // Erase the first knot if needed
        m  = m_mult_sum.begin();
        it = m_knots   .begin();
        if ( *m == 0 )
        {
            m_knots.erase(it);
            m_mult_sum.erase(m);
        }
    }

    /// Increase the multiplicity of the first knot by \a i.
    void increaseMultFirst(int i = 1)
    {
        std::transform(m_mult_sum.begin(), m_mult_sum.end(), m_mult_sum.begin(),
                       std::bind2nd(std::plus<unsigned>(), i) );
    }
    
    /// Increase the multiplicity of the last knot by \a i.
    void increaseMultLast(int i = 1)
    {
        m_mult_sum.back() += i;
    }

    /// Get the multiplicity sum at knot index \a i
    /// \param i knot-index
    int multiplicitySumIndex(size_t const& i) const;

    /// Get the multiplicity of the unique knot indexed i
    /// \param i index of the knot
    unsigned u_multiplicityIndex(size_t const & i) const;

    /// Get the multiplicity of the knot indexed i
    /// \param i index of the knot
    unsigned multiplicityIndex(size_t const & i) const;

    /// Look at the supportIndex function.
    void supportIndex_into(const size_t& i, gsMatrix<unsigned>& result) const;

    /// Get the unique knot index of the beginning and end of support of the
    /// i-th basis function.
    /// \param i index of the basis function
    gsMatrix<unsigned> supportIndex(const size_t& i) const;

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




}// namespace gismo



#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCompactKnotVector.hpp)
#endif
