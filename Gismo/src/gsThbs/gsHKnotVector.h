 
#pragma once

#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>

#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsCompactKnotVector.h>
#include <gsThbs/gsHKnotVectorIter.h>

namespace gismo
{

  /** 
      Class representing a 1D Hierarchical knot vector.

      Template parameter
      \a T is the coefficient type
  */
  
  
template<class T>
class gsHKnotVector : public gsDomain<T>
{
     
public:

  /// Type definitions
  typedef typename std::vector<T>::size_type size_t;

  // Iterators over on all knots of any levels (with/without mults)
  typedef typename gsCompactKnotVector<T>::uiterator       global_uiterator;
  typedef typename gsCompactKnotVector<T>::const_uiterator const_global_uiterator;
  typedef typename gsCompactKnotVector<T>::iterator        global_iterator;
  typedef typename gsCompactKnotVector<T>::const_iterator  const_global_iterator;
  typedef typename gsCompactKnotVector<T>::miterator       global_miterator;
  typedef typename gsCompactKnotVector<T>::const_miterator const_global_miterator;

  // Iterators over on all knots of current level (with/without mults)
  typedef gsHKnotVectorUIter<T, false> uiterator;
  typedef gsHKnotVectorUIter<T, true > const_uiterator;
  typedef gsHKnotVectorIter<T , false> iterator;
  typedef gsHKnotVectorIter<T , true > const_iterator;

  /// Default empty constructor
    gsHKnotVector() : m_knots(),  maxlvl(0), lvl(0) { };

    gsHKnotVector( const gsKnotVector<T> & KV, const unsigned & nlevels = 1) : m_knots(KV), maxlvl(nlevels-1), lvl(0), step(maxlvl) { };

  virtual ~gsHKnotVector() { }; // Constructor

///// LEVEL-related functions

    // to do: rename to: view(i) -- for "viewing" ith-level knot vector
    void setLevel(unsigned const & i) const 
    { assert( i<= maxlvl); lvl=i, step=1 >> (maxlvl-i);};

    int stride() const { return step; };

    void setMaxLevel(unsigned const & i) 
    { assert( i>maxlvl);  maxlvl=i, step=i-lvl;} ; // + refine knots i - maxlvl and mult everything by 2^( i - maxlvl)

    void addLevel() { ++maxlvl;} ;// + refine knots once, and mult everything by 2

    inline T  uGlobalKnot(const size_t & i) { return m_knots[i]; };
    int globalSize() { return m_knots.size(); };
    /// Get the multiplicity of the unique knot indexed i
    /// \param i Global index of the knot
    inline unsigned multiplicityIndexGlobal(size_t const & i) const // OK
        {return m_knots.multiplicityIndex( i ); };
    
    std::vector<T> unique() const 
        { return m_knots.unique(); };

    /// Get degree
    int degree() const {return m_knots.degree(); };
    
    void scale (const T& u0, const T & u1);

    void transform(T const & c, T const & d) 
    { 
        T a = m_knots.front();
        T b = m_knots.back();
        m_knots.front() = c;
        m_knots.back()  = d;
        
        for ( typename gsCompactKnotVector<T>::iterator it= 
           m_knots.begin()+1; it!= m_knots.end()-1; ++it)
            (*it) = c + ((*it) - a )*(d -c) /( b - a) ;
    };

    /// Shift the knot vector so that first knot is equal to t
    void setFirst(const T & t) { addConstant( t - first() ); };

    /// Add a constant to all knots
    void addConstant(const T & t)  
	{ std::transform(m_knots.begin(), m_knots.end(), m_knots.begin(),
			 std::bind2nd(std::plus<T>(), t));
	};

  /// Elevate the degree
    void degreeElevate(int const & i = 1) { } ;


//////////////// All these refer to level lvl;


    inline T    operator [] (const size_t & i) const
        { 
            return this->at(i) ; 
        };

    inline T    at (const size_t & i) const { // OK, but requures uSpanIndex
        return (*this)[ uSpanIndex(i) ]; 
    };


    inline unsigned uSpanIndex(const size_t & i) // TODO transform res to lvl
        { 
            // requires bin search using strind
            return 0;
        };
  
    // returns the value of a unique knot
    inline T  uValue(const size_t & i) { return m_knots.uKnot(i << (maxlvl-lvl) ); }; //OK
    
    inline unsigned firstKnotIndex(const size_t & i) // TODO transform res to lvl
        { 
            return m_knots.firstKnotIndex( i>> (maxlvl-lvl) );
        };
    
    inline unsigned  lastKnotIndex(const size_t & i) // TODO  transform res to lvl
        { return m_knots.lastKnotIndex( i >> (maxlvl-lvl) ); };
    
    inline unsigned  firstInstance(const T & knot) // TODO  transform res to lvl
        { return m_knots.firstInstance(knot) >> (maxlvl-lvl) ; };
    
    inline unsigned  lastInstance(const T & knot) // TODO  transform res to lvl
        { return m_knots.lastInstance(knot) >> (maxlvl-lvl) ; };
    
    inline unsigned findspan (T u) const {return m_knots.findspan(u); }; // = findspan, transform res to lvl
    
    inline unsigned Uniquefindspan (T u) const // OK
        {return m_knots.Uniquefindspan(u) >> (maxlvl-lvl) ; };
    
    /// Get the multiplicity of the knot
    /// \param knot value of a knot
    inline unsigned multiplicity(const T & knot) const // OK .. substruct RR
        { return m_knots.multiplicity(knot);};
    
    /// Get the multiplicity of the unique knot indexed i
    /// \param i index of the knot
    inline unsigned multiplicityIndex(size_t const & i) const // -- substruct RR
        {return m_knots.multiplicityIndex( i >> (maxlvl-lvl) );};
    
    T    first () const { return m_knots.first(); };  // OK
    T    last  () const { return m_knots.last(); };  // OK
    
    int size() const { return m_knots.size() ; } ; // TODO  transform res to lvl 
    
    /// Get the Greville abscissae
    gsMatrix<T> * greville() const { return m_knots.greville(); }; // TODO  reimplement with operator []

//////////////////////////////////////////////////
// Virtual member functions required by the base class
//////////////////////////////////////////////////

  /// Clone function. Used to make a copy of the object
  gsHKnotVector * clone() const
    { return new gsHKnotVector(*this); };

    std::ostream &print(std::ostream &os) const { return m_knots.print(os); }; // TODO  transform res to lvl 

    

private:

    gsCompactKnotVector<T> m_knots;

    unsigned mult_residue;

    unsigned maxlvl;
    mutable int lvl;
    mutable int step;

}; // class gsHKnotVector


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsHKnotVector.hpp)
#endif
