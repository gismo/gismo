
#pragma once

#include <gsCore/gsDomain.h>
#include <gsNurbs/gsCompactKnotVector.h>
#include <gsNurbs/gsCompactKnotVectorIter.h>
#include <gsNurbs/gsKnotVectorIter.h>

namespace gismo
{
    
template<class T> class gsKnotVectorPrivate;

  /** @brief
      A 1D knot vector, i.e., a sequence of non-decreasing knots,
      together with a degree.
      
      \tparam T coefficient type
  */  
template<class T>
class gsKnotVector : public gsDomain<T>
{
     
public:

    /// Type definitions
    typedef typename std::vector<T>::size_type size_t;
    // Iterators over all knots
    typedef typename std::vector<T>::iterator               iterator;
    typedef typename std::vector<T>::const_iterator         const_iterator;
    typedef typename std::vector<T>::reverse_iterator       reverse_iterator;
    typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;

    // Iterators over unique knots
    typedef gsKnotVectorIter<T, false>               uiterator;
    typedef gsKnotVectorIter<T, true >               const_uiterator;
    typedef std::reverse_iterator<uiterator>         reverse_uiterator;
    typedef std::reverse_iterator<const_uiterator>   const_reverse_uiterator;

    // Iterators over multiplicity sums
    //typedef typename ---     miterator;
    //typedef typename ---     const_miterator;


  /// Default empty constructor
  /// Results in empty knot-vector, and degree set to 0
  gsKnotVector();

  /// Copy constructor
  gsKnotVector(gsKnotVector const & other);

  /// Assignment operator
  gsKnotVector & operator=(gsKnotVector other)
    {
      this->swap(other);
      return *this;
    }

  /// Constructor with degree only
  explicit gsKnotVector(int p);

  /// Constructor with degree and size (=number of knots)
  gsKnotVector(int p, unsigned sz );

  /// Construct a knot vector
  /// \param u0 starting parameter
  /// \param u1 end parameter parameter
  /// \param interior number of interior knots
  /// \param mult_ends multiplicity at the two end knots
  /// \param mult_interior multiplicity at the interior knots
  /// \param degree multiplicity of the spline space
  gsKnotVector(T u0, T u1, unsigned interior, unsigned mult_ends=1, unsigned mult_interior=1, int degree=-1);

  /// @brief Construct an open knot vector from the given unique knots.
  /// \param knots sequence of distinct knots
  /// \param degree degree of a spline space
  /// \param regularity of spline space across the knots
  gsKnotVector(std::vector<T> const& knots, int degree, int regularity);

  /// @brief Construct a knot vector with the given degree and knots.
  gsKnotVector(int degree, std::vector<T> const& knots);

  /// Construct a knot vector by a given range
  /// \param deg
  /// \param start iterator pointing the first knot value
  /// \param end iterator pointing the last knot value
  gsKnotVector(int deg, const_iterator start, const_iterator end);

  gsKnotVector(gsCompactKnotVector<T> ckv);

  /// Destructor
  ~gsKnotVector();

public:
    
    void initUniform( T u0, T u1, unsigned interior, unsigned mult_ends, 
                     unsigned mult_interior = 1, int degree = -1);

    void initUniform(unsigned numKnots, unsigned mult_ends,
                    unsigned mult_interior = 1, int degree = -1);

    void initGraded(T u0, T u1, unsigned interior, int degree, 
                   T grading, unsigned mult_interior = 1);

    void initGraded(unsigned numKnots, int degree, 
                   T grading = 0.5, unsigned mult_interior = 1);

    void initClamped(T u0, T u1, int degree, unsigned interior = 0,
                     unsigned mult_interior = 1);

    void initClamped(int degree, unsigned numKnots = 2, unsigned mult_interior = 1);

    // initKnotsDegreeRegularity()
    // initKnotsDegreeRegularity()

    // setInteriorRegularity()
    // setDegree
    // setInteriorMultiplicity()
    // removeKnots()


public:

    /// Get an iterator to the beginning of the knot vector.
    iterator begin();

    reverse_iterator rbegin();
    
    /// Get an iterator to the beginning of the unique knots.
    uiterator ubegin();

    reverse_uiterator urbegin();
    
    /// Get a const-iterator to the beginning of the knot vector.
    const_iterator begin() const;

    const_reverse_iterator rbegin() const;

    /// Get a const-iterator to the beginning of the unique knots.
    const_uiterator ubegin() const;

    const_reverse_uiterator urbegin() const;
    
    /// Get an iterator to the end of the knot vector.
    iterator end();

    reverse_iterator rend();

    /// Get an iterator to the end of the unique knots.
    uiterator uend();

    reverse_uiterator urend();

    /// Get a const-iterator to the end of the knot vector.
    const_iterator end() const;

    const_reverse_iterator rend() const;

    /// Get a const-iterator to the end of the unique knots.
    const_uiterator uend() const;

    const_reverse_uiterator urend() const;


    
  /// \todo implement multiplicity sums iterator
  //const_iterator rbegin() const {  }

  /// Swap contents of this knot vector with \a other.
  void swap(gsKnotVector &other);

  /// Resize the knot vector to the given length.
  void resize(unsigned sz);

  /// Clear all the knots
  void clear();

  /// Print the knot vector to the given stream.
  std::ostream &print(std::ostream &os) const;

  /// Return a string with detailed information on the knot vector.
  std::string detail() const;

  /// Clone function. Used to make a copy of the (derived) geometry
  gsKnotVector * clone() const
    { return new gsKnotVector(*this); }

  /// Access the \a i-th knot.
  T  operator [] (size_t i) const;
  /// Access the \a i-th knot.
  T& operator [] (size_t i);

  /// Access the \a i-th knot, with bounds checking.
  T  at (size_t i) const;
  /// Access the \a i-th knot, with bounds checking.
  T& at (size_t i);
  
  /// Get the first knot.
  T first () const;
  /// Get the last knot.
  T last  () const;

  /// Test if two knot vectors are identical.
  bool operator==(const gsKnotVector<T> &other) const;
  /// Test if two knot vectors are different.
  bool operator!=(const gsKnotVector<T> &other) const;
  
  void push_back( T knot);

  void push_back( T knot, int mult);

  void push_front( T knot);

  void push_front( T knot, int mult);

  /// Return the size of the knot vector.
  int size() const;
  
  /// Get a reference to the underlying std::vector of knots.
  const std::vector<T>& get() const;

  gsVector<T> * getVector() const;
  
  /// Get a vector of all the unique values in the knot vector.
  std::vector<T> unique() const;

  /// Get a vector of all the unique values in the knot-vector
  /// starting from the i-th and up to the j-th knot
  std::vector<T> unique(size_t const i, size_t const & j) const;

  /// Breaks (for integration) \todo same as unique() ?
  std::vector<T> breaks() const;

  /// Find the index of the span in which value \a u lies.
  /// \returns the largest knot-index \c i such that \c knot[i] <= \a u
  unsigned findspan (T u) const;

  /// Find the position of the span in which value \a u lies.
  /// \returns an iterator to the largest knot \c k such that \c k <= \a u
  const_iterator findspanIter (T u) const ;
  
  gsMatrix<unsigned,1> * findspan (const gsMatrix<T,1> & u) const;
  
  //void scale (T u0, T u1);
  //void mirror ();
  
  /// Reverse the knot vector.
  ///
  /// Example [0, 0, 0.2, 1, 1] to [0, 0, 0.8, 1, 1]
  void reverse();

  /// @brief Insert a knot into the knot vector.
  /// \param knot parameter value of the new knot
  /// \param mult multiplicity of the new knot
  void insert(T knot, int mult=1);

  /// @brief Insert knots into the knot vector.
  /// \param knot parameter values of the new knots, stored in a std::vector
  /// \param mult multiplicity of the new knots
  void insert(std::vector<T> const & knots, int mult=1);

  /// True iff the given knot exists in the knot vector.
  bool has(T knot) const;

  /// Transforms the endpoints of the knot vector to [c, d].
  void transform(T c, T d);
  
  /// Insert a knot range into the knot vector.
  void append(const_iterator const & v0, const_iterator const & v1 ) ;

  /// Insert a knot range from another knot vector into the knot vector.
  void merge(gsKnotVector<T> other);
  
  /// Remove last knot from the knot vector.
  void pop_back(int i = 1);
  
  /// Remove first knot from the knot vector.
  void pop_front(int i = 1);
  
  /// Add a constant to all knots.
  void addConstant(T t);

  /// Shift the knot vector so that the first knot is equal to \a t.
  void setFirst(T t);
  
  /// Returns true iff the knot vector has uniform spacing
  bool isUniform() const;

    
  /// Returns the number of knot spans in the knot-vector
  int numKnotSpans() const;  

  /// Returns a vector containing the lenghts of the knot-spans
  std::vector<T> knotSpanLengths() const;

  /// Compute the length of the longest span
  T maxKnotSpanLength() const;

  /// Compute the length of the shortest span
  T minKnotSpanLength() const;

  /// Returns true iff the knot is open (ie. both endpoint
  /// multiplicities equal to m_p+1)
  bool isOpen() const;

  /// Returns the length of the first knot-interval
  T firstInterval() const;

  // Refine uniformly between "start" and "end" by adding \a numKnots
  // knot every two distinct knots
  // void uniformRefine(const T & start , const T & end, int numKnots = 1)

  /// Refine uniformly the knot vector between the interval staring at
  /// interval(0,0) and finishing at interval(0,1), by adding \a
  /// numKnots knot every two distinct knots
  void uniformRefine(gsMatrix<T> const & interval, int numKnots = 1);

  /// Refine uniformly the knot vector
  /// by adding \a numKnots knot every two distinct knots );
  void uniformRefine(int numKnots = 1);

  /// Compute the new knots needed for uniform refinement with the given number of knots per span and return them in \a result.
  void getUniformRefinementKnots(int knotsPerSpan, std::vector<T>& result);

  /// Elevate the degree
  void degreeElevate(int const & i = 1);

  /// Reduce the degree
  void degreeReduce(int const & i = 1);

  /// Reduce the degree keeping interior knots intact
  void degreeDecrease(int const & i = 1);

  /// Trim the knot vector from left and right by \a i knots
  void trim(int i = 1);

  /// Increase the multiplicity of all interior knots by i
  void increaseMultiplicity(int const & i = 1);

  /// Reduce the multiplicity of all interior knots by i
  void reduceMultiplicity(int const & i = 1);

  /// Remove the given knot from the knot vector.
    void remove(T knot, int m = 1);

  /// Get the Greville abscissae.
  gsMatrix<T> * greville() const;

  /// Get the Greville abscissae into \a result.
  void greville_into(gsMatrix<T> & result) const;

  /// Get the \a i-th Greville abscissa.
  T greville(int i) const;

  /// Set degree
  void set_degree(int p);

  /// Get degree
  int degree() const;

  /// Returns vector of multiplicities of the knots.
  std::vector<int> multiplicities() const;

  /// Get the multiplicity of the given knot.
  int multiplicity(T knot) const;

  /// Returns the multiplicity of the first knot
  int multFirst() const;

  /// Returns the multiplicity of the last knot
  int multLast() const;

  /// Increase the multiplicity of the first knot by \a i.
  void increaseMultFirst(int i = 1);

  /// Increase the multiplicity of the last knot by \a i.
  void increaseMultLast(int i = 1);

  /// Look at the supportIndex function.
  inline
  void supportIndex_into(const size_t& i, gsMatrix<unsigned>& result) const
  { GISMO_NO_IMPLEMENTATION }

  /// Get the unique knot index of the beginning and end of support of the
  /// i-th basis function.
  /// \param i index of the basis function
  inline gsMatrix<unsigned> supportIndex(const size_t& i) const
  { GISMO_NO_IMPLEMENTATION }

// Data members
private:
    gsKnotVectorPrivate<T> * my;

}; // class gsKnotVector


}// namespace gismo


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsKnotVector.hpp)
#endif
