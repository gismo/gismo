/** @file gsKnotVector.h

    @brief Knot vector for B-splines.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris, A. Bressan, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsDomain.h>

#include <gsNurbs/gsKnotIterator.h>

namespace gismo
{

/** @brief Class for representing a knot vector.

 Consists of a vector of non-decreasing knots repeated according to
 their multiplicities and a vector of sums of multiplicities up to the
 corresponding unique knot. The repeated knots are stored in `m_repKnots`,
 whereas the multiplicity sum is stored in `m_multSum`.

# Iterators #

There are three means of iterating through the knots:
 + iterator, which is just an STL iterator over the repeated knots;
 + uiterator, which is similar but iterates over the unique knots
   (i.e., without multiplicities) and provides additional functionality;
 + smart_iterator, which iterates over the repeated knots and provides additional functionality.

# Terminology #

 + _repeated knots_ is the non-decreasing sequence, where each knot is repeated according to its multiplicity.

 + _unique knots_ is the increasing sequence of all the distinct values from _repeated knots_.

 + _starting knot_ is the degree-th knot from the beginning, i.e., `m_repKnots[m_deg]`.

 + _ending knot_ is the degree-th knot from the end, i.e., `m_repKnots[ m_repKnots.size() - m_deg - 1 ]`
 or, equivalently, `m_repKnots.end()[-m_deg-1]`.

 + _domain_ is the interval between the _starting knot_ and _ending knot_.

 + _knot interval_ \anchor knotInterval is the interval between two consecutive unique knots.
   It is closed from left and open from right, except for the interval
   [`*(domainUEnd()-1)`, `*domainUEnd()`], which is considered closed from both sides.

# Compatibility notes #

Earlier versions contained several functions that have been removed.

 + `uiterator findElement( const T u ) const` :  use `uFind( const T u )` instead;

 + `unsigned findspan( T u ) const` : use `iFind(u) - begin()` instead;

 + `iterator findspanIter( T u ) const` : use `iFind(u)` instead;

 + `int findElementIndex(T u) const` : use `uFind(u).uIndex()` instead;

 + `unsigned Uniquefindspan (T u) const` : use `uFind(u).uIndex()` instead;

 + `gsMatrix<unsigned,1> * findspan (const gsMatrix<T,1> & u) const` :
  call `iFind(u(0,i)) - begin()` for each `i` from `0` to `u.cols()`.

 + `int numKnotSpans() const` : use `uSize() - 1` instead;

 + `unsigned spans() const` : use `uSize() - 1` instead;

 */
template<typename T>
class gsKnotVector : public gsDomain<T>
{
public: // typedefs

    // IMPORTANT: mult_t must be a signed type, otherwise the
    // knotIterator business would not work.
    typedef index_t mult_t;

    // container types
    typedef std::vector<T>      knotContainer;
    typedef std::vector<mult_t> multContainer;

    // TODO: remove the non-const versions and think about the names.
    // E.g., reverse_uiterator clashes with urbegin().

    /// This iterator is an iterator over the repeated knots.
    /// Can be obtained by calling [r]begin() and [r]end().
    //typedef const T * iterator;
    typedef typename std::vector<T>::const_iterator iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;

    /// This iterator is an iterator over the unique knots. Has additional functions.
    /// Can be obtained by calling [r]ubegin() and [r]uend().
    typedef     internal::gsUKnotIterator<T> uiterator;
    typedef std::reverse_iterator<uiterator> reverse_uiterator;

    /// This iterator is an iterator over repeated knots with additional functions.
    /// Can be obtained by calling [r]sbegin() and [r]send().
    typedef           internal::gsKnotIterator<T>  smart_iterator;
    typedef std::reverse_iterator<smart_iterator>  reverse_smart_iterator;

    //compatibility iterator typedefs
    typedef uiterator        const_uiterator;
    typedef iterator         const_iterator;
    typedef reverse_iterator const_reverse_iterator;
    
public: // iterator ends

    /// Returns iterator pointing to the beginning of the repeated knots.
    iterator               begin()   const;
    /// Returns iterator pointing past the end of the repeated knots.
    iterator               end()     const;    
    /// Returns reverse iterator pointing past the end of the repeated knots.
    reverse_iterator       rbegin()  const;    
    /// Returns reverse iterator pointing to the beginning of the repeated knots.
    reverse_iterator       rend()    const;

    /// Returns an iterator pointing to the first appearance of the
    /// knot with unique index (ie. counted without repetitions, left
    /// ghosts mapped to negatives) equal to \a upos.
    iterator               beginAt(mult_t upos)   const;

    /// Returns an iterator pointing one past the last appearance of
    /// the knot with cardinal index (ie. counted without repetitions,
    /// left ghosts mapped to negatives) equal to \a upos.
    iterator               endAt(mult_t upos)     const;    
    
    /// Returns unique iterator pointing to the beginning of the unique knots.
    uiterator              ubegin()  const;
    /// Returns unique iterator pointing past the end of the unique knots.
    uiterator              uend()    const;
    /// Returns reverse unique iterator pointing past the end of the unique knots.
    reverse_uiterator      urbegin() const;
    /// Returns reverse unique iterator pointing to the beginning of the unique knots.
    reverse_uiterator      urend()   const;

    /// Returns the smart iterator pointing to the beginning of the repeated knots.
    smart_iterator         sbegin()  const;
    /// Returns the smart iterator pointing past the end of the repeated knots.
    smart_iterator         send()    const;    
    /// Returns the reverse smart iterator pointing past the end of the repeated knots.
    reverse_smart_iterator rsbegin() const;
    /// Returns the reverse smart iterator pointing to the beginning of the repeated knots.
    reverse_smart_iterator rsend()   const;

public: // constructors

    /// Empty constructor sets the degree to -1 and leaves the knots empty.
    gsKnotVector() : m_deg(-1)
    { }

    /// Constructs knot vector from the given \a knots (repeated
    /// according to multiplicities) and deduces the degree from the
    /// multiplicity of the endknots.
    explicit gsKnotVector(knotContainer knots, short_t degree = -1 );

    /// Swaps with \a other knot vector.
    void swap( gsKnotVector& other );

    /// Returns the knots as a matrix of size  1 x size()
    gsAsConstMatrix<T> asMatrix() const {return m_repKnots;}
    
public:
    /// Returns a pointer to a copy.
    gsKnotVector* clone() const;

public: // inserting and removing

    /// Inserts knot \a knot into the knot vector with multiplicity \a mult.
    void insert(T knot, mult_t mult = 1 );

    /// Inserts each knot between ibeg and iend with multiplicity 1.
    /// Knots can be repeated in the sequence but they must be sorted.
    /// The complexity is size() + (iend-ibeg) comparisons
    template<typename iterType>
    void insert( iterType ibeg, iterType iend )
    {
        //size_t numKnots = size() + (iend - ibeg);
        // GISMO_ENSURE( numKnots < std::numeric_limits<mult_t>::max(),
        //               "Too many knots." );
        knotContainer temp(m_repKnots.size()+(iend-ibeg) );
        std::merge( begin(), end(), ibeg, iend, temp.begin() );

        m_repKnots.swap(temp);

        rebuildMultSum(); // Possibly sub-optimal but the code is clean.

        GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities.");
    }

    //TODO: insert with uniqIter and multiplicity.

    /// Decreases multiplicity of the knot pointed by \a uit by \a
    /// mult; removes it altogether, should its multiplicity be
    /// decreased to zero.
    void remove( uiterator uit, mult_t mult = 1 );

    /// Decreases multiplicity of the knot \a knot by \a mult.
    void remove( const T knot, mult_t mult = 1 );

    //TODO: remove with two iterators.
    
public: // multiplicities

    /// Returns the multiplicity of the knot-value \a u (or zero if
    /// the value is not a knot)
    mult_t multiplicity( T u ) const;

    /// Returns the multiplicity of the first knot
    mult_t multFirst() const { return m_multSum.front(); }

    /// Returns the multiplicity of the last knot
    mult_t multLast() const { return m_multSum.back() - m_multSum.end()[-2]; }

    /// Returns the multiplicity of the knot number \a i (counter with
    /// repetitions).
    mult_t multiplicityIndex( mult_t i ) const;

public: // modifiers

    /// Scales and translates the knot vector so that its new
    /// endpoints are \a newBeg and \a newEnd.
    void affineTransformTo( T newBeg, T newEnd );

public: // queries

    /// Number of knots (including repetitions).
    inline size_t size() const { return m_repKnots.size(); }

    /// Number of unique knots (i.e., without repetitions).
    inline size_t uSize() const { return m_multSum.size(); }

    /// Provides the i-th knot (numbered including repetitions).
    const T& operator[](const mult_t i) const
    {
        GISMO_ASSERT( static_cast<size_t>(i) < m_repKnots.size(),
                      "Index "<<i<<" not in the knot vector.");
        return m_repKnots[i];
    }

    /// Provides the knot with unique index \a i
    const T& operator()(const mult_t i) const
    {
        return *( this->ubegin()+(numLeftGhosts()+i) );
    }

    /// Provides the knot with unique index \a i
    inline T  uValue(const size_t & i) const
    { return this->operator()(i); }

    /// Number of knot intervals inside domain.
    inline size_t numElements() const { return (domainUEnd() - domainUBegin()); }

public: // getters

    /// Returns unique knots.
    knotContainer unique () const
    {
        return knotContainer(this->ubegin(),this->uend());
    }

    /// Returns vector of multiplicities of the knots.
    multContainer multiplicities() const;

    /// Cast to the full vector of knot values (with repetitions).
    operator const knotContainer& () const {return m_repKnots;}

    /// Returns a pointer to the beginning of the vector of knots.
    const T * data() const {return m_repKnots.data(); }

    /// Returns a pointer to the beginning of the vector of running
    /// multiplicities.
    const mult_t * multSumData() const {return m_multSum.data(); }

public: // Findspan and value query

    /** \brief Returns the uiterator pointing to the knot at the
     * beginning of the _knot interval_ containing \a u.
     * Note that if `u == *domainUEnd()`, it returns the uiterator
     * `domainUEnd() - 1`. Cf. \ref knotInterval "knot interval". */
    uiterator uFind( const T u ) const;

    /** \brief Returns an iterator to the last occurrence of the knot
     * at the beginning of the _knot interval_ containing \a u.
     * Note that if `u == *domainEnd()`, it returns the iterator
     * `domainEnd() - 1`. Cf. \ref knotInterval "knot interval". */
    iterator iFind( const T u ) const;

    /** \brief Returns an iterator pointing to the first knot which
     * compares greater than \a u.
     *
     *  Unlike uLowerBound, the value pointed by the iterator returned
     *  by this function cannot be equal to \a u, only greater.
     */
    uiterator uUpperBound( const T u ) const;

    /** \brief Returns a uiterator pointing to the first knot which
     * does not compare less than \a u.
     * 
     *  Unlike uUpperBound, the value pointed by the iterator returned
     *  by this function may also be equivalent to \a u, and not only
     *  greater.
     */
    uiterator uLowerBound( const T u ) const;
    
public: // miscellaneous

    /// Print the knot vector to the given stream.
    /// TODO: Improve.
    std::ostream &print(std::ostream &os = gsInfo ) const;

    /// Checks whether the knot vector is in a consistent state
    bool check() const;

    typedef typename knotContainer::iterator nonConstIterator    ;
    typedef typename multContainer::iterator nonConstMultIterator;

    /// Returns an iterator pointing to the starting knot of the domain.
    iterator domainBegin() const
    {
        GISMO_ASSERT( size() > static_cast<size_t>(2*m_deg+1), "Not enough knots.");
        return begin() + m_deg;
    }

    /// Returns an iterator pointing to the end-knot of the domain.
    iterator domainEnd() const
    {
        GISMO_ASSERT( size() > static_cast<size_t>(2*m_deg+1), "Not enough knots.");
        return end() - (m_deg + 1);
    }

    /// Returns a smart iterator pointing to the starting knot of the
    /// domain
    smart_iterator domainSBegin() const
    { return sbegin() + m_deg; }

    /// Returns a smart iterator pointing to the ending knot of the
    /// domain
    smart_iterator domainSEnd() const
    { return send() - (m_deg + 1); }

    /// Returns a unique iterator pointing to the starting knot of the domain.
    uiterator domainUBegin() const
    {
        return domainSBegin().uIterator();
        // equivalent:
        //return ubegin() + domainSBegin().uIndex();
    }

    /// Returns a unique iterator pointing to the ending knot of the domain.
    uiterator domainUEnd() const
    {
        return domainSEnd().uIterator();
        // equivalent:
        //return ubegin() + domainSEnd().uIndex();
    }

    /// Checks, whether the given value is inside the domain.
    inline bool inDomain(const T u) const
    {
        return u >= *domainBegin() && u <= *domainEnd();
        // equivalent:
        // return u >= m_repKnots[m_deg] && u <= m_repKnots.end()[-m_deg-1];        
    }

    /// Removes the knots in the range [first,last)
    void erase(const mult_t first, const mult_t last);

    /// Removes the left-most \a numKnots from the knot-vector
    void trimLeft (const mult_t numKnots);
    
    /// Removes the right-most \a numKnots from the knot-vector
    void trimRight(const mult_t numKnots);

    /// Computes the number of left ghosts, i.e., of the knots to the
    /// left of the domain beginning.
    index_t numLeftGhosts() const
    {
        smart_iterator it(*this,0,0);
        it += math::min( (size_t)m_deg, size() );
        return std::distance( uiterator(*this,0,0), it.uIterator() );
    }

    /// Computes the number of right ghosts, i.e., of the knots to the
    /// right of the domain end.
    index_t numRightGhosts() const
    {
        return std::distance(domainUEnd(), uend()) - 1;
    }

public:

    /// Sanity check.
    static bool isConsistent(const knotContainer & repKnots,
                             const multContainer & multSums);

    /// Compare with another knot vector.
    inline bool operator== (const gsKnotVector<T>& other) const
    {
        return (m_repKnots == other.m_repKnots) &&
               (m_multSum  == other.m_multSum ) &&
               (m_deg == other.m_deg);
    }

private: // internals

    // Creates a vector of multiplicity sums from a scratch and m_repKnots.
    void rebuildMultSum();

private: // members

    // Knots including repetitions.
    knotContainer m_repKnots;

    // m_multSum[i] = cardinality of { knots <= unique knot [i] }.
    multContainer m_multSum;




//=======================================================================//
//+++++++ Here follow the functions I would like to get rid of. +++++++++//
//=======================================================================//



public: // Deprecated functions required by gsKnotVector.

    /// Sets the degree and leaves the knots uninitialized.
    explicit gsKnotVector(short_t degree)
    {
        m_deg = degree;
    }

    /// \param first starting parameter
    /// \param last end parameter parameter
    /// \param interior number of interior knots
    /// \param mult_ends multiplicity at the two end knots
    /// \param mult_interior multiplicity at the interior knots
    /// \param degree Degree of the knot. (default: -1)
    gsKnotVector( T first,
                  T last,
                  unsigned interior,
                  mult_t mult_ends=1,
                  mult_t mult_interior=1,
                  short_t degree = -1 );

    /// Constructs knot vector from the given degree and iterators marking its endpoints.
    /// \param deg Degree of the knot.
    /// \param begOfKnots iterator pointing to the beginning of the knots,
    /// \param endOfKnots iterator pointing past the end of the knots.
    /// The knots between \a begOfKnots and \a endOfKnots are assumed to be repeated
    /// according to their multiplicities and sorted.
    template<typename iterType>
    gsKnotVector(short_t deg, const iterType begOfKnots, const iterType endOfKnots)
    {
        insert(begOfKnots,endOfKnots);
        m_deg = deg;
    }

public:
    // TODO If stays, make it private.
    
    /// Resets the knot vector so that it has \a interior knots
    /// between \a first and \a last, each of them repeated \a
    /// mult_interior times (and the endpoints repeated \a mult_ends
    /// times) and the degree is equal to \a degree
    void initUniform( T first,
                      T last,
                      unsigned interior,
                      unsigned mult_ends,
                      unsigned mult_interior,
                      short_t degree=-1);

    /// Resets the knot vector so that it has \a numKnots knots
    /// between 0 and 1, each repeated \a mult_interior times (whereas
    /// 0 and 1 are repeated \a mult_ends times each) and the degree
    /// is equal to \a degree.
    void initUniform( unsigned numKnots,
                      unsigned mult_ends,
                      unsigned mult_interior = 1,
                      short_t degree = - 1)
    {
        initUniform(0.0, 1.0, numKnots - 2, mult_ends, mult_interior, degree );
    }


    /// Resets the knot vector so that its knots are graded towards 0
    /// according to the \a grading parameter and the endpoints are 0
    /// and 1.
    void initGraded(unsigned numKnots, int degree,
                    T grading = 0.5, unsigned mult_interior = 1)
    {
        initGraded( 0.0, 1.0, numKnots - 2, degree, grading, mult_interior);
    }

    /// Resets the knot vector so that its knots are graded.
    void initGraded(T u0, T u1, unsigned interior, int degree,
                    T grading, unsigned mult_interior = 1)
    {
        GISMO_ASSERT(u0<u1,"Knot vector must be an interval.");

        m_deg = degree;
        m_repKnots.reserve( 2*(m_deg+1) + interior*mult_interior );
        m_multSum .reserve(interior+2);

        const T h = (u1-u0) / (interior+1);

        m_repKnots.insert(m_repKnots.begin(), m_deg+1, u0);
        m_multSum .push_back(m_deg+1);

        for ( unsigned i=1; i<=interior; i++ )
        {
            m_repKnots.insert(m_repKnots.end(), mult_interior,
                              math::pow(i*h, 1.0/grading) );
            m_multSum .push_back( mult_interior + m_multSum.back() );
        }
        m_repKnots.insert(m_repKnots.end(), m_deg+1, u1);
        m_multSum .push_back( m_deg+1 + m_multSum.back() );
    }

    /// Returns the greville points of the B-splines defined on this
    /// knot vector.
    gsMatrix<T> * greville() const
    {
        gsMatrix<T> * gr;
        gr = new gsMatrix<T>( 1,this->size() - m_deg - 1 );
        this->greville_into(*gr);
        return gr;
    }

    // TODO Write properly.
     
    /// Compresses the knot-vector by making the knot-sequence strictly
    /// increasing
/*
    void makeCompressed(const T tol = 1e-7)
    {
        // TO DO: Check and improve, use tmp storage
        nonConstIterator     k = m_repKnots.begin();
        nonConstMultIterator m = m_multSum  .begin();
        T prev  = *k;
        while( k+1 != end() )
        {
            if ( math::abs( *k - *(k+1) ) <= tol )
            {
                *m = *(m+1) ;
                m_repKnots.erase(k+1);//all app.
                m_multSum.erase(m+1);
            }
            else
            {
                k++;
                m++;
            }
        }
        m_repKnots.resize( k - m_repKnots.begin() );
        m_multSum .resize( m - m_multSum .begin() );
    }
//*/



    /** \brief Returns true if every knot of \a other (counted with
      repetitions) is found within this knot-vector. Also returns
      true if \a other is empty.
    */
    bool includes(const gsKnotVector<T> & other) const;

    /** \brief Computes the difference between \a this knot-vector and \a other.
        
        All knots in \a this  which are not found in \a other (counted
        with repetitions), are stored in \a result
    */
    void difference(const gsKnotVector<T> & other,
                   std::vector<T>& result) const;

    /** \brief Computes the symmetric difference between \a this knot-vector and \a other.

        All knots which do not exist in both knot-vectors, \a this and
        \a other (counted with repetitions), are stored in \a result
    */
    void symDifference(const gsKnotVector<T> & other,
                       std::vector<T>& result) const;

    /** 
        \brief Computes the union of knot-vectors \a this and \a b.

        Example:
        gsKnotVector<> testKv1(0, 1, 2, 2, 1);		// 0 0 1/3 2/3 1 1
        gsKnotVector<> testKv2(0, 1, 1, 2, 2);		// 0 0 1/2 1/2 1 1 
        gsKnotVector<T> kv = kv1.knotUnion(kv2);

        Results in:  0 0 1/3 1/2 1/2 2/3 1 1
    */
    gsKnotVector knotUnion(const gsKnotVector & b) const;

    /** 
        \brief Computes the intersection of knot-vectors \a this and \a b.
     */
    gsKnotVector knotIntersection(const gsKnotVector & b) const;

    /// Better directly use affineTransformTo.
    void reverse();

    /// @brief Insert knots into the knot vector.
    /// \param knots parameter values of the new knots, stored in a std::vector
    /// \param mult multiplicity of the new knots     
    void insert( const knotContainer &knots, int mult = 1 )
    {
        for( int i = 0; i < mult; ++i)
            insert( knots.begin(), knots.end() ); // inefficient
    }

    /// Sets the degree to \a p.
    void set_degree(short_t p)
    {
        m_deg = p;
    }

    /// Writes unique indices of the knots of the endpoints of \a
    /// i - th B-spline defined on this knot vector to \a result.
    void supportIndex_into(const mult_t &i, gsMatrix<index_t>& result) const;
    // TODO: If the function stays, make a unit test from what is in unifiedKnotVector.cpp.

    /// Inserts \a numKnots (and each of them \a mult - times) between
    /// each two knots.
    void uniformRefine( mult_t numKnots = 1, mult_t mult = 1);

    /// Inserts \a knotsPerSpan knots between each two knots between \a begin and \a end.
    template<typename iterType>
    void refineSpans( iterType begin, iterType end, mult_t knotsPerSpan )
    {
        // We sort the input to be on the safe side.
        multContainer input;
        input.insert(input.begin(), begin, end); // Pitfall: using std::copy requires reserve or somesuch beforehand.
        std::sort(input.begin(), input.end());

        knotContainer newKnots;
        T segmentsPerSpan = knotsPerSpan + 1;
        T newKnot;
        T spanBegin;
        T spanEnd;
        for( typename multContainer::const_iterator it= input.begin();
             it != input.end();
             ++it )
            for( mult_t k = 1; k <= knotsPerSpan; ++k )
            {
                spanBegin =*(this->ubegin()+*it);
                spanEnd  =*(this->ubegin()+*it+1);
                newKnot = ( (segmentsPerSpan-k) * spanBegin + k * spanEnd ) / segmentsPerSpan;
                newKnots.push_back( newKnot );
            }

        insert(newKnots.begin(), newKnots.end());
    }

    /// Adds \a amount to all the knots.
    void addConstant( T amount );

public: // things required by gsKnotVector

    /// \param uKnots unique knots (assumed to be sorted),
    /// \param degree -> endknots have multiplicity \a degree + 1,
    /// \param regularity -> internal knots have multiplicity \a degree - \a regularity
    gsKnotVector( const knotContainer& uKnots,
                         int degree,
                         int regularity);

    /// Resets the knot vector so that it is uniform from 0 to 1 and
    /// the multiplicities of the endpoints are chosen according to
    /// the \a degree.
    void initClamped(int degree, unsigned numKnots = 2, unsigned mult_interior = 1);

    /// Resets the knot vector so that it is uniform and has clamped endknots.
    void initClamped(T u0, T u1, int degree, unsigned interior = 0,
                     unsigned mult_interior = 1)
    {
        return initUniform(u0, u1, interior, degree + 1, mult_interior, degree );
    }

    /// Returns the degree of the knot vector.
    int degree () const;


    /// Writes Greville abscissae of the B-splines defined on this
    /// knot vector to \a result.
    void greville_into(gsMatrix<T> & result) const;

    /// Writes the center points of each element (knot-span inside the
    /// domain) to \a result.
    void centers_into(gsMatrix<T> & result) const;

    /// Attempts to deduce the degree from the multiplicities of the
    /// endpoints (assuming clamped knots).
    int deduceDegree() const;

    /// Returns Greville abscissa of the \a i - the B-spline defined
    /// on the knot vector.
    T greville(int i) const;

    /// Get the first knot.
    T first () const
    {
        GISMO_ASSERT(this->size()>=1, "I need at least one knot.");
        return m_repKnots.front();
    }
    
    /// Get the last knot.
    T last  () const
    {
        GISMO_ASSERT(this->size()>=1, "I need at least one knot.");
        return m_repKnots.back();
    }

    /// See affineTransform().
    void transform(T c, T d)
    {
        affineTransformTo(c,d);
    }

    /// Because of type compatibility, cf. the other version.
    void refineSpans( const std::vector<unsigned> & spanIndices, mult_t knotsPerSpan = 1)
    {
        multContainer transformedIndices;
        transformedIndices.reserve(spanIndices.size());
        for( std::vector<unsigned>::const_iterator it = spanIndices.begin();
             it != spanIndices.end();
             ++it )
            transformedIndices.push_back(static_cast<mult_t>(*it));

        return refineSpans( transformedIndices, knotsPerSpan );
    }

    /// Inserts \a knotsPerSpan knots into the spans corresponding to
    /// indices listed in \a span Indices.
    void refineSpans( const multContainer & spanIndices, mult_t knotsPerSpan = 1);

    /// Increase the degree keeping interior knots intact (add clamped knots only).
    /// Cf. degreeElevate.
    void degreeIncrease(int const & i = 1)
    {
        // update knots
        m_repKnots.reserve(size()+2*i);
        m_repKnots.insert(m_repKnots.begin(), i, m_repKnots.front() );
        m_repKnots.insert(m_repKnots.end()  , i, m_repKnots.back()  );

        // update multiplicity sum
        std::transform(m_multSum.begin(), m_multSum.end(), m_multSum.begin(),
                       GS_BIND2ND(std::plus<mult_t>(), i) );
        m_multSum.back() += i;

        m_deg += i;
    }

    /// Inverse of degreeIncrease.
    void degreeDecrease(int const & i = 1 )
    {
        // note: multiplicities are already updated after each call
        remove( ubegin()  , i );
        remove( uend() - 1, i );
        m_deg -= i;
    }

    /// Increase the multiplicity of all the knots by \a i. If \a
    /// boundary is set to \em false, the boundary knots are not adjusted
    void increaseMultiplicity(const mult_t i = 1, bool boundary = false);

    /// Reduce the multiplicity of all knots by \a i. If \a boundary
    /// is set to \em false, the boundary knots are not adjusted
    void reduceMultiplicity(const mult_t i = 1, bool boundary = false);

    /// Adds the knots between \a begin and \a end to the knot vector.
    template<typename iterType>
    void append ( iterType begin, iterType end )
    {
        insert( begin, end );
    }

    //void _stretchEndKnots()

    /// Get a reference to the underlying std::vector of knots.
    const knotContainer& get() const // to be removed since we have implicit cast
    {
        return m_repKnots;
    }

public: // Deprecated functions required by gsCompactKnotVector.

    /// True iff the knot exists in the vector
    /// \param knot parameter value.
    inline bool has(T knot) const
    {
        return std::binary_search( ubegin(), uend(), knot);
        // equivalent:
        // return 0 != multiplicity(knot);
    }

    /// Get the multiplicity of the unique knot indexed \a i
    /// \param i index of the knot (without repetitions)
    unsigned u_multiplicityIndex(size_t const & i) const
    {
        return (ubegin() + i).multiplicity();
    }

    /// Returns the first knot-index of cardinal index \a i
    /// (i.e. counted without repetitions)
    inline unsigned firstKnotIndex(const size_t & i) const
    {
        return (ubegin()+i).firstAppearance();
        // equivalent:
        // return 0 == i ? 0 : m_multSum[i-1];
    }
     
    /// Returns the last knot-index of cardinal index \a i
    /// (i.e. counted without repetitions)
    inline unsigned lastKnotIndex(const size_t & i) const
    {
        return (ubegin()+i).lastAppearance();
        // equivalent:
        // return m_multSum[i] - 1;
    }
     
    /// Returns the multiplicity sum at unique index \a i
    /// (i.e. counted without repetitions)
    inline unsigned knotsUntilSpan(const size_t & i) const
    {
        return (ubegin()+i).multSum();
    }

    /// Compares with another knot vector.
    bool operator != (const gsKnotVector<T>& other) const
    {
        return ! ((*this)==other);
    }    

    /// Returns the value of the \a i - th knot (counted with repetitions).
    inline T at (const size_t & i) const
    {
        return m_repKnots.at(i);
    }

    /// Checks whether the knot vector is uniform.
    bool isUniform(T tol = 1e-9) const
    {
        const T df = *(ubegin() + 1) - *ubegin();
        for( uiterator uit = ubegin() + 1; uit != uend(); ++uit )
            if( math::abs(*uit - (*uit-1) - df) > tol )
                return false;
        return true;
    }
     
    /// Returns true iff the knot is open (i.e., both endpoints
    /// have multiplicities equal to degree+1).
    bool isOpen() const
    {
        const int dp1 = m_deg + 1;
        return (multFirst() == dp1 &&
                multLast () == dp1 );
        // equivalent
        //return ( ubegin  .multiplicity() == dp1 &&
        //         (--uend).multiplicity() == dp1 );
        // equivalent
        //return m_multSum.front() == dp1 && 
        //    m_multSum.back() - m_multSum.end()[-2] == dp1;
    }

    /// Returns unique knots of the domain (i.e., including the endpoints of the domain).
    // \sa gsDomain
    virtual knotContainer breaks() const
    {
        return knotContainer(domainUBegin(), domainUEnd()+1);
    }

public: // others

    /// Compute the new knots needed for uniform refinement with the
    /// given number of knots per span and return them in \a result.
    /// TODO: Think who is in charge of uniform refinement: basis or knot vector?
    void getUniformRefinementKnots(mult_t knotsPerSpan, knotContainer& result, 
                                   mult_t mult = 1) const;

    /// Return a string with detailed information on the knot vector.
    std::string detail() const;

    /// Returns the maximum interval length of the knot sequence
    T maxIntervalLength() const;

    /// Returns the minimum interval length of the knot sequence
    T minIntervalLength() const;

    /// Elevate the degree. I.e., increase the multiplicity of all the
    /// knots by one and increment the degree.
    void degreeElevate(const short_t & i = 1);

    /// Converse to degreeElevate.
    void degreeReduce(const short_t & i);

    /// Coarsen the knot vector: for each group of \a knotRemove
    /// consecutive interior unique knots, reduce their multiplicity
    /// by \a mul (or remove completely, if \a mul = -1), then skip
    /// \a knotSkip consecutive unique knots and go on
    ///
    /// @param knotRemove    Number of consecutive knots to remove
    /// @param knotSkip      Number of consecutive knots which stay untouched in between
    /// @param mul           Multiplicity drop for each processed knot (or -1, full remove)
    std::vector<T> coarsen(index_t knotRemove = 1, index_t knotSkip = 1, mult_t  mul = -1);

public: // members

    // TODO remove!
    short_t m_deg;
};

template<typename T>
std::ostream& operator << (std::ostream& out, const gsKnotVector<T> KV )
{
    KV.print(out);
    return out;
}

} // namespace gismo


// *****************************************************************
#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsKnotVector.hpp)
#else
#ifdef gsKnotVector_EXPORT
#include GISMO_HPP_HEADER(gsKnotVector.hpp)
#undef  EXTERN_CLASS_TEMPLATE
#define EXTERN_CLASS_TEMPLATE CLASS_TEMPLATE_INST
#endif
namespace gismo
{
EXTERN_CLASS_TEMPLATE gsKnotVector<real_t>;
//EXTERN_CLASS_TEMPLATE internal::gsXml<gsKnotVector<real_t> >;
}
#endif
// *****************************************************************
