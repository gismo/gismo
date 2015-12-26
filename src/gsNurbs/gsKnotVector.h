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

// TODO: Formatting
/** @brief Class for representing a knot vector.

 Consists of a vector of non-decreasing knots repeated according to
 their multiplicities and a vector of sums of multiplicities up to the
 corresponding unique knot.

 There are three means of iterating through the knots:
 1) iterator, which is just an STL iterator over the repeated knots;
 2) uiterator, which is similar but iterates over the unique knots (i.e., without multiplicities) and provides additional functionality;
 3) smart_iterator, which iterates over the repeated knots and provides additional functionality.

Terminology: TODO update

Domain is the interval between the degree-th knot from the beginning
and degree-th knot from the end (note that we number from zero).

Span is the interval between two consecutive unique knots. It is
closed on the left and open on the right, except for the last one in
the domain, which is closed. Span index is the knot index of the last
knot (including repetitions) at the beginning of the span.

Element is a span inside the domain. Their numbering is shifted in
comparison to that of spans so that the element with index 0 starts in
the beginning of the domain.

Example:

Consider the knot vector (0, 0, 0, .2, .3, .7, 1, 1, 1) with degree 2.
The domain is [0, 1]. The spans are [0, .2), [.2, .3), [.3, .7) and
[.7, 1]. Their span indices are 0, 1, 2 and 3, respectively. On a
clamped knot vector the spans and the elements are the same.

Another example:

Consider knot vector (0, .1, .2, .3, .4, .5, .6, .7) with degree 2.
The domain is [.2, .5]. The spans are [0, .1), [.1, .2), [.2, .3),
[.3, .4), [.4, .5], (.5, .6) and [.6, .7] with span indices 0 to
6. Notice that [.4, .5] and [.6, .7] are closed. The elements are [.2,
.3), [.3, .4) and [.4, .5] with element indices 0, 1 and 2.

 
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
    typedef typename knotContainer::const_iterator           iterator;
    typedef typename knotContainer::const_reverse_iterator   reverse_iterator;

    /// This iterator is an iterator over the unique knots. Has additional functions.
    /// Can be obtained by calling [r]ubegin() and [r]uend().
    typedef     internal::gsUKnotIterator<T> uiterator;
    typedef std::reverse_iterator<uiterator> reverse_uiterator;

    /// This iterator is an iterator over repeated knots with additional functions.
    /// Can be obtained by calling [r]sbegin() and [r]send().
    typedef           internal::gsKnotIterator<T>  smart_iterator;
    typedef std::reverse_iterator<smart_iterator>  reverse_smart_iterator;

    //compatibility iterator typedefs
    typedef uiterator const_uiterator;
    typedef iterator  const_iterator;
    
public: // iterator ends

    /// Returns iterator pointing to the beginning of the repeated knots.
    iterator               begin()   const;
    /// Returns iterator pointing past the end of the repeated knots.
    iterator               end()     const;    
    /// Returns reverse iterator pointing past the end of the repeated knots.
    reverse_iterator       rbegin()  const;    
    /// Returns reverse iterator pointing to the beginning of the repeated knots.
    reverse_iterator       rend()    const;
    
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
    explicit gsKnotVector( gsMovable<knotContainer > knots, int degree = -1 );

    /// Swaps with \a other knot vector.
    void swap( gsKnotVector& other );

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
    inline size_t size() const {return m_repKnots.size(); }

    /// Number of unique knots (without repetitins).
    inline size_t uSize() const { return m_multSum.size(); }

    /// Provides i-th knot (numbered including repetitions).
    const T& operator[]( const mult_t i ) const
    {
        GISMO_ASSERT( static_cast<size_t>(i) < m_repKnots.size(),
                      "Index " << i << " not in the knot vector." );
        return m_repKnots[i];
    }

    /// Cast to the full vector of knot values (with repetitions).
    operator const knotContainer& () const {return m_repKnots;}

    /// Returns a pointer to the beginning of the vector of knots.
    const T * data() const {return m_repKnots.data(); }

    /// Returns a pointer to the beginning of the vector of running
    /// multiplicities.
    const mult_t * multSumData() const {return m_multSum.data(); }

public: // miscellaneous

    /// Print the knot vector to the given stream.
    /// TODO: Improve.
    std::ostream &print(std::ostream &os = std::cout ) const;

    /// Checks whether the knot vector is in a consistent state
    bool check() const;

private: // iterator typdefs

    typedef typename knotContainer::iterator nonConstIterator    ;
    typedef typename multContainer::iterator nonConstMultIterator;

    /// Returns a smart iterator pointing to the starting knot of the
    /// domain
    iterator domainBegin() const
    {
        GISMO_ASSERT( size() > static_cast<size_t>(2*m_deg+1), "Not enough knots.");
        return begin() + m_deg;
    }

    /// Returns a smart iterator pointing to the end-knot of the
    /// domain
    iterator domainEnd() const
    {
        GISMO_ASSERT( size() > static_cast<size_t>(2*m_deg+1), "Not enough knots.");
        return end() - (m_deg + 1);
    }

    /// Returns a smart iterator pointing to the starting knot of the
    /// domain
    smart_iterator domainSBegin() const
    { return sbegin() + m_deg; }

    /// Returns a smart iterator pointing to the end-knot of the
    /// domain
    smart_iterator domainSEnd() const
    { return send() - (m_deg + 1); }

    /**
       \brief Returns a unique iterator pointing to the strarting knot of the
    domain. 
    
    Definition: The starting knot is the knot m_repKnots[m_deg]

    The knot-vector represents the domain which is the closed interval
    starting at this knot
    */
    uiterator domainUBegin() const
    {
        return domainSBegin().uIterator();
        // equivalent:
        //return ubegin() + domainSBegin().uIndex();
    }

    /**
    \brief Returns a unique iterator pointing to the end-knot of the domain

    Definition: The ending knot is the knot m_repKnots.end()[-m_deg-1]

    The knot-vector represents the domain closed interval with
    end-point equal to this knot
    */
    uiterator domainUEnd() const
    {
        return domainSEnd().uIterator();
        // equivalent:
        //return ubegin() + domainSEnd().uIndex();
    }

    inline bool inDomain(const T u) const
    {
        return u >= *domainBegin() && u <= *domainEnd();
        // equivalent:
        // return u >= m_repKnots[m_deg] && u <= m_repKnots.end()[-m_deg-1];
        
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
    explicit gsKnotVector(int degree)
    {
        m_deg = degree;
    }

    /// \param first starting parameter
    /// \param last end parameter parameter
    /// \param interior number of interior knots
    /// \param mult_ends multiplicity at the two end knots
    /// \param mult_interior multiplicity at the interior knots
    gsKnotVector( T first,
                  T last,
                  unsigned interior,
                  mult_t mult_ends=1,
                  mult_t mult_interior=1,
                  int degree = -1 );
    
    /// \param knots knots (sorted) including repetitions.
    /// Copies \a knots to the internal storage.
    /// To avoid the copy use gsKnotVector(give(knots)) and keep
    /// in mind that this will destroy knots.
    gsKnotVector( int degree, const knotContainer& knots );

    /// Constructs knot vector from the given degree and iterators marking its endpoints.
    /// \param begOfKnots iterator pointing to the beginning of the knots,
    /// \param endOfKnots iterator pointing past the end of the knots.
    /// The knots between \a begOfKnots and \endOfKnots are assumed to be repeated
    /// according to their multiplicities and sorted.
    template<typename iterType>
    gsKnotVector(int deg, const iterType begOfKnots, const iterType endOfKnots)
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
                      int degree=-1);

    /// Resets the knot vector so that it has \a numKnots knots
    /// between 0 and 1, each repeated \a mult_interior times (whereas
    /// 0 and 1 are repeated \a mult_ends times each) and the degree
    /// is equal to \a degree.
    void initUniform( unsigned numKnots,
                      unsigned mult_ends,
                      unsigned mult_interior = 1,
                      int degree = - 1)
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
    void set_degree(int p)
    {
        m_deg = p;
    }

    /// Writes unique indices of the knots of the endpoints of \a
    /// i - th B-spline defined on this knot vector to \a result.
    void supportIndex_into(const mult_t &i, gsMatrix<unsigned>& result) const;
    // TODO: If the function stays, make a unit test from what is in unifiedKnotVector.cpp.

    /// Returns unique knots.
    knotContainer unique () const
    {
        return knotContainer(this->ubegin(),this->uend());
    }

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

public: // findspan stuff

    /// Returns the uiterator pointing to the knot at the beginning of
    /// the interval (closed from left) containing \a u. Exception:
    /// interval *(domainUEnd()-1, *domainUEnd) is considered closed
    /// from both sides. I.e., if u == *domainUEnd(), it returns the
    /// uiterator domainUEnd() - 1.  Functionality of findElement can
    /// be obtained by calling uIndex() on the result.
    uiterator uFind( const T u ) const;

    /// Returns an iterator to the largest knot (including
    /// repetitions) that is smaller or equal to \a u and is strictly
    /// smaller than *domainEnd().  Functionality of findSpan can be
    /// easily obtained by subtracting begin() from the result.
    iterator iFind( const T u ) const;
    
    /// See uFind().
    // TODO: Remove.
    uiterator findElement( const T u ) const
    {
        return uFind( u );
    }
    
    /// Returns the index of the interval containing \a u. \sa findSpan.
    // TODO: Remove.
    unsigned findspan( T u ) const
    {
        return iFind(u) - begin();
        // equivalent
        // return findElement(u).lastAppearance();
    }

    /// Returns an iterator pointing to the beginning of the span
    /// containing the point \a u.
    // TODO Remove, use iFind() instead.
    iterator findspanIter( T u ) const
    {
        return iFind( u );
    }
    
    /// Returns the index of the "element" containing the point \a u.
    // TODO: Remove
    int findElementIndex(T u) const
    {
        return findElement(u).uIndex();
    }

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

    /// Deduces and sets the degree from the multiplicities of the
    /// endpoints.
    int deduceDegree()
    {
        return uSize() == 0 ? -1 :
            std::max(( ubegin() ).multiplicity(),
                     ( uend()-1 ).multiplicity()) - 1;
    }

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

    /// Returns the number of knot spans in the knot-vector
    int numKnotSpans() const
    {
        return uSize() - 1;
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
                       std::bind2nd(std::plus<unsigned>(), i) );
        m_multSum.back() += i;

        m_deg += i;
    }

    /// Inverse of degreeIncrease.
    void degreeDecrease(int const & i = 1 )
    {
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

      void _stretchEndKnots()
    {
        // Not necessary, I hope.
    }

    /// Get a reference to the underlying std::vector of knots.
    const knotContainer& get() const // to be removed since we have implicit cast
    {
        return m_repKnots;
    }

public: // Deprecated functions required by gsCompactKnotVector.

    /// Returns the index of the span containing the point \a u.
    unsigned Uniquefindspan (T u) const
    {
        return findElementIndex(u);
    }

    /// True iff the knot exists in the vector
    /// \param knot parameter value.
    inline bool has(T knot) const
    {
        return std::binary_search( ubegin(), uend(), knot);
        // equivalent:
        // return 0 != multiplicity(knot);
    }
     
    /// Returns the value of the \a i - th unique index
    inline T  uValue(const size_t & i) const
    { return *(ubegin()+i); }
     
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
    bool operator != (const gsKnotVector<real_t>& other) const
    {
        return ! ((*this)==other);
    }
     
    /// Returns vector of multiplicities of the knots.
    const std::vector<mult_t> multiplicities() const
    {
        std::vector<mult_t> result;
        result.reserve(uSize());
        for( uiterator uit = ubegin(); uit != uend(); ++uit )
            result.push_back( uit.multiplicity() );
        return result;
    }

    /// Returns the value of the \a i - th knot (counted with repetitions).
    inline T at (const size_t & i) const
    {
        return m_repKnots.at(i);
    }

    /// Returns the index (packed in a gsMatrix)
    gsMatrix<unsigned,1> * findspan (const gsMatrix<T,1> & u) const
    {
        // Where is this deleted then? The user is required to do so.
        gsMatrix<unsigned,1> * fs = new gsMatrix<unsigned,1>(1, u.cols() );

        for( index_t i = 0; i < u.cols(); i++ )
            (*fs)(0,i) = findspan( u(0,i) );

        return fs;
    }

    /// Returns the number of knot spans.
    unsigned spans() const
    {
        return this->uSize() - 1;
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
     
    /// Returns true iff the knot is open (ie. both endpoint
    /// multiplicities equal to degree+1)
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

    /// Returns unique knots.
    virtual knotContainer breaks() const
    {
        return knotContainer(ubegin(), uend());
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
    void degreeElevate(int const & i = 1);

    /// Converse to degreeElevate.
    void degreeReduce(int const & i);

public: // members

    // TODO remove!
    int m_deg;
};

template<typename T>
std::ostream& operator << (std::ostream& out, const gsKnotVector<T> KV )
{
    KV.print(out);
    return out;
}

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsKnotVector.hpp)
#endif
