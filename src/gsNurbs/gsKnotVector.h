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
#include <gsDomain/gsDomain.h>

#include <gsNurbs/gsKnotIterator.h>
#include <gsDomain/gsKnotDomainIterator.h>

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

# Periodic knot vectors #

A knot vector may carry a periodic flag (setPeriodic(), isPeriodic()).
It then stores the stretched knots of a periodic spline space. With
 + p the degree, [a,b] the domain and L = b-a the period,
 + the _period part_ s_0..s_{n-1}, the knots t with a <= t < b, with
   multiplicity (a appears mu times first),
 + m = p+1-mu the number of crossing functions (numCrossingFunctions()),

the stretched knots are t_0..t_{K-1}, K = n+m+p+1, with
t_k = s_{(k-m) mod n} + L floor((k-m)/n). The copies of b are stored
as b. Two invariants hold:
 + I1: the stretched knots are the ones generated from their own period part;
 + I2: 1 <= mu <= p and n >= p+1.

Read-only members (size(), the iterators, operator[], domainBegin(),
domainEnd(), greville(), uFind(), ...) keep their raw meaning on the
stretched knots. Mutating members follow one of four rules:
 + _regenerate_: the member acts on the period part and regenerates the
   ghost knots once per call, O(K). These are insert(), remove(),
   uniformRefine(), refineSpans() (domain element indices),
   degreeElevate(), degreeReduce(), degreeIncrease(), degreeDecrease(),
   increaseMultiplicity(), reduceMultiplicity(), coarsen(), set_degree(),
   reverse() (reflects the period part about a+b, so [a,b] is kept),
   affineTransformTo(), transform() and addConstant(T);
 + _refuse_: the member has no periodic meaning. These are trimLeft(),
   trimRight(), append(), erase() and addConstant(T,T) with
   start > first();
 + _reset_: the init*() members build a non-periodic vector;
 + _carry_: copy, move, swap() and operator== include the flag.

A refused call throws std::runtime_error and leaves the vector unchanged.

A knot passed to insert() has one insertion domain: a knot in (a-L,a) or
(b,b+L) is inserted as its translate in (a,b). A knot on the seam (a, b
or a congruent value), or farther than one period away, is refused.

A clamped vector given to setPeriodic() is wrapped with mu = p: a C^0
seam with one crossing function. setPeriodic() validates the seam
condition up to 100 eps max(|a|,|b|,L). Ghost knots are recomputed as
s+qL, so two mathematically equal periodic vectors built by different
routes can differ by round-off, and the exact operator== then reports
them unequal.

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

    // Unique and shared pointers
    typedef memory::unique_ptr< gsKnotVector<T> > uPtr;
    typedef memory::shared_ptr< gsKnotVector<T> > Ptr;

    typedef gsDomainIteratorWrapper<T> domainIter;

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

    // iterator

    domainIter beginAll() const override
    {
        return domainIter(new gsKnotDomainIterator<T>(*this));
    }

    domainIter beginBdr(const boxSide   /* bs */) const override
    {
        return domainIter(new gsKnotDomainIterator<T>(*this));
    }

    domainIter endAll() const override
    {
        return domainIter(new gsKnotDomainIterator<T>(*this,false));
    }

    domainIter endBdr(const boxSide   /* bs */) const override
    {
        return domainIter(new gsKnotDomainIterator<T>(*this,false));
    }

    short_t dim() const override { return 1; }

    typename gsDomain<T>::Ptr component(index_t i) const override
    {
        GISMO_UNUSED(i);
        GISMO_ASSERT(i==0, "gsKnotVector has only one component (i==0)");
        return memory::make_shared_not_owned(const_cast<gsKnotVector<T>*>(this));
    }

    gsMatrix<T> boundingBox() const override
    {
        gsMatrix<T> box(1,2);
        box << m_repKnots.front(), m_repKnots.back();
        return box;
    }

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
    ///
    /// On a periodic knot vector, \a knot is first reduced into the
    /// period: values in (a-L,a) and (b,b+L) are shifted by +-L, the
    /// period part is refined and the stretched knots are regenerated.
    /// A \a knot on the seam (a, b or a congruent value) or farther than
    /// one period away is refused, and so is a negative \a mult; the
    /// vector is then unchanged.
    /// Complexity O(K) (K = number of knots), one regeneration.
    void insert(T knot, mult_t mult = 1 );

    /// Inserts each knot between ibeg and iend with multiplicity 1.
    /// Knots can be repeated in the sequence but they must be sorted.
    /// The complexity is size() + (iend-ibeg) comparisons
    ///
    /// On a periodic knot vector, each value is reduced into the period
    /// as in insert(T,mult_t) (the first invalid value refuses the whole
    /// call), then the values are sorted and merged into the period part
    /// and the stretched knots are regenerated once.
    /// Complexity O(K + r log r) for r inserted knots.
    template<typename iterType>
    void insert( iterType ibeg, iterType iend )
    {
        if (m_periodic)
        {
            knotContainer reduced;
            reduced.reserve(iend-ibeg);
            for (iterType it = ibeg; it != iend; ++it)
                reduced.push_back( _reduceToPeriod(*it) );
            _periodicInsert( give(reduced) );
            return;
        }

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
    ///
    /// On a periodic knot vector, the period part is reduced and the
    /// stretched knots are regenerated. A ghost iterator whose first
    /// stored index is k removes its twin s[(k-m) mod n] in the period
    /// part (n periodic knots s, m = p+1-mu ghosts per side),
    /// independently of round-off. The seam and
    /// a negative \a mult are refused; the result must stay a valid
    /// periodic vector (at least p+1 periodic knots).
    /// Complexity O(K).
    void remove( uiterator uit, mult_t mult = 1 );

    /// Decreases multiplicity of the knot \a knot by \a mult.
    ///
    /// On a periodic knot vector: as remove(uiterator,mult_t), i.e. a
    /// \a knot on the seam or farther than one period away is refused;
    /// a value inside the period that is not a knot is ignored.
    /// Complexity O(K).
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

    /// Returns the maximum multiplicity in the interior
    mult_t maxInteriorMultiplicity() const;

    /// Returns the minimum multiplicity in the interior
    mult_t minInteriorMultiplicity() const;

    /// Returns the multiplicity of the knot number \a i (counter with
    /// repetitions).
    mult_t multiplicityIndex( mult_t i ) const;

public: // modifiers

    /// Scales and translates the knot vector so that its new
    /// endpoints are \a newBeg and \a newEnd.
    ///
    /// On a periodic knot vector, the period part is transformed and the
    /// stretched knots are regenerated; the period becomes
    /// \a newEnd - \a newBeg.
    /// Complexity O(K).
    void affineTransformTo( T newBeg, T newEnd );

public: // periodicity

    /// True iff this is a periodic (stretched) knot vector. O(1).
    ///
    /// A periodic knot vector of degree p stores the stretched knots
    /// t_0..t_{K-1}, K = n+m+p+1, of domain [a,b] = [t_p, t_{K-p-1}] and
    /// period L = b-a. The period part s_0..s_{n-1} are the knots t with
    /// a <= t < b (counted with multiplicity); a appears mu times first
    /// and m = p+1-mu is the number of crossing functions. The stretched
    /// knots are generated from the period part by
    /// t_k = s_{(k-m) mod n} + L floor((k-m)/n), and the copies of b are
    /// stored as b. The vector must satisfy 1 <= mu <= p and n >= p+1.
    /// Read-only members keep their raw meaning on the stretched knots.
    bool isPeriodic() const { return 0 != m_periodic; }

    /// Number m = p+1-mu of crossing functions; 0 when not periodic. O(1).
    mult_t numCrossingFunctions() const { return m_periodic; }

    /// Number of B-splines modulo periodicity:
    /// size()-degree()-1-numCrossingFunctions()
    /// (n on a periodic vector, the ordinary function count otherwise). O(1).
    index_t periodicSize() const
    { return static_cast<index_t>(size()) - m_deg - 1 - m_periodic; }

    /// Length b-a of the domain (the period L on a periodic vector). O(1).
    T period() const { return *domainEnd() - *domainBegin(); }

    /// Converts to a periodic knot vector (\a flag = true) or drops the
    /// periodic flag, keeping the stretched knots (\a flag = false).
    ///
    /// A clamped (open) vector is wrapped with seam multiplicity mu = p
    /// (one crossing function). Otherwise the stretched knots are
    /// validated as they are: 1 <= mu <= p, size() = n+m+p+1, the last
    /// copy of a at index p, equal multiplicities of a and b, and the seam
    /// condition t_{k+n} = t_k + L (k = 0..m+p) up to
    /// 100 eps max(|a|,|b|,L), eps the machine epsilon of T. A vector that fails (also n < p+1 or p < 1)
    /// is refused and left unchanged. The ghost knots are regenerated from
    /// the period part. Calling it on a periodic vector with
    /// \a flag = true does nothing.
    /// Complexity O(K).
    void setPeriodic(bool flag = true);

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
        GISMO_ASSERT(i+numLeftGhosts() >=0  && static_cast<size_t>(i+numLeftGhosts()) < uSize(),
                              "Unique index "<<i<<" not in the knot vector.");
        return *( this->ubegin()+(numLeftGhosts()+i) );
    }

    /// Provides the knot with unique index \a i
    inline T  uValue(const size_t & i) const
    { return this->operator()(i); }

    /// Number of knot intervals inside domain.
    inline size_t numElements() const override
    { return (domainUEnd() - domainUBegin()); }

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
    std::ostream &print(std::ostream &os = gsInfo ) const override;

    /// Checks whether the knot vector is in a consistent state
    bool check() const;

    typedef typename knotContainer::iterator nonConstIterator    ;
    typedef typename multContainer::iterator nonConstMultIterator;

    /// Returns an iterator pointing to the starting knot of the domain.
    iterator domainBegin() const
    {
        GISMO_ASSERT( size() > static_cast<size_t>(2*m_deg+1), "Not enough knots, expecting more than " << 2*m_deg+1);
        return begin() + m_deg;
    }

    /// Returns an iterator pointing to the end-knot of the domain.
    iterator domainEnd() const
    {
        GISMO_ASSERT( size() > static_cast<size_t>(2*m_deg+1), "Not enough knots, expecting more than " << 2*m_deg+1);
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

    /// Removes the knots in the range [first,last).
    /// Not defined on a periodic knot vector (refused).
    void erase(const mult_t first, const mult_t last);

    /// Removes the left-most \a numKnots from the knot-vector.
    /// Not defined on a periodic knot vector (refused).
    void trimLeft (const mult_t numKnots);

    /// Removes the right-most \a numKnots from the knot-vector.
    /// Not defined on a periodic knot vector (refused).
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

    /// Returns true iff the knotvector is clamped (extremal knots are of multiplicities m_deg+1)
    bool clamped() const { return 0==numLeftGhosts() && 0==numRightGhosts(); }
public:

    /// Sanity check.
    static bool isConsistent(const knotContainer & repKnots,
                             const multContainer & multSums);

    /// Compare with another knot vector (knots, degree and periodicity).
    ///
    /// Equality is exact. Two mathematically equal periodic vectors
    /// built by different routes can compare unequal, because
    /// regeneration recomputes the ghost knots as s+qL (round-off).
    inline bool operator== (const gsKnotVector<T>& other) const
    {
        return (m_repKnots == other.m_repKnots) &&
               (m_multSum  == other.m_multSum ) &&
               (m_deg == other.m_deg) &&
               (m_periodic == other.m_periodic);
    }

private: // internals

    // Creates a vector of multiplicity sums from a scratch and m_repKnots.
    void rebuildMultSum();

    // Knots t with a <= t < b (a, b = domain ends), i.e. the period part
    // of a periodic knot vector. O(K).
    knotContainer _periodPart() const;

    // The closed period of a periodic knot vector: its period part
    // followed by mu copies of b, as a clamped knot vector of degree mu-1
    // on exactly [a,b]. Having no ghosts, the non-periodic members
    // (which trim ghosts) act on [a,b] only; its unique span e is domain
    // element e of *this. The result is not periodic. O(K).
    gsKnotVector _closedPeriod() const;

    // Knots of the closed period W that lie below its last knot, i.e.
    // the new period part after W has been modified. O(K).
    static knotContainer _periodOf(const gsKnotVector & W);

    // Maps x in (a-L,a) u (a,b) u (b,b+L) into (a,b) by one shift of
    // +-L; refuses (GISMO_ERROR) everything else, in particular the seam.
    T _reduceToPeriod(T x) const;

    // Single commit point of every periodic mutation: rebuilds the
    // stretched knots of degree p and domain [a,b] from the period part
    // s (sorted, s.front()==a, all < b) by
    //   t_k = s_{(k-m) mod n} + L floor((k-m)/n),  k = 0..n+m+p,
    // with floored division. Validates 1 <= mu <= p, n >= p+1 and that no
    // ghost knot rounds onto a or b, and modifies *this only after all of
    // it succeeded. O(K).
    void _regeneratePeriodic(const knotContainer & s, short_t p, T a, T b);

    // Inserts the already reduced values (see _reduceToPeriod) into the
    // period part and regenerates once. O(K + r log r).
    void _periodicInsert(knotContainer reduced);

    // Decreases the multiplicity of the (reduced) knot in the period
    // part by mult and regenerates. O(K).
    void _periodicRemove(T knot, mult_t mult);

    // Translates the period part by amount and regenerates. O(K).
    void _periodicTranslate(T amount);

private: // members

    // Knots including repetitions.
    knotContainer m_repKnots;

    // m_multSum[i] = cardinality of { knots <= unique knot [i] }.
    multContainer m_multSum;

    // Number m of crossing functions of a periodic knot vector; 0 = not periodic.
    mult_t m_periodic = 0;

    // Degree.
    short_t m_deg;




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
        m_deg = (deg == - 1 ? deduceDegree() : deg);
    }

public:
    // TODO If stays, make it private.

    /// Resets the knot vector so that it has \a interior knots
    /// between \a first and \a last, each of them repeated \a
    /// mult_interior times (and the endpoints repeated \a mult_ends
    /// times) and the degree is equal to \a degree.
    /// The result is not periodic.
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

    /// Resets the knot vector so that its knots on [\a u0, \a u1] are
    /// graded: interior knot i (1 <= i <= \a interior) is
    /// u0 + (u1-u0) * (i/(interior+1))^(1/grading), so \a grading < 1
    /// clusters the knots towards \a u0. The result is not periodic.
    void initGraded(T u0, T u1, unsigned interior, int degree,
                    T grading, unsigned mult_interior = 1)
    {
        m_periodic = 0;
        m_repKnots.clear();
        m_multSum.clear();
        GISMO_ASSERT(u0<u1,"Knot vector must be an interval.");

        m_deg = degree;
        m_repKnots.reserve( 2*(m_deg+1) + interior*mult_interior );
        m_multSum .reserve(interior+2);

        const T h = T(1) / T(interior+1);

        m_repKnots.insert(m_repKnots.begin(), m_deg+1, u0);
        m_multSum .push_back(m_deg+1);

        for ( unsigned i=1; i<=interior; i++ )
        {
            m_repKnots.insert(m_repKnots.end(), mult_interior,
                              u0 + (u1-u0) * math::pow(T(i)*h, T(1)/grading) );
            m_multSum .push_back( mult_interior + m_multSum.back() );
        }
        m_repKnots.insert(m_repKnots.end(), m_deg+1, u1);
        m_multSum .push_back( m_deg+1 + m_multSum.back() );
    }

    /// Returns the greville points of the B-splines defined on this
    /// knot vector. Points are clamped if \a clamp is true
    gsMatrix<T> greville(bool clamp = true) const
    {
        gsMatrix<T> gr( 1,this->size() - m_deg - 1 );
        this->greville_into(gr);
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
*/



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
    ///
    /// On a periodic knot vector the domain [a,b] is kept: the period
    /// part is reflected, s -> a+b-s, and the stretched knots are
    /// regenerated. Complexity O(K).
    void reverse();

    /// @brief Insert knots into the knot vector.
    /// \param knots parameter values of the new knots, stored in a std::vector
    /// \param mult multiplicity of the new knots
    ///
    /// On a periodic knot vector, every knot is reduced into the period
    /// as in insert(T,mult_t) and the stretched knots are regenerated
    /// once (not \a mult times); a negative \a mult is refused and
    /// \a mult = 0 does nothing, as in insert(T,mult_t).
    /// Complexity O(K + r mult log(r mult)) for r knots.
    void insert( const knotContainer &knots, int mult = 1 )
    {
        if (m_periodic)
        {
            GISMO_ENSURE(mult >= 0, "insert: negative multiplicity.");
            if (0 == mult) return;
            knotContainer reduced;
            reduced.reserve(knots.size()*mult);
            for (typename knotContainer::const_iterator it = knots.begin();
                 it != knots.end(); ++it)
                reduced.insert(reduced.end(), mult, _reduceToPeriod(*it));
            _periodicInsert( give(reduced) );
            return;
        }

        for( int i = 0; i < mult; ++i)
            insert( knots.begin(), knots.end() ); // inefficient
    }

    /// Sets the degree to \a p.
    ///
    /// On a periodic knot vector, the period part is kept and the
    /// stretched knots are regenerated for degree \a p (the number of
    /// crossing functions becomes p+1-mu); a degree with mu > p or
    /// n < p+1 is refused. Complexity O(K).
    void set_degree(short_t p)
    {
        if (m_periodic)
        {
            _regeneratePeriodic(_periodPart(), p, *domainBegin(), *domainEnd());
            return;
        }
        m_deg = p;
    }

    /// Writes unique indices of the knots of the endpoints of \a
    /// i - th B-spline defined on this knot vector to \a result.
    void supportIndex_into(const mult_t &i, gsMatrix<index_t>& result) const;
    // TODO: If the function stays, make a unit test from what is in unifiedKnotVector.cpp.

    /// Inserts \a numKnots (and each of them \a mult - times) between
    /// each two knots.
    ///
    /// On a periodic knot vector, the spans of the domain [a,b] are
    /// refined (ghost spans are regenerated) and \a mult must not
    /// exceed p+1. Complexity O(K) plus the number of inserted knots,
    /// one regeneration.
    void uniformRefine( mult_t numKnots = 1, mult_t mult = 1);

    /// Inserts \a knotsPerSpan knots between each two knots between \a begin and \a end.
    ///
    /// On a periodic knot vector the indices are domain element
    /// indices (element e is the e-th span of [a,b]), each must lie in
    /// [0, numElements()), and the stretched knots are regenerated once.
    /// Complexity O(K + r k log(r k)) for r indices, k = \a knotsPerSpan.
    template<typename iterType>
    void refineSpans( iterType begin, iterType end, mult_t knotsPerSpan )
    {
        if (m_periodic)
        {
            multContainer elems;
            elems.insert(elems.begin(), begin, end);
            std::sort(elems.begin(), elems.end());
            GISMO_ENSURE( elems.empty() ||
                          (0 <= elems.front() &&
                           static_cast<size_t>(elems.back()) < numElements()),
                          "refineSpans: element index out of range.");
            gsKnotVector W = _closedPeriod();
            W.refineSpans(elems.begin(), elems.end(), knotsPerSpan);
            _regeneratePeriodic(_periodOf(W), m_deg, W.first(), W.last());
            return;
        }

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
                newKnot = ( (segmentsPerSpan-(T)(k)) * spanBegin + (T)(k) * spanEnd ) / segmentsPerSpan;
                newKnots.push_back( newKnot );
            }

        insert(newKnots.begin(), newKnots.end());
    }

    /// Adds \a amount to all the knots.
    ///
    /// On a periodic knot vector, the period part is translated and the
    /// stretched knots are regenerated. Complexity O(K).
    void addConstant( T amount );

    /// Adds \a amount to all the knots, starting at knot \a start.
    ///
    /// On a periodic knot vector only \a start <= first() (a translation
    /// of all the knots, as addConstant(T)) is defined; a larger \a start
    /// is refused. Complexity O(K).
    void addConstant( T start, T amount );

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
    short_t degree(short_t i = 0) const override //overload from gsDomain
    {
        GISMO_UNUSED(i);
        return m_deg;
    }

    /// Writes Greville abscissae of the B-splines defined on this
    /// knot vector to \a result. Points are clamped if \a clamp is true
    void greville_into(gsMatrix<T> & result, bool clamp = true) const;

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
    ///
    /// On a periodic knot vector, the period part is kept and the stretched
    /// knots are regenerated for degree p+i (crossing functions m+i); a
    /// negative \a i, or p+i >= n, is refused. Complexity O(K).
    void degreeIncrease(int const & i = 1)
    {
        if (m_periodic)
        {
            GISMO_ENSURE(i >= 0, "degreeIncrease: negative increment.");
            if (0 == i) return;
            _regeneratePeriodic(_periodPart(), static_cast<short_t>(m_deg+i),
                                *domainBegin(), *domainEnd());
            return;
        }

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
    ///
    /// On a periodic knot vector, the seam multiplicity mu is kept and
    /// the stretched knots are regenerated for degree p-i (crossing
    /// functions m-i); if \a updateInterior, interior breakpoints of the
    /// period with multiplicity above p-i are reduced to p-i. A negative
    /// \a i, or p-i < mu, is refused. Complexity O(K).
    void degreeDecrease(int const & i = 1, bool updateInterior = false)
    {
        if (m_periodic)
        {
            GISMO_ENSURE(i >= 0, "degreeDecrease: negative decrement.");
            if (0 == i) return;
            const short_t pn = static_cast<short_t>(m_deg - i);
            const T a = *domainBegin(), b = *domainEnd();
            knotContainer s = _periodPart();
            if (updateInterior)
            {
                knotContainer r;
                r.reserve(s.size());
                for (size_t j = 0; j < s.size(); )
                {
                    size_t e = j;
                    while (e < s.size() && s[e] == s[j]) ++e;
                    const size_t cnt = (s[j] == a) ? e-j
                        : (std::min)(e-j, static_cast<size_t>(pn < 0 ? 0 : pn));
                    r.insert(r.end(), cnt, s[j]);
                    j = e;
                }
                s.swap(r);
            }
            _regeneratePeriodic(s, pn, a, b);
            return;
        }

        remove( ubegin()  , i );
        remove( uend() - 1, i );
        m_deg -= i;

        if (updateInterior)
        {
            if ( m_deg <= 0 )
                initUniform(first(), last(), 0, 1, 0, 0);
            else
                for (uiterator itr = ubegin()+1; itr != uend()-1; ++itr)
                    if ( itr.multiplicity() > m_deg )
                        remove( itr, itr.multiplicity() - m_deg );
        }
    }

    /// Increase the multiplicity of all the knots by \a i. If \a
    /// boundary is set to \em false, the boundary knots are not adjusted
    ///
    /// On a periodic knot vector, \a boundary is ignored (the seam is an
    /// interior knot of the period): interior breakpoints of the period
    /// get multiplicity min(p+1, mult+i), the seam min(p, mu+i), and the
    /// stretched knots are regenerated. A negative \a i is refused.
    /// Complexity O(K).
    void increaseMultiplicity(const mult_t i = 1, bool boundary = false);

    /// Reduce the multiplicity of all knots by \a i. If \a boundary
    /// is set to \em false, the boundary knots are not adjusted
    ///
    /// On a periodic knot vector, \a boundary is ignored: interior
    /// breakpoints of the period lose \a i (removed at multiplicity <= 0),
    /// the seam gets max(1, mu-i), and the stretched knots are regenerated.
    /// A negative \a i is refused. Complexity O(K).
    void reduceMultiplicity(const mult_t i = 1, bool boundary = false);

    /// Adds the knots between \a begin and \a end to the knot vector.
    /// Not defined on a periodic knot vector (refused): use insert.
    template<typename iterType>
    void append ( iterType begin, iterType end )
    {
        GISMO_ENSURE(!isPeriodic(), "append: not defined on a periodic knot vector");
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
            if( math::abs(*uit - (*uit-(T)(1)) - df) > tol )
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
    ///
    /// On a periodic knot vector, every breakpoint of [a,b), the seam
    /// included, gains \a i and the stretched knots are regenerated
    /// (mu+i, crossing functions m unchanged). A negative \a i is refused.
    /// Complexity O(K).
    void degreeElevate(const short_t & i = 1);

    /// Converse to degreeElevate.
    ///
    /// On a periodic knot vector, interior breakpoints lose \a i (removed
    /// at multiplicity <= 0), the seam becomes max(1, mu-i), the degree
    /// becomes p-i and the stretched knots are regenerated; a result
    /// with mu > p-i is refused, as is a negative \a i. Complexity O(K).
    void degreeReduce(const short_t & i);

    /// Coarsen the knot vector: for each group of \a knotRemove
    /// consecutive interior unique knots, reduce their multiplicity
    /// by \a mul (or remove completely, if \a mul = -1), then skip
    /// \a knotSkip consecutive unique knots and go on
    ///
    /// @param knotRemove    Number of consecutive knots to remove
    /// @param knotSkip      Number of consecutive knots which stay untouched in between
    /// @param mul           Multiplicity drop for each processed knot (or -1, full remove)
    ///
    /// On a periodic knot vector the interior breakpoints of the period
    /// are processed, the stretched knots are regenerated, and the
    /// returned knots are those removed from the period part (no
    /// translates). Complexity O(K).
    std::vector<T> coarsen(index_t knotRemove = 1, index_t knotSkip = 1, mult_t  mul = -1);

};

template<typename T>
std::ostream& operator << (std::ostream& out, const gsKnotVector<T> KV )
{
    KV.print(out);
    return out;
}


/**
 * @brief      Checks if \a SubSpace is a subspace of \a Space,
 *             i.e. if the knots of \a Space are contained in \a SubSpace
 *
 * @param[in]  SubSpace  The sub space
 * @param[in]  Space     The space
 *
 * @tparam     T         Real tyoe
 *
 * @return     True if nested, False otherwise.
 */
template<typename T>
static bool isNested(const gsKnotVector<T> & SubSpace, const gsKnotVector<T> & Space)
{
    std::vector<T> difference, intersection;
    std::vector<T> knots1, knots2;
    // The unique knots of \a SubSpace contain the ones of \a Space
    std::set_intersection(  SubSpace.ubegin(), SubSpace.uend(),
                            Space.ubegin(), Space.uend(),
                            std::back_inserter(intersection) );
    // The difference of the two must not contain any knot in knots1!
    std::set_difference(    Space.ubegin(), Space.uend(),
                            intersection.begin(), intersection.end(),
                            std::back_inserter(difference) );

    return (difference.size()==0);
}

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsKnotVector
   */
  void pybind11_init_gsKnotVector(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

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
