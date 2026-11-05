/** @file gsKnotVector.hpp

    @brief Definitions of functions from gsKnotVector.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris, A. Bressan, A. Mantzaflaris
*/

#pragma once

#include <gsIO/gsBase64.h>
#include <gsIO/gsXml.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo
{

namespace internal
{

/** \brief Read a KnotVector from XML data
    \ingroup Nurbs

    The optional attribute <c>periodic="1"</c> marks a periodic knot vector
    (see gsKnotVector::setPeriodic). It is written for a periodic vector,
    together with the full stored (stretched) knot list, and applied after
    the vector has been read. With a knot list, setPeriodic(true) is
    applied to the list as read: a clamped (open) list takes its clamped
    route (seam multiplicity p); any other list must satisfy
    t_{k+n} = t_k + L up to 100 eps max(|a|,|b|,L), among setPeriodic's
    other checks. A list that setPeriodic refuses (e.g. p < 1, or fewer
    than p+1 period knots) makes reading throw std::runtime_error and
    leaves the target unchanged. With <c>mode="uniform"</c> the clamped
    uniform vector on [0,1] with <c>csize</c> coefficients is made
    periodic, which gives <c>csize-1</c> periodic functions (refused if
    csize-1 < p+1). Without the attribute the vector is not periodic.
*/
template<class T>
class gsXml< gsKnotVector<T> >
{
private:
    gsXml() { }

public:
    GSXML_COMMON_FUNCTIONS(gsKnotVector<T>)
    GSXML_GET_POINTER(gsKnotVector<T>)
    static std::string tag () { return "KnotVector"; }
    static std::string type() { return ""; }

    static void get_into(gsXmlNode * node, gsKnotVector<T> & result)
    {
        // TODO: make it unused when possible.
        int p = atoi(node->first_attribute("degree")->value() );

        gsXmlAttribute * per = node->first_attribute("periodic");
        const bool periodic = per && 0 != atoi(per->value());

        typename gsKnotVector<T>::knotContainer knotValues;

        gsXmlAttribute * mode = node->first_attribute("mode");

        // Generate Knotvectors based on key
        // mode: uniform, graded, ..
        if (mode)
        {
            if ( !strcmp( mode->value(),"uniform") )
            {
                gsXmlAttribute * szc = node->first_attribute("csize");
                GISMO_ENSURE(szc, "size of knot-vector coefficients is missing (csize attribute).");
                index_t sz = atoi(szc->value());

                //gsXmlAttribute * mlt = node->first_attribute("mult");
                //if mlt

                gsKnotVector<T> kv(0.0, 1.0, sz-p-1, p+1, 1, p);
                if (periodic) kv.setPeriodic(true);
                result.swap(kv);
                return;
            }
        }

        // Read knots from list
        gsXmlAttribute* format = node->first_attribute("format");
        std::string format_flag = format ? format->value() : "ascii";
        std::transform(format_flag.cbegin(), format_flag.cend(),
                       format_flag.begin(),
                       [](unsigned char c) { return ::tolower(c); });

        if (format_flag == "ascii") {
            // Case: mode: none/default
            std::istringstream str;
            str.str(node->value());
            for (T knot; gsGetReal(str, knot);) knotValues.push_back(knot);
        } else {
            Base64::DecodeIntoGsType(node->value(), format_flag, knotValues);
        }
        gsKnotVector<T> kv(give(knotValues), p);
        if (periodic) kv.setPeriodic(true);
        result.swap(kv);
    }

    static gsXmlNode * put (const gsKnotVector<T> & obj, gsXmlTree & data)
    {
        // Write the knot values (for now WITH multiplicities)
        std::ostringstream str;
        str << std::setprecision(REAL_DIG+1);

        for ( typename gsKnotVector<T>::iterator it = obj.begin();
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
        if (obj.isPeriodic())
            tmp->append_attribute( makeAttribute("periodic", "1", data) );

        return tmp;
    }
};

}// namespace internal


//===============//
// iterator ends //
//===============//

template<typename T>
typename gsKnotVector<T>::iterator gsKnotVector<T>::begin()  const
{
    return m_repKnots.begin();
    //return m_repKnots.data();
}

template<typename T>
typename gsKnotVector<T>::iterator gsKnotVector<T>::end()    const
{
    return m_repKnots.end();
    //return m_repKnots.data() + m_repKnots.size();
}

template<typename T>
typename gsKnotVector<T>::iterator gsKnotVector<T>::beginAt(mult_t upos)  const
{
    upos += numLeftGhosts();
    return m_repKnots.begin() + (0 == upos ? 0 : m_multSum[upos-1]);
    //return m_repKnots.data() + (0 == upos ? 0 : m_multSum[upos-1]);
}

template<typename T>
typename gsKnotVector<T>::iterator gsKnotVector<T>::endAt(mult_t upos)    const
{
    upos += numLeftGhosts();
    return m_repKnots.begin() + m_multSum[upos];
    //return m_repKnots.data() + m_multSum[upos];
}

template<typename T>
typename gsKnotVector<T>::reverse_iterator gsKnotVector<T>::rbegin()  const
{
    return std::reverse_iterator<iterator>(end());
}

template<typename T>
bool gsKnotVector<T>::includes(const gsKnotVector<T> & other) const
{
    return std::includes(begin(), end(), other.begin(), other.end() );
}

template<typename T>
void gsKnotVector<T>::difference(const gsKnotVector<T> & other,
                                std::vector<T>& result) const
{
    result.clear();
    const int sz = other.size() - size();
    result.reserve( std::abs(sz) );

    std::set_difference(begin(), end(),
                        other.begin(), other.end(),
                        std::back_inserter(result));
}

template<typename T>
void gsKnotVector<T>::symDifference(const gsKnotVector<T> & other,
                                    std::vector<T>& result) const
{
    result.clear();
    // Next line is ambiguous on MSVC (std does not overload "abs" for size_t)
    // result.reserve(std::abs(other.size()-size()));
    const int sz = other.size() - size();
    result.reserve( std::abs(sz) );

    std::set_symmetric_difference(begin(), end(),
                                  other.begin(), other.end(),
                                  std::back_inserter(result));
}

template<typename T>
gsKnotVector<T> gsKnotVector<T>::knotUnion(const gsKnotVector<T> & b) const
{
    const gsKnotVector<T> & a = *this;
    knotContainer kv;
    kv.reserve( (std::max)(a.size(),b.size()) );
    std::set_union(a.m_repKnots.begin(), a.m_repKnots.end(),
                   b.m_repKnots.begin(), b.m_repKnots.end(), std::back_inserter(kv) );

    // const T newStart = math::min(*a.domainBegin(), *b.domainBegin() );
    // const T newEnd   = math::max(*a.domainEnd()  , *b.domainEnd()   );
    return gsKnotVector<T>( give(kv), (std::max)(a.m_deg, b.m_deg) );
}

template<typename T>
gsKnotVector<T> gsKnotVector<T>::knotIntersection(const gsKnotVector<T> & b) const
{
    const gsKnotVector<T> & a = *this;
    knotContainer kv;
    kv.reserve( (std::min)(a.size(),b.size()) );
    std::set_intersection(a.m_repKnots.begin(), a.m_repKnots.end(),
                          b.m_repKnots.begin(), b.m_repKnots.begin(), std::back_inserter(kv) );
    return gsKnotVector<T>( give(kv), (std::min)(a.m_deg, b.m_deg) );
}

/*
// trim to the minimal domain such that \a dbegin and \a dend are contained
template<typename T>
const gsKnotVector<T> & gsKnotVector<T>::trimDomain(const T dbegin, const T dend) const
{
    iterator lpos = std::upper_bound(begin(), end(), dbegin) - 1; // *lpos<=dbegin
    iterator rpos = std::lower_bound(begin(), end(), dend)      ; // *rpos>=dend
    diffptr_t l = lpos  - begin() - m_deg    ;
    diffptr_t r = end() - rpos    - m_deg - 1;
    gsDebugVar(l);
    gsDebugVar(r);

    return *this;
}
*/


template<typename T>
typename gsKnotVector<T>::reverse_iterator gsKnotVector<T>::rend()    const
{
    return std::reverse_iterator<iterator>(begin());
}

template<typename T>
typename gsKnotVector<T>::uiterator gsKnotVector<T>::ubegin() const
{
    return uiterator(*this,0,numLeftGhosts());
}

template<typename T>
typename gsKnotVector<T>::uiterator gsKnotVector<T>::uend()   const
{
    return uiterator::End(*this);
}

template<typename T>
typename gsKnotVector<T>::reverse_uiterator gsKnotVector<T>::urbegin() const
{
    return reverse_uiterator(uend());
}

template<typename T>
typename gsKnotVector<T>::reverse_uiterator gsKnotVector<T>::urend()   const
{
    return reverse_uiterator(ubegin());
}

template<typename T>
typename gsKnotVector<T>::smart_iterator gsKnotVector<T>::sbegin() const
{
    return smart_iterator(*this,0,numLeftGhosts());
}

template<typename T>
typename gsKnotVector<T>::smart_iterator gsKnotVector<T>::send()   const
{
    return smart_iterator::End(*this);
}

template<typename T>
typename gsKnotVector<T>::reverse_smart_iterator gsKnotVector<T>::rsbegin() const
{
    return reverse_smart_iterator(send());
}

template<typename T>
typename gsKnotVector<T>::reverse_smart_iterator gsKnotVector<T>::rsend()   const
{
    return reverse_smart_iterator(sbegin());
}


//==============//
// constructors //
//==============//


template<typename T>
gsKnotVector<T>::gsKnotVector( knotContainer knots, short_t degree)
{
    knots.swap(m_repKnots);
    rebuildMultSum();

    m_deg = (degree == - 1 ? deduceDegree() : degree);

    GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );
}



template<typename T>
void gsKnotVector<T>::swap( gsKnotVector& other )
{
    m_repKnots.swap( other.m_repKnots );
    m_multSum.swap( other.m_multSum );
    std::swap( m_deg, other.m_deg );
    std::swap( m_periodic, other.m_periodic );

    GISMO_ASSERT(check(), "Unsorted knots or invalid multiplicities.");
}

template<typename T>
void gsKnotVector<T>::setPeriodic(bool flag)
{
    if (!flag)
    {
        m_periodic = 0;
        return;
    }
    if (m_periodic) return;

    const short_t p = m_deg;
    if (p < 1 || size() <= static_cast<size_t>(2*p+1))
        GISMO_ERROR("setPeriodic: need degree >= 1 and more than 2p+1 knots.");
    const T a = *domainBegin(), b = *domainEnd();
    if (!(a < b))
        GISMO_ERROR("setPeriodic: empty domain.");

    knotContainer s = _periodPart();
    if (isOpen())
    {
        // Clamped route: p+1 copies of a become mu = p (m = 1).
        s.erase(s.begin());
    }
    else
    {
        const index_t n  = static_cast<index_t>(s.size());
        index_t mu = 0;
        while (mu < n && s[mu] == a) ++mu;
        const index_t m = p + 1 - mu;
        if (mu < 1 || mu > p)
            GISMO_ERROR("setPeriodic: seam multiplicity must be in [1,p].");
        if (static_cast<index_t>(size()) != n + m + p + 1)
            GISMO_ERROR("setPeriodic: knot count is not n+m+p+1.");
        if (m_repKnots[p+1] == a)
            GISMO_ERROR("setPeriodic: the last copy of a must sit at index p.");
        if (multiplicity(b) != mu)
            GISMO_ERROR("setPeriodic: the seam multiplicities of a and b differ.");
        const T L   = b - a;
        const T tol = 100 * std::numeric_limits<T>::epsilon() * (std::max)( (std::max)(math::abs(a), math::abs(b)), L );
        for (index_t k = 0; k <= m + p; ++k)
            if (math::abs(m_repKnots[k+n] - (m_repKnots[k] + L)) > tol)
                GISMO_ERROR("setPeriodic: knots violate t_{k+n} = t_k + L.");
    }
    _regeneratePeriodic(s, p, a, b);
}

template<typename T>
typename gsKnotVector<T>::knotContainer gsKnotVector<T>::_periodPart() const
{
    const T a = *domainBegin(), b = *domainEnd();
    return knotContainer(std::lower_bound(m_repKnots.begin(), m_repKnots.end(), a),
                         std::lower_bound(m_repKnots.begin(), m_repKnots.end(), b));
}

template<typename T>
gsKnotVector<T> gsKnotVector<T>::_closedPeriod() const
{
    GISMO_ASSERT(m_periodic, "The closed period needs a periodic knot vector.");
    const index_t mu = m_deg + 1 - m_periodic;
    knotContainer w = _periodPart();
    w.insert(w.end(), mu, *domainEnd());
    return gsKnotVector( give(w), static_cast<short_t>(mu-1) );
}

template<typename T>
typename gsKnotVector<T>::knotContainer
gsKnotVector<T>::_periodOf(const gsKnotVector & W)
{
    return knotContainer(W.m_repKnots.begin(),
                         std::lower_bound(W.m_repKnots.begin(), W.m_repKnots.end(),
                                          W.m_repKnots.back()));
}

template<typename T>
T gsKnotVector<T>::_reduceToPeriod(T x) const
{
    const T a = *domainBegin(), b = *domainEnd();
    if (x < a)
        x += b - a;
    else if (x > b)
        x -= b - a;
    if (!(a < x && x < b))
        GISMO_ERROR("Knot " << x << " is not in the open period (" << a << ","
                    << b << ") after reduction (seam or out of range).");
    return x;
}

template<typename T>
void gsKnotVector<T>::_regeneratePeriodic(const knotContainer & s, short_t p, T a, T b)
{
    const index_t n = static_cast<index_t>(s.size());
    GISMO_ENSURE(a < b && n > 0 && s.front() == a && s.back() < b &&
                 std::is_sorted(s.begin(), s.end()),
                 "Periodic knots: the period part must be sorted, start at a and lie in [a,b).");
    index_t mu = 0;
    while (mu < n && s[mu] == a) ++mu;
    GISMO_ENSURE(1 <= mu && mu <= p,
                 "Periodic knots: seam multiplicity "<<mu<<" must be in [1,"<<p<<"].");
    GISMO_ENSURE(n >= p + 1,
                 "Periodic knots: need at least p+1 = "<<p+1<<" periodic knots, got "<<n<<".");

    const index_t m = p + 1 - mu;
    const index_t K = n + m + p + 1;
    const T L = b - a;

    // Floored division: k-m is negative for the left ghosts and C++ '/'
    // and '%' truncate toward zero. Since n > m, q is -1, 0 or 1.
    knotContainer t(K);
    for (index_t k = 0; k < K; ++k)
    {
        const index_t j = k - m;
        index_t q = j / n;
        index_t r = j - q * n;
        if (r < 0) { r += n; --q; }
        // The copies of b are stored as b itself: a+(b-a) need not
        // equal b, and the domain end must stay bit-exact.
        t[k] = (1 == q && r < mu) ? b : s[r] + T(q) * L;
    }
    GISMO_ENSURE(std::is_sorted(t.begin(), t.end()),
                 "Periodic knots: regenerated knots are not sorted (round-off).");
    // The copies of a sit at t_m..t_p and those of b at t_{n+m}..t_{n+p};
    // a ghost rounding onto a or b would change the seam multiplicity.
    GISMO_ENSURE(t[m-1] < a && b < t[n+p+1],
                 "Periodic knots: a ghost knot rounds onto the seam (round-off).");

    multContainer ms;
    for (index_t k = 1; k <= K; ++k)
        if (K == k || t[k] != t[k-1])
            ms.push_back(k);

    // All checks passed: commit.
    m_repKnots.swap(t);
    m_multSum.swap(ms);
    m_deg      = p;
    m_periodic = m;

    GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );
}

template<typename T>
void gsKnotVector<T>::_periodicInsert(knotContainer reduced)
{
    if (reduced.empty()) return;
    std::sort(reduced.begin(), reduced.end());
    gsKnotVector W = _closedPeriod();
    W.insert(reduced.begin(), reduced.end());
    for (typename knotContainer::const_iterator it = reduced.begin();
         it != reduced.end(); ++it)
        if (W.multiplicity(*it) > m_deg + 1)
            GISMO_ERROR("insert: multiplicity of knot " << *it
                        << " would exceed p+1 = " << m_deg + 1 << ".");
    _regeneratePeriodic(_periodOf(W), m_deg, W.first(), W.last());
}

template<typename T>
void gsKnotVector<T>::_periodicRemove(T knot, mult_t mult)
{
    const T x = _reduceToPeriod(knot);
    GISMO_ENSURE(mult >= 0, "remove: negative multiplicity.");
    if (0 == mult) return;
    gsKnotVector W = _closedPeriod();
    if (0 == W.multiplicity(x)) return; // not a knot, nothing to remove
    W.remove(x, mult);
    _regeneratePeriodic(_periodOf(W), m_deg, W.first(), W.last());
}

template<typename T>
void gsKnotVector<T>::_periodicTranslate(T amount)
{
    gsKnotVector W = _closedPeriod();
    W.addConstant(amount);
    _regeneratePeriodic(_periodOf(W), m_deg, W.first(), W.last());
}

template<typename T>
gsKnotVector<T>* gsKnotVector<T>::clone() const
{
    return new gsKnotVector(*this);
}

//========================//
// inserting and removing //
//========================//

// TODO: Possibly insert(uniqIter) alias multiplicity increase.
// Then refactor this insert to use it.
// Maybe not worth the effort.
template<typename T>
void gsKnotVector<T>::insert( T knot, mult_t mult )
{
    if (m_periodic)
    {
        GISMO_ENSURE(mult >= 0, "insert: negative multiplicity.");
        if (0 == mult) return;
        _periodicInsert( knotContainer(static_cast<size_t>(mult),
                                       _reduceToPeriod(knot)) );
        return;
    }

    //size_t numKnots = size()+mult;
    // GISMO_ENSURE( numKnots < std::numeric_limits<mult_t>::max(),
    //               "Too many knots." );

    uiterator uit = std::lower_bound(ubegin(), uend(), knot);
    const mult_t fa = uit.firstAppearance();

    const mult_t lg = numLeftGhosts();
    const mult_t rg = numRightGhosts();

    // update multiplicity sums
    nonConstMultIterator upos = m_multSum.begin() + uit.uCardinalIndex();
    if (upos==m_multSum.end() || *uit != knot) // knot value does not exist ?
        upos = m_multSum.insert(upos, fa );
    std::transform(upos, m_multSum.end(), upos, GS_BIND1ST(std::plus<mult_t>(), mult));

    // insert repeated knots
    m_repKnots.insert(m_repKnots.begin() + fa, mult, knot);

    if (m_deg!=-1 && size() > static_cast<size_t>(2*m_deg+1) )
    {
        if ((0 != lg) && (knot <= *domainBegin()))
            trimLeft(mult);

        if ((0 != rg) && (knot >= *domainEnd()))
            trimRight(mult);
    }

    GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );
}

template<typename T>
void gsKnotVector<T>::remove( uiterator uit, mult_t mult )
{
    GISMO_ASSERT( uit.m_mlt == multSumData(),
                  "The iterator is invalid for this knot vector." );

    if (m_periodic)
    {
        const index_t n = periodicSize();
        index_t r = (static_cast<index_t>(uit.firstAppearance()) - m_periodic) % n;
        if (r < 0) r += n;
        _periodicRemove(_periodPart()[r], mult);
        return;
    }

    mult_t knotMult = uit.multiplicity();
    mult_t toRemove = (std::min<mult_t>)(mult, knotMult);
    nonConstMultIterator upos = m_multSum.begin()  + uit.uCardinalIndex();

    nonConstIterator pos = m_repKnots.begin() + uit.firstAppearance();
    m_repKnots.erase(pos, pos+toRemove);

    if( toRemove ==  knotMult )
        upos = m_multSum.erase( upos );

    std::transform(upos, m_multSum.end(), upos, GS_BIND2ND(std::minus<mult_t>(),toRemove));
}

template<typename T>
void gsKnotVector<T>::remove( const T knot, mult_t mult )
{
    if (m_periodic)
    {
        _periodicRemove(knot, mult);
        return;
    }

    uiterator uit = std::lower_bound( ubegin(), uend(), knot );
    if( uit != uend() && *uit == knot ) // value found ?
        remove( uit, mult );
    // Otherwise the knot is not present and we cannot remove it.
}


template<typename T>
void gsKnotVector<T>::erase(const mult_t first, const mult_t last)
{
    GISMO_ENSURE(!isPeriodic(), "erase: not defined on a periodic knot vector");
    m_repKnots.erase(m_repKnots.begin()+first, m_repKnots.begin()+last);
    nonConstMultIterator fpos =
        std::lower_bound(m_multSum.begin(), m_multSum.end(), first);
    nonConstMultIterator lpos =
        std::upper_bound(m_multSum.begin(), m_multSum.end(), last);
    const mult_t numKnots = last - first;
    *fpos = m_multSum.back() - numKnots;
    lpos  = m_multSum.erase(fpos + 1, lpos);
    std::transform(lpos, m_multSum.end(), lpos, GS_BIND2ND(std::minus<mult_t>(),numKnots));
}

template<typename T>
void gsKnotVector<T>::trimLeft(const mult_t numKnots)
{
    GISMO_ENSURE(!isPeriodic(), "trimLeft: not defined on a periodic knot vector");
    // equiv:
    //erase(0, numKnots);
    //return;
    m_repKnots.erase(m_repKnots.begin(), m_repKnots.begin()+numKnots);
    nonConstMultIterator upos =
        std::upper_bound(m_multSum.begin(), m_multSum.end(), numKnots);
    upos = m_multSum.erase(m_multSum.begin(), upos);
    std::transform(upos, m_multSum.end(), upos, GS_BIND2ND(std::minus<mult_t>(),numKnots));
}

template<typename T>
void gsKnotVector<T>::trimRight(const mult_t numKnots)
{
    GISMO_ENSURE(!isPeriodic(), "trimRight: not defined on a periodic knot vector");
    // equiv:
    //erase(m_multSum.back()-numKnots, m_multSum.back());
    //return;
    m_repKnots.resize(m_repKnots.size()-numKnots);
    const mult_t newSum = m_multSum.back()-numKnots;
    nonConstMultIterator upos =
        std::lower_bound(m_multSum.begin(), m_multSum.end(), newSum) + 1;
    m_multSum.erase(upos, m_multSum.end() );
    m_multSum.back() = newSum;
}

//================//
// multiplicities //
//================//

template<typename T>
typename gsKnotVector<T>::mult_t gsKnotVector<T>::multiplicity( T u ) const
{
    uiterator uit = std::lower_bound( ubegin(), uend(), u );
    if( uit != uend() && *uit == u ) // value found ?
        return uit.multiplicity();
    return 0;
}

template<typename T>
typename gsKnotVector<T>::mult_t gsKnotVector<T>::multiplicityIndex( mult_t knotIndex ) const
{
    GISMO_ASSERT( knotIndex>=0 && static_cast<size_t>(knotIndex)<m_repKnots.size(),
                  "knotIndex " << knotIndex << "out of bounds [0,"
                  << m_repKnots.size() << ")." );

    iterator it = begin() + knotIndex;
    iterator L  = std::find_if(it,end(),GS_BIND1ST(std::not_equal_to<T>(),*it));
    reverse_iterator F = std::find_if(reverse_iterator(it),rend(),
                                      GS_BIND1ST(std::not_equal_to<T>(),*it));
    return L-F.base();
    // equivalent:
    //return (sbegin() + knotIndex).multiplicity();
}

//===========//
// modifiers //
//===========//

template<typename T>
void gsKnotVector<T>::affineTransformTo(T newBeg, T newEnd)
{
    GISMO_ASSERT(newEnd > newBeg+0.00001,
                 "Cannot transform the knot-vector to invalid interval ["<<newBeg<<","<<newEnd<<"].\n");

    if (m_periodic)
    {
        gsKnotVector W = _closedPeriod();
        W.affineTransformTo(newBeg, newEnd);
        _regeneratePeriodic(_periodOf(W), m_deg, W.first(), W.last());
        return;
    }

    const T obeg = *domainBegin();
    const T oend = *domainEnd();
    const T rr = (newEnd - newBeg) / (oend - obeg);

    for (uiterator uit = ubegin(); uit != uend(); ++uit)
        uit.setValue(newBeg + (*uit - obeg) * rr);

    domainUBegin().setValue(newBeg);
    domainUEnd().setValue(newEnd);
    GISMO_ASSERT( check(), "affineTransformTo() has produced an invalid knot vector.");
}

template<typename T>
void gsKnotVector<T>::reverse()
{
    if (m_periodic)
    {
        const T a = *domainBegin(), b = *domainEnd();
        const knotContainer s = _periodPart();
        const index_t mu = m_deg + 1 - m_periodic;
        // mu copies of a, then the interior breakpoints reflected about
        // the midpoint of [a,b]; the ends are the stored a and b since
        // (a+b)-b need not equal a in floating point.
        knotContainer r;
        r.reserve(s.size());
        r.insert(r.end(), mu, a);
        for (typename knotContainer::const_reverse_iterator it = s.rbegin();
             it != s.rend() - mu; ++it)
            r.push_back(a + b - *it);
        _regeneratePeriodic(r, m_deg, a, b);
        return;
    }

    // Not implemented using affineTransformTo() because of efficiency
    // and also to prevent accidental reversing.

    // reverse the multiplicity
    std::reverse  (m_multSum.begin(), m_multSum.end()-1);
    std::transform(m_multSum.begin(), m_multSum.end()-1, m_multSum.begin(),
                   GS_BIND1ST(std::minus<mult_t>(), m_multSum.back() ) );

    // reverse the knots
    std::reverse(m_repKnots.begin(), m_repKnots.end());
    const T ab = m_repKnots.back() + m_repKnots.front();
    for (uiterator uit = ubegin(); uit != uend(); ++uit)
        uit.setValue( ab - uit.value() );

    GISMO_ASSERT( check(), "reverse() produced an invalid knot vector.");
}


//===============//
// miscellaneous //
//===============//

template <typename T>
std::ostream & gsKnotVector<T>::print(std::ostream &os) const
{
    os << "[ " ;
    if ( size() > static_cast<size_t>(2*m_deg+8) )
    {
        for (iterator itr = begin(); itr != begin()+m_deg+3; ++itr)
            os << *itr << " ";
        os << "... ";
        for (iterator itr = end()-m_deg-3; itr != end(); ++itr)
            os << *itr << " ";
    }
    else
    {
        for (iterator itr = begin(); itr != end(); ++itr)
            os << *itr << " ";
    }
    os << "] (deg=" << degree()
       << ", size=" << size()
       << ", minSpan=" << minIntervalLength()
       << ", maxSpan=" << maxIntervalLength()
       << ")";
    return os;
}

template <class T>
T gsKnotVector<T>::maxIntervalLength() const
{
    T hmax = 0.0;
    for (uiterator it = ubegin(); it + 1 < uend(); ++it)
        hmax = math::max(hmax, (T)(*(it+1) - *it) );
    return hmax;
}

template <class T>
T gsKnotVector<T>::minIntervalLength() const
{
    T hmin = std::numeric_limits<T>::max();
    for (uiterator it = ubegin(); it + 1 < uend(); ++it)
        hmin = math::min(hmin, (T)(*(it+1) - *it) );
    return hmin;
}

template<typename T>
bool gsKnotVector<T>::check() const
{
    return isConsistent(m_repKnots,m_multSum);
}

template<typename T>
bool gsKnotVector<T>::isConsistent( const knotContainer & repKnots,
                                    const multContainer & multSum )
{
    // check size
    if (repKnots.size()==0)
    {
        if(multSum.size()==0)
            return true;
        else
            return false;
    }
    if (repKnots.size()!= static_cast<size_t>(multSum.back()) )
        return false;
    // check order and multiplicities
    T prev = repKnots.front();
    mult_t uniqPos = 0;
    for( typename knotContainer::const_iterator kit = repKnots.begin();
         kit != repKnots.end();
         ++kit )
    {
        if( *kit < prev ) // m_repKnots is locally decreasing.
            return false;
        else if( *kit > prev )
        {
            if( multSum[uniqPos] != static_cast<mult_t>(kit - repKnots.begin()) )
                return false;
            ++uniqPos;
            prev = *kit;
        }
    }
    return true;
}

template<typename T>
void gsKnotVector<T>::rebuildMultSum()
{
    m_multSum.clear();

    iterator bb=begin();
    iterator it=begin();
    iterator ee=end();

    while(it!=ee)
    {
        it = std::find_if(it,ee,GS_BIND1ST(std::not_equal_to<T>(),*it));
        m_multSum.push_back(it-bb);
    }
}

template<typename T>
gsKnotVector<T>::gsKnotVector( T first,
                               T last,
                               unsigned interior,
                               mult_t mult_ends,
                               mult_t mult_interior,
                               short_t degree)
{
    initUniform( first, last, interior, mult_ends, mult_interior, degree );
}

template<typename T>
gsKnotVector<T>::gsKnotVector( const knotContainer& uKnots,
                               int degree,
                               int regularity )
{
    // The code is very similar to that of the constructor with first, last, ints, mult, mult.
    GISMO_ASSERT( uKnots.front() < uKnots.back(),
                  "The first element in uknots has to be smaller than the last one." );

    mult_t mult_ends = degree + 1;
    mult_t mult_interior = degree - regularity;

    m_repKnots.clear();
    const size_t nKnots = uKnots.size() * mult_interior + 2 * (mult_ends-mult_interior);
    //GISMO_ENSURE( nKnots < std::numeric_limits<mult_t>::max(),
    //              "Knot vector too big." );
    m_repKnots.reserve( nKnots );
    m_multSum.clear();
    m_multSum.reserve( uKnots.size() );

    m_repKnots.insert( m_repKnots.end(), mult_ends, uKnots.front() );
    m_multSum.push_back( mult_ends );

    // We iterate from the one past begin() to one before end().
    typename knotContainer::const_iterator it = uKnots.begin();
    for( it += 1; it != uKnots.end() - 1; ++it)
    {
        m_repKnots.insert( m_repKnots.end(), mult_interior, *it );
        m_multSum.push_back( mult_interior + m_multSum.back() );
    }

    m_repKnots.insert( m_repKnots.end(), mult_ends, uKnots.back() );
    m_multSum.push_back( mult_ends + m_multSum.back() );

    //GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );

    // TODO remove
    m_deg = (degree == - 1 ? deduceDegree() : degree);
}

template <typename T>
void gsKnotVector<T>::initUniform( T first,
                                   T last,
                                   unsigned interior,
                                   unsigned mult_ends,
                                   unsigned mult_interior,
                                   short_t degree)
{
    m_periodic = 0;
    m_deg = (degree == - 1 ? mult_ends-1 : degree);

    const size_t nKnots = 2 * mult_ends + interior*mult_interior;
    // GISMO_ENSURE( nKnots < std::numeric_limits<mult_t>::max(),
    //               "Knot vector too big." );
    GISMO_ASSERT( first<last, //&& mult_ends<m_deg+2 && mult_interior<deg+2,
                  "The variable first has to be smaller than the variable last." );

    m_repKnots.clear();
    m_repKnots.reserve( nKnots );
    m_multSum .clear();
    m_multSum .reserve(interior+2);

    const T h = (last-first) / T(interior+1);

    for(unsigned i = m_deg - mult_ends + 1, j=1; i!= 0; --i, ++j)
    {   // add left ghost knots
        m_repKnots.push_back(first-(T)(i)*h);
        m_multSum .push_back(j);
    }

    m_repKnots.insert(m_repKnots.end(), mult_ends, first);
    m_multSum .push_back(mult_ends + (m_multSum.empty() ? 0 : m_multSum.back()));

    for( unsigned i=1; i<=interior; ++i)
    {
        m_repKnots.insert( m_repKnots.end(), mult_interior, first + (T)(i)*h );
        m_multSum .push_back( mult_interior + m_multSum.back() );
    }

    m_repKnots.insert( m_repKnots.end(), mult_ends, last );
    m_multSum .push_back( mult_ends + m_multSum.back() );

    for(unsigned i = 1; i!=m_deg - mult_ends + 2; ++i)
    {   // add right ghost knots
        m_repKnots.push_back(last+(T)(i)*h);
        m_multSum .push_back(m_multSum.back() + 1);
    }

    GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );
}

template<typename T>
void gsKnotVector<T>::initClamped(int degree, unsigned numKnots, unsigned mult_interior)
{
    GISMO_ASSERT( numKnots > 1 , "Not enough knots.");
    initUniform(0.0, 1.0, numKnots - 2, degree + 1, mult_interior, degree );
}

template<typename T>
int gsKnotVector<T>::deduceDegree() const
{
    return uSize() == 0 ? -1 :
        (std::max)(( ubegin() ).multiplicity(),
                   ( uend()-1 ).multiplicity()) - 1;
}

template<typename T>
typename gsKnotVector<T>::multContainer gsKnotVector<T>::multiplicities() const
{
    multContainer result;
    result.reserve(uSize());
    for( uiterator uit = ubegin(); uit != uend(); ++uit )
        result.push_back( uit.multiplicity() );
    return result;
}

template<typename T>
typename gsKnotVector<T>::mult_t gsKnotVector<T>::maxInteriorMultiplicity() const
{
    mult_t result = 0;
    for( uiterator uit = std::next(ubegin()); uit != std::prev(uend()); ++uit )
        result = std::max(result,uit.multiplicity());
    return result;
}

template<typename T>
typename gsKnotVector<T>::mult_t gsKnotVector<T>::minInteriorMultiplicity() const
{
    // No interior knots
    if (uSize()==2)
        return 0;

    mult_t result = m_deg+2;
    for( uiterator uit = std::next(ubegin()); uit != std::prev(uend()); ++uit )
        result = std::min(result,uit.multiplicity());
    return result;
}

template<typename T>
typename gsKnotVector<T>::uiterator
gsKnotVector<T>::uFind( const T u ) const
{
    GISMO_ASSERT(size()>1,"Not enough knots."); // todo: check() --> size() > 2*m_deg+1
    GISMO_ASSERT(inDomain(u), "Point "<< u <<" outside active area of the knot vector");

    // The last element is closed from both sides.
    uiterator dend = domainUEnd();

    if (u==*dend) // knot at domain end ?
        return --dend;
    else
        return std::upper_bound( domainUBegin(), dend, u ) - 1;
}

template<typename T>
typename gsKnotVector<T>::uiterator
gsKnotVector<T>::uUpperBound( const T u ) const
{
    GISMO_ASSERT(inDomain(u), "Point outside active area of the knot vector");
    return std::upper_bound( ubegin(), uend(), u );
}

template<typename T>
typename gsKnotVector<T>::uiterator
gsKnotVector<T>::uLowerBound( const T u ) const
{
    GISMO_ASSERT(inDomain(u), "Point outside active area of the knot vector");
    return std::lower_bound( ubegin(), uend(), u );
}

template<typename T>
typename gsKnotVector<T>::iterator
gsKnotVector<T>::iFind( const T u ) const
{
    // GISMO_ASSERT done in uFind().
    return begin() + uFind(u).lastAppearance();

    // equivalent
    /*GISMO_ASSERT(inDomain(u), "Point outside active area of the knot vector");
      iterator dend = domainEnd();
      if ( u == *dend )
      return --dend;
      else
      return std::upper_bound( domainBegin(), dend, u) - 1; */
}

template<typename T>
void gsKnotVector<T>::refineSpans( const multContainer & spanIndices, mult_t knotsPerSpan )
{
    refineSpans(spanIndices.begin(), spanIndices.end(), knotsPerSpan);
}

template <typename T>
std::string gsKnotVector<T>::detail() const
{
    std::stringstream osk, osm;
    osk << "[ " ;
    osm << "( " ;
    for (uiterator itr = ubegin(); itr != uend(); ++itr)
    {
        osk << itr.value()        << " ";
        osm << itr.multiplicity() << " ";
    }
    osm <<")\n";
    osk << "] ~ ";
    osm << "deg="<<m_deg<< ", size="<<size() << ", uSize="<< uSize()
        <<", minSpan="<< minIntervalLength()
        <<", maxSpan="<< maxIntervalLength() <<"\n";

    return osk.str() + osm.str();
}

template<typename T>
void gsKnotVector<T>::uniformRefine( mult_t numKnots, mult_t mult )
{
    //GISMO_ASSERT( numKnots>=0, "Expecting non-negative number");
    if( numKnots <= 0 ) return;

    if (m_periodic)
    {
        GISMO_ENSURE(0 <= mult && mult <= m_deg+1,
                     "uniformRefine: multiplicity out of range [0,p+1].");
        if (0 == mult) return;
        gsKnotVector W = _closedPeriod();
        W.uniformRefine(numKnots, mult);
        _regeneratePeriodic(_periodOf(W), m_deg, W.first(), W.last());
        return;
    }

    const mult_t l = ( domainUBegin() - ubegin() ) * numKnots * mult;
    const mult_t r = (uend() - domainUEnd() - 1  ) * numKnots * mult;

    knotContainer newKnots;
    getUniformRefinementKnots(numKnots,newKnots,mult); // iterate instead of store ?
    insert(newKnots.begin(),newKnots.end()); // complexity: newKnots.size() + size()

    if (0!=l) trimLeft (l);
    if (0!=r) trimRight(r);
}


// TODO use affineTransformTo.
template<typename T>
void gsKnotVector<T>::addConstant( T amount )
{
    if (m_periodic)
    {
        _periodicTranslate(amount);
        return;
    }

    // TODO add test
    // With C++11, you can use:
    // std::for_each( m_repKnots.begin(), m_repKnots.end(),  [amount](T& k){ k+= amount;} );

    std::transform( m_repKnots.begin(), m_repKnots.end(), m_repKnots.begin(),
                    GS_BIND1ST(std::plus<T>(),amount) );
}

template<typename T>
void gsKnotVector<T>::addConstant(T start, T amount)
{
    if (m_periodic)
    {
        GISMO_ENSURE(start <= first(),
                     "addConstant(start,amount) with start > first(): "
                     "not defined on a periodic knot vector");
        _periodicTranslate(amount);
        return;
    }

    typename std::vector<T>::iterator beg   = m_repKnots.begin();
    while (*beg < start)
        beg++;

    // The last element is closed from both sides.
    uiterator dend = domainUEnd();

    // If statement is needed because uFind gives back an iterator to the n-1 entry while the nth is required
    if (start!=*dend)
    {
        uiterator uit = uFind(start);
        for (; uit != uend(); ++uit)
            uit.setValue(*uit + amount);
    }
    else
        dend.setValue(*dend + amount);


    GISMO_ASSERT( check(), "addConstant() has produced an invalid knot vector.");
}

template<typename T>
void gsKnotVector<T>::increaseMultiplicity(const mult_t i, bool boundary)
{
    if (m_periodic)
    {
        GISMO_ENSURE(i >= 0, "increaseMultiplicity: negative increment.");
        const T a = *domainBegin(), b = *domainEnd();
        const mult_t p = m_deg;
        const knotContainer s = _periodPart();
        knotContainer r;
        r.reserve(s.size());
        for (size_t j = 0; j < s.size(); )
        {
            size_t e = j;
            while (e < s.size() && s[e] == s[j]) ++e;
            const mult_t cap = (s[j] == a) ? p : p+1; // seam: mu <= p
            r.insert(r.end(), (std::min<mult_t>)(cap, static_cast<mult_t>(e-j)+i), s[j]);
            j = e;
        }
        _regeneratePeriodic(r, m_deg, a, b);
        return;
    }

    GISMO_ASSERT( i>=0, "Expecting non-negative number");
    const mult_t lg = (domainUBegin() - ubegin()) * i;
    const mult_t rg = (uend() - domainUEnd() - 1) * i;
    size_t newSize = size() + i*(uSize()-2);
    knotContainer tmp;
    tmp.reserve(newSize);

    // fixme: should treat all boundary knots
    tmp.insert(tmp.end(), m_multSum.front() + (boundary ? i : 0),
               m_repKnots.front());

    // for all interior knots
    uiterator uit = ubegin()+1;
    for (; uit != uend()-1; ++uit)
        tmp.insert(tmp.end(), i + uit.multiplicity(), *uit);

    // last knot
    tmp.insert(tmp.end(), uit.multiplicity() + (boundary ? i : 0) , *uit);

    m_repKnots.swap(tmp);

    // update multiplicity sums
    mult_t r = ( boundary ?  1 : 0 );
    for (nonConstMultIterator m = m_multSum.begin(); m != m_multSum.end()-1; ++m)
        *m += i * r++;
    m_multSum.back() += i * (boundary ? r : r-1 );

    if (0 != lg) trimLeft(lg);
    if (0 != rg) trimRight(rg);
}

template<typename T>
void gsKnotVector<T>::reduceMultiplicity(const mult_t i, bool boundary)
{
    if (m_periodic)
    {
        GISMO_ENSURE(i >= 0, "reduceMultiplicity: negative decrement.");
        gsKnotVector W = _closedPeriod();
        W.reduceMultiplicity(i, true);
        _regeneratePeriodic(_periodOf(W), m_deg, W.first(), W.last());
        return;
    }

    GISMO_ASSERT( i>=0, "Expecting non-negative number");
    knotContainer ktmp;
    multContainer mtmp;
    ktmp.reserve(size());
    ktmp.reserve(uSize());

    // fixme: should treat all boundary knots
    mult_t bm = boundary ? (std::max<mult_t>)(1,m_multSum.front()-i) : m_multSum.front();
    mtmp.push_back(bm); // first knot
    ktmp.insert(ktmp.end(), bm, m_repKnots.front());

    uiterator uit = ubegin() + 1;
    for (; uit != uend()-1; ++uit) // for all interior knots
    {
        const mult_t m = uit.multiplicity() - i;
        if ( m > 0 )
        {
            mtmp.push_back( m + mtmp.back() );
            ktmp.insert(ktmp.end(), m, *uit);
        }
    }

    // last knot
    bm = boundary ? (std::max<mult_t>)(1,uit.multiplicity()-i) : uit.multiplicity();
    mtmp.push_back( bm + mtmp.back() );
    ktmp.insert(ktmp.end(), bm, *uit);

    m_multSum .swap(mtmp);
    m_repKnots.swap(ktmp);
}

template<typename T>
void gsKnotVector<T>::degreeElevate(const short_t & i)
{
    if (m_periodic)
    {
        GISMO_ENSURE(i >= 0, "degreeElevate: negative increment.");
        if (0 == i) return;
        // W.increaseMultiplicity (not this->, which caps) raises the
        // multiplicity of every breakpoint of [a,b), the seam included.
        gsKnotVector W = _closedPeriod();
        W.increaseMultiplicity(i, true);
        _regeneratePeriodic(_periodOf(W), static_cast<short_t>(m_deg+i),
                            W.first(), W.last());
        return;
    }

    increaseMultiplicity(i,true);
    m_deg += i;
}

template<typename T>
void gsKnotVector<T>::degreeReduce(const short_t & i)
{
    if (m_periodic)
    {
        GISMO_ENSURE(i >= 0, "degreeReduce: negative decrement.");
        if (0 == i) return;
        gsKnotVector W = _closedPeriod();
        W.reduceMultiplicity(i, true);
        _regeneratePeriodic(_periodOf(W), static_cast<short_t>(m_deg-i),
                            W.first(), W.last());
        return;
    }

    reduceMultiplicity(i,true);
    m_deg -= i;
}

template <typename T>
std::vector<T> gsKnotVector<T>::
coarsen(index_t knotRemove, index_t knotSkip, mult_t mul)
{
    GISMO_ASSERT(knotRemove>=0 && knotSkip>0, "Invalid parameters to knot-coarsening.");

    if (m_periodic)
    {
        if (0 == knotRemove) return std::vector<T>();
        // The default -1 means "remove completely", i.e. p+1 of the
        // periodic vector; on the closed period W it would be mu.
        if (-1 == mul) mul = m_deg + 1;
        gsKnotVector W = _closedPeriod();
        std::vector<T> removed = W.coarsen(knotRemove, knotSkip, mul);
        _regeneratePeriodic(_periodOf(W), m_deg, W.first(), W.last());
        return removed;
    }

    // Special value -1
    if (-1==mul) mul = m_deg+1;

    std::vector<T> coarseKnots, removedKnots;
    if (0==knotRemove) return removedKnots;

    // knots to be removed
    removedKnots.reserve( knotRemove * this->uSize() / (knotRemove+knotSkip) );

    uiterator it   = domainUBegin() + 1;
    uiterator last = domainUEnd();

    for(; it<last; it += knotSkip)
        for(index_t c = 0; it!=last && c!=knotRemove; ++c, ++it)
            removedKnots.insert(
                removedKnots.end(),
                (std::min<mult_t>)( mul, it.multiplicity() ),
                it.value()
            );

    // copy non-removed knots into coarseKnots
    coarseKnots.reserve(m_repKnots.size()-removedKnots.size());
    std::set_difference( m_repKnots.begin(), m_repKnots.end(),
                         removedKnots.begin(), removedKnots.end(),
                         std::back_inserter(coarseKnots) );

    coarseKnots.swap(m_repKnots);
    rebuildMultSum();

    GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );

    return removedKnots;
}

template<typename T>
void gsKnotVector<T>::greville_into(gsMatrix<T> & result, bool clamp) const
{
    iterator itr = begin() + 1;
    result.resize(1, size()-m_deg-1 ) ;
    unsigned i = 1;
    const mult_t lg = numLeftGhosts();
    const mult_t rg = numRightGhosts();
    
    if (clamp && ((lg != 0) || (rg != 0)))
    {
        // The ghost-clamping inserts hit the seam, which a periodic
        // vector refuses; the copy is a plain stretched knot vector.
        gsKnotVector tmp(*this);
        tmp.m_periodic = 0;
        tmp.insert(*domainBegin(), lg);
        tmp.insert(*domainEnd(),rg);
        tmp.greville_into(result);
        return;
    }
    
    if ( 0 != m_deg )
    {   
        result(0,0) = std::accumulate( itr, itr+m_deg, (T)(0.0) ) / (T)(m_deg) ;

        if ( result(0,0) < *(itr-1) )// Ensure that the point is in range
            result(0,0) = *(itr-1);

        for (++itr; itr != end()-m_deg; ++itr, ++i )
        {
            result(0,i) = std::accumulate( itr, itr+m_deg, (T)(0.0) ) / (T)(m_deg) ;
            if ( result(0,i) == result(0,i-1) )
            {
                // perturbe point to remain inside the needed support
                result(0,i-1) -= (T)(1)/(T)(1000000000); // =mpq_class(1,1000000000);
                //to try: result(0,i-1) = math::nextafter(result(0,i-1), *result.data() );
            }
        }

        itr = end()-1;
        i -= 1;
        if ( result(0,i) > *(itr) )// Ensure that the point is in range
            result(0,i) = *(itr);
    }
    else
        copy_n(m_repKnots.data(), result.size(), result.data() );
}

template<typename T>
void gsKnotVector<T>::centers_into(gsMatrix<T> & result) const
{
    result.resize(1, numElements());
    index_t i = 0;
    for (uiterator it = domainUBegin(); it != domainUEnd(); ++it, ++i)
        result.at(i) = (*(it+1) + *it) / (T)(2);
}

template <class T>
T gsKnotVector<T>::greville(int i) const
{
    GISMO_ASSERT(i>=0 && static_cast<size_t>(i) < size()-m_deg-1,
                 "Index of Greville point is out of range.");
    iterator itr = begin() + 1;
    return ( m_deg==0 ? *(itr+i-1) :
             (T)(std::accumulate( itr+i, itr+i+m_deg, (T)(0.0) ) / (T)(m_deg)
             // Special case C^{-1}
             - (*(itr+i) == *(itr+i+m_deg) ? 1e-10 : 0 ))
        );
}

template <class T>
void gsKnotVector<T>::getUniformRefinementKnots(mult_t knotsPerSpan, knotContainer& result, mult_t mult) const
{
    // forward iterator ?

    result.clear();
    result.reserve( (uSize()-1) * knotsPerSpan * mult );

    T prev = m_repKnots.front();
    for( uiterator uit = ubegin()+1; uit != uend(); ++uit )
    {
        const T step = (*uit - prev)/(T)(knotsPerSpan+1);
        for( mult_t i = 1; i <= knotsPerSpan; ++ i)
            result.insert( result.end(), mult, prev + (T)(i)*step );
        prev=*uit;
    }
}

template< typename T>
void gsKnotVector<T>::supportIndex_into(const mult_t& i,
                                        gsMatrix<index_t>& result) const
{
    result.resize(1,2);
    smart_iterator it = sbegin() + i;
    result.at(0) = it.uIndex();
    it += m_deg+1;
    result.at(1) = it.uIndex();
}

} // namespace gismo
