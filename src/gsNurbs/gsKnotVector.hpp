/** @file gsKnotVector.hpp

    @brief Definitions of functions from gsKnotVector.h

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris, A. Bressan, A. Mantzaflaris
*/

#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsXml.h>

#include <sstream>
#include <numeric>
#include <limits>

#pragma once
namespace gismo
{

namespace internal
{


/// Get a KnotVector from XML data
///
/// \ingroup Nurbs
template<class T>
class gsXml< gsKnotVector<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS( gsKnotVector<T>)
    static std::string tag () { return "KnotVector"; }
    static std::string type() { return ""; }

    static gsKnotVector<T> * get (gsXmlNode * node)
    {
        // TODO: make it unused when possible.
        int p = atoi(node->first_attribute("degree")->value() );
        //GISMO_UNUSED(p);

        typename gsKnotVector<T>::knotContainer knotValues;

        std::istringstream str;
        str.str( node->value() );
        for (T knot; str >> knot;)
            knotValues.push_back(knot);

        gsKnotVector<T> * kv = new gsKnotVector<T>(p,knotValues);

        return kv;
    }

    static gsXmlNode * put (const gsKnotVector<T> & obj, gsXmlTree & data)
    {
        // Write the knot values (for now WITH multiplicities)
        std::ostringstream str;
        str << std::setprecision(std::numeric_limits<T>::digits10+1);

        for ( typename gsKnotVector<T>::iterator it = obj.begin();
              it != obj.end(); ++it )
        {
            str << *it <<" ";
        }

        // Make a new XML KnotVector node
        gsXmlNode * tmp = internal::makeNode("KnotVector", str.str(), data);
        // Append the degree attribure
        str.str(std::string());// clean the ostream
        str<< obj.m_deg;
        tmp->append_attribute( makeAttribute("degree", str.str(),data) );

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
}

template<typename T>
typename gsKnotVector<T>::iterator gsKnotVector<T>::end()    const
{
    return m_repKnots.end();
}

template<typename T>
typename gsKnotVector<T>::reverse_iterator gsKnotVector<T>::rbegin()  const
{
    return reverse_iterator(end());
}

template<typename T>
typename gsKnotVector<T>::reverse_iterator gsKnotVector<T>::rend()    const
{
    return reverse_iterator(begin());
}



template<typename T>
typename gsKnotVector<T>::uiterator gsKnotVector<T>::ubegin() const
{
    //return uiterator( m_repKnots.begin(), m_multSum.begin(), m_multSum.begin() );
    return uiterator(*this);
}

template<typename T>
typename gsKnotVector<T>::uiterator gsKnotVector<T>::uend()   const
{
    //return uiterator( m_repKnots.begin(), m_multSum.begin(), m_multSum.end() );
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
    return smart_iterator(*this);
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
gsKnotVector<T>::gsKnotVector( gsMovable< knotContainer > knots )
{
    m_repKnots.swap( knots.ref() );
    rebuildMultSum();
    deduceDegree();
    GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );
}

template<typename T>
void gsKnotVector<T>::swap( gsKnotVector& other )
{
    m_repKnots.swap( other.m_repKnots );
    m_multSum.swap( other.m_multSum );
    std::swap( m_deg, other.m_deg );

    GISMO_ASSERT(check(), "Unsorted knots or invalid multiplicities.");
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
    //size_t numKnots = size()+mult;
    // GISMO_ENSURE( numKnots < std::numeric_limits<mult_t>::max(),
    //               "Too many knots." );

    uiterator uit = std::lower_bound(ubegin(), uend(), knot);
    const mult_t fa = uit.firstAppearance();

    // update multiplicity sums
    fwMIt upos = m_multSum.begin() + uit.uIndex();
    if (upos==m_multSum.end() || *uit != knot) // knot value does not exist ?
        upos = m_multSum.insert(upos, fa );
    std::transform(upos, m_multSum.end(), upos, std::bind1st(std::plus<mult_t>(), mult));

    // insert repeated knots
    m_repKnots.insert(m_repKnots.begin() + fa, mult, knot);

    GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );
}

template<typename T>
void gsKnotVector<T>::remove( uiterator uit, mult_t mult )
{
    GISMO_ASSERT( uit.m_mlt == multSumData(), 
                  "The iterator is invalid for this knot vector." );

    mult_t knotMult = uit.multiplicity();
    mult_t toRemove = std::min(mult, knotMult);
    fwMIt upos = m_multSum.begin()  + uit.uIndex();

    fwKIt pos  = m_repKnots.begin() + uit.firstAppearance(); 
    m_repKnots.erase(pos, pos+toRemove);

    if( toRemove ==  knotMult )
        upos = m_multSum.erase( upos );

    std::transform(upos, m_multSum.end(), upos, std::bind2nd(std::minus<mult_t>(),toRemove));
}

template<typename T>
void gsKnotVector<T>::remove( const T knot, mult_t mult )
{
    uiterator uit = std::lower_bound( ubegin(), uend(), knot );
    if( *uit == knot )
        remove( uit, mult );
    // Otherwise the knot is not present and we cannot remove it.
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

    cFwKIt it = begin() + knotIndex;
    cFwKIt L  = std::find_if(it,m_repKnots.end(),std::bind1st(std::not_equal_to<T>(),*it));
    cBwKIt F  = std::find_if(cBwKIt(it),m_repKnots.rend(),std::bind1st(std::not_equal_to<T>(),*it));
    return L-F.base();
}

//===========//
// modifiers //
//===========//

template<typename T>
void gsKnotVector<T>::affineTransformTo(T newBeg, T newEnd)
{
    GISMO_ASSERT(newEnd > newBeg+0.00001, 
    "Cannot transform the knot-vector to invalid interval ["<<newBeg<<","<<newEnd<<"].\n");

    const T beg   = m_repKnots.front();
    const T rr    = (newEnd - newBeg) / (m_repKnots.back() - beg);
    uiterator uit = ubegin();
    uit.setValue(newBeg);
    ++uit;
    for (; uit != uend()-1; ++uit)
        uit.setValue(newBeg + (*uit - beg) * rr);
    uit.setValue(newEnd);

    GISMO_ASSERT( check(), "affineTransformTo() has produced an invalid knot vector.");
}

template<typename T>
void gsKnotVector<T>::reverse()
{
    // Not implemented using affineTransformTo() because of efficiency
    // and also to prevent accidental reversing.
    
    // reverse the multiplicity
    std::reverse  (m_multSum.begin(), m_multSum.end()-1);
    std::transform(m_multSum.begin(), m_multSum.end()-1, m_multSum.begin(), 
                   std::bind1st(std::minus<T>(), m_multSum.back() ) );

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
    if( size() < 2 ) return 0.0;

    T hmax = 0.0;
    for (uiterator it = ubegin(); it < uend()-2; ++it)
        hmax = math::max(hmax, *(it+1) - *it );
    return hmax;
}

template <class T>
T gsKnotVector<T>::minIntervalLength() const
{
    if( size() < 2 ) return 0.0;

    T hmin = std::numeric_limits<T>::max();
    for (uiterator it = ubegin(); it < uend()-2; ++it)
        hmin = math::min(hmin, *(it+1) - *it );
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
    for( cFwKIt kit = repKnots.begin();
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
        it = std::find_if(it,ee,std::bind1st(std::not_equal_to<T>(),*it));
        m_multSum.push_back(it-bb);
    }
}

template<typename T>
gsKnotVector<T>::gsKnotVector( T first,
                                             T last,
                                             unsigned interior,
                                             mult_t mult_ends,
                                             mult_t mult_interior,
                                             int degree)
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
    cFwKIt it = uKnots.begin();
    for( it += 1; it != uKnots.end() - 1; ++it)
    {
        m_repKnots.insert( m_repKnots.end(), mult_interior, *it );
        m_multSum.push_back( mult_interior + m_multSum.back() );
    }

    m_repKnots.insert( m_repKnots.end(), mult_ends, uKnots.back() );
    m_multSum.push_back( mult_ends + m_multSum.back() );

    //GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );

    // TODO remove
    if( degree == - 1 )
        deduceDegree();
    else
        m_deg = degree;
}

template<typename T>
gsKnotVector<T>::gsKnotVector(int degree, const knotContainer& knots )
{
    m_repKnots = knots;
    m_deg = degree;
    rebuildMultSum();
    GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );
}

template <typename T>
void gsKnotVector<T>::initUniform( T first,
                                          T last,
                                          unsigned interior,
                                          unsigned mult_ends,
                                          unsigned mult_interior,
                                          int degree)
{
    const size_t nKnots = 2 * mult_ends + interior*mult_interior;
    // GISMO_ENSURE( nKnots < std::numeric_limits<mult_t>::max(),
    //               "Knot vector too big." );
    GISMO_ASSERT( first<last,
                  "The variable first has to be smaller than the variable last." );

    m_repKnots.clear();
    m_repKnots.reserve( nKnots );
    m_multSum .clear();
    m_multSum .reserve(interior+2);

    const T h = (last-first) / (interior+1);

    m_repKnots.insert(m_repKnots.begin(), mult_ends, first);
    m_multSum .push_back(mult_ends);

    for( unsigned i=1; i<=interior; i++ )
    {
        m_repKnots.insert( m_repKnots.end(), mult_interior, first + i*h );
        m_multSum .push_back( mult_interior + m_multSum.back() );
    }

    m_repKnots.insert( m_repKnots.end(), mult_ends, last );
    m_multSum .push_back( mult_ends + m_multSum.back() );

    if( degree == - 1 )
        deduceDegree();
    else
        m_deg = degree;
        
    GISMO_ASSERT( check(), "Unsorted knots or invalid multiplicities." );
}

template<typename T>
void gsKnotVector<T>::initClamped(int degree, unsigned numKnots, unsigned mult_interior)
{
    GISMO_ASSERT( numKnots > 1 , "Not enough knots.");
    initUniform(0.0, 1.0, numKnots - 2, degree + 1, mult_interior, degree );
}

template<typename T>
int gsKnotVector<T>::degree() const
{
    return m_deg;
}

template<typename T>
typename gsKnotVector<T>::uiterator
gsKnotVector<T>::uFind( const T u ) const
{
    GISMO_ASSERT(size()>1,"Not enough knots."); // todo: check() --> size() > 2*m_deg+1
    GISMO_ASSERT(inDomain(u), "Point outside active area of the knot vector");

    // The last element is closed from both sides.
    uiterator dend = domainUEnd();

    if (u==*dend) // knot at domain end ? 
        return --dend;
    else 
        return std::upper_bound( domainUBegin(), dend, u ) - 1;
}

template<typename T>
typename gsKnotVector<T>::iterator 
gsKnotVector<T>::iFind( const T u ) const
{ 
    // GISMO_ASSERT done in uFind().
    return  m_repKnots.begin() + uFind(u).lastAppearance();

    // equivalentg
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
    osk << "]\n";
    osm << "deg="<<m_deg<< ", size="<<size() << ", uSize="<< uSize()
        <<", minSpan="<< minIntervalLength() 
        <<", maxSpan="<< maxIntervalLength() <<"\n";

    return osk.str() + osm.str();
}

template<typename T>
void gsKnotVector<T>::uniformRefine( mult_t numKnots, mult_t mult )
{
    //GISMO_ASSERT( numKnots>=0, "Expecting non-negative number");
    if( numKnots < 0 )
        return;
        
    knotContainer newKnots;
    getUniformRefinementKnots(numKnots,newKnots,mult); // iterate instead of store ?
    insert(newKnots.begin(),newKnots.end()); // complexity: newKnots.size() + size()
    // optimal: 
    // newKnots, reserve

    // fixme: trim extra knots
}


// TODO use affineTransformTo.
template<typename T>
void gsKnotVector<T>::addConstant( T amount )
{
    // TODO add test
    // With C++11, you can use:
    // std::for_each( m_repKnots.begin(), m_repKnots.end(),  [amount](T& k){ k+= amount;} );

    std::transform( m_repKnots.begin(), m_repKnots.end(), m_repKnots.begin(), 
                    std::bind1st(std::plus<T>(),amount) );
}

template<typename T>
void gsKnotVector<T>::increaseMultiplicity(const mult_t i, bool boundary)
{
    GISMO_ASSERT( i>=0, "Expecting non-negative number");
    size_t newSize = size() + i*(uSize()-2);
    GISMO_ASSERT( newSize>=0, "Invalid input to adjustMultiplicity");
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
    for (fwMIt m = m_multSum.begin(); m != m_multSum.end()-1; ++m)
        *m += i * r++;
    m_multSum.back() += i * (boundary ? r : r-1 );
}

template<typename T>
void gsKnotVector<T>::reduceMultiplicity(const mult_t i, bool boundary)
{
    GISMO_ASSERT( i>=0, "Expecting non-negative number");
    knotContainer ktmp;
    multContainer mtmp;
    ktmp.reserve(size());
    ktmp.reserve(uSize());

    // fixme: should treat all boundary knots
    mult_t bm = boundary ? math::max(1,m_multSum.front()-i) : m_multSum.front();
    mtmp.push_back(bm); // first knot
    ktmp.insert(ktmp.end(), bm, m_repKnots.front());

    uiterator uit = ubegin() + 1;
    for (; uit != uend()-1; ++uit) // for all interior knots
    {
        const mult_t m = uit.multiplicity();
        if ( m > i )
        {
            mtmp.push_back( m + mtmp.back() );
            ktmp.insert(ktmp.end(), m - i, *uit);
        }
    }
    // last knot
    bm = boundary ? math::max(1,uit.multiplicity()-i) : uit.multiplicity();
    mtmp.push_back( bm + mtmp.back() );
    ktmp.insert(ktmp.end(), bm, *uit);

    m_multSum .swap(mtmp);
    m_repKnots.swap(ktmp);
}

template<typename T>
void gsKnotVector<T>::degreeElevate(int const & i)
{
    increaseMultiplicity(i,true);
    m_deg += i;
}

template<typename T>
void gsKnotVector<T>::degreeReduce(int const & i)
{
    reduceMultiplicity(i,true);
    m_deg -= i;
}

template<typename T>
void gsKnotVector<T>::greville_into(gsMatrix<T> & result) const
{
    // Copy pasted from gsKnotVector.hpp and updated to our needs.
    iterator itr = this->begin() + 1;
    const int p = m_deg;
    result.resize(1, this->size() -p-1 ) ;
    unsigned i = 1;

    if ( m_deg!=0)
    {
        result(0,0) = std::accumulate( itr, itr+p, T(0.0) ) / m_deg ;

        if ( result(0,0) < *(itr-1) )// Ensure that the point is in range
            result(0,0) = *(itr-1);

        for (++itr; itr != this->end()-p; ++itr, ++i )
        {
            result(0,i) = std::accumulate( itr, itr+p, T(0.0) ) / m_deg ;
            if ( result(0,i) == result(0,i-1) )
                // perturbe point to remain inside the needed support
                result(0,i-1) -= 1e-10;
        }

        itr = this->end()-1;
        i -= 1;
        if ( result(0,i) > *(itr) )// Ensure that the point is in range
            result(0,i) = *(itr);
    }
    else
        std::copy(this->begin(), this->end()-1, result.data() );
}

template <class T>
T gsKnotVector<T>::greville(int i) const
{
    // Copy pasted from gsKnotVector::greville(int).
    GISMO_ASSERT(i>=0 && static_cast<size_t>(i) < size()-m_deg-1,
                 "Index of Greville point is out of range.");
    iterator itr = begin() + 1;
    return ( m_deg==0 ? *(itr+i-1) :
             std::accumulate( itr+i, itr+i+m_deg, T(0.0) ) / m_deg
             // Special case C^{-1}
             - (*(itr+i) == *(itr+i+m_deg) ? 1e-10 : 0 )
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
        const T step = (*uit - prev)/T(knotsPerSpan+1);
        for( mult_t i = 1; i <= knotsPerSpan; ++ i)
            result.insert( result.end(), mult, prev + i*step );
        prev=*uit;
    }
}

template< typename T>
void gsKnotVector<T>::supportIndex_into(const mult_t& i,
                                               gsMatrix<unsigned>& result) const
{
    T suppBeg=*(this->begin()+i);
    T suppEnd=*(this->begin()+i+m_deg+1);
    uiterator ubeg   = this->ubegin();
    uiterator indBeg = findElement(suppBeg);
    uiterator indEnd = std::find_if(indBeg, this->uend(), std::bind2nd(std::greater_equal<T>(), suppEnd));
    result.resize(1,2);
    result<<indBeg-ubeg,indEnd-ubeg;
}

} // namespace gismo
