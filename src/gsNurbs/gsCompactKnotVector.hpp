/** @file gsCompactKnotVector.hpp

    @brief Provides implementation of the CompactKnotVector class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsIO/gsXml.h>


namespace gismo
{

template <class T>
gsCompactKnotVector<T>::gsCompactKnotVector(T const u0, T const u1, unsigned const interior,
                                            unsigned const mult_ends, 
                                            unsigned const mult_interior,
                                            int const degree)
{
    initUniform( u0, u1, interior, mult_ends, mult_interior, degree );
}

template <class T>
void gsCompactKnotVector<T>::initUniform(T const u0, T const u1, 
                                         unsigned const interior,
                                         unsigned const mult_ends, 
                                         unsigned const mult_interior, 
                                         int const degree)
{
    T h= (u1-u0)/(interior+1);

    m_knots   .clear();
    m_mult_sum.clear();
    m_knots   .reserve(interior+2);
    m_mult_sum.reserve(interior+2);
    m_knots   .push_back(u0);
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
    GISMO_ASSERT(knots.size() >= 2, "We need at least two knots to define a meaningful knot vector.");
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
    // Also ASSUMES that start + 1 <= end.
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
    return gsCompactKnotVector<T>::const_iterator(*this,static_cast<size_t>(findElementIndex(u)));
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
    size_t c = 1;
    for ( typename  std::vector<unsigned>::iterator it= m_mult_sum.begin();
          it != m_mult_sum.end(); ++it )
        *it += i * c++;
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
void gsCompactKnotVector<T>::degreeDecrease(int const & i)
{
    m_p -= i;
    trim(i);
}

template <class T>
bool gsCompactKnotVector<T>::contains(gsCompactKnotVector<T> & other)
{
    GISMO_NO_IMPLEMENTATION
}

template <class T>
void gsCompactKnotVector<T>::remove(T const& knot)
{
    GISMO_NO_IMPLEMENTATION
/*    typename std::vector<T>::iterator itr =
        std::lower_bound(m_knots.begin(), m_knots.end(), knot);
    m_mult_sum.erase(m_mult_sum.begin() + itr-m_knots.begin() ) ;
    //to do: update mults tail
    m_knots.erase(itr);
*/
}

template <class T>
void gsCompactKnotVector<T>::trim(int i)
{
    std::transform(m_mult_sum.begin(), m_mult_sum.end(), m_mult_sum.begin(),
                   std::bind2nd(std::minus<unsigned>(), i) );
    m_mult_sum.back() -= i;
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
        result_t result = std::equal_range(m_knots.begin(), m_knots.end(), knot);
        if ( result.first == result.second )
            return 0;
        const size_t j = result.first-m_knots.begin();
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
    typename std::vector<T>::const_iterator itr;
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













namespace internal
{

/// Get a KnotVector from XML data
///
/// \ingroup Nurbs
template<class T>
class gsXml< gsCompactKnotVector<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsCompactKnotVector<T>);
    static std::string tag () { return "KnotVector"; }
    static std::string type() { return ""; } // "Compact" ?
    
    static gsCompactKnotVector<T> * get (gsXmlNode * node)
    {
        GISMO_ASSERT( ! strcmp( node->name(), "KnotVector"), "Invalid tag" );
        // && node->first_attribute("type")->value() is "Plain" or "Compact");
        
        int p = atoi(node->first_attribute("degree")->value() );
        gsCompactKnotVector<T> * kv = new gsCompactKnotVector<T>(p);        
        
        std::istringstream str;
        str.str( node->value() );
        for (T knot; str >> knot;) 
            kv->push_back(knot);
        
        return kv;
    }
    
    static gsXmlNode * put (const gsCompactKnotVector<T> & obj, gsXmlTree & data)
    {
        // Write the knot values (for now WITH multiplicities)            
        std::ostringstream str;
        str << std::setprecision(FILE_PRECISION);

        for ( typename gsCompactKnotVector<T>::const_iterator it = obj.begin();
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
        
        // todo : append "Compact" and write out as compact ?
        
        return tmp;
    }
};

}// namespace internal

}// namespace gismo
 
