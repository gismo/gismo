/** @file gsCombinatorics.h

    @brief Provides combinatorial unitilies.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Bressan, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMath.h>

namespace gismo
{

/** \brief Returns the factorial of \a n i.e. \a n!
 * Remember that factorial grow too fast and only n! with n<=13 can be
 * stored in a 32bit that is an unsigned.
 * \ingroup combinatorics
 */
// \ingroup Utils
// could also be in Utils, but doxygen allows only one group for free functions
inline unsigned factorial( unsigned n)
{
    static const unsigned precomputed[]= {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600};
    GISMO_ASSERT(n < util::size(precomputed), "Overflow when computing factorial.");
    return precomputed[n];
}

/*
// \brief Returns the factorial of \a n
// Remember that factorial grow too fast and only n! with n<=21 can
// be stored in a 64bit that is an unsigned long long.  \ingroup
// combinatorics
inline unsigned long long factorial( unsigned long long n)
{
static const unsigned long long precomputed[]= {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000 };
GISMO_ASSERT(n < util::size(precomputed), "Overflow when computing factorial.");
return precomputed[n];
}
*/

/**
   \brief Computes the binomial expansion coefficient binomial(n,r)

   The binomial coefficient indexed by \a n and \a k is is the
   coefficient of the $x^k$ term in the polynomial expansion of the
   binomial power $(1 + x)^n$.

   This functions computes a single binomial coefficient, if many
   binomial coefficients with fixed \a n are needed then it is probably
   faster use the function binomial_into.  This function uses a loop
   implementation.  

   \param n binomial power
   \param r term of the binomial expansion

   \ingroup combinatorics
*/
// \ingroup Utils
// Note: could also be in Utils, but doxygen allows only one group for free functions
template <typename Z>
inline Z binomial(const Z n, const Z r)
{
    GISMO_ASSERT(r>=0, "binomial coefficient (n,r) exists only for n,r positive");
    //GISMO_ASSERT(n>=r, "binomial coefficient (n,r) exists only for n,r such that n>=r");

    const Z diff = math::min( n-r, r );
    int result = 1;
    for (Z i=0;i < diff;)
    {
        result *= n-i;
        result /= ++i;
    }
    return result;
}

/** \brief Returns a vector containing all the binomials (n,r) with n fixed.
 *
 * This functions compute all the binomial coefficients with fixed \a n
 * and store them in the supplied vector \a v.
 * Uses the Pacal triangle rule.
 *
 * \param[in] n
 * \param[out] v a vector with n+1 components each equal to
 * \ingroup combinatorics
 */
//\ingroup Utils
inline void binomial_into( unsigned n, gsVector<unsigned>& v)
{
    v.resize (n+1);
    v[0]=1;

    for ( unsigned i = 1; i <= n; ++i)
    {
        v[i]= 1;
        for ( unsigned j= i-1; j>0; --j)
            v[j] += v[j-1];
    }
}

// Used by binomial<n,r>()
template<int n, int r> class binomialT
{public:enum { value= binomialT<n-1,r-1>::value+binomialT<n-1,r>::value};};
template<int n> class binomialT<n,n> {public:enum { value= 1};};
template<>      class binomialT<0,0> {public:enum { value= 1};};
template<int n> class binomialT<n,0> {public:enum { value= 1};};
/**
   \brief Returns binomial(n,r) as a compile time constant

   This is done using template recursion and can be accessed either by
   b=binomialT<n,r>::value;
   or
   b=binomial<n,r>();
   The second form relies on the compiler optimizations to avoid function call.
   \ingroup combinatorics
*/
template <unsigned n, unsigned r>
inline unsigned binomial () {return binomialT<n,r>::value;}

/// \brief Computes the first \a r-combination of {0,..,n-1}
/// \ingroup combinatorics
template<class Vec>
void firstCombination(const unsigned n, const unsigned r, Vec & res)
{
    if (r<=n) 
        res= Vec::LinSpaced(r,0,r-1);
    else
        std::cerr << "Error: r>n combination requested. r="<< r<<", n="<< n<<"\n";
}
  
/// \brief Computes the next r-combination of {0,..,n-1}, where r = \a v.size().
/// The input \a v is expected to be a valid combination
/// \ingroup combinatorics
template<class Vec>
bool nextCombination (Vec & v, const unsigned n)
{
    const index_t r = v.rows() ;            
    if  (v == Vec::LinSpaced(r,n-r,n-1)) return false;
    int i = r-1;
    while (v[i] == n-r+i) --i;
    v[i] += 1;
    for (index_t j=i+1; j<r; ++j)
        v[j] = v[i] +j-i;
    return true;
}

/**
 * \brief changes current to the first permutation of 0 ... size(current)-1
 * note that you must resize the vector to specify the number of elements
 * \ingroup combinatorics
 */
template<class Vec>
void firstPermutation (Vec &current)
{
    const index_t n=current.size();
    current=Vec::LinSpaced(n,0,n-1);
}

/**
 * \brief Changes current to the next lexicographically ordered permutation
 * \return false when the lexicographically last permutation is given
 * \ingroup combinatorics
 */
template<class Vec>
bool nextPermutation (Vec &current)
{
    const index_t n=current.size();
    return std::next_permutation(current.data(), current.data()+n);
}


/// \brief Iterates through a tensor lattice with the given size. Updates cur
/// and returns true if another entry was available End values (\a
/// size) are not included in the enumerated points, as with
/// iterators.
/// \ingroup combinatorics
template<class Vec>
bool nextLexicographic(Vec& cur, const Vec& size)
{
    const index_t d = cur.size();
    GISMO_ASSERT( d == size.size(), "Vector sizes don't match in nextLexicographic" );

    for (index_t i = 0; i < d; ++i)
    {
        // increase current dimension
        if (++cur[i] == size[i])    // current dimension exhausted ?
        {
            if (i == d - 1)         // was it the last one?
                return false;       // then all elements exhausted
            else
                cur[i] = 0;         // otherwise, reset this to 0 and increase the next dimension
        }
        else
            return true;            // current dimension not yet exhausted, return current vector
    }
    GISMO_ERROR("Something went wrong in nextLexicographic, wrong input?");
}


/// \brief Iterate through a tensor lattice with the given start and end
/// points.  \a end coordinates are not included in the enumerated
/// points, as with iterators.  Updates cur and returns true if
/// another entry was available.
/// \ingroup combinatorics
template<class Vec>
bool nextLexicographic(Vec& cur, const Vec& start, const Vec& end)
{
    const index_t d = cur.size();
    GISMO_ASSERT( d == start.size() && d == end.size(),
                  "Vector sizes don't match in nextLexicographic");

    for (index_t i = 0; i < d; ++i)
    {
        // increase current dimension
        if (++cur[i] == end[i])     // current dimension exhausted ?
        {
            if (i == d - 1)         // was it the last one?
                return false;       // then all elements exhausted
            else
                cur[i] = start[i];  // otherwise, reset this and increase the next dimension
        }
        else
            return true;            // current dimension not yet exhausted, return current vector
    }
    GISMO_ERROR("Something went wrong in nextLexicographic, wrong input?");
}


/// \brief Iterate in lexicographic order through the vertices of the cube
/// [start,end]. Updates cur with the current vertex and returns true
/// if another vertex is available. Cube may be degenerate.
/// \ingroup combinatorics
template<class Vec>
bool nextCubeVertex(Vec& cur, const Vec& start, const Vec& end)
{
    const int d = cur.size();
    GISMO_ASSERT( d == start.size() && d == end.size(),
                  "Vector sizes don't match in nextCubeVertex");

    for (int i = 0; i != d; ++i)
    {
        if ( cur[i] != end[i] )
        {
            cur[i] = end[i];
            return true;
        }
        else
            cur[i] = start[i];
    }
    return false;
}

/// \brief Iterate in lexicographic order through the vertices of the cube
/// [0,end]. Updates cur with the current vertex and returns true
/// if another vertex is available. Cube may be degenerate.
/// \ingroup combinatorics
template<class Vec>
bool nextCubeVertex(Vec& cur, const Vec& end)
{
    const index_t d = cur.size();
    GISMO_ASSERT( d == end.size(),
                  "Vector sizes don't match in nextCubeVertex");

    for (index_t i = 0; i != d; ++i)
    {
        if ( cur[i] != end[i] )
        {
            cur[i] = end[i];
            return true;
        }
        else
            cur[i] = 0;
    }
    return false;
}

/// \brief Iterate in lexigographic order through the vertices of the
/// unit cube. Updates \a cur with the lexicographically next vertex
/// and returns true if another point is available. This is equivalent
/// to iterating over all possible binary sequences of length \em cur.size().
/// The input \a cur is expected to contain only zeros and ones (or true/false).
/// \ingroup combinatorics
template<class Vec>
bool nextCubeVertex(Vec& cur)
{
    const index_t d = cur.size();
    GISMO_ASSERT( (cur.array() >= 0).all() && (cur.array() <= 1).all(),
                  "Input must be a vector of zeros and ones, got: "<<cur.transpose() );

    for (index_t i = 0; i != d; ++i)
    {
        if ( cur[i] == 0 )
        {
            ++cur[i];
            return true;
        }
        else
            cur[i] = 0;
    }
    return false;
}

/// \brief Iterate in lexigographic order through the points of the integer
/// lattice contained in the cube [0,end]. Updates cur with the
/// current point and returns true if another point is available. Cube
/// may be degenerate.
/// \ingroup combinatorics
template<class Vec>
bool nextCubePoint(Vec& cur, const Vec& end)
{
    const index_t d = cur.size();
    GISMO_ASSERT(d == static_cast<int>(end.size()),
                 "Vector sizes don't match in nextCubePoint");

    for (index_t i = 0; i != d; ++i)
    {
        if ( cur[i] != end[i] )
        {
            ++cur[i];
            return true;
        }
        else
            cur[i] = 0;
    }
    return false;
}

/// \brief Iterates in lexigographic order through the points of the integer
/// lattice contained in the cube [start,end]. Updates cur with the
/// current point and returns true if another point is available. Cube
/// may be degenerate.
/// \ingroup combinatorics
template<class Vec>
bool nextCubePoint(Vec& cur, const Vec& start, const Vec& end)
{
    const index_t d = cur.size();
    GISMO_ASSERT( d == static_cast<int>(start.size()) &&
                  d == static_cast<int>(end.size()),
                  "Vector sizes don't match in nextCubePoint");

    for (index_t i = 0; i != d; ++i)
    {
        if ( cur[i] != end[i] )
        {
            ++cur[i];
            return true;
        }
        else
            cur[i] = start[i];

    }
    return false;
}

/// \brief Iterates in lex-order through the boundary points of the cube
/// [start,end]. Updates cur with the current point and returns true
/// if another point is available. Cube may be degenerate.
/// \ingroup combinatorics
template<class Vec>
bool nextCubeBoundary(Vec& cur, const Vec& start, const Vec& end)
{
    const index_t d = cur.size();
    GISMO_ASSERT( d == start.size() && d == end.size(),
                  "Vector sizes don't match in nextCubeBoundary");

    for (index_t i = 0; i != d; ++i)
    {        
        if ( cur[i] != end[i] )
        {
            if ( cur[i] == start[i] && ( i!=d-1 || d==1) ) // || d==1 to treat 1D
            {
                int c=i+1;
                for (int j = c; j!=d; ++j)
                    if ( (cur[j] == start[j]) || 
                         (cur[j] == end[j]) )
                        c++;
                
                if ( c==1 )
                    cur[i] = end[i];
                else
                    cur[i]++;
            }
            else
                cur[i]++;

            for (int k = i-1; k!=-1; --k)
                cur[k] = start[k];
            return true;
        }
    }
    return false;
}

/// \brief Iterates in lex-order through the boundary points of the cube
/// [start,end], with an \ offset to the interior. Updates cur with
/// the current point and returns true if another point is
/// available. Cube may be degenerate.
/// \ingroup combinatorics
template<class Vec>
bool nextCubeBoundaryOffset(Vec& cur, const Vec& start, const Vec& end, Vec & offset)
{
    const index_t d = cur.size();
    GISMO_ASSERT( d == start.size() && d == end.size(),
                  "Vector sizes don't match in nextCubeBoundaryOffset");

    for (index_t i = 0; i != d; ++i)
    {        
        if ( cur[i] != end[i] )
        {
            if ( cur[i] == start[i]+offset[i] && ( i!=d-1 || d==1) ) // || d==1 to treat 1D
            {
                int c = i+1;
                for (int j = c; j!=d; ++j)
                    if (cur[j] <=  start[j] + offset[j] ||
                        cur[j] >=    end[j] - offset[j] )
                        c++;
                
                if ( c==1 )
                    cur[i] = end[i] - offset[i];
                else
                    cur[i]++;
            }
            else
                cur[i]++;

            for (int k = i-1; k!=-1; --k)
                cur[k] = start[k];
            return true;
        }
    }

    return false;
}

/// \brief Iterates in lex-order through the boundary points of the cube
/// [start,end], with offset \a loffset from \ start and \a roffset
/// .from the \a end. Updates cur with the current point and returns
/// true if another point is available. Cube may be degenerate.
/// \ingroup combinatorics
template<class Vec>
bool nextCubeBoundaryOffset(Vec& cur, const Vec& start, const Vec& end, 
                            Vec & loffset, Vec & uoffset)
{
    const index_t d = cur.size();
    GISMO_ASSERT( d == start.size() && d == end.size(),
                  "Vector sizes don't match in nextCubeBoundaryOffset");

    for (index_t i = 0; i != d; ++i)
    {        
        if ( cur[i] != end[i] )
        {
            if ( cur[i] == start[i]+loffset[i] && ( i!=d-1 || d==1) ) // || d==1 to treat 1D
            {
                int c = i+1;
                for (int j = c; j!=d; ++j)
                    if (cur[j] <=  start[j] + loffset[j] ||
                        cur[j] >=    end[j] - uoffset[j] )
                        c++;
                
                if ( c==1 )
                    cur[i] = end[i] - uoffset[i];
                else
                    cur[i]++;
            }
            else
                cur[i]++;

            for (int k = i-1; k!=-1; --k)
                cur[k] = start[k];
            return true;
        }
    }

    return false;
}

/// \brief Returns the number of elements (faces) of dimension \a k 
/// of a \a d-cube
/// \ingroup combinatorics
inline index_t numCubeElements(const index_t k, const index_t d)
{
    GISMO_ASSERT(k >= 0 && k<=d, "Invalid arguments.");
    return binomial(d,k) * (1<<(d-k));
}

/// \brief Returns the dimension of an element (face) of the \a d-cube
/// [0,1]^d. The element is expected to contain 0,1 (corresponding to
/// cube extrema) and the special value 2 at the position of "free"
/// coordinates
/// \ingroup combinatorics
template<class Vec>
inline index_t dimCubeElement(const Vec & cur)
{
    return (cur.array() == 2).count();
}

/// \brief Updates \a cur to contain the lexicographically first
/// element (face) of the cube [0,1]^d of dimension \a k. For k==d the
/// face (2..2) is returned, corresponding to the cube itself.
/// \ingroup combinatorics
template<class Vec>
void firstCubeElement(Vec & cur, const index_t k = 0)
{
    const index_t d = cur.size();
    GISMO_ASSERT(k >= 0 && k<=d, "Invalid arguments.");
    cur.topRows   (k  ).setConstant(2);
    cur.bottomRows(d-k).setConstant(0);
}

/// \brief Iterates in lexicographic order through the elements
/// (faces) of dimension \a k of the cube [0,1]^d. Updates \a cur with
/// the current element (face) and returns true if another element
/// (face) of dimension \a k is available.  Coordinates with value 2
/// indicate free/not-fixed dimensions.
/// \ingroup combinatorics
template<class Vec>
bool nextCubeElement(Vec & cur, const index_t k)
{
    const index_t d = cur.size();
    GISMO_ASSERT(k >= 0 && k<=cur.size(), "Invalid arguments.");

    index_t i;
    do
    {
        for (i = 0; i != d; ++i)
        {
            if ( cur[i] != 2 )
            {
                ++cur[i];
                if ( (cur.array() == 2).count() == k ) // dimCubeElement(cur)==k ?
                    return true;
                else
                    break;// skip face 
            }
            else
                cur[i] = 0;
        }
    }
    while ( i!=d );

    return false;
}

/// \brief Computes the isometry of the unit d-cube 
/// implied by a permutation \a perm of the cube directions
/// plus a relocation \a flip of the cube vertices
///
/// \param[in] flip the relocation of the cube vertices
/// flip[k]==true  : the coordinate \em k of the vertex is not relocated
/// flip[k]==false : the coordinate \em k of the vertex is relocated
/// \param[in] perm the permutation of the cube directions (0,..,d-1)
/// \param[out] result A permutation of the vertices (0,..,2^d-1)
/// \ingroup combinatorics
template <typename Z, int d>
void cubeIsometry( const gsVector<bool,d>    & flip,
                   const gsVector<index_t,d> & perm, 
                   gsVector<Z> & result)
{
    const index_t dd = flip.size(); //binary sequence of length d
    GISMO_ASSERT( dd == perm.size(), "Dimensions do not match in cubeIsometry");
    GISMO_ASSERT( perm.sum() == dd*(dd-1)/2, "Error in the permutation: "<< perm.transpose());

    gsVector<index_t,d> pstr(dd);
    for (index_t k=0; k!=dd; ++k)
        pstr[k] = (1<<perm[k]);

    result.resize(1<<dd);
    index_t r  = 0;
    index_t i;
    gsVector<bool,d> v = gsVector<bool,d>::Zero(dd);
    do
    {
        Z & c = result[r++];
        c = 0;
        for (index_t k=0; k!=dd; ++k)
            c += ( flip[perm[k]] == v[k] ) * pstr[k];
        
        for (i = 0; i != dd; ++i)
        {
            if ( !v[i] )
            {
                v[i] = true;
                break;//for
            }
            else
                v[i] = false;
        }
    }
    while (i!=dd);
}

/// \brief Computes the rotation matrix
/// implied by a permutation \a perm of the cube directions
/// plus a relocation \a flip
///
/// \param[in] flip the relocation of the cube vertices
/// flip[k]==true  : the coordinate \em k is not relocated
/// flip[k]==false : the coordinate \em k is relocated
/// \param[in] perm the permutation of the directions (0,..,d-1)
/// \param[out] result A rotation matrix
/// \ingroup combinatorics
template <int d>
void cubeIsometryMatrix ( const gsVector<bool,d>    & flip,
                          const gsVector<index_t,d> & perm,
                          gsMatrix<int,d,d> & result)
{
    result.setZero(d,d);
    for(index_t i = 0; i < d; ++i)
        result(perm(i),i) = (flip(perm(i)) ? 1 : -1);
}

/// \brief Construct first composition of \a sum into \a dim integers
/// \ingroup combinatorics
template<class Vec>
void firstComposition( typename Vec::Scalar sum, index_t dim, Vec & res)
{  
    res.derived().resize(dim);
    res.setZero();
    res[0] = sum;
}

/// \brief Returns (inplace) the next composition in lexicographic order
/// \ingroup combinatorics
template<class Vec>
inline bool nextComposition(Vec & v)
{
    const index_t k = v.size() - 1;
    
    if (v[k] != v.sum())
    {
        for (index_t i = 0; i <= k; ++i)
        {
            if ( v[i]!=0 )
            {
                const typename Vec::Scalar t = v[i];
                v[i]    = 0  ;
                v[0]    = t-1;
                v[i+1] += 1  ;
                return true;
            }
        }
    }
    return false;
}


/// \brief Number of compositions of \a sum into \a dim integers
/// \ingroup combinatorics
inline unsigned numCompositions(int sum, short_t dim)
{
    return binomial(sum+dim-1,dim-1);
}


/** \brief Constructs first multi-composition of \a a = (a_1,..,a_d) into \a k integers
 \ingroup combinatorics
*/
template<class Vec, class Mat>
void firstMultiComposition(const Vec & a, index_t k, Mat & res)
{  
    const index_t d = a.size();
    res.setZero(k, d);
    res.row(0) = a.transpose();
}

/**
   \brief Returns (inplace) the next multi-composition 
   in lexicographic order

   \f$ m \in \mathbb N^{k\times d} \f$

   \ingroup combinatorics
*/ 
template<class Mat>
inline bool nextMultiComposition(Mat & m)
{
    const index_t k = m.rows();
    const index_t d = m.cols();
           
    for (index_t j = 0; j != d; ++j)
    {
        typename Mat::ColXpr c_j = m.col(j);
        if ( nextComposition(c_j) )
            return true;
        else
        {
            const index_t n_j = c_j.sum();
            firstComposition(n_j, k, c_j);
        }
    }
    return false;
}


/** \brief Number of multi-composition of \a a = (a_1,..,a_d) into \a k integers
    \ingroup combinatorics
*/
template<class Vec>
unsigned numMultiCompositions(const Vec & a, index_t k)
{
    unsigned result = 1;
    const index_t d = a.size();
    for (index_t j = 0; j != d; ++j)
        result *= binomial<unsigned>(a[j]+k-1, k-1);
    return result;
}


// Lexicographic comparison for gsVectors
template<class Z, unsigned d>
struct lex_less
{
    bool operator() (const gsVector<Z,d> & lhs, const gsVector<Z,d> & rhs) const
    {
        unsigned i = 0;
        while( (i<d) && (lhs[i] == rhs[i++]) ) ;
        return ( --i==d ? false : lhs[i] < rhs[i] );
    }
};


} // namespace gismo
