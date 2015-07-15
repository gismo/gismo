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

/// \brief Returns the factorial of \a n i.e. \a n!
/// Remember that factorial grow too fast and only n! with n<=13 can be stored in a 32bit that is an unsigned.
/// \ingroup combinatorics
/// \ingroup Utils
inline unsigned factorial( unsigned n)
{
    GISMO_ASSERT(n<14, "Overflow when computing factorial.");
    static const unsigned precomputed[]= {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600};
    return precomputed[n];
}

/*
/// \brief Returns the factorial of \a n
/// Remember that factorial grow too fast and only n! with n<=21 can
/// be stored in a 64bit that is an unsigned long long.  \ingroup
/// combinatorics
inline unsigned long long factorial( unsigned long long n)
{
    static const unsigned long long precomputed[]= {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600, 6227020800, 87178291200, 1307674368000, 20922789888000, 355687428096000, 6402373705728000, 121645100408832000, 2432902008176640000 };
    return precomputed[n];
}
*/

/**
 \brief Returns binomial(n,r)

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
  
  \ingroup Utils
*/
template <typename T>
inline T binomial(T n, T r)
{
    GISMO_ASSERT(r>=0, "binomial coefficient (n,r) exists only for n,r positive");
    //GISMO_ASSERT(n>=r, "binomial coefficient (n,r) exists only for n,r such that n>=r");

    const T diff = math::min( n-r, r );
    int result = 1;
    for (T i=0;i < diff;)
    {
        result *= n-i;
        result /= ++i;
    }
    return result;
}

// forward declaration
template<int n, int r> class binomialT;

/**
   \brief Returns binomial(n,r) as a compile time constant

   This is done using template recursion and can be accessed either by
   b=binomialT<n,r>::value;
   or
   b=binomial<n,r>();
   The second form relies on the compiler optimizations to avoid function call.
   \ingroup combinatorics
   \ingroup Utils
*/
template <unsigned n, unsigned r>
inline unsigned binomial () {return binomialT<n,r>::value;}

// TEMPLATE IMPLEMENTATION NOT DOCUMENTED, BASED ON PASCAL TRIANGLE
template<int n, int r>
class binomialT
{
public:
   enum { value= binomialT<n-1,r-1>::value+binomialT<n-1,r>::value};
};

template<int n>
class binomialT<n,n>
{
public:
   enum { value= 1};
};

template< >
class binomialT<0,0>
{
public:
   enum { value= 1};
};

template<int n>
class binomialT<n,0>
{
public:
   enum { value= 1};
};



/** \brief Returns a vector containing all the binomials (n,r) with n fixed.
 *
 * This functions compute all the binomial coefficients with fixed \a n
 * and store them in the supplied vector \a v.
 * Uses the Pacal triangle rule.
 *
 * \param[in] n
 * \param[out] v a vector with n+1 components each equal to
 * \ingroup combinatorics
 * \ingroup Utils
 */
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





  /** 
      @brief Class for combinatorics. Generates ...
      
      \ingroup Utils
  */

template<class Vect = gsVector<unsigned> >
class gsCombinat
{

public:

  /// Default empty constructor
  gsCombinat() : n(0) { }

  //destructor
  ~gsCombinat() { } 

public:

  /// Construct first combination of {0,..,n-1}
  void first_combination ( unsigned const& n_, unsigned const & r, Vect & res)
  {
    n=n_;
    if (r<=n) 
      res= Vect::LinSpaced(r,0,r-1);
    else
      std::cerr << "gsCombinat: r>n combination requested. r="<< r<<", n="<< n<<"\n";
  }
  
  /// Construct NEXT combination of {0,..,n-1}
  /// first_combination should have been called before this
  bool next_combination ( Vect & v)
  {
    int r= v.rows() ;            
    if  (v == Vect::LinSpaced(r,n-r,n-1)) return false;
    int i=r-1;
    while (v[i] == n-r+i ) --i;
    v[i] += 1;
    for (int j=i+1; j<r; ++j)
      v[j] = v[i] +j-i;
    return true;
  }

  /// Number of combinations (binomial coefficient)
  /// first_combination should have been called before this
  unsigned combinations(unsigned const & r)
  {
    return binomial(n,r);
  }
  
  /// Number of combinations (binomial coefficient)
  static unsigned combinations( unsigned const& n_, unsigned const & r_)
  {
    return binomial(n_,r_);
  }
  

  /// Number of permutations (binomial coefficient)
  /// first_permutation should have been called before this
  unsigned permutations()
  {
    return factorial(n);
  }

  /// Number of permutations (binomial coefficient)
  static unsigned permutations( unsigned const& n_ )
  {
    //return exp(lgamma(n + 1));// computation using the Gamma function
    return factorial(n_);
  }
  
  // To do:
  //first_sum
  //next_sum
  //sums
  
  // To do:
  //first_binary
  //next_binary
  //binaries
  
  
  /// Construct first lattice point in the rectangle [l, u]
  void first_lattice_point( const Vect& l, const Vect& u, Vect & res)
  { 
    GISMO_ASSERT(l.size()==u.size(), 
		 "Lattice endpoints have invalid dimension");
    low=l; upp=u; res=low;     
}
  
  /// Construct first lattice point in the rectangle [0, u]
  void first_lattice_point( const Vect& u, Vect & res)
  { low.setZero(u.size()); upp=u; res=low; }
  
  /// Next lattice point in [low,upp] after v in lexicographic order
  /// first_lattice_point should have been called before
  inline bool next_lattice_point( Vect & v) const 
  {
    if ( v==upp ) return false;// finished
    //note: May exist a better implementation with: std::mismatch,  std::pair

    int i(0);
    while ( v[i] == upp[i]) { ++i; }
    v[i]+=1;

    --i;
    while ( i>=0 )
      {
        v[i]= low[i];
        --i;
      }
    
    return true;
  }

  /// Number of lattice points in [low,upp]
  /// first_lattice_point should have been called before
  unsigned lattice_points() const 
  {
    return (upp-low+Vect::Ones(upp.size()) ).prod();
  }
  
  /// Number of lattice points in [l,u]
  static unsigned lattice_points( const Vect& l, const Vect& u)
  {
    return (u-l+Vect::Ones(u.size()) ).prod();
  }


  /// Construct first composition
  void first_composition( const Vect& u, Vect & res)
  {  

  }
  
  /// Next composition
  /// first_composition should have been called before
  inline bool nextComposition( Vect & v) const 
  {

    return true;
  }

  /// Number of compositions
  /// first_lattice_point should have been called before
  unsigned numCompositions() const 
  {
      return 1;
  }

  
// Data members
private:
  Vect low, upp;

  unsigned n;

}; // class gsCombinat


/**
 * \brief changes current to the first permutation of 0 ... size(current)-1
 * note that you must resize the vector to specify the number of elements
 */
template<class Vec>
void firstPermutation (Vec &current)
{
    const index_t n=current.size();
    current=Vec::LinSpaced(n,0,n-1);
}

/**
 * \brief changes current to the next lexicographically ordered permutation
 * \return false when the lexicographically last permutation is given
 */
template<class Vec>
bool nextPermutation (Vec &current)
{
    const index_t n=current.size();
    return std::next_permutation(current.data(), current.data()+n);
}


/// iterate through a tensor lattice with the given size. Updates cur
/// and returns true if another entry was available End values (\a
/// size) are not included in the enumerated points, as with
/// iterators.
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


/// Iterate through a tensor lattice with the given start and end
/// points.  \a end coordinates are not included in the enumerated
/// points, as with iterators.  Updates cur and returns true if
/// another entry was available.
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


/// Iterate in lexicographic order through the vertices of the cube
/// [start,end]. Updates cur with the current vertex and returns true
/// if another vertex is available. Cube may be degenerate.
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

/// Iterate in lexicographic order through the vertices of the cube
/// [0,end]. Updates cur with the current vertex and returns true
/// if another vertex is available. Cube may be degenerate.
template<class Vec>
bool nextCubeVertex(Vec& cur, const Vec& end)
{
    const int d = cur.size();
    GISMO_ASSERT( d == end.size(),
        "Vector sizes don't match in nextCubeVertex");

    for (int i = 0; i != d; ++i)
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

/// Iterate in lexigographic order through the points of the integer
/// lattice contained in the cube [0,end]. Updates cur with the
/// current point and returns true if another point is available. Cube
/// may be degenerate.
template<class Vec>
bool nextCubePoint(Vec& cur, const Vec& end)
{
    const int d = cur.size();
    GISMO_ASSERT(d == static_cast<int>(end.size()),
        "Vector sizes don't match in nextCubePoint");

    for (int i = 0; i != d; ++i)
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

/// Iterate in lexigographic order through the points of the integer
/// lattice contained in the cube [start,end]. Updates cur with the
/// current point and returns true if another point is available. Cube
/// may be degenerate.
template<class Vec>
bool nextCubePoint(Vec& cur, const Vec& start, const Vec& end)
{
    const int d = cur.size();
    GISMO_ASSERT( d == static_cast<int>(start.size()) &&
                  d == static_cast<int>(end.size()),
        "Vector sizes don't match in nextCubePoint");

    for (int i = 0; i != d; ++i)
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

/// iterate in lex-order through the boundary points of the cube
/// [start,end]. Updates cur with the current point and returns true
/// if another point is available. Cube may be degenerate.
template<class Vec>
bool nextCubeBoundary(Vec& cur, const Vec& start, const Vec& end)
{
    const int d = cur.size();
    GISMO_ASSERT( d == start.size() && d == end.size(),
        "Vector sizes don't match in nextCubeBoundary");

    for (int i = 0; i != d; ++i)
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

/// iterate in lex-order through the boundary points of the cube
/// [start,end], with an \ offset to the interior. Updates cur with
/// the current point and returns true if another point is
/// available. Cube may be degenerate.
template<class Vec>
bool nextCubeBoundaryOffset(Vec& cur, const Vec& start, const Vec& end, Vec & offset)
{
    const int d = cur.size();
    GISMO_ASSERT( d == start.size() && d == end.size(),
        "Vector sizes don't match in nextCubeBoundaryOffset");

    for (int i = 0; i != d; ++i)
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

/// iterate in lex-order through the boundary points of the cube
/// [start,end], with offset \a loffset from \ start and \a roffset
/// .from the \a end. Updates cur with the current point and returns
/// true if another point is available. Cube may be degenerate.
template<class Vec>
bool nextCubeBoundaryOffset(Vec& cur, const Vec& start, const Vec& end, 
                            Vec & loffset, Vec & uoffset)
{
    const int d = cur.size();
    GISMO_ASSERT( d == start.size() && d == end.size(),
        "Vector sizes don't match in nextCubeBoundaryOffset");

    for (int i = 0; i != d; ++i)
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


/// Construct first composition of \a sum into \a dim integers
template<class Vec>
void firstComposition( typename Vec::Scalar_t sum, index_t dim, Vec & res)
{  
    res.setZero(dim);
    res[0] = sum;
}

/// Next composition in lexicographic order
template<class Vec>
inline bool nextComposition(Vec & v)
{
    const typename Vec::Scalar_t sum = v.sum();
    const index_t k   = v.size() - 1;
    
    if (v[k] != sum)
    {
        for (index_t i = 0; i <= k; i++)
        {
            if ( v[i]!=0 )
            {
                const typename Vec::Scalar_t t = v[i];
                v[i]    = 0  ;
                v[0]    = t-1;
                v[i+1] += 1  ;
                return true;
            }
        }
    }
    return false;
}

inline
unsigned numCompositions(int sum, int dim)
{
    return binomial(sum+dim-1,dim-1);
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


}; // namespace gismo
