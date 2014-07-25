
#pragma once

#include <gsCore/gsLinearAlgebra.h>
# include <gsCore/gsForwardDeclarations.h>


namespace gismo {
  
/** 
    A univariate root solver for B-spline curves.
*/  
template<class T>
class gsBSplineSolver
{

public:
    /// Default empty constructor
    gsBSplineSolver() : m_n(1),m_k(1), m_d(1), eps(1e-7) { }

    /// Destructor
    ~gsBSplineSolver() { } 
  
public:
    /// Return true if the first root exist, the value of the root is in this->value()
    bool firstRoot(gsBSpline<T,gsKnotVector<T> > const & bsp, int const & coord = 0, 
                   T const & tr = 0, T const & tol = 1e-7, unsigned const & N=100);

    /// Return true if the first root exist, the value of the root is in this->value()
    bool firstRoot(gsBSpline<T, gsCompactKnotVector<T> > const & bsp, int const & coord = 0, 
                   T const & tr = 0, T const & tol = 1e-7, unsigned const & N=100);

    /// Next root (requires that first root has been called before)
    bool nextRoot ();

    /// The value of the current root
    inline T value () { return x;};

    // Return a vector with all the roots
    void allRoots (gsBSpline<T,gsKnotVector<T> > const & bsp, std::vector<T> & result, 
                   int const & coord = 0, 
                   T const & tr = 0, T const & tol = 1e-7, unsigned const & N=100);

    // Return a vector with all the roots
    void allRoots (gsBSpline<T,gsCompactKnotVector<T> > const & bsp, std::vector<T> & result, 
                   int const & coord = 0, 
                   T const & tr = 0, T const & tol = 1e-7, unsigned const & N=100);

private:
    /// Initialize the solver with B-spline data
    template<class KnotVectorType>
    void initSolver(gsBSpline<T,KnotVectorType> const & bsp , int const & coord, 
                    T const & tr, T const & tol, unsigned const &N);

    /// insert knot x in interval mu by Boehms algorithm
    /// Note: t,c must have size at least n+1, n+d+2 respectively
    /// require that x>=t[mu]
    int  insertKnot (int mu) ;

// Data members
private:

    // m_t: knot vector, m_c: coefficients
    std::vector<T> m_c, m_t;
    //gsVector<T> m_c, m_t;
    
    // m_n: coefficient size, maxn: maximum coeff. size // m_k span of root
    unsigned m_n, maxn, m_k ;

    int m_d;  // m_d: degree

    T eps; // tolerance

    T x; // root value
    
}; // class gsClass

} //namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////

#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsBSplineSolver.hpp)
#endif

