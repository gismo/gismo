/** @file gsBSplineSolver.h

    @brief Provides classes and functions to solve equations involving B-splines.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, A. Bressan
*/



#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsForwardDeclarations.h>

namespace gismo {


enum RootType
{
    odd,
    even,
    odd_interval,
    even_interval,
    invalid
};




template <typename T=real_t>
struct Root
{
    // requires cleanup
    RootType type;
    T parameter;
    T begParameter;
    T endParameter;
    gsVector<T> point;
    gsVector<T> begPoint;
    gsVector<T> endPoint;

    Root(bool is_odd, T par, const gsVector<T> &mpoint)
        : type(is_odd?odd:even)
        {
            parameter=par;
            point=mpoint;
        }
    Root(bool is_odd, T parB, T parE, const gsVector<T> &pointB, const gsVector<T> &pointE)
        : type(is_odd?odd_interval:even_interval)
        {
            begParameter=parB;
            endParameter=parE;
            begPoint=pointB;
            endPoint=pointE;
        }
    Root() : type(invalid){}
};


/**
    \brief find intersections of a BSpline curve with an hyperplane

    This function tries to be robust and to report correctly intersections
    of the curve with the given hyperplane of thickness 2*tolerance.
    Intersections are roots of a B-Spline curve and as such they can be of
    four types:

    -odd points, i.e. the position with respect of the hyperplane changes
     at the intersection

    -even points, i.e. the position with respect of the hyperplane is the same
     before and after the intersection

    -odd intervals, i.e. the curve stays for a full parametric interval in the
     hyperplane and then exit on the other side of it

    -even intervals, i.e. the curve stays for a full parametric interval in the
     hyperplane and then exit on the side it comes from

    The intersections are reported as encountered while following the curve
    in the direction given by increasing parameter.

    This function assumes an open knot vector: the first and last control points
    are in the curve.
    If the start or the end of the the curve is on the hyperplane it is reported
    as an odd intersection.
    If necessary this special case must be handled outside of this function, for
    instance when used to determine if a point is inside a 2D area bounded by a
    closed curve.
    
    \ingroup Nurbs
**/
template <typename T>
unsigned findHyperPlaneIntersections (
        const gsBSpline<T>    &curve,
        const gsVector<T>     &normal,
        T                      reference,
        T                      tolerance,
        std::vector<Root<T> > &roots
        );


  
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
    inline T value () { return x;}

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

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBSplineSolver.hpp)
#endif

