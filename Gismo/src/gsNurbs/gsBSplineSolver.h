
#pragma once

#include <gsCore/gsLinearAlgebra.h>
# include <gsCore/gsForwardDeclarations.h>


namespace gismo {


enum RootType
{
    odd,
    even,
    odd_interval,
    even_interval
};

enum Position
{
    // new scheme
    greater= 1,
    smaller=-greater,
    equal  = 1<<5,
    undef  =-1<<5
};


template <typename T>
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
            begPoint=pointE;
        }
};


template <typename T>
Position relativePosition (T pos, T ref)
{
    if ( pos < ref )
        return smaller;
    if ( pos > ref )
        return greater;
    if (ref==pos)
        return equal;
    return undef;
}

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

    If the beginning of the curve is on the hyperplane it is reported as odd.
    If necessary this special case must be handled outside of this function, for
    instance when used to determine if a point is inside a 2D area bounded by a
    closed curve.
**/
template <typename T>
unsigned findHyperPlaneIntersections (
        const gsBSpline<T>    &curve,
        const gsVector<T>     &normal,
        T                      reference,
        T                      tolerance,
        std::vector<Root<T> > &roots
        )
{
    // argument check
    GISMO_ASSERT( curve.coefDim() == normal.rows(),
        "cannot intersect a curve in R^"<< curve.coefDim()<<" with a hyperplane in R^"<< normal.rows() );
    if ( !math::isfinite(reference + normal.sum()) )
    {
        gsWarn<<"No intersection reportet between curve and invalid hyperplane."<<std::endl;
    }

    // environment description
    const unsigned deg = curve.degree();
    const T        tol = tolerance;
    const T        ref = reference;

    // internal copy of the curve used for knot insertion
    gsBSpline<T> crv(curve);

    // position status
    Position curP = relativePosition(crv.coef(0).dot(normal), ref) ; // current Position
    if ( curP == undef )
    {
        gsWarn<<__FUNCTION__<<": first coefficient of curve is NAN, exit without results."<<std::endl;
        return 0;
    }
    Position newP = undef;
    Position oldP = undef;
    bool     is_odd = true;

    // knots in the parameter domain
    T oldK = NAN; // last knot inserted
    T newK = NAN; // knot to be insert

    // counters
    unsigned curC = 1; // current coefficient (Control Point)
    unsigned lstC = 1; // coefficient in which the position became the current one
    unsigned rootCount = 0; // total number of roots found

    for ( ; curC<crv.coefsSize(); ++curC )
    {
        newP = relativePosition(crv.coef(curC).dot(normal), ref) ; // current Position
        if ( newP == undef )
        {
            gsWarn<<__FUNCTION__<<": stopping the search for intersection because the curve contains nan coefficients"<<std::endl;
            return rootCount;
        }
        if(curP==newP)
            continue;

        switch (curP)
        {
        case greater:
        case smaller:
            if (curP+newP!=0)
            {
                // from one side of the hyperplane to inside the hyperplane
                // we need to scan till the end of the intersection to determine
                // the type so we only change the state
                oldP=curP;
                curP=newP;
                lstC=curC;
                continue;
            }
            // otherwise follow the code to knot insertion
            is_odd = true;
            break;
        case equal:
            if (curC-lstC>= deg)
            {
                is_odd = newP!=oldP;
                // report even interval
                gsMatrix<T> params;
                params.resize(1,2);
                params(0,0) = crv.knots()[lstC+deg];
                params(0,1) = crv.knots()[curC+1];
                gsMatrix<T> points=crv.eval(params);
                roots.push_back(Root<T>(is_odd, params(0,0),params(0,1),points.col(0),points.col(1)));
                ++rootCount;
                curP=newP;
                continue;
            }
            else if (newP==oldP)
            {
                curP=newP;
                continue;
            }
            else
            {
                curC=lstC;
            }
            break;
        default:
            GISMO_ASSERT(false, __FUNCTION__<<"The impossible happened! Check for B..S!!" );
        }
        T grev1 = crv.knots()[curC];
        for (unsigned i=1;i<deg;++i)
            grev1 += (crv.knots()[curC+i]);
        T grev2 = grev1-(crv.knots()[curC])+(crv.knots()[curC+deg]);

        const T coeff1 = ref-crv.coef(curC-1).dot(normal);
        const T coeff2 = crv.coef(curC).dot(normal)-ref;
        const T denom  = (crv.coef(curC)-crv.coef(curC-1)).dot(normal)*deg;
        newK = (coeff2*grev1+coeff1*grev2)/denom;
        // in rare but real cases we go out of the interval due to numerical approximation
        newK = math::max(newK, curve.knots().first());
        newK = math::min(newK, curve.knots().last());

        if ( math::abs(newK-oldK) < tol )
        {
            gsMatrix<T> point=crv.eval(gsMatrix<T>::Constant(1,1,newK));
            if (math::abs(point.col(0).dot(normal)-ref)<tol)
            {
                roots.push_back(Root<T>(is_odd, newK,point.col(0)));
                ++rootCount;
                curP=newP;
                oldK=NAN;
                lstC=curC;
                continue;
            }
        }
        crv.insertKnot(newK);
        curC = math::max(curC-deg, lstC-1); // incremented by one on end of loop
        oldK = newK;
    }
    return rootCount;
}


  
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

#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsBSplineSolver.hpp)
#endif

