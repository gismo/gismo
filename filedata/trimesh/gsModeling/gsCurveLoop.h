/** @file gsCurveLoop.h

    @brief Interface for gsCurveLoop class, representing a loop of curves, in
    anticlockwise order.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Falini, A. Mantzaflaris, D.-M. Nguyen, M. Pauley
*/

#pragma once

#include <iostream>

# include <gsCore/gsCurve.h>
# include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

/** @brief
    A closed loop given by a collection of curves.
    
    The loop may be oriented clockwise (CW) or counterclockwise (CCW).

    \tparam T coefficient type
    
    \ingroup Modeling
*/

template<class T>
class gsCurveLoop
{

public:
    typedef typename std::vector< gsCurve<T> * >::iterator       iterator;
    typedef typename std::vector< gsCurve<T> * >::const_iterator const_iterator;

    /// Shared pointer for gsCurveLoop
    typedef memory::shared_ptr< gsCurveLoop > Ptr;

    /// Unique pointer for gsCurveLoop
    typedef memory::unique_ptr< gsCurveLoop > uPtr;

public:
    /// Default empty constructor
    gsCurveLoop() { }

    /// Constructor by one curve
    explicit gsCurveLoop( gsCurve<T> * cc ) { this->insertCurve(cc); }

    /// Constructor by a list of curves forming a loop
    explicit gsCurveLoop(std::vector< gsCurve<T> *> const & curves ) : m_curves(curves) { }

    /// Constructor from a sequence of points in 3D. Produces a polygon in 2D
    /// whose angles are in proportion to the corresponding angles in space.
    /// You need to pass in a sequence of bools indicating which angles are
    /// convex.
    gsCurveLoop(const std::vector<gsVector3d<T> *> points3D, const std::vector<bool> isConvex, T margin, gsVector3d<T> *outNormal);

    /// Constructor from a sequence of angles and lengths (measured in 3D).
    /// If \a unitSquare==true, when there are four vertices and isConvex are all true, then return the spline curve
    /// loop of degree 1 of the unit square with each edge being one spline curve.
    gsCurveLoop(const std::vector<T> angles3D, const std::vector<T> lengths3D, const std::vector<bool> isConvex, T margin, const bool unitSquare=false);

    gsCurveLoop(const gsCurveLoop & other)
    : m_curves( other.m_curves.size() )
    {
        cloneAll( other.m_curves.begin(), other.m_curves.end(),
                  this->m_curves.begin() );
    }

    gsCurveLoop& operator=(const gsCurveLoop & other)
    {
        freeAll( m_curves );
        m_curves.resize( other.m_curves.size() );
        cloneAll( other.m_curves.begin(), other.m_curves.end(),
                  this->m_curves.begin() );
        return *this;
    }

    ~gsCurveLoop()
    {
        freeAll( m_curves );
    }

    //GISMO_CLONE_FUNCTION(gsCurveLoop)
    uPtr clone() const
    {
        return uPtr(new gsCurveLoop(*this));
    }

public:

    bool isInterior ( gsVector<T> const & p, const T& tol);

    ///@name isOn
    /// returns true if the points in \param u are ON the curve
    /// ideally it should store in \param paramResult the parametric values
    /// associated to those points.
    bool isOn(gsMatrix<T> const &u,T & paramResult, T tol = 1e-3);

    ///@name is_ccw
    ///returns true if the curve is oriented counterclockwise
    bool is_ccw();

    ///@name reverse
    ///changes the orientation of the curve
    void reverse();

    void translate(gsVector<T> const & v);

    void insertCurve( gsCurve<T> * c )
    {
        m_curves.push_back( c ) ;
    }

    /// Return the curve-loop as a single new B-Spline curve
    typename gsCurve<T>::uPtr singleCurve() const;

    int size() const { return m_curves.size(); }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    { os << "gsCurveLoop with "<< size() <<" curves :\n";
        for (typename std::vector< gsCurve<T> *>::const_iterator it= m_curves.begin();
             it!= m_curves.end(); it++ )
        {
            os <<"("<<*it<<") "<< **it ;
        };
        return os;
    }

    /// Computes the normal vector of curve i at parameter u (is
    /// outer in case of CCW orientation)
    gsMatrix<T> normal( int const & c, gsMatrix<T> const & u );

    /// Split a second curve loop off from this one by inserting two new
    /// curves. Curve \a startIndex will be included in the new loop, while
    /// curve \a endIndex will be kept in the original loop.
    uPtr split(int startIndex, int endIndex,
                           gsCurve<T> * newCurveThisFace, gsCurve<T> * newCurveNewFace);


    /// Computes the intersections with the axis-aligned line with
    /// constant value \a abscissa in \a direction.
    // to do: replace return type with std::multimap<int,T>
    std::vector<T> lineIntersections(int const & direction , T const & abscissa);

    std::vector< gsCurve<T> *> & curves() { return m_curves; }

    gsCurve<T> & curve(int i)
    {
        GISMO_ASSERT( i<int(m_curves.size()), "Curve does not exist" );
        return *m_curves[i];
    }

    const gsCurve<T> & curve(int i) const
    {
        GISMO_ASSERT( i<int(m_curves.size()), "Curve does not exist" );
        return *m_curves[i];
    }

    std::vector< gsCurve<T> *> const & curves() const
    { return m_curves; }

    /// get Number of curves
    int numCurves() const     { return m_curves.size(); }

    /// Sample \a npoints uniformly distributed (in parameter domain) points on the loop
    gsMatrix<T> sample(int npoints = 50, int numEndPoints=2) const;

    gsMatrix<T> getBoundingBox();

    /// Computes the corners of a polygon that matches the specified angles and
    /// relative side lengths as well as possible. Fits the resulting polygon
    /// inside the square [0, 1] x [0, 1]. Fills a gsmatrix whose rows are the
    /// corners of the polygon. Returns true on success.
    static bool approximatingPolygon(const std::vector<T> &signedAngles, const std::vector<T> &lengths, T margin, gsMatrix<T> &result);

    /// return a vector containing whose nth entry is the length of the domin of the nth curve.
    std::vector<T> domainSizes() const;

    /// Scale and shift a polygon (represented by a matrix whose rows are its
    /// corners so that it will be contained in the square
    ///   [0, 1]x[0, 1].
    /// Modifies the polygon in place. Specifying margin > 0 will make the
    /// polygon fit in the interior of the square.
    static void adjustPolygonToUnitSquare(gsMatrix<T> &corners, T const margin);

    /// split the \a curveId^th curve in the loop into two curves. Return the point where it splits.
    /// \param curveId
    /// \param lengthRatio    the ratio between the lengths of the first new curve and that of the original curve
    gsMatrix<T> splitCurve(size_t curveId, T lengthRatio=.5);

    /// Initialize a curve loop from some 3d vertices by projecting them onto the plane of best fit and constructing a polygon
    gsVector3d<T> initFrom3DPlaneFit(const std::vector<gsVector3d<T> *> points3D, T margin);

    /// Initialize a curve loop from some 3d vertices by trying to match the angles as well as possible.
    /// Return true on success.
    bool initFrom3DByAngles(const std::vector<gsVector3d<T> *> points3D, const std::vector<bool> isConvex, T margin);

    /// Initialize a curve loop from a list of bools indicating whether the corners are convex or not
    void initFromIsConvex(const std::vector<bool> isConvex, T margin);

    /// flip a curve loop in the u direction
    void flip1(T minu = 0, T maxu = 1);

private:
    bool initFrom3DByAngles(const std::vector<T>& angles3D, const std::vector<T>& lengths3D, const std::vector<bool>& isConvex, T margin, bool unitsquare=false);

    bool parameterOf(gsMatrix<T> const &u, int i, T & result, T tol = 1e-5);

    void removeCurves(iterator begin, iterator end);

    void removeCurve(iterator it)
    { removeCurves(it, it+1); }

    // Data members

    std::vector< gsCurve<T> *> m_curves; // CCW

}; // class gsCurveLoop

template<class T>
bool gsCurveLoop<T>::isInterior ( gsVector<T> const &, const T&)
{
    GISMO_NO_IMPLEMENTATION
}

/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsCurveLoop<T>& b)
{return b.print(os); }


} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCurveLoop.hpp)
#endif
