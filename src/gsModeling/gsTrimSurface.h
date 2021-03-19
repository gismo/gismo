/** @file gsTrimSurface.h

    @brief Provides interface of gsTrimSurface class. Represents a
    trimmed surface (= spline "master surface" in 3d + a planar domain)

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D.-M. Nguyen, M. Pauley, J. Speh
*/

#pragma once

#include <iostream>
#include <gsCore/gsSurface.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsUtils/gsPointGrid.h>

namespace gismo
{

  /** 
      @brief Class for a trim surface
      
      \ingroup Modeling
  */

  
template<class T>
class gsTrimSurface
{

public:

    typedef memory::shared_ptr<gsTrimSurface> Ptr;
    typedef memory::unique_ptr<gsTrimSurface> uPtr;

    /// Default empty constructor
    gsTrimSurface() : m_domain(NULL) { }
    
    /// Constructor by a shared pointer 
    gsTrimSurface( typename gsSurface<T>::Ptr g, gsPlanarDomain<T> * d) 
    : m_surface(g), m_domain(d) 
    { }

    gsTrimSurface( gsSurface<T>* g, gsPlanarDomain<T> * d) : m_surface(g), m_domain(d) 
    { }
    
    /// Construct a trivial trimmed surface from 4 vertices
    /// Inputs: corner (the 4 corner vertices), patchDeg1 and patchDeg2 (spline degree in dim 1 and 2), curveDeg (spline degree of all trimming curves)
    gsTrimSurface(gsMatrix<T> const & corner, int patchDeg1, int patchDeg2, int curveDeg);

    gsTrimSurface(const gsTrimSurface& other)
    {
        m_surface = other.m_surface;
        m_domain = other.m_domain->clone().release();
    }

    gsTrimSurface& operator=(const gsTrimSurface& other)
    {
        m_surface = other.m_surface;
        delete m_domain;
        m_domain = other.m_domain->clone().release();
        return *this;
    }

    ~gsTrimSurface() //destructor
    {
      // don't delete m_surface since it is now a shared_ptr - it will take care of itself
      // delete m_surface;
      delete m_domain;
    }
    
    /// Clone function. Used to make a copy of the (derived) geometry
    gsTrimSurface<T> * clone() const { return new gsTrimSurface(*this); }
    
    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const;

    friend std::ostream& operator<< (std::ostream& os, const gsTrimSurface& ts)
    { return ts.print(os); }
    
    gsBasis<T> & basis() const { return m_surface->basis(); }

    short_t geoDim() const { return m_surface->geoDim(); }

    //bool isProjective() const { return m_surface->isProjective(); }
    
    typename gsSurface<T>::Ptr getTP() const { return m_surface; }

    int nTrims() const { return m_domain->numLoops(); }
    
    /// Get a boundary curve loop.
    gsCurveLoop<T>& boundaryLoop() { return m_domain->loop(0); }

    /// Look at (non-const) boundaryLoop.
    const gsCurveLoop<T>& boundaryLoop() const { return m_domain->loop(0); }

    /// Returns curveNumber-th curve in loopNumber-th loop.
    gsCurve<T>& getCurve(int loopNumber, int curveNumber) const
    {
        return m_domain->curve(loopNumber, curveNumber);
    }

    void sampleLoop_into(int loopNumber, int npoints, gsMatrix<T> & u) const;

    gsMatrix<T> sampleLoop(int loopNumber, int npoints = 100) const
    {
        gsMatrix<T> u;
        sampleLoop_into(loopNumber, npoints, u);
        return u;
    }
    
    gsMatrix<T> sampleBoundary( int npoints = 100) const
    { return sampleLoop(0, npoints); }

    void sampleCurve_into( int loopNumber, int curveNumber, int npoints, gsMatrix<T> & u) const;

    gsMatrix<T> sampleCurve( int loopNumber, int curveNumber, int npoints = 100) const
    {
        gsMatrix<T> u;
        sampleCurve_into(loopNumber, curveNumber, npoints, u);
        return u;
    }

    /// Evaluates curveNumber-th curve from loopNumber-th loop.
    ///
    /// Curve is evaluated at parameter u. Parameter u is in curve domain and
    /// resulting points are in physical space ( result = surf(curve(u)) ).
    void evalCurve_into(int loopNumber,
                        int curveNumber,
                        const gsMatrix<T>& u,
                        gsMatrix<T>& result) const
    {
        GISMO_ASSERT(0 <= loopNumber && loopNumber < m_domain->numLoops(),
                     "Loop number is out of range!");
        GISMO_ASSERT(0 <= curveNumber &&
                     curveNumber < m_domain->loop(loopNumber).size(),
                     "Curve number is out of range!");

        const gsCurve<T> & curve = getCurve(loopNumber, curveNumber);
        evalCurve_into(curve, u, result);
    }

    /// Look at evalCurve_into
    gsMatrix<T> evalCurve(int loopNumber,
                          int curveNumber,
                          const gsMatrix<T>& u) const
    {
        gsMatrix<T> result;
        evalCurve_into(loopNumber, curveNumber, u, result);
        return result;
    }

    gsMatrix<T> sampleBoundaryCurve( unsigned curveNumber, int npoints = 100) const
    { return sampleCurve(0, curveNumber, npoints); }

    /// Evaluates surface.
    ///
    /// Note that this is NOT trimmed surface. (You can evaluate at parameter
    /// which is trimmed away.)
    ///
    /// \param u parameters
    /// \param result points on surface at parameter u
    void evalSurface_into(const gsMatrix<T>& u,
                     gsMatrix<T>& result) const
    {
        m_surface->eval_into(u, result);
    }

    /// Look at evalSurface_into
    ///
    /// \param u parameters
    ///
    /// \return points on surface at parameter u
    gsMatrix<T> evalSurface(const gsMatrix<T>& u) const
    {
        gsMatrix<T> result;
        evalSurface_into(u, result);
        return result;
    }


    gsPlanarDomain<T> & domain()             { return *m_domain; }
    const gsPlanarDomain<T> & domain() const { return *m_domain; }

    /// split the \a curveId^th curve in the \a loopId^th loop of the planar domain into two curves
    /// \param loopId specifies the loop
    /// \param curveId specifies the curve in the loop
    /// \param lengthRatio   the ratio of the lengths of the first new curve and of the original curve
    gsMatrix<T> splitCurve(size_t loopId, size_t curveId, T lengthRatio=.5)
    {return m_domain->splitCurve(loopId,curveId,lengthRatio);}
    
    /// Compute the partial derivatives of the surface parametrization at a corner point of a trimmed patch 
    /// specified by the index of the trimming curve \a sourceID emanating from the corner point
    gsMatrix<T> derivatives(int sourceID) const;
    
    /// Compute the surface normal at a corner.
    gsVector3d<T> cornerNormal(int const & sourceID) const;           
        
    /// Compute angle discrepancy between the angles defined by the edge *source-target* with the sides of the angle wrt *source*
    void cuttingAngles(int const & sourceID,int const & targetID,T* angle, T* angle1, T* angle2,bool const & isCCWviewFromNormal=true) const;      
    
    /// Compute the (scale of) tangent of the curve emanateing from a vertex in the parameter whose image bisects the corresponding angle
    gsVector<T> TangentCoefs_bisect(int const & sourceID) const;
    
    /// Compute the (scale of) tangent of the curve emanateing from a vertex in the parameter whose image bisects the corresponding angle
    /// In addition, accounting for the bisecting vector defined according to the normal
    gsVector<T> TangentCoefs_bisect(int const & sourceID,gsVector3d<T> normal) const;    
    
    /// Define a spline curve connecting *source-target*
    gsBSpline<T> cuttingCurve(int const & sourceID,int const & targetID) const;          

    /// Return a triangulation of the trimmed surface
    memory::unique_ptr<gsMesh<T> > toMesh(int npoints = 50) const;

    /// Return the coefficients of the representation of the unit tangent of the edge ENAMATING from vertex *sourceID* in terms of the standard basis of the tangent space
    gsMatrix<T> UnitTangentCoefs_next(int const & sourceID,gsMatrix<T> const & corJacobian) const; 
    
    /// Return the coefficients of the representation of the unit tangent of the edge COMING to vertex *sourceID* in terms of the standard basis of the tangent space
    gsMatrix<T> UnitTangentCoefs_prev(int const & sourceID,gsMatrix<T> const & corJacobian) const;
    
      /// Return the coefficients of the representation of the unit tangent of the edge ENAMATING from vertex *sourceID* in terms of the standard basis of the tangent space
    gsMatrix<T> TangentCoefs_next(int const & sourceID) const;   
    
    /// Return the coefficients of the representation of the unit tangent of the edge COMING to vertex *sourceID* in terms of the standard basis of the tangent space
    gsMatrix<T> TangentCoefs_prev(int const & sourceID) const;       
    
    gsMatrix<T> vertexCoord(int const & loopID, int const & curveID) const
    {
      gsMatrix<T> cp = m_domain->curve(loopID,curveID).coefs();
      gsMatrix<T> vert(2,1);
      gsMatrix<T> vert3D;
      vert << cp(0,0), cp(0,1);      
      (*m_surface).eval_into(vert,vert3D);
      return vert3D;
    }
    
    /// Compute the unit normal vector of the trimmed surface at a point in the parameter domain
    gsMatrix<T> unitNormal(gsMatrix<T> point) const
    {
      gsMatrix<T> Jacobian = m_surface->jacobian(point);
      return Jacobian.col(0).template head<3>().cross(
          Jacobian.col(1).template head<3>() ).normalized();
    }
    
    /// sample standard unit normals along a trimming curve
    gsMatrix<T> sampleNormal(int loopNumber, int curveNumber, size_t npoints) const
    {
      assert( (loopNumber>=0) && (loopNumber < m_domain->numLoops()) );
      assert( (curveNumber>=0) && (curveNumber < m_domain->loop(loopNumber).size() ) );
      //gsMatrix<T> u( this->geoDim(), npoints );

      gsMatrix<T> u(3, npoints);
      
      gsMatrix<T> pts = m_domain->sampleCurve(loopNumber, curveNumber, npoints);
      gsMatrix<T> nm(3,1);
      
      for (size_t i=0; i < npoints; i++)
      {
          u.col(i) = unitNormal(pts.col(i));
      }
      
      return u;
    }      
    
    /// Return the tangent vectors of the trimming curve \a curveNumber in trimming loop \a loopNumber    
    gsMatrix<T> trimCurTangents(int loopN, int curveN, size_t npoints) const
    {      
      gsMatrix<T> interval = m_domain->curve(loopN,curveN).parameterRange();
      // sample parameter points 
      gsMatrix<T> tval = gsPointGrid(interval(0,0), interval(0,1), npoints);
      gsMatrix<T> trimCur = m_domain->curve(loopN,curveN).eval(tval);
      gsMatrix<T> trimCurDev = m_domain->curve(loopN,curveN).jacobian(tval);
      gsMatrix<T> trimCurJac = m_surface->jacobian( trimCur );
      //gsMatrix<T> tangents(this->geoDim(),npoints);
      gsMatrix<T> tangents(3,npoints);
      
      for (size_t i=0; i<=npoints-1; i++)
      {
        tangents.col(i) = trimCurJac.middleCols( 2*i,2 )*trimCurDev.col(i);
      }
      return tangents;
    }


    void cleanEndpoints(T eps)
    {
      GISMO_UNUSED(eps);
      gsMatrix<T> supp = this->m_surface->support();
      size_t n = domain().loop(0).curves().size();
      for(size_t i = 0; i < n; i++)
      {
        gsCurve<T> &thisCurve = domain().loop(0).curve(i);
        for(size_t dim = 0; dim < 2; dim++)
        {
          gsMatrix<T> &coefs = thisCurve.coefs();
          GISMO_ASSERT(coefs(0, dim) >= supp(dim, 0) - eps, "invalid curve endpoint");
          coefs(0, dim) = std::max(coefs(0, dim), supp(dim, 0));
          GISMO_ASSERT(coefs(coefs.rows() - 1, dim) <= supp(dim, 1) + eps, "invalid curve endpoint");
          coefs(coefs.rows() - 1, dim) = std::min(coefs(coefs.rows() - 1, dim), supp(dim, 1));
        }
      }
    }

    /// find the parameter of the nearest point on a curve to a given point in space
    T nearestPoint(int loopNumber, int curveNumber, int nTrialPoints, int nIterations, const gsMatrix<T> &spacePoint)
    {
        GISMO_ASSERT(spacePoint.rows() == 3 && spacePoint.cols() == 1, "Invalid dimensions");
        const gsCurve<T> & c = m_domain->curve(loopNumber, curveNumber);
        gsMatrix<T> supp = c.support();
        gsVector<T> start = supp.col(0), end = supp.col(1);
        gsMatrix<T> trialPoints = uniformPointGrid(start, end, nTrialPoints);
        gsMatrix<T> curveVal, curveDeriv, curveDeriv2, surfVal, surfDeriv, surfDeriv2;
        T closestParam(10e100), closestSqDist(10e100);

        for(int idxTrial = 0; idxTrial < nTrialPoints; idxTrial++)
        {
            gsMatrix<T> u = trialPoints.col(idxTrial);
            // apply Newton's method to the function
            // f(u) = ||p - x(u)||^2
            // where x(u) is the parametrisation of the curve in space.
            // (todo - also check the distances at the endpoints of the curve?
            // although this is for splitting the curve and we would never want
            // to split at the endpoint.)
            for(int iteration = 0; iteration < nIterations; iteration++)
            {
                c.eval_into(u, curveVal);
                c.jacobian_into(u, curveDeriv);
                c.deriv2_into(u, curveDeriv2);

                m_surface->eval_into(curveVal, surfVal);
                m_surface->jacobian_into(curveVal, surfDeriv);
                m_surface->deriv2_into(curveVal, surfDeriv2);

                // evaluate derivative of f
                gsMatrix<T> sqDistDeriv = -2 * (spacePoint - surfVal).transpose() *
                        surfDeriv * curveDeriv;
                GISMO_ASSERT(sqDistDeriv.rows() == 1 && sqDistDeriv.cols() == 1, "Derivative should be 1x1");

                // evaluate second derivative of f
                // from comment for gsBasis::deriv2_into, the second deriv of a surface takes the form:
                // ( dxxf_1 dyyf_1 dxyf_1 dxxf_2 dyyf_2 dxy2f_2 dxxf_3 dyyf_3 dxyf_3 )^T

                gsMatrix<T> termFromSurfaceCurvature(3, 1);
                for(int idxJ = 0; idxJ < 3; idxJ++)
                {
                    termFromSurfaceCurvature(idxJ, 0) = surfDeriv2(idxJ * 3) * curveDeriv(0, 0) * curveDeriv(0, 0) +
                            surfDeriv2(idxJ * 3 + 1) * curveDeriv(1, 0) * curveDeriv(1, 0) +
                            2 * surfDeriv2(idxJ * 3 + 2) * curveDeriv(0, 0) * curveDeriv(1, 0);
                }
                gsMatrix<T> sqDistDeriv2Term1 = (surfDeriv * curveDeriv).transpose();
                sqDistDeriv2Term1 *= (surfDeriv * curveDeriv);
                gsMatrix<T> sqDistDeriv2Term2 = (spacePoint - surfVal).transpose() * (termFromSurfaceCurvature + surfDeriv * curveDeriv2);
                gsMatrix<T> sqDistDeriv2 = 2 * (sqDistDeriv2Term1 - sqDistDeriv2Term2);
                GISMO_ASSERT(sqDistDeriv2.rows() == 1 && sqDistDeriv2.cols() == 1, "Second derivative should be 1x1");

                u -= sqDistDeriv / sqDistDeriv2(0, 0);
                u(0, 0) = (u(0, 0) < supp(0, 0))? supp(0, 0): ((u(0, 0) > supp(0, 1)? supp(0, 1): u(0, 0)));
            }
            // compute sqDist for the point found by the last iteration, and compare against the best seen so far
            c.eval_into(u, curveVal);
            m_surface->eval_into(curveVal, surfVal);
            T sqDist = (spacePoint - surfVal).squaredNorm();
            if(idxTrial == 0 || sqDist < closestSqDist)
            {
                closestParam = u(0, 0);
                closestSqDist = sqDist;
            }
        }
        return closestParam;
    }

    /// Computes length of curveNumber-th curve in loopNumber-th loop.
    T getLengthOfCurve(const int loopNumber,
                       const int curveNumber,
                       const T eps = 1e-6,
                       const int nmbSegments = 50) const;


    /// Returns parameters t_i, i = 1, 2, ..., nmbParams, such that
    /// arc between surface( curve (t_i) ) and surface( curve (t_{i + 1}) )
    /// is equal for all i.
    ///
    /// With other words:
    /// Parameters are maped into points that are uniform in physical space.
    ///
    /// curve is the curveNumber-th curve in loopNumber-th loop.
    void getPhysicalyUniformCurveParameters(const int loopNumber,
                                            const int curveNumber,
                                            const int nmbParams,
                                            gsVector<T>& parameters,
                                            const T eps = 1e-6)
    {
        T length = getLengthOfCurve(loopNumber, curveNumber, eps);

        // divide the length in nmbParams equally spaced arcs
        gsVector<T> arcs(nmbParams);
        for (int pt = 0; pt != nmbParams; pt++)
        {
            arcs(pt) = length * (pt * 1.0 / (nmbParams - 1));
        }

        // compute parameters from arcs
        fromArcsToParams(loopNumber, curveNumber, arcs, parameters, eps);

    }

// private member functions
private:
    /// Evaluates curve (in physical space) at parameter u.
    ///
    /// result = surface( curve(u) )
    void evalCurve_into(const gsCurve<T>& curve,
                        const gsMatrix<T>& u,
                        gsMatrix<T> &result) const;


    /// Computes the arc length between
    /// surface( curve(a) ) and surface( curve(b) ).
    ///
    /// Integral is evaluated with the number of quadPoints.
    T arcLength(const gsCurve<T>& curve,
                const T a,
                const T b,
                const int quadPoints = 4) const;

    /// Computes the length of curve: surface( curve(t) ).
    ///
    /// \param curve
    /// \param params precomputed parameters (we evaluate length of the curve
    ///        from the begining till the end of params values)
    /// \param linear true - we do linear approximation
    ///               false - we compute integral of derivative
    ///                       surface( curve(t) )'
    T getLengthOfCurve(const gsCurve<T>& curve,
                       gsMatrix<T>& params,
                       bool linear = false) const;



    /// Let curve be a curveNumber-th curve in loopNumber-th loop.
    /// Let be C( t ) = surface( curve( t ) ), t in [parameter domain]
    /// Function returns parameters t_i such that
    /// C( t_i ) = P_i
    /// and
    /// L( C( 0 ), C( t_i ) ) = arcs( i )
    /// where
    /// L( x, y ) is an arc length of a curve C between x and y.
    ///
    /// We return t_i in a vector parameters.
    void fromArcsToParams(const int loopNumber,
                          const int curveNumber,
                          const gsVector<T>& arcs,
                          gsVector<T>& parameters,
                          const T eps);


    /// Let be C( t ) = surface( curve( t ) ).
    /// Function finds a parameter t0 of the curve C, where C( t0 ) lies arc
    /// away from the C( 0 ).
    ///
    /// Let P1 lies curArc away from C( 0 ) on the curve C.
    /// Let P2 lies less than arc away from C( 0 ) on the curve.
    ///
    /// uppParam corresponds to P1
    /// lowParam corresponts to P2
    T findParameter(const gsCurve<T>& curve,
                    const T arc,
                    T curArc,
                    T lowParam,
                    T uppParam,
                    const T eps);
// Data members
private:

   typename gsSurface<T>::Ptr m_surface;
    
   gsPlanarDomain<T> * m_domain;      

}; // class gsTrimSurface


}// namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTrimSurface.hpp)
#endif


