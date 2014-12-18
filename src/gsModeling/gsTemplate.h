/** @file gsTemplate.h

    @brief Provides definition of gsTemplate class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Falini, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsGeometry.h>
#include <gsCore/gsBoxTopology.h>
#include <gsCore/gsBoundary.h>
#include <gsModeling/gsCurveLoop.h>
#include <gsUtils/gsPointGrid.h>

#include <gsNurbs/gsNurbsCreator.h>

#include <eiquadprog.hpp> // external file


namespace gismo
{
    
/** 
    Class gsTemplate object.

    A template is a structure defining a topological arrangement of a
    collection of parameter domains (that map to physical-domain
    patches).
    Additionally, the template stores the images of the corner points
    on physical space (but not the surface patch itself).


    For ordering the box vertices:

          6---7
         /|  /|
        2---3 |
        | 4-|-5
        |/  |/ 
        0---1
    
    Sides:

        2---3  3---7  7---6  6---2
        |   |  |   |  |   |  |   |
        |   |  |   |  |   |  |   |
        0---1  1---5  5---4  4---0
    
    Bottom/Top:

        0---1  6---7
        |   |  |   |
        |   |  |   |
        4---5  2---3

    Then we describe quads or triangle pairs in
    the right counter-clockwise order:
    
        2---3               3         2--3 
        |   |  becomes     /|   and   | / 
        |   |             / |         |/ 
        0---1            0--1         0 
    
    Triangles 0-1-3 and 0-3-2
    Quad 0-1-3-2
*/    
template<class T>
class gsTemplate  : public gsBoxTopology
{
public: 

    /// Shared pointer for gsTemplate
    typedef gsBoxTopology Base;
    typedef memory::shared_ptr< gsTemplate > Ptr;
    //typedef memory::unique_ptr< gsTemplate > LocalPtr;

public:

    /// Default empty constructor
    gsTemplate() : Base()  { }

  gsTemplate( int const & i ) : Base()  
  { 
    GISMO_ASSERT(i<=2 , "Template with more than one hole not implemented yet.");
    gsMatrix<T> xi(2,4);	
    xi << 1, 0, -1,  0,
          0, 1,  0, -1;
         
    if ( i == 0)
    {
        
      // TO DO: add copy constructor and assignmnet operator in gsPlanarDomain
      m_pdomain.insertHole( new gsCurveLoop<T>(gsNurbsCreator<T>::BSplineFatCircle()) );  
      
     m_pdomain.outer().reverse();
      // gsDebug<<"\n Does your template has the right orientation ? ";
      // gsDebug<< "m_domain:\n"<< m_pdomain;
	//m_pdomain.check() ;

     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(0),0.5*xi.col(0))  );
     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(1),0.5*xi.col(1))  );
     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(2),0.5*xi.col(2))  );
     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(3),0.5*xi.col(3))  );
    addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(0.5*xi.col(0), 0.5*xi.col(1))  );
    addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(0.5*xi.col(1), 0.5*xi.col(2))  );
    addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(0.5*xi.col(2), 0.5*xi.col(3))  );
    addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(0.5*xi.col(3), 0.5*xi.col(0))  );

    }
    else if ( i == 1 )//template with one point removed
    {
       
        
        // TO DO: add copy constructor and assignmnet operator in gsPlanarDomain
      
     m_pdomain = gsPlanarDomain<T>( gsNurbsCreator<T>::BSplineFatCircle() );
      
//     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(0),0.0*xi.col(0))  );
//     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(1),0.0*xi.col(1))  );
//     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(2),0.0*xi.col(2))  );
//     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(3),0.0*xi.col(3))  );
  
     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(0),0.167231*xi.col(0)+0.0*xi.col(1))  );
     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(1),0.167231*xi.col(0)+0.0*xi.col(1))  );
     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(2),0.167231*xi.col(0)+0.0*xi.col(1))  );
     addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(xi.col(3),0.167231*xi.col(0)+0.0*xi.col(1))  );
    }
    else
    {
	gsWarn<< "Empty template.\n";
    }
  }


    /// square template constructed providing the 4 corners Bb matrix (2,4) (x and y coord.)
    gsTemplate(int r, bool square) : Base()
  {

      gsMatrix<T> Bb(2,4);
      Bb<< 0, 1, 0, 1,
           0, 0, 1, 1;


      //construct a square
      //boundaryB= bottom side of the square

      gsBSpline<T>* boundaryB = gsNurbsCreator<T>::BSplineLineSegment(r*Bb.col(0),r*Bb.col(1));

      //boundaryU = upper side of the square
      gsBSpline<T>* boundaryU = gsNurbsCreator<T>::BSplineLineSegment(r*Bb.col(3),r*Bb.col(2));

      //boundaryL = left side of the square
      gsBSpline<T>* boundaryL = gsNurbsCreator<T>::BSplineLineSegment(r*Bb.col(2),r*Bb.col(0));

      //boundaryR = right side of the square
      gsBSpline<T>* boundaryR = gsNurbsCreator<T>::BSplineLineSegment(r*Bb.col(1),r*Bb.col(3));

      //generating points on the four boundaries
      //label 0 : boundaryB
      //label 1 : boundaryU
      //label 2 : boundaryR
      //label 3 : boundaryL

      gsMatrix<T> uniformPoints = gsPointGrid<T>(0,1,102);
      int numCols = uniformPoints.cols();

      gsMatrix<T>storage0(2, numCols),storage1(2,numCols),storage2(2,numCols),storage3(2,numCols);
      storage0.setZero();
      storage0.row(0)=uniformPoints.row(0);

      storage1.setOnes();
      storage1.row(0)=uniformPoints.row(0);

      storage2.setOnes();
      storage2.row(1)=uniformPoints.row(0);

      storage3.setZero();
      storage3.row(1)=uniformPoints.row(0);


      //creating the skeleton

      for(int i=1; i<storage0.cols()-1;++i)
      {
       addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(storage0.col(i),storage1.col(i)));
       addSkeleton( gsNurbsCreator<T>::BSplineLineSegment(storage2.col(i),storage3.col(i)));

      }
        //join boundary curves in one single curve
        std::vector< gsCurve<T> *> BoundaryCurves;
        BoundaryCurves.push_back(boundaryB);
        BoundaryCurves.push_back(boundaryR);
        BoundaryCurves.push_back(boundaryU);
        BoundaryCurves.push_back(boundaryL);

        gsCurveLoop<T>* BoundaryLoop = new gsCurveLoop<T>(BoundaryCurves);
       // BoundaryLoop.singleCurve();
        m_pdomain = gsPlanarDomain<T>(BoundaryLoop);

  }

  /// construct a template for a convex polygon, with a single skeleton curve
  /// between two specified vertices.
  /// \param intervalWidths A vector of widths of the parametrization intervals of each edge
  /// \param verts A gsMatrix whose rows are the vertices of the convex polygon
  /// \param startIdx Row number of the vertex for the bisecting curve to start at
  /// \param endIdx Row number for the end vertex
  /// \param curveProportion middle fraction of the bisecting curve that will actually become part of the skeleton
  gsTemplate(const std::vector<T> & intervalWidths,
             const gsMatrix<T> & verts, size_t startIdx, size_t endIdx,
             T curveProportion) : Base()
  {
    index_t n = verts.rows();

    // create an outer loop of straight lines
    gsCurveLoop<T> *loop = new gsCurveLoop<T>();
    for(index_t i = 0; i < n; i++)
    {
      gsMatrix<T> tcp(2, 2);
      tcp << verts(i, 0), verts(i, 1), verts((i + 1) % n, 0), verts((i + 1) % n, 1);
      gsBSpline<T> * tcurve = new gsBSpline<T>( 0, intervalWidths[i], 0, 1, give(tcp) );
      loop->insertCurve(tcurve);
      // TODO: assert if the vertices do not define a convex polygon.
    }
    m_pdomain.insertHole(loop);


    // compute a cubic through these points. minimize mean squared acceleration subject to
    // fixed initial/final points and normals
    gsMatrix<T> startPoint = verts.row(startIdx);
    gsMatrix<T> endPoint = verts.row(endIdx);
    // bisector must be normal to the difference between the normalized outgoing vectors
    gsMatrix<T> Tang1i = (verts.row((startIdx + n - 1) % n) - startPoint).normalized();
    gsMatrix<T> Tang1o = (verts.row((startIdx + 1) % n) - startPoint).normalized();
    gsMatrix<T> Tang2i = (verts.row((endIdx + n - 1) % n) - endPoint).normalized();
    gsMatrix<T> Tang2o = (verts.row((endIdx + 1) % n) - endPoint).normalized();
    gsMatrix<T> startNormal = Tang1o - Tang1i;
    gsMatrix<T> endNormal = Tang2o - Tang2i;
    gsMatrix<T> normals(2, 2);
    normals.col(0) = startNormal.transpose();
    normals.col(1) = endNormal.transpose();

    // matrices representing linear equality constraints
    gsMatrix<T> constraints(6, 8);
    constraints.setZero();
    constraints.block(0, 0, 2, 2).setIdentity();
    constraints.block(2, 6, 2, 2).setIdentity();
    constraints.block(4, 0, 1, 2) = - startNormal;
    constraints.block(4, 2, 1, 2) = startNormal;
    constraints.block(5, 4, 1, 2) = - endNormal;
    constraints.block(5, 6, 1, 2) = endNormal;

    gsVector<T> constraintsRHS(6);
    constraintsRHS.segment(0, 2) = startPoint.transpose();
    constraintsRHS.segment(2, 2) = endPoint.transpose();
    constraintsRHS.segment(4, 2).setZero();

    // matrices representing linear inequalities. each constraint
    // forces the control points to be inside the convex polygon.
    gsMatrix<T> inequalities(n * 4, 8);
    gsVector<T> ineqRHS(n * 4);
    inequalities.setZero();
    ineqRHS.setZero();
//    inequalities.block(0, 0, 1, 4) << startNormal(0, 1), -startNormal(0, 0), -startNormal(0, 1), startNormal(0, 0);
//    inequalities.block(1, 4, 1, 4) << endNormal(0, 1), -endNormal(0, 0), -endNormal(0, 1), endNormal(0, 0);
    for(int idxE = 0; idxE < n; idxE++)
    {
        gsMatrix<T> edgeDiff = verts.row((idxE + 1) % n) - verts.row(idxE);
        for(int idxCP = 0; idxCP < 4; idxCP++)
        {
            inequalities.block(idxE * 4 + idxCP, 2 * idxCP, 1, 2) << -edgeDiff(0, 1), edgeDiff(0, 0);
            ineqRHS(idxE * 4 + idxCP) = -edgeDiff(0, 1) * verts(idxE, 0) + edgeDiff(0, 0) * verts(idxE, 1);
        }
    }

    // matrix for computing the second differences
    gsMatrix<T> secondDifference(4, 8);
    secondDifference.setZero();
    for(size_t i = 0; i < 4; i++)
    {
      secondDifference(i, i) = 1;
      secondDifference(i, i + 2) = -2;
      secondDifference(i, i + 4) = 1;
    }

    // matrix for computing the quadratic cost in terms of the second differences
    gsMatrix<T> costFromSec(4, 4);
    costFromSec.setIdentity();
    costFromSec(0, 2) = costFromSec(1, 3) = costFromSec(2, 0) = costFromSec(3, 1) = 0.5;

    // matrix for computing the total cost as a function of the original control points
    gsMatrix<T> totalCost = secondDifference.transpose() * costFromSec * secondDifference;

    // minimise quadratic subject to linear equalities and inequalities
    gsVector<T> imageCoefsConcat;
    gsVector<T> zeroCost(8);
    zeroCost.setZero();
    constraints.transposeInPlace();
    inequalities.transposeInPlace();

    constraintsRHS *= -1;
    ineqRHS        *= -1;

    solve_quadprog(totalCost, zeroCost,
                   constraints , constraintsRHS,
                   inequalities, ineqRHS,
                   imageCoefsConcat);
    gsMatrix<T> imageCoefs(2, 4);
    for(size_t i = 0; i < 4; i++)
    {
      imageCoefs.col(i) = imageCoefsConcat.segment(2 * i, 2);
    }
    GISMO_ASSERT((imageCoefs.col(0) - startPoint.transpose()).norm() < 0.0001 &&
                 (imageCoefs.col(imageCoefs.cols() - 1) - endPoint.transpose()).norm() < 0.0001,
                 "Quadratic optimization failed to satisfy constraints");

    // transform the coefficients to restrict to a subset of the original interval
    // g: matrix to transform basis functions to monomials
    gsMatrix<T> g(4, 4);
    g << 1,-3,3,-1,0,3,-6,3,0,0,3,-3,0,0,0,1;
    // k: matrix to transform monomials in t into monomials in u, where t = k0 + k1 u
    T k0 = (T(1) - curveProportion) / 2;
    T k1 = curveProportion;
    gsMatrix<T> k(4, 4);
    k.setZero();
    k(0, 0) = 1;
    for(index_t i = 1; i < 4; i++) // first row done already
    {
      for(index_t j = 0; j < 4; j++)
      {
        k(i, j) = k0 * k(i - 1, j);
        if(j > 0) k(i, j) += k1 * k(i - 1, j - 1);
      }
    }
    // g^(-1): transforms monomials in u into basis functions in u
    gsMatrix<T> transformedCoefs = imageCoefs * g * k * g.inverse();

    // TODO: if any control points are outside the polygon, project them to the nearest point on the boundary. that would
    // ensure that the spline stays inside the polygon.

    // construct the spline from the coefficients.
    transformedCoefs.transposeInPlace();
    addSkeleton(new gsBSpline<T>(0, 1, 0, 3, give(transformedCoefs)));
  }
    
    ~gsTemplate() 
      { 
          freeAll(m_skeleton);
      }
    
    /// Clone function. Used to make a copy of the object
    gsTemplate * clone() const 
        {
            return new gsTemplate(*this);
        };
    
public:

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
        { 
            if ( this->size() > 0 )
            {
                os << "gsTemplate (" << this->size() << ").\n"; 
                //os << "First patch: " << this->patch(0) << ").\n"; 
            }
            else
                os << "gsTemplate ( empty! ).\n"; 
            return os; 
        }

    int parDim() const { return dim; }   

    /// Number of patches    
    int nPatches() const {return this->size() ; } 

    void addPatch( gsMatrix<T>* c) 
    { 
        this->m_corners.push_back( c ) ;        
    }

    void addSkeleton(  gsBSpline<T> * bs) 
    { 
        m_skeleton.push_back( bs ) ;        
    }

    std::vector< gsBSpline<T> * > skeleton() 
    { 
        return m_skeleton;        
    }
    
    gsBSpline<T> * skeleton(std::size_t const & i) 
    { 
        return m_skeleton[i];        
    }

    // to do: add outer()
    gsPlanarDomain<T> const & domain() const
    {
      return m_pdomain;
    }
    
    const gsCurveLoop<T> & outer()
    {
      return m_pdomain.outer();
    }
    
    const gsCurveLoop<T> & loop(int i) const
    {
      return m_pdomain.loop(i);
    }
    
    int skeletonSize() 
    { 
        return m_skeleton.size();        
    }
// Data members
private:

    /// List of corner points
    gsMatrix<T> * m_corners;

    std::vector< gsBSpline<T> * >  m_skeleton; // skeleton represents "segments" inside the template
    
    gsPlanarDomain<T> m_pdomain;
    
    /// 
    gsMatrix<T> ponctuals;
    
public:
    // Needed since m_pdomain is 16B aligned
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

}; // class gsTemplate


//////////////////////////////////////////////////
//////////////////////////////////////////////////

/// Print (as string) a tamplate
template<class T>
std::ostream &operator<<(std::ostream &os, const gsTemplate<T>& b)
{return b.print(os); }

}; // namespace gismo
