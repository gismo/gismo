/** @file gsCurveLoop.hpp

    @brief Implementation of gsCurveLoop class.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Falini, A. Mantzaflaris, D.-M. Nguyen, M. Pauley
*/

#pragma once

# include <gsModeling/gsCurveLoop.h>
# include <gsNurbs/gsKnotVector.h>
# include <gsNurbs/gsBSpline.h>
# include <gsNurbs/gsBSplineSolver.h>
# include <gsModeling/gsModelingUtils.hpp>

namespace gismo
{
  
template<class T>
gsCurveLoop<T>::gsCurveLoop(const std::vector<gsVector3d<T> *> points3D, const std::vector<bool> , T margin, gsVector3d<T> *outNormal)
{
    // attempt to choose a curve loop by using the plane of best fit
    gsVector3d<T> resultNormal = initFrom3DPlaneFit(points3D, margin);
    if(outNormal != NULL) *outNormal = resultNormal;
}

template<class T>
gsCurveLoop<T>::gsCurveLoop(const std::vector<T> angles3D, const std::vector<T> lengths3D, const std::vector<bool> isConvex, T margin, const bool unitSquare)
{
    this->initFrom3DByAngles(angles3D, lengths3D, isConvex, margin, unitSquare);
}
      
template<class T>
bool gsCurveLoop<T>::parameterOf(gsMatrix<T> const &u, int i, T & result, T tol)
{
    gsBSplineSolver<T> slv;
    gsMatrix<T> e;
    gsBSpline<T> * c = dynamic_cast<gsBSpline<T> *>(m_curves[i]);
     
    if ( c !=0 )
    {
        std::vector<T> roots;
        slv.allRoots(*c, roots, 0, u(0,0) ) ;

        if(roots.size()==0)
            return false;

        gsAsMatrix<T> xx(roots); 
        m_curves[i]->eval_into( xx, e ); 
    
        for(index_t r=0; r!=e.cols(); r++)
        {     
            //gsDebug<<"matrix e (1,r): "<< e(1,r) <<"\n";
            if( math::abs(e( 1, r )-u(1,0))< tol )
            {
                result = roots[r];
                return true;
            }
        }
        return false;
    }
    else
    {
        gsWarn <<"isOnCurve only works for B-splines.\n";
        return false;
    }
}
  
template<class T>
bool gsCurveLoop<T>::isOn(gsMatrix<T> const &u, T & paramResult, T tol )
{
    for(int i=0; i<this->numCurves(); i++)
        if( m_curves[i]->isOn(u, tol) )
        {
            this->parameterOf(u,i,paramResult,tol);
            return true;
          }
            return false;
}

template<class T>
bool gsCurveLoop<T>::is_ccw()
{
    T result(0);
    
    for ( typename std::vector< gsCurve<T> *>::const_iterator it=
              m_curves.begin(); it!= m_curves.end() ; it++ )
    {
        index_t nrows = (*it)->coefs().rows();
        for(index_t i=0; i!=nrows-1;++i)
            result +=(*it)->coefs()(i,0)* (*it)->coefs()(i+1,1)- 
                (*it)->coefs()(i,1)* (*it)->coefs()(i+1,0);
    }
    return ( result > 0 );
}

template<class T>
void gsCurveLoop<T>::reverse()
{
    for ( typename std::vector< gsCurve<T> *>::const_iterator it=
              m_curves.begin(); it!= m_curves.end() ; it++ )
        (*it)->reverse();
    gsCurve<T>* swap;
    size_t size=m_curves.size();
    for ( size_t i=0;i<size/2;i++)
    {
        swap = m_curves[i];
        m_curves[i]=m_curves[size-i-1];
        m_curves[size-i-1]=swap;
    }   
}
  
template<class T>
void gsCurveLoop<T>::translate(gsVector<T> const & v)
{
    for ( iterator it = m_curves.begin(); it != m_curves.end(); ++it)
        (*it)->translate(v);
}

template<class T>
typename gsCurve<T>::uPtr gsCurveLoop<T>::singleCurve() const
{
    typename std::vector< gsCurve<T> *>::const_iterator it;
    
    GISMO_ASSERT( m_curves.size(), "CurveLoop is empty.\n");
    GISMO_ASSERT( m_curves.front(), "Something went wrong. Invalid pointer in gsCurveLoop member.\n");
    
    typename gsCurve<T>::uPtr loop = m_curves.front()->clone();

    for ( it= m_curves.begin()+1; it!= m_curves.end() ; it++ )
    {
        loop->merge( *it ) ;
    }
    return loop;
}

template<class T>
gsMatrix<T> gsCurveLoop<T>::normal( int const & c, gsMatrix<T> const & u )
{
    int n = u.cols();
    gsMatrix<T> result(2,n);
    
    for ( int i=0; i<n; ++i )
    {
        gsMatrix<T> b = m_curves[c]->deriv( u.col(i) ) ;
        b.normalize(); // unit tangent vector
        result(0,i) = b(1);
        result(1,i) = -b(0);
    }
    return result;
}

template<class T>
typename gsCurveLoop<T>::uPtr gsCurveLoop<T>::split(int startIndex, int endIndex,
                                       gsCurve<T> * newCurveThisFace, 
                                       gsCurve<T> * newCurveNewFace)
{
    int n = m_curves.size();
    uPtr result(new gsCurveLoop<T>(newCurveNewFace));
    for(int i = startIndex; i != endIndex; i = (i + 1) % n)
    {
        result->insertCurve(m_curves[i]);
    }
    if(startIndex < endIndex)
    {
        m_curves.erase(m_curves.begin() + startIndex, m_curves.begin() + endIndex);
        if(startIndex == 0)
        {
            m_curves.push_back(newCurveThisFace);
        }
        else
        {
            m_curves.insert(m_curves.begin() + startIndex, newCurveThisFace);
        }
    }
    else
    {
        m_curves.erase(m_curves.begin() + startIndex, m_curves.end());
        m_curves.erase(m_curves.begin(), m_curves.begin() + endIndex);
        m_curves.push_back(newCurveThisFace);
    }
    return result;
}

template <class T>
std::vector<T> gsCurveLoop<T>::lineIntersections(int const & direction , T const & abscissa)
{
    // For now only Curve loops with ONE curve are supported
    gsBSplineSolver<T> slv;
    typename gsBSpline<T>::uPtr c = memory::convert_ptr<gsBSpline<T> >(singleCurve());
    // = dynamic_cast<gsBSpline<T> *>(m_curves[0])) 
    
    if ( c )
    {
        std::vector<T> result;
        slv.allRoots(*c, result, direction, abscissa) ;
        /*slv.allRoots(*c, result, direction, abscissa,1e-7,1000) ;*/
        return result;
    }
    gsWarn<<"Could not get intersection for this type of curve!\n";
    return std::vector<T>();
}

template <class T>
gsMatrix<T> gsCurveLoop<T>::sample(int npoints, int numEndPoints) const
{
    GISMO_ASSERT(npoints>=2, "");
    GISMO_ASSERT(numEndPoints>=0 && numEndPoints<=2, "");
    int np; // new number of points
    switch (numEndPoints)
    {
    case (0):
        np = npoints-2;
        break;
    case (1):
        np = npoints-1;
        break;	
    case (2):
        np = npoints;
        break;
    default:
        np = 0;
        break;
    }   
    
    gsMatrix<T> u(2, m_curves.size() * np);
    int i=0;
    gsMatrix<T> interval;
    gsMatrix<T> pts(1,np);
    gsMatrix<T> uCols(2,np);
    T a,b;
    int firstInd=0;
    int secondInd=np-1;
    if (numEndPoints==0) {firstInd=1;secondInd=npoints-2;}
    typename std::vector< gsCurve<T> *>::const_iterator it;
    for ( it= m_curves.begin(); it!= m_curves.end() ; it++ )
    {
        interval = (*it)->parameterRange();
        a = interval(0,0);
        b = interval(0,1);
        for (int ii=firstInd;ii<=secondInd;ii++) pts(0,ii-firstInd)= a + T(ii)*(b-a)/(npoints-1);
        (*it)->eval_into( pts, uCols );
        u.middleCols( i * np,np ) = uCols;
        i++ ;
    };
    return u;
}
  
template <class T>
gsMatrix<T> gsCurveLoop<T>::getBoundingBox()
{
    gsMatrix<T, 2, 2> result;
    int numCurves = m_curves.size();
    T offset = 0.0001;
    assert(numCurves != 0); // bounding box does not exist if there are no curves
    gsMatrix<T> *coefs = &(m_curves[0]->coefs());
    result(0, 0) = coefs->col(0).minCoeff() - offset;
    result(1, 0) = coefs->col(1).minCoeff() - offset;
    result(0, 1) = coefs->col(0).maxCoeff() + offset;
    result(1, 1) = coefs->col(1).maxCoeff() + offset;
    
    for(int i = 1; i < numCurves; i++)
    {
        coefs = &(m_curves[i]->coefs());
        result(0, 0) = math::min(result(0, 0), coefs->col(0).minCoeff()-offset );
        result(1, 0) = math::min(result(1, 0), coefs->col(1).minCoeff()-offset );
        result(0, 1) = math::max(result(0, 1), coefs->col(0).maxCoeff()+offset );
        result(1, 1) = math::max(result(1, 1), coefs->col(1).maxCoeff()+offset );
    }
    return result;
}
  
template<class T>
bool gsCurveLoop<T>::approximatingPolygon(const std::vector<T> &signedAngles, const std::vector<T> &lengths, T margin, gsMatrix<T> &result)
{
    size_t n = signedAngles.size();
    assert(lengths.size() == n);
    
    std::vector<T> scaledAngles(signedAngles); // copy
    T totalAngle(0);
    for(size_t i = 0; i < n; i++)
    {
        totalAngle += scaledAngles[i];
    }
    if(totalAngle <= 0)
    {
        gsWarn << "Treatment of non-positive total turning angle is not implemented.\n";
        return false;
    }
    // Scale the turning angles so they add to 2 * pi. These will be the
    // angles we use at each vertex in the domain.
    T angleScale = T(2.0 * EIGEN_PI / totalAngle);
    for(size_t i = 0; i < n; i++)
    {
        scaledAngles[i] *= angleScale;
        if(math::abs(scaledAngles[i]) >= EIGEN_PI)
        {
            gsWarn << "Scaled turning angle exceeded pi, treatment of this has not been implemented.\n";
            return false;
        }
    }
    
    // Compute the cumulative turning angle (i.e. the total angle turned
    // by the time we reach a given vertex).
    std::vector<T> cumulativeAngles(n);
    cumulativeAngles[0] = 0;
    for(size_t i = 1; i < n; i++)
    {
        cumulativeAngles[i] = cumulativeAngles[i - 1] + scaledAngles[i];
    }
    
    // Build a matrix to provide the coefficients of the following
    // constraints (RHS of equations is built later):
    //   total change of x coordinate (as we go round the polygon) == 0;
    //   total change of y coordinate (as we go round the polygon) == 0;
    gsMatrix<T> constraintCoefs(2, n);
    for(size_t j = 0; j < n; j++)
    {
        constraintCoefs(0, j) = math::cos(cumulativeAngles[j]);
        constraintCoefs(1, j) = math::sin(cumulativeAngles[j]);
    }
    
    // Find a critical point of the following cost function subject to the above constraints:
    //   (1 / 2) * ( sum of (  (length i) / (original length i)  ) ^ 2 - 1 ).
    gsMatrix<T> invL(n, n);
    invL.setZero();
    for(size_t i = 0; i < n; i++)
    {
        invL(i, i) = 1 / lengths[i];
    }
    gsMatrix<T> ones(n, 1);
    ones.setOnes();
    gsMatrix<T> constraintRhs(2, 1);
    constraintRhs.setZero();
    
    gsMatrix<T> resultLengths = optQuadratic(invL, ones, constraintCoefs, constraintRhs);
    
    for(size_t j = 0; j < n; j++)
    {
        if(resultLengths(j, 0) <= 0)
        {
            gsWarn << "Could not compute satisfactory lengths for approximating polygon.\n";
            return false;
        }
    }
    
    // Work out the positions of the corners of this polygon.
    // (Corner n == corner 0). Initially, we put the first vertex at (0,0). Then
    // we apply a scale & shift to fit it all inside the box [0,1]x[0,1].
    result.resize(n, 2);
    result(0, 0) = 0; // set up the first corner
    result(0, 1) = 0;
    for(size_t i = 1; i < n; i++) // set up corners other than the first
    {
        result(i, 0) = result(i - 1, 0) + resultLengths(i - 1) * math::cos(cumulativeAngles[i - 1]);
        result(i, 1) = result(i - 1, 1) + resultLengths(i - 1) * math::sin(cumulativeAngles[i - 1]);
    }
    adjustPolygonToUnitSquare(result, margin);
    return true;
}
  
template<class T>
std::vector<T> gsCurveLoop<T>::domainSizes() const
{
    std::vector<T> result;
    size_t n = m_curves.size();
    for(size_t i = 0; i < n; i++)
    {
        gsMatrix<T> range = m_curves[i]->parameterRange();
        GISMO_ASSERT(range.rows() == 1 && range.cols() == 2, "Expected 1x2 matrix for parameter range of curve");
        result.push_back(range(0, 1) - range(0, 0));
    }
    
    return result;
}
  
template <class T>
void gsCurveLoop<T>::adjustPolygonToUnitSquare(gsMatrix<T> &corners, T const margin)
{
    assert(corners.cols() == 2);
    size_t n = corners.rows();
    
    T minx = corners.col(0).minCoeff();
    T maxx = corners.col(0).maxCoeff();
    T miny = corners.col(1).minCoeff();
    T maxy = corners.col(1).maxCoeff();
    
    T scaleFactor(1);
    scaleFactor += margin * 2;
    scaleFactor = 1 / scaleFactor / std::max(maxx - minx, maxy - miny);
    gsMatrix<T> shiftVector(1, 2);
    shiftVector << 0.5 - 0.5 * scaleFactor * (maxx + minx), 0.5 - 0.5 * scaleFactor * (maxy + miny);
    for(size_t i = 0; i < n; i++)
    {
        corners.row(i) *= scaleFactor;
        corners.row(i) += shiftVector;
    }
}

template<class T>
gsVector3d<T> gsCurveLoop<T>::initFrom3DPlaneFit(const std::vector<gsVector3d<T> *> points3D, T margin)
{
    // attempt to construct a curve loop using the plane of best fit

    size_t n = points3D.size(); // number of corners
    GISMO_ASSERT(n >= 3, "Must have at least 3 points to construct a polygon");

    freeAll(m_curves);
    m_curves.clear();

    // find the center of mass of the points
    gsVector3d<T> com;
    com << 0, 0, 0;
    for(size_t i = 0; i < n; i++)
    {
        com += *(points3D[i]);
    }
    com /= n;

    // construct a matrix with the shifted points
    gsMatrix<T> shiftedPoints(3, n);
    for(size_t i = 0; i < n; i++)
    {
        shiftedPoints.col(i) = *(points3D[i]) - com;
    }

    // compute a singular value decomposition
    Eigen::JacobiSVD< Eigen::Matrix<T, Dynamic, Dynamic> > svd(
        shiftedPoints * shiftedPoints.transpose(), Eigen::ComputeFullU);

    // extract plane projection matrix from the svd
    gsMatrix<T> svd_u(svd.matrixU());
    GISMO_ASSERT(svd_u.rows() == 3 && svd_u.cols() == 3, "Unexpected svd matrix result");
    gsMatrix<T> projMat = svd_u.block(0, 0, 3, 2).transpose();

    // project all the points down onto the plane, keep track of min, max of coords
    std::vector< gsVector<T> > projPts;
    gsMatrix<T> bbox(2, 2); // first column is min of each dimension, second column is max
    for(size_t i = 0; i < n; i++)
    {
        gsVector<T> x = projMat * shiftedPoints.col(i);
        projPts.push_back(x);
        for(size_t d = 0; d < 2; d++)
        {
            if(i == 0 || x(d) < bbox(d, 0)) bbox(d, 0) = x(d);
            if(i == 0 || x(d) > bbox(d, 1)) bbox(d, 1) = x(d);
        }
    }

    // scale and shift so everything fits inside the bounding box minus margin.
    T scale = (1 - 2 * margin) / std::max((bbox(0, 1) - bbox(0, 0)), bbox(1, 1) - bbox(1, 0));
    gsVector<T> shift(2);
    shift << margin - scale * bbox(0, 0), margin - scale * bbox(1, 0);
    for(size_t i = 0; i < n; i++)
    {
        projPts[i] = projPts[i] * scale + shift;
    }

    // interpolate linearly
    for(size_t i = 0; i < n; i++)
    {
        size_t i1 = (i + 1) % n;
        gsMatrix<T> tcp(2, 2);
        tcp << projPts[i](0), projPts[i](1), projPts[i1](0), projPts[i1](1);
        this->insertCurve(new gsBSpline<T>( 0, 1, 0, 1, give(tcp) ));
    }

    // return the normal to the plane of best fit. we could just grab the last
    // column from svd_u but then we'd have to worry about the orientation.
    gsVector3d<T> p1 = projMat.row(0);
    gsVector3d<T> p2 = projMat.row(1);
    return p1.cross(p2).normalized().transpose();

}

template<class T>
bool gsCurveLoop<T>::initFrom3DByAngles(const std::vector<gsVector3d<T> *> points3D, const std::vector<bool> isConvex, T margin)
{
    size_t n = isConvex.size();
    assert(n == points3D.size());

    std::vector<T> angles;
    std::vector<T> lengths;

    for(size_t i = 0; i < n; i++)
    {
        gsVector3d<T> *prev = points3D[(i + n - 1) % n];
        gsVector3d<T> *cur = points3D[i];
        gsVector3d<T> *next = points3D[(i + 1) % n];

        gsVector3d<T> change0 = *cur - *prev;
        gsVector3d<T> change1 = *next - *cur;
        T length0 = change0.norm();
        T length1 = change1.norm();
        T acosvalue= change0.dot(change1) / length0 / length1;
        if(acosvalue>1)
            acosvalue=1;
        else if(acosvalue<-1)
            acosvalue=-1;
        angles.push_back( math::acos( acosvalue));
        lengths.push_back(length1);
    }

    return initFrom3DByAngles(angles, lengths, isConvex, margin);
}
  
template<class T>
void gsCurveLoop<T>::initFromIsConvex(const std::vector<bool> isConvex, T margin)
{
    // find the corners of a regular polygon
    size_t np = isConvex.size();
    gsMatrix<T> corners(np, 2);
    for(size_t i = 0; i < np; i++)
    {
        T angle = (T)i * (T)EIGEN_PI * 2 / np;
        corners(i, 0) = math::cos(angle);
        corners(i, 1) = math::sin(angle);
    }
    // choose control points of cubic splines which ensure the correct angle signs
    gsMatrix<T> cps(np * 4, 2);
    for(size_t i = 0; i < np; i++)
    {
        size_t iprev = (i + np - 1) % np;
        gsMatrix<T> prevPt = corners.row(iprev);
        gsMatrix<T> u = corners.row(i);
        gsMatrix<T> nextPt = corners.row((i + 1) % np);
        cps.row(iprev * 4 + 3) = u;
        cps.row(i * 4) = u;
        if(isConvex[i])
        {
            cps.row(iprev * 4 + 2) = u + (prevPt - u) / 3;
            cps.row(i * 4 + 1) = u + (nextPt - u) / 3;
        }
        else
        {
            cps.row(iprev * 4 + 2) = u - (nextPt - u) / 3;
            cps.row(i * 4 + 1) = u - (prevPt - u) / 3;
        }
    }
    // put all the control points inside a box
    adjustPolygonToUnitSquare(cps, margin);
    // construct the splines
    curves().clear();
    for(size_t i = 0; i < np; i++)
    {
        gsMatrix<T> tcp = cps.block(i * 4, 0, 4, 2);
        this->insertCurve(new gsBSpline<T>( 0, 1, 0, 3, give(tcp) ));
    }
}

  
template<class T>
bool gsCurveLoop<T>::initFrom3DByAngles(const std::vector<T>& angles3D, const std::vector<T>& lengths3D, const std::vector<bool>& isConvex, T margin, bool unitSquare)
{
    size_t n = isConvex.size();
    assert(angles3D.size() == n);
    assert(lengths3D.size() == n);

    freeAll(m_curves);
    m_curves.clear();

    gsMatrix<T> corners4(4, 2);
    if (unitSquare==true && n==4)
    {
        bool isAllConvex=true;
        for (size_t i=0;i<4;i++) {if (isConvex.at(i)==false) {isAllConvex=false;break;}}
        if (isAllConvex==true)
        {
            corners4<<0,0,1,0,1,1,0,1;
            for(size_t i = 0; i < n; i++)
            {
                gsMatrix<T> tcp(2, 2);
                tcp << corners4(i, 0), corners4(i, 1), corners4((i + 1) % n, 0), corners4((i + 1) % n, 1);
                this->insertCurve(new gsBSpline<T>( 0, 1, 0, 1, give(tcp) ));
            }
            return true;
        }        
    }    
    
    // compute the 2D angles corresponding to the 3D ones
    std::vector<T> signedAngles(n);
    for(size_t i = 0; i < n; i++)
    {
        signedAngles[i] = (isConvex[i]? 1: -1) * angles3D[i];
    }
    
    gsMatrix<T> corners;
    bool success = gsCurveLoop<T>::approximatingPolygon(signedAngles, lengths3D, margin, corners);
    if(!success) return false;
    
    // create a loop of B-splines that are all straight lines.
    for(size_t i = 0; i < n; i++)
    {
        gsMatrix<T> tcp(2, 2);
        tcp << corners(i, 0), corners(i, 1), corners((i + 1) % n, 0), corners((i + 1) % n, 1);
        this->insertCurve(new gsBSpline<T>( 0, 1, 0, 1, give(tcp) ));
    }
    return true;
}

template<class T>
void gsCurveLoop<T>::flip1(T minu, T maxu)
{
    size_t nCur = m_curves.size();
    T offset = minu + maxu;
    for(size_t i = 0; i < nCur; i++)
    {
        gsMatrix<T> &coefs = m_curves[i]->coefs();
        size_t ncp = coefs.rows();
        for(size_t j = 0; j < ncp; j++)
        {
            coefs(j, 0) = offset - coefs(j, 0);
        }
    }
}
  
template<class T>
gsMatrix<T> gsCurveLoop<T>::splitCurve(size_t curveId, T lengthRatio)
{
    GISMO_UNUSED(lengthRatio);
    GISMO_ASSERT(lengthRatio>0 && lengthRatio<1, "the second parameter *lengthRatio* must be between 0 and 1");
    GISMO_ASSERT(curveId < m_curves.size(), "the first parameter *curveID* exceeds the number of curves of the loop");
    // the two vertices a and b of the curve curveId, counting counterclockwise from z+
    gsMatrix<T> a(2,1),b(2,1),n(2,1);
    gsMatrix<T> cp,tcp0(2,2),tcp1(2,2);
    cp = curve(curveId).coefs();
    //int deg = curve(curveId).degree();
    int deg = curve(curveId).basis().degree(0);
    //GISMO_ASSERT(deg==1,"Pls extend to code to curves of degree more than 1");
    if (deg!=1)
        gsWarn<<"Pls extend to code to curves of degree more than 1, degree 1 is temporarily used.\n";
    a = cp.row(0);
    b = cp.row(cp.rows() - 1);
    //n = (1-lengthRatio)*a + lengthRatio*b; // new vertex
    gsMatrix<T> supp = curve(curveId).support();
    GISMO_ASSERT(supp.rows() == 1 && supp.cols() == 2, "Unexpected support matrix");
    gsMatrix<T> u(1, 1);
    u << (supp(0, 0) + supp(0, 1)) / 2;
    curve(curveId).eval_into (u, n);
    n.transposeInPlace();
    // create the two new curves
    tcp0 << a(0,0),a(0,1),n(0,0),n(0,1);
    gsBSpline<T> * tcurve0 = new gsBSpline<T>( 0, 1, 0, 1, give(tcp0) );
    tcp1 << n(0,0),n(0,1),b(0,0),b(0,1);
    gsBSpline<T> * tcurve1 = new gsBSpline<T>( 0, 1, 0, 1, give(tcp1) );
    // remove the original curve
    removeCurve(m_curves.begin()+curveId);
    // add the two new curves
    m_curves.insert(m_curves.begin()+curveId,tcurve0);
    m_curves.insert(m_curves.begin()+curveId+1,tcurve1);

    return n;
}

template <class T>
void gsCurveLoop<T>::removeCurves(iterator begin, iterator end)
{
    // free the curves before removing them
    freeAll(begin, end);
    m_curves.erase(begin, end);
}

}
