/** @file gsTrimSurface.hpp

    @brief Provides implementation of gsTrimSurface class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D.-M. Nguyen, M. Pauley, J. Speh
*/


#pragma once

#include <gsUtils/gsMesh/gsMesh.h>
#include <gsModeling/gsModelingUtils.hpp>

#include <gsModeling/gsPlanarDomain.h>

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsKnotVector.h>

namespace gismo
{

template <class T>
std::ostream &gsTrimSurface<T>::print(std::ostream &os) const
{
	os << "Trimmed surface with "<< nTrims()<<" trim loop(s).\n";
	os << "Master surface: "<< *m_surface <<"\n";
	os << "Domain: "<< *m_domain <<"\n";
	return os;
}

template <class T>
gsTrimSurface<T>::gsTrimSurface(gsMatrix<T> const & corner, int patchDeg1, int patchDeg2, int curveDeg)
{
    // Construct the Bezier knot vector for all trimming curves
    gsKnotVector<T> kv(0, 1, 0, curveDeg+1);
    unsigned int ntcp = kv.size() - kv.degree() - 1;

    // Initialize a trimming curve loop
    gsCurveLoop<T> * tloop = new gsCurveLoop<T>();

    // run over the 4 corner vertices of the unit square (the fifth one is the first one)
    gsMatrix<T> DomCor = gsMatrix<T>(5,2);
    DomCor << 0, 0, 1, 0, 1, 1, 0, 1, 0, 0;
    for (int ic=0; ic < 4; ic++)
    {
        // Define a spline curve from the current vertex to the next one
        gsMatrix<T> tcp (ntcp, 2);
        for (unsigned int i=0; i < ntcp; i++)
        {
            for (unsigned int xi=0; xi < 2; xi++)
            {
                tcp(i,xi) = DomCor(ic,xi) + T(i) / (T(ntcp)-1) * (DomCor(ic+1,xi)-DomCor(ic,xi));
            }
        }
        tloop->insertCurve( new gsBSpline<T>(0,1,0,curveDeg, give(tcp)) );
    }

    // Then construct a planar domain with only an outer loop *tloop*
    gsPlanarDomain<T> * domain1= new gsPlanarDomain<T>(tloop);

    // Define the master NURBS surface: degree 2, Bezier patch: todo: derive cps from corners
    gsKnotVector<T> KV1 = gsKnotVector<T>(0, 1, 0, patchDeg1+1);
    gsKnotVector<T> KV2 = gsKnotVector<T>(0, 1, 0, patchDeg2+1);
    typename gsTensorBSpline<2,T>::Ptr tp1(new gsTensorBSpline<2,T>(corner, KV1, KV2));

    this->m_domain = domain1;
    this->m_surface = tp1;
}

template <class T>
gsMatrix<T> gsTrimSurface<T>::derivatives(int sourceID) const
{
    std::vector< gsCurve<T>* >  trimLoop = m_domain->outer().curves();
    //edges of the angle in the parameter domain
    gsMatrix<T> CPside2 = (trimLoop[sourceID]->coefs());
    gsMatrix<T> AngleVertex(2,1); // vertex of the trimmed surface angle
    AngleVertex << CPside2(0,0), CPside2(0,1);
    return m_surface->jacobian(AngleVertex);
}


template <class T>
gsVector3d<T> gsTrimSurface<T>::cornerNormal(int const & sourceID) const
{
    gsMatrix<T> cj = this->derivatives(sourceID);
    gsVector3d<T> dx = cj.col(0);
    gsVector3d<T> dy = cj.col(1);
    gsVector3d<T> n = dx.cross(dy);
    n.normalize();
    return n;
}


template <class T>
void gsTrimSurface<T>::cuttingAngles(int const & sourceID,int const & targetID,T* angle, T* angle1, T* angle2, bool const & isCCWviewFromNormal) const
{
    std::vector< gsCurve<T>* >  trimLoop = m_domain->outer().curves();
    // source
    gsMatrix<T> CPside2 = (trimLoop[sourceID]->coefs());
    gsMatrix<T> AngleVertex(2,1); // vertex of the trimmed surface angle
    AngleVertex << CPside2(0,0), CPside2(0,1);
    // target
    gsMatrix<T> targetCP = (trimLoop[targetID]->coefs());
    gsMatrix<T> TargetVertex(2,1);
    TargetVertex << targetCP(0,0), targetCP(0,1);
    // cutting edge
    gsMatrix<T> cuttingEdge = TargetVertex-AngleVertex;

    //edges of the angle in space
    gsMatrix<T> tangent_side1 = TangentCoefs_prev(sourceID);
    gsMatrix<T> tangent_side2 = TangentCoefs_next(sourceID);
    gsVector3d<T> tangent_side1_space;
    gsVector3d<T> tangent_side2_space;

    gsMatrix<T> corJacobian = derivatives(sourceID);
    tangent_side1_space = corJacobian*tangent_side1;
    tangent_side2_space = corJacobian*tangent_side2;

    gsVector3d<T> cuttingEdge_space = corJacobian*cuttingEdge;
    //
    gsVector3d<T> normal = unitNormal(AngleVertex);
    if (isCCWviewFromNormal==false)
    {
        normal = - normal;
    };

    *angle1 = conditionedAngle<T>(cuttingEdge_space, tangent_side1_space, normal);
    *angle2 = conditionedAngle<T>(tangent_side2_space, cuttingEdge_space, normal);
    *angle  = conditionedAngle<T>(tangent_side2_space, tangent_side1_space, normal);
}

template <class T>
void gsTrimSurface<T>::sampleLoop_into( int loopNumber, int npoints, gsMatrix<T> & u) const
{
    assert( (loopNumber>=0) && (loopNumber < m_domain->numLoops()) );

    gsMatrix<T> pts = m_domain->sampleLoop(loopNumber, npoints);
    return m_surface->eval_into(pts, u); // Compute points on curve of the surface
}

template <class T>
void gsTrimSurface<T>::sampleCurve_into( int loopNumber, int curveNumber, int npoints, gsMatrix<T> & u ) const
{
    assert( (loopNumber>=0) && (loopNumber < m_domain->numLoops()) );
    assert( (curveNumber>=0) && (curveNumber < m_domain->loop(loopNumber).size() ) );

    gsMatrix<T> pts = m_domain->sampleCurve(loopNumber, curveNumber, npoints);
    //m_domain->sample(pts,loopNumber,npoints);
    return m_surface->eval_into(pts, u); // Compute points on curve of the surface
}


template <class T>
gsVector<T> gsTrimSurface<T>::TangentCoefs_bisect(int const & sourceID) const
{
    gsMatrix<T> corJacobian = derivatives(sourceID);
    gsMatrix<T> tangent_side1 = UnitTangentCoefs_prev(sourceID, corJacobian);
    gsMatrix<T> tangent_side2 = UnitTangentCoefs_next(sourceID, corJacobian);
    gsVector<T> coefs(2);
    coefs(0) = tangent_side1(0) + tangent_side2(0);
    coefs(1) = tangent_side1(1) + tangent_side2(1);

    return coefs;
}

template <class T>
gsVector<T> gsTrimSurface<T>::TangentCoefs_bisect(int const & sourceID, gsVector3d<T> normal) const
{
    gsMatrix<T> corJacobian = derivatives(sourceID);
    gsVector3d<T> tangent_side1 = UnitTangentCoefs_prev(sourceID, corJacobian);
    gsVector3d<T> tangent_side2 = UnitTangentCoefs_next(sourceID, corJacobian);
    gsVector<T> coefs(2);
    coefs = TangentCoefs_bisect(sourceID);
    return ( normal.dot( tangent_side1.cross( tangent_side2 ) ) >= 0 ) ? coefs: (-coefs).eval();
}

template <class T>
gsBSpline<T> gsTrimSurface<T>::cuttingCurve(int const & sourceID,int const & targetID) const
{
    std::vector< gsCurve<T>* >  trimLoop = m_domain->outer().curves();
    int curveDeg = trimLoop[sourceID]->degree();
    if (curveDeg<=1)
    {
        gsWarn<<"gsTrimSurface: degree of trimming curve is less than 2, this will fail to work in most cases. The degree is set to 3 instead.\n";
        curveDeg=3;
    }
    gsKnotVector<T> kv(0, 1, 0, curveDeg+1); // todo: get kv from curve00.basis().knots(),curve00.basis().degree()

    gsMatrix<T> CPside2 = (trimLoop[sourceID]->coefs());
    gsMatrix<T> targetCP = (trimLoop[targetID]->coefs());

    gsMatrix<T> tg1 = TangentCoefs_bisect(sourceID); //tangent 1
    gsMatrix<T> tg2 = TangentCoefs_bisect(targetID); // tangent 2

    // weights
    T w_reg = 1;
    T w_app = 0;

    // Exact constraints: point interpolation
    short_t dimPI = 1; // dimension of space of preImage
    short_t dimI = 2;  // dimension of space of image
    int nip=2; // number of interpolating points
    int nn=2; // number of prescribed normals
    gsMatrix<T> image(dimI,nip);
    gsMatrix<T> preImage(dimPI,nip);
    gsMatrix<T> normal(dimI,nn);
    gsMatrix<T> preNormal(dimPI,nn);
    preImage  << kv[0],kv[kv.size()-1];
    preNormal << kv[0],kv[kv.size()-1];
    image.col(0) = CPside2.row(0);
    image.col(1) = targetCP.row(0);
    normal(0,0) = -tg1(1);
    normal(1,0) =  tg1(0);
    normal(0,1) = -tg2(1);
    normal(1,1) =  tg2(0);

    // Approximate constraints
    int nipApp = 1;
    gsMatrix<T> preImageApp(dimPI,nipApp);
    gsMatrix<T> imageApp(dimI,nipApp);
    preImageApp << .5*kv[0]+.5*kv[kv.size()-1];
    imageApp = .5*CPside2.row(0).transpose() + .5*targetCP.row(0).transpose();

    gsMatrix<T> pointResiduals, normalResiduals;
    return gsInterpolate(kv,preImage,image,preNormal,normal,preImageApp,imageApp,w_reg,w_app, pointResiduals, normalResiduals);
}


template <class T>
memory::unique_ptr<gsMesh<T> > gsTrimSurface<T>::toMesh(int npoints) const
{
    typename gsMesh<T>::uPtr msh = m_domain->toMesh(npoints);
    gsMatrix<T> tmp;

    // For all vertices of the msh, push forward the value by m_surface
    for (size_t i = 0; i!= msh->numVertices(); ++i)
    {
        m_surface->eval_into( msh->vertex(i).topRows(2), tmp );
        msh->vertex(i).topRows(m_surface->geoDim() ) = tmp;
    }

    return msh;
}

template <class T>
gsMatrix<T> gsTrimSurface<T>::UnitTangentCoefs_next(int const & sourceID,gsMatrix<T> const & corJacobian) const
{
    gsMatrix<T> tangent_side2 = TangentCoefs_next(sourceID);
    gsMatrix<T> spatialTangent = tangent_side2(0)*corJacobian.col(0)+tangent_side2(1)*corJacobian.col(1);
    tangent_side2 = tangent_side2*(1/spatialTangent.norm());

    return tangent_side2;
}

template <class T>
gsMatrix<T> gsTrimSurface<T>::UnitTangentCoefs_prev(int const & sourceID,gsMatrix<T> const & corJacobian) const
{
    gsMatrix<T> tangent_side1 = TangentCoefs_prev(sourceID);
    gsMatrix<T> spatialTangent = tangent_side1(0)*corJacobian.col(0)+tangent_side1(1)*corJacobian.col(1);
    tangent_side1 = tangent_side1*(1/spatialTangent.norm());

    return tangent_side1;
}

template <class T>
gsMatrix<T> gsTrimSurface<T>::TangentCoefs_next(int const & sourceID) const
{
    std::vector< gsCurve<T>* >  trimLoop = m_domain->outer().curves();
    gsMatrix<T> CPside2 = (trimLoop[sourceID]->coefs());
    gsMatrix<T> tangent_side2(2,1); // side 2 of the trimmed surface angle
    tangent_side2 << CPside2(1,0)-CPside2(0,0),CPside2(1,1)-CPside2(0,1);

    return tangent_side2;
}

template <class T>
gsMatrix<T> gsTrimSurface<T>::TangentCoefs_prev(int const & sourceID) const
{
    std::vector< gsCurve<T>* >  trimLoop = m_domain->outer().curves();
    int priorToSource(0);
    if ( sourceID>0 ) priorToSource = sourceID-1; else priorToSource = trimLoop.size()-1;
    gsMatrix<T> CPside1 = (trimLoop[priorToSource]->coefs());
    int const nrow = CPside1.rows();
    gsMatrix<T> tangent_side1(2,1); // side 1 of the trimmed surface angle
    tangent_side1 << CPside1(nrow-2,0)-CPside1(nrow-1,0),CPside1(nrow-2,1)-CPside1(nrow-1,1);
    return tangent_side1;
}


template <typename T>
T gsTrimSurface<T>::getLengthOfCurve(const int loopNumber,
                                     const int curveNumber,
                                     const T eps,
                                     const int nmbSegments) const
{
    const gsCurve<T>& curve = getCurve(loopNumber, curveNumber);
    gsMatrix<T> support  = curve.support();

    gsVector<T> supStart = support.col(0);
    gsVector<T> supEnd = support.col(1);

    gsMatrix<T> params;

    int n = nmbSegments;
    T length = 1;
    T newLength = 0;

    while (eps < math::abs(length - newLength))
    {
        params = uniformPointGrid(supStart, supEnd, n);
        length = newLength;

        newLength = getLengthOfCurve(curve, params, false);
        n += 10;
    }

    return newLength;
}


// =============================================================================
// private member functions
// =============================================================================

template <class T>
void gsTrimSurface<T>::evalCurve_into(const gsCurve<T>& curve,
                                      const gsMatrix<T>& u,
                                      gsMatrix<T>& result) const
{
    gsMatrix<T> curveResult;
    curve.eval_into(u, curveResult);
    m_surface->eval_into(curveResult, result);
}


template <typename T>
T gsTrimSurface<T>::arcLength(const gsCurve<T>& curve,
                              const T a,
                              const T b,
                              const int quadPoints) const
{
    gsVector<index_t> tmp(1);
    tmp(0) = quadPoints;

    gsGaussRule<T> gaussRule;
    gaussRule.setNodes(tmp);

    gsMatrix<T> nodes;
    gsVector<T> weights;

    gsVector<T> lower(1);
    gsVector<T> upper(1);
    lower(0) = a;
    upper(0) = b;

    gaussRule.mapTo(lower, upper, nodes, weights);

    T length = 0;
    for (int col = 0; col != nodes.cols(); col++)
    {
        gsMatrix<T> param(1, 1);
        param(0, 0) = nodes(0, col);

        gsMatrix<T> grad;
        gsMatrix<T> jacobian;
        gsMatrix<T> pointOnCurve;

        curve.eval_into(param, pointOnCurve);
        curve.deriv_into(param, grad);
        m_surface->jacobian_into(pointOnCurve, jacobian);

        gsMatrix<T> derivative = jacobian * grad;

        T value = 0;
        for (int row = 0; row != derivative.rows(); row++)
        {
            value += derivative(row, 0) * derivative(row, 0);
        }

        value = math::sqrt(value);
        length += weights(col) * value;
    }

    return length;
}


template <typename T>
T gsTrimSurface<T>::getLengthOfCurve(const gsCurve<T>& curve,
                                     gsMatrix<T>& params,
                                     bool linear) const
{
    T length = 0;

    if (linear)
    {
        gsMatrix<T> points;
        evalCurve_into(curve, params, points);

        for (int col = 0; col < points.cols() - 1; col++)
        {
            gsVector<T> p1 = points.col(col);
            gsVector<T> p2 = points.col(col + 1);

            // compute distance between this two points
            T dst = 0;
            for (int row = 0; row < p1.rows(); row++)
            {
                T tmp = (p1(row) - p2(row));
                tmp *= tmp;
                dst += tmp;
            }

            length += math::sqrt(dst);
        }
    }
    else
    {
        for (int col = 0; col != params.cols() - 1; col++)
        {
            const T a = params(0, col);
            const T b = params(0, col + 1);

            T arc = arcLength(curve, a, b);

            length += arc;
        }
    }

    return length;
}


template <typename T>
void gsTrimSurface<T>::fromArcsToParams(const int loopNumber,
                                        const int curveNumber,
                                        const gsVector<T>& arcs,
                                        gsVector<T>& parameters,
                                        const T eps)
{
    const gsCurve<T>& curve = getCurve(loopNumber, curveNumber);

    gsMatrix<T> support  = curve.support();

    gsVector<T> supStart = support.col(0);
    gsVector<T> supEnd = support.col(1);

    gsMatrix<T> par = uniformPointGrid(supStart, supEnd, arcs.size() * 10);

    parameters.resize(arcs.size());


    T curArc = 0; // current arc
    T oldParam = supStart(0); // variable to remember old parameter
    T newParam = oldParam;
    int parCounter = 1; // how many parameters (par!) we already used

    for (int p = 0; p != arcs.size(); p++)
    {
        const T arc = arcs(p);

        // we build curent arc untill we reach arc
        while (curArc < arc)
        {
            oldParam = newParam;
            newParam = par(0, parCounter);
            parCounter++;

            T dst = arcLength(curve, oldParam, newParam);
            curArc += dst;

            if (parCounter == par.cols())
                break;
        }

        // in case we just missed last point -------------------------------
        if (parCounter == par.cols() && p == arcs.size() - 1)
        {
            parameters(p) = supEnd(0);
            break;
        }

        if (parCounter == par.cols() && p != arcs.size() - 1)
        {
            GISMO_ERROR("This should not happen. Last arc is too long...\n");
        } // ---------------------------------------------------------------


        // here we know that arc <= curArc,
        // we want to find parameter t, that arc == curArc
        // we do this with binary search (bisection) technique
        if (arc == curArc)
        {
            parameters(p) = newParam;
        }
        else
        {
            // find parameter via binary search
            T t = findParameter(curve,
                                arc, curArc, oldParam, newParam,
                                eps);
            parameters(p) = t;
        }
    }
}


template <typename T>
T gsTrimSurface<T>::findParameter(const gsCurve<T>& curve,
                                  const T arc,
                                  T curArc,
                                  T lowParam,
                                  T uppParam,
                                  const T eps)
{
    // setting upper arc and lower arc
    T arcUpp = curArc;
    T arcLow = curArc - arcLength(curve, lowParam, uppParam);

    T midParam;

    // do the binary search
    while (eps < math::abs(arcUpp - arcLow))
    {
        midParam = lowParam + (uppParam - lowParam) / 2;

        // "arc" distance
        T dst = arcLength(curve, midParam, uppParam);

        if (arc == arcUpp - dst)
        {
            return midParam;
        }
        else if (arc < arcUpp - dst)
        {
            arcUpp -= dst;
            uppParam = midParam;
        }
        else
        {
            arcLow = arcUpp - dst;
            lowParam = midParam;
        }
    }

    midParam = lowParam + (uppParam - lowParam) / 2;

    return midParam;
}


} // namespace gismo
