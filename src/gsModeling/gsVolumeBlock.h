/** @file gsVolumeBlock.h

    @brief Provides gsVolumeBlock class.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Speh
*/


#pragma once

#include <gsModeling/gsSolid.h>
#include <gsModeling/gsSolidElement.h>
#include <gsModeling/gsFitting.h>

#include <gsCore/gsVolume.h>
#include <gsCore/gsBoundary.h>



namespace gismo {

template <class T>
class gsVolumeBlock : public gsSolidElement<T>
{
public:
    // to do change all types in gsModeling

    typedef gsSolidElement<T> SolidElement;
    typedef typename SolidElement::scalar_t scalar_t;
    typedef gsSolidHalfFace<T> HalfFace;
    typedef gsSolidHalfEdge<T> HalfEdge;
    typedef gsSolidHeVertex<T> Vertex;

public:
    gsVolumeBlock() : SolidElement() { }

    explicit gsVolumeBlock(int const& i) : SolidElement(i) { }
    
    explicit gsVolumeBlock(std::vector<HalfFace*> const & hfaces)
    {
        face = hfaces;
        typename std::vector<HalfFace*>::const_iterator it;
        for (it = hfaces.begin();it!=hfaces.end();++it)
        {
            (*it)->vol = this;
        }
    }

    virtual ~gsVolumeBlock() 
    { }

public:

    /// Function returns the boxSide that corresponds to given face
    ///
    /// This function only works if the VolumeBlock is hexahedron
    ///
    /// \param faceId id of the face
    /// \return side of the hexahedron
    boxSide getSideOfHexahedron(const unsigned faceId)
    {
        if (!isHexahedron())
        {
            GISMO_ERROR("VolumeBlock is not hexahedron.\n");
        }

        HalfEdge* edge = face[0]->loop[0];
        std::vector<HalfEdge* > hex; // hexahedrons initial HalfEdges
        getHexEdges(hex, edge); // each side is presented with leading edge

        HalfFace* myFace = face[faceId];

        int sideIndex = -1; // resulting index of side - 1

        for (size_t edgeIndx = 0; edgeIndx < hex.size(); edgeIndx++)
        {
            HalfEdge* faceEdge = hex[edgeIndx];

            // checks if the edge e is in the face
            HalfEdge* e = myFace->loop[0];
            for (unsigned i = 0; i < 4; i++)
            {
                if (e == faceEdge)
                {
                    sideIndex = static_cast<int>(edgeIndx);
                }

                e = e->next;
            }
        }

        if (sideIndex == -1)
        {
            gsWarn << "Didn't transform face with id = " << faceId << "\n";
        }

        return boxSide(sideIndex + 1);
    }



    /// Function checks if volume presents hexahedron.
    ///
    /// Volume is valid hexahedron if it has 6 faces and each face has
    /// 4 edges.
    ///
    /// This function ignores inner holes.
    bool isHexahedron() const
    {

        if (face.size() != 6)
        {
            return false;
        }
        else
        {
            for (unsigned idFace = 0; idFace != face.size(); idFace++)
            {
                HalfFace* myFace = face[idFace];

                if (1 < myFace->loopN())
                {
                    gsWarn << "Input volume has a hole in its face. "
                              "Ignoring it...\n";
                }


                gsPlanarDomain<T>& domain = myFace->surf->domain();

                // first loop (loop(0)) describes outer boundary
                gsCurveLoop<T>& curveLoop = domain.loop(0);

                if (curveLoop.numCurves() != 4)
                    return false;
            }

            return true;
        }
    }

    /// Function computes volume parameterization of gsVolumeBlock
    gsTensorBSpline<3, T> getTensorBSplineVolume(int numSamplePoints,
                                                 int numInternalKnots,
                                                 int degree,
                                                 T lambda,
                                                 T eps)
    {
        if (!isHexahedron())
        {
            GISMO_ERROR("VolumeBlock is not hexahedron.\n"
                        "Function getVolume can not tranform VolumeBlock into "
                        "gsVolume if VolumeBlock is not hexahedron.\n");

        }

        gsMatrix<T> points;
        gsMatrix<T> params;
        //getUniformBoundaryPoints(numSamplePoints, params, points);
        getUniformPoints(numSamplePoints, params, points, eps);

        gsKnotVector<T> kv1(0, 1, numInternalKnots, degree + 1, 1);
        gsKnotVector<T> kv2(0, 1, numInternalKnots, degree + 1, 1);
        gsKnotVector<T> kv3(0, 1, numInternalKnots, degree + 1, 1);

        gsTensorBSplineBasis<3, T> basis(kv1, kv2, kv3);

        gsDebugVar(kv1.detail());
        gsDebugVar(basis);
        gsFitting<T> fitting(params, points, basis);
        fitting.compute(lambda);

        // maybee delete
        std::vector<T> errors;
        fitting.get_Error(errors, 1);
        gsDebug << "Maximum error is: " << *std::max_element(errors.begin(),
                                                               errors.end())
                  << std::endl;

        // note gsFitting deletes the geomety in destructor
        gsGeometry<T>* geo = fitting.result();

        gsTensorBSpline<3, T>* bspline = dynamic_cast<gsTensorBSpline<3, T>* > (geo);

        if (bspline == NULL)
        {
            GISMO_ERROR("Problems with casting gsGeometry to gsTensorBspline");
        }

        return *bspline;
    }


    /// Returns the points that are uniformly spread on the whole volume.
    ///
    /// \param n number of output points on a adge (a volume then has n^3 points)
    /// \param params output parameters
    /// \param points output uniform points
    void getUniformPoints(const int n,
                          gsMatrix<T>& params,
                          gsMatrix<T>& points,
                          T eps)
    {
        gsMatrix<T> bParams, bPoints;
        getUniformBoundaryPoints(n, bParams, bPoints, eps);

        typedef std::pair<std::pair<T, T>, T> tuple;

        // make a map for fast access
        std::map<tuple, int> map;
        for (int col = 0; col < bParams.cols(); col++)
        {
            tuple  t = makeTuple(bParams(0, col), bParams(1, col), bParams(2, col));
            map[t] = col;
        }

        // save points at corners because they are used often
        const gsVector<T> p000 = bPoints.col(map[makeTuple(0., 0., 0.)]);
        const gsVector<T> p001 = bPoints.col(map[makeTuple(0., 0., 1.)]);
        const gsVector<T> p010 = bPoints.col(map[makeTuple(0., 1., 0.)]);
        const gsVector<T> p011 = bPoints.col(map[makeTuple(0., 1., 1.)]);
        const gsVector<T> p100 = bPoints.col(map[makeTuple(1., 0., 0.)]);
        const gsVector<T> p101 = bPoints.col(map[makeTuple(1., 0., 1.)]);
        const gsVector<T> p110 = bPoints.col(map[makeTuple(1., 1., 0.)]);
        const gsVector<T> p111 = bPoints.col(map[makeTuple(1., 1., 1.)]);

        const int numPoints = (n - 2) * (n - 2) * (n - 2) + bParams.cols();
        params.resize(3, numPoints);
        points.resize(3, numPoints);

        // Coons patch algorithm
        int column = 0;
        for (int indx = 1; indx != n - 1; indx++)
        {
            const T x = (indx * 1.0) / (n - 1);

            for (int indy = 1; indy != n - 1; indy++)
            {
                const T y = (indy * 1.0) / (n - 1);

                for (int indz = 1; indz != n - 1; indz++)
                {
                    const T z = (indz * 1.0) / (n - 1);

                    // parameters
                    params(0, column) = x;
                    params(1, column) = y;
                    params(2, column) = z;

                    // points 3D Coons Patch Algorithm

                    // interpolation of faces
                    const gsVector<T> A =
                            (1 - x) * bPoints.col(map[makeTuple(0, y, z)]) +
                            x * bPoints.col(map[makeTuple(1, y, z)]);

                    const gsVector<T> B =
                            (1 - y) * bPoints.col(map[makeTuple(x, 0, z)]) +
                            y * bPoints.col(map[makeTuple(x, 1, z)]);

                    const gsVector<T> C =
                            (1 - z) * bPoints.col(map[makeTuple(x, y, 0)]) +
                            z * bPoints.col(map[makeTuple(x, y, 1)]);

                    // interpolation of boundary curves
                    const gsVector<T> AB =
                            (1 - x) * (1 - y) * bPoints.col(map[makeTuple(0, 0, z)]) +
                            (1 - x) * y * bPoints.col(map[makeTuple(0, 1, z)]) +
                            x * (1 - y) * bPoints.col(map[makeTuple(1, 0, z)]) +
                            x * y * bPoints.col(map[makeTuple(1, 1, z)]);

                    const gsVector<T> BC =
                            (1 - y) * (1 - z) * bPoints.col(map[makeTuple(x, 0, 0)]) +
                            (1 - y) * z * bPoints.col(map[makeTuple(x, 0, 1)]) +
                            y * (1 - z) * bPoints.col(map[makeTuple(x, 1, 0)]) +
                            y * z * bPoints.col(map[makeTuple(x, 1, 1)]);

                    const gsVector<T> AC =
                            (1 - x) * (1 - z) * bPoints.col(map[makeTuple(0, y, 0)]) +
                            (1 - x) * z * bPoints.col(map[makeTuple(0, y, 1)]) +
                            x * (1 - z) * bPoints.col(map[makeTuple(1, y, 0)]) +
                            x * z * bPoints.col(map[makeTuple(1, y, 1)]);

                    // interpolation of corner points
                    const gsVector<T> ABC =
                            (1 - x) * ((1 - y) * ((1 - z) * p000 + z * p001) +
                                            y  * ((1 - z) * p010 + z * p011)) +
                                  x * ((1 - y) * ((1 - z) * p100 + z * p101) +
                                            y  * ((1 - z) * p110 + z * p111));

                    points.col(column) = A + B + C - AB - BC - AC + ABC;

                    column++;
                }
            }
        }


        // append boundary points to
        for (int col = 0; col < bPoints.cols(); ++col)
        {
            points.col(column) = bPoints.col(col);
            params.col(column) = bParams.col(col);
            column++;
        }
    }


    /// Returns points that are uniformly spread on the boundary. For each edge
    /// we get double points (from each curve)
    ///
    /// \param n number of points in one curve (boundary)
    void getUniformBoundaryPoints(const int n,
                                  gsMatrix<T>& params3D,
                                  gsMatrix<T>& points,
                                  T eps = 1e-6)
    {
        HalfEdge* edge = face[0]->loop[0];
        std::vector<HalfEdge* > hex; // hexahedrons initial HalfEdges
        getHexEdges(hex, edge);

        int nmbOfPts = 6 * (n * n + 4); // 6 faces * number of points on faces
        points  .resize(3, nmbOfPts);
        params3D.resize(3, nmbOfPts);

        gsVector<T> params1D(n);
        for (int index = 0; index != n; index++)
        {
            params1D(index) = index * 1.0 / (n - 1);
        }

        int column = 0;

        for (boxSide side = boxSide::getFirst(3); side<boxSide::getEnd(3); ++side )
        {

            // this variable is true if we must turn around a curve loop in
            // this side
            bool turnAround = turnCurveLoopAround(side);
            T fixedConstant = side.parameter() ? 1.0 : 0.0;

            // fixed - index of coordinate where parameter for this side has
            //         fixedConstant value
            // first - index of coordinate where first curve in hexahedron has
            //         its range
            // second - index of cooridnate where first curve in hexahedron is
            //         constant
            int fixed = side.direction();
            int first = (fixed == 0) ? 1 : 0;
            int second = (fixed == 2) ? 1 : 2;

            if (fixed == 0) // east and west
            {
                int t = first;
                first = second;
                second = t;
            }

//            gsDebug << "side: " << side << "\n"
//                      << "directions side: " << fixed << "\n"
//                      << "parameter: " << fixedConstant << "\n"
//                      << "first: " << first << "\n"
//                      << "seconde: " << second << "\n" << std::endl;

            std::vector<gsVector<T> > faceBoundaryParams;
            HalfEdge* faceEdge = hex[side - 1];



            // compute uniform points on each curve
            for (int curveId = 0; curveId < 4; curveId++)
            {
                HalfEdge* edge2 = faceEdge->moveAlongEdge(curveId);

                // gets index of a curve in Trimmed surface
                int index = edge2->face->indexOfEdge(edge2);

                gsVector<T> params;
                edge2->face->surf->getPhysicalyUniformCurveParameters(
                            0, index, n, params, eps);


                gsMatrix<T> curvePoints;
                edge2->face->surf->evalCurve_into(
                            0, index, params.transpose(), curvePoints);


                bool turnCurve = turnCurveAround(turnAround, curveId);

                // curve has constant parameter along two dimensions
                // (in parameter domain)
                // - first constant is set by position of the face
                // - second constant is set by position of the curve in
                //   curve loop
                T curveConstant = 0.0;
                if ( (turnAround && (curveId == 2 || curveId == 3)) ||
                     (!turnAround && (curveId == 1 || curveId == 2)) )
                {
                    curveConstant = 1.0;
                }

                // curveConstantIndex - index of parameter domain where curve
                //                      is constant
                // curveNonConstIndex - index of parameter domain where curve
                //                      is not constant
                int curveConstantIndex = first;
                int curveNonConstIndex = second;
                if (curveId == 0 || curveId == 2)
                {
                    curveConstantIndex = second;
                    curveNonConstIndex = first;
                }


                for (int col = 0; col < n; col++)
                {
                    params3D(fixed, column) = fixedConstant;
                    params3D(curveConstantIndex, column) = curveConstant;
                    params3D(curveNonConstIndex, column) = params1D(col);

                    if (turnCurve)
                    {
                        points.col(column) = curvePoints.col(n - 1 - col);
                    }
                    else
                    {
                        points.col(column) = curvePoints.col(col);
                    }
                    column++;


//                    gsDebug << "side: " << side << "  curveID: " << curveId
//                              << "  turnCurve: " << turnCurve
//                              << "  turnAround: " << turnAround << "\n"
//                              << "params: " << params3D.col(column - 1).transpose()
//                              << "\n"
//                              << "points: " << points.col(column - 1).transpose()
//                              << "\n" << std::endl;
                }

                if (turnCurve)
                {
                    params.reverseInPlace();
                }

                faceBoundaryParams.push_back(params);
            }

            // compute points in the interior of face

            // sometimes (depends on the orientation of the curve loop
            // we must exchange 1 and 3 curve in coons patch algorithm
            bool switch1and3curve = false;
            if (side == boundary::east ||
                side == boundary::front ||
                side == boundary::north)
            {
                switch1and3curve = true;
            }

            // first we compute parameter for surface via conns patch algorithm
            // for given 4 boundary curves

            gsMatrix<T> cpParams; // coons patch parameters
            modifiedCoonsPatch(hex[side - 1], faceBoundaryParams,
                    params1D, cpParams, switch1and3curve);

            gsMatrix<T> cpPoints;
            faceEdge->face->surf->evalSurface_into(cpParams, cpPoints);

            // set global parameters and points
            for (int indv = 0; indv != n - 2; indv++)
            {
                for (int indu = 0; indu != n - 2; indu++)
                {
                    params3D(fixed, column) = fixedConstant;
                    params3D(first, column) = params1D(indu + 1);
                    params3D(second, column) = params1D(indv + 1);

                    points.col(column) = cpPoints.col(indv * (n - 2) + indu);
                    column++;
                }
            }

        }
    }


private:
    static std::pair<std::pair<T, T>, T> makeTuple(T x, T y, T z)
    {
        return std::make_pair(std::make_pair(x, y), z);
    }

//    static void printTuple(std::pair<std::pair<T, T>, T> tuple)
//    {
//        gsDebug << tuple.first.first << " " << tuple.first.second << " "
//                  << tuple.second;
//    }

    /// Computes parameters in parameter domain via coons patch algorithm.
    ///
    /// There are 4 curves, 1st curve is defined with edge parameter, next
    /// curve is defined with edge->next, ...
    ///
    /// Each curve has stored curve parameters in vector boundaryParams.
    ///
    /// Vector params1D is vector of desired parameters (we want to
    /// reparameterize the curves having uniform parameters, so we use this
    /// new uniform parameters, but values are computed with old parameters).
    ///
    /// switch1and3curve - in coons patch algorithm, we must switch first and
    /// third curve (because curve loop is oriented in opposite direction
    void  modifiedCoonsPatch(HalfEdge* edge,
                             const std::vector<gsVector<T> >& boundaryParams,
                             const gsVector<T>& params1D,
                             gsMatrix<T>& cpPoints,
                             const bool switch1and3curve)
    {
        int n = boundaryParams[0].size();
        int m = boundaryParams[1].size();

        cpPoints.resize(2, (n - 2) * (m - 2));

        // eval all curves...
        std::vector<gsMatrix<T> > values;
        for (int curveId = 0; curveId != 4; curveId++)
        {
            HalfEdge* curveEdge = edge->moveAlongEdge(curveId);
            int internalIndex = getInternalIndexOfCurve(curveEdge);
            gsCurve<T>& curve = edge->face->surf->getCurve(0, internalIndex);

            gsMatrix<T> curvePoints;
            curve.eval_into(boundaryParams[curveId].transpose(), curvePoints);
            values.push_back(curvePoints);
        }

        // save boundary points
        gsVector<T> c00 = values[0].col(0);
        gsVector<T> c01 = values[0].col(n - 1);

        gsVector<T> c20 = values[2].col(0);
        gsVector<T> c21 = values[2].col(n - 1);

        // do the coons patch

        int column = 0;
        // indv is index in "v direction", similar indu
        for (int indv = 1; indv != n - 1; indv++) // ommit boundary
        {
            T v = params1D[indv];


            for (int indu = 1; indu != m - 1; indu++)
            {
                T u = params1D[indu];
                cpPoints.col(column) = (1 - v) * values[0].col(indu) +
                        v * values[2].col(indu)
                        - (1 - v) * (1 - u) * c00
                        - (1 - v) * u * c01
                        - v * (1 - u) * c20
                        - v * u * c21;

                gsVector<T> vec;
                if (switch1and3curve)
                {
                    vec = (1 - u) * values[3].col(indv) +
                          u * values[1].col(indv);
                }
                else
                {
                    vec = (1 - u) * values[1].col(indv) +
                          u * values[3].col(indv);
                }

                cpPoints.col(column) += vec;

//                gsDebug << "u: " << u << "    v: " << v << "\n"
//                          << "indu: " << indu << "    indv: " << indv << "\n"
//                          << "change1and2: " << switch1and3 << "\n";

//                gsMatrix<T> mat(2, 1);
//                mat(0, 0) = cpPoints(0, column);
//                mat(1, 0) = cpPoints(1, column);
//                gsMatrix<T> res;
//                edge->face->surf->evalSurface_into(mat, res);

//                gsDebug << "coons patch: " << res.transpose() << std::endl;

//                mat(0, 0) = values[0](0, indu);
//                mat(1, 0) = values[0](1, indu);
//                edge->face->surf->evalSurface_into(mat, res);

//                gsDebug << "curve0(indu): " << res.transpose() << std::endl;

//                mat(0, 0) = values[2](0, indu);
//                mat(1, 0) = values[2](1, indu);
//                edge->face->surf->evalSurface_into(mat, res);

//                gsDebug << "curve2(indu): " << res.transpose() << std::endl;

//                mat(0, 0) = values[1](0, indv);
//                mat(1, 0) = values[1](1, indv);
//                edge->face->surf->evalSurface_into(mat, res);

//                gsDebug << "curve1(indv): " << res.transpose() << std::endl;

//                mat(0, 0) = values[3](0, indv);
//                mat(1, 0) = values[3](1, indv);
//                edge->face->surf->evalSurface_into(mat, res);

//                gsDebug << "curve3(indv): " << res.transpose() << std::endl;

                column++;
            }
        }
    }


    /// Returns true if the curve is reversed for computation of volume.
    ///
    /// \param curveLoopAround if curve loop is reversed
    /// \param curveId nnumber of curve in curve loop (according to hexahedron)
    static
    bool turnCurveAround(const bool curveLoopAround, const int curveId)
    {
        if ( (!curveLoopAround && (curveId == 2 || curveId == 3)) ||
             (curveLoopAround && (curveId == 0 || curveId == 3))     )
        {
            return true;
        }
        else
        {
            return false;
        }
    }


    /// Returns true if curve Loop is reversed (if it is eather
    /// west, south or back)
    ///
    /// \param side boundary side
    static
    bool turnCurveLoopAround(boxSide side)
    {
        if (side == boundary::west ||
            side == boundary::south ||
            side == boundary::back)
        {
            return true;
        }
        else
        {
            return false;
        }
    }




    /// Returns index of curve in CurveLoop in trimmed surface.
    int getInternalIndexOfCurve(HalfEdge* edge) const
    {
        return edge->face->indexOfEdge(edge);
    }

    /// Hexahedron is defined with 6 faces:
    /// west, east, south, north, front, back
    /// First edge in the face is
    ///  * low edge if face is west, east, front, back
    ///  * front edge if face is south, north
    void getHexEdges(std::vector<HalfEdge*>& hex, HalfEdge* edge)
    {
        // set the west face as the first face of the volume
        // and set the first edge of the west face as the first edge of the
        // first face
        hex.push_back(edge);

        hex.push_back(getEastEdge(edge));

        hex.push_back(getSouthEdge(edge));

        hex.push_back(getNorthEdge(edge));

        hex.push_back(getFrontEdge(edge));

        hex.push_back(getBackEdge(edge));

    }

    /// Compute first edge on the east face from the first edge on the west face
    HalfEdge* getEastEdge(HalfEdge* edge)
    {
        edge = edge->next;
        edge = edge->mate;
        edge = edge->moveAlongEdge(2);
        edge = edge->mate;
        edge = edge->next;
        return edge;
    }

    /// Compute first edge on the south face from the first edge on the west face
    HalfEdge* getSouthEdge(HalfEdge* edge)
    {
        edge = edge->mate;
        edge = edge->prev;
        return edge;
    }

    /// Compute first edge on the north face from the first edge on the west face
    HalfEdge* getNorthEdge(HalfEdge* edge)
    {
        edge = edge->moveAlongEdge(2);
        edge = edge->mate;
        edge = edge->next;
        return edge;
    }

    /// Compute first edge on the front face from the first edge on the west face
    HalfEdge* getFrontEdge(HalfEdge* edge)
    {
        edge = edge->next;
        edge = edge->mate;
        edge = edge->next;
        return edge;
    }

    /// Compute first edge on the back face from the first edge on the west face
    HalfEdge* getBackEdge(HalfEdge* edge)
    {
        edge = edge->prev;
        edge = edge->mate;
        edge = edge->prev;
        return edge;
    }



public:
    std::vector<HalfFace*>  face;
};

} // namespace gismo
