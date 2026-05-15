/** @file example_cutcell_levelset_occ.cpp

    @brief Local Bernstein level-set functions on cut cells using native OCC.

    This is the OCC-native counterpart of example_cutcell_levelset.cpp from
    the gsGmsh module.  Instead of routing through the gmsh API to access the
    OCC kernel, this example uses the OpenCASCADE C++ API directly:

      - Disk geometry  : gp_Circ + BRepBuilderAPI_MakeEdge/Wire/Face
      - Boundary edge  : TopExp_Explorer (TopAbs_EDGE)
      - Point-in-face  : BRepClass_FaceClassifier (2-D parameter space)
      - Closest point  : BRepExtrema_DistShapeShape (vertex vs. edge)

    A circular disk is defined in native OCC and the parameter domain of a
    tensor B-spline patch is cut by the disk boundary.  Cut cells are
    identified by the same Lobatto-sampling strategy used in
    example_cutcell_levelset.cpp.

    For each cut element a **local Bernstein (Bézier) basis** is constructed:
    a gsTensorBSplineBasis<2> over the element box [u0,u1]×[v0,v1] with zero
    interior knots (i.e., full multiplicity at each end), giving one single
    polynomial element of degree p in each direction.

    The signed distance to the trimming geometry is evaluated at the Greville
    anchor points of the Bernstein basis:
      - Physical coordinates are obtained by evaluating the patch at the
        anchors.
      - BRepExtrema_DistShapeShape returns the nearest point on the OCC
        circle edge and its distance.
      - The sign is determined by BRepClass_FaceClassifier: φ < 0 inside the
        geometry (standard level-set convention).

    The level-set data are then interpolated through the anchors with
    basis.interpolateAtAnchors(vals), yielding a gsTensorBSpline<2> over the
    element box.  Each per-element level-set is written to Paraview.

    Usage:
      ./example_cutcell_levelset_occ [-x nX] [-y nY] [-p degree]

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsIO/gsWriteParaview.h>

#include <cmath>
#include <memory>
#include <vector>

// OCC geometry primitives
#include <gp_Ax2.hxx>
#include <gp_Ax3.hxx>
#include <gp_Circ.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Vec.hxx>

// OCC surface/plane handle
#include <Geom_Plane.hxx>

// OCC topology builders
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>

// OCC topology exploration and classification
#include <BRepClass_FaceClassifier.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <TopAbs_State.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Vertex.hxx>

using namespace gismo;

int main(int argc, char** argv)
{
    index_t nX = 8;
    index_t nY = 8;
    index_t p  = 2;    // Bernstein degree for the local level-set

    gsCmdLine cmd("Cut-cell level-set (OCC-native): local Bernstein approximation of signed distance.");
    cmd.addInt("x", "numX",   "Number of elements in X direction", nX);
    cmd.addInt("y", "numY",   "Number of elements in Y direction", nY);
    cmd.addInt("p", "degree", "Bernstein degree for level-set",    p);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // -------------------------------------------------------------------------
    // 1. Background Cartesian mesh and identity patch
    // -------------------------------------------------------------------------
    gsKnotVector<> kvX(0.0, 1.0, nX - 1, 2, 1, 1);
    gsKnotVector<> kvY(0.0, 1.0, nY - 1, 2, 1, 1);
    gsTensorBSplineBasis<2, real_t> bkgBasis(kvX, kvY);
    auto domain = bkgBasis.domain();

    // Bilinear identity patch [0,1]^2 → [0,1]^2 (parameter space = physical space)
    gsMatrix<real_t> coefs(4, 2);
    coefs << 0,0, 1,0, 0,1, 1,1;
    gsKnotVector<> kv1(0,1,0,2);
    gsTensorBSpline<2, real_t> patch(kv1, kv1, give(coefs));

    // -------------------------------------------------------------------------
    // 2. Define the OCC geometry: disk centred at (0.5,0.5), radius 0.4
    //    Built using native OCC primitives — no gmsh wrapper needed.
    // -------------------------------------------------------------------------
    const double cx     = 0.5;
    const double cy     = 0.5;
    const double rDisk  = 0.4;

    // Full circle edge in the xy-plane
    gp_Circ     circ(gp_Ax2(gp_Pnt(cx, cy, 0.0), gp_Dir(0.0, 0.0, 1.0)), rDisk);
    TopoDS_Edge circleEdge = BRepBuilderAPI_MakeEdge(circ).Edge();
    TopoDS_Wire circleWire = BRepBuilderAPI_MakeWire(circleEdge).Wire();

    // Planar face bounded by the circle (disk interior)
    gp_Pln      diskPlane(gp_Pnt(cx, cy, 0.0), gp_Dir(0.0, 0.0, 1.0));
    TopoDS_Face diskFace  = BRepBuilderAPI_MakeFace(diskPlane, circleWire,
                                                    /*OnlyPlane=*/true).Face();

    // Extract the plane's Ax3 to project 3D points into 2D face-parameter space.
    // For BRepClass_FaceClassifier we need (u,v) coordinates on the plane:
    //   u = dot(P - origin, xDir)
    //   v = dot(P - origin, yDir)
    Handle(Geom_Surface) occSurf  = BRep_Tool::Surface(diskFace);
    Handle(Geom_Plane)   occPlane = Handle(Geom_Plane)::DownCast(occSurf);
    GISMO_ENSURE(!occPlane.IsNull(), "Expected a planar face from MakeFace.");
    const gp_Ax3& planeAx3 = occPlane->Position();

    // Helper: 3D point (in xy-plane) → 2D face parameter coordinates
    auto toUV = [&](double x, double y) -> gp_Pnt2d
    {
        gp_Vec vec(planeAx3.Location(), gp_Pnt(x, y, 0.0));
        return gp_Pnt2d(vec.Dot(gp_Vec(planeAx3.XDirection())),
                        vec.Dot(gp_Vec(planeAx3.YDirection())));
    };

    // Helper: true iff the 2D point (x,y) lies inside or on the disk face
    auto isInsideDisk = [&](double x, double y) -> bool
    {
        BRepClass_FaceClassifier clf(diskFace, toUV(x, y), 1e-7);
        const TopAbs_State s = clf.State();
        return s == TopAbs_IN || s == TopAbs_ON;
    };

    // Helper: closest point on the circle edge + unsigned distance
    BRep_Builder brepBuilder;
    auto closestOnCircle = [&](double x, double y,
                                double& closestX, double& closestY) -> double
    {
        TopoDS_Vertex v;
        brepBuilder.MakeVertex(v, gp_Pnt(x, y, 0.0), 1e-7);
        BRepExtrema_DistShapeShape dist(v, circleEdge);
        GISMO_ENSURE(dist.IsDone(), "BRepExtrema_DistShapeShape computation failed.");
        gp_Pnt cp = dist.PointOnShape2(1);
        closestX  = cp.X();
        closestY  = cp.Y();
        return dist.Value();
    };

    // -------------------------------------------------------------------------
    // 3. Classify elements using Lobatto sample points and build level-sets
    // -------------------------------------------------------------------------
    gsVector<index_t> lobDim(2);
    lobDim << 4, 4;
    gsLobattoRule<real_t> lobRule(lobDim);

    gsMatrix<> samplePts;
    gsVector<> sampleWts;

    std::vector<std::unique_ptr<gsGeometry<real_t>>> levelSets;
    std::vector<gsTensorBSpline<2, real_t>>           bezGeoms;
    std::vector<index_t> cutElemIds;

    index_t nCut = 0;

    for (auto& elem : domain->allElements())
    {
        const auto& lo = elem.lowerCorner();
        const auto& hi = elem.upperCorner();

        lobRule.mapTo(lo, hi, samplePts, sampleWts);

        int nIn = 0;
        for (int i = 0; i < samplePts.cols(); ++i)
            if (isInsideDisk(samplePts(0,i), samplePts(1,i))) ++nIn;

        // Only process cut elements (mixed inside/outside)
        const int nTotal = static_cast<int>(samplePts.cols());
        if (nIn == 0 || nIn == nTotal) continue;

        ++nCut;

        // ----------------------------------------------------------------
        // 4. Build a local Bernstein (Bézier) basis over this element
        //    degree p, zero interior knots → full-multiplicity end knots
        // ----------------------------------------------------------------
        const real_t u0 = lo(0), u1 = hi(0);
        const real_t v0 = lo(1), v1 = hi(1);

        gsKnotVector<real_t> kvU(u0, u1, 0, p + 1);
        gsKnotVector<real_t> kvV(v0, v1, 0, p + 1);
        gsTensorBSplineBasis<2, real_t> bezBasis(kvU, kvV);

        // ----------------------------------------------------------------
        // 5. Evaluate signed distance at the Greville anchor points
        // ----------------------------------------------------------------
        gsMatrix<real_t> anchors    = bezBasis.anchors(); // 2 × nDofs
        const index_t   nDofs       = anchors.cols();

        // Map anchors from parameter space to physical space via the patch
        gsMatrix<real_t> physAnchors;
        patch.eval_into(anchors, physAnchors); // 2 × nDofs (identity patch here)

        gsMatrix<real_t> vals(1, nDofs);

        for (index_t k = 0; k < nDofs; ++k)
        {
            const double x = static_cast<double>(physAnchors(0, k));
            const double y = static_cast<double>(physAnchors(1, k));

            // Closest point on the circle boundary + unsigned distance
            double closestX, closestY;
            const double dist = closestOnCircle(x, y, closestX, closestY);

            // Sign: φ < 0 inside the geometry (standard level-set convention)
            const bool inside = isInsideDisk(x, y);
            vals(0, k) = inside ? -static_cast<real_t>(dist)
                                :  static_cast<real_t>(dist);
        }

        // ----------------------------------------------------------------
        // 6. Interpolate to obtain the local level-set geometry
        // ----------------------------------------------------------------
        auto levelSet = bezBasis.interpolateAtAnchors(vals);

        // Build a Bezier geometry patch for this element (physical positions
        // of Greville anchors — for the identity patch these equal the anchors).
        gsMatrix<real_t> geomCoefs = physAnchors.transpose(); // nDofs × geoDim
        gsTensorBSpline<2, real_t> bezGeom(kvU, kvV, give(geomCoefs));

        levelSets.push_back(std::move(levelSet));
        bezGeoms.push_back(bezGeom);
        cutElemIds.push_back(elem.id());
    }

    // -------------------------------------------------------------------------
    // 7. Write each per-element level-set to Paraview
    // -------------------------------------------------------------------------
    gsInfo << "Cut-cell level-set OCC-native (degree=" << p << "):\n";
    gsInfo << "  Number of cut elements: " << nCut << "\n";

    for (std::size_t i = 0; i < levelSets.size(); ++i)
    {
        const std::string fn = "levelset_occ_elem_" + std::to_string(cutElemIds[i]);
        gsWriteParaview(bezGeoms[i], *levelSets[i], fn, 32);
    }

    gsInfo << "  Wrote " << levelSets.size()
           << " level-set patches (levelset_occ_elem_*.vts)\n";

    // Sanity check: evaluate the level-set at the element centre and compare
    // with exact distance to the circle (r = 0.4, centre (0.5,0.5)).
    if (!levelSets.empty())
    {
        gsInfo << "\nSanity check (first cut element):\n";
        const auto& ls = *levelSets[0];
        gsMatrix<real_t> centre(2, 1);
        centre.col(0) = (ls.basis().support().col(0) +
                         ls.basis().support().col(1)) * 0.5;
        gsMatrix<real_t> lsVal;
        ls.eval_into(centre, lsVal);
        const real_t R  = static_cast<real_t>(rDisk);
        const real_t x  = centre(0,0), y = centre(1,0);
        const real_t exactDist = std::sqrt((x - static_cast<real_t>(cx))*(x - static_cast<real_t>(cx))
                                         + (y - static_cast<real_t>(cy))*(y - static_cast<real_t>(cy))) - R;
        gsInfo << "  Level-set at elem centre = " << lsVal(0,0)
               << "  exact = " << exactDist << "\n";
    }

    return 0;
}
