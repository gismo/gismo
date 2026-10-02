/** @file calculate_volume_immersed_example.cpp

    @brief Volume of a 3D triangle mesh (cow.obj) via the Finite Cell Method.

    Pipeline (3D analogue of poisson2_nitsche_immersed_example.cpp):

      1. Read a triangle mesh from a Wavefront .obj file (the cow).
      2. Rescale/center the mesh into the unit parametric box [0,1]^3.
      3. Build an analytic level set phi(x,y,z) from the triangle soup:
           - distance : minimal point-to-triangle distance to the surface.
           - sign     : generalized (solid-angle) winding number; phi < 0
                        inside, phi > 0 outside (the convention used by
                        gsImplicitTrimmedDomain).
      4. Classify the background cells of [0,1]^3 with
         gsImplicitTrimmedDomain<3> (interior / cut / exterior).
      5. Integrate the constant 1 over the immersed domain using the Algoim
         quadrature rule on cut cells (phi < 0 volumetric rule) and a full
         Gauss rule on interior cells.  The result is the volume of the cow
         in the parametric box; dividing by the rescaling factor recovers the
         physical volume.
      6. (--plot) Export the rescaled cow surface together with the Algoim
         volumetric quadrature points to ParaView.

    Usage:
    ./bin/calculate_volume_immersed_example [-f cow.obj] [-o output_cow3d] [-r 3] [-e 0] [--fill 0.9] [--plot]
    ./build/bin/calculate_volume_immersed_example -f obj/spot.obj -o output_spot -r 3 --plot

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;


// =============================================================================
//  main
// =============================================================================
int main(int argc, char * argv[])
{
    std::string filename =
        "obj/spot.obj";
    std::string out = gsFileManager::getCanonicRepresentation(
        gsFileManager::getPath(std::string(__FILE__), true) + "../output_cow3d");
    index_t numRefine  = 3;
    index_t numElevate = 0;
    real_t  fill       = 0.9;  // fraction of [0,1]^3 spanned by the largest extent
    bool    plot       = false;

    gsCmdLine cmd("Volume of a 3D triangle mesh from an .obj file via the Finite Cell Method.");
    cmd.addString("f", "file",            "Wavefront .obj mesh file",            filename);
    cmd.addString("o", "output",          "Output folder for ParaView files",    out);
    cmd.addInt   ("r", "uniformRefine",   "Number of uniform refinement steps",  numRefine);
    cmd.addInt   ("e", "degreeElevation", "Number of degree elevation steps",    numElevate);
    cmd.addReal  ("",  "fill",            "Fill fraction of [0,1]^3 (0..1)",      fill);
    cmd.addSwitch("plot", "Create ParaView output",                              plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    out = gsFileManager::getCanonicRepresentation(out);

    // -------------------------------------------------------------------------
    // 1. Load the triangle mesh
    // -------------------------------------------------------------------------
    gsSurfMesh mesh;
    if (!gsReadSurfMesh(filename, mesh))
    {
        gsWarn << "Failed to read a triangle mesh from: " << filename << "\n";
        return EXIT_FAILURE;
    }
    const std::string outputStem = gsFileManager::getBasename(filename);
    const std::size_t nVert = mesh.n_vertices();
    const std::size_t nTri  = mesh.n_faces();
    gsInfo << "Loaded mesh '" << filename << "': "
           << nVert << " vertices, " << nTri << " triangles.\n";

    // -------------------------------------------------------------------------
    // 2. Rescale / center the mesh into the unit parametric box [0,1]^3
    //    mapped = (p - center) * scale + 0.5,  scale = fill / maxExtent
    // -------------------------------------------------------------------------
    const gsMatrix<real_t> physBox = gsSurfMeshBoundingBox(mesh);
    gsInfo << "Physical bounding box: [" << physBox.col(0).transpose()
           << "] - [" << physBox.col(1).transpose() << "]\n";

    const real_t scale = gsNormalizeToUnitBox(mesh, fill);
    gsInfo << "Rescaling factor (param/phys): " << scale << "\n";

    // -------------------------------------------------------------------------
    // 3. Analytic level set on the rescaled mesh
    // -------------------------------------------------------------------------
    gsMatrix<real_t> bbox(3, 2);
    bbox.col(0).setZero();
    bbox.col(1).setOnes();
    gsMeshSignedDist<real_t> impl_fun(mesh, bbox);

    // -------------------------------------------------------------------------
    // 4. Background box [0,1]^3 (identity map: parameter space == physical space)
    // -------------------------------------------------------------------------
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineCube((real_t)1));

    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);

    // -------------------------------------------------------------------------
    // 5. Output directory
    // -------------------------------------------------------------------------
    gsFileManager::mkdir(out);

    gsInfo << "\nDegree of basis: " << dbasis.maxCwiseDegree() << "\n";
    gsInfo << "Computing the volume of the immersed domain (cow) in [0,1]^3.\n\n";
    gsInfo << std::setw(8) << "refine" << std::setw(10) << "cells/dim"
           << std::setw(18) << "vol (param)" << std::setw(18) << "vol (phys)" << "\n";

    const double invScale3 = 1.0 / (scale * scale * scale);

    // Number of uniform sample points per direction used as fallback when Algoim
    // returns 0 nodes (re-entrant / multi-crossing cut cells).  n^3 = 64 pts.
    const index_t nFallback = 4;

    // -------------------------------------------------------------------------
    // 6. Refinement loop: classify cells and integrate 1 over the cow
    // -------------------------------------------------------------------------
    for (int r = 0; r <= numRefine; ++r)
    {
        dbasis.uniformRefine();

        gsTensorBSplineBasis<3, real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<3, real_t> *>(&dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Expected a tensor B-spline basis.");

        // Cut-cell classification of the background grid.
        gsImplicitTrimmedDomain<3, real_t> tr_domain(impl_fun, *tbsPtr);

        // Algoim volumetric rule for {phi < 0} on cut cells (default dim = -1).
        gsAlgoimGenericRule<real_t> algoimRule(impl_fun, *tbsPtr);

        real_t volume = 0;
        gsMatrix<real_t> pts;
        gsVector<real_t> wts;

        // Interior cells (sign < 0) contribute their full (unit-Jacobian) volume.
        for (auto it = tr_domain.beginInterior();
             it != tr_domain.end<InteriorSign>(); ++it)
            volume += (it.upperCorner() - it.lowerCorner()).prod();

        // Cut cells (sign == 0): Algoim volumetric quadrature of {phi < 0}.
        // When Algoim returns no nodes (re-entrant / multi-crossing cell),
        // fall back to a uniform grid fraction estimate.
        for (auto it = tr_domain.beginBdr(boundary::none);
             it != tr_domain.endBdr(boundary::none); ++it)
        {
            algoimRule.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
            if (wts.size() > 0)
            {
                volume += wts.sum();
            }
            else
            {
                const gsVector<real_t> & lo_c = it.lowerCorner();
                const gsVector<real_t> & hi_c = it.upperCorner();
                const real_t cellVol = (hi_c - lo_c).prod();
                const index_t nTotal = nFallback * nFallback * nFallback;
                gsMatrix<real_t> spts(3, nTotal);
                index_t col = 0;
                for (index_t ii = 0; ii < nFallback; ++ii)
                for (index_t jj = 0; jj < nFallback; ++jj)
                for (index_t kk = 0; kk < nFallback; ++kk, ++col)
                {
                    spts(0,col) = lo_c[0] + (ii+0.5)/nFallback*(hi_c[0]-lo_c[0]);
                    spts(1,col) = lo_c[1] + (jj+0.5)/nFallback*(hi_c[1]-lo_c[1]);
                    spts(2,col) = lo_c[2] + (kk+0.5)/nFallback*(hi_c[2]-lo_c[2]);
                }
                gsMatrix<real_t> fbVals;
                impl_fun.eval_into(spts, fbVals);
                const index_t nInside = static_cast<index_t>((fbVals.array() < 0).count());
                volume += cellVol * static_cast<real_t>(nInside) / nTotal;
            }
        }

        const index_t cellsPerDim = tbsPtr->component(0).numElements();
        gsInfo << std::setw(8) << r << std::setw(10) << cellsPerDim
               << std::setw(18) << std::fixed << std::setprecision(6) << volume
               << std::setw(18) << std::fixed << std::setprecision(6)
               << volume * invScale3 << "\n";

        if (r == numRefine)
        {
            const real_t bboxVol = (physBox.col(1) - physBox.col(0)).prod();
            gsInfo << "\n=== Result (finest level) ===\n";
            gsInfo << "Volume of object (parametric [0,1]^3): "
                   << std::fixed << std::setprecision(6) << volume << "\n";
            gsInfo << "Volume of object (physical units):    "
                   << std::fixed << std::setprecision(6) << volume * invScale3 << "\n";
            gsInfo << "Percentage of [0,1]^3 box occupied:   "
                   << std::fixed << std::setprecision(2) << 100.0 * volume << " %\n";
            gsInfo << "Percentage of object bounding box:    "
                   << std::fixed << std::setprecision(2)
                   << 100.0 * (volume * invScale3) / bboxVol << " %\n\n";
        }

        // ---------------------------------------------------------------------
        // 7. ParaView output (final refinement level only)
        // ---------------------------------------------------------------------
        if (plot && r == numRefine)
        {
            // 7a. Rescaled input geometry as a triangle mesh.
            gsWriteParaview(mesh, out + "/" + outputStem + "_geometry");

            // 7b. Background grid.
            gsMesh<real_t> bgMesh(dbasis.basis(0));
            gsWriteParaview(bgMesh, out + "/" + outputStem + "_background_mesh");

            // 7c. Algoim volumetric quadrature points covering the WHOLE cow.
            //     Cut cells (sign == 0) use the Algoim {phi < 0} rule; interior
            //     cells (sign < 0) use a full Gauss rule so the solid body is
            //     filled with points and no region is left empty.
            //     Row 3 stores the quadrature weight.
            gsMatrix<real_t> quad(4, 0), phys;

            for (auto it = tr_domain.beginBdr(boundary::none);
                 it != tr_domain.endBdr(boundary::none); ++it)
            {
                algoimRule.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
                if (pts.cols() == 0)
                {
                    // Algoim failure (re-entrant cell): fall back to uniform grid.
                    const gsVector<real_t> & lo_c = it.lowerCorner();
                    const gsVector<real_t> & hi_c = it.upperCorner();
                    const real_t cellVol = (hi_c - lo_c).prod();
                    const index_t nTotal = nFallback * nFallback * nFallback;
                    gsMatrix<real_t> spts(3, nTotal);
                    index_t scol = 0;
                    for (index_t ii = 0; ii < nFallback; ++ii)
                    for (index_t jj = 0; jj < nFallback; ++jj)
                    for (index_t kk = 0; kk < nFallback; ++kk, ++scol)
                    {
                        spts(0,scol) = lo_c[0] + (ii+0.5)/nFallback*(hi_c[0]-lo_c[0]);
                        spts(1,scol) = lo_c[1] + (jj+0.5)/nFallback*(hi_c[1]-lo_c[1]);
                        spts(2,scol) = lo_c[2] + (kk+0.5)/nFallback*(hi_c[2]-lo_c[2]);
                    }
                    gsMatrix<real_t> fbVals;
                    impl_fun.eval_into(spts, fbVals);
                    std::vector<index_t> insideCols;
                    for (index_t ci = 0; ci < nTotal; ++ci)
                        if (fbVals(0,ci) < 0) insideCols.push_back(ci);
                    if (insideCols.empty()) continue;
                    pts.resize(3, static_cast<index_t>(insideCols.size()));
                    wts.resize(static_cast<index_t>(insideCols.size()));
                    // const real_t w = cellVol / static_cast<real_t>(insideCols.size());
                    const real_t w = cellVol / static_cast<real_t>(nTotal);
                    for (index_t ci = 0; ci < static_cast<index_t>(insideCols.size()); ++ci)
                    {
                        pts.col(ci) = spts.col(insideCols[ci]);
                        wts[ci] = w;
                    }
                }
                mp.patch(0).eval_into(pts, phys);
                const index_t c = quad.cols();
                quad.conservativeResize(4, c + pts.cols());
                quad.block(0, c, 3, pts.cols()) = phys;
                quad.row(3).segment(c, pts.cols()) = wts.transpose();
            }

            gsGaussRule<real_t> gauss(gsVector<index_t, 3>::Constant(
                dbasis.maxCwiseDegree() + 1));
            for (auto it = tr_domain.beginInterior();
                 it != tr_domain.end<InteriorSign>(); ++it)
            {
                gauss.mapTo(it.lowerCorner(), it.upperCorner(), pts, wts);
                if (pts.cols() == 0) continue;
                mp.patch(0).eval_into(pts, phys);
                const index_t c = quad.cols();
                quad.conservativeResize(4, c + pts.cols());
                quad.block(0, c, 3, pts.cols()) = phys;
                quad.row(3).segment(c, pts.cols()) = wts.transpose();
            }

            if (quad.cols() > 0)
            {
                const std::string base = out + "/" + outputStem + "_quadrature_points";
                gsInfo << "\nExporting " << quad.cols()
                       << " quadrature points to " << base << ".vtp\n";
                gsWriteParaviewPoints(quad, base);
            }

            gsInfo << "ParaView files written to: " << out << "\n";
        }
    }

    return EXIT_SUCCESS;
}
