/** @file commandLineArg_example.cpp

    @brief Batch-approximate G+Smo XML objects to multipatch and export.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

using namespace gismo;

namespace
{
    struct ApproximationSettings
    {
        unsigned samples;
        unsigned meshSamples;
        index_t refinementLevel;
        std::string suffix;
    };

    struct MirroredCheckResult
    {
        std::vector<bool> mirrored;
        std::vector<size_t> usedPerPatch;
        std::vector<size_t> posPerPatch;
        std::vector<size_t> negPerPatch;
    };

    enum class ProjectionPlane
    {
        None,
        XY,
        XZ,
        YZ
    };

    gsMatrix<> projectTo2D(const gsMatrix<>& xyz, ProjectionPlane plane)
    {
        if (plane == ProjectionPlane::None)
            return xyz;

        GISMO_ENSURE(xyz.rows() >= 3, "Projection requested but geometry has < 3 dimensions.");

        gsMatrix<> out(2, xyz.cols());
        if (plane == ProjectionPlane::XY)
        {
            out.row(0) = xyz.row(0);
            out.row(1) = xyz.row(1);
        }
        else if (plane == ProjectionPlane::XZ)
        {
            out.row(0) = xyz.row(0);
            out.row(1) = xyz.row(2);
        }
        else // YZ
        {
            out.row(0) = xyz.row(1);
            out.row(1) = xyz.row(2);
        }
        return out;
    }
    gsMatrix<> mapToSupport(const gsGeometry<>& geo, const gsMatrix<>& uvUnit)
    {
        const gsMatrix<> ab = geo.support();
        const gsVector<> a = ab.col(0);
        const gsVector<> b = ab.col(1);

        gsMatrix<> uvMapped = uvUnit;
        for (index_t k = 0; k < uvUnit.rows(); ++k)
        {
            uvMapped.row(k) = uvUnit.row(k) * (b(k) - a(k));
            uvMapped.row(k).array() += a(k);
        }

        return uvMapped;
    }

    gsMatrix<> samplePatchGrid(const gsBasis<>& basis, index_t samplesPerDirection)
    {
        const gsMatrix<> support = basis.support();
        gsVector<> low = support.col(0);
        gsVector<> upp = support.col(1);

        index_t totalSamples = 1;
        for (index_t dim = 0; dim < support.rows(); ++dim)
            totalSamples *= samplesPerDirection;

        return uniformPointGrid(low, upp, totalSamples);
    }

    MirroredCheckResult computeMirroredCheckPerPatch(const gsMultiPatch<>& mp,
                                                     index_t samplesPerDirection,
                                                     real_t detEps)
    {
        const index_t numPatches = static_cast<index_t>(mp.nPatches());
        MirroredCheckResult result;
        result.mirrored.assign(numPatches, false);
        result.usedPerPatch.assign(numPatches, 0);
        result.posPerPatch.assign(numPatches, 0);
        result.negPerPatch.assign(numPatches, 0);

        for (index_t patch = 0; patch < numPatches; ++patch)
        {
            const gsGeometry<>& geometry = mp.patch(patch);
            if (geometry.parDim() != geometry.geoDim())
                continue;

            const gsBasis<>& basis = geometry.basis();
            const gsMatrix<> uvPatch = samplePatchGrid(basis, samplesPerDirection);
            if (uvPatch.cols() == 0)
                continue;

            const gsMatrix<> support = basis.support();
            std::vector<real_t> margin(static_cast<size_t>(support.rows()), 0.0);
            for (index_t dim = 0; dim < support.rows(); ++dim)
                margin[static_cast<size_t>(dim)] = 0.05 * (support(dim, 1) - support(dim, 0));

            for (index_t pt = 0; pt < uvPatch.cols(); ++pt)
            {
                const gsVector<> param = uvPatch.col(pt);

                bool nearBoundary = false;
                for (index_t dim = 0; dim < support.rows(); ++dim)
                {
                    if (param(dim) <= support(dim, 0) + margin[static_cast<size_t>(dim)] ||
                        param(dim) >= support(dim, 1) - margin[static_cast<size_t>(dim)])
                    {
                        nearBoundary = true;
                        break;
                    }
                }
                if (nearBoundary)
                    continue;

                const real_t det = geometry.jacobian(param).determinant();
                if (!std::isfinite(det) || std::abs(det) <= detEps)
                    continue;

                ++result.usedPerPatch[patch];
                if (det > 0)
                    ++result.posPerPatch[patch];
                else
                    ++result.negPerPatch[patch];
            }
        }

        for (index_t patch = 0; patch < numPatches; ++patch)
        {
            if (result.usedPerPatch[patch] == 0)
                continue;

            if (result.posPerPatch[patch] == 0 && result.negPerPatch[patch] > 0)
            {
                result.mirrored[patch] = true;
                continue;
            }

            const real_t negRatio = static_cast<real_t>(result.negPerPatch[patch]) /
                                    static_cast<real_t>(result.usedPerPatch[patch]);
            result.mirrored[patch] =
                (result.negPerPatch[patch] > result.posPerPatch[patch]) && (negRatio > 0.9);
        }

        return result;
    }

    void logMirroredCheck(std::ostream& log,
                          const MirroredCheckResult& report,
                          const std::string& label,
                          bool logToInfo)
    {
        std::ostringstream buffer;
        buffer << "MIRROREDCHECK_BEGIN [" << label << "]\n";
        buffer << "patches=" << report.mirrored.size() << "\n";

        index_t mirroredCount = 0;
        for (index_t patch = 0; patch < static_cast<index_t>(report.mirrored.size()); ++patch)
        {
            if (report.mirrored[patch])
                ++mirroredCount;

            buffer << "patch " << patch
                   << ": used=" << report.usedPerPatch[patch]
                   << ", pos=" << report.posPerPatch[patch]
                   << ", neg=" << report.negPerPatch[patch];

            if (report.usedPerPatch[patch] > 0)
            {
                const real_t negRatio = static_cast<real_t>(report.negPerPatch[patch]) /
                                        static_cast<real_t>(report.usedPerPatch[patch]);
                buffer << ", negRatio=" << negRatio;
            }
            else
            {
                buffer << ", negRatio=n/a";
            }

            buffer << ", mirrored=" << (report.mirrored[patch] ? "true" : "false") << "\n";
        }

        buffer << "mirroredCount=" << mirroredCount << "\n";
        buffer << "MIRROREDCHECK_END\n";

        log << buffer.str();
        if (logToInfo)
            gsInfo << buffer.str();
    }

    // -----------------------------------------------------------------------
    // Wigglify: apply a sinusoidal displacement field to a multi-patch using
    // *physical* (x,y) coordinates as the wave argument.  Because two control
    // points that share an interface have identical (x,y) positions they
    // receive exactly the same (Δx,Δy) → C0 continuity is preserved without
    // any explicit topology bookkeeping.
    //
    // Displacement field:
    //   Δx = amplitude * sin(freqY * π * ty)   ty = (y-ymin)/(ymax-ymin)
    //   Δy = amplitude * sin(freqX * π * tx)   tx = (x-xmin)/(xmax-xmin)
    //
    // When keepOuterBoundary=true the CPs on sides that are NOT shared with
    // another patch (plus one interior row) are frozen so the outer shape of
    // the domain is unchanged.
    // -----------------------------------------------------------------------
    void wigglifyMultiPatch(gsMultiPatch<>& mp,
                            double amplitude,
                            double freqX,
                            double freqY,
                            bool keepOuterBoundary)
    {
        if (mp.nPatches() == 0) return;

        // --- Global bounding box of all control points ---
        double xmin =  1e100, xmax = -1e100;
        double ymin =  1e100, ymax = -1e100;
        for (index_t p = 0; p < static_cast<index_t>(mp.nPatches()); ++p)
        {
            const gsMatrix<>& C = mp.patch(p).coefs();
            for (index_t k = 0; k < C.rows(); ++k)
            {
                xmin = std::min(xmin, static_cast<double>(C(k,0)));
                xmax = std::max(xmax, static_cast<double>(C(k,0)));
                if (C.cols() >= 2)
                {
                    ymin = std::min(ymin, static_cast<double>(C(k,1)));
                    ymax = std::max(ymax, static_cast<double>(C(k,1)));
                }
            }
        }
        const double wx = (xmax > xmin) ? (xmax - xmin) : 1.0;
        const double wy = (ymax > ymin) ? (ymax - ymin) : 1.0;

        // --- Identify outer-boundary CPs that should be frozen ---
        std::vector<std::set<index_t>> frozen(mp.nPatches());
        if (keepOuterBoundary)
        {
            for (index_t p = 0; p < static_cast<index_t>(mp.nPatches()); ++p)
            {
                const gsBasis<>& basis = mp.patch(p).basis();
                for (int s = 1; s <= 4; ++s)
                {
                    boxSide side(s);
                    bool isIface = false;
                    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
                    {
                        if ((it->first().patch  == p && it->first().side()  == side) ||
                            (it->second().patch == p && it->second().side() == side))
                        { isIface = true; break; }
                    }
                    if (!isIface)
                    {
                        // Freeze boundary row + one interior row for a clean outer edge.
                        for (int offset = 0; offset <= 1; ++offset)
                        {
                            gsMatrix<index_t> ind = basis.boundaryOffset(side, offset);
                            for (index_t k = 0; k < ind.rows(); ++k)
                                frozen[p].insert(ind(k, 0));
                        }
                    }
                }
            }
        }

        // --- Apply displacement ---
        const double kPi = std::acos(-1.0);
        for (index_t p = 0; p < static_cast<index_t>(mp.nPatches()); ++p)
        {
            gsMatrix<>& C = mp.patch(p).coefs();
            for (index_t k = 0; k < C.rows(); ++k)
            {
                if (keepOuterBoundary && frozen[p].count(k)) continue;

                const double x  = static_cast<double>(C(k, 0));
                const double y  = static_cast<double>(C(k, 1));
                const double tx = (x - xmin) / wx;
                const double ty = (y - ymin) / wy;

                C(k, 0) += static_cast<real_t>(amplitude * std::sin(freqY * kPi * ty));
                C(k, 1) += static_cast<real_t>(amplitude * std::sin(freqX * kPi * tx));
            }
        }

        gsInfo << "[wigglify] amplitude=" << amplitude
               << "  freqX=" << freqX << "  freqY=" << freqY
               << "  keepOuterBoundary=" << (keepOuterBoundary ? "yes" : "no") << "\n";
    }

    // -----------------------------------------------------------------------
    // Wigglify interfaces: for every patch-to-patch interface, apply a
    // sinusoidal displacement *perpendicular* to the interface and
    // *sinusoidal along* the interface arc-length parameter t ∈ [0,1].
    //
    //   Δ(t) = amplitude · sin(freq · π · t)   applied in direction n̂
    //
    // where n̂ is the unit normal to the interface (first → last CP direction
    // rotated 90° CCW).  sin(0)=sin(integer·π)=0 so the endpoints (corners
    // shared with other interfaces / outer boundary) are never moved.
    //
    // Both patches receive the *same* displacement on the shared row → C0.
    // The adjacent interior row (offset=1) gets 0.5·Δ for a smooth fade-in.
    // -----------------------------------------------------------------------
    void wigglifyInterfacesMultiPatch(gsMultiPatch<>& mp,
                                      double amplitude,
                                      double freq)
    {
        if (mp.nPatches() == 0) return;
        const double kPi = std::acos(-1.0);

        for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
        {
            const index_t p1 = it->first().patch;
            const index_t p2 = it->second().patch;
            const boxSide  s1 = it->first().side();
            const boxSide  s2 = it->second().side();

            gsMatrix<>& C1 = mp.patch(p1).coefs();
            gsMatrix<>& C2 = mp.patch(p2).coefs();
            const gsBasis<>& b1 = mp.patch(p1).basis();
            const gsBasis<>& b2 = mp.patch(p2).basis();

            gsMatrix<index_t> ind1_0 = b1.boundaryOffset(s1, 0);
            gsMatrix<index_t> ind1_1 = b1.boundaryOffset(s1, 1);
            gsMatrix<index_t> ind2_0 = b2.boundaryOffset(s2, 0);
            gsMatrix<index_t> ind2_1 = b2.boundaryOffset(s2, 1);

            const index_t N = ind1_0.rows();
            if (N != ind2_0.rows() || N < 2) continue;

            // Orientation: ind2_0 may run opposite to ind1_0 — detect and reverse
            {
                const double d_same =
                    (C1.row(ind1_0(0,0))   - C2.row(ind2_0(0,  0))).squaredNorm() +
                    (C1.row(ind1_0(N-1,0)) - C2.row(ind2_0(N-1,0))).squaredNorm();
                const double d_flip =
                    (C1.row(ind1_0(0,0))   - C2.row(ind2_0(N-1,0))).squaredNorm() +
                    (C1.row(ind1_0(N-1,0)) - C2.row(ind2_0(0,  0))).squaredNorm();
                if (d_flip < d_same)
                    for (index_t k = 0; k < N / 2; ++k)
                    {
                        std::swap(ind2_0(k, 0), ind2_0(N-1-k, 0));
                        std::swap(ind2_1(k, 0), ind2_1(N-1-k, 0));
                    }
            }

            // Unit normal to the interface (tangent rotated 90° CCW)
            const double dx  = C1(ind1_0(N-1, 0), 0) - C1(ind1_0(0, 0), 0);
            const double dy  = C1(ind1_0(N-1, 0), 1) - C1(ind1_0(0, 0), 1);
            const double len = std::sqrt(dx*dx + dy*dy);
            if (len < 1e-12) continue;
            const double nx = -dy / len;
            const double ny =  dx / len;

            // Apply
            for (index_t k = 0; k < N; ++k)
            {
                const double t     = static_cast<double>(k) / static_cast<double>(N - 1);
                const double delta = amplitude * std::sin(freq * kPi * t);

                C1(ind1_0(k,0), 0) += delta * nx;   C1(ind1_0(k,0), 1) += delta * ny;
                C1(ind1_1(k,0), 0) += 0.5*delta*nx; C1(ind1_1(k,0), 1) += 0.5*delta*ny;

                C2(ind2_0(k,0), 0) += delta * nx;   C2(ind2_0(k,0), 1) += delta * ny;
                C2(ind2_1(k,0), 0) += 0.5*delta*nx; C2(ind2_1(k,0), 1) += 0.5*delta*ny;
            }

            gsInfo << "[wigglify-boundaries] iface p" << p1 << " s" << static_cast<int>(s1)
                   << " <-> p" << p2 << " s" << static_cast<int>(s2)
                   << "  N=" << N << "  normal=(" << nx << "," << ny << ")\n";
        }

        gsInfo << "[wigglify-boundaries] amplitude=" << amplitude
               << "  freq=" << freq << "\n";
    }

    bool runMirrorAwareJacobianCheck(const gsMultiPatch<>& inputMp,
                                     const gsMultiPatch<>& outputMp,
                                     index_t fitGrid,
                                     real_t detThreshold,
                                     const std::string& label,
                                     std::ostream& log,
                                     bool logToInfo)
    {
        const MirroredCheckResult inputMirrorReport =
            computeMirroredCheckPerPatch(inputMp, fitGrid, 1e-12);
        const MirroredCheckResult outputMirrorReport =
            computeMirroredCheckPerPatch(outputMp, fitGrid, 1e-12);

        logMirroredCheck(log, inputMirrorReport, label + " input", logToInfo);
        logMirroredCheck(log, outputMirrorReport, label + " output", logToInfo);

        log << "\nMirror baseline comparison [" << label << "]\n";
        index_t changedMirrorCount = 0;
        for (index_t patch = 0; patch < static_cast<index_t>(outputMirrorReport.mirrored.size()); ++patch)
        {
            const bool inputMirrored =
                (patch < static_cast<index_t>(inputMirrorReport.mirrored.size())) ? inputMirrorReport.mirrored[patch] : false;
            if (inputMirrored != outputMirrorReport.mirrored[patch])
            {
                ++changedMirrorCount;
                log << "Patch " << patch << ": mirror status changed from "
                    << (inputMirrored ? "true" : "false") << " to "
                    << (outputMirrorReport.mirrored[patch] ? "true" : "false") << "\n";
                if (logToInfo)
                    gsInfo << "Patch " << patch << ": mirror status changed from "
                           << (inputMirrored ? "true" : "false") << " to "
                           << (outputMirrorReport.mirrored[patch] ? "true" : "false") << "\n";
            }
        }
        if (changedMirrorCount == 0)
        {
            log << "Mirror status unchanged on all patches.\n";
            if (logToInfo)
                gsInfo << "Mirror status unchanged on all patches.\n";
        }

        log << "\nJacobian determinant check [" << label << "] (threshold=" << detThreshold << ")\n";

        bool hasInvalidPatch = false;
        for (index_t patch = 0; patch < static_cast<index_t>(outputMp.nPatches()); ++patch)
        {
            const gsGeometry<>& geometry = outputMp.patch(patch);
            if (geometry.parDim() != geometry.geoDim())
            {
                log << "Patch " << patch << ": skipped determinant check because parDim="
                    << geometry.parDim() << " and geoDim=" << geometry.geoDim() << "\n";
                if (logToInfo)
                    gsInfo << "Patch " << patch << ": skipped determinant check because parDim="
                           << geometry.parDim() << " and geoDim=" << geometry.geoDim() << "\n";
                continue;
            }

            const gsBasis<>& basis = geometry.basis();
            const gsMatrix<> uv = samplePatchGrid(basis, fitGrid);
            if (uv.cols() == 0)
                continue;

            real_t minRawDet = std::numeric_limits<real_t>::max();
            real_t minOrientedDet = std::numeric_limits<real_t>::max();
            index_t belowThresholdCount = 0;
            const bool patchMirrored =
                (patch < static_cast<index_t>(inputMirrorReport.mirrored.size())) ? inputMirrorReport.mirrored[patch] : false;

            for (index_t j = 0; j < uv.cols(); ++j)
            {
                const real_t det = geometry.jacobian(uv.col(j)).determinant();
                minRawDet = std::min(minRawDet, det);

                const real_t orientedDet = patchMirrored ? -det : det;
                minOrientedDet = std::min(minOrientedDet, orientedDet);
                if (orientedDet < detThreshold)
                    ++belowThresholdCount;
            }

            if (belowThresholdCount > 0)
                hasInvalidPatch = true;

            log << "Patch " << patch << ": mirroredBaseline=" << (patchMirrored ? "true" : "false")
                << ", min raw det=" << minRawDet
                << ", min oriented det=" << minOrientedDet
                << ", count below threshold=" << belowThresholdCount << "\n";
            if (logToInfo)
                gsInfo << "Patch " << patch << ": mirroredBaseline=" << (patchMirrored ? "true" : "false")
                       << ", min raw det=" << minRawDet
                       << ", min oriented det=" << minOrientedDet
                       << ", count below threshold=" << belowThresholdCount << "\n";
        }

        return hasInvalidPatch;
    }

    void logRefinementBoxAndBasis(const std::vector<index_t>& box,
                                  const gsBasis<>& basis,
                                  const char* label)
    {
        gsInfo << label << " refinement box: [";
        for (size_t i = 0; i < box.size(); ++i)
        {
            gsInfo << box[i];
            if (i + 1 < box.size())
                gsInfo << ", ";
        }
        gsInfo << "]\n";
        gsInfo << label << " THB basis after refinement:\n" << basis << "\n";
    }

    void approximateGeometryToMultiPatch(const gsGeometry<>& geo,
                                         gsMultiPatch<>& out,
                                         const ApproximationSettings& settings,
                                         ProjectionPlane projection,
                                         const std::vector<boxSide>& fixedSides,
                                         const std::vector<gsGeometry<>*>& fixedCurves)
    {
        const index_t d = geo.parDim();

        const gsMatrix<> ab = geo.support();
        const gsVector<> a = ab.col(0);
        const gsVector<> b = ab.col(1);

        gsVector<> aUnit(d);
        gsVector<> bUnit(d);
        aUnit.setZero();
        bUnit.setOnes();
        gsVector<unsigned> np = uniformSampleCount(aUnit, bUnit, settings.samples);
        gsMatrix<> uvUnit = gsPointGrid(aUnit, bUnit, np);
        gsMatrix<> uv = mapToSupport(geo, uvUnit);
        gsMatrix<> xyz = geo.eval(uv);
        if (projection != ProjectionPlane::None)
            xyz = projectTo2D(xyz, projection);

        const short_t degree = geo.basis().maxDegree();
        gsKnotVector<> kv(a(0), b(0), 0, degree + 1, 1);
        const index_t level = settings.refinementLevel;
        const index_t levelFactor = static_cast<index_t>(1) << level;

        if (d == 1)
        {
            gsBSplineBasis<> base(kv);
            gsTHBSplineBasis<1> thb(base);
            const index_t boxMax = static_cast<index_t>(base.numElements()) * levelFactor;
            std::vector<index_t> box;
            box.push_back(level);
            box.push_back(0);
            box.push_back(boxMax);
            thb.refineElements(box);
            logRefinementBoxAndBasis(box, thb, "1D");
            gsFitting<> fit(uv, xyz, thb);
            if (!fixedSides.empty())
                fit.setConstraints(fixedSides, fixedCurves);
            fit.compute();
            gsGeometry<>* res = fit.result();
            if (!res)
            {
                gsWarn << "Fitting failed for 1D geometry. Skipping.\n";
                return;
            }
            out.addPatch(*res);
            return;
        }
        if (d == 2)
        {
            gsKnotVector<> kvU(a(0), b(0), 0, degree + 1, 1);
            gsKnotVector<> kvV(a(1), b(1), 0, degree + 1, 1);
            gsTensorBSplineBasis<2> base(kvU, kvV);
            gsTHBSplineBasis<2> thb(base);
            gsVector<unsigned> ne(2);
            base.numElements_cwise(ne);
            const index_t boxMaxU = static_cast<index_t>(ne(0)) * levelFactor;
            const index_t boxMaxV = static_cast<index_t>(ne(1)) * levelFactor;
            std::vector<index_t> box;
            box.push_back(level);
            box.push_back(0);
            box.push_back(0);
            box.push_back(boxMaxU);
            box.push_back(boxMaxV);
            thb.refineElements(box);
            logRefinementBoxAndBasis(box, thb, "2D");
            gsFitting<> fit(uv, xyz, thb);
            // if (!fixedSides.empty())
            //     fit.setConstraints(fixedSides, fixedCurves);
            fit.compute();
            gsGeometry<>* res = fit.result();
            if (!res)
            {
                gsWarn << "Fitting failed for 2D geometry. Skipping.\n";
                return;
            }
            out.addPatch(*res);
            return;
        }
        if (d == 3)
        {
            gsKnotVector<> kvU(a(0), b(0), 0, degree + 1, 1);
            gsKnotVector<> kvV(a(1), b(1), 0, degree + 1, 1);
            gsKnotVector<> kvW(a(2), b(2), 0, degree + 1, 1);
            gsTensorBSplineBasis<3> base(kvU, kvV, kvW);
            gsTHBSplineBasis<3> thb(base);
            gsVector<unsigned> ne(3);
            base.numElements_cwise(ne);
            const index_t boxMaxU = static_cast<index_t>(ne(0)) * levelFactor;
            const index_t boxMaxV = static_cast<index_t>(ne(1)) * levelFactor;
            const index_t boxMaxW = static_cast<index_t>(ne(2)) * levelFactor;
            std::vector<index_t> box;
            box.push_back(level);
            box.push_back(0);
            box.push_back(0);
            box.push_back(0);
            box.push_back(boxMaxU);
            box.push_back(boxMaxV);
            box.push_back(boxMaxW);
            thb.refineElements(box);
            logRefinementBoxAndBasis(box, thb, "3D");
            gsFitting<> fit(uv, xyz, thb);
            // if (!fixedSides.empty())
            //     fit.setConstraints(fixedSides, fixedCurves);
            fit.compute();
            gsGeometry<>* res = fit.result();
            if (!res)
            {
                gsWarn << "Fitting failed for 3D geometry. Skipping.\n";
                return;
            }
            out.addPatch(*res);
            return;
        }

        gsWarn << "Unsupported parametric dimension " << d << " for geometry. Skipping.\n";
    }

    void approximateMultiPatchToMultiPatch(const gsMultiPatch<>& mp,
                                           gsMultiPatch<>& out,
                                           const ApproximationSettings& settings,
                                           ProjectionPlane projection)
    {
        for (size_t i = 0; i < mp.nPatches(); ++i)
        {
            std::vector<boxSide> fixedSides;
            std::vector<gsGeometry<>::uPtr> fixedCurvesOwned;
            std::vector<gsGeometry<>*> fixedCurves;

            for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
            {
                if (it->first().patch == static_cast<index_t>(i))
                {
                    boxSide side = it->first().side();
                    typename gsGeometry<>::uPtr curve = mp.patch(i).boundary(side);
                    if (!curve)
                    {
                        gsWarn << "Failed to extract boundary curve for patch " << i << ".\n";
                        continue;
                    }
                    fixedSides.push_back(side);
                    fixedCurves.push_back(curve.get());
                    fixedCurvesOwned.push_back(std::move(curve));
                }
                else if (it->second().patch == static_cast<index_t>(i))
                {
                    boxSide side = it->second().side();
                    typename gsGeometry<>::uPtr curve = mp.patch(i).boundary(side);
                    if (!curve)
                    {
                        gsWarn << "Failed to extract boundary curve for patch " << i << ".\n";
                        continue;
                    }
                    fixedSides.push_back(side);
                    fixedCurves.push_back(curve.get());
                    fixedCurvesOwned.push_back(std::move(curve));
                }
            }

            approximateGeometryToMultiPatch(mp.patch(i), out, settings, projection, fixedSides, fixedCurves);
        }

        /* Post-process: enforce C0 by averaging interface control points (disabled)
        for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
        {
            const index_t p1 = it->first().patch;
            const index_t p2 = it->second().patch;
            if (p1 >= static_cast<index_t>(out.nPatches()) || p2 >= static_cast<index_t>(out.nPatches()))
                continue;

            const gsGeometry<>& g1 = out.patch(p1);
            const gsGeometry<>& g2 = out.patch(p2);
            if (g1.parDim() != 2 || g2.parDim() != 2)
                continue;

            const boxSide s1 = it->first().side();
            const boxSide s2 = it->second().side();
            gsMatrix<index_t> ind1 = g1.basis().boundaryOffset(s1, 0);
            gsMatrix<index_t> ind2 = g2.basis().boundaryOffset(s2, 0);

            std::vector<index_t> ids1(ind1.rows());
            std::vector<index_t> ids2(ind2.rows());
            for (index_t k = 0; k < ind1.rows(); ++k)
                ids1[k] = ind1(k, 0);
            for (index_t k = 0; k < ind2.rows(); ++k)
                ids2[k] = ind2(k, 0);

            bool reverseIds2 = false;
            {
                typename gsGeometry<>::uPtr c1 = g1.boundary(s1);
                typename gsGeometry<>::uPtr c2 = g2.boundary(s2);
                if (c1 && c2)
                {
                    gsMatrix<> supp1 = c1->support();
                    gsMatrix<> supp2 = c2->support();
                    gsMatrix<> u1(1, 1);
                    gsMatrix<> u2(1, 1);
                    u1(0, 0) = supp1(0, 0);
                    u2(0, 0) = supp1(0, 1);
                    gsMatrix<> p1s = c1->eval(u1);
                    gsMatrix<> p1e = c1->eval(u2);

                    u1(0, 0) = supp2(0, 0);
                    u2(0, 0) = supp2(0, 1);
                    gsMatrix<> p2s = c2->eval(u1);
                    gsMatrix<> p2e = c2->eval(u2);

                    const real_t distSame = (p1s - p2s).squaredNorm() + (p1e - p2e).squaredNorm();
                    const real_t distFlip = (p1s - p2e).squaredNorm() + (p1e - p2s).squaredNorm();
                    reverseIds2 = (distFlip < distSame);
                }
            }
            if (reverseIds2)
                std::reverse(ids2.begin(), ids2.end());

            if (ids1.size() != ids2.size())
            {
                gsWarn << "Interface size mismatch between patches " << p1 << " and " << p2 << ".\n";
                continue;
            }

            gsGeometry<>& gw1 = out.patch(p1);
            gsGeometry<>& gw2 = out.patch(p2);
            const index_t geoDim = gw1.geoDim();
            for (size_t k = 0; k < ids1.size(); ++k)
            {
                for (index_t c = 0; c < geoDim; ++c)
                {
                    const real_t avg = static_cast<real_t>(0.5) * (gw1.coef(ids1[k], c) + gw2.coef(ids2[k], c));
                    gw1.coef(ids1[k], c) = avg;
                    gw2.coef(ids2[k], c) = avg;
                }
            }
        }
        */
    }

    size_t countBasisFunctionsInFile(const gsFileData<>& fileData)
    {
        size_t total = 0;

        if (fileData.has<gsMultiPatch<>>())
        {
            auto mps = fileData.getAll<gsMultiPatch<>>();
            for (auto& mpPtr : mps)
            {
                if (!mpPtr)
                    continue;
                for (size_t i = 0; i < mpPtr->nPatches(); ++i)
                    total += static_cast<size_t>(mpPtr->patch(i).basis().size());
            }
        }

        if (fileData.has<gsGeometry<>>())
        {
            auto geos = fileData.getAll<gsGeometry<>>();
            for (auto& geoPtr : geos)
            {
                if (!geoPtr)
                    continue;
                total += static_cast<size_t>(geoPtr->basis().size());
            }
        }

        if (fileData.has<gsTrimSurface<>>())
        {
            auto trims = fileData.getAll<gsTrimSurface<>>();
            for (auto& trimPtr : trims)
            {
                if (!trimPtr)
                    continue;
                gsSurface<>::Ptr surf = trimPtr->getTP();
                if (!surf)
                    continue;
                total += static_cast<size_t>(surf->basis().size());
            }
        }

        return total;
    }

    bool approximateFile(const std::string& path,
                         const ApproximationSettings& settings,
                         ProjectionPlane projection,
                         const std::string& suffix)
    {
        gsFileData<> fileData;
        if (!fileData.read(path))
        {
            gsWarn << "Failed to read file: " << path << "\n";
            return false;
        }

        gsMultiPatch<> out;
    gsMultiPatch<> inputForCheck;
        bool any = false;

        const size_t basisCount = countBasisFunctionsInFile(fileData);
        gsInfo << "Total basis functions (input): " << basisCount << "\n";

        bool hadMultiPatch = false;
        const gsMultiPatch<>* topoSource = nullptr; // input MP with topology, for fallback
        if (fileData.has<gsMultiPatch<>>())
        {
            auto mps = fileData.getAll<gsMultiPatch<>>();
            for (auto& mpPtr : mps)
            {
                if (!mpPtr)
                    continue;
                for (size_t i = 0; i < mpPtr->nPatches(); ++i)
                    inputForCheck.addPatch(mpPtr->patch(i));
                approximateMultiPatchToMultiPatch(*mpPtr, out, settings, projection);
                if (mpPtr->nInterfaces() > 0)
                    topoSource = mpPtr.get();
                any = true;
                hadMultiPatch = true;
            }
        }

        if (!hadMultiPatch && fileData.has<gsGeometry<>>())
        {
            auto geos = fileData.getAll<gsGeometry<>>();
            for (auto& geoPtr : geos)
            {
                if (!geoPtr)
                    continue;
                const std::vector<boxSide> noSides;
                const std::vector<gsGeometry<>*> noCurves;
                inputForCheck.addPatch(*geoPtr);
                approximateGeometryToMultiPatch(*geoPtr, out, settings, projection, noSides, noCurves);
                any = true;
            }
        }

        if (fileData.has<gsTrimSurface<>>())
        {
            auto trims = fileData.getAll<gsTrimSurface<>>();
            for (auto& trimPtr : trims)
            {
                if (!trimPtr)
                    continue;
                gsSurface<>::Ptr surf = trimPtr->getTP();
                if (!surf)
                    continue;
                const std::vector<boxSide> noSides;
                const std::vector<gsGeometry<>*> noCurves;
                inputForCheck.addPatch(*surf);
                approximateGeometryToMultiPatch(*surf, out, settings, projection, noSides, noCurves);
                any = true;
            }
        }

        if (!any)
        {
            gsWarn << "No geometry or multipatch found in: " << path << "\n";
            return false;
        }

        // Keep topology initialized before mirror-aware checks.
        if (inputForCheck.nPatches() > 0)
            inputForCheck.computeTopology();
        if (out.nPatches() > 0)
        {
            out.computeTopology();
            // computeTopology() may fail to detect interfaces when LS-fitted
            // boundary curves differ by more than the detection tolerance.
            // Fall back to copying the topology from the loaded input, which
            // preserves patch ordering and side connectivity exactly.
            if (out.nInterfaces() == 0 && topoSource && topoSource->nInterfaces() > 0)
            {
                gsInfo << "[topology] computeTopology found no interfaces; "
                          "copying topology from input (" << topoSource->nInterfaces()
                       << " interfaces).\n";
                for (auto it = topoSource->iBegin(); it != topoSource->iEnd(); ++it)
                    out.addInterface(*it);
                for (auto it = topoSource->bBegin(); it != topoSource->bEnd(); ++it)
                    out.addBoundary(*it);
            }
        }

        const std::string outputDir = R"(C:\Users\heydatey\source\repos\gismo\filedata\generatedMPs\)";
        gsFileManager::mkdir(outputDir);
        const std::string outputStem = outputDir +
                      gsFileManager::getBasename(path) +
                      "_approximation" +
                      suffix;
        const std::string outputPath = outputStem + ".xml";
        const std::string paraviewBase = outputStem;

        if (inputForCheck.nPatches() == out.nPatches())
        {
            std::ostringstream detLog;
            const bool hasInvalidDet = runMirrorAwareJacobianCheck(
                inputForCheck,
                out,
                40,
                -1e-8,
                gsFileManager::getBasename(path),
                detLog,
                true);

            if (hasInvalidDet)
                gsWarn << "Mirror-aware determinant check reported invalid oriented Jacobians for "
                       << path << "\n";
        }
        else
        {
            gsWarn << "Skipping mirror-aware determinant check because source/output patch counts differ ("
                   << inputForCheck.nPatches() << " vs " << out.nPatches() << ").\n";
        }

        gsFileData<> outData;
        outData << out;
        outData.dump(outputPath);
        gsInfo << "Saved XML: " << outputPath << "\n";

        gsWriteParaview(out, paraviewBase, settings.meshSamples, true, false);
        gsInfo << "Saved ParaView collection: " << paraviewBase + ".pvd" << "\n";
        gsFileManager::open(paraviewBase + ".pvd");
        return true;
    }
}

int main(int argc, char* argv[])
{
    // CLI flags:
    //   --input  <path>          — file to process (overrides the hardcoded list)
    //   --output-name <stem>     — output filename stem (written to filedata/generatedMPs/)
    //   --wigglify               — skip fitting; apply sinusoidal wiggle directly to CPs
    //   --wigglify-amplitude <A> — peak displacement (default 0.03)
    //   --wigglify-freq <n>      — wave frequency for both axes (default 3)
    std::string cliInput;
    std::string cliOutputName;
    int    cliRefinementLevel   = -1;   // -1 = use hardcoded default (3)
    bool   cliWigglify           = false;
    bool   cliWigglifyBoundaries = false;
    double cliWigglifyAmplitude  = 0.03;
    double cliWigglifyFreqX     = 3.0;
    double cliWigglifyFreqY     = 3.0;
    for (int ai = 1; ai < argc; ++ai)
    {
        if (std::string(argv[ai]) == "--input" && ai + 1 < argc)
        {
            cliInput = std::string(argv[ai + 1]);
            ++ai;
        }
        else if (std::string(argv[ai]) == "--output-name" && ai + 1 < argc)
        {
            cliOutputName = std::string(argv[ai + 1]);
            ++ai;
        }
        else if (std::string(argv[ai]) == "--refinement-level" && ai + 1 < argc)
        {
            cliRefinementLevel = std::stoi(argv[ai + 1]);
            ++ai;
        }
        else if (std::string(argv[ai]) == "--wigglify")
        {
            cliWigglify = true;
        }
        else if (std::string(argv[ai]) == "--wigglify-boundaries")
        {
            cliWigglifyBoundaries = true;
            cliWigglify = true;  // reuse the same main branch
        }
        else if (std::string(argv[ai]) == "--wigglify-amplitude" && ai + 1 < argc)
        {
            cliWigglifyAmplitude = std::stod(argv[ai + 1]);
            ++ai;
        }
        else if (std::string(argv[ai]) == "--wigglify-freq" && ai + 1 < argc)
        {
            cliWigglifyFreqX = cliWigglifyFreqY = std::stod(argv[ai + 1]);
            ++ai;
        }
    }

    const index_t refinementLevel = (cliRefinementLevel >= 0) ? cliRefinementLevel : 3;
    const std::string levelSuffix = "_fine_L" + std::to_string(refinementLevel);
    const ApproximationSettings settings = {2000u, 8u, refinementLevel, levelSuffix};

    // ----------------------------------------------------------------
    // --wigglify mode: load geometry as-is, perturb CPs, save — no LS
    // fitting is performed.  Requires --input and --output-name.
    // ----------------------------------------------------------------
    if (cliWigglify)
    {
        if (cliInput.empty() || cliOutputName.empty())
        {
            gsWarn << "--wigglify requires both --input and --output-name.\n";
            return 1;
        }

        gsInfo << "[wigglify] Input:  " << cliInput       << "\n";
        gsInfo << "[wigglify] Output: " << cliOutputName  << "\n";
        gsInfo << "[wigglify] amplitude=" << cliWigglifyAmplitude
               << "  freqX=" << cliWigglifyFreqX
               << "  freqY=" << cliWigglifyFreqY << "\n";

        gsFileData<> fileData;
        if (!fileData.read(cliInput))
        {
            gsWarn << "Failed to read: " << cliInput << "\n";
            return 1;
        }

        gsMultiPatch<> mp;
        if (fileData.has<gsMultiPatch<>>())
        {
            fileData.getFirst<gsMultiPatch<>>(mp);
        }
        else if (fileData.has<gsGeometry<>>())
        {
            auto geos = fileData.getAll<gsGeometry<>>();
            for (auto& g : geos)
                if (g) mp.addPatch(*g);
        }

        if (mp.nPatches() == 0)
        {
            gsWarn << "No patches found in: " << cliInput << "\n";
            return 1;
        }

        mp.computeTopology();
        gsInfo << "[wigglify] Loaded " << mp.nPatches() << " patches.\n";

        // Jacobian check on the input (baseline)
        {
            std::ostringstream detLog;
            runMirrorAwareJacobianCheck(mp, mp, 200, -1e-8,
                                        cliOutputName + "_input", detLog, true);
        }

        if (cliWigglifyBoundaries)
            wigglifyInterfacesMultiPatch(mp,
                                         cliWigglifyAmplitude,
                                         cliWigglifyFreqX);
        else
            wigglifyMultiPatch(mp,
                               cliWigglifyAmplitude,
                               cliWigglifyFreqX,
                               cliWigglifyFreqY,
                               /*keepOuterBoundary=*/true);

        // Jacobian check on the result
        {
            std::ostringstream detLog;
            const bool bad = runMirrorAwareJacobianCheck(mp, mp, 200, -1e-8,
                                                          cliOutputName, detLog, true);
            if (bad)
            {
                gsWarn << "[wigglify] HARD FAIL: output has invalid Jacobians "
                          "(min det < 0 detected at 200-pt grid) — "
                          "geometry NOT saved. Reduce --wigglify-amplitude.\n";
                return 1;
            }
        }

        const std::string outputDir =
            R"(C:\Users\heydatey\source\repos\gismo\filedata\generatedMPs\)";
        gsFileManager::mkdir(outputDir);
        const std::string outputStem = outputDir + cliOutputName;

        gsFileData<> outData;
        outData << mp;
        outData.dump(outputStem + ".xml");
        gsInfo << "Saved XML: " << outputStem + ".xml\n";

        gsWriteParaview(mp, outputStem, settings.meshSamples, true, false);
        gsInfo << "Saved ParaView: " << outputStem + ".pvd\n";

        return 0;
    }

    // Build the list of files to process.
    // If --input was supplied on the command line it takes priority over the
    // hardcoded list below.
    const std::vector<std::string> hardcodedFiles = {
        // R"(C:\Users\heydatey\source\repos\gismo\filedata\breps\other\TUDflame.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\filedata\domain2d\yeti_mp2.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\filedata\domain2d\yeti_mp2_thb.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\filedata\planar\hexagon_5p.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\unsupported\camembert2d_u3p.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\unsupported\ellipse_hole.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\unsupported\g_cal.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\unsupported\joystick.xml)",
        R"(C:\Users\heydatey\source\repos\gismo\unsupported\mask.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\unsupported\multipatch_tunnel_thb.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\unsupported\thb_shelfSupport.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\unsupported\triangle2d_u3p.xml)",
        // R"(C:\Users\heydatey\source\repos\gismo\unsupported\tv.xml)",
    };

    const std::vector<std::string> files =
        cliInput.empty() ? hardcodedFiles : std::vector<std::string>{cliInput};

    gsInfo << "Batch approximation started. THB basis with input knots"
           << ", samples=" << settings.samples
           << ", meshSamples=" << settings.meshSamples
           << ", refinementLevel=" << settings.refinementLevel
           << ", suffix=" << settings.suffix << "\n";
    if (!cliInput.empty())
        gsInfo << "[flag] --input=" << cliInput << "\n";
    if (!cliOutputName.empty())
        gsInfo << "[flag] --output-name=" << cliOutputName << "\n";

    for (const auto& path : files)
    {
        gsInfo << "\nProcessing: " << path << "\n";
        // When --output-name is given for a single file, pass it directly as
        // the suffix so the final stem becomes exactly <outputName>.
        // The approximateFile function builds:  outputDir + basename(path) + "_approximation" + suffix
        // To get outputDir/<cliOutputName> we construct a suffix that, combined
        // with the basename logic, yields the desired name.
        // Simplest: if --output-name is set, we rewrite the output path manually
        // by calling approximateFile with a sentinel suffix and then renaming —
        // but that couples things.  Instead we pass cliOutputName as the suffix
        // after stripping the basename+_approximation prefix from it, which is
        // fragile.  Cleanest: pass the desired full stem override directly.
        //
        // We therefore extend approximateFile's suffix parameter to support an
        // override.  When cliOutputName is non-empty, pass it as-is and use a
        // special overloadthat skips basename decoration.
        const std::string effectiveSuffix =
            cliOutputName.empty() ? settings.suffix : std::string();
        const std::string effectiveOutputName = cliOutputName;

        // Use a small lambda to keep it self-contained.
        bool ok = false;
        if (effectiveOutputName.empty())
        {
            ok = approximateFile(path, settings, ProjectionPlane::None, settings.suffix);
        }
        else
        {
            // Duplicate approximateFile inline with a fixed output stem.
            gsFileData<> fileData;
            if (!fileData.read(path))
            {
                gsWarn << "Failed to read file: " << path << "\n";
            }
            else
            {
                gsMultiPatch<> out;
                gsMultiPatch<> inputForCheck;
                bool any = false;
                const gsMultiPatch<>* topoSource2 = nullptr; // input MP with topology, for fallback
                if (fileData.has<gsMultiPatch<>>())
                {
                    auto mps = fileData.getAll<gsMultiPatch<>>();
                    for (auto& mpPtr : mps)
                    {
                        if (!mpPtr) continue;
                        for (size_t i = 0; i < mpPtr->nPatches(); ++i)
                            inputForCheck.addPatch(mpPtr->patch(i));
                        approximateMultiPatchToMultiPatch(*mpPtr, out, settings, ProjectionPlane::None);
                        if (mpPtr->nInterfaces() > 0)
                            topoSource2 = mpPtr.get();
                        any = true;
                    }
                }
                if (!any && fileData.has<gsGeometry<>>())
                {
                    auto geos = fileData.getAll<gsGeometry<>>();
                    for (auto& geoPtr : geos)
                    {
                        if (!geoPtr) continue;
                        const std::vector<boxSide> noSides;
                        const std::vector<gsGeometry<>*> noCurves;
                        inputForCheck.addPatch(*geoPtr);
                        approximateGeometryToMultiPatch(*geoPtr, out, settings, ProjectionPlane::None, noSides, noCurves);
                        any = true;
                    }
                }
                if (any)
                {
                    if (inputForCheck.nPatches() > 0)
                        inputForCheck.computeTopology();
                    if (out.nPatches() > 0)
                    {
                        out.computeTopology();
                        if (out.nInterfaces() == 0 && topoSource2 && topoSource2->nInterfaces() > 0)
                        {
                            gsInfo << "[topology] computeTopology found no interfaces; "
                                      "copying topology from input (" << topoSource2->nInterfaces()
                                   << " interfaces).\n";
                            for (auto it = topoSource2->iBegin(); it != topoSource2->iEnd(); ++it)
                                out.addInterface(*it);
                            for (auto it = topoSource2->bBegin(); it != topoSource2->bEnd(); ++it)
                                out.addBoundary(*it);
                        }
                    }

                    if (inputForCheck.nPatches() == out.nPatches())
                    {
                        std::ostringstream detLog;
                        const bool hasInvalidDet = runMirrorAwareJacobianCheck(
                            inputForCheck, out, 40, -1e-8,
                            effectiveOutputName, detLog, true);
                        if (hasInvalidDet)
                            gsWarn << "Mirror-aware determinant check reported invalid oriented Jacobians.\n";
                    }

                    const std::string outputDir =
                        R"(C:\Users\heydatey\source\repos\gismo\filedata\generatedMPs\)";
                    gsFileManager::mkdir(outputDir);
                    const std::string outputStem = outputDir + effectiveOutputName;
                    const std::string outputPath = outputStem + ".xml";

                    gsFileData<> outData;
                    outData << out;
                    outData.dump(outputPath);
                    gsInfo << "Saved XML: " << outputPath << "\n";

                    gsWriteParaview(out, outputStem, settings.meshSamples, true, false);
                    gsInfo << "Saved ParaView collection: " << outputStem + ".pvd" << "\n";
                    ok = true;
                }
            }
        }
        gsInfo << (ok ? "Done." : "Skipped.") << "\n";
    }

    return 0;
}
