/*  @file fitting_mspline.cpp

    @brief Multipatch distortion tool preserving exact C0 continuity.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include <gismo.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <vector>

using namespace gismo;

namespace
{
    struct MirroredCheckResult
    {
        std::vector<bool> mirrored;
        std::vector<size_t> usedPerPatch;
        std::vector<size_t> posPerPatch;
        std::vector<size_t> negPerPatch;
    };

    struct Options
    {
        index_t fitGrid = 20;         // Sampling density for regularity diagnostics.
        real_t detThreshold = -1e-8;  // Jacobian warning threshold for reporting.
        real_t bisectTol = 1e-2;      // Bisection interval tolerance for strength/alpha tuning.
        std::string outputFile = ""; // If empty, use <input_basename>_NLO.xml.
        std::string logFile = "distortion_log.txt";

        // Per-patch overrides for targeted near-singular construction.
        // --patch-det <p>:<val>   : override minOrientedDet threshold for patch p only.
        // --patch-focal <p>:<u>:<v> : override bell focal point (u,v) for patch p only.
        std::map<index_t, real_t> perPatchMinDet;
        std::map<index_t, real_t> perPatchFocalU;
        std::map<index_t, real_t> perPatchFocalV;
    };

    struct nonlinearOptions
    {
        index_t patchIndex = -1;               // -1 selects a default patch automatically.
        real_t strength = 0.25;                // Single distortion strength parameter.
        real_t minOrientedDet = 1e-3;          // Keep perturbation initially regular.
        real_t swirlAlpha = 0.0;               // Swirl rotation angle (radians). 0 = disabled.
        real_t bellSigma = 1.0;                // Gaussian bell width as fraction of patch radius.
        real_t focalU = -1.0;                  // Normalized [0,1] u-parameter for bell center. <0 = centroid.
        real_t focalV = -1.0;                  // Normalized [0,1] v-parameter for bell center. <0 = centroid.
        real_t sfoldAmplitude = 0.0;           // Peak S-fold displacement as fraction of patchRadius. 0 = disabled.
        real_t sfoldHalfWidth = 0.15;          // Half-width of the S-fold zone as fraction of patchRadius.
    };

    constexpr real_t kPi = 3.14159265358979323846;
    constexpr real_t kC0Tolerance = 1e-10;

    bool nudgeInterfaceControlPoint(gsGeometry<>& patch,
                                    const gsMultiPatch<>& mp,
                                    index_t patchIndex,
                                    real_t moveAbs);

    gsMatrix<> linspace(real_t a, real_t b, index_t n)
    {
        gsMatrix<> t(1, n);
        for (index_t i = 0; i < n; ++i)
        {
            real_t s = (n == 1) ? 0.0 : static_cast<real_t>(i) / (n - 1);
            t(0, i) = a + s * (b - a);
        }
        return t;
    }

    gsVector<> curveNormalFromTangent(const gsVector<>& tangent)
    {
        const index_t d = tangent.rows();
        gsVector<> n(d);
        if (d == 2)
        {
            n(0) = -tangent(1);
            n(1) = tangent(0);
        }
        else
        {
            gsVector<> ref(3);
            ref << 0, 0, 1;
            n = tangent.cross(ref);
            if (n.norm() < 1e-12)
            {
                ref << 0, 1, 0;
                n = tangent.cross(ref);
            }
            if (n.norm() < 1e-12)
            {
                ref << 1, 0, 0;
                n = tangent.cross(ref);
            }
        }
        if (n.norm() > 0)
            n.normalize();
        return n;
    }

    gsMatrix<> distortCurvePoints(const gsGeometry<>& curve,
                                  const gsMatrix<>& params,
                                  real_t epsAbs,
                                  index_t waveNumber)
    {
        gsMatrix<> pts = curve.eval(params);
        const index_t N = params.cols();
        if (N < 2 || epsAbs == 0)
            return pts;

        for (index_t i = 0; i < N; ++i)
        {
            index_t i0 = (i == 0) ? 0 : i - 1;
            index_t i1 = (i + 1 >= N) ? N - 1 : i + 1;
            gsVector<> tangent = pts.col(i1) - pts.col(i0);
            if (tangent.norm() == 0)
                continue;

            gsVector<> normal = curveNormalFromTangent(tangent);
            real_t s = static_cast<real_t>(i) / (N - 1);
            real_t amp = epsAbs * std::sin(static_cast<real_t>(waveNumber) * 2.0 * kPi * s);
            pts.col(i) += amp * normal;
        }
        return pts;
    }

    gsGeometry<>::uPtr fitCurveToPoints(const gsBasis<>& basis,
                                        const gsMatrix<>& params,
                                        const gsMatrix<>& points)
    {
        gsFitting<> fit(params, points, const_cast<gsBasis<>&>(basis));
        fit.compute();
        return gsGeometry<>::uPtr(fit.result()->clone());
    }

    gsGeometry<>::uPtr makeDistortedBoundaryCurve(const gsGeometry<>& patch,
                                                  boxSide side,
                                                  real_t epsAbs,
                                                  index_t waveNumber,
                                                  index_t curveSamples)
    {
        gsGeometry<>::uPtr curve = patch.boundary(side);
        const gsMatrix<> support = curve->basis().support();
        gsMatrix<> params = linspace(support(0, 0), support(0, 1), curveSamples);
        gsMatrix<> pts = distortCurvePoints(*curve, params, epsAbs, waveNumber);
        return fitCurveToPoints(curve->basis(), params, pts);
    }

    void applyBoundaryCurveToPatch(gsGeometry<>& patch,
                                   boxSide side,
                                   const gsGeometry<>& curve)
    {
        gsMatrix<index_t> boundaryIds = patch.basis().boundaryOffset(side, 0);
        GISMO_ENSURE(boundaryIds.rows() == curve.coefsSize(),
                     "Boundary curve size does not match patch boundary coefficients.");

        for (index_t k = 0; k < boundaryIds.rows(); ++k)
        {
            const index_t patchCoef = boundaryIds(k, 0);
            for (index_t c = 0; c < patch.geoDim(); ++c)
                patch.coef(patchCoef, c) = curve.coef(k, c);
        }
    }

    void applyWaveToBoundaryControlPoints(gsGeometry<>& patch,
                                          boxSide side,
                                          real_t epsAbs,
                                          index_t waveNumber)
    {
        if (epsAbs == 0)
            return;

        gsMatrix<index_t> boundaryIds = patch.basis().boundaryOffset(side, 0);
        const index_t n = boundaryIds.rows();
        if (n <= 2)
            return;

        for (index_t k = 1; k < n - 1; ++k)
        {
            const index_t idxPrev = boundaryIds(k - 1, 0);
            const index_t idxCur = boundaryIds(k, 0);
            const index_t idxNext = boundaryIds(k + 1, 0);

            gsVector<> tangent = patch.coef(idxNext) - patch.coef(idxPrev);
            if (tangent.norm() == 0)
                continue;

            gsVector<> normal = curveNormalFromTangent(tangent);
            real_t s = static_cast<real_t>(k) / (n - 1);
            real_t amp = epsAbs * std::sin(static_cast<real_t>(waveNumber) * 2.0 * kPi * s);
            for (index_t c = 0; c < patch.geoDim(); ++c)
                patch.coef(idxCur, c) += amp * normal(c);
        }
    }

    void applyWaveToInterfaceControlPoints(gsGeometry<>& patch1,
                                           boxSide side1,
                                           gsGeometry<>& patch2,
                                           boxSide side2,
                                           bool orientSame,
                                           real_t epsAbs,
                                           index_t waveNumber)
    {
        if (epsAbs == 0)
            return;

        gsMatrix<index_t> boundaryIds1 = patch1.basis().boundaryOffset(side1, 0);
        gsMatrix<index_t> boundaryIds2 = patch2.basis().boundaryOffset(side2, 0);
        GISMO_ENSURE(boundaryIds1.rows() == boundaryIds2.rows(),
                     "Interface sides do not have matching coefficient counts.");

        const index_t n = boundaryIds1.rows();
        if (n <= 2)
            return;

        std::vector<index_t> ids2(n);
        for (index_t k = 0; k < n; ++k)
        {
            const index_t mapped = orientSame ? k : (n - 1 - k);
            ids2[k] = boundaryIds2(mapped, 0);
        }

        for (index_t k = 1; k < n - 1; ++k)
        {
            const index_t idxPrev = boundaryIds1(k - 1, 0);
            const index_t idxCur = boundaryIds1(k, 0);
            const index_t idxNext = boundaryIds1(k + 1, 0);

            gsVector<> tangent = patch1.coef(idxNext) - patch1.coef(idxPrev);
            if (tangent.norm() == 0)
                continue;

            gsVector<> normal = curveNormalFromTangent(tangent);
            real_t s = static_cast<real_t>(k) / (n - 1);
            real_t amp = epsAbs * std::sin(static_cast<real_t>(waveNumber) * 2.0 * kPi * s);
            gsVector<> displaced = patch1.coef(idxCur) + amp * normal;

            for (index_t c = 0; c < patch1.geoDim(); ++c)
            {
                patch1.coef(idxCur, c) = displaced(c);
                patch2.coef(ids2[k], c) = displaced(c);
            }
        }
    }

    void perturbInteriorControlPoints(gsGeometry<>& patch,
                                      real_t epsAbs,
                                      std::mt19937& rng)
    {
        if (epsAbs == 0)
            return;

        const gsBasis<>& basis = patch.basis();
        gsMatrix<index_t> bdy = basis.allBoundary();
        std::vector<bool> isBoundary(patch.coefsSize(), false);
        for (index_t i = 0; i < bdy.rows(); ++i)
            isBoundary[bdy(i, 0)] = true;

        std::uniform_real_distribution<real_t> dist(-1.0, 1.0);
        for (index_t i = 0; i < patch.coefsSize(); ++i)
        {
            if (isBoundary[i])
                continue;
            for (index_t d = 0; d < patch.geoDim(); ++d)
                patch.coef(i, d) += epsAbs * dist(rng);
        }
    }

    std::set<index_t> squareBoundaryIndices(index_t coefCount)
    {
        std::set<index_t> ids;
        if (coefCount <= 0)
            return ids;

        const real_t nReal = std::sqrt(static_cast<real_t>(coefCount));
        const index_t n = static_cast<index_t>(std::llround(nReal));
        if (n <= 0 || n * n != coefCount)
            return ids;

        // Row-major indexing on an n x n control net:
        // first row, last row, and start/end index of every intermediate row.
        for (index_t c = 0; c < n; ++c)
        {
            ids.insert(c);
            ids.insert((n - 1) * n + c);
        }

        for (index_t r = 1; r < n - 1; ++r)
        {
            ids.insert(r * n);
            ids.insert(r * n + (n - 1));
        }

        return ids;
    }

    std::vector<bool> boundaryMask(const gsGeometry<>& patch)
    {
        std::vector<bool> isBoundary(patch.coefsSize(), false);
        const std::set<index_t> ids = squareBoundaryIndices(patch.coefsSize());

        if (!ids.empty())
        {
            for (std::set<index_t>::const_iterator it = ids.begin(); it != ids.end(); ++it)
                isBoundary[*it] = true;
            return isBoundary;
        }

        // Fallback for non-square control nets.
        const gsBasis<>& basis = patch.basis();
        gsMatrix<index_t> bdy = basis.allBoundary();
        for (index_t i = 0; i < bdy.rows(); ++i)
            isBoundary[bdy(i, 0)] = true;
        return isBoundary;
    }

    std::string boundaryPointIndicesString(const gsGeometry<>& patch)
    {
        std::set<index_t> uniqueIds;

        const std::set<index_t> squareIds = squareBoundaryIndices(patch.coefsSize());
        if (!squareIds.empty())
        {
            uniqueIds = squareIds;
        }
        else
        {
            const gsBasis<>& basis = patch.basis();
            gsMatrix<index_t> bdy = basis.allBoundary();
            for (index_t i = 0; i < bdy.rows(); ++i)
                uniqueIds.insert(bdy(i, 0));
        }

        std::ostringstream buffer;
        buffer << "[";
        bool first = true;
        for (std::set<index_t>::const_iterator it = uniqueIds.begin(); it != uniqueIds.end(); ++it)
        {
            if (!first)
                buffer << ", ";
            buffer << *it;
            first = false;
        }
        buffer << "]";
        return buffer.str();
    }

    std::vector<bool> fixedPointMask(const gsMultiPatch<>& mp,
                                     index_t patchIndex)
    {
        std::vector<bool> isFixed = boundaryMask(mp.patch(patchIndex));
        const std::vector<boundaryInterface>& interfaces = mp.interfaces();

        for (size_t i = 0; i < interfaces.size(); ++i)
        {
            const boundaryInterface& bi = interfaces[i];
            patchSide ps;
            bool touchesPatch = false;
            if (bi.first().patch == patchIndex)
            {
                ps = bi.first();
                touchesPatch = true;
            }
            else if (bi.second().patch == patchIndex)
            {
                ps = bi.second();
                touchesPatch = true;
            }

            if (!touchesPatch)
                continue;

            gsMatrix<index_t> ids = mp.patch(patchIndex).basis().boundaryOffset(ps.side(), 0);
            for (index_t k = 0; k < ids.rows(); ++k)
                isFixed[ids(k, 0)] = true;
        }

        return isFixed;
    }

    // Returns the physical-space bell center: either the parametric focal point (if
    // focalU/V >= 0) or the control-net centroid.
    gsVector<> computeBellCenter(const gsGeometry<>& patch,
                                 const nonlinearOptions& pOpt)
    {
        if (pOpt.focalU < 0.0 || pOpt.focalV < 0.0)
        {
            // Centroid of control net (default).
            gsVector<> center = gsVector<>::Zero(patch.geoDim());
            for (index_t i = 0; i < patch.coefsSize(); ++i)
                for (index_t d = 0; d < patch.geoDim(); ++d)
                    center(d) += patch.coef(i, d);
            center /= static_cast<real_t>(patch.coefsSize());
            return center;
        }

        // Map normalized [0,1]^2 → parameter support → physical point.
        const gsMatrix<> supp = patch.basis().support();
        gsMatrix<> param(patch.parDim(), 1);
        for (index_t d = 0; d < patch.parDim(); ++d)
        {
            const real_t t = (d == 0) ? pOpt.focalU : pOpt.focalV;
            param(d, 0) = supp(d, 0) + t * (supp(d, 1) - supp(d, 0));
        }
        gsMatrix<> result = patch.eval(param);
        gsVector<> center(patch.geoDim());
        for (index_t d = 0; d < patch.geoDim(); ++d)
            center(d) = result(d, 0);
        return center;
    }

    gsVector<> patchCentroid(const gsGeometry<>& patch)
    {
        gsVector<> center = gsVector<>::Zero(patch.geoDim());
        if (patch.coefsSize() == 0)
            return center;

        for (index_t i = 0; i < patch.coefsSize(); ++i)
        {
            for (index_t d = 0; d < patch.geoDim(); ++d)
                center(d) += patch.coef(i, d);
        }
        center /= static_cast<real_t>(patch.coefsSize());
        return center;
    }

    void restoreBoundaryControlPointsFromReference(const gsMultiPatch<>& referenceMp,
                                                   index_t patchIndex,
                                                   gsGeometry<>& targetPatch)
    {
        const std::vector<bool> isBoundary = boundaryMask(referenceMp.patch(patchIndex));
        const gsGeometry<>& referencePatch = referenceMp.patch(patchIndex);
        GISMO_ENSURE(referencePatch.coefsSize() == targetPatch.coefsSize(),
                     "Reference and target patches have incompatible control-point counts.");

        for (index_t i = 0; i < targetPatch.coefsSize(); ++i)
        {
            if (isBoundary[i])
                targetPatch.coef(i) = referencePatch.coef(i);
        }
    }

    real_t patchControlNetRadius(const gsGeometry<>& patch,
                                 const gsVector<>& center)
    {
        real_t radius = 0;
        for (index_t i = 0; i < patch.coefsSize(); ++i)
        {
            gsVector<> delta = gsVector<>::Zero(patch.geoDim());
            for (index_t d = 0; d < patch.geoDim(); ++d)
                delta(d) = patch.coef(i, d) - center(d);
            radius = std::max(radius, delta.norm());
        }
        return std::max(static_cast<real_t>(1e-12), radius);
    }

    void applyCompactCenterCompression(gsMultiPatch<>& mp,
                                       index_t patchIndex,
                                       const nonlinearOptions& pOpt,
                                       real_t scale)
    {
        GISMO_ENSURE(patchIndex >= 0 && patchIndex < static_cast<index_t>(mp.nPatches()),
                     "Patch index out of range for internal bell distortion.");

        gsGeometry<>& patch = mp.patch(patchIndex);
        const std::vector<bool> isFixed = fixedPointMask(mp, patchIndex);
        const gsVector<> center = computeBellCenter(patch, pOpt);

        const real_t patchRadius = patchControlNetRadius(patch, patchCentroid(patch));
        const real_t effectiveSigma = std::max(static_cast<real_t>(1e-6), pOpt.bellSigma);
        const real_t safeScale = std::max(static_cast<real_t>(0.0), scale);
        const real_t safeStrength = std::max(static_cast<real_t>(0.0), pOpt.strength);

        for (index_t i = 0; i < patch.coefsSize(); ++i)
        {
            if (isFixed[i])
                continue;

            gsVector<> delta = gsVector<>::Zero(patch.geoDim());
            for (index_t d = 0; d < patch.geoDim(); ++d)
                delta(d) = patch.coef(i, d) - center(d);
            const real_t dist = delta.norm();
            if (dist <= 1e-14)
                continue;

            const real_t rho = dist / (effectiveSigma * patchRadius);
            const real_t bell = std::exp(-rho * rho);
            const real_t disp = safeScale * safeStrength * bell;

            const gsVector<> newPos = center + (static_cast<real_t>(1.0) + disp / dist) * delta;
            for (index_t d = 0; d < patch.geoDim(); ++d)
                patch.coef(i, d) = newPos(d);
        }
    }


    // Applies a Gaussian-bell-weighted rotation around a focal point to all free
    // interior control points.  alpha is the peak rotation angle (radians); positive
    // = CCW.  The bell width is controlled by pOpt.bellSigma (fraction of patch
    // radius) and the center by pOpt.focalU/V (< 0 → control-net centroid).
    // A tight bell (sigma ≈ 0.2) concentrated near the first coarsening cell
    // creates angular deformation that LO cannot untwist.
    void applySwirlDistortion(gsMultiPatch<>& mp,
                              index_t patchIndex,
                              real_t alpha,
                              const nonlinearOptions& pOpt)
    {
        GISMO_ENSURE(patchIndex >= 0 && patchIndex < static_cast<index_t>(mp.nPatches()),
                     "Patch index out of range for swirl distortion.");

        gsGeometry<>& patch = mp.patch(patchIndex);
        GISMO_ENSURE(patch.geoDim() == 2, "Swirl distortion is only implemented for 2-D patches.");

        const std::vector<bool> isFixed = fixedPointMask(mp, patchIndex);
        const gsVector<> center = computeBellCenter(patch, pOpt);
        const real_t patchRadius = patchControlNetRadius(patch, patchCentroid(patch));
        const real_t effectiveSigma = std::max(static_cast<real_t>(1e-6), pOpt.bellSigma);

        for (index_t i = 0; i < patch.coefsSize(); ++i)
        {
            if (isFixed[i])
                continue;

            const real_t dx = patch.coef(i, 0) - center(0);
            const real_t dy = patch.coef(i, 1) - center(1);
            const real_t dist = std::sqrt(dx * dx + dy * dy);
            if (dist <= 1e-14)
                continue;

            const real_t rho = dist / (effectiveSigma * patchRadius);
            const real_t bell = std::exp(-rho * rho);
            const real_t angle = alpha * bell;

            const real_t cosA = std::cos(angle);
            const real_t sinA = std::sin(angle);
            patch.coef(i, 0) = center(0) + cosA * dx - sinA * dy;
            patch.coef(i, 1) = center(1) + sinA * dx + cosA * dy;
        }
    }

    // Applies a localized S-fold distortion: opposing v-displacements on either side
    // of a fold axis running through the bell center parallel to the v-direction.
    // Profile: sin(π*du/halfWidth) * exp(-r²/halfWidth²), peak displacement =
    // amplitude * patchRadius.  A halfWidth matching ~1 refined-cell width ensures
    // the fine basis can represent the fold but the next-coarser level cannot.
    void applySfoldDistortion(gsMultiPatch<>& mp,
                              index_t patchIndex,
                              real_t amplitude,
                              real_t foldHalfWidthFrac,
                              const nonlinearOptions& pOpt)
    {
        GISMO_ENSURE(patchIndex >= 0 && patchIndex < static_cast<index_t>(mp.nPatches()),
                     "Patch index out of range for S-fold distortion.");

        gsGeometry<>& patch = mp.patch(patchIndex);
        GISMO_ENSURE(patch.geoDim() == 2, "S-fold distortion is only implemented for 2-D patches.");

        const std::vector<bool> isFixed = fixedPointMask(mp, patchIndex);
        const gsVector<> center = computeBellCenter(patch, pOpt);
        const real_t patchRadius = patchControlNetRadius(patch, patchCentroid(patch));
        const real_t foldHalfWidth =
            std::max(static_cast<real_t>(1e-6), foldHalfWidthFrac) * patchRadius;

        for (index_t i = 0; i < patch.coefsSize(); ++i)
        {
            if (isFixed[i])
                continue;

            const real_t du = patch.coef(i, 0) - center(0);
            const real_t dv = patch.coef(i, 1) - center(1);

            if (std::abs(du) >= foldHalfWidth)
                continue;

            // S-profile: zero at ±halfWidth, ±peak at ±halfWidth/2, zero at center.
            const real_t s_profile = std::sin(kPi * du / foldHalfWidth);
            // Gaussian bell confines the fold to the focal neighborhood.
            const real_t rho = std::sqrt(du * du + dv * dv) / foldHalfWidth;
            const real_t bell = std::exp(-rho * rho);

            patch.coef(i, 1) += amplitude * patchRadius * s_profile * bell;
        }
    }

    gsMatrix<> samplePatchGrid(const gsBasis<>& basis, index_t grid)
    {
        const gsMatrix<> support = basis.support();
        gsVector<> low = support.col(0);
        gsVector<> upp = support.col(1);
        return uniformPointGrid(low, upp, grid * grid);
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

        if (numPatches == 0)
            return result;

        for (index_t patch = 0; patch < numPatches; ++patch)
        {
            const gsBasis<>& basis = mp.patch(patch).basis();
            gsMatrix<> uvPatch = samplePatchGrid(basis, samplesPerDirection);
            if (uvPatch.rows() < 2 || uvPatch.cols() == 0)
                continue;

            const gsMatrix<> support = basis.support();
            const real_t uMin = support(0, 0);
            const real_t uMax = support(0, 1);
            const real_t vMin = support(1, 0);
            const real_t vMax = support(1, 1);
            const real_t uMargin = 0.05 * (uMax - uMin);
            const real_t vMargin = 0.05 * (vMax - vMin);

            for (index_t pt = 0; pt < uvPatch.cols(); ++pt)
            {
                const gsVector<> param = uvPatch.col(pt);
                if (param.rows() < 2)
                    continue;

                if (param(0) <= uMin + uMargin || param(0) >= uMax - uMargin ||
                    param(1) <= vMin + vMargin || param(1) >= vMax - vMargin)
                    continue;

                const real_t det = mp.patch(patch).jacobian(param).determinant();
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
            {
                result.mirrored[patch] = false;
                continue;
            }

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

        bool hasNegative = false;
        for (size_t patch = 0; patch < outputMp.nPatches(); ++patch)
        {
            const gsBasis<>& basis = outputMp.patch(patch).basis();
            gsMatrix<> uv = samplePatchGrid(basis, fitGrid);
            real_t minRawDet = std::numeric_limits<real_t>::max();
            real_t minOrientedDet = std::numeric_limits<real_t>::max();
            index_t negCount = 0;
            const bool patchMirrored =
                (patch < inputMirrorReport.mirrored.size()) ? inputMirrorReport.mirrored[patch] : false;

            for (index_t j = 0; j < uv.cols(); ++j)
            {
                const real_t det = outputMp.patch(patch).jacobian(uv.col(j)).determinant();
                minRawDet = std::min(minRawDet, det);

                // Mirroring only changes the sign convention used for diagnostics.
                const real_t orientedDet = patchMirrored ? -det : det;
                minOrientedDet = std::min(minOrientedDet, orientedDet);
                if (orientedDet < detThreshold)
                    ++negCount;
            }

            if (negCount > 0)
                hasNegative = true;

            log << "Patch " << patch << ": mirroredBaseline=" << (patchMirrored ? "true" : "false")
                << ", min raw det=" << minRawDet
                << ", min oriented det=" << minOrientedDet
                << ", count below threshold=" << negCount << "\n";
            if (logToInfo)
                gsInfo << "Patch " << patch << ": mirroredBaseline=" << (patchMirrored ? "true" : "false")
                       << ", min raw det=" << minRawDet
                       << ", min oriented det=" << minOrientedDet
                       << ", count below threshold=" << negCount << "\n";
        }

        return hasNegative;
    }

    void enforceC0AcrossInterfaces(gsMultiPatch<>& mp,
                                   std::ostream* log = 0)
    {
        if (log)
            *log << "\nManual C0 enforcement across interfaces\n";

        for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
        {
            const index_t p1 = it->first().patch;
            const index_t p2 = it->second().patch;

            gsGeometry<>& g1 = mp.patch(p1);
            gsGeometry<>& g2 = mp.patch(p2);

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
                gsGeometry<>::uPtr c1 = g1.boundary(s1);
                gsGeometry<>::uPtr c2 = g2.boundary(s2);
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
                if (log)
                    *log << "Interface size mismatch between patches " << p1 << " and " << p2 << "\n";
                continue;
            }

            const index_t geoDim = g1.geoDim();
            for (size_t k = 0; k < ids1.size(); ++k)
            {
                for (index_t c = 0; c < geoDim; ++c)
                    g2.coef(ids2[k], c) = g1.coef(ids1[k], c);
            }
        }
    }

    real_t minOrientedDeterminant(const gsMultiPatch<>& inputMp,
                                  const gsMultiPatch<>& outputMp,
                                  index_t fitGrid)
    {
        const MirroredCheckResult inputMirrorReport =
            computeMirroredCheckPerPatch(inputMp, fitGrid, 1e-12);

        real_t minOrientedDet = std::numeric_limits<real_t>::max();
        for (size_t patch = 0; patch < outputMp.nPatches(); ++patch)
        {
            const gsBasis<>& basis = outputMp.patch(patch).basis();
            gsMatrix<> uv = samplePatchGrid(basis, fitGrid);
            const bool patchMirrored =
                (patch < inputMirrorReport.mirrored.size()) ? inputMirrorReport.mirrored[patch] : false;

            for (index_t j = 0; j < uv.cols(); ++j)
            {
                const real_t det = outputMp.patch(patch).jacobian(uv.col(j)).determinant();
                const real_t orientedDet = patchMirrored ? -det : det;
                minOrientedDet = std::min(minOrientedDet, orientedDet);
            }
        }

        return minOrientedDet;
    }

    bool hasC0Violation(const gsMultiPatch<>& mp,
                        real_t tolerance,
                        real_t& maxGapOut,
                        index_t& violatingInterfacesOut)
    {
        maxGapOut = 0;
        violatingInterfacesOut = 0;

        const std::vector<boundaryInterface>& interfaces = mp.interfaces();
        for (size_t i = 0; i < interfaces.size(); ++i)
        {
            const boundaryInterface& bi = interfaces[i];
            const patchSide& ps1 = bi.first();
            const patchSide& ps2 = bi.second();

            const gsGeometry<>& g1 = mp.patch(ps1.patch);
            const gsGeometry<>& g2 = mp.patch(ps2.patch);
            gsMatrix<index_t> ids1 = g1.basis().boundaryOffset(ps1.side(), 0);
            gsMatrix<index_t> ids2 = g2.basis().boundaryOffset(ps2.side(), 0);

            if (ids1.rows() != ids2.rows())
            {
                ++violatingInterfacesOut;
                maxGapOut = std::numeric_limits<real_t>::infinity();
                continue;
            }

            const bool orientSame = bi.dirOrientation(ps1, 1 - ps1.direction());
            real_t interfaceMaxGap = 0;
            for (index_t k = 0; k < ids1.rows(); ++k)
            {
                const index_t mapped = orientSame ? k : (ids1.rows() - 1 - k);
                const real_t gap = (g1.coef(ids1(k, 0)) - g2.coef(ids2(mapped, 0))).norm();
                interfaceMaxGap = std::max(interfaceMaxGap, gap);
                maxGapOut = std::max(maxGapOut, gap);
            }

            if (interfaceMaxGap > tolerance)
                ++violatingInterfacesOut;
        }

        return violatingInterfacesOut > 0;
    }

    bool hasBoundaryPointDrift(const gsMultiPatch<>& referenceMp,
                               const gsMultiPatch<>& candidateMp,
                               real_t tolerance,
                               real_t& maxDriftOut,
                               index_t& movedBoundaryPointCountOut)
    {
        maxDriftOut = 0;
        movedBoundaryPointCountOut = 0;

        const index_t numPatches = static_cast<index_t>(referenceMp.nPatches());
        GISMO_ENSURE(numPatches == static_cast<index_t>(candidateMp.nPatches()),
                     "Reference and candidate multipatches have different patch counts.");

        for (index_t p = 0; p < numPatches; ++p)
        {
            const std::vector<bool> isBoundary = boundaryMask(referenceMp.patch(p));
            const gsGeometry<>& refPatch = referenceMp.patch(p);
            const gsGeometry<>& candPatch = candidateMp.patch(p);

            GISMO_ENSURE(refPatch.coefsSize() == candPatch.coefsSize(),
                         "Reference and candidate patches have incompatible control-point counts.");

            for (index_t i = 0; i < candPatch.coefsSize(); ++i)
            {
                if (!isBoundary[i])
                    continue;

                const real_t drift = (candPatch.coef(i) - refPatch.coef(i)).norm();
                maxDriftOut = std::max(maxDriftOut, drift);
                if (drift > tolerance)
                    ++movedBoundaryPointCountOut;
            }
        }

        return movedBoundaryPointCountOut > 0;
    }

    gsGeometry<>::uPtr fitPatchWithInteriorAmplitude(
        const gsMultiPatch<>& mp,
        index_t patchIndex,
        real_t epsilonIntDist,
        const std::map<patchSide, gsGeometry<>::uPtr>& constraintCurves,
        index_t fitGrid,
        int seed,
        bool doInterior)
    {
        gsGeometry<>::uPtr target = mp.patch(patchIndex).clone();
        if (doInterior)
        {
            std::mt19937 patchRng(static_cast<unsigned>(seed + static_cast<int>(patchIndex)));
            perturbInteriorControlPoints(*target, epsilonIntDist, patchRng);
        }

        const gsBasis<>& basis = mp.patch(patchIndex).basis();
        gsMatrix<> uv = samplePatchGrid(basis, fitGrid);
        gsMatrix<> xy = target->eval(uv);

        gsFitting<> fit(uv, xy, const_cast<gsBasis<>&>(basis));

        std::vector<boxSide> fixedSides;
        std::vector<gsGeometry<>*> fixedCurves;
        for (std::map<patchSide, gsGeometry<>::uPtr>::const_iterator it = constraintCurves.begin();
             it != constraintCurves.end(); ++it)
        {
            if (it->first.patch == patchIndex)
            {
                fixedSides.push_back(it->first.side());
                fixedCurves.push_back(it->second.get());
            }
        }

        if (!fixedSides.empty())
            fit.setConstraints(fixedSides, fixedCurves);

        fit.compute();
        gsGeometry<>::uPtr fittedPatch(fit.result()->clone());
        // Keep boundary control points exactly on the reference patch.
        restoreBoundaryControlPointsFromReference(mp, patchIndex, *fittedPatch);
        return fittedPatch;
    }

    real_t findPatchInteriorAmplitudeByBisection(
        const gsMultiPatch<>& originalMp,
        const std::vector<gsGeometry<>::uPtr>& selectedPatches,
        index_t patchIndex,
        const std::map<patchSide, gsGeometry<>::uPtr>& constraintCurves,
        index_t fitGrid,
        int seed,
        bool doInterior,
        real_t lowerBound,
        real_t upperBound,
        real_t intervalTolerance,
        real_t detTolerance,
        real_t c0Tolerance,
        std::ostream& log,
        bool logToInfo,
        bool& foundRegular)
    {
        foundRegular = false;

        auto isRegular = [&](real_t epsilonValue,
                             real_t& minDetOut,
                             real_t& maxC0GapOut,
                             gsGeometry<>::uPtr* fittedPatchOut) -> bool
        {
            gsGeometry<>::uPtr fittedPatch = fitPatchWithInteriorAmplitude(
                originalMp, patchIndex, epsilonValue, constraintCurves, fitGrid, seed, doInterior);

            gsMultiPatch<> trialMp;
            for (index_t p = 0; p < static_cast<index_t>(originalMp.nPatches()); ++p)
            {
                if (p < patchIndex && p < static_cast<index_t>(selectedPatches.size()) && selectedPatches[p])
                    trialMp.addPatch(*selectedPatches[p]);
                else if (p == patchIndex)
                    trialMp.addPatch(*fittedPatch);
                else
                    trialMp.addPatch(originalMp.patch(p));
            }

            for (index_t p = 0; p < static_cast<index_t>(trialMp.nPatches()); ++p)
                restoreBoundaryControlPointsFromReference(originalMp, p, trialMp.patch(p));

            trialMp.computeTopology();
            minDetOut = minOrientedDeterminant(originalMp, trialMp, fitGrid);
            index_t violatingInterfaces = 0;
            const bool c0Broken = hasC0Violation(trialMp, c0Tolerance, maxC0GapOut, violatingInterfaces);
            if (fittedPatchOut)
                *fittedPatchOut = std::move(fittedPatch);
            return minDetOut > detTolerance && !c0Broken;
        };

        real_t low = lowerBound;
        real_t high = upperBound;
        real_t lowMinDet = -std::numeric_limits<real_t>::max();
        real_t highMinDet = -std::numeric_limits<real_t>::max();
        real_t lowMaxC0Gap = 0;
        real_t highMaxC0Gap = 0;

        const bool lowRegular = isRegular(low, lowMinDet, lowMaxC0Gap, 0);
        const bool highRegular = isRegular(high, highMinDet, highMaxC0Gap, 0);

        if (!lowRegular)
        {
            log << "Patch " << patchIndex << ": epsilonIntDist lower bound " << low
                << " is already irregular (min oriented det=" << lowMinDet
                << ", max C0 gap=" << lowMaxC0Gap << ")\n";
            if (logToInfo)
                gsInfo << "Patch " << patchIndex << ": epsilonIntDist lower bound " << low
                       << " is already irregular (min oriented det=" << lowMinDet
                       << ", max C0 gap=" << lowMaxC0Gap << ")\n";
            real_t zeroMinDet = -std::numeric_limits<real_t>::max();
            real_t zeroMaxC0Gap = 0;
            if (isRegular(0, zeroMinDet, zeroMaxC0Gap, 0))
            {
                foundRegular = true;
                return 0;
            }
            return 0;
        }

        foundRegular = true;
        if (highRegular)
            return high;

        while (high - low > intervalTolerance)
        {
            const real_t mid = 0.5 * (low + high);
            real_t midMinDet = -std::numeric_limits<real_t>::max();
            real_t midMaxC0Gap = 0;
            if (isRegular(mid, midMinDet, midMaxC0Gap, 0))
                low = mid;
            else
                high = mid;
        }

        return low;
    }

    template <typename ApplyFn>
    real_t tuneParameterByBisection(
        const gsMultiPatch<>& inputMp,
        const gsMultiPatch<>& baseMp,
        index_t patchIndex,
        real_t lowerBound,
        real_t upperBound,
        real_t intervalTolerance,
        real_t determinantTolerance,
        real_t c0Tolerance,
        index_t fitGrid,
        const std::string& parameterName,
        const ApplyFn& applyFn,
        std::ostream& log,
        bool logToInfo,
        bool& foundRegular)
    {
        foundRegular = false;

        auto isRegular = [&](real_t value, real_t& minDetOut, real_t& maxC0GapOut) -> bool
        {
            gsMultiPatch<> trialMp = baseMp;
            applyFn(trialMp, patchIndex, value);
            restoreBoundaryControlPointsFromReference(inputMp, patchIndex, trialMp.patch(patchIndex));
            trialMp.computeTopology();
            minDetOut = minOrientedDeterminant(inputMp, trialMp, fitGrid);
            index_t violatingInterfaces = 0;
            const bool c0Broken = hasC0Violation(trialMp, c0Tolerance, maxC0GapOut, violatingInterfaces);
            return minDetOut > determinantTolerance && !c0Broken;
        };

        real_t low = lowerBound;
        real_t high = upperBound;
        real_t lowMinDet = -std::numeric_limits<real_t>::max();
        real_t highMinDet = -std::numeric_limits<real_t>::max();
        real_t lowMaxC0Gap = 0;
        real_t highMaxC0Gap = 0;

        const bool lowRegular = isRegular(low, lowMinDet, lowMaxC0Gap);
        const bool highRegular = isRegular(high, highMinDet, highMaxC0Gap);

        if (!lowRegular)
        {
            log << "Patch " << patchIndex << ": " << parameterName
                << " lower bound " << low
                << " is irregular (min oriented det=" << lowMinDet
                << ", max C0 gap=" << lowMaxC0Gap << ")\n";
            if (logToInfo)
                gsInfo << "Patch " << patchIndex << ": " << parameterName
                       << " lower bound " << low
                       << " is irregular (min oriented det=" << lowMinDet
                       << ", max C0 gap=" << lowMaxC0Gap << ")\n";
            foundRegular = false;
            return 0.0;
        }

        foundRegular = true;
        if (highRegular)
            return high;

        while (high - low > intervalTolerance)
        {
            const real_t mid = static_cast<real_t>(0.5) * (low + high);
            real_t midMinDet = -std::numeric_limits<real_t>::max();
            real_t midMaxC0Gap = 0;
            if (isRegular(mid, midMinDet, midMaxC0Gap))
                low = mid;
            else
                high = mid;
        }

        return low;
    }

    bool nudgeInterfaceControlPoint(gsGeometry<>& patch,
                                    const gsMultiPatch<>& mp,
                                    index_t patchIndex,
                                    real_t moveAbs)
    {
        if (moveAbs == 0)
            return false;

        const std::vector<boundaryInterface>& interfaces = mp.interfaces();
        bool movedAny = false;
        for (size_t i = 0; i < interfaces.size(); ++i)
        {
            const boundaryInterface& bi = interfaces[i];
            patchSide ps;
            if (bi.first().patch == patchIndex)
                ps = bi.first();
            else if (bi.second().patch == patchIndex)
                ps = bi.second();
            else
                continue;

            gsMatrix<index_t> boundaryIds = patch.basis().boundaryOffset(ps.side(), 0);
            const index_t n = boundaryIds.rows();
            if (n <= 2)
                continue;

            for (index_t k = 1; k < n - 1; ++k)
            {
                const index_t idxPrev = boundaryIds(k - 1, 0);
                const index_t idxCur = boundaryIds(k, 0);
                const index_t idxNext = boundaryIds(k + 1, 0);

                gsVector<> tangent = patch.coef(idxNext) - patch.coef(idxPrev);
                if (tangent.norm() == 0)
                    continue;

                gsVector<> normal = curveNormalFromTangent(tangent);
                const real_t sign = (k % 2 == 1) ? static_cast<real_t>(-1) : static_cast<real_t>(1);
                for (index_t c = 0; c < patch.geoDim(); ++c)
                    patch.coef(idxCur, c) += sign * moveAbs * normal(c);

                movedAny = true;
            }
        }

        return movedAny;
    }

    bool endsWith(const std::string& value, const std::string& suffix)
    {
        if (value.size() < suffix.size())
            return false;
        return value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
    }

    std::string resolveInputFilename(int argc, char* argv[], const std::string& parsed)
    {
        if (gsFileManager::fileExists(parsed))
            return parsed;

        for (int i = 1; i < argc; ++i)
        {
            std::string arg = argv[i];
            if (endsWith(arg, ".xml") && gsFileManager::fileExists(arg))
                return arg;
        }

        return parsed;
    }

    std::string makeNLOOutputName(const std::string& inputPath)
    {
        std::string normalized = inputPath;
        for (size_t i = 0; i < normalized.size(); ++i)
        {
            if (normalized[i] == '\\')
                normalized[i] = '/';
        }

        const size_t slashPos = normalized.find_last_of('/');
        std::string fileName = (slashPos == std::string::npos)
                                   ? normalized
                                   : normalized.substr(slashPos + 1);

        const size_t dotPos = fileName.find_last_of('.');
        const std::string stem = (dotPos == std::string::npos)
                                     ? fileName
                                     : fileName.substr(0, dotPos);

        return stem + "_NLO.xml";
    }
}

int main(int argc, char* argv[])
{
    std::string filename = "generatedMPs\\mask_approximation_fine_L3.xml";

    Options opt;
    nonlinearOptions nlo;

    gsCmdLine cmd("Distort multipatch geometry while preserving exact C0 continuity.");
    cmd.addPlainString("filename", "Input multipatch XML file.", filename);
    cmd.addString("o", "output", "Output distorted multipatch XML.", opt.outputFile);
    cmd.addString("L", "log", "Log file.", opt.logFile);
    cmd.addInt("g", "fit-grid", "Grid size per patch for Jacobian sampling (default 20).", opt.fitGrid);
    cmd.addReal("t", "bisect-tol", "Bisection interval tolerance for strength/alpha (default 1e-2).", opt.bisectTol);
    cmd.addReal("j", "det-threshold", "Report if jacobian determinant < threshold.", opt.detThreshold);
    cmd.addInt("p", "stress-patch", "Patch index for stress test (-1: all patches with bisection).", nlo.patchIndex);
    cmd.addReal("a", "stress-strength", "Single outward distortion strength.", nlo.strength);
    cmd.addReal("M", "stress-min-det", "Minimum oriented Jacobian determinant used for reporting.", nlo.minOrientedDet);
    cmd.addReal("w", "swirl-alpha", "Peak swirl rotation angle in radians (0 = disabled). Positive = CCW.", nlo.swirlAlpha);
    cmd.addReal("s", "bell-sigma", "Gaussian bell width as fraction of patch radius (default 1.0, try 0.2 for concentrated distortion).", nlo.bellSigma);
    cmd.addReal("u", "focal-u", "Normalized [0,1] u-parameter of bell center (default -1 = control-net centroid).", nlo.focalU);
    cmd.addReal("v", "focal-v", "Normalized [0,1] v-parameter of bell center (default -1 = control-net centroid).", nlo.focalV);
    cmd.addReal("f", "sfold-amplitude", "Peak S-fold displacement as fraction of patch radius (0 = disabled).", nlo.sfoldAmplitude);
    cmd.addReal("W", "sfold-half-width", "S-fold zone half-width as fraction of patch radius (default 0.15).", nlo.sfoldHalfWidth);

    // Pre-parse per-patch overrides before gsCmdLine sees them (it would fail on unknown args).
    //   --patch-det <p>:<val>       e.g. --patch-det 7:1e-5
    //   --patch-focal <p>:<u>:<v>   e.g. --patch-focal 7:0.0:0.5
    // These allow targeted near-singular construction: push one patch to det≈ε
    // while other patches use the global --stress-min-det floor.
    // We strip these args from argv before passing to gsCmdLine.
    std::vector<char*> filteredArgv;
    filteredArgv.push_back(argv[0]);
    for (int ai = 1; ai < argc; ++ai)
    {
        const std::string arg = argv[ai];
        if ((arg == "--patch-det") && ai + 1 < argc)
        {
            const std::string spec = argv[++ai];
            const size_t c = spec.find(':');
            if (c != std::string::npos)
            {
                const index_t p = static_cast<index_t>(std::stoi(spec.substr(0, c)));
                const real_t  v = static_cast<real_t>(std::stod(spec.substr(c + 1)));
                opt.perPatchMinDet[p] = v;
            }
            else
            {
                gsWarn << "Ignoring malformed --patch-det spec (expected p:val): " << spec << "\n";
            }
        }
        else if ((arg == "--patch-focal") && ai + 1 < argc)
        {
            const std::string spec = argv[++ai];
            const size_t c1 = spec.find(':');
            const size_t c2 = (c1 != std::string::npos) ? spec.find(':', c1 + 1) : std::string::npos;
            if (c1 != std::string::npos && c2 != std::string::npos)
            {
                const index_t p = static_cast<index_t>(std::stoi(spec.substr(0, c1)));
                const real_t  u = static_cast<real_t>(std::stod(spec.substr(c1 + 1, c2 - c1 - 1)));
                const real_t  vv = static_cast<real_t>(std::stod(spec.substr(c2 + 1)));
                opt.perPatchFocalU[p] = u;
                opt.perPatchFocalV[p] = vv;
            }
            else
            {
                gsWarn << "Ignoring malformed --patch-focal spec (expected p:u:v): " << spec << "\n";
            }
        }
        else
        {
            filteredArgv.push_back(argv[ai]);
        }
    }
    int filteredArgc = static_cast<int>(filteredArgv.size());

    try { cmd.getValues(filteredArgc, filteredArgv.data()); } catch (int rv) { return rv; }

    std::string resolvedInput = resolveInputFilename(argc, argv, filename);
    if (opt.outputFile.empty())
        opt.outputFile = makeNLOOutputName(resolvedInput);

    if (resolvedInput != filename)
        gsInfo << "Using raw argv input path: " << resolvedInput << "\n";

    gsMultiPatch<>::uPtr mp = gsReadFile<>(resolvedInput);
    if (!mp)
    {
        gsInfo << "Failed to read input file: " << resolvedInput << "\n";
        gsInfo << "Hint: use forward slashes (C:/path/file.xml) or double backslashes (C:\\path\\file.xml).\n";
        return EXIT_FAILURE;
    }

    mp->computeTopology();

    std::ofstream log(opt.logFile.c_str());
    log << "Multipatch distortion log\n";
    log << "Input (parsed): " << filename << "\n";
    log << "Input (resolved): " << resolvedInput << "\n";
    log << "Output: " << opt.outputFile << "\n";
    log << "fit-grid: " << opt.fitGrid << "\n";
    log << "c0-tolerance: " << kC0Tolerance << "\n";
    log << "stress-patch: " << nlo.patchIndex << "\n";
    log << "stress-strength: " << nlo.strength << "\n";
    log << "stress-min-oriented-det: " << nlo.minOrientedDet << "\n\n";
    log << "stress-test-mode: internal bell distortion with fixed per-patch strengths (2x logfile)\n";
    log << "swirl-alpha: " << nlo.swirlAlpha << "\n";
    log << "bell-sigma: " << nlo.bellSigma << "\n";
    log << "focal-u: " << nlo.focalU << "\n";
    log << "focal-v: " << nlo.focalV << "\n";
    log << "sfold-amplitude: " << nlo.sfoldAmplitude << "\n";
    log << "sfold-half-width: " << nlo.sfoldHalfWidth << "\n";
    if (!opt.perPatchMinDet.empty())
    {
        log << "per-patch-min-det overrides:";
        for (const auto& kv : opt.perPatchMinDet)
            log << " p" << kv.first << "=" << kv.second;
        log << "\n";
        gsInfo << "[per-patch] min-det overrides:";
        for (const auto& kv : opt.perPatchMinDet)
            gsInfo << " p" << kv.first << "=" << kv.second;
        gsInfo << "\n";
    }
    if (!opt.perPatchFocalU.empty() || !opt.perPatchFocalV.empty())
    {
        log << "per-patch-focal overrides:";
        for (const auto& kv : opt.perPatchFocalU)
        {
            const real_t vFocal = opt.perPatchFocalV.count(kv.first)
                                  ? opt.perPatchFocalV.at(kv.first)
                                  : nlo.focalV;
            log << " p" << kv.first << "=(" << kv.second << "," << vFocal << ")";
        }
        log << "\n";
        gsInfo << "[per-patch] focal overrides applied\n";
    }
    log << "\n";

    std::map<patchSide, gsGeometry<>::uPtr> constraintCurves;

    const std::vector<boundaryInterface>& interfaces = mp->interfaces();
    log << "Interfaces: " << interfaces.size() << "\n";
    for (size_t i = 0; i < interfaces.size(); ++i)
    {
        const boundaryInterface& bi = interfaces[i];
        const patchSide& ps1 = bi.first();
        const patchSide& ps2 = bi.second();
        constraintCurves[ps1] = mp->patch(ps1.patch).boundary(ps1.side());
        constraintCurves[ps2] = mp->patch(ps2.patch).boundary(ps2.side());

        log << "Interface " << i << ": patches " << ps1.patch << "-" << ps2.patch
            << " (fit first, distort later)\n";
    }

    const gsBoxTopology::bContainer& boundaries = mp->boundaries();
    log << "Boundaries: " << boundaries.size() << "\n";
    for (size_t i = 0; i < boundaries.size(); ++i)
    {
        const patchSide& ps = boundaries[i];
        if (constraintCurves.find(ps) != constraintCurves.end())
            continue;

        constraintCurves[ps] = mp->patch(ps.patch).boundary(ps.side());

        log << "Boundary patch " << ps.patch << " side " << ps.side().index()
            << " (fit first, distort later)\n";
    }

    gsMultiPatch<> newmp;
    std::vector<gsGeometry<>::uPtr> fittedPatches(mp->nPatches());

    log << "\nBoundary control-point indices per patch:\n";
    for (size_t p = 0; p < mp->nPatches(); ++p)
    {
        const std::string boundaryIndices = boundaryPointIndicesString(mp->patch(static_cast<index_t>(p)));
        log << "Patch " << p << ": " << boundaryIndices << "\n";
        gsInfo << "Patch " << p << " boundary control-point indices: "
               << boundaryIndices << "\n";
    }

    log << "\nPatch fitting:\n";
    for (size_t p = 0; p < mp->nPatches(); ++p)
    {
        gsGeometry<>::uPtr fitted = fitPatchWithInteriorAmplitude(
            *mp,
            static_cast<index_t>(p),
            0,
            constraintCurves,
            opt.fitGrid,
            0,
            false);
        fittedPatches[p] = gsGeometry<>::uPtr(fitted->clone());
        newmp.addPatch(*fitted);
    }

    newmp.computeTopology();

    const bool hasNegativeAfterFitting =
        runMirrorAwareJacobianCheck(*mp, newmp, opt.fitGrid, opt.detThreshold,
                                    "after fitting, before control-point moves", log, true);
    GISMO_UNUSED(hasNegativeAfterFitting);

    // ===== Phase 0.5: S-fold distortion on the clean fitted geometry =====
    // S-fold: opposing v-displacements left/right of a fold axis at the bell center.
    // The fold zone half-width is set to ~1 refined-cell width so the fine basis can
    // represent the fold but the next-coarser level cannot — forcing LS projection to
    // produce an irregular geometry that LO cannot fix (only NLO can).
    if (nlo.sfoldAmplitude != 0.0)
    {
        log << "\n========== S-FOLD DISTORTION (amplitude=" << nlo.sfoldAmplitude
            << ", half-width=" << nlo.sfoldHalfWidth << ") ==========\n";
        gsInfo << "Applying S-fold distortion, amplitude=" << nlo.sfoldAmplitude
               << ", half-width=" << nlo.sfoldHalfWidth << "\n";

        // Collect all patches that share a patch-to-patch interface.
        // S-fold must NOT distort these: testBoundaryAssembly samples 10 points
        // along each interface, so any distortion there drives featureError above epsilon_f.
        std::set<index_t> interfacePatches;
        {
            const std::vector<boundaryInterface>& ifaces = newmp.interfaces();
            for (const auto& iface : ifaces)
            {
                interfacePatches.insert(iface.first().patch);
                interfacePatches.insert(iface.second().patch);
            }
        }
        if (!interfacePatches.empty())
        {
            log << "Skipping S-fold on interface-adjacent patches:";
            for (index_t p : interfacePatches) log << " " << p;
            log << "\n";
        }

        std::vector<real_t> appliedSfold(newmp.nPatches(), 0.0);

        for (index_t pIdx = 0; pIdx < static_cast<index_t>(newmp.nPatches()); ++pIdx)
        {
            if (nlo.patchIndex >= 0 && pIdx != nlo.patchIndex)
                continue;
            if (interfacePatches.count(pIdx))
                continue;

            // Apply per-patch overrides for near-singular targeted construction.
            const real_t patchMinDetS = opt.perPatchMinDet.count(pIdx)
                ? opt.perPatchMinDet.at(pIdx) : nlo.minOrientedDet;
            nonlinearOptions sfoldNlo = nlo;
            if (opt.perPatchFocalU.count(pIdx)) sfoldNlo.focalU = opt.perPatchFocalU.at(pIdx);
            if (opt.perPatchFocalV.count(pIdx)) sfoldNlo.focalV = opt.perPatchFocalV.at(pIdx);
            if (patchMinDetS != nlo.minOrientedDet)
                gsInfo << "Patch " << pIdx << ": using per-patch minDet=" << patchMinDetS << " for s-fold\n";

            bool foundRegular = false;
            const real_t safeAmplitude = tuneParameterByBisection(
                *mp,
                newmp,
                pIdx,
                0.0,
                nlo.sfoldAmplitude,
                opt.bisectTol,
                patchMinDetS,
                kC0Tolerance,
                opt.fitGrid,
                "sfold-amplitude",
                [sfoldNlo](gsMultiPatch<>& trialMp, index_t patchIdx, real_t amp)
                {
                    applySfoldDistortion(trialMp, patchIdx, amp, sfoldNlo.sfoldHalfWidth, sfoldNlo);
                },
                log,
                true,
                foundRegular);

            if (!foundRegular)
            {
                log << "Patch " << pIdx << ": s-fold skipped (irregular at amplitude=0)\n";
                gsInfo << "Patch " << pIdx << ": s-fold skipped\n";
                continue;
            }

            applySfoldDistortion(newmp, pIdx, safeAmplitude, sfoldNlo.sfoldHalfWidth, sfoldNlo);
            enforceC0AcrossInterfaces(newmp);
            newmp.computeTopology();
            appliedSfold[pIdx] = safeAmplitude;

            log << "Patch " << pIdx << ": s-fold applied, safe amplitude=" << safeAmplitude << "\n";
            gsInfo << "Patch " << pIdx << ": s-fold safe amplitude=" << safeAmplitude << "\n";
        }

        log << "\nApplied s-fold amplitudes:\n";
        for (index_t p = 0; p < static_cast<index_t>(newmp.nPatches()); ++p)
        {
            if (appliedSfold[p] != 0.0)
                log << "  patch " << p << " -> " << appliedSfold[p] << "\n";
        }

        // Per-patch det summary after S-fold — critical for verifying near-singular patches.
        runMirrorAwareJacobianCheck(*mp, newmp, opt.fitGrid, opt.detThreshold,
                                    "after s-fold", log, true);
    }

    // ===== Phase 1: Swirl distortion on the clean fitted geometry =====
    // Applied first while min det is still large (34–193), giving bisection room
    // to find meaningful rotation angles before radial compression eats the budget.
    if (nlo.swirlAlpha != 0.0)
    {
        log << "\n========== SWIRL DISTORTION (alpha=" << nlo.swirlAlpha << ") ==========\n";
        gsInfo << "Applying swirl distortion, alpha=" << nlo.swirlAlpha << "\n";

        std::vector<real_t> appliedSwirl(newmp.nPatches(), 0.0);

        for (index_t pIdx = 0; pIdx < static_cast<index_t>(newmp.nPatches()); ++pIdx)
        {
            if (nlo.patchIndex >= 0 && pIdx != nlo.patchIndex)
                continue;

            // Apply per-patch overrides for near-singular targeted construction.
            const real_t patchMinDetW = opt.perPatchMinDet.count(pIdx)
                ? opt.perPatchMinDet.at(pIdx) : nlo.minOrientedDet;
            nonlinearOptions swirlNlo = nlo;
            if (opt.perPatchFocalU.count(pIdx)) swirlNlo.focalU = opt.perPatchFocalU.at(pIdx);
            if (opt.perPatchFocalV.count(pIdx)) swirlNlo.focalV = opt.perPatchFocalV.at(pIdx);
            if (patchMinDetW != nlo.minOrientedDet)
                gsInfo << "Patch " << pIdx << ": using per-patch minDet=" << patchMinDetW << " for swirl\n";

            bool foundRegular = false;
            const real_t safeAlpha = tuneParameterByBisection(
                *mp,
                newmp,
                pIdx,
                0.0,
                nlo.swirlAlpha,
                opt.bisectTol,
                patchMinDetW,
                kC0Tolerance,
                opt.fitGrid,
                "swirl-alpha",
                [swirlNlo](gsMultiPatch<>& trialMp, index_t patchIdx, real_t alpha)
                {
                    applySwirlDistortion(trialMp, patchIdx, alpha, swirlNlo);
                },
                log,
                true,
                foundRegular);

            if (!foundRegular)
            {
                log << "Patch " << pIdx << ": swirl skipped (irregular at alpha=0)\n";
                gsInfo << "Patch " << pIdx << ": swirl skipped\n";
                continue;
            }

            applySwirlDistortion(newmp, pIdx, safeAlpha, swirlNlo);
            enforceC0AcrossInterfaces(newmp);
            newmp.computeTopology();
            appliedSwirl[pIdx] = safeAlpha;

            log << "Patch " << pIdx << ": swirl applied, safe alpha=" << safeAlpha << "\n";
            gsInfo << "Patch " << pIdx << ": swirl safe alpha=" << safeAlpha << "\n";
        }

        log << "\nApplied swirl alphas:\n";
        for (index_t p = 0; p < static_cast<index_t>(newmp.nPatches()); ++p)
        {
            if (appliedSwirl[p] != 0.0)
                log << "  patch " << p << " -> " << appliedSwirl[p] << "\n";
        }
    }

    // ===== Phase 2: Radial distortion on top of the swirled geometry =====
    log << "\n========== RADIAL DISTORTION (bisection, max strength=" << nlo.strength << ") ==========\n";

    std::vector<real_t> chosenStrength(newmp.nPatches(), 0.0);
    std::vector<bool> tunedPatch(newmp.nPatches(), false);

    for (index_t pIdx = 0; pIdx < static_cast<index_t>(newmp.nPatches()); ++pIdx)
    {
        if (nlo.patchIndex >= 0 && pIdx != nlo.patchIndex)
            continue;

        // Apply per-patch overrides for near-singular targeted construction.
        const real_t patchMinDetR = opt.perPatchMinDet.count(pIdx)
            ? opt.perPatchMinDet.at(pIdx) : nlo.minOrientedDet;
        nonlinearOptions radialNlo = nlo;
        if (opt.perPatchFocalU.count(pIdx)) radialNlo.focalU = opt.perPatchFocalU.at(pIdx);
        if (opt.perPatchFocalV.count(pIdx)) radialNlo.focalV = opt.perPatchFocalV.at(pIdx);
        if (patchMinDetR != nlo.minOrientedDet)
            gsInfo << "Patch " << pIdx << ": using per-patch minDet=" << patchMinDetR << " for radial\n";

        bool foundRegular = false;
        const nonlinearOptions patchOpt = radialNlo;
        const real_t safeStrength = tuneParameterByBisection(
            *mp,
            newmp,
            pIdx,
            0.0,
            nlo.strength,
            opt.bisectTol,
            patchMinDetR,
            kC0Tolerance,
            opt.fitGrid,
            "radial-strength",
            [&patchOpt](gsMultiPatch<>& trialMp, index_t patchIdx, real_t s)
            {
                nonlinearOptions applyOpt = patchOpt;
                applyOpt.patchIndex = patchIdx;
                applyOpt.strength = s;
                applyCompactCenterCompression(trialMp, patchIdx, applyOpt, 1.0);
            },
            log,
            true,
            foundRegular);

        if (!foundRegular)
        {
            log << "Patch " << pIdx << ": radial skipped (irregular at strength=0)\n";
            gsInfo << "Patch " << pIdx << ": radial skipped\n";
            continue;
        }

        nonlinearOptions applyOpt = radialNlo;
        applyOpt.patchIndex = pIdx;
        applyOpt.strength = safeStrength;
        applyCompactCenterCompression(newmp, pIdx, applyOpt, 1.0);
        newmp.computeTopology();
        chosenStrength[pIdx] = safeStrength;
        tunedPatch[pIdx] = true;

        log << "Patch " << pIdx << ": radial applied, safe strength=" << safeStrength << "\n";
        gsInfo << "Patch " << pIdx << ": radial safe strength=" << safeStrength << "\n";
    }

    log << "\nApplied radial strengths:\n";
    for (index_t p = 0; p < static_cast<index_t>(newmp.nPatches()); ++p)
    {
        if (tunedPatch[p])
            log << "  patch " << p << " -> " << chosenStrength[p] << "\n";
    }

    real_t finalMaxC0Gap = 0;
    index_t finalViolatingInterfaces = 0;
    const bool finalC0Broken = hasC0Violation(newmp, kC0Tolerance,
                                              finalMaxC0Gap,
                                              finalViolatingInterfaces);
    real_t finalMaxBoundaryDrift = 0;
    index_t finalMovedBoundaryCount = 0;
    const bool boundaryPointsMoved = hasBoundaryPointDrift(*mp, newmp,
                                                           kC0Tolerance,
                                                           finalMaxBoundaryDrift,
                                                           finalMovedBoundaryCount);
    log << "\nFinal C0 check: max gap=" << finalMaxC0Gap
        << ", violating interfaces=" << finalViolatingInterfaces
        << ", tolerance=" << kC0Tolerance << "\n";
    log << "Boundary-point drift check: max drift=" << finalMaxBoundaryDrift
        << ", moved boundary points=" << finalMovedBoundaryCount
        << ", tolerance=" << kC0Tolerance << "\n";
    if (boundaryPointsMoved)
    {
        log << "Result rejected: boundary control points were modified.\n";
        gsInfo << "Result rejected: boundary control points were modified"
               << " (max drift=" << finalMaxBoundaryDrift << ").\n";
        return EXIT_FAILURE;
    }
    if (finalC0Broken)
    {
        log << "Result rejected: C0 continuity violated. Use smaller perturbation values.\n";
        gsInfo << "Result rejected: C0 continuity violated (max gap=" << finalMaxC0Gap
               << "). Use smaller perturbation values.\n";
        return EXIT_FAILURE;
    }

    const bool hasNegative =
        runMirrorAwareJacobianCheck(*mp, newmp, opt.fitGrid, opt.detThreshold,
                                    "final geometry", log, true);

    gsFileData<> fd;
    fd << newmp;
    const std::string outputDir = "C:/Users/heydatey/source/repos/gismo/filedata/generatedMPs/";
    gsFileManager::mkdir(outputDir);
    std::string outputPath = opt.outputFile;
    if (outputPath.find('/') == std::string::npos && outputPath.find('\\') == std::string::npos)
        outputPath = outputDir + outputPath;
    fd.dump(outputPath);

    std::string paraviewBase = outputPath;
    const size_t paraviewDot = paraviewBase.find_last_of('.');
    if (paraviewDot != std::string::npos)
        paraviewBase = paraviewBase.substr(0, paraviewDot);
    gsWriteParaview(newmp, paraviewBase, 1000, true, false);
    gsFileManager::open(paraviewBase + ".pvd");

    log << "\nOutput saved to: " << outputPath << "\n";
    log.close();

    gsInfo << "Distortion completed. Output: " << outputPath << "\n";
    return EXIT_SUCCESS;
}
