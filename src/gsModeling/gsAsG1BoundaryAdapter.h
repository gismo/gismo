/** @file gsAsG1BoundaryAdapter.h
    @brief Arbitrary Boundary Curve Adapter for Analysis-Suitable G1 (AS-G1) Multi-Patch Geometry Generator in G+Smo.

    Implements automatic boundary partitioning, Chebyshev kernel center finding,
    corner speed matching, and degree/knot matching for arbitrary closed planar curves.

    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gismo.h>
#include <gsModeling/gsAsG1GeometryGenerator.h>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>

namespace gismo {

template <class T = real_t>
class gsAsG1BoundaryAdapter {
public:
    typedef std::function<gsVector<T>(T)> CurveFunction;
    typedef std::function<gsVector<T>(T)> CurveDerivFunction;

    gsAsG1BoundaryAdapter()
        : m_valence(5),
          m_degree(3),
          m_numRefine(0),
          m_tMin(T(0)),
          m_tMax(T(2.0 * M_PI)),
          m_hasUserCenter(false),
          m_inwardEps(T(0.01))
    {
        m_x0.resize(2);
        m_x0.setZero();
    }

    // ====================================================================
    // Boundary Curve Setters
    // ====================================================================

    bool setBoundaryCurve(const gsCurve<T>& curve)
    {
        m_curvePtr = curve.clone();
        gsMatrix<T> pr = m_curvePtr->parameterRange();
        m_tMin = pr(0, 0);
        m_tMax = pr(0, 1);

        m_curveFunc = [this](T t) -> gsVector<T> {
            gsMatrix<T> pt(1, 1);
            pt(0, 0) = t;
            return m_curvePtr->eval(pt).col(0);
        };

        m_curveDeriv = [this](T t) -> gsVector<T> {
            gsMatrix<T> pt(1, 1);
            pt(0, 0) = t;
            return m_curvePtr->deriv(pt).col(0);
        };

        return true;
    }

    bool setBoundaryCurve(CurveFunction func, CurveDerivFunction deriv, T tMin = T(0), T tMax = T(2.0 * M_PI))
    {
        m_curveFunc = func;
        m_curveDeriv = deriv;
        m_tMin = tMin;
        m_tMax = tMax;
        m_curvePtr.reset();
        return true;
    }

    bool setBoundarySegments(const std::vector<gsBSpline<T>>& segments)
    {
        if (segments.empty()) return false;
        m_segments = segments;
        return true;
    }

    // ====================================================================
    // Configuration Setters & Getters
    // ====================================================================

    void setValence(int V)      { m_valence = V; }
    int valence() const         { return m_valence; }

    void setDegree(int p)       { m_degree = p; }
    int degree() const          { return m_degree; }

    void setNumRefine(int r)    { m_numRefine = r; }
    int numRefine() const       { return m_numRefine; }

    void setCenterVertex(const gsVector<T>& x0) {
        m_x0 = x0;
        m_hasUserCenter = true;
    }
    const gsVector<T>& centerVertex() const { return m_x0; }

    void setInwardCornerOffset(T eps) { m_inwardEps = eps; }
    T inwardCornerOffset() const { return m_inwardEps; }

    // ====================================================================
    // Step 1: Chebyshev Kernel Center Finder
    // ====================================================================

    gsVector<T> computeOptimalCenter(int numSamples = 300) const
    {
        if (numSamples < 20) numSamples = 300;

        std::vector<gsVector<T>> pts(numSamples);
        std::vector<gsVector<T>> dpts(numSamples);
        T dt = (m_tMax - m_tMin) / numSamples;
        for (int k = 0; k < numSamples; ++k) {
            T t = m_tMin + k * dt;
            pts[k] = evalCurve(t);
            dpts[k] = evalDeriv(t);
        }

        // Compute signed polygon area
        T signedArea = 0;
        for (int k = 0; k < numSamples; ++k) {
            int k1 = (k + 1) % numSamples;
            signedArea += (pts[k](0) * pts[k1](1) - pts[k1](0) * pts[k](1));
        }
        signedArea *= T(0.5);
        bool isCCW = (signedArea > 0);

        // Centroid of the polygon
        gsVector<T> C(2); C.setZero();
        if (std::abs(signedArea) > T(1e-12)) {
            for (int k = 0; k < numSamples; ++k) {
                int k1 = (k + 1) % numSamples;
                T cross = (pts[k](0) * pts[k1](1) - pts[k1](0) * pts[k](1));
                C(0) += (pts[k](0) + pts[k1](0)) * cross;
                C(1) += (pts[k](1) + pts[k1](1)) * cross;
            }
            C /= (T(6.0) * signedArea);
        } else {
            for (int k = 0; k < numSamples; ++k) C += pts[k];
            C /= static_cast<T>(numSamples);
        }

        // Ray monotonicity check: det(pt - x, dpt)
        auto rayMonotonicity = [&](const gsVector<T>& x) -> T {
            T minCross = 1e12;
            for (int k = 0; k < numSamples; ++k) {
                gsVector<T> r = pts[k] - x;
                T cross = r(0) * dpts[k](1) - r(1) * dpts[k](0);
                if (!isCCW) cross = -cross;
                if (cross < minCross) minCross = cross;
            }
            return minCross;
        };

        // Minimum distance to boundary points
        auto clearance = [&](const gsVector<T>& x) -> T {
            T minD = 1e12;
            for (int k = 0; k < numSamples; ++k) {
                T d = (x - pts[k]).norm();
                if (d < minD) minD = d;
            }
            return minD;
        };

        // Bounding box
        T minX = pts[0](0), maxX = pts[0](0);
        T minY = pts[0](1), maxY = pts[0](1);
        for (int k = 1; k < numSamples; ++k) {
            minX = std::min(minX, pts[k](0)); maxX = std::max(maxX, pts[k](0));
            minY = std::min(minY, pts[k](1)); maxY = std::max(maxY, pts[k](1));
        }

        gsVector<T> bestCenter = C;
        T bestScore = -1e12;

        const int nGrid = 30;
        T dx = (maxX - minX) / (nGrid + 1);
        T dy = (maxY - minY) / (nGrid + 1);

        for (int ix = 1; ix <= nGrid; ++ix) {
            for (int iy = 1; iy <= nGrid; ++iy) {
                gsVector<T> cand(2);
                cand << minX + ix * dx, minY + iy * dy;
                T mono = rayMonotonicity(cand);
                if (mono > T(1e-4)) {
                    T score = clearance(cand);
                    if (score > bestScore) {
                        bestScore = score;
                        bestCenter = cand;
                    }
                } else if (bestScore < 0) {
                    if (mono > bestScore) {
                        bestScore = mono;
                        bestCenter = cand;
                    }
                }
            }
        }

        return bestCenter;
    }

    // ====================================================================
    // Step 2: Feature Detection
    // ====================================================================

    std::vector<T> detectFeatureCorners(T angleThresholdRad = T(0.3)) const
    {
        std::vector<T> corners;
        const int nSamples = 300;
        T dt = (m_tMax - m_tMin) / nSamples;

        std::vector<gsVector<T>> tangents(nSamples);
        for (int k = 0; k < nSamples; ++k) {
            T t = m_tMin + k * dt;
            tangents[k] = evalDeriv(t);
            T len = tangents[k].norm();
            if (len > T(1e-12)) tangents[k] /= len;
        }

        for (int k = 0; k < nSamples; ++k) {
            int k1 = (k + 1) % nSamples;
            T cosAng = tangents[k].dot(tangents[k1]);
            T sinAng = std::abs(det2D(tangents[k], tangents[k1]));
            if (sinAng > std::sin(angleThresholdRad) || cosAng < std::cos(angleThresholdRad)) {
                T t_corner = m_tMin + (k + 0.5) * dt;
                corners.push_back(t_corner);
            }
        }

        return corners;
    }

    // ====================================================================
    // Step 3 & 4: Boundary Partitioning & Corner Speed Matching
    // ====================================================================

    bool partitionAndReparameterize(
        int valence,
        std::vector<gsBSpline<T>>& outEastCurves,
        std::vector<gsBSpline<T>>& outNorthCurves,
        std::vector<gsVector<T>>& outSpokes) const
    {
        if (valence <= 0) valence = m_valence;
        const int V = valence;

        gsVector<T> x0 = m_hasUserCenter ? m_x0 : computeOptimalCenter();

        // Sample curve on fine grid to compute unwrapped polar angles
        const int N_GRID = 2000;
        std::vector<T> t_grid(N_GRID + 1);
        std::vector<T> th_cum(N_GRID + 1);
        T dt = (m_tMax - m_tMin) / N_GRID;

        T last_th = 0;
        for (int k = 0; k <= N_GRID; ++k) {
            T t = m_tMin + k * dt;
            t_grid[k] = t;
            gsVector<T> p = evalCurve(t) - x0;
            T th = std::atan2(p(1), p(0));
            if (k == 0) {
                last_th = th;
                th_cum[0] = th;
            } else {
                T diff = th - last_th;
                while (diff >  M_PI) diff -= 2.0 * M_PI;
                while (diff < -M_PI) diff += 2.0 * M_PI;
                last_th += diff;
                th_cum[k] = last_th;
            }
        }

        T totalTheta = th_cum.back() - th_cum.front();
        bool isPositiveAngle = (totalTheta > 0);
        T baseTheta = th_cum.front();

        auto findTForTheta = [&](T targetTheta) -> T {
            while (targetTheta < th_cum.front()) targetTheta += 2.0 * M_PI;
            while (targetTheta > th_cum.back())  targetTheta -= 2.0 * M_PI;
            int low = 0;
            for (int k = 0; k < N_GRID; ++k) {
                if ((th_cum[k] <= targetTheta && targetTheta <= th_cum[k+1]) ||
                    (th_cum[k] >= targetTheta && targetTheta >= th_cum[k+1])) {
                    low = k;
                    break;
                }
            }
            T denom = th_cum[low+1] - th_cum[low];
            T frac = (std::abs(denom) > T(1e-14)) ? (targetTheta - th_cum[low]) / denom : T(0);
            return t_grid[low] + frac * dt;
        };

        // Numerical arc length helper between two parameter values along the curve
        auto arcDist = [&](T tA, T tB) -> T {
            while (tB < tA) tB += (m_tMax - m_tMin);
            T sum = 0;
            const int ns = 30;
            T h = (tB - tA) / ns;
            for (int j = 0; j < ns; ++j) {
                T t = tA + (j + 0.5) * h;
                sum += evalDeriv(t).norm() * h;
            }
            return sum;
        };

        std::vector<T> t_spokes(V);
        std::vector<T> t_mids(V);
        outSpokes.resize(V);
        std::vector<gsVector<T>> mids(V);

        for (int i = 0; i < V; ++i) {
            T angleSign = isPositiveAngle ? T(1) : T(-1);
            T target_spoke = baseTheta + angleSign * (2.0 * M_PI * i / V);
            T target_mid   = baseTheta + angleSign * (2.0 * M_PI * (i + 0.5) / V);

            t_spokes[i] = findTForTheta(target_spoke);
            t_mids[i]   = findTForTheta(target_mid);

            outSpokes[i] = evalCurve(t_spokes[i]);
            mids[i]      = evalCurve(t_mids[i]);
        }

        // 1D boundary B-spline Basis matching the volume patches
        gsKnotVector<T> kv1D(0.0, 1.0, 0, m_degree + 1);
        for (int r = 0; r < m_numRefine; ++r) kv1D.uniformRefine();
        gsBSplineBasis<T> basis(kv1D);
        const int n = basis.size();

        T u1 = basis.knots()[basis.degree() + 1];
        T factor = basis.degree() / u1;

        // Quadrature for least-squares fitting of intermediate control points
        gsGaussRule<T> rule(basis.degree() + 3);
        gsMatrix<T> nodes; gsVector<T> weights;
        rule.mapToAll(kv1D.breaks(), nodes, weights);
        const int nq = nodes.cols();
        gsMatrix<index_t> act = basis.active(nodes);
        gsMatrix<T> basisVals;
        basis.eval_into(nodes, basisVals);

        std::vector<gsMatrix<T>> customEast(V);
        std::vector<gsMatrix<T>> customNorth(V);

        T eps = m_inwardEps;

        // Speed matching at spoke corners
        for (int i = 0; i < V; ++i) {
            int im1 = (i - 1 + V) % V;
            int ip1 = (i + 1) % V;

            // Determinants of virtual skeleton
            gsVector<T> v_im1 = outSpokes[im1] - x0;
            gsVector<T> v_i   = outSpokes[i] - x0;
            gsVector<T> v_ip1 = outSpokes[ip1] - x0;
            T D1 = det2D(v_im1, v_i);
            T D2 = det2D(v_i, v_ip1);
            T R = D1 / D2;

            T L_east  = arcDist(t_spokes[i], t_mids[i]);
            T L_north = arcDist(t_mids[im1], t_spokes[i]);
            T r_bar = L_north / L_east;
            T kappa = std::sqrt(R / r_bar);

            T s_E = L_east / kappa;
            T s_N = L_north * kappa;

            gsVector<T> dGamma_spoke = evalDeriv(t_spokes[i]);
            T spk_len = dGamma_spoke.norm();
            gsVector<T> t_hat = (spk_len > T(1e-12)) ? (dGamma_spoke / spk_len) : dGamma_spoke;

            gsVector<T> t1_star = -s_N * t_hat;
            gsVector<T> t2_star =  s_E * t_hat;

            // East curve of patch i
            customEast[i].resize(n, 2);
            customEast[i].row(0) = outSpokes[i].transpose();
            customEast[i].row(1) = (outSpokes[i] + t2_star / factor).transpose();
            customEast[i].row(n - 1) = mids[i].transpose();

            gsVector<T> dGamma_mid = evalDeriv(t_mids[i]);
            T mid_len = dGamma_mid.norm();
            gsVector<T> t_hat_mid = (mid_len > T(1e-12)) ? (dGamma_mid / mid_len) : dGamma_mid;
            gsVector<T> norm_mid(2); norm_mid << -t_hat_mid(1), t_hat_mid(0);
            if (norm_mid.dot(x0 - mids[i]) < 0) norm_mid = -norm_mid;

            gsVector<T> t_E_mid = L_east * (t_hat_mid - eps * norm_mid);
            customEast[i].row(n - 2) = (mids[i] - t_E_mid / factor).transpose();

            // North curve of patch im1
            customNorth[im1].resize(n, 2);
            customNorth[im1].row(0) = outSpokes[i].transpose();
            customNorth[im1].row(1) = (outSpokes[i] + t1_star / factor).transpose();
            customNorth[im1].row(n - 1) = mids[im1].transpose();

            gsVector<T> dGamma_mim1 = evalDeriv(t_mids[im1]);
            T mim1_len = dGamma_mim1.norm();
            gsVector<T> t_hat_mim1 = (mim1_len > T(1e-12)) ? (dGamma_mim1 / mim1_len) : dGamma_mim1;
            gsVector<T> norm_mim1(2); norm_mim1 << -t_hat_mim1(1), t_hat_mim1(0);
            if (norm_mim1.dot(x0 - mids[im1]) < 0) norm_mim1 = -norm_mim1;

            gsVector<T> t_N_mid = L_north * (-t_hat_mim1 - eps * norm_mim1);
            customNorth[im1].row(n - 2) = (mids[im1] - t_N_mid / factor).transpose();
        }

        // Fit intermediate control points if n > 4
        const int nFree = n - 4;
        if (nFree > 0) {
            for (int i = 0; i < V; ++i) {
                int im1 = (i - 1 + V) % V;
                T dt_east = t_mids[i] - t_spokes[i];
                while (dt_east < 0) dt_east += (m_tMax - m_tMin);
                while (dt_east > (m_tMax - m_tMin)) dt_east -= (m_tMax - m_tMin);

                T dt_north = t_spokes[i] - t_mids[im1];
                while (dt_north < 0) dt_north += (m_tMax - m_tMin);
                while (dt_north > (m_tMax - m_tMin)) dt_north -= (m_tMax - m_tMin);

                // East interior fit
                gsMatrix<T> Mat_E(nFree, nFree); Mat_E.setZero();
                gsMatrix<T> F_E(nFree, 2);       F_E.setZero();
                for (int q = 0; q < nq; ++q) {
                    T w = weights(q);
                    T v = nodes(0, q);
                    T t_val = t_spokes[i] + v * dt_east;
                    gsVector<T> c_t = evalCurve(t_val);

                    T phi_0 = 0, phi_1 = 0, phi_n2 = 0, phi_n1 = 0;
                    for (int k = 0; k < act.rows(); ++k) {
                        if (act(k, q) == 0)     phi_0  = basisVals(k, q);
                        if (act(k, q) == 1)     phi_1  = basisVals(k, q);
                        if (act(k, q) == n - 2) phi_n2 = basisVals(k, q);
                        if (act(k, q) == n - 1) phi_n1 = basisVals(k, q);
                    }
                    gsVector<T> f_eff = c_t - phi_0  * customEast[i].row(0).transpose()
                                            - phi_1  * customEast[i].row(1).transpose()
                                            - phi_n2 * customEast[i].row(n - 2).transpose()
                                            - phi_n1 * customEast[i].row(n - 1).transpose();

                    for (int k1 = 0; k1 < act.rows(); ++k1) {
                        int g1 = act(k1, q);
                        if (g1 < 2 || g1 > n - 3) continue;
                        int r1 = g1 - 2;
                        F_E.row(r1) += w * basisVals(k1, q) * f_eff.transpose();
                        for (int k2 = 0; k2 < act.rows(); ++k2) {
                            int g2 = act(k2, q);
                            if (g2 < 2 || g2 > n - 3) continue;
                            int r2 = g2 - 2;
                            Mat_E(r1, r2) += w * basisVals(k1, q) * basisVals(k2, q);
                        }
                    }
                }
                gsMatrix<T> sol_E = Mat_E.colPivHouseholderQr().solve(F_E);
                for (int r = 0; r < nFree; ++r) customEast[i].row(r + 2) = sol_E.row(r);

                // North interior fit
                gsMatrix<T> Mat_N(nFree, nFree); Mat_N.setZero();
                gsMatrix<T> F_N(nFree, 2);       F_N.setZero();
                for (int q = 0; q < nq; ++q) {
                    T w = weights(q);
                    T u = nodes(0, q);
                    T t_val = t_spokes[i] - u * dt_north;
                    gsVector<T> c_t = evalCurve(t_val);

                    T phi_0 = 0, phi_1 = 0, phi_n2 = 0, phi_n1 = 0;
                    for (int k = 0; k < act.rows(); ++k) {
                        if (act(k, q) == 0)     phi_0  = basisVals(k, q);
                        if (act(k, q) == 1)     phi_1  = basisVals(k, q);
                        if (act(k, q) == n - 2) phi_n2 = basisVals(k, q);
                        if (act(k, q) == n - 1) phi_n1 = basisVals(k, q);
                    }
                    gsVector<T> f_eff = c_t - phi_0  * customNorth[im1].row(0).transpose()
                                            - phi_1  * customNorth[im1].row(1).transpose()
                                            - phi_n2 * customNorth[im1].row(n - 2).transpose()
                                            - phi_n1 * customNorth[im1].row(n - 1).transpose();

                    for (int k1 = 0; k1 < act.rows(); ++k1) {
                        int g1 = act(k1, q);
                        if (g1 < 2 || g1 > n - 3) continue;
                        int r1 = g1 - 2;
                        F_N.row(r1) += w * basisVals(k1, q) * f_eff.transpose();
                        for (int k2 = 0; k2 < act.rows(); ++k2) {
                            int g2 = act(k2, q);
                            if (g2 < 2 || g2 > n - 3) continue;
                            int r2 = g2 - 2;
                            Mat_N(r1, r2) += w * basisVals(k1, q) * basisVals(k2, q);
                        }
                    }
                }
                gsMatrix<T> sol_N = Mat_N.colPivHouseholderQr().solve(F_N);
                for (int r = 0; r < nFree; ++r) customNorth[im1].row(r + 2) = sol_N.row(r);
            }
        }

        // Convert to gsBSpline curves
        outEastCurves.clear();
        outNorthCurves.clear();
        outEastCurves.reserve(V);
        outNorthCurves.reserve(V);

        for (int i = 0; i < V; ++i) {
            outEastCurves.push_back(gsBSpline<T>(basis, customEast[i]));
            outNorthCurves.push_back(gsBSpline<T>(basis, customNorth[i]));
        }

        m_cachedEastCP = customEast;
        m_cachedNorthCP = customNorth;
        m_cachedSpokes = outSpokes;
        m_cachedX0 = x0;

        return true;
    }

    // ====================================================================
    // Step 5: Configure Generator
    // ====================================================================

    bool configureGenerator(gsAsG1GeometryGenerator<T>& generator) const
    {
        std::vector<gsBSpline<T>> eastCurves, northCurves;
        std::vector<gsVector<T>> spokes;
        if (!partitionAndReparameterize(m_valence, eastCurves, northCurves, spokes))
            return false;

        generator.setValence(m_valence);
        generator.setDegree(m_degree);
        generator.setNumRefine(m_numRefine);
        generator.setCenterVertex(m_cachedX0);
        generator.setBoundaryCorners(m_cachedSpokes);
        generator.setUseVirtualParallelogram(true);
        generator.setCustomBoundaryControlPoints(m_cachedEastCP, m_cachedNorthCP);

        return true;
    }

    // ====================================================================
    // Verification & Diagnostics
    // ====================================================================

    T computeMaxBoundaryError(const gsMultiPatch<T>& mp, int numEvalPerEdge = 100) const
    {
        const int V = mp.nPatches();
        T maxDist = 0;

        for (int i = 0; i < V; ++i) {
            for (int s = 0; s <= numEvalPerEdge; ++s) {
                T param = static_cast<T>(s) / numEvalPerEdge;
                gsMatrix<T> ptE(2, 1); ptE << 1.0, param;
                gsVector<T> pE = mp.patch(i).eval(ptE).col(0);
                T distE = distanceToBoundary(pE);
                if (distE > maxDist) maxDist = distE;

                gsMatrix<T> ptN(2, 1); ptN << param, 1.0;
                gsVector<T> pN = mp.patch(i).eval(ptN).col(0);
                T distN = distanceToBoundary(pN);
                if (distN > maxDist) maxDist = distN;
            }
        }
        return maxDist;
    }

    T distanceToBoundary(const gsVector<T>& pt, int searchSamples = 300) const
    {
        T minDist = 1e12;
        T dt = (m_tMax - m_tMin) / searchSamples;
        int bestK = 0;
        for (int k = 0; k < searchSamples; ++k) {
            T t = m_tMin + k * dt;
            T d = (evalCurve(t) - pt).norm();
            if (d < minDist) {
                minDist = d;
                bestK = k;
            }
        }

        // Local Newton refinement with periodic wrap-around
        T t_est = m_tMin + bestK * dt;
        for (int it = 0; it < 10; ++it) {
            gsVector<T> p = evalCurve(t_est);
            gsVector<T> d1 = evalDeriv(t_est);
            T f = (p - pt).dot(d1);
            T df = d1.squaredNorm();
            if (df < T(1e-12)) break;
            T step = f / df;
            t_est -= step;
            if (std::abs(step) < T(1e-10)) break;
        }
        return (evalCurve(t_est) - pt).norm();
    }

    // ====================================================================
    // Standard Test Shape Creators
    // ====================================================================

    template <typename Func>
    static gsBSpline<T> interpolateFunction(Func func, T tMin, T tMax, int numPoints = 64)
    {
        gsKnotVector<T> kv(tMin, tMax, numPoints - 3, 4);
        gsBSplineBasis<T> basis(kv);
        const int n = basis.size();

        gsMatrix<T> greville = basis.anchors();
        gsMatrix<T> B_act;
        basis.eval_into(greville, B_act);
        gsMatrix<index_t> act;
        basis.active_into(greville, act);

        gsMatrix<T> B(n, n);
        B.setZero();
        for (int j = 0; j < n; ++j) {
            for (int r = 0; r < act.rows(); ++r) {
                B(j, act(r, j)) = B_act(r, j);
            }
        }

        gsMatrix<T> targets(n, 2);
        for (int i = 0; i < n; ++i) {
            T t = greville(0, i);
            targets.row(i) = func(t).transpose();
        }
        gsMatrix<T> coefs = B.colPivHouseholderQr().solve(targets);
        return gsBSpline<T>(basis, coefs);
    }

    static gsBSpline<T> createEllipse(T a = T(1.5), T b = T(0.8), int numPoints = 64)
    {
        return interpolateFunction(
            [=](T t) {
                gsVector<T> p(2);
                p << a * std::cos(t), b * std::sin(t);
                return p;
            },
            T(0), T(2.0 * M_PI), numPoints);
    }

    static gsBSpline<T> createTeardrop(T scale = T(1.0), int numPoints = 64)
    {
        return interpolateFunction(
            [=](T t) {
                gsVector<T> p(2);
                T s = std::sin(t * 0.5);
                p << scale * 1.5 * std::cos(t), scale * 1.0 * std::sin(t) * (s * s + 0.3);
                return p;
            },
            T(0), T(2.0 * M_PI), numPoints);
    }

    static gsBSpline<T> createBean(T scale = T(1.0), int numPoints = 64)
    {
        return interpolateFunction(
            [=](T t) {
                gsVector<T> p(2);
                T r = scale * (T(1.0) + T(0.3) * std::cos(t) - T(0.2) * std::cos(2.0 * t) + T(0.1) * std::sin(3.0 * t));
                p << r * std::cos(t), r * std::sin(t);
                return p;
            },
            T(0), T(2.0 * M_PI), numPoints);
    }

    static gsBSpline<T> createArbitraryClosedBSpline(int numPoints = 64)
    {
        std::vector<gsVector<T>> pts(8);
        for (int i = 0; i < 8; ++i) pts[i].resize(2);
        pts[0] << 1.4,  0.2;
        pts[1] << 1.1,  1.0;
        pts[2] << -0.2,  1.3;
        pts[3] << -1.2,  0.7;
        pts[4] << -1.5, -0.3;
        pts[5] << -0.8, -1.1;
        pts[6] << 0.3, -1.3;
        pts[7] << 1.2, -0.7;

        // Centripetal Catmull-Rom spline through non-symmetric points
        auto catmullRom = [pts](T t) -> gsVector<T> {
            while (t < T(0)) t += T(2.0 * M_PI);
            while (t >= T(2.0 * M_PI)) t -= T(2.0 * M_PI);
            T u = (t / T(2.0 * M_PI)) * T(8.0);
            int i1 = static_cast<int>(u);
            if (i1 >= 8) i1 = 7;
            T s = u - i1;
            int i0 = (i1 - 1 + 8) % 8;
            int i2 = (i1 + 1) % 8;
            int i3 = (i1 + 2) % 8;

            T s2 = s * s;
            T s3 = s2 * s;
            T h00 =  T(2.0)*s3 - T(3.0)*s2 + T(1.0);
            T h10 =         s3 - T(2.0)*s2 + s;
            T h01 = -T(2.0)*s3 + T(3.0)*s2;
            T h11 =         s3 -        s2;

            gsVector<T> m1 = T(0.5) * (pts[i2] - pts[i0]);
            gsVector<T> m2 = T(0.5) * (pts[i3] - pts[i1]);

            return h00 * pts[i1] + h10 * m1 + h01 * pts[i2] + h11 * m2;
        };

        return interpolateFunction(catmullRom, T(0), T(2.0 * M_PI), numPoints);
    }

private:
    int m_valence;
    int m_degree;
    int m_numRefine;
    T m_tMin;
    T m_tMax;
    bool m_hasUserCenter;
    T m_inwardEps;
    gsVector<T> m_x0;

    typename gsCurve<T>::uPtr m_curvePtr;
    CurveFunction m_curveFunc;
    CurveDerivFunction m_curveDeriv;
    std::vector<gsBSpline<T>> m_segments;

    mutable std::vector<gsMatrix<T>> m_cachedEastCP;
    mutable std::vector<gsMatrix<T>> m_cachedNorthCP;
    mutable std::vector<gsVector<T>> m_cachedSpokes;
    mutable gsVector<T> m_cachedX0;

    inline gsVector<T> evalCurve(T t) const
    {
        T period = m_tMax - m_tMin;
        if (period > T(1e-12)) {
            while (t < m_tMin) t += period;
            while (t > m_tMax) t -= period;
        }
        if (m_curveFunc) return m_curveFunc(t);
        gsVector<T> p(2); p.setZero();
        return p;
    }

    inline gsVector<T> evalDeriv(T t) const
    {
        T period = m_tMax - m_tMin;
        if (period > T(1e-12)) {
            while (t < m_tMin) t += period;
            while (t > m_tMax) t -= period;
        }
        if (m_curveDeriv) return m_curveDeriv(t);
        T h = T(1e-6);
        return (evalCurve(t + h) - evalCurve(t - h)) / (T(2.0) * h);
    }
};

template <class T>
bool gsAsG1GeometryGenerator<T>::setBoundaryFromCurve(const gsCurve<T>& boundaryCurve, int valence)
{
    gsAsG1BoundaryAdapter<T> adapter;
    if (!adapter.setBoundaryCurve(boundaryCurve)) return false;
    if (valence > 0) adapter.setValence(valence);
    else adapter.setValence(m_valence);
    adapter.setDegree(m_degree);
    adapter.setNumRefine(m_numRefine);
    return adapter.configureGenerator(*this);
}

} // namespace gismo
