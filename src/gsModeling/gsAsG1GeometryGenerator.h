
/** @file gsAsG1GeometryGenerator.h
    @brief Analysis-Suitable G1 (AS-G1) Multi-Patch Geometry Generator in C++ using G+Smo.

    Implements the "Bilinear Skeleton + Linearized Constrained Optimization" pipeline
    following Kapl, Sangalli, Takacs (2018) and Hasanova, Takacs, Takacs (2026).

    This file is part of the G+Smo library.
    Author(s): F. Hasanova, Antigravity
*/

#pragma once

#include <gismo.h>
#include <gsMatrix/gsSparseSolver.h>
#include <vector>
#include <cmath>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>

namespace gismo {

// ====================================================================
// 2D Determinant Helper
// ====================================================================
template <class T>
inline T det2D(const gsVector<T>& a, const gsVector<T>& b) {
    return a(0) * b(1) - a(1) * b(0);
}

// ====================================================================
// Gluing Data Representation
// ====================================================================
template <class T = real_t>
struct AsG1GluingData {
    T D1_0 = 0, D1_1 = 0;
    T D2_0 = 0, D2_1 = 0;
    T b1_0 = 0, b1_1 = 0;
    T b2_0 = 0, b2_1 = 0;

    T alpha1(T xi) const { return D1_0 * (T(1) - xi) + D1_1 * xi; }
    T alpha2(T xi) const { return D2_0 * (T(1) - xi) + D2_1 * xi; }
    T beta1(T xi)  const { return b1_0 * (T(1) - xi) + b1_1 * xi; }
    T beta2(T xi)  const { return b2_0 * (T(1) - xi) + b2_1 * xi; }
    T b_glue(T xi) const { return alpha1(xi) * beta2(xi) + alpha2(xi) * beta1(xi); }
};

// ====================================================================
// 1D L2 Boundary Curve Fitting with Fixed Endpoints
// ====================================================================
template <class T = real_t>
gsMatrix<T> fitBoundaryCurve(
    const gsBSplineBasis<T>& bndBasis,
    std::function<gsVector<T>(T)> curve)
{
    const int n = bndBasis.size();
    gsMatrix<T> P(n, 2);
    P.row(0)     = curve(T(0)).transpose();
    P.row(n - 1) = curve(T(1)).transpose();
    if (n <= 2) return P;

    const int deg = bndBasis.degree();
    gsGaussRule<T> rule(deg + 2);
    gsMatrix<T> nodes; gsVector<T> weights;
    rule.mapToAll(bndBasis.knots().breaks(), nodes, weights);

    const int nQuad = nodes.cols();
    gsMatrix<index_t> act = bndBasis.active(nodes);
    gsMatrix<T> basisVals;
    bndBasis.eval_into(nodes, basisVals);

    const int nInt = n - 2;
    gsMatrix<T> M(nInt, nInt); M.setZero();
    gsMatrix<T> F(nInt, 2);    F.setZero();

    for (int q = 0; q < nQuad; ++q) {
        T w = weights(q);
        T t = nodes(0, q);
        gsVector<T> c_t = curve(t);

        T phi_0 = 0, phi_last = 0;
        for (int k = 0; k < act.rows(); ++k) {
            if (act(k, q) == 0)     phi_0    = basisVals(k, q);
            if (act(k, q) == n - 1) phi_last = basisVals(k, q);
        }

        gsVector<T> f_eff = c_t - phi_0 * P.row(0).transpose()
                                - phi_last * P.row(n - 1).transpose();

        for (int k1 = 0; k1 < act.rows(); ++k1) {
            int g1 = act(k1, q);
            if (g1 < 1 || g1 > nInt) continue;
            T phi_i = basisVals(k1, q);
            F.row(g1 - 1) += w * phi_i * f_eff.transpose();

            for (int k2 = 0; k2 < act.rows(); ++k2) {
                int g2 = act(k2, q);
                if (g2 < 1 || g2 > nInt) continue;
                T phi_j = basisVals(k2, q);
                M(g1 - 1, g2 - 1) += w * phi_i * phi_j;
            }
        }
    }

    gsMatrix<T> P_int = M.colPivHouseholderQr().solve(F);
    for (int i = 0; i < nInt; ++i)
        P.row(i + 1) = P_int.row(i);

    return P;
}

// ====================================================================
// AS-G1 Geometry Generator Class
// ====================================================================
template <class T = real_t>
class gsAsG1GeometryGenerator {
public:
    typedef std::function<gsVector<T>(T)> CurveFunction;

    gsAsG1GeometryGenerator(int valence = 5, int degree = 3, int numRefine = 0)
        : m_valence(valence), m_degree(degree), m_numRefine(numRefine),
          m_radius(T(1)), m_midRadius(T(1.2)),
          m_laplaceWeight(T(1.0)), m_biharmWeight(T(1e-4)),
          m_x0(2)
    {
        m_x0.setZero();
    }

    void setValence(int V)      { m_valence = V; }
    void setDegree(int p)       { m_degree = p; }
    void setNumRefine(int r)    { m_numRefine = r; }
    void setRadius(T r)         { m_radius = r; }
    void setMidRadius(T rmid)   { m_midRadius = rmid; }
    void setLaplaceWeight(T w)  { m_laplaceWeight = w; }
    void setBiharmWeight(T w)   { m_biharmWeight = w; }
    void setCenterVertex(const gsVector<T>& x0) { m_x0 = x0; }

    void setUseVirtualParallelogram(bool enable = true) {
        m_useVirtualParallelogram = enable;
    }
    bool useVirtualParallelogram() const { return m_useVirtualParallelogram; }

    void setRandomize(bool enable, T anglePerturb = T(0.15), T radiusPerturb = T(0.2), unsigned int seed = 42) {
        m_randomize = enable;
        m_anglePerturb = anglePerturb;
        m_radiusPerturb = radiusPerturb;
        m_seed = seed;
    }
    bool isRandomized() const { return m_randomize; }

    void setBoundaryCorners(const std::vector<gsVector<T>>& corners) {
        m_corners = corners;
    }
    void setMidCorners(const std::vector<gsVector<T>>& midCorners) {
        m_midCorners = midCorners;
    }

    // Pipeline Steps
    void setupBilinearSkeleton();
    void extractGluingData();
    void refineGeometry();
    void assembleStiffnessMatrix();
    void assembleConstraints(
        const std::vector<CurveFunction>& eastCurves = {},
        const std::vector<CurveFunction>& northCurves = {});
    bool solveKKT();

    // Full Pipeline Execution
    bool generate(
        const std::vector<CurveFunction>& eastCurves = {},
        const std::vector<CurveFunction>& northCurves = {})
    {
        setupBilinearSkeleton();
        extractGluingData();
        refineGeometry();
        assembleStiffnessMatrix();
        assembleConstraints(eastCurves, northCurves);
        return solveKKT();
    }

    // Verification and Diagnostics
    T verifyAsG1(int numSamples = 100, std::vector<T>* perInterfaceResiduals = nullptr) const;
    bool checkJacobianRegularity(int samplesPerDir = 20, T* outMinDet = nullptr, T* outMaxDet = nullptr) const;

    // Output and Serialization
    void exportParaview(const std::string& basename, int npts = 1000) const {
        gsWriteParaview(m_mp, basename, npts, true);
    }
    void exportXml(const std::string& filename) const {
        gsWrite(m_mp, filename);
    }

    // Custom Boundaries from file or external tools (Mathematica designer)
    void setCustomBoundaryControlPoints(const std::vector<gsMatrix<T>>& eastCP,
                                       const std::vector<gsMatrix<T>>& northCP) {
        m_custom_east_cp = eastCP;
        m_custom_north_cp = northCP;
    }
    bool hasCustomBoundaries() const {
        return !m_custom_east_cp.empty() && !m_custom_north_cp.empty();
    }
    bool loadBoundaryDefinition(const std::string& filename);
    bool setBoundaryFromCurve(const gsCurve<T>& boundaryCurve, int valence = 0);

    // Accessors
    const gsMultiPatch<T>& multiPatch() const { return m_mp; }
    const gsMultiPatch<T>& coarseMultiPatch() const { return m_mp_coarse; }
    const std::vector<AsG1GluingData<T>>& gluingData() const { return m_gd; }
    int valence() const { return m_valence; }
    int degree() const { return m_degree; }

private:
    int m_valence;
    int m_degree;
    int m_numRefine;
    T m_radius;
    T m_midRadius;
    T m_laplaceWeight;
    T m_biharmWeight;
    gsVector<T> m_x0;

    bool m_randomize = false;
    bool m_useVirtualParallelogram = false;
    T m_anglePerturb = T(0.15);
    T m_radiusPerturb = T(0.2);
    unsigned int m_seed = 42;

    std::vector<gsVector<T>> m_corners;
    std::vector<gsVector<T>> m_midCorners;
    std::vector<gsMatrix<T>> m_custom_east_cp;
    std::vector<gsMatrix<T>> m_custom_north_cp;

    gsMultiPatch<T> m_mp_coarse;
    gsMultiPatch<T> m_mp;
    std::vector<AsG1GluingData<T>> m_gd;

    gsBSplineBasis<T> m_basis1d;
    gsTensorBSplineBasis<2, T> m_patchBasis;

    int m_Nu = 0, m_Nv = 0, m_Npatch = 0, m_Ntotal = 0;
    gsSparseMatrix<T> m_K_global;

    std::vector<std::vector<std::pair<int, T>>> m_constr_rows;
    std::vector<T> m_rhs_x;
    std::vector<T> m_rhs_y;

    inline int globalIdx(int patch, int iu, int iv) const {
        return patch * m_Npatch + (iu + iv * m_Nu);
    }

    void addConstraintRow(const std::vector<std::pair<int, T>>& terms, T gx, T gy) {
        m_constr_rows.push_back(terms);
        m_rhs_x.push_back(gx);
        m_rhs_y.push_back(gy);
    }
};

// ====================================================================
// Implementation of Pipeline Methods
// ====================================================================

template <class T>
void gsAsG1GeometryGenerator<T>::setupBilinearSkeleton()
{
    const int V = m_valence;
    if (m_corners.size() != static_cast<size_t>(V)) {
        m_corners.resize(V);
        m_midCorners.resize(V);

        if (m_randomize) {
            std::mt19937 rng(m_seed);
            T max_delta_ang = m_anglePerturb * 2.0 * M_PI / V;
            std::uniform_real_distribution<double> dist_ang(-max_delta_ang, max_delta_ang);
            std::uniform_real_distribution<double> dist_r(1.0 - m_radiusPerturb, 1.0 + m_radiusPerturb);

            std::vector<double> angles(V);
            for (int i = 0; i < V; ++i) {
                angles[i] = 2.0 * M_PI * i / V + dist_ang(rng);
            }
            std::sort(angles.begin(), angles.end());

            for (int i = 0; i < V; ++i) {
                T r = m_radius * dist_r(rng);
                m_corners[i].resize(2);
                m_corners[i] << r * std::cos(angles[i]), r * std::sin(angles[i]);
            }
            for (int i = 0; i < V; ++i) {
                int i1 = (i + 1) % V;
                m_midCorners[i] = m_corners[i] + m_corners[i1];
            }
        } else {
            for (int i = 0; i < V; ++i) {
                T th = 2.0 * M_PI * i / V;
                m_corners[i].resize(2);
                m_corners[i] << m_radius * std::cos(th), m_radius * std::sin(th);
            }
            for (int i = 0; i < V; ++i) {
                T th_mid = 2.0 * M_PI * (i + 0.5) / V;
                m_midCorners[i].resize(2);
                m_midCorners[i] << m_midRadius * std::cos(th_mid), m_midRadius * std::sin(th_mid);
            }
        }
    } else if (m_useVirtualParallelogram) {
        m_midCorners.resize(V);
        for (int i = 0; i < V; ++i) {
            int i1 = (i + 1) % V;
            m_midCorners[i] = m_corners[i] + m_corners[i1] - m_x0;
        }
    } else if (m_midCorners.size() != static_cast<size_t>(V)) {
        m_midCorners.resize(V);
        for (int i = 0; i < V; ++i) {
            T th_mid = 2.0 * M_PI * (i + 0.5) / V;
            m_midCorners[i].resize(2);
            m_midCorners[i] << m_midRadius * std::cos(th_mid), m_midRadius * std::sin(th_mid);
        }
    }

    m_mp_coarse.clear();
    gsKnotVector<T> kv1(0.0, 1.0, 0, 2);
    gsTensorBSplineBasis<2, T> b_bi(kv1, kv1);

    for (int i = 0; i < V; ++i) {
        int i1 = (i + 1) % V;
        gsMatrix<T> coefs(4, 2);
        coefs.row(0) = m_x0.transpose();
        coefs.row(1) = m_corners[i].transpose();
        coefs.row(2) = m_corners[i1].transpose();
        coefs.row(3) = m_midCorners[i].transpose();
        m_mp_coarse.addPatch(gsTensorBSpline<2, T>(b_bi, coefs));
    }
    m_mp_coarse.computeTopology();
}

template <class T>
void gsAsG1GeometryGenerator<T>::extractGluingData()
{
    const int V = m_valence;
    m_gd.resize(V);

    for (int i = 0; i < V; ++i) {
        int i1 = (i + 1) % V;
        int i2 = (i + 2) % V;

        gsVector<T> P00 = m_x0, P10 = m_corners[i], P01 = m_corners[i1], P11 = m_midCorners[i];
        gsVector<T> Q00 = m_x0, Q10 = m_corners[i1], Q01 = m_corners[i2], Q11 = m_midCorners[i1];

        // Derivatives at xi = 0 (central extraordinary vertex)
        gsVector<T> d1G1_0 = P10 - P00;
        gsVector<T> d2G1_0 = P01 - P00;
        gsVector<T> d1G2_0 = Q10 - Q00;
        gsVector<T> d2G2_0 = Q01 - Q00;

        T D1_0 = det2D(d1G1_0, d2G1_0);
        T D2_0 = det2D(d1G2_0, d2G2_0);
        T C_0  = det2D(d1G1_0, d2G2_0);

        // Derivatives at xi = 1 (outer boundary interface corner)
        gsVector<T> d1G1_1 = P11 - P01;
        gsVector<T> d2G1_1 = P01 - P00;
        gsVector<T> d1G2_1 = Q10 - Q00;
        gsVector<T> d2G2_1 = Q11 - Q10;

        T D1_1 = det2D(d1G1_1, d2G1_1);
        T D2_1 = det2D(d1G2_1, d2G2_1);
        T C_1  = det2D(d1G1_1, d2G2_1);

        T denom_0 = D1_0 * D1_0 + D2_0 * D2_0;

        if (m_useVirtualParallelogram) {
            // Parallelogram completion guarantees D1 and D2 are constant on G0,
            // and setting beta(1) = 0 eliminates outer boundary tangent locking.
            m_gd[i].D1_0 = D1_0; m_gd[i].D1_1 = D1_0;
            m_gd[i].D2_0 = D2_0; m_gd[i].D2_1 = D2_0;
            m_gd[i].b1_0 = C_0 * D2_0 / denom_0;
            m_gd[i].b2_0 = C_0 * D1_0 / denom_0;
            m_gd[i].b1_1 = T(0);
            m_gd[i].b2_1 = T(0);
        } else {
            m_gd[i].D1_0 = D1_0; m_gd[i].D1_1 = D1_1;
            m_gd[i].D2_0 = D2_0; m_gd[i].D2_1 = D2_1;
            T denom_1 = D1_1 * D1_1 + D2_1 * D2_1;
            m_gd[i].b1_0 = C_0 * D2_0 / denom_0;
            m_gd[i].b2_0 = C_0 * D1_0 / denom_0;
            m_gd[i].b1_1 = C_1 * D2_1 / denom_1;
            m_gd[i].b2_1 = C_1 * D1_1 / denom_1;
        }
    }
}

template <class T>
void gsAsG1GeometryGenerator<T>::refineGeometry()
{
    gsKnotVector<T> kv(0.0, 1.0, 0, m_degree + 1);
    for (int r = 0; r < m_numRefine; ++r) {
        kv.uniformRefine();
    }
    m_basis1d = gsBSplineBasis<T>(kv);
    m_patchBasis = gsTensorBSplineBasis<2, T>(kv, kv);

    m_Nu = m_basis1d.size();
    m_Nv = m_basis1d.size();
    m_Npatch = m_Nu * m_Nv;
    m_Ntotal = m_valence * m_Npatch;
}

template <class T>
void gsAsG1GeometryGenerator<T>::assembleStiffnessMatrix()
{
    const int V = m_valence;
    gsSparseMatrix<T> K_patch(m_Npatch, m_Npatch);

    gsGaussRule<T> rule1d(m_degree + 2);
    gsMatrix<T> n1d; gsVector<T> w1d;
    rule1d.mapToAll(m_basis1d.knots().breaks(), n1d, w1d);
    const int nq = n1d.cols();

    gsMatrix<T> pts(2, nq * nq);
    gsVector<T> w2d(nq * nq);
    int idx = 0;
    for (int i = 0; i < nq; ++i) {
        for (int j = 0; j < nq; ++j) {
            pts(0, idx) = n1d(0, i);
            pts(1, idx) = n1d(0, j);
            w2d(idx) = w1d(i) * w1d(j);
            idx++;
        }
    }

    gsMatrix<index_t> act2d = m_patchBasis.active(pts);
    gsMatrix<T> d1, d2;
    m_patchBasis.deriv_into(pts, d1);
    m_patchBasis.deriv2_into(pts, d2);

    gsMatrix<T> K_dense(m_Npatch, m_Npatch);
    K_dense.setZero();

    for (int q = 0; q < idx; ++q) {
        T w = w2d(q);
        const int nAct = act2d.rows();
        for (int k1 = 0; k1 < nAct; ++k1) {
            int a = act2d(k1, q);
            T da_u  = d1(2 * k1 + 0, q);
            T da_v  = d1(2 * k1 + 1, q);
            T da_uu = d2(3 * k1 + 0, q);
            T da_vv = d2(3 * k1 + 1, q);
            T lap_a = da_uu + da_vv;

            for (int k2 = 0; k2 < nAct; ++k2) {
                int b_idx = act2d(k2, q);
                T db_u  = d1(2 * k2 + 0, q);
                T db_v  = d1(2 * k2 + 1, q);
                T db_uu = d2(3 * k2 + 0, q);
                T db_vv = d2(3 * k2 + 1, q);
                T lap_b = db_uu + db_vv;

                T val = m_laplaceWeight * (da_u * db_u + da_v * db_v)
                      + m_biharmWeight  * (lap_a * lap_b);
                K_dense(a, b_idx) += w * val;
            }
        }
    }

    for (int a = 0; a < m_Npatch; ++a)
        for (int b_idx = 0; b_idx < m_Npatch; ++b_idx)
            if (std::abs(K_dense(a, b_idx)) > T(1e-15))
                K_patch.insert(a, b_idx) = K_dense(a, b_idx);
    K_patch.makeCompressed();

    m_K_global = gsSparseMatrix<T>(m_Ntotal, m_Ntotal);
    for (int p_idx = 0; p_idx < V; ++p_idx) {
        int offset = p_idx * m_Npatch;
        for (int k = 0; k < K_patch.outerSize(); ++k) {
            for (typename gsSparseMatrix<T>::InnerIterator it(K_patch, k); it; ++it) {
                m_K_global.insert(offset + it.row(), offset + it.col()) = it.value();
            }
        }
    }
    m_K_global.makeCompressed();
}

template <class T>
void gsAsG1GeometryGenerator<T>::assembleConstraints(
    const std::vector<CurveFunction>& eastCurves,
    const std::vector<CurveFunction>& northCurves)
{
    const int V = m_valence;
    m_constr_rows.clear();
    m_rhs_x.clear();
    m_rhs_y.clear();

    // 1. Central Vertex
    for (int i = 0; i < V; ++i) {
        addConstraintRow({{globalIdx(i, 0, 0), T(1)}}, m_x0(0), m_x0(1));
    }

    // 2. Argyris C2 Vertex Decoupling: d_uu G(0,0) = 0, d_vv G(0,0) = 0
    gsMatrix<T> pt0(2, 1); pt0.setZero();
    gsMatrix<index_t> act0 = m_patchBasis.active(pt0);
    gsMatrix<T> d2_0;
    m_patchBasis.deriv2_into(pt0, d2_0);
    for (int i = 0; i < V; ++i) {
        std::vector<std::pair<int, T>> row_uu, row_vv;
        for (int k = 0; k < act0.rows(); ++k) {
            int a = act0(k, 0);
            T val_uu = d2_0(3 * k + 0, 0);
            T val_vv = d2_0(3 * k + 1, 0);
            if (std::abs(val_uu) > T(1e-14))
                row_uu.push_back({i * m_Npatch + a, val_uu});
            if (std::abs(val_vv) > T(1e-14))
                row_vv.push_back({i * m_Npatch + a, val_vv});
        }
        addConstraintRow(row_uu, T(0), T(0));
        addConstraintRow(row_vv, T(0), T(0));
    }

    // 3. Outer Boundary Constraints
    std::mt19937 rng_bump(m_seed + 100);
    std::uniform_real_distribution<double> dist_bump(0.1, 0.25);

    for (int i = 0; i < V; ++i) {
        int i1 = (i + 1) % V;

        CurveFunction eastFunc = nullptr;
        if (eastCurves.size() == static_cast<size_t>(V) && eastCurves[i]) {
            eastFunc = eastCurves[i];
        } else {
            T bump_coeff = m_randomize ? T(dist_bump(rng_bump)) : T(0.4);
            eastFunc = [this, i, bump_coeff](T t) -> gsVector<T> {
                gsVector<T> p_lin = (T(1) - t) * m_corners[i] + t * m_midCorners[i];
                T bump = bump_coeff * t * t * (T(1) - t);
                gsVector<T> norm(2);
                norm << -(m_midCorners[i](1) - m_corners[i](1)), (m_midCorners[i](0) - m_corners[i](0));
                norm.normalize();
                return p_lin + bump * norm;
            };
        }

        CurveFunction northFunc = nullptr;
        if (northCurves.size() == static_cast<size_t>(V) && northCurves[i]) {
            northFunc = northCurves[i];
        } else {
            T bump_coeff = m_randomize ? T(dist_bump(rng_bump)) : T(0.4);
            northFunc = [this, i, i1, bump_coeff](T t) -> gsVector<T> {
                gsVector<T> p_lin = (T(1) - t) * m_corners[i1] + t * m_midCorners[i];
                T bump = bump_coeff * t * t * (T(1) - t);
                gsVector<T> norm(2);
                norm << (m_midCorners[i](1) - m_corners[i1](1)), -(m_midCorners[i](0) - m_corners[i1](0));
                norm.normalize();
                return p_lin + bump * norm;
            };
        }

        gsMatrix<T> P_east;
        if (m_custom_east_cp.size() == static_cast<size_t>(V) && m_custom_east_cp[i].rows() == m_Nv) {
            P_east = m_custom_east_cp[i];
        } else {
            P_east = fitBoundaryCurve(m_basis1d, eastFunc);
        }

        gsMatrix<T> P_north;
        if (m_custom_north_cp.size() == static_cast<size_t>(V) && m_custom_north_cp[i].rows() == m_Nu) {
            P_north = m_custom_north_cp[i];
        } else {
            P_north = fitBoundaryCurve(m_basis1d, northFunc);
        }

        for (int iv = 0; iv < m_Nv; ++iv) {
            addConstraintRow({{globalIdx(i, m_Nu - 1, iv), T(1)}}, P_east(iv, 0), P_east(iv, 1));
        }
        for (int iu = 0; iu < m_Nu - 1; ++iu) {
            addConstraintRow({{globalIdx(i, iu, m_Nv - 1), T(1)}}, P_north(iu, 0), P_north(iu, 1));
        }
    }

    // 4. C0 Interface Matching: G^(i)(0, k) - G^(i1)(k, 0) = 0
    for (int i = 0; i < V; ++i) {
        int i1 = (i + 1) % V;
        for (int k = 1; k < m_Nv - 1; ++k) {
            addConstraintRow({{globalIdx(i, 0, k), T(1)}, {globalIdx(i1, k, 0), T(-1)}}, T(0), T(0));
        }
    }

    // 5. AS-G1 Interface Constraints via test basis projection
    gsBSplineBasis<T> b_test = m_basis1d;
    b_test.degreeElevate(1);
    const int nTest = b_test.size();

    gsGaussRule<T> rule_interf(b_test.degree() + 2);
    gsMatrix<T> qnodes; gsVector<T> qweights;
    rule_interf.mapToAll(b_test.knots().breaks(), qnodes, qweights);
    const int NQ = qnodes.cols();

    gsMatrix<index_t> test_act = b_test.active(qnodes);
    gsMatrix<T> testVals;
    b_test.eval_into(qnodes, testVals);

    gsMatrix<T> pts_p1(2, NQ), pts_p2(2, NQ);
    for (int q = 0; q < NQ; ++q) {
        T xi = qnodes(0, q);
        pts_p1(0, q) = 0.0; pts_p1(1, q) = xi;
        pts_p2(0, q) = xi;  pts_p2(1, q) = 0.0;
    }

    gsMatrix<index_t> act_p1 = m_patchBasis.active(pts_p1);
    gsMatrix<index_t> act_p2 = m_patchBasis.active(pts_p2);
    gsMatrix<T> d1_p1, d1_p2;
    m_patchBasis.deriv_into(pts_p1, d1_p1);
    m_patchBasis.deriv_into(pts_p2, d1_p2);

    for (int i = 0; i < V; ++i) {
        int i1 = (i + 1) % V;
        const AsG1GluingData<T>& g = m_gd[i];
        std::vector<std::map<int, T>> test_rows(nTest);

        for (int q = 0; q < NQ; ++q) {
            T w = qweights(q);
            T xi = qnodes(0, q);

            T a1 = g.alpha1(xi);
            T a2 = g.alpha2(xi);
            T b_glue = g.b_glue(xi);

            for (int k_test = 0; k_test < test_act.rows(); ++k_test) {
                int a = test_act(k_test, q);
                T phi_a = testVals(k_test, q);

                T factor_d1_p1 = w * phi_a * a2;
                T factor_d2_p1 = w * phi_a * (-b_glue);
                T factor_d2_p2 = w * phi_a * a1;

                for (int k = 0; k < act_p1.rows(); ++k) {
                    int j1 = act_p1(k, q);
                    T du1 = d1_p1(2 * k + 0, q);
                    T dv1 = d1_p1(2 * k + 1, q);
                    T term_p1 = factor_d1_p1 * du1 + factor_d2_p1 * dv1;
                    if (std::abs(term_p1) > T(1e-15)) {
                        int gcol = i * m_Npatch + j1;
                        test_rows[a][gcol] += term_p1;
                    }
                }

                for (int k = 0; k < act_p2.rows(); ++k) {
                    int j2 = act_p2(k, q);
                    T dv2 = d1_p2(2 * k + 1, q);
                    T term_p2 = factor_d2_p2 * dv2;
                    if (std::abs(term_p2) > T(1e-15)) {
                        int gcol = i1 * m_Npatch + j2;
                        test_rows[a][gcol] += term_p2;
                    }
                }
            }
        }

        for (int a = 0; a < nTest; ++a) {
            std::vector<std::pair<int, T>> row;
            for (auto const& kv_pair : test_rows[a]) {
                if (std::abs(kv_pair.second) > T(1e-14))
                    row.push_back({kv_pair.first, kv_pair.second});
            }
            if (!row.empty())
                addConstraintRow(row, T(0), T(0));
        }
    }
}

template <class T>
bool gsAsG1GeometryGenerator<T>::solveKKT()
{
    const int M_all = m_constr_rows.size();

    // Dense conversion of constraints to filter redundant rows via QR
    gsMatrix<T> C_all(M_all, m_Ntotal);
    C_all.setZero();
    for (int r = 0; r < M_all; ++r) {
        for (auto const& term : m_constr_rows[r]) {
            C_all(r, term.first) = term.second;
        }
    }

    gsMatrix<T> CT = C_all.transpose();
    auto qr = CT.colPivHouseholderQr();
    int rank = qr.rank();

    auto P_perm = qr.colsPermutation();
    std::vector<int> indep_rows;
    indep_rows.reserve(rank);
    for (int j = 0; j < rank; ++j) {
        indep_rows.push_back(P_perm.indices()(j));
    }

    // Build Exact KKT System
    const int KKT_size = m_Ntotal + rank;
    gsSparseMatrix<T> KKT(KKT_size, KKT_size);

    for (int k = 0; k < m_K_global.outerSize(); ++k) {
        for (typename gsSparseMatrix<T>::InnerIterator it(m_K_global, k); it; ++it) {
            KKT.insert(it.row(), it.col()) = it.value();
        }
    }

    gsMatrix<T> RHS(KKT_size, 2);
    RHS.setZero();

    for (int i = 0; i < rank; ++i) {
        int orig_r = indep_rows[i];
        for (auto const& term : m_constr_rows[orig_r]) {
            int col = term.first;
            T val = term.second;
            KKT.insert(m_Ntotal + i, col) = val;
            KKT.insert(col, m_Ntotal + i) = val;
        }
        RHS(m_Ntotal + i, 0) = m_rhs_x[orig_r];
        RHS(m_Ntotal + i, 1) = m_rhs_y[orig_r];
    }
    KKT.makeCompressed();

    // Sparse LU Solve
    typename gsSparseSolver<T>::LU solver(KKT);
    if (!solver.succeed()) {
        gsWarn << "KKT factorization failed!\n";
        return false;
    }

    gsMatrix<T> SOL = solver.solve(RHS);
    if (!solver.succeed()) {
        gsWarn << "KKT solve failed!\n";
        return false;
    }

    gsMatrix<T> P_sol = SOL.topRows(m_Ntotal);

    // Update MultiPatch
    m_mp.clear();
    for (int i = 0; i < m_valence; ++i) {
        gsMatrix<T> p_coefs(m_Npatch, 2);
        for (int a = 0; a < m_Npatch; ++a) {
            p_coefs(a, 0) = P_sol(i * m_Npatch + a, 0);
            p_coefs(a, 1) = P_sol(i * m_Npatch + a, 1);
        }
        m_mp.addPatch(gsTensorBSpline<2, T>(m_patchBasis, p_coefs));
    }
    m_mp.computeTopology();

    return true;
}

template <class T>
T gsAsG1GeometryGenerator<T>::verifyAsG1(
    int numSamples,
    std::vector<T>* perInterfaceResiduals) const
{
    const int V = m_valence;
    T global_max_res = 0;
    if (perInterfaceResiduals) perInterfaceResiduals->resize(V);

    gsVector<T> s_pts = gsVector<T>::LinSpaced(numSamples, 0.0, 1.0);

    for (int i = 0; i < V; ++i) {
        int i1 = (i + 1) % V;
        const AsG1GluingData<T>& g = m_gd[i];
        T max_res_i = 0;

        for (int s = 0; s < numSamples; ++s) {
            T xi = s_pts(s);
            gsMatrix<T> pt1(2, 1), pt2(2, 1);
            pt1 << 0.0, xi;
            pt2 << xi, 0.0;

            gsMatrix<T> d_p1 = m_mp.patch(i).deriv(pt1);
            gsMatrix<T> d_p2 = m_mp.patch(i1).deriv(pt2);

            gsVector<T> d1G1(2), d2G1(2), d2G2(2);
            d1G1 << d_p1(0, 0), d_p1(2, 0);
            d2G1 << d_p1(1, 0), d_p1(3, 0);
            d2G2 << d_p2(1, 0), d_p2(3, 0);

            T a1 = g.alpha1(xi);
            T a2 = g.alpha2(xi);
            T b_glue = g.b_glue(xi);

            gsVector<T> R_vec = a2 * d1G1 + a1 * d2G2 - b_glue * d2G1;
            T res_norm = R_vec.norm();
            if (res_norm > max_res_i) max_res_i = res_norm;
        }

        if (perInterfaceResiduals) (*perInterfaceResiduals)[i] = max_res_i;
        if (max_res_i > global_max_res) global_max_res = max_res_i;
    }

    return global_max_res;
}

template <class T>
bool gsAsG1GeometryGenerator<T>::checkJacobianRegularity(
    int samplesPerDir,
    T* outMinDet,
    T* outMaxDet) const
{
    const int V = m_valence;
    T min_det_all = 1e12, max_det_all = -1e12;
    bool all_positive = true;

    for (int i = 0; i < V; ++i) {
        T min_det_patch = 1e12;
        T min_u = 0, min_v = 0;
        for (int iu = 0; iu <= samplesPerDir; ++iu) {
            for (int iv = 0; iv <= samplesPerDir; ++iv) {
                T u_val = static_cast<T>(iu) / samplesPerDir;
                T v_val = static_cast<T>(iv) / samplesPerDir;
                gsMatrix<T> pt(2, 1);
                pt << u_val, v_val;
                gsMatrix<T> d = m_mp.patch(i).deriv(pt);
                gsVector<T> du(2), dv(2);
                du << d(0, 0), d(2, 0);
                dv << d(1, 0), d(3, 0);
                T detJ = det2D(du, dv);

                if (std::isnan(detJ)) {
                    all_positive = false;
                } else {
                    if (detJ < min_det_patch) {
                        min_det_patch = detJ;
                        min_u = u_val;
                        min_v = v_val;
                    }
                    if (detJ < min_det_all) min_det_all = detJ;
                    if (detJ > max_det_all) max_det_all = detJ;
                    if (detJ <= T(0)) all_positive = false;
                }
            }
        }
        std::cout << "  Patch " << i << " min det(J) = " << min_det_patch << " at (u,v)=(" << min_u << "," << min_v << ")" << std::endl;
    }

    if (outMinDet) *outMinDet = min_det_all;
    if (outMaxDet) *outMaxDet = max_det_all;
    return all_positive;
}

template <class T>
gsMatrix<T> fitBoundarySamplePoints(const gsBSplineBasis<T>& basis, const gsMatrix<T>& samples)
{
    const int N = samples.rows();
    const int nBasis = basis.size();
    if (N == nBasis) return samples;

    gsMatrix<T> u(1, N);
    for (int j = 0; j < N; ++j) {
        u(0, j) = static_cast<T>(j) / (N - 1);
    }
    gsMatrix<T> B;
    basis.eval_into(u, B);
    gsMatrix<T> A = B * B.transpose();
    gsMatrix<T> RHS = B * samples;
    gsMatrix<T> P = A.colPivHouseholderQr().solve(RHS);
    P.row(0) = samples.row(0);
    P.row(nBasis - 1) = samples.row(N - 1);
    return P;
}

template <class T>
bool gsAsG1GeometryGenerator<T>::loadBoundaryDefinition(const std::string& filename)
{
    std::ifstream in(filename.c_str());
    if (!in.is_open()) {
        gsWarn << "Cannot open boundary definition file: " << filename << "\n";
        return false;
    }

    std::string token;
    int V = 0;
    int deg = m_degree;
    std::vector<gsVector<T>> corners;
    std::vector<gsVector<T>> midCorners;
    std::vector<gsMatrix<T>> customEast;
    std::vector<gsMatrix<T>> customNorth;
    int curPatch = -1;

    while (in >> token) {
        if (token.empty() || token[0] == '#') {
            std::string line;
            std::getline(in, line);
            continue;
        }
        if (token == "VALENCE") {
            in >> V;
            customEast.resize(V);
            customNorth.resize(V);
        } else if (token == "DEGREE") {
            in >> deg;
        } else if (token == "CORNERS") {
            corners.resize(V);
            for (int i = 0; i < V; ++i) {
                corners[i].resize(2);
                in >> corners[i](0) >> corners[i](1);
            }
        } else if (token == "MIDCORNERS") {
            midCorners.resize(V);
            for (int i = 0; i < V; ++i) {
                midCorners[i].resize(2);
                in >> midCorners[i](0) >> midCorners[i](1);
            }
        } else if (token == "PATCH") {
            in >> curPatch;
        } else if (token == "EAST_CONTROL_POINTS" && curPatch >= 0 && curPatch < V) {
            int n;
            in >> n;
            customEast[curPatch].resize(n, 2);
            for (int j = 0; j < n; ++j) {
                in >> customEast[curPatch](j, 0) >> customEast[curPatch](j, 1);
            }
        } else if (token == "NORTH_CONTROL_POINTS" && curPatch >= 0 && curPatch < V) {
            int n;
            in >> n;
            customNorth[curPatch].resize(n, 2);
            for (int j = 0; j < n; ++j) {
                in >> customNorth[curPatch](j, 0) >> customNorth[curPatch](j, 1);
            }
        } else if (token == "EAST_SAMPLES" && curPatch >= 0 && curPatch < V) {
            int n;
            in >> n;
            gsMatrix<T> samples(n, 2);
            for (int j = 0; j < n; ++j) in >> samples(j, 0) >> samples(j, 1);
            gsBSplineBasis<T> b1d(gsKnotVector<T>(0.0, 1.0, 0, deg + 1));
            customEast[curPatch] = fitBoundarySamplePoints(b1d, samples);
        } else if (token == "NORTH_SAMPLES" && curPatch >= 0 && curPatch < V) {
            int n;
            in >> n;
            gsMatrix<T> samples(n, 2);
            for (int j = 0; j < n; ++j) in >> samples(j, 0) >> samples(j, 1);
            gsBSplineBasis<T> b1d(gsKnotVector<T>(0.0, 1.0, 0, deg + 1));
            customNorth[curPatch] = fitBoundarySamplePoints(b1d, samples);
        }
    }

    if (corners.empty() && customEast.size() == static_cast<size_t>(V)) {
        corners.resize(V);
        midCorners.resize(V);
        for (int i = 0; i < V; ++i) {
            if (customEast[i].rows() > 0) {
                corners[i].resize(2);
                corners[i](0) = customEast[i](0, 0);
                corners[i](1) = customEast[i](0, 1);
                midCorners[i].resize(2);
                midCorners[i](0) = customEast[i](customEast[i].rows() - 1, 0);
                midCorners[i](1) = customEast[i](customEast[i].rows() - 1, 1);
            }
        }
    }

    if (V <= 0 || corners.size() != static_cast<size_t>(V)) {
        gsWarn << "Invalid boundary file: valence or corners missing!\n";
        return false;
    }

    m_valence = V;
    m_degree = deg;
    m_corners = corners;
    m_midCorners = midCorners;
    m_custom_east_cp = customEast;
    m_custom_north_cp = customNorth;
    return true;
}

} // namespace gismo
