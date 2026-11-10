/** @file gsASG1LocalFitter.h
    @brief Header for local AS-G1 multi-patch spline surface fitter.
    Based on: Farahat, Kapl, Kosmac, Vitrih (2024), Computers and Mathematics with Applications 168, 46-57.
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMultiPatch.h>
#include <gsNurbs/gsTensorBSpline.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsModeling/gsAsG1Domain.hpp>

namespace gismo
{

/**
 * @brief Local analysis-suitable G^1 (AS-G^1) multi-patch fitting algorithm.
 *
 * Implements the 3-step decoupled local Lagrange multiplier fitting algorithm from:
 * "A locally based construction of analysis-suitable G^1 multi-patch spline surfaces"
 * (Farahat, Kapl, Kosmač, Vitrih, 2024).
 *
 * Given a conforming or G^1-smooth target multi-patch surface S, this fitter
 * constructs a new multi-patch surface F in S^{p,r}_k such that:
 *  - Interfaces satisfy the exact linear AS-G^1 gluing conditions:
 *      \alpha^{(1)} f_1 = \partial_n F^{(1)} + \beta^{(1)} \partial_t F^{(1)}
 *     -\alpha^{(2)} f_1 = \partial_n F^{(2)} - \beta^{(2)} \partial_t F^{(2)}
 *  - Vertex neighborhoods match C^1 / G^1 constraints with valence \nu >= 2.
 *  - Patch interiors and outer boundaries are fitted via local H^1 projections.
 *
 * @tparam T Floating-point scalar type (e.g., real_t, double).
 */
template <class T>
class gsASG1LocalFitter
{
public:
    enum class FittingMode { Global, Local };

    struct SideInfo {
        short_t nDir, tDir;
        bool par;
        T nSign;
    };

    struct GluingData {
        T a0_1, a1_1, b0_1, b1_1;
        T a0_2, a1_2, b0_2, b1_2;
        bool reversed1, reversed2;
        patchSide ps1, ps2;
        gsVector<T> P0, P1;

        T alpha1(T xi) const { return a0_1 * (T(1) - xi) + a1_1 * xi; }
        T alpha2(T xi) const { return a0_2 * (T(1) - xi) + a1_2 * xi; }
        T beta1(T xi)  const { return b0_1 * (T(1) - xi) + b1_1 * xi; }
        T beta2(T xi)  const { return b0_2 * (T(1) - xi) + b1_2 * xi; }
    };

    struct OrderedVertexPatch {
        patchCorner pc;
        patchSide side_prev;
        patchSide side_next;
        bool is_boundary_prev;
        bool is_boundary_next;
        int dir_prev, dir_next;
        T sign_prev, sign_next;
    };

public:
    gsASG1LocalFitter(const gsMultiPatch<T>& target_surface,
                      short_t degree = 3,
                      short_t regularity = 1,
                      index_t num_inner_knots = 2);

    void fit();

    FittingMode mode() const { return m_mode; }
    void setMode(FittingMode mode) { m_mode = mode; }

    const gsMultiPatch<T>& result() const { return m_result; }
    gsMultiPatch<T> giveResult() { return give(m_result); }

    bool checkASG1(T eps = 1e-4) const;

    T maxAlphaResidual() const { return m_max_res_a; }
    T maxBetaResidual()  const { return m_max_res_b; }

    const std::vector<GluingData>& gluingData() const { return m_gluing; }

    static SideInfo getSideInfo(boxSide s);

    static void mapInterfaceParam(const patchSide& ps, bool reversed, T xi,
                                  T& u, T& v, T& tSign);

    static void evalSTangAndNorm(const gsGeometry<T>& patch, const patchSide& ps,
                                 bool reversed, T xi,
                                 gsVector<T>& nVec, gsVector<T>& tVec);

    static gsSparseMatrix<T> collocationMatrix(const gsBSplineBasis<T> &basis,
                                               const gsMatrix<T> &pts,
                                               index_t derivative = 0);

protected:
    void fitGlobal();
    void fitLocal();

    void step0_computeGluingFunctions();
    void step1_solveVertexProblems();
    void step2_solveInterfaceProblems();
    void step2b_solveBoundaryEdges();
    void step3_solvePatchInteriors();

    std::vector<OrderedVertexPatch> orderVertexPatches(const std::vector<patchCorner>& cornerList) const;

private:
    const gsMultiPatch<T>& m_target;
    short_t m_p;
    short_t m_r;
    index_t m_k;
    T m_sigma;

    index_t m_n;  // Basis size along 1D: p + 1 + k*(p - r)
    index_t m_n2; // 2D basis size: m_n * m_n

    gsKnotVector<T> m_kv;
    gsMultiPatch<T> m_result;
    std::vector<GluingData> m_gluing;

    std::vector<gsMatrix<T>> m_patchCP;
    std::vector<std::vector<bool>> m_patchCP_fixed;

    mutable T m_max_res_a;
    mutable T m_max_res_b;

    FittingMode m_mode = FittingMode::Global;
};

} // namespace gismo

#include <gsModeling/gsASG1LocalFitter.hpp>

