/** @file gsASG1LocalFitter.hpp
    @brief Implementation of local AS-G1 multi-patch spline surface fitter.
    Based on: Farahat, Kapl, Kosmac, Vitrih (2024), Computers and Mathematics with Applications 168, 46-57.
*/

#pragma once

#include <gsModeling/gsASG1LocalFitter.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <queue>
#include <algorithm>
#include <cmath>

namespace gismo
{

template <class T>
typename gsASG1LocalFitter<T>::SideInfo gsASG1LocalFitter<T>::getSideInfo(boxSide s)
{
    SideInfo si;
    si.nDir = s.direction();
    si.tDir = 1 - si.nDir;
    si.par = s.parameter();
    si.nSign = si.par ? T(-1) : T(1);
    return si;
}

template <class T>
void gsASG1LocalFitter<T>::mapInterfaceParam(const patchSide& ps, bool reversed, T xi,
                                             T& u, T& v, T& tSign)
{
    SideInfo si = getSideInfo(ps.side());
    T t = reversed ? (T(1) - xi) : xi;
    tSign = reversed ? T(-1) : T(1);
    if (si.nDir == 0) {
        u = si.par ? T(1) : T(0);
        v = t;
    } else {
        u = t;
        v = si.par ? T(1) : T(0);
    }
}

template <class T>
void gsASG1LocalFitter<T>::evalSTangAndNorm(const gsGeometry<T>& patch, const patchSide& ps,
                                          bool reversed, T xi,
                                          gsVector<T>& nVec, gsVector<T>& tVec)
{
    T u, v, tSign;
    mapInterfaceParam(ps, reversed, xi, u, v, tSign);
    SideInfo si = getSideInfo(ps.side());

    gsMatrix<T> pt(2, 1);
    pt << u, v;
    gsMatrix<T> d = patch.deriv(pt);
    const short_t geoDim = patch.targetDim();
    gsVector<T> du(geoDim), dv(geoDim);
    for (short_t c = 0; c < geoDim; ++c) {
        du(c) = d(2 * c + 0, 0);
        dv(c) = d(2 * c + 1, 0);
    }

    if (si.nDir == 0) {
        nVec = si.nSign * du;
        tVec = tSign * dv;
    } else {
        nVec = si.nSign * dv;
        tVec = tSign * du;
    }
}

template <class T>
gsSparseMatrix<T> gsASG1LocalFitter<T>::collocationMatrix(const gsBSplineBasis<T> &basis,
                                                         const gsMatrix<T> &pts,
                                                         index_t derivative)
{
    gsMatrix<T> vals;
    if (derivative == 0)
        vals = basis.eval(pts);
    else if (derivative == 1)
        vals = basis.deriv(pts);
    else if (derivative == 2)
        vals = basis.deriv2(pts);
    else
        GISMO_ERROR("Invalid derivative order.");

    gsMatrix<index_t> idx = basis.active(pts);
    gsSparseEntries<T> entries;
    entries.reserve(vals.rows() * vals.cols());
    for (int i = 0; i < vals.cols(); ++i)
        for (int j = 0; j < vals.rows(); ++j)
            entries.add(i, idx(j, i), vals(j, i));
    gsSparseMatrix<T> C(pts.cols(), basis.size());
    C.setFrom(entries);
    return C;
}

template <class T>
gsASG1LocalFitter<T>::gsASG1LocalFitter(const gsMultiPatch<T>& target_surface,
                                       short_t degree,
                                       short_t regularity,
                                       index_t num_inner_knots)
    : m_target(target_surface),
      m_geoDim(target_surface.targetDim()),
      m_p(degree),
      m_r(regularity),
      m_k(num_inner_knots),
      m_sigma(T(1) / (degree * (num_inner_knots + 1))),
      m_max_res_a(T(0)),
      m_max_res_b(T(0)),
      m_max_res_n(T(0))
{
    GISMO_ASSERT(m_geoDim == 2 || m_geoDim == 3, "gsASG1LocalFitter supports targetDim 2 or 3.");
    GISMO_ASSERT(m_p >= 3, "Degree p must be >= 3.");
    GISMO_ASSERT(m_r >= 1 && m_r <= m_p - 2, "Regularity r must satisfy 1 <= r <= p - 2.");
    GISMO_ASSERT(m_k >= (5 - m_p + (m_p - m_r - 2)) / (m_p - m_r - 1), "Insufficient inner knots.");

    m_kv = gsKnotVector<T>(T(0), T(1), m_k, m_p + 1, m_p - m_r);
    gsBSplineBasis<T> b1D(m_kv);
    m_n = b1D.size();
    m_n2 = m_n * m_n;

    index_t nP = m_target.nPatches();
    m_patchCP.resize(nP, gsMatrix<T>::Zero(m_n2, m_geoDim));
    m_patchCP_fixed.resize(nP, std::vector<bool>(m_n2, false));
}

template <class T>
std::vector<typename gsASG1LocalFitter<T>::OrderedVertexPatch>
gsASG1LocalFitter<T>::orderVertexPatches(const std::vector<patchCorner>& cornerList) const
{
    int val = cornerList.size();
    if (val == 1) {
        std::vector<patchSide> sides;
        cornerList[0].getContainingSides(2, sides);
        OrderedVertexPatch ovp;
        ovp.pc = cornerList[0];
        ovp.side_prev = sides[0]; ovp.is_boundary_prev = true;
        ovp.side_next = sides[1]; ovp.is_boundary_next = true;
        ovp.dir_prev = getSideInfo(sides[0].side()).tDir;
        ovp.dir_next = getSideInfo(sides[1].side()).tDir;
        bool u_par = (ovp.pc.m_index - 1) & 1;
        bool v_par = ((ovp.pc.m_index - 1) >> 1) & 1;
        ovp.sign_prev = (ovp.dir_prev == 0) ? (u_par ? T(-1) : T(1)) : (v_par ? T(-1) : T(1));
        ovp.sign_next = (ovp.dir_next == 0) ? (u_par ? T(-1) : T(1)) : (v_par ? T(-1) : T(1));
        return {ovp};
    }

    std::map<std::pair<int, int>, patchSide> ifMap;
    for (auto it = m_target.iBegin(); it != m_target.iEnd(); ++it) {
        ifMap[{it->first().patch, it->first().side().m_index}] = it->second();
        ifMap[{it->second().patch, it->second().side().m_index}] = it->first();
    }

    struct Node {
        patchCorner pc;
        patchSide s0, s1;
        int neighbor0, neighbor1;
    };
    std::vector<Node> nodes(val);
    for (int i = 0; i < val; ++i) {
        nodes[i].pc = cornerList[i];
        std::vector<patchSide> sides;
        cornerList[i].getContainingSides(2, sides);
        nodes[i].s0 = sides[0];
        nodes[i].s1 = sides[1];
        nodes[i].neighbor0 = -1;
        nodes[i].neighbor1 = -1;

        auto it0 = ifMap.find({sides[0].patch, sides[0].side().m_index});
        if (it0 != ifMap.end()) {
            for (int j = 0; j < val; ++j) {
                if (j != i && cornerList[j].patch == it0->second.patch) {
                    nodes[i].neighbor0 = j;
                    break;
                }
            }
        }
        auto it1 = ifMap.find({sides[1].patch, sides[1].side().m_index});
        if (it1 != ifMap.end()) {
            for (int j = 0; j < val; ++j) {
                if (j != i && cornerList[j].patch == it1->second.patch) {
                    nodes[i].neighbor1 = j;
                    break;
                }
            }
        }
    }

    int startIdx = 0;
    bool is_boundary_vertex = false;
    for (int i = 0; i < val; ++i) {
        if (nodes[i].neighbor0 == -1 || nodes[i].neighbor1 == -1) {
            startIdx = i;
            is_boundary_vertex = true;
            break;
        }
    }

    std::vector<int> chain;
    std::vector<bool> visited(val, false);
    int curr = startIdx;
    chain.push_back(curr);
    visited[curr] = true;

    int lastSideUsed = -1;
    if (is_boundary_vertex) {
        lastSideUsed = (nodes[curr].neighbor0 == -1) ? 0 : 1;
    }

    for (int step = 1; step < val; ++step) {
        int next = -1;
        int nextSideUsed = -1;
        if (lastSideUsed == 0) {
            next = nodes[curr].neighbor1;
            nextSideUsed = 1;
        } else if (lastSideUsed == 1) {
            next = nodes[curr].neighbor0;
            nextSideUsed = 0;
        } else {
            if (nodes[curr].neighbor0 != -1 && !visited[nodes[curr].neighbor0]) {
                next = nodes[curr].neighbor0;
                nextSideUsed = 0;
            } else if (nodes[curr].neighbor1 != -1 && !visited[nodes[curr].neighbor1]) {
                next = nodes[curr].neighbor1;
                nextSideUsed = 1;
            }
        }

        if (next == -1 || visited[next]) break;
        chain.push_back(next);
        visited[next] = true;

        if (nodes[next].neighbor0 == curr)
            lastSideUsed = 0;
        else if (nodes[next].neighbor1 == curr)
            lastSideUsed = 1;
        else
            lastSideUsed = -1;

        curr = next;
    }

    if ((int)chain.size() < val) {
        for (int i = 0; i < val; ++i)
            if (!visited[i]) chain.push_back(i);
    }

    std::vector<OrderedVertexPatch> ordered(val);
    for (int i = 0; i < val; ++i) {
        int idx = chain[i];
        ordered[i].pc = nodes[idx].pc;
        std::vector<patchSide> sides;
        nodes[idx].pc.getContainingSides(2, sides);

        if (i == 0 && is_boundary_vertex) {
            if (nodes[idx].neighbor0 == -1) {
                ordered[i].side_prev = nodes[idx].s0;
                ordered[i].side_next = nodes[idx].s1;
                ordered[i].is_boundary_prev = true;
                ordered[i].is_boundary_next = (nodes[idx].neighbor1 == -1);
            } else {
                ordered[i].side_prev = nodes[idx].s1;
                ordered[i].side_next = nodes[idx].s0;
                ordered[i].is_boundary_prev = true;
                ordered[i].is_boundary_next = (nodes[idx].neighbor0 == -1);
            }
        } else {
            int prevIdx = chain[(i - 1 + val) % val];
            if (nodes[idx].neighbor0 == prevIdx) {
                ordered[i].side_prev = nodes[idx].s0;
                ordered[i].side_next = nodes[idx].s1;
                ordered[i].is_boundary_prev = false;
                ordered[i].is_boundary_next = (nodes[idx].neighbor1 == -1);
            } else {
                ordered[i].side_prev = nodes[idx].s1;
                ordered[i].side_next = nodes[idx].s0;
                ordered[i].is_boundary_prev = false;
                ordered[i].is_boundary_next = (nodes[idx].neighbor0 == -1);
            }
        }

        ordered[i].dir_prev = getSideInfo(ordered[i].side_prev.side()).tDir;
        ordered[i].dir_next = getSideInfo(ordered[i].side_next.side()).tDir;

        bool u_par = (ordered[i].pc.m_index - 1) & 1;
        bool v_par = ((ordered[i].pc.m_index - 1) >> 1) & 1;
        ordered[i].sign_prev = (ordered[i].dir_prev == 0) ? (u_par ? T(-1) : T(1)) : (v_par ? T(-1) : T(1));
        ordered[i].sign_next = (ordered[i].dir_next == 0) ? (u_par ? T(-1) : T(1)) : (v_par ? T(-1) : T(1));
    }

    return ordered;
}

template <class T>
void gsASG1LocalFitter<T>::step0_computeGluingFunctions()
{
    index_t nI = m_target.nInterfaces();
    m_gluing.resize(nI);

    auto areaMag = [this](const gsVector<T>& a, const gsVector<T>& b) -> T {
        if (m_geoDim == 2)
            return std::abs(a(0) * b(1) - a(1) * b(0));
        const T cx = a(1) * b(2) - a(2) * b(1);
        const T cy = a(2) * b(0) - a(0) * b(2);
        const T cz = a(0) * b(1) - a(1) * b(0);
        return std::sqrt(cx * cx + cy * cy + cz * cz);
    };

    index_t interfIdx = 0;
    for (auto it = m_target.iBegin(); it != m_target.iEnd(); ++it, ++interfIdx) {
        const boundaryInterface& bi = *it;
        m_gluing[interfIdx].ps1 = bi.first();
        m_gluing[interfIdx].ps2 = bi.second();

        SideInfo si1 = getSideInfo(bi.first().side());
        bool flipped = !bi.dirOrientation(bi.first(), si1.tDir);
        m_gluing[interfIdx].reversed1 = false;
        m_gluing[interfIdx].reversed2 = flipped;

        gsVector<T> n1_0, t1_0, n1_1, t1_1;
        gsVector<T> n2_0, t2_0, n2_1, t2_1;

        evalSTangAndNorm(m_target.patch(bi.first().patch), bi.first(), false, T(0), n1_0, t1_0);
        evalSTangAndNorm(m_target.patch(bi.first().patch), bi.first(), false, T(1), n1_1, t1_1);
        evalSTangAndNorm(m_target.patch(bi.second().patch), bi.second(), flipped, T(0), n2_0, t2_0);
        evalSTangAndNorm(m_target.patch(bi.second().patch), bi.second(), flipped, T(1), n2_1, t2_1);

        T a0_1 = std::max(areaMag(n1_0, t1_0), T(1e-12));
        T a1_1 = std::max(areaMag(n1_1, t1_1), T(1e-12));
        T a0_2 = std::max(areaMag(n2_0, t2_0), T(1e-12));
        T a1_2 = std::max(areaMag(n2_1, t2_1), T(1e-12));

        // Appendix A normalization: min_{xi in [0,1]} min(alpha_1(xi), alpha_2(xi)) = 1
        const T minA = std::min(std::min(a0_1, a1_1), std::min(a0_2, a1_2));
        if (minA > T(1e-12)) {
            a0_1 /= minA; a1_1 /= minA;
            a0_2 /= minA; a1_2 /= minA;
        }

        m_gluing[interfIdx].a0_1 = a0_1;
        m_gluing[interfIdx].a1_1 = a1_1;
        m_gluing[interfIdx].a0_2 = a0_2;
        m_gluing[interfIdx].a1_2 = a1_2;

        m_gluing[interfIdx].b0_1 = n1_0.dot(t1_0) / std::max(t1_0.squaredNorm(), T(1e-15));
        m_gluing[interfIdx].b1_1 = n1_1.dot(t1_1) / std::max(t1_1.squaredNorm(), T(1e-15));
        m_gluing[interfIdx].b0_2 = n2_0.dot(t2_0) / std::max(t2_0.squaredNorm(), T(1e-15));
        m_gluing[interfIdx].b1_2 = n2_1.dot(t2_1) / std::max(t2_1.squaredNorm(), T(1e-15));

        T u0, v0, ts0, u1, v1, ts1;
        mapInterfaceParam(bi.first(), false, T(0), u0, v0, ts0);
        mapInterfaceParam(bi.first(), false, T(1), u1, v1, ts1);
        gsMatrix<T> pt0(2, 1), pt1(2, 1);
        pt0 << u0, v0;
        pt1 << u1, v1;
        m_gluing[interfIdx].P0 = m_target.patch(bi.first().patch).eval(pt0);
        m_gluing[interfIdx].P1 = m_target.patch(bi.first().patch).eval(pt1);
    }
}

template <class T>
void gsASG1LocalFitter<T>::step1_solveVertexProblems()
{
    auto getGlobalCPIndex = [&](boxCorner corner, int iu, int iv) -> int {
        bool u_par = (corner.m_index - 1) & 1;
        bool v_par = ((corner.m_index - 1) >> 1) & 1;
        int i = (u_par == 0) ? iu : (m_n - 1 - iu);
        int j = (v_par == 0) ? iv : (m_n - 1 - iv);
        return i + j * m_n;
    };

    const auto& vertices = m_target.topology().vertices();
    int nVerts = vertices.size();

    struct VertexResult {
        std::vector<OrderedVertexPatch> ord;
        gsVector<T> x_v;
        gsMatrix<T> sol;
        bool valid;
    };
    std::vector<VertexResult> vResults(nVerts);

    #pragma omp parallel for schedule(dynamic)
    for (int vId = 0; vId < nVerts; ++vId) {
        auto ord = orderVertexPatches(vertices[vId]);
        int val = ord.size();
        vResults[vId].ord = ord;
        vResults[vId].valid = false;

        int p0 = ord[0].pc.patch;
        bool u_par0 = (ord[0].pc.m_index - 1) & 1;
        bool v_par0 = ((ord[0].pc.m_index - 1) >> 1) & 1;
        gsMatrix<T> ptCorner(2, 1);
        ptCorner << (u_par0 ? T(1) : T(0)), (v_par0 ? T(1) : T(0));
        vResults[vId].x_v = m_target.patch(p0).eval(ptCorner);

        if (val == 1) {
            vResults[vId].valid = true;
            continue;
        }

        int nVars = val * 5;
        gsMatrix<T> H(nVars, nVars); H.setZero();
        gsMatrix<T> F(nVars, m_geoDim); F.setZero();

        for (int j = 0; j < val; ++j) {
            int pId = ord[j].pc.patch;
            boxCorner c = ord[j].pc.corner();
            bool u_par = (c.m_index - 1) & 1;
            bool v_par = ((c.m_index - 1) >> 1) & 1;
            gsMatrix<T> ptC(2, 1);
            ptC << (u_par ? T(1) : T(0)), (v_par ? T(1) : T(0));
            gsMatrix<T> dS = m_target.patch(pId).deriv(ptC);
            gsMatrix<T> d2S = m_target.patch(pId).deriv2(ptC);

            gsVector<T> S_du(m_geoDim), S_dv(m_geoDim), S_duu(m_geoDim), S_dvv(m_geoDim), S_duv(m_geoDim);
            for (short_t d = 0; d < m_geoDim; ++d) {
                S_du(d)  = dS(2 * d + 0, 0);
                S_dv(d)  = dS(2 * d + 1, 0);
                S_duu(d) = d2S(3 * d + 0, 0);
                S_dvv(d) = d2S(3 * d + 1, 0);
                S_duv(d) = d2S(3 * d + 2, 0);
            }

            gsVector<T> S_d1  = (ord[j].dir_prev == 0) ? ord[j].sign_prev * S_du : ord[j].sign_prev * S_dv;
            gsVector<T> S_d2  = (ord[j].dir_next == 0) ? ord[j].sign_next * S_du : ord[j].sign_next * S_dv;
            gsVector<T> S_d11 = (ord[j].dir_prev == 0) ? S_duu : S_dvv;
            gsVector<T> S_d22 = (ord[j].dir_next == 0) ? S_duu : S_dvv;
            gsVector<T> S_d12 = ord[j].sign_prev * ord[j].sign_next * S_duv;

            int base = j * 5;
            T s1 = T(m_p) * T(m_k + 1);
            T s2 = T(m_p) * T(m_p - 1) * T(m_k + 1) * T(m_k + 1);
            T s1Sq = s1 * s1;
            T c1 = m_sigma * s1Sq;
            T c2 = m_sigma * m_sigma * (s2 * s2);
            T c12 = m_sigma * m_sigma * (s1Sq * s1Sq);

            gsVector<T> t1 = S_d1 / s1;
            gsVector<T> t2 = S_d2 / s1;
            gsVector<T> t11 = S_d11 / s2;
            gsVector<T> t22 = S_d22 / s2;
            gsVector<T> t12 = S_d12 / s1Sq;

            H(base + 0, base + 0) += T(2) * c1;
            F.row(base + 0) += T(2) * c1 * t1.transpose();

            H(base + 2, base + 2) += T(2) * c1;
            F.row(base + 2) += T(2) * c1 * t2.transpose();

            H(base + 0, base + 0) += T(2) * c2 * T(4);
            H(base + 1, base + 1) += T(2) * c2 * T(1);
            H(base + 0, base + 1) -= T(2) * c2 * T(2);
            H(base + 1, base + 0) -= T(2) * c2 * T(2);
            F.row(base + 0) -= T(2) * c2 * T(2) * t11.transpose();
            F.row(base + 1) += T(2) * c2 * T(1) * t11.transpose();

            H(base + 2, base + 2) += T(2) * c2 * T(4);
            H(base + 3, base + 3) += T(2) * c2 * T(1);
            H(base + 2, base + 3) -= T(2) * c2 * T(2);
            H(base + 3, base + 2) -= T(2) * c2 * T(2);
            F.row(base + 2) -= T(2) * c2 * T(2) * t22.transpose();
            F.row(base + 3) += T(2) * c2 * T(1) * t22.transpose();

            H(base + 0, base + 0) += T(2) * c12;
            H(base + 2, base + 2) += T(2) * c12;
            H(base + 4, base + 4) += T(2) * c12;
            H(base + 0, base + 2) += T(2) * c12;
            H(base + 2, base + 0) += T(2) * c12;
            H(base + 0, base + 4) -= T(2) * c12;
            H(base + 4, base + 0) -= T(2) * c12;
            H(base + 2, base + 4) -= T(2) * c12;
            H(base + 4, base + 2) -= T(2) * c12;
            F.row(base + 0) -= T(2) * c12 * t12.transpose();
            F.row(base + 2) -= T(2) * c12 * t12.transpose();
            F.row(base + 4) += T(2) * c12 * t12.transpose();
        }

        int numInterf = val - (ord[0].is_boundary_prev ? 1 : 0);
        std::vector<std::vector<std::pair<int, T>>> C_list;

        for (int j = 0; j < numInterf; ++j) {
            int j_next = (j + 1) % val;
            int bj = j * 5;
            int bnext = j_next * 5;

            // 1. Tangent 1st deriv matching
            C_list.push_back({ {bj + 2, T(1)}, {bnext + 0, T(-1)} });

            // 2. Tangent 2nd deriv matching
            C_list.push_back({ {bj + 3, T(1)}, {bnext + 1, T(-1)} });

            // 3. Transversal derivative matching
            T aj = T(1), anext = T(1), bj_coeff = T(0), bnext_coeff = T(0);
            for (const auto& gd : m_gluing) {
                if ((gd.ps1.patch == ord[j].pc.patch && gd.ps1.side().m_index == ord[j].side_next.side().m_index) ||
                    (gd.ps2.patch == ord[j].pc.patch && gd.ps2.side().m_index == ord[j].side_next.side().m_index))
                {
                    T d0 = (vResults[vId].x_v - gd.P0).norm();
                    T d1 = (vResults[vId].x_v - gd.P1).norm();
                    int xi_v = (d0 < d1) ? 0 : 1;
                    T sign_factor = (xi_v == 0) ? T(1) : T(-1);

                    if (gd.ps1.patch == ord[j].pc.patch) {
                        aj = (xi_v == 0) ? gd.a0_1 : gd.a1_1;
                        bj_coeff = sign_factor * ((xi_v == 0) ? gd.b0_1 : gd.b1_1);
                        anext = (xi_v == 0) ? gd.a0_2 : gd.a1_2;
                        bnext_coeff = sign_factor * ((xi_v == 0) ? gd.b0_2 : gd.b1_2);
                    } else {
                        aj = (xi_v == 0) ? gd.a0_2 : gd.a1_2;
                        bj_coeff = sign_factor * ((xi_v == 0) ? gd.b0_2 : gd.b1_2);
                        anext = (xi_v == 0) ? gd.a0_1 : gd.a1_1;
                        bnext_coeff = sign_factor * ((xi_v == 0) ? gd.b0_1 : gd.b1_1);
                    }
                    break;
                }
            }

            std::vector<std::pair<int, T>> r3;
            r3.push_back({bj + 0, T(1) / aj});
            r3.push_back({bj + 2, -bj_coeff / aj});
            r3.push_back({bnext + 2, T(1) / anext});
            r3.push_back({bnext + 0, -bnext_coeff / anext});
            C_list.push_back(r3);
        }

        int numC = C_list.size();
        gsEigen::Matrix<T, gsEigen::Dynamic, gsEigen::Dynamic> C_v(numC, nVars);
        C_v.setZero();
        for (int c = 0; c < numC; ++c) {
            for (const auto& item : C_list[c])
                C_v(c, item.first) += item.second;
            T rowNorm = C_v.row(c).norm();
            if (rowNorm > T(1e-14))
                C_v.row(c) /= rowNorm;
        }

        // Closed interior vertex loops (and near-degenerate boundary fans) induce
        // exact rank deficiency (2 redundant rows in 2D); solve the SPD QP on ker(C_v).
        gsEigen::JacobiSVD<gsEigen::Matrix<T, gsEigen::Dynamic, gsEigen::Dynamic>> svdC(C_v, gsEigen::ComputeFullV);
        const auto& sv = svdC.singularValues();
        int rankC = 0;
        if (sv.size() > 0 && sv(0) > T(0)) {
            const T svTol = std::max(T(1e-10) * sv(0), T(1e-12));
            while (rankC < sv.size() && sv(rankC) > svTol)
                ++rankC;
        }
        int nullDim = nVars - rankC;
        if (nullDim > 0) {
            gsMatrix<T> Z(svdC.matrixV().rightCols(nullDim));
            gsMatrix<T> Hred = Z.transpose() * H * Z;
            gsEigen::FullPivLU<gsMatrix<T>> lu(Hred);
            if (lu.isInvertible()) {
                vResults[vId].sol = Z * lu.solve(Z.transpose() * F);
                vResults[vId].valid = true;
            }
        }
    }

    // Apply solved vertex positions and derivatives to m_patchCP
    for (int vId = 0; vId < nVerts; ++vId) {
        const auto& ord = vResults[vId].ord;
        int val = ord.size();
        const gsVector<T>& x_v = vResults[vId].x_v;

        for (int j = 0; j < val; ++j) {
            int pId = ord[j].pc.patch;
            int idx00 = getGlobalCPIndex(ord[j].pc.corner(), 0, 0);
            m_patchCP[pId].row(idx00) = x_v.transpose();
            m_patchCP_fixed[pId][idx00] = true;
        }

        if (!vResults[vId].valid || val == 1) continue;

        const gsMatrix<T>& sol = vResults[vId].sol;

        for (int j = 0; j < val; ++j) {
            int pId = ord[j].pc.patch;
            boxCorner c = ord[j].pc.corner();
            int base = j * 5;

            auto setCP = [&](int i_prev, int i_next, int solRow) {
                int iu = (ord[j].dir_prev == 0) ? i_prev : i_next;
                int iv = (ord[j].dir_prev == 0) ? i_next : i_prev;
                int gIdx = getGlobalCPIndex(c, iu, iv);
                m_patchCP[pId].row(gIdx) = x_v.transpose() + sol.row(solRow);
                m_patchCP_fixed[pId][gIdx] = true;
            };

            setCP(1, 0, base + 0);
            setCP(2, 0, base + 1);
            setCP(0, 1, base + 2);
            setCP(0, 2, base + 3);
            setCP(1, 1, base + 4);
        }
    }
}

template <class T>
void gsASG1LocalFitter<T>::step2_solveInterfaceProblems()
{
    index_t nI = m_target.nInterfaces();

    gsBSplineBasis<T> basis(m_kv);
    gsGaussRule<T> quad1D(m_p + 3);
    gsMatrix<T> qNodes;
    gsVector<T> qWeights;
    quad1D.mapToAll(m_kv.breaks(), qNodes, qWeights);
    int nQ = qNodes.cols();

    gsSparseMatrix<T> B_1D = collocationMatrix(basis, qNodes, 0);
    gsSparseMatrix<T> D_1D = collocationMatrix(basis, qNodes, 1);
    gsMatrix<T> W_1D = qWeights.asDiagonal();
    gsMatrix<T> M1 = B_1D.transpose() * W_1D * B_1D;
    gsMatrix<T> K1 = D_1D.transpose() * W_1D * D_1D;

    const T s1 = T(m_p) * T(m_k + 1);
    const T s1Sq = s1 * s1;
    gsMatrix<T> H_r0_r0 = T(2) * (M1 + m_sigma * K1 + s1Sq * m_sigma * M1);
    gsMatrix<T> H_r1_r1 = T(2) * (s1Sq * m_sigma * M1);
    gsMatrix<T> H_r0_r1 = -T(2) * (s1Sq * m_sigma * M1);

    gsMatrix<T> greville(1, m_n);
    for (int i = 0; i < m_n; ++i) {
        T sum = 0;
        for (int d = 1; d <= m_p; ++d) sum += m_kv[i + d];
        greville(0, i) = sum / m_p;
    }
    gsMatrix<T> B_greville = gsMatrix<T>(collocationMatrix(basis, greville, 0));
    gsMatrix<T> D_greville = gsMatrix<T>(collocationMatrix(basis, greville, 1));

    gsKnotVector<T> kv0(T(0), T(1), m_k, m_p + 1, m_p - (m_r + 1));
    gsBSplineBasis<T> basis0(kv0);
    int n0 = basis0.size();

    gsKnotVector<T> kv1(T(0), T(1), m_k, (m_p - 1) + 1, (m_p - 1) - m_r);
    gsBSplineBasis<T> basis1(kv1);
    int n1 = basis1.size();

    gsMatrix<T> B0_greville = gsMatrix<T>(collocationMatrix(basis0, greville, 0));
    gsMatrix<T> B1_greville = gsMatrix<T>(collocationMatrix(basis1, greville, 0));

    auto getEdgeCPIndex = [&](const patchSide& ps, int layer, int tIdx) -> int {
        SideInfo si = getSideInfo(ps.side());
        int iu = (si.nDir == 0) ? (si.par ? (m_n - 1 - layer) : layer) : tIdx;
        int iv = (si.nDir == 1) ? (si.par ? (m_n - 1 - layer) : layer) : tIdx;
        return iu + iv * m_n;
    };

    for (int m = 0; m < nI; ++m) {
        const auto& gd = m_gluing[m];
        int p1 = gd.ps1.patch;
        int p2 = gd.ps2.patch;

        const T aScale = std::max(std::max(gd.a0_1, gd.a1_1), std::max(gd.a0_2, gd.a1_2));
        const T safeAScale = (aScale > T(1e-12)) ? aScale : T(1e-12);

        gsMatrix<T> ptsS1(2, nQ), ptsS2(2, nQ);
        for (int q = 0; q < nQ; ++q) {
            T xi = qNodes(0, q);
            T u1, v1, tSign1, u2, v2, tSign2;
            mapInterfaceParam(gd.ps1, gd.reversed1, xi, u1, v1, tSign1);
            mapInterfaceParam(gd.ps2, gd.reversed2, xi, u2, v2, tSign2);
            ptsS1(0, q) = u1; ptsS1(1, q) = v1;
            ptsS2(0, q) = u2; ptsS2(1, q) = v2;
        }

        gsMatrix<T> evalS1 = m_target.patch(p1).eval(ptsS1);
        gsMatrix<T> S_bdr(m_geoDim, nQ), S_tan(m_geoDim, nQ), S_norm(m_geoDim, nQ);
        for (int q = 0; q < nQ; ++q) {
            S_bdr.col(q) = evalS1.col(q);
            gsVector<T> nVec, tVec;
            evalSTangAndNorm(m_target.patch(p1), gd.ps1, gd.reversed1, qNodes(0, q), nVec, tVec);
            S_norm.col(q) = nVec;
            S_tan.col(q) = tVec;
        }

        std::vector<int> free_r0_1, free_r1_1, free_r0_2, free_r1_2;
        for (int j = 3; j <= m_n - 4; ++j) free_r0_1.push_back(j);
        for (int j = 2; j <= m_n - 3; ++j) free_r1_1.push_back(j);
        for (int j = 3; j <= m_n - 4; ++j) free_r0_2.push_back(j);
        for (int j = 2; j <= m_n - 3; ++j) free_r1_2.push_back(j);

        int nf_r0_1 = free_r0_1.size();
        int nf_r1_1 = free_r1_1.size();
        int nf_r0_2 = free_r0_2.size();
        int nf_r1_2 = free_r1_2.size();

        int off_r0_1 = 0;
        int off_r1_1 = off_r0_1 + nf_r0_1;
        int off_r0_2 = off_r1_1 + nf_r1_1;
        int off_r1_2 = off_r0_2 + nf_r0_2;
        int off_d0_aux = off_r1_2 + nf_r1_2;
        int off_d1_aux = off_d0_aux + n0;
        int nVars = off_d1_aux + n1;

        gsMatrix<T> H_int(nVars, nVars); H_int.setZero();
        gsMatrix<T> F_int(nVars, m_geoDim); F_int.setZero();

        gsMatrix<T> rhs_full_r0 = T(2) * (B_1D.transpose() * W_1D * S_bdr.transpose() +
                                          m_sigma * D_1D.transpose() * W_1D * S_tan.transpose() -
                                          s1 * m_sigma * B_1D.transpose() * W_1D * S_norm.transpose());
        gsMatrix<T> rhs_full_r1 = T(2) * (s1 * m_sigma * B_1D.transpose() * W_1D * S_norm.transpose());

        for (int j = 0; j < m_n; ++j) {
            int gIdx0 = getEdgeCPIndex(gd.ps1, 0, gd.reversed1 ? (m_n - 1 - j) : j);
            int gIdx1 = getEdgeCPIndex(gd.ps1, 1, gd.reversed1 ? (m_n - 1 - j) : j);
            if (m_patchCP_fixed[p1][gIdx0]) {
                rhs_full_r0 -= H_r0_r0.col(j) * m_patchCP[p1].row(gIdx0);
                rhs_full_r1 -= H_r0_r1.row(j).transpose() * m_patchCP[p1].row(gIdx0);
            }
            if (m_patchCP_fixed[p1][gIdx1]) {
                rhs_full_r0 -= H_r0_r1.col(j) * m_patchCP[p1].row(gIdx1);
                rhs_full_r1 -= H_r1_r1.col(j) * m_patchCP[p1].row(gIdx1);
            }
        }

        for (int a = 0; a < nf_r0_1; ++a) {
            int ja = free_r0_1[a];
            F_int.row(off_r0_1 + a) = rhs_full_r0.row(ja);
            for (int b = 0; b < nf_r0_1; ++b) H_int(off_r0_1 + a, off_r0_1 + b) = H_r0_r0(ja, free_r0_1[b]);
            for (int b = 0; b < nf_r1_1; ++b) {
                int jb = free_r1_1[b];
                H_int(off_r0_1 + a, off_r1_1 + b) = H_r0_r1(ja, jb);
                H_int(off_r1_1 + b, off_r0_1 + a) = H_r0_r1(ja, jb);
            }
        }
        for (int a = 0; a < nf_r1_1; ++a) {
            int ja = free_r1_1[a];
            F_int.row(off_r1_1 + a) = rhs_full_r1.row(ja);
            for (int b = 0; b < nf_r1_1; ++b) H_int(off_r1_1 + a, off_r1_1 + b) = H_r1_r1(ja, free_r1_1[b]);
        }

        gsMatrix<T> S_bdr2(m_geoDim, nQ), S_tan2(m_geoDim, nQ), S_norm2(m_geoDim, nQ);
        gsMatrix<T> evalS2 = m_target.patch(p2).eval(ptsS2);
        for (int q = 0; q < nQ; ++q) {
            S_bdr2.col(q) = evalS2.col(q);
            gsVector<T> nVec2, tVec2;
            evalSTangAndNorm(m_target.patch(p2), gd.ps2, gd.reversed2, qNodes(0, q), nVec2, tVec2);
            S_norm2.col(q) = nVec2;
            S_tan2.col(q) = tVec2;
        }

        gsMatrix<T> rhs_full_r0_2 = T(2) * (B_1D.transpose() * W_1D * S_bdr2.transpose() +
                                            m_sigma * D_1D.transpose() * W_1D * S_tan2.transpose() -
                                            s1 * m_sigma * B_1D.transpose() * W_1D * S_norm2.transpose());
        gsMatrix<T> rhs_full_r1_2 = T(2) * (s1 * m_sigma * B_1D.transpose() * W_1D * S_norm2.transpose());

        for (int j = 0; j < m_n; ++j) {
            int gIdx0 = getEdgeCPIndex(gd.ps2, 0, gd.reversed2 ? (m_n - 1 - j) : j);
            int gIdx1 = getEdgeCPIndex(gd.ps2, 1, gd.reversed2 ? (m_n - 1 - j) : j);
            if (m_patchCP_fixed[p2][gIdx0]) {
                rhs_full_r0_2 -= H_r0_r0.col(j) * m_patchCP[p2].row(gIdx0);
                rhs_full_r1_2 -= H_r0_r1.row(j).transpose() * m_patchCP[p2].row(gIdx0);
            }
            if (m_patchCP_fixed[p2][gIdx1]) {
                rhs_full_r0_2 -= H_r0_r1.col(j) * m_patchCP[p2].row(gIdx1);
                rhs_full_r1_2 -= H_r1_r1.col(j) * m_patchCP[p2].row(gIdx1);
            }
        }

        for (int a = 0; a < nf_r0_2; ++a) {
            int ja = free_r0_2[a];
            F_int.row(off_r0_2 + a) = rhs_full_r0_2.row(ja);
            for (int b = 0; b < nf_r0_2; ++b) H_int(off_r0_2 + a, off_r0_2 + b) = H_r0_r0(ja, free_r0_2[b]);
            for (int b = 0; b < nf_r1_2; ++b) {
                int jb = free_r1_2[b];
                H_int(off_r0_2 + a, off_r1_2 + b) = H_r0_r1(ja, jb);
                H_int(off_r1_2 + b, off_r0_2 + a) = H_r0_r1(ja, jb);
            }
        }
        for (int a = 0; a < nf_r1_2; ++a) {
            int ja = free_r1_2[a];
            F_int.row(off_r1_2 + a) = rhs_full_r1_2.row(ja);
            for (int b = 0; b < nf_r1_2; ++b) H_int(off_r1_2 + a, off_r1_2 + b) = H_r1_r1(ja, free_r1_2[b]);
        }

        int numC = 2 * m_n + nf_r0_2 + nf_r1_2;
        gsMatrix<T> C_int(numC, nVars); C_int.setZero();
        gsMatrix<T> G_int(numC, m_geoDim); G_int.setZero();

        int c_row = 0;
        // Constraint A: Patch 1 boundary match (j = 0..m_n-1)
        for (int j = 0; j < m_n; ++j) {
            gsMatrix<T> fixed_v = gsMatrix<T>::Zero(1, m_geoDim);
            for (int i = 0; i < m_n; ++i) {
                int gIdx = getEdgeCPIndex(gd.ps1, 0, gd.reversed1 ? (m_n - 1 - i) : i);
                if (m_patchCP_fixed[p1][gIdx]) {
                    fixed_v += B_greville(j, i) * m_patchCP[p1].row(gIdx);
                } else {
                    auto it = std::find(free_r0_1.begin(), free_r0_1.end(), i);
                    if (it != free_r0_1.end()) {
                        int pos = off_r0_1 + std::distance(free_r0_1.begin(), it);
                        C_int(c_row, pos) = B_greville(j, i);
                    }
                }
            }
            for (int l = 0; l < n0; ++l) C_int(c_row, off_d0_aux + l) = -B0_greville(j, l);
            G_int.row(c_row) = -fixed_v;
            c_row++;
        }

        // Constraint B: Patch 1 transversal match (j = 0..m_n-1)
        for (int j = 0; j < m_n; ++j) {
            T xi = greville(0, j);
            T a1 = gd.alpha1(xi);
            T b1 = gd.beta1(xi);
            gsMatrix<T> fixed_v = gsMatrix<T>::Zero(1, m_geoDim);

            for (int i = 0; i < m_n; ++i) {
                T coeff_0 = -b1 * D_greville(j, i) - s1 * B_greville(j, i);
                T coeff_1 = s1 * B_greville(j, i);

                int gIdx0 = getEdgeCPIndex(gd.ps1, 0, gd.reversed1 ? (m_n - 1 - i) : i);
                int gIdx1 = getEdgeCPIndex(gd.ps1, 1, gd.reversed1 ? (m_n - 1 - i) : i);

                if (m_patchCP_fixed[p1][gIdx0]) {
                    fixed_v += coeff_0 * m_patchCP[p1].row(gIdx0);
                } else {
                    auto it = std::find(free_r0_1.begin(), free_r0_1.end(), i);
                    if (it != free_r0_1.end()) {
                        int pos = off_r0_1 + std::distance(free_r0_1.begin(), it);
                        C_int(c_row, pos) = coeff_0;
                    }
                }

                if (m_patchCP_fixed[p1][gIdx1]) {
                    fixed_v += coeff_1 * m_patchCP[p1].row(gIdx1);
                } else {
                    auto it = std::find(free_r1_1.begin(), free_r1_1.end(), i);
                    if (it != free_r1_1.end()) {
                        int pos = off_r1_1 + std::distance(free_r1_1.begin(), it);
                        C_int(c_row, pos) = coeff_1;
                    }
                }
            }

            for (int l = 0; l < n1; ++l) C_int(c_row, off_d1_aux + l) = -(a1 / safeAScale) * B1_greville(j, l);
            G_int.row(c_row) = -fixed_v;
            c_row++;
        }

        // Constraint C: Patch 2 boundary match (j = 3..m_n-4)
        for (int j = 3; j <= m_n - 4; ++j) {
            gsMatrix<T> fixed_v = gsMatrix<T>::Zero(1, m_geoDim);
            for (int i = 0; i < m_n; ++i) {
                int gIdx = getEdgeCPIndex(gd.ps2, 0, gd.reversed2 ? (m_n - 1 - i) : i);
                if (m_patchCP_fixed[p2][gIdx]) {
                    fixed_v += B_greville(j, i) * m_patchCP[p2].row(gIdx);
                } else {
                    auto it = std::find(free_r0_2.begin(), free_r0_2.end(), i);
                    if (it != free_r0_2.end()) {
                        int pos = off_r0_2 + std::distance(free_r0_2.begin(), it);
                        C_int(c_row, pos) = B_greville(j, i);
                    }
                }
            }
            for (int l = 0; l < n0; ++l) C_int(c_row, off_d0_aux + l) = -B0_greville(j, l);
            G_int.row(c_row) = -fixed_v;
            c_row++;
        }

        // Constraint D: Patch 2 transversal match (j = 2..m_n-3)
        for (int j = 2; j <= m_n - 3; ++j) {
            T xi = greville(0, j);
            T a2 = gd.alpha2(xi);
            T b2 = gd.beta2(xi);
            gsMatrix<T> fixed_v = gsMatrix<T>::Zero(1, m_geoDim);

            for (int i = 0; i < m_n; ++i) {
                T coeff_0 = s1 * B_greville(j, i) + b2 * D_greville(j, i);
                T coeff_1 = -s1 * B_greville(j, i);

                int gIdx0 = getEdgeCPIndex(gd.ps2, 0, gd.reversed2 ? (m_n - 1 - i) : i);
                int gIdx1 = getEdgeCPIndex(gd.ps2, 1, gd.reversed2 ? (m_n - 1 - i) : i);

                if (m_patchCP_fixed[p2][gIdx0]) {
                    fixed_v += coeff_0 * m_patchCP[p2].row(gIdx0);
                } else {
                    auto it = std::find(free_r0_2.begin(), free_r0_2.end(), i);
                    if (it != free_r0_2.end()) {
                        int pos = off_r0_2 + std::distance(free_r0_2.begin(), it);
                        C_int(c_row, pos) = coeff_0;
                    }
                }

                if (m_patchCP_fixed[p2][gIdx1]) {
                    fixed_v += coeff_1 * m_patchCP[p2].row(gIdx1);
                } else {
                    auto it = std::find(free_r1_2.begin(), free_r1_2.end(), i);
                    if (it != free_r1_2.end()) {
                        int pos = off_r1_2 + std::distance(free_r1_2.begin(), it);
                        C_int(c_row, pos) = coeff_1;
                    }
                }
            }

            for (int l = 0; l < n1; ++l) C_int(c_row, off_d1_aux + l) = -(a2 / safeAScale) * B1_greville(j, l);
            G_int.row(c_row) = -fixed_v;
            c_row++;
        }

        gsMatrix<T> KKT(nVars + numC, nVars + numC); KKT.setZero();
        KKT.block(0, 0, nVars, nVars) = H_int;
        KKT.block(0, nVars, nVars, numC) = C_int.transpose();
        KKT.block(nVars, 0, numC, nVars) = C_int;

        gsMatrix<T> rhs_kkt(nVars + numC, m_geoDim);
        rhs_kkt.topRows(nVars) = F_int;
        rhs_kkt.bottomRows(numC) = G_int;

        gsEigen::FullPivLU<gsMatrix<T>> lu(KKT);
        if (!lu.isInvertible()) {
            std::cerr << "gsASG1LocalFitter: Warning: Interface " << m << " KKT not invertible!\n";
            continue;
        }
        gsMatrix<T> sol = lu.solve(rhs_kkt);

        for (int a = 0; a < nf_r0_1; ++a) {
            int j = free_r0_1[a];
            int gIdx = getEdgeCPIndex(gd.ps1, 0, gd.reversed1 ? (m_n - 1 - j) : j);
            m_patchCP[p1].row(gIdx) = sol.row(off_r0_1 + a);
            m_patchCP_fixed[p1][gIdx] = true;
        }
        for (int a = 0; a < nf_r1_1; ++a) {
            int j = free_r1_1[a];
            int gIdx = getEdgeCPIndex(gd.ps1, 1, gd.reversed1 ? (m_n - 1 - j) : j);
            m_patchCP[p1].row(gIdx) = sol.row(off_r1_1 + a);
            m_patchCP_fixed[p1][gIdx] = true;
        }
        for (int a = 0; a < nf_r0_2; ++a) {
            int j = free_r0_2[a];
            int gIdx = getEdgeCPIndex(gd.ps2, 0, gd.reversed2 ? (m_n - 1 - j) : j);
            m_patchCP[p2].row(gIdx) = sol.row(off_r0_2 + a);
            m_patchCP_fixed[p2][gIdx] = true;
        }
        for (int a = 0; a < nf_r1_2; ++a) {
            int j = free_r1_2[a];
            int gIdx = getEdgeCPIndex(gd.ps2, 1, gd.reversed2 ? (m_n - 1 - j) : j);
            m_patchCP[p2].row(gIdx) = sol.row(off_r1_2 + a);
            m_patchCP_fixed[p2][gIdx] = true;
        }
    }
}

template <class T>
void gsASG1LocalFitter<T>::step2b_solveBoundaryEdges()
{
    const auto& bdrSides = m_target.boundaries();
    int nB = bdrSides.size();

    gsBSplineBasis<T> basis(m_kv);
    gsGaussRule<T> quad1D(m_p + 3);
    gsMatrix<T> qNodes;
    gsVector<T> qWeights;
    quad1D.mapToAll(m_kv.breaks(), qNodes, qWeights);
    int nQ = qNodes.cols();

    gsSparseMatrix<T> B_1D = collocationMatrix(basis, qNodes, 0);
    gsSparseMatrix<T> D_1D = collocationMatrix(basis, qNodes, 1);
    gsMatrix<T> W_1D = qWeights.asDiagonal();
    gsMatrix<T> M1 = B_1D.transpose() * W_1D * B_1D;
    gsMatrix<T> K1 = D_1D.transpose() * W_1D * D_1D;

    const T s1 = T(m_p) * T(m_k + 1);
    const T s1Sq = s1 * s1;
    gsMatrix<T> H_r0_r0 = T(2) * (M1 + m_sigma * K1 + s1Sq * m_sigma * M1);
    gsMatrix<T> H_r1_r1 = T(2) * (s1Sq * m_sigma * M1);
    gsMatrix<T> H_r0_r1 = -T(2) * (s1Sq * m_sigma * M1);

    auto getEdgeCPIndex = [&](const patchSide& ps, int layer, int tIdx) -> int {
        SideInfo si = getSideInfo(ps.side());
        int iu = (si.nDir == 0) ? (si.par ? (m_n - 1 - layer) : layer) : tIdx;
        int iv = (si.nDir == 1) ? (si.par ? (m_n - 1 - layer) : layer) : tIdx;
        return iu + iv * m_n;
    };

    for (int b = 0; b < nB; ++b) {
        const patchSide& ps = bdrSides[b];
        int pId = ps.patch;

        gsMatrix<T> ptsS(2, nQ);
        for (int q = 0; q < nQ; ++q) {
            T xi = qNodes(0, q);
            T u, v, ts;
            mapInterfaceParam(ps, false, xi, u, v, ts);
            ptsS(0, q) = u; ptsS(1, q) = v;
        }
        gsMatrix<T> evalS = m_target.patch(pId).eval(ptsS);
        gsMatrix<T> S_bdr(m_geoDim, nQ), S_tan(m_geoDim, nQ), S_norm(m_geoDim, nQ);
        for (int q = 0; q < nQ; ++q) {
            S_bdr.col(q) = evalS.col(q);
            gsVector<T> nVec, tVec;
            evalSTangAndNorm(m_target.patch(pId), ps, false, qNodes(0, q), nVec, tVec);
            S_norm.col(q) = nVec;
            S_tan.col(q) = tVec;
        }

        std::vector<int> free_r0, free_r1;
        for (int j = 0; j < m_n; ++j) {
            int gIdx0 = getEdgeCPIndex(ps, 0, j);
            int gIdx1 = getEdgeCPIndex(ps, 1, j);
            if (!m_patchCP_fixed[pId][gIdx0]) free_r0.push_back(j);
            if (!m_patchCP_fixed[pId][gIdx1]) free_r1.push_back(j);
        }

        int nf0 = free_r0.size();
        int nf1 = free_r1.size();
        int nFree = nf0 + nf1;
        if (nFree == 0) continue;

        gsMatrix<T> rhs_full_r0 = T(2) * (B_1D.transpose() * W_1D * S_bdr.transpose() +
                                          m_sigma * D_1D.transpose() * W_1D * S_tan.transpose() -
                                          s1 * m_sigma * B_1D.transpose() * W_1D * S_norm.transpose());
        gsMatrix<T> rhs_full_r1 = T(2) * (s1 * m_sigma * B_1D.transpose() * W_1D * S_norm.transpose());

        for (int j = 0; j < m_n; ++j) {
            int gIdx0 = getEdgeCPIndex(ps, 0, j);
            int gIdx1 = getEdgeCPIndex(ps, 1, j);
            if (m_patchCP_fixed[pId][gIdx0]) {
                rhs_full_r0 -= H_r0_r0.col(j) * m_patchCP[pId].row(gIdx0);
                rhs_full_r1 -= H_r0_r1.row(j).transpose() * m_patchCP[pId].row(gIdx0);
            }
            if (m_patchCP_fixed[pId][gIdx1]) {
                rhs_full_r0 -= H_r0_r1.col(j) * m_patchCP[pId].row(gIdx1);
                rhs_full_r1 -= H_r1_r1.col(j) * m_patchCP[pId].row(gIdx1);
            }
        }

        gsMatrix<T> H_bdr(nFree, nFree); H_bdr.setZero();
        gsMatrix<T> rhs(nFree, m_geoDim);
        for (int a = 0; a < nf0; ++a) {
            int ja = free_r0[a];
            rhs.row(a) = rhs_full_r0.row(ja);
            for (int c = 0; c < nf0; ++c) H_bdr(a, c) = H_r0_r0(ja, free_r0[c]);
            for (int c = 0; c < nf1; ++c) {
                H_bdr(a, nf0 + c) = H_r0_r1(ja, free_r1[c]);
                H_bdr(nf0 + c, a) = H_r0_r1(ja, free_r1[c]);
            }
        }
        for (int a = 0; a < nf1; ++a) {
            int ja = free_r1[a];
            rhs.row(nf0 + a) = rhs_full_r1.row(ja);
            for (int c = 0; c < nf1; ++c) H_bdr(nf0 + a, nf0 + c) = H_r1_r1(ja, free_r1[c]);
        }

        gsEigen::FullPivLU<gsMatrix<T>> lu(H_bdr);
        if (lu.isInvertible()) {
            gsMatrix<T> sol = lu.solve(rhs);
            for (int a = 0; a < nf0; ++a) {
                int j = free_r0[a];
                int gIdx = getEdgeCPIndex(ps, 0, j);
                m_patchCP[pId].row(gIdx) = sol.row(a);
                m_patchCP_fixed[pId][gIdx] = true;
            }
            for (int a = 0; a < nf1; ++a) {
                int j = free_r1[a];
                int gIdx = getEdgeCPIndex(ps, 1, j);
                m_patchCP[pId].row(gIdx) = sol.row(nf0 + a);
                m_patchCP_fixed[pId][gIdx] = true;
            }
        }
    }
}

template <class T>
void gsASG1LocalFitter<T>::step3_solvePatchInteriors()
{
    index_t nP = m_target.nPatches();

    gsBSplineBasis<T> basis(m_kv);
    gsGaussRule<T> quad1D(m_p + 3);
    gsMatrix<T> qNodes;
    gsVector<T> qWeights;
    quad1D.mapToAll(m_kv.breaks(), qNodes, qWeights);
    int nQ = qNodes.cols();

    gsSparseMatrix<T> B_1D = collocationMatrix(basis, qNodes, 0);
    gsSparseMatrix<T> D_1D = collocationMatrix(basis, qNodes, 1);
    gsMatrix<T> W_1D = qWeights.asDiagonal();
    gsMatrix<T> M1 = B_1D.transpose() * W_1D * B_1D;
    gsMatrix<T> K1 = D_1D.transpose() * W_1D * D_1D;

    gsMatrix<T> A_2D(m_n2, m_n2);
    for (int j1 = 0; j1 < m_n; ++j1)
        for (int i1 = 0; i1 < m_n; ++i1) {
            int I = i1 + j1 * m_n;
            for (int j2 = 0; j2 < m_n; ++j2)
                for (int i2 = 0; i2 < m_n; ++i2) {
                    int J = i2 + j2 * m_n;
                    A_2D(I, J) = M1(i1, i2) * M1(j1, j2) +
                                 m_sigma * (K1(i1, i2) * M1(j1, j2) + M1(i1, i2) * K1(j1, j2));
                }
        }

    int nQ2 = nQ * nQ;
    gsMatrix<T> qNodes2D(2, nQ2);
    gsVector<T> qWeights2D(nQ2);
    for (int j = 0; j < nQ; ++j)
        for (int i = 0; i < nQ; ++i) {
            int idx = i + j * nQ;
            qNodes2D(0, idx) = qNodes(0, i);
            qNodes2D(1, idx) = qNodes(0, j);
            qWeights2D(idx) = qWeights(i) * qWeights(j);
        }

    #pragma omp parallel for schedule(dynamic)
    for (int pId = 0; pId < nP; ++pId) {
        gsMatrix<T> evalS = m_target.patch(pId).eval(qNodes2D);
        gsMatrix<T> derivS = m_target.patch(pId).deriv(qNodes2D);

        gsMatrix<T> B_rhs(m_n2, m_geoDim);
        B_rhs.setZero();

        for (int j = 0; j < m_n; ++j) {
            for (int i = 0; i < m_n; ++i) {
                int I = i + j * m_n;
                gsMatrix<T> val = gsMatrix<T>::Zero(1, m_geoDim);
                for (int qj = 0; qj < nQ; ++qj) {
                    T Bj_v = B_1D.coeff(qj, j);
                    T Dj_v = D_1D.coeff(qj, j);
                    for (int qi = 0; qi < nQ; ++qi) {
                        int qIdx = qi + qj * nQ;
                        T w = qWeights2D(qIdx);
                        T Bi_u = B_1D.coeff(qi, i);
                        T Di_u = D_1D.coeff(qi, i);

                        T B_val = Bi_u * Bj_v;
                        T Du_val = Di_u * Bj_v;
                        T Dv_val = Bi_u * Dj_v;

                        for (short_t d = 0; d < m_geoDim; ++d) {
                            T sc = evalS(d, qIdx);
                            T dsc_du = derivS(2 * d + 0, qIdx);
                            T dsc_dv = derivS(2 * d + 1, qIdx);
                            val(0, d) += w * (sc * B_val + m_sigma * (dsc_du * Du_val + dsc_dv * Dv_val));
                        }
                    }
                }
                B_rhs.row(I) = val;
            }
        }

        std::vector<int> free_idx, fixed_idx;
        for (int i = 0; i < m_n2; ++i) {
            if (m_patchCP_fixed[pId][i])
                fixed_idx.push_back(i);
            else
                free_idx.push_back(i);
        }

        int nFree = free_idx.size();
        int nFixed = fixed_idx.size();

        if (nFree > 0) {
            gsMatrix<T> A_ff(nFree, nFree);
            gsMatrix<T> A_fc(nFree, nFixed);
            gsMatrix<T> rhs_ff(nFree, m_geoDim);
            gsMatrix<T> d_fixed(nFixed, m_geoDim);

            for (int a = 0; a < nFixed; ++a)
                d_fixed.row(a) = m_patchCP[pId].row(fixed_idx[a]);

            for (int a = 0; a < nFree; ++a) {
                int ia = free_idx[a];
                rhs_ff.row(a) = B_rhs.row(ia);
                for (int b = 0; b < nFree; ++b) A_ff(a, b) = A_2D(ia, free_idx[b]);
                for (int b = 0; b < nFixed; ++b) A_fc(a, b) = A_2D(ia, fixed_idx[b]);
            }

            rhs_ff -= A_fc * d_fixed;

            gsEigen::FullPivLU<gsMatrix<T>> lu(A_ff);
            gsMatrix<T> sol = lu.solve(rhs_ff);

            for (int a = 0; a < nFree; ++a) {
                int ia = free_idx[a];
                m_patchCP[pId].row(ia) = sol.row(a);
                m_patchCP_fixed[pId][ia] = true;
            }
        }
    }
}

template <class T>
void gsASG1LocalFitter<T>::fit()
{
    if (m_mode == FittingMode::Global) {
        fitGlobal();
    } else {
        fitLocal();
    }
}

template <class T>
void gsASG1LocalFitter<T>::fitLocal()
{
    // Step 0: Precompute linear gluing functions
    step0_computeGluingFunctions();

    // Step 1: Vertex quadratic minimization
    step1_solveVertexProblems();

    // Step 2: Interface quadratic minimization
    step2_solveInterfaceProblems();

    // Step 2b: Boundary edge minimization
    step2b_solveBoundaryEdges();

    // Step 3: Patch interior minimization
    step3_solvePatchInteriors();

    // Build output gsMultiPatch
    index_t nP = m_target.nPatches();
    m_result.clear();
    for (int pId = 0; pId < nP; ++pId) {
        gsTensorBSplineBasis<2, T> tBasis(m_kv, m_kv);
        typename gsGeometry<T>::uPtr patchGeom(new gsTensorBSpline<2, T>(tBasis, m_patchCP[pId]));
        m_result.addPatch(give(patchGeom));
    }
    for (size_t b = 0; b < m_target.boundaries().size(); ++b)
        m_result.addBoundary(m_target.boundaries()[b]);
    for (auto it = m_target.iBegin(); it != m_target.iEnd(); ++it)
        m_result.addInterface(*it);

    m_result.computeTopology();
}

template <class T>
void gsASG1LocalFitter<T>::fitGlobal()
{
    // Step 0: Precompute linear gluing functions
    step0_computeGluingFunctions();

    index_t nP = m_target.nPatches();
    index_t nI = m_target.nInterfaces();

    // Patch basis: S^{p, r}_k
    gsBSplineBasis<T> bbasis_p(m_kv);
    gsTensorBSplineBasis<2, T> tbasis(m_kv, m_kv);
    index_t m_n = bbasis_p.size();
    index_t nCP = tbasis.size();

    // Auxiliary f0 basis: S^{p, r+1}_k
    gsKnotVector<T> kv_0(T(0), T(1), m_k, m_p + 1, m_p - (m_r + 1), m_p);
    gsBSplineBasis<T> bbasis_0(kv_0);
    index_t n0 = bbasis_0.size();

    // Auxiliary f1 basis: S^{p-1, r}_k
    gsKnotVector<T> kv_1(T(0), T(1), m_k, (m_p - 1) + 1, (m_p - 1) - m_r, m_p - 1);
    gsBSplineBasis<T> bbasis_1(kv_1);
    index_t n1 = bbasis_1.size();

    T sigma = m_sigma;

    // 1D Gauss quadrature
    gsGaussRule<T> grule(m_p + 3);
    gsMatrix<T> qNodes1D;
    gsVector<T> qWeights1D;
    grule.mapToAll(m_kv.breaks(), qNodes1D, qWeights1D);
    int nQ1 = qNodes1D.cols();

    gsMatrix<T> B_val1D(nQ1, m_n), D_val1D(nQ1, m_n);
    B_val1D.setZero(); D_val1D.setZero();
    for (int q = 0; q < nQ1; ++q) {
        gsMatrix<T> pt(1, 1);
        pt(0, 0) = qNodes1D(0, q);
        gsMatrix<T> val = bbasis_p.eval(pt);
        gsMatrix<T> der = bbasis_p.deriv(pt);
        gsMatrix<index_t> act = bbasis_p.active(pt);
        for (int i = 0; i < act.rows(); ++i) {
            B_val1D(q, act(i, 0)) = val(i, 0);
            D_val1D(q, act(i, 0)) = der(i, 0);
        }
    }

    gsMatrix<T> W1D = qWeights1D.asDiagonal();
    gsMatrix<T> M1D = B_val1D.transpose() * W1D * B_val1D;
    gsMatrix<T> K1D = D_val1D.transpose() * W1D * D_val1D;

    gsMatrix<T> Q(nCP, nCP); Q.setZero();
    for (int iv1 = 0; iv1 < m_n; ++iv1) {
        for (int iu1 = 0; iu1 < m_n; ++iu1) {
            int a = iu1 + iv1 * m_n;
            for (int iv2 = 0; iv2 < m_n; ++iv2) {
                for (int iu2 = 0; iu2 < m_n; ++iu2) {
                    int b = iu2 + iv2 * m_n;
                    T val_m = M1D(iu1, iu2) * M1D(iv1, iv2);
                    T val_ku = K1D(iu1, iu2) * M1D(iv1, iv2);
                    T val_kv = M1D(iu1, iu2) * K1D(iv1, iv2);
                    Q(a, b) = val_m + sigma * (val_ku + val_kv);
                }
            }
        }
    }

    int nQ2 = nQ1 * nQ1;
    gsMatrix<T> pts2(2, nQ2);
    gsVector<T> w2(nQ2);
    for (int iv = 0; iv < nQ1; ++iv) {
        for (int iu = 0; iu < nQ1; ++iu) {
            int q = iu + iv * nQ1;
            pts2(0, q) = qNodes1D(0, iu);
            pts2(1, q) = qNodes1D(0, iv);
            w2(q) = qWeights1D(iu) * qWeights1D(iv);
        }
    }

    index_t nVar = nP * nCP + nI * (n0 + n1);
    std::vector<gsMatrix<T>> Q_patches(nP, Q);
    gsMatrix<T> F_glob(nVar, m_geoDim);
    F_glob.setZero();

    for (index_t p_idx = 0; p_idx < nP; ++p_idx) {
        gsMatrix<T> valS = m_target.patch(p_idx).eval(pts2);
        gsMatrix<T> derS = m_target.patch(p_idx).deriv(pts2);
        int off = p_idx * nCP;

        for (int iv = 0; iv < nQ1; ++iv) {
            for (int iu = 0; iu < nQ1; ++iu) {
                int q = iu + iv * nQ1;
                T w = w2(q);

                for (int jv = 0; jv < m_n; ++jv) {
                    T bv = B_val1D(iv, jv);
                    T dv_b = D_val1D(iv, jv);
                    for (int ju = 0; ju < m_n; ++ju) {
                        T bu = B_val1D(iu, ju);
                        T du_b = D_val1D(iu, ju);
                        int a = off + (ju + jv * m_n);

                        T phi = bu * bv;
                        T dphi_u = du_b * bv;
                        T dphi_v = bu * dv_b;

                        for (short_t d = 0; d < m_geoDim; ++d) {
                            T sc = valS(d, q);
                            T du_sc = derS(2 * d + 0, q);
                            T dv_sc = derS(2 * d + 1, q);
                            F_glob(a, d) += w * (sc * phi + sigma * (du_sc * dphi_u + dv_sc * dphi_v));
                        }
                    }
                }
            }
        }
    }

    auto getCPIdx = [&](short_t nDir, int nIdx, int tIdx) -> int {
        return (nDir == 0) ? (nIdx + tIdx * m_n) : (tIdx + nIdx * m_n);
    };

    // Anchor exterior domain boundaries (L2 position only) so global AS-G1 preserves corners & curves.
    const T wBndGlob = T(500);
    gsMatrix<T> M1_bnd = wBndGlob * M1D;
    for (const auto& bs : m_target.boundaries()) {
        index_t p_idx = bs.patch;
        SideInfo si = getSideInfo(bs.side());
        int nBnd = si.par ? (m_n - 1) : 0;
        int off = p_idx * nCP;

        gsMatrix<T> ptsB(2, nQ1);
        ptsB.row(si.nDir).setConstant(si.par ? T(1) : T(0));
        ptsB.row(si.tDir) = qNodes1D.row(0);
        gsMatrix<T> valB = m_target.patch(p_idx).eval(ptsB);

        for (int j1 = 0; j1 < m_n; ++j1) {
            int loc1 = getCPIdx(si.nDir, nBnd, j1);
            int a1 = off + loc1;
            for (int j2 = 0; j2 < m_n; ++j2) {
                int loc2 = getCPIdx(si.nDir, nBnd, j2);
                Q_patches[p_idx](loc1, loc2) += M1_bnd(j1, j2);
            }
            for (short_t d = 0; d < m_geoDim; ++d) {
                T f_b = 0;
                for (int q = 0; q < nQ1; ++q)
                    f_b += qWeights1D(q) * valB(d, q) * B_val1D(q, j1);
                F_glob(a1, d) += wBndGlob * f_b;
            }
        }
    }

    gsMatrix<T> greville = bbasis_p.anchors();
    std::vector<std::vector<std::pair<int, T>>> constr_rows;

    int if_idx = 0;
    for (auto it = m_target.iBegin(); it != m_target.iEnd(); ++it, ++if_idx) {
        const auto& bi = *it;
        const auto& gd = m_gluing[if_idx];
        patchSide ps1 = bi.first();
        patchSide ps2 = bi.second();
        bool reversed2 = gd.reversed2;

        short_t nDir1 = ps1.direction();
        short_t nDir2 = ps2.direction();

        int bndIdx1 = ps1.parameter() ? (m_n - 1) : 0;
        int inIdx1  = ps1.parameter() ? (m_n - 2) : 1;
        int bndIdx2 = ps2.parameter() ? (m_n - 1) : 0;
        int inIdx2  = ps2.parameter() ? (m_n - 2) : 1;

        T dn_scale1 = m_p / (ps1.parameter() ? (m_kv[m_n + m_p] - m_kv[m_n - 1]) : (m_kv[m_p + 1] - m_kv[1]));
        T dn_scale2 = m_p / (ps2.parameter() ? (m_kv[m_n + m_p] - m_kv[m_n - 1]) : (m_kv[m_p + 1] - m_kv[1]));

        int f0_offset = nP * nCP + if_idx * (n0 + n1);
        int f1_offset = f0_offset + n0;

        int p1_off = ps1.patch * nCP;
        int p2_off = ps2.patch * nCP;

        for (int j = 0; j < m_n; ++j) {
            T xi = greville(0, j);

            gsMatrix<T> pt_xi(1, 1); pt_xi(0, 0) = xi;
            gsMatrix<T> v_p1 = bbasis_p.eval(pt_xi);
            gsMatrix<index_t> act_p1 = bbasis_p.active(pt_xi);

            T xi2 = reversed2 ? (T(1) - xi) : xi;
            gsMatrix<T> pt_xi2(1, 1); pt_xi2(0, 0) = xi2;
            gsMatrix<T> v_p2 = bbasis_p.eval(pt_xi2);
            gsMatrix<index_t> act_p2 = bbasis_p.active(pt_xi2);

            gsMatrix<T> v_0 = bbasis_0.eval(pt_xi);
            gsMatrix<T> d_0 = bbasis_0.deriv(pt_xi);
            gsMatrix<index_t> act_0 = bbasis_0.active(pt_xi);

            gsMatrix<T> v_1 = bbasis_1.eval(pt_xi);
            gsMatrix<index_t> act_1 = bbasis_1.active(pt_xi);

            // Constraint 1: f0^(m, m1)(xi) - f0^(m)(xi) = 0
            std::vector<std::pair<int, T>> r1;
            for (int i = 0; i < act_p1.rows(); ++i) {
                int tIdx1 = act_p1(i, 0);
                int a1 = p1_off + getCPIdx(nDir1, bndIdx1, tIdx1);
                r1.push_back({a1, v_p1(i, 0)});
            }
            for (int i = 0; i < act_0.rows(); ++i) {
                int idx0 = act_0(i, 0);
                r1.push_back({f0_offset + idx0, -v_0(i, 0)});
            }
            constr_rows.push_back(r1);

            // Constraint 2: f0^(m, m2)(xi2) - f0^(m)(xi) = 0
            std::vector<std::pair<int, T>> r2;
            for (int i = 0; i < act_p2.rows(); ++i) {
                int tIdx2 = act_p2(i, 0);
                int a2 = p2_off + getCPIdx(nDir2, bndIdx2, tIdx2);
                r2.push_back({a2, v_p2(i, 0)});
            }
            for (int i = 0; i < act_0.rows(); ++i) {
                int idx0 = act_0(i, 0);
                r2.push_back({f0_offset + idx0, -v_0(i, 0)});
            }
            constr_rows.push_back(r2);

            // Gluing functions at xi
            T a1 = gd.alpha1(xi);
            T a2 = gd.alpha2(xi);
            T b1 = gd.beta1(xi);
            T b2 = gd.beta2(xi);

            // Constraint 3: d_in1 - b1 (f0)' - a1 f1 = 0
            std::vector<std::pair<int, T>> r3;
            for (int i = 0; i < act_p1.rows(); ++i) {
                int tIdx1 = act_p1(i, 0);
                int a1_bnd = p1_off + getCPIdx(nDir1, bndIdx1, tIdx1);
                int a1_in  = p1_off + getCPIdx(nDir1, inIdx1,  tIdx1);
                r3.push_back({a1_bnd, -dn_scale1 * v_p1(i, 0)});
                r3.push_back({a1_in,   dn_scale1 * v_p1(i, 0)});
            }
            for (int i = 0; i < act_0.rows(); ++i) {
                int idx0 = act_0(i, 0);
                r3.push_back({f0_offset + idx0, -b1 * d_0(i, 0)});
            }
            for (int i = 0; i < act_1.rows(); ++i) {
                int idx1 = act_1(i, 0);
                r3.push_back({f1_offset + idx1, -a1 * v_1(i, 0)});
            }
            constr_rows.push_back(r3);

            // Constraint 4: d_in2 - b2 (f0)' + a2 f1 = 0
            std::vector<std::pair<int, T>> r4;
            for (int i = 0; i < act_p2.rows(); ++i) {
                int tIdx2 = act_p2(i, 0);
                int a2_bnd = p2_off + getCPIdx(nDir2, bndIdx2, tIdx2);
                int a2_in  = p2_off + getCPIdx(nDir2, inIdx2,  tIdx2);
                r4.push_back({a2_bnd, -dn_scale2 * v_p2(i, 0)});
                r4.push_back({a2_in,   dn_scale2 * v_p2(i, 0)});
            }
            for (int i = 0; i < act_0.rows(); ++i) {
                int idx0 = act_0(i, 0);
                r4.push_back({f0_offset + idx0, -b2 * d_0(i, 0)});
            }
            for (int i = 0; i < act_1.rows(); ++i) {
                int idx1 = act_1(i, 0);
                r4.push_back({f1_offset + idx1, a2 * v_1(i, 0)});
            }
            constr_rows.push_back(r4);
        }
    }

    int numC = constr_rows.size();
    std::vector<int> indepRows;
    if (numC > 0) {
        gsMatrix<T> Ct(nVar, numC); Ct.setZero();
        for (int c = 0; c < numC; ++c) {
            for (const auto& item : constr_rows[c])
                Ct(item.first, c) += item.second;
            const T colNorm = Ct.col(c).norm();
            if (colNorm > T(1e-14))
                Ct.col(c) /= colNorm;
        }
        // Select a full-row-rank subset of sparse constraints so closed interior
        // vertex loops never render the KKT matrix singular.
        gsEigen::ColPivHouseholderQR<gsMatrix<T>> qr(Ct);
        qr.setThreshold(T(1e-10));
        const index_t rankC = qr.rank();
        const auto& perm = qr.colsPermutation().indices();
        indepRows.reserve(rankC);
        for (index_t k = 0; k < rankC; ++k)
            indepRows.push_back(perm(k));
    }

    const int numCIndep = static_cast<int>(indepRows.size());
    const int nTot = nVar + numCIndep;

    gsSparseEntries<T> kktEntries;
    kktEntries.reserve(nP * nCP * nCP + 2 * numCIndep * 16);
    for (index_t p_idx = 0; p_idx < nP; ++p_idx) {
        const int off = p_idx * nCP;
        for (int a = 0; a < nCP; ++a)
            for (int b = 0; b < nCP; ++b)
                if (std::abs(Q_patches[p_idx](a, b)) > T(1e-16))
                    kktEntries.add(off + a, off + b, Q_patches[p_idx](a, b));
    }
    for (int r = 0; r < numCIndep; ++r) {
        const auto& row = constr_rows[indepRows[r]];
        for (const auto& item : row) {
            if (std::abs(item.second) > T(1e-16)) {
                kktEntries.add(nVar + r, item.first, item.second);
                kktEntries.add(item.first, nVar + r, item.second);
            }
        }
    }

    gsSparseMatrix<T> KKT(nTot, nTot);
    KKT.setFrom(kktEntries);

    gsMatrix<T> rhs(nTot, m_geoDim);
    rhs.setZero();
    rhs.topRows(nVar) = F_glob;

    gsMatrix<T> sol;
    gsEigen::SparseLU<gsSparseMatrix<T>> sparseLu;
    sparseLu.compute(KKT);
    if (sparseLu.info() == gsEigen::Success) {
        sol = sparseLu.solve(rhs);
    } else {
        gsMatrix<T> denseKKT(KKT);
        sol = denseKKT.fullPivLu().solve(rhs);
    }

    m_result.clear();
    for (index_t p_idx = 0; p_idx < nP; ++p_idx) {
        gsMatrix<T> cp(nCP, m_geoDim);
        for (int a = 0; a < nCP; ++a)
            cp.row(a) = sol.row(p_idx * nCP + a);
        m_patchCP[p_idx] = cp;
        gsTensorBSplineBasis<2, T> tBasis(m_kv, m_kv);
        typename gsGeometry<T>::uPtr patchGeom(new gsTensorBSpline<2, T>(tBasis, cp));
        m_result.addPatch(give(patchGeom));
    }

    for (size_t b = 0; b < m_target.boundaries().size(); ++b)
        m_result.addBoundary(m_target.boundaries()[b]);
    for (auto it = m_target.iBegin(); it != m_target.iEnd(); ++it)
        m_result.addInterface(*it);

    m_result.computeTopology();
}

template <class T>
bool gsASG1LocalFitter<T>::checkASG1(T eps) const
{
    m_max_res_a = T(0);
    m_max_res_b = T(0);
    m_max_res_n = T(0);
    bool all_valid = true;

    for (auto it = m_result.iBegin(); it != m_result.iEnd(); ++it) {
        auto S_samp = sampleInterface(m_result, *it, 5);
        m_max_res_n = std::max(m_max_res_n, S_samp.maxNormalAngle);
        if (m_geoDim == 3 && S_samp.maxNormalAngle * S_samp.maxNormalAngle > eps)
            all_valid = false;
        if (S_samp.D1.minCoeff() * S_samp.D1.maxCoeff() <= 0 ||
            S_samp.D2.minCoeff() * S_samp.D2.maxCoeff() <= 0) {
            all_valid = false;
            continue;
        }
        const bool sameSign = (S_samp.D1.minCoeff() * S_samp.D2.minCoeff() > 0);
        if (sameSign) S_samp.D1 = -S_samp.D1;

        index_t N = S_samp.size();
        gsMatrix<T> Phi(N, 4);
        gsVector<T> p0 = T(1) - S_samp.nodes.row(0).array();
        gsVector<T> p1 = S_samp.nodes.row(0).array();
        Phi.col(0) = p0.cwiseProduct(S_samp.D2);
        Phi.col(1) = p1.cwiseProduct(S_samp.D2);
        Phi.col(2) = p0.cwiseProduct(S_samp.D1);
        Phi.col(3) = p1.cwiseProduct(S_samp.D1);
        gsMatrix<T> G = Phi.transpose() * S_samp.weights.asDiagonal() * Phi;
        const T tikh = T(1e-10) * G.diagonal().cwiseAbs().maxCoeff();
        G(0,0) += tikh; G(1,1) += tikh; G(0,1) -= tikh; G(1,0) -= tikh;
        G(2,2) += tikh; G(3,3) += tikh; G(2,3) -= tikh; G(3,2) -= tikh;
        gsVector<T> c(4); c.setConstant(T(0.5));
        gsMatrix<T> K(5,5); K.setZero();
        K.block(0,0,4,4) = T(2) * G;
        K.block(0,4,4,1) = c;
        K.block(4,0,1,4) = c.transpose();
        gsVector<T> r(5); r.setZero(); r(4) = 1;
        gsVector<T> alpha = K.fullPivLu().solve(r).head(4);
        T res_a = S_samp.integrate((Phi * alpha).array().square());

        gsMatrix<T> Psi = Phi;
        Psi.rightCols(2) *= -1;
        gsMatrix<T> H = Psi.transpose() * S_samp.weights.asDiagonal() * Psi;
        gsVector<T> d = Psi.transpose() * S_samp.weights.asDiagonal() * S_samp.D3;
        const T tikhB = T(1e-10) * H.diagonal().cwiseAbs().maxCoeff();
        H(0,0) += tikhB; H(2,2) += tikhB; H(0,2) -= tikhB; H(2,0) -= tikhB;
        H(1,1) += tikhB; H(3,3) += tikhB; H(1,3) -= tikhB; H(3,1) -= tikhB;
        gsVector<T> e(4);
        {
            gsVector<T> alpha1 = alpha[0] * p0 + alpha[1] * p1;
            gsVector<T> alpha2 = alpha[2] * p0 + alpha[3] * p1;
            e(0) =  S_samp.integrate(alpha1.cwiseProduct(p0));
            e(1) =  S_samp.integrate(alpha1.cwiseProduct(p1));
            e(2) = -S_samp.integrate(alpha2.cwiseProduct(p0));
            e(3) = -S_samp.integrate(alpha2.cwiseProduct(p1));
        }
        gsMatrix<T> Kb(5,5); Kb.setZero();
        Kb.block(0,0,4,4) = T(2) * H;
        Kb.block(0,4,4,1) = e;
        Kb.block(4,0,1,4) = e.transpose();
        gsVector<T> rb(5); rb.head(4) = T(2) * d; rb(4) = 0;
        gsVector<T> beta = Kb.fullPivLu().solve(rb).head(4);
        T res_b = S_samp.integrate((Psi * beta - S_samp.D3).array().square());

        m_max_res_a = std::max(m_max_res_a, res_a);
        m_max_res_b = std::max(m_max_res_b, res_b);
        if (res_a > eps || res_b > eps) all_valid = false;
    }

    return all_valid;
}

template <class T>
bool gsASG1LocalFitter<T>::unifyPatchOrientations(gsMultiPatch<T>& mp)
{
    auto flipPatchU = [](gsGeometry<T>& geom) -> bool {
        if (auto* bsp = dynamic_cast<gsTensorBSpline<2, T>*>(&geom)) {
            int nu = bsp->basis().size(0);
            int nv = bsp->basis().size(1);
            gsMatrix<T> coefs = bsp->coefs();
            gsMatrix<T> new_coefs = coefs;
            for (int i = 0; i < nu; ++i)
                for (int j = 0; j < nv; ++j)
                    new_coefs.row(j * nu + (nu - 1 - i)) = coefs.row(j * nu + i);
            bsp->setCoefs(new_coefs);
            return true;
        }
        if (auto* nurbs = dynamic_cast<gsTensorNurbs<2, T>*>(&geom)) {
            int nu = nurbs->basis().size(0);
            int nv = nurbs->basis().size(1);
            gsMatrix<T> coefs = nurbs->coefs();
            gsMatrix<T> weights = nurbs->weights();
            gsMatrix<T> new_coefs = coefs;
            gsMatrix<T> new_weights = weights;
            for (int i = 0; i < nu; ++i) {
                for (int j = 0; j < nv; ++j) {
                    int old_idx = j * nu + i;
                    int new_idx = j * nu + (nu - 1 - i);
                    new_coefs.row(new_idx) = coefs.row(old_idx);
                    new_weights(new_idx, 0) = weights(old_idx, 0);
                }
            }
            nurbs->setCoefs(new_coefs);
            nurbs->weights() = new_weights;
            return true;
        }
        return false;
    };

    int flippedCount = 0;
    if (mp.targetDim() == 3) {
        const index_t nP = mp.nPatches();
        std::vector<bool> visited(nP, false);
        std::vector<bool> flipped(nP, false);

        auto evalSurfaceNormal = [](const gsGeometry<T>& g, T u, T v) -> gsVector<T> {
            gsMatrix<T> pt(2, 1);
            pt << u, v;
            gsMatrix<T> d = g.deriv(pt);
            gsVector<T> du(3), dv(3), n(3);
            for (int c = 0; c < 3; ++c) {
                du(c) = d(2 * c + 0, 0);
                dv(c) = d(2 * c + 1, 0);
            }
            n(0) = du(1) * dv(2) - du(2) * dv(1);
            n(1) = du(2) * dv(0) - du(0) * dv(2);
            n(2) = du(0) * dv(1) - du(1) * dv(0);
            return n;
        };

        auto evalSideMidpoint = [](const patchSide& ps, bool rev, T& u, T& v) {
            T ts;
            mapInterfaceParam(ps, rev, T(0.5), u, v, ts);
        };

        for (index_t root = 0; root < nP; ++root) {
            if (visited[root]) continue;
            visited[root] = true;
            std::queue<index_t> q;
            q.push(root);
            while (!q.empty()) {
                index_t cur = q.front();
                q.pop();
                for (auto it = mp.iBegin(); it != mp.iEnd(); ++it) {
                    const boundaryInterface& bi = *it;
                    patchSide psCur, psNbr;
                    bool curIsFirst = false;
                    if (bi.first().patch == cur) {
                        psCur = bi.first(); psNbr = bi.second(); curIsFirst = true;
                    } else if (bi.second().patch == cur) {
                        psCur = bi.second(); psNbr = bi.first(); curIsFirst = false;
                    } else {
                        continue;
                    }
                    const index_t nbr = psNbr.patch;
                    if (visited[nbr]) continue;

                    SideInfo si1 = getSideInfo(bi.first().side());
                    const bool rev2 = !bi.dirOrientation(bi.first(), si1.tDir);
                    T u1, v1, u2, v2;
                    evalSideMidpoint(bi.first(), false, u1, v1);
                    evalSideMidpoint(bi.second(), rev2, u2, v2);
                    gsVector<T> nCur = curIsFirst
                        ? evalSurfaceNormal(mp.patch(bi.first().patch), u1, v1)
                        : evalSurfaceNormal(mp.patch(bi.second().patch), u2, v2);
                    gsVector<T> nNbr = curIsFirst
                        ? evalSurfaceNormal(mp.patch(bi.second().patch), u2, v2)
                        : evalSurfaceNormal(mp.patch(bi.first().patch), u1, v1);
                    if (flipped[cur]) nCur = -nCur;

                    visited[nbr] = true;
                    if (nCur.dot(nNbr) < T(0))
                        flipped[nbr] = true;
                    q.push(nbr);
                }
            }
        }
        for (index_t p = 0; p < nP; ++p) {
            if (flipped[p] && flipPatchU(mp.patch(p)))
                ++flippedCount;
        }
    } else {
        for (size_t p = 0; p < mp.nPatches(); ++p) {
            auto& geom = mp.patch(p);
            T sumDet = 0;
            int count = 0;
            for (int i = 0; i <= 2; ++i) {
                for (int j = 0; j <= 2; ++j) {
                    gsVector<T> u(2);
                    u << i * T(0.5), j * T(0.5);
                    gsMatrix<T> jac;
                    geom.jacobian_into(u, jac);
                    T detJ = jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0);
                    sumDet += detJ;
                    count++;
                }
            }
            if (sumDet / count < T(0) && flipPatchU(geom))
                flippedCount++;
        }
    }
    if (flippedCount > 0) {
        mp.computeTopology();
    }
    return flippedCount > 0;
}

} // namespace gismo

