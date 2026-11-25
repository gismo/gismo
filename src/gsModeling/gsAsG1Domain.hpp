/** @file gsAsG1Domain.hpp

    @brief Derive gluing data for an AS-G1 geometry.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Hasanova, S. Takacs
*/

#pragma once

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gsCore/gsPiecewiseFunction.h>
#include <gsNurbs/gsTensorNurbsBasis.h>

namespace gismo {

// ====================================================================
// Gluing-data helpers
// ====================================================================

template <class T>
inline T det2(const gsVector<T>& a, const gsVector<T>& b)
{ return a(0)*b(1) - a(1)*b(0); }

template <class T>
inline T det3(const gsVector<T>& a, const gsVector<T>& b, const gsVector<T>& n)
{
    const T cx = a(1)*b(2) - a(2)*b(1);
    const T cy = a(2)*b(0) - a(0)*b(2);
    const T cz = a(0)*b(1) - a(1)*b(0);
    return cx*n(0) + cy*n(1) + cz*n(2);
}

template <class T>
inline gsVector<T> partial(const gsMatrix<T>& d, index_t dir, index_t col, short_t geoDim = 2)
{
    gsVector<T> v(geoDim);
    for (short_t c = 0; c < geoDim; ++c)
        v(c) = d(dir + 2 * c, col);
    return v;
}

template <class T>
gsVector<T> breaksOf(const gsGeometry<T>& geo, short_t dir)
{
    const gsBasis<T>* basisPtr = &geo.basis();
    if (const auto* rBasis = dynamic_cast<const gsRationalBasis<gsTensorBSplineBasis<2, T>>*>(basisPtr))
        basisPtr = &rBasis->source();
    const gsBSplineBasis<T>& bb =
        dynamic_cast<const gsBSplineBasis<T>&>(basisPtr->component(dir));
    const std::vector<T> &brks = bb.knots().breaks();
    gsVector<T> result(brks.size());
    result.assign(brks.begin(), brks.end());
    return result;
}

/// Sampled determinants on the interface (patch-1 tangential parameter
/// mapped to [0,1]).
template <class T>
struct InterfaceSamples
{
    gsVector<T> D1, D2, D3;
    gsMatrix<T> nodes;
    gsVector<T> weights;
    index_t size() const { return D1.size(); }
    T integrate(const gsVector<T>& f) const { return weights.dot(f); }
    bool flipped;
    T maxNormalAngle = T(0); // max |g2n . n_S| / ||g2n|| for 3D surfaces (0 in 2D)
};

/// Sample D1,D2,D3 along the interface.
template <class T>
InterfaceSamples<T> sampleInterface(
    const gsMultiPatch<T>& mp,
    const boundaryInterface& interf,
    index_t numGaussPerSpan = 0)
{
    InterfaceSamples<T> S;
    const patchSide ps1 = interf.first(), ps2 = interf.second();
    const gsGeometry<T>& g1 = mp.patch(ps1.patch);
    const gsGeometry<T>& g2 = mp.patch(ps2.patch);

    const short_t nDir1 = ps1.direction(), tDir1 = 1 - nDir1;
    const short_t nDir2 = ps2.direction(), tDir2 = 1 - nDir2;
    const bool par1 = ps1.parameter(), par2 = ps2.parameter();
    const T s1 = par1 ? T(-1) : T(1);
    const T s2 = par2 ? T(-1) : T(1);

    const gsMatrix<T> sup1 = g1.support();
    const gsMatrix<T> sup2 = g2.support();
    const T n1 = sup1(nDir1, par1 ? 1 : 0);
    const T n2 = sup2(nDir2, par2 ? 1 : 0);

    S.flipped = !interf.dirOrientation(ps1, tDir1);

    const T t1a = sup1(tDir1, 0), t1b = sup1(tDir1, 1);
    T t2a = sup2(tDir2, 0), t2b = sup2(tDir2, 1);
    if (S.flipped)
        std::swap(t2a, t2b);

    // normalized breakpoints
    const gsVector<T> breaks1 = (breaksOf(g1, tDir1).array() - t1a) / (t1b - t1a);
    const gsVector<T> breaks2 = (breaksOf(g2, tDir2).array() - t2a) / (t2b - t2a);

    // merge normalized breakpoints
    std::set<T> merged(breaks1.begin(), breaks1.end());
    merged.insert(breaks2.begin(), breaks2.end());
    const std::vector<T> breaks(merged.begin(), merged.end());

    // derive Gauss nodes
    const short_t deg = std::max(g1.basis().degree(tDir1), g2.basis().degree(tDir2));
    const index_t nGauss = (numGaussPerSpan > 0) ? numGaussPerSpan : 2*deg + 1;
    const gsGaussRule<T> rule(nGauss);
    rule.mapToAll(breaks, S.nodes, S.weights);

    // map nodes back
    const index_t N = S.nodes.cols();
    gsMatrix<T> pts1(2,N), pts2(2,N);
    pts1.row(nDir1).setConstant(n1);
    pts2.row(nDir2).setConstant(n2);
    pts1.row(tDir1) = t1a + S.nodes.array() * (t1b - t1a);
    pts2.row(tDir2) = t2a + S.nodes.array() * (t2b - t2a);

    // Jacobian of geometry
    const gsMatrix<T> d1 = g1.deriv(pts1);
    const gsMatrix<T> d2 = g2.deriv(pts2);
    const short_t geoDim = mp.targetDim();
    GISMO_ENSURE(geoDim == 2 || geoDim == 3,
                 "sampleInterface supports 2D planar domains and 3D surfaces, got geoDim=" << geoDim);

    S.D1.resize(N); S.D2.resize(N); S.D3.resize(N);
    S.maxNormalAngle = T(0);
    for (index_t i = 0; i < N; ++i)
    {
        const gsVector<T> g1n = partial(d1, nDir1, i, geoDim);
        const gsVector<T> g1t = partial(d1, tDir1, i, geoDim);
        const gsVector<T> g2n = partial(d2, nDir2, i, geoDim);
        const gsVector<T> g2t = partial(d2, tDir2, i, geoDim);
        if (geoDim == 2)
        {
            S.D1(i) = s1      * det2(g1n, g1t);
            S.D2(i) = s2      * det2(g2n, g2t);
            S.D3(i) = s1 * s2 * det2(g1n, g2n);
        }
        else
        {
            gsVector<T> nS(3);
            nS(0) = g1n(1)*g1t(2) - g1n(2)*g1t(1);
            nS(1) = g1n(2)*g1t(0) - g1n(0)*g1t(2);
            nS(2) = g1n(0)*g1t(1) - g1n(1)*g1t(0);
            const T nNorm = nS.norm();
            if (nNorm > T(1e-15))
            {
                nS /= nNorm;
                const T g2nNorm = g2n.norm();
                if (g2nNorm > T(1e-15))
                {
                    const T sinAngle = std::abs(g2n.dot(nS)) / g2nNorm;
                    if (sinAngle > S.maxNormalAngle)
                        S.maxNormalAngle = sinAngle;
                }
                S.D1(i) = s1      * det3(g1n, g1t, nS);
                S.D2(i) = s2      * det3(g2n, g2t, nS);
                S.D3(i) = s1 * s2 * det3(g1n, g2n, nS);
            }
            else
            {
                S.D1(i) = T(0);
                S.D2(i) = T(0);
                S.D3(i) = T(0);
            }
        }
    }
    return S;
}

/// Single-solve linear least-squares for (alpha, beta) on the basis
/// {1-t, t}.  alpha is normalised by integral(alpha_1+alpha_2)=1, beta
/// is gauged by integral(alpha_1*beta_1 - alpha_2*beta_2)=0.  Tikhonov
/// bias resolves null directions (toward constant alpha and beta_1=beta_2).
template <class T>
struct SolveResult { gsVector<T> alpha; gsVector<T> beta; };

template <class T>
SolveResult<T> solveLinearGluing(const InterfaceSamples<T>& S, real_t eps)
{
    SolveResult<T> out;
    const index_t N = S.size();

    const gsVector<T> p0 = 1 - S.nodes.transpose().array();
    const gsVector<T> p1 = S.nodes.transpose();

    gsMatrix<T> Phi(N, 4);
    Phi.col(0) = p0.cwiseProduct(S.D2);
    Phi.col(1) = p1.cwiseProduct(S.D2);
    Phi.col(2) = p0.cwiseProduct(S.D1);
    Phi.col(3) = p1.cwiseProduct(S.D1);

    gsMatrix<T> G = Phi.transpose() * S.weights.asDiagonal() * Phi;

    // Tikhonov bias toward constant alpha
    const T tikh = T(1e-10) * G.diagonal().cwiseAbs().maxCoeff();
    G(0,0) += tikh; G(1,1) += tikh; G(0,1) -= tikh; G(1,0) -= tikh;
    G(2,2) += tikh; G(3,3) += tikh; G(2,3) -= tikh; G(3,2) -= tikh;

    gsVector<T> c(4); c.setConstant(T(0.5));
    gsMatrix<T> K(5,5); K.setZero();
    K.block(0,0,4,4) = 2 * G;
    K.block(0,4,4,1) = c;
    K.block(4,0,1,4) = c.transpose();
    gsVector<T> r(5); r.setZero(); r(4) = 1;
    out.alpha = K.fullPivLu().solve(r).head(4);

    GISMO_ENSURE (S.integrate((Phi * out.alpha).array().square()) <= eps, "Not AS-G1");

    // beta solve: Psi = [(1-t)D2, t*D2, -(1-t)D1, -t*D1]
    gsMatrix<T> Psi = Phi;
    Psi.rightCols(2) *= -1;

    gsMatrix<T> H = Psi.transpose() * S.weights.asDiagonal() * Psi;
    gsVector<T> d = Psi.transpose() * S.weights.asDiagonal() * S.D3;

    // Tikhonov bias toward beta_1 = beta_2
    const T tikhB = T(1e-10) * H.diagonal().cwiseAbs().maxCoeff();
    H(0,0) += tikhB; H(2,2) += tikhB; H(0,2) -= tikhB; H(2,0) -= tikhB;
    H(1,1) += tikhB; H(3,3) += tikhB; H(1,3) -= tikhB; H(3,1) -= tikhB;

    gsVector<T> e(4);
    {
        gsVector<T> alpha1 = out.alpha[0] * p0 + out.alpha[1] * p1;
        gsVector<T> alpha2 = out.alpha[2] * p0 + out.alpha[3] * p1;

        e(0) =  S.integrate(alpha1.cwiseProduct(p0));
        e(1) =  S.integrate(alpha1.cwiseProduct(p1));
        e(2) = -S.integrate(alpha2.cwiseProduct(p0));
        e(3) = -S.integrate(alpha2.cwiseProduct(p1));
    }

    gsMatrix<T> Kb(5,5); Kb.setZero();
    Kb.block(0,0,4,4) = 2 * H;
    Kb.block(0,4,4,1) = e;
    Kb.block(4,0,1,4) = e.transpose();
    gsVector<T> rb(5); rb.head(4) = 2*d; rb(4) = 0;
    out.beta = Kb.fullPivLu().solve(rb).head(4);

    GISMO_ENSURE (S.integrate((Psi * out.beta - S.D3).array().square()) <= eps, "Not AS-G1");

    return out;
}

/// Per-interface driver.
///
/// Returns 8 numbers in patch-1 tangential parametrisation:
///   [ a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1 ]
///
template <class T>
gsVector<T> computeGluingDataForInterface(
    const gsMultiPatch<T>& mp,
    const boundaryInterface& interf,
    const T eps = 1e-8,
    index_t numGaussPerSpan = 0)
{
    gsVector<T> result(8);

    InterfaceSamples<T> S = sampleInterface(mp, interf, numGaussPerSpan);

    GISMO_ENSURE (S.D1.minCoeff() * S.D1.maxCoeff() > 0 &&
        S.D2.minCoeff() * S.D2.maxCoeff() > 0, "Not AS-G1");

    const bool sameSign = (S.D1.minCoeff() * S.D2.minCoeff() > 0);
    if (sameSign) S.D1 = -S.D1;

    SolveResult<T> r = solveLinearGluing(S, eps);

    if (sameSign)
    {
        r.alpha.topRows(2) *= -1;
        r.beta.bottomRows(2) *= -1;
    }

    if (S.flipped)
    {
        // Patch-2 tangent is flipped, so we evaluate the patch-2 alpha
        // and beta at the reversed endpoint pairing (a21 at gd-t=0,
        // a20 at gd-t=1).
        result[0] =  r.alpha[0]; result[1] =  r.alpha[1];
        result[2] =  r.beta [0]; result[3] =  r.beta [1];
        result[4] =  r.alpha[3]; result[5] =  r.alpha[2];
        result[6] =  r.beta [3]; result[7] =  r.beta [2];
    }
    else
    {
        result[0] = -r.alpha[0]; result[1] = -r.alpha[1];
        result[2] = -r.beta [0]; result[3] = -r.beta [1];
        result[4] =  r.alpha[2]; result[5] =  r.alpha[3];
        result[6] = -r.beta [2]; result[7] = -r.beta [3];
    }

    // There is another fix to be made: In the paper we have written
    // that the tangential vector is obtained by rotating the normal
    // vector; the gluing data computation code uses always (+1,0) and (0,+1).
    const patchSide ps1 = interf.first(), ps2 = interf.second();
    if (ps1.parameter()!=ps1.direction())
    {
        result[2] *= -1;
        result[3] *= -1;
    }
    if (ps2.parameter()!=ps2.direction())
    {
        result[6] *= -1;
        result[7] *= -1;
    }
    return result;
}

/// Normalize linear gluing data [a1_0, a1_1, b1_0, b1_1, a2_0, a2_1, b2_0, b2_1]
/// according to Appendix A of Hasanova, S. Takacs, T. Takacs (2026):
/// rescale (alpha_1, alpha_2) so min_{xi in [0,1]} min(|alpha_1(xi)|, |alpha_2(xi)|) = 1.
template <class T>
void normalizeGluingDataAppendixA(gsVector<T>& gd)
{
    GISMO_ASSERT(gd.size() == 8, "Expected 8-entry gluing data vector.");
    const T minA = std::min(
        std::min(std::abs(gd[0]), std::abs(gd[1])),
        std::min(std::abs(gd[4]), std::abs(gd[5])));
    if (minA > T(1e-14))
    {
        gd[0] /= minA; gd[1] /= minA;
        gd[4] /= minA; gd[5] /= minA;
    }
}


// ====================================================================
// Multi-patch driver: assemble the (nPatches x 16) matrix
//   columns per side s in {west=0, east=1, south=2, north=3}:
//     4*s+0 = alpha_0    4*s+1 = alpha_1
//     4*s+2 = beta_0     4*s+3 = beta_1
// Boundary sides keep their default (alpha=1, beta=0).
// ====================================================================
template <class T>
gsMatrix<T> computeGluingData(
    const gsMultiPatch<T>& mp,
    T eps = 1e-8,
    index_t numGaussPerSpan = 0)
{
    const index_t nP = mp.nPatches();
    gsMatrix<T> M(nP, 16);
    for (index_t p = 0; p < nP; ++p)
        for (index_t s = 0; s < 4; ++s) {
            M(p, 4*s+0) = T(1); M(p, 4*s+1) = T(1);
            M(p, 4*s+2) = T(0); M(p, 4*s+3) = T(0);
        }

    auto sideCol = [](const patchSide& ps) -> index_t
        { return 4 * (static_cast<index_t>(ps.side()) - 1); };

    for (auto it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const boundaryInterface& interf = *it;

        gsVector<T> gd = computeGluingDataForInterface(
            mp, interf, eps, numGaussPerSpan);

        const index_t c1 = sideCol(interf.first());
        const index_t c2 = sideCol(interf.second());
        M.row(interf.first().patch ).segment(c1, 4) = gd.head(4).transpose();
        M.row(interf.second().patch).segment(c2, 4) = gd.tail(4).transpose();
    }
    return M;
}

// ====================================================================
// Code for handling boundary conditions at the corners
// ====================================================================


template <typename T>
std::vector<std::vector<patchCorner>> getBoundaryVertices(
    const gsBoxTopology& bt,
    const gsBoundaryConditions<T>& bc,
    const std::string& label)
{
    const short_t dim = 2;
    gsSortedVector<std::vector<patchCorner>> result;
    for (auto it = bc.begin(label); it != bc.end(label); ++it)
    {
        std::vector<patchCorner> corners;
        it->ps.getContainedCorners(dim, corners);
        for (size_t i=0; i<corners.size(); ++i)
        {
            std::vector<patchCorner> cornerGroup;
            bt.getCornerList(corners[i], cornerGroup);
            // sort cornerGroup in order to be recognize matching ones
            std::sort(cornerGroup.begin(), cornerGroup.end());
            result.push_sorted_unique(cornerGroup);
        }

    }
    return result;
}

template <typename T>
gsVector<T> getOuterNormalDerivative(
    const gsMultiPatch<T>& mp,
    const std::vector<patchCorner>& corners,
    T eps = 1e-4)
{
    const short_t dim = 2;
    std::vector<gsVector<T>> normals;
    for (size_t i=0; i<corners.size(); ++i)
    {
        const gsGeometry<T>& geo = mp.patch(corners[i].patch);

        // determine coordinates of corner in parameter domain
        gsMatrix<T> pt(dim,1);
        {
            const gsMatrix<T> supp = geo.support();
            gsVector<bool> parameters = corners[i].parameters(dim);
            for (short_t d=0; d<dim; ++d)
                pt(d,0) = parameters[d] ? supp(d,1) : supp(d,0);
        }

        // determine Jacobian
        gsMatrix<T> deriv;
        geo.deriv_into(pt, deriv);

        // determine cofactor matrix of Jacobian
        gsMatrix<T> jacCf(2,2);
        jacCf(0,0) =  deriv(3,0); // =  jacobian(1,1)
        jacCf(1,0) = -deriv(1,0); // = -jacobian(0,1)
        jacCf(0,1) = -deriv(2,0); // = -jacobian(1,0)
        jacCf(1,1) =  deriv(0,0); // =  jacobian(0,0)

        // determine sides
        std::vector<patchSide> sides;
        corners[i].getContainingSides(dim, sides);
        for (size_t j=0; j<sides.size(); ++j)
        {
            if (mp.isBoundary(sides[j]))
            {
                gsVector<T> paramNormal;
                paramNormal.setZero(dim);
                paramNormal[sides[j].direction()] = sides[j].parameter() ? 1 : -1;
                gsVector<T> normal = jacCf * paramNormal * geo.orientation();
                normal.normalize();
                normals.push_back(normal);
            }
        }
    }
    // Now, select one if it fits and return null vector if the normals disagree
    GISMO_ENSURE (normals.size() == 2, "Expected 2 edges, got "<<normals.size());
    if ( normals[0].dot(normals[1]) >= 1-eps )
        return (normals[0]+normals[1])/2;
    else
    {
        gsVector<T> zero(2);
        zero.setZero();
        return zero;
    }
}

/// Boundary corner frame data needed to enforce C^0 (second biharmonic) or C^1 (first biharmonic)
/// boundary conditions cleanly on both sharp corners (Case A) and smooth boundary vertices (Case B).
template <typename T>
struct gsBoundaryCornerData
{
    bool isSmooth = false;
    gsVector<T> dir1;
    gsVector<T> dir2;
    T curvature = T(0);
    // Optional tangent-plane curvature correction weights in the orthonormal basis
    // E = [e1, e2] of span(dir1, dir2), so that:
    //   DOF 3 = Hess_S u(d1, d1) + k0_e1 * (grad_S u . e1) + k0_e2 * (grad_S u . e2)
    //   DOF 4 = Hess_S u(d2, d2) + k1_e1 * (grad_S u . e1) + k1_e2 * (grad_S u . e2)
    T k0_e1 = T(0);
    T k0_e2 = T(0);
    T k1_e1 = T(0);
    T k1_e2 = T(0);
};

template <typename T>
gsBoundaryCornerData<T> getBoundaryCornerData(
    const gsMultiPatch<T>& mp,
    const std::vector<patchCorner>& corners,
    T eps = 1e-4)
{
    const short_t dim = 2;
    std::vector<gsVector<T>> normals;
    std::vector<gsVector<T>> tangents;
    std::vector<gsVector<T>> accs;
    std::vector<T> velSqs;

    for (size_t i = 0; i < corners.size(); ++i)
    {
        const gsGeometry<T>& geo = mp.patch(corners[i].patch);

        gsMatrix<T> pt(dim, 1);
        const gsMatrix<T> supp = geo.support();
        gsVector<bool> parameters = corners[i].parameters(dim);
        for (short_t d = 0; d < dim; ++d)
            pt(d, 0) = parameters[d] ? supp(d, 1) : supp(d, 0);

        gsMatrix<T> deriv, deriv2;
        geo.deriv_into(pt, deriv);
        geo.deriv2_into(pt, deriv2);

        gsMatrix<T> jacCf(2, 2);
        jacCf(0, 0) =  deriv(3, 0);
        jacCf(1, 0) = -deriv(1, 0);
        jacCf(0, 1) = -deriv(2, 0);
        jacCf(1, 1) =  deriv(0, 0);

        std::vector<patchSide> sides;
        corners[i].getContainingSides(dim, sides);
        for (size_t j = 0; j < sides.size(); ++j)
        {
            if (mp.isBoundary(sides[j]))
            {
                gsVector<T> paramNormal(dim);
                paramNormal.setZero();
                paramNormal[sides[j].direction()] = sides[j].parameter() ? T(1) : T(-1);
                gsVector<T> normal = jacCf * paramNormal * geo.orientation();
                normal.normalize();
                normals.push_back(normal);

                const short_t tDir = 1 - sides[j].direction();
                gsVector<T> vel(2);
                vel(0) = deriv(tDir, 0);
                vel(1) = deriv(tDir + 2, 0);
                const T velSq = vel.squaredNorm();
                gsVector<T> tangent = vel * (parameters[tDir] ? T(-1) : T(1)) / std::sqrt(velSq);
                tangents.push_back(tangent);

                gsVector<T> acc(2);
                acc(0) = deriv2(tDir, 0);
                acc(1) = deriv2(tDir + 3, 0);
                accs.push_back(acc);
                velSqs.push_back(velSq);
            }
        }
    }

    GISMO_ENSURE(normals.size() == 2, "Expected 2 boundary edges at boundary vertex, got " << normals.size());

    gsBoundaryCornerData<T> out;
    if (normals[0].dot(normals[1]) >= T(1) - eps)
    {
        // Case B: smooth boundary vertex (unique normal and tangent)
        out.isSmooth = true;
        out.dir1 = (normals[0] + normals[1]).normalized();
        out.dir2.resize(2);
        out.dir2(0) =  out.dir1(1);
        out.dir2(1) = -out.dir1(0);
        const T ck0 = accs[0].dot(out.dir1) / velSqs[0];
        const T ck1 = accs[1].dot(out.dir1) / velSqs[1];
        out.curvature = T(0.5) * (ck0 + ck1);
    }
    else if (std::abs(det2(tangents[0], tangents[1])) >= eps)
    {
        // Case A: sharp boundary corner with linearly independent tangents [t1, t2]
        out.isSmooth = false;
        out.dir1 = tangents[0];
        out.dir2 = tangents[1];
        out.curvature = T(0);
    }
    else
    {
        // Degenerate 0-angle / 360-angle cusp: tangents are collinear
        out.isSmooth = true;
        out.dir1 = normals[0];
        out.dir2.resize(2);
        out.dir2(0) =  out.dir1(1);
        out.dir2(1) = -out.dir1(0);
        const T ck0 = accs[0].dot(out.dir1) / velSqs[0];
        const T ck1 = accs[1].dot(out.dir1) / velSqs[1];
        out.curvature = T(0.5) * (ck0 + ck1);
    }
    return out;
}

/// 3D surface boundary corner frame data in the shared vertex tangent plane T_v S.
/// For smooth boundary vertices (Case B): dir1 = outward conormal n_co in T_v S,
/// dir2 = boundary unit tangent t in T_v S, curvature = geodesic curvature kappa_g.
/// For sharp boundary corners (Case A): dir1 = t1, dir2 = t2 in T_v S, curvature = 0.
template <typename T>
gsBoundaryCornerData<T> getSurfaceBoundaryCornerData(
    const gsMultiPatch<T>& mp,
    const std::vector<patchCorner>& corners,
    T eps = 1e-4)
{
    const short_t dim = 2;
    GISMO_ENSURE(mp.targetDim() == 3, "getSurfaceBoundaryCornerData expects a 3D surface (targetDim == 3).");

    std::vector<gsVector<T>> conormals;
    std::vector<gsVector<T>> tangents;
    std::vector<gsVector<T>> accs;
    std::vector<T> velSqs;

    for (size_t i = 0; i < corners.size(); ++i)
    {
        const gsGeometry<T>& geo = mp.patch(corners[i].patch);
        gsMatrix<T> pt(dim, 1);
        const gsMatrix<T> supp = geo.support();
        const gsVector<bool> parameters = corners[i].parameters(dim);
        for (short_t d = 0; d < dim; ++d)
            pt(d, 0) = parameters[d] ? supp(d, 1) : supp(d, 0);

        gsMatrix<T> deriv, deriv2;
        geo.deriv_into(pt, deriv);
        geo.deriv2_into(pt, deriv2);

        std::vector<patchSide> sides;
        corners[i].getContainingSides(dim, sides);
        for (size_t j = 0; j < sides.size(); ++j)
        {
            if (mp.isBoundary(sides[j]))
            {
                const short_t nDir = sides[j].direction();
                const short_t tDir = 1 - nDir;

                const gsVector<T> vel = partial(deriv, tDir, 0, 3);
                const T velSq = vel.squaredNorm();
                const gsVector<T> tangent = vel * (parameters[tDir] ? T(-1) : T(1)) / std::sqrt(velSq);
                tangents.push_back(tangent);

                const T nOutSign = sides[j].parameter() ? T(1) : T(-1);
                const gsVector<T> gOut = nOutSign * partial(deriv, nDir, 0, 3);
                gsVector<T> conormal = gOut - gOut.dot(tangent) * tangent;
                conormal.normalize();
                conormals.push_back(conormal);

                gsVector<T> acc(3);
                for (short_t c = 0; c < 3; ++c)
                    acc(c) = deriv2(3 * c + tDir, 0);
                accs.push_back(acc);
                velSqs.push_back(velSq);
            }
        }
    }

    GISMO_ENSURE(conormals.size() == 2,
                 "Expected 2 boundary edges at surface boundary vertex, got " << conormals.size());

    gsBoundaryCornerData<T> out;
    auto cross3 = [](const gsVector<T>& a, const gsVector<T>& b) -> gsVector<T> {
        gsVector<T> r(3);
        r(0) = a(1)*b(2) - a(2)*b(1);
        r(1) = a(2)*b(0) - a(0)*b(2);
        r(2) = a(0)*b(1) - a(1)*b(0);
        return r;
    };

    if (conormals[0].dot(conormals[1]) >= T(1) - eps)
    {
        // Case B: smooth boundary vertex on 3D surface
        out.isSmooth = true;
        out.dir1 = (conormals[0] + conormals[1]).normalized();
        gsVector<T> tanDiff = tangents[0] - tangents[1];
        tanDiff -= tanDiff.dot(out.dir1) * out.dir1;
        if (tanDiff.norm() < T(1e-12))
        {
            tanDiff = tangents[0] - tangents[0].dot(out.dir1) * out.dir1;
        }
        out.dir2 = tanDiff.normalized();
        const T ck0 = accs[0].dot(out.dir1) / velSqs[0];
        const T ck1 = accs[1].dot(out.dir1) / velSqs[1];
        out.curvature = T(0.5) * (ck0 + ck1);
        out.k0_e1 = T(0);
        out.k0_e2 = T(0);
        out.k1_e1 = out.curvature;
        out.k1_e2 = T(0);
    }
    else if (cross3(tangents[0], tangents[1]).norm() >= eps)
    {
        // Case A: sharp boundary corner with linearly independent tangents [t1, t2] in T_v S.
        // Include each boundary edge's geodesic curvature correction along its outward conormal
        // so that DOF 3 = d^2(u|_{e_0})/ds_0^2 and DOF 4 = d^2(u|_{e_1})/ds_1^2 exactly even
        // when one edge is curved and the other is free.
        out.isSmooth = false;
        out.dir1 = tangents[0];
        out.dir2 = tangents[1];
        out.curvature = T(0);

        const gsVector<T> e1 = out.dir1.normalized();
        const gsVector<T> e2 = (out.dir2 - out.dir2.dot(e1) * e1).normalized();
        const T ck0 = accs[0].dot(conormals[0]) / velSqs[0];
        const T ck1 = accs[1].dot(conormals[1]) / velSqs[1];
        out.k0_e1 = ck0 * e1.dot(conormals[0]);
        out.k0_e2 = ck0 * e2.dot(conormals[0]);
        out.k1_e1 = ck1 * e1.dot(conormals[1]);
        out.k1_e2 = ck1 * e2.dot(conormals[1]);
    }
    else
    {
        // Degenerate cusp
        out.isSmooth = true;
        out.dir1 = conormals[0];
        gsVector<T> tProj = tangents[0] - tangents[0].dot(out.dir1) * out.dir1;
        out.dir2 = tProj.normalized();
        const T ck0 = accs[0].dot(out.dir1) / velSqs[0];
        const T ck1 = accs[1].dot(out.dir1) / velSqs[1];
        out.curvature = T(0.5) * (ck0 + ck1);
        out.k0_e1 = T(0);
        out.k0_e2 = T(0);
        out.k1_e1 = out.curvature;
        out.k1_e2 = T(0);
    }
    return out;
}

/// Build shared 3D surface tangent-plane frames for all corners of every patch.
/// Returns a (nPatches x 40) matrix with stride 10 per corner:
///   [d1_x, d1_y, d1_z, d2_x, d2_y, d2_z, k0_e1, k0_e2, k1_e1, k1_e2].
/// Every topological vertex (both interior and boundary) receives a single shared
/// tangent-plane frame across all incident patches, guaranteeing exact 6-jet matching
/// and exact arc-length 2-jets along both smooth and sharp curved boundary edges.
/// Also fills cornerIsSmooth(patch, cornerIdx) (1 for smooth boundary / interior, 0 for sharp corner).
template <typename T>
gsMatrix<T> computeSurfaceCornerFrames(
    const gsMultiPatch<T>& mp,
    gsMatrix<index_t>* cornerIsSmooth = nullptr,
    T eps = 1e-4)
{
    GISMO_ENSURE(mp.targetDim() == 3, "computeSurfaceCornerFrames expects a 3D surface (targetDim == 3).");
    const index_t nP = mp.nPatches();
    gsMatrix<T> frames(nP, 4 * 10);
    frames.setZero();
    if (cornerIsSmooth)
        cornerIsSmooth->setOnes(nP, 4);

    auto cross3 = [](const gsVector<T>& a, const gsVector<T>& b) -> gsVector<T> {
        gsVector<T> r(3);
        r(0) = a(1)*b(2) - a(2)*b(1);
        r(1) = a(2)*b(0) - a(0)*b(2);
        r(2) = a(0)*b(1) - a(1)*b(0);
        return r;
    };

    const auto& vertices = mp.topology().vertices();
    for (size_t v = 0; v < vertices.size(); ++v)
    {
        const std::vector<patchCorner>& grp = vertices[v];
        if (grp.empty()) continue;

        // Count boundary sides incident to this vertex
        index_t nBndSides = 0;
        for (size_t i = 0; i < grp.size(); ++i)
        {
            std::vector<patchSide> sides;
            grp[i].getContainingSides(2, sides);
            for (size_t j = 0; j < sides.size(); ++j)
                if (mp.isBoundary(sides[j]))
                    ++nBndSides;
        }

        gsVector<T> d1(3), d2(3);
        T k0_e1 = T(0), k0_e2 = T(0), k1_e1 = T(0), k1_e2 = T(0);
        bool isSmooth = true;

        if (nBndSides == 2)
        {
            const gsBoundaryCornerData<T> cd = getSurfaceBoundaryCornerData(mp, grp, eps);
            d1 = cd.dir1;
            d2 = cd.dir2;
            k0_e1 = cd.k0_e1;
            k0_e2 = cd.k0_e2;
            k1_e1 = cd.k1_e1;
            k1_e2 = cd.k1_e2;
            isSmooth = cd.isSmooth;
        }
        else
        {
            // Interior vertex (or non-manifold boundary): build deterministic orthonormal basis of T_v S
            gsVector<T> avgN = gsVector<T>::Zero(3);
            gsVector<T> refN = gsVector<T>::Zero(3);
            for (size_t i = 0; i < grp.size(); ++i)
            {
                const gsGeometry<T>& geo = mp.patch(grp[i].patch);
                gsMatrix<T> pt(2, 1);
                const gsMatrix<T> supp = geo.support();
                const gsVector<bool> params = grp[i].parameters(2);
                pt(0, 0) = params[0] ? supp(0, 1) : supp(0, 0);
                pt(1, 0) = params[1] ? supp(1, 1) : supp(1, 0);
                gsMatrix<T> deriv;
                geo.deriv_into(pt, deriv);
                const gsVector<T> gu = partial(deriv, 0, 0, 3);
                const gsVector<T> gv = partial(deriv, 1, 0, 3);
                gsVector<T> n = cross3(gu, gv);
                const T nn = n.norm();
                if (nn > T(1e-14))
                {
                    n /= nn;
                    if (refN.norm() < T(0.5)) refN = n;
                    else if (n.dot(refN) < T(0)) n = -n;
                    avgN += n;
                }
            }
            if (avgN.norm() > T(1e-14))
                avgN.normalize();
            else
                avgN << T(0), T(0), T(1);

            // Deterministic tangent vector e1 orthogonal to avgN (matches [e_x, e_y] when avgN = +-e_z)
            gsVector<T> trial(3);
            if (std::abs(avgN(0)) < T(0.9))
                trial << T(1), T(0), T(0);
            else
                trial << T(0), T(1), T(0);
            d1 = trial - trial.dot(avgN) * avgN;
            d1.normalize();
            d2 = cross3(avgN, d1);
            d2.normalize();
            isSmooth = true;
        }

        for (size_t i = 0; i < grp.size(); ++i)
        {
            const index_t p = grp[i].patch;
            const index_t cIdx = grp[i].m_index - 1;
            frames(p, 10 * cIdx + 0) = d1(0);
            frames(p, 10 * cIdx + 1) = d1(1);
            frames(p, 10 * cIdx + 2) = d1(2);
            frames(p, 10 * cIdx + 3) = d2(0);
            frames(p, 10 * cIdx + 4) = d2(1);
            frames(p, 10 * cIdx + 5) = d2(2);
            frames(p, 10 * cIdx + 6) = k0_e1;
            frames(p, 10 * cIdx + 7) = k0_e2;
            frames(p, 10 * cIdx + 8) = k1_e1;
            frames(p, 10 * cIdx + 9) = k1_e2;
            if (cornerIsSmooth)
                (*cornerIsSmooth)(p, cIdx) = isSmooth ? 1 : 0;
        }
    }
    return frames;
}

// ====================================================================
// Homogeneous-BC manufactured solutions (B-spline bubble & helpers)
// ====================================================================

/**
 * @brief Scalar function returning the physical Laplacian Delta f = sum_d d^2f/dx_d^2
 *        of a wrapped scalar gsFunction, used to assemble the energy-form RHS
 *        l(v) = (Delta u_ex, Delta v) + beta (u_ex, v).
 */
template <typename T>
class gsLaplacianFunction : public gsFunction<T>
{
    std::shared_ptr<const gsFunction<T>> m_f;

public:
    explicit gsLaplacianFunction(const gsFunction<T> &f)
        : m_f(f.clone().release()) {}

    GISMO_CLONE_FUNCTION(gsLaplacianFunction)

    short_t domainDim() const override { return m_f->domainDim(); }
    short_t targetDim() const override { return 1; }

    void eval_into(const gsMatrix<T> &u, gsMatrix<T> &result) const override
    {
        gsMatrix<T> d2;
        m_f->deriv2_into(u, d2);
        const short_t d = m_f->domainDim();
        result = d2.topRows(d).colwise().sum();
    }
};

/**
 * @brief Smooth C^{q-1} trigonometric bubble supported on [a_x, b_x] x [a_y, b_y]:
 *        u(x, y) = sin^q(pi * (x - a_x) / (b_x - a_x)) * sin^q(pi * (y - a_y) / (b_y - a_y))
 *        inside the rectangle and 0 outside.
 */
template <typename T>
class gsTrigBubbleFunction : public gsFunction<T>
{
    T m_ax, m_bx, m_ay, m_by;
    index_t m_q;

    static void eval1D(T x, T a, T b, index_t q, T &val, T &d1, T &d2)
    {
        if (x <= a || x >= b)
        {
            val = T(0);
            d1  = T(0);
            d2  = T(0);
            return;
        }
        const T L = b - a;
        const T pi = T(3.1415926535897932384626433832795);
        const T theta = pi * (x - a) / L;
        const T s = std::sin(theta);
        const T c = std::cos(theta);
        const T s_qm2 = (q >= 2) ? std::pow(s, static_cast<int>(q - 2)) : T(1);
        const T s_qm1 = s_qm2 * s;
        const T s_q   = s_qm1 * s;

        val = s_q;
        d1  = T(q) * (pi / L) * s_qm1 * c;
        d2  = T(q) * (pi / L) * (pi / L) * s_qm2 * (T(q) * c * c - T(1));
    }

public:
    gsTrigBubbleFunction(T ax, T bx, T ay, T by, index_t q = 6)
        : m_ax(ax), m_bx(bx), m_ay(ay), m_by(by), m_q(std::max<index_t>(q, 4)) {}

    GISMO_CLONE_FUNCTION(gsTrigBubbleFunction)

    short_t domainDim() const override { return 2; }
    short_t targetDim() const override { return 1; }

    void eval_into(const gsMatrix<T> &u, gsMatrix<T> &result) const override
    {
        const index_t n = u.cols();
        result.resize(1, n);
        for (index_t k = 0; k < n; ++k)
        {
            T vx, dx1, dx2, vy, dy1, dy2;
            eval1D(u(0, k), m_ax, m_bx, m_q, vx, dx1, dx2);
            eval1D(u(1, k), m_ay, m_by, m_q, vy, dy1, dy2);
            result(0, k) = vx * vy;
        }
    }

    void deriv_into(const gsMatrix<T> &u, gsMatrix<T> &result) const override
    {
        const index_t n = u.cols();
        result.resize(2, n);
        for (index_t k = 0; k < n; ++k)
        {
            T vx, dx1, dx2, vy, dy1, dy2;
            eval1D(u(0, k), m_ax, m_bx, m_q, vx, dx1, dx2);
            eval1D(u(1, k), m_ay, m_by, m_q, vy, dy1, dy2);
            result(0, k) = dx1 * vy;
            result(1, k) = vx * dy1;
        }
    }

    void deriv2_into(const gsMatrix<T> &u, gsMatrix<T> &result) const override
    {
        const index_t n = u.cols();
        result.resize(3, n);
        for (index_t k = 0; k < n; ++k)
        {
            T vx, dx1, dx2, vy, dy1, dy2;
            eval1D(u(0, k), m_ax, m_bx, m_q, vx, dx1, dx2);
            eval1D(u(1, k), m_ay, m_by, m_q, vy, dy1, dy2);
            result(0, k) = dx2 * vy;
            result(1, k) = vx * dy2;
            result(2, k) = dx1 * dy1;
        }
    }
};

/**
 * @brief Computes a large axis-aligned rectangle [rect_xmin, rect_xmax] x [rect_ymin, rect_ymax]
 *        strictly inscribed inside the planar multi-patch domain Omega.
 *
 * Handles simply connected, non-convex, curved, and multiply connected (perforated with holes)
 * domains by chaining external/hole boundary curves into closed loops, performing multi-loop
 * even-odd ray-casting, rasterizing boundary segments to enforce clearance, and running the
 * O(N^2) maximal-rectangle-in-histogram search maximizing area * min(width, height).
 */
template <typename T>
void computeInscribedRectangle(const gsMultiPatch<T> &mp,
                               T &rect_xmin, T &rect_xmax,
                               T &rect_ymin, T &rect_ymax,
                               index_t gridN = 160)
{
    struct BdSeg { std::vector<std::pair<T, T>> pts; };
    std::vector<BdSeg> segs;
    const int K = 48;

    for (auto it = mp.bBegin(); it != mp.bEnd(); ++it)
    {
        BdSeg seg;
        const size_t pi = it->patch;
        const boxSide side = it->side();
        const short_t dir  = side.direction();
        const short_t tang = 1 - dir;
        const bool par     = side.parameter();

        const gsMatrix<T> sup = mp.patch(pi).support();
        const T fixedVal = sup(dir, par ? 1 : 0);
        const T t0 = sup(tang, 0), t1 = sup(tang, 1);

        gsMatrix<T> uv(2, K + 1);
        for (int k = 0; k <= K; ++k)
        {
            uv(dir,  k) = fixedVal;
            uv(tang, k) = t0 + (T(k) / T(K)) * (t1 - t0);
        }
        gsMatrix<T> xy;
        mp.patch(pi).eval_into(uv, xy);
        seg.pts.reserve(K + 1);
        for (int k = 0; k <= K; ++k)
            seg.pts.emplace_back(xy(0, k), xy(1, k));
        segs.push_back(give(seg));
    }

    gsMatrix<T> bbox;
    mp.boundingBox(bbox);
    const T xmin = bbox(0, 0), xmax = bbox(0, 1);
    const T ymin = bbox(1, 0), ymax = bbox(1, 1);
    const T Lx = (xmax - xmin > T(1e-12)) ? (xmax - xmin) : T(1);
    const T Ly = (ymax - ymin > T(1e-12)) ? (ymax - ymin) : T(1);
    const T diagSq = Lx * Lx + Ly * Ly;

    if (segs.empty())
    {
        rect_xmin = xmin; rect_xmax = xmax;
        rect_ymin = ymin; rect_ymax = ymax;
        return;
    }

    auto dsq = [](const std::pair<T, T> &a, const std::pair<T, T> &b) -> T {
        const T dx = a.first - b.first, dy = a.second - b.second;
        return dx * dx + dy * dy;
    };

    // Chain boundary segments into multiple closed loops (outer + any inner holes)
    std::vector<std::vector<std::pair<T, T>>> loops;
    std::vector<bool> used(segs.size(), false);
    const T closeTolSq = T(1e-10) * diagSq;
    const T matchTolSq = T(1e-4)  * diagSq;

    for (size_t start = 0; start < segs.size(); ++start)
    {
        if (used[start])
            continue;
        used[start] = true;
        std::vector<std::pair<T, T>> poly = segs[start].pts;

        while (true)
        {
            const auto &tail = poly.back();
            const auto &head = poly.front();
            const T dClose = dsq(tail, head);
            if (poly.size() > 2 && dClose <= closeTolSq)
                break;

            int best = -1;
            T bestD = T(1e30);
            bool rev = false;
            for (size_t k = 0; k < segs.size(); ++k)
            {
                if (used[k])
                    continue;
                const T d0 = dsq(tail, segs[k].pts.front());
                const T d1 = dsq(tail, segs[k].pts.back());
                if (d0 < bestD) { bestD = d0; best = static_cast<int>(k); rev = false; }
                if (d1 < bestD) { bestD = d1; best = static_cast<int>(k); rev = true; }
            }

            if (best < 0)
                break;
            if (bestD > matchTolSq && dClose <= bestD)
                break;

            used[best] = true;
            const auto &pts = segs[best].pts;
            if (rev)
            {
                for (int k = static_cast<int>(pts.size()) - 2; k >= 0; --k)
                    poly.push_back(pts[k]);
            }
            else
            {
                for (size_t k = 1; k < pts.size(); ++k)
                    poly.push_back(pts[k]);
            }
        }
        if (poly.size() >= 3)
            loops.push_back(give(poly));
    }

    auto inDomain = [&](T px, T py) -> bool {
        bool in = false;
        for (const auto &poly : loops)
        {
            for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
            {
                const T xi = poly[i].first,  yi = poly[i].second;
                const T xj = poly[j].first,  yj = poly[j].second;
                if (((yi > py) != (yj > py)) &&
                    (px < (xj - xi) * (py - yi) / (yj - yi) + xi))
                    in = !in;
            }
        }
        return in;
    };

    auto searchGrid = [&](index_t N, bool excludeBdrSegs,
                          T &outX0, T &outX1, T &outY0, T &outY1) -> bool {
        const T dx = Lx / T(N), dy = Ly / T(N);

        std::vector<std::vector<bool>> vIn(N + 1, std::vector<bool>(N + 1, false));
        for (index_t i = 0; i <= N; ++i)
            for (index_t j = 0; j <= N; ++j)
                vIn[i][j] = inDomain(xmin + T(i) * dx, ymin + T(j) * dy);

        std::vector<std::vector<bool>> cellValid(N, std::vector<bool>(N, false));
        for (index_t i = 0; i < N; ++i)
            for (index_t j = 0; j < N; ++j)
                cellValid[i][j] = vIn[i][j] && vIn[i + 1][j] &&
                                  vIn[i][j + 1] && vIn[i + 1][j + 1] &&
                                  inDomain(xmin + (T(i) + T(0.5)) * dx,
                                           ymin + (T(j) + T(0.5)) * dy);

        if (excludeBdrSegs)
        {
            for (const auto &poly : loops)
            {
                for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
                {
                    const T x0 = std::min(poly[i].first,  poly[j].first);
                    const T x1 = std::max(poly[i].first,  poly[j].first);
                    const T y0 = std::min(poly[i].second, poly[j].second);
                    const T y1 = std::max(poly[i].second, poly[j].second);
                    const index_t iMin = std::max<index_t>(0, static_cast<index_t>(std::floor((x0 - xmin) / dx)));
                    const index_t iMax = std::min<index_t>(N - 1, static_cast<index_t>(std::floor((x1 - xmin) / dx)));
                    const index_t jMin = std::max<index_t>(0, static_cast<index_t>(std::floor((y0 - ymin) / dy)));
                    const index_t jMax = std::min<index_t>(N - 1, static_cast<index_t>(std::floor((y1 - ymin) / dy)));
                    for (index_t ci = iMin; ci <= iMax; ++ci)
                        for (index_t cj = jMin; cj <= jMax; ++cj)
                            cellValid[ci][cj] = false;
                }
            }
        }

        // Exact O(N^2) maximal rectangle in histogram maximizing area * min(w, h)
        std::vector<index_t> h(N, 0);
        T bestScore = T(-1);
        index_t bestI0 = 0, bestI1 = 0, bestJ0 = 0, bestJ1 = 0;
        for (index_t j = 0; j < N; ++j)
        {
            for (index_t i = 0; i < N; ++i)
                h[i] = cellValid[i][j] ? (h[i] + 1) : 0;

            std::vector<index_t> st;
            for (index_t i = 0; i <= N; ++i)
            {
                const index_t curH = (i < N) ? h[i] : 0;
                while (!st.empty() && h[st.back()] >= curH)
                {
                    const index_t hh = h[st.back()];
                    st.pop_back();
                    const index_t i0 = st.empty() ? 0 : (st.back() + 1);
                    const index_t i1 = i;
                    if (hh > 0 && i1 > i0)
                    {
                        const T pw = T(i1 - i0) * dx;
                        const T ph = T(hh) * dy;
                        const T score = pw * ph * std::min(pw, ph);
                        if (score > bestScore)
                        {
                            bestScore = score;
                            bestI0 = i0;
                            bestI1 = i1;
                            bestJ0 = j + 1 - hh;
                            bestJ1 = j + 1;
                        }
                    }
                }
                st.push_back(i);
            }
        }

        if (bestScore <= T(0))
            return false;

        outX0 = xmin + T(bestI0) * dx;
        outX1 = xmin + T(bestI1) * dx;
        outY0 = ymin + T(bestJ0) * dy;
        outY1 = ymin + T(bestJ1) * dy;
        return true;
    };

    if (searchGrid(gridN, true, rect_xmin, rect_xmax, rect_ymin, rect_ymax))
        return;

    if (searchGrid(gridN * 2, false, rect_xmin, rect_xmax, rect_ymin, rect_ymax))
        return;

    // Fallback for extremely thin domains: find interior point farthest from boundary
    auto distSqToBoundary = [&](T px, T py) -> T {
        T minD2 = T(1e30);
        for (const auto &poly : loops)
        {
            for (size_t i = 0, j = poly.size() - 1; i < poly.size(); j = i++)
            {
                const T ax = poly[i].first, ay = poly[i].second;
                const T bx = poly[j].first, by = poly[j].second;
                const T vx = bx - ax, vy = by - ay;
                const T len2 = vx * vx + vy * vy;
                T t = (len2 > T(0)) ? ((px - ax) * vx + (py - ay) * vy) / len2 : T(0);
                t = std::max(T(0), std::min(T(1), t));
                const T cx = ax + t * vx - px;
                const T cy = ay + t * vy - py;
                minD2 = std::min(minD2, cx * cx + cy * cy);
            }
        }
        return minD2;
    };

    const index_t M = 200;
    const T dx = Lx / T(M), dy = Ly / T(M);
    T bestX = xmin + T(0.5) * Lx, bestY = ymin + T(0.5) * Ly, bestR2 = T(0);
    for (index_t i = 1; i < M; ++i)
    {
        const T px = xmin + T(i) * dx;
        for (index_t j = 1; j < M; ++j)
        {
            const T py = ymin + T(j) * dy;
            if (!inDomain(px, py))
                continue;
            const T r2 = distSqToBoundary(px, py);
            if (r2 > bestR2)
            {
                bestR2 = r2;
                bestX = px;
                bestY = py;
            }
        }
    }
    const T halfSide = T(0.65) * std::sqrt(bestR2);
    rect_xmin = bestX - halfSide;
    rect_xmax = bestX + halfSide;
    rect_ymin = bestY - halfSide;
    rect_ymax = bestY + halfSide;
}

template <typename T>
gsKnotVector<T> makeBubbleKnotVector(T lo, T hi, T a, T b, index_t p)
{
    std::vector<T> kn;
    kn.reserve(3 * p + 4);
    for (index_t i = 0; i <= p; ++i)
        kn.push_back(lo);
    for (index_t i = 0; i <= p + 1; ++i)
        kn.push_back(a * T(p + 1 - i) / T(p + 1) + b * T(i) / T(p + 1));
    for (index_t i = 0; i <= p; ++i)
        kn.push_back(hi);
    return gsKnotVector<T>(static_cast<short_t>(p), kn.begin(), kn.end());
}

template <typename T>
struct gsBubbleData
{
    T rect_xmin = 0, rect_xmax = 1, rect_ymin = 0, rect_ymax = 1;
    T a_x = 0, b_x = 1, a_y = 0, b_y = 1;
    gsTensorBSpline<2, T> bubble;
    gsPiecewiseFunction<T> field;
    gsPiecewiseFunction<T> lapField;
};

/**
 * @brief Constructs a normalized single tensor-product B-spline bubble of degree bubbleDeg
 *        supported on an interior rectangle [a_x, b_x] x [a_y, b_y] strictly inside mp.
 */
template <typename T>
gsBubbleData<T> makeBSplineBubble(const gsMultiPatch<T> &mp,
                                  index_t bubbleDeg,
                                  T margin = T(0.15))
{
    gsBubbleData<T> out;
    if (margin <= T(0) || margin >= T(0.5))
        margin = T(0.15);
    if (bubbleDeg < 2)
        bubbleDeg = 2;

    gsMatrix<T> bbox;
    mp.boundingBox(bbox);
    const T x_min = bbox(0, 0), x_max = bbox(0, 1);
    const T y_min = bbox(1, 0), y_max = bbox(1, 1);
    const T Lx = (x_max - x_min > T(1e-12)) ? (x_max - x_min) : T(1);
    const T Ly = (y_max - y_min > T(1e-12)) ? (y_max - y_min) : T(1);

    computeInscribedRectangle(mp, out.rect_xmin, out.rect_xmax,
                              out.rect_ymin, out.rect_ymax);

    const T iLx = out.rect_xmax - out.rect_xmin;
    const T iLy = out.rect_ymax - out.rect_ymin;
    out.a_x = out.rect_xmin + margin * iLx;
    out.b_x = out.rect_xmax - margin * iLx;
    out.a_y = out.rect_ymin + margin * iLy;
    out.b_y = out.rect_ymax - margin * iLy;

    const T pad_x = T(0.05) * Lx;
    const T pad_y = T(0.05) * Ly;

    gsKnotVector<T> kvx = makeBubbleKnotVector(x_min - pad_x, x_max + pad_x, out.a_x, out.b_x, bubbleDeg);
    gsKnotVector<T> kvy = makeBubbleKnotVector(y_min - pad_y, y_max + pad_y, out.a_y, out.b_y, bubbleDeg);
    gsTensorBSplineBasis<2, T> bubbleBasis(kvx, kvy);

    const index_t size0 = kvx.size() - bubbleDeg - 1;
    const index_t bumpIdx = (bubbleDeg + 1) + (bubbleDeg + 1) * size0;

    gsMatrix<T> bubbleCoefs(bubbleBasis.size(), 1);
    bubbleCoefs.setZero();
    bubbleCoefs(bumpIdx, 0) = T(1);

    out.bubble = gsTensorBSpline<2, T>(bubbleBasis, give(bubbleCoefs));

    // Normalize peak amplitude to 1.0 at the bubble center
    gsMatrix<T> center(2, 1);
    center(0, 0) = T(0.5) * (out.a_x + out.b_x);
    center(1, 0) = T(0.5) * (out.a_y + out.b_y);
    const T peak = out.bubble.eval(center)(0, 0);
    if (std::abs(peak) > T(1e-14))
        out.bubble.coefs() /= peak;

    gsLaplacianFunction<T> bubbleLap(out.bubble);
    for (size_t i = 0; i < mp.nPatches(); ++i)
    {
        out.field.addPiece(out.bubble);
        out.lapField.addPiece(bubbleLap);
    }
    return out;
}

/**
 * @brief Constructs a smooth trigonometric bubble of power trigPower supported on
 *        [a_x, b_x] x [a_y, b_y] strictly inside mp.
 */
template <typename T>
gsBubbleData<T> makeTrigBubble(const gsMultiPatch<T> &mp,
                               index_t trigPower,
                               T margin = T(0.15))
{
    gsBubbleData<T> out;
    if (margin <= T(0) || margin >= T(0.5))
        margin = T(0.15);

    computeInscribedRectangle(mp, out.rect_xmin, out.rect_xmax,
                              out.rect_ymin, out.rect_ymax);

    const T iLx = out.rect_xmax - out.rect_xmin;
    const T iLy = out.rect_ymax - out.rect_ymin;
    out.a_x = out.rect_xmin + margin * iLx;
    out.b_x = out.rect_xmax - margin * iLx;
    out.a_y = out.rect_ymin + margin * iLy;
    out.b_y = out.rect_ymax - margin * iLy;

    gsTrigBubbleFunction<T> trig(out.a_x, out.b_x, out.a_y, out.b_y, trigPower);
    gsLaplacianFunction<T> trigLap(trig);
    for (size_t i = 0; i < mp.nPatches(); ++i)
    {
        out.field.addPiece(trig);
        out.lapField.addPiece(trigLap);
    }
    return out;
}

struct gsBoxExactExprs
{
    std::string u_expr;
    std::string grad_x_expr;
    std::string grad_y_expr;
    std::string hess_xx;
    std::string hess_xy;
    std::string hess_yy;
    std::string lapl_expr;
    std::string rhs_expr;
};

/**
 * @brief Builds analytic expressions for box-conforming homogeneous-BC exact solutions
 *        on axis-aligned rectangular bounding boxes [x0, x1] x [y0, y1]:
 *        - First biharmonic:  u = sin^2(a*pi*(x-x0)/Lx) * sin^2(a*pi*(y-y0)/Ly)
 *        - Second biharmonic: u = sin(a*pi*(x-x0)/Lx)   * sin(a*pi*(y-y0)/Ly)
 */
template <typename T>
gsBoxExactExprs makeBoxHomogeneousExpressions(const gsMultiPatch<T> &mp,
                                              bool isSecondBiharmonic,
                                              T freqA,
                                              T beta)
{
    gsMatrix<T> bbox;
    mp.boundingBox(bbox);
    const T x0 = bbox(0, 0), x1 = bbox(0, 1);
    const T y0 = bbox(1, 0), y1 = bbox(1, 1);
    const T Lx = (x1 - x0 > T(1e-12)) ? (x1 - x0) : T(1);
    const T Ly = (y1 - y0 > T(1e-12)) ? (y1 - y0) : T(1);

    auto fmt = [](T v) -> std::string {
        std::ostringstream oss;
        oss << std::setprecision(16) << v;
        return oss.str();
    };

    const std::string sx0 = fmt(x0), sy0 = fmt(y0);
    const std::string sLx = fmt(Lx), sLy = fmt(Ly);
    const std::string sa  = fmt(freqA);
    const std::string sb  = fmt(beta);

    const std::string kx = "(" + sa + "*pi/" + sLx + ")";
    const std::string ky = "(" + sa + "*pi/" + sLy + ")";
    const std::string ax = "(" + kx + "*(x-(" + sx0 + ")))";
    const std::string ay = "(" + ky + "*(y-(" + sy0 + ")))";

    gsBoxExactExprs ex;
    if (isSecondBiharmonic)
    {
        const std::string uStr = "sin(" + ax + ")*sin(" + ay + ")";
        ex.u_expr      = uStr;
        ex.grad_x_expr = kx + "*cos(" + ax + ")*sin(" + ay + ")";
        ex.grad_y_expr = ky + "*sin(" + ax + ")*cos(" + ay + ")";
        ex.hess_xx     = "-(" + kx + ")^2*" + uStr;
        ex.hess_xy     = kx + "*" + ky + "*cos(" + ax + ")*cos(" + ay + ")";
        ex.hess_yy     = "-(" + ky + ")^2*" + uStr;
        ex.lapl_expr   = "0";
        ex.rhs_expr    = "(" + sb + "+((" + kx + ")^2+(" + ky + ")^2)^2)*" + uStr;
    }
    else
    {
        const std::string sx2 = "(sin(" + ax + ")^2)";
        const std::string sy2 = "(sin(" + ay + ")^2)";
        const std::string uStr = sx2 + "*" + sy2;
        const std::string s2x = "sin(2*" + ax + ")";
        const std::string s2y = "sin(2*" + ay + ")";
        const std::string c2x = "cos(2*" + ax + ")";
        const std::string c2y = "cos(2*" + ay + ")";

        ex.u_expr      = uStr;
        ex.grad_x_expr = kx + "*" + s2x + "*" + sy2;
        ex.grad_y_expr = ky + "*" + sx2 + "*" + s2y;
        ex.hess_xx     = "2*(" + kx + ")^2*" + c2x + "*" + sy2;
        ex.hess_xy     = kx + "*" + ky + "*" + s2x + "*" + s2y;
        ex.hess_yy     = "2*(" + ky + ")^2*" + sx2 + "*" + c2y;
        ex.lapl_expr   = "2*(" + kx + ")^2*" + c2x + "*" + sy2 + "+2*(" + ky + ")^2*" + sx2 + "*" + c2y;
        ex.rhs_expr    = sb + "*" + uStr
                       + "-8*(" + kx + ")^4*" + c2x + "*" + sy2
                       + "+8*(" + kx + ")^2*(" + ky + ")^2*" + c2x + "*" + c2y
                       + "-8*(" + ky + ")^4*" + sx2 + "*" + c2y;
    }
    return ex;
}

} // namespace gismo
