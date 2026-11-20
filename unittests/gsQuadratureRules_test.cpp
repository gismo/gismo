/** @file gsQuadratureRules_test.cpp

    @brief Tests available gsQuadratureRules (containing gsGaussRules 1-21)

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Vogl, A. Mantzaflaris
 **/

#include "gismo_unittest.h"       // Brings in G+Smo and the UnitTest++ framework
#include <gsHSplines/gsHElement.h>
#include <gsHSplines/gsHElementHelper.h>
#include <random>

SUITE(gsQuadratureRules_test)                 // The suite should have the same name as the file
{

void testWork(const index_t nodes[], size_t dim);

real_t calcAntiDerivative(gsVector<index_t> const &deg, const index_t dim);

real_t calcPoly(gsVector<index_t> const &deg,
                gsQuadRule<real_t> const &gr,
                const index_t dim);

const char *addPlus(const index_t d, index_t i);

TEST(tensor_quad_3)
{
    index_t array[] = {3};
    testWork(array, 1);
}

TEST(tensor_quad_5)
{
    index_t array[] = {5};
    testWork(array, 1);
}

TEST(tensor_quad_1_5_9)
{
    index_t array[] = {1, 5, 9};
    testWork(array, 3);
}

TEST(tensor_quad_2_6_7)
{
    index_t array[] = {2, 6, 7};
    testWork(array, 3);
}

TEST(tensor_quad_3_4_8)
{
    index_t array[] = {3, 4, 8};
    testWork(array, 3);
}

TEST(tensor_quad_10_17)
{
    index_t array[] = {10, 17};
    testWork(array, 2);
}

TEST(tensor_quad_11_16)
{
    index_t array[] = {11, 16};
    testWork(array, 2);
}

TEST(tensor_quad_12_15)
{
    index_t array[] = {12, 15};
    testWork(array, 2);
}

TEST(tensor_quad_13_14)
{
    index_t array[] = {13, 14};
    testWork(array, 2);
}

TEST(tensor_quad_18)
{
    index_t array[] = {18};
    testWork(array, 1);
}

TEST(tensor_quad_19)
{
    index_t array[] = {19};
    testWork(array, 1);
}

TEST(tensor_quad_20)
{
    index_t array[] = {20};
    testWork(array, 1);
}

TEST(tensor_quad_21)
{
    index_t array[] = {21};
    testWork(array, 1);
}

TEST(tensor_quad_22)
{
    index_t array[] = {22};
    testWork(array, 1);
}

void testWork(const index_t nodes[], const size_t dim)
{
    gsVector<index_t> numNodes = gsAsConstVector<index_t>(nodes, dim);

    // ---------------- Lobatto

    gsGaussRule<real_t> legendreRule(numNodes);
    gsGaussRule<real_t> legendreRuleComp(numNodes, REAL_DIG);
    gsQuadRule<real_t> quadLeg_lookup = gsQuadrature::get<real_t>(gsQuadrature::GaussLegendre, numNodes, 0);
    gsQuadRule<real_t> quadLeg_compute = gsQuadrature::get<real_t>(gsQuadrature::GaussLegendre, numNodes, REAL_DIG);
    CHECK(legendreRule.referenceNodes() == quadLeg_lookup.referenceNodes());
    CHECK_MATRIX_CLOSE(quadLeg_compute.referenceNodes(), quadLeg_lookup.referenceNodes(), EPSILON);
    CHECK_MATRIX_CLOSE(quadLeg_compute.referenceWeights(), quadLeg_lookup.referenceWeights(), EPSILON);

    gsVector<index_t> legVec = (2 * numNodes - 1 * gsVector<index_t>::Ones(dim)).cwiseMax(gsVector<index_t>::Zero(dim));
    real_t expectedLeg = calcAntiDerivative(legVec, dim);
    real_t lookupLeg = calcPoly(legVec, quadLeg_lookup, dim);
    real_t computeLeg = calcPoly(legVec, quadLeg_compute, dim);
    CHECK_CLOSE(expectedLeg, lookupLeg, EPSILON);
    CHECK_CLOSE(expectedLeg, computeLeg, EPSILON);

    // ---------------- Lobatto

    gsLobattoRule<real_t> lobattoRule(numNodes);
    gsLobattoRule<real_t> lobattoRuleComp(numNodes, REAL_DIG);
    gsQuadRule<real_t> quadLob_lookup = gsQuadrature::get<real_t>(gsQuadrature::GaussLobatto, numNodes, 0);
    gsQuadRule<real_t> quadLob_compute = gsQuadrature::get<real_t>(gsQuadrature::GaussLobatto, numNodes, REAL_DIG);
    CHECK(lobattoRule.referenceNodes() == quadLob_lookup.referenceNodes());
    CHECK_MATRIX_CLOSE(quadLob_lookup.referenceNodes(), quadLob_compute.referenceNodes(), EPSILON);
    CHECK_MATRIX_CLOSE(quadLob_compute.referenceWeights(), quadLob_lookup.referenceWeights(), EPSILON);

    gsVector<index_t> lobVec = (2 * numNodes - 3 * gsVector<index_t>::Ones(dim)).cwiseMax(gsVector<index_t>::Zero(dim));
    real_t expectedLob = calcAntiDerivative(lobVec, dim);
    real_t lookupLob = calcPoly(lobVec, quadLob_lookup, dim);
    real_t computeLob = calcPoly(lobVec, quadLob_compute, dim);
    CHECK_CLOSE(expectedLob, lookupLob, EPSILON);
    CHECK_CLOSE(expectedLob, computeLob, EPSILON);
}

real_t calcAntiDerivative(gsVector<index_t> const &deg, const index_t dim)
{
    // Test index_tegration
    gsVector<real_t> u;
    u.setConstant(dim, 1.0123);

    // Construct polynomial
    std::string var = std::string("xyzwuv").substr(0, dim);    // defining variable names
    std::stringstream poly;
    for (index_t i = 0; i < dim; ++i)
    {
        // Construct the anti-derivative
        std::string tmp(var);
        tmp.erase(i, 1); // cut out current variable, that will be to the power of deg[i]
        for (index_t j = dim - 1; j; --j) tmp.insert(j, "*");   // add * between variables
        poly << "(1.0/" << deg[i] + 1 << ")*" << tmp << var[i] << "^" << deg[i] + 1 << addPlus(dim, i);
    }

    gsFunctionExpr<real_t> antideriv(poly.str(), dim);
    return antideriv.eval(u).at(0);
}

real_t calcPoly(gsVector<index_t> const &deg,
                gsQuadRule<real_t> const &gr,
                const index_t dim)
{
    // Test index_tegration
    gsVector<real_t> u;
    u.setConstant(dim, 1.0123);

    const index_t d = gr.dim();
    CHECK_EQUAL(d, dim);

    const gsVector<real_t> l = gsVector<real_t>::Zero(d);

    gsMatrix<real_t> ngrid;
    gsVector<real_t> wgrid;

    // Map rule to index_tegration domain
    gr.mapTo(l, u, ngrid, wgrid);

    // Construct polynomial
    std::string var = std::string("xyzwuv").substr(0, d);    // defining variable names
    std::stringstream poly;
    for (index_t i = 0; i < d; ++i)
    {
        // Make a polynomial of requested degree
        poly << var[i] << "^" << deg[i] << addPlus(d, i);
    }

    gsFunctionExpr<real_t> index_tegrand(poly.str(), d);

    ngrid = index_tegrand.eval(ngrid);
    return wgrid.dot(ngrid.row(0));
}

const char *addPlus(const index_t d, index_t i)
{ return (i == d - 1 ? "" : "+"); }

// Integrate f element-by-element over domain using gsPatchRule
real_t patchIntegrate(const gsDomain<real_t>& domain,
                      gsPatchRule<real_t>& rule,
                      const gsFunctionExpr<real_t>& f)
{
    gsMatrix<real_t> nodes;
    gsVector<real_t> weights;
    real_t result = 0;
    auto it    = domain.beginAll();
    auto itEnd = domain.endAll();
    for (; it < itEnd; ++it)
    {
        rule.mapTo(it.lowerCorner(), it.upperCorner(), nodes, weights);
        if (weights.size() == 0) continue;
        gsMatrix<real_t> vals = f.eval(nodes);
        result += weights.dot(vals.row(0).transpose());
    }
    return result;
}

// Integrate f element-by-element over domain using a Gauss-Legendre rule
real_t gaussIntegrate(const gsDomain<real_t>& domain,
                      gsGaussRule<real_t>& rule,
                      const gsFunctionExpr<real_t>& f)
{
    gsMatrix<real_t> nodes;
    gsVector<real_t> weights;
    real_t result = 0;
    auto it    = domain.beginAll();
    auto itEnd = domain.endAll();
    for (; it < itEnd; ++it)
    {
        rule.mapTo(it.lowerCorner(), it.upperCorner(), nodes, weights);
        gsMatrix<real_t> vals = f.eval(nodes);
        result += weights.dot(vals.row(0).transpose());
    }
    return result;
}

TEST(patchrule_1d_weight_sum)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1); // degree 3, 4 elements
    gsTensorBSplineBasis<1,real_t> basis(kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("1", 1);
    CHECK_CLOSE(1.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_1d_polynomial_exactness)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1);
    gsTensorBSplineBasis<1,real_t> basis(kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("x^3", 1);
    CHECK_CLOSE(0.25, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_2d_weight_sum)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("1", 2);
    CHECK_CLOSE(1.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_2d_polynomial_exactness)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("x^2*y^2", 2);
    // integral of x^2*y^2 over [0,1]^2 = (1/3)*(1/3) = 1/9
    CHECK_CLOSE(1.0/9.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_3d_weight_sum)
{
    gsKnotVector<real_t> kv(0, 1, 2, 3, 1); // degree 3, 3 elements
    gsTensorBSplineBasis<3,real_t> basis(kv, kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("1", 3);
    CHECK_CLOSE(1.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_3d_polynomial_exactness)
{
    gsKnotVector<real_t> kv(0, 1, 2, 3, 1);
    gsTensorBSplineBasis<3,real_t> basis(kv, kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);
    gsFunctionExpr<real_t> f("x^2*y^2*z^2", 3);
    // integral of x^2*y^2*z^2 over [0,1]^3 = (1/3)^3 = 1/27
    CHECK_CLOSE(1.0/27.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_overintegration_2d)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, true);
    gsFunctionExpr<real_t> f("x^3*y^3", 2);
    // integral of x^3*y^3 over [0,1]^2 = (1/4)*(1/4) = 1/16
    CHECK_CLOSE(1.0/16.0, patchIntegrate(*basis.domain(), rule, f), EPSILON);
}

TEST(patchrule_matches_gauss_on_polynomial)
{
    gsKnotVector<real_t> kv(0, 1, 3, 3, 1);
    gsTensorBSplineBasis<2,real_t> basis(kv, kv);
    gsPatchRule<real_t> rule(*basis.domain(), 3, 2, false);

    // Gauss rule with 4 pts/direction is exact for degree 7, well above degree 3
    gsVector<index_t> numPts(2); numPts << 4, 4;
    gsGaussRule<real_t> gauss(numPts);

    gsFunctionExpr<real_t> f("x^3*y", 2);
    // analytic: (1/4)*(1/2) = 1/8
    real_t patchResult = patchIntegrate(*basis.domain(), rule, f);
    real_t gaussResult = gaussIntegrate(*basis.domain(), gauss, f);
    CHECK_CLOSE(1.0/8.0, patchResult, EPSILON);
    CHECK_CLOSE(gaussResult, patchResult, EPSILON);
}

// ============================================================
// gsHPatchRule tests
// ============================================================

// Build a 2D truncated THB with a corner refinement (creates truncation).
gsTHBSplineBasis<2,real_t> makeRefinedTHB(index_t p)
{
    gsKnotVector<real_t> kv(0, 1, 0, p + 1, 1);
    kv.uniformRefine();
    kv.uniformRefine();
    gsTensorBSplineBasis<2,real_t> tb(kv, kv);
    gsTHBSplineBasis<2,real_t> thb(tb);
    thb.uniformRefine();                                  // level-1 everywhere
    // Refine the lower-left corner to level 2 (4x4 cells in level-2 indices)
    std::vector<index_t> boxes = {2, 0,0, 4,4};
    thb.refineElements(boxes);
    return thb;
}

bool hpatchIndexBoxesOverlap(const gsMatrix<index_t,2,2> & a,
                             const gsMatrix<index_t,2,2> & b)
{
    for (short_t k = 0; k != 2; ++k)
        if (a(k,0) >= b(k,1) || a(k,1) <= b(k,0))
            return false;
    return true;
}

index_t hpatchMaxTruncatedRepLevelOnElement(const gsTHBSplineBasis<2,real_t> & basis,
                                            const gsVector<real_t> & lower,
                                            const gsVector<real_t> & upper)
{
    gsHElementHelper<2,real_t> helper(basis);
    gsVector<real_t,2> low  = lower;
    gsVector<real_t,2> upp  = upper;
    index_t targetLevel = 0;

    for (std::map<index_t, gsSparseVector<real_t> >::const_iterator it =
             basis.truncatedBegin(); it != basis.truncatedEnd(); ++it)
    {
        const index_t funIndex = it->first;
        const index_t repLevel = basis.repLevel(funIndex);
        if (repLevel <= targetLevel)
            continue;

        const gsHElement<2,real_t> repElement =
            helper.toElement(low, upp, repLevel);
        gsMatrix<index_t,2,2> repElementBox;
        repElementBox.col(0) = repElement.lowerCorner();
        repElementBox.col(1) = repElement.upperCorner();

        const gsTensorBSplineBasis<2,real_t> & repBasis =
            basis.tensorLevel(repLevel);
        const gsSparseVector<real_t> & coefs = it->second;

        gsMatrix<index_t,2,2> supp;
        for (gsSparseVector<real_t>::InnerIterator coefIt(coefs);
             coefIt; ++coefIt)
        {
            if (coefIt.value() == 0)
                continue;

            repBasis.elementSupport_into(coefIt.index(), supp);
            if (hpatchIndexBoxesOverlap(supp, repElementBox))
            {
                targetLevel = repLevel;
                break;
            }
        }
    }

    return targetLevel;
}

gsTHBSplineBasis<2,real_t>
hpatchSupportWideIntegrationBasis(const gsTHBSplineBasis<2,real_t> & basis)
{
    gsTHBSplineBasis<2,real_t> result(basis);
    gsHElementHelper<2,real_t> helper(result);

    std::vector<index_t> refboxes;
    for (index_t i = 0; i != result.size(); ++i)
    {
        if (!result.isTruncated(i)) continue;

        const index_t repLevel = result.repLevel(i);
        const gsMatrix<real_t> support = result.support(i);
        const gsVector<real_t,2> low = support.col(0);
        const gsVector<real_t,2> upp = support.col(1);

        const gsHElement<2,real_t> element =
            helper.toElement(low, upp, repLevel);
        const gsVector<index_t,2> & lowCorner = element.lowerCorner();
        const gsVector<index_t,2> & uppCorner = element.upperCorner();

        refboxes.push_back(repLevel);
        for (short_t k = 0; k != 2; ++k) refboxes.push_back(lowCorner(k));
        for (short_t k = 0; k != 2; ++k) refboxes.push_back(uppCorner(k));
    }

    if (!refboxes.empty())
        result.refineElements(refboxes);
    return result;
}

// Integrate a scalar function over the THB domain using a quadrature rule.
real_t hpatchIntegrate(const gsHTensorBasis<2,real_t> & basis,
                       gsQuadRule<real_t> & rule,
                       const gsFunctionExpr<real_t> & f)
{
    gsMatrix<real_t> nodes; gsVector<real_t> weights;
    real_t result = 0;
    for (auto & elem : basis.domain()->allElements())
    {
        rule.mapTo(elem.lowerCorner(), elem.upperCorner(), nodes, weights);
        if (weights.size() == 0) continue;
        gsMatrix<real_t> vals = f.eval(nodes);
        result += weights.dot(vals.row(0).transpose());
    }
    return result;
}

// Integrals of all functions of `target`, accumulated over the elements of `domain`
// with `rule`; entry i is the integral of basis function i. A rule is exact on the
// span of `target` exactly when this vector equals the one of a rule known to be
// exact there.
gsVector<real_t> hpatchIntegrateBasis(const gsBasis<real_t> & target,
                                      const gsDomain<real_t> & domain,
                                      const gsQuadRule<real_t> & rule)
{
    gsVector<real_t> result = gsVector<real_t>::Zero(target.size());
    gsMatrix<real_t> nodes, vals;
    gsVector<real_t> weights;
    gsMatrix<index_t> actives;
    for (auto & elem : domain.allElements())
    {
        rule.mapTo(elem.lowerCorner(), elem.upperCorner(), nodes, weights);
        if (weights.size() == 0) continue;
        target.active_into(nodes, actives);
        target.eval_into(nodes, vals);
        // Per column: a node on a knot may have other active functions than the
        // interior nodes of its element.
        for (index_t k = 0; k != nodes.cols(); ++k)
            for (index_t a = 0; a != actives.rows(); ++a)
                result(actives(a,k)) += weights[k] * vals(a,k);
    }
    return result;
}

// Mass matrix sum_k w_k N_i(x_k) N_j(x_k) of `basis` on the boundary side `side`, with the
// nodes of `rule` on the side's elements. Evaluated node by node, since a node on a leaf
// boundary can have other active functions than the element's interior nodes.
gsMatrix<real_t> hpatchBoundaryMass(const gsTHBSplineBasis<2,real_t> & basis,
                                    boxSide side, const gsQuadRule<real_t> & rule)
{
    gsMatrix<real_t> M = gsMatrix<real_t>::Zero(basis.size(), basis.size());
    gsMatrix<real_t> nodes, pt, val;
    gsVector<real_t> weights;
    gsMatrix<index_t> act;
    auto dom = basis.domain();
    gsBasis<real_t>::domainIter it = dom->beginBdr(side), end = dom->endBdr(side);
    for (; it < end; ++it)
    {
        rule.mapTo(it.lowerCorner(), it.upperCorner(), nodes, weights);
        for (index_t k = 0; k != nodes.cols(); ++k)
        {
            pt = nodes.col(k);
            basis.active_into(pt, act);
            basis.eval_into(pt, val);
            for (index_t a = 0; a != act.rows(); ++a)
                for (index_t b = 0; b != act.rows(); ++b)
                    M(act(a,0), act(b,0)) += weights[k] * val(a,0) * val(b,0);
        }
    }
    return M;
}

// Assemble a form (0=mass, 1=stiffness) on an identity geometry with the rule
// selected by quRule, and return the dense matrix.
gsMatrix<real_t> hpatchAssemble(gsTHBSplineBasis<2,real_t> thb, int form,
                                index_t quRule, index_t quA, index_t quB)
{
    gsTensorBSplineBasis<2,real_t> & l0 = thb.tensorLevel(0);
    gsMatrix<real_t> coefs = l0.anchors().transpose();
    gsTensorBSpline<2,real_t> geom(l0, coefs);
    gsMultiPatch<real_t> mp; mp.addPatch(geom);
    gsMultiBasis<real_t> mb; mb.addBasis(thb.clone().release()); mb.setTopology(mp);

    gsExprAssembler<real_t> A(1,1);
    gsOptionList opt = gsExprAssembler<real_t>::defaultOptions();
    opt.setInt("quRule", quRule);
    opt.setReal("quA", static_cast<real_t>(quA));
    opt.setInt("quB", quB);
    A.setOptions(opt);
    A.setIntegrationElements(mb);
    auto G = A.getMap(mp);
    auto u = A.getSpace(mb);
    A.initSystem();
    if (form == 0) A.assemble( u * u.tr() * meas(G) );
    else           A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G) );
    return A.matrix().toDense();
}

gsMatrix<real_t> hpatchAssembleOnIntegrationBasis(
                                const gsTHBSplineBasis<2,real_t> & space,
                                const gsTHBSplineBasis<2,real_t> & integration,
                                int form, index_t quRule,
                                index_t quA, index_t quB)
{
    gsTensorBSplineBasis<2,real_t> & l0 = space.tensorLevel(0);
    gsMatrix<real_t> coefs = l0.anchors().transpose();
    gsTensorBSpline<2,real_t> geom(l0, coefs);
    gsMultiPatch<real_t> mp; mp.addPatch(geom);
    gsMultiBasis<real_t> spaceMB;
    spaceMB.addBasis(space.clone().release());
    spaceMB.setTopology(mp);
    gsMultiBasis<real_t> integMB;
    integMB.addBasis(integration.clone().release());
    integMB.setTopology(mp);

    gsExprAssembler<real_t> A(1,1);
    gsOptionList opt = gsExprAssembler<real_t>::defaultOptions();
    opt.setInt("quRule", quRule);
    opt.setReal("quA", static_cast<real_t>(quA));
    opt.setInt("quB", quB);
    A.setOptions(opt);
    A.setIntegrationElements(integMB);
    auto G = A.getMap(mp);
    auto u = A.getSpace(spaceMB);
    A.initSystem();
    if (form == 0) A.assemble( u * u.tr() * meas(G) );
    else           A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G) );
    return A.matrix().toDense();
}

TEST(hpatchrule_2d_partition_of_unity)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, p-1, false);
    gsFunctionExpr<real_t> one("1", 2);
    // integral of 1 over [0,1]^2 == area == 1
    CHECK_CLOSE(1.0, hpatchIntegrate(thb, rule, one), 1e-12);
}

TEST(hpatchrule_2d_polynomial_exactness)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, p-1, false);
    gsFunctionExpr<real_t> f("x^2*y^2", 2);
    // integral over [0,1]^2 = (1/3)*(1/3) = 1/9
    CHECK_CLOSE(1.0/9.0, hpatchIntegrate(thb, rule, f), 1e-12);
}

TEST(hpatchrule_2d_mass_vs_gauss)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    gsMatrix<real_t> Mg = hpatchAssemble(thb, 0, gsQuadrature::GaussLegendre, 2*p, p-1);
    gsMatrix<real_t> Mh = hpatchAssemble(thb, 0, gsQuadrature::HPatchRule,    2*p, p-1);
    CHECK_CLOSE(0.0, (Mg - Mh).cwiseAbs().maxCoeff(), 1e-11);
}

TEST(hpatchrule_2d_stiffness_vs_gauss)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    gsMatrix<real_t> Kg = hpatchAssemble(thb, 1, gsQuadrature::GaussLegendre, 2*p, p-2);
    gsMatrix<real_t> Kh = hpatchAssemble(thb, 1, gsQuadrature::HPatchRule,    2*p, p-2);
    CHECK_CLOSE(0.0, (Kg - Kh).cwiseAbs().maxCoeff(), 1e-11);
}

TEST(hpatchrule_integration_basis_respects_representation_levels)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    gsTHBSplineBasis<2,real_t> integ =
        gsHPatchRule<real_t>::makeIntegrationBasis<2>(thb);

    for (auto & elem : integ.domain()->allElements())
    {
        const index_t leafLevel =
            static_cast<const gsHDomainIterator<real_t,2,index_t> *>(&elem)->getLevel();
        const index_t requiredLevel =
            hpatchMaxTruncatedRepLevelOnElement(thb,
                                                elem.lowerCorner(),
                                                elem.upperCorner());
        CHECK(leafLevel >= requiredLevel);
    }
}

TEST(hpatchrule_integration_basis_not_support_wide)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    gsTHBSplineBasis<2,real_t> integ =
        gsHPatchRule<real_t>::makeIntegrationBasis<2>(thb);
    gsTHBSplineBasis<2,real_t> supportWide =
        hpatchSupportWideIntegrationBasis(thb);

    gsHPatchRule<real_t> ruleNew(*integ.domain(), 2*p, p-1, false);
    gsHPatchRule<real_t> ruleOld(*supportWide.domain(), 2*p, p-1, false);

    CHECK(integ.numElements() <= supportWide.numElements());
    CHECK(ruleNew.numPointsTotal() <= ruleOld.numPointsTotal());
}

TEST(hpatchrule_integration_basis_assembly_vs_gauss)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    gsTHBSplineBasis<2,real_t> integ =
        gsHPatchRule<real_t>::makeIntegrationBasis<2>(thb);

    gsMatrix<real_t> Mg = hpatchAssemble(thb, 0, gsQuadrature::GaussLegendre, 2*p, p-1);
    gsMatrix<real_t> Mh = hpatchAssembleOnIntegrationBasis(
        thb, integ, 0, gsQuadrature::HPatchRule, 2*p, p-1);
    CHECK_CLOSE(0.0, (Mg - Mh).cwiseAbs().maxCoeff(), 1e-11);

    gsMatrix<real_t> Kg = hpatchAssemble(thb, 1, gsQuadrature::GaussLegendre, 2*p, p-2);
    gsMatrix<real_t> Kh = hpatchAssembleOnIntegrationBasis(
        thb, integ, 1, gsQuadrature::HPatchRule, 2*p, p-2);
    CHECK_CLOSE(0.0, (Kg - Kh).cwiseAbs().maxCoeff(), 1e-11);
}

TEST(hpatchrule_finest_fallback_on_hierarchical)
{
    // gsPatchRule (quRule=3) on a hierarchical basis must default to the finest
    // level and still assemble a correct mass matrix (vs Gauss).
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    gsMatrix<real_t> Mg = hpatchAssemble(thb, 0, gsQuadrature::GaussLegendre, 2*p, p-1);
    gsMatrix<real_t> Mp = hpatchAssemble(thb, 0, gsQuadrature::PatchRule,     2*p, p-1);
    CHECK_CLOSE(0.0, (Mg - Mp).cwiseAbs().maxCoeff(), 1e-11);
}

TEST(hpatchrule_2d_fewer_points_than_gauss)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, p-1, false);

    // Gauss reference count: (p+1)^2 per element
    index_t nGauss = 0;
    gsMatrix<real_t> nodes; gsVector<real_t> weights;
    gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(2, p+1));
    for (auto & elem : thb.domain()->allElements())
    {
        gauss.mapTo(elem.lowerCorner(), elem.upperCorner(), nodes, weights);
        nGauss += weights.size();
    }
    CHECK(rule.numPointsTotal() < nGauss);
}

// Level-0 tensor B-spline basis of degree p on [0,1]^2 with a knot of multiplicity
// `mult` at `pos` in direction 0 (simple knots in direction 1), refined the way
// makeRefinedTHB does: uniform to level 1, then the given boxes.
gsTHBSplineBasis<2,real_t> makeRepeatedKnotTHB(index_t p, index_t mult, real_t pos,
                                               const std::vector<index_t> & boxes)
{
    gsKnotVector<real_t> kv0(0, 1, 0, p + 1, 1);
    kv0.uniformRefine();
    kv0.uniformRefine();          // breakpoints 0, .25, .5, .75, 1

    gsKnotVector<real_t> kv1 = kv0;
    if (mult > 1) kv0.insert(pos, mult - 1);   // direction 0 carries the repeated knot

    gsTensorBSplineBasis<2,real_t> tb(kv0, kv1);
    gsTHBSplineBasis<2,real_t> thb(tb);
    thb.uniformRefine();                                  // level-1 everywhere
    std::vector<index_t> b = boxes;
    thb.refineElements(b);
    return thb;
}

// The target knot vector currently built by gsHPatchRule::_build
// (src/gsAssembler/gsHPatchRule.hpp:124-129): unique breakpoints of the leaf,
// uniform end multiplicity deg+1 and interior multiplicity deg-reg.
gsKnotVector<real_t> hpatchTargetKvCurrent(const gsKnotVector<real_t> & kvLev,
                                           index_t lowerIdx, index_t upperIdx,
                                           index_t deg, index_t reg)
{
    std::vector<real_t> uKnots;
    for (index_t i = lowerIdx; i <= upperIdx; ++i)
        uKnots.push_back(kvLev(i));
    return gsKnotVector<real_t>(uKnots, static_cast<int>(deg), static_cast<int>(reg));
}

// The per-breakpoint construction Task 03 will install: interior breakpoint j gets
// multiplicity min(deg - reg - 1 + mu_j, deg + 1), with mu_j its multiplicity in the
// level knot vector; end breakpoints get deg + 1 (open knot vector on the leaf).
gsKnotVector<real_t> hpatchTargetKvPerBreakpoint(const gsKnotVector<real_t> & kvLev,
                                                 index_t lowerIdx, index_t upperIdx,
                                                 index_t deg, index_t reg)
{
    std::vector<real_t> repKnots;
    for (index_t i = lowerIdx; i <= upperIdx; ++i)
    {
        const real_t z    = kvLev(i);
        const bool   endK = (i == lowerIdx || i == upperIdx);
        index_t mult;
        if (endK)
            mult = deg + 1;
        else
        {
            const index_t muj = kvLev.multiplicity(z);
            mult = math::min(deg - reg - 1 + muj, deg + 1);
        }
        repKnots.insert(repKnots.end(), mult, z);
    }
    return gsKnotVector<real_t>(repKnots, static_cast<short_t>(deg));
}

TEST(hpatchrule_target_kv_simple_knot_characterization)
{
    const index_t ps[] = {2,3,4};
    const index_t ns[] = {1,2,3,4,8};

    for (index_t p : ps)
    {
        const index_t deg = 2*p;
        const index_t regs[] = {p-1, p-2};
        for (index_t reg : regs)
        {
            for (index_t n : ns)
            {
                gsKnotVector<real_t> kvLev(0, 1, n - 1, p + 1, 1);
                const index_t lowerIdx = 0;
                const index_t upperIdx = n;

                gsKnotVector<real_t> kvCur =
                    hpatchTargetKvCurrent(kvLev, lowerIdx, upperIdx, deg, reg);
                gsKnotVector<real_t> kvNew =
                    hpatchTargetKvPerBreakpoint(kvLev, lowerIdx, upperIdx, deg, reg);

                CHECK_EQUAL(kvCur.size(), kvNew.size());

                gsKnotVector<real_t>::const_iterator itCur = kvCur.begin();
                gsKnotVector<real_t>::const_iterator itNew = kvNew.begin();
                for ( ; itCur != kvCur.end() && itNew != kvNew.end(); ++itCur, ++itNew)
                    CHECK_CLOSE(*itCur, *itNew, 1e-14);

                CHECK_EQUAL(static_cast<index_t>(deg), static_cast<index_t>(kvNew.degree()));
                CHECK_EQUAL(static_cast<index_t>(deg), static_cast<index_t>(kvCur.degree()));

                CHECK_EQUAL(gsBSplineBasis<real_t>(kvCur).size(),
                            gsBSplineBasis<real_t>(kvNew).size());
            }
        }
    }
}

TEST(hpatchrule_simple_knot_rule_baseline)
{
    // Golden values of the real gsHPatchRule::_build on simple-knot meshes. Every
    // interior breakpoint has multiplicity 1 there, so any per-breakpoint refinement of
    // the target space must leave these counts bit-identical. No row has a fallback
    // cell: every direction solves its 1D moment equations, except direction 0 of the
    // one-element-wide level-0 leaf, where the patch rule ties with (p+1)-point Gauss
    // and Gauss is used uncounted; p = 4, reg = 2 needs the damped
    // Levenberg-Marquardt phase of _compute1D, since the plain Newton step diverges.
    struct Row { index_t p; index_t reg; index_t nPts; index_t nFallback; };
    const Row rows[] = { {2,1,241,0}, {2,0,381,0}, {3,2,452,0}, {3,1,596,0}, {4,3,684,0}, {4,2,909,0} };

    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(rows[k].p);
        gsHPatchRule<real_t> rule(*thb.domain(), 2*rows[k].p, rows[k].reg, false);
        CHECK_EQUAL(rows[k].nPts,      rule.numPointsTotal());
        CHECK_EQUAL(rows[k].nFallback, rule.numFallbackCells());
    }
}

TEST(hpatchrule_2d_mass_vs_gauss_C0_knot)
{
    const index_t p = 3;
    // C^0 line (knot multiplicity p) at x = 0.5 in direction 0; simple knots in direction 1.
    gsTHBSplineBasis<2,real_t> thb = makeRepeatedKnotTHB(p, 3, 0.5, {2,0,0,4,4});

    // The rule must be the Johannessen rule on every cell: a Gauss-fallback cell is
    // integrated exactly by construction, so agreement with Gauss there would verify
    // nothing about the target space this test exists to probe.
    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, p-1, false);
    CHECK_EQUAL(index_t(0), rule.numFallbackCells());

    gsMatrix<real_t> Mg = hpatchAssemble(thb, 0, gsQuadrature::GaussLegendre, 2*p, p-1);
    gsMatrix<real_t> Mh = hpatchAssemble(thb, 0, gsQuadrature::HPatchRule,    2*p, p-1);
    CHECK_CLOSE(0.0, (Mg - Mh).cwiseAbs().maxCoeff(), 1e-11);
}

TEST(hpatchrule_2d_stiffness_vs_gauss_C0_knot_offcentre)
{
    const index_t p = 3;
    // C^0 line (knot multiplicity p) at x = 0.75: same target-multiplicity cap
    // (mu=3, m_deg+1=7, target regularity r_j-s = -1) as the centred fixture at
    // x = 0.5. numFallbackCells()==0: no Newton solve fails, so every leaf direction
    // wider than one element keeps the patch rule (a direction one element wide ties
    // with element-wise Gauss and uses it). Agreement with Gauss below therefore
    // verifies the min(...,m_deg+1) cap branch of the target space: the C^0 line at
    // x = 0.75 is an interior breakpoint of a leaf at least two elements wide in
    // direction 0. The centred fixture at x = 0.5 follows in
    // hpatchrule_2d_stiffness_vs_gauss_C0_knot_centre.
    gsTHBSplineBasis<2,real_t> thb = makeRepeatedKnotTHB(p, 3, 0.75, {2,0,0,4,4});

    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, p-2, false);
    CHECK_EQUAL(index_t(0), rule.numFallbackCells());

    gsMatrix<real_t> Kg = hpatchAssemble(thb, 1, gsQuadrature::GaussLegendre, 2*p, p-2);
    gsMatrix<real_t> Kh = hpatchAssemble(thb, 1, gsQuadrature::HPatchRule,    2*p, p-2);
    CHECK_CLOSE(0.0, (Kg - Kh).cwiseAbs().maxCoeff(), 1e-11);
}

TEST(hpatchrule_2d_stiffness_vs_gauss_C0_knot_centre)
{
    const index_t p = 3;
    // Same cap as the off-centre fixture (mu=3, target mult = m_deg+1 = 7,
    // target regularity r_j-s = -1) but at x = 0.5. The C^{-1} knot at the centre
    // splits the leaf's target space into two blocks; _ensureEven1D makes each
    // block even-dimensional, so the leaf keeps its patch rule (no Gauss fallback
    // cell). Its 617 points in total are a measured value, pinned so that the count guards the
    // target space of this leaf, the Gauss comparison its exactness.
    gsTHBSplineBasis<2,real_t> thb = makeRepeatedKnotTHB(p, 3, 0.5, {2,0,0,4,4});

    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, p-2, false);
    CHECK_EQUAL(index_t(0), rule.numFallbackCells());
    CHECK_EQUAL(index_t(617), rule.numPointsTotal());

    gsMatrix<real_t> Kg = hpatchAssemble(thb, 1, gsQuadrature::GaussLegendre, 2*p, p-2);
    gsMatrix<real_t> Kh = hpatchAssemble(thb, 1, gsQuadrature::HPatchRule,    2*p, p-2);
    CHECK_CLOSE(0.0, (Kg - Kh).cwiseAbs().maxCoeff(), 1e-11);
}

TEST(hpatchrule_2d_mass_vs_gauss_double_knot)
{
    const index_t p = 3;
    gsTHBSplineBasis<2,real_t> thb = makeRepeatedKnotTHB(p, 2, 0.5, {2,0,0,4,4});

    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, p-1, false);
    CHECK_EQUAL(index_t(0), rule.numFallbackCells());

    gsMatrix<real_t> Mg = hpatchAssemble(thb, 0, gsQuadrature::GaussLegendre, 2*p, p-1);
    gsMatrix<real_t> Mh = hpatchAssemble(thb, 0, gsQuadrature::HPatchRule,    2*p, p-1);
    CHECK_CLOSE(0.0, (Mg - Mh).cwiseAbs().maxCoeff(), 1e-11);
}

TEST(hpatchrule_2d_stiffness_vs_gauss_double_knot)
{
    const index_t p = 3;
    gsTHBSplineBasis<2,real_t> thb = makeRepeatedKnotTHB(p, 2, 0.5, {2,0,0,4,4});

    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, p-2, false);
    CHECK_EQUAL(index_t(0), rule.numFallbackCells());

    gsMatrix<real_t> Kg = hpatchAssemble(thb, 1, gsQuadrature::GaussLegendre, 2*p, p-2);
    gsMatrix<real_t> Kh = hpatchAssemble(thb, 1, gsQuadrature::HPatchRule,    2*p, p-2);
    CHECK_CLOSE(0.0, (Kg - Kh).cwiseAbs().maxCoeff(), 1e-11);
}

// ------------------------------------------------------------
// Single-leaf Newton solve of gsHPatchRule::_compute1D
//
// An unrefined THB basis has one level and one leaf, so the rule solves exactly one
// moment system per direction. For a leaf of n elements with simple interior knots and
// target (q, r) = (2p, r) the target space has
//     dim = (q+1) + (n-1)(q-r)
// functions; the rule has ceil(dim/2) nodes per direction (one knot is inserted when
// dim is odd). Targets: mass r = p-1, stiffness r = p-2, biharmonic r = p-3.
// ------------------------------------------------------------

// High degree (p = 4, 5, 6), 1D leaf of n = 4 and 16 elements: the rule keeps its
// patch rule (no Gauss fallback), has ceil(dim/2) nodes, and integrates every function
// of the target space exactly. nPts is ceil(dim/2) with dim from the formula above,
// e.g. p = 4, mass (8,3), n = 4: dim = 9 + 3*5 = 24, nPts = 12.
// The reference for exactness is (p+1)-point Gauss per element, exact to degree 2p+1.
TEST(hpatchrule_newton_1d_single_leaf_high_degree)
{
    struct Row { index_t p; index_t reg; index_t n; index_t nPts; };
    const Row rows[] = {
        {4,3, 4,12}, {4,3,16,42}, {4,2, 4,14}, {4,2,16,50}, {4,1, 4,15}, {4,1,16,57},
        {5,4, 4,15}, {5,4,16,51}, {5,3, 4,16}, {5,3,16,58}, {5,2, 4,18}, {5,2,16,66},
        {6,5, 4,17}, {6,5,16,59}, {6,4, 4,19}, {6,4,16,67}, {6,3, 4,20}, {6,3,16,74} };

    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        const index_t p = rows[k].p, reg = rows[k].reg, n = rows[k].n;
        gsInfo << "row p=" << p << " reg=" << reg << " n=" << n << "\n";

        gsKnotVector<real_t> kv(0, 1, n - 1, p + 1, 1);
        gsTensorBSplineBasis<1,real_t> tb(kv);
        gsTHBSplineBasis<1,real_t> thb(tb);
        gsHPatchRule<real_t> rule(*thb.domain(), 2*p, reg, false);
        CHECK_EQUAL(index_t(0),       rule.numFallbackCells());
        CHECK_EQUAL(rows[k].nPts,     rule.numPointsTotal());

        // Target space of the leaf before even-ing; even-ing only enlarges it.
        gsKnotVector<real_t> kvT(0, 1, n - 1, 2*p + 1, 2*p - reg);
        gsBSplineBasis<real_t> target(kvT);
        gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(1, p + 1));
        gsVector<real_t> Ig = hpatchIntegrateBasis(target, *tb.domain(), gauss);
        gsVector<real_t> Ih = hpatchIntegrateBasis(target, *tb.domain(), rule);
        // 1e-12 leaves at least 50 times headroom over the Newton stop ||res||_2 < 1e-14:
        // the residual is taken on the integrals of the even-ized basis, and each function
        // of this target space is one of them, or (when dim is odd and a knot is inserted)
        // a combination of at most two of them with coefficients in [0,1].
        CHECK_CLOSE(0.0, (Ih - Ig).cwiseAbs().maxCoeff(), 1e-12);
    }
}

// Same rows in 2D: both directions carry the same 1D system, so the node count is the
// square of the 1D one (e.g. 12*12 = 144) and the rule is exact on the tensor target.
TEST(hpatchrule_newton_2d_single_leaf_high_degree)
{
    struct Row { index_t p; index_t reg; index_t n; index_t nPts; };
    const Row rows[] = {
        {4,3, 4,144}, {4,3,16,1764}, {4,2, 4,196}, {4,2,16,2500}, {4,1, 4,225}, {4,1,16,3249},
        {5,4, 4,225}, {5,4,16,2601}, {5,3, 4,256}, {5,3,16,3364}, {5,2, 4,324}, {5,2,16,4356},
        {6,5, 4,289}, {6,5,16,3481}, {6,4, 4,361}, {6,4,16,4489}, {6,3, 4,400}, {6,3,16,5476} };

    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        const index_t p = rows[k].p, reg = rows[k].reg, n = rows[k].n;
        gsInfo << "row p=" << p << " reg=" << reg << " n=" << n << "\n";

        gsKnotVector<real_t> kv(0, 1, n - 1, p + 1, 1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> thb(tb);
        gsHPatchRule<real_t> rule(*thb.domain(), 2*p, reg, false);
        CHECK_EQUAL(index_t(0),       rule.numFallbackCells());
        CHECK_EQUAL(rows[k].nPts,     rule.numPointsTotal());

        gsKnotVector<real_t> kvT(0, 1, n - 1, 2*p + 1, 2*p - reg);
        gsTensorBSplineBasis<2,real_t> target(kvT, kvT);
        gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(2, p + 1));
        gsVector<real_t> Ig = hpatchIntegrateBasis(target, *tb.domain(), gauss);
        gsVector<real_t> Ih = hpatchIntegrateBasis(target, *tb.domain(), rule);
        CHECK_CLOSE(0.0, (Ih - Ig).cwiseAbs().maxCoeff(), 1e-12);
    }
}

// On single leaves where the Newton solve converges without any safeguard
// (p = 2, 3; mass and stiffness targets), gsHPatchRule and the tensor gsPatchRule solve
// the same system from the same starting guess (pairwise Greville midpoints) with the
// same even-ing, so their nodes and weights must agree to round-off. Both rules are
// first checked exact on the target space, otherwise agreement would not show that the
// solved rule is the right one. The whole leaf is compared with one mapTo on [0,1]
// (not element by element: a node at an element boundary such as 0.5 may fall on
// either side of the half-open [lower, upper) filter in the two independent solves).
TEST(hpatchrule_newton_1d_single_leaf_matches_patchrule)
{
    struct Row { index_t p; index_t reg; index_t n; };
    const Row rows[] = { {2,1, 4}, {2,1,16}, {2,0, 4}, {2,0,16},
                         {3,2, 4}, {3,2,16}, {3,1, 4}, {3,1,16} };

    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        const index_t p = rows[k].p, reg = rows[k].reg, n = rows[k].n;
        gsInfo << "row p=" << p << " reg=" << reg << " n=" << n << "\n";

        gsKnotVector<real_t> kv(0, 1, n - 1, p + 1, 1);
        gsTensorBSplineBasis<1,real_t> tb(kv);
        gsTHBSplineBasis<1,real_t> thb(tb);
        gsHPatchRule<real_t> rule(*thb.domain(), 2*p, reg, false);
        gsPatchRule<real_t> ref(*tb.domain(), 2*p, reg, false);

        gsKnotVector<real_t> kvT(0, 1, n - 1, 2*p + 1, 2*p - reg);
        gsBSplineBasis<real_t> target(kvT);
        gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(1, p + 1));
        gsVector<real_t> Ig = hpatchIntegrateBasis(target, *tb.domain(), gauss);
        CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target, *tb.domain(), ref)  - Ig).cwiseAbs().maxCoeff(), 1e-12);
        CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target, *tb.domain(), rule) - Ig).cwiseAbs().maxCoeff(), 1e-12);

        const gsVector<real_t> lower = gsVector<real_t>::Zero(1), upper = gsVector<real_t>::Ones(1);
        gsMatrix<real_t> nodesH, nodesP;
        gsVector<real_t> weightsH, weightsP;
        rule.mapTo(lower, upper, nodesH, weightsH);
        ref.mapTo(lower, upper, nodesP, weightsP);
        // Different counts mean different target spaces: fail loudly, then skip the
        // value comparison, which would otherwise be undefined.
        CHECK_EQUAL(weightsP.size(), weightsH.size());
        if (weightsP.size() == weightsH.size())
        {
            CHECK_CLOSE(0.0, (nodesH   - nodesP  ).cwiseAbs().maxCoeff(), 1e-12);
            CHECK_CLOSE(0.0, (weightsH - weightsP).cwiseAbs().maxCoeff(), 1e-12);
        }
    }
}

// Same comparison in 2D: the rule is the tensor product of the 1D rules and is
// compared in the tensor order of mapTo (direction 0 fastest).
TEST(hpatchrule_newton_2d_single_leaf_matches_patchrule)
{
    struct Row { index_t p; index_t reg; index_t n; };
    const Row rows[] = { {3,1, 4}, {2,1,16} };

    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        const index_t p = rows[k].p, reg = rows[k].reg, n = rows[k].n;
        gsInfo << "row p=" << p << " reg=" << reg << " n=" << n << "\n";

        gsKnotVector<real_t> kv(0, 1, n - 1, p + 1, 1);
        gsTensorBSplineBasis<2,real_t> tb(kv, kv);
        gsTHBSplineBasis<2,real_t> thb(tb);
        gsHPatchRule<real_t> rule(*thb.domain(), 2*p, reg, false);
        gsPatchRule<real_t> ref(*tb.domain(), 2*p, reg, false);

        gsKnotVector<real_t> kvT(0, 1, n - 1, 2*p + 1, 2*p - reg);
        gsTensorBSplineBasis<2,real_t> target(kvT, kvT);
        gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(2, p + 1));
        gsVector<real_t> Ig = hpatchIntegrateBasis(target, *tb.domain(), gauss);
        CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target, *tb.domain(), ref)  - Ig).cwiseAbs().maxCoeff(), 1e-12);
        CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target, *tb.domain(), rule) - Ig).cwiseAbs().maxCoeff(), 1e-12);

        const gsVector<real_t> lower = gsVector<real_t>::Zero(2), upper = gsVector<real_t>::Ones(2);
        gsMatrix<real_t> nodesH, nodesP;
        gsVector<real_t> weightsH, weightsP;
        rule.mapTo(lower, upper, nodesH, weightsH);
        ref.mapTo(lower, upper, nodesP, weightsP);
        // Different counts mean different target spaces: fail loudly, then skip the
        // value comparison, which would otherwise be undefined.
        CHECK_EQUAL(weightsP.size(), weightsH.size());
        if (weightsP.size() == weightsH.size())
        {
            CHECK_CLOSE(0.0, (nodesH   - nodesP  ).cwiseAbs().maxCoeff(), 1e-12);
            CHECK_CLOSE(0.0, (weightsH - weightsP).cwiseAbs().maxCoeff(), 1e-12);
        }
    }
}

// A 1D leaf of 2000 elements, target (6,2) (mass, p = 3): the shape of a long leaf of a
// hierarchical mesh. The target space has dim = 7 + 1999*4 = 8003 functions (odd, so one
// knot is inserted), hence ceil(8003/2) = 4002 nodes. The rule must build without
// falling back and integrate the target space exactly. The memory footprint of the
// banded solve is measured outside the test (peak RSS of this test alone against a small
// baseline test), since the test also passes with a dense Jacobian.
TEST(hpatchrule_newton_1d_long_leaf)
{
    const index_t p = 3, reg = 2, n = 2000;
    gsKnotVector<real_t> kv(0, 1, n - 1, p + 1, 1);
    gsTensorBSplineBasis<1,real_t> tb(kv);
    gsTHBSplineBasis<1,real_t> thb(tb);
    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, reg, false);
    CHECK_EQUAL(index_t(0),    rule.numFallbackCells());
    CHECK_EQUAL(index_t(4002), rule.numPointsTotal());

    gsKnotVector<real_t> kvT(0, 1, n - 1, 2*p + 1, 2*p - reg);
    gsBSplineBasis<real_t> target(kvT);
    gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(1, p + 1));
    gsVector<real_t> Ig = hpatchIntegrateBasis(target, *tb.domain(), gauss);
    gsVector<real_t> Ih = hpatchIntegrateBasis(target, *tb.domain(), rule);
    CHECK_CLOSE(0.0, (Ih - Ig).cwiseAbs().maxCoeff(), 1e-12);
}

// ------------------------------------------------------------
// Lookup of gsHPatchRule::mapTo on a hierarchical mesh with many leaves
//
// The fixture is makeRefinedTHB(2) with three more level-2 boxes (a centre square, the
// upper-right corner and a strip on the west side), so that the cell lookup has to find
// its leaf among many (10 leaves, counted through the tree's leaf iterator and pinned in
// the shuffled-order and boundary tests, so a change of the fixture is noticed). The rule (2p, p-1) is the mass target.
// ------------------------------------------------------------

// A rule returns, for an element, exactly what it returns when the elements are queried
// in volume order, whatever the query order: the lookup finds the one leaf that contains
// the element (each element lies in exactly one leaf), so nodes and weights are
// bitwise identical.
TEST(hpatchrule_lookup_shuffled_order)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    std::vector<index_t> more = {2, 12,12,20,20,   2, 24,24,32,32,   2, 0,24,4,32};
    thb.refineElements(more);

    index_t nLeaves = 0;
    for (auto leaf = thb.tree().beginLeafIterator(); leaf.good(); leaf.next())
        ++nLeaves;
    CHECK_EQUAL(index_t(10), nLeaves);

    auto dom = thb.domain();
    std::vector<gsVector<real_t> > lows, upps;
    for (auto & elem : dom->allElements())
    {
        lows.push_back(elem.lowerCorner());
        upps.push_back(elem.upperCorner());
    }
    const size_t nElem = lows.size();

    // Reference: one rule, elements in volume order.
    gsHPatchRule<real_t> first(*dom, 2*p, p-1, false);
    std::vector<gsMatrix<real_t> > refNodes(nElem);
    std::vector<gsVector<real_t> > refWeights(nElem);
    for (size_t e = 0; e != nElem; ++e)
        first.mapTo(lows[e], upps[e], refNodes[e], refWeights[e]);

    // A second rule on the same domain, queried in reversed, shuffled and volume order.
    // The orders chain on one rule, so each starts from wherever the last one stopped.
    std::vector<size_t> volume(nElem);
    for (size_t e = 0; e != nElem; ++e) volume[e] = e;
    std::vector<std::vector<size_t> > orders(3, volume);
    std::reverse(orders[0].begin(), orders[0].end());
    std::mt19937 rng(20240607);
    std::shuffle(orders[1].begin(), orders[1].end(), rng);

    gsHPatchRule<real_t> second(*dom, 2*p, p-1, false);
    gsMatrix<real_t> nodes;
    gsVector<real_t> weights;
    for (size_t o = 0; o != orders.size(); ++o)
        for (size_t i = 0; i != nElem; ++i)
        {
            const size_t e = orders[o][i];
            second.mapTo(lows[e], upps[e], nodes, weights);
            CHECK_EQUAL(refWeights[e].size(), weights.size());
            if (refWeights[e].size() == weights.size() && weights.size() > 0)
            {
                CHECK_EQUAL(real_t(0), (nodes   - refNodes[e]  ).cwiseAbs().maxCoeff());
                CHECK_EQUAL(real_t(0), (weights - refWeights[e]).cwiseAbs().maxCoeff());
            }
        }
}

// Boundary iteration with a fixed direction: the elements of a side are degenerate in
// the fixed direction (one node on the side) and the boundary iterator visits only the
// leaves touching the side. The boundary mass matrix of the THB basis along a side is
// a product of two basis functions, which lies in the target space (2p, p-1) of the
// leaf, so the rule must reproduce the (p+1)-point Gauss matrix (exact to degree 2p+1
// along the side; one point in the degenerate normal direction, weight 1). One rule per
// direction serves both opposite sides, so the second side starts the lookup from the
// cell the first one ended on. The functions are a partition of unity, hence the
// entries of the mass matrix sum to the side length 1.
TEST(hpatchrule_lookup_boundary_sides)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    std::vector<index_t> more = {2, 12,12,20,20,   2, 24,24,32,32,   2, 0,24,4,32};
    thb.refineElements(more);

    index_t nLeaves = 0;
    for (auto leaf = thb.tree().beginLeafIterator(); leaf.good(); leaf.next())
        ++nLeaves;
    CHECK_EQUAL(index_t(10), nLeaves);

    const boxSide sides[2][2] = { {boundary::west,  boundary::east},
                                  {boundary::south, boundary::north} };
    for (short_t dir = 0; dir != 2; ++dir)
    {
        gsHPatchRule<real_t> rule(*thb.domain(), 2*p, p-1, false, dir);

        // 1 point normal to the side, p+1 points along it.
        gsVector<index_t> nn(2);
        nn[dir] = 1;
        nn[1 - dir] = p + 1;
        gsGaussRule<real_t> gauss(nn);

        for (index_t s = 0; s != 2; ++s)
        {
            gsInfo << "row dir=" << dir << " side=" << sides[dir][s] << "\n";
            gsMatrix<real_t> Mh = hpatchBoundaryMass(thb, sides[dir][s], rule);
            gsMatrix<real_t> Mg = hpatchBoundaryMass(thb, sides[dir][s], gauss);
            CHECK_CLOSE(0.0, (Mh - Mg).cwiseAbs().maxCoeff(), 1e-12);
            CHECK_CLOSE(1.0, Mh.sum(), 1e-12);
        }
    }
}

// A query that no cell contains raises an error in every build (the lookup cannot return
// nodes for it); the rule stays usable afterwards. The whole domain straddles several
// leaves, and a box outside the domain has no leaf at all. Afterwards one real element
// must give what a fresh rule gives.
TEST(hpatchrule_lookup_miss_throws)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    std::vector<index_t> more = {2, 12,12,20,20,   2, 24,24,32,32,   2, 0,24,4,32};
    thb.refineElements(more);

    auto dom = thb.domain();
    gsHPatchRule<real_t> rule(*dom, 2*p, p-1, false);
    gsMatrix<real_t> nodes;
    gsVector<real_t> weights;

    const gsVector<real_t> zero = gsVector<real_t>::Zero(2);
    const gsVector<real_t> one  = gsVector<real_t>::Ones(2);
    CHECK_THROW(rule.mapTo(zero, one, nodes, weights), std::runtime_error);

    const gsVector<real_t> outLow = gsVector<real_t>::Constant(2, 1.25);
    const gsVector<real_t> outUpp = gsVector<real_t>::Constant(2, 1.5);
    CHECK_THROW(rule.mapTo(outLow, outUpp, nodes, weights), std::runtime_error);

    // First element of the volume iteration, after the two failed queries.
    auto elem = dom->beginAll();
    const gsVector<real_t> low = elem.lowerCorner(), upp = elem.upperCorner();
    rule.mapTo(low, upp, nodes, weights);

    gsHPatchRule<real_t> fresh(*dom, 2*p, p-1, false);
    gsMatrix<real_t> refNodes;
    gsVector<real_t> refWeights;
    fresh.mapTo(low, upp, refNodes, refWeights);
    CHECK_EQUAL(refWeights.size(), weights.size());
    if (refWeights.size() == weights.size() && weights.size() > 0)
    {
        CHECK_EQUAL(real_t(0), (nodes   - refNodes  ).cwiseAbs().maxCoeff());
        CHECK_EQUAL(real_t(0), (weights - refWeights).cwiseAbs().maxCoeff());
    }
}

// ------------------------------------------------------------
// Per-direction fallback of gsHPatchRule::_build
//
// For a leaf direction with n elements and nq = m_deg/2 + 1 Gauss points per element, the
// patch rule has n1d nodes. The direction uses element-wise Gauss (n*nq nodes) when
// n1d >= n*nq, or when its moment solve fails; the other directions keep their patch rule.
// A leaf is counted by numFallbackCells() once, when a solve failed; a tie is not counted.
// ------------------------------------------------------------

// Each direction is solved on [0,1] and scaled to the leaf, with the stop of a leaf of
// length 1 for every length of at least 1. A direction of length 1e4 therefore gets the
// patch rule of a direction of length 1 (scaled), and the leaf does not fall back. Single
// leaf, p = 2, 4 elements per direction, target (4,1): the patch rule has
// dim = 5 + 3*3 = 14 -> 7 nodes per direction, Gauss 4*nq = 4*3 = 12.
// Exactness is the relative error, since the integrals scale with the length.
//
// The last block has a direction with 3 elements of lengths in ratio 1 : 100 : 1e4 whose
// moment solve fails: it uses 3*3 = 9 Gauss nodes, while the other direction (4 uniform
// elements) keeps its 7 patch nodes. The leaf is counted once.
TEST(hpatchrule_fallback_newton_fails_in_one_direction_2d)
{
    struct Row { real_t L0; real_t L1; index_t nPts; index_t nFb; };
    const Row rows[] = { {1,1,    49,0}, {1e4,1,    49,0}, {1,1e4,    49,0} };

    const index_t p = 2, reg = 1;
    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        gsInfo << "row L0=" << rows[k].L0 << " L1=" << rows[k].L1 << "\n";
        gsKnotVector<real_t> kv0(0, rows[k].L0, 3, p + 1, 1);
        gsKnotVector<real_t> kv1(0, rows[k].L1, 3, p + 1, 1);
        gsTensorBSplineBasis<2,real_t> tb(kv0, kv1);
        gsTHBSplineBasis<2,real_t> thb(tb);
        gsHPatchRule<real_t> rule(*thb.domain(), 2*p, reg, false);
        CHECK_EQUAL(rows[k].nFb,  rule.numFallbackCells());
        CHECK_EQUAL(rows[k].nPts, rule.numPointsTotal());

        gsKnotVector<real_t> kvT0(0, rows[k].L0, 3, 2*p + 1, 2*p - reg);
        gsKnotVector<real_t> kvT1(0, rows[k].L1, 3, 2*p + 1, 2*p - reg);
        gsTensorBSplineBasis<2,real_t> target(kvT0, kvT1);
        gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(2, p + 1));
        gsVector<real_t> Ig = hpatchIntegrateBasis(target, *tb.domain(), gauss);
        gsVector<real_t> Ih = hpatchIntegrateBasis(target, *tb.domain(), rule);
        CHECK_CLOSE(0.0, (Ih - Ig).cwiseAbs().maxCoeff() / Ig.cwiseAbs().maxCoeff(), 1e-12);
    }

    // Direction 0: breakpoints 0, 1/s, 101/s, 1 with s = 1 + 100 + 1e4, simple interior
    // knots in the basis and multiplicity 2p - reg = 3 in the target.
    const real_t s = 1 + 100 + 1e4, b1 = 1 / s, b2 = 101 / s;
    const real_t knots0[]  = {0,0,0,    b1,    b2,    1,1,1};
    const real_t knotsT0[] = {0,0,0,0,0, b1,b1,b1, b2,b2,b2, 1,1,1,1,1};
    gsKnotVector<real_t> kv0(std::vector<real_t>(knots0, knots0 + 8), p);
    gsKnotVector<real_t> kvT0(std::vector<real_t>(knotsT0, knotsT0 + 16), 2*p);
    gsKnotVector<real_t> kv1(0, 1, 3, p + 1, 1);
    gsKnotVector<real_t> kvT1(0, 1, 3, 2*p + 1, 2*p - reg);
    gsTensorBSplineBasis<2,real_t> tb(kv0, kv1);
    gsTHBSplineBasis<2,real_t> thb(tb);
    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, reg, false);
    CHECK_EQUAL(index_t(1),  rule.numFallbackCells());
    CHECK_EQUAL(index_t(63), rule.numPointsTotal());

    gsTensorBSplineBasis<2,real_t> target(kvT0, kvT1);
    gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(2, p + 1));
    gsVector<real_t> Ig = hpatchIntegrateBasis(target, *tb.domain(), gauss);
    gsVector<real_t> Ih = hpatchIntegrateBasis(target, *tb.domain(), rule);
    CHECK_CLOSE(0.0, (Ih - Ig).cwiseAbs().maxCoeff() / Ig.cwiseAbs().maxCoeff(), 1e-12);
}

// 3D: with two directions of length 1e4 the leaf keeps its patch rule in all directions
// (7 nodes each). The last block has two directions that fail, graded as in the 2D test:
// they use 9 Gauss nodes each and the third keeps 7 patch nodes; the leaf is counted
// once, not once per direction.
TEST(hpatchrule_fallback_newton_fails_in_one_direction_3d)
{
    struct Row { real_t L0; real_t L1; real_t L2; index_t nPts; index_t nFb; };
    const Row rows[] = { {1,1,1,    343,0}, {1e4,1e4,1,    343,0} };

    const index_t p = 2, reg = 1;
    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        gsInfo << "row L0=" << rows[k].L0 << " L1=" << rows[k].L1 << " L2=" << rows[k].L2 << "\n";
        gsKnotVector<real_t> kv0(0, rows[k].L0, 3, p + 1, 1);
        gsKnotVector<real_t> kv1(0, rows[k].L1, 3, p + 1, 1);
        gsKnotVector<real_t> kv2(0, rows[k].L2, 3, p + 1, 1);
        gsTensorBSplineBasis<3,real_t> tb(kv0, kv1, kv2);
        gsTHBSplineBasis<3,real_t> thb(tb);
        gsHPatchRule<real_t> rule(*thb.domain(), 2*p, reg, false);
        CHECK_EQUAL(rows[k].nFb,  rule.numFallbackCells());
        CHECK_EQUAL(rows[k].nPts, rule.numPointsTotal());

        gsKnotVector<real_t> kvT0(0, rows[k].L0, 3, 2*p + 1, 2*p - reg);
        gsKnotVector<real_t> kvT1(0, rows[k].L1, 3, 2*p + 1, 2*p - reg);
        gsKnotVector<real_t> kvT2(0, rows[k].L2, 3, 2*p + 1, 2*p - reg);
        gsTensorBSplineBasis<3,real_t> target(kvT0, kvT1, kvT2);
        gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(3, p + 1));
        gsVector<real_t> Ig = hpatchIntegrateBasis(target, *tb.domain(), gauss);
        gsVector<real_t> Ih = hpatchIntegrateBasis(target, *tb.domain(), rule);
        CHECK_CLOSE(0.0, (Ih - Ig).cwiseAbs().maxCoeff() / Ig.cwiseAbs().maxCoeff(), 1e-12);
    }

    // Directions 0 and 1: the graded knot vector of the 2D test (breakpoints 0, 1/s, 101/s, 1),
    // whose solve fails, so 3*3 = 9 Gauss nodes each; direction 2: 4 uniform elements with 7
    // patch nodes. Expected 9*9*7 = 567 nodes and one counted leaf.
    const real_t s = 1 + 100 + 1e4, b1 = 1 / s, b2 = 101 / s;
    const real_t knots0[]  = {0,0,0,    b1,    b2,    1,1,1};
    const real_t knotsT0[] = {0,0,0,0,0, b1,b1,b1, b2,b2,b2, 1,1,1,1,1};
    gsKnotVector<real_t> kv0(std::vector<real_t>(knots0, knots0 + 8), p);
    gsKnotVector<real_t> kvT0(std::vector<real_t>(knotsT0, knotsT0 + 16), 2*p);
    gsKnotVector<real_t> kv2(0, 1, 3, p + 1, 1);
    gsKnotVector<real_t> kvT2(0, 1, 3, 2*p + 1, 2*p - reg);
    gsTensorBSplineBasis<3,real_t> tb(kv0, kv0, kv2);
    gsTHBSplineBasis<3,real_t> thb(tb);
    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, reg, false);
    CHECK_EQUAL(index_t(1),   rule.numFallbackCells());
    CHECK_EQUAL(index_t(567), rule.numPointsTotal());

    gsTensorBSplineBasis<3,real_t> target(kvT0, kvT0, kvT2);
    gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(3, p + 1));
    gsVector<real_t> Ig = hpatchIntegrateBasis(target, *tb.domain(), gauss);
    gsVector<real_t> Ih = hpatchIntegrateBasis(target, *tb.domain(), rule);
    CHECK_CLOSE(0.0, (Ih - Ig).cwiseAbs().maxCoeff() / Ig.cwiseAbs().maxCoeff(), 1e-12);
}

// A leaf one element wide in a direction ties there: the target space of a single
// element of degree 2p has dim = 2p+1 (odd), hence p+1 nodes after even-ing, which is
// the p+1 = nq Gauss points of the element. The direction uses Gauss, and the tie is
// not counted as a fallback. The mixed row has direction 0 tying (3 nodes) and direction
// 1 with 4 elements on the patch rule (dim 5 + 3*3 = 14 -> 7 nodes), 3*7 = 21.
// Reference: Gauss with p+1 points per element, exact to degree 2p+1.
TEST(hpatchrule_fallback_tie_not_counted)
{
    struct Row { index_t p; index_t n1; index_t n2; index_t n3; };
    const Row rows[] = { {2,3,9,27}, {3,4,16,64} };

    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        const index_t p = rows[k].p;
        gsInfo << "row p=" << p << "\n";
        gsKnotVector<real_t> kv(0, 1, 0, p + 1, 1);
        gsKnotVector<real_t> kvT(0, 1, 0, 2*p + 1, p + 1);
        gsGaussRule<real_t> gauss1(gsVector<index_t>::Constant(1, p + 1));
        gsGaussRule<real_t> gauss2(gsVector<index_t>::Constant(2, p + 1));
        gsGaussRule<real_t> gauss3(gsVector<index_t>::Constant(3, p + 1));

        gsTensorBSplineBasis<1,real_t> tb1(kv);
        gsTHBSplineBasis<1,real_t> thb1(tb1);
        gsHPatchRule<real_t> rule1(*thb1.domain(), 2*p, p - 1, false);
        CHECK_EQUAL(index_t(0),  rule1.numFallbackCells());
        CHECK_EQUAL(rows[k].n1,  rule1.numPointsTotal());
        gsBSplineBasis<real_t> target1(kvT);
        CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target1, *tb1.domain(), rule1)
                          - hpatchIntegrateBasis(target1, *tb1.domain(), gauss1)).cwiseAbs().maxCoeff(), 1e-12);

        gsTensorBSplineBasis<2,real_t> tb2(kv, kv);
        gsTHBSplineBasis<2,real_t> thb2(tb2);
        gsHPatchRule<real_t> rule2(*thb2.domain(), 2*p, p - 1, false);
        CHECK_EQUAL(index_t(0),  rule2.numFallbackCells());
        CHECK_EQUAL(rows[k].n2,  rule2.numPointsTotal());
        gsTensorBSplineBasis<2,real_t> target2(kvT, kvT);
        CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target2, *tb2.domain(), rule2)
                          - hpatchIntegrateBasis(target2, *tb2.domain(), gauss2)).cwiseAbs().maxCoeff(), 1e-12);

        gsTensorBSplineBasis<3,real_t> tb3(kv, kv, kv);
        gsTHBSplineBasis<3,real_t> thb3(tb3);
        gsHPatchRule<real_t> rule3(*thb3.domain(), 2*p, p - 1, false);
        CHECK_EQUAL(index_t(0),  rule3.numFallbackCells());
        CHECK_EQUAL(rows[k].n3,  rule3.numPointsTotal());
        gsTensorBSplineBasis<3,real_t> target3(kvT, kvT, kvT);
        CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target3, *tb3.domain(), rule3)
                          - hpatchIntegrateBasis(target3, *tb3.domain(), gauss3)).cwiseAbs().maxCoeff(), 1e-12);
    }

    // Mixed: direction 0 one element (tie, Gauss), direction 1 four elements (patch).
    const index_t p = 2;
    gsKnotVector<real_t> kv0(0, 1, 0, p + 1, 1), kv1(0, 1, 3, p + 1, 1);
    gsTensorBSplineBasis<2,real_t> tb(kv0, kv1);
    gsTHBSplineBasis<2,real_t> thb(tb);
    gsHPatchRule<real_t> rule(*thb.domain(), 2*p, p - 1, false);
    CHECK_EQUAL(index_t(0),  rule.numFallbackCells());
    CHECK_EQUAL(index_t(21), rule.numPointsTotal());
    gsKnotVector<real_t> kvT0(0, 1, 0, 2*p + 1, p + 1), kvT1(0, 1, 3, 2*p + 1, p + 1);
    gsTensorBSplineBasis<2,real_t> target(kvT0, kvT1);
    gsGaussRule<real_t> gauss(gsVector<index_t>::Constant(2, p + 1));
    CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target, *tb.domain(), rule)
                      - hpatchIntegrateBasis(target, *tb.domain(), gauss)).cwiseAbs().maxCoeff(), 1e-12);
}

// An odd target degree builds an exact rule: p = 2 basis, target degree 5, so the
// element-wise Gauss reference has nq = 5/2 + 1 = 3 points per element (exact to degree 5).
// The count is ceil(dim/2) with dim = 6 + (n-1)*(5-reg) for simple interior knots: e.g.
// reg 4, n 4: 6 + 3*1 = 9 -> 10 (even) -> 5 nodes. The last row (n = 1) has dim 6 -> 3 = nq, a
// tie. No row counts as a fallback.
TEST(hpatchrule_fallback_odd_degree)
{
    struct Row { index_t reg; index_t n; index_t nPts1; index_t nPts2; };
    const Row rows[] = { {4,4, 5,25}, {3,4, 6,36}, {1,4, 9,81}, {3,1, 3,9} };

    const index_t p = 2, deg = 5;
    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        const index_t reg = rows[k].reg, n = rows[k].n;
        gsInfo << "row reg=" << reg << " n=" << n << "\n";
        gsKnotVector<real_t> kv(0, 1, n - 1, p + 1, 1);
        gsKnotVector<real_t> kvT(0, 1, n - 1, deg + 1, deg - reg);

        gsTensorBSplineBasis<1,real_t> tb1(kv);
        gsTHBSplineBasis<1,real_t> thb1(tb1);
        gsHPatchRule<real_t> rule1(*thb1.domain(), deg, reg, false);
        CHECK_EQUAL(index_t(0),     rule1.numFallbackCells());
        CHECK_EQUAL(rows[k].nPts1,  rule1.numPointsTotal());
        gsBSplineBasis<real_t> target1(kvT);
        gsGaussRule<real_t> gauss1(gsVector<index_t>::Constant(1, deg/2 + 1));
        CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target1, *tb1.domain(), rule1)
                          - hpatchIntegrateBasis(target1, *tb1.domain(), gauss1)).cwiseAbs().maxCoeff(), 1e-12);

        gsTensorBSplineBasis<2,real_t> tb2(kv, kv);
        gsTHBSplineBasis<2,real_t> thb2(tb2);
        gsHPatchRule<real_t> rule2(*thb2.domain(), deg, reg, false);
        CHECK_EQUAL(index_t(0),     rule2.numFallbackCells());
        CHECK_EQUAL(rows[k].nPts2,  rule2.numPointsTotal());
        gsTensorBSplineBasis<2,real_t> target2(kvT, kvT);
        gsGaussRule<real_t> gauss2(gsVector<index_t>::Constant(2, deg/2 + 1));
        CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target2, *tb2.domain(), rule2)
                          - hpatchIntegrateBasis(target2, *tb2.domain(), gauss2)).cwiseAbs().maxCoeff(), 1e-12);
    }
}

// Target degree 6 above 2p+1 = 3 for p = 1 (reg 0, nq = 6/2 + 1 = 4). (p+1)-point Gauss
// is exact only to degree 2p+1 = 3, so element-wise Gauss in the target space needs
// deg/2 + 1 = 4 points per element, and the reference uses 4. A single element ties
// (dim 7 -> 8 -> 4 nodes = nq) in 1D and 2D and uses Gauss without being counted. With
// 4 elements per direction (dim 7 + 3*6 = 25 -> 26) each direction has 13 patch nodes, also
// for a direction of length 1e4, which is solved on [0,1] like a direction of length 1:
// 13*13 = 169 nodes, no counted leaf. The 1e4 row is compared in the relative error.
TEST(hpatchrule_fallback_degree_above_2p_plus_1)
{
    const index_t p = 1, deg = 6, reg = 0;
    gsGaussRule<real_t> gauss1(gsVector<index_t>::Constant(1, deg/2 + 1));
    gsGaussRule<real_t> gauss2(gsVector<index_t>::Constant(2, deg/2 + 1));

    gsKnotVector<real_t> kv(0, 1, 0, p + 1, 1);
    gsKnotVector<real_t> kvT(0, 1, 0, deg + 1, deg - reg);

    gsTensorBSplineBasis<1,real_t> tb1(kv);
    gsTHBSplineBasis<1,real_t> thb1(tb1);
    gsHPatchRule<real_t> rule1(*thb1.domain(), deg, reg, false);
    CHECK_EQUAL(index_t(0), rule1.numFallbackCells());
    CHECK_EQUAL(index_t(4), rule1.numPointsTotal());
    gsBSplineBasis<real_t> target1(kvT);
    CHECK_CLOSE(0.0, (hpatchIntegrateBasis(target1, *tb1.domain(), rule1)
                      - hpatchIntegrateBasis(target1, *tb1.domain(), gauss1)).cwiseAbs().maxCoeff(), 1e-12);

    struct Row { index_t n0; index_t n1; real_t L0; index_t nPts; index_t nFb; };
    const Row rows[] = { {1,1,1,    16,0}, {4,4,1e4,    169,0} };
    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        gsInfo << "row n0=" << rows[k].n0 << " n1=" << rows[k].n1 << " L0=" << rows[k].L0 << "\n";
        gsKnotVector<real_t> kv0(0, rows[k].L0, rows[k].n0 - 1, p + 1, 1);
        gsKnotVector<real_t> kv1(0, 1,          rows[k].n1 - 1, p + 1, 1);
        gsTensorBSplineBasis<2,real_t> tb2(kv0, kv1);
        gsTHBSplineBasis<2,real_t> thb2(tb2);
        gsHPatchRule<real_t> rule2(*thb2.domain(), deg, reg, false);
        CHECK_EQUAL(rows[k].nFb,  rule2.numFallbackCells());
        CHECK_EQUAL(rows[k].nPts, rule2.numPointsTotal());

        gsKnotVector<real_t> kvT0(0, rows[k].L0, rows[k].n0 - 1, deg + 1, deg - reg);
        gsKnotVector<real_t> kvT1(0, 1,          rows[k].n1 - 1, deg + 1, deg - reg);
        gsTensorBSplineBasis<2,real_t> target2(kvT0, kvT1);
        gsVector<real_t> Ig = hpatchIntegrateBasis(target2, *tb2.domain(), gauss2);
        gsVector<real_t> Ih = hpatchIntegrateBasis(target2, *tb2.domain(), rule2);
        CHECK_CLOSE(0.0, (Ih - Ig).cwiseAbs().maxCoeff() / Ig.cwiseAbs().maxCoeff(), 1e-12);
    }
}

// ============================================================
// gsHPatchRuleCache, rule handles and gsHPatchRule::makeFactory
// ============================================================

// Largest |difference| of the nodes and weights that `a` and `b` give on every element of
// `basis`: the volume elements for fd = -1, else the elements of both sides normal to
// direction fd (degenerate in fd). Returns +infinity if the two give a different number of
// points on some element. `nElem` returns the number of elements visited, so that a caller
// can rule out a vacuous zero.
real_t hpatchMaxDiff(const gsQuadRule<real_t> & a, const gsQuadRule<real_t> & b,
                     const gsBasis<real_t> & basis, short_t fd, index_t & nElem)
{
    auto dom = basis.domain();
    real_t diff = 0;
    nElem = 0;
    gsMatrix<real_t> na, nb;
    gsVector<real_t> wa, wb;
    // Sides 2*fd+1 and 2*fd+2 are the two sides normal to direction fd.
    for (short_t k = 0; k != (fd < 0 ? 1 : 2); ++k)
    {
        const boxSide side(static_cast<short_t>(fd < 0 ? 1 : 2*fd + 1 + k));
        gsDomain<real_t>::iterator it  = (fd < 0 ? dom->beginAll() : dom->beginBdr(side));
        gsDomain<real_t>::iterator end = (fd < 0 ? dom->endAll()   : dom->endBdr(side));
        for (; it < end; ++it)
        {
            ++nElem;
            a.mapTo(it.lowerCorner(), it.upperCorner(), na, wa);
            b.mapTo(it.lowerCorner(), it.upperCorner(), nb, wb);
            if (na.rows() != nb.rows() || na.cols() != nb.cols() || wa.size() != wb.size())
                return std::numeric_limits<real_t>::infinity();
            if (wa.size() > 0)
                diff = (std::max)(diff, (std::max)((na - nb).cwiseAbs().maxCoeff(),
                                                   (wa - wb).cwiseAbs().maxCoeff()));
        }
    }
    return diff;
}

// Mass (form 0) or stiffness (form 1) matrix of `thb` on the identity geometry, assembled
// with the leaf rule (deg, reg): through the options (quRule = HPatchRule), or through a
// factory from gsHPatchRule::makeFactory installed on the same assembler.
template<short_t d>
gsMatrix<real_t> hpatchAssembleLeafRule(gsTHBSplineBasis<d,real_t> thb, int form,
                                        index_t deg, index_t reg, bool useFactory)
{
    gsTensorBSplineBasis<d,real_t> & l0 = thb.tensorLevel(0);
    gsMatrix<real_t> coefs = l0.anchors().transpose();
    gsTensorBSpline<d,real_t> geom(l0, coefs);
    gsMultiPatch<real_t> mp; mp.addPatch(geom);
    gsMultiBasis<real_t> mb; mb.addBasis(thb.clone().release()); mb.setTopology(mp);

    gsExprAssembler<real_t> A(1,1);
    gsOptionList opt = gsExprAssembler<real_t>::defaultOptions();
    opt.setInt("quRule", gsQuadrature::HPatchRule);
    opt.setReal("quA", static_cast<real_t>(deg));
    opt.setInt("quB", reg);
    A.setOptions(opt);
    A.setIntegrationElements(mb);
    if (useFactory)
        A.setQuadratureFactory(gsHPatchRule<real_t>::makeFactory(mb, deg, reg, false));
    auto G = A.getMap(mp);
    auto u = A.getSpace(mb);
    A.initSystem();
    if (form == 0) A.assemble( u * u.tr() * meas(G) );
    else           A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G) );
    return A.matrix().toDense();
}

// The cache stores the solved shape of every (leaf, direction) lookup under the key
// (min(L,1), normalised even target knot vector); a lookup either hits an entry or solves
// and inserts one, so hits + solves is the number of lookups and size == solves. A leaf
// direction that ties with Gauss never touches the cache, so for the target (2p, p-1) with
// simple knots the lookups of one build are the (leaf, direction) pairs at least two
// elements wide. Here these are counted from outside on makeRefinedTHB(p): the 1x7 leaf has
// one (direction 1), the 7x8 leaf two, the 4x4 corner leaf two. The 1x7 leaf's direction
// and the 7x8 leaf's direction 0 are 7 elements over [0.125,1] (L = 0.875) and share a key,
// the corner leaf's two directions share a key, and the 8-element direction has its own.
// Rules built through the cache must equal rules built without it bit for bit (a hit
// returns what a fresh solve of that key returns), a rule with another regularity has other
// knot multiplicities and so needs entries of its own, and clear() forgets everything.
TEST(hpatchrule_cache_counts_over_two_builds)
{
    struct Row { index_t p; index_t patchDirs; index_t solvesA; index_t hitsA;
                 index_t solvesD; index_t hitsD; };
    const Row rows[] = { {2, 5, 3, 2, 6, 9},
                         {3, 5, 3, 2, 6, 9},
                         {4, 5, 3, 2, 6, 9} };

    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        const index_t p = rows[k].p;
        gsInfo << "row p=" << p << "\n";
        gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
        auto dom = thb.domain();

        index_t patchDirs = 0;
        for (auto leaf = thb.tree().beginLeafIterator(); leaf.good(); leaf.next())
            for (short_t dd = 0; dd != 2; ++dd)
                if (leaf.upperCorner()[dd] - leaf.lowerCorner()[dd] >= 2)
                    ++patchDirs;
        CHECK_EQUAL(rows[k].patchDirs, patchDirs);

        // Build A: every lookup is new except for the repeated keys.
        gsHPatchRuleCache<real_t> cache;
        gsHPatchRule<real_t> ruleA(*dom, 2*p, p-1, false, -1, &cache);
        CHECK_EQUAL(patchDirs, cache.numHits() + cache.numSolves());
        CHECK_EQUAL(rows[k].solvesA, cache.numSolves());
        CHECK_EQUAL(rows[k].hitsA,   cache.numHits());
        CHECK_EQUAL(cache.numSolves(), cache.size());
        CHECK_EQUAL(index_t(0), cache.numFailures());
        CHECK(cache.failedKeys().empty());
        const index_t solvesA = cache.numSolves();
        const index_t hitsA   = cache.numHits();

        // Build B on the same cache: all keys are present, so nothing is solved.
        gsHPatchRule<real_t> ruleB(*dom, 2*p, p-1, false, -1, &cache);
        CHECK_EQUAL(solvesA, cache.numSolves());
        CHECK_EQUAL(hitsA + patchDirs, cache.numHits());
        CHECK_EQUAL(cache.numSolves(), cache.size());
        CHECK_EQUAL(index_t(0), cache.numFailures());
        CHECK(cache.failedKeys().empty());

        // Rule C has its own cache, local to the build. A is built from fresh solves
        // (and hits on shapes solved in the same build), B only from stored shapes.
        gsHPatchRule<real_t> ruleC(*dom, 2*p, p-1, false);
        index_t nElem = 0;
        CHECK_EQUAL(real_t(0), hpatchMaxDiff(ruleA, ruleC, thb, -1, nElem));
        CHECK(nElem > 0);
        CHECK_EQUAL(real_t(0), hpatchMaxDiff(ruleB, ruleC, thb, -1, nElem));
        CHECK(nElem > 0);

        // Build D with regularity p-2: other knot multiplicities in the target, so other keys.
        gsHPatchRule<real_t> ruleD(*dom, 2*p, p-2, false, -1, &cache);
        CHECK(cache.numSolves() > solvesA);
        CHECK_EQUAL(rows[k].solvesD, cache.numSolves());
        CHECK_EQUAL(rows[k].hitsD,   cache.numHits());
        CHECK_EQUAL(cache.numSolves(), cache.size());
        CHECK_EQUAL(index_t(0), cache.numFailures());
        CHECK(cache.failedKeys().empty());
        gsHPatchRule<real_t> ruleD0(*dom, 2*p, p-2, false);
        CHECK_EQUAL(real_t(0), hpatchMaxDiff(ruleD, ruleD0, thb, -1, nElem));
        CHECK(nElem > 0);

        // clear() drops the shapes: a rebuild solves exactly what A solved.
        cache.clear();
        CHECK_EQUAL(index_t(0), cache.numHits());
        CHECK_EQUAL(index_t(0), cache.numSolves());
        CHECK_EQUAL(index_t(0), cache.numFailures());
        CHECK_EQUAL(index_t(0), cache.size());
        CHECK(cache.failedKeys().empty());
        gsHPatchRule<real_t> ruleE(*dom, 2*p, p-1, false, -1, &cache);
        CHECK_EQUAL(solvesA, cache.numSolves());
        CHECK_EQUAL(hitsA,   cache.numHits());
        CHECK_EQUAL(real_t(0), hpatchMaxDiff(ruleE, ruleC, thb, -1, nElem));
        CHECK(nElem > 0);
    }
}

// The key holds min(L,1), not L: a direction of length 1 or more is solved on the same
// reference interval as one of length 1, so equal normalised shapes share one entry,
// while a shorter direction has a different stopping tolerance (1e-14/L) and so is
// keyed by its own length. A single 1D leaf with 4 elements and the target (4,1) has
// dim = 5 + 3*3 = 14, i.e. 7 nodes, below the 4*3 = 12 Gauss nodes, so each rule makes
// exactly one lookup. The lengths are powers of two, so that mapping the knots to [0,1]
// is exact and no rounding can split a key.
TEST(hpatchrule_cache_key_clamped_length)
{
    struct Row { real_t L; index_t solves; index_t hits; };
    // L = 2: first shape. L = 4 and L = 1: min(L,1) = 1 and the same normalised knots,
    // so hits. L = 0.5 and L = 0.25: other keys, so solves.
    const Row rows[] = { {2,    1, 0},
                         {4,    1, 1},
                         {0.5,  2, 1},
                         {0.25, 3, 1},
                         {1,    3, 2} };

    const index_t p = 2;
    gsHPatchRuleCache<real_t> cache;
    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        gsInfo << "row L=" << rows[k].L << "\n";
        gsKnotVector<real_t> kv(0, rows[k].L, 3, p + 1, 1);
        gsTensorBSplineBasis<1,real_t> tb(kv);
        gsTHBSplineBasis<1,real_t> thb(tb);
        gsHPatchRule<real_t> rule(*thb.domain(), 2*p, 1, false, -1, &cache);
        CHECK_EQUAL(rows[k].solves, cache.numSolves());
        CHECK_EQUAL(rows[k].hits,   cache.numHits());
        CHECK_EQUAL(cache.numSolves(), cache.size());
    }
    CHECK_EQUAL(index_t(3), cache.size());
    CHECK_EQUAL(index_t(0), cache.numFailures());

    // The rule that hit the entry of L = 2 equals one solved without a cache.
    gsKnotVector<real_t> kv4(0, 4, 3, p + 1, 1);
    gsTensorBSplineBasis<1,real_t> tb4(kv4);
    gsTHBSplineBasis<1,real_t> thb4(tb4);
    gsHPatchRule<real_t> cached4(*thb4.domain(), 2*p, 1, false, -1, &cache);
    gsHPatchRule<real_t> plain4(*thb4.domain(), 2*p, 1, false);
    index_t nElem = 0;
    CHECK_EQUAL(real_t(0), hpatchMaxDiff(cached4, plain4, thb4, -1, nElem));
    CHECK_EQUAL(index_t(4), nElem);
}

// A direction whose moment solve fails is recorded once: the key (normalised even target
// knot vector on [0,1], min(L,1)) goes into failedKeys() and numFailures() counts it. A
// failed entry stays in the cache, so a second build hits it without a new failure and the
// leaf again falls back to Gauss in that direction, which numFallbackCells() counts. The
// fixture is the graded leaf of hpatchrule_fallback_newton_fails_in_one_direction_2d:
// direction 0 has breakpoints 0, 1/s, 101/s, 1 with s = 1 + 100 + 1e4 and fails, direction
// 1 is uniform on [0,1] and converges.
TEST(hpatchrule_cache_failed_keys)
{
    const index_t p = 2, reg = 1;
    const real_t s = 1 + 100 + 1e4, b1 = 1 / s, b2 = 101 / s;
    const real_t knots0[] = {0,0,0, b1, b2, 1,1,1};
    gsKnotVector<real_t> kv0(std::vector<real_t>(knots0, knots0 + 8), p);
    gsKnotVector<real_t> kv1(0, 1, 3, p + 1, 1);
    gsTensorBSplineBasis<2,real_t> tb(kv0, kv1);
    gsTHBSplineBasis<2,real_t> thb(tb);
    auto dom = thb.domain();

    gsHPatchRuleCache<real_t> cache;
    gsHPatchRule<real_t> rule1(*dom, 2*p, reg, false, -1, &cache);
    CHECK_EQUAL(index_t(1),  rule1.numFallbackCells());
    CHECK_EQUAL(index_t(63), rule1.numPointsTotal());
    CHECK_EQUAL(index_t(2), cache.numSolves());
    CHECK_EQUAL(index_t(0), cache.numHits());
    CHECK_EQUAL(index_t(1), cache.numFailures());
    CHECK_EQUAL(cache.numSolves(), cache.size());
    CHECK_EQUAL(size_t(1), cache.failedKeys().size());
    if (cache.failedKeys().size() == 1)
    {
        const gsKnotVector<real_t> & key = cache.failedKeys()[0].first;
        CHECK_EQUAL(real_t(1), cache.failedKeys()[0].second);
        CHECK_EQUAL(2*p, key.degree());
        CHECK_EQUAL(real_t(0), key[0]);
        CHECK_EQUAL(real_t(1), key[key.size() - 1]);
    }
    const index_t failures1 = cache.numFailures();

    // The failed entry is hit, not solved again, and the fallback is the same.
    gsHPatchRule<real_t> rule2(*dom, 2*p, reg, false, -1, &cache);
    CHECK_EQUAL(failures1, cache.numFailures());
    CHECK_EQUAL(size_t(1), cache.failedKeys().size());
    CHECK_EQUAL(index_t(2), cache.numHits());
    CHECK_EQUAL(index_t(1), rule2.numFallbackCells());
    index_t nElem = 0;
    CHECK_EQUAL(real_t(0), hpatchMaxDiff(rule1, rule2, thb, -1, nElem));
    CHECK(nElem > 0);

    cache.clear();
    CHECK(cache.failedKeys().empty());
    CHECK_EQUAL(index_t(0), cache.numFailures());
}

// A copy of a rule gives what a freshly built rule gives, for every element and in any
// query order, also when the copy is assigned onto a rule of another domain and when the
// original is gone before the first query. Each rule is queried in a different order and
// the queries are interleaved element by element. What this can show is correct results
// under interleaving, after assignment from another domain and after the original is
// destroyed (the lookup is right from any start cell). That the leaf data is shared and
// each handle keeps its own lookup position is not observable through the public API.
TEST(hpatchrule_handle_copies_interleaved)
{
    const index_t p = 2;
    gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p);
    std::vector<index_t> more = {2, 12,12,20,20,   2, 24,24,32,32,   2, 0,24,4,32};
    thb.refineElements(more);
    auto dom = thb.domain();

    std::vector<gsVector<real_t> > lows, upps;
    for (auto & elem : dom->allElements())
    {
        lows.push_back(elem.lowerCorner());
        upps.push_back(elem.upperCorner());
    }
    const size_t nElem = lows.size();

    // Reference result of every element from a rule that was only ever queried in volume order.
    gsHPatchRule<real_t> fresh(*dom, 2*p, p-1, false);
    std::vector<gsMatrix<real_t> > refNodes(nElem);
    std::vector<gsVector<real_t> > refWeights(nElem);
    for (size_t e = 0; e != nElem; ++e)
        fresh.mapTo(lows[e], upps[e], refNodes[e], refWeights[e]);

    gsHPatchRule<real_t> orig(*dom, 2*p, p-1, false);
    gsHPatchRule<real_t> byCopy(orig);

    // The other domain has other leaves, hence another number of points.
    gsTHBSplineBasis<2,real_t> otherThb = makeRefinedTHB(p);
    auto otherDom = otherThb.domain();
    gsHPatchRule<real_t> byAssign(*otherDom, 2*p, p-1, false);
    CHECK(byAssign.numPointsTotal() != orig.numPointsTotal());
    byAssign = orig;

    // The copy is taken from a rule that is destroyed before anything is queried.
    memory::unique_ptr<gsHPatchRule<real_t> > holder(
        new gsHPatchRule<real_t>(*dom, 2*p, p-1, false));
    gsHPatchRule<real_t> outliving(*holder);
    holder.reset();

    CHECK_EQUAL(orig.numPointsTotal(),   byCopy.numPointsTotal());
    CHECK_EQUAL(orig.numPointsTotal(),   byAssign.numPointsTotal());
    CHECK_EQUAL(orig.numPointsTotal(),   outliving.numPointsTotal());
    CHECK_EQUAL(orig.numFallbackCells(), byCopy.numFallbackCells());
    CHECK_EQUAL(orig.numFallbackCells(), byAssign.numFallbackCells());
    CHECK_EQUAL(orig.numFallbackCells(), outliving.numFallbackCells());

    std::vector<size_t> volume(nElem), reversed(nElem), shuffled(nElem);
    for (size_t e = 0; e != nElem; ++e) { volume[e] = e; reversed[e] = nElem - 1 - e; }
    shuffled = volume;
    std::mt19937 rng(20240607);
    std::shuffle(shuffled.begin(), shuffled.end(), rng);

    gsHPatchRule<real_t> * handles[4] = { &orig, &byCopy, &byAssign, &outliving };
    const std::vector<size_t> * orders[4] = { &volume, &reversed, &shuffled, &volume };
    gsMatrix<real_t> nodes;
    gsVector<real_t> weights;
    for (size_t i = 0; i != nElem; ++i)
        for (size_t h = 0; h != 4; ++h)
        {
            const size_t e = (*orders[h])[i];
            handles[h]->mapTo(lows[e], upps[e], nodes, weights);
            CHECK_EQUAL(refWeights[e].size(), weights.size());
            if (refWeights[e].size() == weights.size() && weights.size() > 0)
            {
                CHECK_EQUAL(real_t(0), (nodes   - refNodes[e]  ).cwiseAbs().maxCoeff());
                CHECK_EQUAL(real_t(0), (weights - refWeights[e]).cwiseAbs().maxCoeff());
            }
        }
}

// The factory returns, for every patch and every slot, the rule that the constructor builds
// for the same (patch, fixDir): fixDir -1 is the volume rule, fixDir d the rule of the two
// sides normal to direction d. Equality is bitwise on the nodes and weights of every volume
// or boundary element, with a caller cache (shared across rows) and without one (the
// factory then owns a cache). The two patches of mb2 have different leaves, so a mixed-up
// patch index changes the point counts. A patch or direction outside the multi-basis is
// rejected by a GISMO_ENSURE, which throws in every build.
TEST(hpatchrule_factory_slots_equal_constructor)
{
    // mb2: patch 0 is the corner-refined THB, patch 1 the 10-leaf one.
    gsTHBSplineBasis<2,real_t> ml = makeRefinedTHB(2);
    std::vector<index_t> more = {2, 12,12,20,20,   2, 24,24,32,32,   2, 0,24,4,32};
    ml.refineElements(more);
    gsMultiBasis<real_t> mb2;
    mb2.addBasis(makeRefinedTHB(2).clone().release());
    mb2.addBasis(ml.clone().release());

    // mb3: one 3D patch, degree 2, refined in a corner.
    gsKnotVector<real_t> kv(0, 1, 0, 3, 1);
    kv.uniformRefine(); kv.uniformRefine();
    gsTensorBSplineBasis<3,real_t> tb(kv, kv, kv);
    gsTHBSplineBasis<3,real_t> thb3(tb);
    thb3.uniformRefine();
    std::vector<index_t> box = {2, 0,0,0, 4,4,4};
    thb3.refineElements(box);
    gsMultiBasis<real_t> mb3;
    mb3.addBasis(thb3.clone().release());

    gsHPatchRuleCache<real_t> callerCache;
    struct Row { const gsMultiBasis<real_t> * mb; index_t deg; index_t reg; bool useCache; };
    const Row rows[] = { {&mb2, 4, 1, true},
                         {&mb2, 4, 0, false},
                         {&mb3, 4, 1, true},
                         {&mb3, 4, 0, false} };

    gsOptionList opt;
    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        const gsMultiBasis<real_t> & mb = *rows[k].mb;
        const short_t d = mb.dim();
        gsHPatchRule<real_t>::QuadratureFactory f = gsHPatchRule<real_t>::makeFactory(
            mb, rows[k].deg, rows[k].reg, false, rows[k].useCache ? &callerCache : nullptr);
        for (index_t patch = 0; patch != static_cast<index_t>(mb.nBases()); ++patch)
            for (short_t fd = -1; fd != d; ++fd)
            {
                auto dom = mb.basis(patch).domain();
                gsHPatchRule<real_t> ref(*dom, rows[k].deg, rows[k].reg, false, fd);
                gsQuadRule<real_t>::uPtr h = f(mb.basis(patch), opt, patch, fd);
                index_t nElem = 0;
                const real_t diff = hpatchMaxDiff(*h, ref, mb.basis(patch), fd, nElem);
                gsInfo << "dim=" << d << " deg=" << rows[k].deg << " reg=" << rows[k].reg
                       << " patch=" << patch << " fd=" << fd << " nElem=" << nElem
                       << (rows[k].useCache ? " cache" : " nocache") << "\n";
                CHECK(nElem > 0);
                CHECK_EQUAL(real_t(0), diff);
            }
    }

    gsHPatchRule<real_t>::QuadratureFactory g =
        gsHPatchRule<real_t>::makeFactory(mb2, 4, 1, false);
    CHECK_THROW(g(mb2.basis(0), opt, 2, -1), std::runtime_error);
    CHECK_THROW(g(mb2.basis(0), opt, 0, 2), std::runtime_error);
}

// A caller's cache sees every build the factory makes: the volume rule of each patch is
// built when the factory is created, a boundary slot on its first request, and later
// requests are served from the stored rule. The lookups of the volume builds are counted
// from outside over both patches: (leaf, direction) pairs at least two elements wide, which
// is where the target (4,1) with simple knots does not tie with Gauss. The boundary slot
// (patch 0, fixDir 0) builds a rule that looks up only direction 1; the three leaves of
// makeRefinedTHB(2) (1x7, 7x8 and 4x4 elements) are all at least 2 wide in direction 1, so
// that adds 3 lookups.
TEST(hpatchrule_factory_cache_counters)
{
    gsTHBSplineBasis<2,real_t> ml = makeRefinedTHB(2);
    std::vector<index_t> more = {2, 12,12,20,20,   2, 24,24,32,32,   2, 0,24,4,32};
    ml.refineElements(more);
    gsTHBSplineBasis<2,real_t> corner = makeRefinedTHB(2);
    gsMultiBasis<real_t> mb2;
    mb2.addBasis(corner.clone().release());
    mb2.addBasis(ml.clone().release());

    index_t lookups = 0;
    const gsTHBSplineBasis<2,real_t> * patches[2] = { &corner, &ml };
    for (size_t b = 0; b != 2; ++b)
        for (auto leaf = patches[b]->tree().beginLeafIterator(); leaf.good(); leaf.next())
            for (short_t dd = 0; dd != 2; ++dd)
                if (leaf.upperCorner()[dd] - leaf.lowerCorner()[dd] >= 2)
                    ++lookups;
    CHECK_EQUAL(index_t(23), lookups);

    gsOptionList opt;
    gsHPatchRuleCache<real_t> c;
    gsHPatchRule<real_t>::QuadratureFactory f =
        gsHPatchRule<real_t>::makeFactory(mb2, 4, 1, false, &c,
                                          gsHPatchRule<real_t>::MergeMode::none);

    // Step 1: the volume rules of both patches.
    CHECK_EQUAL(lookups, c.numHits() + c.numSolves());
    CHECK_EQUAL(index_t(15), c.numHits());
    CHECK_EQUAL(index_t(8),  c.numSolves());
    const index_t afterCreate = c.numHits() + c.numSolves();

    // Step 2: the first request of a boundary slot builds it.
    f(mb2.basis(0), opt, 0, 0);
    CHECK_EQUAL(afterCreate + 3, c.numHits() + c.numSolves());
    CHECK_EQUAL(index_t(18), c.numHits());
    CHECK_EQUAL(index_t(8),  c.numSolves());
    const index_t hits = c.numHits(), solves = c.numSolves(), size = c.size();

    // Step 3: the same slot again and a volume slot are copies of stored rules.
    f(mb2.basis(0), opt, 0, 0);
    f(mb2.basis(0), opt, 0, -1);
    CHECK_EQUAL(hits,   c.numHits());
    CHECK_EQUAL(solves, c.numSolves());
    CHECK_EQUAL(size,   c.size());
    CHECK_EQUAL(index_t(0), c.numFailures());

    // Step 4: invariants.
    CHECK_EQUAL(c.numSolves(), c.size());
    CHECK_EQUAL(index_t(0), c.numFailures());
    CHECK(c.failedKeys().empty());
}

// The factory owns a clone of the multi-basis it was made from: refining the caller's basis
// after makeFactory does not change a boundary slot that is built later, which equals the
// rule of an unrefined twin. A handle stays valid when the factory is released.
TEST(hpatchrule_factory_owns_clone)
{
    gsMultiBasis<real_t> mb;
    mb.addBasis(makeRefinedTHB(2).clone().release());
    gsHPatchRule<real_t>::QuadratureFactory f =
        gsHPatchRule<real_t>::makeFactory(mb, 4, 1, false);

    mb.basis(0).uniformRefine();

    // Slot (0, 1) was never requested, so it is built from the factory's own basis now.
    // A factory that read the caller's refined basis would build a different rule.
    gsOptionList opt;
    gsQuadRule<real_t>::uPtr h = f(mb.basis(0), opt, 0, 1);
    gsTHBSplineBasis<2,real_t> twin = makeRefinedTHB(2);
    auto twinDom = twin.domain();
    gsHPatchRule<real_t> ref(*twinDom, 4, 1, false, 1);
    index_t nElem = 0;
    CHECK_EQUAL(real_t(0), hpatchMaxDiff(*h, ref, twin, 1, nElem));
    CHECK(nElem > 0);

    f = gsHPatchRule<real_t>::QuadratureFactory();
    CHECK_EQUAL(real_t(0), hpatchMaxDiff(*h, ref, twin, 1, nElem));
    CHECK(nElem > 0);
}

// The factory installed on a gsExprAssembler assembles the same matrix as the option path
// (quRule = 4, quA, quB), for the mass matrix with the rule (4,1) and the stiffness matrix
// with (4,0), on two 2D hierarchical bases and a 3D one. On one thread the two matrices are
// bitwise equal (the same rules in the same element order); with several threads the
// assembler adds element contributions in a run-dependent order, which changes the sums by
// rounding only (a few units of the machine epsilon relative to the largest entry), and the
// suite may run threaded. So the bound is a small multiple of the epsilon, about 1e-13 in
// double. Per-element identity of the factory rules is checked by
// hpatchrule_factory_slots_equal_constructor.
TEST(hpatchrule_factory_assembly_equals_option_path)
{
    struct Row { int form; index_t deg; index_t reg; };
    const Row rows[] = { {0, 4, 1}, {1, 4, 0} };
    const real_t tol = 500 * math::limits::epsilon();

    gsTHBSplineBasis<2,real_t> corner = makeRefinedTHB(2);
    gsTHBSplineBasis<2,real_t> ml = makeRefinedTHB(2);
    std::vector<index_t> more = {2, 12,12,20,20,   2, 24,24,32,32,   2, 0,24,4,32};
    ml.refineElements(more);
    gsKnotVector<real_t> kv(0, 1, 0, 3, 1);
    kv.uniformRefine(); kv.uniformRefine();
    gsTensorBSplineBasis<3,real_t> tb(kv, kv, kv);
    gsTHBSplineBasis<3,real_t> thb3(tb);
    thb3.uniformRefine();
    std::vector<index_t> box = {2, 0,0,0, 4,4,4};
    thb3.refineElements(box);

    for (size_t k = 0; k != sizeof(rows)/sizeof(rows[0]); ++k)
    {
        gsInfo << "form=" << rows[k].form << " deg=" << rows[k].deg
               << " reg=" << rows[k].reg << "\n";

        gsMatrix<real_t> Ko = hpatchAssembleLeafRule<2>(corner, rows[k].form, rows[k].deg, rows[k].reg, false);
        gsMatrix<real_t> Kf = hpatchAssembleLeafRule<2>(corner, rows[k].form, rows[k].deg, rows[k].reg, true);
        real_t maxAbs = Ko.cwiseAbs().maxCoeff();
        CHECK(maxAbs > 0);
        CHECK_EQUAL(Ko.rows(), Kf.rows());
        CHECK_EQUAL(Ko.cols(), Kf.cols());
        CHECK((Kf - Ko).cwiseAbs().maxCoeff() <= tol * maxAbs);

        Ko = hpatchAssembleLeafRule<2>(ml, rows[k].form, rows[k].deg, rows[k].reg, false);
        Kf = hpatchAssembleLeafRule<2>(ml, rows[k].form, rows[k].deg, rows[k].reg, true);
        maxAbs = Ko.cwiseAbs().maxCoeff();
        CHECK(maxAbs > 0);
        CHECK_EQUAL(Ko.rows(), Kf.rows());
        CHECK_EQUAL(Ko.cols(), Kf.cols());
        CHECK((Kf - Ko).cwiseAbs().maxCoeff() <= tol * maxAbs);

        Ko = hpatchAssembleLeafRule<3>(thb3, rows[k].form, rows[k].deg, rows[k].reg, false);
        Kf = hpatchAssembleLeafRule<3>(thb3, rows[k].form, rows[k].deg, rows[k].reg, true);
        maxAbs = Ko.cwiseAbs().maxCoeff();
        CHECK(maxAbs > 0);
        CHECK_EQUAL(Ko.rows(), Kf.rows());
        CHECK_EQUAL(Ko.cols(), Kf.cols());
        CHECK((Kf - Ko).cwiseAbs().maxCoeff() <= tol * maxAbs);
    }
}

} // SUITE gsQuadratureRules_test
