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
    // x = 0.5. With numFallbackCells()==0 no leaf falls back to exact Gauss, so
    // agreement with Gauss below verifies the min(...,m_deg+1) cap branch of the
    // target space. The centred fixture at x = 0.5 follows in
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

} // SUITE gsQuadratureRules_test
