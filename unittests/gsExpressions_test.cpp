/** @file gsExpressions_test.cpp

    @brief Tests for gsExpressions

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include "gismo_unittest.h"


namespace
{

// Coefficient read of the reference evaluation loops below.
inline real_t refCoef(const gsMatrix<real_t> & Sv, index_t ii) { return Sv.at(ii); }

// Agreement up to a relative 1e-12, not bitwise: where the compiler contracts
// a*b+c into an FMA (clang on arm64 by default) the batched and the per-dof
// loops may round differently. A NaN never agrees.
const real_t lookupTol = 1e-12;

bool nearlyEqual(real_t a, real_t b)
{
    const real_t scale = math::max((real_t)1, math::max(math::abs(a), math::abs(b)));
    return math::abs(a - b) <= lookupTol * scale;
}

// Shape included
bool nearlyEqual(const gsMatrix<real_t> & a, const gsMatrix<real_t> & b)
{
    if (a.rows()!=b.rows() || a.cols()!=b.cols())
        return false;
    for (index_t i = 0; i != a.size(); ++i)
        if (!nearlyEqual(a.at(i), b.at(i)))
            return false;
    return true;
}

// Reference evaluations of the solution expressions: one gsDofMapper::index
// lookup per active function, evaluated on the data of the space u.

gsMatrix<real_t> refValue(const expr::gsFeSpace<real_t> & u, const gsMatrix<real_t> & Sv, index_t k)
{
    gsMatrix<real_t> res;
    const gsDofMapper & map = u.mapper();
    auto & act = u.data().actives.col(1 == u.data().actives.cols() ? 0:k );
    res.setZero(u.dim(), 1);
    for (index_t c = 0; c!=u.dim(); c++) // for all components
    {
        for (index_t i = 0; i!=u.data().actives.rows(); ++i)
        {
            const index_t ii = map.index( act[i], u.data().patchId, c);
            if ( map.is_free_index(ii) ) // DoF value is in the solVector
                res.at(c) += refCoef(Sv, ii) * u.data().values[0](i,k);
            else
                res.at(c) += u.data().values[0](i,k) *
                    u.fixedPart().at( map.global_to_bindex(ii) );
        }
    }
    return res;
}

gsMatrix<real_t> refGrad(const expr::gsFeSpace<real_t> & u, const gsMatrix<real_t> & Sv, index_t k)
{
    gsMatrix<real_t> res;
    const index_t parDim = u.source().domainDim();
    const gsDofMapper & map = u.mapper();
    auto & act = u.data().actives.col(1 == u.data().actives.cols() ? 0:k );
    res.setZero(u.dim(), parDim);
    for (index_t c = 0; c!= u.dim(); c++)
    {
        for (index_t i = 0; i!=u.data().actives.rows(); ++i)
        {
            const index_t ii = map.index(act[i], u.data().patchId, c);
            if ( map.is_free_index(ii) ) // DoF value is in the solVector
            {
                res.row(c) += refCoef(Sv, ii) *
                    u.data().values[1].col(k).segment(i*parDim, parDim).transpose();
            }
            else
            {
                res.row(c) +=
                    u.fixedPart().at( map.global_to_bindex(ii) ) *
                    u.data().values[1].col(k).segment(i*parDim, parDim).transpose();
            }
        }
    }
    return res;
}

gsMatrix<real_t> refLapl(const expr::gsFeSpace<real_t> & u, const gsMatrix<real_t> & Sv, index_t k)
{
    gsMatrix<real_t> res;
    const index_t parDim = u.source().domainDim();
    const gsDofMapper & map = u.mapper();
    res.setZero(u.dim(), 1); //  scalar, but per component

    index_t numActs = u.data().values[0].rows();
    index_t numDers = parDim * (parDim + 1) / 2;
    gsMatrix<real_t> deriv2;

    auto & act = u.data().actives.col(1 == u.data().actives.cols() ? 0:k );
    for (index_t c = 0; c!= u.dim(); c++)
        for (index_t i = 0; i!=numActs; ++i)
        {
            const index_t ii = map.index(act[i], u.data().patchId, c);
            deriv2 = u.data().values[2].block(i*numDers,k,parDim,1);
            if ( map.is_free_index(ii) ) // DoF value is in the solVector
                res.at(c) += refCoef(Sv, ii) * deriv2.sum();
            else
                res.at(c) +=u.fixedPart().at( map.global_to_bindex(ii) ) * deriv2.sum();
        }
    return res;
}

gsMatrix<real_t> refHess(const expr::gsFeSpace<real_t> & u, const gsMatrix<real_t> & Sv, index_t k)
{
    gsMatrix<real_t> deriv2, res;
    const gsDofMapper & map = u.mapper();
    const index_t numActs = u.data().values[0].rows();
    const index_t pdim = u.source().domainDim();
    index_t numDers = pdim*(pdim+1)/2;
    auto & act = u.data().actives.col(1 == u.data().actives.cols() ? 0:k );

    if (1==u.dim())
    {
        res.setZero(numDers,1);
        for (index_t i = 0; i!=numActs; ++i)
        {
            const index_t ii = map.index(act[i], u.data().patchId, 0);
            deriv2 = u.data().values[2].block(i*numDers,k,numDers,1);
            if ( map.is_free_index(ii) ) // DoF value is in the solVector
                res += refCoef(Sv, ii) * deriv2;
            else
                res +=u.fixedPart().at( map.global_to_bindex(ii) ) * deriv2;
        }
        expr::secDerToHessian(res, pdim, deriv2);
        res.swap(deriv2);
        res.resize(pdim,pdim);
    }
    else
    {
        res.setZero(u.dim(), numDers);
        for (index_t c = 0; c != u.dim(); c++)
            for (index_t i = 0; i != numActs; ++i)
            {
                const index_t ii = map.index(act[i], u.data().patchId, c);
                deriv2 = u.data().values[2].block(i * numDers, k, numDers,
                                                  1).transpose(); // start row, start col, rows, cols
                if (map.is_free_index(ii)) // DoF value is in the solVector
                    res.row(c) += refCoef(Sv, ii) * deriv2;
                else
                    res.row(c) += u.fixedPart().at(map.global_to_bindex(ii)) * deriv2;
            }
    }
    return res;
}

// Counters of one lookup case
struct LookupCounts
{
    bool storageOk = false;
    index_t nCompared = 0, nMismatch = 0;
    index_t nEqFirstDiffAct = 0;  ///< consecutive elements: same patch, equal first active, different actives
    index_t nSameActDiffPatch = 0;///< consecutive elements: identical actives on different patches
    index_t nFixedUsed = 0;       ///< (active, component) pairs that refer to an eliminated dof
    index_t nElCompared = 0, nElMismatch = 0;
};

// Quadrature of one element, as in gsExprEvaluator::compute_impl
template<class Fn>
void elementWiseReference(gsExprHelper<real_t> & h, gsExprEvaluator<real_t> & ev, Fn ref,
                          std::vector<real_t> & out)
{
    out.assign(h.domain().numElements(), 0.0);
    gsQuadRule<real_t>::uPtr rule;
    index_t patch = -1;
    for (auto & elem : h.domain().allElements())
    {
        if (patch != elem.patchIndex())
        {
            patch = elem.patchIndex();
            rule = gsQuadrature::getPtr(*h.domain().subdomain(patch), ev.options());
        }
        rule->mapTo(elem.lowerCorner(), elem.upperCorner(), h.points(), h.weights());
        h.precompute(patch);
        real_t elVal = 0;
        for (index_t k = 0; k != h.weights().rows(); ++k)
            elVal += h.weights()[k] * ref(k);
        out[elem.id()] = elVal;
    }
}

// Compares ev.integralElWise(e) (OpenMP) with the serial reference, element by
// element. The reference reads the data of s.space(); e must be built from s.
template<class E, class Fn>
void compareElementWise(gsExprAssembler<real_t> & A, gsExprEvaluator<real_t> & ev,
                        const expr::_expr<E> & e, const expr::gsFeSolution<real_t> & s,
                        Fn ref, LookupCounts & c)
{
    ev.integralElWise(e);
    const std::vector<real_t> par = ev.elementwise();
    gsExprHelper<real_t> & h = *A.exprData();
    h.parse(e, s); // registers s.space(), which the reference reads
    h.activateFlags(SAME_ELEMENT);
    std::vector<real_t> ser;
    elementWiseReference(h, ev, ref, ser);
    c.nElCompared += (index_t)ser.size();
    if (par.size() != ser.size())
    {
        ++c.nElMismatch;
        return;
    }
    for (size_t i = 0; i != ser.size(); ++i)
        if (!nearlyEqual(par[i], ser[i]))
            ++c.nElMismatch;
}

// Serial loop over all elements through the shared helper; every
// evaluation point of s, grad(s), lapl(s) and hess(s) is compared with the
// reference. With sameElement the actives have one column per element.
template<class S, class G, class L, class H>
void comparePointWise(gsExprAssembler<real_t> & A, gsExprEvaluator<real_t> & ev,
                      const expr::gsFeSpace<real_t> & u, const gsMatrix<real_t> & Sv,
                      const S & s, const G & gs, const L & ls, const H & hs,
                      bool sameElement, LookupCounts & c)
{
    gsExprHelper<real_t> & h = *A.exprData();
    h.parse(s, gs, ls, hs);
    if (sameElement) h.activateFlags(SAME_ELEMENT);

    gsQuadRule<real_t>::uPtr rule;
    index_t patch = -1;
    gsMatrix<index_t> prevAct;
    index_t prevPatch = -1;
    for (auto & elem : h.domain().allElements())
    {
        if (patch != elem.patchIndex())
        {
            patch = elem.patchIndex();
            rule = gsQuadrature::getPtr(*h.domain().subdomain(patch), ev.options());
        }
        rule->mapTo(elem.lowerCorner(), elem.upperCorner(), h.points(), h.weights());
        h.precompute(patch);

        const gsMatrix<index_t> & act = u.data().actives;
        if (sameElement)
        {
            if (prevPatch >= 0)
            {
                const bool same = prevAct.rows() == act.rows() &&
                    (prevAct.col(0).array() == act.col(0).array()).all();
                if (prevPatch == patch && prevAct(0,0) == act(0,0) && !same)
                    ++c.nEqFirstDiffAct;
                if (prevPatch != patch && same)
                    ++c.nSameActDiffPatch;
            }
            prevAct = act.col(0);
            prevPatch = patch;
        }

        const gsDofMapper & map = u.mapper();
        for (index_t k = 0; k != h.weights().rows(); ++k)
        {
            const index_t col = (1 == act.cols() ? 0 : k);
            for (index_t i = 0; i != act.rows(); ++i)
                for (index_t d = 0; d != u.dim(); ++d)
                    if (!map.is_free_index(map.index(act(i,col), patch, d)))
                        ++c.nFixedUsed;

            c.nCompared += 4;
            if (!nearlyEqual(s.eval(k),  refValue(u, Sv, k))) ++c.nMismatch;
            if (!nearlyEqual(gs.eval(k), refGrad(u, Sv, k)))  ++c.nMismatch;
            if (!nearlyEqual(ls.eval(k), refLapl(u, Sv, k)))  ++c.nMismatch;
            if (!nearlyEqual(hs.eval(k), refHess(u, Sv, k)))  ++c.nMismatch;
        }
    }
}

void fillSolution(const gsDofMapper & map, gsMatrix<real_t> & Sv)
{
    Sv.resize(map.freeSize(), 1);
    for (index_t i = 0; i != Sv.rows(); ++i)
        Sv(i,0) = std::sin(1.0 + i);
}

// One case: Dirichlet setup with the given mapper storage, then the
// point-wise comparison (with and without SAME_ELEMENT) and, for scalar
// spaces, the element-wise integrals of s and lapl(s).
LookupCounts runLookupCase(const gsMultiBasis<real_t> & mb,
                           const gsBoundaryConditions<real_t> & bc,
                           index_t ncomp, gsDofMapper::storage st)
{
    LookupCounts c;
    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(mb);
    gsExprEvaluator<real_t> ev(A);
    auto u = A.getSpace(mb, ncomp);
    u.setMapperStorage(st);
    u.setup(bc, dirichlet::l2Projection, 0);
    c.storageOk = (u.mapper().storageMode() == st);

    gsMatrix<real_t> Sv;
    fillSolution(u.mapper(), Sv);
    auto s  = A.getSolution(u, Sv);
    auto gs = grad(s);
    auto ls = lapl(s);
    auto hs = hess(s);

    comparePointWise(A, ev, s.space(), Sv, s, gs, ls, hs, true,  c);
    comparePointWise(A, ev, s.space(), Sv, s, gs, ls, hs, false, c);

    if (1 == ncomp)
    {
        compareElementWise(A, ev, s.val(), s,
                           [&](index_t k) { return refValue(s.space(), Sv, k)(0,0); }, c);
        compareElementWise(A, ev, ls.val(), s,
                           [&](index_t k) { return refLapl(s.space(), Sv, k)(0,0); }, c);
    }
    return c;
}

struct ResetupCounts
{
    bool storageOk = false;
    index_t nCompared = 0, nMismatch = 0;
    bool sameActives = false;    ///< actives at the point identical before/after
    bool indexChanged = false;   ///< some global index differs between the mappers
    bool valueChanged = false;   ///< the reference value differs between the mappers
    index_t nElCompared = 0, nElMismatch = 0;
};

// The same expression objects are evaluated, the space is set up again with
// other boundary conditions, and they are evaluated again at the same point.
// With reparse the second evaluation goes through gsExprEvaluator::eval
// (which parses); otherwise through eval(0) directly, without any parse.
ResetupCounts runResetupCase(const gsMultiBasis<real_t> & mb,
                             const gsBoundaryConditions<real_t> & bc1,
                             const gsBoundaryConditions<real_t> & bc2,
                             gsDofMapper::storage st, bool reparse)
{
    ResetupCounts c;
    gsExprAssembler<real_t> A(1,1);
    A.setIntegrationElements(mb);
    gsExprEvaluator<real_t> ev(A);
    auto u = A.getSpace(mb, 1);
    u.setMapperStorage(st);
    u.setup(bc1, dirichlet::l2Projection, 0);
    c.storageOk = (u.mapper().storageMode() == st);

    gsMatrix<real_t> Sv;
    fillSolution(u.mapper(), Sv);
    auto s  = A.getSolution(u, Sv);
    auto gs = grad(s);
    auto ls = lapl(s);
    auto hs = hess(s);

    const expr::gsFeSpace<real_t> & us = s.space();
    gsVector<real_t> pt(2);
    pt << 0.05, 0.05;

    gsMatrix<real_t> v = ev.eval(s, pt, 0);
    gsMatrix<real_t> g = ev.eval(gs, pt, 0);
    gsMatrix<real_t> l = ev.eval(ls, pt, 0);
    gsMatrix<real_t> h = ev.eval(hs, pt, 0);
    const gsMatrix<index_t> act0 = us.data().actives.col(0);
    const gsMatrix<real_t> refV0 = refValue(us, Sv, 0);
    gsMatrix<index_t> idx0(act0.rows(), 1);
    for (index_t i = 0; i != act0.rows(); ++i)
        idx0(i,0) = u.mapper().index(act0(i,0), 0, 0);

    c.nCompared += 4;
    c.nMismatch += !nearlyEqual(v, refV0);
    c.nMismatch += !nearlyEqual(g, refGrad(us, Sv, 0));
    c.nMismatch += !nearlyEqual(l, refLapl(us, Sv, 0));
    c.nMismatch += !nearlyEqual(h, refHess(us, Sv, 0));

    u.setup(bc2, dirichlet::l2Projection, 0);
    fillSolution(u.mapper(), Sv);

    if (reparse)
    {
        v = ev.eval(s, pt, 0);
        g = ev.eval(gs, pt, 0);
        l = ev.eval(ls, pt, 0);
        h = ev.eval(hs, pt, 0);
    }
    else
    {
        v = s.eval(0);
        g = gs.eval(0);
        l = ls.eval(0);
        h = hs.eval(0);
    }

    const gsMatrix<index_t> & act1 = us.data().actives;
    c.sameActives = act1.cols() == 1 && act1.rows() == act0.rows() &&
        (act1.col(0).array() == act0.col(0).array()).all();
    for (index_t i = 0; i != act0.rows() && i != act1.rows(); ++i)
        if (u.mapper().index(act1(i,0), 0, 0) != idx0(i,0))
            c.indexChanged = true;
    c.valueChanged = !nearlyEqual(refV0, refValue(us, Sv, 0));

    c.nCompared += 4;
    c.nMismatch += !nearlyEqual(v, refValue(us, Sv, 0));
    c.nMismatch += !nearlyEqual(g, refGrad(us, Sv, 0));
    c.nMismatch += !nearlyEqual(l, refLapl(us, Sv, 0));
    c.nMismatch += !nearlyEqual(h, refHess(us, Sv, 0));

    LookupCounts el;
    compareElementWise(A, ev, s.val(), s,
                       [&](index_t k) { return refValue(s.space(), Sv, k)(0,0); }, el);
    c.nElCompared = el.nElCompared;
    c.nElMismatch = el.nElMismatch;
    return c;
}

gsMultiBasis<real_t> squareBasis()
{
    gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare());
    gsMultiBasis<real_t> mb(mp);
    mb.degreeElevate(1);
    mb.uniformRefine(2);
    return mb;
}

} // namespace

SUITE(gsExpressions_test)
{
    gsVector<real_t,2> pt = gsVector<real_t,2>::Constant(0.5);
    gsVector<real_t,2> zero=gsVector<real_t,2>::Zero();

    gsFunctionExpr<real_t>     func1D("x^2 + y^2", 2);
    gsFunctionExpr<real_t>     func2D("x^2","y^2", 2);
    gsMatrix<real_t> ev_func1D = func1D.eval(pt);
    gsMatrix<real_t> ev_func2D = func2D.eval(pt);
    gsMatrix<real_t> der_func1D = func1D.deriv(pt);
    gsMatrix<real_t> der_func2D = func2D.deriv(pt);

    gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare());
    gsMatrix<real_t> ev_mp = mp.patch(0).eval(pt);
    gsMatrix<real_t> der_mp = mp.patch(0).deriv(pt);

    gsMultiBasis<real_t> mb(mp);
    gsMatrix<real_t> ev_mb = mb.basis(0).eval(pt);
    gsMatrix<real_t> der_mb = mb.basis(0).deriv(pt);

    gsMatrix<real_t> solVector = mp.patch(0).coefs().reshape(mp.patch(0).coefs().size(),1);
    gsExprAssembler<real_t> A(1,1);
    auto G = A.getMap(mp);
    auto u = A.getSpace(mb);
    // auto s = A.getSolution(u,solVector);
    auto f1D = A.getCoeff(func1D);
    auto f2D = A.getCoeff(func2D);
    gsExprEvaluator<real_t> ev(A);

    // u.setup();

    TEST(sanity_check)
    {
        CHECK_EQUAL(ev.eval(G,pt), ev_mp);
        CHECK_EQUAL(ev.eval(f1D,pt), ev_func1D);
        CHECK_ARRAY_EQUAL(ev.eval(f2D,pt).col(0), ev_func2D.col(0), 2);
    }
/*
    OPERATORS
*/
    TEST(add_expr)
    {
        // Scalar addition
        CHECK_EQUAL(ev.eval(f1D.val()+0.0,pt), ev_func1D);
        CHECK_EQUAL(ev.eval(0.0+f1D.val(),pt), ev_func1D);

        // Vector addition
        CHECK_EQUAL(ev.eval(f2D+G,pt), ev_func2D+ev_mp);
        CHECK_EQUAL(ev.eval(G+f2D,pt), ev_func2D+ev_mp);

        // NOT WORKING:
        // CHECK_EQUAL(ev.eval(f1D+0.0,pt), pt);
        // CHECK_EQUAL(ev.eval(f2D+zero,pt), pt);
        // CHECK_EQUAL(ev.eval(u+0.0,pt), ev_func1D);
        // CHECK_EQUAL(ev.eval(0.0+u,pt), ev_func1D);
    }

    TEST(sub_expr)
    {
        // Scalar subtraction
        CHECK_EQUAL(ev.eval(f1D.val()-0.0,pt), ev_func1D);
        CHECK_EQUAL(ev.eval(0.0-f1D.val(),pt), -ev_func1D);

        // Vector subtraction
        CHECK_EQUAL(ev.eval(f2D-G,pt), ev_func2D-ev_mp);
        CHECK_EQUAL(ev.eval(G-f2D,pt), ev_mp-ev_func2D);

        // NOT WORKING:
        // CHECK_EQUAL(ev.eval(f1D-0.0,pt), pt);
        // CHECK_EQUAL(ev.eval(f2D-zero,pt), pt);
        // CHECK_EQUAL(ev.eval(u-0.0,pt), ev_func1D);
        // CHECK_EQUAL(ev.eval(0.0-u,pt), -ev_func1D);
    }

    TEST(mult_expr)
    {
        // Scalar multiplication
        CHECK_EQUAL(ev.eval(f1D.val()*0.0,pt), ev_func1D*0.0);
        CHECK_EQUAL(ev.eval(0.0*f1D.val(),pt), 0.0*ev_func1D);

        // Scalar-Vector multiplication
        CHECK_EQUAL((ev.eval(f1D.val()*f2D,pt)-ev_func2D*ev_func1D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(f2D*f1D.val(),pt)-ev_func2D*ev_func1D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(f1D.val()*u,pt)-ev_func1D.value()*ev_mb).norm(), 0.0);
        CHECK_EQUAL((ev.eval(u*f1D.val(),pt)-ev_mb*ev_func1D.value()).norm(), 0.0);

        CHECK_EQUAL((ev.eval(f2D*0.0,pt)-ev_func2D*0.0).norm(), 0.0);
        CHECK_EQUAL((ev.eval(0.0*f2D,pt)-0.0*ev_func2D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(u*0.0,pt)-ev_mb*0.0).norm(), 0.0);
        CHECK_EQUAL((ev.eval(0.0*u,pt)-0.0*ev_mb).norm(), 0.0);

        CHECK_EQUAL((ev.eval(G*0.0,pt)-ev_mp*0.0).norm(), 0.0);
        CHECK_EQUAL((ev.eval(0.0*G,pt)-0.0*ev_mp).norm(), 0.0);
    }

    TEST(div_expr)
    {
        // Scalar division
        CHECK_EQUAL(ev.eval(f1D.val()/1.0,pt).value(), ev_func1D.value()/1.0);
        CHECK_EQUAL(ev.eval(1.0/f1D.val(),pt).value(), 1.0/ev_func1D.value());

        CHECK_EQUAL((ev.eval(u/1.0,pt)-ev_mb/1.0).norm(), 0.0);
        // CHECK_EQUAL((ev.eval(1.0/u,pt)-1.0*ev_mb).norm(), 0.0);

        // Vector division
        CHECK_EQUAL((ev.eval(f2D/1.0,pt)-ev_func2D/1.0).norm(), 0.0);
        CHECK_EQUAL((ev.eval(f2D/f1D.val(),pt)-ev_func2D/ev_func1D.value()).norm(), 0.0);
    }

/*
    MATHEMATICAL FUNCTIONS
*/
    TEST(pow_expr)
    {
        // Scalar power
        CHECK_EQUAL(ev.eval(pow(f1D,2),pt).value(), math::pow(ev_func1D.value(),2.0));
        // CHECK_EQUAL((ev.eval(pow(u,2),pt)-ev_mb.cwiseProduct(ev_mb)).norm(), 0.0);
    }

/*
    DIFFERENTIAL OPERATORS
*/
    TEST(grad_expr)
    {
        const index_t dim = 2;
        const index_t nAct = ev_mb.rows();
        gsMatrix<> grad_func1D = der_func1D.reshape(dim,1).transpose();
        gsMatrix<> grad_func2D = der_func2D.reshape(dim,2).transpose();
        gsMatrix<> grad_mb     = der_mb.reshape(dim,nAct).transpose();

        // Scalar gradient
        CHECK_EQUAL((ev.eval(grad(f1D),pt)    -grad_func1D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(f1D+f1D),pt)-2*grad_func1D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(f1D-f1D),pt)-0*grad_func1D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(2*f1D),pt)  -2*grad_func1D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(f1D*2),pt)  -2*grad_func1D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(f1D/2),pt)  -grad_func1D/2).norm(), 0.0);

        // Vector gradient
        CHECK_EQUAL((ev.eval(grad(f2D),pt)    -grad_func2D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(f2D+f2D),pt)-2*grad_func2D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(f2D-f2D),pt)-0*grad_func2D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(2*f2D),pt)  -2*grad_func2D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(f2D*2),pt)  -2*grad_func2D).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(f2D/2),pt)  -grad_func2D/2).norm(), 0.0);

        // Space gradient
        CHECK_EQUAL((ev.eval(grad(u),pt)-grad_mb).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(u+u),pt)-2*grad_mb).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(u-u),pt)-0*grad_mb).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(2*u),pt)-2*grad_mb).norm(), 0.0);
        CHECK_EQUAL((ev.eval(grad(u*2),pt)-2*grad_mb).norm(), 0.0);

        CHECK_EQUAL((ev.eval(grad(u/2),pt)-grad_mb/2).norm(), 0.0);
    }

    TEST(jac_expr)
    {
        const index_t dim = 2;
        const index_t nAct = ev_mb.rows();
        gsMatrix<> grad_func1D = der_func1D.reshape(dim,1).transpose();
        gsMatrix<> grad_func2D = der_func2D.reshape(dim,2).transpose();
        gsMatrix<> grad_mb     = der_mb.reshape(dim,nAct).transpose();

        // Scalar gradient
        // CHECK_EQUAL((ev.eval(jac(f1D),pt)    -grad_func1D).norm(), 0.0);

        // CHECK_EQUAL((ev.eval(jac(f1D+f1D),pt)-2*grad_func1D).norm(), 0.0);
        // CHECK_EQUAL((ev.eval(jac(f1D-f1D),pt)-0*grad_func1D).norm(), 0.0);

        // CHECK_EQUAL((ev.eval(jac(2*f1D),pt)  -2*grad_func1D).norm(), 0.0);
        // CHECK_EQUAL((ev.eval(jac(f1D*2),pt)  -2*grad_func1D).norm(), 0.0);

        // CHECK_EQUAL((ev.eval(jac(f1D/2),pt)  -grad_func1D/2).norm(), 0.0);

    }

/*
    MIN/MAX REDUCTIONS
*/
    // Pins gsExprEvaluator::max()/min() against an everywhere-negative /
    // everywhere-positive constant field, whose true extremum is exact
    // regardless of quadrature (a constant integrand needs no sampling
    // accuracy). max_op::init() must start below every attainable value;
    // starting at the smallest positive normal double instead (the historic
    // defect) would leave that spurious positive seed as the reported
    // maximum of an always-negative field, so a positive returned value is
    // the unmistakable fingerprint of the bug, checked here via CHECK(<0)
    // in addition to the value.
    TEST(minmax_reduction_seed)
    {
        gsFunctionExpr<real_t> negConst("-1.0", 2);
        gsFunctionExpr<real_t> posConst("1.0", 2);

        gsExprAssembler<real_t> Aloc(1,1);
        // max()/min() iterate the integration domain, unlike the pointwise
        // ev.eval() used elsewhere in this file, so the elements must be set.
        Aloc.setIntegrationElements(mb);
        Aloc.getMap(mp);
        Aloc.getSpace(mb);
        auto negC = Aloc.getCoeff(negConst);
        auto posC = Aloc.getCoeff(posConst);
        gsExprEvaluator<real_t> evloc(Aloc);

        const real_t tol = 1e2 * std::numeric_limits<real_t>::epsilon();

        real_t maxNeg = evloc.max(negC.val());
        CHECK(maxNeg < 0.0);
        CHECK_CLOSE(-1.0, maxNeg, tol);

        // mirror case: min_op::init() must start above every attainable
        // value, pinning the symmetric defect for min() on a positive field
        real_t minPos = evloc.min(posC.val());
        CHECK(minPos > 0.0);
        CHECK_CLOSE(1.0, minPos, tol);
    }

    // Only meaningful with OpenMP: without it there is one accumulator, the
    // merge under test never happens, and the test would pass while covering
    // nothing. Omitted outright rather than left to degenerate silently.
#ifdef _OPENMP

    // Exercises the merge of the per-thread accumulators in acc_global.
    //
    // compute_impl gives every thread its own thValue and merges each one
    // into m_value exactly once, at the end. A CONSTANT field therefore has
    // no power here at all: every thread arrives with the same number, so
    // any interleaving of the merge -- protected or not -- still yields the
    // correct answer. The field below varies over the domain so that the
    // thread-local extrema genuinely differ, which is the precondition for a
    // lost merge to change the result.
    //
    // The oracle is the same reduction on a single thread, where acc_global
    // runs once and no interleaving exists. Both runs sample identical
    // quadrature points, and min/max is a selection rather than an
    // arithmetic reduction -- it is associative and returns one of the
    // sampled values unchanged -- so the comparison is an exact equality.
    // (That is why byte-identity is legitimate here and is not for a sum,
    // whose floating point result depends on the order of accumulation.)
    //
    // This still cannot prove the absence of a race; only ThreadSanitizer
    // with libarcher measures that. Unlike a constant field, it can fail
    // when one occurs.
    TEST(minmax_reduction_threads)
    {
        // Sign-definite on the unit square (ranges [-2,-1] and [1,2]) so the
        // seed defect above is still caught, but varying across elements.
        gsFunctionExpr<real_t> negRamp("x - 2.0", 2);
        gsFunctionExpr<real_t> posRamp("x + 1.0", 2);

        gsMultiBasis<real_t> mbFine(mp);
        mbFine.uniformRefine(4); // many elements -> more work per thread

        gsExprAssembler<real_t> Aloc(1,1);
        Aloc.setIntegrationElements(mbFine);
        Aloc.getMap(mp);
        Aloc.getSpace(mbFine);
        auto negC = Aloc.getCoeff(negRamp);
        auto posC = Aloc.getCoeff(posRamp);
        gsExprEvaluator<real_t> evloc(Aloc);

        // The reference is taken on one thread, where acc_global runs once
        // and no interleaving exists.
        const int nThreads = omp_get_max_threads();

        omp_set_num_threads(1);
        const real_t maxRef = evloc.max(negC.val());
        const real_t minRef = evloc.min(posC.val());

        CHECK(maxRef < 0.0);
        CHECK(minRef > 0.0);

        // The loop runs on as many threads as the host offers, and never
        // fewer than four. A host or CI runner pinned to one thread would
        // otherwise turn this into a serial check that passes while
        // exercising none of the merge it exists to cover, so the thread
        // count is asserted rather than assumed.
        omp_set_num_threads(nThreads > 4 ? nThreads : 4);
        CHECK(omp_get_max_threads() > 1);

        for (int i = 0; i != 20; ++i)
        {
            CHECK_EQUAL(maxRef, evloc.max(negC.val()));
            CHECK_EQUAL(minRef, evloc.min(posC.val()));
        }

        // Restore: the thread count is process-wide and would otherwise
        // change how every later test in this binary runs.
        omp_set_num_threads(nThreads);
    }

#endif // _OPENMP

    // Test for the positive part of an expression ( Macaulay bracket )
    TEST(ppart_expr)
    {

        gsMatrix<real_t> zero1D = gsMatrix<real_t>::Zero(1,1);
        CHECK_EQUAL( ev.eval(f2D.ppart(), zero) , zero);
        CHECK_EQUAL( ev.eval(f1D.ppart(), zero) , zero1D);

        CHECK_EQUAL( ev.eval((-f2D).ppart(), zero) , zero);
        CHECK_EQUAL( ev.eval((-f1D).ppart(), zero) , zero1D);

        // Test positive part with positive values
        gsMatrix<real_t> pos_result2D = ev_func2D.cwiseMax(0.0);
        gsMatrix<real_t> pos_result1D = gsMatrix<real_t>::Constant(1,1,std::max(ev_func1D.value(),0.0));
        CHECK_EQUAL((ev.eval(f2D.ppart(), pt) - pos_result2D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(f1D.ppart(), pt) - pos_result1D).norm(), 0.0);

        // Test for the negative part of an expression
        CHECK_EQUAL( ev.eval(f2D.npart(), zero) , zero);
        CHECK_EQUAL( ev.eval(f1D.npart(), zero) , zero1D);

        CHECK_EQUAL( ev.eval((-f2D).npart(), pt) , -ev_func2D);
        CHECK_EQUAL( ev.eval((-f1D).npart(), pt) , -ev_func1D);

        // Test negative part with positive values (should be zero)
        gsMatrix<real_t> neg_result2D = (-ev_func2D).cwiseMax(0.0);
        gsMatrix<real_t> neg_result1D = gsMatrix<real_t>::Constant(1,1,std::max(-ev_func1D.value(),0.0));
        CHECK_EQUAL((ev.eval(f2D.npart(), pt) - neg_result2D).norm(), 0.0);
        CHECK_EQUAL((ev.eval(f1D.npart(), pt) - neg_result1D).norm(), 0.0);

    }

/*
    SOLUTION EXPRESSIONS: dof lookup
*/

    TEST(feSolution_lookup_thb)
    {
        gsKnotVector<real_t> kv(0, 1, 3, 3);
        gsTensorBSplineBasis<2,real_t> tbasis(kv, kv);
        gsTHBSplineBasis<2,real_t> basis(tbasis);
        basis.refineElements({1, 2, 2, 6, 6});
        gsMultiBasis<real_t> tmb(basis);

        gsMultiPatch<real_t> geo(*gsNurbsCreator<real_t>::BSplineSquare());
        gsFunctionExpr<real_t> g("sin(x)+y", 2);
        gsBoundaryConditions<real_t> bc;
        bc.setGeoMap(geo);
        bc.addCondition(0, boundary::west,  condition_type::dirichlet, &g, 0, false, -1);
        bc.addCondition(0, boundary::south, condition_type::dirichlet, &g, 0, false, -1);

        for (int st = 0; st != 2; ++st)
        {
            const LookupCounts c = runLookupCase(tmb, bc, 1,
                st ? gsDofMapper::storage::sparse : gsDofMapper::storage::dense);
            CHECK(c.storageOk);
            CHECK(c.nCompared > 0);
            CHECK_EQUAL(0, c.nMismatch);
            CHECK(c.nEqFirstDiffAct > 0);
            CHECK(c.nFixedUsed > 0);
            CHECK(c.nElCompared > 0);
            CHECK_EQUAL(0, c.nElMismatch);
        }
    }

    TEST(feSolution_lookup_multipatch_single_element_patches)
    {
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(4, 4);
        gsMultiBasis<real_t> mpb(mp);
        mpb.degreeElevate(1);

        gsFunctionExpr<real_t> g("sin(x)+y", 2);
        gsBoundaryConditions<real_t> bc;
        bc.setGeoMap(mp);
        bc.addCondition(0,  boundary::west,  condition_type::dirichlet, &g, 0, false, -1);
        bc.addCondition(0,  boundary::south, condition_type::dirichlet, &g, 0, false, -1);
        bc.addCondition(15, boundary::east,  condition_type::dirichlet, &g, 0, false, -1);
        bc.addCondition(15, boundary::north, condition_type::dirichlet, &g, 0, false, -1);

        for (int st = 0; st != 2; ++st)
        {
            const LookupCounts c = runLookupCase(mpb, bc, 1,
                st ? gsDofMapper::storage::sparse : gsDofMapper::storage::dense);
            CHECK(c.storageOk);
            CHECK(c.nCompared > 0);
            CHECK_EQUAL(0, c.nMismatch);
            CHECK(c.nSameActDiffPatch > 0);
            CHECK(c.nFixedUsed > 0);
            CHECK(c.nElCompared > 0);
            CHECK_EQUAL(0, c.nElMismatch);
        }
    }

    TEST(feSolution_lookup_multipatch_refined)
    {
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(4, 4);
        gsMultiBasis<real_t> mpb(mp);
        mpb.degreeElevate(1);
        mpb.uniformRefine(1);

        gsFunctionExpr<real_t> g("sin(x)+y", 2);
        gsBoundaryConditions<real_t> bc;
        bc.setGeoMap(mp);
        bc.addCondition(0,  boundary::west,  condition_type::dirichlet, &g, 0, false, -1);
        bc.addCondition(0,  boundary::south, condition_type::dirichlet, &g, 0, false, -1);
        bc.addCondition(15, boundary::east,  condition_type::dirichlet, &g, 0, false, -1);
        bc.addCondition(15, boundary::north, condition_type::dirichlet, &g, 0, false, -1);

        for (int st = 0; st != 2; ++st)
        {
            const LookupCounts c = runLookupCase(mpb, bc, 1,
                st ? gsDofMapper::storage::sparse : gsDofMapper::storage::dense);
            CHECK(c.storageOk);
            CHECK(c.nCompared > 0);
            CHECK_EQUAL(0, c.nMismatch);
            CHECK(c.nFixedUsed > 0);
            CHECK(c.nElCompared > 0);
            CHECK_EQUAL(0, c.nElMismatch);
        }
    }

    TEST(feSolution_lookup_two_components_dirichlet)
    {
        gsMultiBasis<real_t> smb = squareBasis();
        gsMultiPatch<real_t> geo(*gsNurbsCreator<real_t>::BSplineSquare());
        gsFunctionExpr<real_t> g("sin(x)+y", 2);
        gsBoundaryConditions<real_t> bc;
        bc.setGeoMap(geo);
        bc.addCondition(0, boundary::west,  condition_type::dirichlet, &g, 0, false, 0);
        bc.addCondition(0, boundary::south, condition_type::dirichlet, &g, 0, false, 1);

        for (int st = 0; st != 2; ++st)
        {
            const LookupCounts c = runLookupCase(smb, bc, 2,
                st ? gsDofMapper::storage::sparse : gsDofMapper::storage::dense);
            CHECK(c.storageOk);
            CHECK(c.nCompared > 0);
            CHECK_EQUAL(0, c.nMismatch);
            CHECK(c.nFixedUsed > 0);
        }
    }

    TEST(feSolution_lookup_resetup_with_reparse)
    {
        gsMultiBasis<real_t> smb = squareBasis();
        gsMultiPatch<real_t> geo(*gsNurbsCreator<real_t>::BSplineSquare());
        gsFunctionExpr<real_t> g("sin(x)+y", 2);
        gsBoundaryConditions<real_t> bc1, bc2;
        bc1.setGeoMap(geo);
        bc2.setGeoMap(geo);
        bc1.addCondition(0, boundary::west,  condition_type::dirichlet, &g, 0, false, -1);
        bc2.addCondition(0, boundary::west,  condition_type::dirichlet, &g, 0, false, -1);
        bc2.addCondition(0, boundary::south, condition_type::dirichlet, &g, 0, false, -1);

        for (int st = 0; st != 2; ++st)
        {
            const ResetupCounts c = runResetupCase(smb, bc1, bc2,
                st ? gsDofMapper::storage::sparse : gsDofMapper::storage::dense, true);
            CHECK(c.storageOk);
            CHECK(c.nCompared > 0);
            CHECK_EQUAL(0, c.nMismatch);
            CHECK(c.sameActives);
            CHECK(c.indexChanged);
            CHECK(c.valueChanged);
            CHECK(c.nElCompared > 0);
            CHECK_EQUAL(0, c.nElMismatch);
        }
    }

    TEST(feSolution_lookup_resetup_without_reparse)
    {
        gsMultiBasis<real_t> smb = squareBasis();
        gsMultiPatch<real_t> geo(*gsNurbsCreator<real_t>::BSplineSquare());
        gsFunctionExpr<real_t> g("sin(x)+y", 2);
        gsBoundaryConditions<real_t> bc1, bc2;
        bc1.setGeoMap(geo);
        bc2.setGeoMap(geo);
        bc1.addCondition(0, boundary::west,  condition_type::dirichlet, &g, 0, false, -1);
        bc2.addCondition(0, boundary::west,  condition_type::dirichlet, &g, 0, false, -1);
        bc2.addCondition(0, boundary::south, condition_type::dirichlet, &g, 0, false, -1);

        for (int st = 0; st != 2; ++st)
        {
            const ResetupCounts c = runResetupCase(smb, bc1, bc2,
                st ? gsDofMapper::storage::sparse : gsDofMapper::storage::dense, false);
            CHECK(c.storageOk);
            CHECK(c.nCompared > 0);
            CHECK_EQUAL(0, c.nMismatch);
            CHECK(c.sameActives);
            CHECK(c.indexChanged);
            CHECK(c.valueChanged);
            CHECK(c.nElCompared > 0);
            CHECK_EQUAL(0, c.nElMismatch);
        }
    }

}
