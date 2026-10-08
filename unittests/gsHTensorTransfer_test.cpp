/** @file gsHTensorTransfer_test.cpp

    @brief Exactness of the hierarchical (HB and THB) coefficient transfer
    under refinement and unrefinement.

    If the domains Omega^k of the old basis are contained in those of the
    new basis for every level k, the old spline space is nested in the new
    one and the transfer matrix T is the exact change of basis:
    sum_i c_i B_i^old == sum_j (T c)_j B_j^new pointwise.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
**/

#include "gismo_unittest.h"

using namespace gismo;

SUITE(gsHTensorTransfer_test)
{

typedef std::vector<index_t> Boxes;
typedef gsSparseMatrix<real_t,RowMajor> SpMat;

// Number of level-0 elements per direction (degree 2 everywhere)
template<short_t d> index_t n0() { return d == 2 ? 8 : 4; }

template<short_t d, bool Trunc>
gsTHBSplineBasis<d,real_t,Trunc> tensorStart()
{
    gsKnotVector<> kv(0, 1, n0<d>()-1, 3);
    gsTensorBSplineBasis<d> tens(std::vector< gsKnotVector<> >(d, kv));
    return gsTHBSplineBasis<d,real_t,Trunc>(tens);
}

// Box (level, lower corner, upper corner) covering [0,t]^d, indices in
// the unique-knot space of the given level
template<short_t d>
Boxes cornerBox(index_t level, index_t upper)
{
    Boxes b(1, level);
    b.insert(b.end(), d, 0);
    b.insert(b.end(), d, upper);
    return b;
}

// S0: level-1 box [0,1/2]^d
template<short_t d>
Boxes boxS0() { return cornerBox<d>(1, n0<d>()); }

// S1 = S0 + level-2 box [0,1/4]^d
template<short_t d>
Boxes boxS1() { return cornerBox<d>(2, n0<d>()); }

// OVL: level-1 box [1/4,3/4] x [0,1]^(d-1)
template<short_t d>
Boxes boxOVL()
{
    const index_t n = n0<d>();
    Boxes b(1, 1);
    b.push_back(n/2); b.insert(b.end(), d-1, 0);
    b.push_back(3*n/2); b.insert(b.end(), d-1, 2*n);
    return b;
}

// CAS: level-2 box [1/4,3/4] x [0,1/2]^(d-1) and level-3 box [0,3/8] x [0,1/4]^(d-1)
template<short_t d>
Boxes boxCAS()
{
    const index_t n = n0<d>();
    Boxes b(1, 2);
    b.push_back(n); b.insert(b.end(), d-1, 0);
    b.push_back(3*n); b.insert(b.end(), d-1, 2*n);
    b.push_back(3);
    b.insert(b.end(), d, 0);
    b.push_back(3*n); b.insert(b.end(), d-1, 2*n);
    return b;
}

// UNR: level-0 boxes clearing S0+ back to the S0 domain
template<short_t d> Boxes boxUNR();

template<> Boxes boxUNR<2>()
{
    const index_t n = n0<2>();
    index_t v[] = {0, n/2, 0, 3*n/4, n,   0, n/4, n/2, n/2, n};
    return Boxes(v, v + 10);
}

template<> Boxes boxUNR<3>()
{
    const index_t n = n0<3>();
    index_t v[] = {0, n/2, 0,   0,   3*n/4, n,   n,
                   0, n/4, n/2, 0,   n/2,   n,   n,
                   0, n/4, 0,   n/2, n/2,   n/2, n};
    return Boxes(v, v + 21);
}

template<short_t d, bool Trunc>
gsTHBSplineBasis<d,real_t,Trunc> startS0()
{
    gsTHBSplineBasis<d,real_t,Trunc> b = tensorStart<d,Trunc>();
    b.refineElements(boxS0<d>());
    return b;
}

template<short_t d, bool Trunc>
gsTHBSplineBasis<d,real_t,Trunc> startS1()
{
    gsTHBSplineBasis<d,real_t,Trunc> b = startS0<d,Trunc>();
    b.refineElements(boxS1<d>());
    return b;
}

template<short_t d, bool Trunc>
gsTHBSplineBasis<d,real_t,Trunc> startS0plus()
{
    gsTHBSplineBasis<d,real_t,Trunc> b = startS0<d,Trunc>();
    b.refineElements(boxOVL<d>());
    return b;
}

template<short_t d, bool Trunc>
gsTHBSplineBasis<d,real_t,Trunc> startS3()
{
    gsTHBSplineBasis<d,real_t,Trunc> b = startS1<d,Trunc>();
    b.refineElements(boxCAS<d>());
    return b;
}

// Values of sum_i coefs(i,:) B_i at each column of pts -> (coefs.cols() x pts.cols())
template<short_t d, bool Trunc>
gsMatrix<real_t> evalSpline(const gsTHBSplineBasis<d,real_t,Trunc> & basis,
                            const gsMatrix<real_t> & coefs,
                            const gsMatrix<real_t> & pts)
{
    gsMatrix<real_t> vals;
    gsMatrix<index_t> act;
    basis.eval_into(pts, vals);
    basis.active_into(pts, act);
    gsMatrix<real_t> out = gsMatrix<real_t>::Zero(coefs.cols(), pts.cols());
    for (index_t k = 0; k < pts.cols(); ++k)
        for (index_t j = 0; j < act.rows(); ++j)
            out.col(k) += vals(j,k) * coefs.row(act(j,k)).transpose();
    return out;
}

template<short_t d>
gsMatrix<real_t> randomPoints()
{
    gsMatrix<real_t> pts = gsMatrix<real_t>::Random(d, 256);
    pts.array() = 0.5 + 0.49 * pts.array();
    return pts;
}

template<short_t d, bool Trunc>
void checkReproduction(const gsTHBSplineBasis<d,real_t,Trunc> & oldB,
                       const gsMatrix<real_t> & cOld,
                       const gsTHBSplineBasis<d,real_t,Trunc> & newB,
                       const gsMatrix<real_t> & cNew,
                       const gsMatrix<real_t> & pts)
{
    const real_t maxErr = (evalSpline(oldB, cOld, pts)
                           - evalSpline(newB, cNew, pts)).cwiseAbs().maxCoeff();
    CHECK_CLOSE(0.0, maxErr, 1e-12);
}

// Dimension guards for a transfer matrix; false if it is unusable
template<short_t d, bool Trunc>
bool guardMatrix(const gsTHBSplineBasis<d,real_t,Trunc> & oldB,
                 const gsTHBSplineBasis<d,real_t,Trunc> & newB,
                 const SpMat & T)
{
    CHECK_EQUAL(oldB.size(), T.cols());
    CHECK_EQUAL(newB.size(), T.rows());
    return T.cols() == oldB.size() && T.rows() == newB.size();
}

// Refinement through an entry point returning the transfer matrix
template<short_t d, bool Trunc, class Op>
void refineWithMatrix(const gsTHBSplineBasis<d,real_t,Trunc> & start, Op op)
{
    std::srand(2026);
    gsTHBSplineBasis<d,real_t,Trunc> oldB = start;
    gsTHBSplineBasis<d,real_t,Trunc> newB = start;
    const gsMatrix<real_t> pts = randomPoints<d>();
    const gsMatrix<real_t> c = gsMatrix<real_t>::Random(oldB.size(), 2);
    SpMat T;
    op(newB, T);
    CHECK(newB.size() > oldB.size());
    if (!guardMatrix(oldB, newB, T)) return;
    CHECK(T.nonZeros() > T.cols());
    checkReproduction(oldB, c, newB, gsMatrix<real_t>(T * c), pts);
}

// Refinement through an entry point overwriting the coefficients
template<short_t d, bool Trunc, class Op>
void refineWithCoefs(const gsTHBSplineBasis<d,real_t,Trunc> & start, Op op)
{
    std::srand(2026);
    gsTHBSplineBasis<d,real_t,Trunc> oldB = start;
    gsTHBSplineBasis<d,real_t,Trunc> newB = start;
    const gsMatrix<real_t> pts = randomPoints<d>();
    const gsMatrix<real_t> c = gsMatrix<real_t>::Random(oldB.size(), 2);
    gsMatrix<real_t> cNew = c;
    op(newB, cNew);
    CHECK(newB.size() > oldB.size());
    CHECK_EQUAL(newB.size(), cNew.rows());
    if (newB.size() != cNew.rows()) return;
    checkReproduction(oldB, c, newB, cNew, pts);
}

// Unrefinement: the transfer maps coarse coefficients to fine ones
template<short_t d, bool Trunc, class Op>
void unrefineWithMatrix(const gsTHBSplineBasis<d,real_t,Trunc> & start, Op op)
{
    std::srand(2026);
    gsTHBSplineBasis<d,real_t,Trunc> fine = start;
    gsTHBSplineBasis<d,real_t,Trunc> coarse = start;
    const gsMatrix<real_t> pts = randomPoints<d>();
    SpMat T;
    op(coarse, T);
    CHECK(coarse.size() < fine.size());
    if (!guardMatrix(coarse, fine, T)) return;
    CHECK(T.nonZeros() > T.cols());
    const gsMatrix<real_t> c = gsMatrix<real_t>::Random(coarse.size(), 2);
    checkReproduction(coarse, c, fine, gsMatrix<real_t>(T * c), pts);
}

template<short_t d, bool Trunc>
void case01() { refineWithMatrix<d,Trunc>(startS0<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, SpMat & T)
    { b.refineElements_withTransfer(boxOVL<d>(), T); }); }

template<short_t d, bool Trunc>
void case02() { refineWithMatrix<d,Trunc>(startS0<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, SpMat & T)
    { b.refineElements_withTransfer2(boxOVL<d>(), T); }); }

template<short_t d, bool Trunc>
void case03() { refineWithCoefs<d,Trunc>(startS0<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, gsMatrix<real_t> & c)
    { b.refineElements_withCoefs(c, boxOVL<d>()); }); }

template<short_t d, bool Trunc>
void case04() { refineWithCoefs<d,Trunc>(startS0<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, gsMatrix<real_t> & c)
    { b.refineElements_withCoefs2(c, boxOVL<d>()); }); }

template<short_t d, bool Trunc>
void case05() { refineWithCoefs<d,Trunc>(startS0<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, gsMatrix<real_t> & c)
    {
        gsMatrix<real_t> B(d, 2);
        B.col(0).setZero();  B(0,0) = 0.25;
        B.col(1).setOnes();  B(0,1) = 0.75;
        b.refine_withCoefs(c, B);
    }); }

template<short_t d, bool Trunc>
void case06() { refineWithMatrix<d,Trunc>(startS1<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, SpMat & T)
    { b.refineToLevel_withTransfer(2, T); }); }

template<short_t d, bool Trunc>
void case07() { refineWithMatrix<d,Trunc>(startS0<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, SpMat & T)
    { b.refineCoarsestLevel_withTransfer(T); }); }

template<short_t d, bool Trunc>
void case08() { refineWithCoefs<d,Trunc>(startS0<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, gsMatrix<real_t> & c)
    { b.uniformRefine_withCoefs(c); }); }

template<short_t d, bool Trunc>
void case09() { refineWithMatrix<d,Trunc>(startS1<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, SpMat & T)
    { b.refineElements_withTransfer(boxCAS<d>(), T); }); }

template<short_t d, bool Trunc>
void case10() { refineWithMatrix<d,Trunc>(startS1<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, SpMat & T)
    { b.refineElements_withTransfer2(boxCAS<d>(), T); }); }

template<short_t d, bool Trunc>
void case11() { unrefineWithMatrix<d,Trunc>(startS0plus<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, SpMat & T)
    { b.unrefineElements_withTransfer(boxUNR<d>(), T); }); }

template<short_t d, bool Trunc>
void case12() { unrefineWithMatrix<d,Trunc>(startS3<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, SpMat & T)
    { b.unrefineToLevel_withTransfer(1, T); }); }

template<short_t d, bool Trunc>
void case13() { unrefineWithMatrix<d,Trunc>(startS3<d,Trunc>(),
    [](gsTHBSplineBasis<d,real_t,Trunc> & b, SpMat & T)
    { b.unrefineFinestLevel_withTransfer(T); }); }

#define TRANSFER_TESTS(NAME, D, TR) \
    TEST(NAME##_d##D##_01_refineElements_withTransfer)       { case01<D,TR>(); } \
    TEST(NAME##_d##D##_02_refineElements_withTransfer2)      { case02<D,TR>(); } \
    TEST(NAME##_d##D##_03_refineElements_withCoefs)          { case03<D,TR>(); } \
    TEST(NAME##_d##D##_04_refineElements_withCoefs2)         { case04<D,TR>(); } \
    TEST(NAME##_d##D##_05_refine_withCoefs)                  { case05<D,TR>(); } \
    TEST(NAME##_d##D##_06_refineToLevel_withTransfer)        { case06<D,TR>(); } \
    TEST(NAME##_d##D##_07_refineCoarsestLevel_withTransfer)  { case07<D,TR>(); } \
    TEST(NAME##_d##D##_08_uniformRefine_withCoefs)           { case08<D,TR>(); } \
    TEST(NAME##_d##D##_09_cascade_withTransfer)              { case09<D,TR>(); } \
    TEST(NAME##_d##D##_10_cascade_withTransfer2)             { case10<D,TR>(); } \
    TEST(NAME##_d##D##_11_unrefineElements_withTransfer)     { case11<D,TR>(); } \
    TEST(NAME##_d##D##_12_unrefineToLevel_withTransfer)      { case12<D,TR>(); } \
    TEST(NAME##_d##D##_13_unrefineFinestLevel_withTransfer)  { case13<D,TR>(); }

TRANSFER_TESTS(HB,  2, false)
TRANSFER_TESTS(HB,  3, false)
TRANSFER_TESTS(THB, 2, true)
TRANSFER_TESTS(THB, 3, true)

#undef TRANSFER_TESTS

} // SUITE(gsHTensorTransfer_test)
