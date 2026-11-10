/** @file gsDofMapper_test.cpp

    @brief Tests for gismo::gsDofMapper

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M.M. Hodzelmans
**/

#include "gismo_unittest.h"

using namespace gismo;

namespace {

typedef std::pair<index_t,index_t> dofPair;

gsMatrix<index_t> idx(std::initializer_list<index_t> l)
{
    gsMatrix<index_t> r(static_cast<index_t>(l.size()), 1);
    index_t i = 0;
    for (index_t v : l) r(i++, 0) = v;
    return r;
}

std::string join(const std::vector<dofPair> & v)
{
    std::ostringstream os;
    for (const dofPair & p : v)
        os << "(" << p.first << "," << p.second << ")";
    return os.str();
}

// "off=offset(k)... ps=patchSize(k,0)...;patchSize(k,1)..."
std::string layout(const gsDofMapper & m)
{
    std::ostringstream os;
    os << "off=";
    for (size_t k = 0; k != m.numPatches(); ++k)
        os << (k ? " " : "") << m.offset(k);
    os << " ps=";
    for (index_t c = 0; c != m.numComponents(); ++c)
        for (size_t k = 0; k != m.numPatches(); ++k)
            os << (k ? " " : (c ? ";" : "")) << m.patchSize(k, c);
    return os.str();
}

// index(i,k,c) of every local dof, as "c0k0|c0k1;c1k0|c1k1"
std::string indices(const gsDofMapper & m)
{
    std::ostringstream os;
    for (index_t c = 0; c != m.numComponents(); ++c)
        for (size_t k = 0; k != m.numPatches(); ++k)
        {
            if (c || k) os << (k ? "|" : ";");
            for (index_t i = 0; i != static_cast<index_t>(m.patchSize(k, c)); ++i)
                os << (i ? "," : "") << m.index(i, k, c);
        }
    return os.str();
}

// Counters, layout, index/bindex/coupling of every local dof, and preimages
std::string state(const gsDofMapper & m)
{
    std::ostringstream os;
    m.print(os);
    os << layout(m) << " map=" << m.mapSize() << "\n";
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        os << "c" << c << ":" << m.size(c) << "," << m.freeSize(c) << "," << m.firstIndex(c);
        for (size_t k = 0; k != m.numPatches(); ++k)
            for (index_t i = 0; i != static_cast<index_t>(m.patchSize(k, c)); ++i)
            {
                os << " " << m.index(i, k, c);
                if (m.is_boundary(i, k, c)) os << "b" << m.bindex(i, k, c);
                if (m.is_coupled(i, k, c)) os << "c";
            }
        os << "\n";
    }
    std::vector<dofPair> pre;
    for (index_t gl = 0; gl != m.size(); ++gl)
    {
        m.preImage(gl, pre);
        os << join(pre) << " ";
    }
    return os.str();
}

// indexOnPatch(gl,k) per patch, '-' where gl is not on patch k
std::string onPatches(const gsDofMapper & m, index_t gl)
{
    std::ostringstream os;
    for (size_t k = 0; k != m.numPatches(); ++k)
    {
        index_t loc = -1;
        os << (k ? " " : "") << (m.indexOnPatch(gl, k, loc) ? std::to_string(loc) : std::string("-"));
    }
    return os.str();
}

// Matching and elimination on the component with the given index
typedef std::function<void(gsDofMapper &, index_t)> compOps;

// Builds gsDofMapper(sizes) with ops[c] on component c and compares each
// component with the scalar mapper of column c with ops[c]. Its free dofs
// come after those of the preceding components; its eliminated dofs come
// after all free dofs and the eliminated dofs of the preceding components.
gsDofMapper checkAgainstScalars(const gsMatrix<index_t> & sizes, const std::vector<compOps> & ops)
{
    const index_t nc = sizes.cols(), np = sizes.rows();
    gsDofMapper m(sizes);
    for (index_t c = 0; c != nc; ++c)
        ops[c](m, c);
    m.finalize();

    std::vector<gsDofMapper> s;
    index_t F = 0, E = 0, cpld = 0;
    size_t map = 0;
    for (index_t c = 0; c != nc; ++c)
    {
        const gsVector<index_t> col = sizes.col(c);
        s.push_back(gsDofMapper(col));
        ops[c](s.back(), 0);
        s.back().finalize();
        F += s.back().freeSize();
        E += s.back().boundarySize();
        cpld += s.back().coupledSize();
        map += s.back().mapSize();
    }
    std::ostringstream expected, actual;
    expected << "free=" << F << " elim=" << E << " size=" << F + E << " cpld=" << cpld << " map=" << map;
    actual << "free=" << m.freeSize() << " elim=" << m.boundarySize() << " size=" << m.size()
           << " cpld=" << m.coupledSize() << " map=" << m.mapSize();
    std::vector<dofPair> pre;
    index_t Fb = 0, Eb = 0, n = 0;
    for (index_t c = 0; c != nc; ++c)
    {
        const gsDofMapper & sc = s[c];
        const gsVector<index_t> asVec = m.asVector(c);
        expected << "\nc" << c << " size=" << sc.size() << " free=" << sc.freeSize() << " total=" << sc.totalSize();
        actual   << "\nc" << c << " size=" << m.size(c) << " free=" << m.freeSize(c) << " total=" << m.totalSize(c);
        for (index_t k = 0, j = 0; k != np; ++k)
        {
            expected << "\n k" << k << " ps=" << sc.patchSize(k);
            actual   << "\n k" << k << " ps=" << m.patchSize(k, c);
            if (sc.patchSize(k) != m.patchSize(k, c))
                continue;
            for (index_t i = 0; i != static_cast<index_t>(sc.patchSize(k)); ++i, ++j, ++n)
            {
                const index_t g = sc.index(i, k), gl = m.index(i, k, c);
                const index_t e = sc.is_free_index(g) ? Fb + g : F + Eb + g - sc.freeSize();
                sc.preImage(g, pre);
                expected << " " << e << "c" << c << "/" << e << "/" << e << join(pre) << onPatches(sc, g);
                m.preImage(gl, pre);
                actual << " " << gl << "c" << m.componentOf(gl) << "/" << m.mapIndex(n) << "/" << asVec[j]
                       << join(pre) << onPatches(m, gl);
            }
        }
        Fb += sc.freeSize();
        Eb += sc.boundarySize();
    }
    CHECK_EQUAL(expected.str(), actual.str());
    return m;
}

// Raviart-Thomas type pair of degree 2 on nPatches unit squares side by side,
// nx x ny elements each: component c has degree 3 in direction c, 2 in the
// other. Component 0 is eliminated on the west side of patch 0, component 1
// on the south side of every patch.
struct rtFixture
{
    gsMatrix<index_t> sizes;
    std::vector<compOps> ops;

    rtFixture(unsigned nx, unsigned ny, index_t nPatches) : sizes(nPatches, 2)
    {
        const std::vector<boundaryInterface> ifs =
            gsNurbsCreator<real_t>::BSplineSquareGrid(nPatches, 1, 1.0).interfaces();
        for (index_t c = 0; c != 2; ++c)
        {
            const gsTensorBSplineBasis<2,real_t> bc(gsKnotVector<real_t>(0, 1, nx - 1, c == 0 ? 4 : 3),
                                                    gsKnotVector<real_t>(0, 1, ny - 1, c == 1 ? 4 : 3));
            sizes.col(c).setConstant(bc.size());
            ops.push_back([bc, ifs, c, nPatches](gsDofMapper & m, index_t comp) {
                gsMatrix<index_t> b1, b2;
                for (const boundaryInterface & bi : ifs)
                {
                    bc.matchWith(bi, bc, b1, b2);
                    m.matchDofs(bi.first().patch, b1, bi.second().patch, b2, comp);
                }
                if (0 == c)
                    m.markBoundary(0, bc.boundary(boundary::west), comp);
                else
                    for (index_t k = 0; k != nPatches; ++k)
                        m.markBoundary(k, bc.boundary(boundary::south), comp);
            });
        }
    }
};

} // anonymous namespace


SUITE(gsDofMapper_test)
{

// Equal columns give the mapper of gsDofMapper(column, nComp), also after
// coupling and elimination
TEST(matrix_equal_columns_as_vector)
{
    gsVector<index_t> sizes(3);
    sizes << 3, 0, 4;
    for (index_t nc = 1; nc <= 3; nc += 2)
    {
        gsDofMapper a(sizes, nc), b(gsMatrix<index_t>(sizes.replicate(1, nc)));
        for (gsDofMapper * m : {&a, &b})
        {
            m->matchDof(0, 2, 2, 0, -1);
            m->matchDof(0, 1, 2, 1, nc - 1);
            m->matchDof(2, 1, 2, 3, 0);
            m->eliminateDof(0, 0, -1);
            m->markBoundary(2, idx({0, 3}), nc - 1);
            m->markCoupled(2, 2, 0);
            m->finalize();
        }
        CHECK_EQUAL(state(a), state(b));
    }
}

// Each component numbers its dofs like a scalar mapper of its column.
// Component 1 has an empty patch, component 3 is left untouched.
TEST(matrix_components_as_scalar_mappers)
{
    gsMatrix<index_t> sizes(3, 4);
    sizes << 5, 4, 3, 2,
             7, 0, 3, 2,
             2, 6, 3, 0;
    // Array of mapper-operations to be applied to each component
    std::vector<compOps> ops;
    ops.push_back([](gsDofMapper & m, index_t c) {
        m.matchDofs(0, idx({1, 3}), 1, idx({0, 2}), c);
        m.eliminateDof(0, 0, c);
        m.matchDof(0, 0, 1, 1, c);
        m.markCoupled(1, 2, c);
        m.matchDof(1, 6, 2, 0, c); });
    ops.push_back([](gsDofMapper & m, index_t c) {
        m.matchDof(0, 3, 2, 0, c);
        m.eliminateDof(5, 2, c);
        m.matchDof(0, 2, 2, 5, c);
        m.markBoundary(2, idx({1, 2}), c); });
    ops.push_back([](gsDofMapper & m, index_t c) {
        for (index_t k = 0; k != 3; ++k)
            m.eliminateDof(1, k, c);
        m.matchDof(0, 2, 1, 0, c);
        m.matchDof(1, 2, 2, 0, c); });
    ops.push_back([](gsDofMapper &, index_t) { });
    const gsDofMapper m = checkAgainstScalars(sizes, ops);
    CHECK_EQUAL("off=0 5 12 ps=5 7 2;4 0 6;3 3 3;2 2 0", layout(m));
}

// Raviart-Thomas type pair on two patches of 4x2 elements: 28 and 30 dofs
// per patch. Compared component-wise to scalar mappers.
TEST(matrix_raviart_thomas_bases)
{
    const rtFixture rt(4, 2, 2);
    const gsDofMapper m = checkAgainstScalars(rt.sizes, rt.ops);
    CHECK(m.coupledSize() > 0);
}

// The per-patch sizes decide the layout, not the component totals (58 each)
TEST(matrix_patch_sizes_differ_with_equal_totals)
{
    gsMatrix<index_t> sizes(2, 2);
    sizes << 28, 31,
             30, 27;
    gsDofMapper m(sizes);
    m.finalize();
    CHECK_EQUAL("off=0 28 ps=28 30;31 27", layout(m));

    // Flat local dof 29 is on patch 1 for component 0, on patch 0 for component 1
    std::vector<dofPair> pre;
    m.preImage(29, pre);
    CHECK_EQUAL("(1,1)", join(pre));
    m.preImage(58 + 29, pre);
    CHECK_EQUAL("(0,29)", join(pre));
    CHECK(dofPair(1, 0) == m.anyPreImage(58 + 31));
    CHECK_EQUAL("29 -", onPatches(m, 58 + 29));
    CHECK_EQUAL("- 26", onPatches(m, 58 + 57));
}

// comp = -1 and setIdentity() throw on unequal components and leave the
// mapper unchanged; debug builds reject a local dof beyond its component
TEST(matrix_unequal_components_reject_shared_calls)
{
    gsMatrix<index_t> sizes(2, 2);
    sizes << 4, 2,
             4, 6;
    gsDofMapper m(sizes);
    CHECK_THROW(m.matchDof(0, 3, 1, 0, -1), std::runtime_error);
    CHECK_THROW(m.eliminateDof(0, 0, -1), std::runtime_error);
    CHECK_THROW(m.setIdentity(2, 6, 2), std::runtime_error);
#ifndef NDEBUG
    CHECK_THROW(m.eliminateDof(6, 1, 1), std::logic_error);
    CHECK_THROW(m.matchDof(1, 6, 0, 0, 1), std::logic_error);
#endif
    m.finalize();
    CHECK_EQUAL("off=0 4 ps=4 4;2 6 0,1,2,3|4,5,6,7;8,9|10,11,12,13,14,15", layout(m) + " " + indices(m));
}

TEST(matrix_invalid_sizes_throw)
{
    const index_t big = std::numeric_limits<index_t>::max() - 1;
    gsMatrix<index_t> noRows(0, 2), noCols(2, 0), negative(2, 2), compTotal(2, 1), total(1, 2);
    negative << 3, 4, -1, 2;
    compTotal << big, 1;
    total << big, 1;
    CHECK_THROW(gsDofMapper m(noRows), std::runtime_error);
    CHECK_THROW(gsDofMapper m(noCols), std::runtime_error);
    CHECK_THROW(gsDofMapper m(negative), std::runtime_error);
    CHECK_THROW(gsDofMapper m(compTotal), std::runtime_error);
    CHECK_THROW(gsDofMapper m(total), std::runtime_error);
}

// The matrix constructor is explicit; Eigen expressions still select the
// vector constructor (one component, one patch per entry)
static_assert(!std::is_convertible<gsMatrix<index_t>, gsDofMapper>::value,
              "gsDofMapper(const gsMatrix<index_t> &) must be explicit");
TEST(matrix_constructor_overloads)
{
    gsVector<index_t> v(2);
    v << 3, 4;
    gsMatrix<index_t> M(2, 2);
    M << 3, 5,
         4, 6;
    const gsDofMapper cst(gsVector<index_t>::Constant(2, 5)), plus(v.array() + 1), col(M.col(1)),
                      sum(v + v), row(v.transpose()), mat(M);
    for (const gsDofMapper * m : {&cst, &plus, &col, &sum, &row})
    {
        CHECK_EQUAL(1, m->numComponents());
        CHECK_EQUAL(2u, m->numPatches());
    }
    CHECK_EQUAL(2, mat.numComponents());
    CHECK_EQUAL(2u, mat.numPatches());
}

// setIdentity() on a moved-from mapper, which has no offsets left
TEST(identity_after_move)
{
    gsVector<index_t> sizes(2);
    sizes << 5, 7;
    gsDofMapper m(sizes);
    const gsDofMapper moved(std::move(m));
    m.setIdentity(2, 4, 1);
    m.finalize();
    CHECK_EQUAL("off=0 0 ps=0 4 |0,1,2,3", layout(m) + " " + indices(m));
}

// After setIdentity() the dofs are addressed through any patch, as callers
// with mapped bases do
TEST(identity_then_writes)
{
    gsDofMapper m;
    m.setIdentity(3, 6, 2);
    m.markBoundary(1, idx({4}), 0);
    m.matchDof(0, 1, 2, 3, 1);
    m.finalize();
    CHECK_EQUAL("off=0 0 0 ps=0 0 6;0 0 6 ||0,1,2,3,10,4;||5,9,6,9,7,8", layout(m) + " " + indices(m));
}

}
