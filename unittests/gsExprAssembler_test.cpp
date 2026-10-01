/** @file gsExprAssembler_test.cpp

    @brief Tests for gsExprAssembler

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): R. Schneckenleitner
*/

#include "gismo_unittest.h"

#include <gsAssembler/gsDofMapperCreator.h>
#include <gsMSplines/gsMappedBasis.h>

#include <atomic>
#include <functional>
#include <sstream>

namespace
{

/// Deliberately simple non-G+Smo rule used to verify that the assembler accepts
/// arbitrary gsQuadRule implementations, rather than merely another quRule
/// option.
template <class T>
class MidpointRule : public gismo::gsQuadRule<T>
{
public:
    explicit MidpointRule(index_t dim)
    {
        this->m_nodes.setZero(dim, 1);
        this->m_weights.resize(1);
        this->m_weights[0] = T(1);
        for (index_t i = 0; i < dim; ++i)
            this->m_weights[0] *= T(2);
    }
};

/// 2 unit squares side by side, one interface (patch0 east <-> patch1 west),
/// biquadratic and refined so that the interface carries several dofs.
/// Same construction as twoPatchBasis() in gsDofMapperCreator_test.cpp.
gsMultiBasis<real_t> twoPatchBasis(const gsMultiPatch<real_t> & mp)
{
    gsMultiBasis<real_t> mb(mp);
    mb.degreeElevate(1);   // bilinear -> biquadratic
    mb.uniformRefine(2);
    return mb;
}

/// The two components of a 2D Raviart-Thomas pair on the patches of
/// BSplineSquareGrid(2,1,1.0): component 0 = S^{3,2}_{2,1}, component 1 =
/// S^{2,3}_{1,2}, maximal regularity, e0 x e1 elements per patch.  Per patch
/// the components have (e0+3)(e1+2) and (e0+2)(e1+3) dofs: equal on an
/// isotropic mesh (4x4: 42 and 42), unequal on an anisotropic one (4x2: 28
/// and 30).  Same construction as rtPair() in gsDofMapperCreator_test.cpp.
std::vector<gsMultiBasis<real_t> > rtPair(index_t e0, index_t e1)
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2,1,1.0);
    std::vector<gsMultiBasis<real_t> > r;
    for (index_t c = 0; c != 2; ++c)
    {
        const short_t p0 = (0 == c ? 3 : 2), p1 = (0 == c ? 2 : 3);
        gsKnotVector<real_t> kv0(0.0, 1.0, e0-1, p0+1, 1, p0);
        gsKnotVector<real_t> kv1(0.0, 1.0, e1-1, p1+1, 1, p1);
        gsMultiBasis<real_t>::BasisContainer bases;
        bases.push_back(new gsTensorBSplineBasis<2,real_t>(kv0, kv1));
        bases.push_back(new gsTensorBSplineBasis<2,real_t>(kv0, kv1));
        r.push_back(gsMultiBasis<real_t>(bases, mp.topology()));
    }
    return r;
}

/// Collects what is written to std::cerr while in scope; GISMO_ENSURE
/// reports its message there and throws a generic std::runtime_error.
class CerrCapture
{
public:
    CerrCapture() : m_old(std::cerr.rdbuf(m_buf.rdbuf())) { }
    ~CerrCapture() { std::cerr.rdbuf(m_old); }
    bool contains(const std::string & s) const
    { return std::string::npos != m_buf.str().find(s); }
private:
    std::ostringstream m_buf;
    std::streambuf * m_old;
};

template <class F> bool throws(F f)
{
    try { f(); }
    catch (...) { return true; }
    return false;
}

} // anonymous namespace


SUITE(gsExprAssembler_test)
{
    TEST(CustomQuadratureFactory)
    {
        gsBSplineBasis<real_t> bb(0.0, 1.0, 3, 3);
        gsMultiBasis<real_t> mb(bb);
        gsBoundaryConditions<real_t> bcs;

        gsExprAssembler<real_t> standard(1, 1);
        standard.setIntegrationElements(mb);
        auto uStandard = standard.getSpace(mb);
        uStandard.setup(bcs, dirichlet::homogeneous, 0);
        standard.initSystem();
        standard.assemble(uStandard * uStandard.tr());
        const gsSparseMatrix<real_t> Mstandard = standard.matrix();

        gsExprAssembler<real_t> custom(1, 1);
        custom.setIntegrationElements(mb);
        auto uCustom = custom.getSpace(mb);
        uCustom.setup(bcs, dirichlet::homogeneous, 0);
        custom.initSystem();

        std::atomic<index_t> calls(0);
        std::atomic<bool> contextIsCorrect(true);
        custom.setQuadratureFactory(
            [&calls, &contextIsCorrect](const gsBasis<real_t> & basis,
                                        const gsOptionList &,
                                        index_t patch,
                                        short_t fixedDirection)
                -> gsExprAssembler<real_t>::QuadratureRulePtr
            {
                ++calls;
                if (patch != 0 || fixedDirection != -1)
                    contextIsCorrect = false;
                return gsExprAssembler<real_t>::QuadratureRulePtr(
                    new MidpointRule<real_t>(basis.dim()));
            });

        CHECK(custom.hasCustomQuadrature());
        custom.assemble(uCustom * uCustom.tr());
        const gsSparseMatrix<real_t> Mcustom = custom.matrix();
        
        const index_t current_calls = calls.load();
        CHECK(current_calls > 0);
        CHECK(contextIsCorrect.load());
        CHECK((Mstandard - Mcustom).norm() > 1e-8);

        custom.clearQuadratureFactory();
        CHECK(!custom.hasCustomQuadrature());
        custom.clearMatrix();
        custom.assemble(uCustom * uCustom.tr());
        const gsSparseMatrix<real_t> Mrestored = custom.matrix();

        CHECK_EQUAL(current_calls, calls.load());
        CHECK_EQUAL(Mstandard.rows(), Mrestored.rows());
        CHECK_EQUAL(Mstandard.cols(), Mrestored.cols());
        CHECK((Mstandard - Mrestored).norm() < 1e-14);
    }

    TEST(MultiSpaceBlockDims)
    {
        // _blockDims sizes the row/column blocks and resetDimensions sets the
        // shift applied to later blocks. mapper.freeSize() already reports the
        // component-inclusive dof count of a space, so multiplying it by that
        // space's dim counts the dimension twice. A single space of dim 1
        // cannot expose that -- the erroneous factor is 1 -- so this needs two
        // spaces of different, non-trivial dimension in one assembler.
        //
        // A 2D geometry is required: the assemble() arm below needs meas(G)
        // to build a genuine bilinear form per space, rather than only
        // inspecting mappers.
        gsMultiPatch<real_t> mp(*gsNurbsCreator<real_t>::BSplineSquare());
        gsMultiBasis<real_t> dbasis(mp);
        dbasis.degreeElevate(2);
        dbasis.uniformRefine(5);
        gsBoundaryConditions<real_t> bcs;
        const index_t n = dbasis.basis(0).size();

        // Control arm: two SCALAR spaces (dim 1 each). Since dim==1 makes
        // the erroneous factor in the fixed formula equal to 1, this arm
        // cannot itself fail on the bug -- it exists to pin that the shift
        // between blocks is otherwise sane, so a failure in the main arm
        // below can be attributed to the dim>1 handling rather than to the
        // fixture or to block bookkeeping in general.
        {
            gsExprAssembler<real_t> A(2, 2);
            A.setIntegrationElements(dbasis);
            auto u0 = A.getSpace(dbasis, 1, 0);
            auto u1 = A.getSpace(dbasis, 1, 1);
            u0.setup(bcs, dirichlet::homogeneous, 0);
            u1.setup(bcs, dirichlet::homogeneous, 0);
            A.initSystem();
            CHECK_EQUAL(2*n, A.numDofs());
            CHECK_EQUAL(n,   u1.mapper().firstIndex());
        }

        // Main arm: a vector-valued space of dimension d sharing an
        // assembler with a scalar space. Looping d = 2 and d = 3 is
        // essential -- the erroneous shift was dim()*(dim()*freeSize())
        // versus the correct dim()*freeSize(), i.e. an extra factor of
        // dim(). At a single d, a formula quadratic in dim() is
        // indistinguishable from a linear one with a different constant;
        // the pair of values pins the exponent.
        for (index_t d = 2; d <= 3; ++d)
        {
            gsExprAssembler<real_t> A(2, 2);
            A.setIntegrationElements(dbasis);
            auto G = A.getMap(mp);
            auto v = A.getSpace(dbasis, d, 0); // vector-valued space, dim d
            auto p = A.getSpace(dbasis, 1, 1); // scalar space, dim 1
            v.setup(bcs, dirichlet::homogeneous, 0);
            p.setup(bcs, dirichlet::homogeneous, 0);
            A.initSystem();

            CHECK_EQUAL((d+1)*n, A.numDofs());
            CHECK_EQUAL(0,       v.mapper().firstIndex());
            CHECK_EQUAL(d*n,     v.mapper().freeSize());
            CHECK_EQUAL(d*n,     p.mapper().firstIndex());
            CHECK_EQUAL(n,       p.mapper().freeSize());

            // Two distinct-block terms in one assemble() call: this is what
            // actually writes into both diagonal blocks of the system
            // matrix, so a wrong offset for p's block (computed pre-fix as
            // d*(d*n) instead of d*n) would either write out of range or
            // leave part of v's block untouched.
            A.assemble(v*v.tr()*meas(G), p*p.tr()*meas(G));

            // Counting entirely-zero rows separates "wrong total" from
            // "structurally broken": pre-fix, each block was internally
            // self-consistent and merely sat at the wrong offset, so
            // numDofs() alone does not reveal that (d^2-d)*n rows exist
            // that nothing ever writes to.
            const gsSparseMatrix<real_t> & M = A.matrix();
            gsVector<bool> rowTouched(M.rows());
            rowTouched.setZero();
            for (index_t c = 0; c != M.cols(); ++c)
                for (gsSparseMatrix<real_t>::InnerIterator it(M, c); it; ++it)
                    rowTouched(it.row()) = true;
            const index_t zeroRows = M.rows() - rowTouched.array().count();
            CHECK_EQUAL(0, zeroRows);

            // matrixBlockView() is the only assertion here reaching
            // _blockDims directly (numDofs()/firstIndex() above reach the
            // same defect only through resetDimensions): the block sizes
            // must match the per-space free dof counts exactly.
            auto view = A.matrixBlockView();
            CHECK_EQUAL(d*n, view(0,0).rows());
            CHECK_EQUAL(d*n, view(0,0).cols());
            CHECK_EQUAL(n,   view(1,1).rows());
            CHECK_EQUAL(n,   view(1,1).cols());
        }
    }

    // matrix() is `m_modified ? makeMatrix() : m_matrix`, and m_matrix is only
    // populated from the fiber matrix by makeMatrix(). clearMatrix() must
    // therefore invalidate the cache on every path, including the one that
    // resizes rather than zeroes: initSystem() reaches exactly that path, so
    // without the flag matrix() reports the default-constructed 0x0 m_matrix
    // for a system whose dimensions are already known.
    TEST(MatrixSizedAfterInitSystem)
    {
        gsBSplineBasis<real_t> bb(0.0, 1.0, 3, 3);
        gsMultiBasis<real_t> mb(bb);
        gsBoundaryConditions<real_t> bcs;

        gsExprAssembler<real_t> A(1, 1);
        A.setIntegrationElements(mb);
        auto u = A.getSpace(mb, 1, 0);
        u.setup(bcs, dirichlet::homogeneous, 0);
        A.initSystem();

        // No assemble() here: sizing must hold from initSystem() alone.
        CHECK_EQUAL(A.numTestDofs(), A.matrix().rows());
        CHECK_EQUAL(A.numDofs(),     A.matrix().cols());
        CHECK(A.matrix().rows() > 0);
    }

    TEST(InterfaceExpression)
    {
        const index_t numRef = 2;
        gsVector<> translation(2);
        translation << -0.5, -1;
        gsFunctionExpr<> ff("if(x>0,1,-1)", 2); // function with jump
        gsMultiPatch<> patches = gsNurbsCreator<>::BSplineSquareGrid(1,2,1);
        patches.patch(0).translate(translation);
        patches.patch(1).translate(translation);

        gsMultiBasis<> mb(patches);

        for(index_t i = 0; i < numRef; i++)
            mb.uniformRefine();

        gsExprEvaluator<> ev;
        ev.setIntegrationDomain(mb.domain());
        gsExprEvaluator<>::geometryMap G = ev.getMap(patches);
        auto f = ev.getVariable(ff, G);

        ev.integral(f);
        const real_t v = ev.value();
        CHECK( v*v < 1e-10 );

        CHECK(1==patches.interfaces().size());
        ev.integralInterface(f.left() + f.right() , patches.interfaces());
        const real_t w = ev.value();
        CHECK( w*w < 1e-10 );
    }

    TEST(BoundaryIntegral)
    {
        // Create a circle
        gsMultiPatch<> mp;
        mp.addPatch(gsNurbsCreator<>::NurbsDisk(1.0));
        mp.computeTopology();
        mp.embed(3);
        mp.uniformRefine(1);
        mp.uniformRefine(1);

        // Rotate it 30 degrees
        gsVector<real_t,3> rv;
        rv.setZero(); rv[0]=1;
        real_t angle = EIGEN_PI/3;
        mp.patch(0).rotate(angle, rv);

        // Create evaluator
        gsMultiBasis<> mb(mp);
        gsExprEvaluator<> ev;
        ev.setIntegrationDomain(mb.domain());
        ev.options().addReal("quA","",2); // added precision to approx pi
        ev.options().addInt("quB","", 2); // added precision to approx pi
        auto G = ev.getMap(mp);
        typedef gsExprAssembler<>::element element;
        element el = ev.getElement();

        // Test measure of a curve in 3D
        CHECK_CLOSE(ev.integralBdr(meas(G)), 2*EIGEN_PI, EPSILON);
        // Test tangent of a curve in 3D
        CHECK(math::abs(ev.integralBdr(tv(G).norm())-2*EIGEN_PI) < 1e-10);

        // Test measure of a surface in 3D
        CHECK(math::abs(ev.integral(meas(G))-EIGEN_PI) < 1e-10);
        // Test surface normal of a surface in 3D
        CHECK(math::abs(ev.integral(sn(G).norm())-EIGEN_PI) < 1e-10);
        // XXXX
        CHECK(math::abs(ev.integral(el.area(G))-2*EIGEN_PI/32) < 1e-10);

        mp.patch(0).rotate(-angle, rv);
        mp.embed(2);

        // Test measure of a curve in 3D
        CHECK(math::abs(ev.integralBdr(meas(G))-2*EIGEN_PI) < 1e-10);
        // Test tangent of a curve in 3D
        CHECK(math::abs(ev.integralBdr(tv(G).norm())-2*EIGEN_PI) < 1e-10);
        // Test outer normal of a curve in 3D
        CHECK(math::abs(ev.integralBdr(nv(G).norm())-2*EIGEN_PI) < 1e-10);
        //
        CHECK(math::abs(ev.integral(el.area(G))-2*EIGEN_PI/32) < 1e-10);
    }

    // Pins: a custom, FINALIZED, patch-concatenated mapper handed to
    // gsFeSpace::setupMapper survives initSystem() and drives the assembly.
    //
    // The mapper used here is *conforming* while the one the assembler builds
    // for itself (gsFeSpaceData::init) is non-conforming, so the two differ in
    // size: numDofs()==custom.freeSize() < mb.totalSize() can only hold if the
    // caller's mapper was really retained.  coupledSize()>0 pins that its
    // interface identification survived as well.  An accepted uniform mapper
    // must be neither rejected nor silently replaced.
    TEST(CustomUniformMapperRetained)
    {
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2,1,1.0);
        gsMultiBasis<real_t> mb = twoPatchBasis(mp);

        gsDofMapper custom = createMapper(mb, 1, /*conforming=*/true, /*finalize=*/true);
        const index_t customFree    = custom.freeSize();
        const index_t customCoupled = custom.coupledSize();
        const index_t customBdr     = custom.boundarySize();

        // The custom mapper is genuinely different from the assembler default
        CHECK(customCoupled > 0);
        CHECK(customFree < static_cast<index_t>(mb.totalSize()));

        gsExprAssembler<real_t> A(1,1);
        A.setIntegrationElements(mb);
        auto G = A.getMap(mp);
        auto u = A.getSpace(mb, 1);

        u.setupMapper(give(custom)); // takes the mapper by value and moves from it
        const_cast<expr::gsFeSpace<real_t>&>(u).fixedPart()
            .setZero(u.mapper().boundarySize(), 1);

        A.initSystem();

        // Retained, not rebuilt
        CHECK_EQUAL(customFree   , A.numDofs()              );
        CHECK_EQUAL(customFree   , u.mapper().freeSize()    );
        CHECK_EQUAL(customCoupled, u.mapper().coupledSize() );
        CHECK_EQUAL(customBdr    , u.mapper().boundarySize());
        CHECK_EQUAL(1            , u.mapper().numComponents());

        A.assemble(u * u.tr() * meas(G));
        const gsSparseMatrix<real_t> M = A.matrix();

        CHECK_EQUAL(customFree, M.rows());
        CHECK_EQUAL(customFree, M.cols());
        // partition of unity => sum of all mass-matrix entries is the area
        CHECK_CLOSE(2.0, M.sum(), 1e-10);
    }

    // Pins: a mapper produced from a gsMappedBasis -- i.e. one built through
    // gsDofMapper::setIdentity, whose local indices are already global basis
    // indices rather than patch-concatenated ones -- is accepted by
    // setupMapper, survives initSystem() and assembles correctly.
    //
    // Matters because the global-identity layout is a second, structurally
    // different mapper layout that must keep working.
    TEST(MappedBasisIdentityMapperRetained)
    {
        gsNurbsCreator<real_t>::TensorBSpline2Ptr geom =
            gsNurbsCreator<real_t>::BSplineSquare(2); // [0,2]x[0,2], area 4
        gsMultiPatch<real_t> mp(*geom);

        gsMultiBasis<real_t> mb(mp);
        mb.basis(0).uniformRefine(2);
        const index_t sz = mb.basis(0).size();

        // trivial (identity) mapped basis
        gsSparseMatrix<real_t> ident(sz, sz);
        ident.setIdentity();
        gsMappedBasis<2,real_t> mapB(mb, ident);
        CHECK_EQUAL(sz, static_cast<index_t>(mapB.size()));

        // goes through gsDofMapper::setIdentity
        gsDofMapper custom = createMapper(mapB, 1, /*conforming=*/false, /*finalize=*/true);
        CHECK_EQUAL(sz, custom.freeSize());
        CHECK_EQUAL(0 , custom.boundarySize());

        gsExprAssembler<real_t> A(1,1);
        A.setIntegrationElements(mb);
        auto G = A.getMap(mp);
        auto u = A.getSpace(mapB, 1);

        u.setupMapper(give(custom));
        const_cast<expr::gsFeSpace<real_t>&>(u).fixedPart()
            .setZero(u.mapper().boundarySize(), 1);

        A.initSystem();

        CHECK_EQUAL(sz, A.numDofs()          );
        CHECK_EQUAL(sz, u.mapper().freeSize());

        A.assemble(u * u.tr() * meas(G));
        const gsSparseMatrix<real_t> M = A.matrix();

        CHECK_EQUAL(sz, M.rows());
        CHECK_EQUAL(sz, M.cols());
        CHECK_CLOSE(4.0, M.sum(), 1e-10);
    }

    // TODO: verify if this is really desirable behavior.
    //
    // Pins accepted behaviour: installing a mapper whose numComponents()
    // differs from the space dimension is *accepted* rather than an error,
    // and the mapper is then silently replaced.
    //
    // gsFeSpace::setupMapper only asserts
    //     mapSize() == source().size()*dofsMapper.numComponents()
    // which a 3-component mapper over the same basis satisfies, so the install
    // does not throw.  gsFeSpaceData::valid() however is
    //     fs->size()*dim == mapper.mapSize()
    // which is false here, so resetDimensions() calls init() and rebuilds a
    // default 2-component NON-conforming mapper, discarding the caller's.
    //
    // This is exactly the gsBarrierPatch/gsBarrierCore caller convention
    // (createMapper(mb, targetDim) installed into getSpace(mb, d)).  Both
    // halves are pinned: no throw, and the silent rebuild.
    TEST(MapperComponentCountMismatchAccepted)
    {
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2,1,1.0);
        gsMultiBasis<real_t> mb = twoPatchBasis(mp);

        gsExprAssembler<real_t> A(1,1);
        A.setIntegrationElements(mb);
        auto G = A.getMap(mp);
        auto u = A.getSpace(mb, /*dim=*/2);

        // 3 components against a 2-dimensional space
        gsDofMapper custom = createMapper(mb, /*nComp=*/3, /*conforming=*/true,
                                          /*finalize=*/true);
        CHECK_EQUAL(3, custom.numComponents());
        CHECK(custom.coupledSize() > 0);

        bool threw = false;
        try { u.setupMapper(give(custom)); }
        catch (...) { threw = true; }
        CHECK(!threw);                                  // accepted today
        CHECK_EQUAL(3, u.mapper().numComponents());     // and really installed

        A.initSystem();

        // ... but silently discarded and rebuilt as the 2-component default
        CHECK_EQUAL(2, u.mapper().numComponents());
        CHECK_EQUAL(0, u.mapper().coupledSize());       // default is non-conforming
        CHECK_EQUAL(2*static_cast<index_t>(mb.totalSize()), A.numDofs());
        CHECK_EQUAL(2*static_cast<index_t>(mb.totalSize()), u.mapper().freeSize());

        A.assemble(u * u.tr() * meas(G));
        const gsSparseMatrix<real_t> M = A.matrix();
        CHECK_EQUAL(2*static_cast<index_t>(mb.totalSize()), M.rows());
        CHECK_EQUAL(2*static_cast<index_t>(mb.totalSize()), M.cols());
        // one scalar mass matrix per component
        CHECK_CLOSE(2.0*2.0, M.sum(), 1e-10);
    }

    // TODO: check last remark here
    //
    // Pins refine-and-reassemble: a basis refined IN PLACE after a custom
    // mapper was installed invalidates gsFeSpaceData::valid(), so the space's
    // mapper is rebuilt from scratch and fixedDofs is cleared.  The second
    // assembly must be identical to a from-scratch assembly on the refined
    // basis with no custom mapper at all.
    //
    // The custom mapper's conforming (coupling) choices are LOST by this --
    // that is existing behaviour and is pinned here deliberately.
    TEST(RefineAfterMapperInstallRebuilds)
    {
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2,1,1.0);
        gsMultiBasis<real_t> mb = twoPatchBasis(mp);

        // conforming AND with an eliminated Dirichlet boundary, so that both
        // kinds of choice a custom mapper can carry are present and can be
        // observed to disappear
        gsFunctionExpr<real_t> g("0", 2);
        gsBoundaryConditions<real_t> bc;
        bc.addCondition(0, boundary::west, condition_type::dirichlet, &g);

        gsDofMapper custom = createMapper(mb, bc, 1, 0, /*conforming=*/true,
                                          /*finalize=*/true);
        const index_t customFree    = custom.freeSize();
        const index_t customCoupled = custom.coupledSize();
        const index_t customBdr     = custom.boundarySize();
        CHECK(customCoupled > 0);
        CHECK(customBdr     > 0);

        gsExprAssembler<real_t> A(1,1);
        A.setIntegrationElements(mb);
        auto G = A.getMap(mp);
        auto u = A.getSpace(mb, 1);

        u.setupMapper(give(custom));
        const_cast<expr::gsFeSpace<real_t>&>(u).fixedPart()
            .setZero(u.mapper().boundarySize(), 1);
        A.initSystem();
        CHECK_EQUAL(customFree, A.numDofs());
        CHECK_EQUAL(customBdr , u.fixedPart().size());

        // refine the very same gsMultiBasis object the space points at
        mb.uniformRefine();
        A.setIntegrationElements(mb); // the integration domain must see it too
        A.initSystem();

        // the custom mapper is gone: rebuilt, non-conforming, no elimination
        CHECK_EQUAL(static_cast<index_t>(mb.totalSize()), A.numDofs()             );
        CHECK_EQUAL(static_cast<index_t>(mb.totalSize()), u.mapper().freeSize()   );
        CHECK_EQUAL(0, u.mapper().coupledSize()  );
        CHECK_EQUAL(0, u.mapper().boundarySize() );
        CHECK_EQUAL(0, u.fixedPart().size()      );

        A.assemble(u * u.tr() * meas(G));
        const gsSparseMatrix<real_t> M = A.matrix();

        // from-scratch reference on the refined basis, no custom mapper
        gsExprAssembler<real_t> B(1,1);
        B.setIntegrationElements(mb);
        auto G2 = B.getMap(mp);
        auto v = B.getSpace(mb, 1);
        B.initSystem();
        B.assemble(v * v.tr() * meas(G2));
        const gsSparseMatrix<real_t> Mref = B.matrix();

        CHECK_EQUAL(Mref.rows()    , M.rows()    );
        CHECK_EQUAL(Mref.cols()    , M.cols()    );
        CHECK_EQUAL(Mref.nonZeros(), M.nonZeros());
        CHECK((M - Mref).norm() < 1e-14);
        CHECK_CLOSE(2.0, M.sum(), 1e-10);
    }
    // A Raviart-Thomas mapper -- one basis per component, declared distinct
    // by the per-component creator -- is rejected by setupMapper in every
    // build type, on the anisotropic mesh (unequal component sizes) and on
    // the isotropic one.  On the isotropic mesh every size the expression
    // layer could compare agrees with a uniform 2-component mapper over
    // component 0's basis, so only the declared flag can reject it.
    TEST(DistinctComponentMapperRejectedBySetupMapper)
    {
        const index_t meshes[2][2] = { {4,2}, {4,4} };
        for (index_t mesh = 0; mesh != 2; ++mesh)
        {
            const std::vector<gsMultiBasis<real_t> > rt =
                rtPair(meshes[mesh][0], meshes[mesh][1]);
            const gsDofMapper rtMapper =
                createMapper(rt, gsBoundaryConditions<real_t>(), 0, true, true);
            CHECK(rtMapper.hasDistinctComponentSpaces());

            if (meshes[mesh][0] == meshes[mesh][1])
            {
                CHECK_EQUAL(rtMapper.patchSize(0,0), rtMapper.patchSize(0,1));
                CHECK_EQUAL(rtMapper.patchSize(1,0), rtMapper.patchSize(1,1));
                CHECK_EQUAL(2*rt[0].totalSize(), rtMapper.mapSize());
            }

            gsExprAssembler<real_t> A(1,1);
            A.setIntegrationElements(rt[0]);
            auto u = A.getSpace(rt[0], 2);
            const size_t before = u.mapper().mapSize();

            CerrCapture err;
            CHECK_THROW(u.setupMapper(rtMapper), std::runtime_error);
            CHECK(err.contains("distinct per-component bases"));
            CHECK_EQUAL(before, u.mapper().mapSize()); // not installed
            CHECK(!u.mapper().hasDistinctComponentSpaces());
        }
    }

    // The same mappers installed through the mutable gsFeSpace::mapper()
    // reference, which bypasses setupMapper: initSystem() rejects them
    // before its valid()/init() rebuild instead of silently replacing
    // (anisotropic) or assembling (isotropic) them.  gsFeSolution::check()
    // rejects them as well.
    TEST(DistinctComponentMapperRejectedThroughMutableMapper)
    {
        const index_t meshes[2][2] = { {4,2}, {4,4} };
        for (index_t mesh = 0; mesh != 2; ++mesh)
        {
            const std::vector<gsMultiBasis<real_t> > rt =
                rtPair(meshes[mesh][0], meshes[mesh][1]);
            const gsDofMapper rtMapper =
                createMapper(rt, gsBoundaryConditions<real_t>(), 0, true, true);

            gsExprAssembler<real_t> A(1,1);
            A.setIntegrationElements(rt[0]);
            auto u = A.getSpace(rt[0], 2);
            u.mapper() = rtMapper;

            {
                CerrCapture err;
                CHECK_THROW(A.initSystem(), std::runtime_error);
                CHECK(err.contains("distinct per-component bases"));
            }
            CHECK(u.mapper().hasDistinctComponentSpaces()); // not replaced

            gsMatrix<real_t> sol(rtMapper.freeSize(), 1);
            sol.setZero();
            const expr::gsFeSolution<real_t> s(u, sol);
            {
                CerrCapture err;
                CHECK_THROW(s.check(), std::runtime_error);
                CHECK(err.contains("distinct per-component bases"));
            }

            // the same through a separate test space: its mapper is the
            // only one that is unusable
            gsExprAssembler<real_t> B(1,1);
            B.setIntegrationElements(rt[0]);
            auto w = B.getSpace(rt[0], 2);
            auto v = B.getTestSpace(w, rt[0]);
            v.mapper() = rtMapper;
            {
                CerrCapture err;
                CHECK_THROW(B.initSystem(), std::runtime_error);
                CHECK(err.contains("distinct per-component bases"));
            }
            CHECK(v.mapper().hasDistinctComponentSpaces()); // not replaced
        }
    }

    // A mapper replaced through gsFeSpace::mapper() AFTER initSystem() is
    // rejected by every assembly and pattern call, in every build type.  On
    // the isotropic mesh the system dimensions still agree, so nothing else
    // would notice that component 1 is indexed with component 0's basis.
    TEST(DistinctComponentMapperRejectedAfterInitialization)
    {
        const index_t meshes[2][2] = { {4,2}, {4,4} };
        for (index_t mesh = 0; mesh != 2; ++mesh)
        {
            const std::vector<gsMultiBasis<real_t> > rt =
                rtPair(meshes[mesh][0], meshes[mesh][1]);
            const gsDofMapper rtMapper =
                createMapper(rt, gsBoundaryConditions<real_t>(), 0, true, true);

            gsExprAssembler<real_t> A(1,1);
            A.setIntegrationElements(rt[0]);
            auto u = A.getSpace(rt[0], 2);
            A.initSystem();
            u.mapper() = rtMapper;

            // Each entry point checks before it looks at its arguments, so
            // empty boundary and interface containers suffice.
            typedef gsExprAssembler<real_t>::bcRefList bcRefList;
            typedef gsExprAssembler<real_t>::bContainer bContainer;
            typedef gsExprAssembler<real_t>::ifContainer ifContainer;
            auto uu = u * u.tr();
            gsMatrix<real_t> sol(rtMapper.freeSize(), 1);
            sol.setZero();
            expr::gsFeSolution<real_t> s(u, sol);

            const std::vector<std::pair<const char *, std::function<void()> > > entries = {
                { "computePattern",      [&]{ A.computePattern(uu); } },
                { "computePatternBdr",   [&]{ A.computePatternBdr(bcRefList(), uu); } },
                { "computePatternIfc",   [&]{ A.computePatternIfc(ifContainer(), uu); } },
                { "assemble",            [&]{ A.assemble(uu); } },
                { "assembleBdr(bc)",     [&]{ A.assembleBdr(bcRefList(), uu); } },
                { "assembleBdr(sides)",  [&]{ A.assembleBdr(bContainer(), uu); } },
                { "assembleIfc",         [&]{ A.assembleIfc(ifContainer(), uu); } },
                { "assembleJacobian",    [&]{ A.assembleJacobian(u, s); } },
                // assembleJacobianIfc is checked the same way but cannot be
                // instantiated here: it calls gsDomain::beginAll(boxSide),
                // which does not exist.
            };
            for (size_t e = 0; e != entries.size(); ++e)
            {
                CerrCapture err;
                const bool threw = throws(entries[e].second);
                CHECK(threw);
                CHECK(err.contains("distinct per-component bases"));
                if (!threw || !err.contains("distinct per-component bases"))
                    gsInfo << "not rejected by " << entries[e].first << "\n";
            }
        }
    }

    // A mapper that was not declared distinct but whose components differ in
    // size is rejected with the offending component and patch named --
    // including equal component totals split differently over the patches,
    // which the mapSize() consistency assert in setupMapper cannot see.
    TEST(UnequalComponentSizesRejected)
    {
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2,1,1.0);
        gsMultiBasis<real_t> mb = twoPatchBasis(mp);

        gsVector<index_t> sizes0(2), sizes1(2);
        sizes0 << 4, 6;
        sizes1 << 6, 4;
        std::vector<gsVector<index_t> > sizes;
        sizes.push_back(sizes0);
        sizes.push_back(sizes1);
        gsDofMapper split(sizes, /*hasDistinctComponentSpaces=*/false);
        split.finalize();
        CHECK_EQUAL(split.totalSize(0), split.totalSize(1));

        gsDofMapper identity;
        std::vector<size_t> totals(2);
        totals[0] = 5;
        totals[1] = 7;
        identity.setIdentity(2, totals);
        identity.finalize();

        gsExprAssembler<real_t> A(1,1);
        A.setIntegrationElements(mb);
        auto u = A.getSpace(mb, 2);
        {
            CerrCapture err;
            CHECK_THROW(u.setupMapper(split), std::runtime_error);
            CHECK(err.contains("component 1 has 6 dofs on patch 0, component 0 has 4"));
        }
        {
            CerrCapture err;
            CHECK_THROW(u.setupMapper(identity), std::runtime_error);
            CHECK(err.contains("component 1 has 7 dofs, component 0 has 5"));
        }
        {
            u.mapper() = split;
            CerrCapture err;
            CHECK_THROW(A.initSystem(), std::runtime_error);
            CHECK(err.contains("component 1 has 6 dofs on patch 0"));
        }
    }

    // The rejection leaves every mapper the single-basis creators produce
    // alone, including a default-constructed one (the state before the
    // first initSystem()) and a uniform one whose component count differs
    // from the space dimension.
    TEST(UniformMappersNotRejected)
    {
        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2,1,1.0);
        gsMultiBasis<real_t> mb = twoPatchBasis(mp);

        typedef expr::gsFeSpaceData<real_t> Data;
        CHECK(!throws([&]{ Data::ensureUsableByUniformEvaluator(gsDofMapper()); }));
        for (index_t nComp = 1; nComp != 4; ++nComp)
            CHECK(!throws([&]{ Data::ensureUsableByUniformEvaluator(
                                   createMapper(mb, nComp, true, true)); }));

        gsExprAssembler<real_t> A(1,1);
        A.setIntegrationElements(mb);
        auto u = A.getSpace(mb, 2);
        CHECK(!throws([&]{ A.initSystem(); }));
        CHECK(!throws([&]{ u.setupMapper(createMapper(mb, 2, true, true)); }));
        CHECK(!throws([&]{ A.initSystem(); }));
    }
}
