/** @file gsAdaptiveMultiPatchBuilder_test.cpp
    @brief Regression tests for hierarchical 3D moving mesh updates.
    This file is part of the G+Smo library, under the Mozilla Public License v2.0.

    Author(s): M. BAHARI
    */
#include "gismo_unittest.h"
#include <gsAssembler/gsAdaptiveMultiPatchBuilder.h>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace
{
void checkParameterConversion(const gsMultiPatch<>& original,
                              const gsAdaptiveMultiPatchBuilder& builder)
{
    const index_t dim = original.dim();
    gsMatrix<> box(dim,2);
    box.col(0).setZero(); box.col(1).setOnes();
    gsVector<index_t> counts(dim); counts.setConstant(5);
    const gsMatrix<> points = gsGridIterator<real_t,CUBE>(box,counts).toMatrix();
    for (size_t p = 0; p < original.nPatches(); ++p)
    {
        const gsMatrix<> bounds = original.basis(p).support();
        gsMatrix<> native = points;
        for (index_t d = 0; d < dim; ++d)
            native.row(d) = (bounds(d,0)+(bounds(d,1)-bounds(d,0))*points.row(d).array())
                .max(bounds(d,0)).min(bounds(d,1));
        CHECK((builder.parameterConversion().patch(p).eval(points)-native).norm() < 1e-10);
        const auto& normalized = builder.normalizedGeometry().patch(p);
        CHECK((normalized.eval(points)-original.patch(p).eval(native)).norm() < 1e-9);
        CHECK((normalized.coefs()-original.patch(p).coefs()).norm() == 0.);
        CHECK((normalized.basis().weights()-original.basis(p).weights()).norm() == 0.);
        CHECK(normalized.basis().support().col(0).norm() == 0.);
        CHECK((normalized.basis().support().col(1).array()-1.).matrix().norm() == 0.);
        const gsMatrix<> derivative = normalized.jacobian(points);
        gsMatrix<> expected = original.patch(p).jacobian(native);
        for (index_t q = 0; q < points.cols(); ++q)
            for (index_t d = 0; d < dim; ++d)
                expected.col(dim*q+d) *= bounds(d,1)-bounds(d,0);
        CHECK((derivative-expected).norm() < 1e-8);
    }
}

template <short_t d>
void checkDifferentParameterBoxes()
{
    const gsMultiPatch<> unit = d == 2
        ? gsNurbsCreator<>::BSplineSquareGrid(2,1,1.)
        : gsNurbsCreator<>::BSplineCubeGrid(2,1,1,1.,-.5,-.5,-.5);
    gsMultiPatch<> original = unit;
    for (size_t p = 0; p < original.nPatches(); ++p)
        for (index_t dir = 0; dir < d; ++dir)
            dynamic_cast<gsTensorBSplineBasis<d>&>(original.patch(p).basis()).knots(dir)
                .affineTransformTo(-3.-p-dir, 2.+3.*p+2.*dir);
    gsAdaptiveMultiPatchBuilder expected(unit,1,2,4.,0,0,true);
    gsAdaptiveMultiPatchBuilder actual(original,1,2,4.,0,0,true);
    checkParameterConversion(original,actual);
    gsFunctionExpr<> f(d == 2 ? "1+x+2*y" : "1+x+2*y+3*z", d);
    const auto rho = actual.buildAnalyticDensity(f);
    const auto expectedRho = expected.buildAnalyticDensity(f);
    for (size_t p = 0; p < original.nPatches(); ++p)
        CHECK((rho.patch(p).coefs()-expectedRho.patch(p).coefs()).norm() < 1e-10);
    actual.buildMultiPatch(rho);
    expected.buildMultiPatch(expectedRho);
    const auto composition = actual.buildColCompMultiPatch();
    const auto expectedComposition = expected.buildColCompMultiPatch();
    for (size_t p = 0; p < original.nPatches(); ++p)
    {
        CHECK((actual.MAmapping.patch(p).coefs()-expected.MAmapping.patch(p).coefs()).norm() < 1e-10);
        CHECK((composition.patch(p).coefs()-expectedComposition.patch(p).coefs()).norm() < 1e-10);
    }
    gsMultiBasis<> basis(composition,true);
    std::vector<bool> marked(basis.totalElements(),false);
    for (size_t i = 0; i < marked.size(); i += 3) marked[i] = true;
    const auto markedRho = actual.buildDensity(basis,marked);
    const auto expectedMarkedRho = expected.buildDensity(basis,marked);
    for (size_t p = 0; p < original.nPatches(); ++p)
        CHECK((markedRho.patch(p).coefs()-expectedMarkedRho.patch(p).coefs()).norm() < 1e-10);
}

void checkIdentity(const gsMultiPatch<>& mapping)
{
    for (size_t p = 0; p < mapping.nPatches(); ++p)
    {
        const gsMatrix<> points = mapping.basis(p).anchors();
        CHECK((mapping.patch(p).eval(points)-points).norm() < 1e-10);
    }
}

void checkExterior(const gsMultiPatch<>& reference, const gsMultiPatch<>& mapping)
{
    for (auto it = reference.bBegin(); it != reference.bEnd(); ++it)
    {
        const auto& B = mapping.basis(it->patch);
        const gsMatrix<> points = B.anchors();
        const gsMatrix<> identity = B.interpolateData(points, points)->coefs();
        const gsMatrix<index_t> boundary = B.boundary(it->side());
        for (index_t i = 0; i < boundary.size(); ++i)
            CHECK((mapping.patch(it->patch).coefs().row(boundary(i))-identity.row(boundary(i))).norm() < 1e-11);
    }
}

void checkInterfaces(const gsMultiPatch<>& reference, const gsMultiPatch<>& mapping)
{
    for (const auto& face : reference.interfaces())
    {
        const auto first = face.first(), second = face.second();
        const index_t dim = mapping.dim();
        gsVector<index_t> counts(dim);
        counts.setConstant(5);
        counts(first.direction()) = 1;
        gsMatrix<> box(dim,2);
        box.col(0).setZero();
        box.col(1).setOnes();
        box.row(first.direction()).setConstant(first.parameter());
        const gsMatrix<> points = gsGridIterator<real_t,CUBE>(box, counts).toMatrix();
        gsMatrix<> target(dim, points.cols());
        target.row(second.direction()).setConstant(second.parameter());
        for (index_t d = 0; d < dim; ++d)
            if (d != first.direction())
            {
                const index_t mapped = face.dirMap(first,d);
                target.row(mapped) = points.row(d);
                if (!face.dirOrientation(first,d)) target.row(mapped).array() = 1.-target.row(mapped).array();
            }
        const gsMatrix<> left = reference.patch(first.patch).eval(mapping.patch(first.patch).eval(points));
        const gsMatrix<> right = reference.patch(second.patch).eval(mapping.patch(second.patch).eval(target));
        CHECK((left-right).norm() < 1e-10);
    }
}

void checkPositiveJacobian(const gsMultiPatch<>& mapping)
{
    const index_t dim = mapping.dim();
    gsVector<index_t> counts(dim);
    counts.setConstant(7);
    for (size_t p = 0; p < mapping.nPatches(); ++p)
    {
        const gsMatrix<> points = gsGridIterator<real_t,CUBE>(mapping.basis(p).support(), counts).toMatrix();
        const gsMatrix<> J = mapping.patch(p).jacobian(points);
        for (index_t q = 0; q < points.cols(); ++q)
            CHECK(J.middleCols(dim*q,dim).determinant() > 0.);
    }
}
struct InteriorThreadSettings
{
#ifdef _OPENMP
    int threads = omp_get_max_threads();
    int dynamic = omp_get_dynamic();
    InteriorThreadSettings() { omp_set_dynamic(0); omp_set_num_threads(1); }
    ~InteriorThreadSettings() { omp_set_num_threads(threads); omp_set_dynamic(dynamic); }
#endif
};

void checkParallelInteriors(index_t dim)
{
    InteriorThreadSettings settings;
#ifdef _OPENMP
    const std::vector<int> counts = {1,2,4};
#else
    const std::vector<int> counts = {1};
#endif
    for (bool fixed : {false,true})
    {
        gsMultiPatch<> reference;
        if (dim == 2)
            reference = gsNurbsCreator<>::BSplineSquareGrid(2,2,.5);
        else
            reference = gsNurbsCreator<>::BSplineCubeGrid(2,2,1,1.,-.5,-.5,-.5);
        reference.degreeElevate(1);
        gsMultiPatch<> density, expected, expectedWarm;
        for (int threads : counts)
        {
#ifdef _OPENMP
            omp_set_num_threads(threads);
#else
            (void)threads;
#endif
            gsAdaptiveMultiPatchBuilder builder(reference,2,4,4.,0,0,fixed);
            if (density.empty())
                density = builder.buildAnalyticDensity(gsFunctionExpr<>(dim==2 ? "1+x+2*y" : "1+x+2*y+3*z",dim));
            builder.buildMultiPatch(density);
            if (expected.empty()) expected = builder.MAmapping;
            for (size_t p=0; p<expected.nPatches(); ++p)
                CHECK((expected.patch(p).coefs()-builder.MAmapping.patch(p).coefs()).norm() < 2e-11);
            checkInterfaces(reference,builder.MAmapping);
            checkPositiveJacobian(builder.MAmapping);
            if (fixed) checkExterior(reference,builder.MAmapping);

            builder.m_maxIter = 1;
            builder.buildMultiPatch(density);
            if (expectedWarm.empty()) expectedWarm = builder.MAmapping;
            for (size_t p=0; p<expectedWarm.nPatches(); ++p)
                CHECK((expectedWarm.patch(p).coefs()-builder.MAmapping.patch(p).coefs()).norm() < 2e-11);
            checkInterfaces(reference,builder.MAmapping);
            checkPositiveJacobian(builder.MAmapping);
            if (fixed) checkExterior(reference,builder.MAmapping);
        }
    }
}

}

SUITE(gsAdaptiveMultiPatchBuilder_test)
{
    TEST(NonUnitCylinderMatchesUnitParameterization)
    {
        gsMultiPatch<> unit;
        gsFileData<>("surfaces/cylinder.xml").getId(1,unit);
        auto& unitBasis = dynamic_cast<gsTensorNurbsBasis<2>&>(unit.patch(0).basis());
        for (index_t d = 0; d < 2; ++d)
            unitBasis.source().knots(d).affineTransformTo(0.,1.);
        gsMultiPatch<> original = unit;
        auto& nativeBasis = dynamic_cast<gsTensorNurbsBasis<2>&>(original.patch(0).basis());
        nativeBasis.source().knots(0).affineTransformTo(-2.,46.43999263688814);
        nativeBasis.source().knots(1).affineTransformTo(3.,23.74557318168491);
        const gsMatrix<> originalBounds = original.basis(0).support();
        for (bool fixed : {false,true})
        {
            gsAdaptiveMultiPatchBuilder expected(unit,1,4,3.,0,0,fixed);
            gsAdaptiveMultiPatchBuilder actual(original,1,4,3.,0,0,fixed);
            CHECK((original.basis(0).support()-originalBounds).norm() == 0.);
            checkParameterConversion(original,actual);
            gsFunctionExpr<> f("100+x+2*y+3*z",3);
            const auto rho = actual.buildAnalyticDensity(f);
            const auto expectedRho = expected.buildAnalyticDensity(f);
            CHECK((rho.patch(0).coefs()-expectedRho.patch(0).coefs()).norm() < 1e-9);
            actual.buildMultiPatch(rho);
            expected.buildMultiPatch(expectedRho);
            CHECK((actual.MAmapping.patch(0).coefs()-expected.MAmapping.patch(0).coefs()).norm() < 1e-9);
            CHECK((actual.buildColCompMultiPatch().patch(0).coefs()
                  -expected.buildColCompMultiPatch().patch(0).coefs()).norm() < 1e-8);
            CHECK((actual.buildH1CompMultiPatch(.2).patch(0).coefs()
                  -expected.buildH1CompMultiPatch(.2).patch(0).coefs()).norm() < 1e-8);
        }
    }

    TEST(DifferentPatchParameterBoxes2DAnd3D)
    {
        checkDifferentParameterBoxes<2>();
        checkDifferentParameterBoxes<3>();
    }

    TEST(HierarchicalParameterConversionPreservesRefinement)
    {
        gsMultiPatch<> cylinder;
        gsFileData<>("surfaces/cylinder.xml").getId(1,cylinder);
        auto& tensor = dynamic_cast<gsTensorNurbs<2>&>(cylinder.patch(0));
        tensor.knots(0).affineTransformTo(-4.,19.);
        tensor.knots(1).affineTransformTo(7.,43.);
        gsMultiPatch<> original;
        original.addPatch(gsRationalTHBSpline<2>(tensor));
        original.uniformRefine();
        // Refine a corner at the next level, retaining both coarse and fine functions.
        original.patch(0).refineElements({2,0,0,1,1});
        original.computeTopology();
        gsAdaptiveMultiPatchBuilder builder(original,1,0,3.,0,0,true);
        checkParameterConversion(original,builder);
        gsMultiPatch<> normalized = builder.normalizedGeometry();
        original.uniformRefine();
        normalized.uniformRefine();
        gsAdaptiveMultiPatchBuilder refined(original,0,0,3.,0,0,true);
        CHECK((normalized.patch(0).coefs()-refined.normalizedGeometry().patch(0).coefs()).norm() < 1e-9);
        checkParameterConversion(original,refined);
    }

    TEST(SinglePatchNormalBoundaryCorrectionIsThreadSafe)
    {
        InteriorThreadSettings settings;
        for (index_t dim : {2,3})
        {
            const gsMultiPatch<> reference = dim == 2
                ? gsNurbsCreator<>::BSplineSquareGrid(1,1,1.)
                : gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,-.5,-.5,-.5);
            gsAdaptiveMultiPatchBuilder builder(reference,1,0);
            const gsBasis<>& basis = builder.mapping_basis.basis(0);
            const gsMatrix<> points = basis.anchors();
            gsMatrix<> perturbed = basis.interpolateData(points,points)->coefs();
            perturbed.array() += static_cast<real_t>(0.125);
            gsMatrix<> expected = perturbed;
            for (index_t d=0; d<dim; ++d)
                for (index_t end=0; end<2; ++end)
                {
                    const gsMatrix<index_t> boundary = basis.boundary(boxSide(2*d+end+1));
                    for (index_t i=0; i<boundary.size(); ++i)
                        expected(boundary(i),d) = end;
                }
            for (int threads : {1,4})
            {
#ifdef _OPENMP
                omp_set_num_threads(threads);
#else
                (void)threads;
#endif
                gsMultiPatch<> mapping;
                mapping.addPatch(basis.makeGeometry(perturbed));
                builder.setBoundaryControlPointsAlongNormal(mapping);
                CHECK((mapping.patch(0).coefs()-expected).norm() < 1e-11);
            }
        }
    }

    // Compare with independently assembled operators; exercise OpenMP slices when enabled.
    TEST(FastDiagonalization3DMatchesAssembledOperators)
    {
        gsMultiPatch<> cube = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,-.5,-.5,-.5);
        cube.degreeElevate(2);
        cube.uniformRefine();
        cube.uniformRefine();
        const gsBasis<>& B = cube.basis(0);
        gsBoundaryConditions<> bc, empty;
        gsFunctionExpr<> zero("0",3);
        for (index_t s=1; s<=6; ++s)
            bc.addCondition(0,boxSide(s),condition_type::dirichlet,&zero);
        gsExprAssembler<> A(1,1);
        A.setIntegrationDomain(B.domain());
        auto v = A.getSpace(B);
        A.initSystem();
        A.assemble(grad(v)*grad(v).tr());
        const gsSparseMatrix<> K = A.matrix();
        A.initSystem();
        A.assemble(v*v.tr());
        gsSparseSolver<>::SimplicialLDLT mass;
        mass.compute(A.matrix());
        std::vector<index_t> free;
        for (index_t i=0; i<B.size(); ++i)
        {
            index_t remainder=i;
            bool boundary=false;
            for (index_t d=0; d<3; ++d)
            {
                const index_t size=B.component(d).size(), j=remainder%size;
                remainder/=size;
                boundary=boundary||j==0||j==size-1;
            }
            if (!boundary) free.push_back(i);
        }
        gsMatrix<> rhs(free.size(),1), density(B.size(),1);
        for (index_t i=0; i<rhs.rows(); ++i) rhs(i)=math::sin(static_cast<real_t>(i+1));
        for (index_t i=0; i<density.rows(); ++i) density(i)=math::cos(static_cast<real_t>(i+1));
        gsSparseEntries<> entries;
        for (size_t i=0; i<free.size(); ++i)
            for (size_t j=0; j<free.size(); ++j)
                if (K.coeff(free[i],free[j])!=0.) entries.add(i,j,K.coeff(free[i],free[j]));
        gsSparseMatrix<> Kff(free.size(),free.size());
        Kff.setFrom(entries);
        gsSparseSolver<>::SimplicialLDLT stiffness;
        stiffness.compute(Kff);
        const gsMatrix<> expected = stiffness.solve(rhs), expectedMass = mass.solve(density);
        gsPatchPreconditionersCreator<>::Poisson_FastDiag poisson(B,bc,A.options(),0.);
        gsPatchPreconditionersCreator<>::Poisson_FastDiag projector(B,empty,A.options(),0.);
#ifdef _OPENMP
        const int originalThreads = omp_get_max_threads();
        const int originalDynamic = omp_get_dynamic();
        omp_set_dynamic(0);
        const std::vector<int> threadCounts = {1,4};
#else
        const std::vector<int> threadCounts = {1};
#endif
        for (int threads : threadCounts)
        {
#ifdef _OPENMP
            omp_set_num_threads(threads);
#else
            (void)threads;
#endif
            for (int repeat=0; repeat<10; ++repeat)
            {
                const gsMatrix<> solution = poisson.solve(rhs);
                const gsMatrix<> projected = projector.L2ProjectScalar(density);
                CHECK((solution-expected).norm()<1e-10*std::max(static_cast<real_t>(1),expected.norm()));
                CHECK((projected-expectedMass).norm()<1e-10*std::max(static_cast<real_t>(1),expectedMass.norm()));
            }
        }
        std::vector<gsMatrix<> > solutions(16), projections(16);
        #pragma omp parallel for num_threads(4)
        for (index_t i=0; i<16; ++i)
        {
            const real_t factor = 1.+static_cast<real_t>(i)/10.;
            solutions[i] = poisson.solve(gsMatrix<>(factor*rhs));
            projections[i] = projector.L2ProjectScalar(gsMatrix<>(factor*density));
        }
        for (index_t i=0; i<16; ++i)
        {
            const real_t factor = 1.+static_cast<real_t>(i)/10.;
            CHECK((solutions[i]-factor*expected).norm()<1e-10*std::max(static_cast<real_t>(1),factor*expected.norm()));
            CHECK((projections[i]-factor*expectedMass).norm()<1e-10*std::max(static_cast<real_t>(1),factor*expectedMass.norm()));
        }
#ifdef _OPENMP
        omp_set_num_threads(originalThreads);
        omp_set_dynamic(originalDynamic);
#endif
    }

    TEST(ParallelPatchInteriors2DMatchSerial)
    {
        checkParallelInteriors(2);
    }

    TEST(ParallelPatchInteriors3DMatchSerial)
    {
        checkParallelInteriors(3);
    }

    TEST(InteriorFailurePropagatesToCaller)
    {
        InteriorThreadSettings settings;
#ifdef _OPENMP
        omp_set_num_threads(2);
#endif
        // Disconnected cubes keep every face fixed, so the invalid monitor is
        // first used by a volume solve inside the parallel interior stage.
        gsMultiPatch<> reference = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,-.5,-.5,-.5);
        gsGeometry<>::uPtr other = reference.patch(0).clone();
        other->coefs().col(0).array() += 2.;
        reference.addPatch(give(other));
        reference.computeTopology();
        gsAdaptiveMultiPatchBuilder builder(reference,1,1,2.,0,0,true);
        gsMultiPatch<> density = builder.buildAnalyticDensity(gsFunctionExpr<>("1",3));
        density.patch(1).coefs().setConstant(std::numeric_limits<real_t>::quiet_NaN());
        bool failed = false;
        try { builder.buildMultiPatch(density); }
        catch (const std::exception&) { failed = true; }
        CHECK(failed);
        CHECK(builder.MAmapping.empty());
    }

    TEST(ConstantVolumeMonitorPreservesIdentity)
    {
        const gsMultiPatch<> reference = gsNurbsCreator<>::BSplineCubeGrid(2,1,1,1.,-.5,-.5,-.5);
        gsAdaptiveMultiPatchBuilder builder(reference, 2, 3, 8.);
        builder.buildMultiPatch(builder.buildAnalyticDensity(gsFunctionExpr<>("1",3)));
        checkIdentity(builder.MAmapping);
        checkInterfaces(reference, builder.MAmapping);
    }

    TEST(FixedVolumeFacesAndThreeCoordinateUpdates)
    {
        gsMultiPatch<> reference = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,-.5,-.5,-.5);
        // Unequal degrees and span counts exercise tensor ordering in fast diagonalization.
        reference.patch(0).degreeElevate(1,1);
        reference.patch(0).degreeElevate(2,2);
        gsAdaptiveMultiPatchBuilder builder(reference, 2, 1, 4., 0, 0, true);
        builder.buildMultiPatch(builder.buildAnalyticDensity(gsFunctionExpr<>("1+x+2*y+3*z",3)));
        checkExterior(reference, builder.MAmapping);
        const auto& patch = builder.MAmapping.patch(0);
        const gsMatrix<> points = patch.basis().anchors();
        const gsMatrix<> identity = patch.basis().interpolateData(points,points)->coefs();
        for (index_t d = 0; d < 3; ++d)
            CHECK((patch.coefs().col(d)-identity.col(d)).norm() > 1e-6);
        // Independently assemble the stated Picard equation on the free rows.
        // Raw rho=1+x+2y+3z scales from [1,7] to [1,4], and C=4+5.
        gsExprAssembler<> A(1,1);
        A.setIntegrationDomain(patch.basis().domain());
        auto v = A.getSpace(patch.basis());
        gsFunctionExpr<> monitor("1+0.5*x+y+1.5*z",3);
        auto rho = A.getCoeff(monitor);
        std::vector<bool> boundary(patch.coefsSize(),false);
        for (index_t s = 1; s <= 6; ++s)
        {
            const gsMatrix<index_t> indices = patch.basis().boundary(boxSide(s));
            for (index_t i = 0; i < indices.size(); ++i) boundary[indices(i)] = true;
        }
        const std::array<std::string,3> coordinates = {{"x","y","z"}};
        for (index_t d = 0; d < 3; ++d)
        {
            gsFunctionExpr<> coordinate(coordinates[d],3);
            auto old = A.getCoeff(coordinate);
            A.initSystem();
            A.assemble(grad(v)*grad(v).tr(), (1.-rho.val()/9.)*grad(v)*grad(old).tr());
            const gsMatrix<> residual = A.matrix()*patch.coefs().col(d)-A.rhs();
            for (index_t i = 0; i < residual.rows(); ++i)
                if (!boundary[i]) CHECK(math::abs(residual(i)) < 1e-11);
        }
        checkPositiveJacobian(builder.MAmapping);
    }

    TEST(OneSweepMovesEdgesFacesAndVolumesWithOrientedInterface)
    {
        gsMultiPatch<> reference = gsNurbsCreator<>::BSplineCubeGrid(1,1,1,1.,-.5,-.5,-.5);
        gsMultiPatch<> neighbor = reference;
        const gsMatrix<> original = neighbor.patch(0).coefs();
        neighbor.patch(0).coefs().col(0) = original.col(2).array()+1.;
        neighbor.patch(0).coefs().col(1) = 1.-original.col(1).array();
        neighbor.patch(0).coefs().col(2) = original.col(0);
        reference.addPatch(neighbor.patch(0));
        reference.computeTopology();
        CHECK_EQUAL(size_t(1), reference.nInterfaces());
        gsAdaptiveMultiPatchBuilder builder(reference, 2, 1, 4.);
        builder.buildMultiPatch(builder.buildAnalyticDensity(gsFunctionExpr<>("1+x+2*y+3*z",3)));
        checkInterfaces(reference,builder.MAmapping);
        const auto& patch = builder.MAmapping.patch(0);
        const gsMatrix<> points = patch.basis().anchors();
        const gsMatrix<> displacement = patch.eval(points)-points;
        std::array<real_t,3> moved = {{0.,0.,0.}};
        for (index_t q = 0; q < points.cols(); ++q)
        {
            index_t fixedCoordinates = 0;
            for (index_t d = 0; d < 3; ++d)
                if (points(d,q) < 1e-12 || points(d,q) > 1.-1e-12) ++fixedCoordinates;
            if (fixedCoordinates < 3)
                moved[fixedCoordinates] = std::max(moved[fixedCoordinates], displacement.col(q).norm());
        }
        for (index_t d = 0; d < 3; ++d) CHECK(moved[d] > 1e-6);
        checkPositiveJacobian(builder.MAmapping);
    }

    TEST(FourPatchSharedEdgesAndMarkedVolumeDensity)
    {
        const gsMultiPatch<> reference = gsNurbsCreator<>::BSplineCubeGrid(2,2,1,1.,-.5,-.5,-.5);
        gsAdaptiveMultiPatchBuilder builder(reference, 2, 2, 3.);
        builder.buildMultiPatch(builder.buildAnalyticDensity(gsFunctionExpr<>("1+x+2*y+3*z",3)));
        checkInterfaces(reference,builder.MAmapping);
        checkPositiveJacobian(builder.MAmapping);
        gsMultiBasis<> basis(builder.MAmapping);
        std::vector<bool> marked(basis.totalElements(),false);
        for (size_t i = 0; i < marked.size(); i += 3) marked[i] = true;
        const gsMultiPatch<> density = builder.buildDensity(basis,marked);
        CHECK_EQUAL(reference.nPatches(),density.nPatches());
        builder.buildMultiPatch(density);
        checkInterfaces(reference,builder.MAmapping);
    }

    TEST(ExactExteriorWithMovingInternalFace)
    {
        const gsMultiPatch<> reference = gsNurbsCreator<>::BSplineCubeGrid(2,1,1,1.,-.5,-.5,-.5);
        gsAdaptiveMultiPatchBuilder builder(reference, 2, 2, 4., 0, 0, true);
        builder.buildMultiPatch(builder.buildAnalyticDensity(gsFunctionExpr<>("1+x+2*y+3*z",3)));
        checkExterior(reference,builder.MAmapping);
        checkInterfaces(reference,builder.MAmapping);
        checkPositiveJacobian(builder.MAmapping);
        const auto& face = reference.interfaces().front().first();
        const auto& patch = builder.MAmapping.patch(face.patch);
        const gsMatrix<> points = patch.basis().anchors();
        const gsMatrix<> identity = patch.basis().interpolateData(points,points)->coefs();
        const gsMatrix<index_t> indices = patch.basis().boundary(face.side());
        real_t change = 0.;
        for (index_t i = 0; i < indices.size(); ++i)
            change += (patch.coefs().row(indices(i))-identity.row(indices(i))).squaredNorm();
        CHECK(change > 1e-8);
    }

    TEST(VolumeWarmStartAndZeroIterationBudget)
    {
        const gsMultiPatch<> reference = gsNurbsCreator<>::BSplineCubeGrid(2,1,1,1.,-.5,-.5,-.5);
        gsAdaptiveMultiPatchBuilder incremental(reference, 2, 1, 4.);
        const gsMultiPatch<> density = incremental.buildAnalyticDensity(gsFunctionExpr<>("1+x+2*y+3*z",3));
        incremental.buildMultiPatch(density);
        const gsMultiPatch<> first = incremental.MAmapping;
        incremental.m_maxIter = 0;
        incremental.buildMultiPatch(density);
        for (size_t p = 0; p < first.nPatches(); ++p)
            CHECK((first.patch(p).coefs()-incremental.MAmapping.patch(p).coefs()).norm() == 0.);
        incremental.m_maxIter = 1;
        incremental.buildMultiPatch(density);
        gsAdaptiveMultiPatchBuilder together(reference, 2, 2, 4.);
        together.buildMultiPatch(density);
        for (size_t p = 0; p < first.nPatches(); ++p)
            CHECK((together.MAmapping.patch(p).coefs()-incremental.MAmapping.patch(p).coefs()).norm() < 1e-10);
    }
}
