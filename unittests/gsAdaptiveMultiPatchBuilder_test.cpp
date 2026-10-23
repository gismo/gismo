/** @file gsAdaptiveMultiPatchBuilder_test.cpp
    @brief Regression tests for hierarchical 3D moving mesh updates.
    This file is part of the G+Smo library, under the Mozilla Public License v2.0.

    Author(s): M. BAHARI
    */
#include "gismo_unittest.h"
#include <gsAssembler/gsAdaptiveMultiPatchBuilder.h>

namespace
{
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
        gsVector<index_t> counts(3);
        counts.setConstant(5);
        counts(first.direction()) = 1;
        gsMatrix<> box(3,2);
        box.col(0).setZero();
        box.col(1).setOnes();
        box.row(first.direction()).setConstant(first.parameter());
        const gsMatrix<> points = gsGridIterator<real_t,CUBE>(box, counts).toMatrix();
        gsMatrix<> target(3, points.cols());
        target.row(second.direction()).setConstant(second.parameter());
        for (index_t d = 0; d < 3; ++d)
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
    gsVector<index_t> counts(3);
    counts.setConstant(7);
    for (size_t p = 0; p < mapping.nPatches(); ++p)
    {
        const gsMatrix<> points = gsGridIterator<real_t,CUBE>(mapping.basis(p).support(), counts).toMatrix();
        const gsMatrix<> J = mapping.patch(p).jacobian(points);
        for (index_t q = 0; q < points.cols(); ++q)
            CHECK(J.middleCols(3*q,3).determinant() > 0.);
    }
}
}

SUITE(gsAdaptiveMultiPatchBuilder_test)
{
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
        CHECK_EQUAL(1, reference.nInterfaces());
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
