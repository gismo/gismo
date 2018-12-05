/** @file gsPreconditioner_test.cpp

    @brief Tests for various preconditioners for Poisson problems.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn, S. Takacs
*/

#include "gismo_unittest.h"
#include <gsSolver/gsSolverUtils.h>

using namespace gismo;

void runPreconditionerTest( index_t testcase )
{
    index_t numRefine = 4;

    // Define Geometry (Unit square with 4 patches)
    gsMultiPatch<> patches( *gsNurbsCreator<>::NurbsQuarterAnnulus() );

    gsConstantFunction<> one(1,2);

    // Define Boundary conditions
    gsBoundaryConditions<> bc;
    bc.addCondition( boundary::west,  condition_type::neumann,   &one );
    bc.addCondition( boundary::east,  condition_type::neumann,   &one );
    bc.addCondition( boundary::south, condition_type::neumann,   &one );
    bc.addCondition( boundary::north, condition_type::dirichlet, &one );

    // Copy bases for refinement
    gsMultiBasis<> bases( patches );

    // Define discretization space by initial refining the basis of the geometry
    for (int i = 0; i < numRefine; ++i)
        bases.uniformRefine();

    // Initilize Assembler and assemble
    gsOptionList opt = gsAssembler<>::defaultOptions();
    gsPoissonAssembler<> assembler(
        patches,
        bases,
        bc,
        one,
        (dirichlet::strategy) opt.getInt("DirichletStrategy"),
        (iFace::strategy) opt.getInt("InterfaceStrategy")
        );

    assembler.assemble();
    gsSparseMatrix<> mat = assembler.matrix();
    gsMatrix<> rhs = assembler.rhs();
    gsMatrix<> sol;
    sol.setRandom(rhs.rows(), rhs.cols());

    if (testcase==0)
    {
        gsConjugateGradient<> solver(mat, makeJacobiOp(mat));
        solver.setTolerance( 1.e-8 );
        solver.setMaxIterations( 110 );
        solver.solve(rhs,sol);
        CHECK ( solver.error() <= solver.tolerance() );
    }
    else if (testcase==1)
    {
        gsConjugateGradient<> solver(mat, makeSymmetricGaussSeidelOp(mat));
        solver.setTolerance( 1.e-8 );
        solver.setMaxIterations( 35 );
        solver.solve(rhs,sol);
        CHECK ( solver.error() <= solver.tolerance() );
    }
    else if (testcase==2)
    {
        gsConjugateGradient<> solver(mat, gsPatchPreconditionersCreator<>::fastDiagonalizationOp(bases[0],bc));
        solver.setTolerance( 1.e-8 );
        solver.setMaxIterations( 95 );
        solver.solve(rhs,sol);
        CHECK ( solver.error() <= solver.tolerance() );
    }
    else if (testcase==3)
    {
        const real_t h = bases[0].getMinCellLength();
        gsGenericAssembler<> gAssembler(patches,bases,opt,&bc);
        mat += (1/(h*h)) * gAssembler.assembleMass();
        gsConjugateGradient<> solver(mat, gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(bases[0],bc));
        solver.setTolerance( 1.e-8 );
        solver.setMaxIterations( 50 );
        solver.solve(rhs,sol);
        CHECK ( solver.error() <= solver.tolerance() );
    }
}


SUITE(gsPreconditioner_test)
{

    TEST(gsJacobiPreconditioner_test)
    {
        runPreconditionerTest(0);
    }
    TEST(gsSymmetricGaussSeidelPreconditioner_test)
    {
        runPreconditionerTest(1);
    }
    TEST(gsFastDiagonalizationPreconditioner_test)
    {
        runPreconditionerTest(2);
    }
    TEST(gsSubspaceCorrectedMassPreconditioner_test)
    {
        runPreconditionerTest(3);
    }

    TEST(gsPatchPreconditioner_stiff_test)
    {
        gsGeometry<>::uPtr geo = gsNurbsCreator<>::BSplineSquare();
        gsMultiPatch<> mp(*geo);
        gsMultiBasis<> mb(mp);
        mb.uniformRefine();
        mb[0].setDegreePreservingMultiplicity(3);

        gsBoundaryConditions<> bc;
        gsConstantFunction<> one(1,mp.geoDim());
        bc.addCondition( boundary::west,  condition_type::neumann,   &one );
        bc.addCondition( boundary::east,  condition_type::neumann,   &one );
        bc.addCondition( boundary::south, condition_type::neumann,   &one );
        bc.addCondition( boundary::north, condition_type::dirichlet, &one );

        gsOptionList opt = gsAssembler<>::defaultOptions();

        gsSparseMatrix<> stiff1 = gsPatchPreconditionersCreator<>::stiffnessMatrix(mb[0],bc,opt);
        gsMatrix<> stiff2;
        gsPatchPreconditionersCreator<>::stiffnessMatrixOp(mb[0],bc,opt)->toMatrix(stiff2);


        gsPoissonAssembler<> assembler(
            mp,
            mb,
            bc,
            one,
            (dirichlet::strategy) opt.getInt("DirichletStrategy"),
            (iFace::strategy) opt.getInt("InterfaceStrategy")
            );
        assembler.assemble();

        CHECK ( (stiff1-assembler.matrix() ).norm() < 1/real_t(10000) );
        CHECK ( (stiff2-assembler.matrix() ).norm() < 1/real_t(10000) );
    }

    TEST(gsPatchPreconditioner_mass_test)
    {
        gsGeometry<>::uPtr geo = gsNurbsCreator<>::BSplineSquare();
        gsMultiPatch<> mp(*geo);
        gsMultiBasis<> mb(mp);
        mb.uniformRefine();
        mb[0].setDegreePreservingMultiplicity(3);

        gsBoundaryConditions<> bc;
        gsConstantFunction<> one(1,mp.geoDim());
        bc.addCondition( boundary::west,  condition_type::neumann,   &one );
        bc.addCondition( boundary::east,  condition_type::neumann,   &one );
        bc.addCondition( boundary::south, condition_type::neumann,   &one );
        bc.addCondition( boundary::north, condition_type::dirichlet, &one );

        gsOptionList opt = gsAssembler<>::defaultOptions();

        gsGenericAssembler<> assembler(
            mp,
            mb,
            opt,
            &bc
            );
        gsSparseMatrix<> mass0 = assembler.assembleMass();


        gsSparseMatrix<> mass1 = gsPatchPreconditionersCreator<>::massMatrix(mb[0],bc,opt);
        CHECK ( ( mass0-mass1 ).norm() < 1/real_t(10000) );

        gsMatrix<> mass2;
        gsPatchPreconditionersCreator<>::massMatrixOp(mb[0],bc,opt)->toMatrix(mass2);
        CHECK ( ( mass0-mass2 ).norm() < 1/real_t(10000) );


        gsLinearOperator<>::Ptr massInv = gsPatchPreconditionersCreator<>::massMatrixInvOp(mb[0],bc,opt);
        gsMatrix<> result;
        massInv->apply(mass2,result);
        for (index_t i=0; i<result.rows(); ++i) result(i,i)-=1;
        CHECK ( result.norm() < 1/real_t(10000) );
    }


}
