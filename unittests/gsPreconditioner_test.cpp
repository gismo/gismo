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
    // Define Geometry
    gsMultiPatch<> mp( *gsNurbsCreator<>::NurbsQuarterAnnulus() );

    // Create mulibasis
    gsMultiBasis<> mb(mp);

    // Refine multibais
    const index_t numRefine = 4;
    for (int i = 0; i < numRefine; ++i)
        mb.uniformRefine();

    // Define Boundary conditions
    gsConstantFunction<> one(1,mp.geoDim());
    gsBoundaryConditions<> bc;
    bc.addCondition( boundary::west,  condition_type::neumann,   &one );
    bc.addCondition( boundary::east,  condition_type::neumann,   &one );
    bc.addCondition( boundary::south, condition_type::neumann,   &one );
    bc.addCondition( boundary::north, condition_type::dirichlet, &one );

    // Initilize Assembler and assemble
    gsOptionList opt = gsAssembler<>::defaultOptions();
    gsPoissonAssembler<> assembler(
        mp,
        mb,
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
        gsConjugateGradient<> solver(mat, gsPatchPreconditionersCreator<>::fastDiagonalizationOp(mb[0],bc));
        solver.setTolerance( 1.e-8 );
        solver.setMaxIterations( 35 );
        solver.solve(rhs,sol);
        CHECK ( solver.error() <= solver.tolerance() );
    }
    else if (testcase==3)
    {
        const real_t h = mb[0].getMinCellLength();
        gsGenericAssembler<> gAssembler(mp,mb,opt,&bc);
        mat += (1/(h*h)) * gAssembler.assembleMass();
        gsConjugateGradient<> solver(mat, gsPatchPreconditionersCreator<>::subspaceCorrectedMassSmootherOp(mb[0],bc));
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
        // Define Geometry
        gsMultiPatch<> mp( *gsNurbsCreator<>::BSplineSquare() );

        // Create mulibasis
        gsMultiBasis<> mb(mp);

        // Refine multibasis
        mb.uniformRefine();
        dynamic_cast< gsTensorBSplineBasis<2>& >(mb[0]).component(0).uniformRefine();

        // Set degree
        mb[0].setDegreePreservingMultiplicity(3);

        // Define Boundary conditions
        gsConstantFunction<> one(1,mp.geoDim());
        gsBoundaryConditions<> bc;
        bc.addCondition( boundary::west,  condition_type::neumann,   &one );
        bc.addCondition( boundary::east,  condition_type::neumann,   &one );
        bc.addCondition( boundary::south, condition_type::neumann,   &one );
        bc.addCondition( boundary::north, condition_type::dirichlet, &one );

        // Initilize Assembler and assemble
        gsOptionList opt = gsAssembler<>::defaultOptions();
        gsPoissonAssembler<> assembler(
            mp,
            mb,
            bc,
            one,
            (dirichlet::strategy) opt.getInt("DirichletStrategy"),
            (iFace::strategy) opt.getInt("InterfaceStrategy")
            );
        assembler.assemble();
        const gsSparseMatrix<>& stiff0 = assembler.matrix();

        // Get stiffness matrix from gsPatchPreconditionersCreator
        gsSparseMatrix<> stiff1 = gsPatchPreconditionersCreator<>::stiffnessMatrix(mb[0],bc,opt);
        CHECK ( (stiff1-stiff0 ).norm() < 1/(real_t)(10000) );

        // Get stiffness matrix operator from gsPatchPreconditionersCreator
        gsMatrix<> stiff2;
        gsPatchPreconditionersCreator<>::stiffnessMatrixOp(mb[0],bc,opt)->toMatrix(stiff2);
        CHECK ( (stiff2-stiff0 ).norm() < 1/(real_t)(10000) );

        // Get inverse of stiffness matrix from gsPatchPreconditionersCreator
        gsLinearOperator<>::Ptr stiffInv = gsPatchPreconditionersCreator<>::fastDiagonalizationOp(mb[0],bc,opt);
        gsMatrix<> result;
        stiffInv->apply(stiff0,result);
        for (index_t i=0; i<result.rows(); ++i) result(i,i)-=1;
        CHECK ( result.norm() < 1/(real_t)(10000) );
    }

    TEST(gsPatchPreconditioner_mass_test)
    {
        // Define Geometry
        gsMultiPatch<> mp( *gsNurbsCreator<>::BSplineSquare() );

        // Create mulibasis
        gsMultiBasis<> mb(mp);

        // Refine multibasis
        mb.uniformRefine();
        dynamic_cast< gsTensorBSplineBasis<2>& >(mb[0]).component(0).uniformRefine();

        // Set degree
        mb[0].setDegreePreservingMultiplicity(3);

        // Define Boundary conditions
        gsConstantFunction<> one(1,mp.geoDim());
        gsBoundaryConditions<> bc;
        bc.addCondition( boundary::west,  condition_type::neumann,   &one );
        bc.addCondition( boundary::east,  condition_type::neumann,   &one );
        bc.addCondition( boundary::south, condition_type::neumann,   &one );
        bc.addCondition( boundary::north, condition_type::dirichlet, &one );

        // Initilize Assembler and assemble
        gsOptionList opt = gsAssembler<>::defaultOptions();
        gsGenericAssembler<> assembler(
            mp,
            mb,
            opt,
            &bc
            );
        gsSparseMatrix<> mass0 = assembler.assembleMass();


        // Get mass matrix from gsPatchPreconditionersCreator
        gsSparseMatrix<> mass1 = gsPatchPreconditionersCreator<>::massMatrix(mb[0],bc,opt);
        CHECK ( ( mass0-mass1 ).norm() < 1/(real_t)(10000) );

        // Get mass matrix operator from gsPatchPreconditionersCreator
        gsMatrix<> mass2;
        gsPatchPreconditionersCreator<>::massMatrixOp(mb[0],bc,opt)->toMatrix(mass2);
        CHECK ( ( mass0-mass2 ).norm() < 1/(real_t)(10000) );

        // Get inverse of mass matrix from gsPatchPreconditionersCreator
        gsLinearOperator<>::Ptr massInv = gsPatchPreconditionersCreator<>::massMatrixInvOp(mb[0],bc,opt);
        gsMatrix<> result;
        massInv->apply(mass0,result);
        for (index_t i=0; i<result.rows(); ++i) result(i,i)-=1;
        CHECK ( result.norm() < 1/(real_t)(10000) );
    }

    TEST(gsAdditiveOp_test)
    {
        gsSparseMatrix<real_t,RowMajor> t1(3,2);
        t1(0,0)=1;
        t1(1,1)=1;
        gsSparseMatrix<real_t,RowMajor> t2(3,1);
        t2(2,0)=1;

        std::vector< gsSparseMatrix<real_t,RowMajor> > t;
        t.reserve(2);
        t.push_back(t1);
        t.push_back(t2);

        gsMatrix<> o1(2,2);
        o1 << 1,0,   0,2;

        gsMatrix<> o2(1,1);
        o2 << 3;

        std::vector< gsLinearOperator<>::Ptr > o;
        o.reserve(2);
        o.push_back(makeMatrixOp(o1));
        o.push_back(makeMatrixOp(o2));

        gsMatrix<> in(3,1);
        in << 9,8,7;
        gsMatrix<> out(3,1);
        out << 9,16,21;

        {
            gsSumOp<> s;
            for (size_t i=0; i<t.size(); ++i)
            s.addOperator(
                gsProductOp<>::make(
                    makeMatrixOp(t[i].transpose()),
                    o[i],
                    makeMatrixOp(t[i])
                )
            );
            gsMatrix<> res;
            s.apply( in, res );
            CHECK ( (res-out).norm() < 1/(real_t)(10000) );
        }

        {
            gsAdditiveOp<> a(t,o);
            gsMatrix<> res;
            a.apply( in, res );
            CHECK ( (res-out).norm() < 1/(real_t)(10000) );
        }
    }


}
