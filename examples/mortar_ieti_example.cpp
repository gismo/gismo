/** @file mortaring_example.cpp

    @brief Provides an example for mortar.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <ctime>

#include <gismo.h>
#include <gsAssembler/gsJumpAssembler.h>

using namespace gismo;

gsMatrix<> constructGlobalSolutionFromLocalSolutions(const std::vector<gsMatrix<>>& local);

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/annulus_dg2.xml");
    index_t refinements = 3;
    index_t degree = 2;
    bool nitsche = false;
    std::string boundaryConditions("d");
    real_t tolerance = 1.e-8;
    index_t maxIterations = 100;
    bool calcEigenvalues = false;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addString("g", "Geometry",              "Geometry file", geometry);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addSwitch("",  "Nitsche",               "Use Nitsche method for Dirichlet boundary conditions", nitsche);
    cmd.addString("b", "BoundaryConditions",    "Boundary conditions", boundaryConditions);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Maximum iterations for linear solver", maxIterations);
    cmd.addSwitch("",  "Solver.CalcEigenvalues","Estimate eigenvalues based on Lanczos", calcEigenvalues);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Run mortaring_example with options:\n" << cmd << "\n";

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    //! [Define Geometry]
    gsMultiPatch<>::uPtr mpPtr = gsReadFile<>(geometry);
    //! [Define Geometry]
    if (!mpPtr)
    {
        gsInfo << "No geometry found in file " << geometry << ".\n";
        return EXIT_FAILURE;
    }
    gsMultiPatch<>& mp = *mpPtr;

    //! [Define Geometry2]
    //for (index_t i=0; i<splitPatches; ++i)
    //{
    //    gsInfo << "split patches uniformly... " << std::flush;
    //    mp = mp.uniformSplit();
    //}
    //! [Define Geometry2]

    gsInfo << "done.\n";

    /************** Define boundary conditions **************/

    gsInfo << "Define boundary conditions... " << std::flush;

    //! [Define Source]
    // Right-hand-side
    gsFunctionExpr<> f( "2*sin(x)*cos(y)", mp.geoDim() );

    // Dirichlet function
    gsFunctionExpr<> gD( "sin(x)*cos(y)", mp.geoDim() );

    // Neumann
    gsConstantFunction<> gN( 1.0, mp.geoDim() );

    gsBoundaryConditions<> bc;
    //! [Define Source]
    {
        const index_t len = boundaryConditions.length();
        index_t i = 0;
        for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
        {
            char b_local;
            if ( len == 1 )
                b_local = boundaryConditions[0];
            else if ( i < len )
                b_local = boundaryConditions[i];
            else
            {
                gsInfo << "\nNot enough boundary conditions given.\n";
                return EXIT_FAILURE;
            }

            if ( b_local == 'd' )
                bc.addCondition( *it, condition_type::dirichlet, &gD );
            else if ( b_local == 'n' )
                bc.addCondition( *it, condition_type::neumann, &gN );
            else
            {
                gsInfo << "\nInvalid boundary condition given; only 'd' (Dirichlet) and 'n' (Neumann) are supported.\n";
                return EXIT_FAILURE;
            }

            ++i;
        }
        if ( len > i )
            gsInfo << "\nToo many boundary conditions have been specified. Ignoring the remaining ones.\n";
        gsInfo << "done. "<<i<<" boundary conditions set.\n";
    }


    /************ Setup bases and adjust degree *************/

    //! [Define Basis]
    gsMultiBasis<> mb(mp);
    //! [Define Basis]

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    //! [Set degree and refine]
    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();
    //! [Set degree and refine]

    gsInfo << "done.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    //! [Assemble]
    gsPoissonAssembler<> assembler(
        mp,
        mb,
        bc,
        f,
        nitsche ? dirichlet::nitsche : dirichlet::elimination,
        iFace::none // iFace::dg
    );
    assembler.assemble();
    //! [Assemble]

    gsSparseMatrix<> mat = assembler.matrix();

    std::vector<gsSparseMatrix<>> constraints;
    index_t nConstraints = 0;
    for ( typename gsMultiPatch<>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it )
    {
        gsJumpAssembler<real_t> ifAssembler(
            mp,
            mb,
            bc,
            f,
            nitsche ? dirichlet::nitsche : dirichlet::elimination,
            it
        );
        ifAssembler.assemble();
        
        // TODO: assume all patches to have same number of dofs!
        size_t dofsPerPatch = mat.rows()/8;
        size_t k1 = it->first().patch;
        gsSparseMatrix<> tmp = ifAssembler.matrix().middleRows(k1*dofsPerPatch,dofsPerPatch);
        std::vector<index_t> rowsToBeChosen;
        for (index_t i=0; i<tmp.rows(); ++i)
            for (index_t j=0; j<tmp.cols(); ++j)
                if (tmp(i,j)!=0)
                {
                    rowsToBeChosen.push_back(i);
                    break;
                }
        gsSparseMatrix<> chser(rowsToBeChosen.size(), tmp.rows());
        for (size_t i=0; i<rowsToBeChosen.size(); ++i)
            chser(i, rowsToBeChosen[i]) = 1;
        
        //gsInfo << "chser " << chser.rows() << "x" << chser.cols() << "\n";
        //gsInfo << "tmp " << tmp.rows() << "x" << tmp.cols() << "\n";
        //gsInfo << "ifAssembler.matrix()  " << ifAssembler.matrix() .rows() << "x" << ifAssembler.matrix() .cols() << "\n";
        
        gsSparseMatrix<> loc = chser * tmp.middleCols(k1*dofsPerPatch,dofsPerPatch) * chser.transpose();
        //gsInfo << "loc=" << loc << "\n\n";
        gsMatrix<> locinv;
        makeSparseCholeskySolver(loc)->toMatrix(locinv);
        gsSparseMatrix<> locinv2(locinv.rows(), locinv.cols());
        for (index_t i=0; i<locinv.rows(); ++i)
            for (index_t j=0; j<locinv.cols(); ++j)
                locinv2(i,j) = locinv(i,j);
        
        constraints.push_back(locinv2 * chser * tmp);
        nConstraints += rowsToBeChosen.size();
        
    }
        
    gsSparseMatrix<> allConstraints(nConstraints,assembler.matrix().cols());
    {
        index_t r = 0;
        for (size_t i=0; i<constraints.size(); ++i)
            for (index_t j=0; j<constraints[i].rows(); ++j)
            {
                for (index_t c=0; c<constraints[i].cols(); ++c)
                    if (constraints[i](j,c)*constraints[i](j,c)>1e-6)
                        allConstraints(r,c) = constraints[i](j,c);
                ++r;
            }
    }
    
    //gsInfo << "Constraints\n" << allConstraints << "\n\n";
    
    gsInfo << "done: "<<assembler.matrix().rows()<<" dofs.\n";

    /**************** Setup solver and solve ****************/

    gsIetiSystem<> ietiSystem;
    gsScaledDirichletPrec<> prec;
    ietiSystem.reserve(mp.nPatches());
    prec.reserve(mp.nPatches());
    for (size_t k=0; k<mp.nPatches(); ++k)
    {
        size_t dofsPerPatch = mat.rows()/8;
        gsSparseMatrix<> localMat = assembler.matrix().middleRows(k*dofsPerPatch,dofsPerPatch).middleCols(k*dofsPerPatch,dofsPerPatch);
        gsMatrix<> localRhs = assembler.rhs().middleRows(k*dofsPerPatch,dofsPerPatch);
        gsSparseMatrix<real_t,RowMajor> localConstr = allConstraints.middleCols(k*dofsPerPatch,dofsPerPatch);
    
        prec.addSubdomain(gsScaledDirichletPrec<>::restrictToSkeleton(localConstr,localMat,gsScaledDirichletPrec<>::skeletonDofs(localConstr)));

        auto solver = makeSparseCholeskySolver(localMat);
        ietiSystem.addSubdomain(localConstr.moveToPtr(), makeMatrixOp(localMat.moveToPtr()), localRhs, give(solver));
        // TODO: primalDofs
    }
   
    prec.setupMultiplicityScaling(); //TODO: scaling

    gsInfo << "Setup solver and solve... " << std::flush;

    gsMatrix<> lambda;
    lambda.setRandom( ietiSystem.nLagrangeMultipliers(), 1 );
    //! [Define initial guess]

    gsMatrix<> errorHistory;

    // This is the main cg iteration
    //! [Solve]
    gsConjugateGradient<> PCG( ietiSystem.schurComplement(), prec.preconditioner() );
    PCG.setOptions( cmd.getGroup("Solver") ).solveDetailed( ietiSystem.rhsForSchurComplement(), lambda, errorHistory );

    gsMatrix<> x = constructGlobalSolutionFromLocalSolutions(
      // primal.distributePrimalSolution(
        ietiSystem.constructSolutionFromLagrangeMultipliers(lambda)
      //)
    );
    
    gsInfo << "done.\n\n";


    const index_t iter = errorHistory.rows()-1;
    const bool success = errorHistory(iter,0) < tolerance;
    if (success)
        gsInfo << "Reached desired tolerance after " << iter << " iterations:\n";
    else
        gsInfo << "Did not reach desired tolerance after " << iter << " iterations:\n";

    if (errorHistory.rows() < 20)
        gsInfo << errorHistory.transpose() << "\n\n";
    else
        gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

    if (calcEigenvalues)
        gsInfo << "Estimated condition number: " << PCG.getConditionNumber() << "\n";


    /******************** Print end Exit ********************/

    if (!out.empty())
    {
        std::ofstream outfile (out, std::ios_base::app);
        outfile << "mortar_ieti_example\t"
                << geometry << "\t"
                << degree << "\t"
                << refinements << "\t"
                << PCG.getConditionNumber() << "\t"
                << iter << "\t"
                /*<< "***" << "\t"
                << sPatchesX << "\t"
                << sPatchesY << "\t"
                << rhsType << "\t"
                << robin << "\t"
                << alpha << "\t"
                << bdyConds << "\t"
                << primals << "\t"
                << dualPreconder << "\t"
                << extremelyDeluxeParameter << "\t"
                << solverType <<*/ "\n";
        gsInfo << "Write solution to file " << out << "\n";
    }

    if (plot)
    {
        // Construct the solution as a scalar field
        gsMultiPatch<> mpsol;
        assembler.constructSolution(x, mpsol);
        gsField<> sol( assembler.patches(), mpsol );

        // Write solution to paraview files
        gsInfo << "Write Paraview data to file mortaring_result.pvd\n";
        gsWriteParaview<>(sol, "mortaring_result", 1000);
        //gsFileManager::open("mortaring_result.pvd");
    }
    if (!plot&&out.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --out to write solution to xml file.\n";
    }
    return EXIT_SUCCESS;
}

gsMatrix<> constructGlobalSolutionFromLocalSolutions(const std::vector<gsMatrix<>>& local)
{
    gsMatrix<> result(local.size() * local[0].rows(), 1);//TODO
    for (size_t k=0; k<local.size(); ++k)
        result.middleRows(k*local[0].rows(),local[0].rows()) = local[k];
    return result;
}