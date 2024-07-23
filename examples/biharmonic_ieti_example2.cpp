/** @file biharmonic_ieti_example.cpp

    @brief 

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <ctime>
#include <gismo.h>

using namespace gismo;


int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    index_t nPatchesX = 1;
    index_t nPatchesY = 1;
    index_t degree = 2;
    index_t refinements = 1;
    
    real_t tolerance = 1.e-8;
    index_t maxIterations = 100;
    std::string out;
    real_t alpha = 0;
    bool plot = false;
    bool calcEigenvalues = ! true; // TODO

    gsCmdLine cmd("Solves a PDE with an isogeometric discretization using a multigrid solver.");
    cmd.addInt   ("x", "PatchesX",              "Number of patches (coordinate direction x)", nPatchesX);
    cmd.addInt   ("y", "PatchesY",              "Number of patches (coordinate direction y)", nPatchesY);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Stopping criterion for linear solver", maxIterations);
    cmd.addReal  ("a", "Alpha",                 "", alpha);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Run biharmonic_ieti_example with options:\n" << cmd << "\n";

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    //! [Define Geometry]    
    gsMultiPatch<> mp;
    for (index_t i=0; i<nPatchesX; ++i)
        for (index_t j=0; j<nPatchesY; ++j)
            mp.addPatch(gsNurbsCreator<>::BSplineRectangle(i,j,i+1,j+1));
    mp.computeTopology(); 

    gsInfo << "done.\n";

    const index_t nPatches = mp.nPatches();
    GISMO_ASSERT( nPatches == nPatchesX*nPatchesY, "");

    /************ Setup bases and adjust degree *************/

    gsMultiBasis<> mb(mp);
    gsInfo << "Setup bases and adjust degree... " << std::flush;
    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();
    gsInfo << "done.\n";

    /******************* Setup dofmapper ********************/
    gsInfo << "Setup dofmapper... " << std::flush;
    gsVector<index_t> patchDofSizes(nPatches);
    for (index_t k=0; k<nPatches; ++k)
        patchDofSizes[k] = mb[k].size();
    
    gsInfo << "\n\npatchDofSizes="<<patchDofSizes.transpose()<<"\n\n";
    
    gsDofMapper dm(patchDofSizes);
    for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        gsInfo << *it << "\n";
        const index_t k1 = it->first().patch;
        const index_t k2 = it->second().patch;

        gsVector<index_t> s1 = mb.basis(k1).boundary(it->first().side()),
                          s2 = mb.basis(k2).boundary(it->second().side()),
                          s1o = mb.basis(k1).boundaryOffset(it->first().side(),1),
                          s2o = mb.basis(k2).boundaryOffset(it->second().side(),1);
        
        gsInfo << "firstside: " << s1.transpose() << "\n";
        gsInfo << "firstoffset: " << s1o.transpose() << "\n";
        gsInfo << "secondside: " << s2.transpose() << "\n";
        gsInfo << "secondoffset: " << s2o.transpose() << "\n";
        gsInfo << "orient: " << it->dirOrientation(it->first())[1-it->first().direction()] << "\n"; ///??
        // TODO: assume for now that the orientation matches!
        
        GISMO_ASSERT( s1.rows() == s2.rows() && s1.rows() == s1o.rows() && s2.rows() == s2o.rows(), ""); 
        
        
        for (index_t i=0;i<s1.rows();++i)
        {
            dm.matchDof(k1,s1[i],k2,s2[i]);
            dm.matchDof(k1,s1o[i],k2,s2o[i]);
        }
        
        //mb.matchInterface(*it, mapper);
    }

    dm.finalize();

    gsInfo << "done.\n";
    /****************** Setup ietimapper ********************/
    gsInfo << "Setup ietimapper... " << std::flush;
    gsMatrix<> fixedPart;
    gsIetiMapper<> ietiMapper(mb,dm,fixedPart);
    ietiMapper.computeJumpMatrices(false,false); // nothing special for corners
    // TODO: mapper -> primals
    
    gsIetiSystem<> ieti;
    ieti.reserve(nPatches+1);

    gsScaledDirichletPrec<> prec;
    prec.reserve(nPatches);

    gsPrimalSystem<> primal(ietiMapper.nPrimalDofs());
    
    gsInfo << "done.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    for (index_t k=0; k<nPatches; ++k)
    {
        gsInfo << "[" << k << "] " << std::flush; 
        gsMultiPatch<> mp_local(mp[k]);
        gsMultiBasis<> mb_local(mb[k]);
        
        //! [Problem setup]
        gsExprAssembler<real_t> A(1,1);
        //gsInfo<<"Active options:\n"<< A.options() <<"\n";
    
        // Elements used for numerical integration
        A.setIntegrationElements(mb_local);
        gsExprEvaluator<real_t> ev(A);
    
        const index_t dim = 2;
        gsFunctionExpr<> f("sin(2*pi*x)*sin(2*pi*y)",dim);
    
        // Set the geometry map
        auto G = A.getMap(mp_local);
        auto ff = A.getCoeff(f, G);
        auto u = A.getSpace(mb_local);
    
        A.initSystem();
        A.assemble(ilapl(u, G) * ilapl(u, G).tr() * meas(G),u * ff * meas(G));
        A.assemble(u * u.tr() * meas(G)); //TODO
    
        //gsInfo << "\n\nmatrix:\n" << A.matrix().toDense() << "\n\n";
        //gsInfo << "\n\nvector:\n" << A.rhs().transpose() << "\n\n";
        
        const index_t ndofs = A.matrix().rows();
        const index_t ndofs1D = math::sqrt(ndofs);
        GISMO_ASSERT( ndofs==ndofs1D*ndofs1D, ndofs<<"=="<<ndofs1D<<"*"<<ndofs1D );
        GISMO_ASSERT( ndofs1D>3, "" );
        
        gsSparseMatrix<> transformer1D(ndofs1D,ndofs1D);
        transformer1D.setIdentity();
        const real_t h=1./ndofs1D;
        transformer1D(0,1)=1;
        transformer1D(ndofs1D-1,ndofs1D-2)=1;
	if (alpha==0) alpha = h/degree;
        transformer1D(1,1)=alpha;//h/degree;
        transformer1D(ndofs1D-2,ndofs1D-2)=-alpha;//-h/degree;
        gsInfo << "\n\ntransformer1D:\n" << transformer1D.toDense() << "\n\n";
        
        gsSparseMatrix<> transformer=transformer1D.kron(transformer1D);
        //gsInfo << "\n\ntransformer:\n" << transformer.toDense() << "\n\n";
        
        //gsInfo << "\n\ntransformed matrix:\n" << (transformer*A.matrix()*transformer.transpose()).toDense() << "\n\n";
        
        // Fetch data
        gsSparseMatrix<real_t, RowMajor> jumpMatrix  = ietiMapper.jumpMatrix(k);
        gsSparseMatrix<>                 localMatrix = (transformer*A.matrix()*transformer.transpose());
        gsMatrix<>                       localRhs    = transformer*A.rhs();
        //! [Assemble]

        gsInfo << "\njump\n" << jumpMatrix.toDense() <<"\n\n";

        gsMatrix<> schurAsMat, schurInvAsMat;
        {
        gsMatrix<> tmp; makeSparseLUSolver(localMatrix)->apply(jumpMatrix.transpose().toDense(), tmp);
        schurAsMat = jumpMatrix*tmp;
        gsInfo << "\nschurF\n" << schurAsMat << "\n\n";
        makePartialPivLUSolver(schurAsMat)->toMatrix(schurInvAsMat);
	gsInfo << "\nschurFinv\n" << schurInvAsMat << "\n\n";

	}

        gsInfo << "dofs: "<< localMatrix.rows() <<"\n";
        gsInfo << "skeleton dofs: "<< ietiMapper.skeletonDofs(k).size() <<"\n";
        for (index_t i=0; i<ietiMapper.skeletonDofs(k).size(); ++i) gsInfo << ietiMapper.skeletonDofs(k)[i] << " ";
	gsInfo << "\n";

        // Add the patch to the scaled Dirichlet preconditioner
        //! [Patch to preconditioner]
        auto sdp =
            gsScaledDirichletPrec<>::restrictToSkeleton(
                jumpMatrix,
                localMatrix,
                ietiMapper.skeletonDofs(k)
            );
	{ gsInfo << "\npreconder embed\n" << sdp.first << "\n\n"; }
	{ gsMatrix<> tmp; sdp.second->toMatrix(tmp); gsInfo << "\npreconder\n" << tmp << "\n\n"; }
	{ gsMatrix<> tmp; sdp.second->apply(sdp.first*schurAsMat*sdp.first.transpose(),tmp); gsInfo<<"\npreconder2schur\n" << tmp << "\n\n";}
        prec.addSubdomain(give(sdp));
        //! [Patch to preconditioner]

        // This function writes back to jumpMatrix, localMatrix, and localRhs,
        // so it must be called after prec.addSubdomain().
        //! [Patch to primals]
        primal.handleConstraints(
            ietiMapper.primalConstraints(k),
            ietiMapper.primalDofIndices(k),
            jumpMatrix,
            localMatrix,
            localRhs
        );
        //! [Patch to primals]

        // Add the patch to the Ieti system
        //! [Patch to system]
        ieti.addSubdomain(
            jumpMatrix.moveToPtr(),
            makeMatrixOp(localMatrix.moveToPtr()),
            give(localRhs)
        );
        //! [Patch to system]
    //! [End of assembling loop]
    } // end for
    //! [End of assembling loop]

    // Add the primal problem if there are primal constraints
    //! [Primal to system]
    if (ietiMapper.nPrimalDofs()>0)
    {
        gsInfo << "Primal.\n";
        // It is not required to provide a local solver to .addSubdomain,
        // since a sparse LU solver would be set up on the fly if required.
        // Here, we make use of the fact that we can use a Cholesky solver
        // because the primal problem is symmetric and positive definite:
        gsLinearOperator<>::Ptr localSolver
            = makeSparseCholeskySolver(primal.localMatrix());

        ieti.addSubdomain(
            primal.jumpMatrix().moveToPtr(),
            makeMatrixOp(primal.localMatrix().moveToPtr()),
            give(primal.localRhs()),
            localSolver
        );
    }
    //! [Primal to system]

    gsInfo << "done. " << ietiMapper.nPrimalDofs() << " primal dofs.\n";

    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... \n"
        "    Setup multiplicity scaling... " << std::flush;

    // Tell the preconditioner to set up the scaling
    //! [Setup scaling]
    prec.setupMultiplicityScaling();
    //prec.setupIdentityScaling();
    //! [Setup scaling]

    gsInfo << "done.\n    Setup rhs... " << std::flush;
    // Compute the Schur-complement contribution for the right-hand-side
    //! [Setup rhs]
    gsMatrix<> rhsForSchur = ieti.rhsForSchurComplement();
    //! [Setup rhs]

    gsInfo << "done.\n    Setup cg solver for Lagrange multipliers and solve... " << std::flush;
    // Initial guess
    //! [Define initial guess]
    gsMatrix<> lambda;
    lambda.setRandom( ieti.nLagrangeMultipliers(), 1 );
    //! [Define initial guess]

    gsMatrix<> errorHistory;
 gsMatrix<> tmpii; ieti.schurComplement()->toMatrix(tmpii); gsInfo << "\nGlobalSchur=\n"<<tmpii<<"\n\n"; 
{ gsMatrix<> tmp; prec.preconditioner()->toMatrix(tmp); gsInfo << "\nGlobalPrec=\n"<<tmp<<"\n\n"; }
{ gsMatrix<> tmp; prec.preconditioner()->apply(tmpii,tmp);  gsInfo << "\nGlobalPrec2GlobalSchur=\n"<<tmp<<"\n\n"; }

    // This is the main cg iteration
    //! [Solve]
    gsConjugateGradient<> PCG( ieti.schurComplement(), prec.preconditioner() );
    PCG.setOptions( cmd.getGroup("Solver") ).solveDetailed( rhsForSchur, lambda, errorHistory );
    //! [Solve]

    gsInfo << "done.\n    Reconstruct solution from Lagrange multipliers... " << std::flush;
    // Now, we want to have the global solution for u
    //! [Recover]
    gsMatrix<> uVec = ietiMapper.constructGlobalSolutionFromLocalSolutions(
        primal.distributePrimalSolution(
            ieti.constructSolutionFromLagrangeMultipliers(lambda)
        )
    );
    //! [Recover]
    gsInfo << "done.\n\n";

    /******************** Print end Exit ********************/

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

    if (!out.empty())
    {
        gsFileData<> fd;
        std::time_t time = std::time(NULL);
        fd.add(cmd);
        fd.add(uVec);
        fd.addComment(std::string("ieti_example   Timestamp:")+std::ctime(&time));
        fd.save(out);
        gsInfo << "Write solution to file " << out << "\n";
    }

    if (plot)
    {
         //TODO: take care of basis transformations!
    }
    if (!plot&&out.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --out to write solution to xml file.\n";
    }
 
    return 0;
}


