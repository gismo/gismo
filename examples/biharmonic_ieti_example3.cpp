/** @file biharmonic_ieti_example.cpp

    @brief Biharmonic example for an extremely simple multipatch domain

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <ctime>
#include <gismo.h>

using namespace gismo;

gsDofMapper setupTwoLayerDofMapper(const gsMultiPatch<>& mp, const gsMultiBasis<>& mb)
{
    gsVector<index_t> patchDofSizes(mp.nPatches());
    for (size_t k=0; k<mp.nPatches(); ++k)
        patchDofSizes[k] = mb[k].size();
    
    gsDofMapper dm(patchDofSizes);
    for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const index_t k1 = it->first().patch;
        const index_t k2 = it->second().patch;
        gsVector<index_t> s1 = mb.basis(k1).boundary(it->first().side()),
                          s2 = mb.basis(k2).boundary(it->second().side()),
                          s1o = mb.basis(k1).boundaryOffset(it->first().side(),1),
                          s2o = mb.basis(k2).boundaryOffset(it->second().side(),1);
        
        // We assume for now that the orientation matches!
        GISMO_ASSERT( s1.rows() == s2.rows() && s1.rows() == s1o.rows() && s2.rows() == s2o.rows(), ""); 
        for (index_t i=0;i<s1.rows();++i)
        {
            dm.matchDof(k1,s1[i],k2,s2[i]);
            dm.matchDof(k1,s1o[i],k2,s2o[i]);
        }
    }
    dm.finalize();
    return dm;
}

std::pair<std::vector<index_t>,std::vector<index_t>> commonSkeletonDofs(const gsSparseMatrix<real_t,RowMajor>& jump1, const gsSparseMatrix<real_t,RowMajor>& jump2 )
{
    GISMO_ASSERT(jump1.rows()==jump2.rows(),"");
    std::vector<std::pair<index_t,index_t>> tmp(jump1.rows());
    for (index_t i=0; i<jump1.outerSize(); ++i)
        for ( gsSparseMatrix<real_t,RowMajor>::InnerIterator it(jump1, i); it; ++it)
             if (it.value())
                 tmp[it.row()].first = 1+it.col();
    for (index_t i=0; i<jump2.outerSize(); ++i)
        for ( gsSparseMatrix<real_t,RowMajor>::InnerIterator it(jump2, i); it; ++it)
             if (it.value())
                 tmp[it.row()].second = 1+it.col();         
    std::pair<std::vector<index_t>,std::vector<index_t>> result;
    for (index_t i=0; i<jump1.rows(); ++i)
        if (tmp[i].first&&tmp[i].second)
        {
            result.first .push_back(tmp[i].first -1);
            result.second.push_back(tmp[i].second-1);
        }
    return result;
}

gsSparseMatrix<> embeddingForChosen( const std::vector<index_t>& idxes, index_t cols )
{
    // Assuming no duplicates in idxes!
    gsSparseMatrix<> result(idxes.size(), cols);

    // TODO: reserve gsSparseMatrix<>
    for (size_t i=0; i<idxes.size(); ++i)
    {
        GISMO_ASSERT (idxes[i]>=0 && idxes[i]<cols, "Out of range");
        result(i,idxes[i])=1;
    }
    result.makeCompressed();
    return result;
}

gsSparseMatrix<> embeddingForNotChosen( std::vector<index_t> idxes, index_t cols )
{
    std::sort(idxes.begin(), idxes.end());
    
    for (size_t i=1; i<idxes.size(); ++i)
    {
        GISMO_ASSERT( idxes[i-1]<idxes[i], "Assume to be unique." );
    }
    
    std::vector<index_t> remaining_idxes;
    remaining_idxes.reserve(cols-idxes.size());
    idxes.push_back(cols);
    size_t i=0;
    for (index_t j=0; j<cols; ++j)
    {
        while (idxes[i]<j)
            ++i;
        if (idxes[i]>j)
            remaining_idxes.push_back(j);
    }
    return embeddingForChosen(remaining_idxes, cols);
}



std::vector<index_t> commonLMultpliersDofs(const gsSparseMatrix<real_t,RowMajor>& jump1, const gsSparseMatrix<real_t,RowMajor>& jump2 )
{
    GISMO_ASSERT(jump1.rows()==jump2.rows(),"");
    gsVector<char> tmp(jump1.rows());
    tmp.setZero();
    for (index_t i=0; i<jump1.outerSize(); ++i)
        for ( gsSparseMatrix<real_t,RowMajor>::InnerIterator it(jump1, i); it; ++it)
             if (it.value())
                 tmp[it.row()] = 1;
    for (index_t i=0; i<jump2.outerSize(); ++i)
        for ( gsSparseMatrix<real_t,RowMajor>::InnerIterator it(jump2, i); it; ++it)
             if (it.value())
                 tmp[it.row()] |= 2;
    std::vector<index_t> result;
    for (index_t i=0; i<jump1.rows(); ++i)
        if (tmp[i]==3)
            result.push_back(i);
    return result;
}

void add(gsSparseEntries<real_t>& se, index_t row, index_t col, const gsSparseMatrix<real_t>& mat )
{
    for (index_t i=0; i<mat.outerSize(); ++i)
        for ( gsSparseMatrix<real_t>::InnerIterator it(mat, i); it; ++it)
             se.add( it.row()+row, it.col()+col, it.value() );
}


std::vector<std::vector<std::pair<index_t,gsSparseVector<>>>> cornersFromJumpMatrices( const gsDofMapper& dm, const std::vector<gsSparseMatrix<real_t,RowMajor>>& sms )
{
    gsVector<index_t> tmp( dm.freeSize() );
    tmp.setZero();
    for (size_t k=0; k<sms.size(); ++k)
        for (index_t i=0; i<sms[k].outerSize(); ++i)
            for (gsSparseMatrix<real_t,RowMajor>::InnerIterator it(sms[k],i); it; ++it)
                if (dm.is_free(it.col(), k))
                    tmp[dm.index(it.col(), k)] += 1;
    index_t j=0;
    for (index_t i=0; i<tmp.rows(); ++i)
        if (tmp[i]>2)
        {
                tmp[i]=j;
            ++j;
        }
        else
            tmp[i]=-1;
    std::vector<std::vector<std::pair<index_t,gsSparseVector<>>>> result(j);
    for (size_t k=0; k<sms.size(); ++k)
        for (index_t i=0; i<sms[k].outerSize(); ++i)
            for (gsSparseMatrix<real_t,RowMajor>::InnerIterator it(sms[k],i); it; ++it)
                if (tmp[dm.index(it.col(),k)]>=0)
                {
                    std::pair<index_t,gsSparseVector<>> constr(k, gsSparseVector<>(sms[k].cols()));
                    constr.second[it.col()] = 1;
                    GISMO_ASSERT (it.col()<sms[k].cols(), "");
                    std::vector<std::pair<index_t,gsSparseVector<>>>& where
                        = result[tmp[dm.index(it.col(),k)]];
                    bool found = false;
                    for (size_t l=0; l<where.size(); ++l)
                        found |= (where[l].first == constr.first /*&& where[l].second == constr.second*/);
                    if (!found)
                        where.push_back(constr);
                }
    return result;
}

gsSparseMatrix<> makeTransformer(const gsBasis<>& basis)
{
    const index_t d = basis.dim();
    gsSparseMatrix<> result;
    for (index_t i=0; i<d; ++i)
    {
        const index_t ndofs1D = basis.component(d-1-i).size();
        GISMO_ASSERT( ndofs1D>5, "" );
        
        gsSparseMatrix<> transformer1D(ndofs1D,ndofs1D);
        transformer1D.setIdentity();
        
        transformer1D(0,0)=1;
        transformer1D(0,1)=1;

        transformer1D(1,1)=1;
        
        transformer1D(ndofs1D-2,ndofs1D-2)=-1;

        transformer1D(ndofs1D-1,ndofs1D-2)=1;
        transformer1D(ndofs1D-1,ndofs1D-1)=1;

        if (i==0)
            result = give(transformer1D);
        else
            result = result.kron(transformer1D);
    }
    return result;
}


int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    index_t nPatchesX = 4;
    index_t nPatchesY = 4;
    index_t degree = 2;
    index_t refinements = 1;
    index_t rhsType = 2;
    std::string out;
    real_t robin = 100;
    //real_t factor = 0;
    bool plot = false;

    real_t tolerance = 1.e-8;
    index_t maxIterations = 1000;
    std::string primals("");

    gsCmdLine cmd("Biharmonic IETI example for an extremely simple multipatch domain.");
    cmd.addInt   ("t", "RhsType",               "Chosen right-hand side", rhsType);
    cmd.addInt   ("x", "PatchesX",              "Number of patches (coordinate direction x)", nPatchesX);
    cmd.addInt   ("y", "PatchesY",              "Number of patches (coordinate direction y)", nPatchesY);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addReal  ("b", "Robin",                 "Penalty parameter for Robin boundary conditions", robin);
    //cmd.addReal  ("q", "Factor",                "Factor for bdy derivs", factor);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);
    
    cmd.addString("",  "Primals",               "Chosen primal dofs; empty is default", primals);
    cmd.addReal  ("",  "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Stopping criterion for linear solver", maxIterations);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Run biharmonic_ieti_example with options:\n" << cmd << "\n";

    /********************** Define rhs **********************/
    const index_t dim = 2;
    const char *rhsTypes[] =
        {
            "32*pi^4*sin(2*pi*x)*sin(2*pi*y)",
            "2*pi^4*sin(pi*x)*sin(pi*y)",
            "1/8*pi^4*sin(pi*x/2)*sin(pi*y/2)"
        };
    
    gsInfo << "Rhs function is " << rhsTypes[rhsType] << "\n";
    gsFunctionExpr<> f(rhsTypes[rhsType],dim);

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;
    gsMultiPatch<> mp;
    for (index_t i=0; i<nPatchesX; ++i)
        for (index_t j=0; j<nPatchesY; ++j)
            mp.addPatch(gsNurbsCreator<>::BSplineRectangle(i,j,i+1,j+1));
    mp.computeTopology(); 

    gsInfo << "done.\n";
    const index_t nPatches = mp.nPatches();

    /************ Setup bases and adjust degree *************/

    gsMultiBasis<> mb(mp);
    gsInfo << "Setup bases and adjust degree... " << std::flush;
    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();
    gsInfo << "done.\n";

    /******************* Setup dofMapper ********************/
    gsInfo << "Setup dofMapper... " << std::flush;
    gsDofMapper dm = setupTwoLayerDofMapper(mp, mb);
    gsInfo << "done:\n" << dm << "\n";
    
    /****************** Setup ietimapper ********************/
    gsInfo << "Setup ietimapper... " << std::flush;
    gsMatrix<> fixedPart;
    gsIetiMapper<> ietiMapper(mb,dm,fixedPart);
    ietiMapper.computeJumpMatrices(/*fullyRedundant=*/true,/*excludeCorners=*/false);
    if (primals == "c")
        ietiMapper.cornersAsPrimals();
    else if (primals == "x")
    {
        std::vector<gsSparseMatrix<real_t,RowMajor>> jm;
        jm.reserve(nPatches);
        for (index_t k=0; k<nPatches; ++k)
            jm.push_back(ietiMapper.jumpMatrix(k));
        std::vector<std::vector<std::pair<index_t,gsSparseVector<>>>> data = cornersFromJumpMatrices(dm, jm);

        if ( true )
        {
            gsInfo << "[";
            for (size_t i=0; i<data.size(); ++i)
            {
                gsInfo << "[";
                for (size_t j=0; j<data[i].size(); ++j)
                {
                    gsInfo << " {" << data[i][j].first << "," << gsSparseVector<>::InnerIterator(data[i][j].second,0).row() << "}";
                }
                gsInfo << " ]\n";
            }
            gsInfo << "]\n\n";
        }

        for (size_t i=0; i<data.size(); ++i)
            ietiMapper.customPrimalConstraints(data[i]);
        gsInfo << "[" << data.size() << " eXtended cornerdofs added]";
    }
    else
        GISMO_ENSURE(primals=="", "Unknown choice \""<<primals<<"\"");
    
    gsIetiSystem<> ieti;
    ieti.reserve(nPatches+1);

    gsScaledDirichletPrec<> prec;
    prec.reserve(nPatches);

    gsPrimalSystem<> primal(ietiMapper.nPrimalDofs());
    
    gsInfo << "done.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    std::vector<gsSparseMatrix<>> localBasisTransforms; localBasisTransforms.reserve(nPatches);
    std::vector<gsSparseMatrix<>> localStiffnessMatrices; localStiffnessMatrices.reserve(nPatches);
    std::vector<gsMatrix<>> localRhsVectors; localRhsVectors.reserve(nPatches);
    std::vector<gsLinearOperator<>::Ptr> localSchurs; localSchurs.reserve(nPatches);
    
    for (index_t k=0; k<nPatches; ++k)
    {
        gsInfo << "[" << k << "] " << std::flush; 
        gsMultiPatch<> mp_local(mp[k]);
        gsMultiBasis<> mb_local(mb[k]);
        
        //! [Problem setup]
        gsExprAssembler<real_t> A(1,1);

        // Elements used for numerical integration
        A.setIntegrationElements(mb_local);
        gsExprEvaluator<real_t> ev(A);
    
        // Set the geometry map
        auto G = A.getMap(mp_local);
        auto ff = A.getCoeff(f, G);
        auto u = A.getSpace(mb_local);
    
        A.initSystem();
        A.assemble(ilapl(u, G) * ilapl(u, G).tr() * meas(G),u * ff * meas(G));
        A.assemble(u*u.tr()*meas(G)); //TODO
                       
        gsSparseMatrix<> transformer = makeTransformer(mb[k]);
        
        // Fetch data
        gsSparseMatrix<real_t, RowMajor> jumpMatrix  = ietiMapper.jumpMatrix(k);
        gsSparseMatrix<>                 localMatrix = transformer*A.matrix()*transformer.transpose();
        gsMatrix<>                       localRhs    = transformer*A.rhs();
        
        GISMO_ASSERT(jumpMatrix.cols() == localMatrix.rows(), "");
        
        // Penalize Dirichlet boundary
        /*for (gsBoxTopology::const_biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
        {
            if (it->patchIndex()==k)
            {
                // gsInfo << "Found biterator: " << *it << "\n";
                gsVector<index_t> s1 = mb_local.basis(0).boundary(it->side());
                for (index_t i=0; i<s1.rows(); ++i)
                    localMatrix(s1[i],s1[i]) += robin;
            }
        }*/
        
        // Store
        localStiffnessMatrices.push_back(localMatrix);
        localRhsVectors.push_back(localRhs);
        localBasisTransforms.push_back(transformer);
        localSchurs.push_back(gsScaledDirichletPrec<>::schurComplement(localMatrix,ietiMapper.skeletonDofs(k)));
        
        prec.addSubdomain(
            gsScaledDirichletPrec<>::restrictToSkeleton(
                jumpMatrix,
                localMatrix,
                ietiMapper.skeletonDofs(k)
            )
        );
        

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
        gsInfo << "[P] " << std::flush;
        gsLinearOperator<>::Ptr localSolver
            = makeSparseLUSolver(primal.localMatrix());

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

    // This is the main cg iteration
    //! [Solve]
    gsConjugateGradient<> PCG( ieti.schurComplement(), prec.preconditioner() );
    PCG.setOptions( cmd.getGroup("Solver") );
    PCG.setCalcEigenvalues(true);
    PCG.solveDetailed( rhsForSchur, lambda, errorHistory );
    //! [Solve]

    gsInfo << "done.\n    Reconstruct solution from Lagrange multipliers... " << std::flush;
    // Now, we want to have the global solution for u
    //! [Recover]
    std::vector<gsMatrix<>> localSol = primal.distributePrimalSolution(
            ieti.constructSolutionFromLagrangeMultipliers(lambda)
        );
    gsMatrix<> globalSol = ietiMapper.constructGlobalSolutionFromLocalSolutions(localSol);
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

    gsInfo << "Estimated condition number: " << PCG.getConditionNumber() << "\n";

    /********************** Output **************************/    
    if (!out.empty())
    {
        gsFileData<> fd;
        std::time_t time = std::time(NULL);
        fd.add(cmd);
        for (index_t k=0; k<nPatches; ++k)
            fd.add(gsMatrix<>(localBasisTransforms[k].transpose() * localSol[k]));
        
        for (index_t k=0; k<nPatches; ++k)
            fd.add(gsSparseMatrix<>(ietiMapper.jumpMatrix(k)));
        
        for (index_t k=0; k<nPatches; ++k)
            fd.add(localStiffnessMatrices[k]);

        for (index_t k=0; k<nPatches; ++k)
        {
            gsMatrix<> tmp;
            localSchurs[k]->toMatrix(tmp);
            fd.add(tmp);
        }
        
        gsMatrix<> sc;
        ieti.schurComplement()->toMatrix(sc);
        fd.add(sc);
        
        gsMatrix<> pc;
        prec.preconditioner()->toMatrix(pc);
        fd.add(pc);
        
        fd.addComment(std::string("biharmonic_ieti_example   Timestamp:")+std::ctime(&time));
        fd.save(out);
        gsInfo << "Write solution to file " << out << "\n";
    }

    if (plot)
    {
        gsInfo << "Write Paraview data to file ieti_result.pvd\n";
        gsMultiPatch<> mpsol;
        for (index_t k=0; k<nPatches; ++k)
            mpsol.addPatch( mb[k].makeGeometry( localBasisTransforms[k].transpose() * localSol[k] ) );
        gsWriteParaview<>( gsField<>( mp, mpsol ), "ieti_result", 1000);
    }
    if (!plot&&out.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --out to write solution to xml file.\n";
    }
 
    return 0;
}



