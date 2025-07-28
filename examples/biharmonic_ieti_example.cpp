/** @file biharmonic_ieti_example.cpp

    @brief Biharmonic example for an extremely simple multipatch domain

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#include <gismo.h>

using namespace gismo;

// Helper functions, defined at the end of the file
gsMultiPatch<>::uPtr tryGetRectangularGeometry(const char *s);
gsDofMapper setupC1DofMapper(const gsMultiPatch<>& mp, const gsMultiBasis<>& mb, unsigned bdyConds);
std::vector<index_t> setupC1interfaceDofsVector(const gsBasis<>& basis, const gsDofMapper& dm, index_t type);
gsSparseMatrix<> setupC1basisTransformation(const gsBasis<>& basis, const gsDofMapper& dm);

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/fat_quarter_annulus.xml");
    index_t patchSplitsX = 1;
    index_t patchSplitsY = 1;
    std::string rhs("1/8*pi^4*sin(pi*x/2)*sin(pi*y/2)");
    int bdyConds = 2;
    index_t degree = 2;
    index_t multiplicity = 1;
    index_t refinements = 2;
    real_t tolerance = 1.e-8;
    index_t maxIterations = 1000;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("Biharmonic IETI example for an extremely simple multipatch domain.");

    cmd.addString("g", "Geometry",              "Chosen geometry file", geometry);
    cmd.addInt   ("x", "PatchSplitsX",          "Number of splits (coordinate direction x)", patchSplitsX);
    cmd.addInt   ("y", "PatchSplitsY",          "Number of splits (coordinate direction y)", patchSplitsY);
    cmd.addString("f", "RhsFunction",           "Chosen right-hand function", rhs);
    cmd.addInt   ("b", "BdyConds",              "Bounday conditions: (1) second biharmonic (u, \u0394u); (2) first biharmonic (u, \u2202nu)", bdyConds);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addInt   ("m", "Multiplicity",          "Multiplicity of knots for B-spline discretization space", multiplicity);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addReal  ("",  "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Stopping criterion for linear solver", maxIterations);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo << "Run biharmonic_ieti_example with options:\n" << cmd << "\n";

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;
    gsMultiPatch<>::uPtr mpPtr;
    if (!mpPtr)
        mpPtr = tryGetRectangularGeometry(geometry.c_str());
    if (!mpPtr)
        mpPtr = gsReadFile<>(geometry);
    if (!mpPtr)
    {
        gsInfo << "No geometry found in file " << geometry << ".\n";
        return EXIT_FAILURE;
    }
    gsMultiPatch<>& mp = *mpPtr;

    for (index_t i=0; i<patchSplitsX; ++i)
        mp = mp.uniformSplit(0);
    for (index_t j=0; j<patchSplitsY; ++j)
        mp = mp.uniformSplit(1);

    const index_t nPatches = mp.nPatches();
    gsInfo << "done: " << nPatches << " patches, " << mp.interfaces().size() << " interfaces.\n";

    /* Setup bases and adjust grid size, degree and smoothness */
    gsInfo << "Setup bases and adjust grid size, degree and smoothness... " << std::flush;
    gsMultiBasis<> mb(mp);
    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinements; ++i )
        mb.uniformRefine();

    for ( size_t i = 0; i < mb.nBases(); ++ i )
    {
        if (multiplicity>1 && multiplicity<degree)
            mb[i].reduceContinuity(multiplicity-1);
        else if (multiplicity!=1)
        {
            gsInfo << "Multiplicity must be at least 1 and at most degree-1\n";
            return EXIT_FAILURE;
        }
    }
    gsInfo << "done.\n";

    /*** Setup dofMapper and local basis transformations ****/

    if (bdyConds < 1 || bdyConds > 2)
    {
        gsInfo << "Invalid choice for --BdyCond (-b).\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Setup dofMapper and local basis transformations... " << std::flush;
    gsDofMapper dm = setupC1DofMapper(mp, mb, bdyConds);
    gsInfo << "done:\n" << dm << "\n";

    /********************** Define rhs **********************/

    gsInfo << "Rhs function is " << rhs << "\n";
    gsFunctionExpr<> f(rhs,mp.geoDim());

    /****************** Setup ietimapper ********************/

    gsInfo << "Setup ietimapper... " << std::flush;
    gsMatrix<> fixedPart;
    fixedPart.setZero(dm.boundarySize(),1);
    gsIetiMapper<> ietiMapper(mb,dm,fixedPart);
    ietiMapper.computeJumpMatrices(false, false);

    for (index_t k=0; k<nPatches; ++k)
    {
        std::vector<index_t> cornerDofs = setupC1interfaceDofsVector(mb[k], ietiMapper.dofMapperLocal(k),2);
        for (size_t i=0; i<cornerDofs.size(); ++i)
            ietiMapper.declareDofAsPrimal(k, cornerDofs[i], true);
    }

    gsIetiSystem<> ieti;
    ieti.reserve(nPatches+1);

    gsScaledDirichletPrec<> prec;
    prec.reserve(nPatches);

    gsPrimalSystem<> primal(ietiMapper.nPrimalDofs());

    gsInfo << "done: " << ietiMapper.nPrimalDofs() << " primals.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrix... " << std::flush;

    std::vector<gsSparseMatrix<>> localBasisTransforms; localBasisTransforms.reserve(nPatches);

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

        u.setup();
        ietiMapper.initFeSpace(u,k);

        A.initSystem();
        //A.assemble(ilapl(u, G) * ilapl(u, G).tr() * meas(G), u * ff * meas(G));
        A.assemble(ihess(u, G) % ihess(u, G).tr() * meas(G), u * ff * meas(G));

        gsSparseMatrix<> transformer = setupC1basisTransformation(mb[k],ietiMapper.dofMapperLocal(k));

        // Fetch data
        gsSparseMatrix<real_t, RowMajor> jumpMatrix  = ietiMapper.jumpMatrix(k);
        gsSparseMatrix<>                 localMatrix = transformer*A.matrix()*transformer.transpose();
        gsMatrix<>                       localRhs    = transformer*A.rhs();

        GISMO_ASSERT(jumpMatrix.cols() == localMatrix.rows(), "");

        // Store
        localBasisTransforms.push_back(transformer);

        for (index_t j=0; j<2; ++j)
        {
            prec.addSubdomain(
                gsScaledDirichletPrec<>::restrictToSkeleton(
                    jumpMatrix,
                    localMatrix,
                    setupC1interfaceDofsVector(mb_local[0], ietiMapper.dofMapperLocal(k), j)
                )
            );
        }

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
        gsLinearOperator<>::Ptr localSolver = makeSparseLUSolver(primal.localMatrix());

        ieti.addSubdomain(
            gsSparseMatrix<real_t,RowMajor>(primal.jumpMatrix()).moveToPtr(),
            makeMatrixOp(gsSparseMatrix<>(primal.localMatrix()).moveToPtr()),
            give(primal.localRhs()),
            localSolver
        );
    }
    //! [Primal to system]

    gsInfo << "done. " << ietiMapper.nPrimalDofs() << " primal dofs.\n";


    /**************** Setup solver and solve ****************/

    gsInfo << "Setup solver and solve... \n";

    // Tell the preconditioner to set up the scaling
    //! [Setup scaling]
    gsLinearOperator<>::Ptr preconder;
    gsInfo << "    Setup sclaled Dirichlet preconder... " << std::flush;
    prec.setupMultiplicityScaling();
    preconder = prec.preconditioner();

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
    gsConjugateGradient<> solver( ieti.schurComplement(), preconder );
    solver.setOptions( cmd.getGroup("Solver") );
    solver.setCalcEigenvalues(true);
    solver.solveDetailed( rhsForSchur, lambda, errorHistory );
    real_t conditionNumber = solver.getConditionNumber();
    gsMatrix<real_t> eigenvalues;
    solver.getEigenvalues(eigenvalues);
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

    gsInfo << "Estimated condition number: " << conditionNumber << "\n";
    gsInfo << "Eigenvalues: " << eigenvalues.transpose() << "\n";

    /********************** Output **************************/
    if (!out.empty())
    {
        const bool exists = gsFileManager::fileExists(out);
        std::ofstream outfile(out, std::ios_base::app);
        if (!exists)
        {
            outfile << "biharmonic_ieti_example\t"
                    << "geometry" << "\t"
                    << "patchSplitsX" << "\t"
                    << "patchSplitsY" << "\t"
                    << "rhs" << "\t"
                    << "bdyConds" << "\t"
                    << "degree" << "\t"
                    << "mutiplicity" << "\t"
                    << "refinements" << "\t"
                    << "conditionNumber" << "\t"
                    << "iter" << "\n";
         }
        outfile << "biharmonic_ieti_example\t"
                << geometry << "\t"
                << patchSplitsX << "\t"
                << patchSplitsY << "\t"
                << rhs << "\t"
                << bdyConds << "\t"
                << degree << "\t"
                << multiplicity << "\t"
                << refinements << "\t"
                << conditionNumber << "\t"
                << iter << "\n";

        gsInfo << "Write solution data to file " << out << "\n";
    }

    if (plot)
    {
        gsInfo << "Write Paraview data to file ieti_result.pvd\n";
        gsMultiPatch<> mpsol;
        for (index_t k=0; k<nPatches; ++k)
        {
            gsMatrix<> localSolutionStdBasis;
            makeSparseLUSolver(localBasisTransforms[k])->apply(localSol[k], localSolutionStdBasis);
            mpsol.addPatch( mb[k].makeGeometry(localSolutionStdBasis) );
        }
        gsWriteParaview<>(gsField<>( mp, mpsol ), "ieti_result", 1000);
    }
    if (!plot&&out.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --out to write solution to xml file.\n";
    }

    return EXIT_SUCCESS;
}



/// Helper functions

gsMultiPatch<>::uPtr tryGetRectangularGeometry(const char *s)
{
    index_t nSplitsX=0, nSplitsY=0, geoExample=0;

    if (*(s++)!='r'||*(s++)!='_') return nullptr;
    for (char c; c=*s, '0'<=c && c<='9'; ++s)
        nSplitsX = 10*nSplitsX+(c-'0');

    if (*s!='_') return nullptr;
    ++s;

    for (char c; c=*s, '0'<=c && c<='9'; ++s)
        nSplitsY = 10*nSplitsY+(c-'0');

    if (*s=='_')
    {
        ++s;
        for (char c; c=*s, '0'<=c && c<='9'; ++s)
            geoExample = 10*geoExample+(c-'0');
    }
    if (*s!='\0') return nullptr;
    if (nSplitsX == 0 || nSplitsY == 0) return nullptr;

    gsMultiPatch<>::uPtr mp = memory::make_unique(new gsMultiPatch<>());

    if (geoExample==0)
        mp->addPatch(gsNurbsCreator<>::BSplineRectangle(0,0,1,1));
    else if (geoExample==1)
        mp->addPatch(gsNurbsCreator<>::BSplineRectangle(0,-1,1,0,90));
    else if (geoExample==2)
        mp->addPatch(gsNurbsCreator<>::BSplineRectangle(-1,-1,0,0,180));
    else if (geoExample==3)
        mp->addPatch(gsNurbsCreator<>::BSplineRectangle(-1,0,0,1,270));
    else if (geoExample==4)
        mp->addPatch(gsNurbsCreator<>::BSplineRectangle(1,0,0,1));
    else if (geoExample==5)
        mp->addPatch(gsNurbsCreator<>::BSplineRectangle(1,1,0,0));
    else
        return nullptr;

    /// end
    for (index_t i=0; i<nSplitsX; ++i)
        for (index_t j=0; j<nSplitsY; ++j)
            if (i+j)
                mp->addPatch(gsNurbsCreator<>::BSplineRectangle(i,j,i+1,j+1));
    mp->computeTopology();
    gsInfo << "Initial geometry is grid of " << nSplitsX << "x" << nSplitsY << " unit squares.\n";
    return mp;
}

std::pair< gsMatrix<>, gsMatrix<> >
sampleBoundary(const gsGeometry<>& geo, boxSide s, bool inverted, index_t numberSamples)
{
    const short_t dim = geo.geoDim();
    GISMO_ASSERT (dim==2, "Only available for two dimensions.");
    const short_t dir = (s.index()-1)/2;
    const short_t prm = (s.index()-1)%2;
    const gsMatrix<>& parameterRange = geo.parameterRange();
    gsMatrix<> sample(dim,math::ipow(numberSamples,dim-1));

    for (index_t i=0; i<numberSamples; ++i)
    {
        index_t j = inverted ? (numberSamples-1-i) : i;
        sample(  dir, j) = parameterRange(  dir,0) + prm                     * (parameterRange(  dir,1)-parameterRange(  dir,0));
        sample(1-dir, j) = parameterRange(1-dir,0) + i/(numberSamples - 1.0) * (parameterRange(1-dir,1)-parameterRange(1-dir,0));
    }

    gsMatrix<> selector(2,4);
    selector.setZero();
    selector(0,dir)   = prm?1.0:-1.0;
    selector(1,dir+2) = prm?1.0:-1.0;

    return std::pair< gsMatrix<>, gsMatrix<> >( geo.eval(sample), selector*geo.deriv(sample) );
}

// Returns indices of edges, edges with 1 offset, corner dofs (4 per corner)
std::vector<index_t>
setupC1interfaceDofsVector(const gsBasis<>& basis, const gsDofMapper& dm, index_t type)
{
    GISMO_ASSERT (type>=0 && type<3, "Unknown type");
    gsMatrix<index_t> count;
    count.setZero(basis.size(),2);
    for (boxSide s=boxSide::getFirst(basis.dim()); s<boxSide::getEnd(basis.dim()); ++s)
        for (index_t offset = 0; offset < 2; ++offset)
        {
            gsVector<index_t> bdyDofs = basis.boundaryOffset(s,offset);
            for (index_t i=0; i<bdyDofs.rows(); ++i)
                count(bdyDofs[i],offset) += 1;
        }

    std::vector<index_t> result;
    for (index_t i=0; i<count.rows(); ++i)
        if (dm.is_free(i, 0))
        {
            if (type==0 && count(i,0)>0)
                result.push_back(dm.index(i, 0));
            if (type==1 && count(i,0)==0 && count(i,1)>0)
                result.push_back(dm.index(i, 0));
            if (type==2 && count(i,0)+count(i,1)>1)
                result.push_back(dm.index(i, 0));
        }

    return result;
}

real_t
getKnotSpan(const gsKnotVector<>& kv, index_t side, index_t p)
{
    GISMO_ASSERT (side==0||side==1, "Unknown side given");
    GISMO_ASSERT (kv.isOpen(), "Not a p-open knot vector");
    if (side==0)
        return kv[p+1]-kv[0];
    else
        return kv[kv.size()-1]-kv[kv.size()-p-2];
}

gsSparseMatrix<>
setupC1basisTransformation(const gsBasis<>& basis, const gsDofMapper& dm)
{
    // This function creates a sparse matix that changes the sign of the
    // the n-1 st row and column (1D). This is tensorized. This function's
    // result is put into applyDofMapperTwoSided to respect the boundary
    // conditions.
    const index_t d = basis.dim();
    gsSparseMatrix<> sm;
    for (index_t i=0; i<d; ++i)
    {
        const gsBSplineBasis<>* bsp = dynamic_cast<const gsBSplineBasis<>*>(&basis.component(d-1-i));
        GISMO_ENSURE (bsp, "Not a gsBSplineBasis given.");


        const short_t p = bsp->degree(0);
        const index_t ndofs1D = bsp->size();
        GISMO_ASSERT( ndofs1D>3, "" );

        gsSparseMatrix<> transformer1D(ndofs1D,ndofs1D);
        transformer1D.setIdentity();

        transformer1D(0,0)=1;
        transformer1D(0,1)=1;

        transformer1D(1,1)=-getKnotSpan(bsp->knots(0),0,p) / p;

        transformer1D(ndofs1D-2,ndofs1D-2)=getKnotSpan(bsp->knots(0),1,p) / p;

        transformer1D(ndofs1D-1,ndofs1D-2)=1;
        transformer1D(ndofs1D-1,ndofs1D-1)=1;

        if (i==0)
            sm = give(transformer1D);
        else
            sm = sm.kron(transformer1D);
    }

    gsSparseEntries<> se;
    se.reserve(sm.nonZeros());
    for ( index_t i=0; i<sm.outerSize(); ++i )
        for ( gsSparseMatrix<>::InnerIterator it(sm,i); it; ++it )
            if (dm.is_free(it.row(), 0) && dm.is_free(it.col(), 0))
                se.add( dm.index(it.row(), 0), dm.index(it.col(), 0), it.value() );

    gsSparseMatrix<> result(dm.freeSize(), dm.freeSize());
    result.setFrom(se);
    return result;
}


gsDofMapper
setupC1DofMapper(const gsMultiPatch<>& mp, const gsMultiBasis<>& mb, unsigned bdyConds)
{
    // Construct the object
    gsVector<index_t> patchDofSizes(mp.nPatches());
    for (size_t k=0; k<mp.nPatches(); ++k)
        patchDofSizes[k] = mb[k].size();
    gsDofMapper dm(patchDofSizes);

    // Iterate over all interfaces
    for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const index_t k1 = it->first().patch;
        const index_t k2 = it->second().patch;

        // Is the edge parameterized with the same direction (increasing/decreasing) for both patches?
        gsVector<index_t> cm;
        it->cornerMap(cm);
        GISMO_ASSERT (cm.size()==2 && cm[1] == 1-cm[0] && (cm[0] == 0 || cm[0] == 1), "Unexcpected corner map:\n"<<cm);
        const bool inverted = cm[0];

        // Check if the geometry is really C0 and C1
        std::pair< gsMatrix<>, gsMatrix<> > data1 = sampleBoundary(mp[k1], it->first(),  inverted, 11);
        std::pair< gsMatrix<>, gsMatrix<> > data2 = sampleBoundary(mp[k2], it->second(), false,    11);
        GISMO_ENSURE ((data1.first-data2.first).cwiseAbs().maxCoeff()<1e-6, "No C0 geometry between patches " << k1 << " and " << k2);
        GISMO_ENSURE ((data1.second+data2.second).cwiseAbs().maxCoeff()<1e-6, "No C1 geometry between patches " << k1 << " and " << k2);

        // Construct the dofmapper
        gsVector<index_t> s1 = mb.basis(k1).boundary(it->first().side());
        gsVector<index_t> s2 = mb.basis(k2).boundary(it->second().side());
        gsVector<index_t> s1o = mb.basis(k1).boundaryOffset(it->first().side(),1);
        gsVector<index_t> s2o = mb.basis(k2).boundaryOffset(it->second().side(),1);

        GISMO_ASSERT (s1.rows() == s2.rows() && s1.rows() == s1o.rows() && s2.rows() == s2o.rows(), "Bases dimensions not agree.");

        for (index_t i=0;i<s1.rows();++i)
        {
            const index_t j = inverted?(s1.rows()-1-i):i;
            dm.matchDof(k1,s1[i],k2,s2[j]);
            dm.matchDof(k1,s1o[i],k2,s2o[j]);
        }
    }

    // Set boundary conditions:
    // bdyConds == 0: no essential boundary conditions
    // bdyConds == 1: second biharmonic problem: essential boundary conditions only for the function values (one layer)
    // bdyConds == 2: first biharmonic problem: essential boundary conditions for the function values and normals (two layers)

    if (bdyConds)
    {
        for (gsBoxTopology::const_biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
        {
            const index_t k = it->patchIndex();
            gsVector<index_t> s = mb.basis(k).boundary(it->side()),
                              so = mb.basis(k).boundaryOffset(it->side(),1);

            GISMO_ASSERT( s.rows() == so.rows(), "");
            for (index_t i=0;i<s.rows();++i)
            {
                dm.eliminateDof(s[i],k);
                if (bdyConds == 2)
                    dm.eliminateDof(so[i],k);
            }
        }
    }

    dm.finalize();

    return dm;
}
