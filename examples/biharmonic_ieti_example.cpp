/** @file biharmonic_ieti_example.cpp

    @brief Biharmonic example for C1 smooth multipatch domains

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

// TODO: This requires the fix in boxComponent

#include <algorithm> // std::sort

#include <gismo.h>


using namespace gismo;

// Helper functions, defined at the end of the file
gsMultiPatch<>::uPtr tryGetRectangularGeometry(const char *s);
gsDofMapper setupC1DofMapper(const gsMultiPatch<>& mp, const gsMultiBasis<>& mb, index_t problemType);
std::vector<index_t> setupC1interfaceDofsVector(const gsBasis<>& basis, const gsDofMapper& dm, index_t type);
std::vector<gsSparseMatrix<>> setupC1basisTransformation(const gsMultiPatch<>& mp, const gsMultiBasis<>& mb);
gsSparseMatrix<> restrictToFreeDofs(const gsSparseMatrix<>& sm, const gsDofMapper& dm);

gsMatrix<> evalField(const gsField<>& field, const gsMatrix<>& coordinates);

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    std::string geometry("domain2d/c1_spline_lamella.xml");
    index_t patchSplitsX = 2;
    index_t patchSplitsY = 2;
    std::string rhs("1/8*pi^4*sin(pi*x/2)*sin(pi*y/2)");
    index_t problemType = 1;
    index_t degree = 2;
    index_t multiplicity = 1;
    index_t refinements = 2;
    bool standardDirichletPreconder = false;
    real_t tolerance = 1.e-6;
    index_t maxIterations = 100;
    std::string out;
    std::string mat;
    std::string log;
    bool plot = false;

    gsCmdLine cmd("Biharmonic example for C1 smooth multipatch domains.");

    cmd.addString("g", "Geometry",              "Chosen geometry file", geometry);
    cmd.addInt   ("x", "PatchSplitsX",          "Number of splits (coordinate direction x)", patchSplitsX);
    cmd.addInt   ("y", "PatchSplitsY",          "Number of splits (coordinate direction y)", patchSplitsY);
    cmd.addString("f", "Rhs",                   "Chosen right-hand function", rhs);
    cmd.addInt   ("t", "ProblemType",           "Boundary conds: (1) first (u, \u2202nu) or (2) second (u, \u0394u) biharmonic problem", problemType);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addInt   ("m", "Multiplicity",          "Multiplicity of knots for B-spline discretization space", multiplicity);
    cmd.addInt   ("r", "Refinements",           "Number of uniform h-refinement steps to perform before solving", refinements);
    cmd.addSwitch(     "StandardDirichletPreconder", "Do not use component wise Dirichlet preconder", standardDirichletPreconder);
    cmd.addReal  ("",  "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Stopping criterion for linear solver", maxIterations);
    cmd.addString("",  "out",                   "Write solution to file", out);
    cmd.addString("",  "mat",                   "Write matrix to file", mat);
    cmd.addString("",  "log",                   "Write solution data to log file", log);
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

    if (problemType < 1 || problemType > 2)
    {
        gsInfo << "Invalid choice for --Type (-t).\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Setup dofMapper and local basis transformations... " << std::flush;
    gsDofMapper dm = setupC1DofMapper(mp, mb, problemType);
    std::vector<gsSparseMatrix<>> localBasisTransformsFullBasis = setupC1basisTransformation(mp, mb);
    std::vector<gsSparseMatrix<>> localBasisTransforms;
    localBasisTransforms.reserve(nPatches);
    gsInfo << "done:\n" << dm << "\n";

    /********************** Define rhs **********************/

    gsInfo << "Rhs function is " << rhs << "\n";
    gsFunctionExpr<> f(rhs,mp.geoDim());

    /****************** Setup ietimapper ********************/

    gsInfo << "Setup ietimapper... " << std::flush;
    gsMatrix<> fixedPart;
    fixedPart.setZero(dm.boundarySize(),1);
    gsIetiMapper<> ietiMapper(mb,dm,fixedPart);

    std::vector<index_t> cornerDofsGlobal;
    cornerDofsGlobal.reserve(nPatches*16);
    for (index_t k=0; k<nPatches; ++k)
    {
        std::vector<index_t> cornerDofs = setupC1interfaceDofsVector(mb[k], ietiMapper.dofMapperLocal(k), 2);
        for (size_t i=0; i<cornerDofs.size(); ++i)
        {
            cornerDofsGlobal.push_back(dm.index(cornerDofs[i],k));
            ietiMapper.declareDofAsPrimal(k, cornerDofs[i], true);
        }
    }

    ietiMapper.computeJumpMatrices(/*fullyRedundant*/false, /*excludeCorners*/false, /*exlude*/cornerDofsGlobal);

    gsIetiSystem<> ieti;
    ieti.reserve(nPatches+1);

    gsScaledDirichletPrec<> prec;
    prec.reserve(nPatches);

    gsPrimalSystem<> primal(ietiMapper.nPrimalDofs());
    primal.setEliminatePointwiseConstraints(true);

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

        gsSparseMatrix<> localBasisTransform = restrictToFreeDofs(localBasisTransformsFullBasis[k], ietiMapper.dofMapperLocal(k));
        localBasisTransforms.push_back(localBasisTransform);

        // Fetch data
        gsSparseMatrix<real_t, RowMajor> jumpMatrix  = ietiMapper.jumpMatrix(k);
        gsSparseMatrix<>                 localMatrix = localBasisTransform*A.matrix()*localBasisTransform.transpose();
        gsMatrix<>                       localRhs    = localBasisTransform*A.rhs();

        GISMO_ASSERT(jumpMatrix.cols() == localMatrix.rows(), "");

        if (standardDirichletPreconder)
        {
            gsSparseMatrix<> localBasisTransform2 = localBasisTransform.cwiseAbs();
            gsSparseMatrix<> localMatrix2 = localBasisTransform2*A.matrix()*localBasisTransform2.transpose();
            std::vector<index_t> vec0 = setupC1interfaceDofsVector(mb_local[0], ietiMapper.dofMapperLocal(k), 0);
            std::vector<index_t> vec1 = setupC1interfaceDofsVector(mb_local[0], ietiMapper.dofMapperLocal(k), 1);
            vec0.reserve(vec0.size()+vec1.size());
            vec0.insert(vec0.end(), vec1.begin(), vec1.end());
            std::sort(vec0.begin(), vec0.end());
            vec0.erase(std::unique(vec0.begin(), vec0.end()), vec0.end());
            prec.addSubdomain(
                gsScaledDirichletPrec<>::restrictToSkeleton(
                    jumpMatrix,
                    localMatrix2,
                    vec0
                )
            );
        }
        else
        {
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

    // This is the main cg iteration
    //! [Solve]
    gsConjugateGradient<> solver( ieti.schurComplement(), preconder );
    solver.setOptions( cmd.getGroup("Solver") );
    solver.setCalcEigenvalues(true);

    gsMatrix<> errorHistory;
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

    // transform solution back to standard basis
    for (index_t k=0; k<nPatches; ++k)
        localSol[k] = localBasisTransforms[k].transpose() * localSol[k];

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
    if (!log.empty())
    {
        const bool exists = gsFileManager::fileExists(log);
        std::ofstream outfile(log, std::ios_base::app);
        if (!exists)
        {
            outfile << "biharmonic_ieti_example\t"
                    << "Geometry" << "\t"
                    << "PatchSplitsX" << "\t"
                    << "PatchSplitsY" << "\t"
                    << "Rhs" << "\t"
                    << "ProblemType" << "\t"
                    << "Degree" << "\t"
                    << "Mutiplicity" << "\t"
                    << "Refinements" << "\t"
                    << "ConditionNumber" << "\t"
                    << "Iter" << "\n";
         }
        outfile << "biharmonic_ieti_example\t"
                << geometry << "\t"
                << patchSplitsX << "\t"
                << patchSplitsY << "\t"
                << rhs << "\t"
                << problemType << "\t"
                << degree << "\t"
                << multiplicity << "\t"
                << refinements << "\t"
                << conditionNumber << "\t"
                << iter << "\n";

        gsInfo << "Write solution data to file " << log << "\n";
    }

    if (plot)
    {
        gsInfo << "Write Paraview data to file ieti_result.pvd\n";
        gsMultiPatch<> mpsol;
        for (index_t k=0; k<nPatches; ++k)
            mpsol.addPatch(mb[k].makeGeometry(ietiMapper.incorporateFixedPart(k,localSol[k])));
        gsWriteParaview<>(gsField<>(mp, mpsol), "ieti_result", 1000);

    }

    if (!out.empty())
    {
        gsInfo << "Write solution to file " << out << "\n";
        gsMatrix<> coordinates;
        coordinates.setZero(2,2);
        coordinates(0,1)=1;
        coordinates(1,1)=1;
        gsMultiPatch<> mpsol;
        for (index_t k=0; k<nPatches; ++k)
            mpsol.addPatch(mb[k].makeGeometry(ietiMapper.incorporateFixedPart(k,localSol[k])));
        gsMatrix<> result = evalField(gsField<>(mp, mpsol), coordinates);
        std::ofstream outfile(out, std::ios_base::app);
        //outfile << result << "\n\n";
        outfile << "(* created with G+smo *)\ndata = {\n";
        for (index_t i=0; i<result.rows(); ++i)
        {
            outfile << "{";
            for (index_t j=0; j<result.cols(); ++j)
                outfile << std::fixed << std::showpoint << std::setprecision(10) << result(i,j) << ((j<result.cols()-1) ? ", " : "");
            outfile << ((i<result.rows()-1) ? "},\n" : "}");
        }
        outfile << "};\n\n";
    }

    if (!mat.empty())
    {
        gsInfo << "Write mat to file " << mat << "\n";

        gsMatrix<> matrix;
        ieti.saddlePointProblem()->toMatrix(matrix);

        std::ofstream outfile(mat, std::ios_base::app);
        //outfile << matrix << "\n\n";
        outfile << "(* created with G+smo *)\nmat = {\n";
        for (index_t i=0; i<matrix.rows(); ++i)
        {
            outfile << "{";
            for (index_t j=0; j<matrix.cols(); ++j)
                outfile << std::fixed << std::showpoint << std::setprecision(10) << matrix(i,j) << ((j<matrix.cols()-1) ? ", " : "");
            outfile << ((i<matrix.rows()-1) ? "},\n" : "}");
        }
        outfile << "};\n\n";
    }

    if (!plot&&out.empty()&&mat.empty()&&log.empty())
    {
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution or --out to write solution data to xml file.\n";
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
            if (type==0 && (count(i,0)==1 || count(i,0)+count(i,1)>1))
                result.push_back(dm.index(i, 0));
            if (type==1 && (count(i,1)==1 || count(i,0)+count(i,1)>1))
                result.push_back(dm.index(i, 0));
            if (type==2 && count(i,0)+count(i,1)>1)
                result.push_back(dm.index(i, 0));
        }

    return result;
}


gsSparseMatrix<>
setupC1localBasisTransformation(const gsBasis<>& basis, const gsMatrix<>& scaling)
{
    // This function creates a sparse matix that changes the sign of the
    // the n-1 st row and column (1D). This is tensorized. This function's
    // result is put into applyDofMapperTwoSided to respect the boundary
    // conditions.
    const index_t d = basis.dim();
    gsSparseMatrix<> transformation;

    for (index_t i=0; i<d; ++i)
    {
        const gsBSplineBasis<>* bsp = dynamic_cast<const gsBSplineBasis<>*>(&basis.component(d-1-i));
        GISMO_ENSURE (bsp, "Not a gsBSplineBasis given.");

        const short_t p = bsp->degree(0);
        const index_t ndofs1D = bsp->size();
        const gsKnotVector<>& kv = bsp->knots(0);
        const real_t h0 = kv[p+1]-kv[0], hN = kv[kv.size()-1]-kv[kv.size()-p-2];

        GISMO_ENSURE( ndofs1D>3, "Not enough knots to perform basis transformation" );
        GISMO_ENSURE( kv.isOpen(), "Not a p-open knot vector" );

        gsSparseMatrix<> transformation1D(ndofs1D,ndofs1D);
        transformation1D.setIdentity();

        //transformation1D(0,0) = 1;
        transformation1D(0,1) = 1;
        transformation1D(1,1) = h0 / p / scaling(0,2*i);

        transformation1D(ndofs1D-2,ndofs1D-2) = hN / p / scaling(0,2*i+1);
        transformation1D(ndofs1D-1,ndofs1D-2) = 1;
        //transformation1D(ndofs1D-1,ndofs1D-1) = 1;

        if (i==0)
            transformation = give(transformation1D);
        else
            transformation = transformation.kron(transformation1D);
    }
    return transformation;
}

index_t
getCorner(boxSide s, bool parameter)
{
    // parameter == 0  begin
    // parameter == 1  end
    const index_t par_x = (s.direction()==0) ? s.parameter() : parameter;
    const index_t par_y = (s.direction()==1) ? s.parameter() : parameter;
    return par_x + 2*par_y + 1;
}


void
updateBasisTransformationSigns(const gsMultiPatch<>& mp, const gsMultiBasis<>& mb, const gsMatrix<index_t>& cornerRanks, std::vector<gsSparseMatrix<>>& transformations)
{
    for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const index_t k1 = it->first().patch;
        const index_t k2 = it->second().patch;

        const gsVector<index_t> s1 = mb.basis(k1).boundaryOffset(it->first().side(),1);
        const gsVector<index_t> s2 = mb.basis(k2).boundaryOffset(it->second().side(),1);
        GISMO_ASSERT(s1.rows()==s2.rows(), "Bases do not match.");

        gsVector<index_t> cm;
        it->cornerMap(cm);
        GISMO_ASSERT (cm.size()==2 && cm[1] == 1-cm[0] && (cm[0] == 0 || cm[0] == 1), "Unexcpected corner map:\n"<<cm);
        const bool inverted = cm[0];

        const index_t sz = s1.rows();
        for (index_t i=0; i<sz; ++i)
        {
            index_t diff = k1 - k2;

            if (i<2)
                diff = cornerRanks(k1,getCorner(it->first().side(),0)-1)
                     - cornerRanks(k2,getCorner(it->second().side(),inverted)-1);
            if (i>=sz-2)
                diff = cornerRanks(k1,getCorner(it->first().side(),1)-1)
                     - cornerRanks(k2,getCorner(it->second().side(),!inverted)-1);
            if (diff==0)
                diff = k1 - k2;

            const index_t j = inverted ? sz-1-i : i;

            if (diff>0)
                transformations[k1].row(s1[i]) *= -1;
            else
                transformations[k2].row(s2[j]) *= -1;
        }

    }

}


gsMatrix<index_t>
cornerRanker(const gsMultiPatch<>& mp)
{
    std::vector<std::vector<patchComponent>> comps = mp.allComponents();

    gsMatrix<index_t> assignedPatchPerCorner;
    assignedPatchPerCorner.setZero(mp.nPatches(), 4);

    for (size_t i=0; i<comps.size(); ++i)
    {
        // We are only interested in corners in the interior
        if (comps[i].size()!=4||comps[i][0].dim()!=0)
            continue;
        index_t kmin = comps[i][0].patch();
        for (size_t j=0; j<comps[i].size(); ++j)
            if (comps[i][j].patch()<kmin)
                kmin = comps[i][j].patch();
        for (size_t j=0; j<comps[i].size(); ++j)
        {
            const index_t k = comps[i][j].patch();
            const index_t idx = comps[i][j].asCorner().m_index-1;
            assignedPatchPerCorner(k,idx) = kmin;
        }

    }

    gsMatrix<index_t> cornerRanks;
    cornerRanks.setConstant(mp.nPatches(), 4, -1);
    for (size_t i=0; i<comps.size(); ++i)
    {
        // We are only interested in edges between two patches
        if (comps[i].size()!=2||comps[i][0].dim()!=1)
            continue;

        for (index_t r=0; r<2; ++r)
        {
            const index_t k0 = comps[i][r].patch();
            const index_t k1 = comps[i][1-r].patch();
            for (index_t j=0; j<4; ++j)
            {
                if (assignedPatchPerCorner(k0,j)==k0)
                    cornerRanks(k0,j) = -3;
                else if (assignedPatchPerCorner(k0,j)==k1)
                    cornerRanks(k0,j) = -2;
            }
        }

    }
    return cornerRanks;
}

gsMatrix<>
getScaling(const gsMultiPatch<>& mp)
{
    gsMatrix<> scaling;
    scaling.setOnes(mp.nPatches(), math::ipow(2,mp.geoDim()));
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
        real_t scalingFactor = data1.second.cwiseAbs().maxCoeff() / data2.second.cwiseAbs().maxCoeff();
        scaling(k2, it->second().index()-1) = math::sqrt(scalingFactor);
        scaling(k1, it->first().index()-1) = 1/math::sqrt(scalingFactor);
        GISMO_ENSURE ((data1.first-data2.first).cwiseAbs().maxCoeff()<1e-6, "No C0 geometry between patches " << k1 << " and " << k2);
        GISMO_ENSURE ((data1.second+scalingFactor*data2.second).cwiseAbs().maxCoeff()<1e-6, "No C1 geometry between patches " << k1 << " and " << k2);

    }

    return scaling;
}


std::vector<gsSparseMatrix<>>
setupC1basisTransformation(const gsMultiPatch<>& mp, const gsMultiBasis<>& mb)
{

    std::vector<gsSparseMatrix<>> transformations;
    transformations.reserve(mp.nPatches());

    gsMatrix<> scaling = getScaling(mp);

    for (size_t k=0; k<mp.nPatches(); ++k)
        transformations.push_back(setupC1localBasisTransformation(mb[k], scaling.row(k)));

    updateBasisTransformationSigns(mp, mb, cornerRanker(mp), transformations);

    return transformations;
}


gsSparseMatrix<>
restrictToFreeDofs(const gsSparseMatrix<>& sm, const gsDofMapper& dm)
{
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
setupC1DofMapper(const gsMultiPatch<>& mp, const gsMultiBasis<>& mb, index_t problemType)
{
    // Construct the object
    gsVector<index_t> patchDofSizes(mp.nPatches());
    for (size_t k=0; k<mp.nPatches(); ++k)
        patchDofSizes[k] = mb[k].size();
    gsDofMapper dm(patchDofSizes);

    // Set interface conditions
    for (gsBoxTopology::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        const index_t k1 = it->first().patch;
        const index_t k2 = it->second().patch;

        // Is the edge parameterized with the same direction (increasing/decreasing) for both patches?
        gsVector<index_t> cm;
        it->cornerMap(cm);
        GISMO_ASSERT(cm.size()==2 && cm[1] == 1-cm[0] && (cm[0] == 0 || cm[0] == 1), "Unexcpected corner map.");
        const bool inverted = cm[0];

        // Register conditions
        for (index_t offset=0; offset<2; ++offset)
        {
            gsVector<index_t> s1 = mb.basis(k1).boundaryOffset(it->first().side(), offset);
            gsVector<index_t> s2 = mb.basis(k2).boundaryOffset(it->second().side(), offset);

            GISMO_ASSERT(s1.rows() == s2.rows(), "Bases dimensions not agree.");

            for (index_t i=0;i<s1.rows();++i)
            {
                const index_t j = inverted?(s1.rows()-1-i):i;
                dm.matchDof(k1,s1[i],k2,s2[j]);
            }
        }
    }

    // Set boundary conditions
    for (gsBoxTopology::const_biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
    {
        const index_t k = it->patchIndex();
        // For first biharmonic (problemType==1), set two layers for boundary conditons
        // For first biharmonic (problemType==2), set only layer for boundary conditons
        for (index_t offset=0; offset<(problemType==1?2:1); ++offset)
        {
            gsVector<index_t> s = mb.basis(k).boundaryOffset(it->side(), offset);
            for (index_t i=0;i<s.rows();++i)
                dm.eliminateDof(s[i],k);
        }
    }

    dm.finalize();

    return dm;
}



gsMatrix<>
evalField(const gsField<>& field, const gsMatrix<>& coordinates)
{
    GISMO_ASSERT (coordinates.cols() == 2, "Assume minima and maxima");
    GISMO_ASSERT (coordinates.rows() == 2, "Only 2d");
    const index_t n = 100;
    gsMatrix<> result; result.setZero(n,n);
    gsMatrix<> mask; mask.setZero(n,n);

    gsMatrix<> par_point(2,1), ph_point(2,1), value(1,1);
    for (index_t k=0; k<field.nPieces(); ++k)
    {
        for (index_t i=0; i<n; ++i)
            for (index_t j=0; j<n; ++j)
            {
                par_point(0,0) = (real_t)i / (n-1);
                par_point(1,0) = (real_t)j / (n-1);
                ph_point = field.point(par_point,k);
                const real_t a = n * (ph_point(0,0) - coordinates(0,0)) / (coordinates(0,1)-coordinates(0,0));
                const real_t b = n * (ph_point(1,0) - coordinates(1,0)) / (coordinates(1,1)-coordinates(1,0));

                const index_t aa = a;
                const index_t bb = b;

                const real_t value = field.value(par_point,k)(0,0);

                for (index_t o1 = -1; o1<3; ++o1)
                    for (index_t o2 = -1; o2<3; ++o2)
                    {
                        const index_t a0 = aa + o1;
                        const index_t b0 = bb + o2;
                        if (a0>0 && a0<n && b0>0 && b0<n)
                        {
                            const real_t distance = (a-a0)*(a-a0) + (b-b0)*(b-b0);
                            if (mask(a0,b0)==0 || mask(a0,b0)>distance)
                            {
                                mask(a0,b0) = distance;
                                result(a0,b0) = value;
                            }
                        }
                    }

            }
    }

    return result;
}
