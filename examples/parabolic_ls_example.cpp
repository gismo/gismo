/** @file parabolic_ls_example.cpp

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

//#define printMat(m) gsInfo << #m << " (" << m.rows() << "x" << m.cols() << "):\n\n"; for (int i=0; i<m.rows(); ++i) {for (int j=0; j<m.cols(); ++j) gsInfo << m(i,j) << "\t"; gsInfo << "\n";} gsInfo << "\n\n";
#define printMat(m)

template<typename S>
gsLinearOperator<>::Ptr fastDiagnonalization(const gsSparseMatrix<>& A1, const gsSparseMatrix<>& B1, const gsSparseMatrix<>& A2, const gsSparseMatrix<>& B2, const S& makeSolver)
{
    GISMO_ASSERT(A1.rows() == A1.cols() && A1.rows() == A2.rows() && A2.rows() == A2.cols(), "");
    GISMO_ASSERT(B1.rows() == B1.cols() && B1.rows() == B2.rows() && B2.rows() == B2.cols(), "");

    typedef gsMatrix<>::GenSelfAdjEigenSolver EVSolver;
    EVSolver ges;
    ges.compute(A2, A1, gsEigen::ComputeEigenvectors);
    gsSparseMatrix<> D1(A1.rows(), A1.rows()), D2(A1.rows(), A1.rows());
    GISMO_ASSERT (ges.eigenvalues().rows() == A1.rows() && ges.eigenvalues().cols() == 1, "");
    GISMO_ASSERT (ges.eigenvectors().rows() == A1.rows() && ges.eigenvectors().cols() == A1.rows(), "");
    for (index_t i=0; i<A1.rows(); ++i)
    {
        D1(i,i) = 1;
        D2(i,i) = ges.eigenvalues()(i,0);
    }
    gsSparseMatrix<> system = D1.kron(B1)+D2.kron(B2); //TODO: is this orderd such that a direct solver would do this efficiently?
    gsLinearOperator<>::Ptr systemSolver = makeSolver(system);
    gsMatrix<> eigs = ges.eigenvectors();
    gsMatrix<> eigsT = ges.eigenvectors().transpose();

    gsLinearOperator<>::Ptr transform  = gsKroneckerOp<>::make( makeMatrixOp(eigs.moveToPtr()), gsIdentityOp<>::make(B1.rows()) );
    gsLinearOperator<>::Ptr transformT = gsKroneckerOp<>::make( makeMatrixOp(eigsT.moveToPtr()), gsIdentityOp<>::make(B1.rows()) );
    return gsProductOp<>::make( transformT, systemSolver, transform );

}

gsLinearOperator<>::Ptr mkSparseLUSolver(const gsSparseMatrix<>& m) { return makeSparseLUSolver(m); }


int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    index_t geoIdx = 1;
    index_t refinementsX = 2;
    index_t refinementsT = 2;
    index_t degree = 1;
    real_t kappa = 1.;
    real_t h0 = 1.;
    real_t tau0 = 1.;
    index_t maxIterations = 100;
    real_t tolerance = 1.e-6;
    index_t preSmooth = 1;
    index_t postSmooth = 1;
    index_t cycles = 1;
    index_t exactPreconder = 1;
    index_t fdPreconder = 1;
    index_t mgPreconder = 0;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("parabolic_ls_example");
    cmd.addInt   ("g", "Geometry",              "0=Rectangle, 1=Quarter Annulus", geoIdx);
    cmd.addInt   ("r", "RefinementsX",          "Number of uniform h-refinement steps to perform before solving", refinementsX);
    cmd.addInt   ("s", "RefinementsT",          "Number of uniform tau-refinement steps to perform before solving", refinementsT);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addReal  ("k", "Kappa",                 "Diffusion parameter", kappa);
    cmd.addReal  ("",  "HNull",                 "h0 for multilevel", h0);
    cmd.addReal  ("",  "TauNull",               "tau0 for multilevel", tau0);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Maximum iterations for linear solver", maxIterations);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "MG.NumPreSmooth",       "Number of pre smoothing steps (only for mg)", preSmooth);
    cmd.addInt   ("",  "MG.NumPostSmooth",      "Number of post smoothing steps (only for mg)", postSmooth);
    cmd.addInt   ("c", "MG.NumCycles",          "Number of multi-grid cycles for coarse-grid correction, i.e., 1=V, 2=W cycle", cycles);
    cmd.addInt   ("",  "ExactPreconder",        "Use that scheme", exactPreconder);
    cmd.addInt   ("",  "FdPreconder",           "Use that scheme", fdPreconder);
    cmd.addInt   ("",  "MgPreconder",           "Use that scheme", mgPreconder);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    bool ok = true;

    gsMultiPatch<> geos[] = {
        *gsNurbsCreator<>::BSplineRectangle(),
        *gsNurbsCreator<>::BSplineQuarterAnnulus()
    };

    if (geoIdx<0          || geoIdx>=(int)(util::size(geos))         ) { gsInfo << "Unfeasible choice for --Geometry (-g).\n";        ok=false; }
    if (!ok) return -1;

    gsInfo << "Run parabolic_oc_example with options:\n" << cmd << std::endl;

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    const gsMultiPatch<>& mp = geos[geoIdx];
    gsMultiPatch<> tp(*gsNurbsCreator<>::BSplineUnitInterval(2));

    gsInfo << "done.\n";

    /************** Define boundary and initial conditions **************/

    gsInfo << "Define boundary and initial conditions... " << std::flush;

    gsConstantFunction<> zero( 0, mp.geoDim() );

    gsFunctionExpr<> y0( "if(((x-3/2*cos(3*pi/8))^2+(y-3/2*sin(3*pi/8))^2)<0.04,1,0) + "
                         "if(((x-3/2*cos(2*pi/8))^2+(y-3/2*sin(2*pi/8))^2)<0.04,1,0) + "
                         "if(((x-3/2*cos(1*pi/8))^2+(y-3/2*sin(1*pi/8))^2)<0.04,1,0) ", mp.geoDim() );

    gsBoundaryConditions<> bc;
    for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
    {
        bc.addCondition( *it, condition_type::dirichlet, zero );
    }

    gsBoundaryConditions<> ic; // Handle ic separately. Here: none

    gsInfo << "done.\n";

    /************ Setup bases and adjust degree *************/
    gsMultiBasis<> mb(mp);
    gsMultiBasis<> tb1(tp);
    gsMultiBasis<> tb2(tp);

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    for ( size_t i = 0; i < mb.nBases(); ++ i )
        mb[i].setDegreePreservingMultiplicity(degree);

    for ( size_t i = 0; i < tb1.nBases(); ++ i )
        tb1[i].setDegreePreservingMultiplicity(degree);

    for ( size_t i = 0; i < tb2.nBases(); ++ i )
        tb2[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinementsX; ++i )
        mb.uniformRefine();

    for ( index_t i = 0; i < refinementsT; ++i )
        tb1.uniformRefine();

    for ( index_t i = 0; i < refinementsT; ++i )
        tb2.uniformRefine();

    for ( size_t i = 0; i < tb2.nBases(); ++ i )
        tb2[i].reduceContinuity(1);

    gsInfo << "done.\n";
    gsInfo << "Knots in tb1: " << dynamic_cast<gsTensorBSplineBasis<1,real_t>&>(tb1[0]).component(0).knots() << "\n";
    gsInfo << "Knots in tb2: " << dynamic_cast<gsTensorBSplineBasis<1,real_t>&>(tb2[0]).component(0).knots() << "\n";



    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrices... " << std::flush;


    // Assemble space matrices
    gsSparseMatrix<> space_mass;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(mb);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space       u = assembler.getSpace(mb,1,0);
        bc.setGeoMap(mp);
        u.setup(bc, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( u * u.tr() * meas(G) );
        space_mass = assembler.matrix();
    }

    gsSparseMatrix<> space_stiff;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(mb);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space       u = assembler.getSpace(mb,1,0);
        bc.setGeoMap(mp);
        u.setup(bc, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G) );
        space_stiff = assembler.matrix();
    }

    // Assemble in time

    gsSparseMatrix<> time_stiff1;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(tb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(tp);
        gsExprAssembler<>::space u = assembler.getSpace(tb1,1,0);
        ic.setGeoMap(tp);
        u.setup(ic, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G) );
        time_stiff1 = assembler.matrix();
        const index_t n = time_stiff1.rows();
        time_stiff1 = time_stiff1.block(1, 1, n-1, n-1);
    }



    gsSparseMatrix<> time_mass1;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(tb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(tp);
        gsExprAssembler<>::space u = assembler.getSpace(tb1,1,0);
        ic.setGeoMap(tp);
        u.setup(ic, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( u * u.tr() * meas(G) );
        time_mass1 = assembler.matrix();
        const index_t n = time_mass1.rows();
        time_mass1 = time_mass1.block(1, 1, n-1, n-1);
    }


    gsSparseMatrix<> time_mass2;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(tb2);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(tp);
        gsExprAssembler<>::space u = assembler.getSpace(tb2,1,0);
        ic.setGeoMap(tp);
        u.setup(ic, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( u * u.tr() * meas(G) );
        time_mass2 = assembler.matrix();
    }

    gsSparseMatrix<> time_massL;
    {
        gsExprAssembler<> assembler(2,2);
        assembler.setIntegrationElements(tb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(tp);
        gsExprAssembler<>::space u = assembler.getSpace(tb1,1,0);
        gsExprAssembler<>::space v = assembler.getSpace(tb2,1,1);
        ic.setGeoMap(tp);
        u.setup(ic, dirichlet::interpolation, 0);
        v.setup(ic, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( u * v.tr() * meas(G) );
        time_massL = assembler.matrix();

        const index_t n1 = time_mass1.rows()+1;
        const index_t n2 = time_mass2.rows();
        GISMO_ENSURE(time_massL.rows() == n1+n2, time_massL.rows()<<"=="<<n1<<"+"<<n2);
        time_massL = time_massL.block(1, n1, n1-1, n2).transpose(); //TODO
    }


    gsSparseMatrix<> time_gradL;
    {
        gsExprAssembler<> assembler(2,2);
        assembler.setIntegrationElements(tb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(tp);
        gsExprAssembler<>::space u = assembler.getSpace(tb1,1,0);
        gsExprAssembler<>::space v = assembler.getSpace(tb2,1,1);
        ic.setGeoMap(tp);
        u.setup(ic, dirichlet::interpolation, 0);
        v.setup(ic, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( igrad(u,G) * v.tr() * meas(G) );
        time_gradL = assembler.matrix();

        const index_t n1 = time_mass1.rows()+1;
        const index_t n2 = time_mass2.rows();
        GISMO_ENSURE(time_gradL.rows() == n1+n2, "");
        time_gradL = time_gradL.block(1, n1, n1-1, n2).transpose(); //TODO

    }


    if (0)
    {
        gsInfo << "\nspace_mass=\n" << space_mass << "\n";
        gsInfo << "\nspace_stiff=\n" << space_stiff << "\n";
        gsInfo << "\ntime_massL=\n" << time_massL << "\n";
        gsInfo << "\ntime_gradL=\n" << time_gradL << "\n";
        gsInfo << "\ntime_stiff1=\n" << time_stiff1 << "\n";
        gsInfo << "\ntime_mass1=\n" << time_mass1 << "\n";
        gsInfo << "\ntime_mass2=\n" << time_mass2 << "\n";
    }
    gsInfo << "done.\n";
    gsInfo << "Setup of linear operators... " << std::flush;


    gsLinearOperator<>::Ptr Lh
        = gsSumOp<>::make(
            gsKroneckerOp<>::make( makeMatrixOp(time_gradL), makeMatrixOp(space_mass) ),
            gsScaledOp<>::make( gsKroneckerOp<>::make( makeMatrixOp(time_massL), makeMatrixOp(space_stiff) ), kappa )
        );
    gsLinearOperator<>::Ptr LhT
        = gsSumOp<>::make(
            gsKroneckerOp<>::make( makeMatrixOp(time_gradL.transpose()), makeMatrixOp(space_mass) ),
            gsScaledOp<>::make( gsKroneckerOp<>::make( makeMatrixOp(time_massL.transpose()), makeMatrixOp(space_stiff) ), kappa )
        );
    gsLinearOperator<>::Ptr dualPc
        = gsKroneckerOp<>::make( makeSparseCholeskySolver(time_mass2), makeSparseCholeskySolver(space_stiff) );

    gsLinearOperator<>::Ptr leastSquares = gsProductOp<>::make( Lh, dualPc, LhT );

    gsInfo << "done; " << leastSquares->rows() << " dofs.\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Somewhat exact preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    index_t iter1;
    real_t cond1;
    gsInfo << "Setup of somewhat exact preconder... " << std::flush;
    if (exactPreconder)
    {
        const index_t primalDim = time_stiff1.rows() * space_stiff.rows();

        gsSparseMatrix<> preconderSelector(primalDim,2*primalDim);
        for (index_t i=0; i<primalDim; ++i)
            preconderSelector(i,i) = 1;

        gsSparseMatrix<> preconderMatrix(2*primalDim,2*primalDim);
        {
            gsSparseMatrix<> preconderBlock11 =  time_mass1.kron(space_stiff);
            gsSparseMatrix<> preconderBlock12 =  time_stiff1.kron(space_mass);
            gsSparseMatrix<> preconderBlock22 = -time_stiff1.kron(space_stiff);

            gsSparseEntries<> se;
            se.reserve( preconderBlock11.nonZeros() + 2*preconderBlock12.nonZeros() + preconderBlock22.nonZeros() );

            for (index_t i = 0; i < preconderBlock11.outerSize(); ++i)
                for (gsSparseMatrix<>::InnerIterator it(preconderBlock11,i); it; ++it)
                    se.add(it.row(), it.col(), it.value());

            for (index_t i = 0; i < preconderBlock12.outerSize(); ++i)
                for (gsSparseMatrix<>::InnerIterator it(preconderBlock12,i); it; ++it)
                {
                    se.add(it.row(), primalDim+it.col(), it.value());
                    se.add(primalDim+it.col(), it.row(), it.value());
                }

            for (index_t i = 0; i < preconderBlock22.outerSize(); ++i)
                for (gsSparseMatrix<>::InnerIterator it(preconderBlock22,i); it; ++it)
                    se.add(primalDim+it.row(), primalDim+it.col(), it.value());

            preconderMatrix.setFrom(se);

        }

        //gsInfo << "\npreconderMatrix=\n" << preconderMatrix << "\n";


        gsLinearOperator<>::Ptr preconder
            = gsProductOp<>::make( makeMatrixOp(preconderSelector.transpose()), makeSparseLUSolver(preconderMatrix), makeMatrixOp(preconderSelector) );


        gsInfo << "done: " << preconder->rows() << " dofs.\n";

        gsInfo << "Setup cg solver and solve... " << std::flush;

        gsMatrix<> x;
        x.setRandom( leastSquares->rows(), 1 );
        gsMatrix<> rhs;
        rhs.setRandom( leastSquares->rows(), 1 ); // TODO
        gsMatrix<> errorHistory;
        gsConjugateGradient<> solver( leastSquares, preconder );
        solver.setCalcEigenvalues(true);
        solver.setOptions( cmd.getGroup("Solver") ).solveDetailed( rhs, x, errorHistory );

        gsInfo << "done.\n\n";

        iter1 = errorHistory.rows()-1;
        const bool success = errorHistory(iter1,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter1 << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter1 << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond1 = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond1 << "\n";
    }
    else
    {
        gsInfo << "skip.\n";
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  FD preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    index_t iter2;
    real_t cond2;
    gsInfo << "Setup of FD preconder... " << std::flush;
    if (fdPreconder)
    {
        const index_t primalDim = space_stiff.rows();
        const index_t primalDimTime = time_stiff1.rows();
        gsSparseMatrix<> preconderSelector(primalDim*primalDimTime,2*primalDim*primalDimTime);
        for (index_t i=0; i<primalDimTime; ++i)
            for (index_t j=0; j<primalDim; ++j)
                preconderSelector(primalDim*i+j,2*primalDim*i+j) = 1;

        gsLinearOperator<>::Ptr preconder;
        {
            gsSparseMatrix<> preconderMatrix1(2*primalDim,2*primalDim);
            gsSparseMatrix<> preconderMatrix2(2*primalDim,2*primalDim);
            gsSparseEntries<> se1;
            gsSparseEntries<> se2;
            se1.reserve( 2*space_mass.nonZeros() + space_stiff.nonZeros() );
            se2.reserve( space_stiff.nonZeros() );

            for (index_t i = 0; i < space_stiff.outerSize(); ++i)
                for (gsSparseMatrix<>::InnerIterator it(space_stiff,i); it; ++it)
                    se2.add(it.row(), it.col(), it.value());

            for (index_t i = 0; i < space_mass.outerSize(); ++i)
                for (gsSparseMatrix<>::InnerIterator it(space_mass,i); it; ++it)
                {
                    se1.add(it.row(), primalDim+it.col(), it.value());
                    se1.add(primalDim+it.col(), it.row(), it.value());
                }

            for (index_t i = 0; i < space_stiff.outerSize(); ++i)
                for (gsSparseMatrix<>::InnerIterator it(space_stiff,i); it; ++it)
                    se1.add(primalDim+it.row(), primalDim+it.col(), -it.value());

            preconderMatrix1.setFrom(se1);
            preconderMatrix2.setFrom(se2);

            preconder = gsProductOp<>::make(
                makeMatrixOp(preconderSelector.transpose()),
                fastDiagnonalization(time_stiff1, preconderMatrix1, time_mass1, preconderMatrix2, mkSparseLUSolver),
                makeMatrixOp(preconderSelector)
            );

        }

        //gsInfo << "\npreconderMatrix=\n" << preconderMatrix << "\n";
        gsInfo << "done: " << preconder->rows() << " dofs.\n";

        gsInfo << "Setup cg solver and solve... " << std::flush;

        gsMatrix<> x;
        x.setRandom( leastSquares->rows(), 1 );
        gsMatrix<> rhs;
        rhs.setRandom( leastSquares->rows(), 1 ); // TODO
        gsMatrix<> errorHistory;
        gsConjugateGradient<> solver( leastSquares, preconder );
        solver.setCalcEigenvalues(true);
        solver.setOptions( cmd.getGroup("Solver") ).solveDetailed( rhs, x, errorHistory );

        gsInfo << "done.\n\n";

        iter2 = errorHistory.rows()-1;
        const bool success = errorHistory(iter2,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter2 << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter2 << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond2 = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond2 << "\n";
    }
    else
    {
        gsInfo << "skip.\n";
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Multigrid preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    gsInfo << "Setup of multigrid preconder... " << std::flush;
    if (mgPreconder)
    {

        std::vector<gsSparseMatrix<real_t, RowMajor>> embeddingSpace;
        std::vector<gsSparseMatrix<>> massesSpace;
        {
            gsGridHierarchy<> ghSpace = gsGridHierarchy<>::buildByCoarsening(mb, bc, cmd);
            index_t szGhSpace = ghSpace.getTransferMatrices().size();

            embeddingSpace.reserve(szGhSpace+1);
            massesSpace.reserve(szGhSpace+1);

            gsSparseMatrix<real_t, RowMajor> id(space_mass.rows(), space_mass.cols());
            id.setIdentity();
            embeddingSpace.push_back(id);
            massesSpace.push_back(space_mass);

            for (index_t i=0; i<szGhSpace; ++i)
            {
                gsSparseMatrix<real_t, RowMajor> e = ghSpace.getTransferMatrices()[szGhSpace-1-i];
                if (e.cols() < 2) break;
                id = id * e;
                embeddingSpace.push_back(id);
                massesSpace.push_back(id.transpose()*space_mass*id);
            }
        }

        std::vector<gsSparseMatrix<real_t, RowMajor>> embeddingTime;
        std::vector<gsSparseMatrix<>> massesTime;
        {
            gsGridHierarchy<> ghTime  = gsGridHierarchy<>::buildByCoarsening(tb1, ic, cmd);
            index_t szGhTime = ghTime.getTransferMatrices().size();

            embeddingTime.reserve(szGhTime+1);
            massesTime.reserve(szGhTime+1);

            gsSparseMatrix<real_t, RowMajor> id(time_mass1.rows(), time_mass1.cols());
            id.setIdentity();
            embeddingTime.push_back(id);
            massesTime.push_back(time_mass1);

            for (index_t i=0; i<szGhTime; ++i)
            {
                gsSparseMatrix<real_t, RowMajor> e = ghTime.getTransferMatrices()[szGhTime-1-i];
                // incorporate initial condition
                const index_t rr = e.rows(), cc = e.cols();
                if (cc < 3) break;
                id = id * e.block(1,1,rr-1,cc-1);
                embeddingTime.push_back(id);
                massesTime.push_back(id.transpose()*time_mass1*id);
            }
        }

        {
            gsInfo << "Space embedding matrices:";
            for (size_t i=0; i<embeddingSpace.size(); ++i)
                gsInfo << " " << embeddingSpace[i].rows() << "x" << embeddingSpace[i].cols();
            gsInfo << "\n";
            gsInfo << "Space mass matrices:";
            for (size_t i=0; i<massesSpace.size(); ++i)
                gsInfo << " " << massesSpace[i].rows() << "x" << massesSpace[i].cols();
            gsInfo << "\n";
            gsInfo << "Time embedding matrices:";
            for (size_t i=0; i<embeddingTime.size(); ++i)
                gsInfo << " " << embeddingTime[i].rows() << "x" << embeddingTime[i].cols();
            gsInfo << "\n";
            gsInfo << "Time mass matrices:";
            for (size_t i=0; i<massesTime.size(); ++i)
                gsInfo << " " << massesTime[i].rows() << "x" << massesTime[i].cols();
            gsInfo << "\n";
        }

        gsMatrix<> scaling(embeddingSpace.size()+1,embeddingTime.size()+1);
        scaling.setZero();
        for (size_t i=0; i<embeddingSpace.size(); ++i)
        {
            const real_t h = h0 * pow(0.5, embeddingSpace.size()-i);
            for (size_t j=0; j<embeddingTime.size(); ++j)
            {
                const real_t tau = tau0 * pow(0.5, embeddingTime.size()-j);
                scaling(i+1,j+1) = 1./( 1./(h*h) + h*h/(tau*tau) );
                gsInfo << "i=" << i << "; j=" << j << "; h=" << h << "; tau=" << tau << "; scaling=" << scaling(i,j) << "\n";
            }
        }

        gsInfo << "Actual scaling: ";
        gsAdditiveOp<>::Ptr mlPreconder = gsAdditiveOp<>::make();
        for (size_t i=0; i<embeddingSpace.size(); ++i)
            for (size_t j=0; j<embeddingTime.size(); ++j)
            {
                gsSparseMatrix<real_t,RowMajor> transfer = embeddingTime[j].kron(embeddingSpace[i]);
                gsLinearOperator<>::Ptr op = gsKroneckerOp<>::make( makeSparseCholeskySolver(massesTime[j]), makeSparseCholeskySolver(massesSpace[i]) );

                real_t sc = scaling(i+1,j+1) - scaling(i,j+1) - scaling(i+1,j) + scaling(i,j);
                gsInfo << sc << " ";

                op = gsScaledOp<>::make(give(op), sc);
                mlPreconder->addOperator(give(transfer), give(op));
            }
        gsInfo << "\n";

        gsInfo << "Setup cg solver and solve... " << std::flush;
        {
            gsMatrix<> x;
            x.setRandom( leastSquares->rows(), 1 );
            gsMatrix<> rhs;
            rhs.setRandom( leastSquares->rows(), 1 ); // TODO
            gsMatrix<> errorHistory;
            gsConjugateGradient<> solver( leastSquares, mlPreconder );
            solver.setCalcEigenvalues(true);
            solver.setOptions( cmd.getGroup("Solver") ).solveDetailed( rhs, x, errorHistory );

            gsInfo << "done.\n\n";

            iter1 = errorHistory.rows()-1;
            const bool success = errorHistory(iter1,0) < tolerance;
            if (success)
                gsInfo << "Reached desired tolerance after " << iter1 << " iterations:\n";
            else
                gsInfo << "Did not reach desired tolerance after " << iter1 << " iterations:\n";

            if (errorHistory.rows() < 20)
                gsInfo << errorHistory.transpose() << "\n\n";
            else
                gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

            cond1 = solver.getConditionNumber();
            gsInfo << "Estimated condition number: " << cond1 << "\n";
        }
    }
    else
    {
        gsInfo << "skip.\n";
    }

    gsInfo << "fin.\n";

    return 0;

}
