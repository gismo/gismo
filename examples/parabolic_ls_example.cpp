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
/*

struct makeMultiGridSolver {
    const gsMultiBasis<>& mbY;
    const gsBoundaryConditions<>& bcY;
    const index_t nTimeDofs;
    const gsOptionList& opt;
    gsLinearOperator<>::Ptr operator()( const gsSparseMatrix<>& m ) const
    {
        gsInfo << "Setup multigrid solver... " << std::flush;
        gsOptionList cmd;
        cmd.addInt( "InterfaceStrategy", "", iFace::conforming      );
        cmd.addInt( "DirichletStrategy", "", dirichlet::elimination );
        gsGridHierarchy<> gh = gsGridHierarchy<>::buildByCoarsening(mbY, bcY, cmd, 10000, 10);
        const index_t lv = gh.getTransferMatrices().size()+1;

        std::vector<gsSparseMatrix<real_t,RowMajor>> transferMatrices;
        transferMatrices.reserve(lv-1);
        gsSparseMatrix<> id(nTimeDofs,nTimeDofs);
        id.setIdentity();
        // consider using gsKroneckerOp<>
        for (index_t i=0; i<lv-1; ++i)
            transferMatrices.push_back( id.kron(gh.getTransferMatrices()[i]) );

        gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make( m, transferMatrices );
        mg->setOptions(opt);

        for (index_t i = 1; i < mg->numLevels(); ++i)
            mg->setSmoother(i, makeGaussSeidelOp(mg->matrix(i)));

        gsInfo << "done (" << lv << " grid levels).\n";

        return mg;
    }
};

gsSparseMatrix<> embeddingMatrix( const gsVector<index_t>& vec, index_t rows, index_t cols )
{
    gsSparseEntries<> se;
    se.reserve(cols);
    GISMO_ASSERT (vec.rows()>=cols, "not having "<<vec.rows()<<">="<<cols);
    for (index_t i=0; i<cols; ++i)
    {
        GISMO_ASSERT (vec(i,0)<rows, "not having "<<vec(i,0)<<"<"<<rows);
        se.add(vec(i,0),i,1);
    }
    gsSparseMatrix<> result(rows, cols);
    result.setFrom(se);
    return result;
}

gsMatrix<> tensorCoefs( const gsMatrix<>& A, const gsMatrix<>& B)
{
    gsMatrix<> result(A.rows()*B.rows(),A.cols()+B.cols());

    for (index_t i=0; i<B.rows(); ++i)
    {
        result.block( i*A.rows(), 0, A.rows(), A.cols() ) = A;
        for (index_t j=0; j<A.rows(); ++j)
            for (index_t k=0; k<B.cols(); ++k)
                result( i*A.rows()+j, A.cols()+k ) = B(i,k);
    }
    return result;
}

gsTensorBSplineBasis<3> tensorBasis(const gsMultiBasis<>& A, const gsMultiBasis<>& B)
{
    GISMO_ENSURE (A.dim()==1&&B.dim()==2,"Only implemented for 1D x 2D.");
    return gsTensorBSplineBasis<3>(
            dynamic_cast<gsBSplineBasis<>*>(A[0].component(0).clone().release()),
            dynamic_cast<gsBSplineBasis<>*>(B[0].component(0).clone().release()),
            dynamic_cast<gsBSplineBasis<>*>(B[0].component(1).clone().release())
        );
}

real_t getOpNorm(const gsSparseMatrix<>& mat)
{
    real_t max = 0;
    for (index_t i=0; i<mat.rows(); ++i)
       if (mat(i,i)>max) max = mat(i,i);
    return max;
}

gsSparseMatrix<real_t, RowMajor> mkTimeEmbedding(const gsSparseMatrix<real_t, RowMajor>& mat)
{
    // Compensate for initial conditions
    const index_t sz1 = mat.rows();
    gsSparseMatrix<> timeEmbedding1(sz1-1,sz1);
    for (index_t j=0; j<sz1-1; ++j)
        timeEmbedding1(j,j+1)=1;

    const index_t sz2 = mat.cols();
    gsSparseMatrix<> timeEmbedding2(sz2, sz2-1);
    for (index_t j=0; j<sz2-1; ++j)
        timeEmbedding2(j+1,j)=1;

    return timeEmbedding1 * mat * timeEmbedding2;
}

gsLinearOperator<>::Ptr makeSpaceTimeMultiGridSolver(
        gsSparseMatrix<> matrixA, gsSparseMatrix<> matrixB, gsSparseMatrix<> matrixC,
        const gsMultiBasis<>& mbY, const gsBoundaryConditions<>& bcY,
        const gsMultiBasis<>& tiY, const gsBoundaryConditions<>& icY,
        gsOptionList opt,
        std::string& info)
{

    const gsSparseMatrix<> fineMatrix = matrixA + matrixB + matrixC;
    gsInfo << "Setup space-time multigrid solver...\n" << std::flush;
    gsOptionList cmd;
    cmd.addInt( "InterfaceStrategy", "", iFace::conforming      );
    cmd.addInt( "DirichletStrategy", "", dirichlet::elimination );
    gsGridHierarchy<> ghSpace = gsGridHierarchy<>::buildByCoarsening(mbY, bcY, cmd);
    gsGridHierarchy<> ghTime  = gsGridHierarchy<>::buildByCoarsening(tiY, icY, cmd);
    index_t szGhSpace = ghSpace.getTransferMatrices().size();
    index_t szGhTime = ghTime.getTransferMatrices().size();

    std::vector<gsSparseMatrix<real_t, RowMajor>> transfers;

    GISMO_ENSURE( ghSpace.getTransferMatrices().size() > 0, "Need at least one transfer." );
    GISMO_ENSURE( ghTime.getTransferMatrices().size() > 0, "Need at least one transfer." );

    while (true)
    {
        const real_t B_norm = getOpNorm(matrixB);
        const real_t C_norm = getOpNorm(matrixC);
        gsInfo << "  B_norm: " << B_norm << ", C_norm: " << C_norm << ", ratio: " << B_norm/C_norm << "; ";
        if (B_norm < C_norm)
        {
            if (szGhTime<1) {
                gsInfo << "cannot coarsen in time.\n";
                break;
            }
            else
                gsInfo << "coarsen in time... " << std::flush;
            info.append("t");
            const index_t sz = szGhSpace>0 ? ghSpace.getTransferMatrices()[szGhSpace-1].rows() : ghSpace.getTransferMatrices()[szGhSpace].cols();
            gsSparseMatrix<real_t,RowMajor> transferSpace(sz,sz); transferSpace.setIdentity();
            gsSparseMatrix<real_t,RowMajor> transferTime = mkTimeEmbedding(ghTime.getTransferMatrices()[szGhTime-1]);
            transfers.push_back( transferTime.kron(transferSpace) );
            szGhTime -= 1;
        }
        else
        {
            if (szGhSpace<1) {
                gsInfo << "cannot coarsen in space.\n";
                break;
            }
            else
                gsInfo << "coarsen in space... " << std::flush;
            info.append("s");
            gsSparseMatrix<real_t,RowMajor> transferSpace = ghSpace.getTransferMatrices()[szGhSpace-1];
            const index_t sz = szGhTime>0 ? ghTime.getTransferMatrices()[szGhTime-1].rows() : ghTime.getTransferMatrices()[szGhTime].cols();
            gsSparseMatrix<real_t,RowMajor> transferTime(sz-1,sz-1); transferTime.setIdentity();
            transfers.push_back( transferTime.kron(transferSpace) );
            szGhSpace -= 1;
        }
        // update
        matrixB = transfers.back().transpose() * matrixB * transfers.back();
        matrixC = transfers.back().transpose() * matrixC * transfers.back();
        gsInfo << "done.\n";
    }

    std::reverse(transfers.begin(), transfers.end());

    gsInfo << "  Have " << transfers.size() << " transfer matrices.\n";
    GISMO_ENSURE( transfers.size() > 0, "Need at least one transfer." );

    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make( fineMatrix, transfers );
    mg->setOptions(opt);
    for (index_t i=1; i<mg->numLevels(); ++i)
        mg->setSmoother(i, makeGaussSeidelOp(mg->matrix(i)));
    gsInfo << "done.\n";

    return mg;

}

*/

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
    
    index_t iter1;
    real_t cond1;
    gsInfo << "Setup of somewhat exact preconder... " << std::flush;
    if (0)
    {
        const index_t primalDim = time_stiff1.rows() * space_stiff.rows();
    
        gsSparseMatrix<> preconderMatrix(2*primalDim,2*primalDim);
        gsSparseMatrix<> preconderSelector(primalDim,2*primalDim);
        for (index_t i=0; i<primalDim; ++i)
            preconderSelector(i,i) = 1;
        {
            gsSparseMatrix<> preconderBlock11 =  time_mass1.kron(space_stiff);
            gsSparseMatrix<> preconderBlock12 =  time_stiff1.kron(space_mass);
            gsSparseMatrix<> preconderBlock22 = -time_stiff1.kron(space_stiff);
        
            for (index_t i = 0; i < preconderBlock11.outerSize(); ++i)
                for (gsSparseMatrix<>::InnerIterator it(preconderBlock11,i); it; ++it)
                    preconderMatrix(it.row(), it.col()) = it.value();
    
            for (index_t i = 0; i < preconderBlock12.outerSize(); ++i)
                for (gsSparseMatrix<>::InnerIterator it(preconderBlock12,i); it; ++it)
                {
                    preconderMatrix(it.row(), primalDim+it.col()) = it.value();
                    preconderMatrix(primalDim+it.col(), it.row()) = it.value();
                }

            for (index_t i = 0; i < preconderBlock22.outerSize(); ++i)
                for (gsSparseMatrix<>::InnerIterator it(preconderBlock22,i); it; ++it)
                    preconderMatrix(primalDim+it.row(), primalDim+it.col()) = it.value();
            
        }

        //gsInfo << "\npreconderMatrix=\n" << preconderMatrix << "\n";

    
        gsLinearOperator<>::Ptr preconder
            = gsProductOp<>::make( makeMatrixOp(preconderSelector.transpose()), makeSparseLUSolver(preconderMatrix), makeMatrixOp(preconderSelector) );

    
        gsInfo << "done.\n";
    
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
    

    gsInfo << "fin.\n";
    
    return 0;

}
