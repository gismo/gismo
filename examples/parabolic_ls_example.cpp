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

template<typename S>
gsLinearOperator<>::Ptr mkFdpfPreconder(
    const gsSparseMatrix<>& time_stiff, const gsSparseMatrix<>& space_stiff,
    const gsSparseMatrix<>& time_mass, const gsSparseMatrix<>& space_mass, const S& makeSolver)
{
    GISMO_ASSERT(time_mass.rows() == time_mass.cols() && time_mass.rows() == time_stiff.rows() && time_stiff.rows() == time_stiff.cols(), "");
    GISMO_ASSERT(space_stiff.rows() == space_stiff.cols() && space_stiff.rows() == space_mass.rows() && space_mass.rows() == space_mass.cols(), "");

    typedef gsMatrix<>::GenSelfAdjEigenSolver EVSolver;
    EVSolver ges;
    ges.compute(time_stiff, time_mass, gsEigen::ComputeEigenvectors);
    gsSparseMatrix<> D1(time_mass.rows(), time_mass.rows()), D2(time_mass.rows(), time_mass.rows());
    GISMO_ASSERT (ges.eigenvalues().rows() == time_mass.rows() && ges.eigenvalues().cols() == 1, "");
    GISMO_ASSERT (ges.eigenvectors().rows() == time_mass.rows() && ges.eigenvectors().cols() == time_mass.rows(), "");
    for (index_t i=0; i<time_mass.rows(); ++i)
    {
        D1(i,i) = 1;
        D2(i,i) = sqrt(ges.eigenvalues()(i,0));
    }
    gsSparseMatrix<> system = D1.kron(space_stiff)+D2.kron(space_mass); //TODO: is this orderd such that a direct solver would do this efficiently?
    gsLinearOperator<>::Ptr systemSolver = makeSolver(system);
    gsMatrix<> eigs = ges.eigenvectors();
    gsMatrix<> eigsT = ges.eigenvectors().transpose();

    gsSparseMatrix<> multiplierMatrix = D1.kron(space_stiff);
    gsLinearOperator<>::Ptr multiplier = makeMatrixOp(multiplierMatrix.moveToPtr());

    gsLinearOperator<>::Ptr transform  = gsKroneckerOp<>::make( makeMatrixOp(eigs.moveToPtr()), gsIdentityOp<>::make(space_stiff.rows()) );
    gsLinearOperator<>::Ptr transformT = gsKroneckerOp<>::make( makeMatrixOp(eigsT.moveToPtr()), gsIdentityOp<>::make(space_stiff.rows()) );
    gsProductOp<>::Ptr result = gsProductOp<>::make();
    result->addOperator(transformT);
    result->addOperator(systemSolver);
    result->addOperator(multiplier);
    result->addOperator(systemSolver);
    result->addOperator(transform);
    return result;
}



gsLinearOperator<>::Ptr mkSparseLUSolver(const gsSparseMatrix<>& m) { return makeSparseLUSolver(m); }


struct mkMultiGridSolver {
    const gsMultiBasis<>& mb_space;
    const gsBoundaryConditions<>& bc_space;
    const index_t nTimeDofs;
    const gsOptionList& opt;
    gsLinearOperator<>::Ptr operator()( const gsSparseMatrix<>& m ) const
    {
        gsInfo << "Setup multigrid solver... " << std::flush;
        gsOptionList cmd;
        cmd.addInt( "InterfaceStrategy", "", iFace::conforming      );
        cmd.addInt( "DirichletStrategy", "", dirichlet::elimination );
        gsGridHierarchy<> gh = gsGridHierarchy<>::buildByCoarsening(mb_space, bc_space, cmd, 10000, 10);
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

template<typename S>
gsLinearOperator<>::Ptr mkTimeMultiLevelPreconder(
    const gsSparseMatrix<>& time_stiff, const gsSparseMatrix<>& space_stiff,
    const gsSparseMatrix<>& time_mass, const gsSparseMatrix<>& space_mass,
    const gsMultiBasis<>& mb_time, const gsBoundaryConditions<>& bc_time, const gsOptionList& opt,
    const S& makeSolver)
{
        gsInfo << "Setup multilevel solver... " << std::flush;
        gsOptionList cmd;
        cmd.addInt( "InterfaceStrategy", "", iFace::conforming      );
        cmd.addInt( "DirichletStrategy", "", dirichlet::elimination );
        gsGridHierarchy<> gh = gsGridHierarchy<>::buildByCoarsening(mb_time, bc_time, cmd, 10000, 2);
        const index_t lv = gh.getTransferMatrices().size()+1;

        // consider using gsKroneckerOp<>
        std::vector<gsSparseMatrix<real_t,RowMajor>> transferMatrices;
        {
            transferMatrices.reserve(lv-1);
            const index_t nSpaceDofs = space_stiff.rows();
            gsSparseMatrix<> id(nSpaceDofs,nSpaceDofs);
            id.setIdentity();
            for (index_t i=0; i<lv-1; ++i)
            {
                gsSparseMatrix<real_t,RowMajor> tmp = gh.getTransferMatrices()[i];
                gsSparseMatrix<real_t,RowMajor> tmp2 = tmp.block(1, 1, tmp.rows()-1, tmp.cols()-1);
                transferMatrices.push_back( tmp2.kron(id) );
            }
        }

        // accumulated transfer matrices
        std::vector<gsSparseMatrix<real_t,RowMajor>> accumulatedTransferMatrices;
        {
            gsSparseMatrix<real_t,RowMajor> tr(time_stiff.rows() * space_stiff.rows(), time_stiff.rows() * space_stiff.rows());
            tr.setIdentity();
            accumulatedTransferMatrices.push_back(tr);
            for (index_t i=0; i<lv-1; ++i)
            {
                tr = tr*transferMatrices[lv-2-i];
                accumulatedTransferMatrices.push_back( tr );
            }
            //for (index_t i=0; i<lv; ++i)
            //    gsInfo << "accumulatedTransferMatrices["<<i<<"] = "<<accumulatedTransferMatrices[i].rows()<<" x "<<accumulatedTransferMatrices[i].cols() << std::endl;
        }

        gsSumOp<>::Ptr preconder = gsSumOp<>::make();
        real_t tau = pow(.5,lv) * opt.askReal("Sigma", 1); // TODO
        for (index_t i=0; i<lv; ++i)
        {

            const index_t nSpaceDofs = space_stiff.rows();
            const index_t nTimeDofs = accumulatedTransferMatrices[i].cols() / nSpaceDofs;
            gsSparseMatrix<> id(nTimeDofs, nTimeDofs); id.setIdentity();

            gsSparseMatrix<> system = ( gsSparseMatrix<>(sqrt(tau)*id).kron(space_stiff) ) + ( gsSparseMatrix<>(sqrt(1/tau)*id).kron(space_mass) );
            gsLinearOperator<>::Ptr systemSolver = makeSolver(system);

            gsSparseMatrix<> multiplierMatrix = id.kron(space_stiff);
            gsLinearOperator<>::Ptr multiplier = makeMatrixOp(multiplierMatrix.moveToPtr());


            gsSparseMatrix<real_t,RowMajor> transformMatrix  = accumulatedTransferMatrices[i];
            gsSparseMatrix<real_t>          transformTMatrix = accumulatedTransferMatrices[i].transpose();
            gsLinearOperator<>::Ptr transform  = makeMatrixOp( transformMatrix .moveToPtr() );
            gsLinearOperator<>::Ptr transformT = makeMatrixOp( transformTMatrix.moveToPtr() );

            gsProductOp<>::Ptr lvPreconder = gsProductOp<>::make();
            lvPreconder->addOperator(transformT);
            lvPreconder->addOperator(systemSolver);
            lvPreconder->addOperator(multiplier);
            lvPreconder->addOperator(systemSolver);
            lvPreconder->addOperator(transform);
            preconder->addOperator(lvPreconder);

            tau *= 2; // time grid size
        }

        gsInfo << "done (" << lv << " grid levels).\n";

        return preconder;
}


int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    index_t geoIdx = 1;
    index_t refinementsX = 2;
    index_t refinementsT = 2;
    index_t degree = 1;
    real_t kappa = 1.;
    real_t sigma = 1.;
    index_t maxIterations = 100;
    real_t tolerance = 1.e-6;
    index_t preSmooth = 1;
    index_t postSmooth = 1;
    index_t cycles = 1;
    index_t exactPreconder = 1;
    index_t fdPreconder = 1;
    index_t fdpfPreconder = 1;
    index_t fdmgPreconder = 1;
    index_t mlluPreconder = 1;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("parabolic_ls_example");
    cmd.addInt   ("g", "Geometry",              "0=Rectangle, 1=Quarter Annulus", geoIdx);
    cmd.addInt   ("r", "RefinementsX",          "Number of uniform h-refinement steps to perform before solving", refinementsX);
    cmd.addInt   ("s", "RefinementsT",          "Number of uniform tau-refinement steps to perform before solving", refinementsT);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addReal  ("k", "Kappa",                 "Diffusion parameter", kappa);
    cmd.addReal  ("",  "Sigma",                 "Sigma for time-multilevel", sigma);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Maximum iterations for linear solver", maxIterations);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "MG.NumPreSmooth",       "Number of pre smoothing steps (only for mg)", preSmooth);
    cmd.addInt   ("",  "MG.NumPostSmooth",      "Number of post smoothing steps (only for mg)", postSmooth);
    cmd.addInt   ("c", "MG.NumCycles",          "Number of multi-grid cycles for coarse-grid correction, i.e., 1=V, 2=W cycle", cycles);
    cmd.addInt   ("",  "ExactPreconder",        "Use that scheme", exactPreconder);
    cmd.addInt   ("",  "FdPreconder",           "Use that scheme", fdPreconder);
    cmd.addInt   ("",  "FdPfPreconder",         "Use that scheme", fdpfPreconder);
    cmd.addInt   ("",  "FdMgPreconder",         "Use that scheme", fdmgPreconder);
    cmd.addInt   ("",  "MlLuPreconder",         "Use that scheme", mlluPreconder);
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
    gsInfo << "Setup of FD preconder... " << std::flush; // FD in time, not in space
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
    //  FD+PF preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    index_t iter3;
    real_t cond3;
    gsInfo << "Setup of FD+PF preconder... " << std::flush; // FD in time, Pearson factorzation in space
    if (fdpfPreconder)
    {
        gsLinearOperator<>::Ptr preconder = mkFdpfPreconder(time_stiff1, space_stiff, time_mass1, space_mass, mkSparseLUSolver);

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

        iter3 = errorHistory.rows()-1;
        const bool success = errorHistory(iter3,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter3 << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter3 << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond3 = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond3 << "\n";
    }
    else
    {
        gsInfo << "skip.\n";
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  FD+MG preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    index_t iter4;
    real_t cond4;
    gsInfo << "Setup of FD+MG preconder... " << std::flush; // FD in time, mg (using Pearson factorzation) in space
    if (fdmgPreconder)
    {
        const index_t primalDimTime = time_stiff1.rows();

        // Todo: make multigrid configurable

        gsLinearOperator<>::Ptr preconder = mkFdpfPreconder(time_stiff1, space_stiff, time_mass1, space_mass, mkMultiGridSolver{mb,bc,primalDimTime,cmd} );

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

        iter4 = errorHistory.rows()-1;
        const bool success = errorHistory(iter4,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter4 << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter4 << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond4 = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond4 << "\n";
    }
    else
    {
        gsInfo << "skip.\n";
    }



    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  ML+LU preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    index_t iter5;
    real_t cond5;
    gsInfo << "Setup of ML+LU preconder... " << std::flush; // ML in time, mg (using Pearson factorzation) in space
    if (mlluPreconder)
    {

        gsLinearOperator<>::Ptr preconder = mkTimeMultiLevelPreconder(time_stiff1, space_stiff, time_mass1, space_mass, tb1, ic, cmd, mkSparseLUSolver /**mg**/);

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

        iter5 = errorHistory.rows()-1;
        const bool success = errorHistory(iter5,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter5 << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter5 << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond5 = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond5 << "\n";
    }
    else
    {
        gsInfo << "skip.\n";
    }


    return 0;

}
