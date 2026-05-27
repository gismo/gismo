/** @file parabolic_l2ls_example.cpp

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

real_t getOpNorm(const gsSparseMatrix<>& mat) // used by makeSpaceTimeMultiGridSolver
{
    real_t max = 0;
    for (index_t i=0; i<mat.rows(); ++i)
       if (mat(i,i)>max) max = mat(i,i);
    return max;
}

gsSparseMatrix<real_t, RowMajor> mkTimeEmbedding(const gsSparseMatrix<real_t, RowMajor>& mat) // used by makeSpaceTimeMultiGridSolver
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

gsSparseMatrix<> forLevel(gsSparseMatrix<> matrix, index_t lv, const std::vector<gsSparseMatrix<real_t, RowMajor>>& transfers)
{
    for (index_t i=transfers.size()-1; i>=lv; --i)
    {
        GISMO_ASSERT (transfers[i].rows() == matrix.rows(), transfers[i].cols()<<","<<transfers[i].rows() <<"=="<< matrix.rows());
        matrix = transfers[i].transpose() * matrix * transfers[i];
    }
    return matrix;
}

gsLinearOperator<>::Ptr makeSpaceTimeMultiGridSolver(
        gsSparseMatrix<> matrixA, gsSparseMatrix<> matrixB, gsSparseMatrix<> massMatrix,
        const gsMultiBasis<>& mbY, const gsBoundaryConditions<>& bcY,
        const gsMultiBasis<>& tiY, const gsBoundaryConditions<>& icY,
        gsOptionList opt,
        std::string& info,
        char smoother_type
        )
{

    const gsSparseMatrix<> fineMatrix = matrixA + matrixB; //+ matrixC;
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
        const real_t A_norm = getOpNorm(matrixA);
        const real_t B_norm = getOpNorm(matrixB);
        gsInfo << "  A_norm: " << A_norm << ", B_norm: " << B_norm << ", ratio: " << A_norm/B_norm << "; ";
        if (A_norm > B_norm)
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
        matrixA = transfers.back().transpose() * matrixA * transfers.back();
        matrixB = transfers.back().transpose() * matrixB * transfers.back();
        gsInfo << "done.\n";
    }

    std::reverse(transfers.begin(), transfers.end());

    gsInfo << "  Have " << transfers.size() << " transfer matrices.\n";
    GISMO_ENSURE( transfers.size() > 0, "Need at least one transfer." );

    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make( fineMatrix, transfers );
    mg->setOptions(opt);
    for (index_t i=1; i<mg->numLevels(); ++i)
    {
        if (smoother_type=='g')
            mg->setSmoother(i, makeGaussSeidelOp(mg->matrix(i)));
        else if (smoother_type=='m')
        {
            gsSparseMatrix<> mass_lv = forLevel(massMatrix, i, transfers);
            GISMO_ASSERT (mass_lv.rows() == mg->matrix(i).rows(), mass_lv.rows()<<" == "<<mg->matrix(i).rows());
            mass_lv *= getOpNorm(mg->matrix(i)) / getOpNorm(mass_lv) * opt.getReal("SigmaMass");
            mass_lv += mg->matrix(i) * opt.getReal("SigmaStiff");
            gsLinearOperator<>::Ptr smootherOp = makeSparseCholeskySolver(mass_lv);
            mg->setSmoother(i, gsPreconditionerFromOp<>::make(mg->underlyingOp(i),smootherOp));
        }
        else
            GISMO_ENSURE(false, "Unknown smoother type.");
    }
    gsInfo << "done.\n";

    return mg;

}




int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    index_t geoIdx = 1;
    index_t refinementsX = 2;
    index_t refinementsT = 2;
    index_t degree = 2;
    index_t dualAddMultiplicityX = 0;
    index_t dualAddMultiplicityT = 1;
    real_t kappa = 1.;
    index_t maxIterations = 100;
    real_t tolerance = 1.e-6;
    index_t preSmooth = 1;
    index_t postSmooth = 1;
    index_t cycles = 1;

    index_t exactPreconder = 1;
    index_t fdPreconder = 1;
    index_t mggsPreconder = 0;
    index_t mgmsPreconder = 0;
    real_t sigmaMass = 1.;
    real_t sigmaStiff = 1.;

    std::string out;

    gsCmdLine cmd("parabolic_l2ls_example");
    cmd.addInt   ("g", "Geometry",              "0=Rectangle, 1=Quarter Annulus", geoIdx);
    cmd.addInt   ("r", "RefinementsX",          "Number of uniform h-refinement steps to perform before solving", refinementsX);
    cmd.addInt   ("s", "RefinementsT",          "Number of uniform tau-refinement steps to perform before solving", refinementsT);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addInt   ("",  "DualAddMultiplicityX",  "Additional multiplicity in space for dual space", dualAddMultiplicityX);
    cmd.addInt   ("",  "DualAddMultiplicityT",  "Additional multiplicity in time for dual space", dualAddMultiplicityT);
    cmd.addReal  ("k", "Kappa",                 "Diffusion parameter", kappa);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Maximum iterations for linear solver", maxIterations);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "MG.NumPreSmooth",       "Number of pre smoothing steps (only for mg)", preSmooth);
    cmd.addInt   ("",  "MG.NumPostSmooth",      "Number of post smoothing steps (only for mg)", postSmooth);
    cmd.addInt   ("",  "MG.NumCycles",          "Number of multi-grid cycles for coarse-grid correction, i.e., 1=V, 2=W cycle", cycles);

    cmd.addInt   ("",  "useExactPreconder",     "Use that scheme", exactPreconder);
    cmd.addInt   ("",  "useFdPreconder",        "Use that scheme", fdPreconder);
    cmd.addInt   ("",  "useMgGsPreconder",      "Use that scheme", mggsPreconder);
    cmd.addInt   ("",  "useMgMsPreconder",      "Use that scheme", mgmsPreconder);
    cmd.addReal  ("",  "MG.SigmaMass",          "Sigma for mass", sigmaMass);
    cmd.addReal  ("",  "MG.SigmaStiff",         "Sigma for stiff", sigmaStiff);

    cmd.addString("",  "out",                   "Write solution and used options to file", out);

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
    gsMultiBasis<> mb1(mp);
    gsMultiBasis<> mb2(mp);
    gsMultiBasis<> tb1(tp);
    gsMultiBasis<> tb2(tp);

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    for ( size_t i = 0; i < mb1.nBases(); ++ i )
        mb1[i].setDegreePreservingMultiplicity(degree);

    for ( size_t i = 0; i < mb2.nBases(); ++ i )
        mb2[i].setDegreePreservingMultiplicity(degree);

    for ( size_t i = 0; i < tb1.nBases(); ++ i )
        tb1[i].setDegreePreservingMultiplicity(degree);

    for ( size_t i = 0; i < tb2.nBases(); ++ i )
        tb2[i].setDegreePreservingMultiplicity(degree);

    for ( index_t i = 0; i < refinementsX; ++i )
        mb1.uniformRefine();

    for ( index_t i = 0; i < refinementsX; ++i )
        mb2.uniformRefine();

    for ( index_t i = 0; i < refinementsT; ++i )
        tb1.uniformRefine();

    for ( index_t i = 0; i < refinementsT; ++i )
        tb2.uniformRefine();

    GISMO_ENSURE (dualAddMultiplicityX==0 || dualAddMultiplicityX==2, "Non-valid value "<<dualAddMultiplicityX);
    if (dualAddMultiplicityX)
        for ( size_t i = 0; i < mb2.nBases(); ++ i )
            mb2[i].reduceContinuity(dualAddMultiplicityX);


    GISMO_ENSURE (dualAddMultiplicityT==0 || dualAddMultiplicityT==1, "Non-valid value "<<dualAddMultiplicityT);
    if (dualAddMultiplicityT)
        for ( size_t i = 0; i < tb2.nBases(); ++ i )
            tb2[i].reduceContinuity(dualAddMultiplicityT);

    gsInfo << "done.\n";
    gsInfo << "Knots in tb1: " << dynamic_cast<gsTensorBSplineBasis<1,real_t>&>(tb1[0]).component(0).knots() << "\n";
    gsInfo << "Knots in tb2: " << dynamic_cast<gsTensorBSplineBasis<1,real_t>&>(tb2[0]).component(0).knots() << "\n";



    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrices... " << std::flush;


    // Assemble space matrices
    gsSparseMatrix<> space_mass1;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(mb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space       u = assembler.getSpace(mb1,1,0);
        bc.setGeoMap(mp);
        u.setup(bc, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( u * u.tr() * meas(G) );
        space_mass1 = assembler.matrix();
    }

    gsSparseMatrix<> space_mass2;
    if (dualAddMultiplicityX==0)
        space_mass2 = space_mass1;
    else
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(mb2);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space       u = assembler.getSpace(mb2,1,0);
        bc.setGeoMap(mp);
        //u.setup(bc, dirichlet::interpolation, 0); // no boundary conditions
        assembler.initSystem();
        assembler.assemble( u * u.tr() * meas(G) );
        space_mass2 = assembler.matrix();
    }

    gsSparseMatrix<> space_massL;
    if (dualAddMultiplicityX==0)
        space_massL = space_mass1;
    else
    {
        gsExprAssembler<> assembler(2,2);
        assembler.setIntegrationElements(mb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space u = assembler.getSpace(mb1,1,0);
        gsExprAssembler<>::space v = assembler.getSpace(mb2,1,1);
        bc.setGeoMap(mp);
        u.setup(bc, dirichlet::interpolation, 0);
        //v.setup(bc, dirichlet::interpolation, 0); // no boundary conditions
        assembler.initSystem();
        assembler.assemble( u * v.tr() * meas(G) );
        space_massL = assembler.matrix();

        const index_t n1 = space_mass1.rows();
        const index_t n2 = space_mass2.rows();
        GISMO_ENSURE(space_massL.rows() == n1+n2, "");
        space_massL = space_massL.block(0, n1, n1, n2).transpose(); //TODO

    }

    gsSparseMatrix<> space_stiffL;
    if (dualAddMultiplicityX==0)
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(mb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space u = assembler.getSpace(mb1,1,0);
        bc.setGeoMap(mp);
        u.setup(bc, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G) );
        space_stiffL = assembler.matrix();
    }
    else
    {
        gsExprAssembler<> assembler(2,2);
        assembler.setIntegrationElements(mb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space u = assembler.getSpace(mb1,1,0);
        gsExprAssembler<>::space v = assembler.getSpace(mb2,1,1);
        bc.setGeoMap(mp);
        u.setup(bc, dirichlet::interpolation, 0);
        //v.setup(bc, dirichlet::interpolation, 0); // no boundary conditoons
        assembler.initSystem();
        assembler.assemble( igrad(u,G) * igrad(v,G).tr() * meas(G) );
        space_stiffL = assembler.matrix();

        const index_t n1 = space_mass1.rows();
        const index_t n2 = space_mass2.rows();
        GISMO_ENSURE(space_stiffL.rows() == n1+n2, "");
        space_stiffL = space_stiffL.block(0, n1, n1, n2).transpose(); //TODO

    }

    gsSparseMatrix<> space_biharm1;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(mb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space       u = assembler.getSpace(mb1,1,0);
        bc.setGeoMap(mp);
        u.setup(bc, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( ilapl(u,G) * ilapl(u,G).tr() * meas(G) );
        space_biharm1 = assembler.matrix();
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
    if (dualAddMultiplicityT==0)
        time_mass2 = time_mass1;
    else
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(tb2);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(tp);
        gsExprAssembler<>::space u = assembler.getSpace(tb2,1,0);
        ic.setGeoMap(tp);
        u.setup(ic, dirichlet::interpolation, 0); // no initial condition
        assembler.initSystem();
        assembler.assemble( u * u.tr() * meas(G) );
        time_mass2 = assembler.matrix();
    }

    gsSparseMatrix<> time_massL;
    if (dualAddMultiplicityT==0)
        time_massL = time_mass1;
    else
    {
        gsExprAssembler<> assembler(2,2);
        assembler.setIntegrationElements(tb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(tp);
        gsExprAssembler<>::space u = assembler.getSpace(tb1,1,0);
        gsExprAssembler<>::space v = assembler.getSpace(tb2,1,1);
        ic.setGeoMap(tp);
        u.setup(ic, dirichlet::interpolation, 0); // initial condition below
        v.setup(ic, dirichlet::interpolation, 0); // no initial condition
        assembler.initSystem();
        assembler.assemble( u * v.tr() * meas(G) );
        time_massL = assembler.matrix();

        const index_t n1 = time_mass1.rows()+1;
        const index_t n2 = time_mass2.rows();
        GISMO_ENSURE(time_massL.rows() == n1+n2, time_massL.rows()<<"=="<<n1<<"+"<<n2);
        time_massL = time_massL.block(1, n1, n1-1, n2).transpose(); //TODO
    }


    gsSparseMatrix<> time_gradL;
    if (dualAddMultiplicityT==0)
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(tb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(tp);
        gsExprAssembler<>::space u = assembler.getSpace(tb1,1,0);
        ic.setGeoMap(tp);
        u.setup(ic, dirichlet::interpolation, 0); // initial condition below
        assembler.initSystem();
        assembler.assemble( igrad(u,G) * u.tr() * meas(G) );
        time_gradL = assembler.matrix();
        const index_t n1 = time_mass1.rows();
        time_gradL = time_gradL.block(1, 1, n1, n1).transpose(); //TODO
    }
    else
    {
        gsExprAssembler<> assembler(2,2);
        assembler.setIntegrationElements(tb1);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(tp);
        gsExprAssembler<>::space u = assembler.getSpace(tb1,1,0);
        gsExprAssembler<>::space v = assembler.getSpace(tb2,1,1);
        ic.setGeoMap(tp);
        u.setup(ic, dirichlet::interpolation, 0); // initial condition below
        v.setup(ic, dirichlet::interpolation, 0); // no initial condition
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
        gsInfo << "\nspace_mass1=\n" << space_mass1 << "\n";
        gsInfo << "\nspace_mass2=\n" << space_mass2 << "\n";
        gsInfo << "\nspace_massL=\n" << space_massL << "\n";
        gsInfo << "\nspace_stiffL=\n" << space_stiffL << "\n";
        gsInfo << "\nspace_biharm1=\n" << space_biharm1 << "\n";

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
            gsKroneckerOp<>::make( makeMatrixOp(time_gradL), makeMatrixOp(space_massL) ),
            gsScaledOp<>::make( gsKroneckerOp<>::make( makeMatrixOp(time_massL), makeMatrixOp(space_stiffL) ), kappa )
        );
    gsLinearOperator<>::Ptr LhT
        = gsSumOp<>::make(
            gsKroneckerOp<>::make( makeMatrixOp(time_gradL.transpose()), makeMatrixOp(space_massL.transpose()) ),
            gsScaledOp<>::make( gsKroneckerOp<>::make( makeMatrixOp(time_massL.transpose()), makeMatrixOp(space_stiffL.transpose()) ), kappa )
        );
    gsLinearOperator<>::Ptr dualPc
        = gsKroneckerOp<>::make( makeSparseCholeskySolver(time_mass2), makeSparseCholeskySolver(space_mass2) );

    gsLinearOperator<>::Ptr leastSquares = gsProductOp<>::make( Lh, dualPc, LhT );

    gsInfo << "done; " << leastSquares->rows() << " dofs.\n";

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Somewhat exact preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    std::array<index_t,4> iters = {-1,-1,-1,-1};
    std::array<real_t,4> conds  = {-1,-1,-1,-1};


    gsInfo << "Setup of somewhat exact preconder... " << std::flush;
    if (exactPreconder)
    {
        index_t &iter = iters[0];
        real_t &cond = conds[0];

        gsSparseMatrix<> preconderMatrix = time_stiff1.kron(space_mass1) + (kappa*kappa)*time_mass1.kron(space_biharm1);

        gsLinearOperator<>::Ptr preconder = makeSparseCholeskySolver(preconderMatrix);


        gsInfo << "done: " << preconder->rows() << " dofs.\n";

        gsInfo << "Setup cg solver and solve... " << std::flush;

        gsMatrix<> x;
        x.setRandom( leastSquares->rows(), 1 );
        gsMatrix<> rhs;
        rhs.setRandom( leastSquares->rows(), 1 );
        gsMatrix<> errorHistory;
        gsConjugateGradient<> solver( leastSquares, preconder );
        solver.setCalcEigenvalues(true);
        solver.setOptions( cmd.getGroup("Solver") ).solveDetailed( rhs, x, errorHistory );

        gsInfo << "done.\n\n";

        iter = errorHistory.rows()-1;
        const bool success = errorHistory(iter,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond << "\n";
    }
    else
    {
        gsInfo << "skip.\n";
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  FD preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    gsInfo << "Setup of FD preconder... " << std::flush; // FD in time, direct in space
    if (fdPreconder)
    {
        index_t &iter = iters[1];
        real_t &cond = conds[1];

        gsLinearOperator<>::Ptr preconder = fastDiagnonalization(time_stiff1, space_mass1, (kappa*kappa)*time_mass1, space_biharm1, mkSparseLUSolver);

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

        iter = errorHistory.rows()-1;
        const bool success = errorHistory(iter,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond << "\n";
    }
    else
    {
        gsInfo << "skip.\n";
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  MG+GS preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo << "Setup of MG+GS preconder... " << std::flush; // MG in space-time, Gauss-Seidel smoother
    if (mggsPreconder)
    {

        index_t &iter = iters[2];
        real_t &cond = conds[2];

        gsSparseMatrix<> massMatrix;

        std::string info;
        gsLinearOperator<>::Ptr preconder = makeSpaceTimeMultiGridSolver(
            time_stiff1.kron(space_mass1), (kappa*kappa)*time_mass1.kron(space_biharm1), massMatrix,
            mb1, bc,
            tb1, ic,
            cmd.getGroup("MG"),
            info, 'g');


        //gsInfo << "\npreconderMatrix=\n" << preconderMatrix << "\n";
        gsInfo << "done: " << preconder->rows() << " dofs; type: " << info << "\n";

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

        iter = errorHistory.rows()-1;
        const bool success = errorHistory(iter,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond << "\n";
    }
    else
    {
        gsInfo << "skip.\n";
    }


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  MG+MS preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    gsInfo << "Setup of MG+MS preconder... " << std::flush; // MG in space-time, mass smoother (no subspace correction)
    if (mgmsPreconder)
    {

        index_t &iter = iters[3];
        real_t &cond = conds[3];

        std::string info;
        gsLinearOperator<>::Ptr preconder = makeSpaceTimeMultiGridSolver(
            time_stiff1.kron(space_mass1), (kappa*kappa)*time_mass1.kron(space_biharm1), time_mass1.kron(space_mass1),
            mb1, bc,
            tb1, ic,
            cmd.getGroup("MG"),
            info, 'm');


        //gsInfo << "\npreconderMatrix=\n" << preconderMatrix << "\n";
        gsInfo << "done: " << preconder->rows() << " dofs; type: " << info << "\n";

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

        iter = errorHistory.rows()-1;
        const bool success = errorHistory(iter,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond << "\n";
    }
    else
    {
        gsInfo << "skip.\n";
    }


    if (!out.empty())
    {
        const bool exists = gsFileManager::fileExists(out);
        std::ofstream outfile;
        outfile.open(out.c_str(), std::ios_base::app);


        if (!exists)
            outfile << "parabolic_l2ls_example\t"
                "geoIdx\t"
                "refinementsX\t"
                "refinementsT\t"
                "degree\t"
                "dualAddMultiplicityX\t"
                "dualAddMultiplicityT\t"
                "kappa\t"
                "sigmaMass\t"
                "sigmaStiff\t"
                "preSmooth\t"
                "postSmooth\t"
                "cycles\t"
                "iter0\t"
                "cond0\t"
                "iter1\t"
                "cond1\t"
                "iter2\t"
                "cond2\t"
                "iter3\t"
                "cond3\n";

        outfile << "parabolic_l2ls_example\t"
            << geoIdx << "\t"
            << refinementsX << "\t"
            << refinementsT << "\t"
            << degree << "\t"
            << dualAddMultiplicityX << "\t"
            << dualAddMultiplicityT << "\t"
            << kappa << "\t"
            << sigmaMass << "\t"
            << sigmaStiff << "\t"
            << preSmooth << "\t"
            << postSmooth << "\t"
            << cycles << "\t"
            << iters[0] << "\t"
            << conds[0] << "\t"
            << iters[1] << "\t"
            << conds[1] << "\t"
            << iters[2] << "\t"
            << conds[2] << "\t"
            << iters[3] << "\t"
            << conds[3] << "\n";
    }


    return 0;

}
