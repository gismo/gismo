/** @file parabolic_oc_example.cpp

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



template<typename S>
gsLinearOperator<>::Ptr fd(const gsSparseMatrix<>& A1, const gsSparseMatrix<>& B1, const gsSparseMatrix<>& A2, const gsSparseMatrix<>& B2, const S& makeSolver)
{
    GISMO_ASSERT(A1.rows() == A1.cols() && A1.rows() == A2.rows() && A2.rows() == A2.cols(), "");
    GISMO_ASSERT(B1.rows() == B1.cols() && B1.rows() == B2.rows() && B2.rows() == B2.cols(), "");

    typedef gsMatrix<>::GenSelfAdjEigenSolver EVSolver;
    EVSolver ges;
    ges.compute(A2, A1, gsEigen::ComputeEigenvectors);
    gsSparseMatrix<> D1(A1.rows(), A1.rows()), D2(A1.rows(), A1.rows());
    GISMO_ASSERT( ges.eigenvalues().rows() == A1.rows() && ges.eigenvalues().cols() == 1, "");
    GISMO_ASSERT( ges.eigenvectors().rows() == A1.rows() && ges.eigenvectors().cols() == A1.rows(), "");
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

gsLinearOperator<>::Ptr mkSparseCholeskySolver(const gsSparseMatrix<>& m) { return makeSparseCholeskySolver(m); }

struct mkMGSolver {
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


int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    index_t geoIdx = 2;
    index_t refinementsX = 1;
    index_t refinementsT = 1;
    index_t degree = 3;
    //real_t timeHorizon = 1.;
    real_t tolerance = 1.e-6;
    real_t alpha = 1.;
    real_t kappa = 0.01;
    index_t obs = 1;
    index_t maxIterations = 1000;
    index_t type = 1;
    index_t preSmooth = 1;
    index_t postSmooth = 1;
    index_t cycles = 1;
    std::string out;
    bool plot = false;

    gsCmdLine cmd("parabolic_oc_example");
    cmd.addInt   ("g", "Geometry",              "1=Rectangle, 2=Quarter Annulus", geoIdx);
    cmd.addInt   ("r", "RefinementsX",          "Number of uniform h-refinement steps to perform before solving", refinementsX);
    cmd.addInt   ("s", "RefinementsT",          "Number of uniform tau-refinement steps to perform before solving", refinementsT);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addReal  ("a", "Alpha",                 "Regularization parameter", alpha);
    cmd.addReal  ("k", "Kappa",                 "Diffusion parameter", kappa);
    cmd.addInt   ("o", "ObservationType",       "0=no observation, 1=limited observation, 2=full observation", obs);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Maximum iterations for linear solver", maxIterations);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "MG.NumPreSmooth",       "Number of pre smoothing steps (only for mg)", preSmooth);
    cmd.addInt   ("",  "MG.NumPostSmooth",      "Number of post smoothing steps (only for mg)", postSmooth);
    cmd.addInt   ("c", "MG.NumCycles",          "Number of multi-grid cycles for coarse-grid correction, i.e., 1=V, 2=W cycle", cycles);
    cmd.addInt   ("y", "PreconderType",         "0=Direct, 1=FD in time and direct in space, 2=FD in time and multigrid in space", type);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    bool ok = true;
    if (geoIdx<1 || geoIdx>2 ) { gsInfo << "Unfeasible choice for --Geometry (-g).\n"; ok=false; }
    if (obs   <0 || obs   >2 ) { gsInfo << "Unfeasible choice for --ObservationType (-o).\n"; ok=false; }
    if (type  <0 || type  >2 ) { gsInfo << "Unfeasible choice for --PreconderType (-y).\n"; ok=false; }
    if (!ok) return -1;

    gsInfo << "Run parabolic_oc_example with options:\n" << cmd << std::endl;

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    GISMO_ENSURE (geoIdx==1 || geoIdx==2, "Unknown geometry.");
    gsMultiPatch<> mp( geoIdx == 1 ?
        *gsNurbsCreator<>::BSplineRectangle() :
        *gsNurbsCreator<>::BSplineQuarterAnnulus()
    );
    gsMultiPatch<> ti(*gsNurbsCreator<>::BSplineUnitInterval(2));

    gsInfo << "done.\n";

    /************** Define boundary and initial conditions **************/

    gsInfo << "Define boundary and initial conditions... " << std::flush;

    gsFunctionExpr<> zero( "0", mp.geoDim() );

    gsFunctionExpr<> y0( "if(((x-3/2*cos(3*pi/8))^2+(y-3/2*sin(3*pi/8))^2)<0.04,1,0) + "
                         "if(((x-3/2*cos(2*pi/8))^2+(y-3/2*sin(2*pi/8))^2)<0.04,1,0) + "
                         "if(((x-3/2*cos(1*pi/8))^2+(y-3/2*sin(1*pi/8))^2)<0.04,1,0) ",1 );

    gsBoundaryConditions<> bcY;
    gsBoundaryConditions<> bcU;
    for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
    {
        bcY.addCondition( *it, condition_type::dirichlet, zero );
    }

    gsBoundaryConditions<> icY;
    gsBoundaryConditions<> icU;
    for (gsMultiPatch<>::const_biterator it = ti.bBegin(); it < ti.bEnd(); ++it)
    {

        if (it->side() == boundary::left)
        {
            gsInfo << "Set ic for " << *it << "\n";
            icY.addCondition( *it, condition_type::dirichlet, y0 );
        }
        else
        {
            gsInfo << "Set no condition for " << *it << "\n";
        }
    }

    gsInfo << "done.\n";

    /************ Setup bases and adjust degree *************/

    //! [Define Basis]
    gsMultiBasis<> mbY(mp);
    gsMultiBasis<> mbU(mp);
    gsMultiBasis<> tiY(ti);
    gsMultiBasis<> tiU(ti);
    //! [Define Basis]

    gsInfo << "Setup bases and adjust degree... " << std::flush;

    //! [Set degree and refine]
    for ( size_t i = 0; i < mbY.nBases(); ++ i )
    {
        mbY[i].setDegreePreservingMultiplicity(degree);
        mbU[i].setDegreePreservingMultiplicity(degree);
    }

    for ( index_t i = 0; i < refinementsX; ++i )
    {
        mbY.uniformRefine();
        mbU.uniformRefine();
    }

    for ( size_t i = 0; i < mbU.nBases(); ++ i )
        mbU[i].reduceContinuity(2);


    for ( size_t i = 0; i < tiY.nBases(); ++ i )
    {
        tiY[i].setDegreePreservingMultiplicity(degree);
        tiU[i].setDegreePreservingMultiplicity(degree);
    }

    for ( index_t i = 0; i < refinementsT; ++i )
    {
        tiY.uniformRefine();
        tiU.uniformRefine();
    }

    for ( size_t i = 0; i < tiU.nBases(); ++ i )
        tiU[i].reduceContinuity(1);

    gsInfo << "done.\n";

    /********* Setup assembler and assemble matrix **********/

    gsInfo << "Setup assembler and assemble matrices... " << std::flush;

    gsSparseMatrix<> space_massU;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        // We set up the assembler
        gsExprAssembler<> assembler(1,1);
        // Elements used for numerical integration
        assembler.setIntegrationElements(mbU);
        gsExprEvaluator<> ev(assembler);
        // Set the geometry map
        geometryMap G = assembler.getMap(mp);
        // Set the discretization space
        space UU = assembler.getSpace(mbU,1,0);
        // Incorporate Dirichlet BC
        bcU.setGeoMap(mp);
        UU.setup(bcU, dirichlet::interpolation, 0);
        // Initialize the system
        assembler.initSystem();
        assembler.assemble( UU * UU.tr() * meas(G) );
        space_massU = assembler.matrix();
    }
    gsSparseMatrix<> space_massY;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        // We set up the assembler
        gsExprAssembler<> assembler(1,1);
        // Elements used for numerical integration
        assembler.setIntegrationElements(mbY);
        gsExprEvaluator<> ev(assembler);
        // Set the geometry map
        geometryMap G = assembler.getMap(mp);
        // Set the discretization space
        space YY = assembler.getSpace(mbY,1,0);
        // Incorporate Dirichlet BC
        bcY.setGeoMap(mp);
        YY.setup(bcY, dirichlet::interpolation, 0);
        // Initialize the system
        assembler.initSystem();
        assembler.assemble( YY * YY.tr() * meas(G) );
        space_massY = assembler.matrix();
    }
    gsSparseMatrix<> space_massYU;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        // We set up the assembler
        gsExprAssembler<> assembler(2,2);
        // Elements used for numerical integration
        assembler.setIntegrationElements(mbY);
        gsExprEvaluator<> ev(assembler);
        // Set the geometry map
        geometryMap G = assembler.getMap(mp);
        // Set the discretization space
        space YY = assembler.getSpace(mbY,1,0);
        space UU = assembler.getSpace(mbU,1,1);
        // Incorporate Dirichlet BC
        bcY.setGeoMap(mp);
        YY.setup(bcY, dirichlet::interpolation, 0);
        bcU.setGeoMap(mp);
        UU.setup(bcU, dirichlet::interpolation, 0);
        // Initialize the system
        assembler.initSystem();
        assembler.assemble( YY * UU.tr() * meas(G) );
        space_massYU = assembler.matrix().block( 0,             space_massY.rows(),
                                           space_massY.rows(),   space_massU.rows() );
    }

    gsSparseMatrix<> space_laplaceYU;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;

        // We set up the assembler
        gsExprAssembler<> assembler(2,2);

        // Elements used for numerical integration
        assembler.setIntegrationElements(mbY);
        gsExprEvaluator<> ev(assembler);

        // Set the geometry map
        geometryMap G = assembler.getMap(mp);

        // Set the discretization space
        space YY = assembler.getSpace(mbY,1,0);
        space UU = assembler.getSpace(mbU,1,1);

        // Incorporate Dirichlet BC
        bcY.setGeoMap(mp);
        YY.setup(bcY, dirichlet::interpolation, 0);
        bcU.setGeoMap(mp);
        UU.setup(bcU, dirichlet::interpolation, 0);
        //PP.setup(bcU, dirichlet::interpolation, 0);

        // Initialize the system
        assembler.initSystem();

        assembler.assemble( -ilapl(YY,G) * UU.tr() * meas(G) );

        space_laplaceYU = assembler.matrix().block( 0,             space_massY.rows(),
                                           space_massY.rows(),   space_massU.rows() );
        //rhs    = assembler.rhs();
    }
    gsSparseMatrix<> space_biharmY;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        // We set up the assembler
        gsExprAssembler<> assembler(1,1);
        // Elements used for numerical integration
        assembler.setIntegrationElements(mbY);
        gsExprEvaluator<> ev(assembler);
        // Set the geometry map
        geometryMap G = assembler.getMap(mp);
        // Set the discretization space
        space YY = assembler.getSpace(mbY,1,0);
        // Incorporate Dirichlet BC
        bcY.setGeoMap(mp);
        YY.setup(bcY, dirichlet::interpolation, 0);
        // Initialize the system
        assembler.initSystem();

        assembler.assemble( ilapl(YY,G) * ilapl(YY,G).tr() * meas(G) );
        space_biharmY = assembler.matrix();
    }
    gsSparseMatrix<> time_massU;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        // We set up the assembler
        gsExprAssembler<> assembler(1,1);
        // Elements used for numerical integration
        assembler.setIntegrationElements(tiU);
        gsExprEvaluator<> ev(assembler);
        // Set the geometry map
        geometryMap G = assembler.getMap(ti);
        // Set the discretization space
        space UU = assembler.getSpace(tiU,1,0);
        // Incorporate Dirichlet BC
        icU.setGeoMap(ti);
        UU.setup(icU, dirichlet::interpolation, 0);
        // Initialize the system
        assembler.initSystem();
        assembler.assemble( UU * UU.tr() * meas(G) );
        time_massU = assembler.matrix();
    }

    gsSparseMatrix<> time_massY;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        // We set up the assembler
        gsExprAssembler<> assembler(1,1);
        // Elements used for numerical integration
        assembler.setIntegrationElements(tiY);
        gsExprEvaluator<> ev(assembler);
        // Set the geometry map
        geometryMap G = assembler.getMap(ti);
        // Set the discretization space
        space YY = assembler.getSpace(tiY,1,0);
        // Incorporate Dirichlet BC
        icY.setGeoMap(ti);
        YY.setup(icY, dirichlet::interpolation, 0);
        // Initialize the system
        assembler.initSystem();

        assembler.assemble( YY * YY.tr() * meas(G) );
        time_massY = assembler.matrix();
    }

    gsSparseMatrix<> time_massYo;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        // We set up the assembler
        gsExprAssembler<> assembler(1,1);
        // Elements used for numerical integration
        assembler.setIntegrationElements(tiY);
        gsExprEvaluator<> ev(assembler);
        // Set the geometry map
        geometryMap G = assembler.getMap(ti);
        // Set the discretization space
        space YY = assembler.getSpace(tiY,1,0);
        // Incorporate Dirichlet BC
        icY.setGeoMap(ti);
        YY.setup(icY, dirichlet::interpolation, 0);
        // Initialize the system
        assembler.initSystem();

        const char* obstypes[] = {"0",
            "if(x<1/16.,1,0) + if((x>4/16.)&(x<5/16.),1,0) + if((x>12/16.)&(x<13/16.),1,0) + if(x>15/16.,1,0) ",
            "1"
        };

        gsFunctionExpr<> limobs( obstypes[obs], 1);

        variable vlimobs = assembler.getCoeff(limobs);
        assembler.assemble( YY * vlimobs.asDiag() * YY.tr() * meas(G) );
        time_massYo = assembler.matrix();
    }

    gsSparseMatrix<> time_massYU;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        // We set up the assembler
        gsExprAssembler<> assembler(2,2);
        // Elements used for numerical integration
        assembler.setIntegrationElements(tiY);
        gsExprEvaluator<> ev(assembler);
        // Set the geometry map
        geometryMap G = assembler.getMap(ti);
        // Set the discretization space
        space YY = assembler.getSpace(tiY,1,0);
        space UU = assembler.getSpace(tiU,1,1);
        // Incorporate Dirichlet BC
        icY.setGeoMap(ti);
        YY.setup(icY, dirichlet::interpolation, 0);
        icU.setGeoMap(ti);
        UU.setup(icU, dirichlet::interpolation, 0);
        // Initialize the system
        assembler.initSystem();

        assembler.assemble( YY * UU.tr() * meas(G) );
        time_massYU = assembler.matrix().block( 0,                   time_massY.rows(),
                                                time_massY.rows(),   time_massU.rows() );
    }

    gsSparseMatrix<> time_gradYU;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        // We set up the assembler
        gsExprAssembler<> assembler(2,2);
        // Elements used for numerical integration
        assembler.setIntegrationElements(tiY);
        gsExprEvaluator<> ev(assembler);
        // Set the geometry map
        geometryMap G = assembler.getMap(ti);
        // Set the discretization space
        space YY = assembler.getSpace(tiY,1,0);
        space UU = assembler.getSpace(tiU,1,1);
        // Incorporate Dirichlet BC
        icY.setGeoMap(ti);
        YY.setup(icY, dirichlet::interpolation, 0);
        icU.setGeoMap(ti);
        UU.setup(icU, dirichlet::interpolation, 0);
        // Initialize the system
        assembler.initSystem();

        assembler.assemble( igrad(YY,G)[0] * UU.tr() * meas(G) );
        time_gradYU = assembler.matrix().block( 0,                   time_massY.rows(),
                                                time_massY.rows(),   time_massU.rows() );
    }

    gsSparseMatrix<> time_stiffY;
    {
        // The usual stuff for the expression assembler
        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::variable    variable;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;
        // We set up the assembler
        gsExprAssembler<> assembler(1,1);
        // Elements used for numerical integration
        assembler.setIntegrationElements(tiY);
        gsExprEvaluator<> ev(assembler);
        // Set the geometry map
        geometryMap G = assembler.getMap(ti);
        // Set the discretization space
        space YY = assembler.getSpace(tiY,1,0);
        // Incorporate Dirichlet BC
        icY.setGeoMap(ti);
        YY.setup(icY, dirichlet::interpolation, 0);
        // Initialize the system
        assembler.initSystem();
        assembler.assemble( igrad(YY,G) * igrad(YY,G).tr() * meas(G) );
        time_stiffY = assembler.matrix();
    }
    gsInfo << "done.\n";
    gsInfo << "Y has " << space_massY.rows()*time_massY.rows() << " dofs "
        "(" << space_massY.rows() << "x" << space_massY.cols() << ").\n";
    gsInfo << "U has " << space_massU.rows()*time_massU.rows() << " dofs "
        "(" << space_massU.rows() << "x" << space_massU.cols() << ").\n";

    gsLinearOperator<>::Ptr Lh
        = gsSumOp<>::make(
        gsKroneckerOp<>::make( makeMatrixOp(time_gradYU), makeMatrixOp(space_massYU) ),
        gsScaledOp<>::make( gsKroneckerOp<>::make( makeMatrixOp(time_massYU), makeMatrixOp(space_laplaceYU) ), kappa )
        );
    gsLinearOperator<>::Ptr LhT
        = gsSumOp<>::make(
        gsKroneckerOp<>::make( makeMatrixOp(time_gradYU.transpose()), makeMatrixOp(space_massYU.transpose()) ),
        gsScaledOp<>::make( gsKroneckerOp<>::make( makeMatrixOp(time_massYU.transpose()), makeMatrixOp(space_laplaceYU.transpose()) ), kappa )
        );
    gsLinearOperator<>::Ptr Mo
        = gsKroneckerOp<>::make(makeMatrixOp(time_massYo), makeMatrixOp(space_massY));
    gsLinearOperator<>::Ptr Mu
        = gsKroneckerOp<>::make( makeMatrixOp(time_massU), makeMatrixOp(space_massU) );
    gsLinearOperator<>::Ptr MuInv
        = gsKroneckerOp<>::make( makeSparseCholeskySolver(time_massU), makeSparseCholeskySolver(space_massU) );

    gsLinearOperator<>::Ptr schur = gsSumOp<>::make( Mo, gsScaledOp<>::make( gsProductOp<>::make( LhT, MuInv, Lh ), alpha ) );

    gsLinearOperator<>::Ptr schurPreconder;
    if (type==0)
        schurPreconder = makeSparseCholeskySolver(gsSparseMatrix<>(gsSparseMatrix<>(time_massYo+alpha*time_stiffY).kron( space_massY ) + (kappa*kappa*alpha)*time_massY.kron(space_biharmY)));
    else if (type==1)
        schurPreconder = fd( time_massYo+alpha*time_stiffY, space_massY, (kappa*kappa*alpha)*time_massY, space_biharmY, mkSparseCholeskySolver);
    else // if (type==2)
        schurPreconder = fd( time_massYo+alpha*time_stiffY, space_massY, (kappa*kappa*alpha)*time_massY, space_biharmY, mkMGSolver{mbY, bcY, time_massY.rows(), cmd.getGroup("MG")});

    //////

    index_t iter1;
    real_t cond1;

    gsInfo << "Setup cg solver for Schur complement and solve... " << std::flush;
    {
        gsMatrix<> x;
        x.setRandom( schur->rows(), 1 );
        gsMatrix<> rhs;
        rhs.setRandom( schur->rows(), 1 );
        gsMatrix<> errorHistory;
        gsConjugateGradient<> solver( schur, schurPreconder );
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

    ///////

    index_t iter2;
    gsInfo << "\nSetup minres solver for 2x2 system and solve... " << std::flush;
    {

        gsBlockOp<>::Ptr saddle = gsBlockOp<>::make(2,2);
        saddle->addOperator(0,0, Mo);
        saddle->addOperator(0,1, Lh);
        saddle->addOperator(1,0, LhT);
        saddle->addOperator(1,1, gsScaledOp<>::make(Mu,-1/alpha));

        gsBlockOp<>::Ptr saddlePreconder = gsBlockOp<>::make(2,2);
        saddlePreconder->addOperator(0,0, schurPreconder);
        saddlePreconder->addOperator(1,1, gsScaledOp<>::make(MuInv,1/alpha));

        gsMatrix<> x;
        x.setRandom( saddle->rows(), 1 );
        gsMatrix<> rhs;
        rhs.setRandom( saddle->rows(), 1 );
        gsMatrix<> errorHistory;
        gsMinimalResidual<> solver( saddle, saddlePreconder );
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

    }

    ///////

    index_t iter3;
    gsInfo << "\nSetup minres solver for 3x3 system and solve... " << std::flush;
    {

        gsBlockOp<>::Ptr saddle = gsBlockOp<>::make(3,3);
        saddle->addOperator(0,0, Mo);
        saddle->addOperator(0,2, Lh);
        saddle->addOperator(2,0, LhT);
        saddle->addOperator(1,1, gsScaledOp<>::make(Mu,alpha));
        saddle->addOperator(1,2, gsScaledOp<>::make(Mu,1));
        saddle->addOperator(2,1, gsScaledOp<>::make(Mu,1));

        gsBlockOp<>::Ptr saddlePreconder = gsBlockOp<>::make(3,3);
        saddlePreconder->addOperator(0,0, schurPreconder);
        saddlePreconder->addOperator(1,1, gsScaledOp<>::make(MuInv,1/alpha));
        saddlePreconder->addOperator(2,2, gsScaledOp<>::make(MuInv,alpha));

        gsMatrix<> x;
        x.setRandom( saddle->rows(), 1 );
        gsMatrix<> rhs;
        rhs.setRandom( saddle->rows(), 1 );
        gsMatrix<> errorHistory;
        gsMinimalResidual<> solver( saddle, saddlePreconder );
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

    }


    if (!out.empty())
    {
        const bool exists = gsFileManager::fileExists(out);
        std::ofstream outfile;
        outfile.open(out.c_str(), std::ios_base::app);

        if (!exists)
            outfile << "parabolic_oc_example\t"
                "geoIdx\t"
                "obs\t"
                "alpha\t"
                "kappa\t"
                "refinementsX\t"
                "refinementsT\t"
                "degree\t"
                "type\t"
                "iter1\t"
                "cond1\t"
                "iter2\t"
                "iter3\t"
                "Y_spaceDofs\t"
                "Y_timeDofs\t"
                "U_spaceDofs\t"
                "U_timeDofs\n";

        outfile << "parabolic_oc_example\t"
            << geoIdx << "\t"
            << obs << "\t"
            << alpha << "\t"
            << kappa << "\t"
            << refinementsX << "\t"
            << refinementsT << "\t"
            << degree << "\t";
        if (type==0)
            outfile << "direct\t";
        else if (type==1)
            outfile << "fd-d\t";
        else if (type==2)
            outfile << "fd-mg-" << cycles << "-" << preSmooth << "-" << postSmooth << "\t";
        outfile
            << iter1 << "\t"
            << cond1 << "\t"
            << iter2 << "\t"
            << iter3 << "\t"
            << space_massY.rows() << "\t"
            << time_massY.rows() << "\t"
            << space_massU.rows() << "\t"
            << time_massU.rows() << "\n";
    }
    return 0;

}
