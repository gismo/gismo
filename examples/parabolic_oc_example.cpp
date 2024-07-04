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

#define debugMat(m) gsInfo<<#m<<":"<<(m).rows()<<"x"<<(m).cols()<<"\n"
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

gsLinearOperator<>::Ptr mkSparseCholeskySolver(const gsSparseMatrix<>& m) { return makeSparseCholeskySolver(m); }

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

int main(int argc, char *argv[])
{
    /************** Define command line options *************/

    index_t geoIdx = 1;
    index_t refinementsX = 1;
    index_t refinementsT = 1;
    index_t degree = 3;
    //real_t timeHorizon = 1.;
    real_t tolerance = 1.e-6;
    real_t alpha = 1.;
    real_t kappa = 0.01;
    index_t obs = 1;
    index_t maxIterations = 100;
    index_t pcType = 1;
    index_t preSmooth = 1;
    index_t postSmooth = 1;
    index_t cycles = 1;
    std::string out;
    bool plot = false;
    bool randRhs = false;
    index_t desiredStateIdx = 0;

    gsCmdLine cmd("parabolic_oc_example");
    cmd.addInt   ("g", "Geometry",              "0=Rectangle, 1=Quarter Annulus", geoIdx);
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
    cmd.addInt   ("y", "PreconderType",         "0=Direct, 1=FD in time and direct in space, 2=FD in time and multigrid in space", pcType);
    cmd.addInt   ("d", "DesiredState",          "0=const (like i.c), 1=follows parabolic proces", desiredStateIdx);
    cmd.addString("",  "out",                   "Write solution and used options to file", out);
    cmd.addSwitch(     "plot",                  "Plot the result with Paraview", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    bool ok = true;

    const char* obstypes[] = {
        "0",
        "if(x<1/16.,1,0) + if((x>4/16.)&(x<5/16.),1,0) + if((x>12/16.)&(x<13/16.),1,0) + if(x>15/16.,1,0) ",
        "1",
        "if(x<1/16.,1,0) + if((x>12/16.)&(x<13/16.),1,0) ",
    };

    gsMultiPatch<> geos[] = {
        *gsNurbsCreator<>::BSplineRectangle(),
        *gsNurbsCreator<>::BSplineQuarterAnnulus()
    };

    if (geoIdx<0          || geoIdx>(int)(util::size(geos)-1)    ) { gsInfo << "Unfeasible choice for --Geometry (-g).\n";        ok=false; }
    if (obs   <0          || obs   >(int)(util::size(obstypes)-1)) { gsInfo << "Unfeasible choice for --ObservationType (-o).\n"; ok=false; }
    if (pcType<0          || pcType>2                            ) { gsInfo << "Unfeasible choice for --PreconderType (-y).\n";   ok=false; }
    if (desiredStateIdx<0 || desiredStateIdx>1                   ) { gsInfo << "Unfeasible choice for --DesiredState (-d).\n";    ok=false; }
    if (!ok) return -1;

    gsInfo << "Run parabolic_oc_example with options:\n" << cmd << std::endl;

    /******************* Define geometry ********************/

    gsInfo << "Define geometry... " << std::flush;

    const gsMultiPatch<>& mp = geos[geoIdx];
    gsMultiPatch<> ti(*gsNurbsCreator<>::BSplineUnitInterval(2));

    gsInfo << "done.\n";

    /************** Define boundary and initial conditions **************/

    gsInfo << "Define boundary and initial conditions... " << std::flush;

    gsConstantFunction<> zero( 0, mp.geoDim() );

    gsFunctionExpr<> y0( "if(((x-3/2*cos(3*pi/8))^2+(y-3/2*sin(3*pi/8))^2)<0.04,1,0) + "
                         "if(((x-3/2*cos(2*pi/8))^2+(y-3/2*sin(2*pi/8))^2)<0.04,1,0) + "
                         "if(((x-3/2*cos(1*pi/8))^2+(y-3/2*sin(1*pi/8))^2)<0.04,1,0) ", mp.geoDim() );

    gsBoundaryConditions<> bcY;
    gsBoundaryConditions<> bcU; // Space U does not have boundary conditions
    for (gsMultiPatch<>::const_biterator it = mp.bBegin(); it < mp.bEnd(); ++it)
    {
        bcY.addCondition( *it, condition_type::dirichlet, zero );
    }

    gsBoundaryConditions<> icY; // Have an initial conditoion (side::left). Do no set conditions here, but handle them later.
    gsBoundaryConditions<> icU; // Space U does not have initial conditions

    gsInfo << "done.\n";

    /************ Setup bases and adjust degree *************/
    gsMultiBasis<> mbY(mp);
    gsMultiBasis<> mbU(mp);
    gsMultiBasis<> tiY(ti);
    gsMultiBasis<> tiU(ti);

    gsInfo << "Setup bases and adjust degree... " << std::flush;

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

    // Assemble space matrices
    gsSparseMatrix<> space_massY, space_embeddingY;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(mbY);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space YY = assembler.getSpace(mbY,1,0);
        bcY.setGeoMap(mp);
        YY.setup(bcY, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( YY * YY.tr() * meas(G) );
        space_massY = assembler.matrix();
        space_embeddingY = embeddingMatrix( YY.mapper().inverseAsVector(), mbY.size(), space_massY.rows() );
    }
    gsSparseMatrix<> space_massU;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(mbU);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space UU = assembler.getSpace(mbU,1,0);
        bcU.setGeoMap(mp);
        UU.setup(bcU, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( UU * UU.tr() * meas(G) );
        space_massU = assembler.matrix();
    }
    gsSparseMatrix<> space_massYU;
    {
        gsExprAssembler<> assembler(2,2);
        assembler.setIntegrationElements(mbY);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space YY = assembler.getSpace(mbY,1,0);
        gsExprAssembler<>::space UU = assembler.getSpace(mbU,1,1);
        bcY.setGeoMap(mp);
        YY.setup(bcY, dirichlet::interpolation, 0);
        bcU.setGeoMap(mp);
        UU.setup(bcU, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( YY * UU.tr() * meas(G) );
        space_massYU = assembler.matrix().block( 0,                    space_massY.rows(),
                                                 space_massY.rows(),   space_massU.rows() );
    }
    gsSparseMatrix<> space_laplaceYU;
    {
        gsExprAssembler<> assembler(2,2);
        assembler.setIntegrationElements(mbY);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space YY = assembler.getSpace(mbY,1,0);
        gsExprAssembler<>::space UU = assembler.getSpace(mbU,1,1);
        bcY.setGeoMap(mp);
        YY.setup(bcY, dirichlet::interpolation, 0);
        bcU.setGeoMap(mp);
        UU.setup(bcU, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( -ilapl(YY,G) * UU.tr() * meas(G) );
        space_laplaceYU = assembler.matrix().block( 0,                    space_massY.rows(),
                                                    space_massY.rows(),   space_massU.rows() );
    }
    gsSparseMatrix<> space_biharmY;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(mbY);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space YY = assembler.getSpace(mbY,1,0);
        bcY.setGeoMap(mp);
        YY.setup(bcY, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( ilapl(YY,G) * ilapl(YY,G).tr() * meas(G) );
        space_biharmY = assembler.matrix();
    }

    // Assemble in time
    gsSparseMatrix<> time_massY;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(tiY);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(ti);
        gsExprAssembler<>::space YY = assembler.getSpace(tiY,1,0);
        icY.setGeoMap(ti);
        YY.setup(icY, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( YY * YY.tr() * meas(G) );
        time_massY = assembler.matrix();
    }
    gsSparseMatrix<> time_massU;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(tiU);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(ti);
        gsExprAssembler<>::space UU = assembler.getSpace(tiU,1,0);
        icU.setGeoMap(ti);
        UU.setup(icU, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( UU * UU.tr() * meas(G) );
        time_massU = assembler.matrix();
    }
    gsSparseMatrix<> time_massYo;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(tiY);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(ti);
        gsExprAssembler<>::space YY = assembler.getSpace(tiY,1,0);
        icY.setGeoMap(ti);
        YY.setup(icY, dirichlet::interpolation, 0);
        assembler.initSystem();
        gsFunctionExpr<> limobs( obstypes[obs], 1);
        gsExprAssembler<>::variable vlimobs = assembler.getCoeff(limobs);
        assembler.assemble( YY * vlimobs.asDiag() * YY.tr() * meas(G) );
        time_massYo = assembler.matrix();
    }
    gsSparseMatrix<> time_massYU;
    {
        gsExprAssembler<> assembler(2,2);
        assembler.setIntegrationElements(tiY);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(ti);
        gsExprAssembler<>::space YY = assembler.getSpace(tiY,1,0);
        gsExprAssembler<>::space UU = assembler.getSpace(tiU,1,1);
        icY.setGeoMap(ti);
        YY.setup(icY, dirichlet::interpolation, 0);
        icU.setGeoMap(ti);
        UU.setup(icU, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( YY * UU.tr() * meas(G) );
        time_massYU = assembler.matrix().block( 0,                   time_massY.rows(),
                                                time_massY.rows(),   time_massU.rows() );
    }
    gsSparseMatrix<> time_gradYU;
    {
        gsExprAssembler<> assembler(2,2);
        assembler.setIntegrationElements(tiY);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(ti);
        gsExprAssembler<>::space YY = assembler.getSpace(tiY,1,0);
        gsExprAssembler<>::space UU = assembler.getSpace(tiU,1,1);
        icY.setGeoMap(ti);
        YY.setup(icY, dirichlet::interpolation, 0);
        icU.setGeoMap(ti);
        UU.setup(icU, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( igrad(YY,G)[0] * UU.tr() * meas(G) );
        time_gradYU = assembler.matrix().block( 0,                   time_massY.rows(),
                                                time_massY.rows(),   time_massU.rows() );
    }
    gsSparseMatrix<> time_stiffY;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(tiY);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(ti);
        gsExprAssembler<>::space YY = assembler.getSpace(tiY,1,0);
        icY.setGeoMap(ti);
        YY.setup(icY, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( igrad(YY,G) * igrad(YY,G).tr() * meas(G) );
        time_stiffY = assembler.matrix();
    }
    // Construct embedding for time and initial conditon vector manually
    gsMatrix<> initialConditionY;
    gsSparseMatrix<> time_embeddingY, time_embeddingIC;
    {
        gsExprAssembler<> assembler(1,1);
        assembler.setIntegrationElements(mbY);
        gsExprEvaluator<> ev(assembler);
        gsExprAssembler<>::geometryMap G = assembler.getMap(mp);
        gsExprAssembler<>::space YY = assembler.getSpace(mbY,1,0);
        auto ff = assembler.getCoeff(y0, G);
        bcY.setGeoMap(mp);
        YY.setup(bcY, dirichlet::interpolation, 0);
        assembler.initSystem();
        assembler.assemble( YY * ff * meas(G) );
        makeSparseCholeskySolver(space_massY)->apply(assembler.rhs(), initialConditionY);
        const index_t sz = time_massY.rows();
        time_embeddingY .resize(sz, sz-1);
        for (index_t i=0; i<sz-1; ++i)
            time_embeddingY(i+1,i) = 1;
        time_embeddingIC.resize(sz, 1);
        time_embeddingIC(0,0)=1;
    }
    // Handle initial conditions for Y
    gsSparseMatrix<> time_massYo_ic, time_massYU_ic, time_gradYU_ic;
    {
        // Matrices required to incorporate i.c. into rhs
        time_massYo_ic = time_embeddingIC.transpose() * time_massYo  * time_embeddingY;
        time_massYU_ic = time_embeddingIC.transpose() * time_massYU;
        time_gradYU_ic = time_embeddingIC.transpose() * time_gradYU;

        // Incorporate i.c. in matrices, i.e., remove first row+col
        time_massY     = time_embeddingY.transpose()  * time_massY   * time_embeddingY;
        time_massYo    = time_embeddingY.transpose()  * time_massYo  * time_embeddingY;
        time_massYU    = time_embeddingY.transpose()  * time_massYU;
        time_gradYU    = time_embeddingY.transpose()  * time_gradYU;
        time_stiffY    = time_embeddingY.transpose()  * time_stiffY  * time_embeddingY;
    }

    gsInfo << "done.\n";

    {
    printMat(space_massY);
    printMat(space_embeddingY);
    printMat(space_massU);
    printMat(space_massYU);
    printMat(space_laplaceYU);
    printMat(space_biharmY);
    printMat(time_massY);
    printMat(time_massU);
    printMat(time_massYo);
    printMat(time_massYU);
    printMat(time_gradYU);
    printMat(time_stiffY);
    printMat(time_embeddingY);
    printMat(time_embeddingIC);
    printMat(time_massYo_ic);
    printMat(time_massYU_ic);
    printMat(time_gradYU_ic);
    }



    gsInfo << "Y has " << space_massY.rows()*time_massY.rows() << " dofs "
        "(" << space_massY.rows() << "x" << time_massY.cols() << ").\n";
    gsInfo << "U has " << space_massU.rows()*time_massU.rows() << " dofs "
        "(" << space_massU.rows() << "x" << time_massU.cols() << ").\n";

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
    if (pcType==0)
        schurPreconder = makeSparseCholeskySolver(gsSparseMatrix<>(gsSparseMatrix<>(time_massYo+alpha*time_stiffY).kron( space_massY ) + (kappa*kappa*alpha)*time_massY.kron(space_biharmY)));
    else if (pcType==1)
        schurPreconder = fastDiagnonalization( time_massYo+alpha*time_stiffY, space_massY, (kappa*kappa*alpha)*time_massY, space_biharmY, gsSolverOp<gsSparseSolver<>::SimplicialLDLT>::make);
    else // if (pcType==2)
        schurPreconder = fastDiagnonalization( time_massYo+alpha*time_stiffY, space_massY, (kappa*kappa*alpha)*time_massY, space_biharmY, makeMultiGridSolver{mbY, bcY, time_massY.rows(), cmd.getGroup("MG")});


    // For initial conditions (as above)
    gsSparseMatrix<> Mo_ic = time_massYo_ic.kron(space_massY);
    gsSparseMatrix<> Lh_ic = time_gradYU_ic.kron(space_massYU) + kappa * time_massYU_ic.kron(space_laplaceYU);

    //////
    gsMatrix<> desiredState, Mo_desiredState;
    if (!randRhs)
    {
        if (desiredStateIdx==0)
        {
            // constant desired state is constant
            desiredState.resize(space_massY.rows()*time_massY.rows(), 1);
            for (index_t i=0; i<time_massY.rows(); ++i)
                desiredState.block(i*space_massY.rows(),0,space_massY.rows(),1) = initialConditionY;
        }
        else //if (desiredStateIdx==1)
        {
            // desired state follows parabolic process with 0 rhs
            gsLinearOperator<>::Ptr fwPreconder;
            if (pcType==0)
                fwPreconder = makeSparseCholeskySolver(gsSparseMatrix<>(gsSparseMatrix<>(time_stiffY).kron( space_massY ) + (kappa*kappa)*time_massY.kron(space_biharmY)));
            else if (pcType==1)
                fwPreconder = fastDiagnonalization( time_stiffY, space_massY, (kappa*kappa)*time_massY, space_biharmY, gsSolverOp<gsSparseSolver<>::SimplicialLDLT>::make);
            else // if (pcType==2)
                fwPreconder = fastDiagnonalization( time_stiffY, space_massY, (kappa*kappa)*time_massY, space_biharmY, makeMultiGridSolver{mbY, bcY, time_massY.rows(), cmd.getGroup("MG")});
            gsInfo << "Compute parabolic desired state..." << std::flush;
            gsMatrix<> rhs;
            Lh->apply(-Lh_ic.transpose()*initialConditionY, rhs);
            gsGMRes<> gmres(gsProductOp<>::Ptr(gsProductOp<>::make(LhT,Lh)),fwPreconder);
            real_t tol = 1e-10;
            gmres.setTolerance(tol);
            gsMatrix<> errorHistory;
            gmres.solveDetailed(rhs, desiredState, errorHistory);
            gsInfo << "done.\n";

            index_t iter1 = errorHistory.rows()-1;
            const bool success = errorHistory(iter1,0) < tol;
            if (success)
                gsInfo << "Reached desired tolerance after " << iter1 << " iterations:\n";
            else
                gsInfo << "Did not reach desired tolerance after " << iter1 << " iterations:\n";

            if (errorHistory.rows() < 20)
                gsInfo << errorHistory.transpose() << "\n\n";
            else
                gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";
        }
        Mo->apply(desiredState,Mo_desiredState);
    }

    //////
    gsMatrix<> solution_Y, solution_U;
    //////

    index_t iter1;
    real_t cond1;

    gsInfo << "Setup cg solver for Schur complement and solve... " << std::flush;
    {
        gsMatrix<> x;
        x.setRandom( schur->rows(), 1 );
        gsMatrix<> rhs;
        rhs.setRandom( schur->rows(), 1 ); // TODO
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
        if (randRhs)
            rhs.setRandom( saddle->rows(), 1 );
        else
        {
            rhs.setZero( saddle->rows(), 1 );
            // desired state
            rhs.topRows(Mo->rows()) = Mo_desiredState;
            // incorporate initial conditions
            //rhs.topRows(Mo->rows()) -= Mo_ic.transpose() * (initialConditionY-initialConditionForDesiredState);
            rhs.bottomRows(Mu->rows()) -= Lh_ic.transpose() * initialConditionY;
        }
        gsMatrix<> errorHistory;
        gsMinimalResidual<> solver( saddle, saddlePreconder );
        solver.setOptions( cmd.getGroup("Solver") ).solveDetailed( rhs, x, errorHistory );

        gsInfo << "done.\n\n";
        solution_Y = x.topRows(Mo->rows());
        solution_U = (-1./alpha) * x.bottomRows(Mu->rows());

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
        if (randRhs)
            rhs.setRandom( saddle->rows(), 1 );
        else
        {
            rhs.setZero( saddle->rows(), 1 );
            // desired state
            rhs.topRows(Mo->rows()) = Mo_desiredState;
            // incorporate initial conditions
            //rhs.topRows(Mo->rows()) -= Mo_ic.transpose() * (initialConditionY-initialConditionForDesiredState);
            rhs.bottomRows(Mu->rows()) -= Lh_ic.transpose() * initialConditionY;
        }

        gsMatrix<> errorHistory;
        gsMinimalResidual<> solver( saddle, saddlePreconder );
        solver.setOptions( cmd.getGroup("Solver") ).solveDetailed( rhs, x, errorHistory );

        gsInfo << "done.\n\n";

        //solution_Y = x.topRows(Mo->rows());
        //solution_U = x.middleRows(Mo->rows(), Mu->rows());

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
                "pcType\t"
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
        if (pcType==0)
            outfile << "direct\t";
        else if (pcType==1)
            outfile << "fd-d\t";
        else //if (pcType==2)
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

    {
        printMat(initialConditionY.transpose());
        gsMatrix<> solution_Y_reshaped = solution_Y.reshaped( space_massY.rows(), time_massY.rows() ).transpose();
        printMat(solution_Y_reshaped);
        gsMatrix<> solution_U_reshaped = solution_U.reshaped( space_massU.rows(), time_massU.rows() ).transpose();
        printMat(solution_U_reshaped);
    }

    if (plot)
    {
        // Plot desired state Y
        {
            gsMatrix<> solution_Y_reshaped = desiredState.reshaped( space_massY.rows(), time_massY.rows() ).transpose();
            gsMatrix<> solution_Y_embedded = time_embeddingY * solution_Y_reshaped * space_embeddingY.transpose();
            solution_Y_embedded.row(0) = (space_embeddingY * initialConditionY).transpose(); // incorporate ic
            gsMatrix<> solution_Y_full = solution_Y_embedded.reshaped( solution_Y_embedded.rows() * solution_Y_embedded.cols(), 1 );
            gsMultiPatch<real_t> mp_geometry( *tensorBasis( gsMultiBasis<>(ti), gsMultiBasis<>(mp) ).makeGeometry( tensorCoefs(ti[0].coefs(),mp[0].coefs() ) ) );
            gsMultiPatch<real_t> mp_solution_Y( *tensorBasis(tiY, mbY).makeGeometry( solution_Y_full ) );
            gsField<real_t> field_Y(mp_geometry, mp_solution_Y, true);
            gsWriteParaview<>(field_Y, "parabolic_oc_example_des_state", 10000);
            gsInfo << "Write solution to \"parabolic_oc_example_des_state.pvd\".\n";
        }


        // Plot state Y
        {
            gsMatrix<> solution_Y_reshaped = solution_Y.reshaped( space_massY.rows(), time_massY.rows() ).transpose();
            gsMatrix<> solution_Y_embedded = time_embeddingY * solution_Y_reshaped * space_embeddingY.transpose();
            solution_Y_embedded.row(0) = (space_embeddingY * initialConditionY).transpose(); // incorporate ic
            gsMatrix<> solution_Y_full = solution_Y_embedded.reshaped( solution_Y_embedded.rows() * solution_Y_embedded.cols(), 1 );
            gsMultiPatch<real_t> mp_geometry( *tensorBasis( gsMultiBasis<>(ti), gsMultiBasis<>(mp) ).makeGeometry( tensorCoefs(ti[0].coefs(),mp[0].coefs() ) ) );
            gsMultiPatch<real_t> mp_solution_Y( *tensorBasis(tiY, mbY).makeGeometry( solution_Y_full ) );
            gsField<real_t> field_Y(mp_geometry, mp_solution_Y, true);
            gsWriteParaview<>(field_Y, "parabolic_oc_example_sol_state", 10000);
            gsInfo << "Write solution to \"parabolic_oc_example_sol_state.pvd\".\n";
        }

        // Plot state Y-desiredState
        {
            gsMatrix<> solution_Y_reshaped = (solution_Y-desiredState).reshaped( space_massY.rows(), time_massY.rows() ).transpose();
            gsMatrix<> solution_Y_embedded = time_embeddingY * solution_Y_reshaped * space_embeddingY.transpose();
            gsMatrix<> solution_Y_full = solution_Y_embedded.reshaped( solution_Y_embedded.rows() * solution_Y_embedded.cols(), 1 );
            gsMultiPatch<real_t> mp_geometry( *tensorBasis( gsMultiBasis<>(ti), gsMultiBasis<>(mp) ).makeGeometry( tensorCoefs(ti[0].coefs(),mp[0].coefs() ) ) );
            gsMultiPatch<real_t> mp_solution_Y( *tensorBasis(tiY, mbY).makeGeometry( solution_Y_full ) );
            gsField<real_t> field_Y(mp_geometry, mp_solution_Y, true);
            gsWriteParaview<>(field_Y, "parabolic_oc_example_diff_state", 10000);
            gsInfo << "Write solution to \"parabolic_oc_example_diff_state.pvd\".\n";
        }

        // Plot state Y for t=0 (ic)
        {
            gsMatrix<> solution_Y_full = space_embeddingY * initialConditionY;
            gsMultiPatch<real_t> mp_solution_Y( *mbY[0].makeGeometry( solution_Y_full ) );
            gsField<real_t> field_Y(mp, mp_solution_Y, true);
            gsWriteParaview<>(field_Y, "parabolic_oc_example_sol_state_0", 10000);
            gsInfo << "Write solution to \"parabolic_oc_example_sol_state_0.pvd\".\n";
        }

        // Plot state Y for t=T
        {
            gsMatrix<> solution_Y_reshaped = solution_Y.reshaped( space_massY.rows(), time_massY.rows() ).transpose();
            gsMatrix<> solution_Y_embedded = time_embeddingY * solution_Y_reshaped * space_embeddingY.transpose();
            gsMatrix<> solution_Y_full = solution_Y_embedded.row(solution_Y_embedded.rows()-1).transpose();
            gsMultiPatch<real_t> mp_solution_Y( *mbY[0].makeGeometry( solution_Y_full ) );
            gsField<real_t> field_Y(mp, mp_solution_Y, true);
            gsWriteParaview<>(field_Y, "parabolic_oc_example_sol_state_N", 10000);
            gsInfo << "Write solution to \"parabolic_oc_example_sol_state_N.pvd\".\n";
        }

        // Plot control U
        {
            gsMatrix<> solution_U_reshaped = solution_U.reshaped( space_massU.rows(), time_massU.rows() ).transpose();
            gsMatrix<> solution_U_full = solution_U_reshaped.reshaped( solution_U_reshaped.rows() * solution_U_reshaped.cols(), 1 );
            gsMultiPatch<real_t> mp_geometry( *tensorBasis( gsMultiBasis<>(ti), gsMultiBasis<>(mp) ).makeGeometry( tensorCoefs(ti[0].coefs(),mp[0].coefs() ) ) );
            gsMultiPatch<real_t> mp_solution_U( *tensorBasis(tiU, mbU).makeGeometry( solution_U_full ) );
            gsField<real_t> field_U(mp_geometry, mp_solution_U, true);
            gsWriteParaview<>(field_U, "parabolic_oc_example_sol_control", 10000);
            gsInfo << "Write solution to \"parabolic_oc_example_sol_control.pvd\".\n";
        }

    }


    return 0;

}
