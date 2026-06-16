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
        if (lv>1)
        {
            transferMatrices.reserve(lv-1);
            index_t nTimeDofs = m.rows() / gh.getTransferMatrices()[lv-2].rows();
            gsSparseMatrix<> id(nTimeDofs,nTimeDofs);
            id.setIdentity();
            // consider using gsKroneckerOp<>
            for (index_t i=0; i<lv-1; ++i)
                transferMatrices.push_back( id.kron(gh.getTransferMatrices()[i]) );
        }

        gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make( m, transferMatrices );
        mg->setOptions(opt);

        for (index_t i = 1; i < mg->numLevels(); ++i)
        {
            mg->setSmoother(i, makeGaussSeidelOp(mg->matrix(i)));
        }

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
                GISMO_ENSURE (tmp.rows()-1==tmp2.rows(),"");
                GISMO_ENSURE (tmp.cols()-1==tmp2.cols(),"");
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
        real_t tau = pow(.5,lv);
        for (index_t i=0; i<lv; ++i)
        {

            const index_t nSpaceDofs = space_stiff.rows();
            const index_t nTimeDofs = accumulatedTransferMatrices[i].cols() / nSpaceDofs;

#if 1
            gsSparseMatrix<> id(nTimeDofs, nTimeDofs); id.setIdentity();
            gsSparseMatrix<> system = id.kron( opt.askReal("Sigma", 1)*sqrt(tau)*space_stiff + sqrt(1/tau)*space_mass );
            gsLinearOperator<>::Ptr systemSolver = makeSolver(system);
            gsSparseMatrix<> multiplierMatrix = id.kron(space_stiff);
            gsLinearOperator<>::Ptr multiplier = makeMatrixOp(multiplierMatrix.moveToPtr());
#else
            gsLinearOperator<>::Ptr idOp = gsIdentityOp<>::make(nTimeDofs);
            gsSparseMatrix<> system = opt.askReal("Sigma", 1)*sqrt(tau)*space_stiff + sqrt(1/tau)*space_mass;
            gsLinearOperator<>::Ptr systemSolver = gsKroneckerOp<>::make(idOp, makeSolver(system));
            gsLinearOperator<>::Ptr multiplier = gsKroneckerOp<>::make(idOp, makeMatrixOp(gsSparseMatrix<>(space_stiff).moveToPtr()));
#endif
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

real_t
trace(const gsSparseMatrix<>& mat)
{
    real_t coef = 0;
    for (index_t i=0; i<mat.rows(); ++i)
        coef += mat(i,i);
    return coef/mat.rows();
}


gsLinearOperator<>::Ptr
mkOpsForMg(gsSparseMatrix<> time_stiff, gsSparseMatrix<> space_stiff, gsSparseMatrix<> time_mass, gsSparseMatrix<> space_mass)
{
    gsLinearOperator<>::Ptr spmass = makeMatrixOp(space_mass.moveToPtr());
    gsLinearOperator<>::Ptr schur = gsProductOp<>::make(spmass, makeSparseCholeskySolver(space_stiff), spmass);
    return gsSumOp<>::make
        (
            gsKroneckerOp<>::make(
                    makeMatrixOp(time_stiff.moveToPtr()),
                    schur
            ),
            gsKroneckerOp<>::make(
                    makeMatrixOp(time_mass.moveToPtr()),
                    makeMatrixOp(space_stiff.moveToPtr())
            )
        );
}

gsLinearOperator<>::Ptr
mkSmsForMg(const gsSparseMatrix<>& time_stiff, const gsSparseMatrix<>& space_stiff, const gsSparseMatrix<>& time_mass, const gsSparseMatrix<>& space_mass)
{
    const real_t hSq = 1./trace(space_stiff);
    const real_t tauSq = 1./trace(time_stiff);
    const real_t factor = hSq/(tauSq*tauSq) + 1./hSq;
    // TODO: Better scaling
    return gsScaledOp<>::make(gsIdentityOp<>::make(time_stiff.rows() * space_stiff.rows()), 1./factor); // 1/(...) since we want to solve...
}

gsLinearOperator<>::Ptr
makeSpaceTimeMultiGridSolver(
    gsSparseMatrix<> time_stiff,
    gsSparseMatrix<> space_stiff,
    gsSparseMatrix<> time_mass,
    gsSparseMatrix<> space_mass,
    const gsMultiBasis<>& mb_time,
    const gsBoundaryConditions<>& bc_time,
    const gsMultiBasis<>& mb_space,
    const gsBoundaryConditions<>& bc_space,
    const gsOptionList& opt)
{

    // TODO: kappa, sigma,...!?

    gsInfo << "Setup space-time multigrid solver... " << std::flush;

    gsOptionList cmd_time;
    cmd_time.addInt( "InterfaceStrategy", "", iFace::conforming      );
    cmd_time.addInt( "DirichletStrategy", "", dirichlet::elimination );
    gsGridHierarchy<> gh_time = gsGridHierarchy<>::buildByCoarsening(mb_time, bc_time, cmd_time, 10000, 4);
    const index_t lv_time = gh_time.getTransferMatrices().size()+1;
    {
        gsInfo << "\n" << lv_time << " levels in time:";
        for (index_t i=lv_time-2; i>-1; --i)
        {
            gsInfo << " " << (gh_time.getTransferMatrices()[i].rows()-1) << "x" << (gh_time.getTransferMatrices()[i].cols()-1);
        }
        gsInfo << "\n";
    }

    gsOptionList cmd_space;
    cmd_space.addInt( "InterfaceStrategy", "", iFace::conforming      );
    cmd_space.addInt( "DirichletStrategy", "", dirichlet::elimination );
    gsGridHierarchy<> gh_space = gsGridHierarchy<>::buildByCoarsening(mb_space, bc_space, cmd_space, 10000, 8);
    const index_t lv_space = gh_space.getTransferMatrices().size()+1;
    {
        gsInfo << lv_space << " levels in space:";
        for (index_t i=lv_space-2; i>-1; --i)
        {
            gsInfo << " " << gh_space.getTransferMatrices()[i].rows() << "x" << gh_space.getTransferMatrices()[i].cols();
        }
        gsInfo << "\n";
    }


    // Provide transfers, ops and smoother_ops for all levels
    gsInfo << "Setup of combined space-time grid hierarchy..." << std::flush;
    std::vector<gsLinearOperator<>::Ptr> prolongations, restrictions, ops, smoother_ops;
    gsLinearOperator<>::Ptr coarseSolver;
    {
        index_t space_dofs = space_mass.rows();
        index_t time_dofs  = time_mass .rows();

        ops         .push_back(mkOpsForMg(time_stiff,space_stiff,time_mass,space_mass));
        smoother_ops.push_back(mkSmsForMg(time_stiff,space_stiff,time_mass,space_mass));

        for (index_t l_time = lv_time-2, l_space = lv_space-2;;)
        {

            const real_t hSq = trace(space_mass)/trace(space_stiff);
            const real_t tauSq = trace(time_mass)/trace(time_stiff);

            bool coarsenInTime = hSq*hSq > tauSq; // TODO: a bit more elaborate...

            gsSparseMatrix<real_t,RowMajor> prolmat;
            bool doBreak;

            gsInfo << (coarsenInTime?"T":"S") << " (tau=" << sqrt(tauSq) << "; h=" << sqrt(hSq) << "); " << std::flush;

            if (coarsenInTime)
            {
                prolmat = gh_time.getTransferMatrices()[l_time];
                prolmat = prolmat.block(1,1, prolmat.rows()-1, prolmat.cols()-1);
                doBreak = l_time==0;
                --l_time;
            }
            else
            {
                prolmat = gh_space.getTransferMatrices()[l_space];
                doBreak = l_space==1; // TODO: the coarset space level has 0 dofs ...
                --l_space;
            }

            gsSparseMatrix<real_t> restmat = prolmat.transpose();
            if (coarsenInTime)
            {
                GISMO_ENSURE (time_dofs == prolmat.rows(), "Dimension missmatch: "<<time_dofs<<"=="<<prolmat.rows());
                time_dofs = prolmat.cols();

                GISMO_ENSURE(restmat.cols() == time_mass.rows(), "Dimension missmatch: "<<restmat.cols()<<"=="<<time_mass.rows());
                time_mass  = restmat * time_mass  * prolmat;
                time_stiff = restmat * time_stiff * prolmat;

                prolongations.push_back(gsKroneckerOp<>::make(
                        makeMatrixOp(prolmat.moveToPtr()),
                        gsIdentityOp<>::make(space_dofs)
                ));
                restrictions .push_back(gsKroneckerOp<>::make(
                        makeMatrixOp(restmat.moveToPtr()),
                        gsIdentityOp<>::make(space_dofs)
                ));
            }
            else //if (!coarsenInTime)
            {
                GISMO_ENSURE (space_dofs == prolmat.rows(), "Dimension missmatch: "<<space_dofs<<"=="<<prolmat.rows());
                space_dofs = prolmat.cols();

                GISMO_ENSURE(restmat.cols() == space_mass.rows(), "Dimension missmatch: "<<restmat.cols()<<"=="<<space_mass.rows());
                space_mass  = restmat * space_mass  * prolmat;
                space_stiff = restmat * space_stiff * prolmat;

                prolongations.push_back(gsKroneckerOp<>::make(
                        gsIdentityOp<>::make(time_dofs),
                        makeMatrixOp(prolmat.moveToPtr())
                ));
                restrictions .push_back(gsKroneckerOp<>::make(
                        gsIdentityOp<>::make(time_dofs),
                        makeMatrixOp(restmat.moveToPtr())
                ));
            }

            ops         .push_back(mkOpsForMg(time_stiff,space_stiff,time_mass,space_mass));
            smoother_ops.push_back(mkSmsForMg(time_stiff,space_stiff,time_mass,space_mass));

            if (doBreak)
                break;
        }
        // Construct coarseSolver
        {
            gsLinearOperator<>::Ptr op = mkOpsForMg(time_stiff,space_stiff,time_mass,space_mass);
            gsMatrix<> mat;
            op->toMatrix(mat);
            gsSparseMatrix<> sm = mat.sparseView(1,1e-8);
            coarseSolver = makeSparseCholeskySolver(sm);
        }
    }

    const index_t lv_total = prolongations.size()+1;
    {
        gsInfo << "\n" << lv_total << " levels in space-time:";
        for (index_t i=0; i<lv_total-1; ++i)
        {
            gsInfo << " " << prolongations[i]->rows() << "x" << prolongations[i]->cols();
        }
        gsInfo << "\n";
    }

    GISMO_ENSURE (ops.size() == prolongations.size()+1, "Dimension missmatch: "<<ops.size()<<"=="<<prolongations.size()+1);
    GISMO_ENSURE (restrictions.size() == prolongations.size(), "Dimension missmatch: "<<restrictions.size()<<"=="<<prolongations.size());


    std::reverse(prolongations.begin(), prolongations.end());
    std::reverse(restrictions.begin(), restrictions.end());
    std::reverse(ops.begin(), ops.end());
    gsMultiGridOp<>::Ptr mg = gsMultiGridOp<>::make( ops, prolongations, restrictions, coarseSolver );
    mg->setOptions(opt);
    for (index_t i=1; i<mg->numLevels(); ++i)
    {
        gsLinearOperator<>::Ptr smop = smoother_ops[smoother_ops.size()-i-1];
        GISMO_ENSURE (smop->rows() == mg->underlyingOp(i)->rows(), "Dimension missmatch: " << smop->rows()<<"=="<<mg->underlyingOp(i)->rows());
        gsPreconditionerOp<>::Ptr smootherOp = gsPreconditionerFromOp<>::make(mg->underlyingOp(i),smop);
        smootherOp->setOptions(opt); //TODO: How to choose damping?
        mg->setSmoother(i, smootherOp);
    }

    return mg;
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
    real_t damping = 1.;
    index_t exactPreconder = 0;
    index_t fdPreconder    = 0;
    index_t fdpfPreconder  = 0;
    index_t fdmgPreconder  = 0;
    index_t mlluPreconder  = 0;
    index_t mlmgPreconder  = 0;
    index_t stmgPreconder  = 1;
    std::string out;

    gsCmdLine cmd("parabolic_ls_example");
    cmd.addInt   ("g", "Geometry",              "0=Rectangle, 1=Quarter Annulus", geoIdx);
    cmd.addInt   ("r", "RefinementsX",          "Number of uniform h-refinement steps to perform before solving", refinementsX);
    cmd.addInt   ("s", "RefinementsT",          "Number of uniform tau-refinement steps to perform before solving", refinementsT);
    cmd.addInt   ("p", "Degree",                "Degree of the B-spline discretization space", degree);
    cmd.addReal  ("k", "Kappa",                 "Diffusion parameter", kappa);
    cmd.addReal  ("",  "ML.Sigma",              "Sigma for time-multilevel", sigma);
    cmd.addInt   ("",  "Solver.MaxIterations",  "Maximum iterations for linear solver", maxIterations);
    cmd.addReal  ("t", "Solver.Tolerance",      "Stopping criterion for linear solver", tolerance);
    cmd.addInt   ("",  "MG.NumPreSmooth",       "Number of pre smoothing steps (only for mg)", preSmooth);
    cmd.addInt   ("",  "MG.NumPostSmooth",      "Number of post smoothing steps (only for mg)", postSmooth);
    cmd.addInt   ("",  "MG.NumCycles",          "Number of multi-grid cycles for coarse-grid correction, i.e., 1=V, 2=W cycle", cycles);
    cmd.addReal  ("",  "MG.Damping",            "Damping factor for the smoother", damping);
    cmd.addInt   ("",  "useExactPreconder",     "Use that scheme", exactPreconder);
    cmd.addInt   ("",  "useFdPreconder",        "Use that scheme", fdPreconder);
    cmd.addInt   ("",  "useFdPfPreconder",      "Use that scheme", fdpfPreconder);
    cmd.addInt   ("",  "useFdMgPreconder",      "Use that scheme", fdmgPreconder);
    cmd.addInt   ("",  "useMlLuPreconder",      "Use that scheme", mlluPreconder);
    cmd.addInt   ("",  "useMlMgPreconder",      "Use that scheme", mlmgPreconder);
    cmd.addInt   ("",  "useStMgPreconder",      "Use that scheme", stmgPreconder);
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

    index_t iter1 = -1;
    real_t cond1  = -1;
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


    index_t iter2 = -1;
    real_t cond2  = -1;
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


    index_t iter3 = -1;
    real_t cond3  = -1;
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


    index_t iter4 = -1;
    real_t cond4  = -1;
    gsInfo << "Setup of FD+MG preconder... " << std::flush; // FD in time, mg (using Pearson factorzation) in space
    if (fdmgPreconder)
    {

        // Todo: make multigrid configurable

        gsLinearOperator<>::Ptr preconder = mkFdpfPreconder(time_stiff1, space_stiff, time_mass1, space_mass, mkMultiGridSolver{mb,bc,cmd.getGroup("MG")} );

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


    index_t iter5 = -1;
    real_t cond5  = -1;
    gsInfo << "Setup of ML+LU preconder... " << std::flush; // ML in time, lu (using Pearson factorzation) in space
    if (mlluPreconder)
    {

        gsLinearOperator<>::Ptr preconder = mkTimeMultiLevelPreconder(time_stiff1, space_stiff, time_mass1, space_mass, tb1, ic, cmd.getGroup("ML"), mkSparseLUSolver /** TODO Cholesky **/);

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

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  ML+MG preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    index_t iter6 = -1;
    real_t cond6  = -1;
    gsInfo << "Setup of ML+MG preconder... " << std::flush; // ML in time, mg (using Pearson factorzation) in space
    if (mlmgPreconder)
    {
        gsLinearOperator<>::Ptr preconder = mkTimeMultiLevelPreconder(time_stiff1, space_stiff, time_mass1, space_mass, tb1, ic, cmd.getGroup("ML"), mkMultiGridSolver{mb,bc,cmd.getGroup("MG")});

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

        iter6 = errorHistory.rows()-1;
        const bool success = errorHistory(iter6,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter6 << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter6 << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond6 = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond6 << "\n";
    }
    else
    {
        gsInfo << "skip.\n";
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  ST-MG preconder
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    index_t iter7 = -1;
    real_t cond7  = -1;
    gsInfo << "Setup of ST-MG preconder... " << std::flush;
    if (stmgPreconder)
    {
        gsLinearOperator<>::Ptr preconder = makeSpaceTimeMultiGridSolver(time_stiff1, space_stiff, time_mass1, space_mass, tb1, ic, mb, bc, cmd.getGroup("MG"));

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

        iter7 = errorHistory.rows()-1;
        const bool success = errorHistory(iter7,0) < tolerance;
        if (success)
            gsInfo << "Reached desired tolerance after " << iter7 << " iterations:\n";
        else
            gsInfo << "Did not reach desired tolerance after " << iter7 << " iterations:\n";

        if (errorHistory.rows() < 20)
            gsInfo << errorHistory.transpose() << "\n\n";
        else
            gsInfo << errorHistory.topRows(5).transpose() << " ... " << errorHistory.bottomRows(5).transpose()  << "\n\n";

        cond7 = solver.getConditionNumber();
        gsInfo << "Estimated condition number: " << cond7 << "\n";
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
            outfile << "parabolic_ls_example\t"
                "geoIdx\t"
                "refinementsX\t"
                "refinementsT\t"
                "degree\t"
                "kappa\t"
                "sigma\t"
                "preSmooth\t"
                "postSmooth\t"
                "cycles\t"
                "iter1\t"
                "cond1\t"
                "iter2\t"
                "cond2\t"
                "iter3\t"
                "cond3\t"
                "iter4\t"
                "cond4\t"
                "iter5\t"
                "cond5\t"
                "iter6\t"
                "cond6\t"
                "iter7\t"
                "cond7\n";

        outfile << "parabolic_ls_example\t"
            << geoIdx << "\t"
            << refinementsX << "\t"
            << refinementsT << "\t"
            << degree << "\t"
            << kappa << "\t"
            << sigma << "\t"
            << preSmooth << "\t"
            << postSmooth << "\t"
            << cycles << "\t"
            << iter1 << "\t"
            << cond1 << "\t"
            << iter2 << "\t"
            << cond2 << "\t"
            << iter3 << "\t"
            << cond3 << "\t"
            << iter4 << "\t"
            << cond4 << "\t"
            << iter5 << "\t"
            << cond5 << "\t"
            << iter6 << "\t"
            << cond6 << "\t"
            << iter7 << "\t"
            << cond7 << "\n";
    }


    return 0;

}
