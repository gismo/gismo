/** @file flow-over-heated-plate.cpp

    @brief Heat equation participant for the PreCICE example "flow over heated plate"

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsXBraid/gsXBraid.h>

using namespace gismo;

// // Forward declaration of residual
// template<typename T>
// braid_Int my_Residual(braid_App app, braid_Vector ustop, braid_Vector r, braid_StepStatus status);

/**
   \brief Derived class implementing the XBraid wrapper for the heat equation
*/
template<typename T>
class gsXBraidHeatEquation : public gsXBraid< gsMatrix<T> >
{
public:

    gsXBraidHeatEquation(   const gsMpiComm & comm,
                            const T&          tstart,
                            const T&          tstop,
                            const T&          theta,
                            index_t           numSteps,
                            index_t           numRefine,
                            index_t           numElevate,
                            std::string&      fn,
                            bool              residual = false)
    :
    gsXBraid< gsMatrix<T> >(comm, tstart, tstop, (int)numSteps),
    m_numRefine(numRefine),
    m_numElevate(numElevate),
    m_numSteps(numSteps),
    m_tstart(tstart),
    m_tstop(tstop),
    m_theta(theta),
    m_tstep( (tstop-tstart)/numSteps ),
    m_A(1,1)
    {
        /////////////////////////////////////////////////////////////////////////////////////////////
        //                           Code for heat equation starts here                            //
        /////////////////////////////////////////////////////////////////////////////////////////////
        gsFileData<> fd(fn);

        fd.getId(0, m_mp); // id=0: Multipatch domain
        fd.getId(1, m_f); // id=1: source function

        fd.getId(2, m_bc); // id=2: boundary conditions
        m_bc.setGeoMap(m_mp);
        fd.getId(5, m_Aopt); // id=5: assembler options
        fd.getId(6, m_Topt); // id=6: MGRIT options

        m_mb = gsMultiBasis<T>(m_mp);
        m_mb.setDegree( m_mb.maxCwiseDegree() + numElevate);
        for (int r =0; r < numRefine; ++r)
            m_mb.uniformRefine();

        m_A.setIntegrationElements(m_mb);
        auto G = m_A.getMap(m_mp);
        auto u = m_A.getSpace(m_mb);
        auto ff= m_A.getCoeff(m_f,G);
        auto g = m_A.getBdrFunction(G);

        // Assemble mass matrix
        u.setup(m_bc, dirichlet::homogeneous, 0);
        m_A.initSystem();
        m_A.assemble( u * u.tr() * meas(G));
        m_M = m_A.matrix();

        // Assemble F (assume constant forcing) and K (constant BCs)
        u.setup(m_bc, dirichlet::l2Projection, 0);
        m_A.initSystem();
        m_A.assembleBdr(m_bc.get("Neumann"), u * g.val() * nv(G).norm() );
        m_A.assemble( igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G) );
        m_F = m_A.rhs();
        m_K = m_A.matrix();

        //////////////////////////////////////////////////////////////////////////
        // MGRIT options
        //////////////////////////////////////////////////////////////////////////
        fd.getId(6, m_Topt); // id=6: multigrid-in-time options
        if (this->id() == 0) gsInfo << "Multigrid-in-time options:\n" << m_Topt << "\n";

        this->SetCFactor(m_Topt.getInt("CFactor"));
        this->SetMaxIter(m_Topt.getInt("maxIter"));
        this->SetMaxLevels(m_Topt.getInt("maxLevel"));
        this->SetMaxRefinements(m_Topt.getInt("numMaxRef"));
        this->SetMinCoarse(m_Topt.getInt("minCLevel"));
        this->SetNFMG(m_Topt.getInt("numFMG"));
        this->SetNFMGVcyc(m_Topt.getInt("numFMGVcyc"));
        this->SetNRelax(m_Topt.getInt("numRelax"));
        this->SetAccessLevel(m_Topt.getInt("access"));
        this->SetPrintLevel(m_Topt.getInt("print"));
        this->SetStorage(m_Topt.getInt("numStorage"));
        this->SetTemporalNorm(m_Topt.getInt("norm"));

        if (m_Topt.getSwitch("fmg"))           this->SetFMG();
        if (m_Topt.getSwitch("incrMaxLevels")) this->SetIncrMaxLevels();
        if (m_Topt.getSwitch("periodic"))      this->SetPeriodic(1); else this->SetPeriodic(0);
        if (m_Topt.getSwitch("refine"))        this->SetRefine(1);   else this->SetRefine(0);
        if (m_Topt.getSwitch("sequential"))    this->SetSeqSoln(1);  else this->SetSeqSoln(0);
        if (m_Topt.getSwitch("skip"))          this->SetSkip(1);     else this->SetSkip(0);
        if (m_Topt.getSwitch("spatial"))       this->SetSpatialCoarsenAndRefine();
        if (m_Topt.getSwitch("tol"))           this->SetAbsTol(m_Topt.getReal("absTol"));
        else                                 this->SetRelTol(m_Topt.getReal("relTol"));

        // Custom residual
        if (residual) this->SetResidual();
    }

    braid_Int Init(braid_Real    t,
                   braid_Vector *u_ptr)
    {
        gsMatrix<T>* u = new gsMatrix<T>(m_A.numDofs(), 1);

        if (t != m_tstart)
        {
            // Intermediate solution
            u->setZero(m_A.numDofs(),1);
        } else {
            // Initial solution
            u->setZero(m_A.numDofs(),1);
        }

        *u_ptr = (braid_Vector) u;
        return braid_Int(0);
    }

    braid_Int Step(braid_Vector     u,
                   braid_Vector     ustop,
                   braid_Vector     fstop,
                   BraidStepStatus& status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
    {
        gsMatrix<T>* u_ptr = (gsMatrix<T>*) u;
        gsMatrix<T>* ustop_ptr = (gsMatrix<T>*) ustop;

        // Get time step information
        std::pair<braid_Real, braid_Real> time = static_cast<gsXBraidStepStatus&>(status).timeInterval();
        T dt(time.second - time.first);

        // Solve time step
        gsMatrix<T> RHS = m_theta*dt*m_F + (1.0-m_theta)*dt*m_F + (m_M-dt*(1.0-m_theta)*m_K)*(*u_ptr);
        // XBraid forcing
        if (fstop != NULL)
        {
            gsMatrix<T>* fstop_ptr = (gsMatrix<T>*) fstop;
            RHS += *fstop_ptr;
        }

        // Solve the linear system
        (*u_ptr) = m_solver.compute(m_M + dt*m_theta*m_K).solve(RHS);

        // Carry out adaptive refinement in time
        if (static_cast<gsXBraidStepStatus&>(status).level() == 0)
        {
            braid_Real error = static_cast<gsXBraidStepStatus&>(status).error();
            if (error != braid_Real(-1.0))
            {
                braid_Int rfactor = (braid_Int) std::ceil( std::sqrt( error / 1e-3) );
                status.SetRFactor(rfactor);
            }
            else
                status.SetRFactor(1);
        }

        return braid_Int(0);
    }

    braid_Int Residual(braid_Vector     u,
                       braid_Vector     r,
                       BraidStepStatus& status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
    {
        gsMatrix<T>* u_ptr = (gsMatrix<T>*) u;
        gsMatrix<T>* r_ptr = (gsMatrix<T>*) r;

        // Get time step information
        std::pair<braid_Real, braid_Real> time = static_cast<gsXBraidStepStatus&>(status).timeInterval();
        T dt(time.second - time.first);

        // Compute residual
        *r_ptr = (m_M+dt*m_theta*m_K)*(*u_ptr) - (m_theta*dt*m_F + (1.0-m_theta)*dt*m_F + (m_M-dt*(1.0-m_theta)*m_K)*(*r_ptr));
        return braid_Int(0);
    }

    braid_Int BufSize(braid_Int         *size_ptr,
                      BraidBufferStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
    {
        *size_ptr = sizeof(T)*(m_A.numDofs()+2);
        return braid_Int(0);
    }

    braid_Int Access(braid_Vector       u,
                     BraidAccessStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
    {
        if (static_cast<gsXBraidAccessStatus&>(status).done() &&
            static_cast<gsXBraidAccessStatus&>(status).timeIndex() ==
            static_cast<gsXBraidAccessStatus&>(status).times())
        {
            gsMatrix<T>* u_ptr = (gsMatrix<T>*) u;
            gsInfo << "norm of the solution = " << u_ptr->norm() << std::endl;
        }
        return braid_Int(0);
    }

    braid_Int Coarsen(braid_Vector           fu,
                      braid_Vector          *cu_ptr,
                      BraidCoarsenRefStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
    {
        gsMatrix<T> *fu_ptr = (gsMatrix<T>*) fu;
        gsMatrix<T>* cu     = new gsMatrix<T>();
        *cu = *fu_ptr;
        *cu_ptr = (braid_Vector) cu;
        return braid_Int(0);
    }

    braid_Int Refine(braid_Vector           cu,
                     braid_Vector          *fu_ptr,
                     BraidCoarsenRefStatus &status)
#if __cplusplus >= 201103L || _MSC_VER >= 1600
    override
#endif
    {
        gsMatrix<T> *cu_ptr = (gsMatrix<T>*) cu;
        gsMatrix<T>* fu     = new gsMatrix<T>();
        *fu = *cu_ptr;
        *fu_ptr = (braid_Vector) fu;
        return braid_Int(0);
    }

public:
    // getters
    const gsSparseMatrix<T> & M()       const { return m_M; }
    const gsSparseMatrix<T> & K()       const { return m_K; }
    const gsMatrix<T>       & F()       const { return m_F; }
    const T                 & theta()   const { return m_theta; }

protected:
    gsExprAssembler<T> m_A;
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    gsFunctionExpr<T> m_f;
    gsBoundaryConditions<T> m_bc;
    gsOptionList m_Aopt;
    gsOptionList m_Topt;

    gsSparseMatrix<T> m_M;
    gsSparseMatrix<T> m_K;
    gsMatrix<T>       m_F;

    index_t m_numRefine;
    index_t m_numElevate;
    index_t m_numSteps;
    T m_tstart;
    T m_tstop;
    T m_theta;
    T m_tstep;

    typename gsSparseSolver<T>::LU m_solver;
};

// template<typename T>
// braid_Int my_Residual(  braid_App        app,
//                         braid_Vector     ustop,
//                         braid_Vector     r,
//                         BraidStepStatus  status)
// {
//     gsMatrix<T>* ustop_ptr = (gsMatrix<T>*) ustop;
//     gsMatrix<T>* r_ptr = (gsMatrix<T>*) r;

//     // Get time step information
//     std::pair<braid_Real, braid_Real> time = static_cast<gsXBraidStepStatus&>(status).timeInterval();
//     T dt(time.second - time.first);

//     // Get the system from the app
//     gsXBraidHeatEquation<T>* app_ptr = (gsXBraidHeatEquation<T>*) app;
//     const & gsSparseMatrix<T> M = app_ptr->M();
//     const & gsSparseMatrix<T> K = app_ptr->K();
//     const & gsMatrix<T>       F = app_ptr->F();
//     const & T             theta = app_ptr->theta();

//     // Compute residual
//     *r_ptr = (M+dt*theta*K)*(*u_ptr) - (theta*dt*F + (1.0-theta)*dt*F + (M-dt*(1.0-theta)*K)*(*r_ptr));
//     return braid_Int(0);
// }

int main(int argc, char *argv[])
{
#ifdef gsXBraid_ENABLED
    //! [Parse command line]
    bool plot = false;
    bool res  = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    index_t steps = 10;
    real_t theta = 1.0;
    std::string fn("pde/heat2d_square_ibvp1.xml");

    gsCmdLine cmd   ("Tutorial on solving a Poisson problem.");
    cmd.addInt      ("e"    , "degreeElevation" ,"Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt      ("r"    , "uniformRefine"   ,"Number of Uniform h-refinement loops",  numRefine );
    cmd.addInt      ("N"    , "steps"           ,"Number of time steps", steps );
    cmd.addReal     ("t"    , "theta"           ,"Theta value for the time-stepping scheme", theta );
    cmd.addString   ("f"    , "file"            ,"Input XML file", fn );
    cmd.addSwitch   ("plot" ,                    "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch   ("res"  ,                    "Use custom residual", res);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    real_t tstart = 0.0;
    real_t tstop  = 1.0;

    // Initialize the MPI environment and obtain the world communicator
    gsMpiComm comm = gsMpi::init(argc, argv).worldComm();

    // Print MPI/OpenMP configuration
    if (comm.rank() == 0)
    {
        gsInfo << "Number of MPI processes    : " << comm.size() << std::endl;
#ifdef _OPENMP
        gsInfo << "Number of OpenMP processes : " << omp_get_num_procs() << std::endl;
#endif
    }


    gsXBraidHeatEquation<real_t> app(comm, tstart, tstop, theta, steps, numRefine, numElevate, fn,res);

    // Perform parallel-in-time multigrid
    app.solve();

#else

  gsInfo << "\n";

#endif

  return 0;
}
