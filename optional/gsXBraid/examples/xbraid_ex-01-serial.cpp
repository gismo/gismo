/** @file xbraid_example.cpp

    @brief XBraid integration

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, M. Moeller
*/

#include <gismo.h>
#include <gsXBraid/gsTimeIntegratorBase.h>


using namespace gismo;

template <typename T>
class MyTimeOp : public TimeOpBase<T,T,T,T>
{

  using Base = TimeOpBase<T,T,T,T>;

public:

  typedef typename Base::SolType SolType;
  typedef typename Base::ResOp ResOp;
  typedef typename Base::JacOp JacOp;

  // typename typedef SolType T;
  // typename typedef ResOp T;
  // typename typedef JacOp T;

  MyTimeOp()
  {}

  ResOp residual(SolType & u, SolType & du, T t) const
  {
    return -math::sin(u);
  }
  JacOp jacobian(SolType & u, SolType & du, T t) const
  {
    GISMO_NO_IMPLEMENTATION;
  }
  index_t size() const
  {
    return 1;
  }

};

template <typename T, typename OpType>
class ExplicitEuler : public gsTimeIntegratorBase<T,OpType>
{
  typedef typename OpType::SolType SolType;
  typedef typename OpType::ResOp ResOp;
  typedef typename OpType::JacOp JacOp;
public:
  ExplicitEuler(OpType op)
  :
  m_op(op)
  {}

    bool step(SolType & u, SolType & du, T tstart, T tstop) const
    {
      T dt = tstop-tstart;
      u += dt * m_op.residual(u,du,tstart);
      return true;
    }

protected:
  OpType m_op;

};

template <typename T, typename OpType>
class ExplicitEuler : public gsTimeIntegratorBase<T,OpType>
{
  typedef typename OpType::SolType SolType;
  typedef typename OpType::ResOp ResOp;
  typedef typename OpType::JacOp JacOp;
public:
  ExplicitEuler(OpType op)
  :
  m_op(op)
  {}

    bool step(SolType & u, SolType & du, T tstart, T tstop) const
    {
      T dt = tstop-tstart;
      u += dt * m_op.residual(u,du,tstart);
      return true;
    }

protected:
  OpType m_op;

};

template <typename T, typename OpType>
class XBraidIntegrator : public gsXBraid<OpType::SolType>
{
  using Base = gsXBraid<OpType::SolType>;
  typedef typename gsTimeIntegratorBase<T,OpType> StepperType;

public:

  XBraidIntegrator( const StepperType & stepper,
                    const gsMpiComm& comm,
                    const T&         tstart,
                    const T&         tstop,
                    index_t          numSteps)
  :
  Base(comm,tstart,tstop,(int) numSteps)
  {}

  /// Initializes a vector
  braid_Int Init(braid_Real    t, braid_Vector *u_ptr) override
  {
      gsMatrix<T>* u = new gsMatrix<T>(3*m_numDofs, 1);

      // Does this mean zero displacements?
      u->setZero();

      if (m_solver->solutionU().rows()==m_numDofs)
          u->col(0).segment(0          ,m_numDofs) = m_solver->solutionU();
      if (m_solver->solutionV().rows()==m_numDofs)
          u->col(0).segment(m_numDofs  ,m_numDofs) = m_solver->solutionV();
      if (m_solver->solutionA().rows()==m_numDofs)
          u->col(0).segment(2*m_numDofs,m_numDofs) = m_solver->solutionA();

      *u_ptr = (braid_Vector) u;
      return braid_Int(0);
  }

  /// Performs a single step of the parallel-in-time multigrid
  braid_Int Step(braid_Vector    u, braid_Vector    ustop, braid_Vector    fstop, BraidStepStatus &status) override
  {
      gsVector<T>* u_ptr = (gsVector<T>*) u;
      // gsMatrix<T>* ustop_ptr = (gsMatrix<T>*) ustop; // the guess is not used

      // XBraid forcing
      if (fstop != NULL)
      {
          gsVector<T>* fstop_ptr = (gsVector<T>*) fstop;
          *u_ptr += *fstop_ptr;
      }

      // Get time step information
      std::pair<braid_Real, braid_Real> time = static_cast<gsXBraidStepStatus&>(status).timeInterval();
      if (m_options.getSwitch("extraVerbose")) gsInfo<<"Solving interval ["<<time.first<<" , "<<time.second<<"] (level "<<static_cast<gsXBraidStepStatus&>(status).level()<<")\n";
      T t  = time.first;
      T dt = time.second - time.first;

      // Solve time step
      gsVector<T> U = (*u_ptr).segment(0          ,m_numDofs);
      gsVector<T> V = (*u_ptr).segment(m_numDofs  ,m_numDofs);
      gsVector<T> A = (*u_ptr).segment(2*m_numDofs,m_numDofs);

      gsStatus stepStatus = m_solver->step(t,dt,U,V,A);

      u_ptr->segment(0          ,m_numDofs) = U;
      u_ptr->segment(m_numDofs  ,m_numDofs) = V;
      u_ptr->segment(2*m_numDofs,m_numDofs) = A;

      // Carry out adaptive refinement in time
      if (static_cast<gsXBraidStepStatus&>(status).level() == 0)
      {
          if (stepStatus==gsStatus::Success)
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
          // Refine if solution interval failed
          else
          {
              if (m_options.getSwitch("extraVerbose")) gsInfo<<"Step "<<(static_cast<gsXBraidStepStatus&>(status)).timeIndex()<<" did not converge";
              status.SetRFactor((braid_Int)2);
          }
      }
      return braid_Int(0);
  }

  /// Computes the spatial norm of the given vector
  braid_Int SpatialNorm(  braid_Vector  u,
                          braid_Real   *norm_ptr) override
  {
      gsVector<T>* u_ptr = (gsVector<T>*) u;
      *norm_ptr = u_ptr->segment(0,m_numDofs).norm(); // Displacement-based norm
      // *norm_ptr = u_ptr->norm();
      return braid_Int(0);
  }

  /// Sets the size of the MPI communication buffer
  braid_Int BufSize(braid_Int *size_ptr, BraidBufferStatus &status) override
  {
      *size_ptr = sizeof(T)*(m_numDofs*3+2); // +2 comes from rows, cols of the solution vector u.
      return braid_Int(0);
  }

  void setCallback(callback_type callback) const {m_callback = callback;}

  /// Handles access for input/output
  braid_Int Access(braid_Vector u, BraidAccessStatus &status) override
  {
      gsVector<T>* u_ptr = (gsVector<T>*) u;
      m_callback((index_t)    static_cast<gsXBraidAccessStatus&>(status).timeIndex(),
                 (T)          static_cast<gsXBraidAccessStatus&>(status).time(),
                 (gsVector<T>)(*u_ptr).segment(0          ,m_numDofs),
                 (gsVector<T>)(*u_ptr).segment(m_numDofs  ,m_numDofs),
                 (gsVector<T>)(*u_ptr).segment(2*m_numDofs,m_numDofs)
                 );
      return braid_Int(0);
  }

  /// Performs spatial coarsening
  /*
      NOTE: This routine is not implemented. How to do it:
      1. Make the Coarsen and Refine routines virtual
      2. Within the example, overload this class, and define Coarsen and Refine for the problem
      3. Update the m_numDoFs within!!
  */
  braid_Int Coarsen(braid_Vector fu, braid_Vector *cu_ptr, BraidCoarsenRefStatus &status) override
  {
      gsMatrix<T> *fu_ptr = (gsMatrix<T>*) fu;
      gsMatrix<T>* cu     = new gsMatrix<T>();
      *cu = *fu_ptr;
      *cu_ptr = (braid_Vector) cu;
      return braid_Int(0);
  }

  // Performs spatial refinement
  braid_Int Refine(braid_Vector cu, braid_Vector *fu_ptr, BraidCoarsenRefStatus &status) override
  {
      gsMatrix<T> *cu_ptr = (gsMatrix<T>*) cu;
      gsMatrix<T>* fu     = new gsMatrix<T>();
      *fu = *cu_ptr;
      *fu_ptr = (braid_Vector) fu;
      return braid_Int(0);
  }

protected:

  StepperType * m_stepper;

}

int main(int argc, char**argv)
{

  index_t maxSteps = 10;
  real_t t = 0;
  real_t tend = 1;

  gsCmdLine cmd("\nExample 1: Solve a scalar ODE \n\n");
  cmd.addInt("N","maxSteps","Maximum number of steps",maxSteps);
  cmd.addReal("b","begTime","Start time",t);
  cmd.addReal("e","endTime","End time",tend);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  MyTimeOp<real_t> op;
  ExplicitEuler<real_t,MyTimeOp<real_t>> EE(op);
  ExplicitEuler<real_t,MyTimeOp<real_t>> IE(op);
  IMEX<real_t,MyTimeOp<real_t>> IMEX(op);

  real_t dt = (tend-t)/maxSteps;
  real_t u = 1;
  for (index_t step = 0; step!=maxSteps; step++)
  {
    EE.step(u,u,t,t+dt);
    gsDebugVar(u);
  }

  gsXBraidIntegrator<real_t,MyTimeOp<real_t>> xbraid(EE);
  xbraid.solve();

  return 0;

}
