/** @file gsTimeIntegratorBase.h

    @brief

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): H.M. Verhelst
*/

#pragma once

#include <gismo.h>

namespace gismo
{

  template <typename T, typename SolType_t, typename ResOp_t, typename JacOp_t>
  class TimeOpBase
  {
  public:

    using SolType = SolType_t;
    using ResOp = ResOp_t;
    using JacOp = JacOp_t;

    virtual ResOp residual(SolType & u, SolType & du, T t) const = 0;
    virtual JacOp jacobian(SolType & u, SolType & du, T t) const = 0;
    virtual index_t size() const = 0;
  };

  template<typename T, typename OpType>
  class gsTimeIntegratorBase
  {
    typedef typename OpType::SolType SolType;
    typedef typename OpType::ResOp ResOp;
    typedef typename OpType::JacOp JacOp;
  public:
    virtual bool step(SolType & u, SolType & du, T tstart, T tstop) const = 0;
  };

  template<typename T, typename OpType >
  class gsTimeIntegratorExplicitEuler : public gsTimeIntegratorBase<T,OpType>
  {
    typedef typename OpType::SolType SolType;
    typedef typename OpType::ResOp ResOp;
    typedef typename OpType::JacOp JacOp;

  public:
    /// Constructor
    gsTimeIntegratorExplicitEuler(const OpType & op)
    :
    m_op(op)
    {}

    bool step(gsMatrix<T> & u, gsMatrix<T> & du, T tstart, T tstop) const override
    {
      T dt = tstop - tstart;
      // m_op.residual(u, du, tstart);
      // du *= dt;
      // u += du;
      return true;
    }

  protected:
    OpType m_op;

  };

  template<typename T, typename OpType >
  class gsTimeIntegratorImplicitEuler : public gsTimeIntegratorBase<T,OpType>
  {
    typedef typename OpType::SolType SolType;
    typedef typename OpType::ResOp ResOp;
    typedef typename OpType::JacOp JacOp;

  public:
    /// Constructor
    gsTimeIntegratorImplicitEuler(const OpType & op)
    :
    m_op(op)
    {}

    bool step(gsMatrix<T> & u, gsMatrix<T> & du, T tstart, T tstop) const override
    {
      T dt = tstop - tstart;
      ResOp res = m_op.residual(u, du, tstart);
      JacOp jac = m_op.jacobian(u, du, tstart);
      // Solve the system
      typename gsSparseSolver<T>::LU solver(jac);
      u = solver.solve(res);
      // du *= dt;
      // u += du;
      return true;
    }

  protected:
    OpType m_op;

  };

}// namespace gismo




// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsTimeIntegratorBase.hpp)
// #endif
