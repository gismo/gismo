/** @file gsFunctionAdaptor.h

    @brief Class defining an adaptor for using a gsFunction in an
    optimization problem

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsOptimizer/gsOptProblem.h>
#include <gsUtils/gsThreaded.h>

namespace gismo
{


/**
 * @brief Adaptor to see a given gsFunction as (the objective of) an unconstrained optimization problem
 *
 */

template <typename T>
class gsFunctionAdaptor : public gsOptProblem<T>
{
    const gsFunction<T> & m_obj;
public:

    gsFunctionAdaptor(const gsFunction<T> & obj)//, const gsVector<T> &  ,bool withSupport = true)
    : m_obj(obj)
    {
        m_numDesignVars  = obj.domainDim();
        m_numConstraints = 0;
        m_numConJacNonZero = 0;

        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);

        gsMatrix<T> sup = obj.support();
        m_desLowerBounds = sup.col(0); //! bound_relax_factor
        m_desUpperBounds = sup.col(1);

        m_curDesign = 0.5 * (m_desLowerBounds + m_desUpperBounds);
        //gsDebugVar( m_curDesign.transpose() );
    }


public:

    T evalObj( const gsAsConstVector<T> & u ) const
    {
        gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() ); // keep point within bounds
        u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);
        return m_obj.eval(u).value();
    }

    mutable util::gsThreaded<gsMatrix<T> > jac;
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {

        //gsOptProblem<T>::gradObj_into(u,result);
        //gsDebugVar( result.transpose() );
        m_obj.deriv_into(u, jac);
        //gsDebugVar( jac.transpose() );
        result = jac.mine();
    }

    void hessObj_into( const gsAsConstVector<T> & u, gsAsMatrix<T> & result) const
    {
        m_obj.hessian_into(u, jac);
        result = jac.mine();
    }

    void evalCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {

    }

    void jacobCon_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {

    }

private:

    using gsOptProblem<T>::m_numDesignVars;
    using gsOptProblem<T>::m_numConstraints;
    using gsOptProblem<T>::m_numConJacNonZero;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_conLowerBounds;
    using gsOptProblem<T>::m_conUpperBounds;

    using gsOptProblem<T>::m_conJacRows;
    using gsOptProblem<T>::m_conJacCols;

    using gsOptProblem<T>::m_curDesign;
};

} // end namespace gismo
