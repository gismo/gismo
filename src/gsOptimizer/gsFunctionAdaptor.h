/** @file gsMinimizer.h

    @brief Class definiting an optimization problem for a gsFunction

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
 * @brief 
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
        for (index_t i=0; i!=m_numDesignVars; ++i)
            if ( u(i)<m_desLowerBounds.at(i) || u(i)>m_desUpperBounds.at(i) )
                return math::limits::max();
        return m_obj.eval(u).value();
    }

    mutable util::gsThreaded<gsMatrix<T> > jac;
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const  
    {
        for (index_t i=0; i!=m_numDesignVars; ++i)
            if ( u(i)<m_desLowerBounds.at(i) || u(i)>m_desUpperBounds.at(i) )
            {
                result.setZero();
                return;
            }

        //gsOptProblem<T>::gradObj_into(u,result);
        //gsDebugVar( result.transpose() );
        m_obj.deriv_into(u, jac);
        //gsDebugVar( jac.transpose() );
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
