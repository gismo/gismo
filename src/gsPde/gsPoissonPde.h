/** @file gsPoissonPde.h

    @brief Describes a Poisson PDE.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsPde/gsPde.h>
#include <gsCore/gsPiecewiseFunction.h>

namespace gismo
{

/** @brief
    A Poisson PDE.

    This class describes a Poisson PDE, with an arbitrary right-hand side
    function.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsPoissonPde : public gsPde<T>
{

public:
    
    gsPoissonPde( ) { }


    /// Constructor
    gsPoissonPde(const gsMultiPatch<T>         &domain,
                 const gsBoundaryConditions<T> &bc,
                 const gsPiecewiseFunction<T>  &rhs,
                 const gsFunction<T>           * = NULL)
    : gsPde<T>(domain,bc), m_rhs(rhs)
    {
        m_unknownDim.setOnes(1);
    }

    int m_compat_dim;
    GISMO_DEPRECATED
    gsPoissonPde(const gsFunction<T>  &rhs,
                 int                   domdim,
                 const gsFunction<T>  &)
    : m_compat_dim(domdim), m_rhs(rhs)
    {
        m_unknownDim.setOnes(1);

    }
    GISMO_DEPRECATED
    gsPoissonPde(const gsFunction<T>  &rhs,
                 int                   domdim)
       : m_compat_dim(domdim), m_rhs(rhs)

    {
        m_unknownDim.setOnes(1);
    }


    GISMO_DEPRECATED
    gsPoissonPde(void * unused)
    {
        m_rhs=new gsConstantFunction<T>(0);
        m_unknownDim.setOnes(1);
    }



    /**
     * @brief gives the number of rhs functions of the PDEs
     */
    virtual int numRhs() const
    {
        return m_rhs.piece(0).targetDim();
    }

    const gsFunction<T> *    rhs()      const { return &m_rhs.piece(0); }

    virtual int numUnknowns() const     {return 1;}

    virtual bool isSymmetric () const { gsWarn<<"Function is gsPde::isSymmetric should not be used!!"; return true;}

    /// Prints the object as a string.
    virtual std::ostream &print(std::ostream &os) const
    {
	    os<<"Poisson's equation  -\u0394u = f ,  with:\n";
	    os<<"Source function f= "<< m_rhs <<".\n";
	    return os; 
	}

    virtual gsPde<T>* restrictToPatch(unsigned np) const
    {
        gsBoundaryConditions<T> bc;
        m_boundary_conditions.getConditionsForPatch(np,bc);
        return new gsPoissonPde<T>(m_domain.patch(np),bc,m_rhs);
    }

protected:
    using gsPde<T>::m_unknownDim;
    using gsPde<T>::m_domain;
    using gsPde<T>::m_boundary_conditions;

    gsPiecewiseFunction<T> m_rhs;
}; // class gsPoissonPde

} // namespace gismo
