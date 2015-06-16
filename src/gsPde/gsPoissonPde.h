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

namespace gismo
{

template <class T> class gsFunction;

/** @brief
    A Poisson PDE.

    This class describes a Poisson PDE, with an arbitrary right-hand side
    function and optionally a known solution.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsPoissonPde : public gsPde<T>
{
protected:
    gsPoissonPde( ) { }
    using gsPde<T>::m_domain;
public:
    /// Constructor
    gsPoissonPde(
        const gsMultiPatch<T>         &domain,
        const gsBoundaryConditions<T> &bc,
        const gsFunction<T>     &rhs
         )
    : gsPde<T>(domain,bc)
    {
        m_rhs=new gsPiecewiseFunction<T>(rhs,m_domain.size());
        m_unknownDim.push_back(1);
    }

    /// Constructor
    gsPoissonPde(
        const gsMultiPatch<T>         &domain,
        const gsBoundaryConditions<T> &bc,
        const gsPiecewiseFunction<T>  &rhs
         )
    : gsPde<T>(domain,bc), m_rhs(rhs.clone())
    {
        m_unknownDim.push_back(1);
    }

    /// Constructor
    gsPoissonPde(
        const gsMultiPatch<T>         &domain,
        const gsBoundaryConditions<T> &bc,
        const gsPiecewiseFunction<T>  &rhs,
        const gsPiecewiseFunction<T>  &sol
         )
    : gsPde<T>(domain,bc), m_rhs(rhs.clone())
    {
        m_solution.push_back(sol.clone());
        m_unknownDim.push_back(1);
    }

    // FOR COMPATIBILITY WITH OLD STRUCTURE: DO NOT USE
    int m_compat_dim;
    gsPoissonPde(
        const gsFunction<T>  &rhs,
        int                   domdim,
        const gsFunction<T>  &sol
        )
        : m_compat_dim(domdim)
    {
        m_rhs=new gsPiecewiseFunction<T>(rhs);
        m_solution.push_back(new gsPiecewiseFunction<T>(sol));
        m_unknownDim.push_back(1);

    }
    gsPoissonPde(
        const gsFunction<T>  &rhs,
        int                   domdim
        )
       : m_compat_dim(domdim)

    {
        m_rhs=new gsPiecewiseFunction<T>(rhs);
        m_unknownDim.push_back(1);
    }
    gsPoissonPde(
        void * unused
        )
    {
        m_rhs=new gsPiecewiseFunction<T>(gsConstantFunction<T>(0));
        m_unknownDim.push_back(1);
    }

    ~gsPoissonPde( ) 
    { 
        delete m_rhs;
    }

    /**
     * @brief gives the number of rhs functions of the PDEs
     */
    virtual int numRhs() const
    {
        return m_rhs->operator[](0).targetDim();
    }

    gsPiecewiseFunction<T> * rhsField() const { return m_rhs; }
    const gsFunction<T> *          rhs(index_t k=0)      const { return &(m_rhs->operator[](k)); }

    virtual int numUnknowns() const     {return 1;}

    virtual bool isSymmetric () const { gsWarn<<"Function is gsPde::isSymmetric should not be used!!"; return true;}

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
	    os<<"Poisson's equation  -\u0394u = f ,  with:\n";
	    os<<"Source function f= "<< * m_rhs <<".\n";
	    if ( this->solutionGiven() )
            os<<"Exact solution g = "<< * this->m_solution[0] <<".\n";
	    return os; 
	}
protected:
    using gsPde<T>::m_unknownDim;
    using gsPde<T>::m_solution;
    gsPiecewiseFunction<T> *m_rhs;
}; // class gsPoissonPde

} // namespace gismo
