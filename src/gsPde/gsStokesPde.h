
#pragma once

#include <gsPde/gsPde.h>

namespace gismo
{


template <class T> class gsFunction;

/** @brief
    A stationary Stokes PDE.

    This class describes a stationary Stokes PDE, with an arbitrary right-hand side
    function and optionally a known solution.

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsStokesPde : public gsPde<T>
{
protected:
    gsStokesPde( ) { }
    using gsPde<T>::m_domain;
    using gsPde<T>::m_unknownDim;
    using gsPde<T>::m_solution;

public:
    gsStokesPde(
        const gsMultiPatch<T>          &domain,
        const gsBoundaryConditions<T>  &bc,
         gsFunction<T>       *force,
         gsFunction<T>       *velSol = NULL,
         gsFunction<T>       *preSol = NULL,
         gsFunction<T>       *source = NULL,
        const T                    viscosity = 1
        )
        :
            gsPde<T>(domain,bc),    m_viscosity(viscosity)

    {
        m_force  = force  ? force->clone()  : NULL;
        m_source = source ? source->clone() : NULL;
        m_solution.resize(2);
        m_solution[0] = velSol ? velSol->clone() : NULL;
        m_solution[1] = preSol ? preSol->clone() : NULL;

        m_unknownDim.resize(2);
        m_unknownDim[0] = m_domain.dim();
        m_unknownDim[1] = 1;
    }

    ~gsStokesPde( ) 
    { 
        delete m_force;
        delete m_source;
    }

    const gsFunction<T>* rhs() const
    { return m_force; }
    const gsFunction<T>* force() const
    { return m_force; }
    const gsFunction<T>* source() const
    { return m_source; }

    const gsFunction<T>* velocitySolution() const
    { return m_solution[0]; }
    const gsFunction<T>* pressureSolution() const
    { return m_solution[1]; }
    T viscocity() const                 { return m_viscosity; }


    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os<<"Stokes's equation:\n"
          <<"-\u0394u-\u2207p = f,\n"
          <<" \u2207\u00B7u=0"
          <<"with:\n";
        if ( m_force )
        os<<"Force  function f= "<< *m_force <<".\n";
        if ( m_source )
        os<<"Source function g= "<< *m_source <<".\n";
        if ( m_solution[0] )
            os<<"Exact solution u = "<< * m_solution[0] <<".\n";
        if ( m_solution[1] )
            os<<"Exact solution p = "<< * m_solution[1] <<".\n";
	    return os; 
	}
    /// Consistency check
    bool check()
    {
        return true;
    }
protected:
    const gsFunction<T> * m_force;
    const gsFunction<T> * m_source;

    T m_viscosity;
}; // class gsStokesPde

} // namespace gismo
