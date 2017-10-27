
#pragma once

#include <gsPde/gsPde.h>

namespace gismo
{


template <class T> class gsFunction;

/** @brief
    A stationary Stokes PDE.

    This class describes a stationary Stokes PDE, with an arbitrary right-hand side
    function.

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

public:
    gsStokesPde(
        const gsMultiPatch<T>          &domain,
        const gsBoundaryConditions<T>  &bc,
         gsFunction<T>       *force,
        gsFunction<T>       *source = NULL,
        const T                    viscosity = 1
        )
        : gsPde<T>(domain,bc),    m_viscosity(viscosity)
    {
        m_force  = force  ? force->clone().release()  : NULL;
        m_source = source ? source->clone().release() : NULL;

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
