
#pragma once

#include <gsPde/gsPde.h>

namespace gismo
{


// pressure (p):  p=2, 
// velocity (u):  p=3, ( p-refined ) , two components, one for each dimension

// \nu constant = 0.001

template <class T> class gsFunction;

/** @brief
    A stationary Stokes PDE.

    This class describes a stationary Stokes PDE, with an arbitrary right-hand side
    function and optionally a known solution.

    \ingroup Pde
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
    gsStokesPde(const gsMultiPatch<T> &domain, const gsBoundaryConditions<T> &bc, gsFunction<T> * rhs, gsFunction<T> * sol = 0 )
        : gsPde<T>(domain,bc), m_rhs(rhs), m_viscosity(0.001)
    { 
        this->m_unknownDim.resize(2);
        this->m_unknownDim[0] = m_domain.dim();
        this->m_unknownDim[1] = 1;
        this->m_solution.push_back(new gsPiecewiseFunction<T>(sol));
    }

    gsStokesPde(const gsMultiPatch<T> &domain, const gsBoundaryConditions<T> &bc, const gsPiecewiseFunction<T> &rhs, const gsPiecewiseFunction<T> &sol)
        : gsPde<T>(domain,bc), m_rhs(rhs.clone()), m_viscosity(0.001)
    { 
        this->m_unknownDim.resize(2);
        this->m_unknownDim[0] = m_domain.dim();
        this->m_unknownDim[1] = 1;
        this->m_solution.push_back(sol.clone());
    }

    gsStokesPde(const gsMultiPatch<T> &domain, const gsBoundaryConditions<T> &bc, const gsFunction<T> & rhs)
        : gsPde<T>(domain,bc), m_rhs(rhs.clone()), m_viscosity(0.001)
    { 
        this->m_unknownDim.resize(2);
        this->m_unknownDim[0] = m_domain.dim();
        this->m_unknownDim[1] = 1;
    }
    // COMPATIBILITY CONSTRUCTORS, DO NOT USE
    gsStokesPde( const gsFunction<T> &rhs, int domdim)
        :         m_rhs(rhs.clone()), m_viscosity(0.001)
    {
        this->m_unknownDim.resize(2);
        this->m_unknownDim[0] = domdim;
        this->m_unknownDim[1] = 1;
    }
    gsStokesPde( const gsFunction<T> &rhs, int domdim, const gsFunction<T> &sol)
        :         m_rhs(rhs.clone()), m_viscosity(0.001)
    {
        this->m_unknownDim.resize(2);
        this->m_unknownDim[0] = domdim;
        this->m_unknownDim[1] = 1;
        this->m_solution.push_back(new gsPiecewiseFunction<T>(sol));
    }

    ~gsStokesPde( ) 
    { 
        delete m_rhs;
    }

    gsFunction<T> * rhs() const         { return m_rhs; }
    T viscocity() const                 { return m_viscosity; }

    /// Consistency check
    bool check() 
    {
        return true;
    }

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os<<"Stokes's equation:\n"
          <<"-\u0394u-\u2207p = f,\n"
          <<" \u2207\u00B7u=0"
          <<"with:\n";
	    os<<"Source function f= "<< * m_rhs <<".\n";
        if ( this->solutionGiven(0) )
            os<<"Exact solution u = "<< * this->m_solution[0] <<".\n";
        if ( this->solutionGiven(0) )
            os<<"Exact solution p = "<< * this->m_solution[1] <<".\n";
	    return os; 
	}

protected:
    gsFunction<T> * m_rhs;
    T m_viscosity;
}; // class gsStokesPde

} // namespace gismo
