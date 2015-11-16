/** @file gsPde.h

    @brief Base class of descriptions of a PDE problem.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsPiecewiseFunction.h>

namespace gismo
{


/** @brief
    Abstract class representing a PDE (partial differential equation).

    This is the base class for PDEs, like the Poisson equation, Stokes
    equation or convection diffusion equation.
    The instances of these classes contains the problem description:
    - the domain,
    - the boundary conditions,
    - coefficients
    - rhs functions

    The description of the discretization: methods to impose conditions
    or discrete spaces must be provided separately.

    When available, an exact solution is kept as well.
    PDEs can be read and written to xml files.
    
    \ingroup Pde
    \ingroup pdeclass
**/
template<class T>
class gsPde
{
protected:
    gsPde()
    {
    }
public:
    /// @brief Constructor without given exact solution.
    gsPde(const gsMultiPatch<T> &domain, const gsBoundaryConditions<T> &bc)
        : m_domain(domain), m_boundary_conditions(bc)
    {}

    /// @brief Constructor with given exact solution.
    gsPde(const gsMultiPatch<T>               &domain,
          const gsBoundaryConditions<T>       &bc,
          const std::vector<gsFunction<T> *>  &solutions
          )
        : m_domain(domain), m_boundary_conditions(bc), m_solution(solutions)
    {}

    virtual ~gsPde() 
    { 
        freeAll(m_solution);
    }
    
    /**
     * @brief Returns a reference to the Pde domain.
     *
     * There is also a const version returning a const reference.
     */
    gsMultiPatch<T>       & domain() {return m_domain;}

    const gsMultiPatch<T> & domain() const {return m_domain;}

    gsMultiPatch<T>       & patches() {return m_domain;}

    const gsMultiPatch<T> & patches() const {return m_domain;}

    /**
     * @brief Returns a reference to the Pde boundary conditions.
     *
     * There is also a const version returning a const reference.
     */
    gsBoundaryConditions<T> &boundaryConditions() {return m_boundary_conditions;}

    const gsBoundaryConditions<T> &boundaryConditions() const {return m_boundary_conditions;}

    const gsBoundaryConditions<T> & bc() const {return m_boundary_conditions;}

    // Is the associated linear system symmetric?
    // TODO: Remove, because it depends on the method and the used
    // discretization and test spaces whether the resulting linear system
    // is symmetric or not. Hence, the isSymmetric()-Flag does not
    // make too much sense in the specification of the PDE.
    //
    // As of now, the function still remains, because it is called
    // by some other functions.
    virtual bool isSymmetric() const {gsWarn<<"Function is gsPde::isSymmetric should not be used!!"; return false;}

    /// @brief Print a short description of the PDE
    virtual std::ostream &print(std::ostream &os) const = 0;

    /**
     * @brief Indicates whether an exact solution is stored.
     *
     * @param field_id
     * @return true, if an exact solution
     * for the field with index \em field_id was provided.
     */
    bool solutionGiven(index_t field_id = 0) const
    {
        GISMO_ASSERT(field_id<numUnknowns(),"Asked for presence of exact solution for a non existing field of this PDE type");
        return m_solution[field_id] != NULL;
    }
    /**
     * @brief Gives the exact solution for the <em>field_id</em>-th field on the <em>patch_id</em>-th patch.
     *
     * @param field_id
     * @return a pointer to the function or NULL depending if the exact solution is provided
     */
    gsFunction<T>* solution(index_t field_id = 0) const
    {
        GISMO_ASSERT(field_id<numUnknowns(),"Asked for the exact solution for a non existing field of this PDE type");
        return (m_solution[field_id]);
    }

    /**
     * @brief solutions
     * @return a reference to the vector of solutions for each field
     */
    const std::vector<gsFunction<T>*> &solutions() const
    {
    return m_solution;
    }

    /**
     * @brief Gives the exact solution for the \em field_id-th field.
     *
     * @param field_id
     * @return a pointer to the function or NULL depending if the exact solution is provided
     */
    const gsFunction<T>* solutionField(index_t field_id = 0) const
    {
        GISMO_ASSERT(field_id<numUnknowns(),"Asked for size of a non existing field of this PDE type");
        return m_solution[field_id];
    }
    /**
     * @brief Gives the number of unknown fields of the PDEs.
     */
    int numUnknowns() const
    {
        return m_unknownDim.size();
    }
    /**
     * @brief gives the number of rhs functions of the PDEs
     */
    virtual int numRhs() const
    {
        return 1;
    }
    /**
     * @brief gives the dimension of the \em i-th field
     * it returns 1 for scalar fields, 2 for 2d vectors field etc.
     * @param field_id the field index
     */
    GISMO_DEPRECATED
    int fieldDim(index_t field_id = 0)
    {
        GISMO_ASSERT(field_id<numUnknowns(),"Asked for size of an Unknown field for this PDE type");
        return m_unknownDim[field_id];
    }

    /**
     * @brief returns the dimension of the domain
     *
    **/
    GISMO_DEPRECATED
    int dim() const
    {
        return m_domain.dim();
    }

protected:
    /// @brief Description of the unknown fields:
    /// for each one the target dimension.
    std::vector<unsigned>                 m_unknownDim;
    /// @brief Exact solution(s)
    std::vector<gsFunction<T> *>          m_solution;
    /// @brief Computational domain
    gsMultiPatch<T>                       m_domain;
    /// @brief Boundary conditions
    gsBoundaryConditions<T>               m_boundary_conditions;
}; // class gsPde

/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsPde<T>& pde)
{
    return pde.print(os);
}

} // namespace gismo
