/** \file gsPde.h
 *
 *  \brief base class of descriptions of a PDE problem.
 *  Due to the huge difference in data required to describe a Pde it only contains
 *  the domain and the definition of a common interface.
**/

#pragma once

#include <iostream>
#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsPiecewiseFunction.h>
#include <gsCore/gsMultiPatch.h>

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
**/
template<class T>
class gsPde
{
protected:
    gsPde()
    {
    }
public:  
    gsPde(const gsMultiPatch<T> &domain, const gsBoundaryConditions<T> &bc)
        : m_domain(domain), m_boundary_conditions(bc)
    {}
    gsPde(const gsMultiPatch<T>         &domain,
          const gsBoundaryConditions<T> &bc,
          const std::vector<gsPiecewiseFunction<T> *>  &solutions
          )
        : m_domain(domain), m_boundary_conditions(bc), m_solution(solutions)
    {}

    virtual ~gsPde() 
    { 
        freeAll(m_solution);
    }
    
    /**
     * @brief returns a reference to the Pde domain
     *        there is also a const variant returning a const reference
     * @return
     */
    gsMultiPatch<T> &domain() {return m_domain;}
    const gsMultiPatch<T> &domain() const {return m_domain;}

    /**
     * @brief returns a reference to the Pde boundary conditions
     *        there is also a const variant returning a const reference
     * @return
     */
    gsBoundaryConditions<T> &boundaryConditions() {return m_boundary_conditions;}
    const gsBoundaryConditions<T> &boundaryConditions() const {return m_boundary_conditions;}

    /// Is the associated linear system symmetric?
    /// TODO: remove because we cannot know, it depends on the method
    GS_DEPRECATED virtual bool isSymmetric() const {gsWarn<<"Function is gsPde::isSymmetric should not be used!!"; return false;}

    /// Print a short description of the PDE
    virtual std::ostream &print(std::ostream &os) const = 0;

    /**
     * @brief solutionGiven
     * @param field_id
     * @return
     */
    bool solutionGiven(index_t field_id = 0) const
    {
        GISMO_ASSERT(field_id<numUnknowns(),"Asked for presence of exact solution for a non existing field of this PDE type");
        return m_solution[field_id] != NULL;
    }
    /**
     * @brief gives the exact solution for the \em field_id-th field on the \em patch_id-th patch
     * @param field_id
     * @param patch_id
     * @return a pointer to the function or NULL depending if the exact solution is provided
     */
    gsFunction<T>* solution(index_t field_id = 0, index_t patch_id =0) const
    {
        GISMO_ASSERT(field_id<numUnknowns(),"Asked for the exact solution for a non existing field of this PDE type");
        return &(m_solution[field_id]->operator[](patch_id));
    }
    /**
     * @brief gives the exact solution for the \em field_id-th field
     * @param field_id
     * @return a pointer to the function or NULL depending if the exact solution is provided
     */
    const gsPiecewiseFunction<T>* solutionField(index_t field_id = 0) const
    {
        GISMO_ASSERT(field_id<numUnknowns(),"Asked for size of a non existing field of this PDE type");
        return m_solution[field_id];
    }
    /**
     * @brief gives the number of unknown fields of the PDEs
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
    int fieldDim(index_t field_id = 0)
    {
        GISMO_ASSERT(field_id<numUnknowns(),"Asked for size of an Unknown field for this PDE type");
        return m_unknownDim[field_id];
    }
    /**
     * @brief returns the dimension of the domain
     *
    **/
    int dim()
    {
        return m_domain.dim();
    }


protected:
    // description of the unknown fields:
    // for each one the target dimension.
    std::vector<unsigned>                 m_unknownDim;
    // exact solution
    std::vector<gsPiecewiseFunction<T> *> m_solution;
    // description of the problem
    gsMultiPatch<T>                       m_domain;
    gsBoundaryConditions<T>               m_boundary_conditions;
}; // class gsPde

/// Print (as string) operator to be used by all derived classes
template<class T>
std::ostream &operator<<(std::ostream &os, const gsPde<T>& pde)
{
    return pde.print(os);
}

} // namespace gismo
