/** @file gsBiharmonicPde.h

    @brief Describes a Poisson PDE.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/


#pragma once

#include <gsPde/gsPde.h>

namespace gismo
{

/** @brief
    A Biharmonic PDE.

    This class describes a Biharmonic PDE, with an arbitrary right-hand side
    function.

    It has an extra gsBoundaryConditions object since the biharmonic has
    two essential (Dirichlet) and two natural (Neumann). The second
    gsBoundaryConditions contains the second kind of BC.

    The first kind of BCs
    Dirichlet: Enforce \f$v = 0 \f$ on the boundary
    Neumann: Add \f$ (g,v)_\Gamma \f$ on the right-hand side,
    where \f$g = -\nabla \Delta u \cdot \mathbf{n}\f$

    The second kind of BCs
    Dirichlet: Enforce \f$ \nabla v \cdot \mathbf{n} = 0 \f$ on the boundary
    Neumann: Add \f$ (g,\nabla v \cdot \mathbf{n})_\Gamma \f$ on the right-hand side,
    where \f$g = \Delta u \f$

    \ingroup Pde
    \ingroup pdeclass
 */

template<class T>
class gsBiharmonicPde : public gsPde<T>
{
protected:
    gsBiharmonicPde( ) { }
    using gsPde<T>::m_domain;
public:
    /// Constructor
    gsBiharmonicPde(
        const gsMultiPatch<T>         &domain,
        const gsBoundaryConditions<T> &bc1,
        const gsBoundaryConditions<T> &bc2,
        const gsPiecewiseFunction<T>  &rhs
         )
    : gsPde<T>(domain,bc1), m_rhs(rhs), m_boundary_conditions_second(bc2)
    {
        m_unknownDim.setOnes(1);
    }


    ~gsBiharmonicPde( )
    { }

    /**
     * @brief gives the number of rhs functions of the PDEs
     */
    virtual int numRhs() const
    {
        return m_rhs.piece(0).targetDim();
    }

    const gsFunction<T> *    rhs()      const { return &m_rhs.piece(0); }

    const gsBoundaryConditions<T> & bcFirstKind()  const {return this->bc();}

    const gsBoundaryConditions<T> & bcSecondKind() const {return m_boundary_conditions_second;}

    gsBoundaryConditions<T> &boundaryConditionsSecond() {return m_boundary_conditions_second;}

    const gsBoundaryConditions<T> &boundaryConditionsSecond() const {return m_boundary_conditions_second;}

    /// Prints the object as a string.
    std::ostream &print(std::ostream &os) const
    {
        os<<"Biharmonic's equation  -\u0394^2 u = f ,  with:\n";
	    os<<"Source function f= "<< m_rhs.piece(0) <<".\n";
	    return os; 
	}
protected:
    using gsPde<T>::m_unknownDim;

    gsPiecewiseFunction<T> m_rhs;
    /// @brief Boundary conditions of the second kind
    gsBoundaryConditions<T> m_boundary_conditions_second;
}; // class gsBiharmonicPde

} // namespace gismo
