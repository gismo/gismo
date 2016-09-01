/** @file gsHeatEquation.h

    @brief Provides assembler for the heat equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Moore, A. Mantzaflaris
*/

#pragma once


namespace gismo
{

template <class T>
class gsHeatEquation : public gsPoissonAssembler<T>
{
public:
    typedef gsPoissonAssembler<T> Base;

public:

    /** \brief Constructs the assembler for the discretized isogeometric heat equation.

        Spatial solution is discretized using Galerkin's
        approach. Time integration scheme (method of lines) is
        controlled by the \a theta parameter. In particular
        
        - Explicit Euler scheme (theta=0)
        - Crank-Nicolson semi-implicit scheme (theta=0.5)
        - implicit Euler scheme (theta=1)
        
        \ingroup Assembler
     */
    gsHeatEquation( gsMultiPatch<T> const         & patches,
                    gsMultiBasis<T> const         & bases,
                    gsBoundaryConditions<T> const & bconditions,
                    const gsFunction<T>           & rhs,
                    const T                       _theta,
                    dirichlet::strategy           dirStrategy,
                    iFace::strategy               intStrategy = iFace::glue)
    :  Base(patches,bases,bconditions,rhs,dirStrategy,intStrategy), m_theta(_theta)
    { }


public:

    
    /// Initial assembly routine.
    void assemble()
    {
        // Assemble the poisson problem
        Base::assemble();

        // Store the stiffness matrix
        m_system.matrix().swap(m_stiffMat);
        
        // Assemble mass matrix
        assembleMass();
    }

    /** \brief Computes the matrix and right-hand side for the next timestep.

       \param curSolution The solution of the previous timestep

       \param[in,out] curRhs Input is the right-hand side of the
       previous timestep. It is overwritten with the right-hand side
       of the current timestep

       \param Dt Length of time interval of the current time step
    */
    void nextTimeStep(const gsMatrix<T> & curSolution, gsMatrix<T> & curRhs, T Dt);

    /// Returns the system matrix for the lastly computed time step.
    const gsSparseMatrix<T> & matrix() { return m_curMat; }

protected:

    /// Mass assembly routine
    void assembleMass();


protected:

    /// The mass matrix is stored here, as well as the right-hand side
    /// (moment) vector corresponding to the source function
    using Base::m_system;

    /// The stiffness matrix
    gsSparseMatrix<T> m_stiffMat;

    /// The system matrix at the current timestep
    gsSparseMatrix<T> m_curMat;
    
    /// Theta parameter determining the scheme
    T m_theta;
    
    using Base::m_pde_ptr;
    using Base::m_bases;
    //using Base::m_ddof;

};// end class definition


} // namespace gismo


namespace gismo
{
template<class T>
void gsHeatEquation<T>::nextTimeStep(const gsMatrix<T> & curSolution, 
                                     gsMatrix<T> & curRhs, const T Dt)
{
    GISMO_ASSERT( curSolution.rows() == m_system.matrix().rows(),
                  "Wrong size in current solution vector.");
  
    const T c1 = Dt * m_theta;
    m_curMat = m_system.matrix() + c1 * m_stiffMat;

    const T c2 = Dt * (1.0 - m_theta);
    // note: noalias() still works since curRhs is multiplied by scalar only
    curRhs.noalias() = c1 * m_system.rhs() + c2 * curRhs + 
        m_system.matrix() * curSolution - c2 * m_stiffMat * curSolution;    
}

template<class T>
void gsHeatEquation<T>::assembleMass()
{
    const index_t m_dofs = this->numDofs();
    if ( 0  == m_dofs ) // Are there any interior dofs ?
    {
        gsWarn<<"No interior DogFs for mass compuation.\n";
        return;
    }

    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_bases.front().dim(); ++i)
        nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;

    gsSparseMatrix<T> & m_matrix = m_system.matrix();
    m_matrix.resize(m_dofs, m_dofs);
    m_matrix.reservePerColumn(nonZerosPerCol);

    // Assemble mass integrals
    gsVisitorMass<T> mass;
    for (unsigned np=0; np < m_pde_ptr->domain().nPatches(); ++np )
    {
        //Assemble mass matrix for this patch
        this->apply(mass, np);
    }

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();  
}

} // namespace gismo
