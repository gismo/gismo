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

/** \brief Constructs the assembler for the discretized isogeometric heat equation.
    
    Spatial solution is discretized using Galerkin's approach. Time
    integration scheme (method of lines) is controlled by the \a theta
    parameter (to be set in options)t. In particular
    
    - Explicit Euler scheme (theta=0)
    - Crank-Nicolson semi-implicit scheme (theta=0.5)
    - implicit Euler scheme (theta=1)
    
    \ingroup Assembler
*/
template <class T>
class gsHeatEquation : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:

    /// Construction receiving all necessary data
    gsHeatEquation(gsAssembler<T> & stationary,
                   const gsOptionList & opt = Base::defaultOptions() )
    :  Base(stationary),  // note: unnecessary sliced copy here
       m_stationary(&stationary)
    {
        m_options.addReal("theta",
        "Theta parameter determining the time integration scheme[0..1]", 0.5);
    }

public:
    
    /// Initial assembly routine.
    void assemble()
    {
        // Grab theta once and for all
        m_theta = m_options.getReal("theta");
        
        // Assemble the stationary problem
        m_stationary->assemble();

        //copy the Dirichlet values, to enable calling
        // the construct solution functions
        Base::m_ddof = m_stationary->allFixedDofs();
        
        // Assemble mass matrix
        assembleMass();

        GISMO_ASSERT( m_stationary->matrix().rows() == m_mass.rows(),
                      "Something went terribly wrong.");
    }

    /** \brief Computes the matrix and right-hand side for the next timestep.

       \param curSolution The solution of the previous timestep

       \param[in,out] curRhs Input is the right-hand side of the
       previous timestep. It is overwritten with the right-hand side
       of the current timestep

       \param Dt Length of time interval of the current time step
    */
    void nextTimeStep(const gsMatrix<T> & curSolution, gsMatrix<T> & curRhs, T Dt);

protected:

    /// Mass assembly routine
    void assembleMass();


protected:

    using Base::m_options;

    /// The stationary system is stored here
    using Base::m_system;

    gsAssembler<T> * m_stationary;
    
    /// The mass matrix
    gsSparseMatrix<T> m_mass;
    
    /// Theta parameter determining the scheme
    T m_theta;
    
    using Base::m_pde_ptr;
    using Base::m_bases;

};// end class definition


} // namespace gismo


namespace gismo
{
template<class T>
void gsHeatEquation<T>::nextTimeStep(const gsMatrix<T> & curSolution, 
                                     gsMatrix<T> & curRhs, const T Dt)
{
    GISMO_ASSERT( curSolution.rows() == m_mass.cols(),
                  "Wrong size in current solution vector.");

    const T c1 = Dt * m_theta;
    m_system.matrix() = m_mass + c1 * m_stationary->matrix();

    const T c2 = Dt * (1.0 - m_theta);
    // note: noalias() still works since curRhs is multiplied by scalar only
    curRhs.noalias() = c1 * m_stationary->rhs() + c2 * curRhs + 
        m_mass * curSolution - c2 * m_stationary->matrix() * curSolution;
}

template<class T>
void gsHeatEquation<T>::assembleMass()
{
    // the mass visitor used to need empty variable
    // Base::m_ddof.resize(1, gsMatrix<T>() );

    const index_t m_dofs = this->numDofs();
    if ( 0  == m_dofs ) // Are there any interior dofs ?
    {
        gsWarn<<"No interior DoFs for mass compuation.\n";
        return;
    }

    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    m_system.reserve(m_bases[0], m_options, 0);// zero rhs's

    // Assemble mass integrals
    gsVisitorMass<T> mass;
    for (unsigned np=0; np < m_pde_ptr->domain().nPatches(); ++np )
    {
        //Assemble mass matrix for this patch
        this->apply(mass, np);
    }

    // Assembly is done, compress the matrix
    this->finalize();

    // Store the mass matrix once and for all
    m_system.matrix().swap(m_mass);
}

} // namespace gismo
