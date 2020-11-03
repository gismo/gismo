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
    explicit gsHeatEquation(gsAssembler<T> & stationary)
    :  Base(stationary),  // note: unnecessary sliced copy here
       m_stationary(&stationary), m_theta(0.5)
    {
        m_options.addReal("theta",
        "Theta parameter determining the time integration scheme[0..1]", m_theta);
    }

public:

    void setTheta(const T th)
    {
        GISMO_ASSERT(th<=1 && th>=0, "Invalid value");
        m_theta= th;
        m_options.setReal("theta", m_theta);
    }
    
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

        The right-hand side function is assumed constant with respect to time

       \param curSolution The solution of the previous timestep

       \param Dt Length of time interval of the current time step
    */
    void nextTimeStep(const gsMatrix<T> & curSolution, const T Dt);
    
    void nextTimeStep(const gsSparseMatrix<T> & sysMatrix,
                      const gsSparseMatrix<T> & massMatrix,
                      const gsMatrix<T> & rhs0,
                      const gsMatrix<T> & rhs1,                                  
                      const gsMatrix<T> & curSolution,
                      const T Dt);
    
    void nextTimeStepFixedRhs(const gsSparseMatrix<T> & sysMatrix,
                              const gsSparseMatrix<T> & massMatrix,
                              const gsMatrix<T> & rhs,
                              const gsMatrix<T> & curSolution,
                              const T Dt);

    const gsSparseMatrix<T> & mass() const { return m_mass; }
    const gsSparseMatrix<T> & stationaryMatrix() const { return m_stationary->matrix(); }
    const gsSparseMatrix<T> & stationaryRhs() const { return m_stationary->rhs(); }
    
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
void gsHeatEquation<T>::nextTimeStep(const gsMatrix<T> & curSolution, const T Dt)
{
    nextTimeStepFixedRhs(m_stationary->matrix(), m_mass,
                         m_stationary->rhs(), curSolution, Dt);
}

/*
template<class T>
void gsHeatEquation<T>::nextTimeStep(const gsMatrix<T> & curSolution, const T Dt)
{
    // Keep previous time
    m_prevRhs.swap(m_stationary->rhs());
    // Get current time
    m_time += Dt;
    m_stationary->rhs(m_time);

    // note: m_stationary->matrix() assumed constant in time
        
    nextTimeStep(m_stationary->matrix(), m_mass,
                 m_prevRhs, m_stationary->rhs(), curSolution, Dt);
    m_prevRhs.swap(m_stationary->rhs());
}
*/

template<class T>
void gsHeatEquation<T>::nextTimeStep(const gsSparseMatrix<T> & sysMatrix,
                                     const gsSparseMatrix<T> & massMatrix,
                                     const gsMatrix<T> & rhs0,
                                     const gsMatrix<T> & rhs1,
                                     const gsMatrix<T> & curSolution,
                                     const T Dt)
{
    GISMO_ASSERT( curSolution.rows() == massMatrix.cols(),
                  "Wrong size in current solution vector.");

    const T c1 = Dt * m_theta;
    m_system.matrix() = massMatrix + c1 * sysMatrix;

    const T c2 = Dt * (1.0 - m_theta);
    m_system.rhs().noalias() = c1 * rhs1 + c2 * rhs0 + (massMatrix - c2 * sysMatrix) * curSolution;
}

template<class T>
void gsHeatEquation<T>::nextTimeStepFixedRhs(const gsSparseMatrix<T> & sysMatrix,
                                             const gsSparseMatrix<T> & massMatrix,
                                             const gsMatrix<T> & rhs,
                                             const gsMatrix<T> & curSolution,
                                             const T Dt)
{
    GISMO_ASSERT( curSolution.rows() == massMatrix.cols(),
                  "Wrong size in current solution vector.");

    const T c1 = Dt * m_theta;
    m_system.matrix() = massMatrix + c1 * sysMatrix;

    const T c2 = Dt * (1.0 - m_theta);
    m_system.rhs().noalias() = Dt * rhs + (massMatrix - c2 * sysMatrix) * curSolution;
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
    for (size_t np=0; np < m_pde_ptr->domain().nPatches(); ++np )
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
