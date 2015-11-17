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
class gsHeatEquation : public gsPoissonAssembler2<T>
{
public:
    typedef gsPoissonAssembler2<T> Base;

public:

    gsHeatEquation( gsMultiPatch<T> const         & patches,
                    gsMultiBasis<T> const         & bases,
                    gsBoundaryConditions<T> const & bconditions,
                    const gsFunction<T>           & rhs,
                    const T                       _theta,
                    dirichlet::strategy           dirStrategy,
                    iFace::strategy               intStrategy = iFace::glue)
    :  Base(patches,bases,bconditions,rhs,dirStrategy,intStrategy)
    {
        m_theta = _theta ;
    }


public:

    
    /// Initial assembly routine.
    void assemble()
    {
        // Assemble the poisson problem
        Base::assemble();

        // Store the stiffness matrix
        m_matrix.swap(m_stiffMat);

        // Assemble mass matrix
        assembleMass();

        // Initializing rhs to zero
        m_curRhs.setZero(m_dofs,1);
    }
    
    //void nextTimeStep(const gsMatrix<T> & curSolution, T Dt);
    
    void nextTimeStep(const gsMatrix<T> & curSolution, const gsMatrix<T> & curRhs, T Dt);

    const gsMatrix<T> & curRhs() { return m_curRhs; }

    const gsSparseMatrix<T> & matrix() { return m_curMat; }

    const gsMatrix<T> &       rhs() { return m_curRhs; }

protected:

    /// Main assembly routine.
    void assembleMass();


protected:

    gsMatrix<T> m_curRhs;

    gsSparseMatrix<T> m_curMat;

    gsSparseMatrix<T> m_stiffMat;

    T m_theta;
    using Base::m_rhsFun;
    using Base::m_bConditions;
    

    // Members from gsAssemblerBase
    using Base::m_patches;
    using Base::m_bases;
    using Base::m_dofMappers;
    using Base::m_ddof;
    using Base::m_matrix;
    using Base::m_rhs;
    using Base::m_dofs;


};// end class definition


} // namespace gismo



//////////////////////////////////////////////////
//////////////////////////////////////////////////


namespace gismo
{
template<class T>
void gsHeatEquation<T>::nextTimeStep(const gsMatrix<T> & curSolution, 
                                     const gsMatrix<T> & curRhs, T Dt)
{
    GISMO_ASSERT( curSolution.rows() == m_matrix.rows(), "Wrong size in current solution vector.");
  
    // Theta-Scheme
    m_curMat           = m_matrix + Dt * m_theta * m_stiffMat;
    
    m_curRhs.noalias() = Dt * m_theta * m_rhs + m_matrix * curSolution
                       + Dt * (1.0 - m_theta) * curRhs 
                       - Dt * (1.0 - m_theta) * m_stiffMat * curSolution;    
}

template<class T>
void gsHeatEquation<T>::assembleMass()
{
    if (m_dofs == 0 ) // Are there any interior dofs ?
    {
        return;
    }

    // Pre-allocate non-zero elements for each column of the
    // sparse matrix
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_bases.front().dim(); ++i) // to do: improve
        nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1; // need more for DG !

    m_matrix = gsSparseMatrix<T>(m_dofs, m_dofs); // Clean matrices
    m_matrix.reserve( gsVector<int>::Constant(m_dofs, nonZerosPerCol) );
    
    // Assemble mass integrals
    gsVisitorMass<T> mass;
    for (unsigned np=0; np < m_patches.nPatches(); ++np )
    {
        //Assemble mass matrix for this patch
        this->apply(mass, np);
    }

    // Assembly is done, compress the matrix
    m_matrix.makeCompressed();  
}

} // namespace gismo
