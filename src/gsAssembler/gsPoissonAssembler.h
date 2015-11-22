 /** @file gsPoissonAssembler.h

    @brief Provides assembler for the Poisson equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsPde/gsPoissonPde.h>


namespace gismo
{

/** @brief
    Implementation of an (multiple right-hand side) Poisson assembler.

    The Poisson equation: \f$-\Delta\mathbf{u}=\mathbf{f} \f$

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).
    
    \ingroup Assembler
*/
template <class T>
class gsPoissonAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:

    gsPoissonAssembler()
    { }

    /** @brief Main Constructor of the assembler object.

    \param[in] pde A boundary value Poisson problem
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
    */
    gsPoissonAssembler( gsPoissonPde<T> const         & pde,
                        gsMultiBasis<T> const         & bases,
                        dirichlet::strategy           dirStrategy,
                        iFace::strategy               intStrategy = iFace::glue)
    {
        m_options.dirStrategy = dirStrategy;
        m_options.intStrategy = intStrategy;

        Base::initialize(pde, bases, m_options);
    }

    /** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] basis a multi-basis that contains patch-wise bases
    \param[in] bconditions is a gsBoundaryConditions object that holds all boundary conditions.
    \param[in] rhs is the right-hand side of the Poisson equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
    
    \ingroup Assembler
    */
    gsPoissonAssembler( gsMultiPatch<T> const         & patches,
                        gsMultiBasis<T> const         & basis,
                        gsBoundaryConditions<T> const & bconditions,
                        const gsFunction<T>           & rhs,
                        dirichlet::strategy           dirStrategy = dirichlet::elimination,
                        iFace::strategy               intStrategy = iFace::glue)
    : m_ppde(patches,bconditions,rhs)
    {
        m_options.dirStrategy = dirStrategy;
        m_options.intStrategy = intStrategy;
        
        Base::initialize(m_ppde, basis, m_options);
    }

    // Refresh routine
    void refresh();

    // Main assembly routine
    void assemble();

    /// Returns an expression of the "full" assembled sparse
    /// matrix. Note that matrix() might return a lower diagonal
    /// matrix, if we exploit possible symmetry during assembly
    /// (check: m_matrix.symmetry() == true )
    Eigen::SparseSelfAdjointView< typename gsSparseMatrix<T>::Base, Lower> fullMatrix()
    {
        return m_system.matrix().template selfadjointView<Lower>();
    }

protected:

    // fixme: remove this
    gsPoissonPde<T> m_ppde;

protected:

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
};


} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPoissonAssembler.hpp)
#endif
