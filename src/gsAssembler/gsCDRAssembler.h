/** @file gsCDRAssembler.h

    @brief Provides assembler and solver for the Poisson equation, incl. adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris
*/

#pragma once

#include <gsAssembler/gsPoissonAssembler.h>
#include <gsAssembler/gsAssemblerOptions.h>

#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals
#include <gsAssembler/gsVisitorCDR.h> //
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals
//#include <gsAssembler/gsVisitorDg.h>      // Disc. Galerkin interface integrals

#include <gsPde/gsPoissonPde.h>

//#include <gsAssembler/gsAssemblerUtils.h>

namespace gismo
{

/** @brief
    Implementation of an (multiple righ-hand side) Poisson solver.

    The Poisson equation: \f$-\Delta\mathbf{u}=\mathbf{f} \f$

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).
*/
template<class T>
class gsCDRAssembler : public gsPoissonAssembler<T>
{
public:
    typedef gsPoissonAssembler<T> Base;

public:

    /** @brief Main Constructor of the assembler object.
        \param[in] pde A boundary value Poisson problem
        \param[in] bases a multi-basis that contains patch-wise bases
        \param[in] flagStabilization (default noStabilization)
    */
    gsCDRAssembler(const gsConvDiffRePde<T> & pde,
                   const gsMultiBasis<T>    & bases,
                   stabilizerCDR::method      flagStabilization = stabilizerCDR::none)
    {
        // enrich options in constructor, refresh to apply options

        // 0: no stabilization
        // 1: SUPG
        m_options.addInt("Stabilization", "Choice of stabilization method; 0 := no; 1 := SUPG;", flagStabilization);
        Base::initialize(pde, bases, m_options);
    }

    /** @brief Constructor of the assembler object.
        \param[in] pde A boundary value Poisson problem
        \param[in] bases a multi-basis that contains patch-wise bases
        \param[in] dirStrategy option for the treatment of Dirichlet boundary
        \param[in] intStrategy option for the treatment of patch interfaces
        \param[in] flagStabilization (default noStabilization)
    */
    gsCDRAssembler(const gsConvDiffRePde<T> & pde,
                   const gsMultiBasis<T>    & bases,
                   dirichlet::strategy        dirStrategy,
                   iFace::strategy            intStrategy = iFace::glue,
                   stabilizerCDR::method      flagStabilization = stabilizerCDR::none)
    {
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        // 0: no stabilization
        // 1: SUPG
        m_options.addInt("Stabilization", "Choice of stabilization method; 0 := no; 1 := SUPG;", flagStabilization);

        Base::initialize(pde, bases, m_options);
    }

    /** @brief Constructor of the assembler object.
        \param[in] patches is a gsMultiPatch object describing the geometry.
        \param[in] bases a multi-basis that contains patch-wise bases
        \param[in] bconditions is a gsBoundaryConditions object that holds all boundary conditions.
        \param[in] rhs is the right-hand side of the equation, \f$\mathbf{f}\f$.
        \param[in] coeff_A diffusion coefficient.
        \param[in] coeff_b convection velocity.
        \param[in] coeff_c reaction coefficient.
        \param[in] dirStrategy option for the treatment of Dirichlet boundary
        \param[in] intStrategy option for the treatment of patch interfaces
        \param[in] flagStabilization (default noStabilization)
     */
    gsCDRAssembler(gsMultiPatch<T> const         & patches,
                   gsMultiBasis<T> const         & bases,
                   gsBoundaryConditions<T> const & bconditions,
                   const gsFunction<T>           & rhs,
                   const gsFunction<T>           & coeff_A,
                   const gsFunction<T>           & coeff_b,
                   const gsFunction<T>           & coeff_c,
                   dirichlet::strategy             dirStrategy = dirichlet::elimination,
                   iFace::strategy                 intStrategy = iFace::glue,
                   stabilizerCDR::method           flagStabilization = stabilizerCDR::none)
    {
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        // 0: no stabilization
        // 1: SUPG
        m_options.addInt("Stabilization", "Choice of stabilization method; 0 := no; 1 := SUPG;", flagStabilization);

        typename gsPde<T>::Ptr pde(new gsConvDiffRePde<T>
                                       (patches, bconditions, &coeff_A, &coeff_b, &coeff_c, &rhs));
        Base::initialize(pde, bases, m_options);
    }

    /// Main assembly routine
    void assemble()
    {
        GISMO_ASSERT(m_system.initialized(),
                     "Sparse system is not initialized, call initialize() or refresh()");

        // Reserve sparse system
        m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        Base::computeDirichletDofs();

        if (0 == this->numDofs()) // Are there any interior dofs ?
        {
            gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" << "\n";
            return;
        }

        // Assemble volume integrals
        Base::template push<gsVisitorCDR<T> >();

        // Enforce Neumann boundary conditions
        Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides());

        // If requested, enforce Dirichlet boundary conditions by Nitsche's method
        if (m_options.getInt("DirichletStrategy") == dirichlet::nitsche)
            Base::template push<gsVisitorNitsche<T> >(m_pde_ptr->bc().dirichletSides());

            // If requested, enforce Dirichlet boundary conditions by diagonal penalization
        else if (m_options.getInt("DirichletStrategy") == dirichlet::penalize)
            Base::penalizeDirichletDofs();

        // If we are in in dg (Discontinuous Galerkin) mode: add
        // interface contributions
        if (m_options.getInt("InterfaceStrategy") == iFace::dg)
            gsWarn << "DG option is ignored.\n";

        // Assembly is done, compress the matrix
        Base::finalize();
    }

protected:

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_bases;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
};

} // namespace gismo
