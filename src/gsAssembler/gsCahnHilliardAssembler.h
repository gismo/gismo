/** @file gsCahnHilliardAssembler.h

    @brief Provides assembler for a (planar) Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    -----------------------------------------------------------------------
    TODO;
    - Change hmax to a gsExprAssembler<>::element el; el.diam();
    -----------------------------------------------------------------------

    Author(s): H.M. Verhelst, L. Venta Viñuela
*/

#pragma once

#include <gsAssembler/gsExprAssembler.h>
#include <gsAssembler/gsExprEvaluator.h>

namespace gismo
{

template <class T>
class gsCahnHilliardAssembler
{
public:

/**
 * @brief      Constructs a new biharmonic assembler
 *
 * @param[in]  mp     A multi-patch
 * @param[in]  mb     A multi-basis
 * @param[in]  force  The force
 * @param[in]  bcs    The bcs
 */
gsCahnHilliardAssembler(const gsMultiPatch<T> & mp,
                        const gsMultiBasis<T> & mb,
                        const gsBoundaryConditions<T> & bcs
                            );

/// Default empty constructor
gsCahnHilliardAssembler() { }

/// Copy constructor (makes deep copy)
gsCahnHilliardAssembler( const gsCahnHilliardAssembler& other )
{
    operator=(other);
}

/// Move constructor
gsCahnHilliardAssembler( gsCahnHilliardAssembler&& other )
{
    operator=(give(other));
}

/// Assignment operator
gsCahnHilliardAssembler& operator= ( const gsCahnHilliardAssembler& other );

/// Move assignment operator
gsCahnHilliardAssembler& operator= ( gsCahnHilliardAssembler&& other );

protected:

    void _defaultOptions();

    void _getOptions();

public:

    void initialize();

    /**
      * @brief Overwrites the basis to be used as a space,
      *        and keeps the basis defined in the constructor
      *        for integration
      * @note  This option is usually called for assembly on \ref gsMappedBasis
      * @param spaceBasis The basis to be used for the space
      */
    void setSpaceBasis(const gsFunctionSet<T> & spaceBasis);

    /**
     * @brief Overwrites the basis to be used for integration,
     *        and keeps the basis defined in the constructor
     *        for the space
     * @note  This option is usually called for assembly on \ref gsMappedBasis
     * @param integrationBasis The basis to be used for integration
     */
    void setIntegrationBasis(const gsMultiBasis<T> & integrationBasis);

    /**
     * @brief Assembles the mass matrix separately
     */
    void assembleMassMatrix();

    /**
     * @brief Assembles the full system
     * @param C The solution (as a gsFunctionSet, e.g. gsMultiPatch)
     */
    void assemble(const gsFunctionSet<T> & C);

    /**
     * @brief Assembles the full system from coefficient vectors.
     *        Calls constructSolution internally so the gsMultiPatch stays on the C++ stack.
     * @param Cvec The solution coefficient vector
     */
    void assemble(const gsMatrix<T> & Cvec);

    /**
     * @brief Assembles the stiffness matrix (Jacobian) separately.
     *        Does not depend on the time-derivative of the solution.
     * @param C The solution (as a gsFunctionSet, e.g. gsMultiPatch)
     */
    void assembleJacobian(const gsFunctionSet<T> & C);

    /**
     * @brief Assembles the stiffness matrix (Jacobian) from a coefficient vector.
     *        Calls constructSolution internally so the gsMultiPatch stays on the C++ stack.
     * @param Cvec The solution coefficient vector
     */
    void assembleJacobian(const gsMatrix<T> & Cvec);

    /**
     * @brief Assembles the spatial (static) residual — the nonlinear spatial terms only.
     *        The time-derivative contribution M*dC must be added by the caller:
     *            Q = assembler.rhs() + M * dC
     * @param C The solution (as a gsFunctionSet, e.g. gsMultiPatch)
     */
    void assembleResidual(const gsFunctionSet<T> & C);

    /**
     * @brief Assembles the spatial residual from a coefficient vector.
     *        Calls constructSolution internally so the gsMultiPatch stays on the C++ stack.
     * @param Cvec The solution coefficient vector
     */
    void assembleResidual(const gsMatrix<T> & Cvec);

    /**
     * @brief Assembles the Nitsche vector for boundary conditions separately.
     *        Does not depend on the time-derivative of the solution.
     * @param C The solution (as a gsFunctionSet, e.g. gsMultiPatch)
     */
    void assembleNitscheVector(const gsFunctionSet<T> & C);

    /**
     * @brief Assembles the Nitsche vector from a coefficient vector.
     *        Calls constructSolution internally so the gsMultiPatch stays on the C++ stack.
     * @param Cvec The solution coefficient vector
     */
    void assembleNitscheVector(const gsMatrix<T> & Cvec);

    /**
     * @brief Assembles the Nitsche matrix for boundary conditions separately
     * @note  This term does not depend on the solution, hence it can be assembled once
     */
    void assembleNitscheMatrix();

    /**
     * @brief Returns a handle to the latest assembled matrix
     * @return The matrix
     */
    const gsSparseMatrix<T> & matrix() const { return m_assembler.matrix();  }

    /**
     * @brief Moves the latest assembled matrix from the assembler to the output
     * @param out The output matrix
     */
    void matrix_into(gsSparseMatrix<T> & out) { m_assembler.matrix_into(out); }

    /**
     * @brief Returns a handle to the latest assembled right-hand side
     * @return The right-hand side
     */
    const gsMatrix<T>       & rhs()    const { return m_assembler.rhs();     }

    /**
     * @brief Moves the latest assembled right-hand side from the assembler to the output
     * @param out The output right-hand side
     */
    void rhs_into(gsMatrix<T> & out) { m_assembler.rhs_into(out); }

    /**
     * @brief Returns the number of degrees of freedom
     * @return The number of degrees of freedom
     */
    index_t numDofs() const { return m_assembler.numDofs(); };

    /// Returns the finalized trial-space DoF mapper used by this assembler.
    const gsDofMapper & dofMapper() const
    { return m_assembler.trialSpace(0).mapper(); }

    /// Returns values for eliminated trial-space DoFs.
    const gsMatrix<T> & fixedPart() const
    { return m_assembler.trialSpace(0).fixedPart(); }

    /**
     * @brief Returns a handle to the options stored in the class
     * @return The options
     */
    gsOptionList & options() {return m_options;}

    /**
     *  @brief Set the options from an option list. Ignores unknown options
     */
    void setOptions(gsOptionList & options);

    /**
     * @brief Constructs a multi-patch solution from a solution vector
     * @param Cvec The solution vector
     * @param C The solution
     */
    void constructSolution(const gsMatrix<T> & Cvec,
                           gsMultiPatch<T>   & C) const;

    /**
     * @brief Constructs a spline solution from a solution vector
     * @param Cvec The solution vector
     * @param C The solution
     */
    void constructSolution(const gsMatrix<T>   & Cvec,
                           gsMappedSpline<2,T> & C) const;

    /**
      * @brief Constructs a solution vector from a multi-patch solution
      * @param C The solution
      * @param Cvec The solution vector
      */
    void constructSolution(const gsMultiPatch<T> & C,
                                 gsMatrix<T>     & Cvec) const;

    /**
     * @brief Computes the mass M(c) = ∫ c dx
     * @param C The solution (as a gsFunctionSet, e.g. gsMultiPatch)
     * @return The computed mass integral
     */
    T computeMass(const gsFunctionSet<T> & C);

    /**
     * @brief Computes the mass M(c) = ∫ c dx from a coefficient vector.
     *        Calls constructSolution internally so the gsMultiPatch stays on the C++ stack.
     * @param Cvec The solution coefficient vector
     * @return The computed mass integral
     */
    T computeMass(const gsMatrix<T> & Cvec);

    /**
     * @brief Computes the dissipation D(c) = ∫ |∇mu|^2 dx where mu = -c + c^3 - lambda Δc
     * @param C The solution (as a gsFunctionSet, e.g. gsMultiPatch)
     * @return The computed dissipation integral
     */
    T computeDissipation(const gsFunctionSet<T> & C);

    /**
     * @brief Computes the dissipation D(c) = ∫ |∇mu|^2 dx where mu = -c + c^3 - lambda Δc from a coefficient vector.
     *        Calls constructSolution internally so the gsMultiPatch stays on the C++ stack.
     * @param Cvec The solution coefficient vector
     * @return The computed dissipation integral
     */
    T computeDissipation(const gsMatrix<T> & Cvec);

    /**
     * @brief Computes the dissipation D(c) = ∫ M(c) |∇mu|^2 dx where mu = -c + c^3 - lambda Δc
     * @param C The solution (as a gsFunctionSet, e.g. gsMultiPatch)
     * @param mu The chemical potential (as a gsFunctionSet, e.g. gsMultiPatch)
     * @return The computed dissipation integral
     * @note This implementation requires a field for the chemical potential mu, avoiding the requirement for
     */
    T computeDissipation(const gsFunctionSet<T> & C, const gsFunctionSet<T> & mu);

    /**
     * @brief Computes the dissipation D(c) = ∫ M(c) |∇mu|^2 dx where mu = -c + c^3 - lambda Δc from coefficient vectors.
     *        Calls constructSolution internally so the gsMultiPatch stays on the C++ stack.
     * @param Cvec The solution coefficient vector
     * @param muVec The chemical potential coefficient vector
     * @return The computed dissipation integral
     * @note This implementation requires a field for the chemical potential mu, avoiding the requirement for
     */
    T computeDissipation(const gsMatrix<T> & Cvec, const gsMatrix<T> & muVec);

    /**
     * @brief Computes the energy E(c) = ∫ [ 1/4 (c^2 - 1)^2 + lambda/2 |∇c|^2 ] dx
     * @param C The solution (as a gsFunctionSet, e.g. gsMultiPatch)
     * @return The computed energy integral
     */
    T computeEnergy(const gsFunctionSet<T> & C);
    /**
     * @brief Computes the energy E(c) = ∫ [ 1/4 (c^2 - 1)^2 + lambda/2 |∇c|^2 ] dx from a coefficient vector.
     *        Calls constructSolution internally so the gsMultiPatch stays on the C++ stack.
     * @param Cvec The solution coefficient vector
     * @return The computed energy integral
     */
    T computeEnergy(const gsMatrix<T> & Cvec);

protected:

    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::space       space;
    typedef typename gsExprAssembler<T>::solution    solution;
    typedef typename gsExprAssembler<T>::element     element;

    mutable index_t m_continuity;
    mutable T m_lambda;
    // mutable T m_M0;
    mutable T m_penalty;

    gsExprAssembler<T> m_assembler;
    gsExprEvaluator<T> m_evaluator;

    gsMultiPatch<T>           m_patches;
    const gsMultiBasis<T> *   m_integrationBasis;
    const gsFunctionSet<T> *  m_spaceBasis;
    gsBoundaryConditions<T>   m_bcs;
    bool m_initialized;

    mutable gsOptionList m_options;

  }; // class gsCahnHilliardAssembler

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsCahnHilliardAssembler
   */
  void pybind11_init_gsCahnHilliardAssembler(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsCahnHilliardAssembler.hpp)
#endif
