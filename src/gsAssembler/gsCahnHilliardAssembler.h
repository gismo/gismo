/** @file gsCahnHilliardAssembler.h

    @brief Provides assembler for a (planar) Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
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
     * @brief Assembles the mass matrix separately
     */
    void assembleMassMatrix();

    /**
     * @brief Assembles the stiffness matrix separately
     * @param C The solution
     * @param DC The time-derivative of the solution
     */
    void assembleJacobian(const gsMultiPatch<T> & C, const gsMultiPatch<T> & DC);

    /**
     * @brief Assembles the residual
     * @param C The solution
     * @param DC The time-derivative of the solution
     */
    void assembleResidual(const gsMultiPatch<T> & C, const gsMultiPatch<T> & DC);

    /**
     * @brief Assembles the Nitsche vector for boundary conditions separately
     * @param C The solution
     * @param DC The time-derivative of the solution
     */
    void assembleNitscheVector(const gsMultiPatch<T> & C, const gsMultiPatch<T> & DC);

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
    void constructSolution(gsMatrix<T>     & Cvec,
                           gsMultiPatch<T> & C);

    /**
     * @brief Constructs a spline solution from a solution vector
     * @param Cvec The solution vector
     * @param C The solution
     */
    void constructSolution(gsMatrix<T>         & Cvec,
                           gsMappedSpline<2,T> & C);

    /**
     * @brief Constructs a solution vector from a multi-patch solution
     * @param C The solution
     * @param Cvec The solution vector
     */
    void constructSolution(const gsMultiPatch<T> & C,
                                 gsMatrix<T>     & Cvec);

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
    mutable gsMultiBasis<T>   m_basis;
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
