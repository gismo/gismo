/** @file gsBiharmonicExprAssembler.h

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
class gsBiharmonicExprAssembler
{
public:

gsBiharmonicExprAssembler(  const gsMultiPatch<T> & mp,
                            const gsMultiBasis<T> & mb,
                            const gsFunctionSet<T> & force,
                            const gsBoundaryConditions<T> & bcs
                            );

/// Default empty constructor
gsBiharmonicExprAssembler() { }

/// Copy constructor (makes deep copy)
gsBiharmonicExprAssembler( const gsBiharmonicExprAssembler& other )
{
    operator=(other);
}

/// Move constructor
gsBiharmonicExprAssembler( gsBiharmonicExprAssembler&& other )
{
    operator=(give(other));
}

/// Assignment operator
gsBiharmonicExprAssembler& operator= ( const gsBiharmonicExprAssembler& other );

/// Move assignment operator
gsBiharmonicExprAssembler& operator= ( gsBiharmonicExprAssembler&& other );

protected:

    void _defaultOptions();

    void _getOptions();

    void _initialize();

    void _setup(const expr::gsFeSpace<T> & u);

public:

    /// See \ref gsThinShellAssemblerBase for details
    void setSpaceBasis(const gsFunctionSet<T> & spaceBasis);

    /// Assembles the mass matrix
    void assembleMass();

    /// Assembles the full system
    void assemble();

    /// Assembles the LHS
    void assembleLHS();

    /// Assembles the RHS
    void assembleRHS();

    const gsSparseMatrix<T> & matrix() const { return m_assembler.matrix();  }
    const gsMatrix<T>       & rhs()    const { return m_assembler.rhs();     }

    T penalty() const { return m_penalty; };
    index_t numDofs() const { return m_assembler.numDofs(); };

    /// Returns the L2 error given \a solVector and \a exact solution
    T l2error(gsMatrix<T> & solVector, const gsFunctionSet<T> & exact);

    /// Returns the H1 error given \a solVector and \a exact solution
    T h1error(gsMatrix<T> & solVector, const gsFunctionSet<T> & exact);

    /// Returns the H2 error given \a solVector and \a exact solution
    T h2error(gsMatrix<T> & solVector, const gsFunctionSet<T> & exact);

    /// Returns the L2, H1 and H2 errors given \a solVector and \a exact solution
    std::tuple<T,T,T> errors(gsMatrix<T> & solVector, const gsFunctionSet<T> & exact);

    /// Returns the H2 error given \a solVector and \a exact solution
    T interfaceError(gsMatrix<T> & solVector, const gsFunctionSet<T> & exact);

    /// Access a reference to the stored options.
    gsOptionList & options() {return m_options;}

    /// Set the options from an option list. Ignores unknown options
    void setOptions(gsOptionList & options);

    void constructSolution(gsMatrix<T> & solVector);

    typename gsFunctionSet<T>::Ptr getSolution() const;

private:

    // template <typename basisT>
    void _setMapperForBiharmonic(const gsBoundaryConditions<T> & bc,
                                const gsFunctionSet<T> & bb2,
                                gsDofMapper & mapper);

    void _getDirichletNeumannValuesL2Projection(const gsMultiPatch<T> & mp,
                                                const gsMultiBasis<T> & dbasis,
                                                const gsBoundaryConditions<T> & bc,
                                                const gsFunctionSet<T> & bb2,
                                                const expr::gsFeSpace<T> & u);

    void _computeStabilityParameter(const gsMultiPatch<T> & mp,
                                    const gsMultiBasis<T> & dbasis,
                                    gsMatrix<T> & mu_interfaces);

protected:

    typedef typename gsExprAssembler<T>::geometryMap geometryMap;
    typedef typename gsExprAssembler<T>::space       space;
    typedef typename gsExprAssembler<T>::solution    solution;
    typedef typename gsExprAssembler<T>::element     element;

    mutable T m_penalty;
    mutable T m_lambda;
    mutable bool m_second;
    mutable index_t m_continuity;

    gsExprAssembler<T> m_assembler;
    gsExprEvaluator<T> m_evaluator;

    gsMultiPatch<T>           m_patches;
    mutable gsMultiBasis<T>   m_basis;
    const gsFunctionSet<T> *  m_spaceBasis;
    const gsFunctionSet<T> * m_force;
    gsBoundaryConditions<T>   m_bcs;

    mutable gsMappedSpline<2,T> m_mspline;
    mutable gsMultiPatch<T> m_sol;

    mutable gsOptionList m_options;


    // gsFunctionExpr<T> m_ms;

  }; // class gsBiharmonicExprAssembler

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsBiharmonicExprAssembler
   */
  void pybind11_init_gsBiharmonicExprAssembler(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace gismo


//////////////////////////////////////////////////
//////////////////////////////////////////////////


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsBiharmonicExprAssembler.hpp)
#endif
