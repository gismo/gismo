/** @file gsBiharmonicMethods.h

    @brief Compute the biharmonic equation with the Argyris basis functions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsArgyris/gsC1Argyris.h>

#include <gsAssembler/gsBiharmonicArgyrisAssembler.h>
#include <gsAssembler/gsBiharmonicNitscheAssembler.h>
#include <gsAssembler/gsBiharmonicArgyrisDirectAssembler.h>

#include <gsArgyris/gsErrorAnalysis/gsC1ArgyrisNorms.h>
#include <gsArgyris/gsErrorAnalysis/gsC1ArgyrisJumpNorm.h>
#include <gsArgyris/gsErrorAnalysis/gsC1NitscheNorms.h>
#include <gsArgyris/gsErrorAnalysis/gsC1NitscheJumpNorm.h>
#include <gsG1Basis/gsG1Norm.h>

namespace gismo
{
template<class T>
class gsBiharmonic
{

public:
    gsBiharmonic () {};

    //virtual ~gsBiharmonic();

    virtual void init() { GISMO_NO_IMPLEMENTATION };
    virtual void assemble(gsBoundaryConditions<T> const & bconditions,
                          gsBoundaryConditions<T> const & bconditions2,
                          const gsFunction<T>           & rhs) { GISMO_NO_IMPLEMENTATION };

    virtual void constructSolution(gsMatrix<T> & solVector) { GISMO_NO_IMPLEMENTATION };
    virtual void error(const gsFunctionWithDerivatives<T> &exactSolution) { GISMO_NO_IMPLEMENTATION };

    virtual index_t numDofs() const { GISMO_NO_IMPLEMENTATION };
    virtual const gsSparseMatrix<T> & matrix() const { GISMO_NO_IMPLEMENTATION };
    virtual gsMatrix<T> & rhs() { GISMO_NO_IMPLEMENTATION };

    virtual T valueL2() const { GISMO_NO_IMPLEMENTATION };
    virtual T valueH1() const { GISMO_NO_IMPLEMENTATION };
    virtual T valueH2() const { GISMO_NO_IMPLEMENTATION };
    virtual gsVector<T> valueJump() const { GISMO_NO_IMPLEMENTATION };
};

template<class T>
class gsBiharmonicArgyris : public gsBiharmonic<T>
{
public:

    /// Empty constructor
    gsBiharmonicArgyris() { }

    gsBiharmonicArgyris(const gsMultiPatch<T> & mp,
                        const gsMultiBasis<T> & mb,
                        const gsOptionList & optionList)
                        : m_mp(mp), m_mb(mb), m_optionList(optionList)
    {
        c1Argyris = gsC1Argyris<2, real_t>(m_mp, m_mb, m_optionList);
    }

    index_t numDofs() const { return g1BiharmonicAssembler->numDofs(); }
    const gsSparseMatrix<T> & matrix() const { return g1BiharmonicAssembler->matrix(); }
    gsMatrix<T> & rhs() { return g1BiharmonicAssembler->rhs(); }

    void init()
    {
        c1Argyris.init();
        c1Argyris.createArgyrisSpace();
        c1Argyris.getMultiBasis(mb_argyris);
        sparseMatrix_argyris = c1Argyris.getSystem();
        mappedBasis.init(mb_argyris, sparseMatrix_argyris.transpose());
    }

    void assemble(gsBoundaryConditions<T> const & bcInfo,
                  gsBoundaryConditions<T> const & bcInfo2,
                  const gsFunction<T>           & source)
    {
        g1BiharmonicAssembler = new gsBiharmonicArgyrisAssembler<real_t>(m_mp, mappedBasis, bcInfo, bcInfo2, source, m_optionList.getSwitch("twoPatch"));
        g1BiharmonicAssembler->assemble();
    }

    void constructSolution(gsMatrix<T> & solVector)
    {
        gsMatrix<real_t> solFull;
        g1BiharmonicAssembler->constructSolution(solVector, solFull);
        sparseMatrix_argyris = solFull.asDiagonal() * sparseMatrix_argyris;
        c1Argyris.setSystem(sparseMatrix_argyris);
    }

    void error(const gsFunctionWithDerivatives<T> &solution)
    {
        mappedBasis.init(mb_argyris, sparseMatrix_argyris.transpose());
        gsC1ArgyrisNorms<real_t> argyrisNorms(m_mp, mappedBasis, solution);
        argyrisNorms.compute();
        gsC1ArgyrisJumpNorm<real_t> c1ArgyrisJumpNorm(m_mp, mappedBasis, solution);
        c1ArgyrisJumpNorm.compute();

        l2Error = argyrisNorms.valueL2();
        h1Error = argyrisNorms.valueH1();
        h2Error = argyrisNorms.valueH2();
        jumpError = c1ArgyrisJumpNorm.value();
    }

    T valueL2() const { return l2Error; }
    T valueH1() const { return h1Error; }
    T valueH2() const { return h2Error; }
    gsVector<T> valueJump() const { return jumpError; }


protected:
    /// Multipatch
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;

    /// Optionlist
    gsOptionList m_optionList;

protected:
    gsC1Argyris<2, real_t> c1Argyris;
    gsMappedBasis<2,real_t> mappedBasis;

    gsSparseMatrix<> sparseMatrix_argyris;
    gsMultiBasis<> mb_argyris;
protected:
    gsBiharmonicArgyrisAssembler<real_t> *g1BiharmonicAssembler;

    gsVector<T> jumpError;
    T l2Error, h1Error, h2Error;

}; // class gsBiharmonicArgyris


template<class T>
class gsBiharmonicNitsche : public gsBiharmonic<T>
{
public:

    /// Empty constructor
    gsBiharmonicNitsche() { }

    gsBiharmonicNitsche(const gsMultiPatch<T> & mp,
                        const gsMultiBasis<T> & mb,
                        const gsOptionList & optionList)
            : m_mp(mp), m_mb(mb), m_optionList(optionList)
    {

    }

    index_t numDofs() const { return biharmonicNitscheAssembler->numDofs(); }
    const gsSparseMatrix<T> & matrix() const { return biharmonicNitscheAssembler->matrix(); }
    gsMatrix<T> & rhs() { return biharmonicNitscheAssembler->rhs(); }

    void init()
    {

    }

    void assemble(gsBoundaryConditions<T> const & bcInfo,
                  gsBoundaryConditions<T> const & bcInfo2,
                  const gsFunction<T>           & source)
    {
        biharmonicNitscheAssembler = new gsBiharmonicNitscheAssembler<real_t>(m_mp, m_mb, bcInfo, bcInfo2, source, m_optionList);
        biharmonicNitscheAssembler->assemble();
    }

    void constructSolution(gsMatrix<T> & solVector)
    {
        biharmonicNitscheAssembler->constructSolution(solVector, mpsol);
    }

    void error(const gsFunctionWithDerivatives<T> &solution)
    {
        gsC1NitscheNorms<real_t> c1NitscheNorms(m_mp, mpsol, solution);
        c1NitscheNorms.compute();
        gsC1NitscheJumpNorm<real_t> c1NitscheJumpNorm(m_mp, mpsol, solution);
        c1NitscheJumpNorm.compute();

        l2Error = c1NitscheNorms.valueL2();
        h1Error = c1NitscheNorms.valueH1();
        h2Error = c1NitscheNorms.valueH2();
        jumpError = c1NitscheJumpNorm.value();
    }

    T valueL2() const { return l2Error; }
    T valueH1() const { return h1Error; }
    T valueH2() const { return h2Error; }
    gsVector<T> valueJump() const { return jumpError; }


protected:
    /// Multipatch
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;

    /// Optionlist
    gsOptionList m_optionList;

protected:
    gsMultiPatch<> mpsol;

protected:
    gsBiharmonicNitscheAssembler<real_t> *biharmonicNitscheAssembler;

    gsVector<T> jumpError;
    T l2Error, h1Error, h2Error;

}; // class gsBiharmonicNitsche


    template<class T>
    class gsBiharmonicArgyrisDirect : public gsBiharmonic<T>
    {
    public:

        /// Empty constructor
        gsBiharmonicArgyrisDirect() { }

        gsBiharmonicArgyrisDirect(const gsMultiPatch<T> & mp,
                            const gsMultiBasis<T> & mb,
                            const gsOptionList & optionList)
                : m_mp(mp), m_mb(mb), m_optionList(optionList)
        {

        }

        index_t numDofs() const { return biharmonicArgyrisDirectAssembler->numDofs(); }
        const gsSparseMatrix<T> & matrix() const { return biharmonicArgyrisDirectAssembler->matrix(); }
        gsMatrix<T> & rhs() { return biharmonicArgyrisDirectAssembler->rhs(); }

        void init()
        {

        }

        void assemble(gsBoundaryConditions<T> const & bcInfo,
                      gsBoundaryConditions<T> const & bcInfo2,
                      const gsFunction<T>           & source)
        {
            biharmonicArgyrisDirectAssembler = new gsBiharmonicArgyrisDirectAssembler<real_t>(m_mp, m_mb, bcInfo,
                                                          bcInfo2, source);
            biharmonicArgyrisDirectAssembler->assemble();
        }

        void constructSolution(gsMatrix<T> & solVector)
        {
            biharmonicArgyrisDirectAssembler->constructSolution(solVector, mpsol);
            biharmonicArgyrisDirectAssembler->constructG1Solution(solVector, g1Sol);
        }

        void error(const gsFunctionWithDerivatives<T> &solution)
        {
            gsG1Norm<real_t> g1Norm(m_mp, m_mb, mpsol, g1Sol, solution);
            g1Norm.compute();

            l2Error = g1Norm.valueL2();
            h1Error = g1Norm.valueH1();
            h2Error = g1Norm.valueH2();
            //jumpError = c1ArgyrisJumpNorm.value();
            jumpError.setZero(1);
        }

        T valueL2() const { return l2Error; }
        T valueH1() const { return h1Error; }
        T valueH2() const { return h2Error; }
        gsVector<T> valueJump() const { return jumpError; }


    protected:
        /// Multipatch
        gsMultiPatch<T> m_mp;
        gsMultiBasis<T> m_mb;

        /// Optionlist
        gsOptionList m_optionList;

    protected:
        gsMultiPatch<> mpsol;
        gsMatrix<real_t> g1Sol;

    protected:
        gsBiharmonicArgyrisDirectAssembler<real_t> *biharmonicArgyrisDirectAssembler;

        gsVector<T> jumpError;
        T l2Error, h1Error, h2Error;

    }; // class gsBiharmonicArgyris

} // namespace gismo