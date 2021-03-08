/** @file gsGluingDataAssembler.h

    @brief Provides assembler for the gluing data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

#include <gsArgyris/gsGluingData/gsApproxGluingDataVisitor.h>

namespace gismo
{

template <class T, class bhVisitor = gsApproxGluingDataVisitor<T> >
class gsApproxGluingDataAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsApproxGluingDataAssembler(const gsGeometry<> & patch,
                                gsBSplineBasis<T> & bsp_Gd,
                                index_t uv,
                                const gsOptionList & optionList)
        : m_patch(patch), m_bspGD(bsp_Gd), m_uv(uv), m_optionList(optionList)
    {
        refresh();
        assemble();
        solve();
    }

    void refresh();

    void assemble();

    void solve();

    inline void apply(bhVisitor & visitor,
                      index_t patchIndex = 0,
                      boxSide side = boundary::none);

    gsBSpline<T> getAlphaS() { return alpha_S; }
    gsBSpline<T> getBetaS() { return beta_S; }

protected:
    // interface + geometry for computing alpha and beta
    gsMultiPatch<T> m_patch;
    gsBSplineBasis<T> m_bspGD;
    index_t m_uv;

    gsOptionList m_optionList;

    //using Base::m_system;
    using Base::m_ddof;

    gsSparseSystem<T> m_system_alpha_S;
    gsSparseSystem<T> m_system_beta_S;

    gsBSpline<T> alpha_S;
    gsBSpline<T> beta_S;

}; // class gsG1BasisAssembler


template <class T, class bhVisitor>
void gsApproxGluingDataAssembler<T, bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map_alpha_S(m_bspGD);

    gsDofMapper map_beta_S(m_bspGD);

    map_alpha_S.finalize();
    map_beta_S.finalize();

    // 2. Create the sparse system
    m_system_alpha_S = gsSparseSystem<T>(map_alpha_S);
    m_system_beta_S = gsSparseSystem<T>(map_beta_S);

} // refresh()

template <class T, class bhVisitor>
void gsApproxGluingDataAssembler<T, bhVisitor>::assemble()
{
    //GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_bspGD,2,1,0.333333);
    m_system_alpha_S.reserve(nz, 1);
    m_system_beta_S.reserve(nz, 1);


    if(m_ddof.size()==0)
        m_ddof.resize(2); // One for L, one for R, x2 for beta, alpha

    const gsDofMapper & mapper_alpha = m_system_alpha_S.colMapper(0);
    const gsDofMapper & mapper_beta = m_system_beta_S.colMapper(0);

    m_ddof[0].setZero(mapper_alpha.boundarySize(), 1 ); // tilde
    m_ddof[1].setZero(mapper_beta.boundarySize(), 1 ); // bar

    // Compute Boundary of Gluing Data

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor);

    m_system_alpha_S.matrix().makeCompressed();
    m_system_beta_S.matrix().makeCompressed();

}

template <class T, class bhVisitor>
inline void gsApproxGluingDataAssembler<T, bhVisitor>::apply(bhVisitor & visitor,
                                                        index_t patchIndex,
                                                        boxSide side)
{
#pragma omp parallel
    {
        bhVisitor
#ifdef _OPENMP
        // Create thread-private visitor
        visitor_(visitor);
        const int tid = omp_get_thread_num();
        const int nt  = omp_get_num_threads();
#else
            &visitor_ = visitor;
#endif

        gsQuadRule<T> quRule ; // Quadrature rule
        gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
        gsVector<T> quWeights; // Temp variable for mapped weights

        // Initialize reference quadrature rule and visitor data
        visitor_.initialize(m_bspGD,quRule);

        //const gsGeometry<T> & patch = m_geo.patch(patchIndex); // 0 = patchindex

        // Initialize domain element iterator -- using unknown 0
        typename gsBasis<T>::domainIter domIt = m_bspGD.makeDomainIterator(boundary::none);

#ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
        for (; domIt->good(); domIt->next() )
#endif
        {

            // Map the Quadrature rule to the element
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

            // Perform required evaluations on the quadrature nodes
            visitor_.evaluate(m_patch, m_bspGD, quNodes, m_uv, m_optionList);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex, m_ddof, m_system_alpha_S, m_system_beta_S); // omp_locks inside

        }

    }//omp parallel
}

template <class T, class bhVisitor>
void gsApproxGluingDataAssembler<T, bhVisitor>::solve()
{
    gsSparseSolver<real_t>::LU solver;
    gsVector<> sol_a, sol_b;

    // Solve alpha^S
    solver.compute(m_system_alpha_S.matrix());
    sol_a = solver.solve(m_system_alpha_S.rhs());

    gsGeometry<>::uPtr tilde_temp;
    tilde_temp = m_bspGD.makeGeometry(sol_a);
    alpha_S = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    // Solve beta^S
    solver.compute(m_system_beta_S.matrix());
    sol_b = solver.solve(m_system_beta_S.rhs());

    tilde_temp = m_bspGD.makeGeometry(sol_b);
    beta_S = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

} // solve

} // namespace gismo