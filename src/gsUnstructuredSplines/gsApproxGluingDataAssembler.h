/** @file gsGluingDataAssembler.h

    @brief Provides assembler for the gluing data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

#include <gsUnstructuredSplines/gsApproxGluingDataVisitor.h>

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
        interpolation_value = false;
        interpolation_deriv = false;

        if (m_bspGD.size() < 4) // Only if the gd space is too small
            interpolation_deriv = false;

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

    bool interpolation_value, interpolation_deriv;

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

    if (interpolation_value)
    {
        gsMatrix<index_t> act;
        act = m_bspGD.boundaryOffset(1,0); // WEST
        map_beta_S.markBoundary(0, act); // Patch 0
        map_alpha_S.markBoundary(0, act); // Patch 0

        act = m_bspGD.boundaryOffset(2,0); // East
        map_beta_S.markBoundary(0, act); // Patch 0
        map_alpha_S.markBoundary(0, act); // Patch 0

        if (interpolation_deriv)
        {
            act = m_bspGD.boundaryOffset(1,1); // WEST
            map_beta_S.markBoundary(0, act); // Patch 0
            map_alpha_S.markBoundary(0, act); // Patch 0

            act = m_bspGD.boundaryOffset(2,1); // East
            map_beta_S.markBoundary(0, act); // Patch 0
            map_alpha_S.markBoundary(0, act); // Patch 0
        }
    }


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

    m_ddof[0].setZero(mapper_alpha.boundarySize(), 1 ); // alpha
    m_ddof[1].setZero(mapper_beta.boundarySize(), 1 ); // beta

    // Compute Boundary of Gluing Data
    if (mapper_alpha.boundarySize() > 0)
    {
        // alpha^S
        gsMatrix<> points(1,2), uv, ev, D0;
        points.setZero();
        points(0,1) = 1.0;

        uv.setZero(2,points.cols());
        uv.row(m_uv) = points; // u

        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        const gsGeometry<> & P0 = m_patch.patch(0); // Right
        for (index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i), ev);
            uv(0, i) = 1 * ev.determinant();
        }
        //m_ddof[0].block(0,0,2,1) = uv.row(0).transpose();
        m_ddof[0](mapper_alpha.bindex(0),0) = uv(0,0);
        m_ddof[0](mapper_alpha.bindex(m_bspGD.size()-1),0) = uv(0,1);


        // Beta
        uv.setZero(2,points.cols());
        uv.row(m_uv) = points; // u
        for(index_t i = 0; i < uv.cols(); i++)
        {
            P0.jacobian_into(uv.col(i),ev);
            D0 = ev.col(m_uv);
            real_t D1 = 1/ D0.norm();
            uv(0,i) = - 1 * D1 * D1 * ev.col(1).transpose() * ev.col(0);
        }
        //m_ddof[1].block(0,0,2,1)  = uv.row(0).transpose();
        m_ddof[1](mapper_beta.bindex(0),0) = uv(0,0);
        m_ddof[1](mapper_beta.bindex(m_bspGD.size()-1),0) = uv(0,1);

        if (interpolation_deriv)
        {
            gsMatrix<> ev2;
            // alpha^S
            uv.setZero(2,points.cols());
            uv.row(m_uv) = points; // u

            // ======== Determine bar{alpha^(L)} == Patch 0 ========
            for (index_t i = 0; i < uv.cols(); i++)
            {
                P0.jacobian_into(uv.col(i), ev);
                P0.deriv2_into(uv.col(i), ev2);
                if (m_uv == 1)
                    uv(0, i) = 1.0 * (ev2(2,0)*ev(1,1) + ev2(4,0)*ev(0,0) -
                                                ev2(1,0)*ev(1,0) - ev2(5,0)*ev(0,1));
                else if (m_uv == 0)
                    uv(0, i) = 1.0 * (ev2(0,0)*ev(1,1) + ev2(5,0)*ev(0,0) -
                                                ev2(2,0)*ev(1,0) - ev2(3,0)*ev(0,1));
            }
            // Correctur
            //gsInfo << "Alpha deriv 2: " << uv.row(0) << "\n";
            uv(0,0) = (uv(0,0) - m_ddof[0](mapper_alpha.bindex(0),0) * m_bspGD.derivSingle(0,points)(0,0)) *
                    (1 / m_bspGD.derivSingle(1,points)(0,0));
            uv(0,1) = (uv(0,1) - m_ddof[0](mapper_alpha.bindex(m_bspGD.size()-1),0) * m_bspGD.derivSingle(m_bspGD.size()-1,points)(0,1)) *
                      (1 / m_bspGD.derivSingle(m_bspGD.size()-2,points)(0,1));
            //m_ddof[0].block(2,0,2,1)  = uv.row(0).transpose();
            //gsInfo << "Alpha deriv exact 2: " << uv(0,1)*m_bspGD.derivSingle(m_bspGD.size()-2,points)(0,1) + m_ddof[0](1,0) * m_bspGD.derivSingle(m_bspGD.size()-1,points)(0,1) << "\n";
            m_ddof[0](mapper_alpha.bindex(1),0) = uv(0,0);
            m_ddof[0](mapper_alpha.bindex(m_bspGD.size()-2),0) = uv(0,1);

            // beta^S
            uv.setZero(2,points.cols());
            uv.row(m_uv) = points; // u

            const index_t d = m_patch.patch(0).parDim();
            gsVector<> D0(d);

            // ======== Determine bar{beta}^L ========
            for(index_t i = 0; i < uv.cols(); i++)
            {
                P0.jacobian_into(uv.col(i),ev);
                P0.deriv2_into(uv.col(i), ev2);
                D0 = ev.col(m_uv);
                real_t D1 = 1/ D0.squaredNorm();
                real_t D2 = D0.squaredNorm();
                if (m_uv == 1)
                    uv(0,i) = - 1.0 * D1 * D1 * (D2*(ev2(2,0)*ev(0,1) + ev2(1,0)*ev(0,0)+
                                                               ev2(5,0)*ev(1,1) + ev2(4,0)*ev(1,0)) -
                                                           (ev.col(1).transpose() * ev.col(0))(0,0) * 2.0 * (ev2(1,0)*ev(0,1) + ev2(4,0)*ev(1,1)));
                else if (m_uv == 0)
                    uv(0,i) = - 1.0 * D1 * D1 * (D2*(ev2(0,0)*ev(0,1) + ev2(2,0)*ev(0,0)+
                                                               ev2(3,0)*ev(1,1) + ev2(5,0)*ev(1,0)) -
                                                           (ev.col(1).transpose() * ev.col(0))(0,0) * 2.0 * (ev2(0,0)*ev(0,0) + ev2(3,0)*ev(1,0)));

            }
            // Correctur
            uv(0,0) = (uv(0,0) - m_ddof[1](mapper_beta.bindex(0),0) * m_bspGD.derivSingle(0,points)(0,0)) *
                      (1 / m_bspGD.derivSingle(1,points)(0,0));
            uv(0,1) = (uv(0,1) - m_ddof[1](mapper_beta.bindex(m_bspGD.size()-1),0) * m_bspGD.derivSingle(m_bspGD.size()-1,points)(0,1)) *
                      (1 / m_bspGD.derivSingle(m_bspGD.size()-2,points)(0,1));
            //m_ddof[1].block(2,0,2,1)  = uv.row(0).transpose();
            m_ddof[1](mapper_beta.bindex(1),0) = uv(0,0);
            m_ddof[1](mapper_beta.bindex(m_bspGD.size()-2),0) = uv(0,1);
        }
    }

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
    if (m_system_alpha_S.matrix().dim().first > 0)
    {
        solver.compute(m_system_alpha_S.matrix());
        sol_a = solver.solve(m_system_alpha_S.rhs());
    }



    gsVector<> alpha_sol, beta_sol;
    alpha_sol.setZero(m_bspGD.size());
    beta_sol.setZero(m_bspGD.size());

    const gsDofMapper & mapper = m_system_alpha_S.colMapper(0); // unknown = 0
    for (index_t i = 0; i < m_bspGD.size(); ++i)
    {
        if (mapper.is_free(i, 0)) // DoF value is in the solVector // 0 = unitPatch
        {
            alpha_sol[i] = sol_a(mapper.index(i, 0),0);
        }
        else // eliminated DoF: fill with Dirichlet data
        {
            alpha_sol[i] = m_ddof[0]( mapper.bindex(i, 0), 0); // = 0
        }
    }


    gsGeometry<>::uPtr tilde_temp;
    tilde_temp = m_bspGD.makeGeometry(alpha_sol);
    alpha_S = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    // Solve beta^S
    if (m_system_beta_S.matrix().dim().first > 0)
    {
        solver.compute(m_system_beta_S.matrix());
        sol_b = solver.solve(m_system_beta_S.rhs());
    }



    const gsDofMapper & mapper_b = m_system_beta_S.colMapper(0); // unknown = 0
    for (index_t i = 0; i < m_bspGD.size(); ++i)
    {
        if (mapper_b.is_free(i, 0)) // DoF value is in the solVector // 0 = unitPatch
        {
            beta_sol[i] = sol_b(mapper_b.index(i, 0),0);
        }
        else // eliminated DoF: fill with Dirichlet data
        {
            beta_sol[i] = m_ddof[1]( mapper_b.bindex(i, 0), 0); // = 0
        }
    }

    tilde_temp = m_bspGD.makeGeometry(beta_sol);
    beta_S = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    gsMatrix<> points(1,2);
    points.setZero();
    points(0,1) = 1.0;
    /*
    gsInfo << "Alpha: " << alpha_S.eval(points) << "\n";
    gsInfo << "Beta: " << beta_S.eval(points) << "\n";
    gsInfo << "Alpha deriv: " << alpha_S.deriv(points) << "\n";
    gsInfo << "Beta deriv: " << beta_S.deriv(points) << "\n";
     */

} // solve

} // namespace gismo