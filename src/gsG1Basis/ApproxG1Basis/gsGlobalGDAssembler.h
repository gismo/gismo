/** @file gsGluingDataAssembler.h

    @brief Provides assembler for the gluing data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

#include <gsG1Basis/ApproxG1Basis/gsVisitorGlobalApproxGD.h>

namespace gismo
{

template <class T, class bhVisitor = gsVisitorGlobalApproxGD<T> >
class gsGlobalGDAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsGlobalGDAssembler(gsBasis<T> const & basis,
                           index_t const & uv,
                           index_t const & patchID,
                           gsMultiPatch<T> const & mp,
                           real_t const & gamma,
                           bool & isBoundary,
                           bool twoPatch)
        : m_uv(uv), m_patchID(patchID), m_mp(mp), m_gamma(gamma), m_isBoundary(isBoundary), m_twoPatch(twoPatch)
    {

        m_basis.push_back(basis); // Basis for alpha and beta

        refresh();
    }

    void refresh();

    void assemble();

    void computeBoundaryValues();

    inline void apply(bhVisitor & visitor,
                      int patchIndex = 0,
                      boxSide side = boundary::none);

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix_alpha() const { return m_system_alpha.matrix(); }
    const gsMatrix<T> & rhs_alpha() const { return m_system_alpha.rhs(); }
    const gsMatrix<T> & bdy_alpha() const { return m_ddof[0]; }

    const gsSparseMatrix<T> & matrix_beta() const { return m_system_beta_S.matrix(); }
    const gsMatrix<T> & rhs_beta() const { return m_system_beta_S.rhs(); }

protected:
    index_t m_uv, m_patchID;

    // interface + geometry for computing alpha and beta
    gsMultiPatch<T> m_mp;

    real_t m_gamma;
    bool m_isBoundary;
    bool m_twoPatch;

    // Space for phi_0,i, phi_1,j
    std::vector< gsMultiBasis<T> > m_basis;

    //using Base::m_system;
    using Base::m_ddof;

    gsSparseSystem<T> m_system_alpha;
    gsSparseSystem<T> m_system_beta_S;

}; // class gsG1BasisAssembler


template <class T, class bhVisitor>
void gsGlobalGDAssembler<T, bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map_alpha(m_basis[0]);

    gsDofMapper map_beta_S(m_basis[0]);

    gsMatrix<unsigned> act;
    act = m_basis[0].basis(0).boundaryOffset(1,0); // WEST
    if (m_twoPatch)
    {
        map_beta_S.markBoundary(0, act); // Patch 0
        map_alpha.markBoundary(0, act); // Patch 0
    }

    act = m_basis[0].basis(0).boundaryOffset(2,0); // East
    if (m_twoPatch)
    {
        map_beta_S.markBoundary(0, act); // Patch 0
        map_alpha.markBoundary(0, act); // Patch 0
    }
    map_alpha.finalize();
    map_beta_S.finalize();

    // 2. Create the sparse system
    m_system_alpha = gsSparseSystem<T>(map_alpha);
    m_system_beta_S = gsSparseSystem<T>(map_beta_S);

} // refresh()

template <class T, class bhVisitor>
void gsGlobalGDAssembler<T, bhVisitor>::assemble()
{
    //GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis[0][0],2,1,0.333333);
    m_system_alpha.reserve(nz, 1);
    m_system_beta_S.reserve(nz, 1);


    if(m_ddof.size()==0)
        m_ddof.resize(2); // One for L, one for R, x2 for beta, alpha

    const gsDofMapper & mapper_alpha = m_system_alpha.colMapper(0);
    const gsDofMapper & mapper_beta = m_system_beta_S.colMapper(0);

    m_ddof[0].setZero(mapper_alpha.boundarySize(), 1 ); // tilde
    m_ddof[1].setZero(mapper_beta.boundarySize(), 1 ); // bar

    // Boundary
    //computeBoundaryValues(); == 0
    // alpha^S
    gsMatrix<> points(1,2), uv, ev;
    points.setZero();
    points(0,1) = 1.0;

    uv.setZero(2,points.cols());
    uv.row(m_uv) = points; // u

    // ======== Determine bar{alpha^(L)} == Patch 0 ========
    const gsGeometry<> & P0 = m_mp.patch(m_patchID); // Right
    for (index_t i = 0; i < uv.cols(); i++)
    {
        P0.jacobian_into(uv.col(i), ev);
        uv(0, i) = 1 * ev.determinant();

    }
    if (m_isBoundary)
        uv.setOnes();

    if (m_twoPatch)
        m_ddof[0] = uv.row(0).transpose();

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor);

    m_system_alpha.matrix().makeCompressed();
    m_system_beta_S.matrix().makeCompressed();

}

template <class T, class bhVisitor>
inline void gsGlobalGDAssembler<T, bhVisitor>::apply(bhVisitor & visitor,
                                                        int patchIndex,
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

        gsBasis<T> & basis = m_basis[0].basis(0); // = 0

        // Initialize reference quadrature rule and visitor data
        visitor_.initialize(basis,quRule);

        //const gsGeometry<T> & patch = m_geo.patch(patchIndex); // 0 = patchindex

        // Initialize domain element iterator -- using unknown 0
        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator(boundary::none);

#ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
        for (; domIt->good(); domIt->next() )
#endif
        {

            // Map the Quadrature rule to the element
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

            // Perform required evaluations on the quadrature nodes
            visitor_.evaluate(basis, quNodes, m_uv, m_mp, m_patchID, m_gamma, m_isBoundary, m_twoPatch);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex, m_ddof, m_system_alpha, m_system_beta_S); // omp_locks inside

        }

    }//omp parallel
}

template <class T, class bhVisitor>
void gsGlobalGDAssembler<T, bhVisitor>::computeBoundaryValues()
{
    real_t D1, lambda0, lambda1;
    gsMatrix<> uv, ev, D0;

    gsMatrix<> alpha_S, beta_S;

    gsMatrix<> zeroOne(2,2);
    zeroOne.setZero();
    zeroOne(1,1) = 1.0; // v

    gsMatrix<> points(1,2);
    points.setZero();
    points(0,1) = 1.0;

    const gsGeometry<> & PR = m_mp.patch(0); // Right
    const gsGeometry<> & PL = m_mp.patch(1); // Left

    PR.jacobian_into(zeroOne.col(1), ev);
    lambda1 = 1/ev.determinant(); // alpha_R

    D0 = ev.col(1);
    D1 = 1/ D0.norm();
    lambda1 *= - m_gamma * D1 * D1 * ev.col(1).transpose() * ev.col(0);



    PL.jacobian_into(zeroOne.col(0), ev);
    lambda0 = 1/ev.determinant(); // alpha_L

    D0 = ev.col(0);
    D1 = 1/ D0.norm();
    lambda0 *= - m_gamma * D1 * D1 * ev.col(1).transpose() * ev.col(0);

    if (m_uv == 1)
        lambda0 = - lambda0;

    if (m_uv == 0)
        lambda1 = - lambda1;

    // alpha^S
    if (m_uv==1)
    {
        uv.setZero(2,points.cols());
        uv.bottomRows(1) = points; // v
    }
    else if (m_uv==0)
    {
        uv.setZero(2,points.cols());
        uv.topRows(1) = points; // u
    }


    // ======== Determine bar{alpha^(L)} == Patch 0 ========
    const gsGeometry<> & P0 = m_mp.patch(m_patchID); // Right
    for (index_t i = 0; i < uv.cols(); i++)
    {
        P0.jacobian_into(uv.col(i), ev);
        uv(0, i) = 1 * ev.determinant();

    }
    if (m_isBoundary)
        uv.setOnes();
    alpha_S = uv.row(0);

    // beta^S
    if (m_uv==1)
    {
        uv.setZero(2,points.cols());
        uv.bottomRows(1) = points; // v
    }
    else if (m_uv==0)
    {
        uv.setZero(2,points.cols());
        uv.topRows(1) = points; // u
    }

    // ======== Determine bar{beta}^L ========
    for(index_t i = 0; i < uv.cols(); i++)
    {
        P0.jacobian_into(uv.col(i),ev);
        D0 = ev.col(m_uv);
        real_t D1 = 1/ D0.norm();
        uv(0,i) = - m_gamma * D1 * D1 * ev.col(1).transpose() * ev.col(0);
        //uv(0,i) = - ev.col(1).transpose() * ev.col(0);

    }
    if (m_isBoundary)
        uv.setZero();
    beta_S = uv.row(0);


    // ++++++++++++++++++++++++++++++++
    // ================================
    // ++++++++++++++++++++++++++++++++
    gsMatrix<> ones;
    ones.setOnes(1,points.cols());

    m_ddof[1] = (beta_S - lambda0 * (ones - points).cwiseProduct(alpha_S) - lambda1 * (points).cwiseProduct(alpha_S)).transpose();

}


} // namespace gismo