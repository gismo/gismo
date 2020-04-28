/** @file gsApproxBetaSAssembler.h

    @brief Provides assembler for the gluing data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

#include <gsG1Basis/gsVisitorApproxBetaS.h>

namespace gismo
{

template <class T, class bhVisitor = gsVisitorApproxBetaS<T> >
class gsApproxBetaSAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsApproxBetaSAssembler(gsMultiPatch<T> const & mp,
                           gsMultiBasis<T> mb,
                           gsG1OptionList optionList)
        : m_mp(mp), m_mb(mb), m_optionList(optionList)
    {

        index_t p_tilde = m_optionList.getInt("p_tilde");
        index_t r_tilde = m_optionList.getInt("r_tilde");

        gsKnotVector<T> kv(0,1,0,p_tilde+1,p_tilde-r_tilde); // first,last,interior,mult_ends,mult_interior
        gsBSplineBasis<T> bsp_gD(kv);

        const boundaryInterface & item = mp.interfaces()[0];
        index_t uv_L = item.second().index() < 3 ? 1 : 0;
        index_t L = item.second().patch;

        gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(L).component(uv_L));

        index_t degree = temp_basis_first.maxDegree();
        for (size_t i = degree+1; i < temp_basis_first.knots().size() - (degree+1); i += temp_basis_first.knots().multiplicityIndex(i))
            bsp_gD.insertKnot(temp_basis_first.knot(i),p_tilde-r_tilde);

        m_basis_betaS = bsp_gD; // beta_L and beta_R

        refresh();
        assemble();
        solve();
    }

    gsBSpline<> get_beta_0() { return beta_R; }
    gsBSpline<> get_beta_1() { return beta_L; }

    void refresh();

    void assemble();

    inline void apply(bhVisitor & visitor,
                      int patchIndex = 0,
                      boxSide side = boundary::none);

    void solve() {
        gsSparseSolver<real_t>::CGDiagonal solver;
        gsVector<> sol;

        // alpha^S
        solver.compute(m_system.matrix());
        sol = solver.solve(m_system.rhs());
        
        gsInfo << "sol: "<< sol << "\n";

        gsGeometry<>::uPtr tilde_temp;
        tilde_temp = m_basis_betaS.makeGeometry(sol.topRows(m_basis_betaS.size()));
        beta_L = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

        tilde_temp = m_basis_betaS.makeGeometry(sol.bottomRows(m_basis_betaS.size()));
        beta_R = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

        gsWriteParaview(beta_L,"beta_L",5000);
        gsWriteParaview(beta_R,"beta_R",5000);

        // BETA
        // first,last,interior,mult_ends,mult_interior,degree
        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(m_mp.patch(0).basis().component(1)); // 0 -> v, 1 -> u
        index_t m_p = basis_edge.maxDegree(); // Minimum degree at the interface // TODO if interface basis are not the same

        gsKnotVector<> kv(0, 1, basis_edge.numElements()-1, 2 * m_p  + 1, 2 * m_p - 1 );
        gsBSplineBasis<> bsp(kv);

        gsMatrix<> greville = bsp.anchors();
        gsMatrix<> uv1, uv0, ev1, ev0;

        const index_t d = 2;
        gsMatrix<> D0(d,d);

        gsGeometry<>::Ptr beta_temp;

        uv0.setZero(2,greville.cols());
        uv0.bottomRows(1) = greville;

        uv1.setZero(2,greville.cols());
        uv1.topRows(1) = greville;

        const gsGeometry<> & P0 = m_mp.patch(0); // iFace.first().patch = 1
        const gsGeometry<> & P1 = m_mp.patch(1); // iFace.second().patch = 0
        // ======================================

        // ======== Determine bar{beta} ========
        for(index_t i = 0; i < uv1.cols(); i++)
        {
            P0.jacobian_into(uv0.col(i),ev0);
            P1.jacobian_into(uv1.col(i),ev1);

            D0.col(1) = ev0.col(0); // (DuFL, *)
            D0.col(0) = ev1.col(1); // (*,DuFR)

            uv0(0,i) = D0.determinant();
        }

        beta_temp = bsp.interpolateData(uv0.topRows(1), uv0.bottomRows(1));
        gsBSpline<> beta = dynamic_cast<gsBSpline<> &> (*beta_temp);

        uv0.setZero(2,greville.cols());
        uv0.bottomRows(1) = greville;

        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        for (index_t i = 0; i < uv0.cols(); i++)
        {
            P0.jacobian_into(uv0.col(i), ev0);
            uv0(0, i) = 1 * ev0.determinant();

        }

        beta_temp = bsp.interpolateData(uv0.topRows(1), uv0.bottomRows(1));
        gsBSpline<> alpha0 = dynamic_cast<gsBSpline<> &> (*beta_temp);


        uv1.setZero(2,greville.cols());
        uv1.topRows(1) = greville;


        // ======== Determine bar{alpha^(L)} == Patch 0 ========
        for (index_t i = 0; i < uv1.cols(); i++)
        {
            P1.jacobian_into(uv1.col(i), ev1);
            uv1(1, i) = 1 * ev1.determinant();

        }

        beta_temp = bsp.interpolateData(uv1.bottomRows(1), uv1.topRows(1));
        gsBSpline<> alpha1 = dynamic_cast<gsBSpline<> &> (*beta_temp);

        index_t p_size = 8;
        gsMatrix<> points(1, p_size);
        points.setRandom();
        points = points.array().abs();

        gsVector<> vec;
        vec.setLinSpaced(p_size,0,1);
        points = vec.transpose();

        gsInfo << "Points: " << points << " \n";
        gsInfo << "Beta: " << beta.eval(points) << " \n";
        gsInfo << "alpha1: " << alpha1.eval(points) << " \n";
        gsInfo << "alpha0: " << alpha0.eval(points) << " \n";
        gsInfo << "beta_L: " << beta_L.eval(points) << " \n";
        gsInfo << "beta_R: " << beta_R.eval(points) << " \n";

        gsMatrix<> temp;
        temp = alpha1.eval(points).cwiseProduct(beta_L.eval(points))
            + alpha0.eval(points).cwiseProduct(beta_R.eval(points))
            - beta.eval(points);


        gsInfo << "Conditiontest Gluing data: \n" << temp.array().abs().maxCoeff() << "\n\n";

    }


protected:
    // interface + geometry for computing alpha and beta
    gsMultiPatch<T> m_mp;
    // space for beta, alpha
    gsMultiBasis<T>  m_mb;

    gsG1OptionList m_optionList;

    gsBSplineBasis<> m_basis_betaS;

    //using Base::m_system;
    using Base::m_ddof;

    gsSparseSystem<T> m_system;

    gsBSpline<> beta_L, beta_R;

}; // class gsG1BasisAssembler


template <class T, class bhVisitor>
void gsApproxBetaSAssembler<T, bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsVector<index_t> size(1);
    size.at(0) = 2* m_basis_betaS.size(); // 4 * beta_S

    gsDofMapper map(size);

    map.finalize();


    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(map);

} // refresh()

template <class T, class bhVisitor>
void gsApproxBetaSAssembler<T, bhVisitor>::assemble()
{
    //GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis_betaS,2,1,0.333333);
    m_system.reserve(nz, 1);


    if(m_ddof.size()==0)
        m_ddof.resize(1); // One for L, one for R, x2 for beta, alpha

    const gsDofMapper & mapper = m_system.colMapper(0);

    m_ddof[0].setZero(mapper.boundarySize(), 1 );

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor);

    m_system.matrix().makeCompressed();

}

template <class T, class bhVisitor>
inline void gsApproxBetaSAssembler<T, bhVisitor>::apply(bhVisitor & visitor,
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

        gsBasis<T> & basis = m_basis_betaS; // = 0

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
            visitor_.evaluate(m_basis_betaS, quNodes, m_mp);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex, m_ddof, m_system); // omp_locks inside

        }

    }//omp parallel
}


} // namespace gismo