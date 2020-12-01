/** @file gsG1BasisVertex.h

    @brief Provides assembler for a G1 Basis for multiPatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

#include <gsG1Basis/ApproxG1Basis/gsApproxGluingData.h>

#include <gsG1Basis/ApproxG1Basis/gsVisitorApproxG1BasisVertex.h>
#include <gsG1Basis/gsG1OptionList.h>

namespace gismo
{
template<class T, class bhVisitor = gsVisitorG1BasisVertex<T>>
class gsApproxG1BasisVertex : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsApproxG1BasisVertex(gsMultiPatch<> mp, // Single Patch
                          gsMultiBasis<> basis, // Single Basis
                          std::vector<bool> isBoundary,
                          real_t sigma,
                          gsG1OptionList & g1OptionList)
        : m_mp(mp), m_basis(basis), m_isBoundary(isBoundary), m_sigma(sigma), m_g1OptionList(g1OptionList)
    {

        for (index_t dir = 0; dir < m_mp.parDim(); dir++) // For the TWO directions
        {
            // Computing the G1 - basis function at the edge
            // Spaces for computing the g1 basis
            gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(dir)); // 0 -> u, 1 -> v
            index_t m_p = basis_edge.maxDegree(); // Minimum degree at the interface // TODO if interface basis are not the same

            index_t m_r = g1OptionList.getInt("regularity") > m_p - 2 ? (m_p-2) : g1OptionList.getInt("regularity"); // TODO CHANGE IF DIFFERENT REGULARITY IS NECESSARY


            // first,last,interior,mult_ends,mult_interior
            gsKnotVector<T> kv_plus(0,1,0,m_p+1,m_p-1-m_r); // p,r+1 //-1 bc r+1
            gsBSplineBasis<> basis_plus(kv_plus);

            for (size_t i = m_p+1; i < basis_edge.knots().size() - (m_p+1); i = i+(m_p-m_r))
                basis_plus.insertKnot(basis_edge.knot(i),m_p-1-m_r);

            m_basis_plus.push_back(basis_plus);

            gsKnotVector<T> kv_minus(0,1,0,m_p+1-1,m_p-1-m_r); // p-1,r //-1 bc p-1
            gsBSplineBasis<> basis_minus(kv_minus);

            for (size_t i = m_p+1; i < basis_edge.knots().size() - (m_p+1); i = i+(m_p-m_r))
                basis_minus.insertKnot(basis_edge.knot(i),m_p-1-m_r);

            m_basis_minus.push_back(basis_minus);

            // Computing the gluing data
            gsApproxGluingData<T> gluingData(m_mp, m_basis, dir, m_isBoundary[dir], m_g1OptionList);
            if (g1OptionList.getInt("gluingData") == gluingData::local)
                gluingData.setLocalGluingData(basis_plus, basis_minus, "vertex");
            else if (g1OptionList.getInt("gluingData") == gluingData::global)
                gluingData.setGlobalGluingData(0,dir);

            m_gD.push_back(gluingData);
        }

        // Basis for the G1 basis
        m_basis_g1 = m_basis.basis(0);
/*
        if ( m_isBoundary[0] == 1 && m_isBoundary[1] == 0 )
        {
            m_basis_g1.basis(0).degreeElevate(3);
            gsInfo << "hab hier \n";
        }
        else if ( m_isBoundary[0] == 0 && m_isBoundary[1] == 1 )
        {
            m_basis_g1.basis(0).degreeElevate(3);
            gsInfo << "hab hier 2\n";
        }
*/
    }


    void refresh(index_t kindOfVertex);
    void assemble();
    inline void apply(bhVisitor & visitor, int patchIndex);
    void solve();
    void computeDofsIntpl();

    void constructSolution(gsMultiPatch<T> & result);

    void setG1BasisVertex(gsMultiPatch<T> & result, index_t kindOfVertex)
    {

        gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mp.basis(0).component(0)); // u
        index_t degree = temp_basis_first.maxDegree();

        gsMatrix<T> ab = m_basis_plus[0].support(2);
/*        if (kindOfVertex == 0)
        {
            gsMatrix<T> ab_temp = ab;
            for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
            {
                gsMatrix<T> xy = temp_basis_first.support(i);
                if ((xy(0, 0) < ab(0, 0)) && (xy(0, 1) > ab(0, 0)))
                    ab_temp(0, 0) = xy(0, 0);
                if ((xy(0, 0) < ab(0, 1)) && (xy(0, 1) > ab(0, 1)))
                    ab_temp(0, 1) = xy(0, 1);
            }
            ab = ab_temp;
            ab_temp = ab;
            for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
            {
                gsMatrix<T> xy = temp_basis_first.support(i);
                if ((xy(0, 0) < ab(0, 0)) && (xy(0, 1) > ab(0, 0)))
                    ab_temp(0, 0) = xy(0, 0);
                if ((xy(0, 0) < ab(0, 1)) && (xy(0, 1) > ab(0, 1)))
                    ab_temp(0, 1) = xy(0, 1);
            }
            ab = ab_temp;
        }
*/
        gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, 1);
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i += temp_basis_first.knots().multiplicityIndex(i))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv.insert(temp_basis_first.knot(i), 1);

        // #### v ####

        temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mp.basis(0).component(1)); // v
        degree = temp_basis_first.maxDegree();

        ab = m_basis_plus[1].support(2);
/*        if (kindOfVertex == 0)
        {
            gsMatrix<T> ab_temp = ab;
            for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
            {
                gsMatrix<T> xy = temp_basis_first.support(i);
                if ((xy(0, 0) < ab(0, 0)) && (xy(0, 1) > ab(0, 0)))
                    ab_temp(0, 0) = xy(0, 0);
                if ((xy(0, 0) < ab(0, 1)) && (xy(0, 1) > ab(0, 1)))
                    ab_temp(0, 1) = xy(0, 1);
            }
            ab = ab_temp;
            for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
            {
                gsMatrix<T> xy = temp_basis_first.support(i);
                if ((xy(0, 0) < ab(0, 0)) && (xy(0, 1) > ab(0, 0)))
                    ab_temp(0, 0) = xy(0, 0);
                if ((xy(0, 0) < ab(0, 1)) && (xy(0, 1) > ab(0, 1)))
                    ab_temp(0, 1) = xy(0, 1);
            }
            ab = ab_temp;
        }
*/
        gsKnotVector<T> kv2(ab.at(0), ab.at(1), 0, 1);
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i += temp_basis_first.knots().multiplicityIndex(i))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv2.insert(temp_basis_first.knot(i), 1);

        gsTensorBSplineBasis<2, T> bsp_geo_local(kv, kv2);
        if (m_g1OptionList.getInt("g1BasisVertex") == g1BasisVertex::local && kindOfVertex == -1)
            m_geo = bsp_geo_local; // Basis for Integration
        else
            m_geo = m_basis_g1;

        refresh(kindOfVertex);
        assemble();
        solve();

        constructSolution(result);
    }


protected:

    // Input
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_basis;
    std::vector<bool> m_isBoundary;
    real_t m_sigma;
    gsG1OptionList m_g1OptionList;

    // Gluing data
    std::vector<gsApproxGluingData<T>> m_gD;

    // Basis plus minus
    std::vector<gsBSplineBasis<>> m_basis_plus;
    std::vector<gsBSplineBasis<>> m_basis_minus;

    // Basis for the G1 Basis
    gsMultiBasis<T> m_basis_g1;

    // Basis for Integration
    gsMultiBasis<T> m_geo;

    // System
    std::vector<gsSparseSystem<T> > m_f;

    // For Dirichlet boundary
    using Base::m_ddof;

    std::vector<gsMatrix<>> solVec;
}; // class gsG1BasisEdge


template <class T, class bhVisitor>
void gsApproxG1BasisVertex<T, bhVisitor>::constructSolution(gsMultiPatch<T> & result)
{

    result.clear();

    // Dim is the same for all basis functions
    const index_t dim = ( 0!=solVec.at(0).cols() ? solVec.at(0).cols() :  m_ddof[0].cols() );

    gsMatrix<T> coeffs;
    for (index_t p = 0; p < 6; ++p)
    {

        const gsDofMapper & mapper = m_f.at(p).colMapper(0); // unknown = 0

        // Reconstruct solution coefficients on patch p
        index_t sz;
        sz = m_basis_g1.basis(0).size();

        coeffs.resize(sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if (mapper.is_free(i, 0)) // DoF value is in the solVector // 0 = unitPatch
            {
                coeffs.row(i) = solVec.at(p).row(mapper.index(i, 0));
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                coeffs.row(i) = m_ddof[p].row( mapper.bindex(i, 0) ).head(dim); // = 0
            }
        }
        result.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));
    }
}

template <class T, class bhVisitor>
void gsApproxG1BasisVertex<T, bhVisitor>::refresh(index_t kindOfVertex)
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map(m_basis_g1.basis(0));
    gsMatrix<unsigned> act;
    if ((m_g1OptionList.getInt("g1BasisVertex") == g1BasisVertex::local) && kindOfVertex == -1)
    {
        for (index_t dir = 0; dir < 2; dir++)
        {
            gsMatrix<T> ab = m_basis_plus[dir].support(2);
/*            if (kindOfVertex == 0)
            {
                gsMatrix<T> ab_temp = ab;
                for (index_t i = 0; i < m_basis.basis(0).component(dir).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    gsMatrix<T> xy = m_basis.basis(0).component(dir).support(i);
                    if ((xy(0, 0) < ab(0, 0)) && (xy(0, 1) > ab(0, 0)))
                        ab_temp(0, 0) = xy(0, 0);
                    if ((xy(0, 0) < ab(0, 1)) && (xy(0, 1) > ab(0, 1)))
                        ab_temp(0, 1) = xy(0, 1);
                }
                ab = ab_temp;

                for (index_t i = 0; i < m_basis.basis(0).component(dir).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    gsMatrix<T> xy = m_basis.basis(0).component(dir).support(i);
                    if ((xy(0, 0) < ab(0, 0)) && (xy(0, 1) > ab(0, 0)))
                        ab_temp(0, 0) = xy(0, 0);
                    if ((xy(0, 0) < ab(0, 1)) && (xy(0, 1) > ab(0, 1)))
                        ab_temp(0, 1) = xy(0, 1);
                }
                ab = ab_temp;
            }
*/
            for (index_t i = 0; i < m_basis_g1.basis(0).component(dir).size();
                 i++) // only the first two u/v-columns are Dofs (0/1)
            {
                gsMatrix<T> xy = m_basis_g1.basis(0).component(dir).support(i);
                if ((xy(0, 0) < ab(0, 0)) || (xy(0, 1) > ab(0, 1)))
                {
                    act = m_basis_g1.basis(0).boundaryOffset(dir == 0 ? 1 : 3, i); // WEST
                    map.markBoundary(0, act); // Patch 0
                }
            }
        }
    }

    /*
 * ...
 * o o o
 * x x o
 * x x o ...
 */
    act.resize(4, 1);
    act(0, 0) = (unsigned) 0;
    act(1, 0) = (unsigned) 1;
    act(2,0) = (unsigned) m_basis_g1.basis(0).component(0).size();
    act(3,0) = (unsigned) m_basis_g1.basis(0).component(0).size()+1;
    //map.markBoundary(0, act);

    map.finalize();
    //gsInfo << "map : " << map.asVector() << "\n";

    // 2. Create the sparse system
    gsSparseSystem<T> m_system = gsSparseSystem<T>(map);
    for (index_t i = 0; i < 6; i++)
        m_f.push_back(m_system);

} // refresh()

template <class T, class bhVisitor>
void gsApproxG1BasisVertex<T, bhVisitor>::assemble()
{
    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis_g1[0],2,1,0.333333);
    for (unsigned i = 0; i < m_f.size(); i++)
        m_f.at(i).reserve(nz, 1);


    const gsDofMapper & map = m_f.at(0).colMapper(0); // Map same for every

    m_ddof.resize(6); // 0,1
    for (size_t i = 0; i < m_ddof.size(); i++)
        m_ddof[i].setZero(map.boundarySize(), 1 ); // plus

    //computeDofsIntpl();

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor,0); // patch 0

    for (unsigned i = 0; i < m_f.size(); i++)
        m_f.at(i).matrix().makeCompressed();

} // assemble()

template <class T, class bhVisitor>
void gsApproxG1BasisVertex<T, bhVisitor>::apply(bhVisitor & visitor, int patchIndex)
{
#pragma omp parallel
    {

        gsQuadRule<T> quRule ; // Quadrature rule
        gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
        gsVector<T> quWeights; // Temp variable for mapped weights

        bhVisitor
#ifdef _OPENMP
        // Create thread-private visitor
    visitor_(visitor);
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#else
            &visitor_ = visitor;
#endif

        gsBasis<T> & basis_g1 = m_basis_g1.basis(0); // basis for construction

        // Same for all patches
        gsBasis<T> & basis_geo = m_basis.basis(0); // 2D

        // Initialize reference quadrature rule and visitor data
        visitor_.initialize(basis_g1, quRule);

        const gsGeometry<T> & patch = m_mp.patch(0);

        // Initialize domain element iterator
        typename gsBasis<T>::domainIter domIt = m_geo.basis(0).makeDomainIterator(boundary::none);

#ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
        for (; domIt->good(); domIt->next() )
#endif
        {
            // Map the Quadrature rule to the element
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );
#pragma omp critical(evaluate)
            // Perform required evaluations on the quadrature nodes
            visitor_.evaluate(basis_g1, basis_geo, m_basis_plus, m_basis_minus, patch, quNodes, m_gD, m_isBoundary, m_sigma, m_g1OptionList);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex, m_ddof, m_f); // omp_locks inside // patchIndex == 0
        }
    }//omp parallel
} // apply

template <class T, class bhVisitor>
void gsApproxG1BasisVertex<T, bhVisitor>::solve()
{
    gsSparseSolver<real_t>::LU solver;


    for (index_t i = 0; i < 6; i++) // Tilde
    {
        solver.compute(m_f.at(i).matrix());
        solVec.push_back(solver.solve(m_f.at(i).rhs()));
    }
} // solve

template <class T, class bhVisitor>
void gsApproxG1BasisVertex<T,bhVisitor>::computeDofsIntpl()

{
    const gsGeometry<T>    & geo = m_mp.patch(0); // patch
    gsBasis<T> & basis_geo = m_basis.basis(0); // 2D

    gsMatrix<> points(2,1);
    points.setZero();

    // Computing the basis functions at the vertex
    gsMatrix<> Phi(6,6);
    Phi.setIdentity();

    Phi.row(1) *= m_sigma;
    Phi.row(2) *= m_sigma;
    Phi.row(3) *= m_sigma * m_sigma;
    Phi.row(4) *= m_sigma * m_sigma;
    Phi.row(5) *= m_sigma * m_sigma;

    // Computing c, c+ and c-
    std::vector<gsMatrix<>> c_0, c_1;
    std::vector<gsMatrix < >> c_0_plus, c_1_plus, c_2_plus;
    std::vector<gsMatrix < >> c_0_plus_deriv, c_1_plus_deriv, c_2_plus_deriv;
    std::vector<gsMatrix < >> c_0_minus, c_1_minus;
    for (index_t i = 0; i < 2; i++) // i == 0 == u , i == 1 == v
    {
        gsMatrix<> b_0, b_1;
        gsMatrix<> b_0_plus, b_1_plus, b_2_plus;
        gsMatrix<> b_0_plus_deriv, b_1_plus_deriv, b_2_plus_deriv;
        gsMatrix<> b_0_minus, b_1_minus;

        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(basis_geo.component(i));
        real_t p = bsp_temp.maxDegree();
        real_t h_geo = bsp_temp.knots().at(p + 2);

        basis_geo.component(i).evalSingle_into(0, points.row(i),b_0); // first
        basis_geo.component(i).evalSingle_into(1, points.row(i),b_1); // second

        m_basis_plus[i].evalSingle_into(0, points.row(i),b_0_plus);
        m_basis_plus[i].evalSingle_into(1, points.row(i),b_1_plus);
        m_basis_plus[i].evalSingle_into(2, points.row(i),b_2_plus);

        m_basis_plus[i].derivSingle_into(0, points.row(i),b_0_plus_deriv);
        m_basis_plus[i].derivSingle_into(1, points.row(i),b_1_plus_deriv);
        m_basis_plus[i].derivSingle_into(2, points.row(i),b_2_plus_deriv);

        m_basis_minus[i].evalSingle_into(0, points.row(i),b_0_minus);
        m_basis_minus[i].evalSingle_into(1, points.row(i),b_1_minus);

        c_0.push_back(b_0 + b_1);
        c_1.push_back((h_geo / p) * b_1);

        c_0_minus.push_back(b_0_minus + b_1_minus);
        c_1_minus.push_back(h_geo/ (p-1) * b_1_minus);

        // TODO IF CASE
        if ( p == 3)
        {
            // WORKS ONLY FOR p=3 AND r=1
            c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
            c_1_plus.push_back((h_geo / p) * (b_1_plus + 3 * b_2_plus));
            c_2_plus.push_back((h_geo * h_geo / (p * (p - 1))) * 2 * b_2_plus);

            c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
            c_1_plus_deriv.push_back((h_geo / p) * (b_1_plus_deriv + 3 * b_2_plus_deriv));
            c_2_plus_deriv.push_back((h_geo * h_geo / (p * (p - 1))) * 2 * b_2_plus_deriv);
        }
        else
        {
            c_0_plus.push_back(b_0_plus + b_1_plus + b_2_plus);
            c_1_plus.push_back((h_geo / p) * (b_1_plus + 2 * b_2_plus));
            c_2_plus.push_back((h_geo * h_geo / (p * (p - 1))) * b_2_plus);

            c_0_plus_deriv.push_back(b_0_plus_deriv + b_1_plus_deriv + b_2_plus_deriv);
            c_1_plus_deriv.push_back((h_geo / p) * (b_1_plus_deriv + 2 * b_2_plus_deriv));
            c_2_plus_deriv.push_back((h_geo * h_geo / (p * (p - 1))) * b_2_plus_deriv);

        }
    }

    // Point zero
    gsMatrix<> zero;
    zero.setZero(2,1);

    std::vector<gsMatrix<>> alpha, beta, alpha_0, beta_0, alpha_deriv, beta_deriv;
/*
        if (g1OptionList.getInt("gluingData") == gluingData::global)
        {
            //gsInfo << "Evtl kann beta falsch sein!! \n";
            alpha.push_back(gluingData[0].get_alpha_S_tilde(0).eval(md.points.row(0))); // u
            alpha.push_back(gluingData[1].get_alpha_S_tilde(0).eval(md.points.row(1))); // v
            alpha_0.push_back(gluingData[0].get_alpha_S_tilde(0).eval(zero.row(0))); // u
            alpha_0.push_back(gluingData[1].get_alpha_S_tilde(0).eval(zero.row(0))); // v
            alpha_deriv.push_back(gluingData[0].get_alpha_S_tilde(0).deriv(zero.row(0))); // u
            alpha_deriv.push_back(gluingData[1].get_alpha_S_tilde(0).deriv(zero.row(0))); // v

            beta.push_back(gluingData[0].get_beta_S_tilde(0).eval(md.points.row(0))); // u
            beta.push_back(gluingData[1].get_beta_S_tilde(0).eval(md.points.row(1))); // v
            beta_0.push_back(gluingData[0].get_beta_S_tilde(0).eval(zero.row(0))); // u
            beta_0.push_back(gluingData[1].get_beta_S_tilde(0).eval(zero.row(0))); // v
            beta_deriv.push_back(gluingData[0].get_beta_S_tilde(0).deriv(zero.row(0))); // u
            beta_deriv.push_back(gluingData[1].get_beta_S_tilde(0).deriv(zero.row(0))); // v

        }
        else if (g1OptionList.getInt("gluingData") == gluingData::local)
        {
            alpha.push_back(gluingData[0].get_alpha_S_tilde(0).eval(md.points.row(0))); // u
            alpha.push_back(gluingData[1].get_alpha_S_tilde(0).eval(md.points.row(1))); // v
            alpha_0.push_back(gluingData[0].get_alpha_S_tilde(0).eval(zero.row(0))); // u
            alpha_0.push_back(gluingData[1].get_alpha_S_tilde(0).eval(zero.row(0))); // v
            alpha_deriv.push_back(gluingData[0].get_alpha_S_tilde(0).deriv(zero.row(0))); // u
            alpha_deriv.push_back(gluingData[1].get_alpha_S_tilde(0).deriv(zero.row(0))); // v

            beta.push_back(gluingData[0].get_beta_S_tilde(0).eval(md.points.row(0))); // u
            beta.push_back(gluingData[1].get_beta_S_tilde(0).eval(md.points.row(1))); // v
            beta_0.push_back(gluingData[0].get_beta_S_tilde(0).eval(zero.row(0))); // u
            beta_0.push_back(gluingData[1].get_beta_S_tilde(0).eval(zero.row(0))); // v
            beta_deriv.push_back(gluingData[0].get_beta_S_tilde(0).deriv(zero.row(0))); // u
            beta_deriv.push_back(gluingData[1].get_beta_S_tilde(0).deriv(zero.row(0))); // v
        }
        else if (g1OptionList.getInt("gluingData") == gluingData::exact)
*/        {
        gsMatrix < T > temp_mat;
        m_gD[0].eval_alpha_into(points.row(0), temp_mat);
        alpha.push_back(temp_mat); // u
        m_gD[1].eval_alpha_into(points.row(1), temp_mat);
        alpha.push_back(temp_mat); // v

        m_gD[0].eval_alpha_into(zero.row(0), temp_mat);
        alpha_0.push_back(temp_mat); // u
        m_gD[1].eval_alpha_into(zero.row(0), temp_mat);
        alpha_0.push_back(temp_mat); // v

        m_gD[0].deriv_alpha_into(zero.row(0), temp_mat);
        alpha_deriv.push_back(temp_mat); // u
        m_gD[1].deriv_alpha_into(zero.row(0), temp_mat);
        alpha_deriv.push_back(temp_mat); // v

        m_gD[0].eval_beta_into(points.row(0), temp_mat);
        beta.push_back(temp_mat); // u
        m_gD[1].eval_beta_into(points.row(1), temp_mat);
        beta.push_back(temp_mat); // v

        m_gD[0].eval_beta_into(zero.row(0), temp_mat);
        beta_0.push_back(temp_mat); // u
        m_gD[1].eval_beta_into(zero.row(0), temp_mat);
        beta_0.push_back(temp_mat); // v

        m_gD[0].deriv_beta_into(zero.row(0), temp_mat);
        beta_deriv.push_back(temp_mat); // u
        m_gD[1].deriv_beta_into(zero.row(0), temp_mat);
        beta_deriv.push_back(temp_mat); // v
    }

    // Compute dd^^(i_k) and dd^^(i_k-1)
    gsMatrix<> dd_ik_plus, dd_ik_minus;
    gsMatrix<> dd_ik_minus_deriv, dd_ik_plus_deriv;
    dd_ik_minus = -1/(alpha_0[0](0,0)) * (geo.jacobian(zero).col(1) +
        beta_0[0](0,0) * geo.jacobian(zero).col(0));

    dd_ik_plus = 1/(alpha_0[1](0,0)) * (geo.jacobian(zero).col(0) +
        beta_0[1](0,0) * geo.jacobian(zero).col(1));

    gsMatrix<> geo_deriv2_12(2,1), geo_deriv2_11(2,1), geo_deriv2_22(2,1);
    geo_deriv2_12.row(0) = geo.deriv2(zero).row(2);
    geo_deriv2_12.row(1) = geo.deriv2(zero).row(5);
    geo_deriv2_11.row(0) = geo.deriv2(zero).row(0);
    geo_deriv2_11.row(1) = geo.deriv2(zero).row(3);
    geo_deriv2_22.row(0) = geo.deriv2(zero).row(1);
    geo_deriv2_22.row(1) = geo.deriv2(zero).row(4);
    gsMatrix<> alpha_squared_u = alpha_0[0]*alpha_0[0];
    gsMatrix<> alpha_squared_v = alpha_0[1]*alpha_0[1];

    dd_ik_minus_deriv = -1/(alpha_squared_u(0,0)) * // N^2
        ((geo_deriv2_12 + (beta_deriv[0](0,0) * geo.jacobian(zero).col(0) +
            beta_0[0](0,0) * geo_deriv2_11))*alpha_0[0](0,0) -
            (geo.jacobian(zero).col(1) + beta_0[0](0,0) * geo.jacobian(zero).col(0)) *
                alpha_deriv[0](0,0));

    dd_ik_plus_deriv = 1/(alpha_squared_v(0,0)) *
        ((geo_deriv2_12 + (beta_deriv[1](0,0) * geo.jacobian(zero).col(1) +
            beta_0[1](0,0) * geo_deriv2_22))*alpha_0[1](0,0) -
            (geo.jacobian(zero).col(0) + beta_0[1](0,0) * geo.jacobian(zero).col(1)) *
                alpha_deriv[1](0,0));

    // Comupute d_(0,0)^(i_k), d_(1,0)^(i_k), d_(0,1)^(i_k), d_(1,1)^(i_k) ; i_k == 2
    std::vector<gsMatrix<>> d_ik;
    d_ik.push_back(Phi.col(0));
    d_ik.push_back(Phi.block(0,1,6,2) * geo.jacobian(zero).col(0) ); // deriv into u
    d_ik.push_back(Phi.block(0,1,6,2) * geo.jacobian(zero).col(1) ); // deriv into v
    d_ik.push_back((geo.jacobian(zero)(0,0) * Phi.col(3) + geo.jacobian(zero)(1,0) * Phi.col(4))*geo.jacobian(zero)(0,1) +
        (geo.jacobian(zero)(0,0) * Phi.col(4) + geo.jacobian(zero)(1,0) * Phi.col(5))*geo.jacobian(zero)(1,1) +
        Phi.block(0,1,6,1) * geo.deriv2(zero).row(2) +
        Phi.block(0,2,6,1) * geo.deriv2(zero).row(5)); // Hessian

    // Compute d_(*,*)^(il,ik)
    std::vector<gsMatrix<>> d_ilik_minus, d_ilik_plus;
    d_ilik_minus.push_back(Phi.col(0));
    d_ilik_minus.push_back(Phi.block(0,1,6,2) * geo.jacobian(zero).col(0));
    d_ilik_minus.push_back((geo.jacobian(zero)(0,0) * Phi.col(3) + geo.jacobian(zero)(1,0) * Phi.col(4))*geo.jacobian(zero)(0,0) +
        (geo.jacobian(zero)(0,0) * Phi.col(4) + geo.jacobian(zero)(1,0) * Phi.col(5))*geo.jacobian(zero)(1,0) +
        Phi.block(0,1,6,1) * geo.deriv2(zero).row(0) +
        Phi.block(0,2,6,1) * geo.deriv2(zero).row(3));
    d_ilik_minus.push_back(Phi.block(0,1,6,2) * dd_ik_minus);
    d_ilik_minus.push_back((geo.jacobian(zero)(0,0) * Phi.col(3) + geo.jacobian(zero)(1,0) * Phi.col(4))*dd_ik_minus(0,0) +
        (geo.jacobian(zero)(0,0) * Phi.col(4) + geo.jacobian(zero)(1,0) * Phi.col(5))*dd_ik_minus(1,0) +
        Phi.block(0,1,6,1) * dd_ik_minus_deriv.row(0) +
        Phi.block(0,2,6,1) * dd_ik_minus_deriv.row(1));

    d_ilik_plus.push_back(Phi.col(0));
    d_ilik_plus.push_back(Phi.block(0,1,6,2) * geo.jacobian(zero).col(1));
    d_ilik_plus.push_back((geo.jacobian(zero)(0,1) * Phi.col(3) + geo.jacobian(zero)(1,1) * Phi.col(4))*geo.jacobian(zero)(0,1) +
        (geo.jacobian(zero)(0,1) * Phi.col(4) + geo.jacobian(zero)(1,1) * Phi.col(5))*geo.jacobian(zero)(1,1) +
        Phi.block(0,1,6,1) * geo.deriv2(zero).row(1) +
        Phi.block(0,2,6,1) * geo.deriv2(zero).row(4));
    d_ilik_plus.push_back(Phi.block(0,1,6,2) * dd_ik_plus);
    d_ilik_plus.push_back((geo.jacobian(zero)(0,1) * Phi.col(3) + geo.jacobian(zero)(1,1) * Phi.col(4))*dd_ik_plus(0,0) +
        (geo.jacobian(zero)(0,1) * Phi.col(4) + geo.jacobian(zero)(1,1) * Phi.col(5))*dd_ik_plus(1,0) +
        Phi.block(0,1,6,1) * dd_ik_plus_deriv.row(0) +
        Phi.block(0,2,6,1) * dd_ik_plus_deriv.row(1));

    gsMatrix<T> fpts(4,6);
    for (index_t i = 0; i < 6; i++)
    {
        fpts(0,i) = d_ilik_minus.at(0)(i,0); // f*_(ik-1,ik)
        fpts(0,i) += d_ilik_plus.at(0)(i,0) ; // f*_(ik+1,ik)
        fpts(0,i) -= d_ik.at(0)(i,0); // f*_(ik)

        // partial_u
        fpts(1,i) = d_ilik_minus.at(1)(i,0); // f*_(ik-1,ik)
        fpts(1,i) += -1 * d_ilik_plus.at(1)(i,0) * beta[1](0,0) +  d_ilik_plus.at(3)(i,0) * alpha[1](0,0) ; // f*_(ik+1,ik)
        fpts(1,i) -= d_ik.at(1)(i,0); // f*_(ik)

        // partial_v
        fpts(2,i) = -1 * d_ilik_minus.at(1)(i,0) * beta[0](0,0) - d_ilik_minus.at(3)(i,0) * alpha[0](0,0) ; // f*_(ik-1,ik)
        fpts(2,i) += d_ilik_plus.at(1)(i,0) ; // f*_(ik+1,ik)
        fpts(2,i) -= d_ik.at(2)(i,0); // f*_(ik)

        // partial_u partial_v
        fpts(3,i) = -1 * d_ilik_minus.at(1)(i,0) * beta_deriv[0](0,0) - d_ilik_minus.at(2)(i,0) * beta[0](0,0) -
            d_ilik_minus.at(3)(i,0) * alpha_deriv[0](0,0) - d_ilik_minus.at(4)(i,0) * alpha[0](0,0); // f*_(ik-1,ik)
        fpts(3,i) += -1 * d_ilik_plus.at(1)(i,0) * beta_deriv[1](0,0) -  d_ilik_plus.at(2)(i,0) * beta[1](0,0) +
            d_ilik_plus.at(3)(i,0) * alpha_deriv[1](0,0) + d_ilik_plus.at(4)(i,0) * alpha[1](0,0); // f*_(ik+1,ik)
        fpts(3,i) -= d_ik.at(3)(i,0); // f*_(ik)
    }

    gsSparseMatrix<real_t> mat(4,4);
    mat.setZero();
    gsVector<index_t> bindex(4);
    bindex.at(0) = 0;
    bindex.at(1) = 1;
    bindex.at(2) = m_basis_g1.basis(0).component(0).size();
    bindex.at(3) = m_basis_g1.basis(0).component(0).size()+1;

    for (index_t i = 0; i < 4; i++)
    {
        const int ii = bindex.at(i);
        mat(0,i) = basis_geo.evalSingle(ii,points)(0,0);
        mat(1,i) = basis_geo.derivSingle(ii,points)(0,0);
        mat(2,i) = basis_geo.derivSingle(ii,points)(1,0);
        mat(3,i) = basis_geo.deriv2Single(ii,points)(2,0);
    }

    // Save corresponding boundary dofs
    for (index_t l = 0; l < 6; ++l)
    {
        gsSparseSolver<real_t>::LU solver;

        solver.compute(mat);
        gsVector<real_t> sol = solver.solve(fpts.col(l));
        for (index_t i = 0; i < 4; i++)
        {
            const int ii = m_f[l].colMapper(0).bindex(bindex.at(i), 0);
            m_ddof[l](ii, 0) = sol(i);
        }
    }


    }

} // namespace gismo