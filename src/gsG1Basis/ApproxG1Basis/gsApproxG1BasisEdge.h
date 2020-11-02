/** @file gsApproxG1BasisEdge.h

    @brief Provides assembler for a G1 Basis for multiPatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once


#include <gsG1Basis/ApproxG1Basis/gsApproxGluingData.h>
#include <gsG1Basis/ApproxG1Basis/gsVisitorApproxG1BasisEdge.h>

# include <gsG1Basis/gsG1OptionList.h>

//# include <gsG1Basis/gsApproxSingleEdgeAssembler.h>

# include "gsG1Basis/Norm/gsNormL2BasisFunction.h"

namespace gismo
{
template<class T, class bhVisitor = gsVisitorApproxG1BasisEdge<T>>
class gsApproxG1BasisEdge : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsApproxG1BasisEdge(gsMultiPatch<> mp, // single patch
                  gsMultiBasis<> basis, // single basis
                  std::vector<gsBSplineBasis<>> basis_pm,
                  gsApproxGluingData<T> gluingData,
                  index_t uv, // !!! 0 == u; 1 == v !!!
                  bool isBoundary,
                  gsG1OptionList & g1OptionList)
        : m_mp(mp), m_basis(basis), m_uv(uv), m_isBoundary(isBoundary), m_g1OptionList(g1OptionList)
    {

        m_basis_plus = basis_pm[0];
        m_basis_minus = basis_pm[1];

        // Computing the gluing data
/*        gsApproxGluingData<T> gluingData(m_mp, m_basis, m_uv, m_isBoundary, m_g1OptionList);
        if (g1OptionList.getInt("gluingData") == gluingData::local)
            gluingData.setLocalGluingData(m_basis_plus, m_basis_minus, "edge");
        else if (g1OptionList.getInt("gluingData") == gluingData::global)
            gluingData.setGlobalGluingData();
*/
        m_gD.push_back(gluingData);

        // Basis for the G1 basis
        m_basis_g1 = m_basis.basis(0);
        //m_basis_g1.degreeReduce(1);


    }

    // Computed the gluing data globally
    void setG1BasisEdge(gsMultiPatch<T> & result);

    void refresh(index_t bfID, std::string typeBf);
    void assemble(index_t i, std::string typeBf); // i == number of bf
    inline void apply(bhVisitor & visitor, index_t i, std::string typeBf); // i == number of bf

    void computeDofsIntpl(index_t bfID, std::string typeBf);

    void constructSolution(const gsMatrix<> & solVector, gsMultiPatch<T> & result);

    gsBSpline<> get_alpha() { return m_gD[0].get_alpha_tilde(); }
    gsBSpline<> get_beta() { return m_gD[0].get_beta_tilde(); }

    gsMatrix<> eval_beta(gsMatrix<> points) {
        gsMatrix<> result;
        m_gD[0].eval_beta_into(points,result);
        return result;
    }

    gsMatrix<> eval_alpha(gsMatrix<> points) {
        gsMatrix<> result;
        m_gD[0].eval_alpha_into(points,result);
        return result;
    }

    void set_beta_tilde(gsBSpline<T> beta_t) { m_gD[0].set_beta_tilde(beta_t); }

    void plotGluingData(index_t numGd) { m_gD[0].plotGluingData(numGd); }

    index_t get_plus() { return m_basis_plus.size(); }

    real_t get_error() { return m_error; }

protected:

    // Input
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_basis;
    index_t m_uv;
    bool m_isBoundary;
    gsG1OptionList m_g1OptionList;

    // Gluing data
    std::vector<gsApproxGluingData<T>> m_gD;

    // Basis for getting the G1 Basis
    gsBSplineBasis<> m_basis_plus;
    gsBSplineBasis<> m_basis_minus;

    // Basis for the G1 Basis
    gsMultiBasis<T> m_basis_g1;

    // Basis for Integration
    gsMultiBasis<T> m_geo;

    // For Dirichlet boundary
    using Base::m_ddof;
    using Base::m_system;

    // For special projection
    gsBSpline<> result_singleEdge;

    real_t m_error;

}; // class gsG1BasisEdge

template <class T, class bhVisitor>
void gsApproxG1BasisEdge<T,bhVisitor>::setG1BasisEdge(gsMultiPatch<T> & result)
{
    result.clear();

    index_t n_plus = m_basis_plus.size();
    index_t n_minus = m_basis_minus.size();

    float progress = 0.0;
    int barWidth = 50;

    gsMultiPatch<> g1EdgeBasis, g1EdgeBasis2;
    index_t bfID_init = 3;
    if (m_g1OptionList.getSwitch("twoPatch"))
        bfID_init = 0;

    for (index_t bfID = bfID_init; bfID < n_plus - bfID_init; bfID++) // first 3 and last 3 bf are eliminated
    {
        if (m_g1OptionList.getSwitch("info"))
        {
            std::cout << "[";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i) {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " %\r";
            std::cout.flush();
        }


        gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(m_uv)); // u
        index_t degree = temp_basis_first.maxDegree();

        gsMatrix<T> ab = m_basis_plus.support(bfID);

        gsMatrix<T> ab_temp = ab;
        for (index_t pp = 0; pp < degree; pp++)
        {
            for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
            {
                gsMatrix<T> xy = temp_basis_first.support(i);
                if ( (xy(0,0) < ab(0,0)) && (xy(0,1) > ab(0,0)))
                    ab_temp(0,0) = xy(0,0);
                if ( (xy(0,0) < ab(0,1)) && (xy(0,1) > ab(0,1)))
                    ab_temp(0,1) = xy(0,1);
            }
            ab = ab_temp;
        }

        gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, 1);
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i += temp_basis_first.knots().multiplicityIndex(i))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv.insert(temp_basis_first.knot(i), 1);

        temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(1 - m_uv)); // v
        ab = temp_basis_first.support(1);
        gsKnotVector<T> kv2(ab.at(0), ab.at(1), 0, 1);
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i += temp_basis_first.knots().multiplicityIndex(i))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv2.insert(temp_basis_first.knot(i), 1);

        gsTensorBSplineBasis<2, T> bsp_geo_local(kv, kv2);
        gsTensorBSplineBasis<2, T> temp_basis(kv2, kv);
        if (m_uv == 1)
            bsp_geo_local.swap(temp_basis);

        if (m_g1OptionList.getInt("g1BasisEdge") == g1BasisEdge::local)
            m_geo = bsp_geo_local; // Basis for Integration
        else
            m_geo = m_basis_g1;

        //m_geo = m_basis_g1;
//// NEUNEUNEU
/*
        // Compute the assembler
        gsApproxSingleEdgeAssembler<real_t> approxSingleEdgeAssembler(m_mp, m_gD[0], m_g1OptionList, m_uv, bfID);

        gsSparseSolver<real_t>::CGDiagonal solver2;
        gsMatrix<> sol2;
        solver2.compute(approxSingleEdgeAssembler.matrix());
        sol2 = solver2.solve(approxSingleEdgeAssembler.rhs());
        approxSingleEdgeAssembler.constructSolution(sol2, result_singleEdge);
*/
//// NEUNEUNEU


        refresh(bfID,"plus");

        assemble(bfID,"plus"); // i == number of bf

        gsSparseSolver<real_t>::CGDiagonal solver;
        gsMatrix<> sol;
        solver.compute(m_system.matrix());
        sol = solver.solve(m_system.rhs());
/*
        if (bfID == 0)
        {
            //gsInfo << "Mat: " << m_system.matrix().toDense() << "\n";
            gsInfo << "RHS: " << m_system.rhs() << "\n";
            gsInfo << "sol: " << sol << "\n";
        }
*/

        //std::cout << "#iterations:     " << solver.iterations() << std::endl;
        //std::cout << "estimated error: " << solver.error()      << std::endl;

        constructSolution(sol,g1EdgeBasis);

        progress += (float) 1.0 / (n_plus+n_minus); // for demonstration only
    }


    bfID_init = 2;
    if (m_g1OptionList.getSwitch("twoPatch"))
        bfID_init = 0;

    for (index_t bfID = bfID_init; bfID < n_minus-bfID_init; bfID++)  // first 2 and last 2 bf are eliminated
    {
        if (m_g1OptionList.getSwitch("info"))
        {
            std::cout << "[";
            int pos = barWidth * progress;
            for (int i = 0; i < barWidth; ++i)
            {
                if (i < pos) std::cout << "=";
                else if (i == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] " << int(progress * 100.0) << " %\r";
            std::cout.flush();
        }

        gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(m_uv)); // u
        index_t degree = temp_basis_first.maxDegree();

        gsMatrix<T> ab = m_basis_minus.support(bfID);

        gsMatrix<T> ab_temp = ab;
        for (index_t pp = 0; pp < degree; pp++)
        {
            for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
            {
                gsMatrix<T> xy = temp_basis_first.support(i);
                if ( (xy(0,0) < ab(0,0)) && (xy(0,1) > ab(0,0)))
                    ab_temp(0,0) = xy(0,0);
                if ( (xy(0,0) < ab(0,1)) && (xy(0,1) > ab(0,1)))
                    ab_temp(0,1) = xy(0,1);
            }
            ab = ab_temp;
        }
/*
        if (!m_isBoundary)
        {
            gsMatrix<T> ab_temp = ab;
            for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
            {
                gsMatrix<T> xy = temp_basis_first.support(i);
                if ( (xy(0,0) < ab(0,0)) && (xy(0,1) > ab(0,0)))
                    ab_temp(0,0) = xy(0,0);
                if ( (xy(0,0) < ab(0,1)) && (xy(0,1) > ab(0,1)))
                    ab_temp(0,1) = xy(0,1);
            }
            ab = ab_temp;

            for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
            {
                gsMatrix<T> xy = temp_basis_first.support(i);
                if ( (xy(0,0) < ab(0,0)) && (xy(0,1) > ab(0,0)))
                    ab_temp(0,0) = xy(0,0);
                if ( (xy(0,0) < ab(0,1)) && (xy(0,1) > ab(0,1)))
                    ab_temp(0,1) = xy(0,1);
            }
            ab = ab_temp;
        }
*/
        gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, 1);
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i += temp_basis_first.knots().multiplicityIndex(i))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv.insert(temp_basis_first.knot(i), 1);


        temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(1 - m_uv)); // v
        ab = temp_basis_first.support(1);
        gsKnotVector<T> kv2(ab.at(0), ab.at(1), 0, 1);
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i += temp_basis_first.knots().multiplicityIndex(i))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv2.insert(temp_basis_first.knot(i), 1);
        gsTensorBSplineBasis<2, T> bsp_geo_local(kv, kv2);
        gsTensorBSplineBasis<2, T> temp_basis(kv2, kv);
        if (m_uv == 1)
            bsp_geo_local.swap(temp_basis);

        if (m_g1OptionList.getInt("g1BasisEdge") == g1BasisEdge::local)
            m_geo = bsp_geo_local; // Basis for Integration
        else
            m_geo = m_basis_g1;

        refresh(bfID,"minus");

        assemble(bfID,"minus"); // i == number of bf

        gsSparseSolver<real_t>::CGDiagonal solver;
        gsMatrix<> sol;
        solver.compute(m_system.matrix());
        sol = solver.solve(m_system.rhs());

        constructSolution(sol,g1EdgeBasis2);
        g1EdgeBasis.addPatch(g1EdgeBasis2.patch(bfID));

        progress += (float) 1.0 / (n_plus+n_minus); // for demonstration only
    }


    result = g1EdgeBasis;
    /*
    gsInfo << "BEfore: " << result.basis(0) << "\n";

    gsMultiPatch<T> result_refined;
    for ( size_t i = 0; i < result.nPatches(); i++ )
    {

        gsMatrix<T> iVals, iPts = m_basis.basis(0).anchors();
        result.patch(i).eval_into(iPts, iVals);
        typename gsGeometry<T>::uPtr g = m_basis.basis(0).interpolateData(iVals, iPts);
        result_refined.addPatch(m_basis.basis(0).interpolateData(iVals, iPts));
    }
    result.swap(result_refined);
    gsInfo << "result: " << result.basis(0) << "\n";
    */

    if (m_g1OptionList.getSwitch("info"))
        std::cout << std::endl;


    gsNormL2BasisFunction<T> normL2BasisFunction(g1EdgeBasis2, m_basis, m_basis_plus, m_basis_minus, m_gD[0], m_uv);
    normL2BasisFunction.compute();

    m_error = normL2BasisFunction.value();
    if (m_g1OptionList.getSwitch("info"))
        gsInfo << "\nVALUE: " << normL2BasisFunction.value() << "\n\n";
} // setG1BasisEdge

template <class T, class bhVisitor>
void gsApproxG1BasisEdge<T,bhVisitor>::constructSolution(const gsMatrix<> & solVector, gsMultiPatch<T> & result)
{
    // Dim is the same for all basis functions
    const index_t dim = ( 0!=solVector.cols() ? solVector.cols() :  m_ddof[0].cols() );

    gsMatrix<T> coeffs;
    const gsDofMapper & mapper = m_system.colMapper(0); // unknown = 0

    // Reconstruct solution coefficients on patch p
    index_t sz;
    sz = m_basis_g1.basis(0).size();

    coeffs.resize(sz, dim);

    for (index_t i = 0; i < sz; ++i)
    {
        if (mapper.is_free(i, 0)) // DoF value is in the solVector // 0 = unitPatch
        {
            coeffs.row(i) = solVector.row(mapper.index(i, 0));
        }
        else // eliminated DoF: fill with Dirichlet data
        {
            coeffs.row(i) = m_ddof[0].row( mapper.bindex(i, 0) ).head(dim); // = 0
        }
    }

    result.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));

}

template <class T, class bhVisitor>
void gsApproxG1BasisEdge<T,bhVisitor>::refresh(index_t bfID, std::string typeBf)
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map(m_basis_g1.basis(0));

    gsMatrix<unsigned> act;

    //index_t n_plus = m_basis_plus.size();
    //index_t n_minus = m_basis_minus.size();

    for (index_t i = 2; i < m_basis_g1.basis(0).component(1 - m_uv).size();
         i++) // only the first two u/v-columns are Dofs (0/1)
    {
        act = m_basis_g1.basis(0).boundaryOffset(m_uv == 0 ? 3 : 1, i); // WEST
        map.markBoundary(0, act); // Patch 0
    }

    if (m_g1OptionList.getInt("g1BasisEdge") == g1BasisEdge::local)
    {
        for (index_t i = 2; i < m_basis_g1.basis(0).component(1 - m_uv).size();
             i++) // only the first two u/v-columns are Dofs (0/1)
        {
            act = m_basis_g1.basis(0).boundaryOffset(m_uv == 0 ? 3 : 1, i); // WEST
            map.markBoundary(0, act); // Patch 0
        }

//        if (m_isBoundary)
//        {
            if (typeBf == "plus")
            {
                gsMatrix<T> ab = m_basis_plus.support(bfID);


                gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(m_uv)); // u
                index_t degree = temp_basis_first.maxDegree();

                gsMatrix<T> ab_temp = ab;
                for (index_t pp = 0; pp < degree; pp++)
                {
                    for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
                    {
                        gsMatrix<T> xy = temp_basis_first.support(i);
                        if ( (xy(0,0) < ab(0,0)) && (xy(0,1) > ab(0,0)))
                            ab_temp(0,0) = xy(0,0);
                        if ( (xy(0,0) < ab(0,1)) && (xy(0,1) > ab(0,1)))
                            ab_temp(0,1) = xy(0,1);
                    }
                    ab = ab_temp;
                }


                for (index_t i = 0; i < m_basis_g1.basis(0).component(m_uv).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    gsMatrix<T> xy = m_basis_g1.basis(0).component(m_uv).support(i);
                    if ( (xy(0, 0) < ab(0, 0)-1e-10) || (xy(0, 1) > ab(0, 1)+1e-10) ) // only subsets
                    //if ( (xy(0, 1) < ab(0, 0)+1e-10) || (xy(0, 0) > ab(0, 1)-1e-10) ) // all non-empty set
                    {
                        act = m_basis_g1.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, i); // WEST
                        map.markBoundary(0, act); // Patch 0
                    }
                }
/*
                if (bfID > 0 && bfID < n_plus - 1)
                {
                    // set first row to zero
                    act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, 0); // WEST
                    map.markBoundary(0, act); // Patch 0
                    act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 2 : 4, 0); // WEST
                    map.markBoundary(0, act); // Patch 0
                }
*/
            }
            else if (typeBf == "minus")
            {
                gsMatrix<T> ab = m_basis_minus.support(bfID);

                gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(m_uv)); // u
                index_t degree = temp_basis_first.maxDegree();

                gsMatrix<T> ab_temp = ab;
                for (index_t pp = 0; pp < degree; pp++)
                {
                    for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
                    {
                        gsMatrix<T> xy = temp_basis_first.support(i);
                        if ( (xy(0,0) < ab(0,0)) && (xy(0,1) > ab(0,0)))
                            ab_temp(0,0) = xy(0,0);
                        if ( (xy(0,0) < ab(0,1)) && (xy(0,1) > ab(0,1)))
                            ab_temp(0,1) = xy(0,1);
                    }
                    ab = ab_temp;
                }

                for (index_t i = 0; i < m_basis_g1.basis(0).component(m_uv).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    gsMatrix<T> xy = m_basis_g1.basis(0).component(m_uv).support(i);
                    if ( (xy(0, 0) < ab(0, 0)-1e-10) || (xy(0, 1) > ab(0, 1)+1e-10) ) // only subsets
                    {
                        act = m_basis_g1.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, i); // WEST
                        map.markBoundary(0, act); // Patch 0
                    }
                }

                // set first row to zero
                act = m_basis_g1.basis(0).boundaryOffset(m_uv == 0 ? 3 : 1, 0); // WEST
                map.markBoundary(0, act); // Patch 0
/*
                if (bfID > 0 && bfID < n_minus - 1)
                {
                    // set first row to zero
                    act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, 0); // WEST
                    map.markBoundary(0, act); // Patch 0
                    act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 2 : 4, 0); // WEST
                    map.markBoundary(0, act); // Patch 0
                }
*/

            }
/*        }
        else if (!m_isBoundary)
        {
            if (typeBf == "plus")
            {
                gsMatrix<T> ab = m_basis_plus.support(bfID);

                gsMatrix<T> ab_temp = ab;
                for (index_t i = 0; i < m_basis.basis(0).component(m_uv).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    gsMatrix<T> xy = m_basis.basis(0).component(m_uv).support(i);
                    if ((xy(0, 0) < ab(0, 0)) && (xy(0, 1) > ab(0, 0)))
                        ab_temp(0, 0) = xy(0, 0);
                    if ((xy(0, 0) < ab(0, 1)) && (xy(0, 1) > ab(0, 1)))
                        ab_temp(0, 1) = xy(0, 1);
                }
                ab = ab_temp;

                for (index_t i = 0; i < m_basis.basis(0).component(m_uv).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    gsMatrix<T> xy = m_basis.basis(0).component(m_uv).support(i);
                    if ((xy(0, 0) < ab(0, 0)) && (xy(0, 1) > ab(0, 0)))
                        ab_temp(0, 0) = xy(0, 0);
                    if ((xy(0, 0) < ab(0, 1)) && (xy(0, 1) > ab(0, 1)))
                        ab_temp(0, 1) = xy(0, 1);
                }
                ab = ab_temp;

                for (index_t i = 0; i < m_basis.basis(0).component(m_uv).size(); i++)
                {
                    gsMatrix<T> xy = m_basis.basis(0).component(m_uv).support(i);
                    if ((xy(0, 0) < ab(0, 0)) || (xy(0, 1) > ab(0, 1)))
                    {
                        act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, i); // WEST
                        map.markBoundary(0, act); // Patch 0
                    }
                }

                act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, 0); // WEST
                map.markBoundary(0, act); // Patch 0
                act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, m_basis.basis(0).component(m_uv).size() - 1);
                map.markBoundary(0, act); // Patch 0

            }
            else if (typeBf == "minus")
            {
                gsMatrix<T> ab = m_basis_minus.support(bfID);

                gsMatrix<T> ab_temp = ab;
                for (index_t i = 0; i < m_basis.basis(0).component(m_uv).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    gsMatrix<T> xy = m_basis.basis(0).component(m_uv).support(i);
                    if ((xy(0, 0) < ab(0, 0)) && (xy(0, 1) > ab(0, 0)))
                        ab_temp(0, 0) = xy(0, 0);
                    if ((xy(0, 0) < ab(0, 1)) && (xy(0, 1) > ab(0, 1)))
                        ab_temp(0, 1) = xy(0, 1);
                }
                ab = ab_temp;

                for (index_t i = 0; i < m_basis.basis(0).component(m_uv).size();
                     i++) // only the first two u/v-columns are Dofs (0/1)
                {
                    gsMatrix<T> xy = m_basis.basis(0).component(m_uv).support(i);
                    if ((xy(0, 0) < ab(0, 0)) && (xy(0, 1) > ab(0, 0)))
                        ab_temp(0, 0) = xy(0, 0);
                    if ((xy(0, 0) < ab(0, 1)) && (xy(0, 1) > ab(0, 1)))
                        ab_temp(0, 1) = xy(0, 1);
                }
                ab = ab_temp;

                for (index_t i = 0; i < m_basis.basis(0).component(m_uv).size(); i++)
                {
                    gsMatrix<T> xy = m_basis.basis(0).component(m_uv).support(i);
                    if ((xy(0, 0) < ab(0, 0)) || (xy(0, 1) > ab(0, 1)))
                    {
                        act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, i); // WEST
                        map.markBoundary(0, act); // Patch 0
                    }
                }

                act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, 0); // WEST
                map.markBoundary(0, act); // Patch 0
                act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, m_basis.basis(0).component(m_uv).size() - 1);
                map.markBoundary(0, act); // Patch 0

            }
       }
*/    }


    //act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, 0); // WEST
    //map.markBoundary(0, act); // Patch 0
    //act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, m_basis.basis(0).component(m_uv).size() - 1);
    //map.markBoundary(0, act); // Patch 0

    /*
     * x x o
     * x o o
     * ...
     * o o o
     * x o o
     * x x o
     */
    act.resize(6, 1);
    act(0, 0) = (unsigned) 0;
    act(1,0) = (unsigned) (m_uv == 0 ? m_basis.basis(0).component(m_uv).size() : 1);
    act(2,0) = (unsigned) (m_uv == 0 ? 1 : m_basis.basis(0).component(m_uv).size());
    act(3,0) = (unsigned) (m_uv == 0 ? (m_basis.basis(0).component(m_uv).size()* 2)-1: (m_basis.basis(0).component(m_uv).size()-1) * (m_basis.basis(0).component(1-m_uv).size()) +1);
    act(4,0) = (unsigned) (m_uv == 0 ? m_basis.basis(0).component(m_uv).size() - 1 : (m_basis.basis(0).component(m_uv).size()-1) * (m_basis.basis(0).component(1-m_uv).size()));
    act(5,0) = (unsigned) (m_uv == 0 ? m_basis.basis(0).component(m_uv).size() - 2 : (m_basis.basis(0).component(m_uv).size()-2) * (m_basis.basis(0).component(1-m_uv).size()));
    //map.markBoundary(0, act);

    map.finalize();

    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(map);

} // refresh()

template <class T, class bhVisitor>
void gsApproxG1BasisEdge<T,bhVisitor>::assemble(index_t bfID, std::string typeBf)
{
    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis[0],2,1,0.333333);
    m_system.reserve(nz, 1);

    if(m_ddof.size()==0)
        m_ddof.resize(1); // 0,1

    const gsDofMapper & map = m_system.colMapper(0); // Map same for every functions

    m_ddof[0].setZero(map.boundarySize(), 1 );

    //computeDofsIntpl(bfID, typeBf);

    // Compute dirichlet value
    // m_ddof[0].row(map.global_toBindex()) = value of function

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor, bfID, typeBf); // basis function i

    m_system.matrix().makeCompressed();

} // assemble()

template <class T, class bhVisitor>
void gsApproxG1BasisEdge<T,bhVisitor>::apply(bhVisitor & visitor, int bf_index, std::string typeBf)
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
        gsBasis<T> & basis_geo = m_basis.basis(0).component(1-m_uv);
        gsBasis<T> & basis_plus = m_basis_plus;
        gsBasis<T> & basis_minus = m_basis_minus;

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

            // Perform required evaluations on the quadrature nodes
            visitor_.evaluate(bf_index, typeBf, basis_g1, basis_geo, basis_plus, basis_minus, patch, quNodes, m_uv, m_gD[0], m_isBoundary, m_g1OptionList, result_singleEdge);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(0, m_ddof, m_system); // omp_locks inside // patchIndex == 0
        }
    }//omp parallel
} // apply


template <class T, class bhVisitor>
void gsApproxG1BasisEdge<T,bhVisitor>::computeDofsIntpl(index_t bfID, std::string typeBf)

{

    gsVector<index_t> sides(2); // Two boundaries
    sides.at(0) = m_uv == 1 ? 3 : 1;
    sides.at(1) = m_uv == 1 ? 4 : 2;

    const gsBasis<T> & basis = m_basis_g1.basis(0);

    for (index_t side = 0; side < sides.rows(); side++)
    {
        // Get dofs on this boundary
        gsMatrix<unsigned> boundary(3,1); // bottom boundary v = 0 or u = 0
        //const gsMatrix<unsigned> boundary = m_basis_g1.basis(0).boundaryOffset(3,1); // bottom boundary v = 0 or u = 0
        //gsInfo << "Boundary: " << boundary << "\n";
        if (side == 0)
        {
            // Value
            // D_(1-uv)
            // D_uv
            boundary(0, 0) = (unsigned) 0;
            boundary(1, 0) = (unsigned) (m_uv == 0 ? basis.component(m_uv).size() : 1);
            boundary(2, 0) = (unsigned) (m_uv == 0 ? 1 : basis.component(m_uv).size());
        }
        else
        {
            boundary(0, 0) =
                (unsigned) (m_uv == 0 ? basis.component(m_uv).size() - 1 : (basis.component(m_uv).size() - 1)
                    * (basis.component(1 - m_uv).size()));
            boundary(1, 0) =
                    (unsigned) (m_uv == 0 ? (basis.component(m_uv).size()* 2)-1: (basis.component(m_uv).size()-1) * (basis.component(1-m_uv).size()) +1);
            boundary(2, 0) =
                (unsigned) (m_uv == 0 ? basis.component(m_uv).size() - 2: (basis.component(m_uv).size() - 2)
                    * (basis.component(1 - m_uv).size()));
        }

        gsMatrix<> points(2,1);
        points.setZero();
        points(m_uv,0) = side == 0 ? 0 : 1;

        // #### Evaluate the points ####
        // Same for all patches
        gsBasis<T> & basis_geo = m_basis.basis(0).component(1 - m_uv);
        gsBasis<T> & basis_plus = m_basis_plus;
        gsBasis<T> & basis_minus = m_basis_minus;

        // tau/p
        gsBSplineBasis<T> bsp_temp = dynamic_cast<gsBSplineBasis<> & >(basis_geo);

        real_t p = basis_geo.maxDegree();
        real_t tau_1 = bsp_temp.knots().at(p + 1); // p + 2

        gsMatrix<T> alpha, beta,
            N_0, N_1,
            N_j_minus, N_i_plus,
            der_N_i_plus;

        gsMatrix<T> der_N_0, der_N_1, der_N_00, der_N_01, der_N_uv;

        basis.derivSingle_into(boundary(0, 0),points,der_N_00);
        basis.derivSingle_into(boundary(1, 0),points,der_N_01);
        basis.derivSingle_into(boundary(2, 0),points,der_N_uv);



        gsMatrix<T> fpts(3,1);
        if (m_uv == 1) // edge is in v-direction
        {
            if (m_g1OptionList.getInt("gluingData") == gluingData::global)
            {
                m_gD[0].get_alpha_S_tilde(1-m_uv).eval_into(points.bottomRows(1), alpha); // v
                m_gD[0].get_beta_S_tilde(1-m_uv).eval_into(points.bottomRows(1), beta);
            }

            basis_geo.evalSingle_into(0, points.topRows(1), N_0); // u
            basis_geo.evalSingle_into(1, points.topRows(1), N_1); // u

            basis_geo.derivSingle_into(0, points.topRows(1), der_N_0); // u
            basis_geo.derivSingle_into(1, points.topRows(1), der_N_1); // u

            // Initialize local matrix/rhs
            if (typeBf == "plus")
            {
                basis_plus.evalSingle_into(bfID, points.bottomRows(1), N_i_plus); // v
                basis_plus.derivSingle_into(bfID, points.bottomRows(1), der_N_i_plus);

                gsMatrix<> temp = beta.cwiseProduct(der_N_i_plus);

                fpts.row(0) = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p;
                if (der_N_01(1-m_uv,0) == 0)
                    fpts(1,0) = 0;
                else
                    fpts(1,0) = -fpts(0,0)*der_N_00(1-m_uv,0)/der_N_01(1-m_uv,0);

                if (der_N_uv(m_uv,0) == 0)
                    fpts(2,0) = 0;
                else
                    fpts.row(2) = ( der_N_i_plus - fpts.row(0)*der_N_00(m_uv,0) )/der_N_uv(m_uv,0);

            } // n_plus
            else if (typeBf == "minus")
            {
                basis_minus.evalSingle_into(bfID, points.bottomRows(1), N_j_minus); // v


                fpts.row(0) = alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p; // == 0
                if (der_N_01(1-m_uv,0) == 0)
                    fpts(1,0) = 0;
                else
                {
                    fpts.row(1) = (alpha.cwiseProduct(N_j_minus.cwiseProduct(der_N_1)) * tau_1 / p)* 1/der_N_01(1-m_uv,0);
                }
                fpts(2,0) = 0;

            } // n_minus

        } // Patch 0
        else if (m_uv == 0) // edge is in u-direction
        {
            if (m_g1OptionList.getInt("gluingData") == gluingData::global)
            {
                m_gD[0].get_alpha_S_tilde(1-m_uv).eval_into(points.topRows(1), alpha); // u
                m_gD[0].get_beta_S_tilde(1-m_uv).eval_into(points.topRows(1), beta);
            }

            basis_geo.evalSingle_into(0, points.bottomRows(1), N_0); // v
            basis_geo.evalSingle_into(1, points.bottomRows(1), N_1); // v

            basis_geo.derivSingle_into(0, points.bottomRows(1), der_N_0); // v
            basis_geo.derivSingle_into(1, points.bottomRows(1), der_N_1); // v

            // Initialize local matrix/rhs
            if (typeBf == "plus")
            {
                basis_plus.evalSingle_into(bfID, points.topRows(1), N_i_plus); // u
                basis_plus.derivSingle_into(bfID, points.topRows(1), der_N_i_plus);

                gsMatrix<T> temp = beta.cwiseProduct(der_N_i_plus);

                //gsInfo << "uv = 0 : " << temp - result_singleEdge.eval(md.points.topRows(1)) << "\n";
                //temp = result_singleEdge.eval(md.points.topRows(1));
                fpts.row(0) = N_i_plus.cwiseProduct(N_0 + N_1) - temp.cwiseProduct(N_1) * tau_1 / p;
                if (der_N_01(1-m_uv,0) == 0)
                    fpts(1,0) = 0;
                else
                    fpts(1,0) = -fpts(0,0)*der_N_00(1-m_uv,0)/der_N_01(1-m_uv,0);

                if (der_N_uv(m_uv,0) == 0)
                    fpts(2,0) = 0;
                else
                    fpts.row(2) = ( der_N_i_plus - fpts.row(0)*der_N_00(m_uv,0) )/der_N_uv(m_uv,0);

            } // n_tilde
            else if (typeBf == "minus")
            {
                basis_minus.evalSingle_into(bfID, points.topRows(1), N_j_minus); // u

                fpts.row(0) = -alpha.cwiseProduct(N_j_minus.cwiseProduct(N_1)) * tau_1 / p; // == 0
                if (der_N_01(1-m_uv,0) == 0)
                    fpts(1,0) = 0;
                else
                {
                    fpts.row(1) = (-alpha.cwiseProduct(N_j_minus.cwiseProduct(der_N_1)) * tau_1 / p )* 1/der_N_01(1-m_uv,0);

                }
                fpts(2,0) = 0;

            } // n_bar

        } // Patch 1

/*
    // Compute dirichlet values
    gsMatrix<T> fpts;
    if ( it->parametric() )
        fpts = it->function()->eval( gsPointGrid<T>( rr ) );
    else
        fpts = it->function()->eval( m_pde_ptr->domain()[it->patch()].eval(  gsPointGrid<T>( rr ) ) );
*/

        //gsInfo << "typeBf " << typeBf << " bfid " << bfID << " " <<
        //    " points " << points << " fpts: " << fpts << "\n";

        // Interpolate dirichlet boundary
        //typename gsBasis<T>::uPtr h = basis.boundaryBasis(sides.at(side));
        //typename gsGeometry<T>::uPtr geo = h->interpolateAtAnchors(fpts);
        //const gsMatrix<T> & dVals = geo->coefs();

        //gsInfo << fpts << "\n";

        // Save corresponding boundary dofs
        for (index_t l = 0; l != boundary.size(); ++l)
        {
            const int ii = m_system.colMapper(0).bindex(boundary.at(l), 0);
            m_ddof[0].row(ii) = fpts.row(l);
        }
    }
}




//void gluingDataCondition(gsBSpline<> alpha_0, gsBSpline<> alpha_1, gsBSpline<> beta_0, gsBSpline<> beta_1)
//{
//    // BETA
//    // first,last,interior,mult_ends,mult_interior,degree
//    gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(auxGeom[0].getPatch().basis().component(1)); // 0 -> v, 1 -> u
//    index_t m_p = basis_edge.maxDegree(); // Minimum degree at the interface // TODO if interface basis are not the same
//
//    gsKnotVector<> kv(0, 1, basis_edge.numElements()-1, 2 * m_p  + 1, 2 * m_p - 1 );
//    gsBSplineBasis<> bsp(kv);
//
//    gsMatrix<> greville = bsp.anchors();
//    gsMatrix<> uv1, uv0, ev1, ev0;
//
//    const index_t d = 2;
//    gsMatrix<> D0(d,d);
//
//    gsGeometry<>::Ptr beta_temp;
//
//    uv0.setZero(2,greville.cols());
//    uv0.bottomRows(1) = greville;
//
//    uv1.setZero(2,greville.cols());
//    uv1.topRows(1) = greville;
//
//    const gsGeometry<> & P0 = auxGeom[0].getPatch(); // iFace.first().patch = 1
//    const gsGeometry<> & P1 = auxGeom[1].getPatch(); // iFace.second().patch = 0
//    // ======================================
//
//    // ======== Determine bar{beta} ========
//    for(index_t i = 0; i < uv1.cols(); i++)
//    {
//        P0.jacobian_into(uv0.col(i),ev0);
//        P1.jacobian_into(uv1.col(i),ev1);
//
//        D0.col(1) = ev0.col(0); // (DuFL, *)
//        D0.col(0) = ev1.col(1); // (*,DuFR)
//
//        uv0(0,i) = D0.determinant();
//    }
//
//    beta_temp = bsp.interpolateData(uv0.topRows(1), uv0.bottomRows(1));
//    gsBSpline<> beta = dynamic_cast<gsBSpline<> &> (*beta_temp);
//
//
//    index_t p_size = 10000;
//    gsMatrix<> points(1, p_size);
//    points.setRandom();
//    points = points.array().abs();
//
//    gsVector<> vec;
//    vec.setLinSpaced(p_size,0,1);
//    points = vec.transpose();
//
//    gsMatrix<> temp;
//    temp = alpha_1.eval(points).cwiseProduct(beta_0.eval(points))
//        + alpha_0.eval(points).cwiseProduct(beta_1.eval(points))
//        - beta.eval(points);
//
//
//    gsInfo << "Conditiontest Gluing data: \n" << temp.array().abs().maxCoeff() << "\n\n";
//
//
//}
//
//void g1ConditionRep(gsBSpline<> alpha_0, gsBSpline<> alpha_1, gsMultiPatch<> g1Basis_0,  gsMultiPatch<> g1Basis_1)
//{
//    // BETA
//    // first,last,interior,mult_ends,mult_interior,degree
//    gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(auxGeom[0].getPatch().basis().component(1)); // 0 -> v, 1 -> u
//    index_t m_p = basis_edge.maxDegree(); // Minimum degree at the interface // TODO if interface basis are not the same
//
//    gsKnotVector<> kv(0, 1, basis_edge.numElements()-1, 2 * m_p  + 1, 2 * m_p - 1 );
//    gsBSplineBasis<> bsp(kv);
//
//    gsMatrix<> greville = bsp.anchors();
//    gsMatrix<> uv1, uv0, ev1, ev0;
//
//    const index_t d = 2;
//    gsMatrix<> D0(d,d);
//
//    gsGeometry<>::Ptr beta_temp;
//
//    uv0.setZero(2,greville.cols());
//    uv0.bottomRows(1) = greville;
//
//    uv1.setZero(2,greville.cols());
//    uv1.topRows(1) = greville;
//
//    const gsGeometry<> & P0 = auxGeom[0].getPatch(); // iFace.first().patch = 1
//    const gsGeometry<> & P1 = auxGeom[1].getPatch(); // iFace.second().patch = 0
//    // ======================================
//
//    // ======== Determine bar{beta} ========
//    for(index_t i = 0; i < uv1.cols(); i++)
//    {
//        P0.jacobian_into(uv0.col(i),ev0);
//        P1.jacobian_into(uv1.col(i),ev1);
//        D0.col(1) = ev0.col(0); // (DuFL, *)
//        D0.col(0) = ev1.col(1); // (*,DuFR)
//
//        uv0(0,i) = D0.determinant();
//    }
//
//    beta_temp = bsp.interpolateData(uv0.topRows(1), uv0.bottomRows(1));
//    gsBSpline<> beta = dynamic_cast<gsBSpline<> &> (*beta_temp);
//
//
//
//    index_t p_size = 10000;
//    gsMatrix<> points(1, p_size);
//    points.setRandom();
//    points = points.array().abs();
//
//    gsVector<> vec;
//    vec.setLinSpaced(p_size,0,1);
//    points = vec.transpose();
//
//    gsMatrix<> points2d_0(2, p_size);
//    gsMatrix<> points2d_1(2, p_size);
//
//    points2d_0.setZero();
//    points2d_1.setZero();
//    points2d_0.row(1) = points; // v
//    points2d_1.row(0) = points; // u
//
//    real_t g1Error = 0;
//
//    for (size_t i = 0; i < g1Basis_0.nPatches(); i++)
//    {
//        gsMatrix<> temp;
//        temp = alpha_1.eval(points).cwiseProduct(g1Basis_0.patch(i).deriv(points2d_0).topRows(1))
//            + alpha_0.eval(points).cwiseProduct(g1Basis_1.patch(i).deriv(points2d_1).bottomRows(1))
//            + beta.eval(points).cwiseProduct(g1Basis_0.patch(i).deriv(points2d_0).bottomRows(1));
//
//        if (temp.array().abs().maxCoeff() > g1Error)
//            g1Error = temp.array().abs().maxCoeff();
//    }
//
//    gsInfo << "Conditiontest G1 continuity Rep: \n" << g1Error << "\n\n";
//}

} // namespace gismo