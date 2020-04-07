/** @file gsApproxG1BasisEdge.h

    @brief Provides assembler for a G1 Basis for multiPatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

#include <gsG1Basis/gsGluingData.h>
#include <gsG1Basis/gsApproxGluingData.h>
#include <gsG1Basis/gsG1ASGluingData.h>
#include <gsG1Basis/gsVisitorApproxG1BasisEdge.h>
# include <gsAssembler/gsAssembler.h>
# include <gsG1Basis/gsG1OptionList.h>

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
                  index_t uv, // !!! 0 == u; 1 == v !!!
                  bool isBoundary,
                  gsG1OptionList & g1OptionList)
        : m_mp(mp), m_basis(basis), m_uv(uv), m_isBoundary(isBoundary), m_g1OptionList(g1OptionList)
    {

        // Computing the G1 - basis function at the edge
        // Spaces for computing the g1 basis
        index_t m_r = m_g1OptionList.getInt("regularity"); // TODO CHANGE IF DIFFERENT REGULARITY IS NECESSARY

        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(m_uv)); // 0 -> v, 1 -> u
        index_t m_p = basis_edge.maxDegree(); // Minimum degree at the interface // TODO if interface basis are not the same

        // first,last,interior,mult_ends,mult_interior
        gsKnotVector<T> kv_plus(0,1,0,m_p+1,m_p-1-m_r); // p,r+1 //-1 bc r+1
        gsBSplineBasis<> basis_plus(kv_plus);

        // TODO if interface basis are not the same DO NOT DELETE THE FOLLOWING LINES
        /*
        if (basis_1.numElements() <= basis_2.numElements()) //
            for (size_t i = basis_1.degree()+1; i < basis_1.knots().size() - (basis_1.degree()+1); i = i+(basis_1.degree()-m_r))
                basis_plus.insertKnot(basis_1.knot(i),m_p-1-m_r);
        else
            for (size_t i = basis_2.degree()+1; i < basis_2.knots().size() - (basis_2.degree()+1); i = i+(basis_2.degree()-m_r))
                basis_plus.insertKnot(basis_2.knot(i),m_p-1-m_r);
        */
        for (size_t i = m_p+1; i < basis_edge.knots().size() - (m_p+1); i = i+(m_p-m_r))
            basis_plus.insertKnot(basis_edge.knot(i),m_p-1-m_r);

        m_basis_plus = basis_plus;
        //gsInfo << "Basis plus : " << basis_plus << "\n";

        gsKnotVector<T> kv_minus(0,1,0,m_p+1-1,m_p-1-m_r); // p-1,r //-1 bc p-1
        gsBSplineBasis<> basis_minus(kv_minus);

        // TODO if interface basis are not the same DO NOT DELETE THE FOLLOWING LINES
        /*
        if (basis_1.numElements() <= basis_2.numElements()) //
            for (size_t i = basis_1.degree()+1; i < basis_1.knots().size() - (basis_1.degree()+1); i = i+(basis_1.degree()-m_r))
                basis_minus.insertKnot(basis_1.knot(i),m_p-1-m_r);
        else
            for (size_t i = basis_2.degree()+1; i < basis_2.knots().size() - (basis_2.degree()+1); i = i+(basis_2.degree()-m_r))
                basis_minus.insertKnot(basis_2.knot(i),m_p-1-m_r);
        */
        for (size_t i = m_p+1; i < basis_edge.knots().size() - (m_p+1); i = i+(m_p-m_r))
            basis_minus.insertKnot(basis_edge.knot(i),m_p-1-m_r);

        m_basis_minus = basis_minus;

        // Computing the gluing data
        gsApproxGluingData<T> gluingData(m_mp, m_basis, m_uv, m_isBoundary, m_g1OptionList);
        if (g1OptionList.getInt("gluingData") == gluingData::local)
            gluingData.setLocalGluingData(basis_plus, basis_minus, "edge");
        else if (g1OptionList.getInt("gluingData") == gluingData::global)
            gluingData.setGlobalGluingData();

        m_gD.push_back(gluingData);

        // Basis for the G1 basis
        m_basis_g1 = m_basis.basis(0);
    }

    // Computed the gluing data globally
    void setG1BasisEdge(gsMultiPatch<T> & result);

    void refresh(index_t bfID, std::string typeBf);
    void assemble(index_t i, std::string typeBf); // i == number of bf
    inline void apply(bhVisitor & visitor, index_t i, std::string typeBf); // i == number of bf
    void solve();

    void constructSolution(const gsMatrix<> & solVector, gsMultiPatch<T> & result);

    gsBSpline<> get_alpha() { return m_gD[0].get_alpha_tilde(); }
    gsBSpline<> get_beta() { return m_gD[0].get_beta_tilde(); }

    void plotGluingData(index_t numGd) { m_gD[0].plotGluingData(numGd); }

    index_t get_plus() { return m_basis_plus.size(); }

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


}; // class gsG1BasisEdge

template <class T, class bhVisitor>
void gsApproxG1BasisEdge<T,bhVisitor>::setG1BasisEdge(gsMultiPatch<T> & result)
{
    result.clear();

    index_t n_plus = m_basis_plus.size();
    index_t n_minus = m_basis_minus.size();

    gsMultiPatch<> g1EdgeBasis;

    for (index_t bfID = 3; bfID < n_plus - 3; bfID++) // first 3 and last 3 bf are eliminated
    {
        gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mp.basis(0).component(m_uv)); // u
        index_t degree = temp_basis_first.maxDegree();

        gsMatrix<T> ab = m_basis_plus.support(bfID);

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

        ab_temp = ab;
        for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
        {
            gsMatrix<T> xy = temp_basis_first.support(i);
            if ( (xy(0,0) < ab(0,0)) && (xy(0,1) > ab(0,0)))
                ab_temp(0,0) = xy(0,0);
            if ( (xy(0,0) < ab(0,1)) && (xy(0,1) > ab(0,1)))
                ab_temp(0,1) = xy(0,1);
        }
        ab = ab_temp;

        gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, 1);
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i += temp_basis_first.knots().multiplicityIndex(i))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv.insert(temp_basis_first.knot(i), 1);

        temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mp.basis(0).component(1 - m_uv)); // v
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

        refresh(bfID,"plus");

        assemble(bfID,"plus"); // i == number of bf

        gsSparseSolver<real_t>::BiCGSTABILUT solver;
        gsMatrix<> sol;
        solver.compute(m_system.matrix());
        sol = solver.solve(m_system.rhs());

        constructSolution(sol,g1EdgeBasis);
    }
    for (index_t bfID = 2; bfID < n_minus-2; bfID++)  // first 2 and last 2 bf are eliminated
    {
        gsBSplineBasis<> temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mp.basis(0).component(m_uv)); // u
        index_t degree = temp_basis_first.maxDegree();

        gsMatrix<T> ab = m_basis_minus.support(bfID);

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

        ab_temp = ab;
        for (index_t i = 0; i < temp_basis_first.size(); i++) // only the first two u/v-columns are Dofs (0/1)
        {
            gsMatrix<T> xy = temp_basis_first.support(i);
            if ( (xy(0,0) < ab(0,0)) && (xy(0,1) > ab(0,0)))
                ab_temp(0,0) = xy(0,0);
            if ( (xy(0,0) < ab(0,1)) && (xy(0,1) > ab(0,1)))
                ab_temp(0,1) = xy(0,1);
        }

        gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, 1);
        for (size_t i = degree + 1; i < temp_basis_first.knots().size() - (degree + 1); i += temp_basis_first.knots().multiplicityIndex(i))
            if ((temp_basis_first.knot(i) > ab.at(0)) && (temp_basis_first.knot(i) < ab.at(1)))
                kv.insert(temp_basis_first.knot(i), 1);


        temp_basis_first = dynamic_cast<gsBSplineBasis<> &>(m_mp.basis(0).component(1 - m_uv)); // v
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

        constructSolution(sol,g1EdgeBasis);
    }

    result = g1EdgeBasis;
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
    sz = m_basis.basis(0).size();

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
    gsDofMapper map(m_basis.basis(0));

    gsMatrix<unsigned> act;

    for (index_t i = 2; i < m_basis.basis(0).component(1-m_uv).size(); i++) // only the first two u/v-columns are Dofs (0/1)
    {
        act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 3 : 1, i); // WEST
        map.markBoundary(0, act); // Patch 0
    }
/*
    if (typeBf == "plus")
    {
        gsMatrix<T> ab = m_basis_plus.support(bfID);

        gsMatrix<T> ab_temp = ab;
        for (index_t i = 0; i < m_basis.basis(0).component(m_uv).size(); i++) // only the first two u/v-columns are Dofs (0/1)
        {
            gsMatrix<T> xy = m_basis.basis(0).component(m_uv).support(i);
            if ( (xy(0,0) < ab(0,0)) && (xy(0,1) > ab(0,0)))
                ab_temp(0,0) = xy(0,0);
            if ( (xy(0,0) < ab(0,1)) && (xy(0,1) > ab(0,1)))
                ab_temp(0,1) = xy(0,1);
        }
        ab = ab_temp;

        for (index_t i = 0; i < m_basis.basis(0).component(m_uv).size(); i++) // only the first two u/v-columns are Dofs (0/1)
        {
            gsMatrix<T> xy = m_basis.basis(0).component(m_uv).support(i);
            if ( (xy(0,1) < ab(0,0) + 1e-10) || (xy(0,0) > ab(0,1) - 1e-10) ) // || (xy[0] < ab[0] - 1e-10) || (xy[1] > ab[1] + 1e-10))
            {
                act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, i); // WEST
                map.markBoundary(0, act); // Patch 0
            }
        }

        act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, 0); // WEST
        map.markBoundary(0, act); // Patch 0
        act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3,  m_basis.basis(0).component(m_uv).size() -1); // WEST
        map.markBoundary(0, act); // Patch 0
    }
    else if (typeBf == "minus")
    {
        gsMatrix<T> ab = m_basis_minus.support(bfID);

        gsMatrix<T> ab_temp = ab;
        for (index_t i = 0; i < m_basis.basis(0).component(m_uv).size(); i++) // only the first two u/v-columns are Dofs (0/1)
        {
            gsMatrix<T> xy = m_basis.basis(0).component(m_uv).support(i);
            if ( (xy(0,0) < ab(0,0)) && (xy(0,1) > ab(0,0)))
                ab_temp(0,0) = xy(0,0);
            if ( (xy(0,0) < ab(0,1)) && (xy(0,1) > ab(0,1)))
                ab_temp(0,1) = xy(0,1);
        }
        ab = ab_temp;

        for (index_t i = 0; i < m_basis.basis(0).component(m_uv).size(); i++) // only the first two u/v-columns are Dofs (0/1)
        {
            gsMatrix<T> xy = m_basis.basis(0).component(m_uv).support(i);
            if ( (xy(0,1) < ab(0,0) + 1e-10) || (xy(0,0) > ab(0,1) - 1e-10) ) //|| (xy[0] < ab[0] - 1e-10) || (xy[1] > ab[1] + 1e-10))
            {
                act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, i); // WEST
                map.markBoundary(0, act); // Patch 0
            }
        }

        act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3, 0); // WEST
        map.markBoundary(0, act); // Patch 0
        act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 1 : 3,  m_basis.basis(0).component(m_uv).size() -1); // WEST
        map.markBoundary(0, act); // Patch 0

        if (m_isBoundary)
        {
            act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 3 : 1, 0); // WEST
            map.markBoundary(0, act); // Patch 0
        }

    }
*/
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
            visitor_.evaluate(bf_index, typeBf, basis_g1, basis_geo, basis_plus, basis_minus, patch, quNodes, m_uv, m_gD[0], m_isBoundary, m_g1OptionList);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(0, m_ddof, m_system); // omp_locks inside // patchIndex == 0
        }
    }//omp parallel
} // apply

} // namespace gismo