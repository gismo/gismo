/** @file gsG1BasisEdge.h

    @brief Provides assembler for a G1 Basis for multiPatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

#include <gsG1Basis/gsApproxGluingData.h>
#include <gsG1Basis/gsG1ASGluingData.h>
#include <gsG1Basis/gsVisitorG1BasisEdge.h>
# include <gsAssembler/gsAssembler.h>

# include <gsG1Basis/gsG1OptionList.h>

namespace gismo
{
template<class T, class bhVisitor = gsVisitorG1BasisEdge<T>>
class gsG1BasisEdge : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsG1BasisEdge(gsMultiPatch<> geo, // single patch
                 gsMultiBasis<> basis, // single basis
                 index_t uv, // !!! 0 == u; 1 == v !!!
                 bool isBoundary,
                 gsG1OptionList & g1OptionList)
        : m_geo(geo), m_basis(basis), m_uv(uv), m_isBoundary(isBoundary), m_g1OptionList(g1OptionList)
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
        n_plus = m_basis_plus.size();
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
        n_minus = m_basis_minus.size();

        // Computing the gluing data
        gsApproxGluingData<T> gluingData(m_geo, m_basis, m_uv, m_isBoundary, m_g1OptionList);
        if (g1OptionList.getInt("gluingData") == gluingData::local)
            gluingData.setLocalGluingData(basis_plus, basis_minus);
        else if (g1OptionList.getInt("gluingData") == gluingData::l2projection)
            gluingData.setGlobalGluingData();

        m_gD.push_back(gluingData);

        // Basis for the G1 basis
        m_basis_g1 = m_basis.basis(0);

        refresh();
        assemble();
        solve();
    }


    void refresh();
    void assemble();
    inline void apply(bhVisitor & visitor, int patchIndex, boxSide side = boundary::none);
    void solve();

    void constructSolution(gsMultiPatch<T> & result);

    gsBSpline<> get_alpha() { return m_gD[0].get_alpha_tilde(); }
    gsBSpline<> get_beta() { return m_gD[0].get_beta_tilde(); }

protected:

    // Input
    gsMultiPatch<T> m_geo;
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

    // Size of the basis
    index_t n_plus, n_minus;

    // System
    std::vector<gsSparseSystem<T> > m_f_0, m_f_1;

    // For Dirichlet boundary
    using Base::m_ddof;

    std::vector<gsMatrix<>> solVec_t, solVec_b;
    gsMultiPatch<T> g1Basis;

    index_t reduced_minus, reduced_plus;

}; // class gsG1BasisEdge


template <class T, class bhVisitor>
void gsG1BasisEdge<T,bhVisitor>::constructSolution(gsMultiPatch<T> & result)
{

    result.clear();

    // Dim is the same for all basis functions
    const index_t dim = ( 0!=solVec_t.at(0).cols() ? solVec_t.at(0).cols() :  m_ddof[0].cols() );

    gsMatrix<T> coeffs;
    for (index_t p = 0; p < m_basis_plus.size(); ++p)
    {
        const gsDofMapper & mapper = m_f_0.at(p).colMapper(0); // unknown = 0

        // Reconstruct solution coefficients on patch p
        index_t sz;
        sz = m_basis.basis(0).size();

        coeffs.resize(sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if (mapper.is_free(i, 0)) // DoF value is in the solVector // 0 = unitPatch
            {
                coeffs.row(i) = solVec_t.at(p).row(mapper.index(i, 0));
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                coeffs.row(i) = m_ddof[0].row( mapper.bindex(i, 0) ).head(dim); // = 0
            }
        }
        result.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));
    }

    for (index_t p = 0; p < m_basis_minus.size(); ++p)
    {
        const gsDofMapper & mapper = m_f_1.at(p).colMapper(0); // unknown = 0

        // Reconstruct solution coefficients on patch p
        index_t sz;
        sz = m_basis.basis(0).size();

        coeffs.resize(sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if (mapper.is_free(i, 0)) // DoF value is in the solVector // 0 = unitPatch
            {
                //gsInfo << "mapper index: " << mapper.index(i, p) << "\n";
                coeffs.row(i) = solVec_b.at(p).row(mapper.index(i, 0));
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                coeffs.row(i) = m_ddof[0].row( mapper.bindex(i, 0) ).head(dim); // = 0
            }
        }

        result.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));

    }


    g1Basis = result;

    //g1Condition();
}

template <class T, class bhVisitor>
void gsG1BasisEdge<T,bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map(m_basis.basis(0));

    gsMatrix<unsigned> act;

    for (index_t i = 2; i < m_basis.basis(0).component(1-m_uv).size(); i++) // only the first two u/v-columns are Dofs (0/1)
    {
        act = m_basis.basis(0).boundaryOffset(m_uv == 0 ? 3 : 1, i); // WEST
        map.markBoundary(0, act); // Patch 0
    }

    map.finalize();
    //gsInfo << "map : " << map.asVector() << "\n";
    //map.print();

    // 2. Create the sparse system
    gsSparseSystem<T> m_system = gsSparseSystem<T>(map);
    for (index_t i = 0; i < m_basis_plus.size() ; i++)
        m_f_0.push_back(m_system);
    for (index_t i = 0; i < m_basis_minus.size(); i++)
        m_f_1.push_back(m_system);


} // refresh()

template <class T, class bhVisitor>
void gsG1BasisEdge<T,bhVisitor>::assemble()
{
    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis[0],2,1,0.333333);
    for (index_t i = 0; i < m_basis_plus.size(); i++)
        m_f_0.at(i).reserve(nz, 1);

    for (index_t i = 0; i < m_basis_minus.size(); i++)
        m_f_1.at(i).reserve(nz, 1);

    if(m_ddof.size()==0)
        m_ddof.resize(2); // 0,1

    const gsDofMapper & map_0 = m_f_0.at(0).colMapper(0); // Map same for every 0
    const gsDofMapper & map_1 = m_f_1.at(0).colMapper(0); // Map same for every 1

    m_ddof[0].setZero(map_0.boundarySize(), 1 ); // plus
    m_ddof[1].setZero(map_1.boundarySize(), 1 ); // minus

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor,0); // patch 0

    for (index_t i = 0; i < m_basis_plus.size(); i++)
        m_f_0.at(i).matrix().makeCompressed();
    for (index_t i = 0; i < m_basis_minus.size(); i++)
        m_f_1.at(i).matrix().makeCompressed();

} // assemble()

template <class T, class bhVisitor>
void gsG1BasisEdge<T,bhVisitor>::apply(bhVisitor & visitor, int patchIndex, boxSide side)
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
        visitor_.initialize(basis_g1, basis_plus, basis_minus, quRule);

        const gsGeometry<T> & patch = m_geo.patch(0);

        // Initialize domain element iterator
        typename gsBasis<T>::domainIter domIt = basis_g1.makeDomainIterator(boundary::none);

#ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
        for (; domIt->good(); domIt->next() )
#endif
        {

            // Map the Quadrature rule to the element
            quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

            // Perform required evaluations on the quadrature nodes
            visitor_.evaluate(basis_g1, basis_geo, basis_plus, basis_minus, patch, quNodes, m_uv, m_gD[0], m_isBoundary, m_g1OptionList);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex, m_ddof, m_f_0, m_f_1); // omp_locks inside // patchIndex == 0
        }
    }//omp parallel
} // apply

template <class T, class bhVisitor>
void gsG1BasisEdge<T,bhVisitor>::solve()
{
    gsSparseSolver<real_t>::CGDiagonal solver;

    for (index_t i = 0; i < m_basis_plus.size(); i++) // Tilde
    {
        solver.compute(m_f_0.at(i).matrix());
        solVec_t.push_back(solver.solve(m_f_0.at(i).rhs()));
    }
    for (index_t i = 0; i < m_basis_minus.size(); i++)
    {
        solver.compute(m_f_1.at(i).matrix());
        solVec_b.push_back(solver.solve(m_f_1.at(i).rhs()));
    }
} // solve

} // namespace gismo