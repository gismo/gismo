/** @file gsG1BasisVertex.h

    @brief Provides assembler for a G1 Basis for multiPatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

#include <gsG1Basis/gsGluingData.h>
#include <gsG1Basis/gsVisitorG1BasisVertex.h>


namespace gismo
{
template<class T, class bhVisitor = gsVisitorG1BasisVertex<T>>
class gsG1BasisVertex : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsG1BasisVertex(gsMultiPatch<> geo, // Single Patch
                  gsMultiBasis<> basis, // Single Basis
                  std::vector<bool> isBoundary,
                  gsOptionList & optionList)
        : m_geo(geo), m_basis(basis), m_isBoundary(isBoundary), m_optionList(optionList)
    {
        for (index_t uv = 0; uv < 2; uv++) // For the TWO interface
        {
            // Computing the gluing data
            gsGluingData<T> gluingData(m_geo,m_basis,uv,m_isBoundary[uv],m_optionList);
            m_gD.push_back(gluingData);

            // Computing the G1 - basis function at the edge
            // Spaces for computing the g1 basis
            index_t m_r = m_optionList.getInt("regularity"); // TODO CHANGE IF DIFFERENT REGULARITY IS NECESSARY

            gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(uv)); // 0 -> u, 1 -> v
            index_t m_p = basis_edge.maxDegree(); // Minimum degree at the interface // TODO if interface basis are not the same

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
        }

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


    void plotG1Basis(gsMultiPatch<T> & basisG1_L, gsMultiPatch<T> & basisG1_R, gsMultiPatch<T> & mp, std::string baseName)
    {

        const std::string baseName1(baseName + "_0");
        gsParaviewCollection collection1(baseName1);

        const std::string baseName2(baseName + "_1");
        gsParaviewCollection collection2(baseName2);

        std::string fileName, fileName2;
        for (unsigned i = 0; i < basisG1_L.nPatches(); i++)
        {

            fileName = baseName1 + "_" + util::to_string(i);
            gsField<> temp_field_L(mp.patch(0),basisG1_L.patch(i));
            gsWriteParaview(temp_field_L,fileName,5000);
            collection1.addTimestep(fileName,i,"0.vts");

        }
        for (unsigned i = 0; i < basisG1_R.nPatches(); i++)
        {

            fileName2 = baseName2 + "_" + util::to_string(i);
            gsField<> temp_field_R(mp.patch(1),basisG1_R.patch(i));
            gsWriteParaview(temp_field_R,fileName2,5000);
            collection2.addTimestep(fileName2,i,"0.vts");

        }
        collection1.save();
        collection2.save();
    }

    void plotG1BasisBoundary(gsMultiPatch<T> & basisG1_boundary, gsMultiPatch<T> mp, std::string baseName)
    {
        gsParaviewCollection collection1(baseName);
        std::string fileName;
        for (unsigned i = 0; i < basisG1_boundary.nPatches(); i++)
        {

            fileName = baseName + "_" + util::to_string(i);
            gsField<> temp_field_L(mp.patch(0),basisG1_boundary.patch(i));
            gsWriteParaview(temp_field_L,fileName,5000);
            collection1.addTimestep(fileName,i,"0.vts");

        }
        collection1.save();
    }

    void g1Condition()
    {
        gsMatrix<> points(1,1000);
        points.setRandom();
        points = points.array().abs();

        gsMatrix<> points2d_L(2, 1000);
        gsMatrix<> points2d_R(2, 1000);

        points2d_L.setZero();
        points2d_R.setZero();
        points2d_L.row(1) = points; // v
        points2d_R.row(0) = points; // u

        real_t g1Error = 0;
        /*
        for (size_t i = 0; i < g1Basis.nPatches(); i++)
        {
            gsMatrix<> temp;
            temp = m_gD.get_alpha_tilde_1().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).topRows(1))
                + m_gD.get_alpha_tilde_0().eval(points).cwiseProduct(g1Basis_R.patch(i).deriv(points2d_R).bottomRows(1))
                + m_gD.get_beta_bar().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).bottomRows(1));

            if (temp.array().abs().maxCoeff() > g1Error)
                g1Error = temp.array().abs().maxCoeff();
        }
        */
        gsInfo << "Conditiontest G1 continuity: \n" << g1Error << "\n\n";


        //gsInfo << "\nConditiontest G1 continuity: \n" << g1Basis_L.patch(0).coefs() << "\n\n";

    }


protected:

    // Input
    gsMultiPatch<T> m_geo;
    gsMultiBasis<T> m_basis;
    std::vector<bool> m_isBoundary;
    gsOptionList m_optionList;

    // Gluing data
    std::vector<gsGluingData<T>> m_gD;

    // Basis for getting the G1 Basis
    std::vector<gsBSplineBasis<>> m_basis_plus;
    std::vector<gsBSplineBasis<>> m_basis_minus;

    // Basis for the G1 Basis
    gsMultiBasis<T> m_basis_g1;

    // System
    std::vector<gsSparseSystem<T> > m_f;

    // For Dirichlet boundary
    using Base::m_ddof;

    std::vector<gsMatrix<>> solVec;
    gsMultiPatch<T> g1Basis;

}; // class gsG1BasisEdge


template <class T, class bhVisitor>
void gsG1BasisVertex<T,bhVisitor>::constructSolution(gsMultiPatch<T> & result)
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
        sz = m_basis.basis(0).size();

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
                coeffs.row(i) = m_ddof[0].row( mapper.bindex(i, 0) ).head(dim); // = 0
            }
        }
        result.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));
    }

    g1Basis = result;

    //g1Condition();
}

template <class T, class bhVisitor>
void gsG1BasisVertex<T,bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map(m_basis.basis(0));

    gsMatrix<unsigned> act(m_basis.basis(0).size()-6,1);
    gsVector<unsigned> vec;
    index_t dimU = m_basis.basis(0).component(0).size();
    vec.setLinSpaced(m_basis.basis(0).size(),0,m_basis.basis(0).size());
    vec.block(2*dimU,0,vec.size()-2*dimU-1,1) = vec.block(2*dimU +1,0,vec.size()-2*dimU-1,1); // 2*dim U
    vec.block(dimU,0,vec.size()-dimU-2,1) = vec.block(dimU +2,0,vec.size()-dimU-2,1); // dim U, dimU+1
    vec.block(0,0,vec.size()-3,1) = vec.block(3,0,vec.size()-3,1); // 0,1,2
    act = vec.block(0,0,vec.size()-6,1);
    //map.markBoundary(0,act); // TODO TODO TODO TODO is wrong, need bigger support!!!!!!
    map.finalize();
    //gsInfo << "map : " << map.asVector() << "\n";
    //map.print();

    // 2. Create the sparse system
    gsSparseSystem<T> m_system = gsSparseSystem<T>(map);
    for (index_t i = 0; i < 6; i++)
        m_f.push_back(m_system);

} // refresh()

template <class T, class bhVisitor>
void gsG1BasisVertex<T,bhVisitor>::assemble()
{
    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis[0],2,1,0.333333);
    for (unsigned i = 0; i < m_f.size(); i++)
        m_f.at(i).reserve(nz, 1);


    if(m_ddof.size()==0)
        m_ddof.resize(1); // 0,1

    const gsDofMapper & map = m_f.at(0).colMapper(0); // Map same for every

    m_ddof[0].setZero(map.boundarySize(), 1 ); // plus

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor,0); // patch 0

    for (unsigned i = 0; i < m_f.size(); i++)
        m_f.at(i).matrix().makeCompressed();

} // assemble()

template <class T, class bhVisitor>
void gsG1BasisVertex<T,bhVisitor>::apply(bhVisitor & visitor, int patchIndex, boxSide side)
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
#pragma omp critical(evaluate)
            // Perform required evaluations on the quadrature nodes
            visitor_.evaluate(basis_g1, basis_geo, m_basis_plus, m_basis_minus, patch, quNodes, m_gD, m_isBoundary, m_optionList);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex, m_ddof, m_f); // omp_locks inside // patchIndex == 0
        }
    }//omp parallel
} // apply

template <class T, class bhVisitor>
void gsG1BasisVertex<T,bhVisitor>::solve()
{
    gsSparseSolver<real_t>::CGDiagonal solver;

    for (index_t i = 0; i < 6; i++) // Tilde
    {
        solver.compute(m_f.at(i).matrix());
        solVec.push_back(solver.solve(m_f.at(i).rhs()));
    }
} // solve

} // namespace gismo