/** @file gsG1Basis.h

    @brief Provides assembler for a G1 Basis.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

#include <gsG1Basis/gsGluingData.h>
#include <gsG1Basis/gsVisitorG1Basis.h>

namespace gismo
{
template <class T, class bhVisitor = gsVisitorG1Basis<T> >
class gsG1Basis : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsG1Basis(gsMultiPatch<T> const & mp,
              gsMultiBasis<T> & mb,
              index_t const  k,
              index_t const  r,
              index_t const  p_t,
              index_t const  r_t,
              bool const dir,
              bool const loc,
              bool const pl)
              : m_mp(mp), m_mb(mb), m_k(k), p_tilde(p_t), r_tilde(r_t), m_r(r), direct(dir), local(loc), plot(pl)
    {
        const boundaryInterface & iFace = *m_mp.iBegin();// assume one single interface
        gsGluingData<T> gluingData(m_mp,m_mb,iFace,m_k,m_r,p_tilde,r_tilde,direct,plot);

        m_gD = gluingData;

        // Spaces for g1 basis
        // first,last,interior,mult_ends,mult_interior
        gsKnotVector<T> kv_tilde(0,1,m_k,m_mb.basis(0).component(1).maxDegree()+1,m_mb.basis(0).component(1).maxDegree()-1-m_r); // p,r+1 //-1 bc r+1
        gsBSplineBasis<> basis_tilde(kv_tilde);
        m_basis_tilde = basis_tilde;

        n_tilde = m_basis_tilde.size();

        gsKnotVector<T> kv_bar(0,1,m_k,m_mb.basis(0).component(1).maxDegree()+1-1,m_mb.basis(0).component(1).maxDegree()-1-m_r); // p-1,r //-1 bc p-1
        gsBSplineBasis<> basis_bar(kv_bar);
        m_basis_bar = basis_bar;

        n_bar = m_basis_bar.size();

        gsKnotVector<T> kv_geo(0,1,m_k,m_mb.basis(0).component(0).maxDegree()+1,m_mb.basis(0).component(0).maxDegree()-m_r); // p,r
        gsBSplineBasis<> basis_geo(kv_geo);
        m_basis_geo = basis_geo;

        m_basis_g1 = m_mb; // basis for the g1 basis functions

        refresh();
    }

    void refresh();
    void assemble();
    inline void apply(bhVisitor & visitor, int patchIndex = 0, boxSide side = boundary::none);

    void solve();
    void constructSolution(gsMultiPatch<T> & result_L, gsMultiPatch<T> & result_R);

    index_t get_n_tilde() { return n_tilde; }
    index_t get_n_bar() { return n_bar; }

    void g1Condition()
    {
        gsMatrix<> points(1,1000);
        points.setRandom();
        points = points.array().abs();

        gsMatrix<> points2d_R(2,1000);
        points2d_R.setZero();
        points2d_R.row(1) = points;
        gsMatrix<> points2d_L(2,1000);
        points2d_L.setOnes();
        points2d_L.row(1) = points;

        if (direct)
        {
            gsMatrix<T> a_L, b_L, a_R, b_R;
            m_gD.eval_into_alpha_L(points2d_L,a_L);
            m_gD.eval_into_alpha_R(points2d_R,a_R);
            m_gD.eval_into_beta_L(points2d_L,b_L);
            m_gD.eval_into_beta_R(points2d_R,b_R);

            real_t g1Error = 0;
            for (size_t i = 0; i < g1Basis_L.nPatches(); i++)
            {
                gsMatrix<> temp = a_R.cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).topRows(1))
                    + a_L.cwiseProduct(g1Basis_R.patch(i).deriv(points2d_R).topRows(1))
                    + m_gD.get_beta_bar().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).bottomRows(1));

                if (temp.array().abs().maxCoeff() > g1Error)
                    g1Error = temp.array().abs().maxCoeff();
            }

            //gsInfo << "\nConditiontest G1 continuity: \n" << g1Error << "\n\n";
        }
        else
        {
            real_t g1Error = 0;
            for (size_t i = 0; i < g1Basis_L.nPatches(); i++)
            {
                gsMatrix<> temp = m_gD.get_alpha_tilde_R().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).topRows(1))
                    + m_gD.get_alpha_tilde_L().eval(points).cwiseProduct(g1Basis_R.patch(i).deriv(points2d_R).topRows(1))
                    + m_gD.get_beta_bar().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).bottomRows(1));

                if (temp.array().abs().maxCoeff() > g1Error)
                    g1Error = temp.array().abs().maxCoeff();
            }

            gsInfo << "\nConditiontest G1 continuity: \n" << g1Error << "\n\n";
        }

        //gsInfo << "\nConditiontest G1 continuity: \n" << g1Basis_L.patch(0).coefs() << "\n\n";

    }

    void plotG1Basis(gsMultiPatch<T> & basisG1_L, gsMultiPatch<T> & basisG1_R)
    {

        const std::string baseName1("gDBasis_L");
        gsParaviewCollection collection1(baseName1);

        const std::string baseName2("gDBasis_R");
        gsParaviewCollection collection2(baseName2);

        std::string fileName, fileName2;
        for (unsigned i = 0; i < basisG1_L.nPatches(); i++)
        {

            fileName = baseName1 + "_" + util::to_string(i);

            gsField<> temp_field_L(m_mp.patch(0),basisG1_L.patch(i));
            gsWriteParaview(temp_field_L,fileName,5000);
            collection1.addTimestep(fileName,i,"0.vts");

        }
        for (unsigned i = 0; i < basisG1_R.nPatches(); i++)
        {

            fileName2 = baseName2 + "_" + util::to_string(i);

            gsField<> temp_field_R(m_mp.patch(1),basisG1_R.patch(i));
            gsWriteParaview(temp_field_R,fileName2,5000);
            collection2.addTimestep(fileName2,i,"0.vts");

        }
        collection1.save();
        collection2.save();
    }

protected:

    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;

    index_t m_k, p_tilde, r_tilde, m_r;

    bool direct;
    bool local;
    bool plot;

    // Gluing data
    gsGluingData<T> m_gD;

    // Spaces for the g1 basis
    gsMultiBasis<> m_basis_g1;
    gsBSplineBasis<> m_basis_geo;
    gsBSplineBasis<> m_basis_tilde;
    gsBSplineBasis<> m_basis_bar;

    std::vector< gsSparseSystem<T> > m_phi_i, m_phi_j;

    using Base::m_ddof;

    std::vector<gsMatrix<>> solVec_t, solVec_b;
    gsMultiPatch<T> g1Basis_L, g1Basis_R;

    index_t n_tilde, n_bar;
}; // class gsG1Basis


template <class T, class bhVisitor>
void gsG1Basis<T,bhVisitor>::constructSolution(gsMultiPatch<T> & result_L, gsMultiPatch<T> & result_R)
{

    result_L.clear();
    result_R.clear();

    // Dim is the same for all basis functions
    const index_t dim = ( 0!=solVec_t.at(0).cols() ? solVec_t.at(0).cols() :  m_ddof[0].cols() );

    gsMatrix<T> coeffs;
    for (size_t patchNum = 0; patchNum < m_mp.nPatches(); patchNum++)
    {
        for (index_t p = 0; p < m_basis_tilde.size(); ++p)
        {
            const gsDofMapper & mapper = m_phi_i.at(p).colMapper(0); // unknown = 0

            // Reconstruct solution coefficients on patch p
            const int sz = m_basis_g1.basis(0).size();
            coeffs.resize(sz, dim);

            for (index_t i = 0; i < sz; ++i)
            {
                if (mapper.is_free(i, patchNum)) // DoF value is in the solVector // 0 = unitPatch
                {
                    coeffs.row(i) = solVec_t.at(p).row(mapper.index(i, patchNum));
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                    coeffs.row(i) = m_ddof[0].row( mapper.bindex(i, patchNum) ).head(dim); // = 0
                }
            }
            if (patchNum == 0)
                result_L.addPatch(m_basis_g1.basis(patchNum).makeGeometry(give(coeffs)));
            if (patchNum == 1)
                result_R.addPatch(m_basis_g1.basis(patchNum).makeGeometry(give(coeffs)));
        }

        for (index_t p = 0; p < m_basis_bar.size(); ++p)
        {
            const gsDofMapper & mapper = m_phi_j.at(p).colMapper(0); // unknown = 0

            // Reconstruct solution coefficients on patch p
            const int sz = m_basis_g1.basis(0).size();
            coeffs.resize(sz, dim);

            for (index_t i = 0; i < sz; ++i)
            {
                if (mapper.is_free(i, patchNum)) // DoF value is in the solVector // 0 = unitPatch
                {
                    //gsInfo << "mapper index: " << mapper.index(i, p) << "\n";
                    coeffs.row(i) = solVec_b.at(p).row(mapper.index(i, patchNum));
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                    coeffs.row(i) = m_ddof[0].row( mapper.bindex(i, patchNum) ).head(dim); // = 0
                }
            }
            if (patchNum == 0)
                result_L.addPatch(m_basis_g1.basis(patchNum).makeGeometry(give(coeffs)));
            if (patchNum == 1)
                result_R.addPatch(m_basis_g1.basis(patchNum).makeGeometry(give(coeffs)));
        }
    }

    g1Basis_L = result_L;
    g1Basis_R = result_R;

    //g1Condition();
}


template <class T, class bhVisitor>
void gsG1Basis<T,bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map(m_basis_g1);


    const boundaryInterface & iFace = *m_mp.iBegin();// assume one single interface
    gsBasis<> & B_R = m_mb.basis(iFace.first().patch); // Patch Right
    gsBasis<> & B_L = m_mb.basis(iFace.second().patch); // Patch Left

    gsMatrix<unsigned> actR, actL;

    for (index_t i = 2; i < B_L.component(0).size(); i++) // B_L == B_R, i = 0,1 interface
    {
        actR = B_R.boundaryOffset(iFace.first() .side(), i); //first adj. face
        actL = B_L.boundaryOffset(iFace.second().side(), i); //second adj. face

        map.markBoundary(iFace.first().patch, actR);
        map.markBoundary(iFace.second().patch, actL);
    }
    map.finalize();
    //map.print();

    // 2. Create the sparse system
    gsSparseSystem<T> m_system = gsSparseSystem<T>(map);
    for (index_t i = 0; i < m_basis_tilde.size(); i++)
        m_phi_i.push_back(m_system);
    for (index_t i = 0; i < m_basis_bar.size(); i++)
        m_phi_j.push_back(m_system);


} // refresh()

template <class T, class bhVisitor>
void gsG1Basis<T,bhVisitor>::assemble()
{
    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_mb[0],2,1,0.333333);
    for (index_t i = 0; i < m_basis_tilde.size(); i++)
        m_phi_i.at(i).reserve(nz, 1);

    for (index_t i = 0; i < m_basis_bar.size(); i++)
        m_phi_j.at(i).reserve(nz, 1);

    if(m_ddof.size()==0)
        m_ddof.resize(1); // tilde = bar

    const gsDofMapper & map = m_phi_i.at(0).colMapper(0); // Map same for tilde and bar

    m_ddof[0].setZero(map.boundarySize(), 1 ); // tilde

    // Assemble volume integrals
    bhVisitor visitor;
    for (size_t np = 0; np < m_mp.nPatches(); ++np)
        apply(visitor,np);

    for (index_t i = 0; i < m_basis_tilde.size(); i++)
        m_phi_i.at(i).matrix().makeCompressed();
    for (index_t i = 0; i < m_basis_bar.size(); i++)
        m_phi_j.at(i).matrix().makeCompressed();

} // assemble()

template <class T, class bhVisitor>
void gsG1Basis<T,bhVisitor>::apply(bhVisitor & visitor, int patchIndex, boxSide side)
{
    gsQuadRule<T> quRule ; // Quadrature rule
    gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
    gsVector<T> quWeights; // Temp variable for mapped weights

    gsBasis<T> & basis_g1 = m_basis_g1.basis(patchIndex); // basis for construction

    // Same for all patches
    gsBasis<T> & basis_geo = m_basis_geo;
    gsBasis<T> & basis_tilde = m_basis_tilde;
    gsBasis<T> & basis_bar = m_basis_bar;

    // Initialize reference quadrature rule and visitor data
    visitor.initialize(basis_g1, basis_geo, basis_tilde, basis_bar, quRule);

    const gsGeometry<T> & patch = m_mp.patch(patchIndex);

    // Initialize domain element iterator
    typename gsBasis<T>::domainIter domIt = basis_g1.makeDomainIterator(boundary::none);
    for (; domIt->good(); domIt->next() )
    {

        // Map the Quadrature rule to the element
        quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

        // Perform required evaluations on the quadrature nodes
        visitor.evaluate(basis_g1, basis_geo, basis_tilde, basis_bar, patch, quNodes, m_gD, direct, local);

        // Assemble on element
        visitor.assemble(*domIt, quWeights);

        // Push to global matrix and right-hand side vector
        visitor.localToGlobal(patchIndex, m_ddof, m_phi_i, m_phi_j); // omp_locks inside
    }
} // apply

template <class T, class bhVisitor>
void gsG1Basis<T,bhVisitor>::solve()
{
    gsSparseSolver<real_t>::CGDiagonal solver;

    for (index_t i = 0; i < m_basis_tilde.size(); i++) // Tilde
    {
        solver.compute(m_phi_i.at(i).matrix());
        solVec_t.push_back(solver.solve(m_phi_i.at(i).rhs()));
    }
    for (index_t i = 0; i < m_basis_bar.size(); i++)
    {
        solver.compute(m_phi_j.at(i).matrix());
        solVec_b.push_back(solver.solve(m_phi_j.at(i).rhs()));
    }
} // solve

} // namespace gismo