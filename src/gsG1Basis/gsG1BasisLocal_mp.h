/** @file gsG1Basis_mp.h

    @brief Provides assembler for a G1 Basis for multiPatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

#include <gsG1Basis/gsGluingData_mp.h>
#include <gsG1Basis/gsVisitorG1BasisLocal_mp.h>

namespace gismo
{
template<class T, class bhVisitor = gsVisitorG1BasisLocal_mp<T> >
class gsG1BasisLocal_mp : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsG1BasisLocal_mp(gsMultiPatch<T> const & mp,
                 gsMultiBasis<T> & mb,
                 boundaryInterface iFace,
                 index_t const p_t,
                 index_t const r_t,
                 bool const dir,
                 bool const loc,
                 bool const pl)
        : m_mp(mp), m_mb(mb), m_iFace(iFace), p_tilde(p_t), r_tilde(r_t), direct(dir), local(loc), plot(pl)
    {

        index_t max_degree; // Assume that the degree is the same for BOTH patches
        index_t m_k_geo;
        if (m_iFace.first().side().index() >= 3) // north or south of BOTH patches
        {
            gsBasis<T> & basis_first = m_mb.basis(m_iFace.first().patch).component(0);
            gsBasis<T> & basis_second = m_mb.basis(m_iFace.second().patch).component(0);
            if (basis_first.numElements() >= basis_second.numElements())
                m_k = basis_second.numElements()-1;
            else if (basis_first.numElements() < basis_second.numElements())
                m_k = basis_first.numElements()-1;

            max_degree = basis_first.maxDegree();
        }
        else if (m_iFace.first().side().index() <= 2) // west or east of BOTH patches
        {
            gsBasis<T> & basis_first = m_mb.basis(m_iFace.first().patch).component(1);
            gsBasis<T> & basis_second = m_mb.basis(m_iFace.second().patch).component(1);
            if (basis_first.numElements() >= basis_second.numElements())
                m_k = basis_second.numElements()-1;
            else if (basis_first.numElements() < basis_second.numElements())
                m_k = basis_first.numElements()-1;

            max_degree = basis_first.maxDegree();
        }
        m_r = 1; // ASSUME that the regularity is the same on each patch

        //gsInfo << " m_K : " << m_k << "\n";

        gsGluingData_mp<T> gluingData(m_mp,m_mb,m_iFace,m_k,m_r,p_tilde,r_tilde,direct,plot);
        m_gD = gluingData;

        // Spaces for g1 basis
        // first,last,interior,mult_ends,mult_interior
        gsKnotVector<T> kv_tilde(0,1,0,max_degree+1,max_degree-1-m_r); // p,r+1 //-1 bc r+1
        gsTensorBSplineBasis<2,real_t> & temp_Tensorbasis_L = dynamic_cast<gsTensorBSplineBasis<2,real_t> &>(m_mb.basis(0));
        gsBSplineBasis<> basis_temp_L = dynamic_cast<gsBSplineBasis<> &>(temp_Tensorbasis_L.component(1));
        gsTensorBSplineBasis<2,real_t> & temp_Tensorbasis_R = dynamic_cast<gsTensorBSplineBasis<2,real_t> &>(m_mb.basis(1));
        gsBSplineBasis<> basis_temp_R = dynamic_cast<gsBSplineBasis<> &>(temp_Tensorbasis_R.component(1));
        gsBSplineBasis<> BASIS_tilde_R(kv_tilde);
        gsBSplineBasis<> BASIS_tilde_L(kv_tilde);
        gsBSplineBasis<> basis_tilde(kv_tilde);

        typename gsBasis<T>::domainIter it_L = basis_temp_L.makeDomainIterator();
        typename gsBasis<T>::domainIter it_R = basis_temp_R.makeDomainIterator();

        it_L->good();
        it_R->good();
        for (index_t i = max_degree+1; i < basis_temp_L.knots().size() - (max_degree+1); i = i+(max_degree-m_r))
        {
            //if (it_L->getMaxCellLength() >= it_R->getMaxCellLength())
                basis_tilde.insertKnot(basis_temp_L.knot(i),max_degree-1-m_r);
            //else if (it_R->getMaxCellLength() > it_L->getMaxCellLength())
                //basis_tilde.insertKnot(basis_temp_R.knot(i),max_degree-1-m_r);

            BASIS_tilde_L.insertKnot(basis_temp_L.knot(i),max_degree-1-m_r);
            BASIS_tilde_R.insertKnot(basis_temp_R.knot(i),max_degree-1-m_r);

            it_L->next();
            it_R->next();
        }
        m_basis_tilde = basis_tilde;

        //gsInfo << m_basis_tilde << "\n";

        n_tilde = m_basis_tilde.size();

        gsKnotVector<T> kv_bar(0,1,0,max_degree+1-1,max_degree-1-m_r); // p-1,r //-1 bc p-1
        gsBSplineBasis<> basis_bar(kv_bar);

        it_L->good();
        it_R->good();
        for (index_t i = max_degree+1; i < basis_temp_L.knots().size() - (max_degree+1); i = i+(max_degree-m_r))
        {
            //if (it_L->getMaxCellLength() >= it_R->getMaxCellLength())
                basis_bar.insertKnot(basis_temp_L.knot(i),max_degree-1-m_r);
            //else if (it_R->getMaxCellLength() > it_L->getMaxCellLength())
                //basis_bar.insertKnot(basis_temp_R.knot(i),max_degree-1-m_r);

            //basis_bar.insertKnot(basis_temp_L.knot(i),max_degree-1-m_r);
            //basis_bar.insertKnot(basis_temp_R.knot(i),max_degree-1-m_r);

            it_L->next();
            it_R->next();
        }
        m_basis_bar = basis_bar;

        n_bar = m_basis_bar.size();

        if (m_iFace.first().side().index() >= 3) // north or south of BOTH patches
        {
            gsBasis<T> & basis_first = m_mb.basis(m_iFace.first().patch).component(1); // right
            gsBasis<T> & basis_second = m_mb.basis(m_iFace.second().patch).component(1); // left

            m_k = basis_second.numElements() -1;
            max_degree = basis_second.maxDegree();
            gsKnotVector<T> kv_geo(0,1,m_k,max_degree+1,max_degree-m_r); // p,r
            gsBSplineBasis<> basis_geo(kv_geo);
            m_basis_geo.push_back(basis_geo);

            m_k = basis_first.numElements() -1;
            max_degree = basis_first.maxDegree();
            gsKnotVector<T> kv_geo2(0,1,m_k,max_degree+1,max_degree-m_r); // p,r
            gsBSplineBasis<> basis_geo2(kv_geo2);
            m_basis_geo.push_back(basis_geo2);
        }
        else if (m_iFace.first().side().index() <= 2) // west or east of BOTH patches
        {
            gsBasis<T> & basis_first = m_mb.basis(m_iFace.first().patch).component(0);
            gsBasis<T> & basis_second = m_mb.basis(m_iFace.second().patch).component(0);

            m_k = basis_second.numElements() -1;
            max_degree = basis_second.maxDegree();
            gsKnotVector<T> kv_geo(0,1,m_k,max_degree+1,max_degree-m_r); // p,r
            gsBSplineBasis<> basis_geo(kv_geo);
            m_basis_geo.push_back(basis_geo);

            m_k = basis_first.numElements() -1;
            max_degree = basis_first.maxDegree();
            gsKnotVector<T> kv_geo2(0,1,m_k,max_degree+1,max_degree-m_r); // p,r
            gsBSplineBasis<> basis_geo2(kv_geo2);
            m_basis_geo.push_back(basis_geo2);
        }

        m_basis_g1 = m_mb; // basis for the g1 basis functions

        refresh();

        // Make local geometry:
/*
        gsMatrix<T> ab_0 = BASIS_tilde_L.support(0);
        gsMatrix<T> ab_1 = BASIS_tilde_R.support(0);

        gsKnotVector<T> kv_0(ab_0.at(0), ab_0.at(1),0, max_degree+1);
        gsKnotVector<T> kv_1(ab_1.at(0), ab_1.at(1),0, max_degree+1);

        for (index_t i = 0; i < BASIS_tilde_L.knots().size(); i++)
        {
            if (BASIS_tilde_L.knot(i) > ab_0.at(0) && BASIS_tilde_L.knot(i) < ab_0.at(1))
                kv_0.insert(BASIS_tilde_L.knot(i),max_degree-m_r);

            if (BASIS_tilde_R.knot(i) > ab_1.at(0) && BASIS_tilde_R.knot(i) < ab_1.at(1))
                kv_1.insert(BASIS_tilde_R.knot(i),max_degree-m_r);
        }
        gsTensorBSplineBasis<2,T> bsp_local_geo(m_basis_geo.at(0).component(0).knots(),kv_0);
        gsTensorBSplineBasis<2,T> bsp_local_geo_R(m_basis_geo.at(1).component(0).knots(),kv_1);
        m_basis_geo_local.push_back(bsp_local_geo);
        m_basis_geo_local.push_back(bsp_local_geo_R);
*/
        //
        // gsInfo << m_basis_geo_local.at(1) << "\n";
    }

    void refresh();
    void assemble();
    inline void apply(bhVisitor & visitor, int patchIndex = 0, boxSide side = boundary::none);

    void solve();
    void constructSolution(gsMultiPatch<T> & result_L, gsMultiPatch<T> & result_R);

    index_t get_n_tilde() { return n_tilde; }
    index_t get_n_bar() { return n_bar; }

    void plotG1Basis(gsMultiPatch<T> & basisG1_L, gsMultiPatch<T> & basisG1_R, std::string baseName)
    {

        const std::string baseName1(baseName + "_L");
        gsParaviewCollection collection1(baseName1);

        const std::string baseName2(baseName + "_R");
        gsParaviewCollection collection2(baseName2);

        std::string fileName, fileName2;
        for (unsigned i = 0; i < basisG1_L.nPatches(); i++)
        {

            fileName = baseName1 + "_" + util::to_string(i);

            gsField<> temp_field_L(m_mp.patch(m_iFace.second().patch),basisG1_L.patch(i));
            gsWriteParaview(temp_field_L,fileName,5000);
            collection1.addTimestep(fileName,i,"0.vts");

        }
        for (unsigned i = 0; i < basisG1_R.nPatches(); i++)
        {

            fileName2 = baseName2 + "_" + util::to_string(i);

            gsField<> temp_field_R(m_mp.patch(m_iFace.first().patch),basisG1_R.patch(i));
            gsWriteParaview(temp_field_R,fileName2,5000);
            collection2.addTimestep(fileName2,i,"0.vts");

        }
        collection1.save();
        collection2.save();
    }

    void g1Condition()
    {
        gsMatrix<> points(1,1000);
        points.setRandom();
        points = points.array().abs();

        gsMatrix<> points2d_L(2, 1000);
        gsMatrix<> points2d_R(2, 1000);
        if (m_iFace.first().side().index() <= 2) // west or east of BOTH patches
        {
            if (m_iFace.first().side().index() == 2) // East
            {
                points2d_L.setOnes();
                points2d_R.setZero();
            }

            if (m_iFace.first().side().index() == 1) // West
            {
                points2d_L.setZero();
                points2d_R.setOnes();
            }

            points2d_L.row(0) = points;
            points2d_R.row(0) = points;
        }
        else if (m_iFace.first().side().index() >= 3) // north or south of BOTH patches
        {
            if (m_iFace.first().side().index() == 4) // North
            {
                points2d_L.setOnes();
                points2d_R.setZero();
            }
            if (m_iFace.first().side().index() == 3) // South
            {
                points2d_L.setZero();
                points2d_R.setOnes();
            }

            points2d_L.row(1) = points;
            points2d_R.row(1) = points;
        }

        if (direct)
        {

        }
        else
        {
            real_t g1Error = 0;
            for (size_t i = 0; i < g1Basis_L.nPatches(); i++)
            {
                gsMatrix<> temp;
                if (m_iFace.first().side().index() <= 3) // north or south of BOTH patches
                    temp = m_gD.get_alpha_tilde_R().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).bottomRows(1))
                        + m_gD.get_alpha_tilde_L().eval(points).cwiseProduct(g1Basis_R.patch(i).deriv(points2d_R).bottomRows(1))
                        + m_gD.get_beta_bar().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).topRows(1));
                else if (m_iFace.first().side().index() <= 2) // west or east of BOTH patches
                    temp = m_gD.get_alpha_tilde_R().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).topRows(1))
                        + m_gD.get_alpha_tilde_L().eval(points).cwiseProduct(g1Basis_R.patch(i).deriv(points2d_R).topRows(1))
                        + m_gD.get_beta_bar().eval(points).cwiseProduct(g1Basis_L.patch(i).deriv(points2d_L).bottomRows(1));

                if (temp.array().abs().maxCoeff() > g1Error)
                    g1Error = temp.array().abs().maxCoeff();
            }

            gsInfo << "Conditiontest G1 continuity: \n" << g1Error << "\n\n";
        }

        //gsInfo << "\nConditiontest G1 continuity: \n" << g1Basis_L.patch(0).coefs() << "\n\n";

    }

protected:

    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;

    boundaryInterface m_iFace;

    index_t p_tilde, r_tilde;
    index_t m_k, m_r;

    bool direct;
    bool local;
    bool plot;


// Gluing data
    gsGluingData_mp<T> m_gD;
    // Spaces for the g1 basis
    gsMultiBasis<> m_basis_g1;
    std::vector<gsBSplineBasis<>> m_basis_geo;
    gsBSplineBasis<> m_basis_tilde;
    gsBSplineBasis<> m_basis_bar;

    std::vector<gsTensorBSplineBasis<2,T>> m_basis_geo_local;

    std::vector<gsSparseSystem<T> > m_phi_i, m_phi_j;

    using Base::m_ddof;

    std::vector<gsMatrix<>> solVec_t, solVec_b;
    gsMultiPatch<T> g1Basis_L, g1Basis_R;

    index_t n_tilde, n_bar;

}; // class gsG1Basis_mp

template <class T, class bhVisitor>
void gsG1BasisLocal_mp<T,bhVisitor>::constructSolution(gsMultiPatch<T> & result_L, gsMultiPatch<T> & result_R)
{

    result_L.clear();
    result_R.clear();

    // Dim is the same for all basis functions
    const index_t dim = ( 0!=solVec_t.at(0).cols() ? solVec_t.at(0).cols() :  m_ddof[0].cols() );

    gsMatrix<T> coeffs;
    for (size_t patchNum = 0; patchNum < 2; patchNum++)
    {
        for (index_t p = 0; p < m_basis_tilde.size(); ++p)
        {
            const gsDofMapper & mapper = m_phi_i.at(p).colMapper(0); // unknown = 0

            // Reconstruct solution coefficients on patch p
            index_t sz;
            if (patchNum == 0)
                 sz = m_basis_g1.basis(m_iFace.second().patch).size();
            if (patchNum == 1)
                 sz = m_basis_g1.basis(m_iFace.first().patch).size();

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
                result_L.addPatch(m_basis_g1.basis(m_iFace.second().patch).makeGeometry(give(coeffs)));
            if (patchNum == 1)
                result_R.addPatch(m_basis_g1.basis(m_iFace.first().patch).makeGeometry(give(coeffs)));
        }

        for (index_t p = 0; p < m_basis_bar.size(); ++p)
        {
            const gsDofMapper & mapper = m_phi_j.at(p).colMapper(0); // unknown = 0

            // Reconstruct solution coefficients on patch p
            index_t sz;
            if (patchNum == 0)
                sz = m_basis_g1.basis(m_iFace.second().patch).size();
            if (patchNum == 1)
                sz = m_basis_g1.basis(m_iFace.first().patch).size();

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
                result_L.addPatch(m_basis_g1.basis(m_iFace.second().patch).makeGeometry(give(coeffs)));
            if (patchNum == 1)
                result_R.addPatch(m_basis_g1.basis(m_iFace.first().patch).makeGeometry(give(coeffs)));
        }
    }

    g1Basis_L = result_L;
    g1Basis_R = result_R;

    //g1Condition();
}


template <class T, class bhVisitor>
void gsG1BasisLocal_mp<T,bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsBasis<> & B_R = m_mb.basis(m_iFace.first().patch); // Patch Right
    gsBasis<> & B_L = m_mb.basis(m_iFace.second().patch); // Patch Left

    gsVector<index_t> sz(2);
    sz(0) = B_L.size();
    sz(1) = B_R.size();

    gsDofMapper map(sz);

    gsMatrix<unsigned> actR, actL;

    for (index_t i = 2; i < B_L.component(0).size(); i++) // B_L == B_R, i = 0,1 interface
    {
        actL = B_L.boundaryOffset(m_iFace.second().side(), i); //second adj. face
        map.markBoundary(0, actL);
    }
    for (index_t i = 2; i < B_R.component(0).size(); i++) // B_L == B_R, i = 0,1 interface
    {
        actR = B_R.boundaryOffset(m_iFace.first() .side(), i); //first adj. face
        map.markBoundary(1, actR);
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
void gsG1BasisLocal_mp<T,bhVisitor>::assemble()
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
    apply(visitor,m_iFace.second().patch); // left
    apply(visitor,m_iFace.first().patch); // right

    for (index_t i = 0; i < m_basis_tilde.size(); i++)
        m_phi_i.at(i).matrix().makeCompressed();
    for (index_t i = 0; i < m_basis_bar.size(); i++)
        m_phi_j.at(i).matrix().makeCompressed();

} // assemble()

template <class T, class bhVisitor>
void gsG1BasisLocal_mp<T,bhVisitor>::apply(bhVisitor & visitor, int patchIndex, boxSide side)
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

    gsBasis<T> & basis_g1 = m_basis_g1.basis(patchIndex); // basis for construction

    index_t leftright;
    if (patchIndex == m_iFace.first().patch)
        leftright = 1;
    if (patchIndex == m_iFace.second().patch)
        leftright = 0;

    // Same for all patches
    gsBasis<T> & basis_geo = m_basis_geo.at(leftright);
    gsBasis<T> & basis_tilde = m_basis_tilde;
    gsBasis<T> & basis_bar = m_basis_bar;

    // Initialize reference quadrature rule and visitor data
    visitor_.initialize(basis_g1, basis_geo, basis_tilde, basis_bar, quRule);

    const gsGeometry<T> & patch = m_mp.patch(patchIndex);

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
        visitor_.evaluate(basis_g1, basis_geo, basis_tilde, basis_bar, patch, quNodes, m_gD, direct, local);

        // Assemble on element
        visitor_.assemble(*domIt, quWeights);

        // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
        visitor_.localToGlobal(leftright, m_ddof, m_phi_i, m_phi_j); // omp_locks inside
    }
}//omp parallel
} // apply

template <class T, class bhVisitor>
void gsG1BasisLocal_mp<T,bhVisitor>::solve()
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

}