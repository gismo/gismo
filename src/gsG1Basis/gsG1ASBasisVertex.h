//
// Created by afarahat on 4/16/20.
//


#pragma once

#include <gsG1Basis/gsG1ASGluingData.h>
#include <gsG1Basis/gsG1ASVisitorBasisVertex.h>
#include <gsG1Basis/gsG1OptionList.h>

namespace gismo
{
template<class T, class bhVisitor = gsG1ASVisitorBasisVertex<T>>
class gsG1ASBasisVertex : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsG1ASBasisVertex(gsMultiPatch<> mp, // Single Patch
                    gsMultiBasis<> basis, // Single Basis
                    std::vector<bool> isBoundary,
                    real_t sigma,
                    gsG1OptionList & g1OptionList,
                    gsMatrix<> gluingD)
        : m_mp(mp), m_basis(basis), m_isBoundary(isBoundary), m_sigma(sigma), m_g1OptionList(g1OptionList), m_gD(gluingD)
    {

        for (index_t dir = 0; dir < m_mp.parDim(); dir++) // For the TWO directions
        {
            // Computing the G1 - basis function at the edge
            // Spaces for computing the g1 basis
            index_t m_r = m_g1OptionList.getInt("regularity"); // TODO CHANGE IF DIFFERENT REGULARITY IS NECESSARY

            gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(dir)); // 0 -> u, 1 -> v
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


    }


    void refresh();
    void assemble();
    inline void apply(bhVisitor & visitor, int patchIndex);
    void solve();

    void constructSolution(gsMultiPatch<T> & result);

    void setG1BasisVertex(gsMultiPatch<T> & result, index_t kinfOfVertex)
    {
        // SIMPLIFY THE m_geo
        m_geo = m_basis_g1; // INTEGRATION AREA

        refresh();
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
    gsMatrix<> m_gD;

    // Basis for getting the G1 Basis
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
void gsG1ASBasisVertex<T,bhVisitor>::constructSolution(gsMultiPatch<T> & result)
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
}

template <class T, class bhVisitor>
void gsG1ASBasisVertex<T,bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map(m_basis.basis(0));

    // SET THE DOFS


    map.finalize();


    // 2. Create the sparse system
    gsSparseSystem<T> m_system = gsSparseSystem<T>(map);
    for (index_t i = 0; i < 6; i++)
        m_f.push_back(m_system);

} // refresh()

template <class T, class bhVisitor>
void gsG1ASBasisVertex<T,bhVisitor>::assemble()
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
void gsG1ASBasisVertex<T,bhVisitor>::apply(bhVisitor & visitor, int patchIndex)
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
void gsG1ASBasisVertex<T,bhVisitor>::solve()
{
    gsSparseSolver<real_t>::CGDiagonal solver;


    for (index_t i = 0; i < 6; i++) // Tilde
    {
        solver.compute(m_f.at(i).matrix());
        solVec.push_back(solver.solve(m_f.at(i).rhs()));
    }
} // solve

} // namespace gismo
