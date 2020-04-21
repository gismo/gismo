//
// Created by afarahat on 3/25/20.
//


#pragma once

#include <gsG1Basis/gsG1ASGluingDataVisitorGlobal.h>

namespace gismo
{

template <class T, class bhVisitor = gsG1ASGluingDataVisitorGlobal<T> >
class gsG1ASGluingDataAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsG1ASGluingDataAssembler(gsBasis<T> const & basis,
                        index_t const & uv,
                        gsMultiPatch<T> const & mp,
                        index_t const & gamma,
                        bool & isBoundary)
        : m_uv(uv), m_mp(mp), m_gamma(gamma), m_isBoundary(isBoundary)
    {

        m_basis.push_back(basis); // Basis for alpha and beta

        refresh();
    }

    void refresh();
    void refreshBeta();


    void assemble();
    void assembleBeta();


    inline void apply(bhVisitor & visitor,
                      int patchIndex = 0,
                      boxSide side = boundary::none);

    inline void applyBeta(bhVisitor & visitor,
                      int patchIndex = 0,
                      boxSide side = boundary::none);

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix_alpha() const { return m_system_alpha.matrix(); }
    const gsMatrix<T> & rhs_alpha() const { return m_system_alpha.rhs(); }

    const gsSparseMatrix<T> & matrix_beta() const { return m_system_beta_S.matrix(); }
    const gsMatrix<T> & rhs_beta() const { return m_system_beta_S.rhs(); }

protected:
    index_t m_uv;

    // interface + geometry for computing alpha and beta
    gsMultiPatch<T> m_mp;

    index_t m_gamma;
    bool m_isBoundary;

    // Space for phi_0,i, phi_1,j
    std::vector< gsMultiBasis<T> > m_basis;

    //using Base::m_system;
    using Base::m_ddof;

    gsSparseSystem<T> m_system_alpha;
    gsSparseSystem<T> m_system_beta_S;

}; // class gsG1BasisAssembler


template <class T, class bhVisitor>
void gsG1ASGluingDataAssembler<T, bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map_alpha(m_basis[0]);

    gsDofMapper map_beta_S(m_basis[0]);


    map_alpha.finalize();
    map_beta_S.finalize();

    // 2. Create the sparse system
    m_system_alpha = gsSparseSystem<T>(map_alpha);
    m_system_beta_S = gsSparseSystem<T>(map_beta_S);

} // refresh()



template <class T, class bhVisitor>
void gsG1ASGluingDataAssembler<T, bhVisitor>::assemble()
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

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor);

    m_system_alpha.matrix().makeCompressed();
    m_system_beta_S.matrix().makeCompressed();

}



template <class T, class bhVisitor>
inline void gsG1ASGluingDataAssembler<T, bhVisitor>::apply(bhVisitor & visitor,
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

        gsBasis<T> & basis = m_basis[0].basis(patchIndex); // = 0

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
            visitor_.evaluate(basis, quNodes, m_uv, m_mp, m_gamma, m_isBoundary);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(patchIndex, m_ddof, m_system_alpha, m_system_beta_S); // omp_locks inside

        }

    }//omp parallel
}

} // namespace gismo
