/** @file gsApproxGDAssembler.h

    @brief Provides assembler for the gluing data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

#include <gsG1Basis/gsVisitorLocalGD.h>

namespace gismo
{

template <class T, class bhVisitor = gsVisitorLocalGD<T> >
class gsLocalGDAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsLocalGDAssembler(gsBSplineBasis<T> const & basis,
                       gsBSplineBasis<T> const & basis_geo,
                       gsMultiPatch<T> const & mp,
                       boundaryInterface const & iFace,
                       index_t const & gamma,
                       std::string const & alphaBeta)
        : m_basis_geo(basis_geo), m_mp(mp), m_iFace(iFace), m_gamma(gamma), m_alphaBeta(alphaBeta)
    {
        m_basis.push_back(basis);
        refresh();
    }

    void refresh();

    void assemble();

    inline void apply(bhVisitor & visitor,
                      int patchIndex = 0,
                      boxSide side = boundary::none);

    void constructSolution(const std::vector<gsMatrix<T>>& solVector,
                           const std::vector<gsMatrix<T>>& solVector_b,
                           gsMultiPatch<T>& result, int unk = 0);

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix_L() const { return m_system_L.matrix(); }

    /// @brief Returns the left-hand side vector(s)
    const gsMatrix<T> & rhs_L() const { return m_system_L.rhs(); }

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix_R() const { return m_system_R.matrix(); }

    /// @brief Returns the left-hand side vector(s)
    const gsMatrix<T> & rhs_R() const { return m_system_R.rhs(); }

protected:
    // Basis for integration
    gsBSplineBasis<T> m_basis_geo;

    // interface + geometry for computing alpha and beta
    gsMultiPatch<T> m_mp;
    boundaryInterface m_iFace;

    index_t m_gamma;
    std::string m_alphaBeta;

    // Space for phi_0,i, phi_1,j
    std::vector< gsMultiBasis<T> > m_basis;

    //using Base::m_system;
    using Base::m_ddof;

    gsSparseSystem<T> m_system_L;
    gsSparseSystem<T> m_system_R;



}; // class gsLocalGDAssembler


template <class T, class bhVisitor>
void gsLocalGDAssembler<T, bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map_L(m_basis[0]);
    gsDofMapper map_R(m_basis[0]);


    map_L.finalize();
    //map_L.print();

    map_R.finalize();
    //map_R.print();

    // 2. Create the sparse system
    m_system_L = gsSparseSystem<T>(map_L);
    m_system_R = gsSparseSystem<T>(map_R);

} // refresh()

template <class T, class bhVisitor>
void gsLocalGDAssembler<T, bhVisitor>::assemble()
{
    //GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis[0][0],2,1,0.333333);
    m_system_L.reserve(nz, 1);
    m_system_R.reserve(nz, 1);


    if(m_ddof.size()==0)
        m_ddof.resize(2); // One for L, one for R

    const gsDofMapper & mapper_L = m_system_L.colMapper(0);
    const gsDofMapper & mapper_R = m_system_R.colMapper(0);

    m_ddof[0].setZero(mapper_L.boundarySize(), 1 ); // tilde
    m_ddof[1].setZero(mapper_R.boundarySize(), 1 ); // bar

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor);

    m_system_L.matrix().makeCompressed();
    m_system_R.matrix().makeCompressed();
}

template <class T, class bhVisitor>
inline void gsLocalGDAssembler<T, bhVisitor>::apply(bhVisitor & visitor,
                                                    int patchIndex,
                                                    boxSide side)
{
    gsQuadRule<T> quRule ; // Quadrature rule
    gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
    gsVector<T> quWeights; // Temp variable for mapped weights

    gsBasis<T> & basis = m_basis[0].basis(patchIndex);

    // Initialize reference quadrature rule and visitor data
    visitor.initialize(basis,quRule);

    //const gsGeometry<T> & patch = m_geo.patch(patchIndex); // 0 = patchindex

    // Initialize domain element iterator -- using unknown 0
    typename gsBasis<T>::domainIter domIt = m_basis_geo.makeDomainIterator(boundary::none);

    for (; domIt->good(); domIt->next() )
    {
        // Map the Quadrature rule to the element
        quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

        // Perform required evaluations on the quadrature nodes
        visitor.evaluate(basis, quNodes,  m_mp, m_iFace, m_gamma, m_alphaBeta);

        // Assemble on element
        visitor.assemble(*domIt, quWeights);

        // Push to global matrix and right-hand side vector
        visitor.localToGlobal(patchIndex, m_ddof, m_system_L, m_system_R); // omp_locks inside

    }
}


} // namespace gismo
