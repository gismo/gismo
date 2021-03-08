/** @file gsApproxArgyrisEdgeBasis.h

    @brief Provides assembler for a G1 Basis for multiPatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once


#include <gsArgyris/gsGluingData/gsApproxGluingData.h>
#include <gsArgyris/gsApproxArgyrisEdgeBasisVisitor.h>



namespace gismo
{
template<short_t d, class T, class bhVisitor = gsApproxArgyrisEdgeBasisVisitor<d, T>>
class gsApproxArgyrisEdgeBasis : public gsAssembler<T>
{
private:
    typedef typename std::vector<gsC1ArgyrisAuxiliaryPatch<d,T>> ArgyrisAuxPatchContainer;

public:
    typedef gsAssembler<T> Base;

public:
    gsApproxArgyrisEdgeBasis() { };

    gsApproxArgyrisEdgeBasis(ArgyrisAuxPatchContainer & auxPatches,
                        gsApproxGluingData<d, T> & approxGluingData,
                        index_t patchID,
                        const gsOptionList & optionList,
                        const bool & isBoundary = false)
        : m_auxPatches(auxPatches), m_approxGluingData(approxGluingData), m_patchID(patchID), m_optionList(optionList), m_isBoundary(isBoundary)
    {

        m_side = m_auxPatches[patchID].side();
        m_dir = patchID == 0 ? 1 : 0;

        // Collect the needed basis
        m_basis_g1 = m_auxPatches[patchID].getArygrisBasisRotated().getEdgeBasis(m_side);

        m_basis_plus = m_auxPatches[patchID].getArygrisBasisRotated().getBasisPlus(m_side);
        m_basis_minus = m_auxPatches[patchID].getArygrisBasisRotated().getBasisMinus(m_side);
        m_basis_geo = m_auxPatches[patchID].getArygrisBasisRotated().getBasisGeo(m_side);
    }

    gsApproxArgyrisEdgeBasis(ArgyrisAuxPatchContainer & auxPatches,
                    index_t patchID,
                    const gsOptionList & optionList,
                    const bool & isBoundary = true)
    : m_auxPatches(auxPatches), m_patchID(patchID), m_optionList(optionList), m_isBoundary(isBoundary)
    {
        m_side = m_auxPatches[patchID].side();
        m_dir = patchID == 0 ? 1 : 0;

        // Collect the needed basis
        m_basis_g1 = m_auxPatches[patchID].getArygrisBasisRotated().getEdgeBasis(m_side);

        m_basis_plus = m_auxPatches[patchID].getArygrisBasisRotated().getBasisPlus(m_side);
        m_basis_minus = m_auxPatches[patchID].getArygrisBasisRotated().getBasisMinus(m_side);
        m_basis_geo = m_auxPatches[patchID].getArygrisBasisRotated().getBasisGeo(m_side);
    }

    // Computed the gluing data globally
    void setG1BasisEdge(gsMultiPatch<T> & result);

    void refresh(index_t bfID, std::string typeBf);
    void assemble(index_t bfID, std::string typeBf); // i == number of bf
    inline void apply(bhVisitor & visitor, index_t bfID, std::string typeBf); // i == number of bf

    void constructSolution(const gsMatrix<> & solVector, gsMultiPatch<T> & result);

protected:

    // Input
    ArgyrisAuxPatchContainer m_auxPatches;
    gsApproxGluingData<d, T> m_approxGluingData;
    index_t m_patchID;
    const gsOptionList m_optionList;
    const bool m_isBoundary;

    // Side for access the AuxContainer
    index_t m_side, m_dir;

    // Basis for the G1 Basis
    gsTensorBSplineBasis<d, T> m_basis_g1;

    // Basis need for construction
    gsBSplineBasis<T> m_basis_plus;
    gsBSplineBasis<T> m_basis_minus;
    gsBSplineBasis<T> m_basis_geo;

    // For Dirichlet boundary
    using Base::m_ddof;
    using Base::m_system;

}; // class gsG1BasisEdge

template <short_t d, class T, class bhVisitor>
void gsApproxArgyrisEdgeBasis<d, T,bhVisitor>::setG1BasisEdge(gsMultiPatch<T> & result)
{
    result.clear();

    index_t n_plus = m_basis_plus.size();
    index_t n_minus = m_basis_minus.size();

    float progress = 0.0;
    int barWidth = 50;

    gsMultiPatch<> g1EdgeBasis, g1EdgeBasis2;
    index_t bfID_init = 3;
    if (m_optionList.getSwitch("twoPatch"))
        bfID_init = 0;

    if (m_isBoundary && m_optionList.getSwitch("twoPatch"))
        bfID_init = 2;

    for (index_t bfID = bfID_init; bfID < n_plus - bfID_init; bfID++) // first 3 and last 3 bf are eliminated
    {
        if (m_optionList.getSwitch("info"))
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

        refresh(bfID,"plus");

        assemble(bfID,"plus"); // i == number of bf

        //gsSparseSolver<real_t>::CGDiagonal solver;
        gsSparseSolver<real_t>::LU solver;
        gsMatrix<real_t> sol;
        solver.compute(m_system.matrix());
        sol = solver.solve(m_system.rhs());

        constructSolution(sol,g1EdgeBasis);

        progress += (float) 1.0 / (n_plus+n_minus); // for demonstration only
    }


    bfID_init = 2;
    if (m_optionList.getSwitch("twoPatch"))
        bfID_init = 0;

    if (m_isBoundary && m_optionList.getSwitch("twoPatch"))
        bfID_init = 2;

    for (index_t bfID = bfID_init; bfID < n_minus-bfID_init; bfID++)  // first 2 and last 2 bf are eliminated
    {
        if (m_optionList.getSwitch("info"))
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

        refresh(bfID,"minus");

        assemble(bfID,"minus"); // i == number of bf

        //gsSparseSolver<real_t>::CGDiagonal solver;
        gsSparseSolver<real_t>::LU solver;
        gsMatrix<> sol;
        solver.compute(m_system.matrix());
        sol = solver.solve(m_system.rhs());

        constructSolution(sol,g1EdgeBasis2);
        g1EdgeBasis.addPatch(g1EdgeBasis2.patch(bfID-bfID_init));

        progress += (float) 1.0 / (n_plus+n_minus); // for demonstration only
    }


    result = g1EdgeBasis;

    if (m_optionList.getSwitch("info"))
        std::cout << std::endl;

} // setG1BasisEdge

template <short_t d, class T, class bhVisitor>
void gsApproxArgyrisEdgeBasis<d, T,bhVisitor>::constructSolution(const gsMatrix<> & solVector, gsMultiPatch<T> & result)
{
    // Dim is the same for all basis functions
    const index_t dim = ( 0!=solVector.cols() ? solVector.cols() :  m_ddof[0].cols() );

    gsMatrix<T> coeffs;
    const gsDofMapper & mapper = m_system.colMapper(0); // unknown = 0

    // Reconstruct solution coefficients on patch p
    index_t sz;
    sz = m_basis_g1.size();

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

    result.addPatch(m_basis_g1.makeGeometry(give(coeffs)));

}

template <short_t d, class T, class bhVisitor>
void gsApproxArgyrisEdgeBasis<d, T,bhVisitor>::refresh(index_t bfID, std::string typeBf)
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map(m_basis_g1);

    gsMatrix<index_t> act;

    for (index_t i = 2; i < m_basis_g1.component(1 - m_dir).size();
         i++) // only the first two u/v-columns are Dofs (0/1)
    {
        act = m_basis_g1.boundaryOffset(m_dir == 0 ? 3 : 1, i); // WEST
        map.markBoundary(0, act); // Patch 0
    }

    map.finalize();

    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(map);

} // refresh()

template <short_t d, class T, class bhVisitor>
void gsApproxArgyrisEdgeBasis<d, T,bhVisitor>::assemble(index_t bfID, std::string typeBf)
{
    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis_g1,2,1,0.333333);
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

template <short_t d, class T, class bhVisitor>
void gsApproxArgyrisEdgeBasis<d, T,bhVisitor>::apply(bhVisitor & visitor, int bf_index, std::string typeBf)
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

        gsBasis<T> & basis_g1 = m_basis_g1; // basis for construction

        // Initialize reference quadrature rule and visitor data
        visitor_.initialize(basis_g1, quRule);

        const gsGeometry<T> & patch = m_auxPatches[m_patchID].getPatch();

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
            visitor_.evaluate(patch, m_basis_g1, m_basis_plus, m_basis_minus, m_basis_geo, m_approxGluingData,
                              quNodes, bf_index, typeBf, m_dir, m_isBoundary, m_optionList);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(0, m_ddof, m_system); // omp_locks inside // patchIndex == 0
        }
    }//omp parallel
} // apply



} // namespace gismo