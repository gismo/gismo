/** @file gsGluingDataAssembler.h

    @brief Provides assembler for the gluing data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


#pragma once

#include <gsG1Basis/gsVisitorGlobalApproxGD4.h>
#include <gsAssembler/gsAssembler.h>
#include <gsG1Basis/gsG1OptionList.h>

namespace gismo
{

template <class T, class bhVisitor = gsVisitorGlobalApproxGD4<T> >
class gsGlobalGDAssembler4 : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsGlobalGDAssembler4(gsMultiPatch<T> const & mp,
                         gsBSplineBasis<T> const & basis,
                         gsG1OptionList optionList)
        : m_mp(mp), m_optionList(optionList)
    {
        m_basis.push_back(basis); // Alpha^L
        m_basis.push_back(basis); // Alpha^R

        refresh();
    }

    void refresh();

    void assemble();

    inline void apply(bhVisitor & visitor);

protected:
    // interface + geometry for computing alpha
    gsMultiPatch<T> m_mp;

    gsG1OptionList m_optionList;

    // Space for alpha^S, beta^S
    std::vector<gsBSplineBasis<>> m_basis;

    //using Base::m_system;
    using Base::m_ddof;
    using Base::m_system;

}; // class gsG1BasisAssembler


template <class T, class bhVisitor>
void gsGlobalGDAssembler4<T, bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsVector<> size;
    size.setZero(1);
    for (size_t i = 0; i < m_basis.size(); i++)
        size.at(0) += m_basis[i].size(); // alpha^S

    gsDofMapper map(size);

    map.finalize();

    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(map);


} // refresh()

template <class T, class bhVisitor>
void gsGlobalGDAssembler4<T, bhVisitor>::assemble()
{
    //GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis[0],2,1,0.333333);
    m_system.reserve(nz, 1);


    if(m_ddof.size()==0)
        m_ddof.resize(1);

    const gsDofMapper & mapper = m_system.colMapper(0);

    m_ddof[0].setZero(mapper.boundarySize(), 1 );

    // Assemble volume integrals
    bhVisitor visitor;
    apply(visitor);

    m_system.matrix().makeCompressed();

}

template <class T, class bhVisitor>
inline void gsGlobalGDAssembler4<T, bhVisitor>::apply(bhVisitor & visitor)
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

        gsBasis<T> & basis = m_mp.basis(0).component(1); // = 0

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
            visitor_.evaluate(m_basis, quNodes, m_mp, m_optionList);

            // Assemble on element
            visitor_.assemble(*domIt, quWeights);

            // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
            visitor_.localToGlobal(0, m_ddof, m_system); // omp_locks inside

        }
    }//omp parallel
}


} // namespace gismo