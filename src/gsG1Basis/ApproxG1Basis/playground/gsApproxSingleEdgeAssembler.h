/** @file gsApproxSingleEdgeAssembler.h

    @brief Provides assembler for testing the approx spaces.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/


# pragma once

# include <gsAssembler/gsAssembler.h>

# include <gsG1Basis/ApproxG1Basis/playground/gsVisitorApproxSingleEdge.h>
# include <gsG1Basis/ApproxG1Basis/gsApproxGluingData.h>


namespace gismo
{

template <class T, class bhVisitor = gsVisitorApproxSingleEdge<T> >
class gsApproxSingleEdgeAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsApproxSingleEdgeAssembler(gsMultiPatch<T> & mp,
                                gsApproxGluingData<T> gluingData,
                                gsG1OptionList & optionList,
                                index_t uv,
                                index_t bfID)
        : m_mp(mp), m_gluingData(gluingData), m_optionList(optionList), m_uv(uv), m_bfID(bfID)
    {
        // --------------------------------------------- CHANGE ALSO IN L2 NORM IF YOU CHANGE THAT
        index_t m_r = m_optionList.getInt("regularity");

        gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(m_mp.basis(0).component(m_uv)); // 0 -> v, 1 -> u
        index_t m_p = basis_edge.maxDegree();

        // first,last,interior,mult_ends,mult_interior
        gsKnotVector<T> kv_plus(0,1,0,m_p+1,m_p-1-m_r); // p,r+1 //-1 bc r+1
        gsBSplineBasis<> basis_plus(kv_plus);

        for (size_t i = m_p+1; i < basis_edge.knots().size() - (m_p+1); i = i+(m_p-m_r))
            basis_plus.insertKnot(basis_edge.knot(i),m_p-1-m_r);

        m_basis_plus = basis_plus;

        gsMatrix<T> ab = m_basis_plus.support(m_bfID);

        gsKnotVector<T> kv(ab.at(0), ab.at(1), 0, 1);
        for (size_t i = m_p + 1; i < basis_edge.knots().size() - (m_p + 1); i += basis_edge.knots().multiplicityIndex(i))
            if ((basis_edge.knot(i) > ab.at(0)) && (basis_edge.knot(i) < ab.at(1)))
                kv.insert(basis_edge.knot(i), 1);

        gsBSplineBasis<T> bsp_geo_local(kv);


        m_basis_target = basis_edge;
        //m_basis_geo = bsp_geo_local;
        m_basis_geo = basis_edge;

        refresh();
        assemble();
    }

    void refresh();

    void assemble();

    inline void apply(bhVisitor & visitor,
                      int patchIndex = 0,
                      boxSide side = boundary::none);

    void constructSolution(const gsMatrix<> & solVector, gsBSpline<T> & result);

    /// @brief Returns the left-hand global matrix
    const gsSparseMatrix<T> & matrix() const { return m_system.matrix(); }

    /// @brief Returns the left-hand side vector(s)
    const gsMatrix<T> & rhs() const { return m_system.rhs(); }


protected:

    // geometry for computing alpha and beta
    gsMultiPatch<T> m_mp;

    gsApproxGluingData<T> m_gluingData;

    gsG1OptionList m_optionList;

    // Basis for integration
    gsBSplineBasis<T> m_basis_target;
    gsBSplineBasis<T> m_basis_geo;



    index_t m_uv, m_bfID;

    gsBSplineBasis<T> m_basis_plus, m_basis_minus;



    //using Base::m_system;
    using Base::m_ddof;
    using Base::m_system;

}; // class gsApproxProjectionAssembler


template <class T, class bhVisitor>
void gsApproxSingleEdgeAssembler<T,bhVisitor>::constructSolution(const gsMatrix<> & solVector, gsBSpline<T> & result)
{
    // Dim is the same for all basis functions
    const index_t dim = ( 0!=solVector.cols() ? solVector.cols() :  m_ddof[0].cols() );

    const gsDofMapper & mapper = m_system.colMapper(0); // unknown = 0

    // Reconstruct solution coefficients on patch p
    index_t sz;
    sz = m_basis_target.basis(0).size();

    gsMatrix<> coeffs;
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

    result = dynamic_cast<gsBSpline<> &> (*m_basis_target.makeGeometry(coeffs));

}

template <class T, class bhVisitor>
void gsApproxSingleEdgeAssembler<T, bhVisitor>::refresh()
{
    // 1. Obtain a map from basis functions to matrix columns and rows
    gsDofMapper map(m_basis_target);

    gsMatrix<T> ab = m_basis_geo.support();

    gsMatrix<unsigned> act(1,1);

    for (index_t i = 0; i < m_basis_target.size(); i++) // only the first two u/v-columns are Dofs (0/1)
    {
        gsMatrix<T> xy = m_basis_target.support(i);
        // if ( (xy(0, 1) < ab(0, 0)+1e-10) || (xy(0, 0) > ab(0, 1)-1e-10) ) // all non-empty set
        if ( (xy(0, 0) < ab(0, 0)-1e-10) || (xy(0, 1) > ab(0, 1)+1e-10) ) // only subsets
        {
            act << i;
            //map.markBoundary(0, act); // Patch 0
        }
    }

    map.finalize();
    //map_L.print();

    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(map);

} // refresh()

template <class T, class bhVisitor>
void gsApproxSingleEdgeAssembler<T, bhVisitor>::assemble()
{
    //GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    const index_t nz = gsAssemblerOptions::numColNz(m_basis_target,2,1,0.333333);
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
inline void gsApproxSingleEdgeAssembler<T, bhVisitor>::apply(bhVisitor & visitor,
                                                             int patchIndex,
                                                             boxSide side)
{
    gsQuadRule<T> quRule ; // Quadrature rule
    gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
    gsVector<T> quWeights; // Temp variable for mapped weights

    // Initialize reference quadrature rule and visitor data
    visitor.initialize(m_basis_target,quRule);

    // Initialize domain element iterator -- using unknown 0
    typename gsBasis<T>::domainIter domIt = m_basis_geo.makeDomainIterator(boundary::none);

    for (; domIt->good(); domIt->next() )
    {
        // Map the Quadrature rule to the element
        quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

        // Perform required evaluations on the quadrature nodes
        visitor.evaluate(m_basis_target, quNodes, m_mp, m_basis_plus, m_gluingData, m_optionList, m_uv, m_bfID);

        // Assemble on element
        visitor.assemble(*domIt, quWeights);

        // Push to global matrix and right-hand side vector
        visitor.localToGlobal(patchIndex, m_ddof, m_system); // omp_locks inside

    }
}


} // namespace gismo
