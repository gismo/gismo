/** @file gsC1SurfEdge.h

    @brief Creates the (approx) C1 Edge space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/gsC1SurfGluingData.h>
#include <gsUnstructuredSplines/gsC1SurfVisitorBasisEdge.h>

namespace gismo
{
    template<class T, class bhVisitor = gsC1SurfVisitorBasisEdge<T>>
    class gsC1SurfBasisEdge : public gsAssembler<T>
    {
    public:
        typedef gsAssembler<T> Base;

    public:
        gsC1SurfBasisEdge(gsMultiPatch<> mp, // single patch
                        gsMultiBasis<> basis, // single basis
                        index_t uv, // !!! 0 == u; 1 == v !!!
                        bool isBoundary,
                        gsC1SurfGluingData<real_t> gluingD)
                : m_mp(mp), m_basis(basis), m_uv(uv), m_isBoundary(isBoundary), m_gD(gluingD)
        {

            // Computing the G1 - basis function at the edge
            // Spaces for computing the g1 basis
            index_t m_r = 1; // TODO CHANGE IF DIFFERENT REGULARITY IS NECESSARY

            gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(m_basis.basis(0).component(m_uv)); // 0 -> v, 1 -> u
            index_t m_p = basis_edge.maxDegree(); // Minimum degree at the interface // TODO if interface basis are not the same

            // first,last,interior,mult_ends,mult_interior
            gsKnotVector<T> kv_plus(0,1,0,m_p+1,m_p-1-m_r); // p,r+1 //-1 bc r+1
            gsBSplineBasis<> basis_plus(kv_plus);


            for (size_t i = m_p+1; i < basis_edge.knots().size() - (m_p+1); i = i+(m_p-m_r))
                basis_plus.insertKnot(basis_edge.knot(i),m_p-1-m_r);

            m_basis_plus = basis_plus;


            gsKnotVector<T> kv_minus(0,1,0,m_p+1-1,m_p-1-m_r); // p-1,r //-1 bc p-1
            gsBSplineBasis<> basis_minus(kv_minus);

            for (size_t i = m_p+1; i < basis_edge.knots().size() - (m_p+1); i = i+(m_p-m_r))
                basis_minus.insertKnot(basis_edge.knot(i),m_p-1-m_r);

            m_basis_minus = basis_minus;

            // Basis for the G1 basis
            m_basis_g1 = m_basis.basis(0);

//        gsTensorBSplineBasis<2, T> basis_edge_ab = dynamic_cast<gsTensorBSplineBasis<2, real_t> &>(m_basis_g1.basis(0));
//        gsInfo << "Basis edge 0: " << basis_edge_ab.component(0).knots().asMatrix() << "\n";
//        gsInfo << "Basis edge 1: " << basis_edge_ab.component(1).knots().asMatrix() << "\n";
        }

        // Computed the gluing data globally
        void setG1BasisEdge(gsMultiPatch<T> & result);

        void refresh();
        void assemble(index_t i, std::string typeBf); // i == number of bf
        inline void apply(bhVisitor & visitor, index_t i, std::string typeBf); // i == number of bf
        void solve();

        void constructSolution(const gsMatrix<> & solVector, gsMultiPatch<T> & result);

    protected:

        // Input
        gsMultiPatch<T> m_mp;
        gsMultiBasis<T> m_basis;
        index_t m_uv;
        bool m_isBoundary;

        // Gluing data
        gsC1SurfGluingData<T> m_gD;

        // Basis for getting the G1 Basis
        gsBSplineBasis<> m_basis_plus;
        gsBSplineBasis<> m_basis_minus;

        // Basis for the G1 Basis
        gsMultiBasis<T> m_basis_g1;

        // Basis for Integration
        gsMultiBasis<T> m_geo;

        // For Dirichlet boundary
        using Base::m_ddof;
        using Base::m_system;


    }; // class gsG1BasisEdge

    template <class T, class bhVisitor>
    void gsC1SurfBasisEdge<T,bhVisitor>::setG1BasisEdge(gsMultiPatch<T> & result)
    {
        result.clear();

        index_t n_plus = m_basis_plus.size();
        index_t n_minus = m_basis_minus.size();

        gsMultiPatch<> g1EdgeBasis;
        index_t bfID_init = 3;

        for (index_t bfID = bfID_init; bfID < n_plus - bfID_init; bfID++) // first 3 and last 3 bf are eliminated
        {
            m_geo = m_basis_g1; // Basis for Integration

            refresh();

            assemble(bfID,"plus"); // i == number of bf

            gsSparseSolver<real_t>::CGDiagonal solver;
//        gsSparseSolver<real_t>::LU solver;
            gsMatrix<> sol;
            solver.compute(m_system.matrix());
            sol = solver.solve(m_system.rhs());

            constructSolution(sol,g1EdgeBasis);
        }
        bfID_init = 2;
        for (index_t bfID = bfID_init; bfID < n_minus-bfID_init; bfID++)  // first 2 and last 2 bf are eliminated
        {


            m_geo = m_basis_g1; // Basis for Integration

            refresh();

            assemble(bfID,"minus"); // i == number of bf

            gsSparseSolver<real_t>::CGDiagonal solver;
//        gsSparseSolver<real_t>::LU solver;
            gsMatrix<> sol;
            solver.compute(m_system.matrix());
            sol = solver.solve(m_system.rhs());

            constructSolution(sol,g1EdgeBasis);
        }

        result = g1EdgeBasis;
    } // setG1BasisEdge

    template <class T, class bhVisitor>
    void gsC1SurfBasisEdge<T,bhVisitor>::constructSolution(const gsMatrix<> & solVector, gsMultiPatch<T> & result)
    {
        // Dim is the same for all basis functions
        const index_t dim = ( 0!=solVector.cols() ? solVector.cols() :  m_ddof[0].cols() );

        gsMatrix<T> coeffs;
        const gsDofMapper & mapper = m_system.colMapper(0); // unknown = 0

        // Reconstruct solution coefficients on patch p
        index_t sz;
        sz = m_basis.basis(0).size();

        coeffs.resize(sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if (mapper.is_free(i, 0)) // DoF value is in the solVector // 0 = unitPatch
            {
                coeffs.row(i) = solVector.row(mapper.index(i, 0));
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                //gsInfo << "mapper index dirichlet: " << m_ddof[unk].row( mapper.bindex(i, p) ).head(dim) << "\n";
                coeffs.row(i) = m_ddof[0].row( mapper.bindex(i, 0) ).head(dim); // = 0
            }
        }
        result.addPatch(m_basis_g1.basis(0).makeGeometry(give(coeffs)));

    }

    template <class T, class bhVisitor>
    void gsC1SurfBasisEdge<T,bhVisitor>::refresh()
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

        // 2. Create the sparse system
        m_system = gsSparseSystem<T>(map);

    } // refresh()

    template <class T, class bhVisitor>
    void gsC1SurfBasisEdge<T,bhVisitor>::assemble(index_t bfID, std::string typeBf)
    {
        // Reserve sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_basis[0],2,1,0.333333);
        m_system.reserve(nz, 1);

        if(m_ddof.size()==0)
            m_ddof.resize(1); // 0,1

        const gsDofMapper & map = m_system.colMapper(0); // Map same for every functions

        m_ddof[0].setZero(map.boundarySize(), 1 );

        // Assemble volume integrals
        bhVisitor visitor;
        apply(visitor, bfID, typeBf); // basis function i

        m_system.matrix().makeCompressed();

    } // assemble()

    template <class T, class bhVisitor>
    void gsC1SurfBasisEdge<T,bhVisitor>::apply(bhVisitor & visitor, int bf_index, std::string typeBf)
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

                // Perform required evaluations on the quadrature nodes
                visitor_.evaluate(bf_index, typeBf, basis_g1, basis_geo, basis_plus, basis_minus, patch, quNodes, m_uv, m_gD, m_isBoundary);

                // Assemble on element
                visitor_.assemble(*domIt, quWeights);

                // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
                visitor_.localToGlobal(0, m_ddof, m_system); // omp_locks inside // patchIndex == 0
            }
        }//omp parallel
    } // apply

} // namespace gismo
