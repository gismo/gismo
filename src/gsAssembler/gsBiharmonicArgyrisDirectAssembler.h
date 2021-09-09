/** @file gsBiharmonicArgyrisDirectAssembler.h
    @brief Provides assembler for a homogenius Biharmonic equation.
    This file is part of the G+Smo library.
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    Author(s): J. Sogn & P. Weinmüller
*/


#pragma once

#include <gsAssembler/gsAssembler.h>

#include <gsPde/gsBiharmonicPde.h>

#include <gsAssembler/gsVisitorBiharmonicArgyrisDirect.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorLaplaceBoundaryBiharmonic.h>
//#include <gsAssembler/gsVisitorNitscheBiharmonic.h>

#include <gsC1BasisDirect/gsG1MultiBasis.h>

namespace gismo
{

/** @brief
    Implementation of a homogeneous Biharmonic Assembler.
    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can only be enforced strongly (i.e Nitsche is
    not implemented).
*/
    template <class T, class bhVisitor = gsVisitorBiharmonicArgyrisDirect<T> >
    class gsBiharmonicArgyrisDirectAssembler : public gsAssembler<T>
    {
    public:
        typedef gsAssembler<T> Base;

    public:
/** @brief
    Constructor of the assembler object.
    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] bconditions  is a gsBoundaryConditions object that holds boundary conditions on the form:
    \f[ \text{Dirichlet: } u = g \text{ on } \Gamma, \text{ and Neumann: } \nabla \Delta u \cdot \mathbf{n} = h \text{ on } \Gamma\f]
    \param[in] bconditions2 is a gsBoundaryConditions object that holds Neumann boundary conditions on the form:
    \f[\text{Neumann: } \nabla \Delta u \cdot \mathbf{n} = g\, \rightarrow \,(g,\nabla v \cdot \mathbf{n})_\Gamma, \f] where \f$ g \f$ is the Neumann data,
    \f$ v \f$ is the test function and \f$ \Gamma\f$ is the boundary side.
    \param[in] rhs is the right-hand side of the Biharmonic equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary in the \em bconditions object.
    \param[in] intStrategy option for the treatment of patch interfaces
*/
        gsBiharmonicArgyrisDirectAssembler( gsMultiPatch<T> const         & patches,
                               gsMultiBasis<T> const         & bases,
                               gsBoundaryConditions<T> const & bconditions,
                               gsBoundaryConditions<T> const & bconditions2,
                               const gsFunction<T>           & rhs,
                               dirichlet::strategy           dirStrategy = dirichlet::elimination,
                               iFace::strategy               intStrategy = iFace::none)
                : m_ppde(patches,bconditions,bconditions2,rhs)
        {
            m_options.setInt("DirichletStrategy", dirStrategy);
            m_options.setInt("InterfaceStrategy", intStrategy);

            Base::initialize(m_ppde, bases, m_options);
        }

        void refresh();

        void assemble();

        void push();

        void apply(bhVisitor & visitor,
                   size_t patchIndex,
                   boxSide side= boundary::none);

        void constructSolution(const gsMatrix<T>& solVector,
                               gsMultiPatch<T>& result, short_t unk = 0);

        void constructG1Solution(const gsMatrix<T>& solVector,
                                 gsMatrix<T>& result, short_t unk = 0);

        void computeDirichletAndNeumannDofs();

        void plotParaview(gsField<> &solField_interior, gsMatrix<T> &result);

    protected:

        // fixme: add constructor and remove this
        gsBiharmonicPde<T> m_ppde;

        // Members from gsAssembler
        using Base::m_pde_ptr;
        using Base::m_bases;
        using Base::m_ddof;
        using Base::m_options;
        using Base::m_system;
    };

    template <class T, class bhVisitor>
    void gsBiharmonicArgyrisDirectAssembler<T,bhVisitor>::refresh()
    {
        // We use predefined helper which initializes the system matrix
        // rows and columns using the same test and trial space

        // Pascal
        if (m_ppde.domain().nPatches() == 1)
        {
            gsDofMapper map(m_bases[0]);

            gsMatrix<index_t> act;
            for (typename gsBoundaryConditions<T>::bcContainer::const_iterator it
                    = m_ppde.bcFirstKind().dirichletSides().begin();
                 it != m_ppde.bcFirstKind().dirichletSides().end(); ++it) {
                act = m_bases[0].basis(it->patch()).boundaryOffset(it->side().index(), 0); // First
                map.markBoundary(it->patch(), act);
            }

            for (typename gsBoundaryConditions<T>::bcContainer::const_iterator it
                    = m_ppde.bcSecondKind().neumannSides().begin();
                 it != m_ppde.bcSecondKind().neumannSides().end(); ++it) {
                act = m_bases[0].basis(it->patch()).boundaryOffset(it->side().index(), 1); // Second
                // without the first and the last (already marked from dirichlet boundary)
                map.markBoundary(it->patch(), act.block(1, 0, act.rows() - 2, 1));
            }

            map.finalize();

            // 2. Create the sparse system
            m_system = gsSparseSystem<T>(map);
        }
        else if (m_ppde.domain().nPatches() == 2)
        {
            gsVector<index_t> sz(m_ppde.domain().nPatches());
            sz.setZero();
            for (size_t np = 0; np < m_ppde.domain().nPatches(); np++)
                sz[np] = m_bases[0].basis(np).size();  // 0 == one unknown

            gsVector<index_t> sz_int(m_ppde.domain().interfaces().size());
            sz_int.setZero();
            for (size_t numInt = 0; numInt < m_ppde.domain().interfaces().size(); numInt++)
            {
                const boundaryInterface &item = m_ppde.domain().interfaces()[numInt];

                // Get the dimension for the spaces at the patch
                index_t dir = item.first().m_index < 3 ? 1 : 0;
                gsBSplineBasis<> basis_edge = dynamic_cast<gsBSplineBasis<> &>(m_bases[0].basis(item.first().patch).component(dir)); // If the interface matches!!!
                index_t m_p = basis_edge.maxDegree();
                index_t m_r = m_p - basis_edge.knots().multiplicityIndex(m_p + 1) > m_p - 2 ? m_p - 2 : m_p - basis_edge.knots().multiplicityIndex(m_p + 1);
                index_t m_n = basis_edge.numElements();

                sz_int[numInt] = 2*(m_p - m_r - 1) * (m_n - 1) + 2*m_p + 1;
                sz[item.first().patch] += sz_int[numInt];
                sz[item.second().patch] += sz_int[numInt]; // Should be the same!!!
            }

            gsDofMapper map(sz);

            // Boundary
            gsMatrix<index_t> act;
            for (typename gsBoundaryConditions<T>::bcContainer::const_iterator it
                    = m_ppde.bcFirstKind().dirichletSides().begin();
                 it != m_ppde.bcFirstKind().dirichletSides().end(); ++it)
            {
                act = m_bases[0].basis(it->patch()).boundaryOffset(it->side().index(), 0); // Dirichlet
                map.markBoundary(it->patch(), act);
            }

            // Delete all dofs at the interfaces
            for (size_t numInt = 0; numInt < m_ppde.domain().interfaces().size(); numInt++)
            {
                const boundaryInterface &item = m_ppde.domain().interfaces()[numInt];

                gsMatrix<index_t> act;
                act = m_bases[0].basis(item.first().patch).boundaryOffset(item.first().index(), 0); // First
                map.markBoundary(item.first().patch, act);
                act = m_bases[0].basis(item.first().patch).boundaryOffset(item.first().index(), 1); // Second
                map.markBoundary(item.first().patch, act);

                act = m_bases[0].basis(item.second().patch).boundaryOffset(item.second().index(), 0); // First
                map.markBoundary(item.second().patch, act);
                act = m_bases[0].basis(item.second().patch).boundaryOffset(item.second().index(), 1); // Second
                map.markBoundary(item.second().patch, act);

                // Match interface (only for dofs, not for bdy)
                for (index_t i = 0; i < sz_int[numInt]; i++)
                    map.matchDof(item.first().patch,m_bases[0].basis(item.first().patch).size()+i,
                                 item.second().patch,m_bases[0].basis(item.second().patch).size()+i);


                // Mark g1 basis dofs as bdy TODO Only for two patch!!!
                index_t sz_plus = (sz_int[numInt] + 1)/2;
                map.eliminateDof(m_bases[0].basis(item.first().patch).size(), item.first().patch);
                map.eliminateDof(m_bases[0].basis(item.first().patch).size()+sz_plus-1, item.first().patch);
                map.eliminateDof(m_bases[0].basis(item.first().patch).size()+sz_plus, item.first().patch);
                map.eliminateDof(m_bases[0].basis(item.first().patch).size()+sz_int[numInt]-1, item.first().patch);

                map.eliminateDof(m_bases[0].basis(item.second().patch).size(), item.second().patch);
                map.eliminateDof(m_bases[0].basis(item.second().patch).size()+sz_plus-1, item.second().patch);
                map.eliminateDof(m_bases[0].basis(item.second().patch).size()+sz_plus, item.second().patch);
                map.eliminateDof(m_bases[0].basis(item.second().patch).size()+sz_int[numInt]-1, item.second().patch);
            }

            map.finalize();

            // Mark all dofs at the interfaces as tagged
            for (size_t numInt = 0; numInt < m_ppde.domain().interfaces().size(); numInt++)
            {
                const boundaryInterface &item = m_ppde.domain().interfaces()[numInt];

                gsMatrix<index_t> act;
                act = m_bases[0].basis(item.first().patch).boundaryOffset(item.first().index(), 0); // First
                for (index_t i = 0; i < act.rows(); ++i)
                    map.markTagged(act(i,0), item.first().patch);

                act = m_bases[0].basis(item.first().patch).boundaryOffset(item.first().index(), 1); // Second
                for (index_t i = 0; i < act.rows(); ++i)
                    map.markTagged(act(i,0), item.first().patch);

                act = m_bases[0].basis(item.second().patch).boundaryOffset(item.second().index(), 0); // First
                for (index_t i = 0; i < act.rows(); ++i)
                    map.markTagged(act(i,0), item.second().patch);

                act = m_bases[0].basis(item.second().patch).boundaryOffset(item.second().index(), 1); // Second
                for (index_t i = 0; i < act.rows(); ++i)
                    map.markTagged(act(i,0), item.second().patch);

            }
/*
            map.print();
            gsInfo << map.asVector() << "\n";
            for (index_t i = 0; i < map.taggedSize(); i++)
                gsInfo << "Tagged: " << map.getTagged()[i] << "\n";
*/
            // 2. Create the sparse system
            m_system = gsSparseSystem<T>(map);

        }
        else
            gsInfo << "Not implemented for patches > 3 \n";
        // END
    }

    template <class T, class bhVisitor>
    void gsBiharmonicArgyrisDirectAssembler<T,bhVisitor>::assemble()
    {
        GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

        // Reserve sparse system
        const index_t nz = gsAssemblerOptions::numColNz(m_bases[0][0],2,1,0.333333);
        m_system.reserve(nz, this->pde().numRhs());

        // Compute the Dirichlet Degrees of freedom (if needed by m_options)
        if (m_ppde.bcSecondKind().laplaceSides().size() != 0)
            Base::computeDirichletDofs();
        else
            computeDirichletAndNeumannDofs();

        m_ddof[0].setZero(); // Change TODO

        // Assemble volume integrals
        push();

        // Neuman conditions of first kind
        Base::template push<gsVisitorNeumann<T> >(
                m_ppde.bcFirstKind().neumannSides() );


        // Laplace conditions of second kind
        gsG1MultiBasis<T> g1MultiBasis(m_ppde.domain(), m_bases[0]); // Maybe earlier? Only need for #patches>2
        Base::template push<gsVisitorLaplaceBoundaryBiharmonic<T> >(
                m_ppde.bcSecondKind().laplaceSides(), g1MultiBasis );


        if ( m_options.getInt("InterfaceStrategy") == iFace::dg )
            gsWarn <<"DG option ignored.\n";

        /*
        // If requested, force Dirichlet boundary conditions by Nitsche's method
        this->template push<gsVisitorNitscheBiharmonic<T> >(
        m_ppde.bcSecondKind().dirichletSides() );
        */

        // Assembly is done, compress the matrix
        Base::finalize();
    }

// Pascal
    template <class T, class bhVisitor>
    void gsBiharmonicArgyrisDirectAssembler<T,bhVisitor>::push()
    {
        for (size_t np = 0; np < m_ppde.domain().nPatches(); ++np)
        {
            bhVisitor visitor(*m_pde_ptr);
            //Assemble (fill m_matrix and m_rhs) on patch np
            apply(visitor, np);
        }
    }

// Pascal
    template <class T, class bhVisitor>
    void gsBiharmonicArgyrisDirectAssembler<T,bhVisitor>::apply(bhVisitor & visitor,
                                                   size_t patchIndex,
                                                   boxSide side)
    {
        //gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";

        const gsBasisRefs<T> bases(m_bases, patchIndex);

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

            gsG1MultiBasis<T> g1MultiBasis(m_ppde.domain(), m_bases[0]); // Maybe earlier? Only need for #patches>2

            // Initialize reference quadrature rule and visitor data
            visitor_.initialize(bases, patchIndex, m_options, quRule);

            const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIndex];

            // Initialize domain element iterator -- using unknown 0
            typename gsBasis<T>::domainIter domIt = bases[0].makeDomainIterator(side);

            // Start iteration over elements
#ifdef _OPENMP
            for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#else
            for (; domIt->good(); domIt->next() )
#endif
            {
                // Map the Quadrature rule to the element
                quRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

                // Perform required evaluations on the quadrature nodes
                visitor_.evaluate(bases, g1MultiBasis, patch, quNodes);

                // Assemble on element
                visitor_.assemble(*domIt, quWeights);

                // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
                visitor_.localToGlobal(patchIndex, m_ddof, m_system); // omp_locks inside
            }
        }//omp parallel

    }


    // Pascal
    template <class T, class bhVisitor>
    void gsBiharmonicArgyrisDirectAssembler<T,bhVisitor>::constructSolution(const gsMatrix<T>& solVector,
                                                               gsMultiPatch<T>& result, short_t unk)
    {
        // we might need to get a result even without having the system ..
        //GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

        // fixme: based on \a m_options and \a unk choose the right dof mapper
        const gsDofMapper & mapper = m_system.colMapper(unk);

        result.clear(); // result is cleared first

        /*
        GISMO_ASSERT(solVector.rows() == m_dofs,
                     "The provided solution vector does not match the system."
                     " Expected: "<<mapper.freeSize()<<", Got:"<<solVector.rows() );
        */
        // is the situation whtn solVector has more than one columns important?
        const index_t dim = ( 0!=solVector.cols() ? solVector.cols() :  m_ddof[unk].cols() );

        // to do: test unknown_dim == dim

        gsMatrix<T> coeffs;
        for (size_t p = 0; p < m_pde_ptr->domain().nPatches(); ++p)
        {
            // Reconstruct solution coefficients on patch p
            const size_t sz  = m_bases[m_system.colBasis(unk)][p].size();
            coeffs.resize(sz, dim);

            for (size_t i = 0; i < sz; ++i)
            {
                if ( mapper.is_free(i, p) ) // DoF value is in the solVector
                {
                    coeffs.row(i) = solVector.row( mapper.index(i, p) );
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    //coeffs.row(i) = m_ddof[unk].row( mapper.bindex(i, p) ).head(dim);
                    gsMatrix<> zero(1,1);
                    zero.setZero();
                    coeffs.row(i) = zero; // TODO
                }
            }

            result.addPatch( m_bases[m_system.colBasis(unk)][p].makeGeometry( give(coeffs) ) );
        }

        // AM: result topology ?
    }

    template <class T, class bhVisitor>
    void gsBiharmonicArgyrisDirectAssembler<T,bhVisitor>::constructG1Solution(const gsMatrix<T>& solVector,
                                                                 gsMatrix<T>& result, short_t unk)
    {
        // we might need to get a result even without having the system ..
        //GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

        // fixme: based on \a m_options and \a unk choose the right dof mapper
        const gsDofMapper & mapper = m_system.colMapper(unk);

        gsG1MultiBasis<T> g1MultiBasis(m_ppde.domain(), m_bases[0]); // Maybe earlier? Only need for #patches>2
        result.setZero(g1MultiBasis.nBasisFunctions(),1); // result is cleared first

        for (size_t p = 0; p < m_pde_ptr->domain().interfaces().size(); ++p) // TODO == 1
        {
            for (index_t i = 0; i < result.rows(); ++i)
            {
                if ( mapper.is_free(i+m_bases[0].basis(0).size(), p) ) // DoF value is in the solVector
                {
                    result.row(i) = solVector.row( mapper.index(i+m_bases[0].basis(0).size(), p) );
                }
                else // eliminated DoF: fill with Dirichlet data
                {
                    //coeffs.row(i) = m_ddof[unk].row( mapper.bindex(i, p) ).head(dim);
                    gsMatrix<> zero(1,1);
                    zero.setZero();
                    result.row(i) = zero; // TODO
                }
            }
        }

        // AM: result topology ?
    }

// Pascal
    template <class T, class bhVisitor>
    void gsBiharmonicArgyrisDirectAssembler<T,bhVisitor>::computeDirichletAndNeumannDofs()
    {

        const gsMultiBasis<T> & mbasis = m_bases[m_system.colBasis(0)];
        const gsDofMapper & mapper =
                dirichlet::elimination == m_options.getInt("DirichletStrategy") ?
                m_system.colMapper(0) :
                mbasis.getMapper(dirichlet::elimination,
                                 static_cast<iFace::strategy>(m_options.getInt("InterfaceStrategy")),
                                 m_pde_ptr->bc(), 0);

        if(m_ddof.size()==0)
            m_ddof.resize(m_system.numUnknowns());

        m_ddof[0].resize( mapper.boundarySize(), m_system.unkSize(0)*m_system.rhs().cols());  //m_pde_ptr->numRhs() );

        // Set up matrix, right-hand-side and solution vector/matrix for
        // the L2-projection
        gsSparseEntries<T> projMatEntries;
        gsMatrix<T>        globProjRhs;
        globProjRhs.setZero( mapper.boundarySize(), m_system.unkSize(0)*m_system.rhs().cols() );

        // Temporaries
        gsVector<T> quWeights;

        gsMatrix<T> rhsVals, rhsVals2;
        gsMatrix<index_t> globIdxAct;
        gsMatrix<T> basisVals, basisGrads, physBasisGrad;

        gsVector<T> unormal;

        real_t lambda = 1e-2; // TODO not so optimal :/ Don't like that idea

        gsMapData<T> md(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM | NEED_JACOBIAN);


        typename gsBoundaryConditions<T>::const_iterator
                iter_dir = m_pde_ptr->bc().dirichletBegin();

        for ( typename gsBoundaryConditions<T>::const_iterator
                      iter = m_ppde.bcSecondKind().neumannBegin();
              iter != m_ppde.bcSecondKind().neumannEnd(); ++iter )
        {
            if (iter->isHomogeneous() )
                continue;

            GISMO_ASSERT(iter_dir->function()->targetDim() == m_system.unkSize(0)*m_system.rhs().cols(),
                         "Given Dirichlet boundary function does not match problem dimension. "
                                 << iter_dir->function()->targetDim()<<" != "<<m_system.unkSize(0)<<"\n");

            const int unk = iter_dir->unknown();

            const int patchIdx   = iter_dir->patch();
            const gsBasis<T> & basis = (m_bases[unk])[patchIdx];

            GISMO_ASSERT(iter_dir->patch() == patchIdx && iter_dir->side().index() == iter->side().index(),
                         "Given Dirichlet boundary edge does not match the neumann edge."
                                 <<iter_dir->patch()<<" != "<<patchIdx<<" and "
                                 <<iter_dir->side().index()<<" != "<<iter->side().index()<<"\n");

            const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIdx];

            // Set up quadrature to degree+1 Gauss points per direction,
            // all lying on iter->side() except from the direction which
            // is NOT along the element

            gsGaussRule<T> bdQuRule(basis, 1.0, 1, iter->side().direction());

            // Create the iterator along the given part boundary.
            typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

            for(; bdryIter->good(); bdryIter->next() )
            {
                bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                                md.points, quWeights);

                //geoEval->evaluateAt( md.points );
                patch.computeMap(md);

                // the values of the boundary condition are stored
                // to rhsVals. Here, "rhs" refers to the right-hand-side
                // of the L2-projection, not of the PDE.
                rhsVals = iter_dir->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );
                rhsVals2 = iter->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );

                basis.eval_into( md.points, basisVals);
                basis.deriv_into( md.points, basisGrads);

                // Indices involved here:
                // --- Local index:
                // Index of the basis function/DOF on the patch.
                // Does not take into account any boundary or interface conditions.
                // --- Global Index:
                // Each DOF has a unique global index that runs over all patches.
                // This global index includes a re-ordering such that all eliminated
                // DOFs come at the end.
                // The global index also takes care of glued interface, i.e., corresponding
                // DOFs on different patches will have the same global index, if they are
                // glued together.
                // --- Boundary Index (actually, it's a "Dirichlet Boundary Index"):
                // The eliminated DOFs, which come last in the global indexing,
                // have their own numbering starting from zero.

                // Get the global indices (second line) of the local
                // active basis (first line) functions/DOFs:
                basis.active_into(md.points.col(0), globIdxAct );
                gsMatrix<index_t> localIdx = globIdxAct;
                mapper.localToGlobal( globIdxAct, patchIdx, globIdxAct);

                // Out of the active functions/DOFs on this element, collect all those
                // which correspond to a boundary DOF.
                // This is checked by calling mapper.is_boundary_index( global Index )

                // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
                // something like a "element-wise index"
                std::vector<index_t> eltBdryFcts;
                eltBdryFcts.reserve(mapper.boundarySize());
                for( index_t i=0; i < globIdxAct.rows(); i++)
                    if( mapper.is_boundary_index(globIdxAct(i,0)) )
                        eltBdryFcts.push_back( i );

                // Do the actual assembly:
                for (index_t k = 0; k < md.points.cols(); k++)
                {
                    // Compute the outer normal vector on the side
                    outerNormal(md, k, iter->side(), unormal);

                    // Multiply quadrature weight by the measure of normal
                    const T weight_k = quWeights[k] * md.measure(k);
                    unormal.normalize();

                    transformGradients(md, k, basisGrads, physBasisGrad);

                    // Only run through the active boundary functions on the element:
                    for (size_t i0 = 0; i0 < eltBdryFcts.size(); i0++)
                    {
                        // Each active boundary function/DOF in eltBdryFcts has...
                        // ...the above-mentioned "element-wise index"
                        const unsigned i = eltBdryFcts[i0];
                        // ...the boundary index.
                        const unsigned ii = mapper.global_to_bindex(globIdxAct(i));

                        for (size_t j0 = 0; j0 < eltBdryFcts.size(); j0++)
                        {
                            const unsigned j = eltBdryFcts[j0];
                            const unsigned jj = mapper.global_to_bindex(globIdxAct(j));

                            // Use the "element-wise index" to get the needed
                            // function value.
                            // Use the boundary index to put the value in the proper
                            // place in the global projection matrix.
                            projMatEntries.add(ii, jj, weight_k * (basisVals(i, k) * basisVals(j, k) + lambda *
                                                                                                       ((physBasisGrad.col(i).transpose() * unormal)(0, 0)
                                                                                                        * (physBasisGrad.col(j).transpose() * unormal)(0, 0))));
                        } // for j

                        globProjRhs.row(ii) += weight_k * (basisVals(i, k) * rhsVals.col(k).transpose() + lambda *
                                                                                                          (physBasisGrad.col(i).transpose() * unormal) * (rhsVals2.col(k).transpose() * unormal));

                    } // for i
                } // for k

            } // bdryIter
            iter_dir++;
        } // boundaryConditions-Iterator

        gsSparseMatrix<T> globProjMat( mapper.boundarySize(), mapper.boundarySize() );
        globProjMat.setFrom( projMatEntries );
        globProjMat.makeCompressed();

        // Solve the linear system:
        // The position in the solution vector already corresponds to the
        // numbering by the boundary index. Hence, we can simply take them
        // for the values of the eliminated Dirichlet DOFs.
        typename gsSparseSolver<T>::CGDiagonal solver;

        m_ddof[0] = solver.compute( globProjMat ).solve ( globProjRhs );
        //m_ddof[0].setZero();
    }
// End

    template <class T, class bhVisitor>
    void gsBiharmonicArgyrisDirectAssembler<T,bhVisitor>::plotParaview(gsField<> &solField_interior, gsMatrix<T> &result)
    {

        std::string fn = "G1Biharmonic";
        index_t npts = 10000;
        gsParaviewCollection collection2(fn);
        std::string fileName2;

        gsG1MultiBasis<T> g1MultiBasis(m_ppde.domain(), m_bases[0]); // Maybe earlier? Only need for #patches>2


        for ( size_t pp =0; pp < m_pde_ptr->domain().nPatches(); ++pp ) // Patches
        {
            fileName2 = fn + util::to_string(pp);
            //writeSinglePatchField( field, i, fileName, npts );

            const gsFunction<T> & geometry = m_pde_ptr->domain().patch(pp);
            const gsFunction<T> & parField = solField_interior.function(pp);

            const int n = geometry.targetDim();
            const int d = geometry.domainDim();

            gsMatrix<T> ab = geometry.support();
            gsVector<T> a = ab.col(0);
            gsVector<T> b = ab.col(1);

            gsVector<unsigned> np = uniformSampleCount(a, b, npts);
            gsMatrix<T> pts = gsPointGrid(a, b, np);

            gsMatrix<T> eval_geo = geometry.eval(pts);//pts
            gsMatrix<T> eval_field = solField_interior.isParametric() ? parField.eval(pts) : parField.eval(eval_geo);

            // Here add g1 basis
            //eval_field.setZero();
            std::vector<gsMatrix<T>> eval_g1;
            g1MultiBasis.evalAllDers_into(pts, 0, eval_g1, pp);

            gsDofMapper mapper = m_system.rowMapper(0);
            for (index_t i = 0; i < eval_g1[0].rows(); i++)
            {
                if (mapper.is_free(i+ m_bases[0].basis(pp).size(), pp))
                    eval_field += eval_g1[0].row(i) * result(
                            mapper.index(i + m_bases[0].basis(pp).size(), pp),0); // Only 1 dim solution!!!

                // TODO Add here bdy
            }

            if ( 3 - d > 0 )
            {
                np.conservativeResize(3);
                np.bottomRows(3-d).setOnes();
            }
            else if (d > 3)
            {
                gsWarn<< "Cannot plot 4D data.\n";
                return;
            }

            if ( 3 - n > 0 )
            {
                eval_geo.conservativeResize(3,eval_geo.cols() );
                eval_geo.bottomRows(3-n).setZero();
            }
            else if (n > 3)
            {
                gsWarn<< "Data is more than 3 dimensions.\n";
            }

            if ( eval_field.rows() == 2)
            {
                eval_field.conservativeResize(3,eval_geo.cols() );
                eval_field.bottomRows(1).setZero(); // 3-field.dim()
            }

            gsWriteParaviewTPgrid(eval_geo, eval_field, np.template cast<index_t>(), fileName2);


            collection2.addPart(fileName2, ".vts");
        }
        collection2.save();
    }


} // namespace gismo