/** @file gsG1BiharmonicAssembler.h

    @brief Provides assembler for a homogenius Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/


#pragma once

#include <gsAssembler/gsAssembler.h>

#include <gsPde/gsBiharmonicPde.h>

#include <gsAssembler/gsVisitorG1Biharmonic.h>
#include <gsAssembler/gsVisitorNeumann.h>
#include <gsAssembler/gsVisitorLaplaceBiharmonic.h>
//#include <gsAssembler/gsVisitorNitscheBiharmonic.h>

#include <gsArgyris/gsC1ArgyrisBasis.h>

namespace gismo
{

/** @brief
    Implementation of a homogeneous Biharmonic Assembler.

    It sets up an assembler and assembles the system patch wise and
    combines the patch-local stiffness matrices into a global system.
    Dirichlet boundary can only be enforced strongly (i.e Nitsche is
    not implemented).
*/
template <class T, class bhVisitor = gsVisitorG1Biharmonic<T> >
class gsG1BiharmonicAssembler : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;
    typedef typename gsBoundaryConditions<T>::bcContainer bcContainer;

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
    gsG1BiharmonicAssembler( gsMultiPatch<T> const         & patches,
                           gsMappedBasis<2,T> const         & bases,
                           gsBoundaryConditions<T> const & bconditions,
                           gsBoundaryConditions<T> const & bconditions2,
                           const gsFunction<T>           & rhs,
                           const bool                    & twoPatch,
                           dirichlet::strategy           dirStrategy = dirichlet::none,
                           iFace::strategy               intStrategy = iFace::none)
    : m_ppde(patches,bconditions,bconditions2,rhs), m_bases(bases), m_twoPatch(twoPatch)
    {
        m_options.setInt("DirichletStrategy", dirStrategy);
        m_options.setInt("InterfaceStrategy", intStrategy);

        typename gsPde<T>::Ptr _pde = memory::make_shared_not_owned(&m_ppde);
        m_pde_ptr = _pde;

        refresh();
    }

    void refresh();
    
    void assemble();
    void push();
    void apply(bhVisitor & visitor,
               size_t patchIndex,
               boxSide side = boundary::none);

    void push(const bcContainer & bc);
    void apply(gsVisitorLaplaceBiharmonic<T> & visitor,
               size_t patchIndex,
               boxSide side = boundary::none);

    void constructSolution(const gsMatrix<T>& solVector,
                           gsMatrix<T>& result, short_t unk = 0) const;

    void computeDirichletDofsL2Proj(const short_t unk_ = 0);

protected:

    // fixme: add constructor and remove this
    gsBiharmonicPde<T> m_ppde;

    // G1 Basis
    gsMappedBasis<2,T> const m_bases;

    bool m_twoPatch;

    // Members from gsAssembler
    using Base::m_pde_ptr;
    using Base::m_ddof;
    using Base::m_options;
    using Base::m_system;
};

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::refresh()
{
    // We use predefined helper which initializes the system matrix
    // rows and columns using the same test and trial space
    //Base::scalarProblemGalerkinRefresh();
    gsVector<index_t> sz(m_bases.nPatches());
    for (size_t np = 0; np < m_bases.nPatches(); np++)
    {
        const gsC1ArgyrisBasis<2,T> c1ArgyrisBasis = static_cast<const gsC1ArgyrisBasis<2,T>&>(m_bases.getBase(np));
        sz[np] = c1ArgyrisBasis.size_rows();
    }

    gsDofMapper map(sz);

    typedef std::vector< patchSide >::const_iterator b_const_iter;
    for(b_const_iter iter = m_ppde.domain().bBegin();iter!=m_ppde.domain().bEnd();++iter)
    {
        gsMatrix<index_t> boundaryDofs;
        boundaryDofs = m_bases.getBase(iter->patch).boundaryOffset(*iter, 0);

        map.markBoundary(iter->patch, boundaryDofs);
    }

    typedef std::vector<boundaryInterface>::const_iterator i_const_iter;
    for(i_const_iter iter = m_ppde.domain().iBegin();iter!=m_ppde.domain().iEnd();++iter)
    {
        index_t numIntBdy = 0;
        std::vector<bool> kink;
        kink.resize(2);
        kink[0] = false; // u = 0
        kink[1] = false; // u = 1

        gsMatrix<> points_1, points_2, matrix_det(2,2);
        points_1.setZero(2,2);
        points_2.setZero(2,2);

        switch(iter->first().side().index())
        {
            case 1:
                points_1(1,1) = 1;
                break;
            case 2:
                points_1.setOnes();
                points_1(1,0) = 0;
                break;
            case 3:
                points_1(0,1) = 1;
                break;
            case 4:
                points_1.setOnes();
                points_1(0,0) = 0;
                break;
            default:
                gsInfo << "Something went wrong! \n";
                break;
        }

        switch(iter->second().side().index())
        {
            case 1:
                points_2(1,1) = 1;
                break;
            case 2:
                points_2.setOnes();
                points_2(1,0) = 0;
                break;
            case 3:
                points_2(0,1) = 1;
                break;
            case 4:
                points_2.setOnes();
                points_2(0,0) = 0;
                break;
            default:
                gsInfo << "Something went wrong! \n";
                break;
        }



        gsMatrix<index_t> matched_dofs1, matched_dofs2;
        matched_dofs1 = m_bases.getBase(iter->first().patch).boundaryOffset(iter->first().side().index(), 0);
        matched_dofs2 = m_bases.getBase(iter->second().patch).boundaryOffset(iter->second().side().index(), 0);

        if (m_twoPatch)
        {
            const index_t dir_1 = iter->first().side().index() > 2 ? 0 : 1;
            const index_t dir_2 = iter->second().side().index() > 2 ? 0 : 1;

            matrix_det.col(0) = m_ppde.domain().patch(iter->first().patch).jacobian(points_1.col(0)).col(1-dir_1);
            matrix_det.col(1) = m_ppde.domain().patch(iter->second().patch).jacobian(points_2.col(0)).col(1-dir_2);
            if (matrix_det.determinant()*matrix_det.determinant() > 1e-25)
            {
                kink[0] = true;
                numIntBdy += 1;
            }

            matrix_det.col(0) = m_ppde.domain().patch(iter->first().patch).jacobian(points_1.col(1)).col(1-dir_1);
            matrix_det.col(1) = m_ppde.domain().patch(iter->second().patch).jacobian(points_2.col(1)).col(1-dir_2);
            if (matrix_det.determinant()*matrix_det.determinant() > 1e-25)
            {
                kink[1] = true;
                numIntBdy += 1;
            }

            gsMatrix<index_t> bdy_interface(2+numIntBdy,1);
            bdy_interface.row(0) = matched_dofs1.row(0);
            bdy_interface.row(1) = matched_dofs1.bottomRows(1);
            if(kink[0])
                bdy_interface.row(2) = matched_dofs1.row(1);
            if(kink[1])
                bdy_interface(1+numIntBdy,0) = matched_dofs1(-2,0);

            map.markBoundary(iter->first(), bdy_interface); // TODO Shift to vertex
        }



        map.matchDofs(iter->first().patch, matched_dofs1, iter->second().patch, matched_dofs2);

        matched_dofs1 = m_bases.getBase(iter->first().patch).boundaryOffset(iter->first().side().index(), 1);
        matched_dofs2 = m_bases.getBase(iter->second().patch).boundaryOffset(iter->second().side().index(), 1);

        if (m_twoPatch)
        {
            gsMatrix<index_t> bdy_interface(2,1);
            bdy_interface.setZero(2,1);
            bdy_interface.row(0) = matched_dofs1.row(0);
            bdy_interface.row(1) = matched_dofs1.bottomRows(1);


            map.markBoundary(iter->first(), bdy_interface); // TODO Shift to vertex
        }

        map.matchDofs(iter->first().patch, matched_dofs1, iter->second().patch, matched_dofs2);

    }

    for (size_t numVer = 0; numVer < m_ppde.domain().vertices().size(); numVer++)
    {
        std::vector<patchCorner> allcornerLists = m_ppde.domain().vertices()[numVer];
        std::vector<size_t> patchIndex;
        std::vector<size_t> vertIndex;
        for (size_t j = 0; j < allcornerLists.size(); j++)
        {
            patchIndex.push_back(allcornerLists[j].patch);
            vertIndex.push_back(allcornerLists[j].m_index);
        }

        std::vector<gsMatrix<index_t>> bdy_cornerContainer;
        for (size_t i = 0; i < patchIndex.size(); ++i)
        {

            gsMatrix<index_t> bdy_corner;
            bdy_corner = m_bases.getBase(patchIndex[i]).boundaryOffset(vertIndex[i]+4, 0); // + 4 of the edge index // 0 == bdy
            if (bdy_corner(0,0) != -1) // Just for the two Patch case
                map.markBoundary(patchIndex[i], bdy_corner);

            bdy_corner = m_bases.getBase(patchIndex[i]).boundaryOffset(vertIndex[i]+4, 1); // + 4 of the edge index // 1 == coupled

            if (bdy_corner(0,0) != -1) // Just for the two Patch case
                bdy_cornerContainer.push_back(bdy_corner);
        }
        if (bdy_cornerContainer.size() > 1)
            for (size_t i = 1; i < bdy_cornerContainer.size(); ++i)
            {
                gsInfo << "matched: " << bdy_cornerContainer[0] << " with " << bdy_cornerContainer[i] <<"\n";
                map.matchDofs(patchIndex[i], bdy_cornerContainer[i], patchIndex[0], bdy_cornerContainer[0]);
            }



    }


    map.finalize();
    map.print();

    gsInfo << map.asVector() << "\n";

    // 2. Create the sparse system
    m_system = gsSparseSystem<T>(map);
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::assemble()
{
    GISMO_ASSERT(m_system.initialized(), "Sparse system is not initialized, call refresh()");

    // Reserve sparse system
    index_t nz = 1;
    for (short_t i = 0; i != m_bases.domainDim(); ++i) // 2 == dim
        nz *= static_cast<index_t>(2 * m_bases.degree(0,i) + 1 + 0.5); // Patch 0
    nz = static_cast<index_t>(nz*(1.333333));
    m_system.reserve(nz, 1);

    // Compute the Dirichlet Degrees of freedom (if needed by m_options)
    //Base::computeDirichletDofs();
    //m_ddof.resize(m_system.numUnknowns());
    //m_ddof[0].setZero(m_system.colMapper(0).boundarySize(), m_system.unkSize(0) * m_system.rhs().cols());
    computeDirichletDofsL2Proj();

    // Assemble volume integrals
    push();
    
    // Neuman conditions of first kind
    Base::template push<gsVisitorNeumann<T> >(
        m_ppde.bcFirstKind().neumannSides() );
    
    // Laplace conditions of second kind
    push( m_ppde.bcSecondKind().laplaceSides() );

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


template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::push()
{
    for (size_t np = 0; np < m_pde_ptr->domain().nPatches(); ++np)
    {
        bhVisitor visitor(*m_pde_ptr);
        //Assemble (fill m_matrix and m_rhs) on patch np
        apply(visitor, np);
    }
}

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::push(const bcContainer & BCs)
{
    for (typename bcContainer::const_iterator it
         = BCs.begin(); it!= BCs.end(); ++it)
    {
        gsVisitorLaplaceBiharmonic<T> visitor(*m_pde_ptr, *it);
        //Assemble (fill m_matrix and m_rhs) contribution from this BC
        apply(visitor, it->patch(), it->side());
    }
}


template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::apply(bhVisitor & visitor, size_t patchIndex, boxSide side)
{
    //gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";

    //const gsBasisRefs<T> bases(m_bases.getBase(patchIndex), patchIndex);

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

    // Initialize reference quadrature rule and visitor data
    visitor_.initialize(m_bases, patchIndex, m_options, quRule);

    const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIndex];

    // Initialize domain element iterator -- using unknown 0
    typename gsBasis<T>::domainIter domIt = m_bases.getBase(patchIndex).makeDomainIterator(side);

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
        visitor_.evaluate(m_bases, patch, quNodes);

        // Assemble on element
        visitor_.assemble(*domIt, quWeights);

        // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
        visitor_.localToGlobal(patchIndex, m_ddof, m_system); // omp_locks inside
    }
}//omp parallel

} // apply

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::apply(gsVisitorLaplaceBiharmonic<T> & visitor, size_t patchIndex, boxSide side)
{
    //gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";

    //const gsBasisRefs<T> bases(m_bases.getBase(patchIndex), patchIndex);

#pragma omp parallel
{
    gsQuadRule<T> quRule ; // Quadrature rule
    gsMatrix<T> quNodes  ; // Temp variable for mapped nodes
    gsVector<T> quWeights; // Temp variable for mapped weights

    gsVisitorLaplaceBiharmonic<T>
#ifdef _OPENMP
    // Create thread-private visitor
    visitor_(visitor);
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#else
            &visitor_ = visitor;
#endif

    // Initialize reference quadrature rule and visitor data
    visitor_.initialize(m_bases, patchIndex, m_options, quRule);

    const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIndex];

    // Initialize domain element iterator -- using unknown 0
    typename gsBasis<T>::domainIter domIt = m_bases.getBase(patchIndex).makeDomainIterator(side);

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
        visitor_.evaluate(m_bases, patch, quNodes);

        // Assemble on element
        visitor_.assemble(*domIt, quWeights);

        // Push to global matrix and right-hand side vector
#pragma omp critical(localToGlobal)
        visitor_.localToGlobal(patchIndex, m_ddof, m_system); // omp_locks inside
    }
}//omp parallel

} // apply


template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::constructSolution(const gsMatrix<T>& solVector,
                                       gsMatrix<T>& result, short_t unk) const
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

    result.resize(0,1);
    gsMatrix<T> coeffs;
    for (size_t np = 0; np < m_pde_ptr->domain().nPatches(); ++np)
    {
        // Reconstruct solution coefficients on patch p
        const gsC1ArgyrisBasis<2,T> c1ArgyrisBasis = static_cast<const gsC1ArgyrisBasis<2,T>&>(m_bases.getBase(np));
        index_t sz  = c1ArgyrisBasis.size_rows();
        coeffs.resize(sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if ( mapper.is_free(i, np) ) // DoF value is in the solVector
            {
                coeffs.row(i) = solVector.row( mapper.index(i, np) );
                //coeffs.row(i) = 0 * solVector.row( mapper.index(i, np) );
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                coeffs.row(i) = m_ddof[unk].row( mapper.bindex(i, np) ).head(dim);
            }
        }

        result.conservativeResize(result.rows()+coeffs.rows(), 1 );
        result.bottomRows(coeffs.rows()) = coeffs;
    }

    // AM: result topology ?
} // constructSolution

template <class T, class bhVisitor>
void gsG1BiharmonicAssembler<T,bhVisitor>::computeDirichletDofsL2Proj(const short_t unk_)
{
    const gsDofMapper & mapper = m_system.colMapper(unk_);
    m_ddof.resize(m_system.numUnknowns());
    m_ddof[unk_].resize( mapper.boundarySize(), m_system.unkSize(unk_)* m_pde_ptr->numRhs());

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>        globProjRhs;
    globProjRhs.setZero( mapper.boundarySize(), m_system.unkSize(unk_)* m_pde_ptr->numRhs() );

    // Temporaries
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals;
    gsMatrix<index_t> globIdxAct;
    gsMatrix<T> basisVals;

    gsMapData<T> md(NEED_MEASURE);

    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
          iter = m_pde_ptr->bc().dirichletBegin();
          iter != m_pde_ptr->bc().dirichletEnd(); ++iter )
    {
        if (iter->isHomogeneous() )
            continue;

        GISMO_ASSERT(iter->function()->targetDim() == m_system.unkSize(unk_)* m_pde_ptr->numRhs(),
                     "Given Dirichlet boundary function does not match problem dimension."
                     <<iter->function()->targetDim()<<" != "<<m_system.unkSize(unk_)<<"x"<<m_system.rhs().cols()<<"\n");

        const short_t unk = iter->unknown();
        if(unk!=unk_)
            continue;
        const index_t patchIdx   = iter->patch();

        const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIdx];

        // Set up quadrature to degree+1 Gauss points per direction,
        // all lying on iter->side() except from the direction which
        // is NOT along the element

        gsGaussRule<T> bdQuRule(m_bases.getBase(patchIdx), 1.0, 1, iter->side().direction());

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = m_bases.getBase(patchIdx).makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                            md.points, quWeights);

            //geoEval->evaluateAt( md.points );
            patch.computeMap(md);

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = iter->function()->eval( m_pde_ptr->domain()[patchIdx].eval( md.points ) );

            //basis.eval_into( md.points, basisVals);
            std::vector<gsMatrix<T>> basisData;
            m_bases.evalAllDers_into(patchIdx, md.points, 0, basisData);
            basisVals = basisData[0];

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
            m_bases.active_into(patchIdx, md.points.col(0), globIdxAct );

            for (index_t i = 0; i < patchIdx; i++)
                globIdxAct.array() -= mapper.patchSize(i);

            mapper.localToGlobal( globIdxAct, patchIdx, globIdxAct);

            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper.is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/globIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            eltBdryFcts.reserve(mapper.boundarySize());
            for( index_t i=0; i < globIdxAct.rows(); i++)
                if( mapper.is_boundary_index( globIdxAct(i,0)) )
                    eltBdryFcts.push_back( i );

            // Do the actual assembly:
            for( index_t k=0; k < md.points.cols(); k++ )
            {
                const T weight_k = quWeights[k] * md.measure(k);

                // Only run through the active boundary functions on the element:
                for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
                {
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    const unsigned i = eltBdryFcts[i0];
                    // ...the boundary index.
                    const unsigned ii = mapper.global_to_bindex( globIdxAct( i ));

                    for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                    {
                        const unsigned j = eltBdryFcts[j0];
                        const unsigned jj = mapper.global_to_bindex( globIdxAct( j ));

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        projMatEntries.add(ii, jj, weight_k * basisVals(i,k) * basisVals(j,k));
                    } // for j

                    globProjRhs.row(ii) += weight_k *  basisVals(i,k) * rhsVals.col(k).transpose();

                } // for i
            } // for k
        } // bdryIter
    } // boundaryConditions-Iterator

    gsSparseMatrix<T> globProjMat( mapper.boundarySize(), mapper.boundarySize() );
    globProjMat.setFrom( projMatEntries );
    globProjMat.makeCompressed();

    // Solve the linear system:
    // The position in the solution vector already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    typename gsSparseSolver<T>::CGDiagonal solver;
    m_ddof[unk_] = solver.compute( globProjMat ).solve ( globProjRhs );

    gsInfo << "rhs: " << globProjRhs << "\n";
    gsInfo << m_ddof[unk_] << "\n";

} // computeDirichletDofsL2Proj

} // namespace gismo



