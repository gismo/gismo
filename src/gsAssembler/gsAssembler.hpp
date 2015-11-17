/** @file gsAssembler.hpp

    @brief Provides assembler implementation for the Poisson equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss, A. Mantzaflaris, J. Sogn
*/

#include <gsAssembler/gsAssembler.h>
#include <gsAssembler/gsGaussRule.h>
#include <gsCore/gsMultiBasis.h>
#include <gsCore/gsDomainIterator.h>
#include <gsCore/gsField.h>
#include <gsUtils/gsPointGrid.h>

#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals and load vector


namespace gismo
{

template<class T>
void gsAssembler<T>::penalizeDirichlet()
{
    static const T PP = 1e9; // magic number
    
    // Store mapper
    gsDofMapper & elim = m_system.colMapper(0);
    gsDofMapper old    = elim;
    
    // Recompute dof mapper
    const bool conforming = ( m_options.intStrategy == iFace::conforming );
    m_bases.front().getMapper(conforming, m_pde_ptr->bc(), 0, elim); //unk

    // Compute dirichlet values
    computeDirichletDofs(); // fixme: add mapper

    // apply BCs
    for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_pde_ptr->bc().dirichletBegin();
          it != m_pde_ptr->bc().dirichletEnd(); ++it )
    {            
        const gsBasis<T> & basis = (m_bases[0])[it->patch()];

        gsMatrix<unsigned> bnd = safe( basis.boundary(it->side() ) ); 
        for (index_t k=0; k!= bnd.size(); ++k)
        {
            const index_t ii = old            .index ( bnd(k) , it->patch() );
            const index_t bb = elim.bindex( bnd(k) , it->patch() );

            m_system.matrix()(ii,ii) = PP;                  // corr. free indices
            m_system.rhs().row(ii)   = PP * m_ddof.row(bb); // boundary indices
        }
    }

    // restore mapper
    elim = old;
}

template<class T>
void gsAssembler<T>::assemble()
{GISMO_NO_IMPLEMENTATION}

template<class T>
void gsAssembler<T>::assemble(const gsMultiPatch<T> & curSolution)
{GISMO_NO_IMPLEMENTATION}


template<class T>
void gsAssembler<T>::setDirichletDofs(const gsMatrix<T> & vals, int unk)
{
    m_ddof = vals;
    // Assuming that the DoFs are already set by the user
    GISMO_ENSURE( m_ddof.rows() == m_system.colMapper(unk).boundarySize() && 
                  m_ddof.cols() == m_pde_ptr->numRhs(), 
                  "The Dirichlet DoFs were not provided correctly.");
}

template<class T>
void gsAssembler<T>::computeDirichletDofs() // index_t unk, mapper
{
    if ( m_options.dirStrategy != dirichlet::nitsche)
    {
        switch (m_options.dirValues)
        {
        case dirichlet::homogeneous:
            // If we have a homogeneous Dirichlet problem fill boundary
            // DoFs with zeros
            m_ddof.setZero( m_system.colMapper(0).boundarySize(), m_pde_ptr->numRhs() );
            break;
        case dirichlet::interpolation:        
            computeDirichletDofsIntpl();
            break;
        case dirichlet::l2Projection:
            computeDirichletDofsL2Proj();
            break;
        case dirichlet::user:
            // Assuming that the DoFs are already set by the user
            // fixme: more unknowns
            GISMO_ENSURE( m_ddof.rows() == m_system.colMapper(0).boundarySize() && 
                          m_ddof.cols() == m_pde_ptr->numRhs(), 
                          "The Dirichlet DoFs are not correctly provided.");
            break;
        default:
            GISMO_ERROR("Something went wrong with Dirichlet values.");
        }
    }

    // Corner values
    const gsDofMapper & mapper = m_system.colMapper(0);
    for ( typename gsBoundaryConditions<T>::const_citerator
              it = m_pde_ptr->bc().cornerBegin();
          it != m_pde_ptr->bc().cornerEnd(); ++it )
    {
        const int i  = m_bases[it->unknown][it->patch].functionAtCorner(it->corner);
        const int ii = mapper.bindex( i , it->patch );
        m_ddof.row(ii).setConstant(it->value);
    }    
}


// SKleiss: Note that this implementation is not useable for (T)HB-Splines!
//
// 1. Computation of the Dirichlet values explicitly uses gsPointGrid(rr)
// Computing a grid of evaluation points does not make sense for any locally
// refined basis.
// Also, "component(i)" is used.
// I'm afraid this makes sense ONLY FOR TENSOR-PRODUCT bases.
//
// 2. As of now (16.May 2014), the boundaryBasis of (T)HB-spline basis is not
// implemented, as far as I know.
//
// 3. gsInterpolate uses the anchors of the boundary basis.
// With truncated hierarchical B-splines, the use of classical anchors does
// not work, because functions might be truncated to zero at these points.
template<class T> // 
void gsAssembler<T>::computeDirichletDofsIntpl()
{
    const gsDofMapper & mapper = m_system.colMapper(0);
    
    m_ddof.resize( mapper.boundarySize(), m_pde_ptr->numRhs() );
    
    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_pde_ptr->bc().dirichletBegin();
          it != m_pde_ptr->bc().dirichletEnd(); ++it )
    {
        const int unk = it->unknown();
        const int k   = it->patch();
        const gsBasis<T> & basis = (m_bases[unk])[k];

        // Get dofs on this boundary
        gsMatrix<unsigned> * boundary = basis.boundary(it->side()) ;

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t i=0; i!= boundary->size(); ++i)
            {
                const int ii= mapper.bindex( (*boundary)(i) , k );
                m_ddof.row(ii).setZero();
            }
            delete boundary;
            continue;
        }

        // Get the side information
        int dir = it->side().direction( );
        index_t param = (it->side().parameter() ? 1 : 0);

        // Compute grid of points on the face ("face anchors")
        std::vector< gsVector<T> > rr;
        rr.reserve( this->patches().parDim() );

        for ( int i=0; i < this->patches().parDim(); ++i)
        {
            if ( i==dir )
            {
                gsVector<T> b(1); 
                b[0] = ( basis.component(i).support() ) (0, param);
                rr.push_back(b);
            }
            else
            {   
                rr.push_back( basis.component(i).anchors()->transpose() );
            }
        }

        GISMO_ASSERT(it->function()->targetDim() == m_pde_ptr->numRhs(),
                     "Given Dirichlet boundary function does not match problem dimension.\n");

        // Compute dirichlet values
        gsMatrix<T> fpts;
        if ( it->parametric() )
            fpts = it->function()->eval( gsPointGrid( rr ) );
        else
            fpts = it->function()->eval( m_pde_ptr->domain()[it->patch()].eval(  gsPointGrid( rr ) ) );

        // Interpolate dirichlet boundary 
        gsBasis<T> * h = basis.boundaryBasis(it->side());
        gsGeometry<T> * geo = h->interpolateAtAnchors(fpts);
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t l=0; l!= boundary->size(); ++l)
        {
            const int ii= mapper.bindex( (*boundary)(l) , it->patch() );

            m_ddof.row(ii) = dVals.row(l);
        }

        delete h;
        delete geo;
        delete boundary;
    }
}

template<class T>
void gsAssembler<T>::computeDirichletDofsL2Proj()
{
    const gsDofMapper & mapper = m_system.colMapper(0);

    m_ddof.resize( mapper.boundarySize(), m_pde_ptr->numRhs() );

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseEntries<T> projMatEntries;
    gsMatrix<T>       globProjRhs;   
    globProjRhs.setZero( mapper.boundarySize(), m_pde_ptr->numRhs() );

    // Temporaries
    gsMatrix<T> quNodes;
    gsVector<T> quWeights;

    gsMatrix<T> rhsVals;
    gsMatrix<unsigned> globIdxAct;
    gsMatrix<T> basisVals;

    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
              iter = m_pde_ptr->bc().dirichletBegin();
          iter != m_pde_ptr->bc().dirichletEnd(); ++iter )
    {
        const int unk = iter->unknown();
        const int patchIdx   = iter->patch();
        const gsBasis<T> & basis = (m_bases[unk])[patchIdx];

        typename gsGeometry<T>::Evaluator geoEval( m_pde_ptr->domain()[patchIdx].evaluator(NEED_MEASURE));

        // Set up quadrature. The number of integration points in the direction
        // that is NOT along the element has to be set to 1.
        gsVector<index_t> numQuNodes( basis.dim() );
        for( int i=0; i < basis.dim(); i++)
            numQuNodes[i] = (basis.degree(i)+1);
        numQuNodes[ iter->side().direction()] = 1;

        gsGaussRule<T> bdQuRule(numQuNodes);

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            bdQuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),
                          quNodes, quWeights);

            geoEval->evaluateAt( quNodes );

            // the values of the boundary condition are stored
            // to rhsVals. Here, "rhs" refers to the right-hand-side
            // of the L2-projection, not of the PDE.
            rhsVals = iter->function()->eval( m_pde_ptr->domain()[patchIdx].eval( quNodes ) );

            basis.eval_into( quNodes, basisVals);

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
            basis.active_into(quNodes.col(0), globIdxAct );
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
            for( index_t k=0; k < quNodes.cols(); k++ )
            {
                const T weight_k = quWeights[k] * geoEval->measure(k);

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
    m_ddof = solver.compute( globProjMat ).solve ( globProjRhs );
    
} // computeDirichletDofsL2Proj

template<class T>
void gsAssembler<T>::constructSolution(const gsMatrix<T>& solVector, 
                                       gsMultiPatch<T>& result, int unk) const
{
    // we might need to get a result even without having the system ..
    //GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    // fixme: based on \a m_options and \a unk choose the right dof mapper
    const gsDofMapper & mapper = m_system.colMapper(unk);

    GISMO_ASSERT(solVector.rows() == mapper.freeSize(), "Something went wrong, solution vector is not OK.");

    result.clear(); // result is cleared first
    
    const index_t dim = m_pde_ptr->numRhs();
    
    for (size_t p=0; p < m_pde_ptr->domain().nPatches(); ++p )
    {    
        // Reconstruct solution coefficients on patch p
        const int sz  = m_bases[unk][p].size();
        gsMatrix<T> coeffs( sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if ( mapper.is_free(i, p) ) // DoF value is in the solVector
            {
                coeffs.row(i) = solVector.row( mapper.index(i, p) );
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                coeffs.row(i) = m_ddof.row( mapper.bindex(i, p) );
            }
        }
        
        result.addPatch( m_bases[unk][p].makeGeometry( give(coeffs) ) );
    }

    // AM: result topology ?
}

template<class T>
gsField<T> *  gsAssembler<T>::constructSolution(const gsMatrix<T>& solVector,
                                                int unk) const
{
    const gsDofMapper & mapper = m_system.colMapper(unk);

    std::vector<gsFunction<T> * > sols ;

    const index_t dim = m_pde_ptr->numRhs();
    gsMatrix<T> coeffs;

    for (size_t p=0; p < m_pde_ptr->domain().nPatches(); ++p )
    {
        // Reconstruct solution coefficients on patch p
        const int sz  = m_bases[unk][p].size();
        coeffs.resize( sz, dim);

        for (index_t i = 0; i < sz; ++i)
        {
            if ( mapper.is_free(i, p) ) // DoF value is in the solVector
            {
                coeffs.row(i) = solVector.row( mapper.index(i, p) );
            }
            else // eliminated DoF: fill with Dirichlet data
            {
                coeffs.row(i) = m_ddof.row( mapper.bindex(i, p) );
            }
        }

        sols.push_back( m_bases[unk][p].makeGeometry( give(coeffs) ) );
    }

    return new gsField<T>(m_pde_ptr->domain(), sols);
}

}// namespace gismo
