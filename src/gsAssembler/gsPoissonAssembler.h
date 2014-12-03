/** @file gsPoissonAssembler.h

    @brief Provides assembler and solver for the Poisson equation, incl. adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsAssembler/gsAssemblerBase.h>
#include <gsAssembler/gsAssemblerOptions.h>

#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals
#include <gsAssembler/gsVisitorDg.h>      // Disc. Galerkin interface integrals
#include <gsAssembler/gsVisitorResidual.h>// Residual error estimator

#include <gsPde/gsPoissonPde.h>

//#include <gsAssembler/gsAssemblerUtils.h>

namespace gismo
{

/** @brief
    Implementation of an (multiple righ-hand side) Poisson solver.

    The Poisson equation: \f$-\Delta\mathbf{u}=\mathbf{f} \f$

    It sets up an assembler and assembles the system patch wise and combines
    the patch-local stiffness matrices into a global system by various methods
    (see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
    conditions in various ways (see gismo::gsDirichletStrategy).
*/
template <class T>
class gsPoissonAssembler : public gsAssemblerBase<T>
{
public:
    typedef gsAssemblerBase<T> Base;

public:

/** @brief
    Main Constructor of the assembler object.

    \param[in] pde A boundary value poisson problem
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces

    gsPoissonAssembler( gsPoissonPde<T> const         & pde,
                        gsMultiBasis<T> const         & bases,
                        gsDirichletStrategy           dirStrategy,
                        gsInterfaceStrategy           intStrategy)
    :  Base(patches), 
       m_rhsFun(&rhs),
       m_bConditions(bconditions),
       m_dirStrategy(dirStrategy), 
       m_intStrategy(intStrategy)
    {
        m_bases.push_back(bases);

        const bool conforming = ( m_intStrategy == glue );

        if ( m_dirStrategy == elimination || m_dirStrategy == homogeneous)
            m_dofMapper.push_back( bases.makeMapper(conforming, bconditions) );
        else
            m_dofMapper.push_back( bases.makeMapper(conforming) );

        m_dofs = m_dofMapper.front()->freeSize();

        // Resize system matrix and right hand side
        m_matrix.resize(m_dofs, m_dofs);
        m_rhs.resize(m_dofs, rhs.targetDim() );
    }
//*/


/** @brief
    Constructor of the assembler object.

    \param[in] patches is a gsMultiPatch object describing the geometry.
    \param[in] bases a multi-basis that contains patch-wise bases
    \param[in] bconditions is a gsBoundaryConditions object that holds all boundary conditions.
    \param[in] rhs is the right-hand side of the Poisson equation, \f$\mathbf{f}\f$.
    \param[in] dirStrategy option for the treatment of Dirichlet boundary
    \param[in] intStrategy option for the treatment of patch interfaces
*/
    gsPoissonAssembler( gsMultiPatch<T> const         & patches,
                        gsMultiBasis<T> const         & bases,
                        gsBoundaryConditions<T> const & bconditions,
                        const gsFunction<T>           & rhs,
                        dirichlet::strategy           dirStrategy,
                        iFace::strategy               intStrategy = iFace::glue)
    :  Base(patches), 
       m_rhsFun(&rhs),
       m_bConditions(bconditions),
       m_dirStrategy(dirichlet::none), 
       m_intStrategy(iFace::none)
    {
        m_bases.push_back(bases);

        gsAssemblerOptions options(dirStrategy, intStrategy);
        setOptions(options);
        
        m_dofs = m_dofMappers.front().freeSize();

    }

    virtual void setOptions(const gsAssemblerOptions  & options)
    {
        if ( m_dirStrategy != options.dirStrategy ||
             m_intStrategy != options.intStrategy )
        {
            m_dirStrategy = options.dirStrategy;
            m_intStrategy = options.intStrategy;
            
            const bool conforming = ( m_intStrategy == iFace::glue );
            m_dofMappers.resize(1);
            if ( m_dirStrategy == dirichlet::elimination || 
                 m_dirStrategy == dirichlet::homogeneous)
                m_bases.front().getMapper(conforming, m_bConditions, m_dofMappers.front() );
            else
                m_bases.front().getMapper(conforming, m_dofMappers.front() );
            
            m_dofs = m_dofMappers.front().freeSize();
        }
    }

    /// Main assembly routine
    virtual void assemble()
    {
        // If we have a homogeneous Dirichlet problem fill boundary
        // DoFs with zeros
        if ( m_dirStrategy == dirichlet::homogeneous)
            m_ddof.setZero( m_dofMappers[0].boundarySize(), m_rhsFun->targetDim() );

        // If the Dirichlet strategy is elimination then precompute
        // Dirichlet dofs (m_dofMapper excludes these from the system)
        if ( m_dirStrategy == dirichlet::elimination)
            computeDirichletDofsL2Proj(); // replacing computeDirichletDofs();
        

        if (m_dofs == 0 ) // Are there any interior dofs ?
        {
            gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
            return;
        }

        // Pre-allocate non-zero elements for each column of the
        // sparse matrix
        int nonZerosPerCol = 1;
        for (int i = 0; i < m_bases.front().dim(); ++i) // to do: improve
            nonZerosPerCol *= 2 * m_bases.front().maxDegree(i) + 1;

        m_matrix = gsSparseMatrix<T>(m_dofs, m_dofs); // Clean matrices
        m_matrix.reserve( gsVector<int>::Constant(m_dofs, nonZerosPerCol) );
        
        // Resize the load vector
        m_rhs.setZero(m_dofs, m_rhsFun->targetDim() );


        // Assemble volume stiffness and load vector integrals
        gsVisitorPoisson<T> poisson(*m_rhsFun);
        for (unsigned np=0; np < m_patches.nPatches(); ++np )
        {
            //Assemble stiffness matrix and rhs for the local patch
            // with index np and add to m_matrix and m_rhs
            this->apply(poisson, np);
        }

        // If requested, force Dirichlet boundary conditions by Nitsche's method
        if ( m_dirStrategy == dirichlet::nitsche )
            assembleNitsche();

        // Enforce Neumann boundary conditions
        assembleNeumann();

        // If we are in in dg (Discontinuous Galerkin) mode: add
        // interface contributions
        if ( m_intStrategy == iFace::dg )
            assembleDg();
        
        // Assembly is done, compress the matrix
        m_matrix.makeCompressed();   
    }
    
    
    virtual void assembleNitsche()
    {
        for ( typename gsBoundaryConditions<T>::const_iterator
                  it = m_bConditions.dirichletBegin();
              it != m_bConditions.dirichletEnd(); ++it )
        {
            gsVisitorNitsche<T> nitsche(*it->function(), penalty(it->patch()), it->side());
            
            // Note: it->unknown()
            this->apply(nitsche, it->patch(), it->side() );
        }
    }
    
    virtual void assembleNeumann()
    {
        for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_bConditions.neumannBegin();
              it != m_bConditions.neumannEnd(); ++it )
        {
            gsVisitorNeumann<T> neumann(*it->function(), it->side());

            // Note: it->unknown()
            this->apply(neumann, it->patch(), it->side() );
        }
    }
    
    virtual void assembleDg()
    {
        for ( typename gsMultiPatch<T>::iiterator it =
                  m_patches.iBegin(); it != m_patches.iEnd(); ++it )
        {
            if ( m_bases[0][(*it)[0].patch].numElements() < 
                 m_bases[0][(*it)[1].patch].numElements() )
                std::swap( (*it)[0], (*it)[1] );
            
            gsVisitorDg<T> dg(penalty(it->ps1.patch), it->ps1.side());
            this->apply(dg, *it);
        }
    }
    
//    void estimateResidue(const gsField<T> & sol, std::vector<gsMatrix<T> > & errors )
//    {
//        errors.resize( m_patches.nPatches() );
//
//        for (unsigned np=0; np < m_patches.nPatches(); ++np )
//        {
//            gsVisitorResidual<T> resEst(*m_rhsFun, sol);// can move out if not in parallel
//
//            int numEl = m_bases[0][np].numElements();
//
//            m_rhs.resize(numEl); // (!) element index, otherwise coords here !
//            // or the centerpoint of the element
//
//            //Assemble stiffness matrix and rhs for the local patch
//            // with index np and add to m_matrix and m_rhs
//            this->apply(resEst, np);
//
//            std::swap(errors[np], m_rhs);
//        }
//
//        /*
//          // Accumulate Neumann contribution
//        for ( typename gsBoundaryConditions<T>::const_iterator
//              it = m_bConditions.neumannBegin();
//              it != m_bConditions.neumannEnd(); ++it )
//        {
//            gsVisitorNeuResidual<T> neuResEst(*it->function(), sol, it->side() );
//
//            int numEl = m_bases[0][np].numElements();
//            m_rhs.setZero(numEl); // (!) element index, otherwise coords here !
//
//            // Note: it->unknown()
//            this->apply(neuResEst, it->patch(), it->side() );
//
//            errors[np] += m_rhs;
//        }
//        //*/
//
//    }


    /// Penalty constant for patch \a k, used for Nitsche and
    /// Discontinuous Galerkin methods
    T penalty(int k) const
    {
        //return gsAssemblerUtils<T>::getMu(m_bases[0][k]);
        const int deg = m_bases[0][k].maxDegree();
        return (deg + m_bases[0][k].dim()) * (deg + 1) * T(2.0);
    }

    /// Computes the Dirichlet DoF values by interpolation
    void computeDirichletDofs();

    void computeDirichletDofsL2Proj();

    /// Reconstruct solution field from computed solution vector
    gsField<T> * constructSolution(const gsMatrix<T> & solVector) const;

protected:

    // Right hand side function
    const gsFunction<T> * m_rhsFun;

    /// Boundary conditions
    gsBoundaryConditions<T> m_bConditions;

    // Strategy for dealing with Dirichlet dofs
    dirichlet::strategy m_dirStrategy;

    // Strategy for dealing with patch interface
    iFace::strategy m_intStrategy;

    // Members from gsAssemblerBase
    using gsAssemblerBase<T>::m_patches;
    using gsAssemblerBase<T>::m_bases;
    using gsAssemblerBase<T>::m_dofMappers;
    using gsAssemblerBase<T>::m_ddof;
    using gsAssemblerBase<T>::m_matrix;
    using gsAssemblerBase<T>::m_rhs;
    using gsAssemblerBase<T>::m_dofs;
};


//////////////////////////////////////////////////
//////////////////////////////////////////////////


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
void gsPoissonAssembler<T>::computeDirichletDofs()
{
    const gsDofMapper & mapper = m_dofMappers.front();
    
    m_ddof.resize( mapper.boundarySize(), m_rhsFun->targetDim() ); //--mrhs
    
    for ( typename gsBoundaryConditions<T>::const_iterator
              it = m_bConditions.dirichletBegin();
          it != m_bConditions.dirichletEnd(); ++it )
    {
        const int unk = it->unknown();
        const int k   = it->patch();
        const gsBasis<T> & basis = (m_bases[unk])[k];

        // Get dofs on this boundary
        gsMatrix<unsigned> * boundary = basis.boundary(it->side()) ;

        // If the condition is homogeneous then fill with zeros
        if ( it->isHomogeneous() )
        {
            for (index_t k=0; k!= boundary->size(); ++k)
            {
                const int ii= mapper.bindex( (*boundary)(k) , it->patch() );
                m_ddof.row(ii).setZero();
            }
            delete boundary;
            continue;
        }

        // Get the side information
        int dir = direction( it->side() );
        index_t param = (parameter( it->side() ) ? 1 : 0);

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

        // Compute dirichlet values
        gsMatrix<T> fpts = 
            it->function()->eval( m_patches[it->patch()].eval(  gsPointGrid( rr ) ) );

        // Interpolate dirichlet boundary 
        gsBasis<T> * h = basis.boundaryBasis(it->side());
        gsGeometry<T> * geo = gsInterpolate( *h, *h->anchors(), fpts );
        const gsMatrix<T> & dVals =  geo->coefs();

        // Save corresponding boundary dofs
        for (index_t k=0; k!= boundary->size(); ++k)
        {
            const int ii= mapper.bindex( (*boundary)(k) , it->patch() );
            m_ddof.row(ii) = dVals.row(k);
        }
        delete h;
        delete geo;
        delete boundary;
    }
}


// S.Kleiss
/** \brief Incorporates Dirichlet-boundary conditions by L2-projection.
 *
 * ...if the Dirichlet-strategy \em elimination is chosen.\n
 * A global \f$L_2\f$-projection is applied to the given Dirichlet data
 * and the eliminated coefficients are set to the corresponding values.
 * The projection is global in the sense that all Dirichlet-DOFs are
 * computed at once.
 */
template<class T>
void gsPoissonAssembler<T>::computeDirichletDofsL2Proj()
{
    const gsDofMapper & mapper = m_dofMappers.front();

    m_ddof.resize( mapper.boundarySize(), m_rhsFun->targetDim() );

    // Set up matrix, right-hand-side and solution vector/matrix for
    // the L2-projection
    gsSparseMatrix<T> globProjMat( mapper.boundarySize(), mapper.boundarySize() );
    gsMatrix<T> globProjRhs( mapper.boundarySize(), m_rhsFun->targetDim() );
    gsMatrix<T> globProjCoeffs( mapper.boundarySize(), m_rhsFun->targetDim() );

    globProjMat.setZero();
    globProjRhs.setZero();
    globProjCoeffs.setZero();

    // Iterate over all patch-sides with Dirichlet-boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
              iter = m_bConditions.dirichletBegin();
          iter != m_bConditions.dirichletEnd(); ++iter )
    {
        const int unk = iter->unknown();
        const int patchIdx   = iter->patch();
        const gsBasis<T> & basis = (m_bases[unk])[patchIdx];

        // Set up quadrature. The number of integration points in the direction
        // that is NOT along the element has to be set to 1.
        gsVector<index_t> numQuNodes( basis.dim() );
        for( int i=0; i < basis.dim(); i++)
            numQuNodes[i] = (basis.degree(i)+2);
        numQuNodes[ direction( iter->side() )] = 1;

        gsGaussRule<T> QuRule;
        QuRule.setNodes( numQuNodes);
        gsMatrix<T> quNodes;
        gsVector<T> quWeights;

        // Create the iterator along the given part boundary.
        typename gsBasis<T>::domainIter bdryIter = basis.makeDomainIterator(iter->side());

        for(; bdryIter->good(); bdryIter->next() )
        {
            //bdryIter->computeQuadratureRule( numQuNodes );
            QuRule.mapTo( bdryIter->lowerCorner(), bdryIter->upperCorner(),\
                             quNodes, quWeights);

            gsMatrix<T> rhsVals = iter->function()->eval( m_patches[patchIdx].eval( quNodes ) );

            gsMatrix<unsigned> locIdxAct;
            basis.active_into( quNodes, locIdxAct );

            gsMatrix<T> basisVals;
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


            // Get the global indices of the local active basis functions/DOFs:
            gsMatrix<unsigned> globIdxAct;
            mapper.localToGlobal( locIdxAct, patchIdx, globIdxAct);

            // Out of the active functions/DOFs on this element, collect all those
            // which correspond to a boundary DOF.
            // This is checked by calling mapper::is_boundary_index( global Index )

            // eltBdryFcts stores the row in basisVals/locIdxAct, i.e.,
            // something like a "element-wise index"
            std::vector<index_t> eltBdryFcts;
            for( index_t i=0; i < globIdxAct.rows(); i++)
                if( mapper.is_boundary_index( globIdxAct(i,0)) )
                    eltBdryFcts.push_back( i );

            // Do the actual assembly:
            for( index_t k=0; k < quNodes.cols(); k++ )
            {
                T weight_k = quWeights[k];

                // Only run through the active boundary functions on the element:
                for( size_t i0=0; i0 < eltBdryFcts.size(); i0++ )
                {
                    // Each active boundary function/DOF in eltBdryFcts has...
                    // ...the above-mentioned "element-wise index"
                    unsigned i = eltBdryFcts[i0];
                    // ...the boundary index.
                    unsigned ii = mapper.global_to_bindex( globIdxAct( i ));
                    for( size_t j0=0; j0 < eltBdryFcts.size(); j0++ )
                    {
                        unsigned j = eltBdryFcts[j0];
                        unsigned jj = mapper.global_to_bindex( globIdxAct( j ));

                        // Use the "element-wise index" to get the needed
                        // function value.
                        // Use the boundary index to put the value in the proper
                        // place in the global projection matrix.
                        globProjMat.coeffRef(ii,jj) += weight_k * basisVals(i,k) * basisVals(j,k);
                    } // for j

                    globProjRhs.row(ii) += weight_k *  basisVals(i,k) * rhsVals.col(k).transpose();

                } // for i
            } // for k
        } // bdryIter
    } // boundaryConditions-Iterator

    globProjMat.makeCompressed();

    // Solve the linear system
    Eigen::ConjugateGradient< gsSparseMatrix<T> > solver;
    globProjCoeffs = solver.compute( globProjMat ).solve ( globProjRhs );

    // The position in globProjCoeffs already corresponds to the
    // numbering by the boundary index. Hence, we can simply take them
    // for the values of the eliminated Dirichlet DOFs.
    m_ddof = globProjCoeffs;

} // computeDirichletDofsL2Proj


template<class T>
gsField<T> *  gsPoissonAssembler<T>::constructSolution(const gsMatrix<T>& solVector) const
//gsField<T> & result ) const
{
    GISMO_ASSERT(m_dofs == m_rhs.rows(), "Something went wrong, assemble() not called?");

    const gsDofMapper & mapper = m_dofMappers.front();

    std::vector<gsFunction<T> * > sols ;
    gsMatrix<T> coeffs;
    
    for (size_t p=0; p < m_patches.nPatches(); ++p )
    {    
        // Reconstruct solution coefficients on patch p
        const int sz  = m_bases[0][p].size();
        const index_t dim = m_rhsFun->targetDim();
        
        coeffs.resize( sz, dim );
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
        
        sols.push_back( m_bases[0][p].makeGeometry( give(coeffs) ) );
    }

    //result = gsField<T>(m_patches, sols);
    return new gsField<T>(m_patches, sols);
}

} // namespace gismo



