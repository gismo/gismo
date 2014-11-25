/** @file gsPoissonAssembler.hpp

    @brief Provides assembler and solver for the Poisson equation, incl. adaptive refinement.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss A. Mantzaflaris, J. Sogn
*/

#pragma once

#include <gsCore/gsDebug.h>
#include <gsAssembler/gsAssemblerUtils.h>
#include <gsAssembler/gsPoissonAssembler.h>

//#include <gismo.h>

namespace gismo {


//Virtual function
template<class T>
void gsPoissonAssembler<T>::initialize()
{
    //Tells the user that elimination strategy does not work in combination with dg.
    GISMO_ASSERT(!(m_interfaceStrategy == iFace::dg && 
                   m_dirStrategy == dirichlet::elimination),
                 "Must use Nitsche method for Dirichlet BC when using discontinuous Galerkin for patch interface");

    this->m_fixedDofs.resize(1); //One unknown u (either scalar or vector valued)

    // Initialize the DofMapper and finds m_idofs and m_idofs
    this->gsPdeAssembler<T>::initDofMapper(m_interfaceStrategy, m_dirStrategy, true);

    // Initialize the assembler by resize the global matrix and the right-hand side.
    this->m_matrix.resize(m_idofs,m_idofs);
    this->m_rhs.resize(m_idofs,this->m_unknownDim[0]);
    this->m_matrix.setZero();
    this->m_rhs.setZero();
}


// Assembles the final system with all boundary conditions contained
template<class T>
void gsPoissonAssembler<T>::assemble()
{
    // If the Dirichlet strategy is elimination then precompute
    // Dirichlet dofs (m_dofMapper excludes these from the system)
    if ( m_dirStrategy == dirichlet::elimination || 
         m_dirStrategy == dirichlet::homogeneous)
        this->computeDirichletDofs(m_dirStrategy);
    
    if (m_idofs == 0 ) // Are there any interior dofs ?
    {
        gsWarn << " No internal DOFs. Computed Dirichlet boundary only.\n" <<"\n" ;
        return;
    }
    
    // Assemble the system matrix and right-hand side
    if (this->m_patches.nPatches() == 1)
    {
        assemblePatch(0);
        this->m_matrix.makeCompressed();
    }
    
    else
        assembleMultipatch();
    
    // Enforce Dirichlet boundary conditions by Nitsche's method
    if ( m_dirStrategy == dirichlet::nitsche )
    {
        for ( typename gsBoundaryConditions<T>::const_iterator
                  it = this->m_bconditions[0]->dirichletBegin();
              it != this->m_bconditions[0]->dirichletEnd(); ++it )
        {
            // to do remove members from args
            applyBoundary(this->m_bases[0][it->patch()], *it, *this->m_dofMapper[0]);
        }
    }

    // Enforce Neumann boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
          it = this->m_bconditions[0]->neumannBegin();
          it != this->m_bconditions[0]->neumannEnd(); ++it )
    {
        applyBoundary(this->m_bases[0][it->patch()], *it, *this->m_dofMapper[0]);
    }
    

    // If we are in in dg (Discontinuous Galerkin) mode add interface
    // contributions
    if ( m_interfaceStrategy == iFace::dg )
    {
        for ( typename gsMultiPatch<T>::const_iiterator it =
                  m_patches.iBegin(); it != m_patches.iEnd(); ++it )

            applyDG( this->m_bases[0][it->ps1.patch], this->m_bases[0][it->ps2.patch],
                     m_patches.patch( it->ps1.patch ), m_patches.patch( it->ps2.patch ),
                     *it, *this->m_dofMapper[0], m_matrix );
    }
}


template<class T>
void gsPoissonAssembler<T>::assembleMultipatch()
{
    for (unsigned np=0; np < this->m_patches.nPatches(); ++np )
    {
         //Assemble stiffness matrix and rhs for the local patch
         // with index np and add to m_matrix and m_rhs
         assemblePatch(np);
    }
    this->m_matrix.makeCompressed();
}

template<class T>
void gsPoissonAssembler<T>::assemblePatch(int patchIndex)
{
    const gsBasis<T> & basis       = m_bases[0][patchIndex];
    const gsDofMapper & mapper  = *m_dofMapper[0];
    const gsMatrix<T> & ddof       =  m_fixedDofs[0];

    // Quadrature  nodes
    gsVector<int> numCwiseNodes = gsAssemblerUtils<T>::getNumIntNodesFor( basis );
    index_t       numNodes      = numCwiseNodes.prod();
    gsGaussRule<T> QuRule( numCwiseNodes ); // reference rule
    gsMatrix<T> quNodes  ; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights

    // Empty sparse matrix to fill with the contribution of the current patch
    gsSparseMatrix<T> sys_K(m_idofs, m_idofs );

    // Estimate max non-zeros per row (only valid for tensor product bases right now)
    unsigned nzRowsPerCol = 1;
    for (int i = 0; i < basis.dim(); ++i)
        nzRowsPerCol *= 2 * basis.degree(i) + 1;
    // Reserve space
    // \todo get the fragment of dofs in patch numbered patchIndex
    sys_K.reserve( gsVector<int>::Constant(m_idofs, nzRowsPerCol) );

    // Empty sparse rhs vector to fill with the contribution of the current patch
    gsMatrix<T> sys_b(mapper.freeSize(), this->m_unknownDim[0]);
    sys_b.setZero();

    // values of the right-hand side
    gsMatrix<T> rhsVals;

    // basisData contains stacked the values and the gradients of all
    // basis functions at one quadrature node
    // trf_grads_k keeps the (transformed) physical gradients of the
    // basis functions at one quadrature node
    gsMatrix<T> basisData, trf_grads_k;

    // active basis functions at one quadrature node
    gsMatrix<unsigned> actives; 
    
    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell
    gsMatrix<T> localRhs;       // local load vector
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
    this->m_patches.patch(patchIndex).evaluator(NEED_VALUE   |
                                                NEED_MEASURE |
                                                NEED_GRAD_TRANSFORM) );

    int targetDim = this->m_unknownDim[0];

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

        // Compute the active basis functions
        basis.active_into(domIt->centerPoint() , actives);
        const index_t numActive = actives.rows();

        // Local Dofs to global dofs
        mapper.localToGlobal(actives, patchIndex, actives);
 
        // Evaluate basis functions on element
        basis.evalAllDers_into( quNodes, 1, basisData);
        const typename gsMatrix<T>::Block basisValues = 
            basisData.topRows(numActive);
        const typename gsMatrix<T>::Block basisGrads  = 
            basisData.middleRows(numActive, numActive*basis.dim());

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(quNodes);

        // evaluate right-hand side at the geometry points
        if (m_rhs_function)
            m_rhs_function->eval_into( geoEval->values(), rhsVals );

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);
        localRhs.setZero(numActive, targetDim);

        for (index_t k = 0; k < numNodes; ++k)      // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = quWeights[k] * geoEval->measure(k);
                
            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval->transformGradients(k, basisGrads, trf_grads_k);
            
            //  if (m_rhs_function) 
            // assumption: if ,_rhs_function = NULL, then the Laplace problem is solved.
            // {
            localRhs += weight * ( basisValues.col(k) * rhsVals.col(k).transpose() ) ;
            //  }

            localStiffness.noalias() += weight * (trf_grads_k.transpose() * trf_grads_k);
        }  // end loop Gauss nodes

        // add contributions from local stiffness matrix to global
        // stiffness matrix and load vector

	    gsAssemblerUtils<T>::localToGlobal_withBC(localStiffness, localRhs, ddof, mapper,
				 actives.cast<int>(), sys_K, sys_b, m_isSymmetric);

    } //end loop over all domain elements

    this->m_matrix += sys_K;
    this->m_rhs    += sys_b;
}

// Solves the linear system and fills up \a m_sysSolution
template<class T>
void gsPoissonAssembler<T>::solveSystem()
{
    // Initialize solver
    Eigen::ConjugateGradient< gsSparseMatrix<T> > solver;
    // Solve linear system
    this->m_sysSolution = solver.compute(  this->m_matrix ).solve (  this->m_rhs );
    
    gsInfo << "residual error: " << solver.error() << "\n";
    gsInfo << "    iterations: " << solver.iterations() << "\n";
}


template<class T> void 
gsPoissonAssembler<T>::applyBoundary( const gsBasis<T>   & B,
				   const boundary_condition<T> & bc, 
                   const gsDofMapper& mapper)
{
    switch ( bc.type() )
        {
        case boundary::dirichlet:
            this->gsPdeAssembler<T>::boundaryNitsche(B, bc.patch(), *bc.function(), bc.side(), mapper);
            break;
        case boundary::neumann:
            this->gsPdeAssembler<T>::boundaryNeumann(B, bc.patch(), *bc.function(), bc.side(), mapper);
            break;
        default:
            gsWarn<<"Unknown boundary condition.\n";
        }
}

template<class T> void
gsPoissonAssembler<T>::applyDG( const gsBasis<T> & B1, const gsBasis<T> & B2,
                             const gsGeometry<T> & geo1,
                             const gsGeometry<T> & geo2,
                             const boundaryInterface & bi,
                             const gsDofMapper& mapper,
                             gsSparseMatrix<T> & system_matrix )
{
    gsSparseMatrix<T> & lhs    = system_matrix;
    const int d                = B1.dim();
    const T mu                 = gsAssemblerUtils<T>::getMu(B1);

    //gsDebug<<"Apply DG on "<< bi <<"(mu="<<mu<<").\n";

    const int patch1           = bi[0].patch;
    const int patch2           = bi[1].patch;
    const boundary::side side1 = bi[0].side;
    const boundary::side side2 = bi[1].side;

    //const int dir1             = direction(side1);
    //const int dir2             = direction(side2);
    //GISMO_ASSERT( B1.component(!dir1).size() == B2.component(!dir2).size(),
    //              "DG method not implemented yet for non matching interfaces");

    // Quadrature for boundary integrals
    gsVector<int> intNodes1 = 
        gsAssemblerUtils<T>::getNumIntNodesForInterface( B1, B2, bi, true  );
    gsVector<int> intNodes2 = 
        gsAssemblerUtils<T>::getNumIntNodesForInterface( B1, B2, bi, false );
    const index_t numNodes1 = intNodes1.prod();
    //const index_t numNodes2 = intNodes2.prod();

    // Quadrature  nodes
    gsGaussRule<T> QuRule1(intNodes1), QuRule2(intNodes2); // reference rules
    gsMatrix<T> quNodes1  , quNodes2  ; // Mapped nodes
    gsVector<T> quWeights1, quWeights2; // Mapped weights

    // Evaluators for the two patches
    typename gsGeometry<T>::Evaluator geoEval1( 
        geo1.evaluator(NEED_VALUE| NEED_JACOBIAN| NEED_GRAD_TRANSFORM) );
    typename gsGeometry<T>::Evaluator geoEval2( 
        geo2.evaluator(NEED_VALUE| NEED_JACOBIAN| NEED_GRAD_TRANSFORM) );
    
    // Temporaries
    gsMatrix<T> grads_k_1, grads_k_2, basisData1, basisData2;
    gsVector<T> unormal(d);
    // active basis functions at one quadrature node
    gsMatrix<unsigned> actives1, actives2;

    gsMatrix<T> B11, B22, B12, B21,
                E11, E22, E12, E21,
                N1 , N2;

    // iterator on grid cells on the "right"
    typename gsDomainIterator<T>::uPtr domIter2= B2.makeDomainIterator(side2);

    //int count = 0;
    // iterate over all boundary grid cells on the "left"
    for (typename gsDomainIterator<T>::uPtr domIter1 = B1.makeDomainIterator(side1);
         domIter1->good(); domIter1->next())
    {
        //count++;
        // Get the element of the other side in domIter2
        //domIter1->adjacent( bi.orient, *domIter2 );

        // Compute the quadrature rule on both sides
        QuRule1.mapTo( domIter1->lowerCorner(), domIter1->upperCorner(), quNodes1, quWeights1 );
        QuRule2.mapTo( domIter2->lowerCorner(), domIter2->upperCorner(), quNodes2, quWeights2 );

        // Compute the active basis functions
        B1.active_into(domIter1->centerPoint() , actives1);
        B2.active_into(domIter2->centerPoint() , actives2);
        const index_t numActive = actives1.rows(); // assuming numActive1=numActive2

        // Evaluate basis functions and their first derivatives
        B1.evalAllDers_into( quNodes1, 1, basisData1);
        const typename gsMatrix<T>::Block ev1 = basisData1.topRows(numActive);
        const typename gsMatrix<T>::Block grads1 = 
            basisData1.middleRows(numActive, numActive*d );
        B2.evalAllDers_into( quNodes2, 1, basisData2);
        const typename gsMatrix<T>::Block ev2 = basisData2.topRows(numActive);
        const typename gsMatrix<T>::Block grads2 = 
            basisData2.middleRows(numActive, numActive*d );

        // Local Dofs to global dofs
        mapper.localToGlobal(actives1, patch1, actives1);
        mapper.localToGlobal(actives2, patch2, actives2);

        // Push forward the quad-points to the physical domain
        geoEval1->evaluateAt(quNodes1);
        geoEval2->evaluateAt(quNodes2);

        B11.setZero(numActive, numActive); B22.setZero(numActive, numActive);
        B12.setZero(numActive, numActive); B21.setZero(numActive, numActive);
        E11.setZero(numActive, numActive); E22.setZero(numActive, numActive);
        E12.setZero(numActive, numActive); E21.setZero(numActive, numActive);

        // assuming quNodes1.cols() == quNodes2.cols()
        for (index_t k=0; k!= numNodes1; ++k)
        {
            // Compute the outer normal vector from patch1
            geoEval1->outerNormal(k, side1, unormal);

            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T fff = quWeights1[k] *  unormal.norm();

            // Compute the outer unit normal vector from patch1
            unormal.normalize();

            // Transform the basis gradients
            geoEval1->transformGradients(k, grads1, grads_k_1);
            geoEval2->transformGradients(k, grads2, grads_k_2);

            // Compute element matrices
            const gsMatrix<T> & val1 = ev1.col(k);
            const gsMatrix<T> & val2 = ev2.col(k);
            const T c1 = fff * T(0.5);
            N1.noalias() = unormal.transpose() * grads_k_1;
            N2.noalias() = unormal.transpose() * grads_k_2;
            B11.noalias() += c1 * ( val1 * N1 );
            B22.noalias() -= c1 * ( val2 * N2 );
            B12.noalias() += c1 * ( val1 * N2 );
            B21.noalias() -= c1 * ( val2 * N1 );
            const T c2 = fff * mu;
            E11.noalias() += c2 * ( val1 * val1.transpose() );
            E22.noalias() += c2 * ( val2 * val2.transpose() );
            E12.noalias() += c2 * ( val1 * val2.transpose() );
            E21.noalias() += c2 * ( val2 * val1.transpose() );
        }

        // Push element contributions to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const index_t jj1 = actives1(j); // N1_j
            const index_t jj2 = actives2(j); // N2_j
            for (index_t i=0; i!=numActive; ++i)
            {
                const index_t  ii1 = actives1(i); // N1_i
                const index_t  ii2 = actives2(i); // N2_i

                if ( !m_isSymmetric || jj1 <= ii1 )
                    lhs( ii1, jj1 ) -=  B11(i,j) + B11(j,i) - E11(i,j);
                if ( !m_isSymmetric || jj2 <= ii2 )
                    lhs( ii2, jj2 ) -=  B22(i,j) + B22(j,i) - E22(i,j);
                if ( !m_isSymmetric || jj2 <= ii1 )
                    lhs( ii1, jj2)  -=  B12(i,j) + B21(j,i) + E12(i,j);
                if ( !m_isSymmetric || jj1 <= ii2 )
                    lhs( ii2, jj1)  -=  B21(i,j) + B12(j,i) + E21(i,j);
            }
        }

        domIter2->next();
    }
}


// S.Kleiss
template <class T>
void gsPoissonAssembler<T>::boundaryFixDofs( const gsVector<unsigned> & Idx ,
                  const gsMatrix<T> & Val )
{

    // REMARK:
    // Probably less efficient than the implementation of
    // gsForceDirichletConditions, but it does the correction of
    // system-matrix and right-hand-side, as well as fixing the
    // corresponding value in the right-hand-side in one go.

    GISMO_ASSERT( Val.cols() == m_rhs.cols() , "Prescribed values must have same no. of components as RHS of PDE");
    GISMO_ASSERT( Val.rows() == Idx.size() , "Number of entries in Val and Idx do not match.");

    // Get a brutal "inverse mapping" that tells you
    // whether a given dof is one that should be fixed,
    // and if so, which entry of Idx/Val it corresponds to.
    //
    // Obviously, requires O(n) memory, where n is the
    // number of degrees of freedom, but
    // let's just call this cost "negligible" compared to the
    // global stiffness matrix.
    gsVector<int> Flag_Enf( m_matrix.rows() );
    Flag_Enf.setConstant(-1); // -1 corresponds to "do not fix"

    for( int i=0; i < Val.rows(); i++ )
    {
        // Make sure that the indices are within the range of
        // DOF-indices, i.e., that they do not exceed the size
        // of the system matrix.
        GISMO_ASSERT( int( Idx[i] ) < m_matrix.rows(), "Trying to fix value for DOF-index out of range.");

        // Set up the "inverse mapping".
        Flag_Enf[ Idx[i] ] = i;

        // Fix the desired value(s) in the right-hand-side.
        m_rhs.row( Idx[i] ) = Val.row( i );
    }

    int r;
    int c;

    // Treatment of the symmetric and non-symmetric case are almost the same.
    // However, instead of having an if-test in each loop, there is one
    // test in the beginning and some copy-pasted code.
    if( m_isSymmetric )
    {
        // Due to the symmetric structure where only the lower half of the
        // matrix is stored, entries which are really above the diagonal will have to
        // be accessed by using the corresponding ("transposed") below-diagonal entry.
        // Hence, it might not safe to set a value to zero, right away, because its value
        // might still be used later (as a virtual above-diagonal-value, so to say).
        // So, instead of setting any entry to zero right away, the indices are stored
        // and the values set to zero at the end.
        std::vector<int> setZero_r;
        std::vector<int> setZero_c;

        for (int k=0; k< m_matrix.outerSize(); ++k)
          for (typename gsSparseMatrix<T>::InnerIterator it(m_matrix,k); it; ++it)
          {
            r = it.row();
            c = it.col();

            GISMO_ASSERT( r >= c, "boundaryFixDofs: Assumes LOWER part of matrix is stored!");

            if( Flag_Enf[r] >= 0 )
            {
                // Row corresponds to a fix-dof.
                // Right-hand-side has been fixed above.
                if( r == c )
                {
                    // If row == column, we are on the diagonal,
                    // so put a 1 on the diagonal.
                    m_matrix.coeffRef(r,c) = T(1);
                }
                else
                {
                    // Otherwise, it is an off-diagonal entry.
                    // Since only the lower half of the matrix is
                    // stored, this entry corresponds to one entry
                    // of the upper half.
                    // Hence, unless the row of the above-diagonal-entry
                    // (which is the column of this below-diagonal-entry)
                    // corresponds to a fix-dof, this entry has to be
                    // corrected in the right-hand-side.
                    if( Flag_Enf[c] < 0 )
                        m_rhs.row(c) -= it.value() * Val.row( Flag_Enf[r] );

                    // Mark the element for setting it to zero.
                    setZero_r.push_back( r );
                    setZero_c.push_back( c );
                }
            }
            else if( Flag_Enf[c] >= 0 && r != c )
            {
                // The entry is off-diagonal and its column corresponds
                // to a fix-dof. Correct the right-hand-side.
                m_rhs.row(r) -= it.value() * Val.row( Flag_Enf[c] );

                // Mark the element for setting it to zero.
                setZero_r.push_back( r );
                setZero_c.push_back( c );
            }
            // else...
            //
            // otherwise, we are at an entry where the
            // row corresponds to a non-fix-dof and the
            // column also corresponds to a non-fix-dof.
            // In this case, do nothing.

          } // end for inner iterator

        for( int i=0; i < int( setZero_r.size() ); i++)
            m_matrix.coeffRef( setZero_r[i], setZero_c[i] ) = T(0);
    }
    else // m_isSymmetric == false, i.e., the full stiffness matrix is available
    {
        // Iterate through the entries of the sparse matrix.
        // Does not make a difference between ColMajor and RowMajor.
        for (int k=0; k< m_matrix.outerSize(); ++k)
          for (typename gsSparseMatrix<T>::InnerIterator it(m_matrix,k); it; ++it)
          {
            r = it.row();
            c = it.col();

            if( Flag_Enf[r] >= 0 )
            {
                // Row corresponds to a fix-dof.
                // Right-hand-side has been fixed above.
                // If row == column, we are on the diagonal,
                // so put a 1 on the diagonal.
                // Otherwise, it is an off-diagonal entry,
                // so put a 0 there.
                if( r == c )
                    m_matrix.coeffRef(r,c) = T(1);
                else
                    m_matrix.coeffRef(r,c) = T(0);
            }
            else if( Flag_Enf[c] >= 0 )
            {
                // Row corresponds to a non-fix-dof,
                // and the column corresponds to an fix-dof.
                // So, correct the right-hand-side...
                m_rhs.row(r) -= it.value() * Val.row( Flag_Enf[c] );
                // ...and put a zero in the matrix.
                m_matrix.coeffRef(r,c) = T(0);
            }
            // else...
            //
            // otherwise, we are at an entry where the
            // row corresponds to a non-fix-dof and the
            // column also corresponds to a non-fix-dof.
            // In this case, do nothing.

          } // end for inner iterator
    } // end if-else m_isSymmetric

    // Compress sparse matrix.
    m_matrix.makeCompressed();

} // end boundaryFixDofs



// S.Kleiss
// Assembles and returns the raw stiffness matrix without any boundary conditions
// Not really tested yet (07.May 2014).
template<class T>
gsSparseMatrix<T> * gsPoissonAssembler<T>::stiffnessMatrixPatch(int patchIndex)
{
    const gsBasis<T> & basis = m_bases[0][patchIndex];
    const index_t Ndofs = basis.size();

    // Quadrature  nodes
    gsVector<int> numCwiseNodes = gsAssemblerUtils<T>::getNumIntNodesFor( basis );
    index_t       numNodes      = numCwiseNodes.prod();
    gsGaussRule<T> QuRule( numCwiseNodes ); // reference rule
    gsMatrix<T> quNodes  ; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights

    // Empty sparse matrix to fill with the contribution of the current patch
    gsSparseMatrix<T> * Kh = new gsSparseMatrix<T>(Ndofs,Ndofs );

    gsMatrix<T> trf_grads_k;    // transformed (physical) gradients of all basis functions at one quadrature node

    // active basis functions at one quadrature node
    gsMatrix<unsigned> actives; 
    gsMatrix<T> basisDerivs;

    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
        this->m_patches.patch(patchIndex).evaluator(NEED_VALUE   |
                                                    NEED_MEASURE |
                                                    NEED_GRAD_TRANSFORM) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

        // Compute the active basis functions
        basis.active_into(domIt->centerPoint() , actives);
        const index_t numActive = actives.rows();

        // Evaluate basis function derivatives on element
        basis.deriv_into( quNodes, basisDerivs);

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(quNodes);

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);

        for (index_t k = 0; k < numNodes; ++k)      // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = quWeights[k] * geoEval->measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval->transformGradients(k, basisDerivs, trf_grads_k);

            localStiffness.noalias() += weight * (trf_grads_k.transpose() * trf_grads_k);
        }  // end loop Gauss nodes

        // add contributions from local stiffness matrix to global one.
        // Note: assembled as FULL stiffness matrix by setting "symmetric" to false.
        gsAssemblerUtils<T>::localToGlobal( localStiffness, actives, * Kh, false);
    } //end loop over all domain elements

    Kh->makeCompressed();
    return Kh;
}


// S.Kleiss
// Assembles and returns the raw mass matrix without any boundary conditions
// Not really tested yet (07.May 2014).
template<class T>
gsSparseMatrix<T> * gsPoissonAssembler<T>::massMatrixPatch(int patchIndex)
{
    const gsBasis<T> & basis = m_bases[0][patchIndex];
    const index_t Ndofs = basis.size();

    // Quadrature  nodes
    gsVector<int> numCwiseNodes = gsAssemblerUtils<T>::getNumIntNodesFor( basis );
    index_t       numNodes      = numCwiseNodes.prod();
    gsGaussRule<T> QuRule( numCwiseNodes ); // reference rule
    gsMatrix<T> quNodes  ; // Mapped nodes
    gsVector<T> quWeights; // Mapped weights

    // Empty sparse matrix to fill with the contribution of the current patch
    gsSparseMatrix<T> * Mh = new gsSparseMatrix<T>(Ndofs,Ndofs );

    // Estimate max non-zeros per row (only valid for tensor product bases right now)
    unsigned nzRowsPerCol = 1;
    for (int i = 0; i < basis.dim(); ++i)
        nzRowsPerCol *= 2 * basis.degree(i) + 1;
    // Reserve space
    // \todo get the fragment of dofs in patch numbered patchIndex
    Mh->reserve( gsVector<int>::Constant(Ndofs, nzRowsPerCol) );


    gsMatrix<T> localMass; // (dense) mass matrix within one grid cell

    // active basis functions at one quadrature node
    gsMatrix<unsigned> actives; 
    gsMatrix<T>    basisValues;

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
        this->m_patches.patch(patchIndex).evaluator(NEED_VALUE   |
                                                    NEED_MEASURE) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        // Map the Quadrature rule to the element
        QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );

        // Compute the active basis functions
        basis.active_into(domIt->centerPoint() , actives);
        const index_t numActive = actives.rows();

        // Evaluate basis functions on element
        basis.eval_into( quNodes, basisValues);

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(quNodes);

        // initialize local linear system to 0
        localMass.setZero(numActive, numActive);

        for (index_t k = 0; k < numNodes; ++k)      // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = quWeights[k] * geoEval->measure(k);

            localMass.noalias() += weight * ( basisValues.col(k) * basisValues.col(k).transpose());
        }  // end loop Gauss nodes

        // add contributions from local mass matrix to global one.
        // Note: assembled as FULL mass matrix by setting "symmetric" to false.
        gsAssemblerUtils<T>::localToGlobal( localMass, actives, * Mh, false);

    } //end loop over all domain elements

    Mh->makeCompressed();
    return Mh;
}


// S.Kleiss
// Should become main function to call error indication.
// Currently just calling the one with bubble functions
// TODO: Add a way of selecting one out of several error indicators
template<class T>
gsVector< gsMatrix<T> > gsPoissonAssembler<T>::errorIndicator()
{
    gsVector< gsMatrix<T> > errInd;
    errInd.resize( this->m_patches.nPatches() );

    // For the time being, only one error indicator implemented
    // TODO: Add other error indicators and
    // add a way of selecting one of them.
    for( unsigned i=0; i<this->m_patches.nPatches(); i++)
        errInd[i] = errIndPatch_Bubble( i );

    return errInd;
} //end gsPoissonAssembler<T>::errorIndicator()


// S.Kleiss
// Uses error indicator with bubble functions.
// Based on error estimator presented in
// M. R. Doerfel, B. Juettler, B. Simeon. Adaptive isogeometric analysis by local
// h-refinement with T-splines. Computer Methods in Applied Mechanics and
// Engineering, 199: 264-275, 2010.
// with a certin modification that will be documented sometime soon.
template<class T>
gsMatrix<T> gsPoissonAssembler<T>::errIndPatch_Bubble( const index_t patchIndex )
{

    gsDebug << "Start: errIndPatch_Bubble()\n";

    // access the local basis on patch:
    const gsBasis<T> & basis = m_bases[0][patchIndex];
    // construct domain iterator for basis:
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // number of elements on patch:
    const index_t NE = basis.numElements();
    // dimension of the parameter space of the basis:
    const index_t d = basis.dim();

    // access local solution...
    const gsField<T> * sol = m_solutions[patchIndex];
    // ...and extract the "discrete-solution-part", and...
    const gsGeometry<T> & igaFct = sol->igaFunction(patchIndex);
    // ...the "geometry-information-part".
    //const gsGeometry<T> * geo = sol->patch(patchIndex);

    // initialize the matrix for the local errors and
    // element coordinates
    gsMatrix<T> errInd( NE, 1 + 2*d);

    // compute ...
    // ...the number of quadrature nodes
    gsVector<int> numQuNodes(d);
    // ...the degrees of the bubbles
    gsVector<index_t> bubbleDeg(d);
    // ...the number of bubble functions per coordinate direction
    gsVector<index_t> bub_n(d);
    // ...the total number of basis functions
    index_t bubInt_nTotal = 1;
    for( index_t i=0; i < d; i++)
    {
        bubbleDeg[i] = basis.degree(i) + 1;
        //bubbleDeg[i] = 2;
        if( bubbleDeg[i] < 2 )         // at least degree 2, otherwise no bubble
            bubbleDeg[i] = 2;
        bub_n[i] = bubbleDeg[i] + 1;    // number of univariate basis fcts
        bubInt_nTotal *= bub_n[i] - 2;  // number of interior univ. basis fcts

        numQuNodes[i] = bubbleDeg[i] + 1;
    }

    // compute the indices of the interior bubble functions
    gsVector<index_t> bubInt_indices( bubInt_nTotal );
    if( d == 2)
    {
        index_t c = 0;
        for( int j=1; j < (bub_n[1]-1); j++)
            for( int i=1; i < (bub_n[0]-1); i++)
            {
                bubInt_indices[c] = i + bub_n[0] * j;
                c += 1;
            }
    }

    if( d == 3)
    {
        index_t c = 0;
        for( int k=1; k < (bub_n[2]-1); k++)
            for( int j=1; j < (bub_n[1]-1); j++)
                for( int i=1; i < (bub_n[0]-1); i++)
                {
                    bubInt_indices[c] = i + bub_n[0] * j + bub_n[0]*bub_n[1] * k;
                    c += 1;
                }
    }

    GISMO_ASSERT( bubInt_indices.size() == bubInt_nTotal, " bub INT! ");

    gsVector< gsKnotVector<T> > bub_KV(3);

    for( index_t i = 0; i < d; i++)
    {
        bub_KV[i] = gsKnotVector<T>(0,1,0, bubbleDeg[i]+1 );
    }

    // initialize basis for the bubbles
    gsBasis<T>* bubbles;
    if( d == 2 )
    {
        bubbles = new gsTensorBSplineBasis<2,T>( bub_KV[0], bub_KV[1] );
    }
    else
    {
        bubbles = new gsTensorBSplineBasis<3,T>( bub_KV[0], bub_KV[1], bub_KV[2] );
    }

    // initialize the quadrature rule in the domain iterator
    gsGaussRule<T> QuRule( numQuNodes );
    const index_t numNodes = numQuNodes.prod();

    // could be more efficient, but good enough for now:

    // Evaluate the bubbles once on the reference hypercube
    gsVector<T> v0(d);
    gsVector<T> v1(d);
    v0.setZero();
    v1.setOnes();
    gsMatrix<T> Qp;
    gsVector<T> Qw;
    QuRule.mapTo( v0,v1, Qp, Qw);
    gsMatrix<T> bub_v_tmp = bubbles->eval( Qp );
    gsMatrix<T> bub_d_tmp = bubbles->deriv( Qp );

    gsMatrix<T> bub_v( bubInt_nTotal, bub_v_tmp.cols() );
    gsMatrix<T> bub_d( d * bubInt_nTotal, bub_d_tmp.cols() );

    for( index_t i = 0; i < bubInt_indices.size(); i++)
    {
        index_t ii = bubInt_indices[i];
        bub_v.row(i) = bub_v_tmp.row(ii);
        bub_d.row(d*i) = bub_d_tmp.row(d*ii);
        bub_d.row(d*i+1) = bub_d_tmp.row(d*ii+1);
        if( d == 3)
            bub_d.row(d*i+2) = bub_d_tmp.row(d*ii+2);
    }

    bub_v_tmp.resize( bub_v.rows(), bub_v.cols() );
    bub_d_tmp.resize( bub_d.rows(), bub_d.cols() );
    bub_v_tmp.setZero();
    bub_d_tmp.setZero();

    // set up gsGeometryEvaluator
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
                this->m_patches.patch(patchIndex).evaluator(NEED_VALUE   |
                                                            NEED_MEASURE |
                                                            NEED_GRAD_TRANSFORM) );

    gsMatrix<T> rhsVals;        // values of the right-hand side
    gsMatrix<T> uhGrads;        // gradients of discrete solution
    gsMatrix<T> trf_grads_k;    // transformed (physical) gradients of all basis functions at one quadrature node
    gsMatrix<T> trf_uhGradsT;   // transformed (physical) gradients of discrete solution

    Eigen::Matrix<T,Dynamic,Dynamic> Lh( bubInt_nTotal, bubInt_nTotal);     // local system matrix
    Eigen::Matrix<T,Dynamic,Dynamic> Rh( bubInt_nTotal, 1 );     // local right-hand-side

    index_t actElt = 0;
    for (; domIt->good(); domIt->next(), actElt++ )
    {
        Lh.setZero();
        Rh.setZero();

        const gsVector<T> & lowC = domIt->lowerCorner();
        const gsVector<T> & uppC = domIt->upperCorner();
        const gsVector<T> h = uppC - lowC;

        // adjust the derivatives of the bubble functions
        // to the current element.
        bub_d_tmp = bub_d;
        for( index_t i = 0; i < d ; i++)
        {
            for( index_t j = 0; j < bub_v.rows(); j++)
                bub_d_tmp.row( j*d + i) /= h[i];
        }

        // Map the Quadrature rule to the element
        QuRule.mapTo( lowC, uppC, Qp, Qw );

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(Qp);

        // compute gradients of discrete solution and resize result
        igaFct.deriv_into( Qp, uhGrads );
        gsMatrix<T> uhGradsT(d, bub_v.cols() );
        for( int i=0; i < d; i++)
            for( int j=0; j < bub_v.cols(); j++)
                uhGradsT(i,j) = uhGrads(0,j*d+i);

        // evaluate right-hand side at the geometry points
        if (m_rhs_function)
            m_rhs_function->eval_into( geoEval->values(), rhsVals );

        for (index_t k = 0; k < numNodes; ++k)      // loop over quadrature nodes
        {
            // update quadrature weights with determinant of Jacobian.
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = Qw[k] * geoEval->measure(k);

            // compute physical gradients of bubble functions at k as a Dim x NumActive matrix
            geoEval->transformGradients(k, bub_d_tmp, trf_grads_k);

            // compute physical gradients of discrete solution
            geoEval->transformGradients(k, uhGradsT, trf_uhGradsT);

            // add to right-hand-side
            if (m_rhs_function)
                Rh += bub_v.col(k) * (weight * rhsVals.col(k)).transpose();

            Rh -= weight * ( trf_grads_k.transpose() * trf_uhGradsT);

            // add to system matrix
            Lh.noalias() += weight * (trf_grads_k.transpose() * trf_grads_k);

        }  // end loop Gauss nodes

        // apply direct solver to these small local problems
        Eigen::Matrix<T,Dynamic,Dynamic> eh = Lh.fullPivLu().solve( Rh );

        // recover the error indicator
        T errIndLocal = 0.0;
        for( index_t i = 0; i < eh.size(); i++)
            for( index_t j = 0; j < eh.size(); j++)
                errIndLocal += eh(i,0) * eh(j,0) * Lh(i,j);

        // store result and coordinates of element in parameter domain
        errInd( actElt, 0 ) = errIndLocal;
        for( index_t i = 0; i < d; i++)
        {
            errInd( actElt, 1 + i) = lowC[i];
            errInd( actElt, 1 + d + i) = uppC[i];
        }

    } //end for domIt

    // ...that's it...

    gsDebug << "--End: errIndPatch_Bubble()\n";

    delete bubbles;
    return errInd;

} // end gsPoissonAssembler<T>::errorIndicatorBubble()



template<class T>
void gsPoissonAssembler<T>::adaptiveRefine( const gsVector< gsMatrix<T> >& errInd, const int refCriterion, const T a )
{
    gsDebug << "Start: adaptiveRefine()\n";
    gsDebug << "         Criterion " << refCriterion << "\n";

    T Thr = 0;
    // Total number of elements
    int NE = 0;

    if( refCriterion == 1 )
    {
        // First, conduct a brutal search for the maximum local error
        T maxErr = 0;
        for( int patchIndex = 0; patchIndex < errInd.size(); patchIndex++)
        {
            NE += errInd[patchIndex].rows();
            for( int i = 0; i < errInd[patchIndex].rows(); i++ )
                if( maxErr < errInd[patchIndex](i,0) )
                    maxErr = errInd[patchIndex](i,0);
        }

        // Compute the threshold:
        // Note that the element-wise errors estimates are given as
        // estimates of the squared local error. Hence, parameter a is
        // also squared to compute the threshold.
        Thr = a*a * maxErr;
    }
    else if ( refCriterion == 2 || refCriterion == 3 )
    {
        for( int patchIndex = 0; patchIndex < errInd.size(); patchIndex++)
            NE += errInd[patchIndex].rows();

        // Get the list of all local errors estimates and the total error estimate:
        std::vector<T> errList;
        errList.reserve( NE );
        T TotalError = 0;
        for( int patchIndex = 0; patchIndex < errInd.size(); patchIndex++)
            for( int i = 0; i < errInd[patchIndex].rows(); i++ )
            {
                errList.push_back( errInd[patchIndex](i,0) );
                TotalError += errInd[patchIndex](i,0);
            }

        if( refCriterion == 2)
        {

            // For RefCriterion 2:
            // Assume a sorted (from small to large) list of local errors:
            // Compute the index of the list-element from which
            // the refinement should start.
            unsigned IdxRefineStart = static_cast<unsigned>( floor( a * NE) );
            if( IdxRefineStart == errList.size() )
                IdxRefineStart -= 1;

            // Sort the list using bubblesort.
            // After each loop, the largest elements are at the end
            // of the list. Since we are only interested in the largest elements,
            // it is enough to run the sorting until enough "largest" elements
            // have been found.
            unsigned LastSwapDone = errList.size() - 1;
            unsigned LastCheckIdx = LastSwapDone;

            bool DidSwap;
            do{
                DidSwap = false;
                T tmp;                
                LastCheckIdx = LastSwapDone;
                for( unsigned i=0; i < LastCheckIdx; i++)
                    if( errList[i] > errList[i+1] )
                    {
                        tmp = errList[i];
                        errList[i] = errList[i+1];
                        errList[i+1] = tmp;

                        DidSwap = true;
                        LastSwapDone = i;
                    }
            }while( DidSwap && (LastSwapDone+1 >= IdxRefineStart ) );
            Thr = errList[ IdxRefineStart ];

            gsDebug << "         NE        = " << NE << "\n";
            gsDebug << "         IdxStart  = " << IdxRefineStart << "\n";
            gsDebug << "         Thr       = " << Thr << "\n";
        } // end refCrit == 2

        if( refCriterion == 3)
        {
            // For RefCriterion 3:
            T ErrorReduce = (1-a) * TotalError;
            T cummulErrRed = 0;
            unsigned LastSwapDone = errList.size() - 1;

            do{
                T tmp;
                for( unsigned i=0; i < LastSwapDone; i++)
                    if( errList[i] > errList[i+1] )
                    {
                        tmp = errList[i];
                        errList[i] = errList[i+1];
                        errList[i+1] = tmp;
                    }

                cummulErrRed += errList[ LastSwapDone  ];
                LastSwapDone -= 1;

            }while( cummulErrRed < ErrorReduce && LastSwapDone > 0 );

            Thr = errList[ LastSwapDone + 1 ];
        } // end refCrit == 3

    } // end refCrit == 2 || refCrit == 3

    // Initialize the refinement-box
    const int d = m_bases[0][0].dim();

    for( int patchIndex = 0; patchIndex < errInd.size(); patchIndex++)
    {
        // Using std::vectors, because the number of cells to be refined
        // is not known a priori.
        std::vector< std::vector<T> > Ref_lo;
        std::vector< std::vector<T> > Ref_up;

        int RefCount = 0;
        for( int i = 0; i < errInd[patchIndex].rows(); i++ )
            if( Thr <= errInd[patchIndex](i,0) )
            {
                // If the local error is larger than the threshold, refine.

                std::vector<T> Ref_lo_t;
                std::vector<T> Ref_up_t;

                for( int j = 0; j < d; j++)
                {
                    // Write the corners of the element in the RefBox.
                    // The little trick with the "h" is necessary,
                    // because the refine()-method does not work properly,
                    // if coordinates of the refinement-box coincide with knots.
                    T h = errInd[patchIndex](i,1+d+j) - errInd[patchIndex](i,1+j);

                    Ref_lo_t.push_back( errInd[patchIndex](i,1+j) + h/10 );
                    Ref_up_t.push_back( errInd[patchIndex](i,1+d+j) - h/10 );
                }

                Ref_lo.push_back( Ref_lo_t );
                Ref_up.push_back( Ref_up_t );

                RefCount += 1;
            }

        gsMatrix<T> RefBox( d, 2*Ref_lo.size() );
        RefBox.setZero();
        for( int di=0; di < d; di++)
            for( unsigned i=0; i < Ref_lo.size(); i++)
            {
                RefBox(di,2*i) = Ref_lo[i][di];
                RefBox(di,2*i+1) = Ref_up[i][di];
            }

        gsDebug << "         refining " << RefCount << " of " << NE << " Element(s) = " << 100*RefCount/NE << "%\n";
        m_bases[0][ patchIndex ].refine( RefBox );
    }

    gsDebug << "--End: adaptiveRefine()\n";
} // end gsPoissonAssembler<T>::adaptiveRefinement


template<class T>
void gsPoissonAssembler<T>::uniformRefine()
{
    for (unsigned np=0; np < this->m_patches.nPatches(); ++np )
    {
        m_bases[0][np]->uniformRefine();
    }
}



} // namespace gismo


