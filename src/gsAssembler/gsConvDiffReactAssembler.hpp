#pragma once

#include <gsCore/gsDebug.h>
#include <gsAssembler/gsAssemblerUtils.h>
#include <gsAssembler/gsConvDiffReactAssembler.h>

//#include <gismo.h>

namespace gismo {


//Virtual function
template<class T>
void gsConvDiffReactAssembler<T>::initialize()
{
    //Tells the user that elimination strategy does not work in combination with dg.
    GISMO_ASSERT(!(m_interfaceStrategy == dg && m_dirStrategy == elimination),
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
void gsConvDiffReactAssembler<T>::assemble()
{
    // If the Dirichlet strategy is elimination then precompute
    // Dirichlet dofs (m_dofMapper excludes these from the system)
    if ( m_dirStrategy == elimination || m_dirStrategy == homogeneous)
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
    if ( m_dirStrategy == nitsche )
    {
        for ( typename gsBoundaryConditions<T>::const_iterator
                  it = this->m_bconditions[0]->dirichletBegin();
              it != this->m_bconditions[0]->dirichletEnd(); ++it )
        {
            // to do remove members from args
            applyBoundary(*this->m_bases[0][it->patch()], *it, *this->m_dofMapper[0]);
        }
    }

    // Enforce Neumann boundary conditions
    for ( typename gsBoundaryConditions<T>::const_iterator
          it = this->m_bconditions[0]->neumannBegin();
          it != this->m_bconditions[0]->neumannEnd(); ++it )
    {
        applyBoundary(*this->m_bases[0][it->patch()], *it, *this->m_dofMapper[0]);
    }


    // If we are in in dg (Discontinuous Galerkin) mode add interface
    // contributions
    if ( m_interfaceStrategy == dg )
    {
        for ( typename gsMultiPatch<T>::const_iiterator it =
                  m_patches.iBegin(); it != m_patches.iEnd(); ++it )
        {
            applyDG( *this->m_bases[0][it->ps1.patch], *this->m_bases[0][it->ps2.patch],
                     m_patches.patch( it->ps1.patch ), m_patches.patch( it->ps2.patch ),
                     *it, *this->m_dofMapper[0], m_matrix );
        }
    }

}


template<class T>
void gsConvDiffReactAssembler<T>::assembleMultipatch()
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
void gsConvDiffReactAssembler<T>::assemblePatch(int patchIndex)
{
    const gsBasis<T> & basis = * m_bases[0][patchIndex];
    const gsDofMapper  & mapper = *m_dofMapper[0];
    const gsMatrix<T> & ddof = m_fixedDofs[0];

    const int d = basis.dim();
    GISMO_ASSERT( d == 2 || d == 3, "Only d==2 or d==3 considered!");

    if( m_doSUPG)
        GISMO_ASSERT( m_isSymmetric == false, "Cannot do SUPG if you say the stiffness matrix is symmetric");

    gsVector<int> numNodes = gsAssemblerUtils<T>::getNumIntNodesFor( basis );
    for( int i=0; i < numNodes.size(); i++)
        numNodes[i] += 1;

    // Empty sparse matrix to fill with the contribution of the current patch
    gsSparseMatrix<T> sys_K(m_idofs, m_idofs );

    // Estimate max non-zeros per row (only valid for tensor product bases right now)
    unsigned nzRowsPerCol = 1;
    for (int i = 0; i < basis.dim(); ++i)
        nzRowsPerCol *= 2 * basis.component(i).degree() + 1;
    // Reserve space
    // \todo get the fragment of dofs in patch numbered patchIndex
    sys_K.reserve( gsVector<int>::Constant(m_idofs, nzRowsPerCol) );

    // Empty sparse rhs vector to fill with the contribution of the current patch
    gsMatrix<T> sys_b(mapper.freeSize(), this->m_unknownDim[0]);
    sys_b.setZero();

    // Local values of the right-hand side
    gsMatrix<T> rhsVals;
    // Local values of the coefficients
    gsMatrix<T> coeff_A_vals;
    gsMatrix<T> coeff_b_vals;
    gsMatrix<T> coeff_c_vals;
    // Auxiliary variable
    gsMatrix<T> tmp_A;

    // basisData contains stacked the values and the gradients of all
    // basis functions at one quadrature node
    // trf_grads_pp keeps the (transformed) physical gradients of the
    // basis functions at one quadrature node
    // b_grads_pp will be used to pre-compute "b \cdot \nabla v"
    gsMatrix<T> basisData;
    gsMatrix<T> trf_grads_pp;
    gsMatrix<T> b_grads_pp;

    // Active basis functions at one quadrature node
    gsMatrix<unsigned> actives;

    gsMatrix<T> localMat;     // (dense) system matrix within one grid cell
    gsMatrix<T> localMatSUPG; // SUPG-part for local system matrix
    gsMatrix<T> localRhs;     // local load vector
    gsMatrix<T> localRhsSUPG; // SUPG-part for local load vector

    T SUPG_tau = T(0.0); // SUPG stabilization parameter
    gsMatrix<T> SUPG_termA_pp; // auxiliary variable

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
    this->m_patches.patch(patchIndex).evaluator(NEED_VALUE   |
                                                NEED_MEASURE |
                                                NEED_GRAD_TRANSFORM) );

    int targetDim = this->m_unknownDim[0];
    GISMO_ASSERT(targetDim == 1, "AS OF NOW, ONLY targetDim == 1 IMPLEMENTED!");

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        // Compute the active basis functions
        basis.active_into(domIt->center, actives);
        const index_t numActive = actives.rows();
        //domIt->computeActiveDofs(mapper, patchIndex).rows();
        mapper.localToGlobal(actives, patchIndex, actives);

        // compute and "split" up function values and first
        // and second derivatives.
        basis.evalAllDers_into(domIt->quNodes, 2, basisData);

        const typename gsMatrix<T>::Block basisValues = basisData.topRows(numActive);
        const typename gsMatrix<T>::Block basisGrads  =
            basisData.middleRows(numActive, numActive*d);
        const typename gsMatrix<T>::Block basis2ndDerivs  =
            basisData.middleRows(numActive*(1+d), numActive*d*(d+1)/2 );

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(domIt->quNodes);

        // evaluate right-hand side at the geometry points
        if (m_rhs_function)
            m_rhs_function->eval_into( geoEval->values(), rhsVals );

        // evaluate coefficients of the PDE
        m_coeff_A->eval_into( geoEval->values(), coeff_A_vals );
        m_coeff_b->eval_into( geoEval->values(), coeff_b_vals );
        m_coeff_c->eval_into( geoEval->values(), coeff_c_vals );

        // initialize local linear system to 0
        localMat.setZero(numActive, numActive);
        localRhs.setZero(numActive, targetDim);
        localMatSUPG.setZero(numActive, numActive);
        localRhsSUPG.setZero(numActive, targetDim);

        // auxiliary variable, number of second derivatives
        const int d2_offsetter = d*(d+1)/2;

        for (index_t pp = 0; pp < domIt->numQuNodes(); ++pp)      // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = domIt->quWeights[pp] * geoEval->measure(pp);

            // compute physical gradients at quadrature point k
            // d... dimension of parameter space
            // N... number of active functions

            // trf_grads_pp:            d x N
            geoEval->transformGradients(pp, basisGrads, trf_grads_pp);
            int N = trf_grads_pp.cols();

            // coeff_b_vals.col(pp):    d x 1
            // trf_grads_pp:            d x N
            // b_grads_pp:              1 x N
            b_grads_pp = coeff_b_vals.col(pp).transpose() * trf_grads_pp;

            // "Standard" contribution to right-hand-side
            localRhs += basisValues.col(pp) * (weight * rhsVals.col(pp)).transpose();

            // Diffusion coefficient
            tmp_A = coeff_A_vals.col(pp);
            tmp_A.resize(d,d);

            // --- "Standard" 2nd-order term
            // trf_grads_pp.transpose(): N x d
            // A:                        d x d
            // trf_grads_pp:             d x N
            localMat.noalias() += weight * (trf_grads_pp.transpose() * (tmp_A * trf_grads_pp) );

            // --- "Standard" 1st-order term
            // b_grads_pp:              1 x N
            // basisValues.col(pp):     N x 1
            localMat.noalias() += weight * ( basisValues.col(pp) * b_grads_pp );

            // --- "Standard" 0th-order term
            // coeff_c_vals.col(pp):    1 x 1
            // basisValues.col(pp):     N x 1
            localMat.noalias() += weight * coeff_c_vals(0,pp) * ( basisValues.col(pp) * (basisValues.col(pp).transpose() ) );

            if( m_doSUPG )
            {
                // Compute stabilization parameter
                SUPG_tau = get_SUPG_parameter( domIt->lowerCorner(), domIt->upperCorner(), patchIndex);

                // Inverse of Jacobian...
                // --- WARNING:
                // The function gradTransforms() returns the
                // TRANSPOSED Jacobian at the moment (01.Jul.2014)!
                gsMatrix<T> JinvT = geoEval->gradTransforms(); // transposed invJacs

                // SUPG_termA_pp:            d x N
                SUPG_termA_pp.resize( d, N );
                SUPG_termA_pp.setZero();
                for( int ni=0; ni < N ; ni++)
                {
                    for( int dj=0; dj < d; dj++)
                        for( int i=0; i < d; i++ )
                        {
                            T tmp = T(0);

                            for( int k = 0; k < d; k++)
                                for( int l = 0; l < d; l++)
                                {
                                    int kl = 2*ni*d2_offsetter; // default will throw an error
                                    if( k == l )
                                        kl = k;
                                    else if( (k==0 && l==1) || (k==1 && l==0 ) )
                                        kl = d;
                                    else if( (k==0 && l==2) || (k==2 && l==0 ) )
                                        kl = 4;
                                    else if( (k==1 && l==2) || (k==2 && l==1 ) )
                                        kl = 5;

                                    // WORKING WITH THE TRANSPOSED INVERSE JinvT
                                    tmp += basis2ndDerivs(ni*d2_offsetter+kl,pp)
                                            * JinvT( i, pp*d + k ) * JinvT( dj, pp*d + l );

                                    // Version where Jinv really is the inverse of the Jacobian
                                    //tmp += basis2ndDerivs(ni*d2_offsetter+kl,pp)
                                    //        * Jinv( k, ni*d + i ) * Jinv( l, ni*d + dj );
                                }
                            SUPG_termA_pp(dj,ni) += coeff_b_vals(i,pp) * tmp;
                        }
                }

                // System matrix contributions due to SUPG stabilization
                // --- SUPG 2nd-order term
                // SUPG_termA_pp:   d x N
                // tmp_A:           d x d
                // trf_grads_pp:    d x N
                localMatSUPG.noalias() += weight * ( SUPG_termA_pp.transpose() * (tmp_A * trf_grads_pp) );

                // --- SUPG 1st-order term
                localMatSUPG.noalias() += weight * ( b_grads_pp.transpose() * b_grads_pp );

                // --- SUPG 0th-order term
                localMatSUPG.noalias() += (weight * coeff_c_vals(0,pp) )*( b_grads_pp.transpose() * ( basisValues.col(pp).transpose() ) );

                // Right-hand-side contributions due to SUPG stabilization
                localRhsSUPG.noalias() += weight * ( b_grads_pp.transpose() * (rhsVals.col(pp)).transpose() );
            }
        }  // end loop Gauss nodes

        // Add SUPG terms to system matrix and right-hand-side
        if( m_doSUPG )
        {
            localMat.noalias() += SUPG_tau * localMatSUPG;
            localRhs.noalias() += SUPG_tau * localRhsSUPG;
        }

        // Add contributions from local stiffness matrix to global
        // stiffness matrix and load vector
        gsAssemblerUtils<T>::localToGlobal_withBC(localMat, localRhs, ddof, mapper,
                                                 actives.cast<int>(), sys_K, sys_b, m_isSymmetric);
    } //end loop over all domain elements

    this->m_matrix += sys_K;
    this->m_rhs    += sys_b;
}

template<class T>
T gsConvDiffReactAssembler<T>::get_SUPG_parameter( const gsVector<T> & lo, const gsVector<T> & up, const int patchIdx )
{

    const int N = 2;

    const int d = lo.size();
    gsMatrix<T> ctr(d,1);
    gsMatrix<T> phys_ctr;
    gsMatrix<T> b_ctr;

    // compute the center point of the cell...
    for( int i=0; i < d; i++)
        ctr(i,0) = ( lo[i] + up[i] ) / 2;
    // ...map it to the physical space...
    m_patches[ patchIdx ].eval_into( ctr, phys_ctr );
    // ...evaluate the convection coefficient there, ...
    m_coeff_b->eval_into( phys_ctr, b_ctr );
    // ...and get it's norm.
    T b_norm = 0;
    for( int i=0; i < d; i++)
        b_norm += b_ctr(i,0) * b_ctr(i,0);
    b_norm = math::sqrt( b_norm );

    T SUPG_tau = T(0);
    if( b_norm > 0 )
    {
        gsMatrix<T> bdryPts;
        gsMatrix<T> phys_bdryPts;
        if( d == 2 )
        {
            int ij = 0;
            T a;
            bdryPts.resize(d,4*N);

            // west and east
            for( int i=0; i <= N; i++ )
            {
                a = i/N;
                bdryPts(0,ij) = lo[0];
                bdryPts(1,ij) = (1-a) * lo[1] + a * up[1];
                bdryPts(0,ij+1) = up[0];
                bdryPts(1,ij+1) = (1-a) * lo[1] + a * up[1];
                ij += 2;
            }
            // south and north
            for( int i=1; i <= N-1; i++ )
            {
                a = i/N;
                bdryPts(0,ij) = (1-a) * lo[0] + a * up[0];
                bdryPts(1,ij) = lo[1];
                bdryPts(0,ij+1) = (1-a) * lo[0] + a * up[0];
                bdryPts(1,ij+1) = up[1];
                ij += 2;
            }
        }
        else if( d == 3 )
        {
            int ij = 0;
            T a1;
            T a2;

            // redundant points.
            bdryPts.resize(d,6*(N+1)*(N+1));
            for( int i=0; i <= N; i++ )
                for( int j=0; j <= N; j++ )
                {
                    a1 = i/N;
                    a2 = j/N;
                    bdryPts(0,ij) = lo[0];
                    bdryPts(1,ij) = (1-a1) * lo[1] + a1 * up[1];
                    bdryPts(2,ij) = (1-a2) * lo[2] + a2 * up[2];
                    bdryPts(0,ij+1) = up[0];
                    bdryPts(1,ij+1) = (1-a1) * lo[1] + a1 * up[1];
                    bdryPts(2,ij+1) = (1-a2) * lo[2] + a2 * up[2];
                    ij += 2;
                }
            for( int i=0; i <= N; i++ )
                for( int j=0; j <= N; j++ )
                {
                    a1 = i/N;
                    a2 = j/N;
                    bdryPts(0,ij) = (1-a1) * lo[0] + a1 * up[0];
                    bdryPts(1,ij) = lo[1];
                    bdryPts(2,ij) = (1-a2) * lo[2] + a2 * up[2];
                    bdryPts(0,ij+1) = (1-a1) * lo[0] + a1 * up[0];
                    bdryPts(1,ij+1) = up[1];
                    bdryPts(2,ij+1) = (1-a2) * lo[2] + a2 * up[2];
                    ij += 2;
                }
            for( int i=0; i <= N; i++ )
                for( int j=0; j <= N; j++ )
                {
                    a1 = i/N;
                    a2 = j/N;
                    bdryPts(0,ij) = (1-a1) * lo[0] + a1 * up[0];
                    bdryPts(1,ij) = (1-a1) * lo[1] + a1 * up[1];
                    bdryPts(2,ij) = lo[2];
                    bdryPts(0,ij+1) = (1-a1) * lo[0] + a1 * up[0];
                    bdryPts(1,ij+1) = (1-a1) * lo[1] + a1 * up[1];
                    bdryPts(2,ij+1) = up[2];
                    ij += 2;
                }


            GISMO_ERROR(" d == 3 not tested yet (03.jul.2014)");
        }

        m_patches[ patchIdx ].eval_into( bdryPts, phys_bdryPts );

        T vmin = T(0);
        for( int i = 0; i < d; i++)
            vmin += phys_bdryPts(i,0) * b_ctr(i);
        T vmax = vmin;

        for( int k=1; k<phys_bdryPts.cols(); k++)
        {
            T vtmp = T(0);
            for( int i = 0; i < d; i++)
                vtmp += phys_bdryPts(i,k) * b_ctr(i);

            if( vtmp < vmin )
                vmin = vtmp;
            if( vtmp > vmax )
                vmax = vtmp;
        }

        SUPG_tau = ( vmax - vmin ) / ( 2 *  b_norm );
    }
    //else
    //{
    //    SUPG_tau = 0;
    //}
    // ...taken care of by initialization of SUPG_tau

    return SUPG_tau;
}

// Solves the linear system and fills up \a m_sysSolution
template<class T>
void gsConvDiffReactAssembler<T>::solveSystem()
{
    // Initialize solver
    //Eigen::ConjugateGradient< gsSparseMatrix<T> > solver; //
    Eigen::BiCGSTAB< gsSparseMatrix<T>, Eigen::IncompleteLUT<T> > solver;
    // Solve linear system
    this->m_sysSolution = solver.compute(  this->m_matrix ).solve (  this->m_rhs );

    gsInfo << "residual error: " << solver.error() << "\n";
    gsInfo << "    iterations: " << solver.iterations() << "\n";
}


template<class T> void
gsConvDiffReactAssembler<T>::applyBoundary( const gsBasis<T>   & B,
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
gsConvDiffReactAssembler<T>::applyDG( const gsBasis<T> & B1, const gsBasis<T> & B2,
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
    gsVector<int> intNodes1 = gsAssemblerUtils<T>::getNumIntNodesForInterface( B1, B2, bi, true  );
    gsVector<int> intNodes2 = gsAssemblerUtils<T>::getNumIntNodesForInterface( B1, B2, bi, false );

    // Evaluators for the two patches
    std::auto_ptr< gsGeometryEvaluator<T> >
        geoEval1 ( geo1.evaluator(NEED_VALUE| NEED_JACOBIAN| NEED_GRAD_TRANSFORM) );
    std::auto_ptr< gsGeometryEvaluator<T> >
        geoEval2 ( geo2.evaluator(NEED_VALUE| NEED_JACOBIAN| NEED_GRAD_TRANSFORM) );

    // Temporaries
    gsMatrix<T> grads_k_1, grads_k_2;
    gsVector<T> unormal(d);

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
        domIter1->computeQuadratureRule(intNodes1);
        domIter2->computeQuadratureRule(intNodes2);

        // Push forward the quad-points to the physical domain
        geoEval1->evaluateAt(domIter1->quNodes);
        geoEval2->evaluateAt(domIter2->quNodes);

        // Evaluate basis functions and their first derivatives
        // assuming numActive1=numActive2
        const index_t numActive = domIter1->computeActiveDofs(mapper, patch1).rows();
        domIter1->evaluateBasis( 1 );
        const typename gsMatrix<T>::Block ev1  = domIter1->basisValues();
        domIter2->computeActiveDofs(mapper, patch2);
        domIter2->evaluateBasis( 1 );
        const typename gsMatrix<T>::Block ev2  = domIter2->basisValues();

        B11.setZero(numActive, numActive); B22.setZero(numActive, numActive);
        B12.setZero(numActive, numActive); B21.setZero(numActive, numActive);
        E11.setZero(numActive, numActive); E22.setZero(numActive, numActive);
        E12.setZero(numActive, numActive); E21.setZero(numActive, numActive);

        // assuming domIter1->quNodes.cols() == domIter2->quNodes.cols()
        for (index_t k=0; k!= domIter1->numQuNodes(); ++k)
        {
            // Compute the outer normal vector from patch1
            geoEval1->outerNormal(k, side1, unormal);

            // Integral transformation and quadrature weight (patch1)
            // assumed the same on both sides
            const T fff = domIter1->quWeights[k] *  unormal.norm();

            // Compute the outer unit normal vector from patch1
            unormal.normalize();

            // Transform the basis gradients
            geoEval1->transformGradients(k, domIter1->basisDerivs(1), grads_k_1);
            geoEval2->transformGradients(k, domIter2->basisDerivs(1), grads_k_2);

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
            const index_t jj1 = domIter1->activeDofs[j]; // N1_j
            const index_t jj2 = domIter2->activeDofs[j]; // N2_j
            for (index_t i=0; i!=numActive; ++i)
            {
                const index_t  ii1 = domIter1->activeDofs[i]; // N1_i
                const index_t  ii2 = domIter2->activeDofs[i]; // N2_i

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
void gsConvDiffReactAssembler<T>::boundaryFixDofs( const gsVector<unsigned> & Idx ,
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
gsSparseMatrix<T> * gsConvDiffReactAssembler<T>::stiffnessMatrixPatch(int patchIndex)
{
    const gsBasis<T> & basis = * m_bases[0][patchIndex];
    const index_t Ndofs = basis.size();

    gsVector<int> numNodes = gsAssemblerUtils<T>::getNumIntNodesFor( basis );

    // Empty sparse matrix to fill with the contribution of the current patch
    gsSparseMatrix<T> * Kh = new gsSparseMatrix<T>(Ndofs,Ndofs );

    gsMatrix<T> trf_grads_k;    // transformed (physical) gradients of all basis functions at one quadrature node
    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
        this->m_patches.patch(patchIndex).evaluator(NEED_VALUE   |
                                                    NEED_MEASURE |
                                                    NEED_GRAD_TRANSFORM) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->activeFuncs.rows();
        domIt->evaluateBasis( 1 );

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(domIt->quNodes);

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = domIt->quWeights[k] * geoEval->measure(k);

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval->transformGradients(k, domIt->basisDerivs(1), trf_grads_k);

            localStiffness.noalias() += weight * (trf_grads_k.transpose() * trf_grads_k);
        }  // end loop Gauss nodes

        // add contributions from local stiffness matrix to global one.
        // Note: assembled as FULL stiffness matrix by setting "symmetric" to false.
        gsAssemblerUtils<T>::localToGlobal( localStiffness, domIt->activeFuncs, * Kh, false);

    } //end loop over all domain elements

    Kh->makeCompressed();
    return Kh;
}


// S.Kleiss
// Assembles and returns the raw mass matrix without any boundary conditions
// Not really tested yet (07.May 2014).
template<class T>
gsSparseMatrix<T> * gsConvDiffReactAssembler<T>::massMatrixPatch(int patchIndex)
{
    const gsBasis<T> & basis = * m_bases[0][patchIndex];
    const index_t Ndofs = basis.size();

    gsVector<int> numNodes = gsAssemblerUtils<T>::getNumIntNodesFor( basis );

    // Empty sparse matrix to fill with the contribution of the current patch
    gsSparseMatrix<T> * Mh = new gsSparseMatrix<T>(Ndofs,Ndofs );

    gsMatrix<T> localMass; // (dense) mass matrix within one grid cell
    gsMatrix<T> BasisFcts_k;

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
        this->m_patches.patch(patchIndex).evaluator(NEED_VALUE   |
                                                    NEED_MEASURE) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();

    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->activeFuncs.rows();
        domIt->evaluateBasis( 0 );

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(domIt->quNodes);

        // initialize local linear system to 0
        localMass.setZero(numActive, numActive);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = domIt->quWeights[k] * geoEval->measure(k);

            BasisFcts_k = domIt->basisValues().col(k);
            localMass.noalias() += weight * ( BasisFcts_k * BasisFcts_k.transpose());
        }  // end loop Gauss nodes

        // add contributions from local mass matrix to global one.
        // Note: assembled as FULL mass matrix by setting "symmetric" to false.
        gsAssemblerUtils<T>::localToGlobal( localMass, domIt->activeFuncs, * Mh, false);

    } //end loop over all domain elements

    Mh->makeCompressed();
    return Mh;
}


// S.Kleiss
// Should become main function to call error indication.
// Currently just calling the one with bubble functions
// TODO: Add a way of selecting one out of several error indicators
template<class T>
gsVector< gsMatrix<T> > gsConvDiffReactAssembler<T>::errorIndicator()
{
    gsVector< gsMatrix<T> > errInd;
    errInd.resize( this->m_patches.nPatches() );

    // For the time being, only one error indicator implemented
    // TODO: Add other error indicators and
    // add a way of selecting one of them.
    for( unsigned i=0; i<this->m_patches.nPatches(); i++)
    {
        //errInd[i] = errIndPatch_Bubble( i );
        errInd[i] = errIndPatch_Residual( i );
    }

    return errInd;
} //end gsConvDiffReactAssembler<T>::errorIndicator()


template<class T>
gsMatrix<T> gsConvDiffReactAssembler<T>::errIndPatch_Residual( const index_t patchIndex )
{
    gsDebug << "Start: errIndPatch_Residual()\n";

    // Comptes the local residuals in a straight-forward manner on an element K
    // size( K ) \int_K ( f + div (A grad uh) - b (grad uh) - c uh )^2 dx
    // WARNING: THE 2nd-ORDER TERM IS NOT TRANSFORMED TO THE PHYSICAL DOMAIN YET!!!

//    std::cout << "gsConvDiffReactAssembler<T>::errIndPatch_Residual\n";
//    std::cout << "WARNING: THE 2nd-ORDER TERM IS NOT TRANSFORMED ";
//    std::cout << "TO THE PHYSICAL DOMAIN YET!!!\n\n";

    // access the local basis on patch:
    const gsBasis<T> & basis = * m_bases[0][patchIndex];
    // construct domain iterator for basis:
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // number of elements on patch:
    const index_t NE = basis.numElements();
    // dimension of the parameter space of the basis:
    const index_t d = basis.dim();

    // access local solution...
    gsField<T> * sol = m_solutions[0];

    // ...and extract the "discrete-solution-part".
    const gsGeometry<T> * igaFct = & sol->igaFunction(patchIndex);

    // initialize the matrix for the local errors and
    // element coordinates
    gsMatrix<T> errInd( NE, 1 + 2*d);

    // set up gsGeometryEvaluator
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
                this->m_patches.patch(patchIndex).evaluator(NEED_VALUE   |
                                                            NEED_MEASURE |
                                                            NEED_GRAD_TRANSFORM) );

    gsMatrix<T> trf_uhGradsT;   // transformed (physical) gradients of discrete solution

    gsMatrix<T> rhsVals;        // gradients of discrete solution

    gsMatrix<T> uhVals;         // gradients of discrete solution
    gsMatrix<T> uhGrads;        // gradients of discrete solution
    gsMatrix<T> uh2ndDerivs;        // gradients of discrete solution

    gsMatrix<T> coeff_A_vals;
    gsMatrix<T> coeff_b_vals;
    gsMatrix<T> coeff_c_vals;

    // compute ...
    // ...the number of quadrature nodes
    gsVector<int> numQuNodes(d);
    for( index_t i=0; i < d; i++)
        numQuNodes[i] = basis.degree(i) + 1;

    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numQuNodes );

    //std::ofstream out("../../matlab_Gismo/uh140708a.m");
    //out << "D = [...\n";

    index_t actElt = 0;
    for (; domIt->good(); domIt->next(), actElt++ )
    {
        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt( domIt->quNodes );

        gsMatrix<T> Jac = geoEval->jacobians();
        GISMO_ASSERT( Jac.rows() == d ,"m1031, just to check");

        const gsVector<T> lowC = domIt->lowerCorner();
        const gsVector<T> uppC = domIt->upperCorner();
        const gsVector<T> h = uppC - lowC;
        T hh = T(1.0);
        for( int i=0; i < d; i++)
            hh *= h[i];

        m_coeff_A->eval_into( geoEval->values(), coeff_A_vals );
        m_coeff_b->eval_into( geoEval->values(), coeff_b_vals );
        m_coeff_c->eval_into( geoEval->values(), coeff_c_vals );

        m_rhs_function->eval_into( geoEval->values(), rhsVals );

        int quN = rhsVals.cols();

        igaFct->eval_into( domIt->quNodes, uhVals );

        GISMO_ASSERT( uhVals.cols() == quN, "Mark 1119" );
        GISMO_ASSERT( domIt->numQuNodes() == quN, "Mark 1120" );

        // compute gradients of discrete solution and resize result
        igaFct->deriv_into( domIt->quNodes, uhGrads );
        gsMatrix<T> uhGradsT( d, quN );
        for( int i=0; i < d; i++)
            for( int j=0; j < quN; j++)
                uhGradsT(i,j) = uhGrads(0,j*d+i);

        igaFct->deriv2_into( domIt->quNodes, uh2ndDerivs );

        T errIndLocal = T(0.0);

        for (index_t pp = 0; pp < quN; ++pp)      // loop over quadrature nodes
        {
            T contrib_pp = T(0.0);

            gsMatrix<T> Jac_pp( d,d );
            for ( int i = 0; i < d; i++)
                Jac_pp.col(i) = Jac.col( pp*d + i );

            // computing the inverse of the Jacobian here, because
            // at the moment (10.Jul.2014), it is not so clear what
            // is stored at m_Jinv of the gsGeometryEvaluator.
            gsMatrix<T> JacInv_pp;
            JacInv_pp.setZero(d,d);
            if( d == 2 )
            {
                T tmpDet = Jac_pp(0,0) * Jac_pp(1,1) - Jac_pp(0,1) * Jac_pp(1,0);
                JacInv_pp(0,0) = Jac_pp(1,1) / tmpDet;
                JacInv_pp(1,1) = Jac_pp(0,0) / tmpDet;
                JacInv_pp(1,0) = (-1.0) * Jac_pp(1,0) / tmpDet;
                JacInv_pp(0,1) = (-1.0) * Jac_pp(0,1) / tmpDet;
            }
            else
            {
                JacInv_pp = Jac_pp.inverse();
            }

            gsMatrix<T> trf_uh2ndDerivs( d, d );
            trf_uh2ndDerivs.setZero();

            // update quadrature weights with determinant of Jacobian.
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = domIt->quWeights[pp] * geoEval->measure(pp);

            // compute physical gradients of discrete solution
            geoEval->transformGradients(pp, uhGradsT, trf_uhGradsT);

            // getting a map from double-indices to a single index
            // referring to the mixed derivatives, which are stored in
            // uh2ndDerivs in non-lexicographic order.
            gsMatrix<int> klM( d, d );
            if( d == 2 )
            {
                klM(0,0) = 0;
                klM(0,1) = 2;
                klM(1,0) = 2;
                klM(1,1) = 1;
            }
            else if( d == 3 )
            {
                klM(0,0) = 0;
                klM(0,1) = 3;
                klM(0,2) = 4;
                klM(1,0) = 3;
                klM(1,1) = 1;
                klM(1,2) = 5;
                klM(2,0) = 4;
                klM(2,1) = 5;
                klM(2,2) = 2;
            }

            // transforming second derivatives to the physical domain
            // WARNINGS:
            // ! Assuming (piecewise) constant coefficient A
            // ! Neglecting the term with the second derivatives of the inverse geometric mapping!
            for( int i = 0; i < d; i++)
                for( int j = 0; j < d; j++)
                    for( int k = 0; k < d; k++)
                        for( int l = 0; l < d; l++)
                        {
                            int kl = klM(k,l);
                            trf_uh2ndDerivs( i,j ) += uh2ndDerivs(kl,pp) * JacInv_pp(l,j) * JacInv_pp(k,i);
                        }

            // add div A grad uh... doing this manually to be on the safe side
            T divAgradUh = T(0.0);
            if( d == 2 )
            {
                divAgradUh += coeff_A_vals(0,pp) * trf_uh2ndDerivs(0,0); // d_xx
                divAgradUh += coeff_A_vals(1,pp) * trf_uh2ndDerivs(0,1); // d_xy
                divAgradUh += coeff_A_vals(2,pp) * trf_uh2ndDerivs(1,0); // d_yx
                divAgradUh += coeff_A_vals(3,pp) * trf_uh2ndDerivs(1,1); // d_yy
            }
            else if( d == 3 )
            {
                divAgradUh += coeff_A_vals(0,pp) * trf_uh2ndDerivs(0,0); // d_xx
                divAgradUh += coeff_A_vals(1,pp) * trf_uh2ndDerivs(0,1); // d_xy
                divAgradUh += coeff_A_vals(2,pp) * trf_uh2ndDerivs(0,2); // d_xz
                divAgradUh += coeff_A_vals(3,pp) * trf_uh2ndDerivs(1,0); // d_yx
                divAgradUh += coeff_A_vals(4,pp) * trf_uh2ndDerivs(1,1); // d_yy
                divAgradUh += coeff_A_vals(5,pp) * trf_uh2ndDerivs(1,2); // d_yz
                divAgradUh += coeff_A_vals(6,pp) * trf_uh2ndDerivs(2,0); // d_zx
                divAgradUh += coeff_A_vals(7,pp) * trf_uh2ndDerivs(2,1); // d_zy
                divAgradUh += coeff_A_vals(8,pp) * trf_uh2ndDerivs(2,2); // d_zz
            }
            contrib_pp += divAgradUh;

            contrib_pp -= coeff_b_vals(0,pp) * trf_uhGradsT(0,0);
            contrib_pp -= coeff_b_vals(1,pp) * trf_uhGradsT(1,0);
            if( d == 3 )
                contrib_pp -= coeff_b_vals(2,pp) * trf_uhGradsT(2,0);

            contrib_pp -= coeff_c_vals(0,pp) * uhVals(0,pp);

            // add right-hand-side
            contrib_pp += rhsVals(0,pp);

            errIndLocal += weight * contrib_pp * contrib_pp;

        }  // end loop Gauss nodes

        // store result and coordinates of element in parameter domain
        errInd( actElt, 0 ) = hh * errIndLocal;
        for( index_t i = 0; i < d; i++)
        {
            errInd( actElt, 1 + i) = lowC[i];
            errInd( actElt, 1 + d + i) = uppC[i];
        }

    } // domIt

    //out << "];\n";
    //out << "NE = " << errInd.rows() << ";\n";

    gsDebug << "--End: errIndPatch_Residual()\n";

    return errInd;
}

// S.Kleiss
// Uses error indicator with bubble functions.
// Based on error estimator presented in
// M. R. Doerfel, B. Juettler, B. Simeon. Adaptive isogeometric analysis by local
// h-refinement with T-splines.
// Computer Methods in Applied Mechanics and Engineering, 199: 264-275, 2010.
//
// Using a version described in Section 3.8 of
// V. John. A numerical study of a posteriori error estimators for
// convection-diffusion equations.
// Computer Methods in Applied Mechanics and Engineering, 190: 757-781, 2000.
template<class T>
gsMatrix<T> gsConvDiffReactAssembler<T>::errIndPatch_Bubble( const index_t patchIndex )
{

    gsDebug << "Start: errIndPatch_Bubble()\n";

    // access the local basis on patch:
    const gsBasis<T> & basis = * m_bases[0][patchIndex];
    // construct domain iterator for basis:
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // number of elements on patch:
    const index_t NE = basis.numElements();
    // dimension of the parameter space of the basis:
    const index_t d = basis.dim();
    GISMO_ASSERT( d == 2 || d == 3 , "CHECK DIMENSION in errIndPatch_Bubble!" );
    // number of second derivatives
    const index_t d2d = d * ( d+1 ) / 2;

    // access local solution...
    const gsField<T> * sol = m_solutions[patchIndex];
    // ...and extract the "discrete-solution-part", and...
    const gsGeometry<T> & igaFct = sol->igaFunction(patchIndex);
    // ...the "geometry-information-part".
    //const gsGeometry<T> * geo = sol->patch(patchIndex);

    //std::cout << " igaFct.m_coefs\n" << igaFct.m_coefs << "\n";

    // initialize the matrix for the local errors and
    // element coordinates
    gsMatrix<T> errInd( NE, 1 + 2*d );

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
        gsKnotVector<T> bub_KVtmp(0,1,0, bubbleDeg[i]+1 );
        bub_KV[i] = bub_KVtmp;
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
    domIt->computeQuadratureRule( numQuNodes );
    // ...and again separately, because I can't get the
    // quadrature rule out of the domain iterator.
    gsGaussRule<T> QuRule( numQuNodes );

    // could be more efficient, but good enough for now:

    // Evaluate the bubbles once on the reference hypercube,
    // thus obtaining function values, first and second derivatives
    gsVector<T> v0(d);
    gsVector<T> v1(d);
    v0.setZero();
    v1.setOnes();
    gsMatrix<T> Qp;
    gsVector<T> Qw;
    QuRule.mapTo( v0,v1, Qp, Qw);
    gsMatrix<T> bub_v_tmp = bubbles->eval( Qp );
    gsMatrix<T> bub_d_tmp = bubbles->deriv( Qp );
    gsMatrix<T> bub_d2_tmp = bubbles->deriv2( Qp );

    gsMatrix<T> bub_v( bubInt_nTotal, bub_v_tmp.cols() );
    gsMatrix<T> bub_d( d * bubInt_nTotal, bub_d_tmp.cols() );
    gsMatrix<T> bub_d2( d*(d+1)/2 * bubInt_nTotal, bub_d_tmp.cols() );

    for( index_t i = 0; i < bubInt_indices.size(); i++)
    {
        index_t ii = bubInt_indices[i];
        bub_v.row(i) = bub_v_tmp.row(ii);
        bub_d.row(d*i    ) = bub_d_tmp.row(d*ii);
        bub_d.row(d*i + 1) = bub_d_tmp.row(d*ii + 1);

        bub_d2.row(d2d*i     ) = bub_d2_tmp.row( d2d*ii );
        bub_d2.row(d2d*i + 1 ) = bub_d2_tmp.row( d2d*ii + 1 );
        bub_d2.row(d2d*i + 2 ) = bub_d2_tmp.row( d2d*ii + 2 );

        if( d == 3)
        {
            bub_d.row(d*i+2) = bub_d_tmp.row(d*ii+2);

            bub_d2.row(d2d * i+3 ) = bub_d2_tmp.row( d2d * ii+3 );
            bub_d2.row(d2d * i+4 ) = bub_d2_tmp.row( d2d * ii+4 );
            bub_d2.row(d2d * i+5 ) = bub_d2_tmp.row( d2d * ii+5 );
        }
    }

    bub_v_tmp.resize( bub_v.rows(), bub_v.cols() );
    bub_d_tmp.resize( bub_d.rows(), bub_d.cols() );
    bub_d2_tmp.resize( bub_d2.rows(), bub_d2.cols() );
    bub_v_tmp.setZero();
    bub_d_tmp.setZero();
    bub_d2_tmp.setZero();

    // set up gsGeometryEvaluator
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval (
                this->m_patches.patch(patchIndex).evaluator(NEED_VALUE   |
                                                            NEED_MEASURE |
                                                            NEED_GRAD_TRANSFORM) );

    gsMatrix<T> rhsVals;       // values of the right-hand side
    gsMatrix<T> uhVals;        // values of discrete solution
    gsMatrix<T> uhGrads;       // gradients of discrete solution
    gsMatrix<T> trf_grads_pp;  // transformed (physical) gradients of all basis functions at one quadrature node
    gsMatrix<T> trf_uhGradsT;  // transformed (physical) gradients of discrete solution
    gsMatrix<T> b_grads_pp;    // for computation of b.( grad basis-function)
    gsMatrix<T> b_gradsUh_pp;  // for computation of b.( grad uh)

    gsMatrix<T> coeff_A_vals;
    gsMatrix<T> coeff_b_vals;
    gsMatrix<T> coeff_c_vals;
    gsMatrix<T> tmp_A;

    T SUPG_tau = T(0);
    gsMatrix<T> SUPG_termA_pp;

    Eigen::Matrix<T,Dynamic,Dynamic> localMat( bubInt_nTotal, bubInt_nTotal); // local system matrix
    Eigen::Matrix<T,Dynamic,Dynamic> localRhs( bubInt_nTotal, 1 );            // local right-hand-side
    Eigen::Matrix<T,Dynamic,Dynamic> localMatSUPG( bubInt_nTotal, bubInt_nTotal); // local system matrix
    Eigen::Matrix<T,Dynamic,Dynamic> localRhsSUPG( bubInt_nTotal, 1 );        // local right-hand-side
    // matrix for H1-products of basis functions...
    Eigen::Matrix<T,Dynamic,Dynamic> H1Mat( bubInt_nTotal, bubInt_nTotal);

    index_t actElt = 0;
    for (; domIt->good(); domIt->next(), actElt++ )
    {
        localMat.setZero();
        localRhs.setZero();
        localMatSUPG.setZero();
        localRhsSUPG.setZero();
        H1Mat.setZero();

        const gsVector<T> lowC = domIt->lowerCorner();
        const gsVector<T> uppC = domIt->upperCorner();
        const gsVector<T> h = uppC - lowC;

        bub_v_tmp = bub_v;
        // adjust the derivatives of the bubble functions
        // to the current element.
        bub_d_tmp = bub_d;
        for( index_t i = 0; i < d ; i++)
        {
            for( index_t j = 0; j < bub_v.rows(); j++)
                bub_d_tmp.row( j*d + i) /= h[i];
        }

        // adjust the 2nd derivatives of the bubble functions
        // to the current element.
        bub_d2_tmp = bub_d2;
        gsVector<T> hh( d2d );
        hh[0] = h[0]*h[0];
        hh[1] = h[1]*h[1];
        if( d == 2 )
        {
            hh[2] = h[0]*h[1];
        }
        else if( d == 3 )
        {
            hh[2] = h[2]*h[2];
            hh[3] = h[0]*h[1];
            hh[4] = h[0]*h[2];
            hh[5] = h[1]*h[2];
        }

        for( index_t i = 0; i < d2d; i++)
        {
            for( index_t j = 0; j < bub_v.rows(); j++)
                bub_d2_tmp.row( j*d2d + i) /= hh[i];
        }

        // inefficient, but more convenient so that the
        // code from assemblePatch can be copy-pasted.
        // "basis" in those variables now corresponds to bubble functions!
        const gsMatrix<T>basisValues    = bub_v_tmp;
        const gsMatrix<T>basisGrads     = bub_d_tmp;
        const gsMatrix<T>basis2ndDerivs = bub_d2_tmp;

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(domIt->quNodes);

        // compute values of discrete solution
        igaFct.eval_into( domIt->quNodes, uhVals );
        // compute gradients of discrete solution and resize result
        igaFct.deriv_into( domIt->quNodes, uhGrads );
        gsMatrix<T> uhGradsT(d, bub_v.cols() );
        for( int i=0; i < d; i++)
            for( int j=0; j < bub_v.cols(); j++)
                uhGradsT(i,j) = uhGrads(0,j*d+i);

        // evaluate right-hand side at the geometry points
        m_rhs_function->eval_into( geoEval->values(), rhsVals );

        m_coeff_A->eval_into( geoEval->values(), coeff_A_vals );
        m_coeff_b->eval_into( geoEval->values(), coeff_b_vals );
        m_coeff_c->eval_into( geoEval->values(), coeff_c_vals );

        const int d2_offsetter = d*(d+1)/2;

        for (index_t pp = 0; pp < domIt->numQuNodes(); ++pp)      // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = domIt->quWeights[pp] * geoEval->measure(pp);

            // compute physical gradients at quadrature point k
            // d... dimension of parameter space
            // N... number of active functions

            // compute physical gradients of discrete solution
            geoEval->transformGradients(pp, uhGradsT, trf_uhGradsT);

            // trf_grads_pp:            d x N
            geoEval->transformGradients(pp, basisGrads, trf_grads_pp);
            int N = trf_grads_pp.cols();

            // coeff_b_vals.col(pp):    d x 1
            // trf_grads_pp:            d x N
            // b_grads_pp:              1 x N
            b_grads_pp = coeff_b_vals.col(pp).transpose() * trf_grads_pp;

            // b_gradsUh_pp:            1 x 1
            b_gradsUh_pp = coeff_b_vals.col(pp).transpose() * trf_uhGradsT;

            tmp_A = coeff_A_vals.col(pp);
            tmp_A.resize(d,d);

            // "usual" right-hand-side term
            localRhs += basisValues.col(pp) * (weight * rhsVals.col(pp)).transpose();

            // Terms coming from a( uh, bubble )
            // --- 2-nd order term
            // trf_grads_pp:           N x d
            // tmp_A * trf_uhGradsT:   d x 1
            localRhs -= weight * ( trf_grads_pp.transpose() * (tmp_A * trf_uhGradsT));

            // --- 1-st order term
            // basisValues.col(pp):  N x 1
            // b_gradsUh_pp:         1 x 1
            localRhs -= weight * ( basisValues.col(pp) * b_gradsUh_pp );

            // --- 0-th order term
            localRhs -= weight * coeff_c_vals(0,pp) * ( basisValues.col(pp) * (uhVals.col(pp).transpose() ) );

            // --- 2nd-order term
            // trf_grads_pp.transpose(): N x d
            // A:                        d x d
            // trf_grads_pp:             d x N
            localMat.noalias() += weight * (trf_grads_pp.transpose() * (tmp_A * trf_grads_pp) );

            H1Mat.noalias() += weight * ( trf_grads_pp.transpose() * trf_grads_pp );

            // --- 1st-order term
            // b_grads_pp:              1 x N
            // basisValues.col(pp):     N x 1
            localMat.noalias() += weight * ( basisValues.col(pp) * b_grads_pp );

            // --- 0th-order term
            // coeff_c_vals.col(pp):    1 x 1
            // basisValues.col(pp):     N x 1
            localMat.noalias() += weight * coeff_c_vals(0,pp) * ( basisValues.col(pp) * (basisValues.col(pp).transpose() ) );

            if( m_doSUPG )
            {
                // SUPG_termA_pp:            d x N
                SUPG_termA_pp.resize( d, N );
                SUPG_termA_pp.setZero();

                SUPG_tau = get_SUPG_parameter( domIt->lowerCorner(), domIt->upperCorner(), patchIndex);

                // Inverse of Jacobian...
                // --- WARNING:
                // The function gradTransforms() returns the
                // TRANSPOSED Jacobian at the moment (01.Jul.2014)!
                gsMatrix<T> JinvT = geoEval->gradTransforms();

                gsMatrix<int> klMatrix(d,d);
                klMatrix.setZero();
                if( d == 2 )
                    klMatrix << 0,2,2,1;
                else if ( d == 3 )
                    klMatrix << 0,3,4, 3,1,5, 4,5,2;

//                if( d == 2 )
//                {
//                    klMatrix(0,0) = 0;
//                    klMatrix(0,1) = 2;
//                    klMatrix(1,0) = 2;
//                    klMatrix(1,1) = 1;
//                }
//                else if( d == 3 )
//                {
//                    klMatrix(0,0) = 0;
//                    klMatrix(0,1) = 3;
//                    klMatrix(0,2) = 4;
//                    klMatrix(1,0) = 3;
//                    klMatrix(1,1) = 1;
//                    klMatrix(1,2) = 5;
//                    klMatrix(2,0) = 4;
//                    klMatrix(2,1) = 5;
//                    klMatrix(2,2) = 2;
//                }

                SUPG_termA_pp.setZero();
                for( int ni=0; ni < N ; ni++)
                {
                    for( int dj=0; dj < d; dj++)
                        for( int i=0; i < d; i++ )
                        {
                            T tmp = T(0);

                            for( int k = 0; k < d; k++)
                                for( int l = 0; l < d; l++)
                                {
//                                    int kl = 2*ni*d2_offsetter; // default will throw an error
//                                    if( k == l )
//                                        kl = k;
//                                    else if( (k==0 && l==1) || (k==1 && l==0 ) )
//                                        kl = d;
//                                    else if( (k==0 && l==2) || (k==2 && l==0 ) )
//                                        kl = 4;
//                                    else if( (k==1 && l==2) || (k==2 && l==1 ) )
//                                        kl = 5;
                                    int kl = klMatrix(k,l);

                                    // WORKING WITH THE TRANSPOSED INVERSE JinvT
                                    tmp += basis2ndDerivs( ni*d2_offsetter + kl, pp )
                                            * JinvT( i, pp*d + k ) * JinvT( dj, pp*d + l );

                                    // Version where Jinv really is the inverse of the Jacobian
                                    // tmp += basis2ndDerivs(ni*d2_offsetter+kl,pp)
                                    //         * Jinv( k, ni*d + i ) * Jinv( l, ni*d + dj );
                                }

                            SUPG_termA_pp(dj,ni) += coeff_b_vals(i,pp) * tmp;
                        }
                }

                // --- SUP 2nd-order term is NOT added (see [V. John])
                // --- SUPG 1st-order term
                localMatSUPG.noalias() += weight * ( b_grads_pp.transpose() * b_grads_pp );
                // --- SUPG 0th-order term
                localMatSUPG.noalias() += ( weight * coeff_c_vals(0,pp) )*( b_grads_pp.transpose() * ( basisValues.col(pp).transpose() ));

                // --- SUPG right-hand-side term with uh
                localRhsSUPG.noalias() += weight * ( b_grads_pp.transpose() * (rhsVals.col(pp)).transpose() );
                // --- SUPG 2nd-order term with uh
                localRhsSUPG.noalias() -= weight * ( SUPG_termA_pp.transpose() * ( tmp_A * trf_uhGradsT ) );
                // --- SUPG 1st-order term with uh
                localRhsSUPG.noalias() -= weight * ( b_grads_pp.transpose() * b_gradsUh_pp );
                // --- SUPG 0th-order term with uh
                localRhsSUPG.noalias() -= (weight * coeff_c_vals(0,pp) )*( b_grads_pp.transpose() * ( uhVals.col(pp).transpose() ) );

            } // SUPG

        }  // end loop Gauss nodes

        if( m_doSUPG )
        {
            localMat.noalias() += SUPG_tau * localMatSUPG;
            localRhs.noalias() += SUPG_tau * localRhsSUPG;
        }

        // apply direct solver to these small local problems
        Eigen::Matrix<double,Dynamic,Dynamic> eh = localMat.fullPivLu().solve( localRhs );
        //std::cout << "\n" << eh << "\n";

        // recover the error indicator
        T errIndLocal = 0.0;
        for( index_t i = 0; i < eh.size(); i++)
            for( index_t j = 0; j < eh.size(); j++)
            {
                errIndLocal += eh(i,0) * eh(j,0) * localMat(j,i); // see [Doerfel et al]
                //errIndLocal += eh(i,0) * eh(j,0) * H1Mat(j,i); // see [V. John]
            }

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

} // end gsConvDiffReactAssembler<T>::errorIndicatorBubble()


// S.K.
template<class T>
void gsConvDiffReactAssembler<T>::adaptiveRefine( const gsVector< gsMatrix<T> >& errInd, const int refCriterion, const T a, const unsigned RefExtension)
{
    gsDebug << "Start: adaptiveRefine()\n";
    gsDebug << "         Criterion " << refCriterion << "\n";

    std::cout << "Start: adaptiveRefine()\n";
    std::cout << "         Criterion " << refCriterion << std::endl;

    std::cout << "" << std::endl;

    T Thr = 0;
    // Total number of elements
    int NE = 0;

    // Checking whether the underlying
    // basis is an HTensorBasis.
    // Will be used further down below.
    const gsHTensorBasis<2,T> * HTensB2cast = dynamic_cast< const gsHTensorBasis<2,T> * >( & m_bases[0].basis(0) );
    const gsHTensorBasis<3,T> * HTensB3cast = dynamic_cast< const gsHTensorBasis<3,T> * >( & m_bases[0].basis(0) );

    bool Flag_BasisIsHTensB = false;
    if( HTensB2cast != 0 || HTensB3cast != 0 )
        Flag_BasisIsHTensB = true;


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
    const int d = m_bases[0][0]->dim();

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


        if( Flag_BasisIsHTensB == false || RefExtension == 0)
        {
            m_bases[0][ patchIndex ]->refine( RefBox );
        }
        else
        {
            // If the refinement areas is to be extended
            // for a gsHTensorBasis, we will use the refinement
            // function refine( std::vector ) where one can explicitly give the
            // level to which the are should be refined to.
            //
            // If the other refinement function, refine(gsMatrix),
            // is used, then overlapping extensions might lead
            // to cells which are refined more than once, which would
            // mess things up.


            std::vector<unsigned> RefVec;

            // transform every box given by parameter values in the
            // gsMatrix RefBox to levels and unique knot-indices.
            for( int j = 0; j < RefBox.cols(); j+= 2 )
            {
                gsMatrix<double> mid(d,1);

                // find the point in the lower left (front) quarter
                // of the cell.
                for( int di=0; di < d; di++ )
                    mid(di,0) = 3*RefBox(di,j)/4 + RefBox(di,j+1)/4;

                // get the level of the mesh at this particular point.
                gsMatrix<int> level(1,1);
                if( HTensB2cast != 0 )
                {
                    level = HTensB2cast->getLevelAtPoint( mid );
                }
                else if( HTensB3cast != 0 )
                {
                    level = HTensB3cast->getLevelAtPoint( mid );
                }

                std::vector<unsigned> tmpRefVec(1+2*d);
                tmpRefVec[0] = level(0,0)+1;

                // find out, which cell of the finer mesh this point
                // is contained in.
                for( int di =0; di < d; di++)
                {
                    unsigned tmpKnotIdx = HTensB2cast->m_bases[level(0,0)+1]->component(di).knots().Uniquefindspan( mid(di,0) );

                    // lower corner of the refinement area
                    if ( tmpKnotIdx < RefExtension )
                        tmpRefVec[di+1] = 0;
                    else
                        tmpRefVec[di+1] = tmpKnotIdx - RefExtension;

                    // for the upper corner, add 1 to get the next index
                    // and another 1, because the indices already correspond
                    // to the next finer level.
                    // And then, of course, add the RefExtension.
                    tmpRefVec[di+d+1] = tmpKnotIdx + 2 + RefExtension;
                }

                for( unsigned i = 0; i< tmpRefVec.size(); i++)
                    RefVec.push_back( tmpRefVec[i] );

            } // for j

            m_bases[0][ patchIndex ]->refineElements( RefVec );

        } // if is HTensB

    } // patchIndex

    gsDebug << "  End: adaptiveRefine()\n";
} // end gsConvDiffReactAssembler<T>::adaptiveRefinement


template<class T>
void gsConvDiffReactAssembler<T>::uniformRefine()
{
    for (unsigned np=0; np < this->m_patches.nPatches(); ++np )
    {
        m_bases[0][np]->uniformRefine();
    }
}



} // namespace gismo


