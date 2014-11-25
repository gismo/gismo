
#pragma once

#include <gsUtils/gsQuadrature.h>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsDebug.h>
#include <gsCore/gsDomain.h>

#include <gsUtils/gsCombinatorics.h>

#include <gsTensor/gsTensorDomainBoundaryIterator.h>

#include <gsCore/gsBoundary.h>

namespace gismo {


template<class T>
gsSparseMatrix<T>*
gsGaussAssembler<T>::assembleMass( const gsBasis<T>& basis, const gsDofMapper& mapper, int patchIndex)
{
    gsVector<int> numNodes = getNumIntNodesFor( basis );
    gsSparseMatrix<T> *M = this->initMatrix( mapper.freeSize(), basis, true );
    gsMatrix<T> localMass;
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( this->geometry().evaluator(NEED_MEASURE) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->computeActiveDofs(mapper, patchIndex).rows();
        domIt->evaluateBasis( 0 );

        geoEval->evaluateAt(domIt->quNodes);

        localMass.setZero(numActive, numActive);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = domIt->quWeights[k] * geoEval->measure(k); 
            localMass.noalias() += weight * (domIt->basisValues().col(k) * domIt->basisValues().col(k).transpose());
        }

        this->localToGlobal_withBC(localMass, mapper, domIt->activeDofs, *M, true);
        //this->localToGlobal_withBC(localMass, mapper, domIt->activeDofs, *M, false);
    }

    M->makeCompressed();
    return M;
}


template<class T>
gsSparseSystem<T>
gsGaussAssembler<T>::assemblePoisson( const gsBasis<T>& basis, const gsDofMapper& mapper, const gsMatrix<T> & ddof, const gsPoissonPde<T> & pde, int patchIndex)
{
    gsVector<int> numNodes = getNumIntNodesFor( basis );
    gsSparseSystem<T> sys = this->initLinearSystem( mapper.freeSize(), basis, pde.isSymmetric() );
    
    gsMatrix<T> rhsVals;        // values of the right-hand side
    gsMatrix<T> trf_grads_k;    // transformed (physical) gradients of all basis functions at one quadrature node
    
    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell
    gsVector<T> localRhs;       // local load vector

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( this->geometry().evaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->computeActiveDofs(mapper, patchIndex).rows();
        domIt->evaluateBasis( 1 );

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(domIt->quNodes);

        // evaluate right-hand side at the geometry points
        if (pde.rhs())
            pde.rhs()->eval_into( geoEval->values(), rhsVals );

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);
        localRhs.setZero(numActive);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            const T weight = domIt->quWeights[k] * geoEval->measure(k); // weight * abs(det J), where J is geometry Jacobian

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval->transformGradients(k, domIt->basisDerivs(1), trf_grads_k);

            if (pde.rhs())
                localRhs += (weight * rhsVals(0,k)) * domIt->basisValues().col(k);

            localStiffness.noalias() += weight * (trf_grads_k.transpose() * trf_grads_k);
        }  // end loop Gauss nodes

        // add contributions from local stiffness matrix to global
        // stiffness matrix and load vector
        this->localToGlobal_withBC(localStiffness, localRhs, ddof, mapper, 
                                   domIt->activeDofs, *sys.matrix(), *sys.rhs(), true);

    } //end loop over all domain elements

    sys.matrix()->makeCompressed();
    return sys;
}

template<class T>
gsSparseSystem<T>
gsGaussAssembler<T>::assembleCDR( const gsBasis<T>& basis, const gsDofMapper& mapper, const gsMatrix<T> & ddof, const gsConvDiffRePde<T> & pde, int patchIndex)
{
    gsVector<int> numNodes = getNumIntNodesFor( basis );
    gsSparseSystem<T> sys = this->initLinearSystem( mapper.freeSize(), basis, pde.isSymmetric() );
    
    gsMatrix<T> rhsVals;        // values of the right-hand side
    gsMatrix<T> diffVals, convVals, reacVals;   // PDE coefficient values
    gsMatrix<T> trf_grads_k;    // transformed (physical) gradients of all basis functions at one quadrature node

    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell
    gsVector<T> localRhs;       // local load vector

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( this->geometry().evaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->computeActiveDofs(mapper, patchIndex).rows();
        domIt->evaluateBasis( 1 );

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(domIt->quNodes);

        // evaluate right-hand side and PDE coefficients at the geometry points
        if (pde.rhs())              pde.rhs()->eval_into( geoEval->values(), rhsVals );
        if (pde.diffusion())        pde.diffusion()->eval_into( geoEval->values(), diffVals );
        if (pde.convection())       pde.convection()->eval_into( geoEval->values(), convVals );
        if (pde.reaction())         pde.reaction()->eval_into( geoEval->values(), reacVals );

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);
        localRhs.setZero(numActive);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            const T weight = domIt->quWeights[k] * geoEval->measure(k); // weight * abs(det J), where J is geometry Jacobian

            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval->transformGradients(k, domIt->basisDerivs(1), trf_grads_k);

            if (pde.rhs())
                localRhs += (weight * rhsVals(0,k)) * domIt->basisValues().col(k);

            if (pde.diffusion())
                localStiffness.noalias() += (weight * diffVals(0,k)) * (trf_grads_k.transpose() * trf_grads_k);

            if (pde.convection())
                localStiffness.noalias() += weight * (domIt->basisValues().col(k) * (convVals.col(k).transpose() * trf_grads_k));

            if (pde.reaction())
                localStiffness.noalias() += (weight * reacVals(0,k)) * (domIt->basisValues().col(k) * domIt->basisValues().col(k).transpose());
        }  // end loop Gauss nodes

        // add contributions from local stiffness matrix to global stiffness matrix and load vector
        this->localToGlobal_withBC(localStiffness, localRhs, ddof, mapper, domIt->activeDofs, *sys.matrix(), *sys.rhs(), pde.isSymmetric());
    } // loop over all domain elements

    sys.matrix()->makeCompressed();
    return sys;
}


template<class T>
gsSparseSystem<T>
gsGaussAssembler<T>::assembleBeam( const gsBasis<T>& basis, const gsDofMapper& mapper, const gsMatrix<T> & ddof, const gsEulerBernoulliBeamPde<T> & pde, int patchIndex)
{
    const int d = this->geometry().parDim();
    if (d != 1)
        GISMO_ERROR("Beam must be 1D");
    gsVector<int> numNodes = getNumIntNodesFor( basis );

    gsSparseSystem<T> sys = this->initLinearSystem( mapper.freeSize(), basis, pde.isSymmetric() );

    gsMatrix<T> rhsVals;        // values of the right-hand side

    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell
    gsVector<T> localRhs;       // local load vector
    gsVector<T> der2;

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( this->geometry().evaluator(NEED_VALUE | NEED_JACOBIAN | NEED_2ND_DER) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->computeActiveDofs(mapper, patchIndex).rows();
        domIt->evaluateBasis( 2 );

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(domIt->quNodes);

        // evaluate right-hand side at the geometry points
        if (pde.rhs())
            pde.rhs()->eval_into( geoEval->values(), rhsVals );

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);
        localRhs.setZero(numActive);

        const gsMatrix<T>& geoDer1 = geoEval->jacobians();
        const gsMatrix<T>& geoDer2 = geoEval->derivs2();
        const typename gsMatrix<T>::Block basDer1 = domIt->basisDerivs(1);
        const typename gsMatrix<T>::Block basDer2 = domIt->basisDerivs(2);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            const T weight = domIt->quWeights[k] * geoDer1(0,k);

            // transform parametric 2nd basis derivatives into physical ones
            der2 = (basDer2.col(k) - (geoDer2(0,k)/geoDer1(0,k)) * basDer1.col(k)) / (geoDer1(0,k)*geoDer1(0,k));

            if (pde.rhs())
                localRhs += (weight * rhsVals(0,k)) * domIt->basisValues().col(k);

            localStiffness.noalias() += weight * (der2 * der2.transpose());
        }  // end loop Gauss nodes

        // add contributions from local stiffness matrix to global stiffness matrix and load vector
        this->localToGlobal_withBC(localStiffness, localRhs, ddof, mapper, domIt->activeDofs, *sys.matrix(), *sys.rhs(), true);
    } // loop over all domain elements

    sys.matrix()->makeCompressed();
    return sys;
}

template<class T>
gsSparseSystem<T>
gsGaussAssembler<T>::assembleStokes( const gsBasis<T>& basis_u, 
                                     //const gsBasis<T>& basis_p,
                                     const gsDofMapper& mapper, // mapper for u
                                     const gsMatrix<T> & ddof, const gsStokesPde<T> & pde, 
                                     int patchIndex)
{
    gsBasis<T> * basis_p = basis_u.clone();
    basis_p->degreeReduce(1);

    gsDebugVar( basis_u);
    gsDebugVar(*basis_p );
    
    gsVector<int> numNodes_u  = getNumIntNodesFor( basis_u );
    //gsVector<int> numNodes_p  = getNumIntNodesFor(*basis_p );
    //gsVector<int> numNodes_pu = getNumIntNodesForCoupled( *basis_p, basis_u );
    
    int sz_u = mapper.freeSize();
    int sz_p = basis_p->size();
    int sz = 2 * sz_u + sz_p;
    gsSparseMatrix<T> * glM = new gsSparseMatrix<T>(sz,sz);
    gsVector<T> * glB = new gsVector<T>(sz);
    gsSparseMatrix<T> & globalM = *glM;
    gsVector<T>       & globalB = *glB;

    unsigned nzRowsPerCol = 1;
    for (int i = 0; i < basis_u.dim(); ++i)
        nzRowsPerCol *= 4 * basis_u.component(i).degree() + 1;
    glM->reserve( gsVector<int>::Constant(sz, nzRowsPerCol) );

    globalM.setZero();
    globalB.setZero();

    // Get a block view of the matrix
    gsVector<index_t> nblocks(3);
    nblocks << sz_u, sz_u, sz_p ;
    typename gsSparseMatrix<T>::BlockView Mblock = globalM.blockView(nblocks, nblocks);
    typename gsSparseMatrix<T>::Block & A1  = Mblock(0,0); //
    typename gsSparseMatrix<T>::Block & A2  = Mblock(1,1); //
    typename gsSparseMatrix<T>::Block & B1  = Mblock(2,0); //
    typename gsSparseMatrix<T>::Block & B1t = Mblock(0,2); //  
    typename gsSparseMatrix<T>::Block & B2  = Mblock(2,1); //
    typename gsSparseMatrix<T>::Block & B2t = Mblock(1,2); //  
    //gsDebug<< Mblock <<"\n";

    typename gsVector<T>::BlockView Bblock = globalB.blockView(nblocks);
    typename gsVector<T>::Block & rhs_u1  = Bblock(0); //
    typename gsVector<T>::Block & rhs_u2  = Bblock(1); //
    typename gsVector<T>::Block & rhs_p   = Bblock(2); //
    //gsDebug<< Bblock <<"\n";

    const T nu = pde.viscocity();
    
    gsMatrix<T> rhsVals;        // values of the right-hand side -- = 0
    gsMatrix<T> trf_grads_k;    // transformed (physical) gradients of all basis functions at one quadrature node
    
    gsMatrix<T> localB1, localB2, localA;
    gsVector<T> localRhs_u, localRhs_p;

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( 
        this->geometry().evaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM) );
    
    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt_u = basis_u.makeDomainIterator();
    typename gsBasis<T>::domainIter domIt_p = basis_p->makeDomainIterator();

    // Set the number of integration points for each element
    domIt_u->computeQuadratureRule( numNodes_u );
    domIt_p->computeQuadratureRule( numNodes_u );

    // Start iteration over elements
    for (; domIt_u->good(); domIt_u->next())
    {
        const index_t numAct_u  = domIt_u->computeActiveDofs(mapper, patchIndex).rows();
        const index_t numAct_p  = domIt_p->computeActiveFunctions().rows();
        domIt_u->evaluateBasis( 1 );
        domIt_p->evaluateBasis( 0 );
        const typename gsMatrix<T>::Block evp  = domIt_p->basisValues();
        
        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(domIt_u->quNodes);
        
        // evaluate right-hand side at the geometry points
        if (pde.rhs())
            pde.rhs()->eval_into( geoEval->values(), rhsVals );
        
        // initialize local linear system to 0
        localA.setZero (numAct_u, numAct_u);
        localB1.setZero(numAct_p, numAct_u);
        localB2.setZero(numAct_p, numAct_u);
        localRhs_u.setZero(numAct_u);
        localRhs_p.setZero(numAct_p);
        
        //for (index_t k = 0; k < domIt_u->numQuNodes(); ++k)
        //{        }
        
        // loop over quadrature nodes
        for (index_t k = 0; k < domIt_u->numQuNodes(); ++k)      
        {
            // weight * abs(det J), where J is geometry Jacobian
            const T weight = domIt_u->quWeights[k] * geoEval->measure(k);
            
            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval->transformGradients(k, domIt_u->basisDerivs(1), trf_grads_k);
            
            // Right-hand side function
            localRhs_u += (weight * rhsVals(0,k)) * domIt_u->basisValues().col(k);
            localRhs_p += (weight * rhsVals(0,k)) * domIt_p->basisValues().col(k);
            
            localA.noalias()  += weight * nu * (trf_grads_k.transpose() * trf_grads_k);
            localB1.noalias() += weight * ( evp.col(k) * trf_grads_k.row(0) );
            localB2.noalias() += weight * ( evp.col(k) * trf_grads_k.row(1) );
        }  // end loop Gauss nodes
        
        // add contributions from local to global
        for (index_t i=0; i < numAct_u; ++i)
        {
            const int ii = domIt_u->activeDofs[i];
            if ( mapper.is_free_index(ii) )
            {
                //globalB[ii] += localRhs_u[i];
                
                for (index_t k=0; k < numAct_p; ++k)
                {
                    const int kk = domIt_p->activeFuncs(k,0);
                    B1 .coeffRef(kk, ii ) += localB1(k, i);
                    B1t.coeffRef(ii, kk ) -= localB1(k, i);
                    B2 .coeffRef(kk, ii ) += localB2(k, i);
                    B2t.coeffRef(ii, kk ) -= localB2(k, i);
                }
                
                for (index_t j=0; j < numAct_u; ++j)
                {
                    const int jj = domIt_u->activeDofs[j];
                    if ( mapper.is_free_index(jj) )
                    {
                        A1.coeffRef(ii, jj) += localA(i, j);
                        A2.coeffRef(ii, jj) += localA(i, j);                     
                    }
                    else if ( mapper.is_boundary_index(jj) )
                    {
                        const int bb = mapper.global_to_bindex(jj);
                        rhs_u1[ii] -= ddof(bb, 0) * localA(i, j);
                        rhs_u2[ii] -= ddof(bb, 1) * localA(i, j);
                    }
                }
            }
            else if ( mapper.is_boundary_index(ii) )
            {
                const int bb = mapper.global_to_bindex(ii);
                for (index_t k=0; k < numAct_p; ++k)
                {
                    const int kk = domIt_p->activeFuncs(k,0);
                    rhs_p[kk] -= ddof(bb, 0 ) * localB1(k, i)
                               + ddof(bb, 1 ) * localB2(k, i);
                }
            }
        }
        
        domIt_p->next();
    } // loop over all domain elements
    
    delete basis_p;
    
    // There are two unknowns in this problem (velocity "u" and
    // pressure "p"). Currently the velocity solution is returned, and the
    // pressure component is ignored, since we can return only one
    // unknown from the solver for now. The following code is the
    // recontruction of the pressure component:
    // 
    // gsBasis<T> * basis_p = m_bases[p]->clone();
    // basis_p->degreeReduce(1);
    // int szp = basis_p->size();
    // gsMatrix<T> coeffsp = data.bottomRows(szp);
    // return basis_p->makeGeometry( give(coeffsp) );

    glM->makeCompressed();    
    return gsSparseSystem<T>(glM, glB);
}

template<class T>
gsSparseSystem<T>
gsGaussAssembler<T>::assembleSurfacePoisson( const gsBasis<T>& basis, const gsDofMapper& mapper, const gsMatrix<T> & ddof, const gsSurfacePoissonPde<T> & pde, int patchIndex)
{
    gsDebug<< "Running for surface.\n";

    const int d = 2;

    gsVector<int> numNodes = getNumIntNodesFor( basis );
    gsSparseSystem<T> sys = this->initLinearSystem( mapper.freeSize(), basis, pde.isSymmetric() );
    
    gsMatrix<T> rhsVals;        // values of the right-hand side
    gsMatrix<T> trf_grads_k;    // transformed (physical) gradients of all basis functions at one quadrature node
    
    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell
    gsVector<T> localRhs;       // local load vector

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( 
        this->geometry().evaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->computeActiveDofs(mapper, patchIndex).rows();
        domIt->evaluateBasis( 1 );

        // compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(domIt->quNodes);

        // evaluate right-hand side at the geometry points
        if (pde.rhs())
            pde.rhs()->eval_into( geoEval->values(), rhsVals );

        // initialize local linear system to 0
        localStiffness.setZero(numActive, numActive);
        localRhs.setZero(numActive);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            const gsMatrix<T> & Jk = 
                geoEval->jacobians().block(0, k*d, d+1, d);
            gsMatrix<T> FirstFund = Jk.transpose()*Jk;
            const T weight = domIt->quWeights[k] * math::sqrt(FirstFund.determinant());
            FirstFund = FirstFund.inverse();

            trf_grads_k.resize( 2, numActive );
            for (index_t i = 0; i < numActive; ++i)
            {
                trf_grads_k.template block<2, 1>(0, i) =
                    domIt->basisDerivs(1).template block<2, 1>(i * 2, k);
            }
            
            if (pde.rhs())
                localRhs += (weight * rhsVals(0,k)) * domIt->basisValues().col(k);
            
            localStiffness.noalias() += weight * (trf_grads_k.transpose() *FirstFund* trf_grads_k);
        }  // end loop Gauss nodes
        
        // add contributions from local stiffness matrix to global
        // stiffness matrix and load vector
        this->localToGlobal_withBC(localStiffness, localRhs, ddof, mapper, domIt->activeDofs, *sys.matrix(), *sys.rhs(), true);
    } // loop over all domain elements

    sys.matrix()->makeCompressed();
    return sys;
}
/*
template<class T>
gsSparseSystem<T>
gsGaussAssembler<T>::assemblePde( const gsBasis<T>& basis, 
                                  const gsDofMapper& mapper,
                                  const gsMatrix<T> & ddof, 
                                  const gsPoissonPde<T> & pde,
                                  int patchIndex)
{
    gsVector<int> numNodes = getNumIntNodesFor( basis );
    gsSparseSystem<T> sys = this->initLinearSystem( mapper.freeSize(), basis, pde.isSymmetric() );
    
    gsMatrix<T> localMatrix; // (dense) stiffness matrix within one grid cell
    gsMatrix<T> localRhs;       // local load vector

    // Evaluate proper geometry data
    std::auto_ptr< gsGeometryEvaluator<T> > 
        geoEval ( this->geometry().evaluator(NEED_VALUE    |  
                                             NEED_MEASURE  |
                                             NEED_GRAD_TRANSFORM) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();
    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        // initialize local linear system
        const index_t numActive = domIt->computeActiveDofs(mapper, patchIndex).rows();
        localMatrix.setZero(numActive, numActive);
        localRhs.setZero(numActive,1);

        pde.apply( *domIt, *geoEval, localMatrix, localRhs );
         
         // Add contributions to the global system
        this->localToGlobal_withBC(localMatrix, localRhs, ddof, mapper, 
                                   domIt->activeDofs, *sys.matrix(), *sys.rhs(), true);
    }

    sys.matrix()->makeCompressed();
    return sys;
}
*/

template<class T>
gsSparseMatrix<T> * gsGaussAssembler<T>::massMatrix( const gsBasis<T>& B )
{

    const gsVector<int> numNodes = getNumIntNodesFor( B );

    const int bb= B.size();
    const int d = B.dim();

    // Allocate the resulting sparse matrix
    gsSparseMatrix<T> * Mh = new gsSparseMatrix<T>(bb,bb);
    Mh->setZero();
    unsigned nzRowsPerCol=1;
    for (int i=0;i!=d;++i)
        nzRowsPerCol *= 2*B.degree(i)+1;
    Mh->reserve( gsVector<int>::Constant(bb, (nzRowsPerCol+1)/2 ));

    gsMatrix<T> localMassMatrix;
    gsMatrix<T> BasisFcts_k;

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( this->geometry().evaluator(
                                                          NEED_MEASURE | NEED_GRAD_TRANSFORM) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = B.makeDomainIterator();
    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->computeActiveFunctions().rows();
        domIt->evaluateBasis( 0 );

        geoEval->evaluateAt(domIt->quNodes);

        localMassMatrix.setZero(numActive, numActive);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            const T weight = domIt->quWeights[k] * geoEval->measure(k); // weight * abs(det J), where J is geometry Jacobian

            BasisFcts_k = domIt->basisValues().col(k);
            localMassMatrix.noalias() += weight * ( BasisFcts_k * BasisFcts_k.transpose());
        }

        this->localToGlobal( localMassMatrix, domIt->activeFuncs, *Mh, true );
    }

    Mh->makeCompressed();
    return Mh;
}


template<class T>
gsSparseMatrix<T> * gsGaussAssembler<T>::stiffness( const gsBasis<T>& B )
{
    const gsVector<int> numNodes = getNumIntNodesFor( B );
    const int bb= B.size();
    const int d = B.dim();

    // Allocate the resulting sparse matrix
    gsSparseMatrix<T> * K = new gsSparseMatrix<T>(bb,bb);
    K->setZero();
    unsigned nzRowsPerCol=1;
    for (int i=0;i!=d;++i)
        nzRowsPerCol *= 2*B.degree(i)+1;
    K->reserve( gsVector<int>::Constant(bb, (nzRowsPerCol+1)/2 ));

    gsMatrix<T> localStiffness;
    gsMatrix<T> trf_grads_k;    // transformed (physical) gradients of all basis functions at one quadrature node

    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( this->geometry().evaluator(
                                                          NEED_MEASURE | NEED_GRAD_TRANSFORM) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = B.makeDomainIterator();
    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        const index_t numActive = domIt->computeActiveFunctions().rows();
        domIt->evaluateBasis( 1 );

        geoEval->evaluateAt(domIt->quNodes);

        localStiffness.setZero(numActive, numActive);

        for (index_t k = 0; k < domIt->numQuNodes(); ++k)      // loop over quadrature nodes
        {
            const T weight = domIt->quWeights[k] * geoEval->measure(k); // weight * abs(det J), where J is geometry Jacobian
            // compute physical gradients at k as a Dim x NumActive matrix
            geoEval->transformGradients(k, domIt->basisDerivs(1), trf_grads_k);

            localStiffness.noalias() += weight * (trf_grads_k.transpose() * trf_grads_k);
            /* // Equivalent:
            for (index_t i=0; i!= numActive; ++i)
                for (index_t j=0; j!= numActive; ++j)
                    localStiffness(i,j) += weight * trf_grads_k.col(i).adjoint() * trf_grads_k.col(j);
            */
        }

        this->localToGlobal( localStiffness, domIt->activeFuncs, *K, true );
    }

    K->makeCompressed();
    return K;
}


template<class T>
gsVector<T> * gsGaussAssembler<T>::moments( const gsBasis<T>& B, gsFunction<T> const& f )
{
    gsVector<int> numNodes = getNumIntNodesFor( B );
    unsigned bb= B.size();

    gsVector<T> * b = new gsVector<T>(bb);
    b->setZero();

    gsMatrix<T> fvalues;
    
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( 
        this->geometry().evaluator(NEED_MEASURE) );

    // Make domain element iterator
    typename gsBasis<T>::domainIter domIt = B.makeDomainIterator();
    // Set the number of integration points for each element
    domIt->computeQuadratureRule( numNodes );

    // Start iteration over elements
    for (; domIt->good(); domIt->next())
    {
        domIt->computeActiveFunctions();
        const gsMatrix<unsigned> & act = domIt->activeFuncs;
        geoEval->evaluateAt(domIt->quNodes);
        domIt->evaluateBasis( 0 ); // evaluate basis functions
        f.eval_into( geoEval->values(), fvalues ) ; // evaluate the function

        for (index_t k=0; k!= domIt->numQuNodes(); ++k)
        {
            const T weight = domIt->quWeights[k] * geoEval->measure(k);
            const gsVector<T> & ev = domIt->basisValues().col(k);

            for (index_t i=0; i!=act.rows(); ++i)
                (*b)( act(i,0) ) +=  weight * fvalues(0,k) * ev[i] ;
        }
    }

    return b;
}

template<class T>
gsVector<T> * gsGaussAssembler<T>::boundaryMoments( const gsBasis<T>& B, gsFunction<T> const& f , boundary::side const& s )
{
  const int d  = this->m_geometry->parDim() ;
  unsigned bb= B.size();
  
  gsVector<T> * v = new gsVector<T>(bb);
  v->setZero();

  gsGeometry<T> *GSide= this->m_geometry->boundary(s);

  //std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( this->geometry().evaluator() );
  gsMatrix<T> ev, geoPts, fvals, dets;
  gsMatrix<unsigned> act;

  if (d == 1)   // special case for 1D since no quadrature needed/possible
  {
      // get boundary point
      gsMatrix<T> param;
      param.setZero(1, 1);
      GSide->eval_into(param, geoPts);

      // evaluate given function at boundary point
      f.eval_into(geoPts, fvals);

      // get active basis value(s) (only one for boundary interpolatory bases)
      B.eval_into(param, ev);
      B.active_into(param, act);

      for (index_t i = 0; i < act.rows(); ++i)
          (*v)( act(i) ) +=  (ev)(i,0) * fvals(0);

      delete GSide;
      return v;
  }

  std::vector<std::vector<T> > breaks;
  gsVector<unsigned> meshSize(d-1);
    
  for (int i=0; i < d-1; ++i) 
  {
    breaks.push_back( B.component(i).domain()->breaks() ) ;
    meshSize[i] = breaks[i].size() - 1;   // for n breaks, we have n-1 elements (spans)
  }

  gsMatrix<T> ngrid;          // tensor Gauss nodes
  gsVector<T> wgrid;          // tensor Gauss weights
  
  gsVector<T> lower(d-1), upper(d-1), center(d-1);
  
  // Quadrature points for the boundary integral
  gsVector<int> numNodesAll = getNumIntNodesFor( B );
  gsVector<int> intnodes(d-1);
  int dir = direction(s);
  for (int i = 0; i < dir; ++i)
    intnodes[i] = numNodesAll[i];
  for (int i = dir+1; i<d; ++i)
    intnodes[i-1] = numNodesAll[i];
  
  // tensor index of the current grid cell
  gsVector<unsigned> curElement = gsVector<unsigned>::Zero(d-1);
  
  do
  {
    // loop: we are now in grid cell (curElement(0), ..., curElement(d-1))
    for (int i = 0; i < d-1; ++i)
    {
      lower[i] = breaks[i][curElement[i]];
      upper[i] = breaks[i][curElement[i]+1];
      center[i] = T(0.5) * (lower[i] + upper[i]);
    }
    
    tensorGaussRule<T>(ngrid, wgrid, intnodes, lower, upper); 
    
    B.eval_into(ngrid, ev);            // evaluate basis functions
    GSide->eval_into(ngrid, geoPts);    // evaluate geometry points
    f.eval_into(geoPts, fvals);         // evaluate function f
    computeSurfaceDets_into( *GSide, ngrid, dets ); // compute surface determinants
    B.active_into(center, act);        // compute active basis functions

    for (index_t k=0; k < ev.cols(); ++k)
    {
      for (index_t i = 0; i < act.rows(); ++i)
	(*v)( act(i) ) +=  wgrid[k] * dets(k) * (ev)(i,k) * fvals(k);
    }
  } 
  while (nextLexicographic(curElement, meshSize));

  delete GSide;
  
  return v;
}


template<class T> void 
gsGaussAssembler<T>::applyBoundary( const gsBasis<T>   & B,
                                    const boundary_condition<T> & bc,
                                    const gsDofMapper& mapper,
                                    gsSparseSystem<T> & system )
{    
    switch ( bc.type() )
    {
    case boundary::dirichlet:
        boundaryNitsche(B, bc.patch(), bc.side(), *bc.function(), mapper, system);
        break;
    case boundary::neumann:
        boundaryNeumann(B, bc.patch(), bc.side(), *bc.function(), mapper, system);
        break;
    default:
        gsWarn<<"Unknown boundary condition.\n";
    }
}


template<class T> void 
gsGaussAssembler<T>::boundaryNeumann( const gsBasis<T> & B,
                                      const int patch, 
                                      const boundary::side s,
                                      const gsFunction<T> & f,
                                      const gsDofMapper& mapper,
                                      gsSparseSystem<T> & system )
{  
    //gsDebug<<"Neumann boundary: side="<< s<<", patch="<<patch<<"\n";
    const int d   = this->m_geometry->parDim() ;
    gsVector<T> & rhs = *system.rhs();

    // Quadrature for boundary integral: we fix coordinates for
    // direction = dir to be the fixed coordinate on the edge/face
    gsVector<int> bd_intNodes = getNumIntNodesForSide( B, direction(s) );

    std::auto_ptr< gsGeometryEvaluator<T> > 
        geoEval ( this->geometry().evaluator(NEED_VALUE | NEED_JACOBIAN) );

    // Temporaries
    gsMatrix<T> fev;
    gsVector<T> localRhs, unormal(d);
        
    // iterate over all boundary grid cells
    for (typename gsDomainIterator<T>::uPtr domIter = B.makeDomainIterator(s); 
         domIter->good(); domIter->next())
    {
        // Compute the quadrature rule (nodes and weights)
        domIter->computeQuadratureRule(bd_intNodes);

        // Evaluate the geometry on the Gauss points
        geoEval->evaluateAt(domIter->quNodes);

        // Evaluate the basis functions
        const index_t numActive = domIter->computeActiveDofs(mapper, patch).rows();
        domIter->evaluateBasis();

        // Evaluate the Neumann data
        f.eval_into(geoEval->values(), fev);

        localRhs.setZero(numActive);

        for (index_t k=0; k!= domIter->numQuNodes(); ++k) // For all quadrature points
        {
            // Compute the outer normal vector on the side
            geoEval->outerNormal(k, s, unormal);

            // Sum up quadrature evaluations
            const T fff = domIter->quWeights[k] * fev(0,k) *unormal.norm();           
            localRhs.noalias() += fff * domIter->basisValues().col(k);
        }
        
        // Push element contribution to the global load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            // convert local dof index to global dof index
            const unsigned jj = domIter->activeDofs[j];
            if (mapper.is_free_index(jj))
                rhs[jj] += localRhs[j];
        }
    }

}

template<class T> void 
gsGaussAssembler<T>::boundaryNitsche( const gsBasis<T> & B,
                                      const int patch,
                                      const boundary::side s,
                                      const gsFunction<T> & f,
                                      const gsDofMapper& mapper,
                                      gsSparseSystem<T> & system )
{
    gsSparseMatrix<T> & lhs = *system.matrix();
    gsVector<T>       & rhs = *system.rhs();
    const int d   = this->m_geometry->parDim() ;

    const T mu = this->getMu(B);
    // const T h = std::pow( (T) minElementVol(B), 1.0 / B.dim() );
    // const T mu = (this->geometry().basis().degree(0)+B.dim())*
    //              (B.degree(0)+B.dim())* (B.degree(0)+1) * 2.0 / h ;

    gsDebug<<"Nitsche boundary: side="<< s<<", patch="<<patch<<"(mu="<<mu<<").\n";

    // Quadrature for boundary integral: we fix coordinates for
    // direction = dir to be the fixed coordinate on the edge/face
    gsVector<int> bd_intNodes = getNumIntNodesForSide( B, direction(s) );
    
    std::auto_ptr< gsGeometryEvaluator<T> > 
        geoEval ( this->geometry().evaluator(NEED_VALUE    | 
                         NEED_JACOBIAN | NEED_GRAD_TRANSFORM) );

    // Temporaries
    gsMatrix<T> fev, grads_k;
    gsVector<T> unormal(d);

    // Local matrix and load vector
    gsMatrix<T> LM;
    gsVector<T> LB;

    // iterate over all boundary grid cells
    for (typename gsDomainIterator<T>::uPtr domIter = B.makeDomainIterator(s); 
         domIter->good(); domIter->next())
    {
        // Compute the quadrature rule (nodes and weights)
        domIter->computeQuadratureRule(bd_intNodes);

        // Evaluate the geometry on the Gauss points
        geoEval->evaluateAt(domIter->quNodes);

        // Evaluate basis functions and their first derivatives
        const index_t numActive = domIter->computeActiveDofs(mapper, patch).rows();
        domIter->evaluateBasis( 1 );
        const typename gsMatrix<T>::Block ev  = domIter->basisValues();

        // Evaluate the Dirichlet data
        f.eval_into(geoEval->values(), fev);

        LM.setZero(numActive, numActive);
        LB.setZero(numActive);

        for (index_t k=0; k!= domIter->numQuNodes(); ++k) // For all quadrature points
        {
            // Compute the outer normal vector
            geoEval->outerNormal(k, s, unormal);

            // Integral transformation and quadarature weight
            const T fff = domIter->quWeights[k] * unormal.norm();

            // Compute the unit normal vector 
            unormal.normalize();

            // Transform the basis gradients
            geoEval->transformGradients(k, domIter->basisDerivs(1), grads_k);

            // Sum up quadrature point evaluations
            LB.noalias() += fff * fev(0,k) * ( grads_k.transpose() * unormal - mu * ev.col(k) );
            LM.noalias() += fff * ( ev.col(k) * unormal.transpose() * grads_k
                                +  (ev.col(k) * unormal.transpose() * grads_k).transpose()
                                -  mu * ev.col(k) * ev.col(k).transpose() );
        }

        // Push element contribution to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const unsigned jj = domIter->activeDofs[j];
            rhs[jj] -= LB[j];
            for (index_t i=0; i!=numActive; ++i)
            {
                const unsigned ii = domIter->activeDofs[i];
                if ( jj <= ii ) // assuming symmetric problem
                    lhs( ii, jj ) -= LM(i,j);
            }
        }
    }
}


template<class T> void 
gsGaussAssembler<T>::applyDG( const gsBasis<T> & B1, const gsBasis<T> & B2,
                              const gsGeometry<T> & geo2,
                              const boundaryInterface & bi,
                              const gsDofMapper& mapper,
                              gsSparseSystem<T> & system )
{
    const gsGeometry<T> & geo1 = *this->m_geometry;
    gsSparseMatrix<T> & lhs    = *system.matrix();
    const int d                = this->m_geometry->parDim() ;
    
    // const T h = std::pow( (T) minElementVol(B1), 1.0 / B1.dim() );
    // const T mu = (this->geometry().basis().degree(0)+B1.dim())*
    //              (B1.degree(0)+B1.dim())* (B1.degree(0)+1) * 2.0 / h ;
    
    int patch1           = bi[0].patch;
    int patch2           = bi[1].patch;
    boundary::side side1 = bi[0].side;
    boundary::side side2 = bi[1].side;
    
    // iterator on grid cells on the "right"
    typename gsDomainIterator<T>::uPtr domIter1 = B1.makeDomainIterator(side1);
    typename gsDomainIterator<T>::uPtr domIter2= B2.makeDomainIterator(side2);
    
    // Evaluators for the two patches
    std::auto_ptr< gsGeometryEvaluator<T> > 
        geoEval1 ( geo1.evaluator(NEED_VALUE| NEED_JACOBIAN| NEED_GRAD_TRANSFORM) );
    std::auto_ptr< gsGeometryEvaluator<T> > 
        geoEval2 ( geo2.evaluator(NEED_VALUE| NEED_JACOBIAN| NEED_GRAD_TRANSFORM) );

    int bSize1 = domIter1->numElements();
    int bSize2 = domIter2->numElements();

    const gsBasis<T> * Bs1 = &B1;
    const gsBasis<T> * Bs2 = &B2;

    if( bSize1 < bSize2 )
    {
        // switch sides
        std::swap( patch1   , patch2   );
        std::swap( side1    , side2    );
        std::swap( domIter1 , domIter2 );
        std::swap( geoEval1 , geoEval2 );
        std::swap( bSize1   , bSize2   );
        std::swap( Bs1      , Bs2      );
    }
    
    const int ratio = bSize1 / bSize2;
    gsDebugVar(bSize1);
    gsDebugVar(bSize2);
    gsDebugVar(ratio);
    
    const T mu                 = this->getMu(*Bs1);
    gsDebug<<"Apply DG on "<< bi <<"(mu="<<mu<<").\n";

    // const int dir1             = direction(side1);
    // const int dir2             = direction(side2);
    const bool par2               = parameter(side2);
    //GISMO_ASSERT( Bs1->component(!dir1).size() == B2.component(!dir2).size(), 
    //              "DG method not implemented yet for non matching interfaces");
    
       
    // Temporaries
    gsMatrix<T> grads_k_1, grads_k_2, basisData1, basisData2;
    gsVector<T> unormal(d);
    
    // active basis functions at one quadrature node
    gsMatrix<unsigned> actives1, actives2;

    gsMatrix<T> B11, B22, B12, B21, 
                E11, E22, E12, E21,
                N1 , N2;

    // Quadrature for boundary integrals
    gsVector<int> intNodes1 = getNumIntNodesForInterface( *Bs1, *Bs2, side1, side2, true );
    //gsVector<int> intNodes2 = getNumIntNodesForInterface( *Bs1, *Bs2, side1, side2, false );
    const index_t numNodes1 = intNodes1.prod();

    gsGaussRule<T> QuRule1(intNodes1); //, QuRule2(intNodes2);
    gsMatrix<T> quNodes1, quNodes2  ; // Mapped nodes
    gsVector<T> quWeights1; // Mapped weights

    int count = 0;
    // iterate over all boundary grid cells on the "left"
    for (; domIter1->good(); domIter1->next())
    {
        count++;
        // Get the element of the other side in domIter2
        //domIter1->adjacent( bi.orient, *domIter2 );

        // Compute the quadrature rule on both sides
        QuRule1.mapTo( domIter1->lowerCorner(), domIter1->upperCorner(), quNodes1, quWeights1 );

        mapGaussNodes( quNodes1, side1, side2, ( par2 ? 1.0 : 0.0 ), quNodes2 );
        //QuRule2.mapTo( domIter2->lowerCorner(), domIter2->upperCorner(), quNodes2, quWeights1 );

        // Compute the active basis functions
        Bs1->active_into(domIter1->centerPoint() , actives1);
        Bs2->active_into(domIter2->centerPoint() , actives2);
        const index_t numActive1 = actives1.rows();
        const index_t numActive2 = actives2.rows();

        // Evaluate basis functions and their first derivatives
        Bs1->evalAllDers_into( quNodes1, 1, basisData1);
        const typename gsMatrix<T>::Block ev1 = basisData1.topRows(numActive1);
        const typename gsMatrix<T>::Block grads1 = 
            basisData1.middleRows(numActive1, numActive1*d );
        Bs2->evalAllDers_into( quNodes2, 1, basisData2);
        const typename gsMatrix<T>::Block ev2 = basisData2.topRows(numActive2);
        const typename gsMatrix<T>::Block grads2 = 
            basisData2.middleRows(numActive2, numActive2*d );

        // Local Dofs to global dofs
        mapper.localToGlobal(actives1, patch1, actives1);
        mapper.localToGlobal(actives2, patch2, actives2);

        // Push forward the quad-points to the physical domain
        geoEval1->evaluateAt(quNodes1);
        geoEval2->evaluateAt(quNodes2);

        B11.setZero(numActive1, numActive1); B22.setZero(numActive2, numActive2); 
        B12.setZero(numActive1, numActive2); B21.setZero(numActive2, numActive1);
        E11.setZero(numActive1, numActive1); E22.setZero(numActive2, numActive2); 
        E12.setZero(numActive1, numActive2); E21.setZero(numActive2, numActive1);

        // assuming quNodes1.cols() == quNodes2.cols()
        for (index_t k=0; k!= numNodes1; ++k)
        {            
            // Compute the outer normal vector from patch1
            geoEval1->outerNormal(k, side1, unormal);

            // Integral transformation and quadarature weight (patch1)
            // assumed the same on both sides
            const T fff = quWeights1[k] * unormal.norm();
            
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
            B12.noalias() += c1 * ( val1 * N2 );

            B22.noalias() -= c1 * ( val2 * N2 );
            B21.noalias() -= c1 * ( val2 * N1 );

            const T c2 = fff * mu;

            E11.noalias() += c2 * ( val1 * val1.transpose() );
            E12.noalias() += c2 * ( val1 * val2.transpose() );

            E22.noalias() += c2 * ( val2 * val2.transpose() );
            E21.noalias() += c2 * ( val2 * val1.transpose() );
        }

        // Push element contributions 1-2 to the global matrix and load vector
        for (index_t j=0; j!=numActive1; ++j)
        {
            const index_t jj1 = actives1(j); // N1_j
            for (index_t i=0; i!=numActive1; ++i)
            {
                const index_t  ii1 = actives1(i); // N1_i

                if ( jj1 <= ii1 )
                    lhs( ii1, jj1 ) -=  B11(i,j) + B11(j,i) - E11(i,j);
            }

            for (index_t i=0; i!=numActive2; ++i)
            {
                const index_t  ii2 = actives2(i); // N2_i

                if ( jj1 <= ii2 ) 
                    lhs( ii2, jj1)  -=  B21(i,j) + B12(j,i) + E21(i,j);
            }
        }

        // Push element contributions 2-1 to the global matrix and load vector
        for (index_t j=0; j!=numActive2; ++j)
        {
            const index_t jj2 = actives2(j); // N2_j

            for (index_t i=0; i!=numActive2; ++i)
            {
                const index_t  ii2 = actives2(i); // N2_i
                
                if ( jj2 <= ii2 ) 
                    lhs( ii2, jj2 ) -=  B22(i,j) + B22(j,i) - E22(i,j);
            }

            for (index_t i=0; i!=numActive1; ++i)
            {
                const index_t  ii1 = actives1(i); // N1_i

                if ( jj2 <= ii1 ) 
                    lhs( ii1, jj2)  -=  B12(i,j) + B21(j,i) + E12(i,j);
            }

        }        
        if ( count % ratio == 0 ) 
            domIter2->next();
    }
}


template<class T>
gsMatrix<T> * gsGaussAssembler<T>::innerProduct( const gsBasis<T>& B1, const gsBasis<T>& B2)
{
  unsigned b1= B1.size();
  unsigned b2= B2.size();
  
  gsMatrix<T> * K = new gsMatrix<T>(b1,b2) ;
  K->setZero();

  std::vector<T> breaks = B1.domain()->breaks() ;
  gsMatrix<T> ngrid;          // tensor Gauss nodes
  gsVector<T> wgrid;          // tensor Gauss weights
  // Quadrature points
  int nGauss = int( ceil( double(B1.degree() + B2.degree() + 1)/2 ) );
  if (nGauss<1) nGauss=1;
  iteratedGaussRule(ngrid, wgrid, nGauss, breaks ) ;
    
  gsMatrix<T>   ev1  = B1.eval(ngrid); // Evaluate over the grid
  gsMatrix<unsigned> act1 = B1.active(ngrid);
  gsMatrix<T>   ev2  = B2.eval(ngrid); // Evaluate over the grid
  gsMatrix<unsigned> act2 = B2.active(ngrid);
  
  for (index_t k=0; k!= ngrid.cols(); ++k)
  {
    for (index_t i=0; i!=act1.rows(); ++i)
      for (index_t j=0; j!=act2.rows(); ++j)
        (*K)( act1(i,k) , act2(j,k) ) +=  wgrid(k) * ev1(i,k) * ev2(j,k) ;
  }

  return K;
}

template<class T>
gsMatrix<T> * gsGaussAssembler<T>::innerProduct1( const gsBasis<T>& B1, const gsBasis<T>& B2)
{
  unsigned b1= B1.size();
  unsigned b2= B2.size();
  
  gsMatrix<T> * K = new gsMatrix<T>(b1,b2) ;
  K->setZero();

  std::vector<T> breaks = B1.domain()->breaks() ;
  gsMatrix<T> ngrid;          // tensor Gauss nodes
  gsVector<T> wgrid;          // tensor Gauss weights
  // Quadrature points
  int nGauss = int( ceil( double(B1.degree()-1 + B2.degree()-1 + 1)/2 ) );
  if (nGauss<1) nGauss=1;
  iteratedGaussRule(ngrid, wgrid, nGauss, breaks ) ;
    
  gsMatrix<T>   ev1  = B1.deriv(ngrid); // Evaluate over the grid
  gsMatrix<unsigned> act1 = B1.active(ngrid);
  gsMatrix<T>   ev2  = B2.deriv(ngrid); // Evaluate over the grid
  gsMatrix<unsigned> act2 = B2.active(ngrid);
  
  for (index_t k=0; k!= ngrid.cols(); ++k)
  {
    for (index_t i=0; i!=act1.rows(); ++i)
      for (index_t j=0; j!=act2.rows(); ++j)
        (*K)( act1(i,k) , act2(j,k) ) +=  wgrid(k) * ev1(i,k) * ev2(j,k) ;
  }

  return K;
}

template<class T>
gsMatrix<T> * gsGaussAssembler<T>::innerProduct2( const gsBasis<T>& B1, const gsBasis<T>& B2)
{
  unsigned b1= B1.size();
  unsigned b2= B2.size();
  
  gsMatrix<T> * K = new gsMatrix<T>(b1,b2) ;
  K->setZero();

  std::vector<T> breaks = B1.domain()->breaks() ;
  gsMatrix<T> ngrid;          // tensor Gauss nodes
  gsVector<T> wgrid;          // tensor Gauss weights
  // Quadrature points
  int nGauss = int( ceil( double(B1.degree()-2 + B2.degree()-2 + 1)/2 ) );
  if (nGauss<1) nGauss=1;
  iteratedGaussRule(ngrid, wgrid, nGauss, breaks ) ;
    
  gsMatrix<T>   ev1  = B1.deriv2(ngrid); // Evaluate over the grid
  gsMatrix<unsigned> act1 = B1.active(ngrid);
  gsMatrix<T>   ev2  = B2.deriv2(ngrid); // Evaluate over the grid
  gsMatrix<unsigned>  act2 = B2.active(ngrid);
  
  for (index_t k=0; k!= ngrid.cols(); ++k)
  {
    for (index_t i=0; i!=act1.rows(); ++i)
      for (index_t j=0; j!=act2.rows(); ++j)
        (*K)( act1(i,k) , act2(j,k) ) +=  wgrid(k) * ev1(i,k) * ev2(j,k) ;
  }

  return K;
}
 

template<class T>
void gsGaussAssembler<T>::computeSurfaceDets_into( const gsGeometry<T>& surf, 
                                                   const gsMatrix<T>& x, 
                                                   gsMatrix<T>& result )
{
    const int d  = surf.geoDim();
    GISMO_ASSERT( d == surf.parDim() + 1, "Expected surface geometry");

    gsMatrix<T> jacs;
    surf.deriv_into( x, jacs );
    result.resize( 1, x.cols() );

    switch (d)
    {
    case 2:
        for (int i = 0; i < x.cols(); ++i)
            result(i) = jacs.col(i).norm();
        break;
    case 3:
        for (int i = 0; i < x.cols(); ++i)
        {
            result(i) =
                jacs.col(2*i).template segment<3>(0).cross(
                    jacs.col(2*i+1).template segment<3>(0)).norm();
        }
        break;
    default:
        GISMO_ERROR("Invalid geometry dimension in surface integral");
    }
}


template<class T>
T gsGaussAssembler<T>::minElementVol(const gsBasis<T>& b)
{
    T min(-1.0);
    T sum;
    
    std::auto_ptr< gsGeometryEvaluator<T> > 
        geoEval ( this->geometry().evaluator(NEED_VALUE    | 
                                             NEED_MEASURE) );

    gsVector<int> numNodes = gsGaussAssembler<T>::getNumIntNodesFor( this->geometry().basis() );
    
    sum = T(0.0);
    typename gsBasis<T>::domainIter domIt = b.makeDomainIterator();
    domIt->computeQuadratureRule( numNodes );
    for (; domIt->good(); domIt->next())
    {
        geoEval->evaluateAt(domIt->quNodes);
        
        for (index_t k = 0; k < domIt->numQuNodes(); ++k)
            sum += domIt->quWeights[k] * fabs( geoEval->measure(k) );

    if (sum > min)
        min = sum;
    }
    
    //gsDebugVar(min);
    return min;
}



}; // namespace gismo
