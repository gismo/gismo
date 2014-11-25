
#pragma once

#include <iostream>
#include <gsAssembler/gsAssembler.h>
#include <gsPde/gsSurfacePoissonPde.h>

namespace gismo
{    

/** @brief
    Implementation of a surface assembler using Gauss quadrature.
*/
    
template<class T>
class gsGaussSurfaceAssembler : public gsAssembler<T>
{
public:
    /// Default empty constructor
    gsGaussSurfaceAssembler() : gsAssembler<T>() { }
    
    /// Construct by a provided geometry
    gsGaussSurfaceAssembler(const gsGeometry<T> & geom) : gsAssembler<T>(geom)
    {
        this->setGeometry( geom );
    }
    
    ~gsGaussSurfaceAssembler()                 //destructor
    { }

    virtual void setGeometry(const gsGeometry<T> & geom)
    { 
      this->m_geometry = &geom;
      // Caution: integration points should be set wrt the degree of
      // the discretization basis, not the geometry map
    }

public:

    static gsVector<int> getNumIntNodesFor(const gsBasis<T>& b)
    {
        gsVector<int> numNodes( b.dim() );
        for (int i = 0; i < b.dim(); ++i)
            numNodes[i] = b.degree(i) + 1;
        return numNodes;
    }

    virtual gsSparseSystem<T>
    assemble( const gsBasis<T>& basis, const gsDofMapper& mapper, 
              const gsMatrix<T> & ddof, const gsPde<T> & pde, int patchIndex=0)
    {
        GISMO_NO_IMPLEMENTATION
    }

    virtual gsSparseSystem<T>
    assemble( const gsBasis<T>& basis,
              const gsMatrix<T> & ddof, const gsPde<T> & pde, int patchIndex=0)
    {
        const gsSurfacePoissonPde<T> *surfacepoisson = 
            dynamic_cast<const gsSurfacePoissonPde<T>*>(&pde);
        if (surfacepoisson)
            return assembleSurfacePoisson(basis, mapper,ddof, *surfacepoisson, patchIndex);

        GISMO_ERROR("Unknown PDE type in assemble()");
    }

    /// Assembler for single patch Surface Poisson equation
    gsSparseSystem<T> assembleSurfacePoisson( const gsBasis<T>& basis, 
                                              const gsDofMapper& mapper,
                                              const gsMatrix<T> & ddof, 
                                              const gsSurfacePoissonPde<T> & pde, 
                                       int patchIndex=0);

    void applyBoundary( const gsBasis<T>   & B,
                        const boundary_condition<T> & bc,
                        const gsDofMapper& mapper,
                        gsSparseSystem<T> & system );
    
    void applyDG( const gsBasis<T> & B1, const gsBasis<T> & B2,
                  const gsGeometry<T> & geo2,
                  const boundaryInterface & bi,
                  const gsDofMapper& mapper,
                  gsSparseSystem<T> & system );

private:

    /// Add contribution of Nitsche Dirichlet boundary to matrix K
    /// \param B is a boundary basis, \param f is the Dirichlet function
    void boundaryNitsche( const gsBasis<T>   & B,
                          const int patch, 
                          const boundary::side s,
                          const gsFunction<T> & f,
                          gsSparseSystem<T> & system );

    /// Add contribution of Neumann boundary condition to he \a system
    /// \param B is a boundary basis, \param f is the Dirichlet function    
    void boundaryNeumann( const gsBasis<T>   & B,
                          const int patch, 
                          const boundary::side s,
                          const gsFunction<T> & f,

                          gsSparseSystem<T> & system );

    static gsVector<int> getNumIntNodesForSide(const gsBasis<T>& b, int dir)
    {
        gsVector<int> numNodes ( b.dim() );
        for (int i = 0; i < dir; ++i)
            numNodes[i] = b.degree(i) + 1;
        numNodes[dir] = 1;
        for (int i = dir+1; i < b.dim(); ++i)
            numNodes[i] = b.degree(i) + 1;
        return numNodes;
    }

    static gsVector<int> getNumIntNodesForInterface(const gsBasis<T>& b1, const gsBasis<T>& b2,
                           const boundaryInterface & bi, bool left = true)
    {
        // assumes matching orientation
        gsVector<int> numNodes ( b1.dim() );
        const int dir = ( left ? direction( bi.ps1.side ) : direction( bi.ps2.side ) );
        for (int i = 0; i < dir; ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        numNodes[dir] = 1;
        for (int i = dir+1; i < b1.dim(); ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        return numNodes;
    }

    static gsVector<int> getNumIntNodesForCoupled(const gsBasis<T>& b1, const gsBasis<T>& b2)
    {
        gsVector<int> numNodes ( b1.dim() );
        for (int i = 0; i < b1.dim(); ++i)
            numNodes[i] = ( b1.degree(i) + b2.degree(i) + 2 )/ 2 ;
        return numNodes;
    }

public:

    gsDofMapper mapper;

}; // class gsGaussSurfaceAssembler

//////////////////////////////////////////////////
//////////////////////////////////////////////////


template<class T>
gsSparseSystem<T>
gsGaussSurfaceAssembler<T>::assembleSurfacePoisson( const gsBasis<T>& basis, 
                                                    const gsDofMapper& mapper, 
                                                    const gsMatrix<T> & ddof, 
                                                    const gsSurfacePoissonPde<T> & pde, 
                                                    int patchIndex)
{
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
            
            localStiffness.noalias() += weight * (trf_grads_k.transpose() *FirstFund* trf_grads_k
                                       //Plus a mass term
                                       + domIt->basisValues().col(k) * domIt->basisValues().col(k).transpose()
                                     );

;


        }  // end loop Gauss nodes
        
        // add contributions from local stiffness matrix to global
        // stiffness matrix and load vector
        this->localToGlobal_withBC(localStiffness, localRhs, ddof, mapper, 
                                   domIt->activeDofs, *sys.matrix(), *sys.rhs(), true);
    } // end loop over all domain elements

    sys.matrix()->makeCompressed();
    return sys;
}


template<class T> void 
gsGaussSurfaceAssembler<T>::applyDG( const gsBasis<T> & B1, const gsBasis<T> & B2,
                              const gsGeometry<T> & geo2,
                              const boundaryInterface & bi,
                              const gsDofMapper& mapper,
                              gsSparseSystem<T> & system )
{
    const gsGeometry<T> & geo1 = *this->m_geometry;
    gsSparseMatrix<T> & lhs    = *system.matrix();
    const int d                = this->m_geometry->parDim() ;
    const T mu                 = this->getMu(B1);
    
    //gsDebug<<"Apply surf. DG on "<< bi <<"(mu="<<mu<<").\n";
    
    const int patch1           = bi[0].patch;
    const int patch2           = bi[1].patch;
    const boundary::side side1 = bi[0].side;
    const boundary::side side2 = bi[1].side;
    const T orient = sideOrientation(side1);
    //const T orient2 = sideOrientation(side2);
    //gsDebugVar( orient );
    //gsDebugVar( orient2 );

    const int dir1             = direction(side1);
    //const int dir2             = direction(side2);

    //GISMO_ASSERT( B1.component(!dir1).size() == B2.component(!dir2).size(), 
    //              "DG method not implemented yet for non matching interfaces");
    
    // Quadrature for boundary integrals
    gsVector<int> intNodes1 = getNumIntNodesForInterface( B1, B2, bi, true  );
    gsVector<int> intNodes2 = getNumIntNodesForInterface( B1, B2, bi, false );
    
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
        domIter2->computeActiveDofs(mapper, patch2).rows();
        domIter2->evaluateBasis( 1 );
        const typename gsMatrix<T>::Block ev2  = domIter2->basisValues();
        
        B11.setZero(numActive, numActive); B22.setZero(numActive, numActive); 
        B12.setZero(numActive, numActive); B21.setZero(numActive, numActive);
        E11.setZero(numActive, numActive); E22.setZero(numActive, numActive); 
        E12.setZero(numActive, numActive); E21.setZero(numActive, numActive);

        // assuming domIter1->quNodes.cols() == domIter2->quNodes.cols()
        for (index_t k=0; k!= domIter1->numQuNodes(); ++k)
        {            
            // Compute first fund. form
            const gsMatrix<T> & jac1 = geoEval1->jacobians().block(0, k*d, d+1, d);
            const gsMatrix<T> & jac2 = geoEval2->jacobians().block(0, k*d, d+1, d);
            const gsMatrix<T> FirstFund1 = jac1.transpose()*jac1;
            const gsMatrix<T> FirstFund2 = jac2.transpose()*jac2;
            
            // Transform the basis gradients
            grads_k_1 = domIter1->basisDerivs(1).col(k);
            grads_k_2 = domIter2->basisDerivs(1).col(k);
            grads_k_1.resize( 2, numActive );
            grads_k_2.resize( 2, numActive );
            grads_k_1 = jac1 * FirstFund1.inverse() * grads_k_1;
            grads_k_2 = jac2 * FirstFund2.inverse() * grads_k_2;// param. assumed conforming    
            // *** Compute the co-normal vector 
            const gsVector<T,3> tangent  =  orient * jac1.template block<3, 1>(0,!dir1);
            // Check for zero tangent
            if ( tangent.squaredNorm() < 1e-10 ) 
            {
                gsDebug<< "Skip "<< geoEval1->values().col(k).transpose() <<"\n";
                continue;
            }

            // parametrization normal side 1
            const gsVector<T,3> outer_normal = jac1.template block<3, 1>(0,0).cross(
                                               jac1.template block<3, 1>(0,1) ).normalized();
       
            unormal = tangent.cross( outer_normal ).normalized() ;
            //gsDebugVar( unormal.transpose() );

            // Integral transformation and quadarature weight
            const T fff = domIter1->quWeights[k] * tangent.norm(); //unormal.norm();
            //gsDebugVar( tangent2.cross( u_surf_normal ).normalized().transpose() ) ;

/*
            // parametrization normal side 1
            const gsVector<T,3> s_normal1 = jac1.template block<3, 1>(0,0).cross(
                jac1.template block<3, 1>(0,1) ).normalized();

            // parametrization normal side 2
            const gsVector<T,3> s_normal2 = jac2.template block<3, 1>(0,0).cross(
                jac2.template block<3, 1>(0,1) ).normalized();

            // Tangent vector side 2
            gsVector<T,3> tangent2 =  orient2 * jac2.template block<3, 1>(0,!dir2);

            // Check for mirrored parametrization patch 1
            if ( (s_normal1 - outer_normal ).squaredNorm() > 1.0 ) 
            {
                //tangent.array() *= -1.0;
                gsDebug<<"**toggle tangent 1\n";
            }

            gsDebugVar( geoEval1->values().col(k).transpose() );
            gsDebugVar( s_normal1.transpose() );
            gsDebugVar( s_normal2.transpose() );
            gsDebugVar( tangent.transpose()   );
            gsDebugVar( tangent2.transpose()  );
*/

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
        
        //gsDebugVar( N1 );
        
        // Push element contributions to the global matrix and load vector
        for (index_t j=0; j!=numActive; ++j)
        {
            const index_t jj1 = domIter1->activeDofs[j]; // N1_j
            const index_t jj2 = domIter2->activeDofs[j]; // N2_j
            for (index_t i=0; i!=numActive; ++i)
            {
                // assuming symmetric problem
                const index_t  ii1 = domIter1->activeDofs[i]; // N1_i
                const index_t  ii2 = domIter2->activeDofs[i]; // N2_i
                
                if ( jj1 <= ii1 )
                    lhs( ii1, jj1 ) -=  B11(i,j) + B11(j,i) - E11(i,j);
                if ( jj2 <= ii2 ) 
                    lhs( ii2, jj2 ) -=  B22(i,j) + B22(j,i) - E22(i,j);
                if ( jj2 <= ii1 ) 
                    lhs( ii1, jj2)  -=  B12(i,j) + B21(j,i) + E12(i,j);
                if ( jj1 <= ii2 ) 
                    lhs( ii2, jj1)  -=  B21(i,j) + B12(j,i) + E21(i,j);
            }
        }
        
        domIter2->next();
    }
}

template<class T> void 
gsGaussSurfaceAssembler<T>::applyBoundary( const gsBasis<T>   & B,
                                           const boundary_condition<T> & bc,
                                           const gsDofMapper& mapper,
                                           gsSparseSystem<T> & system )
{    
    switch ( bc.type() )
    {
    case boundary::dirichlet:
        boundaryNitsche(B, bc.patch(), bc.side(), *bc.function(), system);
        break;
    case boundary::neumann:
        boundaryNeumann(B, bc.patch(), bc.side(), *bc.function(), system);
        break;
    default:
        gsWarn<<"Unknown boundary condition.\n";
    }
}


template<class T> void 
gsGaussSurfaceAssembler<T>::boundaryNitsche( const gsBasis<T> & B,
                                      const int patch,
                                      const boundary::side s,
                                      const gsFunction<T> & f,
                                      gsSparseSystem<T> & system )
{
    gsSparseMatrix<T> & lhs = *system.matrix();
    gsVector<T>       & rhs = *system.rhs();
    const int d   = this->m_geometry->parDim() ;

    const T mu = this->getMu(B);
    const int dir      = direction(s);

    //gsDebug<<"Nitsche boundary: side="<< s<<", patch="<<patch<<"(mu="<<mu<<").\n";

    // Quadrature for boundary integral: we fix coordinates for
    // direction = dir to be the fixed coordinate on the edge/face
    gsVector<int> bd_intNodes = getNumIntNodesForSide( B, direction(s) );
    
    std::auto_ptr< gsGeometryEvaluator<T> > 
        geoEval ( this->geometry().evaluator(NEED_VALUE    | 
                         NEED_JACOBIAN | NEED_GRAD_TRANSFORM) );
    const T orient = sideOrientation(s);

    // Temporaries
    gsMatrix<T> fev, grads_k;
    gsVector<T> unormal(d);

    // Local matrix and load vector
    gsMatrix<T> LM;
    gsVector<T> LB;

    typename gsDomainIterator<T>::uPtr domIter = B.makeDomainIterator(s); 

    // iterate over all boundary grid cells
    for (; domIter->good(); domIter->next())
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
            //  *** Compute first fund. form
            const gsMatrix<T> & jac = geoEval->jacobians().block(0, k*d, d+1, d);
            const gsMatrix<T> FirstFund = jac.transpose()*jac;
                        
            //  *** Transform the basis gradients
            grads_k = domIter->basisDerivs(1).col(k);
            grads_k.resize( 2, numActive );
            grads_k = jac * FirstFund.inverse()*grads_k;

            // *** Compute the co-normal vector 
            const gsVector<T,3> tangent  = orient * jac.template block<3, 1>(0,!dir);
            // Check for zero tangent
            if ( tangent.squaredNorm() < 1e-10 ) 
            {
                gsDebug<< "Skip "<< geoEval->values().col(k).transpose() <<"\n";
                continue;
            }

            const gsVector<T,3> outer_normal = jac.template block<3, 1>(0,0).cross(
                                                jac.template block<3, 1>(0,1) ).normalized();
            unormal = tangent.cross( outer_normal ).normalized() ;

            // Integral transformation and quadarature weight
            const T fff = domIter->quWeights[k] * tangent.norm();
           
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
gsGaussSurfaceAssembler<T>::boundaryNeumann( const gsBasis<T> & B,
                                      const int patch, 
                                      const boundary::side s,
                                      const gsFunction<T> & f,
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
            rhs[jj] += localRhs[j];
        }
    }

}



} // namespace gismo

