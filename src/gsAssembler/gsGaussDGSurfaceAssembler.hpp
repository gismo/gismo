#pragma once

#include <gsUtils/gsQuadRules.hpp>
#include <gsUtils/gsQuadrature.hpp>
#include <gsCore/gsMemory.hpp>
#include <gsUtils/gsCombinat.hpp>
#include <gsCore/gsDofMapper.h>
#include <gsCore/gsDebug.h>

namespace gismo {
  
//------------------------------------------------------------------------------------//
  
template<class T >
// class gsGaussDGSurfaceAssembler;
gsSparseSystem<T>
gsGaussDGSurfaceAssembler<T>::assembleDGSurfacePoisson( const gsBasis<T>& m_basis, 
					     const gsDofMapper& mapper, 
					     const gsMatrix<T> & ddof, 
					     const gsPoissonPde<T> & pde, 
					     int penalty, 
					     int patchIndex)
// gsGaussDGSurfaceAssembler<T>::assemblePoisson( const gsBasis<T>& m_basis, const gsDofMapper& mapper, const gsMatrix<T> & ddof, const gsFunction<T> & rhs, int patchIndex)
// TODO: plus.. Neumann condtions
  {
  // Eigen provides special optimizations for statically
  // sized small matrices, so take advantage of that by
  // dispatching statically on dimensions 1-4
    switch (this->geometry().parDim())
    {
       case 2:  return assembleDGSurfacePoisson_impl<2> (m_basis, mapper, ddof, pde, patchIndex);
       default: return assembleDGSurfacePoisson_impl<Dynamic> (m_basis, mapper, ddof, pde, patchIndex);
     }
   }
   
//-----------------------------------------------------------------------------------//

template<class T>
template <int D>
gsSparseSystem<T>
gsGaussDGSurfaceAssembler<T>::assembleDGSurfacePoisson_impl( const gsBasis<T>& m_basis, 
							     const gsDofMapper& mapper, 
							     const gsMatrix<T> & ddof, 
							     const gsPoissonPde<T> & pde, 
							     int penalty, 
							     int patchIndex)
  // TODO: plus.. Neumann condtions
  {
    const int d  = this->geometry().parDim() ;
    assert( D == Dynamic || D == d );

    // TODO :  grab m_basis from mapper instead of argument
    //gsBasis<T>& m_basis = mapper.basis();

    std::vector<std::vector<T> > breaks;
    gsVector<unsigned, D> meshSize;
    meshSize.setZero(d);
    for (int i=0; i!=d; ++i) 
    {
        breaks.push_back( m_basis.component(i).domain()->breaks() ) ;
        meshSize[i] = breaks[i].size() - 1;   // for n breaks, we have n-1 elements (spans)
    }
    
    const int nDofs= mapper.freeSize();
    
    // initialize stiffness matrix and reserve space for entries
    gsSparseMatrix<T> * K = new gsSparseMatrix<T>(nDofs,nDofs) ;
    unsigned nzRowsPerCol=1;
    for (int i=0;i!=d;++i)
	nzRowsPerCol *= 2 * m_basis.component(i).degree() + 1;    
    K->reserve( gsVector<T>::Constant(nDofs,nzRowsPerCol));

    // initialize right-hand side vector to zero
    gsVector<T> * b = new gsVector<T>(nDofs) ;
    b->setZero();
    
    gsMatrix<T> ngrid;          // tensor Gauss nodes
    gsVector<T> wgrid;          // tensor Gauss weights

    gsMatrix<unsigned> act;     // active basis functions
    gsVector<int> activeIndex;  // dof indices of the active basis functions
    gsMatrix<T> values;         // values and derivatives of the active basis functions
    gsMatrix<T> rhsVals;        // values of the right-hand side
    gsMatrix<T> diffVals;       // Diffussion coefficient values

    // temporaries
    gsVector<T>     gradients_k;
    gsVector<T,D>   grad_i(d);
    gsMatrix<T,D,D> FirstFund(d,d), invFirstFund(d,d), ggT(d,d);
    gsMatrix<T,D+1,D> jac(d+1,d);

    gsMatrix<T> localStiffness; // (dense) stiffness matrix within one grid cell
//     gsVector<T> localRhs ; // rhs within one grid cell

    // bounding box of the current grid cell
    // initialization needed for silencing Eigen warning
    gsVector<T,D> lower  = gsVector<T, D>::Zero(d);
    gsVector<T,D> upper  = gsVector<T, D>::Zero(d);
    gsVector<T,D> center = gsVector<T, D>::Zero(d);

    // tensor index of the current grid cell
    gsVector<unsigned, D> curElement = gsVector<unsigned, D>::Zero(d);

    // Incorrect results from the evaluator, going back to the standard functions.
    std::auto_ptr< gsGeometryEvaluator<T> > geoEval ( this->geometry().evaluator() );

    do
    {
        // loop: we are now in grid cell (curElement(0), ..., curElement(d-1))

        for (int i = 0; i < d; ++i)
        {
            lower[i] = breaks[i][curElement[i]];
            upper[i] = breaks[i][curElement[i]+1];
            center[i] = T(0.5) * (lower[i] + upper[i]);
        }

        tensorGaussRule<T>(ngrid, wgrid, m_numIntNodes, lower, upper);

        // ngrid is a dxN matrix of Gauss nodes, where N is the number of quadrature points
        // wgrid is a vector of Gauss weights

        m_basis.active_into(center, act);               // basis functions which are active in this element
        const index_t numActive = act.rows();
	

	// precompute global dof indices of the active functions in the current element (optimization)
        activeIndex.resize(numActive);
        for (int i = 0; i < numActive; ++i)
            activeIndex[i] = mapper.index(act(i,0), patchIndex);

        // Problematic?
        // create blocks which refer to the basis functions and to the gradients 
        m_basis.evalAllDers_into(ngrid, 1, values);	
        typename gsMatrix<T>::Block ev  = values.topRows(numActive);
        typename gsMatrix<T>::Block dev = values.bottomRows(values.rows() - numActive);

	// compute image of Gauss nodes under geometry mapping as well as Jacobians
        geoEval->evaluateAt(ngrid);
	
        // evaluate right-hand side at the geometry points
	if (pde.rhs())
            pde.rhs()->eval_into( geoEval->values(), rhsVals );
// 	 pde.diffusion()->eval_into( geoEval->values(), diffVals );
	
       // initialize element stiffness and load matrix to 0	
        localStiffness.setZero(numActive, numActive);

        for (index_t k = 0; k < ev.cols(); ++k)      // loop over Gauss nodes
        {
            // column vector which contains all gradients at point k of the
            // active basis functions stacked on top of each other
            gradients_k = dev.col(k);
            const T weight = wgrid[k];

            // get geometry Jacobian at point k
            jac = geoEval->jacobians().block(0, k*d, d+1, d); //Jacobian is R^{3x2} 
  
	    // noalias helps in optimization of computation
	    FirstFund.noalias() = jac.transpose()*jac;  //First fundamental form is R^{2x2} 
	    invFirstFund = FirstFund.inverse() ;
	    
 	    // functional determinant of geometry mapping
            const T funcDet = sqrt(fabs( FirstFund.determinant() )) ; 
	    
	    // transformation matrix for a grad * grad term (gradients transform by [J*F^(-1)]^(T))
	    ggT.noalias() = (weight * funcDet) * invFirstFund.transpose() ;
	    	               
            for (index_t i = 0; i < numActive; ++i)
            {
                const int ii = activeIndex[i];// index as a row of stiffness matrix without boundary dofs
                if ( mapper.is_free_index(ii) )
                {
		    if (pde.rhs())
                       (*b)[ii] +=  weight * rhsVals(0,k) * ev(i,k) * funcDet;
		    
                    // weighted and transformed gradient of i-th active basis function
                    // grad_i.noalias() = weight * funcDet * invFirstFund * gradients_k.segment(i*d,d);
		    grad_i.noalias() = ggT * gradients_k.segment(i*d,d);

                    for (index_t j=0; j < numActive; ++j)
                    {
//                         // dN(i,p)dN(j,p)|Det|
                                localStiffness(i, j) +=  grad_i.dot(gradients_k.segment(j*d,d));
                      } // end loop j
                  }
                }  // end loop i
            }  // end loop Gauss nodes

        // add contributions from local stiffness matrix to global stiffness matrix
        for (index_t i=0; i < numActive; ++i)
        {
            const int ii = activeIndex[i];
            if ( mapper.is_free_index(ii) )
            { 
                for (index_t j=0; j < numActive; ++j)
                {
                    const int jj = activeIndex[j];
                    if ( mapper.is_free_index(jj) )
		    {
		      // if matrix is symmetric, store only lower triangular part
		      if (jj <= ii)
			K->coeffRef(ii, jj) += localStiffness(i, j);
		    }
                    else if ( mapper.is_boundary_index(jj) ) // subtract Dirichlet boundary
		      {
			(*b)[ii] -= ddof( mapper.global_to_bindex(jj) ) *localStiffness(i,j);
		      }
		}
	    }
	}
    }
    while (nextLexicographic(curElement, meshSize));
    
       K->makeCompressed();
      
    return gsSparseSystem<T>(K,b);
}
  
// //---------------------------------------------------------------------------------------------//
// 
// template<class T>
// gsMatrix<T> * gsGaussDGSurfaceAssembler<T>::dgPoisson ( const gsBasis<T>* B, gsMatrix<T> * const& ngrid, gsMatrix<T> * const& wgrid )
// {
//  
//   unsigned bb= B->size();
//   gsMatrix<T> * K = new gsMatrix<T>(gsMatrix<T>::Zero(bb,bb) ) ;
// 
//   // TO DO : maybe compute inside the loops ??
//   gsMatrix<T> * ev  = B->eval(*ngrid); // Evaluate over the grid
//   gsMatrix<unsigned> * act = B->active(*ngrid);
//   T tmp;
// 
//   for (index_t k=0; k!= ev->cols(); ++k) // for all quad points
//   {
//     // evaluate the determinant of the Jacobian at point k
//     tmp= fabs( this->m_geometry->jac( ngrid->col(k) )->determinant() ) ;
// 
//     for (index_t i=0; i!=act->rows(); ++i)
//       for (index_t j=0; j!=act->rows(); ++j)
//         (*K)( (*act)(i,k) , (*act)(j,k) ) +=  (*wgrid)(k) * (*ev)(i,k) * (*ev)(j,k) * tmp ;
//   }
// 
//   delete ev;
//   delete act;
// 
//   return K;
// }  
//   
// };


}// namespace gismo
