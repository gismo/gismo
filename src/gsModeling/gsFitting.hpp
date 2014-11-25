/** @file gsFitting.hpp

    @brief Provides implementation of data fitting algorithms by least
    squares approximation.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Kapl, G. Kiss, A. Mantzaflaris

*/

#include <gsCore/gsGeometry.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsTensor/gsTensorDomainIterator.h>


namespace gismo
{

template<class T>
gsFitting<T>::~gsFitting() 
{
    if ( m_result )
        delete m_result;
}

template<class T>
gsFitting<T>::  gsFitting(gsMatrix<T> const & param_values, 
                          gsMatrix<T> const & points, 
                          gsBasis<T>  & basis)
{
    m_param_values = param_values;
    m_points = points;
    m_result = NULL;
    m_basis = &basis;
    m_points.transposeInPlace();
}

template<class T>
void gsFitting<T>::compute(T lambda)
{
    // Wipe out previous result
    if ( m_result )
        delete m_result;

    //number of basis functions
    const int num_basis=m_basis->size();
    // number of points
    const int num_points=m_points.rows();

    //for computing the value of the basis function
    gsMatrix<T> value;
    gsMatrix<unsigned> actives;

    const int dimension=m_points.cols();

    //left side matrix
    //gsMatrix<T> A_mat(num_basis,num_basis);
    gsSparseMatrix<T> A_mat(num_basis,num_basis);
    //gsMatrix<T>A_mat(num_basis,num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here
    A_mat.reserve( gsVector<T>::Constant(num_basis, 15 ) );

    A_mat.setZero(); // ensure that all entries are zero in the beginning
    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis, dimension);
    // enusure that all entries are zero in the beginning
    m_B.setZero(); 
    // building the matrix A and the vector b of the system of linear
    // equations A*x==b
    for(index_t k=0;k<num_points;k++)
    {
        const gsMatrix<T> & curr_point = m_param_values.col(k);
        //computing the values of the basis functions at the current point
        m_basis->eval_into(curr_point, value);
        // which functions have been computed i.e. which are active
        m_basis->active_into(curr_point, actives);
        const index_t numActive = actives.rows();

        for (index_t i=0; i!=numActive; ++i)
        {
            const int ii = actives(i,0);
            m_B.row(ii) += value(i,0)*m_points.row(k);
            for (index_t j=0; j!=numActive; ++j)
                A_mat( ii, actives(j,0) ) += value(i,0)*value(j,0);
        }
    }

    // --- Smoothing matrix computation
    //test degree >=3
    if(lambda > 0)
      applySmoothing(lambda, A_mat);

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)

    A_mat.makeCompressed();
    Eigen::BiCGSTAB< gsSparseMatrix<T>,  Eigen::IncompleteLUT<T> > solver( A_mat );

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
        m_result = NULL;
        return;
    }
    // Solves for many right hand side  columns
    gsMatrix<T> x;
    x = solver.solve(m_B); //toDense()
    //gsMatrix<T> x (m_B.rows(), m_B.cols());
    //x=A_mat.fullPivHouseholderQr().solve( m_B);
    // Solves for many right hand side  columns
    // finally generate the B-spline curve
    m_result = m_basis->makeGeometry( give(x) );
}



template<class T>
void gsFitting<T>::applySmoothing(T lambda, gsSparseMatrix<T> & A_mat)
{
    const int dim = m_basis->dim();
    const int stride = dim*(dim+1)/2;
    gsVector<int> numNodes(dim);
    for ( int i = 0; i!= dim; ++i )
        numNodes[i] = this->m_basis->degree(i);//+1;
    
    gsMatrix<T> der2, localA;
    
    typename gsBasis<T>::domainIter domIt = m_basis->makeDomainIterator();
    domIt->computeQuadratureRule( numNodes );
    for (; domIt->good(); domIt->next() )
    {
        m_basis->deriv2_into(domIt->quNodes, der2);
        const index_t numActive = domIt->computeActiveFunctions().rows();
        localA.setZero(numActive, numActive );

        // perform the quadrature
        for (index_t k = 0; k < domIt->numQuNodes(); ++k)
        {
            const T weight = domIt->quWeights[k] * lambda;

            for (index_t i=0; i!=numActive; ++i)
                for (index_t j=0; j!=numActive; ++j)
                {
                    T localAij = 0; // temporary variable

                    for (int s = 0; s < stride; s++)
                    {
                        // The pure second derivatives
                        // d^2u N_i * d^2u N_j + ...
                        if (s < dim)
                        {
                            localAij += der2(i * stride + s, k) *
                                        der2(j * stride + s, k);
                        }
                        // Mixed derivatives 2 * dudv N_i * dudv N_j + ...
                        else
                        {
                            localAij += 2 * der2(i * stride + s, k) *
                                            der2(j * stride + s, k);
                        }
                    }

                    localA(i, j) += weight * localAij;

                    // old code, just for the case if I break something

//                    localA(i,j) += weight * (
//                            // der2.template block<3,1>(i*stride,k)
//                            der2(i*stride  , k) * der2(j*stride  , k) +  // d^2u N_i * d^2u N_j
//                            der2(i*stride+1, k) * der2(j*stride+1, k) +  // d^2v N_i * d^2v N_j
//                            2 * der2(i*stride+2, k) * der2(j*stride+2, k)// dudv N_i * dudv N_j
//                            );
                }
        }
        
        for (index_t i=0; i!=numActive; ++i)
        {
            const int ii = domIt->activeFuncs(i,0);
            for (index_t j=0; j!=numActive; ++j)
                A_mat( ii, domIt->activeFuncs(j,0) ) += localA(i,j);
        }
    }
}
  
template<class T>
void gsFitting<T>::computeApproxError(T& error, int type) const
{
    gsMatrix<T> results;
    m_result->eval_into(m_param_values, results);
    results.transposeInPlace();
    error = 0;

    //computing the approximation error = sum_i ||x(u_i)-p_i||^2

    for (index_t row = 0; row != m_points.rows(); row++)
    {
        T err = 0;
        for (index_t col = 0; col != m_points.cols(); col++)
        {
            err += pow(m_points(row, col) - results(row, col), 2);
        }

        switch (type) {
        case 0:
            error += err;
            break;
        case 1:
            error += sqrt(err);
            break;
        default:
            gsWarn << "Unknown type in computeApproxError(error, type)...\n";
            break;
        }
    }
}

template<class T>
void gsFitting<T>::get_Error(std::vector<T>& errors, int type) const
{ 
    errors.clear();
    gsMatrix<T> results;
    m_result->eval_into(m_param_values,results);
    results.transposeInPlace();

    for (index_t row = 0; row != m_points.rows(); row++)
    {
        T err = 0;
        for (index_t col = 0; col != m_points.cols(); col++)
        {
            err += pow(m_points(row, col) - results(row, col), 2);
        }

        switch (type)
        {
        case 0:
            errors.push_back(err);
            break;
        case 1:
            errors.push_back(sqrt(err));
            break;
        default:
            gsWarn << "Unknown type in get_Error(errors, type)...\n";
            break;
        }
    }
}

} // namespace gismo

/* from svn version 1230-
template<class T>
void gsFitting<T>::compute(T lambda)
{
    //number of basis functions
    int num_basis=m_basis->size();

    gsDebug << "num_basis: "<< num_basis  <<"\n";
    int num_points=m_points.rows();
    //for computing the value of the basis function

    gsMatrix<T> values;
    gsMatrix<unsigned> actives;

    //computing the values of the basis functions at some position
    m_basis->eval_into(m_param_values,values);
    // which functions have been computed i.e. which are active
    m_basis->active_into(m_param_values,actives);
    // dimension of points will be also later dimension of coefficients
    int m_dimension=m_points.cols();
    //left side matrix
    gsMatrix<T> m_A(num_basis,num_basis); // TO DO: change to sparse matrix

    //gsSparseMatrix<T> m_A(num_basis,num_basis);
    // m_A.reserve(...); //to opotimiye sparse matrix an estimation of nonyero elements per column can be given here
    m_A.setZero(); // ensure that all entries are zero in the beginning
    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis,m_dimension);
    m_B.setZero(); // enusure that all entris are zero in the beginning
    for(index_t k=0;k<num_points;k++){
        for(index_t i=0;i<actives.rows();i++){
            m_B.row(actives(i,k)) += values(i,k)*m_points.row(k);
            for(index_t j=0;j<actives.rows();j++){
                m_A(actives(i,k), actives(j,k) ) += values(i,k)*values(j,k);
            }
        }
    }

    if(lambda > 0){
     //test degree >=3
            int dim = m_basis->dim();
            std::vector<std::vector<T> > breaks;
            gsVector<unsigned> meshSize(dim);
            for (int i=0; i!=dim; ++i)
              {
            breaks.push_back( m_basis->component(i).domain()->breaks() ) ;
            meshSize[i] = breaks[i].size() - 1;   // for n breaks, we have n-1 elements (spans)
              }
            gsMatrix<T> ngrid1;          // tensor Gauss nodes
            gsVector<T> wgrid1;          // tensor Gauss weights
            gsVector<T> lower(dim);
            gsVector<T> upper(dim);
            int stride = dim*(dim+1)/2;
            gsVector<int> num_nodes(dim);
            for ( index_t i = 0; i!= num_nodes.size(); ++i )
                num_nodes[i] = 2; // this should be  ceil( (degree in dir i) / 2 )
            gsVector<unsigned> curElement = gsVector<unsigned>::Zero(dim);
            do
              {
            for (int i = 0; i < dim; ++i)
              {
                lower[i] = breaks[i][curElement[i]];
                upper[i] = breaks[i][curElement[i]+1];
              }
            tensorGaussRule<T>(ngrid1, wgrid1, num_nodes, lower, upper);
            gsMatrix<T> der2;
            m_basis->deriv2_into(ngrid1, der2);
            gsMatrix<unsigned> act;
            m_basis->active_into(ngrid1, act);
            for (index_t k=0; k!= der2.cols(); ++k)
              {
                for (index_t i=0; i!=act.rows(); ++i)
                  for (index_t j=0; j!=act.rows(); ++j)
                {
                  m_A( act(i,k) , act(j,k) ) += wgrid1[k] * lambda *
                    (
                        der2(i*stride  , k) * der2(j*stride  , k) +  // d^2u N_i * d^2u N_j
                        der2(i*stride+1, k) * der2(j*stride+1, k) +  // d^2v N_i * d^2v N_j
                    2 * der2(i*stride+2, k) * der2(j*stride+2, k)   // dudv N_i * dudv N_j
                     );
                      }
            }
          }
          while (nextLexicographic(curElement, meshSize));
    gsDebug << "m_A size: "<< m_A.rows() <<" x "<< m_A.cols() <<"\n";
     }

    gsMatrix<T> x (m_B.rows(), m_B.cols());
    x=m_A.fullPivHouseholderQr().solve( m_B);
    //Eigen::BiCGSTAB< gsSparseMatrix<T>,  Eigen::IncompleteLUT<T> > solver( m_A );

    // Solves for many right hand side  columns
    //x =  solver.solve( m_B); //toDense()
    // finally generate the B-spline curve
    m_result = m_basis->makeGeometry(&x);
}

*/
