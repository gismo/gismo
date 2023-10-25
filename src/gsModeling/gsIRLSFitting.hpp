/** @file gsIRLSFitting.hpp

    @brief Provides implementation of data fitting algorithms by iterative re-weighted least
    squares approximation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore

*/

#include <gsCore/gsBasis.h>
#include <gsCore/gsGeometry.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsBSpline.h>
#include <gsTensor/gsTensorDomainIterator.h>
#include <gsModeling/gsIRLSFitting.h>


namespace gismo
{

template<class T>
gsIRLSFitting<T>::~gsIRLSFitting()
{
    if ( m_result )
        delete m_result;
}

template<class T>
gsIRLSFitting<T>::  gsIRLSFitting(gsMatrix<T> const & param_values,
                                  gsMatrix<T> const & points,
                                  gsBasis<T>  & basis,
                                  std::vector<T> & weights)
{
    GISMO_ASSERT(points.cols()==param_values.cols(), "Pointset dimensions problem "<< points.cols() << " != " <<param_values.cols() );
    GISMO_ASSERT(points.cols()==weights.size(), "Weight and observation dimensions problem "<< points.cols() << " != " << weights.size() );
    GISMO_ASSERT(basis.domainDim()==param_values.rows(), "Parameter values are inconsistent: "<< basis.domainDim() << " != " <<param_values.rows() );

    m_param_values = param_values;
    m_points = points;
    m_result = NULL;
    m_basis = &basis;
    m_weights = weights; // weights for the weighted least square fitting: not const since they can be changed within an iterative procedure
    m_points.transposeInPlace(); // points are now on the rows
}

template<class T>
void gsIRLSFitting<T>::compute(T lambda)
{
    m_last_lambda = lambda;

    // Wipe out previous result
    if ( m_result )
        delete m_result;

    const int num_basis=m_basis->size();
    const short_t dimension=m_points.cols();

    //left side matrix
    //gsMatrix<T> A_mat(num_basis,num_basis);
    gsSparseMatrix<T> A_mat(num_basis + m_constraintsLHS.rows(), num_basis + m_constraintsLHS.rows());
    //gsMatrix<T>A_mat(num_basis,num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_basis->dim(); ++i) // to do: improve
        // nonZerosPerCol *= m_basis->degree(i) + 1;
        nonZerosPerCol *= ( 2 * m_basis->degree(i) + 1 ) * 4;
    // TODO: improve by taking constraints nonzeros into account.
    A_mat.reservePerColumn( nonZerosPerCol );

    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis + m_constraintsRHS.rows(), dimension);
    m_B.setZero(); // enusure that all entries are zero in the beginning

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b
    // weighted least square: B'WB * x = B'W * b

    assembleSystem(A_mat, m_weights, m_B); // weighted least square fitting

    // --- Smoothing matrix computation
    //test degree >=3
    if(lambda > 0)
      applySmoothing(lambda, A_mat);

    if(m_constraintsLHS.rows() > 0)
	extendSystem(A_mat, m_B);

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)

    //gsDebugVar( A_mat.nonZerosPerCol().maxCoeff() );
    //gsDebugVar( A_mat.nonZerosPerCol().minCoeff() );
    A_mat.makeCompressed();

    typename gsSparseSolver<T>::BiCGSTABILUT solver( A_mat );

    if ( solver.preconditioner().info() != gsEigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";
        m_result = NULL;
        return;
    }
    // Solves for many right hand side  columns
    gsMatrix<T> x;

    x = solver.solve(m_B); //toDense()

    // If there were constraints, we obtained too many coefficients.
    x.conservativeResize(num_basis, gsEigen::NoChange);

    //gsMatrix<T> x (m_B.rows(), m_B.cols());
    //x=A_mat.fullPivHouseholderQr().solve( m_B);
    // Solves for many right hand side  columns
    // finally generate the B-spline curve
    m_result = m_basis->makeGeometry( give(x) ).release();
}

template <class T>
void gsIRLSFitting<T>::parameterCorrection(T accuracy,
                                       index_t maxIter,
                                       T tolOrth)
{
    if ( !m_result )
        compute(m_last_lambda);

    const index_t d = m_param_values.rows();
    const index_t n = m_points.cols();
    T maxAng, avgAng;
    std::vector<gsMatrix<T> > vals;
    gsMatrix<T> DD, der;
    for (index_t it = 0; it!=maxIter; ++it)
    {
        maxAng = -1;
        avgAng = 0;
        //auto der = Eigen::Map<typename gsMatrix<T>::Base, 0, Eigen::Stride<-1,-1> >
        //(vals[1].data()+k, n, m_points.rows(), Eigen::Stride<-1,-1>(d*n,d) );

#       pragma omp parallel for default(shared) private(der,DD,vals)
        for (index_t s = 0; s<m_points.rows(); ++s)
            //for (index_t s = 1; s<m_points.rows()-1; ++s) //(! curve) skip first and last point
        {
            vals = m_result->evalAllDers(m_param_values.col(s), 1);
            for (index_t k = 0; k!=d; ++k)
            {
                der = vals[1].reshaped(d,n);
                DD = vals[0].transpose() - m_points.row(s);
                const T cv = ( DD.normalized() * der.row(k).transpose().normalized() ).value();
                const T a = math::abs(0.5*EIGEN_PI-math::acos(cv));
#               pragma omp critical (max_avg_ang)
                {
                    maxAng = math::max(maxAng, a );
                    avgAng += a;
                }
            }
            /*
            auto der = Eigen::Map<typename gsMatrix<T>::Base, 0, Eigen::Stride<-1,-1> >
                (vals[1].data()+k, n, m_points.rows(), Eigen::Stride<-1,-1>(d*n,d) );
            maxAng = ( DD.colwise().normalized() *
                       der.colwise().normalized().transpose()
                ).array().acos().maxCoeff();
            */
        }

        avgAng /= d*m_points.rows();
        //gsInfo << "Avg-deviation: "<< avgAng << " / max: "<<maxAng<<"\n";

        // if (math::abs(0.5*EIGEN_PI-maxAng) <= tolOrth ) break;

        gsVector<T> newParam;
#       pragma omp parallel for default(shared) private(newParam)
        for (index_t i = 0; i<m_points.rows(); ++i)
        //for (index_t i = 1; i<m_points.rows()-1; ++i) //(!curve) skip first last pt
        {
            newParam = m_param_values.col(i);
            m_result->closestPointTo(m_points.row(i).transpose(),
                                     newParam, accuracy, true);

            // (!) There might be the same parameter for two points
            // or ordering constraints in the case of structured/grid data
            m_param_values.col(i) = newParam;
        }

        // refit
        compute(m_last_lambda);
    }
}


template <class T>
void gsIRLSFitting<T>::assembleSystem(gsSparseMatrix<T>& A_mat,
                                      std::vector<T>& m_weights,
                                      gsMatrix<T>& m_B)
{
    const int num_points = m_points.rows();

    //for computing the value of the basis function
    gsMatrix<T> value, curr_point;
    gsMatrix<index_t> actives;

#   pragma omp parallel for default(shared) private(curr_point,actives,value)
    for(index_t k = 0; k != num_points; ++k)
    {
        curr_point = m_param_values.col(k);

        //computing the values of the basis functions at the current point
        m_basis->eval_into(curr_point, value);

        //associate the weight to the basis functions:
        // value = curr_weight * value;

        // which functions have been computed i.e. which are active
        m_basis->active_into(curr_point, actives);

        const index_t numActive = actives.rows();

        for (index_t i = 0; i != numActive; ++i)
        {
            const index_t ii = actives.at(i);
#           pragma omp critical (acc_m_B)
            m_B.row(ii) += value.at(i) * m_points.row(k) * m_weights[k];
            for (index_t j = 0; j != numActive; ++j)
#               pragma omp critical (acc_A_mat)
                A_mat(ii, actives.at(j)) += m_weights[k] * value.at(i) * value.at(j);
        }
    }
}

template <class T>
void gsIRLSFitting<T>::extendSystem(gsSparseMatrix<T>& A_mat,
				gsMatrix<T>& m_B)
{
    index_t basisSize = m_basis->size();

    // Idea: maybe we can resize here instead of doing it in compute().

    // This way does not work, as these block operations are read-only for sparse matrices.
    /*A_mat.topRightCorner(m_constraintsLHS.cols(), m_constraintsLHS.rows()) = m_constraintsLHS.transpose();
      A_mat.bottomLeftCorner(m_constraintsLHS.rows(), m_constraintsLHS.cols()) = m_constraintsLHS;*/
    m_B.bottomRows(m_constraintsRHS.rows()) = m_constraintsRHS;

    for (index_t k=0; k<m_constraintsLHS.outerSize(); ++k)
    {
	for (typename gsSparseMatrix<T>::InnerIterator it(m_constraintsLHS,k); it; ++it)
	{
	    A_mat(basisSize + it.row(), it.col()) = it.value();
	    A_mat(it.col(), basisSize + it.row()) = it.value();
	}
    }
}

template<class T>
gsSparseMatrix<T> gsIRLSFitting<T>::smoothingMatrix(T lambda) const
{
    m_last_lambda = lambda;

    const int num_basis=m_basis->size();

    gsSparseMatrix<T> A_mat(num_basis + m_constraintsLHS.rows(), num_basis + m_constraintsLHS.rows());
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_basis->dim(); ++i) // to do: improve
        nonZerosPerCol *= ( 2 * m_basis->degree(i) + 1 );
    A_mat.reservePerColumn( nonZerosPerCol );
    const_cast<gsIRLSFitting*>(this)->applySmoothing(lambda, A_mat);
    return A_mat;
}

template<class T>
void gsIRLSFitting<T>::applySmoothing(T lambda, gsSparseMatrix<T> & A_mat)
{
    m_last_lambda = lambda;

    const short_t dim = m_basis->dim();
    const short_t stride = dim*(dim+1)/2;

    gsVector<index_t> numNodes(dim);
    for ( short_t i = 0; i!= dim; ++i )
        numNodes[i] = this->m_basis->degree(i);//+1;
    gsGaussRule<T> QuRule( numNodes ); // Reference Quadrature rule
    gsMatrix<T> quNodes, der2, localA;
    gsVector<T> quWeights;
    gsMatrix<index_t> actives;

    typename gsBasis<T>::domainIter domIt = m_basis->makeDomainIterator();

#   ifdef _OPENMP
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
    for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#   else
    for (; domIt->good(); domIt->next() )
#   endif
    {
        // Map the Quadrature rule to the element and compute basis derivatives
        QuRule.mapTo( domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights );
        m_basis->deriv2_into(quNodes, der2);
        m_basis->active_into(domIt->center, actives);
        const index_t numActive = actives.rows();
        localA.setZero(numActive, numActive );

        // perform the quadrature
        for (index_t k = 0; k < quWeights.rows(); ++k)
        {
            const T weight = quWeights[k] * lambda;

            for (index_t i=0; i!=numActive; ++i)
                for (index_t j=0; j!=numActive; ++j)
                {
                    T localAij = 0; // temporary variable

                    for (short_t s = 0; s < stride; s++)
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
                            localAij += (T)2 * der2(i * stride + s, k) *
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
            const int ii = actives(i,0);
            for (index_t j=0; j!=numActive; ++j)
                A_mat( ii, actives(j,0) ) += localA(i,j);
        }
    }
}

template<class T>
void gsIRLSFitting<T>::computeErrors()
{
    m_pointErrors.clear();
    m_pointWErrors.clear();

    gsMatrix<T> val_i;
    //m_result->eval_into(m_param_values.col(0), val_i);
    m_result->eval_into(m_param_values, val_i);
    m_pointErrors.push_back(  (m_points.row(0) - val_i.col(0).transpose()).norm() );
    m_pointWErrors.push_back(  ( (m_points.row(0) - val_i.col(0).transpose()).norm() ) * m_weights[0] ); // weighted error computation
    gsInfo << "Initial weight: " << m_weights[0] << "\n";
    m_weights[0] = 1/( (m_points.row(0) - val_i.col(0).transpose()).norm() );
    gsInfo << "Update weight: " << m_weights[0] << "\n";
    m_max_error = m_min_error = m_pointErrors.back();

    gsInfo << "Initial last weighted-error: " << m_pointWErrors.back() << "\n";

    for (index_t i = 1; i < m_points.rows(); i++)
    {
        //m_result->eval_into(m_param_values.col(i), val_i);

        const T err = (m_points.row(i) - val_i.col(i).transpose()).norm();
        const T werr = ( (m_points.row(i) - val_i.col(i).transpose()).norm()) * m_weights[i] ;

        m_pointErrors.push_back(err);
        m_pointWErrors.push_back(werr); // weighted errors


        m_weights[i] = 1/( (m_points.row(i) - val_i.col(i).transpose()).norm() );

        if ( err > m_max_error ) m_max_error = err;
        if ( err < m_min_error ) m_min_error = err;
    }

    gsInfo << "Update last weighted-error: " << m_pointWErrors.back() << "\n";



}


template<class T>
void gsIRLSFitting<T>::computeMaxNormErrors()
{
    m_pointErrors.clear();

    gsMatrix<T> values;
    m_result->eval_into(m_param_values, values);

    for (index_t i = 0; i != m_points.rows(); i++)
    {
        const T err = (m_points.row(i) - values.col(i).transpose()).cwiseAbs().maxCoeff();

        m_pointErrors.push_back(err);

        if ( i == 0 || m_max_error < err ) m_max_error = err;
        if ( i == 0 || err < m_min_error ) m_min_error = err;
    }
}



template<class T>
void gsIRLSFitting<T>::computeApproxError(T& error, int type) const
{
    gsMatrix<T> results;
    m_result->eval_into(m_param_values, results);
    error = 0;

    //computing the approximation error = sum_i ||x(u_i)-p_i||^2

    for (index_t i = 0; i != m_points.rows(); ++i)
    {
        const T err = (m_points.row(i) - results.col(i).transpose()).squaredNorm();

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
void gsIRLSFitting<T>::get_Error(std::vector<T>& errors, int type) const
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
            err += math::pow(m_points(row, col) - results(row, col), 2);
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

template<class T>
void gsIRLSFitting<T>::setConstraints(const std::vector<index_t>& indices,
				  const std::vector<gsMatrix<T> >& coefs)
{
    if(coefs.size() == 0)
	return;

    GISMO_ASSERT(indices.size() == coefs.size(),
		 "Prescribed indices and coefs must have the same length.");

    gsSparseMatrix<T> lhs(indices.size(), m_basis->size());
    gsMatrix<T> rhs(indices.size(), coefs[0].cols());

    index_t duplicates = 0;
    for(size_t r=0; r<indices.size(); r++)
    {
        index_t fix = indices[r];
        lhs(r-duplicates, fix) = 1;
        rhs.row(r-duplicates) = coefs[r];
    }

    setConstraints(lhs, rhs);
}

template<class T>
void gsIRLSFitting<T>::setConstraints(const std::vector<boxSide>& fixedSides)
{
    if(fixedSides.size() == 0)
    return;

    std::vector<index_t> indices;
    std::vector<gsMatrix<T> > coefs;

    for(std::vector<boxSide>::const_iterator it=fixedSides.begin(); it!=fixedSides.end(); ++it)
    {
    gsMatrix<index_t> ind = this->m_basis->boundary(*it);
    for(index_t r=0; r<ind.rows(); r++)
    {
        index_t fix = ind(r,0);
        // If it is a new constraint, add it.
        if(std::find(indices.begin(), indices.end(), fix) == indices.end())
        {
        indices.push_back(fix);
        coefs.push_back(this->m_result->coef(fix));
        }
    }
    }

    setConstraints(indices, coefs);
}

template<class T>
void gsIRLSFitting<T>::setConstraints(const std::vector<boxSide>& fixedSides,
                      const std::vector<gsBSpline<T> >& fixedCurves)
{
    std::vector<gsBSpline<T> > tmp = fixedCurves;
    std::vector<gsGeometry<T> *> fixedCurves_input(tmp.size());
    for (size_t k=0; k!=fixedCurves.size(); k++)
        fixedCurves_input[k] = dynamic_cast<gsGeometry<T> *>(&(tmp[k]));
    setConstraints(fixedSides, fixedCurves_input);
}

template<class T>
void gsIRLSFitting<T>::setConstraints(const std::vector<boxSide>& fixedSides,
                      const std::vector<gsGeometry<T> * >& fixedCurves)
{
    if(fixedSides.size() == 0)
    return;

    GISMO_ASSERT(fixedCurves.size() == fixedSides.size(),
         "fixedCurves and fixedSides are of different sizes.");

    std::vector<index_t> indices;
    std::vector<gsMatrix<T> > coefs;
    for(size_t s=0; s<fixedSides.size(); s++)
    {
    gsMatrix<T> coefsThisSide = fixedCurves[s]->coefs();
    gsMatrix<index_t> indicesThisSide = m_basis->boundaryOffset(fixedSides[s],0);
    GISMO_ASSERT(coefsThisSide.rows() == indicesThisSide.rows(),
             "Coef number mismatch between prescribed curve and basis side.");

    for(index_t r=0; r<indicesThisSide.rows(); r++)
    {
        index_t fix = indicesThisSide(r,0);
        // If it is a new constraint, add it.
        if(std::find(indices.begin(), indices.end(), fix) == indices.end())
        {
        indices.push_back(fix);
        coefs.push_back(coefsThisSide.row(r));
        }
    }
    }

    setConstraints(indices, coefs);
}


} // namespace gismo