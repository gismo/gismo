/** @file gsFitting.hpp

    @brief Provides implementation of data fitting algorithms by least
    squares approximation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Kapl, G. Kiss, A. Mantzaflaris, D. Mokris

*/

#include <gsCore/gsBasis.h>
#include <gsCore/gsBoundary.h>
#include <gsCore/gsGeometry.h>
#include <gsCore/gsLinearAlgebra.h>
#include <gsNurbs/gsBSpline.h>
#include <gsTensor/gsTensorDomainIterator.h>

#ifdef gsParasolid_ENABLED

#include <gsParasolid/gsClosestPoint.h>

// TODO: We assume that m_result is a gsTHBSpline*, which is not always the case.
#include <gsHSplines/gsTHBSpline.h>

#endif  // gsParasolid_ENABLED

#include <gsIO/gsFileData.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

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
    GISMO_ASSERT(points.cols()==param_values.cols(), "Pointset dimensions problem "<< points.cols() << " != " <<param_values.cols() );
    GISMO_ASSERT(basis.domainDim()==param_values.rows(), "Parameter values are inconsistent: "<< basis.domainDim() << " != " <<param_values.rows() );

    m_param_values = param_values;
    m_points = points;
    m_result = nullptr;
    m_basis = &basis;
    m_points.transposeInPlace();

    m_offset.resize(2);
    m_offset[0] = 0;
    m_offset[1] = m_points.rows();
}

template<class T>
gsFitting<T>::gsFitting(gsMatrix<T> const& param_values,
                        gsMatrix<T> const& points,
                        gsVector<index_t>  offset,
                        gsMappedBasis<2,T> & mbasis)
{
    m_param_values = param_values;
    m_points = points;
    m_result = nullptr;
    m_basis = &mbasis;
    m_points.transposeInPlace();
    m_offset = give(offset);
}

template<class T>
void gsFitting<T>::compute(T lambda)
{
    m_last_lambda = lambda;

    // Wipe out previous result
    if ( m_result )
        delete m_result;

    const int num_basis = m_basis->size();
    const short_t dimension = m_points.cols();

    //left side matrix
    //gsMatrix<T> A_mat(num_basis,num_basis);
    gsSparseMatrix<T> A_mat(num_basis + m_constraintsLHS.rows(), num_basis + m_constraintsLHS.rows());
    //gsMatrix<T>A_mat(num_basis,num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_basis->domainDim(); ++i) // to do: improve
        // nonZerosPerCol *= m_basis->degree(i) + 1;
        nonZerosPerCol *= ( 2 * m_basis->basis(0).degree(i) + 1 ) * 4;
    // TODO: improve by taking constraints nonzeros into account.
    A_mat.reservePerColumn( nonZerosPerCol );

    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis + m_constraintsRHS.rows(), dimension);
    m_B.setZero();  // enusure that all entries are zero in the beginning

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b

    assembleSystem(A_mat, m_B);


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

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";

        return;
    }
    // Solves for many right hand side  columns
    gsMatrix<T> x;

    x = solver.solve(m_B); //toDense()

    // If there were constraints, we obtained too many coefficients.
    x.conservativeResize(num_basis, Eigen::NoChange);

    //gsMatrix<T> x (m_B.rows(), m_B.cols());
    //x=A_mat.fullPivHouseholderQr().solve( m_B);
    // Solves for many right hand side  columns
    // finally generate the B-spline curve

    if (const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis))
        m_result = bb->makeGeometry( give(x) ).release();
    else
        m_mresult = gsMappedSpline<2,T> ( *static_cast<gsMappedBasis<2,T>*>(m_basis),give(x));
}

template<class T>
void gsFitting<T>::computePen(T lambda1, index_t d1, T lambda2, index_t d2)
{
    // m_last_lambda = lambda; // what about this member?

    // Wipe out previous result
    if ( m_result )
        delete m_result;

    const int num_basis = m_basis->size();
    const short_t dimension = m_points.cols();

    //left side matrix
    //gsMatrix<T> A_mat(num_basis,num_basis);
    gsSparseMatrix<T> A_mat(num_basis + m_constraintsLHS.rows(), num_basis + m_constraintsLHS.rows());
    //gsMatrix<T>A_mat(num_basis,num_basis);
    //To optimize sparse matrix an estimation of nonzero elements per
    //column can be given here
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_basis->domainDim(); ++i) // to do: improve
        // nonZerosPerCol *= m_basis->degree(i) + 1;
        nonZerosPerCol *= ( 2 * m_basis->basis(0).degree(i) + 1 ) * 4;
    // TODO: improve by taking constraints nonzeros into account.
    A_mat.reservePerColumn( nonZerosPerCol );

    //right side vector (more dimensional!)
    gsMatrix<T> m_B(num_basis + m_constraintsRHS.rows(), dimension);
    m_B.setZero();  // enusure that all entries are zero in the beginning

    // building the matrix A and the vector b of the system of linear
    // equations A*x==b

    assembleSystem(A_mat, m_B); // A_mat is the collocation matrix


    // --- Smoothing matrix computation
    //test degree >=3
    if(lambda1 > 0 || lambda2 > 0)
    {
      // gsInfo << "lambda1 = " << lambda1 << ", lambda2 = " << lambda2 << ".\n";
      applyPenalization(lambda1, d1, lambda2, d2, A_mat);
    }

    if(m_constraintsLHS.rows() > 0)
	extendSystem(A_mat, m_B);

    //Solving the system of linear equations A*x=b (works directly for a right side which has a dimension with higher than 1)

    //gsDebugVar( A_mat.nonZerosPerCol().maxCoeff() );
    //gsDebugVar( A_mat.nonZerosPerCol().minCoeff() );
    A_mat.makeCompressed();

    typename gsSparseSolver<T>::BiCGSTABILUT solver( A_mat );

    if ( solver.preconditioner().info() != Eigen::Success )
    {
        gsWarn<<  "The preconditioner failed. Aborting.\n";

        return;
    }
    // Solves for many right hand side  columns
    gsMatrix<T> x;

    x = solver.solve(m_B); //toDense()

    // If there were constraints, we obtained too many coefficients.
    x.conservativeResize(num_basis, Eigen::NoChange);

    //gsMatrix<T> x (m_B.rows(), m_B.cols());
    //x=A_mat.fullPivHouseholderQr().solve( m_B);
    // Solves for many right hand side  columns
    // finally generate the B-spline curve

    if (const gsBasis<T> * bb = dynamic_cast<const gsBasis<T> *>(m_basis))
        m_result = bb->makeGeometry( give(x) ).release();
    else
        m_mresult = gsMappedSpline<2,T> ( *static_cast<gsMappedBasis<2,T>*>(m_basis),give(x));
}

#ifdef gsParasolid_ENABLED

template <class T>
void gsFitting<T>::parameterCorrection(T accuracy,
                                       index_t maxIter,
                                       T tolOrth)
{
    // Silence warnings, cf. https://stackoverflow.com/a/1486931/10348053
    (void)accuracy;
    (void)maxIter;
    (void)tolOrth;

    if ( !m_result )
        compute(m_last_lambda);

    extensions::gsPKSession::start();
    for (index_t it = 0; it!=maxIter; ++it)
    {
        // convert m_result to B-spline
        gsTHBSpline<2, T> resultTHB  = *static_cast<gsTHBSpline<2>*>(m_result);
        gsTensorBSpline<2, T> result;
        resultTHB.convertToBSpline(result);

        // Less efficient version:

        // gsVector<T, 3> point;
        // gsVector<T, 2> newParam;
        // for(index_t i=0; i<m_points.rows(); i++)
        // {
        //     point = m_points.row(i);
        //     extensions::gsClosestParam(result, point, newParam);
        //     m_param_values.col(i) = newParam;
        // }

        // More efficient version:
        extensions::gsClosestParam(result, m_points, m_param_values);

        // refit
        compute(m_last_lambda);
    }
    extensions::gsPKSession::stop();
}

#else // Parasolid not enabled

template <class T>
void gsFitting<T>::smoothParameterCorrection(T accuracy, index_t maxIter)
{
    if(!m_result)
	compute(m_last_lambda);

    // TODO: Find a better place for this.
    initParametricDomain();

    for (index_t it = 0; it!=maxIter; ++it)
    {
	gsInfo << "Correcting for the " << it << "-th time." << std::endl;

	gsMatrix<T> idealPars = m_param_values;
	gsVector<T> newParam;
	for (index_t i = 0; i<m_points.rows(); ++i)
	{
	    newParam = m_param_values.col(i);
	    m_result->closestPointTo(m_points.row(i).transpose(), newParam, accuracy, true);
	    idealPars.col(i) = newParam;
	}

	index_t deg = 3;
	index_t numKnots = 7;
	gsKnotVector<T> uKnots(m_uMin, m_uMax, numKnots, deg + 1);
	gsKnotVector<T> vKnots(m_vMin, m_vMax, numKnots, deg + 1);
	gsTensorBSplineBasis<2, T> rebasis(uKnots, vKnots);
	gsFitting<T> reparam(m_param_values, idealPars, rebasis);

	// TODO: How to choose proper smoothing?
	// TODO: It seems to work but not to cause that much difference.
	// Test on a more challenging example!
	reparam.compute(1e-5);

	// gsFileData<> fd;
	// fd << *reparam.result();
	// fd.dump("reparam");

	// Evaluate the reparametrization and replace m_param_values with it.
	gsMatrix<T> oldPars = m_param_values;
	reparam.result()->eval_into(oldPars, m_param_values);

	compute(m_last_lambda);
    }
}

template <class T>
void gsFitting<T>::parameterCorrection(T accuracy,
                                       index_t maxIter,
                                       T tolOrth)
{
    if ( !m_result )
        compute(m_last_lambda);

    const index_t d = m_param_values.rows();
    const index_t n = m_points.cols();

     for (index_t it = 0; it!=maxIter; ++it)
     {
        gsVector<T> newParam;
#       pragma omp parallel for default(shared) private(newParam)
        for (index_t i = 0; i<m_points.rows(); ++i)
        //for (index_t i = 1; i<m_points.rows()-1; ++i) //(!curve) skip first last pt
        {
            const auto & curr = m_points.row(i).transpose();
            newParam = m_param_values.col(i);
            m_result->closestPointTo(curr, newParam, accuracy, true);

            // Decide whether to accept the correction or to drop it
            if ((m_result->eval(newParam) - curr).norm()
                    < (m_result->eval(m_param_values.col(i))
                        - curr).norm())
                    m_param_values.col(i) = newParam;

            // (!) There might be the same parameter for two points
            // or ordering constraints in the case of structured/grid data
        }

        // refit
        compute(m_last_lambda);
    }
}


template <class T>
void gsFitting<T>::parameterCorrectionPen(T accuracy,
                                       index_t maxIter,
                                       T tolOrth,
                                       T lambda1,
                                       index_t d1,
                                       T lambda2,
                                       index_t d2)
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
        gsVector<T> newParam;
        gsMatrix<T> supp = this->result()->support();
#       pragma omp parallel for default(shared) private(newParam)
        for (index_t i = 0; i<m_points.rows(); ++i)
        //for (index_t i = 1; i<m_points.rows()-1; ++i) //(!curve) skip first last pt
        {
            const auto & curr = m_points.row(i).transpose();
            newParam = m_param_values.col(i);
            m_result->closestPointTo(curr, newParam, accuracy, true);

            // Decide whether to accept the correction or to drop it
            if ((m_result->eval(newParam) - curr).norm()
                    < (m_result->eval(m_param_values.col(i))
                        - curr).norm())
                    m_param_values.col(i) = newParam;
        }

        // refit
        computePen(lambda1, d1, lambda2, d2);
    }
}


// template<short_t d, class T>
// void gsFitting<d, T>::computeApproximation(index_t maxPcIter)
// {
//     // We run one fitting step and compute the errors
//     this->compute(m_lambda);
//     //parameter correction
//     this->parameterCorrection(1e-7, maxPcIter, 1e-4);//closestPoint accuracy, orthogonality tolerance
//     //Approximation error
//     this->computeErrors();
// }



template <class T>
bool gsFitting<T>::is_corner(gsMatrix<T> & p_domain,
                             gsVector<T> & param)
{
  bool corner_check = false;
  if( (math::abs(param(0) - p_domain(0,0)) < 1e-15) && (math::abs(param(1) - p_domain(1,0)) < 1e-15) ){
    // gsInfo << "param:\n" << param << "\n";
    // gsInfo << param(0,0) << " == " << p_domain(0,0) << "\n";
    // gsInfo << param(0,1) << " == " << p_domain(1,0) << "\n";
    corner_check = true;
  }
  else if( (math::abs(param(0) - p_domain(0,1)) < 1e-15) && (math::abs(param(1) - p_domain(1,0)) < 1e-15) ){
    corner_check = true;
  }
  else if( (math::abs(param(0) - p_domain(0,1)) < 1e-15) && (math::abs(param(1) - p_domain(1,1)) < 1e-15) ){
    corner_check = true;
  }
  else if( (math::abs(param(0) - p_domain(0,0)) < 1e-15) && (math::abs(param(1) - p_domain(1,1)) < 1e-15) ){
    corner_check = true;
  }
  else{
    corner_check = false;
  }
  return corner_check;
}

template <class T>
void gsFitting<T>::parameterCorrectionFixedBoundary(T accuracy,
                                                    index_t maxIter,
                                                    index_t sepIndex)
{
    if ( !m_result )
        compute(m_last_lambda);

    for (index_t it = 0; it!=maxIter; ++it)
    {

      gsVector<T> newParam;
      gsMatrix<T> geoSupport = m_result->support();
#     pragma omp parallel for default(shared) private(newParam)
      for (index_t i = 0; i < sepIndex; ++i)
      {
          // gsInfo << "Interior:\n" << newParam << "\n";
          const auto & curr = m_points.row(i).transpose();
          newParam = m_param_values.col(i);
          m_result->closestPointTo(curr, newParam, accuracy, true);

          // Decide whether to accept the correction or to drop it
          // gsInfo << "Correction:\n" << newParam << "\n";
          // gsInfo << "On the boundary ? " << "\n";
          // gsInfo << geoSupport(0,0) << "<" << newParam(0) << "<" << geoSupport(0,1) << "\n";
          // gsInfo << geoSupport(1,0) << "<" << newParam(1) << "<" << geoSupport(1,1) << "\n";
          // if( ( geoSupport(0,0) < newParam(0) ) && ( newParam(0) < geoSupport(0,1) ) && ( geoSupport(1,0) < newParam(1) ) && ( newParam(1) < geoSupport(1,1) ) )
          // {
            // gsInfo << "And we get here.\n";
          if ((m_result->eval(newParam) - curr).norm() < (m_result->eval(m_param_values.col(i))- curr).norm()){
            m_param_values.col(i) = newParam;
          }

          // (!) There might be the same parameter for two points
          // or ordering constraints in the case of structured/grid data
      }

      }

      // refit
      compute(m_last_lambda);
}

template <class T>
void gsFitting<T>::parameterCorrectionSepBoundary(T accuracy,
                                                  index_t maxIter,
                                                  index_t sepIndex)
{
    if ( !m_result )
        compute(m_last_lambda);

    for (index_t it = 0; it!=maxIter; ++it)
    {

      gsVector<T> newParam;
      gsMatrix<T> geoSupport = m_result->support();
#     pragma omp parallel for default(shared) private(newParam)
      for (index_t i = 0; i < sepIndex; ++i)
      //for (index_t i = 1; i<m_points.rows()-1; ++i) //(!curve) skip first last pt
      {
          // gsInfo << "Interior:\n" << newParam << "\n";
          const auto & curr = m_points.row(i).transpose();
          newParam = m_param_values.col(i);
          m_result->closestPointTo(curr, newParam, accuracy, true);

          // Decide whether to accept the correction or to drop it
          // gsInfo << "Correction:\n" << newParam << "\n";
          // gsInfo << "On the boundary ? " << "\n";
          // gsInfo << geoSupport(0,0) << "<" << newParam(0) << "<" << geoSupport(0,1) << "\n";
          // gsInfo << geoSupport(1,0) << "<" << newParam(1) << "<" << geoSupport(1,1) << "\n";
          // if( ( geoSupport(0,0) < newParam(0) ) && ( newParam(0) < geoSupport(0,1) ) && ( geoSupport(1,0) < newParam(1) ) && ( newParam(1) < geoSupport(1,1) ) )
          // {
            // gsInfo << "And we get here.\n";
          if ((m_result->eval(newParam) - curr).norm() < (m_result->eval(m_param_values.col(i))- curr).norm()){
            m_param_values.col(i) = newParam;
          }
          //}

          // (!) There might be the same parameter for two points
          // or ordering constraints in the case of structured/grid data
      }

      // Correct the parameters, but keep them on the boundary


      for (index_t i = sepIndex; i < m_points.rows(); ++i)
      {
          gsVector<T> newBoundaryParam(1);
          gsVector<T> oldBoundaryParam(1);
          const auto & curr = m_points.row(i).transpose();
          newParam = m_param_values.col(i);

          //if (!is_corner(geoSupport, newParam)){

            // gsInfo << "West  : " << newParam(0) << "==" << geoSupport(0,0) << "\n";
            // gsInfo << "East  : " << newParam(0) << "==" << geoSupport(0,1) << "\n";
            // gsInfo << "South : " << newParam(1) << "==" << geoSupport(1,0) << "\n";
            // gsInfo << "North : " << newParam(1) << "==" << geoSupport(1,1) << "\n";


            gsVector<T> leftParam = m_param_values.col(i);
            gsVector<T> rightParam = m_param_values.col(i);

            if ( (i==sepIndex) && !is_corner(geoSupport, newParam) ){
              leftParam = m_param_values.col(m_points.rows()-1);
              rightParam = m_param_values.col(i+1);
            }

            if ( (i==m_points.rows()-1) && !is_corner(geoSupport, newParam) ){
              leftParam = m_param_values.col(i-1);
              rightParam = m_param_values.col(sepIndex);
            }

            if ( (i>sepIndex) && (i < m_points.rows()-1) && !is_corner(geoSupport, newParam) ){

              leftParam = m_param_values.col(i-1);
              rightParam = m_param_values.col(i+1);
            }

            if( math::abs(newParam(0) - geoSupport(0,0)) < 1e-15 )
            {
              // gsInfo << newParam(0) << " == " << geoSupport(0,0) << "\n";
              // gsInfo << "1 = West-boundary:\n" << newParam << "\n";
              typename gsGeometry<T>::uPtr b = m_result->boundary(1);


              newBoundaryParam << newParam(1);
              oldBoundaryParam << newParam(1);
              // gsInfo << "param to optimize: " << newBoundaryParam << "\n";
              // gsInfo << "West-boundary correction:\n";
              b->closestPointTo(curr, newBoundaryParam, accuracy, true);
              // gsInfo << "optimized: " << newBoundaryParam << "\n";

              // assumption: boundary parameters in anti-clockwise order
              if( (leftParam(1) > newBoundaryParam(0)) && (newBoundaryParam(0) > rightParam(1)) ){
                newParam(1) = newBoundaryParam(0);
              }
              if ((b->eval(newBoundaryParam) - curr).norm() < (b->eval(oldBoundaryParam)- curr).norm()){
                m_param_values.col(i) = newParam;
              }

            }
            else if ( math::abs(newParam(0) - geoSupport(0,1)) < 1e-15 )
            {
              // gsInfo << newParam(0) << " == " << geoSupport(0,1) << "\n";
              // gsInfo << "2 = East-boundary:\n" << newParam << "\n";
              typename gsGeometry<T>::uPtr b = m_result->boundary(2);


              newBoundaryParam << newParam(1);
              oldBoundaryParam << newParam(1);
              // gsInfo << "param to optimize: " << newBoundaryParam << "\n";
              // gsInfo << "East-boundary correction:\n";
              b->closestPointTo(curr, newBoundaryParam, accuracy, true);
              // gsInfo << "optimized: " << newBoundaryParam << "\n";
              if( (leftParam(1) < newBoundaryParam(0)) && (newBoundaryParam(0) < rightParam(1)) ){
                newParam(1) = newBoundaryParam(0);
              }

              if ((b->eval(newBoundaryParam) - curr).norm() < (b->eval(oldBoundaryParam)- curr).norm()){
                m_param_values.col(i) = newParam;
              }

            }
            else if ( math::abs(newParam(1) - geoSupport(1,0)) < 1e-15 )
            {
              // gsInfo << newParam(1) << " == " << geoSupport(1,0) << "\n";
              // gsInfo << "3 = South-boundary:\n" << newParam << "\n";
              typename gsGeometry<T>::uPtr b = m_result->boundary(3);


              newBoundaryParam << newParam(0);
              oldBoundaryParam << newParam(0);
              // gsInfo << "param to optimize: " << newBoundaryParam << "\n";
              // gsInfo << "South-boundary correction:\n";
              b->closestPointTo(curr, newBoundaryParam, accuracy, true);
              // gsInfo << "optimized: " << newBoundaryParam << "\n";
              if( (leftParam(0) < newBoundaryParam(0)) && (newBoundaryParam(0) < rightParam(0)) ){
                newParam(0) = newBoundaryParam(0);
              }

              if ((b->eval(newBoundaryParam) - curr).norm() < (b->eval(oldBoundaryParam)- curr).norm()){
                m_param_values.col(i) = newParam;
              }
            }
            else{
              // gsInfo << newParam(1) << " == " << geoSupport(1,1) << "\n";
              // gsInfo << i << "-th point:\n" << "\n";
              // gsInfo << rightParam(0)<< " < " << newParam(0) << " < " << leftParam(0) << "\n";
              typename gsGeometry<T>::uPtr b = m_result->boundary(4);

              // gsInfo << newParam(1) << "==" << geoSupport(1,1) << "\n";
              newBoundaryParam << newParam(0);
              oldBoundaryParam << newParam(0);
              // gsInfo << "param to optimize: " << newBoundaryParam << "\n";
              // gsInfo << "North-boundary correction:\n";
              b->closestPointTo(curr, newBoundaryParam, accuracy, true);
              // gsInfo << "optimized: " << newBoundaryParam << "\n";
              // assumption: paramters in anti-clockwise order
              // gsInfo << "Projection:\n";
              // gsInfo << rightParam(0)<< " < " << newBoundaryParam(0) << " < " << leftParam(0) << "\n";
              if( (leftParam(0) > newBoundaryParam(0)) && (newBoundaryParam(0) > rightParam(0)) ){
                // gsInfo << "Before:\n" << newParam << "\n";
                newParam(0) = newBoundaryParam(0);
                // gsInfo << "After:\n" << newParam << "\n";
              }
              if ((b->eval(newBoundaryParam) - curr).norm() < (b->eval(oldBoundaryParam)- curr).norm()){
                m_param_values.col(i) = newParam;
              }
            }


          //} // corner: to be moved or not.
      }

      // refit
      compute(m_last_lambda);
    }
}

#endif // gsParasolid_ENABLED


template <class T>
void gsFitting<T>::assembleSystem(gsSparseMatrix<T>& A_mat,
                                  gsMatrix<T>& m_B)
{
    const int num_patches ( m_basis->nPieces() ); //initialize

    //for computing the value of the basis function
    gsMatrix<T> value, curr_point;
    gsMatrix<index_t> actives;

    for (index_t h = 0; h < num_patches; h++ )
    {
        auto & basis = m_basis->basis(h);

//#   pragma omp parallel for default(shared) private(curr_point,actives,value)
        for (index_t k = m_offset[h]; k < m_offset[h+1]; ++k)
        {
            curr_point = m_param_values.col(k);

            //computing the values of the basis functions at the current point
            basis.eval_into(curr_point, value);

            // which functions have been computed i.e. which are active
            basis.active_into(curr_point, actives);

            const index_t numActive = actives.rows();

            for (index_t i = 0; i != numActive; ++i)
            {
                const index_t ii = actives.at(i);
//#           pragma omp critical (acc_m_B)
                m_B.row(ii) += value.at(i) * m_points.row(k);
                for (index_t j = 0; j != numActive; ++j)
//#               pragma omp critical (acc_A_mat)
                    A_mat(ii, actives.at(j)) += value.at(i) * value.at(j);
            }
        }
    }
}

template <class T>
void gsFitting<T>::extendSystem(gsSparseMatrix<T>& A_mat,
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
gsSparseMatrix<T> gsFitting<T>::smoothingMatrix(T lambda) const
{
    m_last_lambda = lambda;

    const int num_basis=m_basis->size();

    gsSparseMatrix<T> A_mat(num_basis + m_constraintsLHS.rows(), num_basis + m_constraintsLHS.rows());
    int nonZerosPerCol = 1;
    for (int i = 0; i < m_basis->domainDim(); ++i) // to do: improve
        nonZerosPerCol *= ( 2 * m_basis->basis(0).degree(i) + 1 );
    A_mat.reservePerColumn( nonZerosPerCol );
    const_cast<gsFitting*>(this)->applySmoothing(lambda, A_mat);
    return A_mat;
}

template<class T>
void gsFitting<T>::applySmoothing(T lambda, gsSparseMatrix<T> & A_mat)
{
    m_last_lambda = lambda;
    const int num_patches(m_basis->nPieces()); //initialize
    const short_t dim(m_basis->domainDim());
    const short_t stride = dim * (dim + 1) / 2;

    gsVector<index_t> numNodes(dim);
    gsMatrix<T> quNodes, der2, localA;
    gsVector<T> quWeights;
    gsMatrix<index_t> actives;

#   ifdef _OPENMP
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#   endif

    for (index_t h = 0; h < num_patches; h++)
    {
        auto & basis = m_basis->basis(h);

        //gsDebugVar(dim);
        //gsDebugVar(stride);

        for (short_t i = 0; i != dim; ++i)
        {
            numNodes[i] = basis.degree(i);//+1;
        }

        gsGaussRule<T> QuRule(numNodes); // Reference Quadrature rule

        typename gsBasis<T>::domainIter domIt = basis.makeDomainIterator();


#       ifdef _OPENMP
        for ( domIt->next(tid); domIt->good(); domIt->next(nt) )
#       else
        for (; domIt->good(); domIt->next() )
#       endif
        {
            // Map the Quadrature rule to the element and compute basis derivatives
            QuRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);
            basis.deriv2_into(quNodes, der2);
            basis.active_into(domIt->center, actives);
            const index_t numActive = actives.rows();
            localA.setZero(numActive, numActive);

            // perform the quadrature
            for (index_t k = 0; k < quWeights.rows(); ++k)
            {
                const T weight = quWeights[k] * lambda;

                for (index_t i = 0; i != numActive; ++i)
                    for (index_t j = 0; j != numActive; ++j)
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
                                localAij += T(2) * der2(i * stride + s, k) *
                                    der2(j * stride + s, k);
                            }
                        }
                        localA(i, j) += weight * localAij;
                    }
            }

            for (index_t i = 0; i != numActive; ++i)
            {
                const int ii = actives(i, 0);
                for (index_t j = 0; j != numActive; ++j)
                    A_mat(ii, actives(j, 0)) += localA(i, j);
            }
        }

    }
}



template<class T>
gsSparseMatrix<T> diff(const gsSparseMatrix<T> & E)
{
  gsSparseMatrix<T> E1 = E.block(0, 0, E.rows()-1, E.cols());
  gsSparseMatrix<T> E2 = E.block(1, 0, E.rows()-1, E.cols());
  return E2 - E1;
}

template<class T>
gsSparseMatrix<T> diff(const gsSparseMatrix<T> & E, const index_t d)
{
  GISMO_ASSERT( E.rows() <= d, "Wrong regulatization amount.");
  if (d == 0)
    return E;
  else
  {
    index_t r = 1;
    gsSparseMatrix<T> D = E;
    while( r <= d)
    {
      D = diff(D);
      r += 1;
    }
    return D;
  }
}

template<class T>
gsSparseMatrix<T> gsKroneckerMatProduct(const gsSparseMatrix<T>& A, const gsSparseMatrix<T>& B)
{
  gsMatrix<T> C(B.rows() * A.rows(), B.cols() * A.cols());
  C.setZero();

  for(index_t i = 0; i < A.rows(); i++)
  {
    for(index_t j = 0; j < A.cols(); j++)
    {
      // gsInfo << A(i,j) << " *\n"<< B << "\n";
      C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i,j) * B;
    }
  }

  // gsInfo << C << "\n";

  return C.sparseView();
}



template<class T>
void gsFitting<T>::applyPenalization(T lambda1, index_t d1, T lambda2, index_t d2, gsSparseMatrix<T> & A_mat)
{

    const int num_patches(m_basis->nPieces()); //initialize
    const short_t dim(m_basis->domainDim());
    const short_t stride = dim * (dim + 1) / 2;

    gsMatrix<T> localP1;
    gsMatrix<T> localP2;

    gsMatrix<index_t> actives;

#   ifdef _OPENMP
    const int tid = omp_get_thread_num();
    const int nt  = omp_get_num_threads();
#   endif

    for (index_t h = 0; h < num_patches; h++)
    {
        auto & basis = m_basis->basis(h);

        index_t n1 = basis.component(0).size(); // n
        if(n1 <= d1)
        {
          gsInfo << "Too high penalization in 0-direction, set " << d1 << " to " << n1-1 << "\n";
          d1 = n1 - 1;
        }

        index_t n2 = basis.component(1).size(); // n-tilde
        if(n2 <= d2)
        {
          gsInfo << "Too high penalization in 1-direction, set " << d2 << " to " << n2-1 << "\n";
          d2 = n2 - 1;
        }

        gsSparseMatrix<T> I1(n1, n1);
        I1.setIdentity();
        // gsInfo << "I1 =\n" << I1 << "\n";
        gsSparseMatrix<T> D1 = diff(I1, d1);
        // gsInfo << "D1 =\n" << D1 << "\n";



        gsSparseMatrix<T> I2(n2, n2);
        I2.setIdentity();
        gsSparseMatrix<T> D2 = diff(I2, d2);

        gsSparseMatrix<T> P1_left = gsKroneckerMatProduct(I2,D1);
        gsSparseMatrix<T> P1 = P1_left.transpose() * P1_left;

        // gsInfo << "P1_left =\n" << P1_left << "\n";

        gsSparseMatrix<T> P2_left = gsKroneckerMatProduct(I1,D2);
        gsSparseMatrix<T> P2 = P2_left.transpose() * P2_left;

        // gsInfo << "Collocation matrix =\n" << A_mat.rows() << " x " << A_mat.cols() << "\n";
        // gsInfo << "Apply penalizaion to patch n. " << h << ".\n";

        A_mat += lambda1 * P1 + lambda2 * P2;


    } // num patches
}




template<class T>
void gsFitting<T>::computeErrors()
{
    m_pointErrors.clear();

    gsMatrix<T> val_i;
    //->eval_into(m_param_values.col(0), val_i);
    m_result->eval_into(m_param_values, val_i);
    m_pointErrors.push_back( (m_points.row(0) - val_i.col(0).transpose()).norm() );
    m_max_error = m_min_error = m_pointErrors.back();

    for (index_t i = 1; i < m_points.rows(); i++)
    {
        //m_result->eval_into(m_param_values.col(i), val_i);

        const T err = (m_points.row(i) - val_i.col(i).transpose()).norm() ;

        m_pointErrors.push_back(err);

        if ( err > m_max_error ) m_max_error = err;
        if ( err < m_min_error ) m_min_error = err;
    }
}


template<class T>
void gsFitting<T>::computeMaxNormErrors()
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
void gsFitting<T>::computeApproxError(T& error, int type) const

{
    gsMatrix<T> curr_point, results;

    const int num_patches(m_basis->nPieces());

    error = 0;

    for (index_t h = 0; h < num_patches; h++)
    {

        for (index_t k = m_offset[h]; k < m_offset[h + 1]; ++k)
        {
            curr_point = m_param_values.col(k);

            if (m_result)
                m_result->eval_into(curr_point, results);
            else
            {
                m_mresult.eval_into(h, curr_point, results);
            }

                const T err = (m_points.row(k) - results.transpose()).squaredNorm();

                switch (type) {
                case 0:
                    error += err;
                    break;
                case 1:
                    error += math::sqrt(err);
                    break;
                default:
                    gsWarn << "Unknown type in computeApproxError(error, type)...\n";
                    break;
                }
        }
    }
}

template<class T>
void gsFitting<T>::get_Error(std::vector<T>& errors, int type) const
{
    errors.clear();
    errors.reserve(m_points.rows());

    gsMatrix<T> curr_point, results;

    T err = 0;

    const int num_patches(m_basis->nPieces());

    for (index_t h = 0; h < num_patches; h++)
    {
        for (index_t k = m_offset[h]; k < m_offset[h + 1]; ++k)
        {
            curr_point = m_param_values.col(k);

            if (m_result)
                m_result->eval_into(curr_point, results);
            else
            {
                m_mresult.eval_into(h, curr_point, results);
            }

            results.transposeInPlace();

            err = (m_points.row(k) - results).template lpNorm<Eigen::Infinity>();

                    switch (type)
                    {
                    case 0:
                        errors.push_back(err);
                        break;
                    case 1:
                        errors.push_back(math::sqrt(err));
                        break;
                    default:
                        gsWarn << "Unknown type in get_Error(errors, type)...\n";
                        break;
                    }
        }
    }
}

template<class T>
void gsFitting<T>::setConstraints(const std::vector<index_t>& indices,
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
void gsFitting<T>::setConstraints(const std::vector<boxSide>& fixedSides)
{
    if(fixedSides.size() == 0)
    return;

    std::vector<index_t> indices;
    std::vector<gsMatrix<T> > coefs;

    for(std::vector<boxSide>::const_iterator it=fixedSides.begin(); it!=fixedSides.end(); ++it)
    {
        gsMatrix<index_t> ind = this->m_basis->basis(0).boundary(*it);
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
void gsFitting<T>::setConstraints(const std::vector<boxSide>& fixedSides,
                      const std::vector<gsBSpline<T> >& fixedCurves)
{
    std::vector<gsBSpline<T> > tmp = fixedCurves;
    std::vector<gsGeometry<T> *> fixedCurves_input(tmp.size());
    for (size_t k=0; k!=fixedCurves.size(); k++)
        fixedCurves_input[k] = dynamic_cast<gsGeometry<T> *>(&(tmp[k]));
    setConstraints(fixedSides, fixedCurves_input);
}

template<class T>
void gsFitting<T>::setConstraints(const std::vector<boxSide>& fixedSides,
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
    gsMatrix<index_t> indicesThisSide = m_basis->basis(0).boundaryOffset(fixedSides[s],0);
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
