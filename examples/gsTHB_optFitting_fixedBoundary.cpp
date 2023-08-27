/** @file gsTHB_optFitting_fixedBoundary.cpp

    @brief Tensor product BSpline surface fitting with HLBFGS optimization.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <iostream>

#include <gismo.h>

#include <gsOptimizer/gsOptProblem.h>
#include <gsOptimizer/gsGradientDescent.h>

#ifdef gsHLBFGS_ENABLED
#include <gsHLBFGS/gsHLBFGS.h>
#include <gsIpOpt/gsIpOpt.h>
#endif


using namespace gismo;

template<class T>
void scalePoints(const gsMatrix<T> & xyz,
                 gsMatrix<T> & points)
{
  T p_min = xyz.minCoeff(),
         p_max = xyz.maxCoeff();
  T den = p_max - p_min;

  points.resize(xyz.rows(), xyz.cols());
  points = (1/den)*(xyz - p_min * gsMatrix<T>::Ones(xyz.rows(), xyz.cols()));
  gsWriteParaviewPoints(points, "scaled_points");
}



/* input : a parametrized point cloud with parameters and points
  output : parameters and points ordered by : interior (parameters/points) and
           boundary (parameters/points) ordered anticlockwise south-east-north-west edges,
           plus the 4 corner domains stored in a vector [c1, c2, c3, c4].
*/
template<class T>
void sortPointCloud(gsMatrix<T> & parameters,
                    gsMatrix<T> & points,
                    std::vector<index_t> & corners)
{
  // The following matrices and vectors store the parameters and points values and indeces.
  // There is no need to store these information, we could also use only one matrix and 1 std::vector and overwirte them each time.
  gsMatrix<T> uv_interiors, uv_south, uv_east, uv_north, uv_west;
  gsMatrix<T> p_interiors, p_south, p_east, p_north, p_west;
  std::vector<index_t> interiors, b_west, b_east, b_south, b_north;

  // Determine the parameter domain by mi/max of parameter values
  T u_min = parameters.row(0).minCoeff(),
    u_max = parameters.row(0).maxCoeff(),
    v_min = parameters.row(1).minCoeff(),
    v_max = parameters.row(1).maxCoeff();

  gsVector<T> curr_point(2,1);
  for(index_t i=0; i < parameters.cols(); i++)
  {
    curr_point = parameters.col(i);
    if( (u_min < curr_point(0)) && (curr_point(0) < u_max) && (v_min < curr_point(1)) && (curr_point(1) < v_max) )
      interiors.push_back(i);
    else // not interior point
    {
      if( (math::abs(curr_point(0) - u_min) < 1e-15) && (curr_point(1) > v_min) )
        b_west.push_back(i);//west edge
      else if( (math::abs(curr_point(0) - u_max) < 1e-15) && curr_point(1) < v_max)
        b_east.push_back(i);// east edge
      else if( (math::abs(curr_point(1) - v_min) < 1e-15) && (curr_point(0) < u_max) )
        b_south.push_back(i);// south edge
      else
        b_north.push_back(i);// north edge
    }
  }

  // gsInfo << "There are " << interiors.size() << " interior points.\n";
  // gsDebugVar(interiors.size());
  corners.push_back(interiors.size()); // c1
  // gsInfo << "There are " << b_south.size() << " south points.\n";
  // gsDebugVar(interiors.size() + b_south.size());
  corners.push_back(interiors.size() + b_south.size()); // c2
  // gsInfo << "There are " << b_east.size() << " east points.\n";
  // gsDebugVar(interiors.size() + b_south.size() + b_east.size());
  corners.push_back(interiors.size() + b_south.size() + b_east.size()); // c3
  // gsInfo << "There are " << b_north.size() << " north points.\n";
  // gsDebugVar(interiors.size() + b_south.size() + b_east.size() + b_north.size());
  corners.push_back(interiors.size() + b_south.size() + b_east.size() + b_north.size()); // c4
  // gsInfo << "There are " << b_west.size() << " west points.\n";

  uv_interiors.resize(2, interiors.size());
  p_interiors.resize(3, interiors.size());
  for( index_t i = 0; i < interiors.size(); i++ )
  {
    uv_interiors.col(i) = parameters.col(interiors[i]);
    p_interiors.col(i) = points.col(interiors[i]);
  }

  uv_west.resize(2, b_west.size());
  gsMatrix<T> tmp_west(3, b_west.size());
  for( index_t i = 0; i < b_west.size(); i++ )
  {
    uv_west.col(i) = parameters.col(b_west[i]);
    tmp_west.col(i) = points.col(b_west[i]);
  }

  uv_east.resize(2, b_east.size());
  gsMatrix<T> tmp_east(3, b_east.size());
  for( index_t i = 0; i < b_east.size(); i++ )
  {
    uv_east.col(i) = parameters.col(b_east[i]);
    tmp_east.col(i) = points.col(b_east[i]);
  }

  uv_south.resize(2, b_south.size());
  gsMatrix<T> tmp_south(3, b_south.size());
  for( index_t i = 0; i < b_south.size(); i++ )
  {
    uv_south.col(i) = parameters.col(b_south[i]);
    tmp_south.col(i) = points.col(b_south[i]);
  }

  uv_north.resize(2, b_north.size());
  gsMatrix<T> tmp_north(3, b_north.size());
  for( index_t i = 0; i < b_north.size(); i++ )
  {
    uv_north.col(i) = parameters.col(b_north[i]);
    tmp_north.col(i) = points.col(b_north[i]);
  }

  uv_south.transposeInPlace();
  uv_east.transposeInPlace();
  uv_north.transposeInPlace();
  uv_west.transposeInPlace();


  std::vector<index_t> tmp = uv_south.idxByColumn(0);
  // gsDebugVar(uv_south);
  p_south.resize(tmp_south.rows(), tmp_south.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_south.col(i) = tmp_south.col(tmp[i]);
  }
  uv_south.transposeInPlace();


  tmp.clear();
  tmp = uv_east.idxByColumn(1);
  // gsDebugVar(uv_east);
  p_east.resize(tmp_east.rows(), tmp_east.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_east.col(i) = tmp_east.col(tmp[i]);
  }
  uv_east.transposeInPlace();


  tmp.clear();
  tmp = uv_north.idxByColumn(0);
  std::reverse(tmp.begin(),tmp.end());
  // gsDebugVar(uv_north);
  gsVector<T> tcol = uv_north.col(0).reverse();
  uv_north.col(0) = tcol;
  tcol = uv_north.col(1).reverse();
  uv_north.col(1) = tcol;
  // gsDebugVar(uv_north);
  // for (std::vector<index_t>::iterator it = tmp.begin(); it != tmp.end(); ++it)
  //   gsInfo << *it <<"\n";
  p_north.resize(tmp_north.rows(), tmp_north.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_north.col(i) = tmp_north.col(tmp[i]);
  }
  uv_north.transposeInPlace();


  tmp.clear();
  tmp = uv_west.idxByColumn(1);
  // gsDebugVar(uv_west);
  tcol = uv_west.col(0).reverse();
  uv_west.col(0) = tcol;
  tcol = uv_west.col(1).reverse();
  uv_west.col(1) = tcol;
  // gsDebugVar(uv_west);
  std::reverse(tmp.begin(),tmp.end());

  p_west.resize(tmp_west.rows(), tmp_west.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_west.col(i) = tmp_west.col(tmp[i]);
  }
  uv_west.transposeInPlace();


  gsWriteParaviewPoints(uv_interiors, "uv_interiors");
  gsWriteParaviewPoints(p_interiors, "p_interiors");

  gsWriteParaviewPoints(uv_west, "uv_west");
  gsWriteParaviewPoints(tmp_west, "p_west");

  gsWriteParaviewPoints(uv_east, "uv_east");
  gsWriteParaviewPoints(tmp_east, "p_east");

  gsWriteParaviewPoints(uv_south, "uv_south");
  gsWriteParaviewPoints(tmp_south, "p_south");

  gsWriteParaviewPoints(uv_north, "uv_north");
  gsWriteParaviewPoints(tmp_north, "p_north");

  // reordering of the input point cloud (parameters and points)
  parameters.resize(uv_interiors.rows(), points.cols());
  parameters << uv_interiors.row(0), uv_south.row(0), uv_east.row(0), uv_north.row(0), uv_west.row(0),
                uv_interiors.row(1), uv_south.row(1), uv_east.row(1), uv_north.row(1), uv_west.row(1);

  points.resize(p_interiors.rows(), parameters.cols());
  points << p_interiors.row(0), p_south.row(0), p_east.row(0), p_north.row(0), p_west.row(0),
                p_interiors.row(1), p_south.row(1), p_east.row(1), p_north.row(1), p_west.row(1),
                p_interiors.row(2), p_south.row(2), p_east.row(2), p_north.row(2), p_west.row(2);

} // end sortPointCloud


//
//
// // To define an optimization problem we inherit from gsOptProblem class
// // and implement the default constructor and few inherited virtual functions
//
// //! [OptProblemExample Class]
// template <typename T>
// class gsOptProblemExample : public gsOptProblem<T>
// //! [OptProblemExample Class]
// {
// public:
//
//     //! [OptProblemExample Constructor]
//     // The constructor defines all properties of our optimization problem
//     gsOptProblemExample(gsHFitting<2, T>& mp, const gsMatrix<T> & params, const gsMatrix<T> & X)
//     :
//     m_mp(&mp),
//     m_params(params),
//     m_X(X)
//     {
//         // Number of design variables: how many variables we optimize, i.e. coefficiets + parametric values
//         // m_numDesignVars  = m_mp[0].coefs().size() + m_params.size(); // dim * spline-dofs + numPts
//         m_numDesignVars  = m_mp->result()->coefs().size() + m_params.size();
//
//         // design bounds
//         m_desLowerBounds.resize(m_numDesignVars);
//         m_desUpperBounds.resize(m_numDesignVars);
//
//         // coefficiets in R^3, namely no bounds
//         for(index_t i = 0; i < m_mp->result()->coefs().size(); i++)
//         {
//           m_desLowerBounds[i] = -1.0e19; // lower bound on the coefficients
//           m_desUpperBounds[i] =  1.0e19; // upper bound on the coefficients
//         }
//
//         // parameters in [0,1]^2
//         currentparams = m_params.transpose(); // m_mp.returnParamValues().transpose();
//         for(index_t i = 0; i < currentparams.size(); i++)
//         {
//           m_desLowerBounds[m_mp->result()->coefs().size() + i] = 0.; // lower bound on the parameters
//           m_desUpperBounds[m_mp->result()->coefs().size() + i] = 1.; // upper bound on the parameters
//         }
//
//         // Initialization of the smoothing matrix that we need to define th objective function.
//         m_G.resize(m_mp->result()->coefs().rows(), m_mp->result()->coefs().rows());
//         m_G.reservePerColumn( cast<T,index_t>( (2 * m_mp->result()->basis().maxDegree() + 1) * 1.333 ) );
//         m_mp->applySmoothing(m_mp->lambda(), m_G);
//
//         // design variables: whant we do optimize.
//         // c_x, c_y, c_z, u, v
//         m_curDesign.resize(m_numDesignVars,1);
//         m_curDesign << m_mp->result()->coefs().reshape(m_mp->result()->coefs().size(),1), currentparams.reshape(currentparams.size(),1);
//     }
//     //! [OptProblemExample Constructor]
//
// public:
//
//     //temporary storage
//     mutable gsMatrix<T> tmp, currentparams;
//     mutable gsSparseMatrix<T> c;
//     mutable std::vector<gsSparseMatrix<T>> c_matrices;
//
//     //! [OptProblemExample evalObj]
//     // The evaluation of the objective function must be implemented
//     T evalObj( const gsAsConstVector<T> & u ) const
//     {
//
//         gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() ); // keep point within design bounds
//         u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);
//
//         gsAsConstMatrix<T> currentcoefs(u.data(),  m_mp->result()->coefs().rows(), m_mp->result()->coefs().cols());
//         currentparams = gsAsConstMatrix<T>(u.data() + m_mp->result()->coefs().size(), m_X.cols(), 2).transpose();
//
//         c = m_mp->result()->basis().collocationMatrix( currentparams ) ;
//         tmp.noalias() = c * currentcoefs - m_X.transpose();
//         return 0.5 * ( (tmp * tmp.transpose()).trace() + (currentcoefs.transpose() * m_G * currentcoefs).trace() * m_mp->lambda() );
//     }
//
//     //! [OptProblemExample gradObj_into]
//     // The gradient of the objective function (resorts to finite differences if left unimplemented)
//     void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
//     {
//       // make sure that everything stays in bounds.
//       gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() );
//       u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);
//
//       // finite differences, as if is left unimplemented; use to make the check if the implementation is correct.
//       // gsOptProblem<T>::gradObj_into(u, check_result);
//
//       gsAsConstMatrix<T> currentcoefs(u.data(),  m_mp->result()->coefs().rows(), m_mp->result()->coefs().cols());
//       currentparams = gsAsConstMatrix<T>(u.data() + m_mp->result()->coefs().size(), m_X.cols(), 2).transpose();
//
//       // c_matrices = collocationMatrix1(m_mp->result()->basis(), currentparams);
//       c_matrices = m_mp->result()->basis().collocationMatrixWithDeriv( currentparams );
//       // c_matrices[0] : collocation matrix
//       // c_matrices[1] : basis partial first derivatives
//       tmp.noalias() = c_matrices[0] * currentcoefs - m_X.transpose();
//
//       result.head(currentcoefs.size()).noalias() = // d_coefs.asVector();
//       ( ( c_matrices[0].transpose()*c_matrices[0] + m_mp->lambda() * m_G ) * currentcoefs - c_matrices[0].transpose() * m_X.transpose() )
//           .reshaped(currentcoefs.size(),1);
//
//       result.middleRows(currentcoefs.size(),m_X.cols()).noalias() = (tmp * (c_matrices[1] * currentcoefs).transpose()).diagonal();
//       result.tail(m_X.cols()).noalias() = (tmp * (c_matrices[2] * currentcoefs).transpose()).diagonal();
//     }
//
//
//
// private:
//
//     gsHFitting<2, T> *m_mp;
//     const gsMatrix<T> m_params;
//     const gsMatrix<T> m_X;
//     gsSparseMatrix<T> m_G;
//     // Lastly, we forward the memebers of the base clase gsOptProblem
//     using gsOptProblem<T>::m_numDesignVars;
//
//     using gsOptProblem<T>::m_desLowerBounds;
//     using gsOptProblem<T>::m_desUpperBounds;
//
//     using gsOptProblem<T>::m_curDesign;
//
// };
// //! [OptProblem]

// This optimization problem addresses the optimization of coefficients and paramters
// for a tensor product bspline fitting model.
// Assumption: the spline space is fixed and does not change withing the optimization loop
//             the basis is encoded in the gsOptProblem member m_mp, i.e. m_mp->result()->basis()

// Input: gsFitting<T>& mp fitting object;
//        const gsMatrix<T> & uv_interiors interior parameters in [0,1]^2;
//        const gsMatrix<T> & uv_south south edge parameters in [0,1]^2;
//        const gsMatrix<T> & uv_east east edge parameters in [0,1]^2;
//        const gsMatrix<T> & uv_north north edge parameters in [0,1]^2;
//        const gsMatrix<T> & uv_west west edge parameters in [0,1]^2;
//        const gsMatrix<T> & X_interiors interior points in [0,1]^3;
//        const gsMatrix<T> & X_south south edge points in [0,1]^3;
//        const gsMatrix<T> & X_east east edge points in [0,1]^3;
//        const gsMatrix<T> & X_north north edge points in [0,1]^3;
//        const gsMatrix<T> & X_west west edge points in [0,1]^3;

// Output: gsVector<T> u containing the optimized coefficients (cx, cy, cz) in R^3,  and parameters (u, v) in [0,1]^2.
//         dimension of the tensor-product spline space: n = m_mp->result()->basis().size(),
//         number of points m = X.cols()
//         then, u = [cx_0, cx_1, ... , cx_{n-1}, cy_1, cy_2, ... , cy_{n-1}, cz_1, cz_2, ... , cz_{n-1},
//                    uinteriors_i, ..., uinteriors_f,
//                    usouth_i, ..., usouth_f,
//                    unorth_i, ..., unorth_f,
//                    veast_i, ..., veast_f,
//                    vwest_i, ..., vwest_f]

//! [OptProblemExample Class]
template <typename T>
class gsOptProblemExample : public gsOptProblem<T>
//! [OptProblemExample Class]
{
public:

    //! [OptProblemExample Constructor]
    // The constructor defines all properties of our optimization problem
    gsOptProblemExample(gsFitting<T>& mp,
                        const gsMatrix<T> & params,
                        const gsMatrix<T> & X,
                        index_t c1, // corner 0, south edge
                        index_t c2, // corner 1, east edge
                        index_t c3, // corner 2, north edge
                        index_t c4) // corner 3, west edge
    :
    m_mp(&mp),
    m_params(params),
    m_X(X),
    m_c1(c1), m_c2(c2), m_c3(c3), m_c4(c4)
    {
        // Number of design variables: how many variables we optimize, i.e. coefficiets + parametric values
        // m_numDesignVars  = m_mp[0].coefs().size() + m_params.size(); // dim * spline-dofs + numPts
        m_numDesignVars  = m_mp->result()->coefs().size() + m_params.size();

        // design bounds
        m_desLowerBounds.resize(m_numDesignVars);
        m_desUpperBounds.resize(m_numDesignVars);

        // coefficiets in R^3, namely no bounds
        for(index_t i = 0; i < m_mp->result()->coefs().size(); i++)
        {
          m_desLowerBounds[i] = -1.0e19; // lower bound on the coefficients
          m_desUpperBounds[i] =  1.0e19; // upper bound on the coefficients
        }

        // parameters in [0,1]^2
        currentparams = m_params.transpose(); // m_mp.returnParamValues().transpose();
        // gsDebugVar(currentparams.rows());
        // gsDebugVar(currentparams.cols());
        // for(index_t i = 0; i < currentparams.size(); i++)
        for(index_t i = 0; i < m_c1; i++) // u_interior parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = 0.; // lower bound on the interior parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = 1.; // upper bound on the interior parameters
        }
        // gsDebugVar(c1);
        // gsDebugVar(m_mp->result()->coefs().size());
        for(index_t i = m_c1; i < m_c2; i++) // u_south parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = 0.; // lower bound on the interior parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = 1.; // upper bound on the interior parameters
        }
        // gsDebugVar(m_mp->result()->coefs().size() + m_c1);
        for(index_t i = m_c2; i < m_c3; i++) // u_east parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = 1.; // lower bound on the interior parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = 1.; // upper bound on the interior parameters
        }
        // gsDebugVar(m_mp->result()->coefs().size() + m_c2);
        for(index_t i = m_c3; i < m_c4; i++) // u_north parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size()+ i] = 0.; // lower bound on the interior parameters
          m_desUpperBounds[m_mp->result()->coefs().size()+ i] = 1.; // upper bound on the interior parameters
        }
        // gsDebugVar(m_mp->result()->coefs().size() + m_c3);
        for(index_t i = m_c4; i < currentparams.rows(); i++) // u_west parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = 0.; // lower bound on the interior parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = 0.; // upper bound on the interior parameters
        }
        // gsDebugVar(m_mp->result()->coefs().size() + m_c4);

        index_t v_shift = m_mp->result()->coefs().size() + currentparams.rows();
        // gsDebugVar(v_shift);
        for(index_t i = 0; i < m_c1; i++) // v_interior parameters
        {
          m_desLowerBounds[v_shift + i] = 0.; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = 1.; // upper bound on the interior parameters
        }
        for(index_t i = m_c1; i < m_c2; i++) // v_south parameters
        {
          m_desLowerBounds[v_shift + i] = 0.; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = 0.; // upper bound on the interior parameters
        }
        for(index_t i = m_c2; i < m_c3; i++) // v_east parameters
        {
          m_desLowerBounds[v_shift + i] = 0.; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = 1.; // upper bound on the interior parameters
        }
        for(index_t i = m_c3; i < m_c4; i++) // u_north parameters
        {
          m_desLowerBounds[v_shift + i] = 1.; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = 1.; // upper bound on the interior parameters
        }
        for(index_t i = m_c4; i < currentparams.size()/2; i++) // u_west parameters
        {
          m_desLowerBounds[v_shift + i] = 0.; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = 1.; // upper bound on the interior parameters
        }


        // Initialization of the smoothing matrix that we need to define the objective function.
        m_G.resize(m_mp->result()->coefs().rows(), m_mp->result()->coefs().rows());
        m_G.reservePerColumn( cast<T,index_t>( (2 * m_mp->result()->basis().maxDegree() + 1) * 1.333 ) );
        m_mp->applySmoothing(m_mp->lambda(), m_G);

        // design variables: whant we do optimize.
        // c_x, c_y, c_z, u, v
        m_curDesign.resize(m_numDesignVars,1);
        m_curDesign << m_mp->result()->coefs().reshape(m_mp->result()->coefs().size(),1), currentparams.reshape(currentparams.size(),1);
    }
    //! [OptProblemExample Constructor]

public:

    //temporary storage
    mutable gsMatrix<T> tmp, currentparams;
    mutable gsSparseMatrix<T> c;
    mutable std::vector<gsSparseMatrix<T>> c_matrices;

    //! [OptProblemExample evalObj]
    // The evaluation of the objective function must be implemented
    // look at the manuscript for its rigurous definition
    // idea: 1/2 * ( (spline_model - points)^2 + lambda * smoothing_term )
    T evalObj( const gsAsConstVector<T> & u ) const
    {

        gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() ); // keep point within design bounds
        u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);

        gsAsConstMatrix<T> currentcoefs(u.data(),  m_mp->result()->coefs().rows(), m_mp->result()->coefs().cols());
        currentparams = gsAsConstMatrix<T>(u.data() + m_mp->result()->coefs().size(), m_X.cols(), 2).transpose();

        c = m_mp->result()->basis().collocationMatrix( currentparams ) ;
        tmp.noalias() = c * currentcoefs - m_X.transpose();
        return 0.5 * ( (tmp * tmp.transpose()).trace() + (currentcoefs.transpose() * m_G * currentcoefs).trace() * m_mp->lambda() );
    }

    //! [OptProblemExample gradObj_into]
    // The gradient of the objective function (resorts to finite differences if left unimplemented)

    //   compute the derivative with respect to the coefficients (cx, cy, cz)
    //   compute the partial derivative with respect to u-parameter
    //   compute the partial derivatice with respect to v-parameter
    void gradObj_into( const gsAsConstVector<T> & u, gsAsVector<T> & result) const
    {
      // make sure that everything stays in bounds.
      gsAsVector<T> u1(const_cast<T*>(u.data()), u.size() );
      u1 = u1.cwiseMax(m_desLowerBounds).cwiseMin(m_desUpperBounds);

      // finite differences, as if is left unimplemented; use to make the check if the implementation is correct.
      // gsOptProblem<T>::gradObj_into(u, check_result);

      gsAsConstMatrix<T> currentcoefs(u.data(),  m_mp->result()->coefs().rows(), m_mp->result()->coefs().cols());
      currentparams = gsAsConstMatrix<T>(u.data() + m_mp->result()->coefs().size(), m_X.cols(), 2).transpose();

      //c_matrices = collocationMatrix1(m_mp->result()->basis(), currentparams);
      c_matrices = m_mp->result()->basis().collocationMatrixWithDeriv( currentparams );
      // c_matrices[0] : collocation matrix
      // c_matrices[1] : basis partial first derivatives
      tmp.noalias() = c_matrices[0] * currentcoefs - m_X.transpose();

      result.head(currentcoefs.size()).noalias() = // d_coefs.asVector();
      ( ( c_matrices[0].transpose()*c_matrices[0] + m_mp->lambda() * m_G ) * currentcoefs - c_matrices[0].transpose() * m_X.transpose() )
          .reshaped(currentcoefs.size(),1);

      result.middleRows(currentcoefs.size(),m_X.cols()).noalias() = (tmp * (c_matrices[1] * currentcoefs).transpose()).diagonal();
      result.tail(m_X.cols()).noalias() = (tmp * (c_matrices[2] * currentcoefs).transpose()).diagonal();
    }



private:

    gsFitting<T> *m_mp;
    const gsMatrix<T> m_params;
    const gsMatrix<T> m_X;
    index_t m_c1;
    index_t m_c2;
    index_t m_c3;
    index_t m_c4;
    gsSparseMatrix<T> m_G;

    // Lastly, we forward the memebers of the base clase gsOptProblem
    using gsOptProblem<T>::m_numDesignVars;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_curDesign;

};
//! [OptProblem]


int main(int argc, char *argv[])
{

    index_t maxIter = 1; // for the optimization algorithm
    index_t refIter = 1; // for the adaptive loop
    index_t maxPcIter = 1; // for the parameter correction
    index_t numURef = 0;
    index_t mupdate = 20;
    real_t deg = 2;
    index_t numKnots = 2;
    real_t lambda = 1e-6;
    real_t gtoll = 1e-7;
    std::string fn = "../filedata/fitting/shiphull_scalePts.xml";
    index_t verbosity = 1;
    index_t extension = 2;
    real_t tolerance = 1e-04;
    bool apdm = false;

    gsCmdLine cmd("Adaptive global THB-spline surface fitting by L-BFGS: http://dx.doi.org/10.1016/j.cagd.2012.03.004");

    // maximum number of iterations.
    cmd.addInt("i", "iter", "number of maximum iterations for the optimization algorithm for CPDM.", maxIter);
    cmd.addInt("l", "level", "number of maximum iterations for the adaptive loop.", refIter);
    cmd.addInt("c", "step", "number of parameter correction steps for APDM.", maxPcIter);

    // hierarchical fitting settings
    cmd.addReal("d", "degree", "bi-degree (d,d).", deg);
    cmd.addReal("s", "smoothing", "smoothing weight", lambda);
    cmd.addInt("n", "interiors", "number of interior knots in each direction.", numKnots);
    cmd.addString("f", "filename", "name of the .xml file containing the data", fn);
    cmd.addReal("e", "tolerance", "error tolerance (desired upper bound for pointwise error)", tolerance);

    //HLBFGS options:
    cmd.addInt("m", "update", "number of LBFGS updates.", mupdate);
    cmd.addInt("b", "print", "set printing verbosity", verbosity);
    cmd.addReal("g", "gtoll", "stopping criteria on ||g||", gtoll);

    // enamble for comparison
    cmd.addSwitch("a", "apdm", "run the A-PDM algorithm.", apdm);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    real_t threshold = tolerance;
    gsStopwatch time;

    std::ofstream file_opt, file_pc;

    gsFileData<> fd_in(fn);
    gsMatrix<> uv, X;
    std::vector<index_t> corners;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, X);
    //! [Read data]

    gsWriteParaviewPoints(uv, "parameters");
    gsWriteParaviewPoints(X, "points");

    GISMO_ENSURE( uv.cols() == X.cols() && uv.rows() == 2 && X.rows() == 3,
                  "Wrong input, check id of matrices in the .xml file");

    gsInfo << "Reordering parameters and points as interiors,\n"
              "and anticlockwise boundaried, i.e. south edge, east edge, north edge, west edge.\n";
    sortPointCloud(uv,X,corners);
    index_t c1 = corners[0];
    index_t c2 = corners[1];
    index_t c3 = corners[2];
    index_t c4 = corners[3];

    gsWriteParaviewPoints(uv, "parameters");
    gsWriteParaviewPoints(X, "points");

    GISMO_ENSURE( uv.cols() == X.cols() && uv.rows() == 2 && X.rows() == 3,
                  "Wrong input");

    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, numKnots, deg+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, numKnots, deg+1 ) ;

    // Create a tensor-basis and apply initial uniform refinement
    gsTensorBSplineBasis<2> T_tbasis( u_knots, v_knots );
    T_tbasis.uniformRefine( (1<<numURef)-1 );

    // Create Initial hierarchical basis
    gsTHBSplineBasis<2>  basis ( T_tbasis ) ;


    std::vector<unsigned> ext;
    ext.push_back(extension);
    ext.push_back(extension);

    gsHFitting<2, real_t> opt_f(uv, X, basis, 0, ext, lambda); // gsHFitting object for C-PDM
    const std::vector<real_t> & errors = opt_f.pointWiseErrors();
    std::vector<real_t> errors2;
    real_t sum_of_errors2;


    file_opt.open("results_adaptive_CPDM.csv");
    file_opt << "it, time, fraction, dofs, rmse\n";
    real_t finaltime_adaptiveLoop = 0;
    for(int i = 0; i <= refIter; i++) // adaptive loop on the spline space
    {
        gsInfo<<"------------------------------------------------\n";
        gsInfo << "C-PDM for the adaptive loop.\n";
        gsInfo<<"Adaptive loop iteration "<<i<<".."<<"\n";
        real_t finaltime_itLoop = 0.;

        time.restart();
        opt_f.nextIteration(tolerance, threshold, 0); // no parameter correction here, but the combined approach.
        finaltime_itLoop += time.stop();
        // the first opt_f is empty, therefore only the fitting and the computation of the errors is applied.
        // the computation of the errors is not needed, but it is by default done.

        // after the least square fitting with THB-splines, we optimize the control points and the parametric values.
        gsMatrix<> coefs(basis.size(), 3);
        coefs = opt_f.result()->coefs();
        gsMatrix<> params(2, X.cols());
        params = opt_f.returnParamValues();

        // gsOptProblemExample<real_t> problem(opt_f, params, X); // FIRST THE fitting_object, THEN THE PARAMETERS AND THEN THE POINTS.
        gsOptProblemExample<real_t> problem(opt_f, params, X, c1, c2, c3, c4); // FIRST THE fitting object, THEN THE PARAMETERS AND THEN THE POINTS.
        // the params here are not the original ones, as in the tensor product approach,
        // but the ones computed in each adaptive loop.
        gsOptimizer<real_t> * optimizer;
        gsGeometry<>::uPtr tmp = (opt_f.result()->basis().makeGeometry(opt_f.result()->coefs()));
        gsTHBSpline<2,real_t> & original = dynamic_cast<gsTHBSpline<2,real_t> &>(*tmp);
        // gsTHBSpline<2, real_t> original(opt_f.result()->basis(), opt_f.result()->coefs());

        gsMatrix<> currentparams(2, X.cols());
        gsMatrix<> currentcoefs( original.coefs().rows(), original.coefs().cols() ); // the original for each iteration of the adaptive loop.

        optimizer = new gsHLBFGS<real_t>(&problem);
        optimizer->options().setInt("Verbose",verbosity);
        optimizer->options().setReal("MinGradLen", gtoll);
        optimizer->options().setInt("LBFGSUpdates", mupdate);

        //! [Solve]
        // Start the optimization
        gsVector<> in(original.coefs().size() + params.size());
        gsMatrix<> t_params(params.cols(), params.rows());
        t_params = params.transpose();
        in << original.coefs().reshape( original.coefs().size() ,1), t_params.reshape( t_params.size() ,1);
        gsInfo << "Optimization loop, for " << maxIter << "iterations.\n";
        optimizer->options().setInt("MaxIterations", maxIter);
        time.restart();
        optimizer->solve(in);
        finaltime_itLoop += time.stop();

        gsMatrix<> finaldesign = optimizer->currentDesign();

        for(index_t i=0; i < original.coefs().rows(); i++)
        {
          currentcoefs(i,0) = finaldesign(i,0);
          currentcoefs(i,1) = finaldesign(i + original.coefs().rows(),0);
          currentcoefs(i,2) = finaldesign(i + 2*original.coefs().rows(),0);
        }

        for(index_t i = original.coefs().size(); i < finaldesign.size(); i++)
        {
          index_t poff = X.cols() + original.coefs().size();
          if ( i < X.cols() + original.coefs().size())
          {
            currentparams(0, i - original.coefs().size() ) = finaldesign(i);
          }
          else
          {
            currentparams(1, i - poff ) = finaldesign(i);
          }

        }

        gsInfo << "max error before optimization: " << opt_f.maxPointError() << "\n";
        opt_f.updateGeometry(currentcoefs, currentparams);
        gsInfo << "max error after optimization: " << opt_f.maxPointError() << "\n";

        gsMesh<> mesh(opt_f.result()->basis());
        gsMatrix<> uv_fitting = opt_f.returnParamValues() ;
        gsWriteParaview(mesh, internal::to_string(i+1) + "_iter_mesh_cpdm");
        gsWriteParaview(*opt_f.result(), internal::to_string(i+1) + "_iter_geo_cpdm", 100000, true);
        gsWriteParaviewPoints(uv_fitting, internal::to_string(i+1) + "_iter_fitting_parameters_pdm");

        // compute mean squared error
        opt_f.get_Error(errors2, 0);
        sum_of_errors2 = std::accumulate(errors2.begin(), errors2.end(), 0.0);

        gsInfo<<"Fitting time: "<< finaltime_itLoop <<"\n";

        real_t rmse = 0.; // fitting error
        gsMatrix<> mtmp = opt_f.result()->eval(opt_f.returnParamValues()) - X;
        real_t pred_eval = (mtmp * mtmp.transpose()).trace();
        rmse += math::pow(pred_eval/X.cols(), 0.5);


        index_t dofs = opt_f.result()->basis().size();
        real_t minPointError = opt_f.minPointError();
        real_t maxPointError = opt_f.maxPointError();
        real_t mseError = sum_of_errors2/errors2.size();
        real_t percentagePoint = 100.0 * opt_f.numPointsBelow(tolerance)/errors.size();

        gsInfo<<"Fitted with "<< opt_f.result()->basis() <<"\n";
        gsInfo    << "DOFs         : "<< dofs <<"\n";
        std::cout << "Min distance : "<< minPointError << std::scientific <<"\n";
        std::cout << "Max distance : "<< maxPointError << std::scientific <<"\n";
        std::cout << "         MSE : "<< mseError << std::scientific <<"\n";
        gsInfo<<"Points below tolerance: "<< percentagePoint <<"%.\n";

        finaltime_adaptiveLoop += finaltime_itLoop;
        //file_opt << "it, time, fraction, rmse\n";
        std::setprecision(12);
        file_opt << std::to_string(i+1) << "," << std::to_string(finaltime_adaptiveLoop) << "," << std::to_string(finaltime_itLoop) << "," << std::to_string(dofs) << "," << std::to_string(rmse) << "\n";

        if ( opt_f.maxPointError() < tolerance )
        {
            gsInfo<<"Error tolerance achieved after "<< i <<" iterations.\n";
            break;
        }

    } // adaptive loop
    file_opt.close();
    gsInfo << "C-PDM total time: " << finaltime_adaptiveLoop << "\n";

    if(apdm)
    {
      finaltime_adaptiveLoop = 0;
      gsInfo << "Running the A-PDM algorithm for comparion.\n";
      gsTHBSplineBasis<2>  refbasis ( T_tbasis ) ;
      gsHFitting<2, real_t> ref(uv, X, refbasis, 0, ext, lambda); // gsHFitting object for A-PDM
      const std::vector<real_t> & adapt_errors = ref.pointWiseErrors();
      std::vector<real_t> adapt_errors2;
      real_t adapt_sum_of_errors2;

      file_pc.open("results_adaptive_APDM.csv");
      file_pc << "it, time, fraction, dofs, rmse\n";
      for(int i = 0; i <= refIter; i++) // adaptive loop on the spline space
      {
          gsInfo<<"------------------------------------------------\n";
          gsInfo << "A-PDM for the adaptive loop.\n";
          gsInfo<<"Adaptive loop iteration "<<i<<".."<<"\n";
          real_t finaltime_itLoop = 0.;

          time.restart();
          ref.nextIteration(tolerance, threshold, maxPcIter);
          finaltime_itLoop += time.stop();

          real_t rmse = 0.; // fitting error
          gsMatrix<> mtmp = ref.result()->eval(ref.returnParamValues()) - X;
          real_t pred_eval = (mtmp * mtmp.transpose()).trace();
          rmse += math::pow(pred_eval/X.cols(), 0.5);

          gsMesh<> mesh(ref.result()->basis());
          gsMatrix<> uv_fitting = ref.returnParamValues() ;
          gsWriteParaview(mesh, internal::to_string(i+1) + "_iter_mesh_apdm");
          gsWriteParaview(*ref.result(), internal::to_string(i+1) + "_iter_geo_apdm");
          gsWriteParaviewPoints(uv_fitting, internal::to_string(i+1) + "_iter_fitting_parameters_apdm");

          ref.get_Error(adapt_errors2, 0);
          adapt_sum_of_errors2 = std::accumulate(errors2.begin(), errors2.end(), 0.0);

          gsInfo<<"Fitting time: "<< finaltime_itLoop <<"\n";
          index_t dofs = ref.result()->basis().size();
          real_t minPointError = ref.minPointError();
          real_t maxPointError = ref.maxPointError();
          real_t mseError = adapt_sum_of_errors2/adapt_errors2.size();
          real_t percentagePoint = 100.0 * ref.numPointsBelow(tolerance)/adapt_errors.size();

          gsInfo<<"Fitted with "<< ref.result()->basis() <<"\n";
          gsInfo    << "DOFs         : "<< dofs <<"\n";
          std::cout << "Min distance : "<< minPointError << std::scientific <<"\n";
          std::cout << "Max distance : "<< maxPointError << std::scientific <<"\n";
          std::cout << "         MSE : "<< mseError << std::scientific <<"\n";
          gsInfo<<"Points below tolerance: "<< percentagePoint <<"%.\n";

          finaltime_adaptiveLoop += finaltime_itLoop;
          std::setprecision(12);
          file_pc << std::to_string(i+1) << "," << std::to_string(finaltime_adaptiveLoop) << "," << std::to_string(finaltime_itLoop) << "," << std::to_string(dofs) << "," << std::to_string(rmse) << "\n";

          if ( ref.maxPointError() < tolerance )
          {
              gsInfo<<"Error tolerance achieved after "<< i <<" iterations.\n";
              break;
          }

      }
      file_pc.close();
      gsInfo << "A-PDM total time: " << finaltime_adaptiveLoop << "\n";
    } // apdm

    return EXIT_SUCCESS;

}
