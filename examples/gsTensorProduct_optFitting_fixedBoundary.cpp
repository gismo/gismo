/** @file gsTensorProduct_optFitting_fixedBoundary.cpp

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

/*
input  : a point cloud in R^N
output : corresponding scaled point cloud in [0,1]^N
*/
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

  gsInfo << "+++++++++++++++++++ POINTS +++++++++++++++++++\n";
  gsInfo << "Interior points: " << p_interiors.rows() << " x " << p_interiors.cols() << "\n";
  gsInfo << "south points: " << p_south.rows() << " x " << p_south.cols() << "\n";
  gsInfo << "east points: " << p_east.rows() << " x " << p_east.cols() << "\n";
  gsInfo << "north points: " << p_north.rows() << " x " << p_north.cols() << "\n";
  gsInfo << "west points: " << p_west.rows() << " x " << p_west.cols() << "\n";

  gsInfo << "+++++++++++++++++++ PARAMS +++++++++++++++++++\n";
  gsInfo << "Interior params: " << uv_interiors.rows() << " x " << uv_interiors.cols() << "\n";
  gsInfo << "south params: " << uv_south.rows() << " x " << uv_south.cols() << "\n";
  gsInfo << "east params: " << uv_east.rows() << " x " << uv_east.cols() << "\n";
  gsInfo << "north params: " << uv_north.rows() << " x " << uv_north.cols() << "\n";
  gsInfo << "west params: " << uv_west.rows() << " x " << uv_west.cols() << "\n";

} // end sortPointCloud



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

        //gsDebugVar(m_c1); //[0,0]
        m_desLowerBounds[m_mp->result()->coefs().size() + m_c1] = 0.; // lower bound on the LEFT SOUTH corner
        m_desUpperBounds[m_mp->result()->coefs().size() + m_c1] = 0.; // upper bound on the LEFT SOUTH corner
        for(index_t i = m_c1+1; i < m_c2; i++) // u_south parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = 0.; // lower bound on the south parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = 1.; // upper bound on the south parameters
        }
        //gsDebugVar(m_mp->result()->coefs().size() + m_c1);

        //gsDebugVar(m_c2); //[1,0]
        for(index_t i = m_c2; i < m_c3; i++) // u_east parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = 1.; // lower bound on the east parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = 1.; // upper bound on the east parameters
        }
        //gsDebugVar(m_mp->result()->coefs().size() + m_c2);

        //gsDebugVar(m_c3); // [1,1]
        m_desLowerBounds[m_mp->result()->coefs().size()+ m_c3] = 1.; // lower bound on the RIGHT NORTH corner
        m_desUpperBounds[m_mp->result()->coefs().size()+ m_c3] = 1.; // upper bound on the RIGHT NORTH corner
        for(index_t i = m_c3+1; i < m_c4; i++) // u_north parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size()+ i] = 0.; // lower bound on the north parameters
          m_desUpperBounds[m_mp->result()->coefs().size()+ i] = 1.; // upper bound on the north parameters
        }
        //gsDebugVar(m_mp->result()->coefs().size() + m_c3);

        //gsDebugVar(m_c4); //[0,1]
        for(index_t i = m_c4; i < currentparams.rows(); i++) // u_west parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = 0.; // lower bound on the west parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = 0.; // upper bound on the west parameters
        }
        // gsDebugVar(m_mp->result()->coefs().size() + m_c4);

        index_t v_shift = m_mp->result()->coefs().size() + currentparams.rows();
        //gsDebugVar(v_shift);
        for(index_t i = 0; i < m_c1; i++) // v_interior parameters
        {
          m_desLowerBounds[v_shift + i] = 0.; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = 1.; // upper bound on the interior parameters
        }

        // m_c1 = [0,0]
        for(index_t i = m_c1; i < m_c2; i++) // v_south parameters
        {
          m_desLowerBounds[v_shift + i] = 0.; // lower bound on the south parameters
          m_desUpperBounds[v_shift + i] = 0.; // upper bound on the south parameters
        }

        // m_c2 = [1,0] -> set v to 0!
        m_desLowerBounds[v_shift + m_c2] = 0.; // lower bound on the east parameters
        m_desUpperBounds[v_shift + m_c2] = 0.; // upper bound on the east parameters
        for(index_t i = m_c2+1; i < m_c3; i++) // v_east parameters
        {
          m_desLowerBounds[v_shift + i] = 0.; // lower bound on the east parameters
          m_desUpperBounds[v_shift + i] = 1.; // upper bound on the east parameters
        }

        // m_c3 = [1,1]
        for(index_t i = m_c3; i < m_c4; i++) // v_north parameters
        {
          m_desLowerBounds[v_shift + i] = 1.; // lower bound on the north parameters
          m_desUpperBounds[v_shift + i] = 1.; // upper bound on the north parameters
        }

        // m_c4 = [0,1]
        m_desLowerBounds[v_shift + m_c4] = 1.; // lower bound on the west parameters
        m_desUpperBounds[v_shift + m_c4] = 1.; // upper bound on the west parameters
        for(index_t i = m_c4+1; i < currentparams.size()/2; i++) // v_west parameters
        {
          m_desLowerBounds[v_shift + i] = 0.; // lower bound on the west parameters
          m_desUpperBounds[v_shift + i] = 1.; // upper bound on the west parameters
        }

        //gsDebugVar(m_desLowerBounds);
        //gsDebugVar(m_desUpperBounds);


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

    bool apdm = false; // a
    index_t verbosity = 0; // b 0 (no videoprint), 1 (some videoprint), 2 (a lot of videoprint)
    real_t funcTol = 1e-4; // c
    real_t deg = 2; // d
    // e
    std::string fn = "../filedata/fitting/shiphull_v200_scalePts.xml"; // f
    real_t gtoll = 1e-7; // g, to decrease to push trough the iterations
    // h
    index_t maxIter = 1; // i
    int maxEval = 100; // j
    index_t plotIt = maxIter; // k
    // l
    index_t mupdate = 20; // m
    index_t numKnots = 2; // n
    // o
    bool ptype = false; // p, keep it false.
    // q, r
    real_t lambda = 1e-6; //s
    // t, u, v, w, x, y,
    bool plotInParaview = false; // z

    gsCmdLine cmd("Tensor product B-spline surface fitting by L-BFGS: http://dx.doi.org/10.1016/j.cagd.2012.03.004");

    cmd.addSwitch("a", "apdm", "run the A-PDM algorithm.", apdm);
    cmd.addInt("b", "print", "set printing verbosity", verbosity);
    cmd.addReal("c", "funcTol", "function tolerance used in line-search", funcTol);
    cmd.addReal("d", "degree", "bi-degree (d,d).", deg);
    // e
    cmd.addString("f", "filename", "name of the .xml file containing the data", fn);
    cmd.addReal("g", "gtoll", "stopping criteria on ||g||", gtoll);
    // h
    cmd.addInt("i", "iter", "number of maximum iterations of the optimization algorithm(s).", maxIter);
    cmd.addInt("j", "maxEval", "the max number of evaluation in line-search", maxEval);
    cmd.addInt("k", "kplot", "iteration of the optimization procedure to be plotted", plotIt);
    // l
    cmd.addInt("m", "update", "number of LBFGS updates.", mupdate);
    cmd.addInt("n", "interiors", "number of interior knots in each direction.", numKnots);
    // o
    cmd.addSwitch("p", "parameters", "input parameters: (0) from .xml file; (1) for foot-point projection;", ptype);
    // q, r
    cmd.addReal("s", "smoothing", "smoothing weight", lambda);
    // t, u, v, w, x, y,
    cmd.addSwitch("z", "plot", "(0): no paraview plot generated.", plotInParaview);


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsStopwatch gsTime;

    time_t now = time(0);

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


    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, numKnots, deg+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, numKnots, deg+1 ) ;

    // Create a tensor-basis and apply initial uniform refinement
    gsTensorBSplineBasis<2> basis( u_knots, v_knots );

    gsFitting<real_t> fitting_object(uv, X, basis);
    fitting_object.compute(lambda);

    gsMatrix<> coefs(basis.size(), 3);
    coefs = fitting_object.result()->coefs();

    gsTensorBSpline<2, real_t> original(basis, coefs);


    // step 2. Find the foot point P(t_k) on P(t) for evry data point X_k.
    // already as initial parameters, the foot point projection are provided.
    // REMARK: we can start with any different parameterization, not necessary with the foot-point projection.
    gsMatrix<> params(2, X.cols());
    for (index_t i = 0; i < X.cols(); ++i)
    {
        if (ptype)
        {
          gsVector<> newParam;
          original.closestPointTo(X.col(i).transpose(), newParam, 1e-10, true);
          params.col(i) = newParam;
        }
        else
        {
          params.col(i) = uv.col(i);
        }
    }

    gsInfo << "Fixed spline space:\n" << original.basis() << "\n";
    gsInfo << "Initial coefficients:\n" << original.coefs().rows() << " x " << original.coefs().cols() << "\n";
    gsInfo << "Initial parameters:\n" << params.rows() << " x " << params.cols() << "\n";
    gsInfo << "x0 size: " << original.coefs().size() + params.size()  << "\n";

    gsInfo << "Starting geometry:\n " << original << "\n";

    if (plotInParaview)
    {
      gsWriteParaview( original, "geo_it0", 1000, true, false);
      gsWriteParaview( original, "cnet_it0", 1000, false, true);
      gsMatrix<> originalCoefsToPlot(original.coefs().cols(), original.coefs().rows());
      originalCoefsToPlot = original.coefs().transpose();
      gsWriteParaviewPoints( originalCoefsToPlot, "coefs_it0");
      gsWriteParaviewPoints( params, "params_it0");
    }

    real_t rmse0 = 0.;
    // compute the fitting error.
    gsMatrix<> tmp0 = original.eval(params) - X;
    real_t pred0_eval = (tmp0 * tmp0.transpose()).trace();
    rmse0 += math::pow(pred0_eval/X.cols(), 0.5);

    std::cout << "Initial fitting error: " << rmse0 << std::scientific << "\n";

    // Initialization of the optimization problem.
    gsOptProblemExample<real_t> problem(fitting_object, params, X, c1, c2, c3, c4); // FIRST THE fitting object, THEN THE PARAMETERS AND THEN THE POINTS.
    gsOptimizer<real_t> * optimizer;

    // to store the information in .csv files.
    std::ofstream file_opt, file_apdm;
    file_opt.open(std::to_string(now)+"_CPDM_results.csv");
    file_opt << "m, deg, pen, dofs, it_opt, min, max, mse, rmse, time\n";


    gsInfo << "Fast fitting with HLBFGS:\n";
    // TODO: store time and fitting error in a proper way.
    // Uncomment the following line to avoid the foor loop on the maximum number of iterations.
    gsMatrix<> currentparams(2, X.cols());
    gsMatrix<> currentcoefs( original.coefs().rows(), original.coefs().cols() );
    index_t it_opt = maxIter;
    index_t maxIterComparison = maxIter;
    std::string prefix;
    for(index_t it_opt = 1; it_opt <= maxIter; it_opt++)
    {
      optimizer = new gsHLBFGS<real_t>(&problem);
      optimizer->options().setInt("Verbose",verbosity);
      gsDebugVar(gtoll);
      optimizer->options().setReal("MinGradLen", gtoll);
      optimizer->options().setInt("LBFGSUpdates", mupdate);

      optimizer->options().setReal("FuncTol", funcTol);
      optimizer->options().setInt("MaxEval", maxEval);
      optimizer->options().setReal("MinStepLen", 1e-12);


      //! [Solve]
      // Start the optimization
      gsVector<> in(original.coefs().size() + params.size());
      gsMatrix<> t_params(params.cols(), params.rows());
      t_params = params.transpose();
      in << original.coefs().reshape( original.coefs().size() ,1), t_params.reshape( t_params.size() ,1);
      optimizer->options().setInt("MaxIterations",it_opt); // set maximum number of iterations
      gsTime.restart(); // start optimization algorithm
      optimizer->solve(in);
      real_t finaltime = gsTime.stop(); // end optimization algorithm
      prefix = "uniform"+internal::to_string(optimizer->iterations())+"optIt_";

      // assemble the new geometry with optimized coefficiets and parameters
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

      if(plotInParaview && it_opt == plotIt) // plot certain output
      {
        gsWriteParaviewPoints( currentparams, prefix + "params_CPDM_it");
        gsMatrix<> currentCoefsToPlot(currentcoefs.cols(), currentcoefs.rows());
        currentCoefsToPlot = currentcoefs.transpose();
        gsWriteParaviewPoints( currentCoefsToPlot, prefix + "coefs_CPDM_it");

        // gsTensorBSpline<2, real_t> final( basis, give(currentcoefs));
        gsTensorBSpline<2, real_t> final( basis, currentcoefs);
        gsWriteParaview( final, prefix + "geo_CPDM_it", 1000, true, false);
        gsWriteParaview( final, prefix + "cnet_CPDM_it", 1000, false, true);

        gsMatrix<> pcolors(4, X.cols());
        pcolors << X.row(0), X.row(1), X.row(2), final.pointWiseErrors(currentparams,X);

        gsWriteParaviewPoints(pcolors, prefix + "colors_it");
      }

      gsTensorBSpline<2, real_t> currentGeo(basis, currentcoefs);
      std::vector<real_t> sol_min_max_mse = currentGeo.MinMaxMseErrors(params,X);

      gsInfo << it_opt << ", " << gsTime << ", " << math::sqrt(sol_min_max_mse[2]) <<"\n";
      if (verbosity > 1)
      {
        gsInfo << "\nNumber of iterations : " << optimizer->iterations() <<"\n";
        gsInfo << "Final objective value: " << optimizer->objective() <<"\n";
        gsInfo<<"Fitting time: "<< gsTime <<"\n";
        // gsInfo << "Final design:\n" << optimizer->currentDesign() <<"\n"; // this plot the whole vector given output from the optimizer.
        gsInfo << "params are moved from the originals by: " << (currentparams - params).norm() << "\n";
        gsInfo << "coefficients are moved from the originals by: " << (currentcoefs - original.coefs()).norm() << "\n";
      }

      //file_opt <<std::setprecision(3)<< std::to_string(it_opt)<<std::setprecision(12) << "," << std::to_string(finaltime) << "," << std::to_string(rmse) << "\n";

      //file_opt << "m, deg, pen, dofs, it_opt, min, max, mse, rmse, time\n";
      file_opt << X.cols() << "," << deg << "," << lambda << ","
                  << basis.size() << ","<< optimizer->iterations() << ","
                  << sol_min_max_mse[0] << std::scientific << ","
  					      << sol_min_max_mse[1] << std::scientific << ","
  					      << sol_min_max_mse[2] << std::scientific << ","
                  << math::sqrt(sol_min_max_mse[2]) << std::scientific << ","
                  << finaltime << "\n";

      gsWriteParaviewPoints( currentparams, prefix + "params_CPDM_it");
      gsMatrix<> currentCoefsToPlot(currentcoefs.cols(), currentcoefs.rows());
      currentCoefsToPlot = currentcoefs.transpose();
      gsWriteParaviewPoints( currentCoefsToPlot, prefix + "coefs_CPDM_it");

      // gsTensorBSpline<2, real_t> final( basis, give(currentcoefs));
      gsTensorBSpline<2, real_t> final( basis, currentcoefs);
      gsWriteParaview( final, prefix + "geo_CPDM_it", 1000, true, false);
      gsWriteParaview( final, prefix + "cnet_CPDM_it", 1000, false, true);

      maxIterComparison = optimizer->iterations();
      gsInfo << "Optimizer iterations performed = " << maxIterComparison << "\n";
    } // maxIter
    file_opt.close();

    if(apdm)
    {
      gsInfo << "Running the A-PDM algorithm for comparion.\n";

      file_apdm.open(std::to_string(now)+"_APDM_results.csv");
      file_apdm << "m, deg, pen, dofs, it_opt, pc, min, max, mse, rmse, time\n";

      gsFitting<real_t> ref(uv, X, basis); // original geometry, same starting point for C-PDM;
      ref.compute(lambda);
      if(verbosity > 1)
        gsInfo << "it     time     rmse\n";
      //index_t step = maxIterComparison; // uncomment to avoid foor loop on maximum number of iterations.
      for (index_t step=1; step <= maxIterComparison; step ++)
      {
        prefix = "uniform" + internal::to_string(step) + "pc_";
        gsTime.restart(); // start optimization procedure: 1 step = points projection + refit to update the control points.
        ref.parameterCorrection(1e-7, step, 1e-4); //closestPoint accuracy, orthogonality tolerance
        real_t finaltime = gsTime.stop(); // end of the optimization algorithm

      //if(plotInParaview && step == plotIt) // plot certain output
      {
        gsWriteParaviewPoints( ref.returnParamValues(), prefix + "params_APDM_it");
        gsMatrix<> currentCoefsToPlot(ref.result()->coefs().cols(), ref.result()->coefs().rows());
        currentCoefsToPlot = ref.result()->coefs().transpose();
        gsWriteParaviewPoints( currentCoefsToPlot, prefix + "coefs_APDM_it");

        gsWriteParaview( *ref.result(), prefix + "geo_APDM_it", 1000, true, false);
        gsWriteParaview( *ref.result(), prefix + "cnet_APDM_it", 1000, false, true);
      }

      // fitting error
      std::vector<real_t> sol_min_max_mse = ref.result()->MinMaxMseErrors(ref.returnParamValues(), X);

      // store data in .csv file
      file_apdm << X.cols() << "," << deg << "," << lambda << ","
              << basis.size() << ","<< step << ","
              << sol_min_max_mse[0] << std::scientific << ","
  					  << sol_min_max_mse[1] << std::scientific << ","
  					  << sol_min_max_mse[2] << std::scientific << ","
              << math::sqrt(sol_min_max_mse[2]) << std::scientific << ","
              << finaltime << "\n";

      } // maxIter
      file_apdm.close();
    } // fi apdm

    return EXIT_SUCCESS;

}
