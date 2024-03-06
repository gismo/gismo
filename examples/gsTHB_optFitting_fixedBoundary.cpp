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

#ifdef gsParasolid_ENABLED
#include <gsParasolid/gsWriteParasolid.h>
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
  for( size_t i = 0; i < interiors.size(); i++ )
  {
    uv_interiors.col(i) = parameters.col(interiors[i]);
    p_interiors.col(i) = points.col(interiors[i]);
  }

  uv_west.resize(2, b_west.size());
  gsMatrix<T> tmp_west(3, b_west.size());
  for( size_t i = 0; i < b_west.size(); i++ )
  {
    uv_west.col(i) = parameters.col(b_west[i]);
    tmp_west.col(i) = points.col(b_west[i]);
  }

  uv_east.resize(2, b_east.size());
  gsMatrix<T> tmp_east(3, b_east.size());
  for( size_t i = 0; i < b_east.size(); i++ )
  {
    uv_east.col(i) = parameters.col(b_east[i]);
    tmp_east.col(i) = points.col(b_east[i]);
  }

  uv_south.resize(2, b_south.size());
  gsMatrix<T> tmp_south(3, b_south.size());
  for( size_t i = 0; i < b_south.size(); i++ )
  {
    uv_south.col(i) = parameters.col(b_south[i]);
    tmp_south.col(i) = points.col(b_south[i]);
  }

  uv_north.resize(2, b_north.size());
  gsMatrix<T> tmp_north(3, b_north.size());
  for( size_t i = 0; i < b_north.size(); i++ )
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
  for(size_t i = 0; i<tmp.size(); i++)
  {
    p_south.col(i) = tmp_south.col(tmp[i]);
  }
  uv_south.transposeInPlace();


  tmp.clear();
  tmp = uv_east.idxByColumn(1);
  // gsDebugVar(uv_east);
  p_east.resize(tmp_east.rows(), tmp_east.cols());
  for(size_t i = 0; i<tmp.size(); i++)
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
  for(size_t i = 0; i<tmp.size(); i++)
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
  for(size_t i = 0; i<tmp.size(); i++)
  {
    p_west.col(i) = tmp_west.col(tmp[i]);
  }
  uv_west.transposeInPlace();


  // gsWriteParaviewPoints(uv_interiors, "uv_interiors");
  // gsWriteParaviewPoints(p_interiors, "p_interiors");
  //
  // gsWriteParaviewPoints(uv_west, "uv_west");
  // gsWriteParaviewPoints(tmp_west, "p_west");
  //
  // gsWriteParaviewPoints(uv_east, "uv_east");
  // gsWriteParaviewPoints(tmp_east, "p_east");
  //
  // gsWriteParaviewPoints(uv_south, "uv_south");
  // gsWriteParaviewPoints(tmp_south, "p_south");
  //
  // gsWriteParaviewPoints(uv_north, "uv_north");
  // gsWriteParaviewPoints(tmp_north, "p_north");

  // reordering of the input point cloud (parameters and points)
  parameters.resize(uv_interiors.rows(), points.cols());
  parameters << uv_interiors.row(0), uv_south.row(0), uv_east.row(0), uv_north.row(0), uv_west.row(0),
                uv_interiors.row(1), uv_south.row(1), uv_east.row(1), uv_north.row(1), uv_west.row(1);

  points.resize(p_interiors.rows(), parameters.cols());
  points << p_interiors.row(0), p_south.row(0), p_east.row(0), p_north.row(0), p_west.row(0),
                p_interiors.row(1), p_south.row(1), p_east.row(1), p_north.row(1), p_west.row(1),
                p_interiors.row(2), p_south.row(2), p_east.row(2), p_north.row(2), p_west.row(2);

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
                        index_t c4, // corner 3, west edge
                        bool constrainCorners)
    :
    m_mp(&mp),
    m_params(params),
    m_X(X),
    m_c1(c1), m_c2(c2), m_c3(c3), m_c4(c4)
    {
        T u_min = m_params.row(0).minCoeff(),
          u_max = m_params.row(0).maxCoeff(),
          v_min = m_params.row(1).minCoeff(),
          v_max = m_params.row(1).maxCoeff();

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

        for(index_t i = 0; i < m_c1; i++) // u_interior parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = u_min; // lower bound on the interior parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = u_max; // upper bound on the interior parameters
        }

        // m_c1 = [0,0]
        m_desLowerBounds[m_mp->result()->coefs().size() + m_c1] = u_min; // lower bound on the LEFT SOUTH corner
        if(constrainCorners)
            m_desUpperBounds[m_mp->result()->coefs().size() + m_c1] = u_min; // upper bound on the LEFT SOUTH corner
        else
            m_desUpperBounds[m_mp->result()->coefs().size() + m_c1] = u_max; // upper bound on the LEFT SOUTH corner

        for(index_t i = m_c1+1; i < m_c2; i++) // u_south parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = u_min; // lower bound on the interior parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = u_max; // upper bound on the interior parameters
        }

        // m_c2 = [1,0]
        for(index_t i = m_c2; i < m_c3; i++) // u_east parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = u_max; // lower bound on the interior parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = u_max; // upper bound on the interior parameters
        }

        // m_c3 = [1,1]
        if(constrainCorners)
            m_desLowerBounds[m_mp->result()->coefs().size()+ m_c3] = u_max; // lower bound on the interior parameters
        else
            m_desLowerBounds[m_mp->result()->coefs().size()+ m_c3] = u_min; // lower bound on the interior parameters
        m_desUpperBounds[m_mp->result()->coefs().size()+ m_c3] = u_max; // upper bound on the interior parameters

        for(index_t i = m_c3+1; i < m_c4; i++) // u_north parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size()+ i] = u_min; // lower bound on the interior parameters
          m_desUpperBounds[m_mp->result()->coefs().size()+ i] = u_max; // upper bound on the interior parameters
        }

        // m_c4 = [0,1]
        for(index_t i = m_c4; i < currentparams.rows(); i++) // u_west parameters
        {
          m_desLowerBounds[m_mp->result()->coefs().size() + i] = u_min; // lower bound on the interior parameters
          m_desUpperBounds[m_mp->result()->coefs().size() + i] = u_min; // upper bound on the interior parameters
        }


        index_t v_shift = m_mp->result()->coefs().size() + currentparams.rows();
        for(index_t i = 0; i < m_c1; i++) // v_interior parameters
        {
          m_desLowerBounds[v_shift + i] = v_min; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = v_max; // upper bound on the interior parameters
        }

        // m_c1 = [0,0]
        for(index_t i = m_c1; i < m_c2; i++) // v_south parameters
        {
          m_desLowerBounds[v_shift + i] = v_min; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = v_min; // upper bound on the interior parameters
        }

        // m_c2 = [1,0]
        m_desLowerBounds[v_shift + m_c2] = v_min; // lower bound on the interior parameters
        if(constrainCorners)
            m_desUpperBounds[v_shift + m_c2] = v_min; // upper bound on the interior parameters
        else
            m_desUpperBounds[v_shift + m_c2] = v_max; // upper bound on the interior parameters

        for(index_t i = m_c2+1; i < m_c3; i++) // v_east parameters
        {
          m_desLowerBounds[v_shift + i] = v_min; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = v_max; // upper bound on the interior parameters
        }
        // m_c3 = [1,1]
        for(index_t i = m_c3; i < m_c4; i++) // u_north parameters
        {
          m_desLowerBounds[v_shift + i] = v_max; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = v_max; // upper bound on the interior parameters
        }
        // m_c4 = [0,1]

        if(constrainCorners)
            m_desLowerBounds[v_shift + m_c4] = v_max; // lower bound on the interior parameters
        else
            m_desLowerBounds[v_shift + m_c4] = v_min; // lower bound on the interior parameters
        m_desUpperBounds[v_shift + m_c4] = v_max; // upper bound on the interior parameters

        for(index_t i = m_c4+1; i < currentparams.size()/2; i++) // u_west parameters
        {
          m_desLowerBounds[v_shift + i] = v_min; // lower bound on the interior parameters
          m_desUpperBounds[v_shift + i] = v_max; // upper bound on the interior parameters
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

    // Last, we forward the members of the base clase gsOptProblem.
    using gsOptProblem<T>::m_numDesignVars;

    using gsOptProblem<T>::m_desLowerBounds;
    using gsOptProblem<T>::m_desUpperBounds;

    using gsOptProblem<T>::m_curDesign;

};
//! [OptProblem]


int main(int argc, char *argv[])
{
    bool apdm = false; // a, run a-pdm
    index_t verbosity = 1; // b, b=2 prints a lot of gsInfo
    index_t maxPcIter = 0; // c, parameter correction steps
    index_t deg = 2; // d, fitting degree
    real_t tolerance = 1e-4; // e hierarchical refinement tolerance
    std::string fn = "../filedata/fitting/shipHullPts55_scale01.xml"; // f, string file input data
    real_t gtoll = 1e-6; // g, HLBFGS stopping criteria
    index_t maxIter = 200; // i, max number of iteration for the hlbfgs optimization algorithm
    bool atdm = false; // j, run a-tdm
    // k
    index_t maxRef = 1; // l for the adaptive loop
    index_t mupdate = 20; // m, HLBFGS hessian updates
    index_t numKnots = 2; // n, fitting initial number of knots in each direction
    bool jopt = false; // o, run joint optimization
    real_t pdm_system = 0.01; // p
    index_t method = 4; // q
	  bool constrainCorners = false; // r
    real_t lambda = 1e-7; // s, fitting smoothing weight
    bool callScalePoints = false; // t
    // u, v
    index_t pc0 = 0; // w
    index_t kx = -1; // x
    index_t ky = -1; // y
    bool admissibleRef = false; // z admissible refinement method.

    // default settings, not to be changed
    real_t tdm_system = 1.;
    index_t extension = 2;
    index_t numURef = 0; // number of initial uniform refinements

    gsCmdLine cmd("Adaptive global THB-spline surface fitting by L-BFGS: http://dx.doi.org/10.1016/j.cagd.2012.03.004");


    cmd.addSwitch("a", "apdm", "run the A-PDM algorithm.", apdm); // enable for comparison
    cmd.addInt("b", "print", "set printing verbosity", verbosity);
    cmd.addInt("c", "step", "number of parameter correction steps for APDM.", maxPcIter);
    cmd.addInt("d", "degree", "bi-degree (d,d).", deg);
    cmd.addReal("e", "tolerance", "error tolerance (desired upper bound for pointwise error)", tolerance);
    cmd.addString("f", "filename", "name of the .xml file containing the data", fn);
    cmd.addReal("g", "gtoll", "stopping criteria on ||g||", gtoll);
    cmd.addInt("i", "iter", "number of maximum iterations for the optimization algorithm for CPDM.", maxIter);
    cmd.addSwitch("j", "atdm", "run the A-TDM algorithm.", atdm); // enable for comparison
    // k
    cmd.addInt("l", "level", "number of maximum iterations for the adaptive loop.", maxRef);
    cmd.addInt("m", "update", "number of LBFGS updates.", mupdate);
    cmd.addInt("n", "interiors", "number of interior knots in each direction.", numKnots);
    cmd.addSwitch("o", "jopt", "run the Joing OPTimization algorithm.", jopt); // enable for comparison
    cmd.addReal("p", "pdm", "use pdm method", pdm_system);
    cmd.addInt("q", "tdm_method", "specify tdm method", method);
    cmd.addSwitch("r", "constrainCorners", "constrain the corners", constrainCorners);
    cmd.addReal("s", "smoothing", "smoothing weight", lambda);
    cmd.addSwitch("t", "scale", "scale input data.", callScalePoints); // enable for comparison
    //u, v
    cmd.addInt("w", "pc0", "number of initial parameter correction steps.", pc0);
    cmd.addInt("x", "xknt", "number of interior knots in x-direction.", kx);
    cmd.addInt("y", "yknt", "number of interior knots in y-direction.", ky);
    cmd.addSwitch("z", "admRef", "Enable admissible refinement.", admissibleRef); // enable for admissibiliy

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    std::string method_name = "";

    gsFitting<real_t>::tdm_method method_enum;
    if(method == 0)
        {
          gsInfo << "Vanilla TDM: do not use.\n";
          return -1;
        }
    else if(method == 1)
      {
        gsInfo << "TDM with boundary PDM: do not use.\n";
        return -1;
      }
    else if(method == 2)
    {
        gsInfo << "TDM with boundary TANGENT: do not use.\n";
        return -1;
      }
    else if(method == 3)
    {
        gsInfo << "PDM: use different implementation.\n";
        return -1;
      }
    else if(method == 4)
    {
        gsInfo << "HDM with CONSTANT blending weights.\n";
        method_enum = gsFitting<real_t>::hybrid_pdm_tdm_boundary_pdm;
        method_name = "constant = " + std::to_string(pdm_system);
      }
    else if(method == 5)
    {
        gsInfo << "HDM with boundary TANGENT: do not use.\n";
        return -1;
      }
    else if(method == 6)
    {
        gsInfo << "HDM with ERROR blending weights.\n";
        method_enum = gsFitting<real_t>::hybrid_error_pdm_tdm_boundary_pdm;
        method_name = "error";
      }
    else if(method == 7)
    {
        gsInfo << "HDM with CURVATURE blending weights.\n";
        method_enum = gsFitting<real_t>::hybrid_curvature_pdm_tdm_boundary_pdm;
        method_name = "curvature";
    }
    else
    {
        gsWarn << "Unknown method, exiting." << std::endl;
        return -1;
    }


    if(kx < 0)
      kx = numKnots;
    if(ky < 0)
      ky = numKnots;

    real_t threshold = tolerance;
    gsStopwatch gsTime;
    time_t now = time(0);

    std::ofstream file_opt, file_apdm, file_atdm;

    gsFileData<> fd_in(fn);
    gsMatrix<> uv, P, X;
    std::vector<index_t> corners;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, P);
    //! [Read data]

    gsWriteParaviewPoints(uv, "uv");
    gsWriteParaviewPoints(P, "points");

    GISMO_ENSURE( uv.cols() == P.cols() && uv.rows() == 2 && P.rows() == 3,
                  "Wrong input, check id of matrices in the .xml file");

    gsInfo << "Reordering parameters and points as interiors,\n"
              "and anticlockwise boundaried, i.e. south edge, east edge, north edge, west edge.\n";

    if(callScalePoints)
        scalePoints(P,X);
    else
        X = P;

    sortPointCloud(uv,X,corners);
    index_t c1 = corners[0];
    index_t c2 = corners[1];
    index_t c3 = corners[2];
    index_t c4 = corners[3];

    gsWriteParaviewPoints(uv, "parameters");
    gsWriteParaviewPoints(X, "points_x");

    GISMO_ENSURE( uv.cols() == X.cols() && uv.rows() == 2 && X.rows() == 3,
                  "Wrong input");

    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    // // Set the default
    // if(numVKnots == -1)
    //     numVKnots = numUKnots;

    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, kx, deg+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, ky, deg+1 ) ;

    // Create a tensor-basis and apply initial uniform refinement
    gsTensorBSplineBasis<2> T_tbasis( u_knots, v_knots );
    T_tbasis.uniformRefine( (1<<numURef)-1 );

    gsFitting<real_t> initObj( uv, X, T_tbasis);
    initObj.compute(0);
    gsGeometry<> * initGeom = initObj.result();
    gsMatrix<> initParam = initObj.returnParamValues();


    // Create Initial hierarchical basis
    gsTHBSplineBasis<2>  basis ( T_tbasis ) ;
    std::vector<unsigned> ext;
    ext.push_back(extension);
    ext.push_back(extension);

    gsHFitting<2, real_t> opt_f(uv, X, basis, 0, ext, lambda); // gsHFitting object for C-PDM
    const std::vector<real_t> & errors = opt_f.pointWiseErrors();
    std::vector<real_t> errors2;

    std::string prefix = "adaptive";

    file_opt.open(std::to_string(now)+"results_adaptive_CPDM.csv");
    file_opt << "m, deg, pen, dofs, optTol, optIt, refIt, min, max, mse, rmse, perc, refTol, time\n";

    file_apdm.open(std::to_string(now)+"results_adaptive_APDM.csv");
    file_apdm << "m, deg, pen, dofs, refIt, pc, min, max, mse, rmse, perc, refTol, time\n";

    file_atdm.open(std::to_string(now)+"results_adaptive_ATDM.csv");
    file_atdm << "m, deg, pen, dofs, refIt, pc, min, max, mse, rmse, perc, refTol, weights, time\n";

    ////////////////////////////////////////////////////////////////////////////////
    //
    // @Sofia: Can we please partition the main function a bit?
    // Something like
    //
    // callJPDM(...);
    //
    // if(apdm)
    //     callAPDM(...);
    //
    // if(tdm)
    //     callTDM(...);
    //
    // Also, there is quite some copy-pasting that could be turned into functions,
    // e.g., writing the CSV files.
    //
    ////////////////////////////////////////////////////////////////////////////////

    tdm_system = 1.-pdm_system;

    real_t finaltime_adaptiveLoop = 0;
    // gtoll = gtoll * 4;
    // maxIter = maxIter / 2;
    if (jopt)
    {
    for(int refIt = 0; refIt <= maxRef; refIt++) // adaptive loop on the spline space
    {
        gsInfo<<"------------------------------------------------\n";
        gsInfo << "C-PDM for the adaptive loop.\n";
        gsInfo<<"Adaptive loop iteration "<< refIt <<".."<<"\n";
        prefix = std::to_string(now) + "adaptive" + internal::to_string(refIt) + "it_";
        real_t finaltime_itLoop = 0.;

        gsTime.restart();
        // opt_f.nextIteration(tolerance, threshold, maxPcIter); // no parameter correction here, but the combined approach.
        opt_f.nextRefinement(tolerance, threshold, 0, admissibleRef); //(tolerance, threshold, maxPcIter = 0); // no parameter correction here, but the combined approach.
        finaltime_itLoop += gsTime.stop();
        // the first opt_f is empty, therefore only the fitting and the computation of the errors is applied.
        // the computation of the errors is not needed, but it is by default done.

        gsWriteParaview(*opt_f.result(), prefix + "cpdm_geo_in", 100000, true);
        gsInfo << opt_f.result()->basis().size() << ", " << opt_f.result()->coefs().rows() << " x " << opt_f.result()->coefs().cols() << "\n";

#ifdef gsParasolid_ENABLED
        gsTHBSpline<2>* result = static_cast<gsTHBSpline<2>*>(opt_f.result());
        extensions::gsWriteParasolid<real_t>(*result, prefix + "cpdm_geo_in");
#endif


        // after the least square fitting with THB-splines, we optimize the control points and the parametric values.
        gsMatrix<> coefs(opt_f.result()->basis().size(), 3);
        coefs = opt_f.result()->coefs();
        gsMatrix<> params(2, X.cols());
        params = opt_f.returnParamValues();


        gsInfo << "MAKE GEOMETRY.\n";
        gsInfo << "basis = \n" << opt_f.result()->basis().size() << "\n";
        gsInfo << "coefs = \n" << opt_f.result()->coefs().rows() << " x " <<  opt_f.result()->coefs().cols() << "\n";
        gsGeometry<>::uPtr tmp = (opt_f.result()->basis().makeGeometry(opt_f.result()->coefs()));
        gsTHBSpline<2,real_t> & original = dynamic_cast<gsTHBSpline<2,real_t> &>(*tmp);
        // gsTHBSpline<2, real_t> original(opt_f.result()->basis(), opt_f.result()->coefs());

        gsInfo << "Input geometry = \n"<< original << "\n";
        gsInfo << "max error before optimization: " << opt_f.maxPointError() << "\n";

        gsMatrix<> currentparams(2, X.cols());
        gsMatrix<> currentcoefs( original.coefs().rows(), original.coefs().cols() ); // the original for each iteration of the adaptive loop.


        // gtoll = gtoll / 4;
        // maxIter = maxIter * 2;

        gsInfo << "Initialization of the optimization problem:\n";

        // gsOptProblemExample<real_t> problem(opt_f, params, X); // FIRST THE fitting_object, THEN THE PARAMETERS AND THEN THE POINTS.
        gsInfo << opt_f.result()->coefs().rows() << " * " << opt_f.result()->coefs().cols() << " + "
               << opt_f.returnParamValues().rows() << " * " << opt_f.returnParamValues().cols() << " = "
               << opt_f.result()->coefs().rows() * opt_f.result()->coefs().cols() + opt_f.returnParamValues().rows() * opt_f.returnParamValues().cols() << "\n";

        gsOptProblemExample<real_t> problem(opt_f, params, X, c1, c2, c3, c4, constrainCorners); // FIRST THE fitting object, THEN THE PARAMETERS AND THEN THE POINTS.
        // the params here are not the original ones, as in the tensor product approach,
        // but the ones computed in each adaptive loop.
        gsOptimizer<real_t> * optimizer;
        optimizer = new gsHLBFGS<real_t>(&problem);
        optimizer->options().setInt("Verbose",verbosity);
        optimizer->options().setReal("MinGradLen", gtoll); // 1e-6 : more or less as refinement tolerance; there should be a balance between the two;
        optimizer->options().setInt("LBFGSUpdates", mupdate);
        optimizer->options().setReal("MinStepLen", 1e-12);
        optimizer->options().setInt("MaxEval", 100);


        //! [Solve]
        // Start the optimization
        gsVector<> in(original.coefs().size() + params.size());
        gsMatrix<> t_params(params.cols(), params.rows());
        t_params = params.transpose();
        in << original.coefs().reshape( original.coefs().size() ,1), t_params.reshape( t_params.size() ,1);
        optimizer->options().setInt("MaxIterations", maxIter);
        gsTime.restart();

        gsInfo << "optimizer = " << optimizer << "\n";
        gsInfo << "input = " << in.rows() << " x " << in.cols() << "\n";
        optimizer->solve(in);
        finaltime_itLoop += gsTime.stop();
        gsInfo << "Optimization loop, for " << optimizer->iterations() << " iterations.\n";

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

        opt_f.updateGeometry(currentcoefs, currentparams); // inside we comput the erros again.

        // gsDebugVar(currentparams);

        gsWriteParaviewPoints(currentparams, prefix + "_adaptive_params");
        gsInfo << "max error after optimization: " << opt_f.maxPointError() << "\n";

        gsMesh<> mesh(opt_f.result()->basis());
        gsMatrix<> uv_fitting = opt_f.returnParamValues() ;
        gsWriteParaview(mesh, prefix + "cpdm_mesh");
        gsWriteParaview(*opt_f.result(), prefix + "cpdm_geo_moved", 100000, true);

        gsMatrix<> ptwise = opt_f.pointWiseErrors(opt_f.returnParamValues(), X);
        gsMatrix<> point_colors(4, X.cols());
        point_colors << X.row(0), X.row(1), X.row(2), ptwise.row(0);

        gsWriteParaviewPoints(point_colors, prefix + "cpdm_colors");

#ifdef gsParasolid_ENABLED
        result = static_cast<gsTHBSpline<2>*>(opt_f.result());
        extensions::gsWriteParasolid<real_t>(*result, prefix + "cpdm_geo_moved");
#endif
        gsWriteParaviewPoints(uv_fitting, prefix + "cpdm_parameters");

        // compute mean squared error
        opt_f.get_Error(errors2, 0);

        gsInfo<<"Fitting time: "<< finaltime_itLoop <<"\n";

        real_t rmse = 0.; // fitting error
        gsMatrix<> mtmp = opt_f.result()->eval(opt_f.returnParamValues()) - X;
        real_t pred_eval = (mtmp * mtmp.transpose()).trace();
        rmse += math::pow(pred_eval/X.cols(), 0.5);


        index_t dofs = opt_f.result()->basis().size();
        std::vector<real_t> sol_min_max_mse = opt_f.result()->MinMaxMseErrors(currentparams,X);
        real_t minPointError = sol_min_max_mse[0];//opt_f.minPointError();
        real_t maxPointError = sol_min_max_mse[1];//opt_f.minPointError();
        real_t mseError = sol_min_max_mse[2];//opt_f.minPointError();
        real_t percentagePoint = 100.0 * opt_f.numPointsBelow(tolerance)/errors.size();

        gsInfo<<"Fitted with "<< opt_f.result()->basis() <<"\n";
        gsInfo    << "DOFs         : "<< dofs <<"\n";
        std::cout << "Min distance : "<< minPointError << std::scientific <<"\n";
        std::cout << "Max distance : "<< maxPointError << std::scientific <<"\n";
        std::cout << "         MSE : "<< mseError << std::scientific <<"\n";
        std::cout << "      myRMSE : "<< math::sqrt(mseError) << std::scientific <<"\n";
        std::cout << "     refRMSE : "<< rmse << std::scientific <<"\n";
        gsInfo<<"Points below tolerance: "<< percentagePoint <<"%.\n";

        finaltime_adaptiveLoop += finaltime_itLoop;
        file_opt << X.cols() << "," << deg << "," << lambda << "," << basis.size()<< ","
                 << gtoll << ", " << optimizer->iterations() << "," << refIt << ","
                    << sol_min_max_mse[0] << std::scientific << ","
    					      << sol_min_max_mse[1] << std::scientific << ","
    					      << sol_min_max_mse[2] << std::scientific << ","
                    << math::sqrt(sol_min_max_mse[2]) << std::scientific << ","
                    << percentagePoint << "," << tolerance << ","
                    << finaltime_itLoop << "\n";



        if ( opt_f.maxPointError() < tolerance )
        {
            gsInfo<<"Error tolerance achieved after "<< refIt <<" iterations.\n";
            break;
        }

    } // adaptive loop
    file_opt.close();
    gsInfo << "C-PDM total time: " << finaltime_adaptiveLoop << "\n";
  } // j-pdm

  if(apdm)
  {
      index_t stepPC = maxPcIter;
      //for(index_t stepPC = 1; stepPC < maxPcIter; stepPC ++)
      {
        finaltime_adaptiveLoop = 0;
        gsInfo << "Running the A-PDM algorithm for comparion.\n";
        gsTHBSplineBasis<2>  refbasis ( T_tbasis ) ;
        gsHFitting<2, real_t> ref(initParam, X, refbasis, 0, ext, lambda); // gsHFitting object for A-PDM
        ref.compute(lambda);
        const std::vector<real_t> & adapt_errors = ref.pointWiseErrors();
        std::vector<real_t> adapt_errors2;

        prefix = "adaptive";
        for(int refIt = 0; refIt <= maxRef; refIt++) // adaptive loop on the spline space
        {
          gsInfo<<"------------------------------------------------\n";
          gsInfo << "A-PDM for the adaptive loop.\n";
          gsInfo<<"Adaptive loop iteration "<< refIt <<".."<<"\n";
          prefix = std::to_string(now) + "adaptive" + internal::to_string(refIt) + "it_";
          real_t finaltime_itLoop = 0.;

          gsTime.restart();
          ref.nextIteration_pdm(tolerance, threshold, stepPC, corners, admissibleRef);
          finaltime_itLoop += gsTime.stop();

          real_t rmse = 0.; // fitting error
          gsMatrix<> mtmp = ref.result()->eval(ref.returnParamValues()) - X;
          real_t pred_eval = (mtmp * mtmp.transpose()).trace();
          rmse += math::pow(pred_eval/X.cols(), 0.5);

          gsMesh<> mesh(ref.result()->basis());
          gsMatrix<> uv_fitting = ref.returnParamValues() ;
          gsWriteParaview(mesh, prefix + "apdm_mesh");
          gsWriteParaview(*ref.result(), prefix + "apdm_geo", 100000, true);
          gsWriteParaviewPoints(uv_fitting, prefix + "apdm_parameters");

          ref.get_Error(adapt_errors2, 0);

          gsInfo<<"Fitting time: "<< finaltime_itLoop <<"\n";

          std::vector<real_t> sol_min_max_mse = ref.result()->MinMaxMseErrors(ref.returnParamValues(), X);
          index_t dofs = ref.result()->basis().size();
          real_t percentagePoint = 100.0 * ref.numPointsBelow(tolerance)/adapt_errors.size();

          gsInfo<<"Fitted with "<< ref.result()->basis() <<"\n";
          gsInfo    << "DOFs         : "<< dofs <<"\n";
          std::cout << "Min distance : "<< sol_min_max_mse[0] << std::scientific <<"\n";
          std::cout << "Max distance : "<< sol_min_max_mse[1] << std::scientific <<"\n";
          std::cout << "         MSE : "<< sol_min_max_mse[2] << std::scientific <<"\n";
          std::cout << "        RMSE : "<< math::sqrt(sol_min_max_mse[2]) << std::scientific <<"\n";
          gsInfo<<"Points below tolerance: "<< percentagePoint <<"%.\n";

          finaltime_adaptiveLoop += finaltime_itLoop;
          file_apdm << X.cols() << "," << deg << "," << lambda << ","
                  << dofs << ","<< refIt << "," << stepPC << ","
                  << sol_min_max_mse[0] << std::scientific << ","
                  << sol_min_max_mse[1] << std::scientific << ","
                  << sol_min_max_mse[2] << std::scientific << ","
                  << math::sqrt(sol_min_max_mse[2]) << std::scientific << ","
                  << percentagePoint << "," << tolerance << ","
                  << finaltime_itLoop << "\n";

          gsMatrix<> ptwise = ref.pointWiseErrors(ref.returnParamValues(), X);
          gsMatrix<> point_colors(4, X.cols());
          point_colors << X.row(0), X.row(1), X.row(2), ptwise.row(0);
          gsWriteParaviewPoints(point_colors, prefix + "APDM_colors");

          if ( ref.maxPointError() < tolerance )
          {
              gsInfo<<"Error tolerance achieved after "<< refIt <<" iterations.\n";
              break;
          }

      }
    } // loop on PC
      file_apdm.close();
    } // apdm

    if(atdm)
    {
      index_t stepPC = maxPcIter;
      //for(index_t stepPC = 3; stepPC < 4; stepPC ++)
        {

          finaltime_adaptiveLoop = 0;
          gsInfo << "Running the A-TDM algorithm for comparion.\n";
          gsTHBSplineBasis<2>  refbasis ( T_tbasis ) ;
          gsHFitting<2, real_t> ref(initParam, X, refbasis, 0, ext, lambda);
          ref.compute(lambda);
          // ref.initializeGeometry(initGeom->coefs(), ref.returnParamValues());
          // ref.parameterProjectionSepBoundary(1e-8, corners);
          const std::vector<real_t> & adapt_errors = ref.pointWiseErrors();
          std::vector<real_t> adapt_errors2;

          prefix = "adaptive";
          // ref.initializeGeometry(initGeom->coefs(), ref.returnParamValues());
          // ref.compute_tdm(lambda, 0., 1., corners);
          for(int refIt = 0; refIt <= maxRef; refIt++) // adaptive loop on the spline space
            {
              gsInfo<<"------------------------------------------------\n";
              gsInfo << "A-TDM for the adaptive loop.\n";
              gsInfo<<"Adaptive loop iteration "<< refIt <<".."<<"\n";
              prefix = std::to_string(now) + "adaptive" + internal::to_string(refIt) + "it_";
              real_t finaltime_itLoop = 0.;

              gsTime.restart();
              ref.nextIteration_tdm(tolerance, threshold, stepPC, pdm_system, tdm_system, corners, method_enum, admissibleRef);
              //(T tolerance, T err_threshold, index_t maxPcIter, T mu, T sigma, const std::vector<index_t> & interpIdx);

              finaltime_itLoop += gsTime.stop();

              real_t rmse = 0.; // fitting error
              gsMatrix<> mtmp = ref.result()->eval(ref.returnParamValues()) - X;
              real_t pred_eval = (mtmp * mtmp.transpose()).trace();
              rmse += math::pow(pred_eval/X.cols(), 0.5);

              gsMesh<> mesh(ref.result()->basis());
              gsMatrix<> uv_fitting = ref.returnParamValues() ;
              gsWriteParaview(mesh, prefix + "aTDM_mesh");
              gsWriteParaview(*ref.result(), prefix + "aTDM_geo", 100000, true);
              gsWriteParaviewPoints(uv_fitting, prefix + "aTDM_parameters");

              ref.get_Error(adapt_errors2, 0);

              gsInfo<<"Fitting time: "<< finaltime_itLoop <<"\n";

              std::vector<real_t> sol_min_max_mse = ref.result()->MinMaxMseErrors(ref.returnParamValues(), X);
              index_t dofs = ref.result()->basis().size();
              real_t percentagePoint = 100.0 * ref.numPointsBelow(tolerance)/adapt_errors.size();

              gsInfo<<"Fitted with "<< ref.result()->basis() <<"\n";
              gsInfo    << "DOFs         : "<< dofs <<"\n";
              std::cout << "Min distance : "<< sol_min_max_mse[0] << std::scientific <<"\n";
              std::cout << "Max distance : "<< sol_min_max_mse[1] << std::scientific <<"\n";
              std::cout << "         MSE : "<< sol_min_max_mse[2] << std::scientific <<"\n";
              std::cout << "        RMSE : "<< math::sqrt(sol_min_max_mse[2]) << std::scientific <<"\n";
              gsInfo<<"Points below tolerance: "<< percentagePoint <<"%.\n";

              finaltime_adaptiveLoop += finaltime_itLoop;
              file_atdm << X.cols() << "," << deg << "," << lambda << ","
                      << dofs << ","<< refIt << "," << stepPC << ","
                      << sol_min_max_mse[0] << std::scientific << ","
          					  << sol_min_max_mse[1] << std::scientific << ","
          					  << sol_min_max_mse[2] << std::scientific << ","
                      << math::sqrt(sol_min_max_mse[2]) << std::scientific << ","
                      << percentagePoint << "," << tolerance << ","
                      << method_name << ","
                      << finaltime_itLoop << "\n";

              gsMatrix<> ptwise = ref.pointWiseErrors(ref.returnParamValues(), X);
              gsMatrix<> point_colors(4, X.cols());
              point_colors << X.row(0), X.row(1), X.row(2), ptwise.row(0);
              gsWriteParaviewPoints(point_colors, prefix + "ATDM_colors");

              if ( ref.maxPointError() < tolerance )
              {
                  gsInfo<<"Error tolerance achieved after "<< refIt <<" iterations.\n";
                  break;
              }

          }
        } // loop on PC
          file_atdm.close();
    } // atdm


    return EXIT_SUCCESS;

}
