/** @file pdm_tdm_surface.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <iostream>
#include <gismo.h>
#include <ctime>

using namespace gismo;

void writeToCSVfile(std::string name, gsMatrix<> matrix)
{
  std::ofstream file(name.c_str());
  for(int  i = 0; i < matrix.rows(); i++)
  {
    for(int j = 0; j < matrix.cols(); j++)
    {
       std::string str = std::to_string(matrix(i,j));
       if(j+1 == matrix.cols())
       {
           file<<std::setprecision(10)<<str;
       }
       else
       {
           file<<std::setprecision(10)<<str<<',';
       }
    }
    file<<'\n';
  }
}


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

} // end sortPointCloud




gsMatrix<> pointWiseErrors(const gsTensorBSpline<2, real_t> & geometry,
                                  const gsMatrix<> & parameters,
                                  const gsMatrix<> & points)
{
  gsMatrix<> eval;
  geometry.eval_into(parameters, eval);
  gsMatrix<> errorMat(1, eval.cols());
  for (index_t col = 0; col != eval.cols(); col++)
  {
      errorMat(0, col) = (eval.col(col) - points.col(col)).norm();
  }
  return errorMat;
}



std::vector<real_t> computeErrors(const gsTensorBSpline<2, real_t> & geometry,
                                  const gsMatrix<> & parameters,
                                  const gsMatrix<> & points)
{
  std::vector<real_t> min_max_mse;
  // gsMatrix<> eval;
  // geometry.eval_into(parameters, eval);
  //
  gsMatrix<> errorMat = pointWiseErrors(geometry, parameters, points);
  //
  // for (index_t col = 0; col != eval.cols(); col++)
  // {
  //     pointWiseErrors(0, col) = (eval.col(col) - points.col(col)).norm();
  // }

  real_t min_error = 1e6;
  real_t max_error = 0;
  real_t mse_error = 0;

  gsMatrix<> err_point(3,1);

  for (index_t i = 1; i < errorMat.cols(); i++)
  {
    const real_t err = errorMat(0,i) ;
    mse_error += err * err ;
    if ( err > max_error )
    {
      max_error = err;
      err_point << points.col(i);
      gsWriteParaviewPoints(err_point, "piMax");
      geometry.eval_into(parameters.col(i), err_point);
      gsWriteParaviewPoints(err_point, "siMax");
    }

    if ( err < min_error )
    {
      min_error = err;
    }
  }

  min_max_mse.push_back(min_error);
  min_max_mse.push_back(max_error);
  min_max_mse.push_back(mse_error/errorMat.cols());

  return min_max_mse;
}


int main(int argc, char *argv[])
{
    // a
    bool corner_conditions = false; // b
    bool boundary_curves = false; // c
    index_t deg = 2; // d
    // e
    std::string fn = "../filedata/fitting/shiphull_scalePts.xml"; // f
    bool uvcorrection = false; // g
    // h, i, j, k, l,
    real_t sigma = 0; // m
    index_t numKnots = 5; // n
    // o
    real_t mu = 1; // p
    index_t deg_x = -1; // q
    index_t deg_y = -1; // r
    real_t lambda = 0; // s
    // t, u, v, w,
    index_t kx = -1; // x
    index_t ky = -1; // y
    // z

    gsCmdLine cmd("Tensor product B-spline surface fitting with Tangent Distance Minimization.");
    // a
    cmd.addSwitch("b", "boundary", "add boundary constraints to the fitting.", boundary_curves);
    cmd.addSwitch("c", "corners", "add conrners constraints to the fitting.", corner_conditions);
    cmd.addInt("d", "degree", "bi-degree (d,d).", deg);
    // e
    cmd.addString("f", "filename", "name of the .xml file containing the data", fn);
    cmd.addSwitch("g", "uvc", "apply parameter correction.", uvcorrection);
    // h, i, j, k, l,
    cmd.addReal("m", "tdm", "add TDM to the system", sigma);
    cmd.addInt("n", "interiors", "number of interior knots in each direction.", numKnots);
    // o
    cmd.addReal("p", "pdm", "add PDM to the system", mu);
    cmd.addInt("q", "xdeg", "x-degree (dx,dy).", deg_x);
    cmd.addInt("r", "ydeg", "y-degree (dy,dy).", deg_y);
    cmd.addReal("s", "smoothing", "smoothing weight", lambda);
    // t, u, v, w,
    cmd.addInt("x", "xknt", "number of interior knots in x-direction.", kx);
    cmd.addInt("y", "yknt", "number of interior knots in y-direction.", ky);
    // z


    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    if(deg_x < 0)
      deg_x = deg;
    if(deg_y < 0)
        deg_y = deg;

    if(kx < 0)
      kx = numKnots;
    if(ky < 0)
      ky = numKnots;


    gsStopwatch gsTime;
    time_t now = time(0);

    gsFileData<> fd_in(fn);
    gsMatrix<> uv, xyz, X;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);
    //! [Read data]

    scalePoints(xyz, X);

    GISMO_ENSURE( uv.cols() == X.cols() && uv.rows() == 2 && X.rows() == 3,
                  "Wrong input, check id of matrices in the .xml file");

    std::vector<index_t> interpIdx;
    sortPointCloud(uv, X, interpIdx);



    std::ofstream pdm_results;
    pdm_results.open(std::to_string(now) + "pdm_results.csv");
    pdm_results << "m, deg, mesh, dofs, pc, pen, min, max, mse\n";

    std::ofstream tdm_results;
    tdm_results.open(std::to_string(now) + "tdm_results.csv");
    // tdm_results << "m, deg, mesh, dofs, pc, pen, min, max, mse, mu*PDM, sigma*TDM, boundary, pos-semi-def, cond, s_min, s_max, s-free\n";
    tdm_results << "m, deg, mesh, dofs, pc, pen, min, max, mse, mu*PDM, sigma*TDM, boundary, cond, s_min, s_max, s-free\n";

    // for loop on space dimension: n = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
      // 0 PC and 1 PC
      // impose boundary conditions
      // combine with PDM: m = 1, p = 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1


    for (std::vector<index_t>::iterator it = interpIdx.begin(); it != interpIdx.end(); ++it)
    {
      gsInfo << "idx = " << *it << ", params = " << uv.col(*it).transpose() << ", points = " << X.col(*it).transpose() <<"\n";
      gsMatrix<> uvcol = uv.col(*it);
      gsMatrix<> xcol = X.col(*it);
      gsWriteParaviewPoints(uvcol, "uv_"+internal::to_string(*it));
      gsWriteParaviewPoints(xcol, "p_"+internal::to_string(*it));
    }

    gsWriteParaviewPoints(uv, "parameters");
    gsWriteParaviewPoints(X, "points");



    // GISMO_ENSURE(lambda == 0, "Smoothing matrix not available (yet).");

    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
        u_max = uv.row(0).maxCoeff(),
        v_min = uv.row(1).minCoeff(),
        v_max = uv.row(1).maxCoeff();

    // Create knot-vectors without interior knots
    gsKnotVector<> u_knots (u_min, u_max, kx, deg_x+1 ) ;
    gsKnotVector<> v_knots (v_min, v_max, ky, deg_y+1 ) ;

    // Create a tensor-basis and apply initial uniform refinement
    gsTensorBSplineBasis<2> basis( u_knots, v_knots );

    gsFitting<real_t> fitting_object(uv, X, basis);
    if(uvcorrection)
    {
      fitting_object.compute(lambda);
      fitting_object.parameterCorrection(1e-10,1,1e-8);
      fitting_object.compute(lambda);
    }
    else
    {
      fitting_object.compute(lambda);
    }


    gsMatrix<> coefs(basis.size(), 3);
    coefs = fitting_object.result()->coefs();

    gsTensorBSpline<2, real_t> original(basis, coefs);

    gsInfo << "Least square fitting, with smoothing:\n " << original << "\n";
    gsWriteParaview( original, "pdm_0", 1000, true, false);
    gsWriteParaview( original, "pdm_0", 1000, false, true);
    // writeSingleControlNet(original, "pdm_0_cnet");

    std::vector<real_t> pdm_min_max_mse = computeErrors(original,uv,X);


    pdm_results << X.cols() << "," << deg << "," << kx << "x" << ky << ","
                << basis.size() << "," << uvcorrection << ","
                << lambda << ","
                << pdm_min_max_mse[0] << std::scientific << ","
					      << pdm_min_max_mse[1] << std::scientific << ","
					      << pdm_min_max_mse[2] << std::scientific << "\n";
    pdm_results.close();

    gsMatrix<> pdmColors(4, X.cols());
    pdmColors << X.row(0), X.row(1), X.row(2),pointWiseErrors(original,uv,X);

    gsWriteParaviewPoints(pdmColors, "pdm_colors");


    // Standard least square approximation is our starting point.

    bool psdMatrix = true;
    // Now we perform 1 step of parameter parameterCorrection
    // 1 step of control point update.
    // These two steps represent one iteration of TDM method.
    gsMatrix<> params(2, X.cols());
    // gsDebugVar(coefs);
    gsMatrix<> sol(coefs.rows(), coefs.cols());
    sol = coefs;
    index_t iter = 0;
    //for (index_t iter = 0; iter < 10; iter ++)
    {
      // Apply parameter correction
      if (uvcorrection)
      {
        gsInfo << "Apply parameter correction:\n";
        for (index_t i = 0; i < X.cols(); ++i)
        {
          gsVector<> newParam;
          original.closestPointTo(X.col(i).transpose(), newParam, 1e-10, true);
          params.col(i) = newParam;
        }
      }
      else
      {
        gsInfo << "NO parameter correction:\n";
        params << uv;
      }



      gsMapData<> md(NEED_VALUE|NEED_NORMAL);
      md.points = params; // parametric values

      original.computeMap(md); // gsGeometry->computeMap(md);
      gsInfo << "Succesful map computation.\n";
      gsMatrix<> values  = md.values[0]; // a 3 x u.cols() matrix
      gsMatrix<> normals = md.normals;   // a 3 x u.cols() matrix

      // gsExprEvaluator<> ev;
      // // auto G = ev.getMap(tdm_surface);
      // auto G = ev.getMap(original);
      // gsMatrix<> normals(3, params.cols());
      // gsInfo << "Normals computation:\n";
      // for(index_t j=0; j< params.cols(); j++)
      // {
      //   normals.col(j) = ev.eval(sn(G).normalized(),params.col(j));
      //   // gsInfo << "n^1 : " << normals(0,j) << "\n";
      //   // gsInfo << "n^2 : " << normals(1,j) << "\n";
      //   // gsInfo << "n^3 : " << normals(2,j) << "\n\n";
      // }

      gsMatrix<> N_x(params.cols(), 1);
      gsMatrix<> N_y(params.cols(), 1);
      gsMatrix<> N_z(params.cols(), 1);

      for (index_t k = 0; k < normals.cols(); ++k)
      {
        normals.col(k).normalize();
        N_x(k,0) = normals(0,k);
        N_y(k,0) = normals(1,k);
        N_z(k,0) = normals(2,k);
        // gsInfo << "nx: " << N_x(k,0) << ", ny: " << N_y(k,0) << ", nz: " << N_z(k,0) << "\n\n";
      }

      gsInfo << "Normals from expression evaluator: " << normals.rows() << " x " << normals.cols() << "\n";
      gsInfo << " ? Unit normals ? " << normals.col(0).norm() << "\n";

      gsInfo << "Points: " << X.rows() << " x " << X.cols() << "\n";
      gsInfo << "Normals: " << normals.rows() << " x " << normals.cols() << "\n";


      gsSparseMatrix<> sparseColloc(X.cols(), basis.size());
      sparseColloc = basis.collocationMatrix( params ) ;
      gsMatrix<> tmp = sparseColloc;
      gsInfo << "Collocation matrix: " << tmp.rows() << " x " << tmp.cols() << "\n";

      gsWriteParaview(basis, "basis");

  //#   pragma omp parallel for default(shared) private(curr_point,actives,value)
      gsInfo << "assemble TDM system:\n";

      gsInfo << "collocation matrix ...\n";
      // gsMatrix<> B_mat(X.cols() * 3, basis.size() * 3);
      gsSparseMatrix<> B_mat(X.cols() * 3, basis.size() * 3);
      gsMatrix<> Bb(X.cols() * 3, basis.size() * 3);
      Bb.setZero();
      Bb.block(0,0,tmp.rows(),tmp.cols()) = tmp;
      Bb.block(tmp.rows(),tmp.cols(),tmp.rows(),tmp.cols()) = tmp;
      Bb.block(2*tmp.rows(),2*tmp.cols(),tmp.rows(),tmp.cols()) = tmp;
      gsInfo << "... assembled.\n";
      B_mat = Bb.sparseView();


      // gsInfo << "~B:\n" << B_mat << "\n";

      gsInfo << "Nomral matrix ...\n";
      // gsMatrix<> N_diag(X.cols() * 3, X.cols());
      gsSparseMatrix<> N_diag(X.cols() * 3, X.cols());
      N_diag.setZero();
      for(index_t i = 0; i < X.cols(); i++)
      {
        N_diag(i,i) = N_x(i,0);
        N_diag(X.cols()+i,i) = N_y(i,0);
        N_diag(2*X.cols()+i,i) = N_z(i,0);
      }
      gsInfo << "... assembled.\n";


      // N_diag(0,0) = 1;
      // N_diag(X.cols()-1, X.cols()-1) = 1;
      // N_diag(X.cols(), X.cols()-1) = 1;
      // N_diag(2*X.cols()-1, X.cols()-1) = 1;

      gsMatrix<> X_tilde(X.cols() * 3, 1);
      X_tilde << X.row(0).transpose(), X.row(1).transpose(), X.row(2).transpose();

      // gsMatrix<> A_tilde = N_diag.transpose() * B_mat;
      // gsMatrix<> rhs = N_diag.transpose() * X_tilde;


      // gsMatrix<> Im(3*X.cols(), 3*X.cols());
      gsSparseMatrix<> Im(3*X.cols(), 3*X.cols());
      Im.setIdentity();

      // N_diag.block(0,0,X.cols(), X.cols()) += Im;
      // N_diag.block(X.cols(),0, X.cols(), X.cols()) += Im;
      // N_diag.block(2*X.cols(),0, X.cols(), X.cols()) += Im;

      // Initialization of the smoothing matrix that we need to define the objective function.
      gsSparseMatrix<> m_G;
      m_G.resize(fitting_object.result()->coefs().rows(), fitting_object.result()->coefs().rows());
      m_G.reservePerColumn( cast<real_t,index_t>( (2 * fitting_object.result()->basis().maxDegree() + 1) * 1.333 ) );
      fitting_object.applySmoothing(lambda, m_G);

      gsInfo << "Smoothing matrix: " << m_G.rows() << " x " << m_G.cols() << "\n";


      gsInfo << "System matrix ...\n";
      gsTime.restart();
      // gsMatrix<> A_tilde = B_mat.transpose() * ( mu * Im + sigma * N_diag * N_diag.transpose()) * B_mat;
      gsSparseMatrix<> A_tilde = B_mat.transpose() * ( mu * Im + sigma * N_diag * N_diag.transpose()) * B_mat;

      gsSparseMatrix<> G_mat(A_tilde.rows(), A_tilde.cols());
      gsMatrix<> Gg(A_tilde.rows(), A_tilde.cols());
      Bb.setZero();
      Gg.block(0,0,m_G.rows(),m_G.cols()) = m_G;
      Gg.block(m_G.rows(),m_G.cols(),m_G.rows(),m_G.cols()) = m_G;
      Gg.block(2*m_G.rows(),2*m_G.cols(),m_G.rows(),m_G.cols()) = m_G;
      gsInfo << "... assembled.\n";
      G_mat = Gg.sparseView();

      A_tilde += lambda * G_mat;
      gsTime.stop();
      gsInfo << "... assembled in "<< gsTime << ".\n";
      gsInfo << "System rhs ...\n";
      gsMatrix<> rhs = (X_tilde.transpose() * (mu * Im + sigma * N_diag * N_diag.transpose()) * B_mat).transpose();
      gsInfo << "... assembled.\n";

      gsInfo << A_tilde.rows() << " x " << A_tilde.cols()  << "\n";
      gsInfo << rhs.rows() << " x " << rhs.cols()  << "\n";

      // gsMatrix<> Im(X.cols(), X.cols());
      // Im.setIdentity();
      // N_diag.block(0,0,X.cols(), X.cols()) += Im;
      // N_diag.block(X.cols(),0, X.cols(), X.cols()) += Im;
      // N_diag.block(2*X.cols(),0, X.cols(), X.cols()) += Im;

      // gsMatrix<> A_tilde =  B_mat;
      // gsMatrix<> rhs = X_tilde;


      // gsMatrix<> b_x(N_x.rows(), N_x.cols());
      // b_x = X.row(0).transpose().cwiseProduct(normals.row(0).transpose());
      // gsMatrix<> b_y(N_y.rows(), N_y.cols());
      // b_y = X.row(1).transpose().cwiseProduct(normals.row(1).transpose());
      // gsMatrix<> b_z(N_z.rows(), N_z.cols());
      // b_z = X.row(2).transpose().cwiseProduct(normals.row(2).transpose());
      // gsMatrix<> m_B(N_x.rows(), 1);
      // m_B = b_x+b_y+b_z;

      // gsInfo << "right hand side difference: " << (rhs - m_B).norm() << "\n";

      gsEigen::LLT<gsMatrix<>::Base> A_llt(A_tilde);
      if (!A_tilde.isApprox(A_tilde.transpose()) || A_llt.info() == gsEigen::NumericalIssue)
      {
        psdMatrix = false;
        gsInfo << "Possibly non semi-positive definitie matrix!\n";
      }
      // // check the system matrix
      // gsEigen::JacobiSVD<gsMatrix<>::Base> svd(A_tilde);
      // real_t sigma_min = svd.singularValues()(svd.singularValues().size()-1);
      // real_t sigma_max = svd.singularValues()(0);
      // real_t cond = sigma_max/sigma_min;
      //
      // gsInfo << "System matrix condition number: " << cond << "\n";
      // gsInfo << "System matrix maximum singular values: " << sigma_max << "\n";
      // gsInfo << "System matrix minimum singular values: " << sigma_min << "\n";
      //
      // gsMatrix<real_t> Rb = (svd.singularValues().array() < 1e-6).cast<real_t>();;
      // // gsInfo << Rb << "\n";
      // gsInfo << "Degrees of freedom: " << Rb.sum() << "\n";


      index_t c1 = 0;
      index_t c2 = basis.component(0).size()-1;
      index_t c3 = (basis.component(1).size()-1) *  basis.component(0).size() ;
      index_t c4 = basis.size()-1;

      gsInfo << "System matrix size: " << A_tilde.rows() << " x " << A_tilde.cols() << "\n";
      gsInfo << "basis at the domain extreme = " << c1 << ", " << c2 << ", "<< c3 << ", " << c4 << "\n";
      gsInfo << "indeces at the domain extreme = " << c1 + basis.size() << ", " << c2 + basis.size()<< ", " << c3 + basis.size() << ", " << c4 + basis.size() << "\n";
      gsInfo << "indeces at the domain extreme = " << c1 + 2*basis.size() << ", " << c2 + 2*basis.size()<< ", " << c3 + 2*basis.size() << ", " << c4 + 2*basis.size() << "\n";

      std::vector<boxSide> fixedSides;

      fixedSides.push_back(boxSide(3));
      fixedSides.push_back(boxSide(2));
      fixedSides.push_back(boxSide(4));
      fixedSides.push_back(boxSide(1));

      for(std::vector<boxSide>::const_iterator it=fixedSides.begin(); it!=fixedSides.end(); ++it)
        gsInfo << basis.boundary(*it).transpose() << "\n";


      std::vector<index_t> idxConstraints;
      std::vector<gsMatrix<> > coefsConstraints;

      for(std::vector<boxSide>::const_iterator it=fixedSides.begin(); it!=fixedSides.end(); ++it)
      {
        gsMatrix<index_t> ind = basis.boundary(*it);
        for(index_t r=0; r<ind.rows(); r++)
        {
          index_t fix = ind(r,0);
          // If it is a new constraint, add it.
          if(std::find(idxConstraints.begin(), idxConstraints.end(), fix) == idxConstraints.end())
          {
            idxConstraints.push_back(fix);
            coefsConstraints.push_back(original.coef(fix));
          }
        }
      }

      gsInfo << "First coefficient to contrain:\n";
      gsInfo << coefsConstraints[0] << "\n";


      if(boundary_curves)
      {
        gsInfo << "Constrains on the boundary curves.\n";
        // for(std::vector<index_t>::const_iterator it=idxConstraints.begin(); it!=idxConstraints.end(); ++it)

      //   for (index_t k=0; k<m_constraintsLHS.outerSize(); ++k)
      //   {
    	//   for (typename gsSparseMatrix<T>::InnerIterator it(m_constraintsLHS,k); it; ++it)
    	//     {
    	//       A_mat(basisSize + it.row(), it.col()) = it.value();
    	//       A_mat(it.col(), basisSize + it.row()) = it.value();
    	//     }
      //   }

        gsMatrix<> A_comp = A_tilde;

        for(index_t el = 0; el < idxConstraints.size(); ++el)
        {
          for (int k=0; k<B_mat.outerSize(); ++k)
          {
            for (gsSparseMatrix<>::InnerIterator it(A_tilde,k); it; ++it)
            {
              if (it.row() == idxConstraints[el])
                A_tilde(it.row(), it.col()) = 0;
              if (it.row() == idxConstraints[el] && it.col() == idxConstraints[el])
                A_tilde(it.row(), it.col()) = 1;

              if (it.row() == idxConstraints[el] + basis.size())
                A_tilde(it.row(), it.col()) = 0;
              if (it.row() == idxConstraints[el] + basis.size() && it.col() == idxConstraints[el] + basis.size())
                A_tilde(it.row(), it.col()) = 1;

              if (it.row() == idxConstraints[el] + 2*basis.size())
                A_tilde(it.row(), it.col()) = 0;
              if (it.row() == idxConstraints[el] + 2*basis.size() && it.col() == idxConstraints[el] + 2*basis.size())
                A_tilde(it.row(), it.col()) = 1;

            }
          }

          rhs(idxConstraints[el], 0) = coefsConstraints[el](0,0); // x-component
          rhs(idxConstraints[el] + basis.size(),0) = coefsConstraints[el](0,1); // y-component
          rhs(idxConstraints[el] + 2*basis.size(),0) = coefsConstraints[el](0,2); // z-component
        }
      // } boundary curves constraints



      for(index_t el = 0; el < idxConstraints.size(); ++el)
      {
        // gsInfo << el << "\n";
        // gsInfo << coefsConstraints[el] << "\n";

        A_comp.row(idxConstraints[el]) = gsVector<>::Zero(A_comp.cols(),1);
        A_comp( idxConstraints[el], idxConstraints[el] ) = 1;

        A_comp.row( idxConstraints[el] + basis.size()) = gsVector<>::Zero(A_comp.cols(),1);
        A_comp( idxConstraints[el] + basis.size(), idxConstraints[el] + basis.size() ) = 1;

        A_comp.row( idxConstraints[el] + 2*basis.size()) = gsVector<>::Zero(A_comp.cols(),1);
        A_comp( idxConstraints[el] + 2*basis.size(), idxConstraints[el] + 2*basis.size() ) = 1;

      }

      gsInfo << "========================================================================================================================\n";
      gsInfo << "Check constraints application: " << (A_tilde - A_comp.sparseView()).norm() << "\n";
      gsInfo << "========================================================================================================================\n";

      gsEigen::LLT<gsMatrix<>::Base> A_bc(A_tilde);
      if (!A_tilde.isApprox(A_tilde.transpose()) || A_bc.info() == gsEigen::NumericalIssue)
      {
        psdMatrix = false;
        gsInfo << "After imposing boundary constrains:\n+++Possibly non semi-positive definitie matrix!+++\n";
      }
      else
      {
        psdMatrix = true;
      }

    }

      if(corner_conditions)
      {
        gsInfo << "Constrains on the domain corners.\n";
        // corners conditions
        // x-component
        A_tilde.row(c1) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c1, c1 ) = 1;

        A_tilde.row( c2 ) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c2 , c2 ) = 1;

        A_tilde.row( c3 ) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c3 , c3 ) = 1;

        A_tilde.row( c4 ) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c4 , c4 ) = 1;
        //

        // y-component
        A_tilde.row(c1 + basis.size()) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c1 + basis.size(), c1 + basis.size() ) = 1;

        A_tilde.row( c2 + basis.size()) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c2 + basis.size(), c2 + basis.size()) = 1;

        A_tilde.row( c3 + basis.size()) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c3 + basis.size(), c3 + basis.size()) = 1;

        A_tilde.row( c4 + basis.size()) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c4 + basis.size() , c4 + basis.size()) = 1;
        //

        // z-component
        A_tilde.row(c1 + 2*basis.size()) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c1 + 2*basis.size(), c1 + 2*basis.size() ) = 1;

        A_tilde.row( c2 + 2*basis.size()) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c2 + 2*basis.size(), c2 + 2*basis.size() ) = 1;

        A_tilde.row( c3 + 2*basis.size()) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c3 + 2*basis.size(), c3 + 2*basis.size() ) = 1;

        A_tilde.row( c4 + 2*basis.size()) = gsVector<>::Zero(A_tilde.cols(),1);
        A_tilde( c4 + 2*basis.size(), c4 + 2*basis.size()) = 1;
        //

        //
        rhs(c1,0) = X(0,interpIdx[0]); // x-component
        rhs(c2,0) = X(0,interpIdx[1]); // x-component
        rhs(c3,0) = X(0,interpIdx[3]); // x-component
        rhs(c4,0) = X(0,interpIdx[2]); // x-component

        rhs(c1 + basis.size(),0) = X(1,interpIdx[0]); // y-component
        rhs(c2 + basis.size(),0) = X(1,interpIdx[1]); // y-component
        rhs(c3 + basis.size(),0) = X(1,interpIdx[3]); // y-component
        rhs(c4 + basis.size(),0) = X(1,interpIdx[2]); // y-component

        rhs(c1 + 2*basis.size(),0) = X(2,interpIdx[0]); // z-component
        rhs(c2 + 2*basis.size(),0) = X(2,interpIdx[1]); // z-component
        rhs(c3 + 2*basis.size(),0) = X(2,interpIdx[3]); // z-component
        rhs(c4 + 2*basis.size(),0) = X(2,interpIdx[2]); // z-component

        gsEigen::LLT<gsMatrix<>::Base> A_cc(A_tilde);
        if (A_cc.info() == gsEigen::NumericalIssue)
        {
          psdMatrix = false;
          gsInfo << "After imposing corner constrains:\n+++Possibly non semi-positive definitie matrix!+++\n";
        }
        else
        {
          psdMatrix = true;
        }
      }

      //
      // for(index_t x=0; x<rhs.size(); x++)
      //   gsInfo << x << ", " << rhs(x,0) << "\n";


      // gsSparseSolver<>::QR qrsolver_tilde(A_tilde.sparseView());
      // gsSparseSolver<>::QR qrsolver_tilde(A_tilde);
      // gsMatrix<> sol_tilde = qrsolver_tilde.solve(rhs);

      // check the system matrix
      gsEigen::JacobiSVD<gsMatrix<>::Base> svd(A_tilde);
      real_t sigma_min = svd.singularValues()(svd.singularValues().size()-1);
      real_t sigma_max = svd.singularValues()(0);
      real_t cond = sigma_max/sigma_min;

      gsInfo << "System matrix condition number: " << cond << "\n";
      gsInfo << "System matrix maximum singular values: " << sigma_max << "\n";
      gsInfo << "System matrix minimum singular values: " << sigma_min << "\n";

      gsMatrix<real_t> Rb = (svd.singularValues().array() < 1e-8).cast<real_t>();;
      // gsInfo << Rb << "\n";
      gsInfo << "Degrees of freedom: " << Rb.sum() << "\n";


      A_tilde.makeCompressed();

      typename gsSparseSolver<real_t>::BiCGSTABILUT solver( A_tilde );

      if ( solver.preconditioner().info() != gsEigen::Success )
      {
          gsWarn<<  "The preconditioner failed. Aborting.\n";

          return 0;
      }
      // Solves for many right hand side  columns

      gsMatrix<> sol_tilde = solver.solve(rhs); //toDense()



      // gsInfo << "Tilde solution:\n";
      // gsInfo << sol_tilde << "\n";

      gsMatrix<> coefs_tilde(basis.size(), 3);
      for(index_t j=0; j<basis.size(); j++)
      {
        coefs_tilde(j,0) = sol_tilde(j);
        coefs_tilde(j,1) = sol_tilde(basis.size()+j);
        coefs_tilde(j,2) = sol_tilde(2*basis.size()+j);
      }

      gsTensorBSpline<2, real_t> surface_tilde(basis, coefs_tilde);
      gsWriteParaview( surface_tilde, "tdm_surf", 1000);
      gsWriteParaview( surface_tilde, "tdm_cnet", 1000, false, true);

      std::vector<real_t> tdm_min_max_mse = computeErrors(surface_tilde,params,X);

      tdm_results << X.cols() << "," << deg << "," << kx << "x" << ky << ","
                  << basis.size() << ","<< uvcorrection << ","
                  << lambda << ","
                  << tdm_min_max_mse[0] << std::scientific << ","
  					      << tdm_min_max_mse[1] << std::scientific << ","
  					      << tdm_min_max_mse[2] << std::scientific << ","
                  << mu << ","
                  << sigma << ","
                  << boundary_curves << ","
                  //<< psdMatrix << ","
                  << cond << "," << sigma_min << "," << sigma_max << ","
                  << Rb.sum() << "\n";
      tdm_results.close();

      gsMatrix<> tdmColors(4, X.cols());
      tdmColors << X.row(0), X.row(1), X.row(2), pointWiseErrors(surface_tilde,params,X);

      gsWriteParaviewPoints(tdmColors, "tdm_colors");

    } // end A-TDM loop

}
