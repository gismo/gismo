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
    real_t sigma = 1; // m
    index_t numKnots = 5; // n
    // o
    real_t mu = 0; // p
    index_t deg_x = -1; // q
    index_t deg_y = -1; // r
    real_t lambda = 0; // s
    // t, u, v, w,
    index_t kx = -1; // x
    index_t ky = -1; // y
    index_t pc_step = 0; // y


    gsCmdLine cmd("Tensor product B-spline surface fitting with Point/Tangent Distance Minimization.");
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
    cmd.addInt("z", "pc", "parameter correction step.", pc_step);


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
    real_t computeTime = 0;
    time_t now = time(0);

    gsFileData<> fd_in(fn);
    gsMatrix<> uv, xyz, X;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);
    //! [Read data]

    scalePoints(xyz, X);

    GISMO_ENSURE( uv.cols() == X.cols() && uv.rows() == 2 && X.rows() == 3,
                  "Wrong input, check id of matrices in the .xml file");

    // Determine the parameter domain by mi/max of parameter values
    real_t u_min = uv.row(0).minCoeff(),
           u_max = uv.row(0).maxCoeff(),
           v_min = uv.row(1).minCoeff(),
           v_max = uv.row(1).maxCoeff();

   // Create knot-vectors without interior knots
   gsKnotVector<> u_knots (u_min, u_max, numKnots, deg_x+1 ) ;
   gsKnotVector<> v_knots (v_min, v_max, numKnots, deg_y+1 ) ;

   // Create a tensor-basis nad apply initial uniform refinement
   gsTensorBSplineBasis<2> tbasis( u_knots, v_knots );
   //tbasis.uniformRefine( (1<<numURef)-1 );

    std::vector<index_t> interpIdx;
    sortPointCloud(uv, X, interpIdx);


    std::ofstream pdm_results;
    pdm_results.open(std::to_string(now) + "pdm_results.csv");
    pdm_results << "m, deg, mesh, dofs, pc, pen, min, max, mse\n";

    std::ofstream tdm_results;
    tdm_results.open(std::to_string(now) + "tdm_results.csv");
    // tdm_results << "m, deg, mesh, dofs, pc, pen, min, max, mse, mu*PDM, sigma*TDM, boundary, pos-semi-def, cond, s_min, s_max, s-free\n";
    tdm_results << "m, deg, mesh, dofs, pc, pen, min, max, mse, mu*PDM, sigma*TDM, boundary, cond, s_min, s_max, s-free\n";

    gsFitting<real_t> pdm_obj( uv, X, tbasis);
    pdm_obj.compute(lambda);
    pdm_obj.parameterCorrection(1e-6,pc_step,1e-6);

    gsWriteParaview(*pdm_obj.result(), "pdm_geo");

    gsFitting<real_t> tdm_obj( uv, X, tbasis);
    tdm_obj.compute_tdm(lambda, mu, sigma, interpIdx);
    tdm_obj.parameterCorrection_tdm(1e-6,pc_step, mu, sigma, interpIdx);

    gsWriteParaview(*tdm_obj.result(), "tdm_geo");


    // gsFitting<2, real_t> tdm_obj( uv, X, tbasis, lambda);


}
