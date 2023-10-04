/** @file gsTensorProduct_opt_decoupled.cpp

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

  corners.push_back(interiors.size()); // c1
  corners.push_back(interiors.size() + b_south.size()); // c2
  corners.push_back(interiors.size() + b_south.size() + b_east.size()); // c3
  corners.push_back(interiors.size() + b_south.size() + b_east.size() + b_north.size()); // c4

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

  p_south.resize(tmp_south.rows(), tmp_south.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_south.col(i) = tmp_south.col(tmp[i]);
  }
  uv_south.transposeInPlace();


  tmp.clear();
  tmp = uv_east.idxByColumn(1);

  p_east.resize(tmp_east.rows(), tmp_east.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_east.col(i) = tmp_east.col(tmp[i]);
  }
  uv_east.transposeInPlace();

  tmp.clear();
  tmp = uv_north.idxByColumn(0);
  std::reverse(tmp.begin(),tmp.end());

  gsVector<T> tcol = uv_north.col(0).reverse();
  uv_north.col(0) = tcol;
  tcol = uv_north.col(1).reverse();
  uv_north.col(1) = tcol;

  p_north.resize(tmp_north.rows(), tmp_north.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_north.col(i) = tmp_north.col(tmp[i]);
  }
  uv_north.transposeInPlace();

  tmp.clear();
  tmp = uv_west.idxByColumn(1);

  tcol = uv_west.col(0).reverse();
  uv_west.col(0) = tcol;
  tcol = uv_west.col(1).reverse();
  uv_west.col(1) = tcol;
  std::reverse(tmp.begin(),tmp.end());

  p_west.resize(tmp_west.rows(), tmp_west.cols());
  for(index_t i = 0; i<tmp.size(); i++)
  {
    p_west.col(i) = tmp_west.col(tmp[i]);
  }
  uv_west.transposeInPlace();


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
    bool boundary_curves = false; // b
    // c
    index_t deg = 2; // d
    // e
    std::string fn = "../filedata/fitting/shiphull_scalePts.xml"; // f
    bool uvcorrection = false; // g
    // h
    index_t maxIter = 0; // i
    //j, k, l,
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
    cmd.addInt("d", "degree", "bi-degree (d,d).", deg);
    // e
    cmd.addString("f", "filename", "name of the .xml file containing the data", fn);
    cmd.addSwitch("g", "uvc", "apply parameter correction.", uvcorrection);
    // h
    cmd.addInt("i", "iter", "number of maximum iterations of the optimization algorithm(s).", maxIter);
    //j, k, l,
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

    std::vector<index_t> interpIdx;
    sortPointCloud(uv, X, interpIdx);

    gsWriteParaviewPoints(uv, "0parameters");
    gsWriteParaviewPoints(X, "points");


    std::ofstream pdm_results;
    pdm_results.open(std::to_string(now) + "pdm_results.csv");
    pdm_results << "m, deg, mesh, dofs, pc, pen, min, max, mse\n";

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
    fitting_object.compute(lambda);

    gsMatrix<> coefs(basis.size(), 3);
    coefs = fitting_object.result()->coefs();

    gsTensorBSpline<2, real_t> original(basis, coefs);

    gsInfo << "Least square fitting, with smoothing:\n " << original << "\n";
    gsWriteParaview( original, "0pdm", 1000, true, false);
    gsWriteParaview( original, "0pdm", 1000, false, true);
    // writeSingleControlNet(original, "pdm_0_cnet");

    std::vector<real_t> pdm_min_max_mse = original.MinMaxMseErrors(uv,X);

    pdm_results << X.cols() << "," << deg << "," << kx << "x" << ky << ","
                << basis.size() << "," << uvcorrection << ","
                << lambda << ","
                << pdm_min_max_mse[0] << std::scientific << ","
					      << pdm_min_max_mse[1] << std::scientific << ","
					      << pdm_min_max_mse[2] << std::scientific << "\n";
    pdm_results.close();

    gsMatrix<> pdmColors(4, X.cols());
    pdmColors << X.row(0), X.row(1), X.row(2), original.pointWiseErrors(uv,X);

    gsWriteParaviewPoints(pdmColors, "pdm_colors");


    // Standard least square approximation is our starting point.
    // Now we perform 1 step of parameter parameterCorrection
    // 1 step of control point update.
    // These two steps represent one iteration of TDM method.
    gsMatrix<> params(2, X.cols());
    // gsDebugVar(coefs);
    gsMatrix<> sol(coefs.rows(), coefs.cols());
    sol = coefs;
    index_t iter = 0;

    gsInfo << "Start alternating loop.\n";
    // alternating loop
    std::ofstream results;
    results.open(std::to_string(now) + "results.csv");
    results << "m, deg, mesh, dofs, it, time, pen, min, max, mse, rmse, mu*PDM, sigma*TDM, boundary, cond, s_min, s_max, s-free\n";

    for (index_t iter = 0; iter < maxIter; iter ++)
    {
      gsInfo << "Iteration ... " << iter << " ...\n";
      // Apply parameter correction
      gsInfo << "Apply parameter correction.\n";
      gsTime.restart();
      for (index_t i = 0; i < X.cols(); ++i)
      {
        gsVector<> newParam;
        original.closestPointTo(X.col(i), newParam, 1e-10, true);
        params.col(i) = newParam;
      }
      computeTime += gsTime.stop();


      gsInfo << "Normals computation.\n";
      gsTime.restart();

      gsExprEvaluator<> ev;
      auto G = ev.getMap(original);
      gsMatrix<> normals(3, params.cols());
      gsSparseMatrix<> N_diag(X.cols() * 3, X.cols());
      N_diag.setZero();
      for(index_t j=0; j< params.cols(); j++)
      {
        normals.col(j) = ev.eval(sn(G).normalized(),params.col(j));
        N_diag(j,j) = normals(0,j);
        N_diag(X.cols()+j,j) = normals(1,j);
        N_diag(2*X.cols()+j,j) = normals(2,j);
      }
      computeTime += gsTime.stop();

      gsInfo << "Block collocation matrix.\n";
      gsTime.restart();

      gsSparseMatrix<> sparseColloc(X.cols(), basis.size());
      sparseColloc = basis.collocationMatrix( params ) ;

      gsMatrix<> tmp = sparseColloc;
      gsMatrix<> Bb(X.cols() * 3, basis.size() * 3);
      Bb.setZero();
      Bb.block(0,0,tmp.rows(),tmp.cols()) = tmp;
      Bb.block(tmp.rows(),tmp.cols(),tmp.rows(),tmp.cols()) = tmp;
      Bb.block(2*tmp.rows(),2*tmp.cols(),tmp.rows(),tmp.cols()) = tmp;

      gsSparseMatrix<> B_mat(X.cols() * 3, basis.size() * 3);
      B_mat = Bb.sparseView();

      gsMatrix<> X_tilde(X.cols() * 3, 1);
      X_tilde << X.row(0).transpose(), X.row(1).transpose(), X.row(2).transpose();
      computeTime += gsTime.stop();

      // gsMatrix<> A_tilde = N_diag.transpose() * B_mat;
      // gsMatrix<> rhs = N_diag.transpose() * X_tilde;

      gsSparseMatrix<> A_tilde(B_mat.cols(), B_mat.cols());
      gsMatrix<> rhs(B_mat.cols(),1);
      if ( mu > 0 && sigma > 0)
      {
        gsInfo << mu << "*PDM + " << sigma << "*TDM.\n";
        gsTime.restart();
        gsSparseMatrix<> Im(3*X.cols(), 3*X.cols());
        Im.setIdentity();
        A_tilde = B_mat.transpose() * ( mu * Im + sigma * N_diag * N_diag.transpose()) * B_mat;
        rhs =  B_mat.transpose() * ( mu * Im + sigma * N_diag * N_diag.transpose() ) * X_tilde ;
        computeTime += gsTime.stop();
      }
      else if (mu > 0 && sigma == 0)
      {
        gsInfo << "PDM.\n";
        gsTime.restart();
        A_tilde = B_mat.transpose() *  B_mat;
        rhs = B_mat.transpose() * X_tilde ;
      }
      else if (sigma > 0 && mu == 0)
      {
        gsInfo << "TDM.\n";
        gsTime.restart();
        gsSparseMatrix<> Im(3*X.cols(), 3*X.cols());
        Im.setIdentity();
        A_tilde = B_mat.transpose() * (  N_diag * N_diag.transpose() ) * B_mat;
        rhs =  B_mat.transpose() * ( N_diag * N_diag.transpose() ) * X_tilde ;
        computeTime += gsTime.stop();
      }
      else
      {
        gsWarn<<  "No available fitting method. Aborting.\n";
        results.close();
        return 0;
      }


      if (lambda > 0)
      {
        gsInfo << "Applying thin-plate energy smoothing...\n";
        gsTime.restart();
        gsSparseMatrix<> m_G;
        m_G.resize(basis.size(), basis.size());
        m_G.reservePerColumn( cast<real_t,index_t>( (2 * basis.maxDegree() + 1) * 1.333 ) );
        fitting_object.applySmoothing(lambda, m_G);

        gsSparseMatrix<> G_mat(A_tilde.rows(), A_tilde.cols());
        gsMatrix<> Gg(A_tilde.rows(), A_tilde.cols());
        Gg.setZero();
        Gg.block(0,0,m_G.rows(),m_G.cols()) = m_G;
        Gg.block(m_G.rows(),m_G.cols(),m_G.rows(),m_G.cols()) = m_G;
        Gg.block(2*m_G.rows(),2*m_G.cols(),m_G.rows(),m_G.cols()) = m_G;
        //gsInfo << "... assembled.\n";
        G_mat = Gg.sparseView();

        A_tilde += lambda * G_mat;
        computeTime += gsTime.stop();
      }


      gsInfo << A_tilde.rows() << " x " << A_tilde.cols()  << "\n";
      gsInfo << rhs.rows() << " x " << rhs.cols()  << "\n";

      gsEigen::LLT<gsMatrix<>::Base> A_llt(A_tilde);
      if (!A_tilde.isApprox(A_tilde.transpose()) || A_llt.info() == gsEigen::NumericalIssue)
        gsInfo << "Possibly non semi-positive definitie matrix!\n";

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


      gsEigen::JacobiSVD<gsMatrix<>::Base> svd(A_tilde);
      real_t sigma_min = svd.singularValues()(svd.singularValues().size()-1);
      real_t sigma_max = svd.singularValues()(0);
      real_t cond = sigma_max/sigma_min;

      gsInfo << "System matrix condition number: " << cond << "\n";
      gsInfo << "System matrix maximum singular values: " << sigma_max << "\n";
      gsInfo << "System matrix minimum singular values: " << sigma_min << "\n";

      gsMatrix<real_t> Rb = (svd.singularValues().array() < 1e-6).cast<real_t>();

      // gsInfo << Rb << "\n";
      gsInfo << "Almost linearly dependent columns: " << Rb.sum() << "\n";

      if(boundary_curves)
      {
        gsInfo << "Constrains on the boundary curves.\n";

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

          gsDebugVar(idxConstraints[el]);
          gsDebugVar(coefsConstraints[el]);
          rhs(idxConstraints[el], 0) = coefsConstraints[el](0,0); // x-component
          rhs(idxConstraints[el] + basis.size(),0) = coefsConstraints[el](0,1); // y-component
          rhs(idxConstraints[el] + 2*basis.size(),0) = coefsConstraints[el](0,2); // z-component
        }
          // check the system matrix
          gsEigen::LLT<gsMatrix<>::Base> A_bc(A_tilde);
          if (!A_tilde.isApprox(A_tilde.transpose()) || A_bc.info() == gsEigen::NumericalIssue)
            gsInfo << "After imposing boundary constrains:\n +++Possibly non semi-positive definitie matrix!+++\n";


          gsEigen::JacobiSVD<gsMatrix<>::Base> A_svd(A_tilde);
          real_t sigma_min = A_svd.singularValues()(svd.singularValues().size()-1);
          real_t sigma_max = A_svd.singularValues()(0);
          real_t cond = sigma_max/sigma_min;

          gsInfo << "System matrix condition number: " << cond << "\n";
          gsInfo << "System matrix maximum singular values: " << sigma_max << "\n";
          gsInfo << "System matrix minimum singular values: " << sigma_min << "\n";
          gsMatrix<real_t> Rb = (A_svd.singularValues().array() < 1e-6).cast<real_t>();
          // gsInfo << Rb << "\n";
          gsInfo << "Almost linearly dependent columns: " << Rb.sum() << "\n";
      }

      gsTime.restart();
      A_tilde.makeCompressed();

      typename gsSparseSolver<real_t>::BiCGSTABILUT solver( A_tilde );

      if ( solver.preconditioner().info() != gsEigen::Success )
      {
          gsWarn<<  "The preconditioner failed. Aborting.\n";
          results.close();
          return 0;
      }
      // Solves for many right hand side  columns

      gsMatrix<> sol_tilde = solver.solve(rhs); //toDense()

      computeTime += gsTime.stop();

      gsInfo << "==============================================================================\n";
      gsInfo << "COMPUTATIONAL TIME " <<  computeTime << "\n";

      // gsInfo << "Tilde solution:\n";
      // gsInfo << sol_tilde << "\n";

      gsMatrix<> coefs_tilde(basis.size(), 3);
      for(index_t j=0; j<basis.size(); j++)
      {
        coefs_tilde(j,0) = sol_tilde(j);
        coefs_tilde(j,1) = sol_tilde(basis.size()+j);
        coefs_tilde(j,2) = sol_tilde(2*basis.size()+j);
      }

      //gsTensorBSpline<2, real_t> surface_tilde(basis, coefs_tilde);
      original.coefs() = coefs_tilde;
      gsWriteParaview( original, internal::to_string(iter)+"surf", 1000);
      gsWriteParaview( original, internal::to_string(iter)+"cnet", 1000, false, true);

      std::vector<real_t> sol_min_max_mse = original.MinMaxMseErrors(params,X);
      // results << "m, deg, mesh, dofs, it, time, pen, min, max, mse, rmse, mu*PDM, sigma*TDM, boundary, cond, s_min, s_max, s-free\n";
      results << X.cols() << "," << deg << "," << kx << "x" << ky << ","
                  << basis.size() << ","<< iter << "," << computeTime << ","
                  << lambda << ","
                  << sol_min_max_mse[0] << std::scientific << ","
  					      << sol_min_max_mse[1] << std::scientific << ","
  					      << sol_min_max_mse[2] << std::scientific << ","
                  << math::sqrt(sol_min_max_mse[2]) << std::scientific << ","
                  << mu << ","
                  << sigma << ","
                  << boundary_curves << ","
                  << cond << "," << sigma_min << "," << sigma_max << ","
                  << Rb.sum() << "\n";


      gsMatrix<> tdmColors(4, X.cols());
      tdmColors << X.row(0), X.row(1), X.row(2), original.pointWiseErrors(params,X);

      gsWriteParaviewPoints(tdmColors, "tdm_colors");

    } // end A-TDM loop
    results.close();
    return 0;

}
