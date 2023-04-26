/** @file fitting_example.cpp

    @brief Demonstrates adaptive fitting of data samples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <gismo.h>

using namespace gismo;

template<class T>
gsMatrix<T> diff(const gsMatrix<T> & E)
{
  gsMatrix<T> E1 = E.block(0, 0, E.rows()-1, E.cols());
  gsMatrix<T> E2 = E.block(1, 0, E.rows()-1, E.cols());
  return E2 - E1;
}

template<class T>
gsMatrix<T> diff(const gsMatrix<T> & E, const index_t d)
{
  GISMO_ASSERT( E.rows() <= d, "Wrong regulatization amount.");
  if (d == 0)
    return E;
  else
  {
    index_t r = 1;
    gsMatrix<T> D = E;
    while( r <= d)
    {
      D = diff(D);
      r += 1;
    }
    return D;
  }
}

template<class T>
gsMatrix<T> gsKroneckerIdxProduct(const gsMatrix<T>& A, const gsMatrix<T>& B)
{
  gsMatrix<T> C(B.rows() * A.rows(), B.cols() * A.cols());

  for (index_t i = 0; i < A.rows(); i++)
  {
    for(index_t j = 0; j < A.cols(); j++)
    {
      for(index_t k = 0; k < B.rows(); k++)
      {
        for(index_t m = 0; m < B.cols(); m++)
        {
          C(i * B.rows()+ k, j * B.cols() + m) = A(i,j) * B(k,m);
        }
      }
    }
  }
  return C;
}

template<class T>
gsMatrix<T> gsKroneckerMatProduct(const gsMatrix<T>& A, const gsMatrix<T>& B)
{
  gsMatrix<T> C(B.rows() * A.rows(), B.cols() * A.cols());
  C.setZero();

  for(index_t i = 0; i < A.rows(); i++)
  {
    for(index_t j = 0; j < A.cols(); j++)
    {
      C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i,j) * B;
    }
  }
  return C;
}


int main(int argc, char *argv[])
{

  index_t maxPcIter = 1;
  index_t deg_x = 3;
  index_t deg_y = 3;
  index_t numInt = 7; // number of interior points in one direction
  index_t d1 = 2;
  index_t d2 = 2;
  real_t lambda1 = 1e-6;
  real_t lambda2 = 1e-6;
  std::ofstream file_results;
  file_results.open("carShell2_fitting_example_results.csv");

  std::vector<std::string> filenames;
  filenames.push_back("../../testdata/carShell2/carShell2_m55_scale_m55_deg3.xml"); // cnn data
  filenames.push_back("../../testdata/carShell2/carShell2_m55_scale.xml"); // standard data

  // Reading options from the command line
  gsCmdLine cmd("Penalized tensor-product b-spline approximation.");
  cmd.addInt("c", "parcor", "Steps of parameter correction", maxPcIter);
  cmd.addInt("x", "deg_x", "degree in x direction", deg_x);
  cmd.addInt("y", "deg_y", "degree in y direction", deg_y);
  // cmd.addReal("l1", "lambda1", "smoothing coefficient", lambda1);
  // cmd.addInt("d1", "penalization1", "penalization order in direction 0", d1);
  // cmd.addReal("l2", "lambda2", "smoothing coefficient", lambda2);
  // cmd.addInt("d2", "penalization2", "penalization order in direction 1", d2);
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  gsInfo << "Fine-tuning of the P-splines tensor product approximation.\n";

  index_t maxP1 = numInt + deg_x + 1;
  index_t maxP2 = numInt + deg_y + 1;

  index_t best_d1;
  index_t best_d2;
  real_t best_l1;
  real_t best_l2;

  real_t last_max = 1e3;

  for(index_t pen1 = 1; pen1 < maxP1-1; pen1++)
  {
    for(index_t pen2 = 1; pen2 < maxP2-1; pen2++)
    {
      real_t lambda1 = 1e-9;
      real_t lambda2 = 1e-9;
      for(index_t tpow1 = 0; tpow1 < 10; tpow1++)
      {
        for(index_t tpow2 = 0; tpow2 < 10; tpow2++)
        {
          std::string fn = filenames[0];
          gsFileData<> fd_in(fn);
          gsMatrix<> uv, xyz;
          fd_in.getId<gsMatrix<> >(0, uv );
          fd_in.getId<gsMatrix<> >(1, xyz);
          real_t u_min = uv.row(0).minCoeff(), u_max = uv.row(0).maxCoeff(),
                 v_min = uv.row(1).minCoeff(), v_max = uv.row(1).maxCoeff();
          gsKnotVector<> ku (u_min, u_max, numInt, deg_x+1 ) ;
          gsKnotVector<> kv (v_min, v_max, numInt, deg_y+1 ) ;
          gsTensorBSplineBasis<2> basis( ku, kv );
          gsFitting<real_t> pls(uv, xyz, basis);

          lambda1 = lambda1 * 10;
          lambda2 = lambda2 * 10;

          pls.computePen(lambda1, pen1, lambda2, pen2); // least square with penalization
          pls.parameterCorrectionPen(1e-7, maxPcIter, 1e-4, lambda1, pen1, lambda2, pen2);
          pls.computeErrors();

          real_t curr_max = pls.maxPointError();
          if (curr_max < last_max)
          {
            last_max = curr_max;
            best_d1 = pen1;
            best_d2 = pen2;
            best_l1 = lambda1;
            best_l2 = lambda2;
          }
        }
      }
    }
  }

  gsInfo << "Best configurations: lambda1 = " << best_l1 << ", d1 = " << best_d1 << ", lambda2 = " << best_l2 << ", d2 = " << d2 <<"\n";

  file_results << "method, step, lambda1, d1, lambda2, d2, dofs, min, max, mse\n";

  for(index_t method = 0; method < filenames.size(); method++)
  {
    std::string fn = filenames[method];
    std::string pname;
    if(method == 0)
    {
      pname = "cnn";
    }
    else
    {
      pname = "std";
    }

    gsFileData<> fd_in(fn);
    gsMatrix<> uv, xyz;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, xyz);

    gsWriteParaviewPoints(xyz, "points");

    real_t u_min = uv.row(0).minCoeff(), u_max = uv.row(0).maxCoeff(),
           v_min = uv.row(1).minCoeff(), v_max = uv.row(1).maxCoeff();


    gsKnotVector<> ku (u_min, u_max, numInt, deg_x+1 ) ;
    gsKnotVector<> kv (v_min, v_max, numInt, deg_y+1 ) ;

    // Create a tensor-basis nad apply initial uniform refinement
    gsTensorBSplineBasis<2> basis( ku, kv );

    gsFitting<real_t> ls(uv, xyz, basis);
    ls.compute(0); // least square with no penalization
    ls.parameterCorrection(1e-7, maxPcIter, 1e-4);
    ls.computeErrors();
    gsWriteParaview(*ls.result(), "ls_geo");


    real_t sum_of_2errors = 0;
    std::vector<real_t> errors2 = ls.pointWiseErrors();
    for(index_t el = 0; el != errors2.size(); el++){
        sum_of_2errors += std::pow(errors2[el], 2);
    }
    sum_of_2errors = sum_of_2errors / errors2.size();

    gsInfo << "----------------------------------------------------------------------------\n";
    gsInfo << "----------------------------------------------------------------------------\n";
    gsInfo << "standard Least Squares:\n";
    index_t dofs = ls.result()->basis().size();
    real_t minPointError = ls.minPointError();
    real_t maxPointError = ls.maxPointError();
    real_t mseError = sum_of_2errors;
    gsInfo<<"DOFs         : "<< dofs <<"\n";
    gsInfo << "PC step = " << maxPcIter << "\n";
    std::cout << "Max distance : "<< maxPointError << std::scientific <<"\n";
    std::cout << "MSE    error : "<< mseError << std::scientific <<"\n";

    file_results << pname + "-ls" << "," << std::to_string(maxPcIter) << "," << std::to_string(0) << "," << "-" << "," << "-" << "," << "-" << "," << std::to_string(dofs) << "," << std::ostringstream(std::to_string(minPointError)).str() << "," << std::ostringstream(std::to_string(maxPointError)).str() << "," << std::ostringstream(std::to_string(mseError)).str() << "\n";
    gsInfo << "----------------------------------------------------------------------------\n";



    errors2.clear();
    sum_of_2errors = 0;
    gsFitting<real_t> sls(uv, xyz, basis);
    real_t lambda = best_l1;
    sls.compute(lambda); // least square with smoothing
    sls.parameterCorrection(1e-7, maxPcIter, 1e-4);
    sls.computeErrors();
    gsWriteParaview(*sls.result(), "sls_geo");

    sum_of_2errors = 0;
    errors2 = sls.pointWiseErrors();
    for(index_t el = 0; el != errors2.size(); el++){
        sum_of_2errors += std::pow(errors2[el], 2);
    }
    sum_of_2errors = sum_of_2errors / errors2.size();

    gsInfo << "Smoothing Least Squares:\n";
    dofs = sls.result()->basis().size();
    minPointError = sls.minPointError();
    maxPointError = sls.maxPointError();
    mseError = sum_of_2errors;
    gsInfo<<"DOFs         : "<< dofs <<"\n";
    gsInfo << "PC step = " << maxPcIter << "\n";
    std::cout << "Max distance : "<< maxPointError << std::scientific <<"\n";
    std::cout << "MSE    error : "<< mseError << std::scientific <<"\n";
    file_results << pname + "-sls" << "," << std::to_string(maxPcIter) << "," << std::to_string(lambda) << "," << "-" << "," << "-" << "," << "-" << "," << std::to_string(dofs) << "," << std::ostringstream(std::to_string(minPointError)).str() << "," << std::ostringstream(std::to_string(maxPointError)).str() << "," << std::ostringstream(std::to_string(mseError)).str() << "\n";
    gsInfo << "----------------------------------------------------------------------------\n";



    errors2.clear();
    sum_of_2errors = 0;
    gsFitting<real_t> pls(uv, xyz, basis);
    real_t lambda1 = best_l1;
    index_t d1 = best_d1;
    real_t lambda2 = best_l2;
    index_t d2 = best_d2;
    pls.computePen(lambda1, d1, lambda2, d2); // least square with penalization
    pls.parameterCorrectionPen(1e-7, maxPcIter, 1e-4, lambda1, d1, lambda2, d2);
    pls.computeErrors();
    gsWriteParaview(*pls.result(), "pls_geo");

    sum_of_2errors = 0;
    errors2 = pls.pointWiseErrors();
    for(index_t el = 0; el != errors2.size(); el++){
        sum_of_2errors += std::pow(errors2[el], 2);
    }
    sum_of_2errors = sum_of_2errors / errors2.size();

    gsInfo << "Penalized Least Squares:\n";
    dofs = pls.result()->basis().size();
    minPointError = pls.minPointError();
    maxPointError = pls.maxPointError();
    mseError = sum_of_2errors;
    gsInfo<<"DOFs         : "<< dofs <<"\n";
    gsInfo << "PC step = " << maxPcIter << "\n";
    std::cout << "Max distance : "<< maxPointError << std::scientific <<"\n";
    std::cout << "MSE    error : "<< mseError << std::scientific <<"\n";
    file_results << pname + "-pls" << "," << std::to_string(maxPcIter) << "," << std::to_string(lambda1) << "," << std::to_string(d1) << "," << std::to_string(lambda2) << "," << std::to_string(d2) << "," << std::to_string(dofs) << "," << std::ostringstream(std::to_string(minPointError)).str() << "," << std::ostringstream(std::to_string(maxPointError)).str() << "," << std::ostringstream(std::to_string(mseError)).str() << "\n";
    gsInfo << "----------------------------------------------------------------------------\n";

  }


  return 0;
}
