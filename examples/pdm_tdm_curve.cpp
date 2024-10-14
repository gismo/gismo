/** @file pdm_tdm_curve.cpp

    @brief Tutorial on gsBSpline class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <iostream>
#include <gismo.h>

using namespace gismo;


int main(int argc, char *argv[])
{
    index_t numPts  = 100;
    index_t numKnt = 5;
    index_t deg = 3;
    real_t lambda = 0;
    real_t noiseFactor = 0;

    gsCmdLine cmd("Planar B-spline curve fitting with Tangent Distance Minimization.");
    cmd.addInt("m", "psize", "number of planar points to sample.", numPts);
    cmd.addInt("n", "interiors", "number of interior knots of the basis.", numKnt);
    cmd.addInt("d", "degree", "degree of the basis.", deg);
    cmd.addReal("e", "noise", "noise perturbation to apply at the parameters", noiseFactor);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // Data generation.
    // Make a BSpline curve

    gsKnotVector<> kv(0, 1, 5, 3+1);//start,end,interior knots, start/end multiplicites of knots
    gsMatrix<> coefs(9, 2);

    coefs << 0, 0,
             0., 0.25,
             0.125, 0.5,
             0.25, 0.55,
             0.375, 0.5,
             0.5, 0.25,
             0.625, 0.25,
             0.8, 0.25,
             1, 0.5;

    gsMatrix<> coefs_plot;
    coefs_plot = coefs.transpose();
    gsWriteParaviewPoints(coefs_plot, "coefs_original");
    gsBSpline<> original( kv, give(coefs));

    gsInfo << "Original curve:\n" << original << "\n";

    gsWriteParaview( original, "originalPlanarBsplineCurve", 1000);
    gsWriteParaview( original, "originalcnet", 1000, false, true);


    gsMatrix<> X, parameters, uv;

    gsMatrix<> uline = gsVector<>::LinSpaced(numPts, 0, 1);
    parameters = uline.transpose();
    original.eval_into(parameters, X);

    uv = parameters;
    gsInfo << "Points:\n" << X << "\n";
    gsInfo << "Parameters:\n" << parameters << "\n";
    gsWriteParaviewPoints(X, "points");
//
    if (noiseFactor > 0)
    {
    real_t umin = uv.minCoeff();
    real_t umax = uv.maxCoeff();

    real_t urange = umax - umin; // umax (=1) - umin (=0)
    gsMatrix<> mu = gsMatrix<>::Random(1,numPts); // 1xnumPts Matrix filled with random numbers between (-1,1)
    mu = (mu + gsMatrix<>::Constant(1,numPts,1.))*urange/2.; // add umax to the matrix to have values between 0 and 2; multiply with range/2
    mu = (mu + gsMatrix<>::Constant(1,numPts,umin)); //set LO (=0) as the lower bound (offset)

    uv += noiseFactor * mu;

    for(auto row : uv.rowwise())
      std::sort(row.begin(), row.end());

    uv = (uv - gsMatrix<>::Constant(uv.rows(), uv.cols(), uv.minCoeff())) * (1/(uv.maxCoeff() - uv.minCoeff()));
    gsInfo << "Here is the scaled, sorted and PERTUBATED parameters:\n" << uv << "\n";
    }


    gsInfo << "Original basis:" << original.basis() << "\n";

    gsKnotVector<> kf(0, 1, numKnt, deg+1);//start,end,interior knots, start/end multiplicites of knots
    gsBSplineBasis<> basis(kf);
    gsInfo << "----------------------------------------------------\n";
    gsInfo << " + Original basis size = " << original.basis().size() << "\n";
    gsInfo << " + Fitting basis size = " << basis.size() << "\n";
    gsInfo << "----------------------------------------------------\n";
    gsWriteParaview( basis, "basis");
    GISMO_ENSURE(basis.size() < numPts, "Too few sampled points.");
    gsFitting<> fitted_curve(uv, X, basis);
    fitted_curve.compute(lambda);
    gsBSpline<>curve(basis, give(fitted_curve.result()->coefs()));
    coefs_plot = curve.coefs().transpose();
    gsWriteParaviewPoints(coefs_plot, "coefs_pdm");
    if(lambda > 0)
      gsInfo << "Least Squares with smoothing:\n " << curve << "\n";
    else
      gsInfo << "Ordinary Least Squares:\n " << curve << "\n";
    gsWriteParaview( curve, "pdm_curve", 1000);
    gsWriteParaview( curve, "pdm_cnet", 1000, false, true);

    //Apply parameter correction
    for (index_t i = 1; i < X.cols()-1; ++i)
    {
      gsVector<> newParam;
      curve.closestPointTo(X.col(i).transpose(), newParam, 1e-10, true);
      uv.col(i) = newParam;
    }

    gsInfo << "TDM - method:\n";


    gsMapData<> md(NEED_VALUE|NEED_NORMAL);
    md.points = uv; // parametric values
    curve.computeMap(md); // gsGeometry->computeMap(md);
    gsInfo << "Succesful map computation.\n";
    gsMatrix<> values  = md.values[0]; // a 3 x u.cols() matrix
    gsMatrix<> normals = md.normals;   // a 3 x u.cols() matrix

    // computing the unit normal of the geometry at the current parameter
    gsMatrix<> N_x(uv.cols(), 1);
    gsMatrix<> N_y(uv.cols(), 1);


    for (index_t k = 0; k < uv.cols(); ++k)
    {
      normals.col(k).normalize();
      N_x(k,0) = normals(0,k);
      N_y(k,0) = normals(1,k);
      GISMO_ENSURE(math::abs(normals.col(k).norm() - 1) < 1e-6, "no unit normal.");
    }



    gsSparseMatrix<> tmp(X.cols(), basis.size());
    gsMatrix<>  tmp_col(tmp.rows(), 1);
    tmp = basis.collocationMatrix( uv ) ;



    gsInfo << "Collocation matrix: " << tmp.rows() << " x " << tmp.cols() << "\n";
    // gsDebugVar(tmp);

    // gsMatrix<> A_mat(X.cols(), basis.size() * 2);


    gsMatrix<> B_mat(X.cols() * 2, basis.size() * 2);
    B_mat.setZero();
    B_mat.block(0,0,tmp.rows(),tmp.cols()) = tmp;
    B_mat.block(tmp.rows(),tmp.cols(),tmp.rows(),tmp.cols()) = tmp;
    // gsInfo << "~B:\n" << B_mat << "\n";

    gsMatrix<> N_diag(X.cols() * 2, X.cols());
    N_diag.setZero();
    for(index_t i = 0; i < X.cols(); i++)
    {
      N_diag(i,i) = N_x(i,0);
      N_diag(X.cols()+i,i) = N_y(i,0);
    }
    // N_diag(0,0) = 1;
    // N_diag(X.cols()-1, X.cols()-1) = 1;
    // N_diag(X.cols(), X.cols()-1) = 1;
    // N_diag(2*X.cols()-1, X.cols()-1) = 1;

    // gsInfo << "N_diag:\n" << N_diag << "\n";

    gsMatrix<> X_tilde(X.cols() * 2, 1);
    X_tilde << X.row(0).transpose(), X.row(1).transpose();

    gsMatrix<> A_tilde = N_diag.transpose() * B_mat;
    gsMatrix<> rhs = N_diag.transpose() * X_tilde;

    gsEigen::FullPivLU<gsMatrix<>::Base> lu_decomp(A_tilde);
    auto rank_tilde = lu_decomp.rank();
    gsInfo << "tilde size = " << A_tilde.rows() << " x " << A_tilde.cols() << "\n";
    gsInfo << "tilde rank = " << rank_tilde << "\n";

    // impose boundary conditions
    // x-component
    A_tilde.row(0) = gsVector<>::Zero(A_tilde.cols(),1);
    A_tilde(0,0) = 1;

    A_tilde.row(basis.size()-1) = gsVector<>::Zero(A_tilde.cols(),1);
    A_tilde(basis.size()-1, basis.size()-1) = 1;
    //
    // y-component
    A_tilde.row(basis.size()) = gsVector<>::Zero(A_tilde.cols(),1);
    A_tilde(basis.size(), basis.size()) = 1;

    A_tilde.row(2*basis.size()-1) = gsVector<>::Zero(A_tilde.cols(),1);
    A_tilde(2*basis.size()-1, 2*basis.size()-1) = 1;

    // X_tilde(X.cols() * 2, 1);
    rhs(0,0) = X_tilde(0,0); // x-component
    rhs(basis.size()-1,0) = X_tilde(X.cols()-1,0); // x-component

    rhs(basis.size(),0) = X_tilde(X.cols(),0); // y-component
    rhs(2*basis.size()-1,0) = X_tilde(2*X.cols()-1,0); // y-component

    gsSparseSolver<>::QR qrsolver_tilde(A_tilde.sparseView());
    gsMatrix<> sol_tilde = qrsolver_tilde.solve(rhs);

    gsInfo << "Tilde solution:\n";
    gsInfo << sol_tilde << "\n";

    gsMatrix<> coefs_tilde(basis.size(), 2);
    for(index_t j=0; j<basis.size(); j++)
    {
      coefs_tilde(j,0) = sol_tilde(j);
      coefs_tilde(j,1) = sol_tilde(basis.size()+j);
    }

    // check the system matrix
    // gsEigen::JacobiSVD<gsMatrix<>::Base> svd(A_tilde);
    // real_t cond = svd.singularValues()(0)/ svd.singularValues()(svd.singularValues().size()-1);

    gsBSpline<> curve_tilde(basis, coefs_tilde);
    gsWriteParaview(curve_tilde, "tdm_curve", 1000);
    gsWriteParaview( original, "tdm_cnet", 1000, false, true);

}
