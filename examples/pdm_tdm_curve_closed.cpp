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
    index_t numPts  = 20;
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

    gsKnotVector<> kv;
    kv.initUniform(0 , 1, 5, 1, 1, 3); //start,end,num_interior knots, interior multiplicity, end multiplicity, degree
    gsInfo << "Cyclic knot vector:\n" << kv << "\n";
    gsMatrix<> coefs(9, 2);

    coefs << 0, 0.25,
             0., 0.75,
             0.25, 1,
             0.5, 0.75,
             0.5, 0.25,
             0.25,0,
             0, 0.25,
             0., 0.75,
             0.25, 1;


    gsMatrix<> coefs_plot;
    coefs_plot = coefs.transpose();
    gsWriteParaviewPoints(coefs_plot, "coefs_original");
    gsBSpline<> original( kv, give(coefs));

    gsInfo << "Original curve:\n" << original << "\n";

    gsWriteParaview( original, "originalPlanarBsplineCurve", 1000);
    gsWriteParaview( original, "originalcnet", 1000, false, true);


    gsMatrix<> X, parameters, uv;

//
    real_t umin = 0.;
    real_t umax = 1.;

    real_t urange = umax - umin; // umax (=1) - umin (=0)
    gsMatrix<> mu = gsMatrix<>::Random(1,numPts); // 1xnumPts Matrix filled with random numbers between (-1,1)
    mu = (mu + gsMatrix<>::Constant(1,numPts,1.))*urange/2.; // add umax to the matrix to have values between 0 and 2; multiply with range/2
    mu = (mu + gsMatrix<>::Constant(1,numPts,umin)); //set LO (=0) as the lower bound (offset)


//
    gsMatrix<> uline = gsVector<>::LinSpaced(numPts, 0, 1);
    parameters = uline.transpose();
    original.eval_into(parameters, X);

    uv = parameters + noiseFactor * mu;

    for(auto row : uv.rowwise())
      std::sort(row.begin(), row.end());

    uv = (uv - gsMatrix<>::Constant(uv.rows(), uv.cols(), uv.minCoeff())) * (1/(uv.maxCoeff() - uv.minCoeff()));


    gsInfo << "Points:\n" << X << "\n";
    gsInfo << "Parameters:\n" << parameters << "\n";
    if(noiseFactor > 0)
      gsInfo << "Here is the scaled, sorted and PERTUBATED parameters:\n" << uv << "\n";

    gsWriteParaviewPoints(X, "points");


    gsInfo << "Original basis:" << original.basis() << "\n";

    gsKnotVector<> kf;
    kf.initUniform(0 , 1, numKnt, 1, 1, deg); //start,end,num_interior knots, interior multiplicity, end multiplicity, degree
    gsBSplineBasis<> basis(kf);
    gsInfo << "----------------------------------------------------\n";
    gsInfo << " + Original basis size = " << original.basis().size() << "\n";
    gsInfo << " + Fitting basis size = " << basis.size() << "\n";
    gsInfo << "----------------------------------------------------\n";
    gsWriteParaview( basis, "basis");
    GISMO_ENSURE(basis.size() < numPts, "Too few sampled points.");
    gsFitting<> fitted_curve(uv, X, basis);
    fitted_curve.compute(lambda);


    gsMatrix<> coefs_pdm;
    coefs_pdm = fitted_curve.result()->coefs();

    gsInfo << "Cyclic mapping:\n";
    for(index_t j=0; j<deg; j++)
    {
      gsInfo << j << "___" << basis.size()-deg+j << "\n";
      coefs_pdm.row(basis.size()-deg+j) = coefs_pdm.row(j);
    }

    gsBSpline<>curve(basis, coefs_pdm);
    coefs_plot = curve.coefs().transpose();
    gsWriteParaviewPoints(coefs_plot, "coefs_pdm");
    if(lambda > 0)
      gsInfo << "Least Squares with smoothing:\n " << curve << "\n";
    else
      gsInfo << "Ordinary Least Squares:\n " << curve << "\n";
    gsWriteParaview( curve, "pdm_curve", 1000);
    gsWriteParaview( curve, "pdm_cnet", 1000, false, true);

    // Apply parameter correction
    // for (index_t i = 1; i < X.cols()-1; ++i)
    // {
    //   gsVector<> newParam;
    //   curve.closestPointTo(X.col(i).transpose(), newParam, 1e-10, true);
    //   uv.col(i) = newParam;
    // }

    gsInfo << "TDM - method:\n";


    gsMapData<> md(NEED_VALUE|NEED_NORMAL|NEED_OUTER_NORMAL);
    md.points = uv; // parametric values
    curve.computeMap(md); // gsGeometry->computeMap(md);
    gsInfo << "Succesful map computation.\n";
    gsMatrix<> values  = md.values[0]; // a 3 x u.cols() matrix
    gsMatrix<> normals = md.normals;   // a 3 x u.cols() matrix
    //gsMatrix<> outers = md.outNormals;   // a 3 x u.cols() matrixs

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

    gsMatrix<> A_mat(X.cols(), basis.size() * 2);
    gsMatrix<> m_B(X.cols(), 1);
    m_B = X.row(0).transpose().cwiseProduct(N_x) + X.row(1).transpose().cwiseProduct(N_y);

    gsInfo << "m_B : " << m_B.rows() << " x " << m_B.cols() << "\n";

    index_t k2 = 0;


    for(index_t j=0; j < basis.size(); j++)
    {
      tmp_col = tmp.col(j);

      A_mat.col(k2) = tmp_col.cwiseProduct(N_x);
      k2 = k2+1;

      A_mat.col(k2) = tmp_col.cwiseProduct(N_y);
      k2 = k2 +1;
    }





    gsInfo << "Computing solution with QR factorization: A * c = b.\n";
    // gsEigen::FullPivHouseholderQR<gsMatrix<>::Base> qrsolver( A_mat );
    // gsDebugVar(A_mat);
    // sol_qr = qrsolver.solve(m_B);
    gsSparseSolver<>::QR qrsolver_tmp(A_mat.sparseView());
    gsMatrix<> sol_qr = qrsolver_tmp.solve(m_B);

    gsMatrix<> coefs_qr(basis.size(), 2);
    for(index_t j=0; j<basis.size(); j++)
    {
      coefs_qr(j,0) = sol_qr(2*j);
      coefs_qr(j,1) = sol_qr(2*j+1);
    }

    gsInfo << "Cyclic mapping:\n";
    for(index_t j=0; j<deg; j++)
    {
      gsInfo << j << "___" << basis.size()-deg+j << "\n";
      coefs_qr.row(j) = coefs_qr.row(basis.size()-deg+j);
    }

    gsInfo << coefs_qr << "\n";
    gsBSpline<> curve_qr( basis, coefs_qr);
    gsDebugVar(curve_qr);
    coefs_plot = coefs_qr.transpose();
    gsWriteParaviewPoints(coefs_plot, "coefs_qr");
    gsWriteParaview( curve_qr, "tdm_curve_qr");
    gsWriteParaview( curve_qr, "tdm_cnet_qr", 1000, false, true);

    gsInfo << "Computing solution with normal equations: A^tA * c = A^T b.\n";
    // sol_aa = (A_mat.transpose() * A_mat).ldlt().solve(A_mat.transpose() * m_B);
    gsMatrix<> AA = A_mat.transpose() * A_mat;
    gsSparseSolver<>::SimplicialLDLT LDLTsolver(AA.sparseView());
    gsMatrix<> sol_aa = LDLTsolver.solve(A_mat.transpose() * m_B);




    gsMatrix<> coefs_aa(basis.size(), 2);
    for(index_t j=0; j<basis.size(); j++)
    {
      coefs_aa(j,0) = sol_aa(2*j);
      coefs_aa(j,1) = sol_aa(2*j+1);
    }

    gsInfo << "Cyclic mapping:\n";
    for(index_t j=0; j<deg; j++)
    {
      gsInfo << j << "___" << basis.size()-deg+j << "\n";
      coefs_aa.row(basis.size()-deg+j) = coefs_aa.row(j);
    }

    gsInfo << coefs_aa << "\n";

    gsBSpline<> curve_aa( basis, coefs_aa);
    coefs_plot = coefs_aa.transpose();
    gsWriteParaviewPoints(coefs_plot, "coefs_aa");
    gsWriteParaview( curve_aa, "tdm_curve_aa");
    gsWriteParaview( curve_aa, "tdm_cnet_aa", 1000, false, true);


}
