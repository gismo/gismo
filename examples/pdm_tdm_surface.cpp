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

using namespace gismo;


int main(int argc, char *argv[])
{
    real_t deg = 2;
    index_t numKnots = 1;
    real_t lambda = 0;
    std::string fn = "../filedata/fitting/floaterPts_out.xml";

    gsCmdLine cmd("Tensor product B-spline surface fitting with Tangent Distance Minimization.");

    cmd.addReal("d", "degree", "bi-degree (q,q).", deg);
    cmd.addReal("s", "smoothing", "smoothing weight", lambda);
    cmd.addInt("n", "interiors", "number of interior knots in each direction.", numKnots);
    cmd.addString("f", "filename", "name of the .xml file containing the data", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }


    gsStopwatch time;

    gsFileData<> fd_in(fn);
    gsMatrix<> uv, X;
    fd_in.getId<gsMatrix<> >(0, uv );
    fd_in.getId<gsMatrix<> >(1, X);
    //! [Read data]

    gsWriteParaviewPoints(uv, "parameters");
    gsWriteParaviewPoints(X, "points");

    GISMO_ENSURE( uv.cols() == X.cols() && uv.rows() == 2 && X.rows() == 3,
                  "Wrong input, check id of matrices in the .xml file");

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

    gsInfo << "Least square fitting, with smoothing:\n " << original << "\n";
    gsWriteParaview( original, "pdm_0", 100000);
    writeSingleControlNet(original, "pdm_0_cnet");

    gsMapData<> md(NEED_VALUE|NEED_NORMAL);
    md.points = uv; // parametric values

    original.computeMap(md); // gsGeometry->computeMap(md);
    gsInfo << "Succesful map computation.\n";
    gsMatrix<> values  = md.values[0]; // a 3 x u.cols() matrix
    gsMatrix<> normals = md.normals;   // a 3 x u.cols() matrix


    // computing the unit normal of the geometry at the current parameter
    gsMatrix<> N_x(uv.cols(), 1);
    gsMatrix<> N_y(uv.cols(), 1);
    gsMatrix<> N_z(uv.cols(), 1);

    for (index_t k = 0; k < uv.cols(); ++k)
    {
      normals.col(k).normalize();
      // gsInfo << normals.col(k).transpose() << "\n";
      N_x(k,0) = normals(0,k);
      N_y(k,0) = normals(1,k);
      N_z(k,0) = normals(2,k);
      // gsInfo << "nx: " << N_x(k,0) << ", ny: " << N_y(k,0) << ", nz: " << N_z(k,0) << "\n\n";
    }

    gsBSplineBasis<> lineBasis(0,1,0,1); // start, end, interior, degree, mult_int, periodic
    gsWriteParaview(lineBasis, "1dbasis");
    gsMatrix<> vlines = values + 0.05 * normals;   // a 3 x u.cols() matrix

    // gsMultiPatch<> mp_line;
    // for (index_t k = 0; k < uv.cols(); ++k)
    // {
    //   gsMatrix<> lineCoefs(2, 3);
    //   lineCoefs.row(0) = values.col(k);
    //   lineCoefs.row(1) = vlines.col(k);
    //   gsBSpline<> line(lineBasis, give(lineCoefs));
    //   gsMatrix<> mv(3, 1);
    //   mv = vlines.col(0);
    //   mp_line.addPatch(line);
    //   // gsWriteParaview(line, "line" + internal::to_string(k));
    // }
    // gsWriteParaview(mp_line, "allNormals");



//#   pragma omp parallel for default(shared) private(curr_point,actives,value)
    gsInfo << "assemble TDM system for x-component:\n";
    gsSparseMatrix<> Ax_mat(basis.size(), basis.size());
    gsSparseMatrix<> Ay_mat(basis.size(), basis.size());
    gsSparseMatrix<> Az_mat(basis.size(), basis.size());
    gsMatrix<> x_B(basis.size(), 1);
    gsMatrix<> y_B(basis.size(), 1);
    gsMatrix<> z_B(basis.size(), 1);

    gsSparseMatrix<> tmp(X.cols(), basis.size());
    tmp = basis.collocationMatrix( uv ) ;

    gsSparseMatrix<> X_mat(tmp.rows(), tmp.cols());
    gsSparseMatrix<> Y_mat(tmp.rows(), tmp.cols());
    gsSparseMatrix<> Z_mat(tmp.rows(), tmp.cols());
    for(index_t j = 0; j < tmp.cols(); j ++)
    {
      // gsInfo << "column:\n" << tmp.col(j) << "\n";
      // gsInfo << "vector:\n" << N_x << "\n";
      X_mat.col(j) = tmp.col(j).cwiseProduct(N_x + N_y + N_z);
      Y_mat.col(j) = tmp.col(j).cwiseProduct(N_x + N_y + N_z);
      Z_mat.col(j) = tmp.col(j).cwiseProduct(N_x + N_y + N_z);
    }

    gsMatrix<> b_x(N_x.rows(), N_x.cols());
    b_x = X.row(0).transpose().cwiseProduct(N_x + N_y + N_z);

    gsMatrix<> b_y(N_y.rows(), N_y.cols());
    b_y = X.row(1).transpose().cwiseProduct(N_x + N_y + N_z);

    gsMatrix<> b_z(N_z.rows(), N_z.cols());
    b_z = X.row(2).transpose().cwiseProduct(N_x + N_y + N_z);

    Ax_mat = X_mat.transpose() * X_mat;
    x_B = X_mat.transpose() * b_x;

    Ax_mat.makeCompressed();

    Eigen::LDLT<gsMatrix<>::Base> xsolver( Ax_mat );
    gsMatrix<> sol_x;
    sol_x = xsolver.solve(x_B);

    // gsInfo << "Initial x-coefs:\n" << original.coefs().col(0) << "\n";
    // gsInfo << "New x-coefs:\n" << sol_x << "\n";

    Ay_mat = Y_mat.transpose() * Y_mat;
    y_B = Y_mat.transpose() * b_y;

    Ay_mat.makeCompressed();

    Eigen::LDLT<gsMatrix<>::Base> ysolver( Ay_mat );
    gsMatrix<> sol_y;
    sol_y = ysolver.solve(y_B);

    // gsInfo << "Initial y-coefs:\n" << original.coefs().col(1) << "\n";
    // gsInfo << "New y-coefs:\n" << sol_y << "\n";

    Az_mat = Z_mat.transpose() * Z_mat;
    z_B = Z_mat.transpose() * b_z;

    Az_mat.makeCompressed();

    Eigen::LDLT<gsMatrix<>::Base> zsolver( Az_mat );
    gsMatrix<> sol_z;
    sol_z = zsolver.solve(z_B);

    gsInfo << "x-coefs diff --- y-coefs diff --- z-coefs diff\n";
    gsInfo << (coefs.col(0) - sol_x).norm() << "        " << (coefs.col(1) - sol_y).norm() << "        " << (coefs.col(2) - sol_z).norm() << "\n";

    gsMatrix<> sol(tmp.cols(), 3);
    sol << sol_x, sol_y, sol_z;
    // gsInfo << "TDM coefs:\n" << sol << "\n";
    gsTensorBSpline<2, real_t> tdm_surface(basis, sol);
    gsWriteParaview( tdm_surface, "tdm_0", 100000);
    writeSingleControlNet(tdm_surface, "tdm_0_cnet");

}
