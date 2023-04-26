/** @file gsScatterGeometry.cpp

    @brief Take a gsGeometry and make a point-cloud out of it.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <gismo.h>
#include <string>
#include <iostream>
#include <Eigen/Dense>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values
    bool save     = false;
    std::string in_name("in_geometry");
    index_t numSamples(10000);
    int numPts = 1000;
    int numPatch = 0;
    real_t umin = 0.;
    real_t umax = 0.;
    real_t vmin = 0.;
    real_t vmax = 0.;
    bool more     = false;
    int morePts = 100;
    bool plot_mesh = true;
    bool plot_net = false;
    std::string fn = "fitting/deepdrawingC.xml";

    // Reading options from the command line
    gsCmdLine cmd("Sample a geometry given in input (or one of its patch) to generate an unorganized point-cloud."
                   "Expected input file is an XML file containing the geometry."
                   "Return in output two matrices (<Matrix>), with \nMatrix id 0 : contains a 2 x N matrix. "
                   "Every column represents a (u,v) parametric coordinate\nMatrix id 1 : contains a "
                   "3 x N matrix. Every column represents a point (x,y,z) in space.");
    cmd.addSwitch("save", "Save result in XML format", save);
    cmd.addString("t", "data", "Input sample data", fn);
    cmd.addInt("m", "psize", "number of sampled points", numPts);
    cmd.addInt("p", "patch", "number of the patch to be gridded", numPatch);
    cmd.addReal("a", "u-min", "Select manually the parametric domain [a,b]x[c,d]: a value.", umin);
    cmd.addReal("b", "u-max", "Select manually the parametric domain [a,b]x[c,d]: b value.", umax);
    cmd.addReal("c", "v-min", "Select manually the parametric domain [a,b]x[c,d]: c value.", vmin);
    cmd.addReal("d", "v-max", "Select manually the parametric domain [a,b]x[c,d]: d value.", vmax);
    cmd.addSwitch("more", "Additional points along parametric isolines.", more);
    cmd.addInt("r", "rsize", "number of additional sampled points", morePts);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //! [Read data]
    gsFileData<> fd_in(fn);

    std::vector<gsGeometry<>::uPtr> geo = fd_in.getAll< gsGeometry<> >();
    if ( ! geo.empty() ){
        gsInfo<< "Got "<< geo.size() <<" patch"<<(geo.size() == 1 ? "." : "es.") <<"\n";
    gsWriteParaview(memory::get_raw(geo), in_name, numSamples, plot_mesh, plot_net);
    }
    else{
        gsInfo<< "----------------------------------------------" <<"\n";
        gsInfo<< "No geometry has been found in input, quitting." <<"\n";
        gsInfo<< "----------------------------------------------" <<"\n";
        return 0;
    }

    if ( numPatch > geo.size()-1 ){
        gsInfo<< "--------------------------------------------------------------------" <<"\n";
        gsInfo<< "Patch number "<< numPatch <<" has not been found in input, quitting." <<"\n";
        gsInfo<< "--------------------------------------------------------------------" <<"\n";
        return 0;
    }

    gsMatrix<> ab(2,2);
    if (umin + umax + vmin + vmax > 0.)
    {
        ab << umin, umax, vmin, vmax;
        gsInfo << "Selected domain diagonal corners:\n" << ab << "\n";
    }
    else
    {
        ab = geo[numPatch]->support();
        gsInfo << "Selected support diagonal corners:\n" << ab << "\n";
    }

    real_t totPts;
    if(more){
      totPts = numPts + 5 * morePts;
    }
    else{
      totPts = numPts;
    }

    gsVector<unsigned> numPtsVec(2);
    numPtsVec<<numPts,numPts;
    gsVector<> a = ab.col(0);
    gsVector<> b = ab.col(1);

    //real_t upper = a(0); // umin
    //real_t lower = a(1); // umax
    real_t urange= b(0)-a(0); // umax - umin
    gsInfo << "u-range sampling: " << urange << "\n";
    gsMatrix<> mu = gsMatrix<>::Random(1,numPts); // 3x3 Matrix filled with random numbers between (-1,1)
    mu = (mu + gsMatrix<>::Constant(1,numPts,1))*urange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
    mu = (mu + gsMatrix<>::Constant(1,numPts,a(0))); //set LO as the lower bound (offset)
    //gsInfo << "mu =\n" << mu << "\n";

    //real_t upper = b(0); // vmin
    //real_t lower = b(1); // vmax
    real_t vrange= b(1)-a(1); // vmax - vmin
    gsInfo << "v-range sampling: " << vrange << "\n";
    gsMatrix<> mv= gsMatrix<>::Random(1,numPts); // 3x3 Matrix filled with random numbers between (-1,1)
    mv = (mv + gsMatrix<>::Constant(1,numPts,1))*vrange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
    mv = (mv + gsMatrix<>::Constant(1,numPts,a(1))); //set LO as the lower bound (offset)
    //gsInfo << "mv =\n" << mv << "\n";


    gsMatrix<> uv_interiors(2, numPts);
    uv_interiors << mu, mv; //= gsPointGrid(a,b, numPtsVec);
    //gsInfo << "Random uniform interior parameters:\n" << uv_interiors.rows() << " x " << uv_interiors.cols() << "\n";
    gsWriteParaviewPoints(uv_interiors, "interior_parameters");

    int numBts = math::ceil(math::sqrt(totPts)) + 2; // b = math.ceil(math.sqrt(n)) + 2
    //gsInfo << "Number of boundary points: " << numBts*4-4 << "\n";

    // Sample the boundaries and the corners.
    gsMatrix<> uv_boundary(2, numBts*4-4);
    gsMatrix<> b_0(1, numBts-1);
    gsMatrix<> b_1(1, numBts-1);
    gsMatrix<> b_2(1, numBts-1);
    gsMatrix<> b_3(1, numBts-1);
    for (index_t pace=0; pace < 4; pace++){
      if (pace == 0){
        gsMatrix<> mu = gsMatrix<>::Random(numBts-2, 1); // 1xnumBts-1 Matrix filled with random numbers between (-1,1)
        mu = (mu + gsMatrix<>::Constant(numBts-2, 1, 1))*urange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
        mu = (mu + gsMatrix<>::Constant(numBts-2,1,a(0))); //set LO as the lower bound (offset)
        mu.sortByColumn(0);
        mu = mu.reshape(1, numBts-2);
        b_0 << a(0), mu;
        //gsInfo << pace << "-boundary paramters: " << b_0.rows() << "x" << b_0.cols() << "\n" << b_0 << "\n";
      }
      if (pace == 1){
        gsMatrix<> mu = gsMatrix<>::Random(numBts-2, 1); // 1xnumBts-1 Matrix filled with random numbers between (-1,1)
        mu = (mu + gsMatrix<>::Constant(numBts-2, 1, 1))*vrange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
        mu = (mu + gsMatrix<>::Constant(numBts-2,1, a(1))); //set LO as the lower bound (offset)
        mu.sortByColumn(0);
        mu = mu.reshape(1, numBts-2);
        b_1 << a(1), mu;
        //gsInfo << pace << "-boundary paramters: " << b_1.rows() << "x" << b_1.cols() << "\n" << b_1 << "\n";
      }
      if (pace == 2){
        gsMatrix<> mu = gsMatrix<>::Random(numBts-2, 1); // 1xnumBts-1 Matrix filled with random numbers between (-1,1)
        mu = (mu + gsMatrix<>::Constant(numBts-2, 1, 1))*urange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
        mu = (mu + gsMatrix<>::Constant(numBts-2,1,a(0))); //set LO as the lower bound (offset)
        mu.sortByColumn(0);
        mu = mu.reshape(1, numBts-2);
        b_2 << b(0), mu.reverse();
        //gsInfo << pace << "-boundary paramters: " << b_2.rows() << "x" << b_2.cols() << "\n" << b_2 << "\n";
      }
      if (pace == 3){
        gsMatrix<> mu = gsMatrix<>::Random(numBts-2, 1); // 1xnumBts-1 Matrix filled with random numbers between (-1,1)
        mu = (mu + gsMatrix<>::Constant(numBts-2, 1, 1))*vrange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
        mu = (mu + gsMatrix<>::Constant(numBts-2,1,a(1))); //set LO as the lower bound (offset)
        mu.sortByColumn(0);
        mu = mu.reshape(1, numBts-2);
        b_3 << b(1), mu.reverse();
        //gsInfo << pace << "-boundary paramters: " << b_3.rows() << "x" << b_3.cols() << "\n" << b_3 << "\n";
      }
    }

    gsMatrix<> u_zeros = gsMatrix<>::Constant(1, numBts-1, a(1));
    gsMatrix<> v_zeros = gsMatrix<>::Constant(1, numBts-1, a(0));

    gsMatrix<> u_ones = gsMatrix<>::Constant(1, numBts-1, b(1));

    gsMatrix<> v_ones = gsMatrix<>::Constant(1, numBts-1, b(0));

    gsMatrix<> zeros = gsMatrix<>::Zero(1, numBts-1);
    gsMatrix<> ones  = gsMatrix<>::Ones(1, numBts-1);

    //gsInfo << "Zeros:\n" << zeros << "\n";
    //gsInfo << "Ones:\n" << ones << "\n";
    uv_boundary << b_0,     v_ones, b_2,    v_zeros,
                   u_zeros, b_1,    u_ones, b_3;
    gsWriteParaviewPoints(uv_boundary, "uv_boundary");

    gsInfo << "South:\n";
    gsInfo << b_0 << "\n" << u_zeros << "\n";
    gsInfo << "East:\n";
    gsInfo << v_ones << "\n" << b_1 << "\n";
    gsInfo << "North:\n";
    gsInfo << b_2 << "\n" << u_ones << "\n";
    gsInfo << "West:\n";
    gsInfo << v_zeros << "\n" << b_3 << "\n";


    //gsInfo << "Boundary points:\n" << uv_boundary << "\n";

    gsMatrix<> uv;
    if (more){
      // gsMatrix<> addPts;
      // gsMatrix<> uv_add(2, morePts);
      // gsMatrix<> uline = gsVector<>::LinSpaced(morePts, a(0), b(0));
      // //gsInfo << "uline:\n" << uline << "\n";
      // gsMatrix<> vline = gsMatrix<>::Random(morePts, 1); // 1xnumBts-1 Matrix filled with random numbers between (-1,1)
      // vline = (vline + gsMatrix<>::Constant(morePts, 1, 1.))*0.25/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
      // vline = (vline + gsMatrix<>::Constant(morePts,1,0.75)); //set LO as the lower bound (offset)
      //
      // uv_add.row(0) = uline.reshape(1,morePts);
      // uv_add.row(1) = vline.col(0);
      // geo[numPatch]->eval_into(uv_add, addPts);
      // gsWriteParaviewPoints(addPts, "addPoints");
      // gsWriteParaviewPoints(uv_add, "addParameters");

      gsInfo << "Filling the sharp features parameters." << "\n";
      gsMatrix<> removeVoids;
      gsMatrix<> uv_more(2, 5*morePts);
      //gsMatrix<> uv_more(2, 5*morePts);
      gsMatrix<> segment0(2, morePts);
      gsMatrix<> segment1(2, morePts);
      gsMatrix<> segment2(2, morePts);
      gsMatrix<> segment3(2, morePts);
      gsMatrix<> segment4(2, morePts);

      for(index_t segment=0; segment < 5; segment++){
        real_t v0 = 0;
        real_t v1 = 0.2;
        real_t u0 = 0.5;
        real_t u1 = 0.75;

        if (segment == 1){
          v0 = 0.2;
          v1 = 0.4;
        }
        if (segment == 2){
          v0 = 0.4;
          v1 = 0.6;
        }
        if (segment == 3){
          v0 = 0.6;
          v1 = 0.8;
        }
        if (segment == 4){
          v0 = 0.8;
          v1 = 1;
        }
        real_t urange = u1 - u0;
        gsMatrix<> temp_uv_more(2, morePts);
        gsMatrix<> vvoid = gsVector<>::LinSpaced(morePts, v0, v1);
        //gsInfo << "uline:\n" << uline << "\n";
        gsMatrix<> uvoid = gsMatrix<>::Random(morePts, 1); // 1xnumBts-1 Matrix filled with random numbers between (-1,1)
        uvoid = (uvoid + gsMatrix<>::Constant(morePts, 1, 1.))*urange/2.; // add 1 to the matrix to have values between 0 and 2; multiply with range/2
        uvoid = (uvoid + gsMatrix<>::Constant(morePts,1,u0)); //set LO as the lower bound (offset)
        temp_uv_more.row(0) = uvoid.reshape(1,morePts);
        temp_uv_more.row(1) = vvoid.col(0);

        gsInfo << "Filling the segments." << "\n";
        if (segment == 0){
          segment0 << temp_uv_more;
          gsInfo << "Segment 0:" << segment0 << "\n";
        }
        if (segment == 1){
          segment1 << temp_uv_more;
          gsInfo << "Segment 1\n";
        }
        if (segment == 2){
          segment2 << temp_uv_more;
          gsInfo << "Segment 2\n";
        }
        if (segment == 3){
          segment3 << temp_uv_more;
          gsInfo << "Segment 3\n";
        }
        if (segment == 4){
          segment4 << temp_uv_more;
          gsInfo << "Segment 4\n";
        }
      }
      gsInfo << "Filling the segments all together:\n" << "\n";
      uv_more << segment0, segment1, segment2, segment3, segment4;
      gsInfo << "Done." << "\n";




      geo[numPatch]->eval_into(uv_more, removeVoids);
      gsWriteParaviewPoints(removeVoids, "voidPoints");
      gsWriteParaviewPoints(uv_more, "voidParameters");

      gsMatrix<> x_I(2, numPts + 5 * morePts);
      x_I << uv_interiors, uv_more;
      Eigen::PermutationMatrix<Dynamic,Dynamic> perm(x_I.cols());
      perm.setIdentity();
      std::random_shuffle(perm.indices().data(), perm.indices().data()+perm.indices().size());
      gsMatrix<> A_perm = x_I * perm; // permute columns
      //A_perm = perm * A; // permute rows


      gsMatrix<> params(2, numPts + numBts*4-4 + 5*morePts);
      params << A_perm, uv_boundary;
      uv = params;
    }
    else{
      gsMatrix<> params(2, numPts + numBts*4-4);
      params << uv_interiors, uv_boundary;
      uv = params;
    }

    gsMatrix<> xyz;
    geo[numPatch]->eval_into(uv, xyz);
    gsWriteParaviewPoints(xyz, "points");


    gsMatrix<> uv_normalized = uv;
    uv_normalized.row(0) = (uv.row(0) - gsMatrix<>::Constant(1, uv.cols(), a(0))) / (b(0)-a(0)) ;
    uv_normalized.row(1) = (uv.row(1) - gsMatrix<>::Constant(1, uv.cols(), a(1))) / (b(1)-a(1)) ;

    gsWriteParaviewPoints(uv_normalized, "normalized_parameters");
    gsWriteParaviewPoints(uv, "parameters");

    gsInfo << "\nPoint cloud dimension:" << "\n";
    gsInfo << "         Total number of points: m = " << xyz.cols() << "\n";
    gsInfo << "      number of interior points: n = " << uv.cols() -  uv_boundary.cols() << "\n";
    gsInfo << "  number of boundary points: b*4-4 = " << uv_boundary.cols() << "\n";


    gsFileData<> fd;
    fd << xyz;
    fd << uv_normalized;

    fd.dump("floaterPts_out");

    for(index_t p=0; p != geo.size(); p++){
    gsMatrix<> pp;
    geo[p]->eval_into(uv, pp);
    gsWriteParaviewPoints(pp, "datapatch" + internal::to_string(p));
    }

    return 0;
}
