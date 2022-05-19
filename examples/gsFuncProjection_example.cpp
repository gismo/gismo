/** @file biharmonic2_example.cpp

    @brief Tutorial on how to use expression assembler and the (approx.) C1 basis function
                to solve the Biharmonic equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

//! [Include namespace]
#include <gismo.h>

#include <gsUtils/gsFuncProjection.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    //! [Parse command line]
    bool plot  = false;
    bool mesh  = false;


    index_t numRefine  = 3;
    index_t degree = 3;
    index_t smoothness = 2;

    index_t gluingDataDegree = -1;
    index_t gluingDataSmoothness = -1;

    std::string output;
    std::string fn;

    gsCmdLine cmd(".");
    // Flags related to the input/geometry
    cmd.addString( "f", "file", "Input geometry file from path (with .xml)", fn );

    // Flags related to the discrete settings
    cmd.addInt( "p", "degree", "Set the polynomial degree of the basis.", degree );
    cmd.addInt( "s", "smoothness", "Set the smoothness of the basis.",  smoothness );
    cmd.addInt( "r", "numRefine", "Number of refinement-loops.",  numRefine );

    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);
    cmd.addSwitch("mesh", "Plot the mesh", mesh);

    cmd.addString("o", "output", "Output in xml (for python)", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    //! [Initialize data]
    gsMultiPatch<real_t> mp;

    //! [Read Argument inputs]
    //! [Read geometry]
    std::string string_geo;
    string_geo = fn;

    gsInfo << "Filedata: " << string_geo << "\n";
    gsReadFile<>(string_geo, mp);
    mp.clearTopology();
    mp.computeTopology();
    gsMultiBasis<real_t> dbasis(mp, true);//true: poly-splines (not NURBS)

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    dbasis.setDegree( degree); // preserve smoothness
    //dbasis.degreeElevate(degree- mp.patch(0).degree(0));

    if (degree-mp.patch(0).degree(0) >0)
        mp.degreeElevate(degree-mp.patch(0).degree(0));
    else if (degree-mp.patch(0).degree(0) <0)
    {
        mp.degreeReduce(mp.patch(0).degree(0)-degree);

        // // Perform
        // gsMultiPatch<> result;
        // gsFuncProjection<real_t>::L2(dbasis,mp,result);
        // gsWriteParaview(result,"mp_corrected",1000,true);
        // gsWrite(result,"mp_corrected");
        // mp = result;
        // gsInfo<<"\nProjected to lower-order basis of degree "<<degree<<"\n\n";
    }

    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        dbasis.uniformRefine(1, degree-smoothness);

    //! [Export visualization in ParaView]
    if (plot)
    {
        gsInfo<<"Plotting in Paraview...\n";
        gsWriteParaview( mp, "geom",1000,true);
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";
    //! [Export visualization in ParaView]

    // //! [Export data to xml]
    // if (!output.empty())
    // {
    //     index_t cols = method == MethodFlags::NITSCHE ? 7+penalty.cols() : 7;
    //     gsMatrix<real_t> error_collection(l2err.rows(), cols);
    //     error_collection.col(0) = meshsize;
    //     error_collection.col(1) = dofs;
    //     error_collection.col(2) = l2err;
    //     error_collection.col(3) = h1err;
    //     error_collection.col(4) = h2err;
    //     error_collection.col(5) = IFaceErr;
    //     error_collection.col(6) = cond_num;
    //     if (method == MethodFlags::NITSCHE)
    //         error_collection.block(0,7,penalty.rows(),penalty.cols()) = penalty;

    //     gsFileData<real_t> xml_out;
    //     xml_out << error_collection;
    //     xml_out.addString("Meshsize, dofs, l2err, h1err, h2err, iFaceErr, cond_num, (penalty)","Label");
    //     xml_out.addString(std::to_string(degree),"Degree");
    //     xml_out.addString(std::to_string(smoothness),"Regularity");
    //     xml_out.addString(std::to_string(numRefine),"NumRefine");
    //     xml_out.addString(std::to_string(method),"Method");
    //     // Add solution
    //     // [...]
    //     xml_out.save(output);
    //     gsInfo << "XML saved to " + output << "\n";
    // }
    // //! [Export data to xml]

    // if (!write.empty())
    // {
    //     std::ofstream file(write.c_str());
    //     file<<"Meshsize, dofs, l2err, h1err, h2err, iFaceErr"<<"\n";
    //     for (index_t k=0; k<meshsize.size(); ++k)
    //     {
    //         file<<meshsize[k]<<","<<dofs[k]<<","<<l2err[k]<<","<<h1err[k]<<","<<h2err[k]<<","<<IFaceErr[k]<<"\n";
    //     }
    //     file.close();
    // }


    return EXIT_SUCCESS;
}// end main
