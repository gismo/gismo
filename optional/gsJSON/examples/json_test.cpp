/** @file

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/

#include <gismo.h>
#include <gsJSON/gsJSON.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    std::string fileName;

    gsCmdLine cmd("Tutorial for reading and writing JSON files in G+Smo.");
    cmd.addString("f","file", "Input/output file name", fileName);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsInfo<<"=========================================================\n";
    gsInfo<<"Creating a JSON file...\n";

    // Create JSON writer
    gsJSON j;

    // Write integer
    j["a"] = 2;
    gsInfo<<"* Wrote an integer.\n";
    // Write double
    j["b"] = 3.0;
    gsInfo<<"* Wrote a double.\n";
    // Write string in a nested entry
    j["c"]["d"] = "four";
    gsInfo<<"* Wrote a string in a nested entry.\n";
    // Write an array of strings
    j["e"] = {"e1", "e2", "e3"};
    gsInfo<<"* Wrote an array of strings.\n";

    // Write a matrix
    gsMatrix<> mat(2,2);
    mat << 1, 2, 3, 4;
    j["A"] = mat;
    gsInfo<<"* Wrote a matrix.\n";

    // Write a vector
    gsVector<index_t> vec(3);
    vec << 1, 2, 3;
    j["v"] = vec;
    gsInfo<<"* Wrote a vector.\n";

    // Write a knot vector
    gsKnotVector<> kv(0,1,0,2);
    j["kv"] = kv;
    gsInfo<<"* Wrote a knot vector.\n";

    // Write a BSplineBasis
    gsBSplineBasis<> basis(kv);
    j["basis"] = basis;
    gsInfo<<"* Wrote a BSplineBasis.\n";

    // Write a TensorBSplineBasis
    gsTensorBSplineBasis<2> tbasis(kv,kv);
    j["tbasis"] = tbasis;
    gsInfo<<"* Wrote a TensorBSplineBasis.\n";

    // Write a TensorBSpline
    gsMatrix<> coefs = gsMatrix<>::Random(tbasis.size(),2);
    gsTensorBSpline<2> tgeom(tbasis, coefs);
    j["tgeom"] = tgeom;
    gsInfo<<"* Wrote a TensorBSpline.\n";

    // Print the file
    gsInfo<<"The contents of the JSON file are:\n"<<j<<"\n";

    // Save the file
    if (!fileName.empty())
    {
        gsInfo<<"Saving JSON file to "<<fileName<<"\n";
        j.save(fileName);
    }

    gsInfo<<"=========================================================\n";
    gsInfo<<"Construction of a JSON file from an option list.\n";
    // Define a gsOptionList
    gsOptionList opt;
    opt.addInt("a", "", 2);
    opt.addReal("b", "", 3.0);
    gsInfo<<"The option list is:\n"<<opt<<"\n";

    // Construct a JSON file from an option list
    // gsJSON j2(opt);
    j = gsJSON(opt); // This is the same as the above line, but using the gsJSON constructor directly
    // Print the file
    gsInfo<<"The contents of the JSON file constructed from an option list are:\n"<<j<<"\n";
    // Obtain the option list from the JSON object
    opt = j.get<gsOptionList>();
    gsInfo<<"The option list obtained from the JSON object is:\n"<<opt<<"\n";

    // Read the file written before
    if (!fileName.empty())
    {
        gsInfo<<"=========================================================\n";
        gsInfo<<"Reading the JSON file from "<<fileName<<"\n";
        j = gsJSON(fileName);
        // Print the knot vector
        kv = j["kv"].get<gsKnotVector<>>();
        gsInfo<<"* The knot vector read from the file is:\n"<<kv<<"\n";
        // Print the basis
        basis = j["basis"].get<gsBSplineBasis<>>();
        gsInfo<<"* The BSplineBasis read from the file is:\n"<<basis<<"\n";
        // Print the tensor basis
        tbasis = j["tbasis"].get<gsTensorBSplineBasis<2>>();
        gsInfo<<"* The TensorBSplineBasis read from the file is:\n"<<tbasis<<"\n";
        // Print the tensor geometry
        tgeom = j["tgeom"].get<gsTensorBSpline<2>>();
        gsInfo<<"* The TensorBSpline read from the file is:\n"<<tgeom<<"\n";
    }

    return EXIT_SUCCESS;

}//main
