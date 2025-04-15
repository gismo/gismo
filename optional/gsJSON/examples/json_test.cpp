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
    // Create JSON writer
    gsJSON j;

    // Write integer
    j["a"] = 2;
    // Write double
    j["b"] = 3.0;
    // Write string in a nested entry
    j["c"]["d"] = "four";
    // Write an array of strings
    j["e"] = {"e1", "e2", "e3"};

    // Write a matrix
    gsMatrix<> mat(2,2);
    mat << 1, 2, 3, 4;
    j["A"] = mat;

    // Write a vector
    gsVector<index_t> vec(3);
    vec << 1, 2, 3;
    j["v"] = vec;

    // Write a knot vector
    gsKnotVector<> kv(0,1,0,2);
    j["kv"] = kv;

    // Write a BSplineBasis
    gsBSplineBasis<> basis(kv);
    j["basis"] = basis;

    // Write a TensorBSplineBasis
    gsTensorBSplineBasis<2> tbasis(kv,kv);
    j["tbasis"] = tbasis;

    // Write a TensorBSpline
    gsMatrix<> coefs = gsMatrix<>::Random(tbasis.size(),2);
    gsTensorBSpline<2> tgeom(tbasis, coefs);
    j["tgeom"] = tgeom;

    // Print the file
    gsInfo<<j<<"\n";

    // Save the file
    j.save("test.json");

    // Define a gsOptionList
    gsOptionList opt;
    opt.addInt("a", "", 2);
    opt.addReal("b", "", 3.0);

    // Construct a JSON file from an option list
    gsJSON j2(opt);
    // Print the file
    gsInfo<<j2<<"\n";
    // Obtain the option list from the JSON object
    gsOptionList opt2 = j2.get<gsOptionList>();
    gsInfo<<opt2;

    // Read the file written before
    gsJSON j3("test.json");
    // Print the knot vector
    gsKnotVector<> kv2 = j3["kv"].get<gsKnotVector<>>();
    gsInfo<<kv2<<"\n";
    // Print the basis
    gsBSplineBasis<> basis2 = j3["basis"].get<gsBSplineBasis<>>();
    gsInfo<<basis2<<"\n";
    // Print the tensor basis
    gsTensorBSplineBasis<2> tbasis2 = j3["tbasis"].get<gsTensorBSplineBasis<2>>();
    gsInfo<<tbasis2<<"\n";
    // Print the tensor geometry
    gsTensorBSpline<2> tgeom2 = j3["tgeom"].get<gsTensorBSpline<2>>();
    gsInfo<<tgeom2<<"\n";

    return EXIT_SUCCESS;

}//main
