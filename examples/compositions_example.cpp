/** @file compositions_example.cpp

    @brief Tutorial on how to generate a composed basis and a composed geometry

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
    // Construct the composition
    gsGeometry<>::uPtr composition = gsNurbsCreator<>::BSplineSquareDeg(2);
    composition->coefs().row(4).setConstant(0.9);
    gsInfo<<"The composition has basis:\n"<<composition->basis()<<"\n";
    gsInfo<<"The control points of the map are:\n"<<composition->coefs()<<"\n";

    // Plot the composition
    gsWriteParaview(*composition,"composition");

    // Construct a basis
    gsKnotVector<> kv1({0,0,0,1./3.,2./3.,1,1,1},2);
    gsKnotVector<> kv2({0,0,0,0.25,0.50,0.75,1,1,1},2);
    gsTensorBSplineBasis<2> tbasis(kv1,kv2);

    // Composte it
    gsComposedBasis<> cbasis(*composition,tbasis);

    // Plot the basis
    gsWriteParaview(cbasis,"basis",1000);
    // .. and its mesh
    gsMesh<> mesh(tbasis);
    cbasis.mapMesh(mesh);
    gsWriteParaview(mesh,"mesh");

    // Construct a random geometry
    gsMatrix<> coefs(cbasis.size(),3);
    coefs.leftCols(2) = cbasis.anchors().transpose();
    coefs.col(2).setRandom();

    // Make the geometries (composed and non-composed)
    gsGeometry<>::uPtr geom = tbasis.makeGeometry(coefs);
    gsGeometry<>::uPtr cgeom= cbasis.makeGeometry(coefs);

    // Plot the geometries (composed and non-composed)
    gsWriteParaview(*geom,"geom",1000,true);
    gsWriteParaview(*cgeom,"cgeom",1000,true);


    return EXIT_SUCCESS;

}// end main
