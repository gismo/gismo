/** @file lofting_example.cpp

    @brief Demonstrates adaptive fitting of data samples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <gismo.h>
#include <gsModeling/gsLofting.h>

using namespace gismo;

gsMultiPatch<> getMPfromGeom(gsFileData<> & data);

int main(int argc, char *argv[])
{
    // Options with default values
    
    std::string fn = "curves3d/loft_exm0.xml";
    real_t lambda = 1e-07;
    index_t numSamples = 1000;
    index_t deg_v = 2;
    real_t cm = 1.;

    // Reading options from the command line
    gsCmdLine cmd("Loft parametrized B-spline cuves. Expected input file is an XML "
            "file containing the curves to be lofted.");
    
    cmd.addString("f", "filedata", "Input sample data", fn);
    cmd.addReal("l", "lambda", "smoothing coefficient", lambda);
    cmd.addInt("s", "samples", "number of plotting samples", numSamples);
    cmd.addInt("v", "vdeg", "degree in v direction", deg_v);
    cmd.addReal("m", "method", "isoparametric curves placement method: (0) unfirom; (0.5) centripetal; (1) chord-length", cm);
    

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    
    // typename std::vector<gsBSpline<real_t>> container;
    gsMultiPatch<real_t> container;
    //! [Read data]

    gsReadFile<>(fn, container);

    gsInfo << container << "\n";
    gsInfo << "number of lofting curves: " << container.nPatches() << "\n";

    gsWriteParaview(container, "section_curves", numSamples);
    
    // choose the degree, parameters and knot-vector in the v-direction.
    // degree must be less than #section curves


    // Create knot-vectors without interior knots
    // gsKnotVector<> u_knots (u_min, u_max, 0, deg_x+1 ) ;
    index_t interiors = container.nPatches() - (deg_v + 1);
    gsKnotVector<> kv (0., 1., interiors, deg_v+1) ;


    // Create lofting object
    // gsLofting<real_t> loft(container);
    // gsLofting<real_t> loft(container, deg_v);
    gsLofting<real_t> loft(container, kv);
    loft.compute();

    // This is for outputing an XML file, if requested
    gsFileData<> fd;
    fd << *loft.result() ;

    fd.dump("lofting_out");


    return 0;
}

gsMultiPatch<> getMPfromGeom(gsFileData<> & data)
{
    gsMultiPatch<> container;
    std::vector< memory::unique_ptr< gsGeometry<> > >  geoContainer = data.getAll< gsGeometry<> >();

    for(index_t el = 0; el < geoContainer.size(); el++)
        container.addPatch(*geoContainer[el]);

    return container;
}