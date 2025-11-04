/** @file paraview_example.cpp

    @brief Tutorial on how to export geometries and solution data to ParaView files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Karampatzakis
*/


#include <iostream>
#include <gismo.h>

using namespace gismo;


int main(int argc, char* argv[])
{

    //! [Parse command line]
    std::string input("surfaces/egg.xml");
    std::string output("");
    bool bezier=false; 

    gsCmdLine cmd("Tutorial on ParaView Output");
    cmd.addPlainString("filename", "G+Smo input geometry file.", input);
    cmd.addString("o", "output", "Base name of the output file", output);
    cmd.addSwitch("b","bezier", "Output using Bezier elements", bezier);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    if (!gsFileManager::fileExists(input))
    {
        gsWarn << "The file cannot be found!\n";
        return EXIT_FAILURE;
    }

    gsInfo << "Read file \"" << input << "\"\n";

    //! [Read geometry])
    gsFileData<> fileData(input);


    gsGeometry<>::uPtr pGeom;
    gsMultiPatch<> mPatch;
    if(!fileData.has< gsGeometry<> >())
    {
        gsWarn << "Input file doesn't contain any geometry.\n";
        return EXIT_FAILURE;
    }
    else if(fileData.has< gsMultiPatch<> >())
    {
        mPatch = * fileData.getFirst< gsMultiPatch<> >();
        gsInfo << "Found a ";
    }
    else
    {
        std::vector<gsGeometry<>::uPtr> geo = fileData.getAll< gsGeometry<> >();
        // Put all geometries in mPatch
        for (index_t g = 0; g < (index_t)geo.size(); ++g)
        {
            mPatch.addPatch(give(*geo[g]));
        }
        gsInfo << "Created a ";
    }
    gsInfo << mPatch << "\n";
    //! [Read geometry]

    if ( output.empty() )
    {
        gsInfo << "Call program with option -o <basename> to write data to files\n";
        return EXIT_SUCCESS;
    }

    //! [gsWriteParaview]
    gsMultiBasis<> mBasis(mPatch);
    gsWriteParaview(mPatch, output+"MultiPatch");
    gsWriteParaview(mPatch, mBasis, output+"MultiBasis");
    gsMatrix<> cPoints = mPatch.coefs().transpose();
    gsWriteParaviewPoints(cPoints, output+"Points");
    gsInfo << "Wrote " << output+"MultiPatch.pvd, " << output+"MultiPatch_*.vts, "<< output+"MultiBasis.pvd, "<< output+"MultiBasis_*.vts"   << "\n\n";

    // Fabricate some data to write to csv file and plot
    gsMatrix<> lineData(1000, 2);
    for (index_t i = 0; i < 1000; ++i)
    {
        lineData(i,0) = i/1000.0;
        lineData(i,1) = std::sin(2*EIGEN_PI*i/1000.0) + 0.25*std::sin(8*EIGEN_PI*i/1000.0);
    }
    // Export to csv with headers
    std::vector<std::string> headers = {"x", "f(x)"};
    gsWriteCsv(output + "LineData.csv" ,lineData, headers);
    //! [gsWriteParaview]


    //! [gsExprEvaluator]
    gsExprEvaluator<> evaluator; 
    gsExprEvaluator<>::geometryMap geoMap = evaluator.getMap(mPatch);
    evaluator.setIntegrationElements(mBasis);

    evaluator.writeParaview( meas(geoMap), geoMap, output+"ExprEval");
    evaluator.writeParaview(usn(geoMap), geoMap, output+"ExprEval");
    gsInfo << "Wrote " << output+"ExprEval.pvd"<<","<<output+"ExprEval.vts\n\n";
    //! [gsExprEvaluator]


    //! [gsParaviewCollection]

    // Initialize the Paraview Collection  with the desired filename of the .pvd file 
    // and the evaluator that will be used to evaluate the expressions (optionally)
    gsParaviewCollection PVCollection(output + "PVCollection.pvd", &evaluator);

    // Number of evaluation points per patch
    PVCollection.options().setInt("numPoints", 1000);
    // Number of decimal points in the output
    PVCollection.options().setInt("precision", 5);
    // Plot the element mesh and set it's resolution
    PVCollection.options().setSwitch("plotElements", false);
    PVCollection.options().setInt("plotElements.resolution", -1);
    // Plot the control net
    PVCollection.options().setSwitch("plotControlNet", false);
    // Export the vtk files to a subfolder
    PVCollection.options().setSwitch("makeSubfolder", true);
    // Name of the subfolder
    PVCollection.options().setString("subfolder", "");
    // Export in base64 binary format
    PVCollection.options().setSwitch("base64", false);

    // Create a new timestep (e.g for the initial state of the problem)
    // This initialises the appropriate files on disk
    PVCollection.newTimeStep(&mPatch);
    // Write the measure of the geometry to the files
    PVCollection.addField( meas(geoMap), "Measure");
    // Write the unit surface normal to the files
    PVCollection.addField(  usn(geoMap), "Surface normal");
    // Save the timestep, so that all files are complete.
    PVCollection.saveTimeStep();

    // Give ellapsed time as seed for random numbers
    std::srand(std::time(0)); 
    // Perturb each patch of the multi patch with random numbers
    // This simulates a deformation of the geometry that might 
    // be the result of a structural simulation
    for (index_t p=0; p < (index_t)mPatch.nPatches(); ++p)
    {
        mPatch.patch(p).coefs() +=  gsEigen::MatrixXd::Random(mPatch.patch(p).coefs().rows(), mPatch.patch(p).coefs().cols());
    }

    // Perform the same steps as before for the deformed geometry
    PVCollection.newTimeStep(&mPatch);
    PVCollection.addField( meas(geoMap), "Measure");
    PVCollection.addField(  usn(geoMap), "Surface normal");
    PVCollection.saveTimeStep();

    // Save the collection, so that the main .pvd file is completed.
    PVCollection.save(); 
    gsInfo << "Wrote " << output + "PVCollection.pvd and all files referenced therein."<< "\n\n";
    //! [gsParaviewCollection]

    //! [Bezier]
    if (bezier) 
    {
        gsWriteParaviewBezier(mPatch, output+"Bezier");
        gsInfo << "Wrote " << output+"Bezier" << ".vtu" << ", " << output+"Bezier" << ".pvd\n";
    }
    //! [Bezier]

    return EXIT_SUCCESS;
}
