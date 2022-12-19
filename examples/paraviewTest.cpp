#include <iostream>
#include <fstream>
#include <gismo.h>
#include <gsIO/gsFileManager.h>


using namespace gismo;

int main(int argc, char const *argv[])
{
    gsTensorBSpline<2>  bottom = *gsNurbsCreator<>::BSplineRectangle(-10,-10,10,0);
    gsTensorBSpline<2>  top    = *gsNurbsCreator<>::BSplineRectangle(-5,0,5,10);

    gsMultiPatch<> mPatch;
    mPatch.addPatch(top);
    mPatch.addPatch(bottom);


    gsExprAssembler<> assembler(1,1);
    gsMultiBasis<> mBasis(mPatch);
    assembler.setIntegrationElements(mBasis);

    // Get handles to the geometry map to be used in expressions.
    gsExprAssembler<>::geometryMap initGeo = assembler.getMap(mPatch); // Geometry map in the initial configuration


    gsExprEvaluator<> evaluator(assembler);

    // Create a Collection object [i.e. .pvd file]
    gsParaviewCollection collection("outputFiles/ParaviewExample", &evaluator);

    // Specify plotting options
    collection.options().setInt("numPoints", 100);
    collection.options().setInt("precision",5);
    collection.options().setSwitch("plotElements",true);
    collection.options().setSwitch("plotControlNet",true);

    // Create a DataSet object [i.e. a group of .vts files]
    collection.newTimeStep(&mPatch);

    // Solution field names
    std::vector<std::string> labels;
    labels.push_back("measure");
    labels.push_back("norm");

    // Evaluate and write expressions to the .vts files
    collection.addFields( labels, meas(initGeo), initGeo.norm() );
    // Reference the .vts files in the .pvd file
    collection.saveTimeStep();

    // Dummy Solution Loop
    for (int i=0;i!=3;i++)
    {
        // Create a clean DataSet [ new .vts files]
        collection.newTimeStep(&mPatch);

        // Dummy solution, just displaces one of the patches
        mPatch.patch(0).coefs().array() += 1;

        // Evaluate and write expressions to the .vts files
        collection.addFields( labels, meas(initGeo), initGeo.norm() );

        // Reference the .vts files in the .pvd file
        collection.saveTimeStep();

    }
    // Solution has finished

    // Finalize and save the main [.pvd] file
    collection.save(); 

    // Just to check I have not messed up the current implementation
    // evaluator.writeParaview(initGeo.norm(),initGeo, "evOutput");
    // gsWriteParaview(mBasis, mPatch, "test", 1000);

    return 0; 
}