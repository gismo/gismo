#include <iostream>
#include <fstream>
#include <gismo.h>


using namespace gismo;

int main(int argc, char const *argv[])
{
    gsTensorBSpline<2>  bottom = *gsNurbsCreator<>::BSplineRectangle(-10,-10,10,0);
    gsTensorBSpline<2>  top    = *gsNurbsCreator<>::BSplineRectangle(-5,0,5,10);

    gsMultiPatch<> mPatch;
    mPatch.addPatch(top); // top - 0 - slave
    mPatch.addPatch(bottom); // bottom - 1 - master


    gsExprAssembler<> assembler(1,1);
    // assembler.options() = assembler.defaultOptions();
    gsMultiBasis<> mBasis(mPatch);
    assembler.setIntegrationElements(mBasis);

    gsMultiPatch<> currentMP = mPatch;  // Make a copy of mPatch for the current configuration
    // Get handles to the geometry map and test/trial functions to be used in expressions.
    gsExprAssembler<>::geometryMap initGeo = assembler.getMap(mPatch); // Geometry map in the initial configuration


    gsBoundaryConditions<> bc;
    bc.setGeoMap(mPatch);

    gsConstantFunction<> zero(0,2);
    gsConstantFunction<> ten(10,2);
    gsFunctionExpr<> ten_conditional("10*( (x<-5) or (x>5) )",2);

    gsExprEvaluator<> evaluator(assembler);

    //std::vector<std::string> out( evaluator.expr2vtk(initGeo, "Geometry") );
    //gsDebugVar( out[0] );

    //gsParaviewCollection collection("collection", &evaluator);

    gsParaviewCollection collection("collect", &evaluator);

    std::vector<std::string> labels;
    labels.push_back("measure");
    labels.push_back("norm");


    gsParaviewDataSet dSet = collection.newTimeStep(&initGeo);
    dSet.addFields( labels, meas(initGeo), initGeo.norm() );
    // for ( index_t i=0; i!= dSet.m_numPatches; i++) gsInfo << dSet.filenames()[i] << "\n";
    // dSet.save();

    collection.addDataSet(dSet);

    mPatch.patch(0).coefs().array() += 1;

    dSet = collection.newTimeStep(&initGeo);
    dSet.addFields( labels, meas(initGeo), initGeo.norm() );

    collection.addDataSet(dSet);

    collection.save(); 



    return 0; 
}