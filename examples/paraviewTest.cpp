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

    gsParaviewDataSet dSet("test_file", &initGeo, &evaluator);
    dSet.addField( meas(initGeo) ,"measure");
    for ( index_t i=0; i!= dSet.m_numPatches; i++) gsInfo << dSet.filenames()[i] << "\n";
    dSet.save();



    return 0; 
}