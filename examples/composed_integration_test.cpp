/** @file composed_integration_test.cpp

    @brief Tutorial on how to use expression assembler to solve the Poisson equation

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

//! [Include namespace]
#include <gismo.h>
#include <gsNurbs/gsSquareDomain.h>
#include <iostream>
#include <numeric>

using namespace gismo;
//! [Include namespace]

// Function to return LCM of two numbers
index_t LCM(index_t a, index_t b)
{
    index_t greater  = std::max(a, b);
    index_t smallest = std::min(a, b);
    for (index_t i = greater; ; i += greater) {
        if (i % smallest  == 0)
            return i;
    }
}




template <short_t _DIM, class T>
gsTensorBSplineBasis<_DIM,T> integrationBasis(const gsTensorBSplineBasis<_DIM,T> & basis1,
                                           const gsTensorBSplineBasis<_DIM,T> & basis2)
{
    gsTensorBSplineBasis<_DIM,T> ibasis(basis1);
    // Integration basis: parent basis with knots of composition basis inserted, and the degree is the sum of the two degrees (?)
    index_t targetDegree;
    for (size_t d = 0; d!=_DIM; d++)
    {
        // 1. Insert interior knots of composition basis
        for (typename gsKnotVector<T>::uiterator it = std::next(basis2.knots(d).ubegin());
                                                    it!= std::prev(basis2.knots(d).uend());
                                                    ++it)
            {
                if (ibasis.knots(d).has(*it))
                    continue;
                ibasis.insertKnot(*it,d);
            }
        // 2. Increase the degree
        index_t lcmdeg = LCM(ibasis.degree(d), basis2.degree(d));
        gsInfo << "LCM of degrees: " << lcmdeg << "\n";
        targetDegree = ibasis.degree(d) * basis2.degree(d);
        ibasis.degreeIncrease(targetDegree-ibasis.degree(d),d);

    }
    return ibasis;
}

int main(int argc, char *argv[])
{
    // Input options
    index_t numKnots    = 3;
    index_t degree      = 3;
    index_t numKnotsMap = 2;
    index_t degreeMap   = 2;
    bool plot           = false;

    std::string options = "options/assembler_options.xml";
    std::string fun     = "x^2*y^2";
    gsCmdLine cmd("Shell modal solver.");

    cmd.addInt("n", "numKnots", "Number of interior knots", numKnots);
    cmd.addInt("p", "degree", "Degree of the basis", degree);
    cmd.addInt("m", "numKnotsMap", "Number of interior knots for the mapping", numKnotsMap);
    cmd.addInt("q", "degreeMap", "Degree of the basis for the mapping", degreeMap);
    cmd.addString("o", "options", "Path to the options file", options);
    cmd.addSwitch("plot", "Plot result in ParaView format", plot);
    cmd.addString("f","fun","Function to integrate",fun);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsOptionList evaluatorOptions;
    gsFileData<> fd(options);
    fd.getFirst<gsOptionList>(evaluatorOptions);

    typedef real_t T;

    gsKnotVector<T> kv1(0,1,numKnotsMap,degreeMap+1);
    gsKnotVector<T> kv2(0,1,numKnots   ,degree+1   );

    gsTensorBSplineBasis<2,T> basisMap(kv1,kv1);
    gsTensorBSplineBasis<2,T> basisAnalysis(kv2,kv2);
    gsTensorBSplineBasis<2,T> basisI = integrationBasis<2,T>(basisMap,basisAnalysis);

    gsInfo<<"Knots of basis 1: "<<kv1<<"\n";
    gsInfo<<"Knots of basis 2: "<<kv2<<"\n";
    gsInfo<<"Knots of integration basis: "<<basisI.knots(0)<<"\n";
    // Construct a square domain
    gsSquareDomain<T> squareDomain(basisMap);
    gsVector<T> controls = squareDomain.getControls();
    controls += basisAnalysis.knots(0).maxIntervalLength()*0.25*gsVector<T>::Random(controls.size());
    squareDomain.setControls(controls);

    // Construct a composed basis
    gsComposedBasis<T> composedBasis(&squareDomain, &basisI);

    // Construct a square (for the geometry map)
    gsGeometry<T>::uPtr domain = gsNurbsCreator<T>::BSplineSquare();
    gsComposedGeometry<T> composedGeometry(squareDomain, *domain);

    // Construct a function
    gsFunctionExpr<T> f(fun,2);
    gsComposedFunction<T> cf(&composedGeometry, &f);

    // Integrate
    gsExprEvaluator<T> ev;
    ev.options().update(evaluatorOptions,gsOptionList::addIfUnknown);
    auto G  = ev.getMap(*domain);
    auto F  = ev.getVariable(f);
    auto cG = ev.getMap(composedGeometry);
    auto cF = ev.getVariable(cf);

    gsInfo<<"Original basis:\n";
    gsMultiBasis<T> ob(basisAnalysis);
    ev.setIntegrationElements(ob);
    gsInfo<<"* int( F*| G|) = "<<ev.integral( F * meas( G))<<"\n";
    gsInfo<<"* int(cF*|cG|) = "<<ev.integral(cF * meas(cG))<<"\n";
    gsInfo<<"Integration basis:\n";
    gsMultiBasis<T> ib(basisI);
    ev.setIntegrationElements(ib);
    gsInfo<<"* int( F*| G|) = "<<ev.integral( F * meas( G))<<"\n";
    gsInfo<<"* int(cF*|cG|) = "<<ev.integral(cF * meas(cG))<<"\n";

    // Export quadrature points
    gsMatrix<T> nodesMap = gsQuadrature::getAllNodes(basisMap,evaluatorOptions);
    gsMatrix<T> nodesAnalysis = gsQuadrature::getAllNodes(basisAnalysis,evaluatorOptions);
    gsMatrix<T> nodesI = gsQuadrature::getAllNodes(basisI,evaluatorOptions);
    gsMatrix<T> cnodesMap, cnodesAnalysis, cnodesI;

    squareDomain.eval_into(nodesMap,cnodesMap);
    squareDomain.eval_into(nodesAnalysis,cnodesAnalysis);
    squareDomain.eval_into(nodesI,cnodesI);

    gsWriteParaviewPoints(nodesMap,"nodesMap");
    gsWriteParaviewPoints(nodesAnalysis,"nodesAnalysis");
    gsWriteParaviewPoints(nodesI,"nodesI");
    gsWriteParaviewPoints(cnodesMap,"cnodesMap");
    gsWriteParaviewPoints(cnodesAnalysis,"cnodesAnalysis");
    gsWriteParaviewPoints(cnodesI,"cnodesI");

    gsMesh<T> meshMap(basisMap,2);
    gsMesh<T> meshAnalysis(basisAnalysis,2);
    gsMesh<T> meshI(basisI,2);

    gsWriteParaview(meshMap,"meshMap",false);
    gsWriteParaview(meshAnalysis,"meshAnalysis",false);
    gsWriteParaview(meshI,"meshI",false);

    gsWriteParaview(composedGeometry,"composedGeometry",1000, true, true);
    gsWriteParaview(squareDomain.domain(),"squareDomain",1000,true,true);
    gsWriteParaview(composedBasis,"composedBasis",1000, true);


    std::vector<std::string> headers = {"X","Y"};
    nodesMap.transposeInPlace();
    nodesAnalysis.transposeInPlace();
    nodesI.transposeInPlace();
    cnodesMap.transposeInPlace();
    cnodesAnalysis.transposeInPlace();
    cnodesI.transposeInPlace();
    gsWriteCsv("nodesMap.csv",nodesMap,headers);
    gsWriteCsv("nodesAnalysis.csv",nodesAnalysis,headers);
    gsWriteCsv("nodesI.csv",nodesI,headers);
    gsWriteCsv("cnodesMap.csv",cnodesMap,headers);
    gsWriteCsv("cnodesAnalysis.csv",cnodesAnalysis,headers);
    gsWriteCsv("cnodesI.csv",cnodesI,headers);


    return EXIT_SUCCESS;

}// end main
