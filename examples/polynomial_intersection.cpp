/** @file polynomial_intersection.cpp

    @brief Demonstrates adaptive fitting of data samples

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore , Y. Ji
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Options with default values

    std::string fn = "fitting/deepdrawingC.xml";

    // Reading options from the command line
    gsCmdLine cmd("Intersection of two polynomial curves.");
    cmd.addString("d", "data", "Input sample data", fn);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    index_t deg = 3;
    index_t order = deg+1;

    gsKnotVector<> ku(0, 1, 5, deg+1);//start,end,interior knots, start/end multiplicites of knots
    gsMatrix<> coefs(9, 2);

    coefs << 0, 0,
             0., 0.25,
             0.125, 0.5,
             0.25, 0.55,
             0.375, 0.5,
             0.5, 0.25,
             0.625, 0.25,
             0.8, 0.25,
             1, 0.5;

    gsMatrix<> coefs_plot;
    coefs_plot = coefs.transpose();
    gsWriteParaviewPoints(coefs_plot, "coefs_original");
    gsBSpline<> curve1( ku, give(coefs));








    gsInfo << "curve:\n" << curve1 << "\n";

    gsWriteParaview( curve1, "c1", 1000);
    gsWriteParaview( curve1, "c1", 1000, false, true);


    // TODO: add bezier extraction, by proper knot insertion algorithm.
    auto multContainer = ku.multiplicities();
    gsInfo << multContainer[0] << "\n";
    //gsDebugVar(multContainer);
    // for(index_t i = 1; i < multContainer.size()-1; i++)
    // {
    //   c1.insertKnot(ku.knot(i), degree - multContainer[i]);
    // }
    //
    ku.increaseMultiplicity(deg-1, false);

    gsMatrix<> X, uv;
    ku.greville_into(uv);
    curve1.eval_into(uv, X);

    gsInfo << "Points:\n" << X << "\n";
    gsInfo << "Parameters:\n" << uv << "\n";
    gsWriteParaviewPoints(X, "points");


    gsBSplineBasis<> basis(ku);
    gsFitting<> bezierExtraction(uv, X, basis);
    bezierExtraction.compute();

    gsMatrix<> newCoefs = bezierExtraction.result()->coefs();
    gsDebugVar(basis);
    gsDebugVar(newCoefs);

    gsBSpline<> bezierExtractionCurve( basis, newCoefs);
    gsWriteParaview(bezierExtractionCurve, "bezierC1", 1000);


    gsKnotVector<> kv(0, 1, 0, 2);//start,end,interior knots, start/end multiplicites of knots
    gsMatrix<> coefs2(2, 2);

    coefs2 << 0.2, 0,
             0.2, 2;

   gsBSpline<> curve2( kv, give(coefs2));

   gsInfo << "curve:\n" << curve2 << "\n";

   gsWriteParaview( curve2, "line", 1000);
   gsWriteParaview( curve2, "line", 1000, false, true);


   //find bounding box for bezier curves segments
   gsInfo << "----------------------------------------------------\n";
   gsInfo << deg <<"\n";
   gsKnotVector<> auxKnt(0, 1, 0, deg+1);//start,end,interior knots, start/end multiplicites of knots
   gsBSplineBasis<> basisAux(auxKnt);

   gsInfo << basisAux << "\n";
   gsMatrix<> coefsPatch(deg+1,2);
   gsInfo << "# Bezier patches =  " << basis.size() / deg << "\n";
   gsInfo << newCoefs.rows() << " x " << newCoefs.cols() << "\n";
   for(index_t i = 0; i < basis.size()/deg-1; i++)
   {
     // gsInfo << newCoefs.middleRows(i*(deg), deg+1) << "\n\n";
     coefsPatch = newCoefs.middleRows(i*(deg), deg+1);
     gsBSpline<> bezierPatch(basisAux, coefsPatch);
   }


   gsMultiPatch<> mp(bezierExtractionCurve);
   gsMatrix<> bb;
   mp.boundingBox(bb);

   gsInfo << "bounding box:\n" << bb << "\n";

    return 0;
}
