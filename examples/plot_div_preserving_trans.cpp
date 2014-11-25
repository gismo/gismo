/** @file plot_div_preserving_trans.cpp

    @brief Example programs that reconstructs a 2D vector field using differnt tarsformations

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#include <iostream>
#include <gismo.h>
#include <gsCore/gsDivConSolution.h>


using std::cout;
using std::endl;
using namespace gismo;

int main(int argc, char *argv[])
{
    // Define Geometry
    gsMultiPatch<> * patches;
    patches = new gsMultiPatch<>( *safe(gsNurbsCreator<>::BSplineFatQuarterAnnulus()) );
    gsMultiPatch<> * patchesSquare;
    patchesSquare = new gsMultiPatch<>(*safe(gsNurbsCreator<>::BSplineSquare(1.0,0,0)));


    // Define and refine discretization space by refining the basis of the geometry
    std::vector< gsMultiBasis<> >  refine_bases;
    refine_bases.push_back(gsMultiBasis<>( *patches ));//Basis for x direction of the vector field
    refine_bases.push_back(gsMultiBasis<>( *patches ));//Basis for y direction of the vector field

    int numRefine = 3;
    for (int i = 0; i < numRefine; ++i)
    {
        refine_bases[0].uniformRefine();
        refine_bases[1].uniformRefine();
    }

    refine_bases[0].degreeElevateComponent(0);
    refine_bases[1].degreeElevateComponent(0);



    cout << "Degrees of basis for x direction: "<< refine_bases[0].degree(0,0)<< " and "<< refine_bases[0].degree(0,1)<<endl;
    cout << "Degrees of basis for y direction: "<< refine_bases[1].degree(0,0)<< " and "<< refine_bases[1].degree(0,1)<<endl;


    int szMax  = refine_bases[0][0].size();
    if (szMax < refine_bases[1][0].size())
        szMax = refine_bases[1][0].size();

    //Coefficients for vector field in the parametric domain
    //In the parametric domain the field is (0,1) everywhere
    gsMatrix<> coeffs(szMax,2);
    for (index_t k=0; k< szMax; ++k)
    {
        coeffs(k,0) = 0;
        coeffs(k,1) = 1;
    }
    gsMatrix<> coeffsCP = coeffs;
    std::vector<gsBasis<> *> bases_vec;
    bases_vec.push_back(&refine_bases[0][0]);
    bases_vec.push_back(&refine_bases[1][0]);


    //Creating a vector field in the physical domain with
    //divergence preserving transformation (Piola).
    gsMultiPatch<> patchesTest = *patches;
    gsFunction<> * funcPiola = new gsDivConSolution<real_t, gsBSplineBasis<real_t,gsKnotVector<> > > (coeffs,patchesTest[0],bases_vec);
    gsField<> fieldPiola( *patches , funcPiola) ;

    //Creating a vector field in the physical domain with
    //inverse composition (standard IGA)
    gsFunction<> * funcIGA= refine_bases[0][0].makeGeometry( give(coeffs) );
    gsField<> fieldIGA( *patches , funcIGA) ;

    //Creating a vector field in the physical domain with
    //no transform (parametric space)
    gsFunction<> * funcNo= refine_bases[1][0].makeGeometry( give(coeffsCP) );
    gsField<> fieldNoTars( *patchesSquare , funcNo) ;


    // Write solution to paraview files
    gsWriteParaview<>( fieldPiola, "transPiola", 500);
    gsWriteParaview<>( fieldIGA, "transIGA", 500);
    gsWriteParaview<>( fieldNoTars, "transNo", 500);

    delete patches;
    delete patchesSquare;


    cout << "Test is done: Exiting" << endl;
    return 0;
}
