/** @file adaptRefinementThb_example.cpp

    @brief Tutorial on how to use G+Smo to solve the Poisson equation,
    see the \ref PoissonTutorial

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Kleiss
*/

//! [Include namespace]
# include <gismo.h>
# include <gsAssembler/gsAdaptiveRefUtils.h>

using namespace gismo;
//! [Include namespace]

int main(int argc, char *argv[])
{
   //! [Parse command line]
    bool plot = false;
    index_t numRefine  = 0;
    index_t numElevate = 0;
    bool last = false;
    int function = 0;

    gsCmdLine cmd("Tutorial on solving a Poisson problem.");
    cmd.addInt( "e", "degreeElevation",
                "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate );
    cmd.addInt( "r", "uniformRefine", "Number of Uniform h-refinement steps to perform before solving",  numRefine );
    cmd.addInt( "i", "function", "Function to refine",  function );
    cmd.addSwitch("last", "Solve solely for the last level of h-refinement", last);
    cmd.addSwitch("plot", "Create a ParaView visualization file with the solution", plot);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    //! [Parse command line]

    gsMultiPatch<> mp0;
    mp0 = *gsNurbsCreator<>::BSplineRectangle(0,0,1,1);

    //! [Read input file]

    //! [Refinement]

    // Elevate and p-refine the basis to order p + numElevate
    // where p is the highest degree in the bases
    mp0.degreeElevate(numElevate);
    // h-refine each basis
    for (int r =0; r < numRefine; ++r)
        mp0.uniformRefine();

    gsInfo  <<"Original : "<< mp0.basis(0)<<"\t"
            <<"size = "<< mp0.basis(0).size()<<"\n";

    // gsWriteParaview<>(mp0, "mp", 1000, true);

    // Cast all patches of the mp object to THB splines
    gsMultiPatch<> mp;
    gsTHBSpline<2,real_t> thb;
    for (index_t k=0; k!=mp0.nPatches(); ++k)
    {
        gsTensorBSpline<2,real_t> *geo = dynamic_cast< gsTensorBSpline<2,real_t> * > (&mp0.patch(k));

        thb = gsTHBSpline<2,real_t>(*geo);
    // gsDebugVar(*geo);
    // gsDebugVar(*basis0);
        mp.addPatch(thb);
    }

    //gsWriteParaview<>(mp, "mp", 1000, true);

    gsMultiBasis<> basis(mp);
    gsWriteParaview<>(basis.basis(0), "oldBasis", 1000, true);
    gsWriteParaview<>(mp, "mp", 1000, true);

    gsInfo<<"Basis Primal: "<<basis.basis(0)<<"\n";
    gsInfo << "Patches: "<< mp.nPatches() <<", degree: "<< basis.minCwiseDegree() <<"\n";


    gsTensorBSplineBasis<2,real_t> *basis0 = dynamic_cast< gsTensorBSplineBasis<2,real_t> * > (&mp0.basis(0));


    gsHBSplineBasis<2,real_t> THBbasis(*basis0);
    THBbasis.refineBasisFunction(function);

    std::vector<bool> funMarked( mp.basis(0).size() );
    for(std::vector<bool>::iterator i = funMarked.begin(); i!=  funMarked.end(); ++i)
        *i = false;
    funMarked[function] = true;

    // PRINT
    for (std::vector<bool>::const_iterator i = funMarked.begin(); i != funMarked.end(); ++i)
        gsInfo << *i << ' ';
    gsInfo<<"\n";
    // ! PRINT

    gsRefineMarkedFunctions(mp,funMarked,0);


    gsInfo<<"Basis After: "<<mp.basis(0)<<"\n";
    gsInfo  <<"Original : "<< mp.basis(0)<<"\t"
            <<"size = "<< mp.basis(0).size()<<"\n";

    gsWriteParaview<>(mp.basis(0), "newBasisL", 1000, true);
    gsWriteParaview<>(mp, "mpL", 1000, true);

    mp.degreeElevate(1);
    gsWriteParaview<>(mp.basis(0), "newBasisH", 1000, true);
    gsWriteParaview<>(mp, "mpH", 1000, true);


    mp0.basis(0).uniformRefine();
    gsWriteParaview<>(mp0.basis(0), "uniformRefine", 1000, true);

   return EXIT_SUCCESS;

}// end main
