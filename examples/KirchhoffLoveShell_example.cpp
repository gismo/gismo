/** @file KirchhoffLoveShell_example.cpp

    @brief A Kirchhoff-Love example.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
*/

# include <gismo.h>
# include <gsAssembler/gsKirchhoffLoveShellAssembler.h>
# include <gsG1Basis/gsG1AuxiliaryEdgeMultiplePatches.h>
# include <gsG1Basis/gsG1AuxiliaryVertexMultiplePatches.h>


using namespace gismo;

int main(int argc, char *argv[])
{
    index_t numRefine = 5;
    index_t numDegree = 1;
    bool plot = true;

    gsCmdLine cmd("Example for solving the Kirchhoff-Love problem.");
    cmd.addInt("r", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "degree", "Polynomial degree", numDegree);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    dirichlet::strategy dirStrategy = dirichlet::elimination;
    iFace::strategy intStrategy = iFace::glue;

    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    gsFunctionExpr<>sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                             "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
    gsFunctionExpr<>sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                             "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                             " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)", 2);
    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der, sol2der);

    //gsFileData<> fileSrc("KirchhoffLoveGeo/geo_fivePatchDiffParam.xml");
    gsFileData<> fileSrc("KirchhoffLoveGeo/square_diffParam1.xml");

    gsInfo << "Loaded file " << fileSrc.lastPath() << "\n";

    gsMultiPatch<> geo;
    gsInfo << "Geometry taken correctly \n";
    fileSrc.getId(2, geo);
    geo.computeTopology();
    gsInfo << "Geometry computed correctly\n";


    gsOptionList optionList;
    optionList.addInt("p_tilde","Grad",1);
    optionList.addInt("r_tilde","Reg",0);
    optionList.addInt("regularity","Regularity of the initial geometry",1);
    optionList.addInt("refine","Plot in Paraview",numRefine);
    optionList.addInt("degree","Degree",numDegree);
    optionList.addSwitch("local","Local projection for gluing data",false);
    optionList.addSwitch("direct","Local projection for gluing data",false);
    optionList.addSwitch("plot","Plot in Paraview",false);

    geo.degreeElevate(optionList.getInt("degree"));
    geo.uniformRefine_withSameRegularity(optionList.getInt("refine"),optionList.getInt("regularity"));
    gsMultiBasis<> basis(geo);
    gsInfo << "Old: " << basis << "\n";




//     Interface loop
//    for (const boundaryInterface &  item : geo.interfaces() )
//    {
//
//
//        gsG1AuxiliaryEdgeMultiplePatches a(geo, item.first().patch, item.second().patch);
//
//        a.computeG1InterfaceBasis(optionList);
//
//    }



    //     Loop over the boundary edges
   // for ( auto & it : geo.boundaries()){
//        gsInfo << "Patch: " << geo.boundaries()[0].patch << "\n";
//        gsInfo << "m_index: " << geo.boundaries()[0].m_index << "\n";
//        gsG1AuxiliaryEdgeMultiplePatches a(geo, geo.boundaries()[0].patch);
//        a.computeG1BoundaryBasis(optionList, geo.boundaries()[0].m_index);
    //}



//     Vertices loop
    std::vector<std::vector<patchCorner>> allcornerLists = geo.vertices();
    size_t i = 1;
//    for(size_t i=0; i < allcornerLists.size(); i++)
    {
        std::vector<size_t> patchIndex;
        std::vector<size_t> vertIndex;
        for(size_t j = 0; j < allcornerLists[i].size(); j++)
        {
            patchIndex.push_back(allcornerLists[i][j].patch);
            vertIndex.push_back(allcornerLists[i][j].m_index);
            gsInfo << "Patch: " << allcornerLists[i][j].patch << "\t Index: " << allcornerLists[i][j].m_index << "\n";

        }
        gsInfo << "\n";

        gsG1AuxiliaryVertexMultiplePatches a(geo, patchIndex, vertIndex);
        a.computeG1InternalVertexBasis(optionList);
    }





//    for (const std::vector<patchCorner> & it : allcornerLists)
//    {
//        gsInfo << "Dimension of the vector: " << it.size() << "\n";
//        gsInfo << "Corner " << it.at(0).m_index << " in Patch " << it.at(0).patch << "\n";
//        for (const patchCorner & it_corner : it)
//        {
//            gsInfo << "Patch : " << it_corner.patch << "\t Corner: " << it_corner.m_index << "\n";
//        }
//    }



//    gsWriteParaview(newgeom1, "Geometry", 1000);

//    // Write file .xml of the new geometry
//    gsFileData<> fd;
//    fd << test;
//    // output is a string. The extention .xml is added automatically
//    fd.save("newGeo");


    //Setting up oundary conditions
    gsBoundaryConditions<> bcInfo;
    gsBoundaryConditions<> bcInfo2;
    for (gsMultiPatch<>::const_biterator
             bit = geo.bBegin(); bit != geo.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &solution );
        bcInfo2.addCondition( *bit,  condition_type::neumann, &laplace);
    }




    //Initilize solver
    gsKirchhoffLoveShellAssembler<real_t> KirchhoffLoveShellAssembler( geo,basis,bcInfo,bcInfo2,source,
                                                       dirStrategy, intStrategy);

    gsInfo<<"Assembling..." << "\n";
    KirchhoffLoveShellAssembler.assemble();

    gsInfo<<"Solving with direct solver, "<< KirchhoffLoveShellAssembler.numDofs()<< " DoFs..."<< "\n";
    gsSparseSolver<real_t>::LU solver;
    solver.analyzePattern(KirchhoffLoveShellAssembler.matrix() );
    solver.factorize(KirchhoffLoveShellAssembler.matrix());
    gsMatrix<> solVector= solver.solve(KirchhoffLoveShellAssembler.rhs());

    //Reconstruct solution
    gsMultiPatch<> mpsol;
    KirchhoffLoveShellAssembler.constructSolution(solVector, mpsol);
    gsField<> solField(KirchhoffLoveShellAssembler.patches(), mpsol);

    //Contruct the H2 norm, part by part.
    real_t errorH2Semi = solField.distanceH2(solution, false);
    real_t errorH1Semi = solField.distanceH1(solution, false);
    real_t errorL2 = solField.distanceL2(solution, false);
    real_t errorH1 = math::sqrt(errorH1Semi*errorH1Semi + errorL2*errorL2);
    real_t errorH2 = math::sqrt(errorH2Semi*errorH2Semi + errorH1Semi*errorH1Semi + errorL2*errorL2);

    gsInfo << "The L2 error of the solution is : " << errorL2 << "\n";
    gsInfo << "The H1 error of the solution is : " << errorH1 << "\n";
    gsInfo << "The H2 error of the solution is : " << errorH2 << "\n";

    // Plot solution in paraview
    if (plot)
    {
        // Write approximate and exact solution to paraview files
        gsInfo<<"Plotting in ParaView...\n";
        gsWriteParaview<>(solField, "KirchhoffLoveShell2d", 5000);
        const gsField<> exact( geo, solution, false );
        gsWriteParaview<>( exact, "KirchhoffLoveShell2d_exact", 5000);
    }
    else
        gsInfo << "Done. No output created, re-run with --plot to get a ParaView "
                  "file containing the solution.\n";

    return  0;
}
