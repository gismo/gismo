/** @file biharmonic_multiPatch.cpp

    @brief A Biharmonic example for ONLY TWO-PATCHES

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
# include <gismo.h>
# include <gsG1Basis/gsG1AuxiliaryEdgeMultiplePatches.h>
# include <gsG1Basis/gsG1AuxiliaryVertexMultiplePatches.h>
# include <gsAssembler/gsG1BiharmonicAssembler.h>
# include <gsG1Basis/gsG1System.h>
# include <gsG1Basis/gsG1Mapper_pascal.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Geometry data
    index_t geometry = 0; // Which geometry

    index_t numRefine = 4;
    index_t numDegree = 0;
    index_t regularity = 1;

    // For the spline space of the gluing data
    index_t p_tilde = 1;
    index_t r_tilde = 0;

    index_t threads = 1;

    bool plot = false;
    bool direct = false;
    bool local = false;
    bool loop = false;
    bool local_g1 = false;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("k", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "p_tilde", "Polynomial degree for tilde{p}", p_tilde);
    cmd.addInt("r", "r_tilde", "Regularity for tilde{r}", r_tilde);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    cmd.addSwitch( "direct", "Construction of the G1 basis functions", direct );
    cmd.addSwitch( "loop", "If you want to solve several levels", loop );
    cmd.addSwitch( "local_g1", "If you want to solve several levels", local_g1 );
    cmd.addInt("g", "geometry", "Geometry", geometry);
    cmd.addInt("t", "threads", "Threads", threads);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======= Solution =========
    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> solVal("x",2);
    gsFunctionExpr<>sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                             "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
    gsFunctionExpr<>sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                             "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                             " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)", 2);
    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der, sol2der);

    // ======= Geometry =========
    std::string string_geo;
    switch(geometry)
    {
        case 0:
            string_geo = "planar/multiPatches/4_square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 1:
            string_geo = "planar/multiPatches/6_square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 2:
            string_geo = "planar/multiPatches/square_curved.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 3:
            string_geo = "planar/twoPatches/square_curved.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 4:
            string_geo = "planar/twoPatches/square_curved_deg_5.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 5:
            string_geo = "planar/twoPatches/square_curved_deg_7.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 6:
            string_geo = "planar/twoPatches/square_non_conform.xml";
            numDegree = 2; // 0 == degree 3
            break;
        case 7:
            string_geo = "planar/twoPatches/square_bent.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 8:
            string_geo = "planar/twoPatches/square_complex_bent.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 9:
            string_geo = "planar/twoPatches/square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        default:
            gsInfo << "No geometry is used! \n";
            break;
    }

    gsFileData<> fd(string_geo);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> multiPatch;
    fd.getId(0, multiPatch); // id=0: Multipatch domain
    multiPatch.computeTopology();


    gsWriteParaview(multiPatch,"geometry",5000,true);

    gsOptionList optionList;
    optionList.addInt("p_tilde","Grad",p_tilde);
    optionList.addInt("r_tilde","Reg",r_tilde);
    optionList.addInt("regularity","Regularity of the initial geometry",regularity);
    optionList.addSwitch("local","Local projection for gluing data",local);
    optionList.addSwitch("direct","Local projection for gluing data",direct);
    optionList.addSwitch("plot","Plot in Paraview",plot);
    optionList.addInt("refine","Refinement",numRefine);
    optionList.addInt("degree","Degree",numDegree);

    //multiPatch.patch(1).degreeElevate(1,0);
    multiPatch.degreeElevate(optionList.getInt("degree"));

    multiPatch.uniformRefine_withSameRegularity(optionList.getInt("refine"), optionList.getInt("regularity"));
    gsMultiBasis<> mb(multiPatch);


    gsWriteParaview(mb.basis(0),"basis",5000);

/*
    // Interface loop
    std::vector<gsG1AuxiliaryPatch> g1_interface;
    for (const boundaryInterface &  item : multiPatch.interfaces() )
    {
        gsG1AuxiliaryEdgeMultiplePatches a(multiPatch, item.first().patch, item.second().patch);
        a.computeG1InterfaceBasis(optionList);
        g1_interface.push_back(a.getSinglePatch(0));
        g1_interface.push_back(a.getSinglePatch(1));
    }

    std::vector<gsG1AuxiliaryPatch> g1_boundaries;
    for (gsMultiPatch<>::const_biterator bit = multiPatch.bBegin(); bit != multiPatch.bEnd(); ++bit)
    {
        gsInfo << "Patch: " << bit->patch << "\n";
        gsInfo << "m_index: " << bit->m_index << "\n";
        gsG1AuxiliaryEdgeMultiplePatches a(multiPatch, bit->patch);
        a.computeG1BoundaryBasis(optionList, bit->m_index);
        g1_boundaries.push_back(a.getSinglePatch(0));
    }
*/
    gsG1Mapper_pascal<real_t> g1Mapper(multiPatch);

    std::vector<gsG1AuxiliaryPatch> g1_edges;
    for (size_t np = 0; np < multiPatch.nPatches(); np++)
    {
        for (index_t side_index = 1; side_index < 5; side_index++)
        {
            gsG1AuxiliaryEdgeMultiplePatches a(multiPatch, np);
            a.computeG1EdgeBasis(optionList,side_index,multiPatch.isBoundary(np,side_index));
            g1Mapper.markGlobalIndex_Edge(g1_edges.size(),side_index,np,multiPatch.isBoundary(np,side_index)); // This case is the same
            g1_edges.push_back(a.getSinglePatch(0));
        }
    }

    g1Mapper.markBoundary_Edge(g1_edges, multiPatch);

//     Vertices loop
    gsG1AuxiliaryPatch init;
    std::vector<gsG1AuxiliaryPatch> g1_vertices(4*multiPatch.nPatches(),init);
    std::vector<std::vector<patchCorner>> allcornerLists = multiPatch.vertices();
    for(size_t i=0; i < allcornerLists.size(); i++)
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

        gsG1AuxiliaryVertexMultiplePatches a(multiPatch, patchIndex, vertIndex);
        a.computeG1InternalVertexBasis(optionList);
        index_t kindBdr = a.kindOfVertex();
        for (size_t j = 0; j < vertIndex.size(); j++)
        {
            g1Mapper.markGlobalIndex_Vertex(patchIndex.at(j)*4 + vertIndex.at(j) -1,vertIndex.at(j),patchIndex.at(j),kindBdr); // TODO
            g1_vertices.at(patchIndex.at(j)*4 + vertIndex.at(j) -1) = a.getSinglePatch(j);
        }
        g1Mapper.coupleVertex(patchIndex, vertIndex);

    }
    g1Mapper.markBoundary_Vertex(multiPatch,g1_vertices);

    gsG1System<real_t> g1System(mb, g1_edges, g1_vertices, g1Mapper);
    //if (plot)
    //g1System.plotParaview(multiPatch,g1_interface,g1_boundaries,g1_vertices,"BasisFunctions");
    g1System.plotParaview(multiPatch,g1_edges,g1_vertices,"BasisFunctions");


// NEW NEW NEW NEW NEW NEW NEW NEW NEW

    gsBoundaryConditions<> bcInfo, bcInfo2;
    for (gsMultiPatch<>::const_biterator bit = multiPatch.bBegin(); bit != multiPatch.bEnd(); ++bit)
    {
        bcInfo.addCondition( *bit, condition_type::dirichlet, &solVal ); // = 0
        bcInfo2.addCondition(*bit, condition_type::neumann, &laplace ); // = 0
    }

    // BiharmonicAssembler
    gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(multiPatch, mb, bcInfo, bcInfo2, source);
    g1BiharmonicAssembler.assemble();

    g1BiharmonicAssembler.computeDirichletDofsL2Proj(g1_edges, g1_vertices, g1Mapper );
    g1BiharmonicAssembler.constructDirichletSolution(g1_edges, g1_vertices, g1Mapper );

    g1System.setGlobalMapper(g1BiharmonicAssembler.get_mapper());
    g1System.initialize();
    g1System.assemble();














} // main