/** @file g1biharmonic_example.cpp

    @brief A Biharmonic example only for AS geometries

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/
# include <omp.h>

# include <gismo.h>
# include <gsG1Basis/gsG1AuxiliaryEdgeMultiplePatches.h>
# include <gsG1Basis/gsG1AuxiliaryVertexMultiplePatches.h>
# include <gsAssembler/gsG1BiharmonicAssembler.h>
# include <gsG1Basis/gsG1System.h>

# include <gsG1Basis/gsG1OptionList.h>

# include <gsG1Basis/gsNormL2.h>
# include <gsG1Basis/gsSeminormH1.h>
# include <gsG1Basis/gsSeminormH2.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    // Geometry data
    index_t geometry = 0; // Which geometry

    index_t numRefine = 4; // Initial refinement
    index_t numDegree = 0; // Degree
    index_t regularity = 1; // Regularity

    // For the spline space of the gluing data
    index_t p_tilde = 1;
    index_t r_tilde = 0;

    index_t threads = 1; // For parallel computing

    index_t loop = 1; // Number of refinement steps

    real_t threshold = 1e-5;

    bool plot = false;
    bool latex = false;
    bool localGd = false;
    bool localEdge = false;
    bool localVertex = false;

    gluingData::strategy gluingData_strategy = gluingData::global;
    g1BasisEdge::strategy g1BasisEdge_strategy = g1BasisEdge::l2projection;
    g1BasisVertex::strategy g1BasisVertex_strategy = g1BasisVertex::global;

    gsCmdLine cmd("Example for solving the biharmonic problem.");
    cmd.addInt("k", "refine", "Number of refinement steps", numRefine);
    cmd.addInt("p", "p_tilde", "Polynomial degree for tilde{p}", p_tilde);
    cmd.addInt("r", "r_tilde", "Regularity for tilde{r}", r_tilde);
    cmd.addInt("g", "geometry", "Geometry", geometry);
    cmd.addInt("t", "threads", "Threads", threads);
    cmd.addInt( "l", "loop", "The number of refinement steps", loop);
    cmd.addSwitch( "plot", "Plot result in ParaView format", plot );
    cmd.addSwitch( "localGd", "To compute the gluing data with local support", localGd );
    cmd.addSwitch( "localEdge", "To compute the G1 edge basis functions with local support", localEdge );
    cmd.addSwitch( "localVertex", "To compute the G1 vertex basis functions with the average dd_ik", localVertex );
    cmd.addSwitch("latex","Print the rate and error latex-ready",latex);
    cmd.addReal("e","threshold", "The threshold for computing the kernel", threshold);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    // ======= Solution =========
/*    gsFunctionExpr<> source  ("4096*pi*pi*pi*pi*(4*cos(8*pi*x)*cos(8*pi*y) - cos(8*pi*x) - cos(8*pi*y))",2);
    gsFunctionExpr<> laplace ("-64*pi*pi*(2*cos(8*pi*x)*cos(8*pi*y) - cos(8*pi*x) - cos(8*pi*y))",2);
    gsFunctionExpr<> solVal("(cos(8*pi*x) - 1) * (cos(8*pi*y) - 1)",2);
    gsFunctionExpr<>sol1der ("-8*pi*(cos(8*pi*y) - 1)*sin(8*pi*x)",
                             "-8*pi*(cos(8*pi*x) - 1)*sin(8*pi*y)",2);
    gsFunctionExpr<>sol2der ("-64*pi^2*(cos(8*pi*y) - 1)*cos(8*pi*x)",
                             "-64*pi^2*(cos(8*pi*x) - 1)*cos(8*pi*y)",
                             " 64*pi^2*sin(8*pi*x)*sin(8*pi*y)", 2);
*/
    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    gsFunctionExpr<>sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                             "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
    gsFunctionExpr<>sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                             "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                             " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)", 2);
/*
    // Solution of Marios paper
    gsFunctionExpr<> source  ("-1*cos(x/2)*sin(y/2)",2);
    gsFunctionExpr<> laplace ("2*cos(x/2)*sin(y/2)",2);
    gsFunctionExpr<> solVal("-4*cos(x/2)*sin(y/2)",2);
    gsFunctionExpr<>sol1der ("2*sin(x/2)*sin(y/2)",
                             "-2*cos(x/2)*cos(y/2)",2);
    gsFunctionExpr<>sol2der ("cos(x/2)*sin(y/2)",
                             "cos(x/2)*sin(y/2)",
                             "cos(x/2)*sin(y/2)", 2);
*/

    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der, sol2der);

    // ======= Geometry =========
    std::string string_geo;
    switch(geometry)
    {
        case 0:
            string_geo = "planar/twoPatches/square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 1:
            string_geo = "planar/multiPatches/4_square_diagonal.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 2:
            string_geo = "planar/twoPatches/square_curved.xml";
            numDegree = 0; // 0 == degree 3
            break;
        case 3:
            string_geo = "KirchhoffLoveGeo/geo_fivePatchDiffParam.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 4:
            string_geo = "planar/multiPatches/4_square_curved.xml";
            numDegree = 0; // 2 == degree 3
            break;
        default:
            gsInfo << "No geometry is used! \n";
            break;
    }

    gsFileData<> fd(string_geo);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> multiPatch_init;
    fd.getId(0, multiPatch_init); // id=0: Multipatch domain
    multiPatch_init.computeTopology();


    gsG1OptionList g1OptionList;
    g1OptionList.addInt("p_tilde","Grad",p_tilde);
    g1OptionList.addInt("r_tilde","Reg",r_tilde);
    g1OptionList.addInt("regularity","Regularity of the initial geometry",regularity);
    g1OptionList.addSwitch("plot","Plot in Paraview",plot);
    g1OptionList.addInt("refine","Refinement",numRefine);
    g1OptionList.addInt("degree","Degree",numDegree);
    g1OptionList.addReal("threshold","Threshold",threshold);

    if (localGd)
        gluingData_strategy = gluingData::local;
    if (localEdge)
        g1BasisEdge_strategy = g1BasisEdge::local;
    if (localVertex)
        g1BasisVertex_strategy = g1BasisVertex::local;


    g1OptionList.addInt("gluingData","The strategy for the gluing data",gluingData_strategy);
    g1OptionList.addInt("g1BasisEdge","The strategy for the g1 basis edge",g1BasisEdge_strategy);
    g1OptionList.addInt("g1BasisVertex","The strategy for the g1 basis vertex",g1BasisVertex_strategy);
    g1OptionList.addInt("user", "User ID", user::pascal); // Set the user

    if (g1OptionList.getInt("user") == user::pascal)
        gsInfo << "User is pascal!\n";
    else if (g1OptionList.getInt("user") == user::andrea)
        gsInfo << "User is andrea!\n";

    user::name pascal_id = user::pascal;
    if (pascal_id == g1OptionList.getInt("user"))
        gsInfo << "User is pascal!\n";


    //multiPatch.patch(1).degreeElevate(1,0);
    multiPatch_init.degreeElevate(g1OptionList.getInt("degree"));

    gsVector<real_t> l2Error_vec(loop + 1);
    gsVector<real_t> h1SemiError_vec(loop + 1);
    gsVector<real_t> h2SemiError_vec(loop + 1);
    l2Error_vec.setZero();
    h1SemiError_vec.setZero();
    h2SemiError_vec.setZero();

    gsVector<index_t> num_knots(loop);
    num_knots[0] = g1OptionList.getInt("refine");
    for (index_t i = 1; i < loop; i++)
        num_knots[i] = num_knots[i-1]*2 + 1;

    for (index_t refinement_level = 0; refinement_level < loop; refinement_level++)
    {
        gsMultiPatch<> multiPatch(multiPatch_init);
        multiPatch.uniformRefine_withSameRegularity(num_knots[refinement_level], g1OptionList.getInt("regularity"));

        gsInfo << "###### Level: " << refinement_level << " with " << num_knots[refinement_level] << " inner knots ###### " << "\n";

        gsMultiBasis<> mb(multiPatch);

        omp_set_num_threads(threads);
        omp_set_nested(1);

        gsG1System<real_t> g1System(multiPatch, mb);

        // ########### EDGE FUNCTIONS ###########
        // Interface loop
        for (size_t numInt = 0; numInt < multiPatch.interfaces().size(); numInt++ )
        {
            const boundaryInterface & item = multiPatch.interfaces()[numInt];

            std::string fileName;
            std::string basename = "InterfaceBasisFunctions" + util::to_string(numInt);
            gsParaviewCollection collection(basename);

            gsG1AuxiliaryEdgeMultiplePatches singleInt(multiPatch, item.first().patch, item.second().patch);
            singleInt.computeG1InterfaceBasis(g1OptionList);
            singleInt.deleteBasisFunctions(0,g1System.sizePlusInterface(numInt));
            singleInt.deleteBasisFunctions(1,g1System.sizePlusInterface(numInt));

            for (size_t i = 0; i < singleInt.getSinglePatch(0).getG1Basis().nPatches(); i++)
            {
                gsMultiPatch<> edgeSingleBF;

                edgeSingleBF.addPatch(singleInt.getSinglePatch(0).getG1Basis().patch(i));
                edgeSingleBF.addPatch(singleInt.getSinglePatch(1).getG1Basis().patch(i));

                g1System.insertInterfaceEdge(edgeSingleBF,item,numInt,i);

                if (plot)
                {
                    // First Interface Side
                    fileName = basename + "_0_" + util::to_string(i);
                    gsField<> temp_field(multiPatch.patch(item.first().patch),edgeSingleBF.patch(0));
                    gsWriteParaview(temp_field,fileName,5000);
                    collection.addTimestep(fileName,i,"0.vts");
                    // Second Interface Side
                    fileName = basename + "_1_" + util::to_string(i);
                    gsField<> temp_field_1(multiPatch.patch(item.second().patch),edgeSingleBF.patch(1));
                    gsWriteParaview(temp_field_1,fileName,5000);
                    collection.addTimestep(fileName,i,"0.vts");
                }
            }
            collection.save();
        }
        // Boundaries loop
        for (size_t numBdy = 0; numBdy < multiPatch.boundaries().size(); numBdy++ )
        {
            const patchSide & bit = multiPatch.boundaries()[numBdy];

            std::string fileName;
            std::string basename = "BoundaryBasisFunctions" + util::to_string(numBdy);
            gsParaviewCollection collection(basename);

            gsG1AuxiliaryEdgeMultiplePatches singleBdy(multiPatch, bit.patch);
            singleBdy.computeG1BoundaryBasis(g1OptionList, bit.m_index);
            singleBdy.deleteBasisFunctions(0,g1System.sizePlusBoundary(numBdy));

            for (size_t i = 0; i < singleBdy.getSinglePatch(0).getG1Basis().nPatches(); i++)
            {
                gsMultiPatch<> edgeSingleBF;

                edgeSingleBF.addPatch(singleBdy.getSinglePatch(0).getG1Basis().patch(i));

                g1System.insertBoundaryEdge(edgeSingleBF,bit,numBdy,i);

                if (plot)
                {
                    fileName = basename + "_0_" + util::to_string(i);
                    gsField<> temp_field(multiPatch.patch(bit.patch),edgeSingleBF.patch(0));
                    gsWriteParaview(temp_field,fileName,5000);
                    collection.addTimestep(fileName,i,"0.vts");
                }
            }
            collection.save();
        }

        // Vertices
        for(size_t numVer=0; numVer < multiPatch.vertices().size(); numVer++)
        {
            std::string fileName;
            std::string basename = "VerticesBasisFunctions" + util::to_string(numVer);
            gsParaviewCollection collection(basename);

            std::vector<patchCorner> allcornerLists = multiPatch.vertices()[numVer];
            std::vector<size_t> patchIndex;
            std::vector<size_t> vertIndex;
            for(size_t j = 0; j < allcornerLists.size(); j++)
            {
                patchIndex.push_back(allcornerLists[j].patch);
                vertIndex.push_back(allcornerLists[j].m_index);
            }

            gsG1AuxiliaryVertexMultiplePatches singleVertex(multiPatch, patchIndex, vertIndex);
            singleVertex.computeG1InternalVertexBasis(g1OptionList);
            for (index_t i = 0; i < 6; i++)
            {
                gsMultiPatch<> singleBasisFunction;
                for (size_t np = 0; np < vertIndex.size(); np++)
                {
                    singleBasisFunction.addPatch(singleVertex.getSinglePatch(np).getG1Basis().patch(i));
                    if (plot)
                    {
                        fileName = basename + "_" + util::to_string(np) + "_" + util::to_string(i);
                        gsField<> temp_field(multiPatch.patch(patchIndex[np]),singleBasisFunction.patch(np));
                        gsWriteParaview(temp_field,fileName,5000);
                        collection.addTimestep(fileName,i,"0.vts");
                    }
                }
                g1System.insertVertex(singleBasisFunction,patchIndex,numVer,singleVertex.get_internalDofs(),i);
            }
            collection.save();
        }

        gsBoundaryConditions<> bcInfo, bcInfo2;
        for (gsMultiPatch<>::const_biterator bit = multiPatch.bBegin(); bit != multiPatch.bEnd(); ++bit)
        {
            bcInfo.addCondition( *bit, condition_type::dirichlet, &solVal );
            bcInfo2.addCondition(*bit, condition_type::neumann, &laplace );
        }

        // BiharmonicAssembler
        gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(multiPatch, mb, bcInfo, bcInfo2, source);
        g1BiharmonicAssembler.assemble();
        g1BiharmonicAssembler.computeDirichletDofsL2Proj(g1System); // Compute boundary values (Type 1)

        g1System.finalize(multiPatch,mb,g1BiharmonicAssembler.get_bValue());

        gsInfo << "Solving system... \n";
        gsMatrix<> solVector = g1System.solve(g1BiharmonicAssembler.matrix(), g1BiharmonicAssembler.rhs());
        gsInfo << "Solving finished! \n";

        // construct solution: INTERIOR
        gsMultiPatch<> mpsol;
        g1BiharmonicAssembler.constructSolution(solVector.bottomRows(g1BiharmonicAssembler.matrix().dim().first),mpsol);
        gsField<> solField(multiPatch, mpsol);

        // construct solution: G1 Basis
        std::vector<gsMultiPatch<>> g1Basis;
        g1System.constructG1Solution(solVector,g1Basis, multiPatch);

        g1BiharmonicAssembler.plotParaview(solField, g1Basis);

        omp_set_num_threads(1); // Set to 1 because of memmory problems :/
        omp_set_nested(1);
#pragma omp parallel for
        for (index_t e = 0; e < 4; ++e)
        {
            if (e == 0)
            {
                gsNormL2<real_t> errorL2(solField, solVal, g1Basis);
                errorL2.compute();
                l2Error_vec[refinement_level] = errorL2.value();
            }
            else if (e == 1)
            {
                gsSeminormH1<real_t> errorSemiH1(solField, solVal, g1Basis);
                errorSemiH1.compute();
                h1SemiError_vec[refinement_level] = errorSemiH1.value();
            }
            else if (e == 2)
            {
                gsSeminormH2<real_t> errorSemiH2(solField, solVal, g1Basis);
                errorSemiH2.compute();
                h2SemiError_vec[refinement_level] = errorSemiH2.value();
            }
        }
    }

    gsInfo << "=====================================================================\n";
    if (loop > 1)
    {
        gsMatrix<> rate(loop + 1,3);
        rate.setZero();
        printf("|%-5s|%-14s|%-5s|%-14s|%-5s|%-14s|%-5s\n", "k","L2-error", "Rate", "Semi-H1-error",
            "Rate", "Semi-H2-error", "Rate");
        printf("|%-5s|%-14s|%-5s|%-14s|%-5s|%-14s|%-5s\n", "-----", "--------------", "-----", "--------------",
            "-----", "--------------", "-----");
        printf("|%-5d|%-14.6e|%-5.2f|%-14.6e|%-5.2f|%-14.6e|%-5.2f\n", num_knots[0], l2Error_vec[0],
            rate(0,0),h1SemiError_vec[0], rate(0,1),h2SemiError_vec[0], rate(0,2));
        for (index_t i = 1; i < loop; i++)
        {
            rate(i,0) = log2(l2Error_vec[i-1] / l2Error_vec[i]);
            rate(i,1) = log2(h1SemiError_vec[i-1] / h1SemiError_vec[i]);
            rate(i,2) = log2(h2SemiError_vec[i-1] / h2SemiError_vec[i]);
            printf("|%-5d|%-14.6e|%-5.2f|%-14.6e|%-5.2f|%-14.6e|%-5.2f\n", num_knots[i], l2Error_vec[i],
                rate(i,0),h1SemiError_vec[i], rate(i,1),h2SemiError_vec[i], rate(i,2));
        }
        if (latex)
        {
            printf("%-5d & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\ \n", num_knots[0],
                l2Error_vec[0], rate(0,0),h1SemiError_vec[0], rate(0,1), h2SemiError_vec[0], rate(0,2));
            for (index_t i = 1; i < loop; i++)
            {
                rate(i,0) = log2(l2Error_vec[i-1] / l2Error_vec[i]);
                rate(i,1) = log2(h1SemiError_vec[i-1] / h1SemiError_vec[i]);
                rate(i,2) = log2(h2SemiError_vec[i-1] / h2SemiError_vec[i]);
                printf("%-5d & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\ \n", num_knots[i],
                    l2Error_vec[i], rate(i,0),h1SemiError_vec[i], rate(i,1),h2SemiError_vec[i], rate(i,2));
            }
        }

    }
    else
    {
        gsInfo << "L2 Error: " << l2Error_vec[0] << "\n";
        gsInfo << "H1 Semi-error: " << h1SemiError_vec[0] << "\n";
        gsInfo << "H2 Semi-error: " << h2SemiError_vec[0] << "\n";
    }
    gsInfo << "=====================================================================\n";


} // main