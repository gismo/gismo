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

# include <gsG1Basis/Norm/gsNormL2.h>
# include <gsG1Basis/Norm/gsSeminormH1.h>
# include <gsG1Basis/Norm/gsSeminormH2.h>
# include <gsG1Basis/Norm/gsH1NormWithJump.h>

#include <iostream>

using namespace gismo;

int main(int argc, char *argv[])
{

    gsG1OptionList g1OptionList;
    g1OptionList.initialize(argc, argv);

    g1OptionList.addInt("user","Pascal",user::pascal);

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
    gsFunctionExpr<> source  ("0",2);
    gsFunctionExpr<> laplace ("0",2);
    gsFunctionExpr<> solVal("1",2);
    gsFunctionExpr<>sol1der ("0",
                             "0",2);
    gsFunctionExpr<>sol2der ("0",
                             "0",
                             " 0", 2);
/*
    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
    gsFunctionExpr<>sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                             "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
    gsFunctionExpr<>sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                             "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                             " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)", 2);

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
    index_t numDegree = 0;
    switch(g1OptionList.getInt("geometry"))
    {
        case 0:
            string_geo = "planar/twoPatches/square_diagonal.xml";
            numDegree = 2;
            break;
        case 1:
            string_geo = "planar/twoPatches/square_curved.xml";
            numDegree = 0;
            break;
        case 2:
            string_geo = "planar/twoPatches/2patch_curved.xml";
            numDegree = 0;
            break;
        case 3:
            string_geo = "planar/twoPatches/2patch_C1curved.xml";
            numDegree = 0;
            break;



        case 10:
            string_geo = "planar/multiPatches/4_square_diagonal.xml";
            numDegree = 2;
            break;
        case 11:
            string_geo = "planar/multiPatches/4_square_curved.xml";
            numDegree = 0;
            break;
        case 12:
            string_geo = "planar/multiPatches/3_patch_curved.xml";
            numDegree = 0;
            break;
        case 13:
            string_geo = "planar/multiPatches/6_square_diagonal.xml";
            numDegree = 2;
            break;
        case 14:
            string_geo = "planar/multiPatches/4_patch_linear.xml";
            numDegree = 2;
            break;
        case 15:
            string_geo = "planar/multiPatches/3_patch_corner.xml";
            numDegree = 2;
            break;


        default:
            gsInfo << "No geometry is used! \n";
            break;
    }
    g1OptionList.addInt("degree","Degree", numDegree);

    gsFileData<> fd(string_geo);
    gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> multiPatch_init;
    fd.getId(0, multiPatch_init); // id=0: Multipatch domain
    multiPatch_init.computeTopology();

    gsWriteParaview(multiPatch_init,"geoemtry_init",2000,true);

    //multiPatch.patch(1).degreeElevate(1,0);
    multiPatch_init.degreeElevate(g1OptionList.getInt("degree"));

    gsVector<real_t> l2Error_vec(g1OptionList.getInt("loop") + 1);
    gsVector<real_t> h1SemiError_vec(g1OptionList.getInt("loop") + 1);
    gsVector<real_t> h2SemiError_vec(g1OptionList.getInt("loop") + 1);
    gsMatrix<real_t> h1SemiError_jump_edge(g1OptionList.getInt("loop") + 1, multiPatch_init.interfaces().size());
    gsMatrix<real_t> h1SemiError_jump_vertex(g1OptionList.getInt("loop") + 1, multiPatch_init.interfaces().size());
    gsMatrix<real_t> h1SemiError_jump_all(g1OptionList.getInt("loop") + 1, multiPatch_init.interfaces().size());
    l2Error_vec.setZero();
    h1SemiError_vec.setZero();
    h2SemiError_vec.setZero();
    h1SemiError_jump_edge.setZero();
    h1SemiError_jump_vertex.setZero();
    h1SemiError_jump_all.setZero();

    gsVector<index_t> num_knots(g1OptionList.getInt("loop"));
    num_knots[0] = g1OptionList.getInt("numRefine");
    for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
        num_knots[i] = num_knots[i-1]*2 + 1;

    for (index_t refinement_level = 0; refinement_level < g1OptionList.getInt("loop"); refinement_level++)
    {
        real_t old_lambda = g1OptionList.getReal("lambda");
        g1OptionList.setReal("lambda",old_lambda/4);

        gsMultiPatch<> multiPatch(multiPatch_init);
        multiPatch.uniformRefine_withSameRegularity(num_knots[refinement_level], g1OptionList.getInt("regularity"));

        gsInfo << "###### Level: " << refinement_level << " with " << num_knots[refinement_level] << " inner knots ###### " << "\n";

        gsMultiBasis<> mb(multiPatch);

#ifdef _OPENMP
        omp_set_num_threads( g1OptionList.getInt("threads"));
        omp_set_nested(1);
#endif
        gsG1System<real_t> g1System(multiPatch, mb, g1OptionList.getSwitch("neumann"));

        // ########### EDGE FUNCTIONS ###########
        // Interface loop
        gsInfo << "Computing Interface basis functions ... \n";
        for (size_t numInt = 0; numInt < multiPatch.interfaces().size(); numInt++ )
        {
            const boundaryInterface & item = multiPatch.interfaces()[numInt];

            std::string fileName;
            std::string basename = "InterfaceBasisFunctions" + util::to_string(numInt);
            gsParaviewCollection collection(basename);

            gsG1AuxiliaryEdgeMultiplePatches singleInt(multiPatch, item.first().patch, item.second().patch);
            singleInt.computeG1InterfaceBasis(g1OptionList);

            for (size_t i = 0; i < singleInt.getSinglePatch(0).getG1Basis().nPatches(); i++)
            {
                gsMultiPatch<> edgeSingleBF;

                edgeSingleBF.addPatch(singleInt.getSinglePatch(0).getG1Basis().patch(i));
                edgeSingleBF.addPatch(singleInt.getSinglePatch(1).getG1Basis().patch(i));

                g1System.insertInterfaceEdge(edgeSingleBF,item,numInt,i);

                if (g1OptionList.getSwitch("plot"))
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
        gsInfo << "Computing Boundary basis functions ... \n";
        for (size_t numBdy = 0; numBdy < multiPatch.boundaries().size(); numBdy++ )
        {
            const patchSide & bit = multiPatch.boundaries()[numBdy];

            std::string fileName;
            std::string basename = "BoundaryBasisFunctions" + util::to_string(numBdy);
            gsParaviewCollection collection(basename);

            gsG1AuxiliaryEdgeMultiplePatches singleBdy(multiPatch, bit.patch);
            singleBdy.computeG1BoundaryBasis(g1OptionList, bit.m_index);

            for (size_t i = 0; i < singleBdy.getSinglePatch(0).getG1Basis().nPatches(); i++)
            {
                gsMultiPatch<> edgeSingleBF;

                edgeSingleBF.addPatch(singleBdy.getSinglePatch(0).getG1Basis().patch(i));

                g1System.insertBoundaryEdge(edgeSingleBF,bit,numBdy,i);

                if (g1OptionList.getSwitch("plot"))
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
        gsInfo << "Computing Vertex basis functions ... \n";
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
                    if (g1OptionList.getSwitch("plot"))
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
            if (!g1OptionList.getSwitch("neumann"))
                bcInfo2.addCondition( *bit, condition_type::laplace, &laplace);
            else
                bcInfo2.addCondition(*bit, condition_type::neumann, &sol1der );
        }

        // BiharmonicAssembler
        gsInfo << "Computing Internal basis functions ... \n";
        gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(multiPatch, mb, bcInfo, bcInfo2, source);
        g1BiharmonicAssembler.assemble();
        gsInfo << "Computing Boundary data ... \n";
        if (!g1OptionList.getSwitch("neumann"))
            g1BiharmonicAssembler.computeDirichletDofsL2Proj(g1System); // Compute boundary values with laplace
        else
            g1BiharmonicAssembler.computeDirichletAndNeumannDofsL2Proj(g1System); // Compute boundary values with neumann

        g1System.finalize(multiPatch,mb,g1BiharmonicAssembler.get_bValue());

        gsInfo << "Solving system... \n";
        gsMatrix<> solVector = g1System.solve(g1BiharmonicAssembler.matrix(), g1BiharmonicAssembler.rhs());
        gsInfo << "Solving finished! \n";

        if (g1OptionList.getSwitch("plot"))
        {
            // construct solution: INTERIOR
            gsMultiPatch<> mpsol;
            g1BiharmonicAssembler.constructSolution(solVector.bottomRows(g1BiharmonicAssembler.matrix().dim().first),mpsol);
            gsField<> solField(multiPatch, mpsol);

            // construct solution for plotting
            std::vector<gsMultiPatch<>> g1Basis;
            g1System.constructG1Solution(solVector,g1Basis, multiPatch);

            g1BiharmonicAssembler.plotParaview(solField, g1Basis);
        }

        // construct solution: G1 Basis
        gsSparseMatrix<real_t> Sol_sparse;
        g1System.constructSparseG1Solution(solVector,Sol_sparse);

#ifdef _OPENMP
        // omp_set_num_threads(g1OptionList.getInt("threads"));
        omp_set_num_threads(1);
        omp_set_nested(1);
#endif

#pragma omp parallel for
        for (index_t e = 0; e < 6; ++e)
        {
            if (e == 0)
            {
                gsNormL2<real_t> errorL2(multiPatch, Sol_sparse, solVal);
                errorL2.compute(g1System.get_numBasisFunctions(), true);
                l2Error_vec[refinement_level] = errorL2.value();
            }

            else if (e == 1)
            {
                gsSeminormH1<real_t> errorSemiH1(multiPatch, Sol_sparse, solVal);
                errorSemiH1.compute(g1System.get_numBasisFunctions());
                h1SemiError_vec[refinement_level] = errorSemiH1.value();
            }
            else if (e == 2)
            {
                gsSeminormH2<real_t> errorSemiH2(multiPatch, Sol_sparse, solVal);
                errorSemiH2.compute(g1System.get_numBasisFunctions());
                h2SemiError_vec[refinement_level] = errorSemiH2.value();
            }
            else if (e == 3)
            {
                gsH1NormWithJump<real_t> errorJump(multiPatch, Sol_sparse);
                errorJump.compute(g1System.get_numBasisFunctions(), g1System.get_numInterfaceFunctions(), "edge");
                h1SemiError_jump_edge.row(refinement_level) = errorJump.value().transpose();
            }
            else if (e == 4)
            {
                gsH1NormWithJump<real_t> errorJump(multiPatch, Sol_sparse);
                errorJump.compute(g1System.get_numBasisFunctions(), g1System.get_numVertexFunctions(), "vertex");
                h1SemiError_jump_vertex.row(refinement_level) = errorJump.value().transpose();
            }
            else if (e == 5)
            {
                gsH1NormWithJump<real_t> errorJump(multiPatch, Sol_sparse);
                errorJump.compute(g1System.get_numBasisFunctions(), g1System.get_numVertexFunctions(), "all");
                h1SemiError_jump_all.row(refinement_level) = errorJump.value().transpose();
            }
        }
    }

    for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
    {
        h2SemiError_vec[i] = math::sqrt(h2SemiError_vec[i]*h2SemiError_vec(i) +
            h1SemiError_vec[i]*h1SemiError_vec[i] + l2Error_vec[i]*l2Error_vec[i]);
        h1SemiError_vec[i] = math::sqrt(h1SemiError_vec[i]*h1SemiError_vec[i] +
            l2Error_vec[i]*l2Error_vec[i]);
    }

    if (g1OptionList.getInt("loop") > 1)
    {
        gsInfo << "=====================================================================\n";

        gsMatrix<> rate(g1OptionList.getInt("loop") + 1,3);
        rate.setZero();
        printf("|%-5s|%-14s|%-5s|%-14s|%-5s|%-14s|%-5s\n", "k","L2-error", "Rate", "H1-error",
            "Rate", "H2-error", "Rate");
        printf("|%-5s|%-14s|%-5s|%-14s|%-5s|%-14s|%-5s\n", "-----", "--------------", "-----", "--------------",
            "-----", "--------------", "-----");
        printf("|%-5d|%-14.6e|%-5.2f|%-14.6e|%-5.2f|%-14.6e|%-5.2f\n", num_knots[0], l2Error_vec[0],
            rate(0,0),h1SemiError_vec[0], rate(0,1),h2SemiError_vec[0], rate(0,2));
        for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
        {
            rate(i,0) = log2(l2Error_vec[i-1] / l2Error_vec[i]);
            rate(i,1) = log2(h1SemiError_vec[i-1] / h1SemiError_vec[i]);
            rate(i,2) = log2(h2SemiError_vec[i-1] / h2SemiError_vec[i]);

            printf("|%-5d|%-14.6e|%-5.2f|%-14.6e|%-5.2f|%-14.6e|%-5.2f\n", num_knots[i], l2Error_vec[i],
                rate(i,0),h1SemiError_vec[i], rate(i,1),h2SemiError_vec[i], rate(i,2));
        }
        if (g1OptionList.getSwitch("latex"))
        {
            printf("%-5d & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\ \n", num_knots[0],
                l2Error_vec[0], rate(0,0),h1SemiError_vec[0], rate(0,1), h2SemiError_vec[0], rate(0,2));
            for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
            {
                printf("%-5d & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\ \n", num_knots[i],
                    l2Error_vec[i], rate(i,0),h1SemiError_vec[i], rate(i,1),h2SemiError_vec[i], rate(i,2));
            }
        }
        gsInfo << "=====================================================================\n\n";

        // JUMP
        rate.setZero(g1OptionList.getInt("loop") + 1,multiPatch_init.interfaces().size());
        gsMatrix<> rate_vertex(g1OptionList.getInt("loop") + 1,multiPatch_init.interfaces().size());
        rate_vertex.setZero();
        gsMatrix<> rate_all(g1OptionList.getInt("loop") + 1,multiPatch_init.interfaces().size());
        rate_all.setZero();
        gsInfo << "======";
        for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
            gsInfo << "===============================================================";

        gsInfo << "\n";
        printf("|%-5s", "k");
        for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
            printf("|%-14s|%-5s|%-14s|%-5s|%-14s|%-5s", ("Single E. " + std::to_string(i)).c_str(), "Rate",
                ("Single V. " + std::to_string(i)).c_str(), "Rate", ("IFace " + std::to_string(i)).c_str(), "Rate");
        gsInfo << "\n";
        printf("|%-5s","-----");
        for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
            printf("|%-14s|%-5s|%-14s|%-5s|%-14s|%-5s", "--------------", "-----", "--------------", "-----", "--------------", "-----");
        gsInfo << "\n";

        printf("|%-5d",num_knots[0]);
        for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
            printf("|%-14.6e|%-5.2f|%-14.6e|%-5.2f|%-14.6e|%-5.2f", h1SemiError_jump_edge(0,i), rate(0,i),
                h1SemiError_jump_vertex(0,i), rate_vertex(0,i), h1SemiError_jump_all(0,i), rate_all(0,i));
        printf("\n");

        for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
        {
            printf("|%-5d",num_knots[i]);
            for (size_t j = 0; j < multiPatch_init.interfaces().size(); j++)
            {
                rate(i,j) = log2(h1SemiError_jump_edge(i-1,j) / h1SemiError_jump_edge(i,j));
                rate_vertex(i,j) = log2(h1SemiError_jump_vertex(i-1,j) / h1SemiError_jump_vertex(i,j));
                rate_all(i,j) = log2(h1SemiError_jump_all(i-1,j) / h1SemiError_jump_all(i,j));
                printf("|%-14.6e|%-5.2f|%-14.6e|%-5.2f|%-14.6e|%-5.2f", h1SemiError_jump_edge(i,j), rate(i,j),
                    h1SemiError_jump_vertex(i,j), rate_vertex(i,j), h1SemiError_jump_all(i,j), rate_all(i,j));
            }
            printf("\n");
        }

        gsInfo << "======";
        for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
            gsInfo << "===============================================================";

        gsInfo << "\n";

        if (g1OptionList.getSwitch("latex"))
        {
            for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
                printf("%-5d & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\", num_knots[0],
                   h1SemiError_jump_edge(0,i), rate(0,i),h1SemiError_jump_vertex(0,i), rate_vertex(0,i), h1SemiError_jump_all(0,i), rate_all(0,i));
            printf("\n");
            for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
            {
                printf("%-5d & ",num_knots[i]);
                for (size_t j = 0; j < multiPatch_init.interfaces().size(); j++)
                {
                    printf("%-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\",
                           h1SemiError_jump_edge(i,j),
                           rate(i, j),
                           h1SemiError_jump_vertex(i,j),
                           rate_vertex(i, j),
                           h1SemiError_jump_all(i,j),
                           rate_all(i, j));
                }
                printf("\n");
            }
        }
    }
    else
    {
        gsInfo << "=====================================================================\n";
        gsInfo << "L2 Error: " << l2Error_vec[0] << "\n";
        gsInfo << "H1 Semi-error: " << h1SemiError_vec[0] << "\n";
        gsInfo << "H2 Semi-error: " << h2SemiError_vec[0] << "\n";
        gsInfo << "Jump error Edge: " << h1SemiError_jump_edge.row(0) << "\n";
        gsInfo << "Jump error Vertex: " << h1SemiError_jump_vertex.row(0) << "\n";
        gsInfo << "Jump error all: " << h1SemiError_jump_all.row(0) << "\n";
        gsInfo << "=====================================================================\n";

    }


} // main
