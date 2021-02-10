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
# include <gsG1Basis/Norm/gsH1NormProof.h>

#include <iostream>
#include <typeinfo>

using namespace gismo;

int main(int argc, char *argv[])
{
    gsG1OptionList g1OptionList;
    g1OptionList.initialize(argc, argv);

    g1OptionList.addInt("user","Pascal",user::pascal);
    g1OptionList.setSwitch("twoPatch",true);

    // ======= Solution =========
/*   gsFunctionExpr<> source  ("4096*pi*pi*pi*pi*(4*cos(8*pi*x)*cos(8*pi*y) - cos(8*pi*x) - cos(8*pi*y))",2);
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
    gsFunctionExpr<> source  ("0",2);
    gsFunctionExpr<> laplace ("0",2);
    gsFunctionExpr<> solVal("1",2);
    gsFunctionExpr<>sol1der ("0",
                             "0",2);
    gsFunctionExpr<>sol2der ("0",
                             "0",
                             "0", 2);

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
        case 4:
            string_geo = "planar/twoPatches/square_complex_bent.xml";
            numDegree = 0;
            break;
        case 5:
            string_geo = "planar/twoPatches/2patch_C1curved_complex.xml";
            numDegree = 0;
            break;
        case 6:
            string_geo = "planar/twoPatches/2patch_puzzle.xml";
            numDegree = 0;
            break;
        case 7:
            string_geo = "planar/twoPatches/square_diagonal_partial_matching.xml";
            numDegree = 2;
            break;
        case 8:
            string_geo = "planar/twoPatches/square_cuttedCornerBSpline.xml";
            numDegree = 0;
            break;
        case 9:
            string_geo = "planar/twoPatches/square_curved_large.xml";
            numDegree = 0;
            break;
        case 10:
            string_geo = "planar/twoPatches/square_curved_shift.xml";
            numDegree = 0;
            break;


        case 11:
            string_geo = "planar/multiPatches/4_square_curved.xml";
            numDegree = 0;
            break;
        case 12:
            string_geo = "planar/multiPatches/3_patch_curved.xml";
            numDegree = 0;
            break;

        default:
            //gsInfo << "No geometry is used! \n";
            break;
    }
    g1OptionList.addInt("degree","Degree", numDegree);

    gsFileData<> fd(string_geo);
    //gsInfo << "Loaded file "<< fd.lastPath() <<"\n";

    gsMultiPatch<> multiPatch_init;
    fd.getId(0, multiPatch_init); // id=0: Multipatch domain
    multiPatch_init.computeTopology();

/*   TODO NURBS
    gsMultiBasis<> mb_test(multiPatch_init, false);

    gsMatrix<> point(2,1);
    point << 0.1, 0.1;
    gsInfo << mb_test.basis(0) << "\n";
    gsInfo << mb_test.basis(1).eval(point) << "\n";
    gsInfo << multiPatch_init.basis(1).eval(point) << "\n";
*/
    multiPatch_init.degreeElevate(g1OptionList.getInt("degree") + g1OptionList.getInt("P_geo"));

    //std::vector<int> mul={multiPatch_init.patch(1).degree(1) - g1OptionList.getInt("regularity"), multiPatch_init.patch(1).degree(1) - g1OptionList.getInt("regularity")};
    //multiPatch_init.patch(1).uniformRefine(1,mul);
    //multiPatch_init.patch(1).degreeElevate(1);

    gsVector<real_t> l2Error_vec(g1OptionList.getInt("loop") + 1);
    gsVector<real_t> h1SemiError_vec(g1OptionList.getInt("loop") + 1);
    gsVector<real_t> h2SemiError_vec(g1OptionList.getInt("loop") + 1);
    gsMatrix<real_t> h1SemiError_jump_edge(g1OptionList.getInt("loop") + 1, multiPatch_init.interfaces().size());
    gsMatrix<real_t> h1SemiError_jump_proof(g1OptionList.getInt("loop") + 1, multiPatch_init.interfaces().size());
    l2Error_vec.setZero();
    h1SemiError_vec.setZero();
    h2SemiError_vec.setZero();
    h1SemiError_jump_edge.setZero();
    h1SemiError_jump_proof.setZero();

    gsVector<index_t> num_knots(g1OptionList.getInt("loop"));
    num_knots[0] = g1OptionList.getInt("numRefine");
    for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
        num_knots[i] = num_knots[i-1]*2 + 1;

    gsVector<real_t> mesh_size(g1OptionList.getInt("loop"));
    gsVector<index_t> dofs_size(g1OptionList.getInt("loop"));

    gsMatrix<real_t> l2Error_basisfunction(g1OptionList.getInt("loop") + 1, multiPatch_init.interfaces().size());
    l2Error_basisfunction.setZero();

    for (index_t refinement_level = 0; refinement_level < g1OptionList.getInt("loop"); refinement_level++)
    {
        // Proof
        if (g1OptionList.getSwitch("h1projectionProof"))
            g1OptionList.setSwitch("isogeometric", false);
        std::vector<gsSparseMatrix<real_t>> sol_proof;
        std::vector<std::vector<gsMultiBasis<>>> mb_proof;
        std::vector<gsG1System<real_t>> g1System_proof;
        std::vector<gsMultiPatch<>> multiPatch_proof;
        for (index_t projectionProof = 0; projectionProof < (g1OptionList.getSwitch("h1projectionProof") ? 2 : 1); projectionProof++)
        {
            if (g1OptionList.getSwitch("h1projectionProof") && projectionProof == 1)
                g1OptionList.setSwitch("isogeometric", !g1OptionList.getSwitch("isogeometric"));

            // Proof end

            gsMultiPatch<> multiPatch(multiPatch_init);

            //if ((multiPatch_init.basis(0).maxDegree() - g1OptionList.getInt("regularity")) > 1) // P - R > 1
            multiPatch.uniformRefine_withSameRegularity(num_knots[refinement_level], g1OptionList.getInt("regularity"));
            //else
            //    multiPatch.uniformRefine_withDifferentRegularity(num_knots[refinement_level], g1OptionList.getInt("regularity"));

            //gsTensorBSpline<2,real_t> bsp_temp = dynamic_cast<gsTensorBSpline<2,real_t>&>(multiPatch.patch(0));
            //gsInfo << "knotvector: " << bsp_temp.knots(0).asMatrix() << "\n";

            gsWriteParaview(multiPatch, "geometry_init", 2000, true);

            //gsInfo << "###### Level: " << refinement_level << " with " << num_knots[refinement_level] << " inner knots ###### " << "\n";

            mesh_size[refinement_level] =
                multiPatch.basis(0).getMinCellLength() > multiPatch.basis(1).getMinCellLength() ?
                multiPatch.basis(1).getMinCellLength() : multiPatch.basis(0).getMinCellLength();

            std::vector<gsMultiBasis<>> mb;
            mb.push_back(gsMultiBasis<>(multiPatch));
            if (!g1OptionList.getSwitch("isogeometric"))
            {
                gsMultiPatch<> multiPatch_temp(multiPatch_init);

                //std::vector<int> mul={0, multiPatch_init.patch(0).degree(1) - g1OptionList.getInt("regularity")};
                //multiPatch_temp.patch(0).uniformRefine(1,mul);

                multiPatch_temp.patch(0).degreeElevate(g1OptionList.getInt("p_tilde") - 1, 1);
                multiPatch_temp.patch(1).degreeElevate(g1OptionList.getInt("p_tilde") - 1, 1);

                // Degree in Patch 0 and 1 the same!!! TODO
                std::vector<std::vector<int>> patch_mul;
                std::vector<int> single_mul;

                single_mul.push_back(g1OptionList.getInt("regularity")); // u
                single_mul.push_back(math::min(g1OptionList.getInt("regularity"), math::min(g1OptionList.getInt("r_tilde"), multiPatch_init.patch(0).degree(1)-2))); // v
                patch_mul.push_back(single_mul);

                single_mul.clear();
                single_mul.push_back(g1OptionList.getInt("regularity")); // u
                single_mul.push_back(math::min(g1OptionList.getInt("regularity"), math::min(g1OptionList.getInt("r_tilde"),multiPatch_init.patch(0).degree(1)-2))); // v
                patch_mul.push_back(single_mul);

                multiPatch_temp.uniformRefine_withSameRegularity(num_knots[refinement_level], patch_mul);

                gsMultiBasis<> mb_g1(multiPatch_temp);
                /*
                if (g1OptionList.getInt("regularity") + g1OptionList.getInt("P_geo") >= g1OptionList.getInt("r_tilde"))
                    mb_g1.degreeIncrease(g1OptionList.getInt("p_tilde") -1,1); // (keeping the multiplicity)
                else
                    mb_g1.degreeElevate(g1OptionList.getInt("p_tilde") -1,1); // (keeping the smoothness)
                */

                //gsInfo << "Basis G1: " << mb_g1.basis(0) << "\n";
                //gsInfo << "Basis G1: " << mb_g1.basis(1) << "\n";
                mb.push_back(mb_g1);

                if (g1OptionList.getSwitch("info"))
                {
                    gsBSplineBasis<> basis_bspline = dynamic_cast<gsBSplineBasis<real_t> &>(mb[0].basis(0).component(1));
                    gsBSplineBasis<> basis_bspline_g1 = dynamic_cast<gsBSplineBasis<real_t> &>(mb[1].basis(0).component(1));
                    gsInfo << "Basis: " << basis_bspline.knots().asMatrix() << "\n";
                    gsInfo << "Basis G1: " << basis_bspline_g1.knots().asMatrix() << "\n";
                }
            }

            //mb.degreeIncrease(1,0);
            //mb.uniformRefine(num_knots[refinement_level], 3-g1OptionList.getInt("regularity"));

            //gsBSplineBasis<> basis_bspline = dynamic_cast<gsBSplineBasis<real_t> &>(mb[0].basis(0).component(1));
            //gsInfo << "Basis: " << basis_bspline.knots().asMatrix() << "\n";
            //gsInfo << "Basis: " << mb.basis(1) << "\n";
            //gsInfo << "Basis: " << mb[0].basis(0) << "\n";
            //gsInfo << "Basis: " << mb[0].basis(1) << "\n";

#ifdef _OPENMP
            omp_set_num_threads(g1OptionList.getInt("threads"));
            omp_set_nested(1);
#endif
            gsG1System<real_t> g1System(multiPatch,
                                        mb,
                                        g1OptionList.getSwitch("neumann"),
                                        g1OptionList.getSwitch("twoPatch"),
                                        g1OptionList.getSwitch("isogeometric"));

            gsStopwatch clock;

            // ########### EDGE FUNCTIONS ###########
            // Interface loop
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Computing Interface basis functions ... \n";
            for (size_t numInt = 0; numInt < multiPatch.interfaces().size(); numInt++)
            {


                const boundaryInterface & item = multiPatch.interfaces()[numInt];

                std::string fileName;
                std::string basename = "InterfaceBasisFunctions" + util::to_string(numInt);
                gsParaviewCollection collection(basename);

                gsG1AuxiliaryEdgeMultiplePatches singleInt(multiPatch,
                                                           mb[g1OptionList.getSwitch("isogeometric") ? 0 : 1],
                                                           item.first().patch,
                                                           item.second().patch);
                singleInt.computeG1InterfaceBasis(g1OptionList);

                l2Error_basisfunction(refinement_level,0) = singleInt.getError();
                //gsInfo << "interfaces: " << singleInt.getSinglePatch(0).getG1Basis().nPatches() << "\n";

                if (g1OptionList.getSwitch("info"))
                    gsInfo << "Adding to the system \n";

                for (size_t i = 0; i < singleInt.getSinglePatch(0).getG1Basis().nPatches(); i++)
                {

                    gsMultiPatch<> edgeSingleBF;

                    edgeSingleBF.addPatch(singleInt.getSinglePatch(0).getG1Basis().patch(i));
                    edgeSingleBF.addPatch(singleInt.getSinglePatch(1).getG1Basis().patch(i));

                    g1System.insertInterfaceEdge(edgeSingleBF, item, numInt, i);

                    if (g1OptionList.getSwitch("plot"))
                    {
                        // First Interface Side
                        fileName = basename + "_0_" + util::to_string(i);
                        gsField<> temp_field(multiPatch.patch(item.first().patch), edgeSingleBF.patch(0));
                        gsWriteParaview(temp_field, fileName, 5000);
                        collection.addTimestep(fileName, i, "0.vts");
                        // Second Interface Side
                        fileName = basename + "_1_" + util::to_string(i);
                        gsField<> temp_field_1(multiPatch.patch(item.second().patch), edgeSingleBF.patch(1));
                        gsWriteParaview(temp_field_1, fileName, 5000);
                        collection.addTimestep(fileName, i, "0.vts");
                    }
                }

                collection.save();

                real_t dimU = singleInt.getSinglePatch(0).getG1Basis().patch(0).basis().component(0).size();

                gsMatrix<> coefs_interface(2*dimU,3);
                coefs_interface.setZero();

                for (size_t i = 0; i < dimU; i++)
                {
                    if (singleInt.getSinglePatch(0).getG1Basis().patch(0).coef(i,0)*singleInt.getSinglePatch(0).getG1Basis().patch(0).coef(i,0) > 1e-25)
                        coefs_interface(i, 0) = singleInt.getSinglePatch(0).getG1Basis().patch(0).coef(i,0);
                    if (singleInt.getSinglePatch(0).getG1Basis().patch(1).coef(i,0)*singleInt.getSinglePatch(0).getG1Basis().patch(1).coef(i,0) > 1e-25)
                        coefs_interface(i, 1) = singleInt.getSinglePatch(0).getG1Basis().patch(1).coef(i,0);
                    if (singleInt.getSinglePatch(0).getG1Basis().patch(singleInt.getPlus()).coef(i,0)*singleInt.getSinglePatch(0).getG1Basis().patch(singleInt.getPlus()).coef(i,0) > 1e-25)
                        coefs_interface(i, 2) = singleInt.getSinglePatch(0).getG1Basis().patch(singleInt.getPlus()).coef(i,0);

                    if (singleInt.getSinglePatch(1).getG1Basis().patch(0).coef(i,0)*singleInt.getSinglePatch(1).getG1Basis().patch(0).coef(i,0) > 1e-25)
                        coefs_interface(dimU + i, 0) = singleInt.getSinglePatch(1).getG1Basis().patch(0).coef(i,0);
                    if (singleInt.getSinglePatch(1).getG1Basis().patch(1).coef(i,0)*singleInt.getSinglePatch(1).getG1Basis().patch(1).coef(i,0) > 1e-25)
                        coefs_interface(dimU + i, 1) = singleInt.getSinglePatch(1).getG1Basis().patch(1).coef(i,0);
                    if (singleInt.getSinglePatch(1).getG1Basis().patch(singleInt.getPlus()).coef(i,0)*singleInt.getSinglePatch(1).getG1Basis().patch(singleInt.getPlus()).coef(i,0) > 1e-25)
                        coefs_interface(dimU + i, 2) = singleInt.getSinglePatch(1).getG1Basis().patch(singleInt.getPlus()).coef(i,0);
                }

                //gsInfo << "MATRIX: " << coefs_interface << "\n";

                Eigen::FullPivLU<gsMatrix<>> SmallLU(coefs_interface);
                SmallLU.setThreshold(1e-5);
                gsInfo << "Kernel: " << SmallLU.kernel() << "\n";
                gsInfo << "Kernel 2: " << math::pow(mesh_size[refinement_level],(g1OptionList.getInt("p_tilde") +1)) << "\n";

            }
            //gsInfo << "Done. Interface computing time is: " << clock.stop() << "\n";
            clock.restart();
            // Boundaries loop
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Computing Boundary basis functions ... \n";
            for (size_t numBdy = 0; numBdy < multiPatch.boundaries().size(); numBdy++)
            {
                const patchSide & bit = multiPatch.boundaries()[numBdy];

                index_t dir = multiPatch.boundaries()[numBdy].m_index < 3 ? 1 : 0;
                gsBSplineBasis<> basis_edge =
                    dynamic_cast<gsBSplineBasis<> &>(mb[0].basis(multiPatch.boundaries()[numBdy].patch)
                                                          .component(dir)); // 0 -> u, 1 -> v

                std::string fileName;
                std::string basename = "BoundaryBasisFunctions" + util::to_string(numBdy);
                gsParaviewCollection collection(basename);

                gsG1AuxiliaryEdgeMultiplePatches singleBdy(multiPatch, mb[0], bit.patch);
                singleBdy.computeG1BoundaryBasis(g1OptionList, bit.m_index);

                for (size_t i = 0; i < singleBdy.getSinglePatch(0).getG1Basis().nPatches(); i++)
                {
                    gsMultiPatch<> edgeSingleBF;

                    edgeSingleBF.addPatch(singleBdy.getSinglePatch(0).getG1Basis().patch(i));

                    g1System.insertBoundaryEdge(edgeSingleBF, bit, numBdy, i);

                    if (g1OptionList.getSwitch("plot"))
                    {
                        fileName = basename + "_0_" + util::to_string(i);
                        gsField<> temp_field(multiPatch.patch(bit.patch), edgeSingleBF.patch(0));
                        gsWriteParaview(temp_field, fileName, 5000);
                        collection.addTimestep(fileName, i, "0.vts");
                    }
                }
                collection.save();
            }
            //gsInfo << clock.stop() << "\n";
            clock.restart();
            // Vertices
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Computing Vertex basis functions ... \n";
            for (size_t numVer = 0; numVer < multiPatch.vertices().size(); numVer++)
            {
                std::string fileName;
                std::string basename = "VerticesBasisFunctions" + util::to_string(numVer);
                gsParaviewCollection collection(basename);

                std::vector<patchCorner> allcornerLists = multiPatch.vertices()[numVer];
                std::vector<size_t> patchIndex;
                std::vector<size_t> vertIndex;
                for (size_t j = 0; j < allcornerLists.size(); j++)
                {
                    patchIndex.push_back(allcornerLists[j].patch);
                    vertIndex.push_back(allcornerLists[j].m_index);
                }
                if (patchIndex.size() == 1)
                {
                    gsG1AuxiliaryVertexMultiplePatches singleVertex(multiPatch, mb[0], patchIndex, vertIndex);
                    singleVertex.computeG1InternalVertexBasis(g1OptionList);

                    for (index_t i = 0; i < 4; i++)
                    {
                        gsMultiPatch<> singleBasisFunction;
                        for (size_t np = 0; np < vertIndex.size(); np++)
                        {
                            singleBasisFunction.addPatch(singleVertex.getSinglePatch(np).getG1Basis().patch(i));
                            if (g1OptionList.getSwitch("plot"))
                            {
                                fileName = basename + "_" + util::to_string(np) + "_" + util::to_string(i);
                                gsField<> temp_field(multiPatch.patch(patchIndex[np]), singleBasisFunction.patch(np));
                                gsWriteParaview(temp_field, fileName, 5000);
                                collection.addTimestep(fileName, i, "0.vts");
                            }
                        }

                        g1System
                            .insertVertex(singleBasisFunction, patchIndex, numVer, singleVertex.get_internalDofs(), i);

                    }
                }
                collection.save();
            }
            //gsInfo << clock.stop() << "\n";
            clock.restart();

            gsBoundaryConditions<> bcInfo, bcInfo2;
            for (gsMultiPatch<>::const_biterator bit = multiPatch.bBegin(); bit != multiPatch.bEnd(); ++bit)
            {
                bcInfo.addCondition(*bit, condition_type::dirichlet, &solVal);
                bcInfo2.addCondition(*bit, condition_type::laplace, &laplace);
            }

            // BiharmonicAssembler
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Computing Internal basis functions ... \n";
            gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(multiPatch, mb[0], bcInfo, bcInfo2, source);
            g1BiharmonicAssembler.assemble();

            clock.restart();
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Computing Boundary dofs ... \n";
            g1BiharmonicAssembler
                .computeDirichletDofsL2Proj(mb,g1System, g1OptionList.getSwitch("isogeometric")); // Compute boundary values (Type 1) // maybe too much memmory!!!
            //gsInfo << clock.stop() << "\n";
            clock.restart();

            if (!g1OptionList.getSwitch("isogeometric"))
            {
                // For interface basis
                gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler_g22(multiPatch, mb[1], bcInfo, bcInfo2, source);
                g1BiharmonicAssembler_g22.assemble();
                //g1BiharmonicAssembler_g22.computeDirichletDofsL2Proj(mb,g1System, false);
/*
            // Mixed
            gsBSplineBasis<> basis_1u = dynamic_cast<gsBSplineBasis<> &>(mb[0].basis(0).component(0)); // 0 -> u, 1 -> v
            gsBSplineBasis<> basis_1v = dynamic_cast<gsBSplineBasis<> &>(mb[1].basis(0).component(1)); // 0 -> u, 1 -> v

            gsTensorBSplineBasis<2,real_t> Tbasis_1(basis_1u.knots(),basis_1v.knots());

            gsBSplineBasis<> basis_2u = dynamic_cast<gsBSplineBasis<> &>(mb[0].basis(1).component(0)); // 0 -> u, 1 -> v
            gsBSplineBasis<> basis_2v = dynamic_cast<gsBSplineBasis<> &>(mb[1].basis(1).component(1)); // 0 -> u, 1 -> v

            gsTensorBSplineBasis<2,real_t> Tbasis_2(basis_2u.knots(),basis_2v.knots());

            gsMultiBasis<> mb_mixed;
            mb_mixed.addBasis(&Tbasis_1);
            mb_mixed.addBasis(&Tbasis_2);

            gsMultiBasis<> mb_mixed2;
            mb_mixed2.addBasis(&Tbasis_2);
            mb_mixed2.addBasis(&Tbasis_1);
*/

                std::vector<std::vector<int>> patch_mul;
                gsMultiPatch<> multiPatch_temp(multiPatch_init);
                //multiPatch_temp.patch(0).degreeElevate(g1OptionList.getInt("p_tilde") -1,1);
                multiPatch_temp.patch(1).degreeElevate(g1OptionList.getInt("p_tilde") - 1, 1);

                std::vector<int> single_mul;

                single_mul.push_back(g1OptionList.getInt("regularity")); // u
                single_mul.push_back(g1OptionList.getInt("regularity")); // v
                patch_mul.push_back(single_mul);

                single_mul.clear();
                single_mul.push_back(g1OptionList.getInt("regularity")); // u
                single_mul.push_back(math::min(g1OptionList.getInt("regularity"), math::min(g1OptionList.getInt("r_tilde"),multiPatch_init.patch(0).degree(1)-2))); // v

                patch_mul.push_back(single_mul);
                multiPatch_temp.uniformRefine_withSameRegularity(num_knots[refinement_level], patch_mul);


                gsMultiBasis<> mb_g1(multiPatch_temp);

                //gsInfo << "Basis: " << mb_g1.basis(0) << "\n";
                //gsInfo << "Basis2 : " << mb_g1.basis(1) << "\n";
/*
            if (g1OptionList.getInt("regularity") + g1OptionList.getInt("P_geo") >= g1OptionList.getInt("r_tilde"))
                mb_g1.basis(1).degreeIncrease(g1OptionList.getInt("p_tilde") -1,1);
            else
                mb_g1.basis(1).degreeElevate(g1OptionList.getInt("p_tilde") -1,1);
*/
                //mb_g1.basis(0).degreeElevate(g1OptionList.getInt("p_tilde") -1,1);
                //mb_g1.basis(1).degreeElevate(g1OptionList.getInt("p_tilde") -1,1);
                if (g1OptionList.getSwitch("info"))
                    gsInfo << "Computing Internal basis functions 2 ... \n";
                gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler_g12(multiPatch, mb_g1, bcInfo, bcInfo2, source);
                g1BiharmonicAssembler_g12.assemble(false, 0);

                //gsMatrix<> tempo = g1BiharmonicAssembler_g12.matrix().toDense().block(mb_g1.basis(0).size(),0,mb_g1.basis(1).size(),mb_g1.basis(0).size());

                //gsMatrix<> tempo2 = g1BiharmonicAssembler_g12.matrix().toDense().block(0,0,mb_g1.basis(1).size(),mb_g1.basis(0).size());
                //gsInfo << g1BiharmonicAssembler_g12.matrix().row(mb_g1.basis(0).size()) << "\n";
                //gsInfo << mb_g1.basis(0).size() << "\n";
                //gsInfo << mb_g1.basis(1).size() << "\n";
                if (g1OptionList.getSwitch("info"))
                    gsInfo << "Computing Internal basis functions 3 ... \n";
                gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler_g21(multiPatch, mb_g1, bcInfo, bcInfo2, source);
                g1BiharmonicAssembler_g21.assemble(false, 1);

                //gsInfo << g1BiharmonicAssembler.matrix().block(0,0,10,mb_g1.basis(0).size()) << "\n";
                //gsInfo << g1BiharmonicAssembler_g21.matrix().bottomRows(2) << "\n";
                //gsInfo << "test: " <<  g1BiharmonicAssembler_g21.matrix().block(mb_g1.basis(0).size(),0,10,mb_g1.basis(0).size()) << "\n";
                //gsInfo << g1BiharmonicAssembler_g21.matrix().row(mb_g1.basis(0).size()) << "\n";
                if (g1OptionList.getSwitch("info"))
                    gsInfo << "Construct system ... \n";
                g1BiharmonicAssembler.constructSystem(g1BiharmonicAssembler_g22.matrix(),
                                                      g1BiharmonicAssembler_g22.rhs(),
                                                      g1BiharmonicAssembler_g12.matrix(),
                                                      g1BiharmonicAssembler_g21.matrix(),
                                                      mb_g1);

                //gsInfo << " test " << g1BiharmonicAssembler.matrix().bottomRows(2) << "\n";
                //gsInfo << " test 3" << g1BiharmonicAssembler.matrix().rightCols(2) << "\n";
            }
            //gsInfo << clock.stop() << "\n";

            //gsInfo << "bdy: " << g1BiharmonicAssembler.get_bValue().transpose() << "\n";
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Finalize ... ";
            g1System.finalize(multiPatch, mb[0], g1BiharmonicAssembler.get_bValue());
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Solving system... \n";
            gsMatrix<> solVector = g1System.solve(g1BiharmonicAssembler.matrix(), g1BiharmonicAssembler.rhs());
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Solving finished! \n";

            dofs_size[refinement_level] = solVector.rows();


            if (g1OptionList.getSwitch("plot"))
            {
                // construct solution: INTERIOR
                gsMultiPatch<> mpsol;
                g1BiharmonicAssembler.constructSolution(solVector
                                                            .block(g1System.get_numBoundaryVertexFunctions().last(),
                                                                   0,
                                                                   mb[0].size(),
                                                                   1), mpsol);
                gsField<> solField(multiPatch, mpsol);

                // construct solution for plotting
                std::vector<gsMultiPatch<>> g1Basis;
                g1System.constructG1Solution(solVector, g1Basis, multiPatch, mb);

                g1BiharmonicAssembler.plotParaview(solField, g1Basis);
            }

            // construct solution: G1 Basis
            if (g1OptionList.getSwitch("info"))
                gsInfo << "Construct Solution ... \n";
            gsSparseMatrix<real_t> Sol_sparse;
            g1System.constructSparseG1Solution(solVector, Sol_sparse);


            sol_proof.push_back(Sol_sparse);
            mb_proof.push_back(mb);
            g1System_proof.push_back(g1System);
            multiPatch_proof.push_back(multiPatch);
        } // Proof


#ifdef _OPENMP
        omp_set_num_threads(g1OptionList.getInt("threads"));
        //omp_set_num_threads(1);
        omp_set_nested(1);
#endif
        if (g1OptionList.getSwitch("info"))
            gsInfo << "Compute Error ... \n";
#pragma omp parallel for
        for (index_t e = 0; e < 6; ++e)
        {
            if (e == 0)
            {
                gsNormL2<real_t> errorL2(multiPatch_proof.back(), mb_proof.back(), sol_proof.back(), solVal);
                errorL2.compute(g1System_proof.back(), g1OptionList.getSwitch("isogeometric"));
                l2Error_vec[refinement_level] = errorL2.value();
            }
            else if (e == 1)
            {
                gsSeminormH1<real_t> errorSemiH1(multiPatch_proof.back(), mb_proof.back(), sol_proof.back(), solVal);
                errorSemiH1.compute(g1System_proof.back(), g1OptionList.getSwitch("isogeometric"));
                h1SemiError_vec[refinement_level] = errorSemiH1.value();
            }
            else if (e == 2)
            {
                gsSeminormH2<real_t> errorSemiH2(multiPatch_proof.back(), mb_proof.back(), sol_proof.back(), solVal);
                errorSemiH2.compute(g1System_proof.back(), g1OptionList.getSwitch("isogeometric"));
                h2SemiError_vec[refinement_level] = errorSemiH2.value();
            }
            else if (e == 3)
            {
                gsH1NormWithJump<real_t> errorJump(multiPatch_proof.back(), mb_proof.back(), sol_proof.back());
                errorJump.compute(g1System_proof.back(), g1OptionList.getSwitch("isogeometric"), "all");
                h1SemiError_jump_edge.row(refinement_level) = errorJump.value().transpose();
            }

        }

        if (g1OptionList.getSwitch("h1projectionProof"))
        {
            gsH1NormProof<real_t> errorProof(multiPatch_proof, mb_proof, sol_proof);
            errorProof.compute(g1System_proof, g1OptionList.getSwitch("isogeometric"), "all");
            h1SemiError_jump_proof.row(refinement_level) = errorProof.value().transpose();
        }

    } // refinement_level

    for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
    {
        h2SemiError_vec[i] = math::sqrt(h2SemiError_vec[i]*h2SemiError_vec(i) +
            h1SemiError_vec[i]*h1SemiError_vec[i] + l2Error_vec[i]*l2Error_vec[i]);
        h1SemiError_vec[i] = math::sqrt(h1SemiError_vec[i]*h1SemiError_vec[i] +
            l2Error_vec[i]*l2Error_vec[i]);
    }

    const char* var1 = "|%-5d|%-14.6e|%-5.2f|%-14.6e|%-5.2f|%-14.6e|%-5.2f\n";
    const char* var2 = "%-5d & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\ \n";
    const char* var3 = "|%-14.6e|%-5.2f";
    const char* var4 = "%-5f %-8d %-14.6e %-5.2f %-14.6e %-5.2f %-14.6e %-5.2f %-14.6e %-5.2f \n";

    if (std::string (typeid(real_t).name())== "e") // long double
    {
        var1 = "|%-5d|%-14.6Le|%-5.2Lf|%-14.6Le|%-5.2Lf|%-14.6Le|%-5.2Lf\n";
        var2 = "%-5d & %-14.6Le & %-5.2Lf & %-14.6Le & %-5.2Lf & %-14.6Le & %-5.2Lf \\\\ \n";
        var3 = "|%-14.6Le|%-5.2Lf";
        var4 = "%-5Lf %-8d %-14.6Le %-5.2Lf %-14.6Le %-5.2Lf %-14.6Le %-5.2Lf %-14.6Le %-5.2Lf \n";
        gsInfo << "Double \n";
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
        printf(var1, num_knots[0], l2Error_vec[0],
               rate(0,0),h1SemiError_vec[0], rate(0,1),h2SemiError_vec[0], rate(0,2));
        for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
        {
            rate(i,0) = log2(l2Error_vec[i-1] / l2Error_vec[i]);
            rate(i,1) = log2(h1SemiError_vec[i-1] / h1SemiError_vec[i]);
            rate(i,2) = log2(h2SemiError_vec[i-1] / h2SemiError_vec[i]);

            printf(var1, num_knots[i], l2Error_vec[i],
                   rate(i,0),h1SemiError_vec[i], rate(i,1),h2SemiError_vec[i], rate(i,2));
        }
        if (g1OptionList.getSwitch("latex"))
        {
            printf(var2, num_knots[0],
                   l2Error_vec[0], rate(0,0),h1SemiError_vec[0], rate(0,1), h2SemiError_vec[0], rate(0,2));
            for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
            {
                printf(var2, num_knots[i],
                       l2Error_vec[i], rate(i,0),h1SemiError_vec[i], rate(i,1),h2SemiError_vec[i], rate(i,2));
            }
        }

        // JUMP
        gsMatrix<> rate_jump;
        rate_jump.setZero(g1OptionList.getInt("loop") + 1,multiPatch_init.interfaces().size());
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
            printf("|%-14s|%-5s", ("Single E. " + std::to_string(i)).c_str(), "Rate");
        gsInfo << "\n";
        printf("|%-5s","-----");
        for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
            printf("|%-14s|%-5s", "--------------", "-----");
        gsInfo << "\n";

        printf("|%-5d",num_knots[0]);
        for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
            printf(var3, h1SemiError_jump_edge(0,i), rate_jump(0,i));
        printf("\n");

        for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
        {
            printf("|%-5d",num_knots[i]);
            for (size_t j = 0; j < multiPatch_init.interfaces().size(); j++)
            {
                rate_jump(i,j) = log2(h1SemiError_jump_edge(i-1,j) / h1SemiError_jump_edge(i,j));
                printf(var3, h1SemiError_jump_edge(i,j), rate_jump(i,j));
            }
            printf("\n");
        }

        if (g1OptionList.getSwitch("latex_plot"))
        {
            gsInfo << "======";
            for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
                gsInfo << "===============================================================";

            gsInfo << "\n";

            gsInfo << "======================= Latex Plot ===================================\n\n";

            printf("%-8s %-8s %-14s %-5s %-14s %-5s %-14s %-5s %-14s %-5s\n", "k","Dofs","L2error", "Rate", "H1error", "Rate", "H2error", "Rate", "Jump", "Rate");
            printf(var4, mesh_size[0], dofs_size[0], l2Error_vec[0], rate(0,0),
                   h1SemiError_vec[0], rate(0,1), h2SemiError_vec[0], rate(0,2), h1SemiError_jump_edge(0,0), rate_jump(0,0));

            for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
            {
                printf(var4, mesh_size[i], dofs_size[i], l2Error_vec[i], rate(i,0),
                       h1SemiError_vec[i], rate(i,1), h2SemiError_vec[i], rate(i,2), h1SemiError_jump_edge(i,0), rate_jump(i,0));
            }

        }


        gsInfo << "=====================================================================\n\n";


        if (g1OptionList.getSwitch("h1projectionProof"))
        {
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
                printf("|%-14s|%-5s", ("Single E. " + std::to_string(i)).c_str(), "Rate");
            gsInfo << "\n";
            printf("|%-5s","-----");
            for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
                printf("|%-14s|%-5s", "--------------", "-----");
            gsInfo << "\n";

            printf("|%-5d",num_knots[0]);
            for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
                printf(var3, h1SemiError_jump_proof(0,i), rate(0,i));
            printf("\n");

            for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
            {
                printf("|%-5d",num_knots[i]);
                for (size_t j = 0; j < multiPatch_init.interfaces().size(); j++)
                {
                    rate(i,j) = log2(h1SemiError_jump_proof(i-1,j) / h1SemiError_jump_proof(i,j));
                    printf(var3, h1SemiError_jump_proof(i,j), rate(i,j));
                }
                printf("\n");
            }

            gsInfo << "======";
            for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
                gsInfo << "===============================================================";

            gsInfo << "\n";
        }
/*
        if (g1OptionList.getSwitch("latex"))
        {
            for (size_t i = 0; i < multiPatch_init.interfaces().size(); i++)
                printf("%-5d & %-14.6e & %-5.2f \\\\", num_knots[0],
                       h1SemiError_jump_edge(0,i), rate(0,i));
            printf("\n");
            for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
            {
                printf("%-5d & ",num_knots[i]);
                for (size_t j = 0; j < multiPatch_init.interfaces().size(); j++)
                {
                    printf("%-14.6e & %-5.2f  \\\\",
                           h1SemiError_jump_edge(i,j),
                           rate(i, j));
                }
                printf("\n");
            }
        }
*/
    }
    else
    {
        gsInfo << "=====================================================================\n";
        gsInfo << "L2 Error: " << l2Error_vec[0] << "\n";
        gsInfo << "H1 Semi-error: " << h1SemiError_vec[0] << "\n";
        gsInfo << "H2 Semi-error: " << h2SemiError_vec[0] << "\n";
        gsInfo << "Jump error Edge: " << h1SemiError_jump_edge.row(0) << "\n";
        gsInfo << "=====================================================================\n";

    }
} // main
