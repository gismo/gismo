/** @file g1BiharmonicSurface_example.cpp

    @brief A Biharmonic Surface example

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat
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

# include <gsG1Basis/Norm/gsG1ASResidualNormL2.h>
# include <gsG1Basis/Norm/gsG1ASResidualSeminormH1.h>
# include <gsG1Basis/Norm/gsG1ASResidualSeminormH2.h>

using namespace gismo;

int main(int argc, char *argv[])
{

    gsG1OptionList g1OptionList;
    g1OptionList.initialize(argc, argv);

    g1OptionList.addInt("user", "User defined gluingData", user::name::andrea);
    g1OptionList.setSwitch("twoPatch",false);



    // ======= Solution =========
//    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",3);
    gsFunctionExpr<> source  ("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",3); // L2 approximation RHS

    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",3);
    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",3);
    gsFunctionExpr<>sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
                             "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",
                             "0",3);
    gsFunctionExpr<>sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
                             "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
                             "0",
                             "16*pi^2*sin(4*pi*x)*sin(4*pi*y)",
                             "0",
                             "0", 3);



//    gsFunctionExpr<> source  ("(1 / ( (1 + x^2 + \n"
//                                              "  y^2)^5) ) * 4 * (8 * x^10 - 12 * x^9 * (1 + y^2)^2 - \n"
//                                              "   6 * x * (1 + y^2)^2 * (2 - 12 * y + 3 * y^2 + 14 * y^3 + y^4) + \n"
//                                              "   x^8 * (35 - 6 * y + 33 * y^2 - 38 * y^3 - 30 * y^5) + (1 + y^2)^2 * (2 - \n"
//                                              "      12 * y + 13 * y^2 - 10 * y^3 + 19 * y^4 - 12 * y^5 + 8 * y^6) + \n"
//                                              "   2 * x^7 * (-17 - 42 * y + 17 * y^2 - 60 * y^3 + 73 * y^4 + 35 * y^6) + \n"
//                                              "   x^6 * (59 - 30 * y + 104 * y^2 - 150 * y^3 + 39 * y^4 - 110 * y^5 + 70 * y^7) - \n"
//                                              "   2 * x^5 * (22 + 48 * y + 21 * y^2 - 120 * y^3 + 39 * y^4 - 216 * y^5 + 55 * y^6 + \n"
//                                              "      15 * y^8) - \n"
//                                              "   2 * x^3 * (17 - 30 * y + 71 * y^2 - 210 * y^3 + 110 * y^4 - 120 * y^5 + 75 * y^6 + \n"
//                                              "      60 * y^7 + 19 * y^8) + \n"
//                                              "   x^2 * (17 - 42 * y + 72 * y^2 - 142 * y^3 + 126 * y^4 - 42 * y^5 + 104 * y^6 + \n"
//                                              "      34 * y^7 + 33 * y^8 - 24 * y^9) + \n"
//                                              "   x^4 * (47 - 54 * y + 126 * y^2 - 220 * y^3 + 118 * y^4 - 78 * y^5 + 39 * y^6 + \n"
//                                              "      146 * y^7 - 12 * y^9))",3);
//
//    gsFunctionExpr<> laplace ("-((2 * (- ( -1 + y)^2 * y^2 + 6 * x * (-1 + y)^2 * y^2 + \n"
//                                              "    x^3 * (2 - 12 * y + 8 * y^2 + 12 * y^3 - 6 * y^4) + \n"
//                                              "    x^2 * (-1 + 6 * y - 10 * y^2 + 8 * y^3 - 5 * y^4) + \n"
//                                              "    x^4 * (-1 + 6 * y - 5 * y^2 - 6 * y^3 + 4 * y^4) ) ) / (1 + x^2 + y^2))",3);
//
//    gsFunctionExpr<> solVal("-( (2 * (- ( -1 + y)^2 * y^2 + 6 * x * (-1 + y)^2 * y^2 + \n"
//                                            "    x^3 * (2 - 12 * y + 8 * y^2 + 12 * y^3 - 6 * y^4) + \n"
//                                            "    x^2 * (-1 + 6 * y - 10 * y^2 + 8 * y^3 - 5 * y^4) + \n"
//                                            "    x^4 * (-1 + 6 * y - 5 * y^2 - 6 * y^3 + 4 * y^4) ) ) / (1 + x^2 + y^2))",3);
//
//    gsFunctionExpr<>sol1der ("(2 * (-1 + x) * x * (-1 + 2 * x) * (1 + x^2) * (-1 + y)^2 * y^2) / (1 + x^2 + y^2) - (\n"
//                                            " 2 * (-1 + x)^2 * x^3 * (-1 + y) * y^2 * (-1 + 2 * y) ) / (1 + x^2 + y^2)",
//                             "-( (2 * (-1 + x) * x^2 * (-1 + 2 * x) * (-1 + y)^2 * y^3) / (1 + x^2 + y^2)) + (\n"
//                                            " 2 * (-1 + x)^2 * x^2 * (-1 + y) * y * (-1 + 2 * y) * (1 + y^2) ) / (1 + x^2 + y^2)",
//                             "2 * (-1 + x) * x * (-1 + 2 * x) * (-1 + \n"
//                                             "    y)^2 * y^2 * (-( ( x^2 * y) / (1 + x^2 + y^2)) + ( (1 + x^2) * y) / (\n"
//                                             "    1 + x^2 + y^2) ) + \n"
//                                             " 2 * (-1 + x)^2 * x^2 * (-1 + y) * y * (-1 + \n"
//                                             "    2 * y) * (-( ( x * y^2) / (1 + x^2 + y^2)) + (x * (1 + y^2) ) / (1 + x^2 + y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( (2 * v^2 (4 * u^2 * (-1 + v) - (-1 + v)^2 + 6 * u * (-1 + v)^2 - \n"
//                                             "    12 * u^3 * (-1 + v) * v - 2 * u^5 * (2 - 6 * v + 3 * v^2) + \n"
//                                             "    u^6 * (1 - 6 * v + 4 * v^2) + u^4 * (-2 - 4 * v + 5 * v^2) ) ) / (1 + u^2 + v^2)^2\n"
//                                             " )",
//                             "-( ( 2 * u^2 * (-1 + 6 * v - 4 * v^2 - 2 * v^4 - 4 * v^5 + v^6 - \n"
//                                             "    2 * u * (-1 + 6 * v - 2 * v^2 - 6 * v^3 + 2 * v^4 - 6 * v^5 + 3 * v^6) + \n"
//                                             "    u^2 * (-1 + 6 * v - 12 * v^3 + 5 * v^4 - 6 * v^5 + 4 * v^6) ) ) / (1 + u^2 + \n"
//                                             "   v^2)^2)",
//                             "(2 * ( (-1 + v)^2 * v^4 - 6 * u * (-1 + v)^2 * v^4 - \n"
//                                             "   12 * u^3 * v^2 * (1 - 3 * v + 2 * v^2) - 2 * u^5 * (1 - 6 * v + 6 * v^2) + \n"
//                                             "   u^6 * (1 - 6 * v + 6 * v^2) + \n"
//                                             "   2 * u^2 * v^2 * (2 - 6 * v + 7 * v^2 - 6 * v^3 + 3 * v^4) + \n"
//                                             "   u^4 * (1 - 6 * v + 14 * v^2 - 24 * v^3 + 16 * v^4) ) ) / (1 + u^2 + v^2)^2",
//                             "(2 * u * v * (2 - 6 * v + 5 * v^2 - 4 * v^3 + 3 * v^4 - \n"
//                                             "   6 * u * (1 - 3 * v + 2 * v^2 - v^3 + v^4) - \n"
//                                             "   2 * u^3 * (2 - 3 * v + 2 * v^2 - 6 * v^3 + 3 * v^4) + \n"
//                                             "   u^2 * (5 - 12 * v + 6 * v^2 - 4 * v^3 + 3 * v^4) + \n"
//                                             "   u^4 * (3 - 6 * v + 3 * v^2 - 6 * v^3 + 4 * v^4) ) ) / (1 + u^2 + v^2)^2",
//                             "(2 * v * ( (-1 + v)^2 * v^2 - 6 * u * (-1 + v)^2 * v^2 + u^5 * (-4 + 6 * v) + \n"
//                                             "   6 * u^3 * (-1 + v)^2 * (-1 + v + v^2) + u^6 * (3 - 6 * v + 2 * v^2) + \n"
//                                             "   u^2 * (-1 + v)^2 * (2 - 2 * v + 3 * v^2) + \n"
//                                             "   u^4 * (5 - 12 * v + 8 * v^2 - 2 * v^4) ) ) / (1 + u^2 + v^2)^2",
//                             "-( (2 * u * (2 * u^3 * (1 - 6 * v + 4 * v^2 + 3 * v^3) + \n"
//                                             "    v^2 * (-2 + 6 * v - 5 * v^2 + 4 * v^3 - 3 * v^4) + \n"
//                                             "    6 * u * v^2 * (1 - 3 * v + 2 * v^2 - v^3 + v^4) + \n"
//                                             "    u^4 * (-1 + 6 * v - 3 * v^2 - 6 * v^3 + 2 * v^4) - \n"
//                                             "    u^2 * (1 - 6 * v + 9 * v^2 - 12 * v^3 + 8 * v^4 + 2 * v^6) ) ) / (1 + u^2 + v^2)^2\n"
//                                             " )", 3);


// ==============================================================================================================================================
/*
 *  QUADRATIC SURFACE (u, v, u * v)
 */


//    gsFunctionExpr<> source  ("(8 * (3 * x^10 + 3 * (1 + y^2)^3 + 12 * x^2 * (1 + y^2)^2 * (1 + 5 y^2) + \n"
//                              "   x^8 * (14 + 15 * y^2) + x^6 * (25 + 36 * y^2 + 5 * y^4) + \n"
//                              "   x^4 * (23 + 96 * y^2 + 178 * y^4 + 105 * y^6))) / (1 + x^2 + y^2)^5 ", 3);
//    gsFunctionExpr<> laplace ("(12 * (x^2 + x^4) ) / (1 + x^2 + y^2)",3);
//    gsFunctionExpr<> solVal("x^4",3);
//    gsFunctionExpr<>sol1der ("(4 * x^3 * (1 + x^2)) / (1 + x^2 + y^2)",
//                             "-((4 * x^4 * y) / (1 + x^2 + y^2))",
//                             "4 * x^3 * ( -( (x^2 * y) / (1 + x^2 + y^2) ) + ( (1 + x^2) * y ) / (1 + x^2 + y^2) )", 3);
//    gsFunctionExpr<>sol2der ("(12 * (x + x^3)^2 ) / (1 + x^2 + y^2)^2",
//                             "(12 * x^4 * y^2) / (1 + x^2 + y^2)^2",
//                             "(12 * x^2 * y^2) / (1 + x^2 + y^2)^2",
//                             "-((12 * x^3 * (1 + x^2) y) / (1 + x^2 + y^2)^2)",
//                             "(12 * x^2 * (1 + x^2) * y) / (1 + x^2 + y^2)^2",
//                             "-((12 * x^3 * y^2) / (1 + x^2 + y^2)^2)", 3);


//    gsFunctionExpr<> source  ("(4 * x * (1 - 15 * y^2 - 10 * y^4 + 6 * y^6 + x^4 * (1 + 15 * y^2) + \n"
//                              "   x^2 * (2 - 35 * y^4))) / (1 + x^2 + y^2)^5 ",3);
//    gsFunctionExpr<> laplace ("0",3);
//    gsFunctionExpr<> solVal("x",3);
//    gsFunctionExpr<>sol1der ("(1 + x^2) / (1 + x^2 + y^2)",
//                             "-( (x * y) / (1 + x^2 + y^2) )",
//                             "-( (x^2 * y) / (1 + x^2 + y^2) ) + ( (1 + x^2) * y) / (1 + x^2 + y^2)",3);
//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);


// PLANAR SOLUTION

//    gsFunctionExpr<> source  ("0 ",3);
//    gsFunctionExpr<> laplace ("0",3);
//    gsFunctionExpr<> solVal("x",3);
//    gsFunctionExpr<>sol1der ("1",
//                             "0",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);


//    gsFunctionExpr<> source  ("x^4 ",3);
//    gsFunctionExpr<> laplace ("12 * x^2",3);
//    gsFunctionExpr<> solVal("x^4",3);
//    gsFunctionExpr<>sol1der ("4 * x^3",
//                             "0",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("12 * x^2",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);


// ==============================================================================================================================================


//    gsFunctionExpr<> source  ("0",3);
//
//    gsFunctionExpr<> laplace ("0",3);
//
//    gsFunctionExpr<> solVal("(x - 1) * (y - 1)",3);
//
//    gsFunctionExpr<>sol1der ("(y - 1)",
//                             "(x - 1)",
//                             "0",3);
//
//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);

//    gsFunctionExpr<> source  ("0",3);
//
//    gsFunctionExpr<> laplace ("0",3);
//
//    gsFunctionExpr<> solVal("(x - 1) * (z - 1)",3);
//
//    gsFunctionExpr<>sol1der ("(z - 1)",
//                             "0",
//                             "(x - 1)",3);
//
//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);



//    gsFunctionExpr<> source  ("0",3);
//
//    gsFunctionExpr<> laplace ("0",3);
//
//    gsFunctionExpr<> solVal("1",3);
//
//    gsFunctionExpr<>sol1der ("0",
//                             "0",
//                             "0",3);
//
//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);


//    gsFunctionExpr<> source  ("(4 * x * y * ( 13 - 10 * x^4 + 3 * y^2 - 10 * y^4 + x^2 * (3 + 36 * y^2) ) ) / (1 + x^2 + y^2)^5",3);
//
//    gsFunctionExpr<> laplace ("-( (2 * x * y) / (1 + x^2 + y^2))",3);
//
//    gsFunctionExpr<> solVal("x * y",3);
//
//    gsFunctionExpr<>sol1der ("-( (x^2 * y) / (1 + x^2 + y^2) ) + ((1 + x^2) * y) / (1 + x^2 + y^2)",
//                             "-( (x * y^2) / (1 + x^2 + y^2) ) + (x * (1 + y^2))/(1 + x^2 + y^2)",
//                             "y * (-( ( x^2 * y) / (1 + x^2 + y^2) ) + ( (1 + x^2) * y) / (1 + x^2 + y^2)) + \n"
//                                             " x * (-( ( x * y^2) / (1 + x^2 + y^2)) + (x * (1 + y^2) ) / (1 + x^2 + y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);



    gsFunctionWithDerivatives<real_t> solution(solVal, sol1der, sol2der);

    // ======= Geometry =========
    std::string string_geo;
    index_t numDegree = 0;
    switch(g1OptionList.getInt("geometry"))
    {
        case 0:
            string_geo = "KirchhoffLoveGeo/parabola_surfaceRoundedBoundary.xml";
            numDegree = 1; // 1 == degree 3
            break;
        case 1:
            string_geo = "KirchhoffLoveGeo/parabola_surfaceSquaredBoundary.xml";
            numDegree = 1; // 1 == degree 3
            break;
        case 2:
            string_geo = "KirchhoffLoveGeo/flag_surface.xml";
            numDegree = 0; // 1 == degree 3
            break;
        case 3:
            string_geo = "KirchhoffLoveGeo/parabola_surfaceTwoPatchRoundBoundary.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 4:
            string_geo = "KirchhoffLoveGeo/square3dPositiveOrientation.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 5:
            string_geo = "KirchhoffLoveGeo/square3dNegativeOrientation.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 6:
            string_geo = "KirchhoffLoveGeo/squareSurface3d.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 7:
            string_geo = "KirchhoffLoveGeo/bentSquareSurface.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 8:
            string_geo = "KirchhoffLoveGeo/surface_fourPatch.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 9:
            string_geo = "KirchhoffLoveGeo/singlePatch_quadraticParamSurf.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 10:
            string_geo = "KirchhoffLoveGeo/singlePatch_firstCoordSquared.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 11:
            string_geo = "KirchhoffLoveGeo/square_TwoPatch3d.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 12:
            string_geo = "KirchhoffLoveGeo/parabola_surfaceCubicParamSquaredBdy.xml";
            numDegree = 0; // 1 == degree 3
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


    std::vector<gsSparseMatrix<>> sol_vector;
    std::vector<gsMultiBasis<>> sol_vec_basis;
    std::vector<gsG1System<real_t>> sol_vec_sys;

    for (index_t refinement_level = 0; refinement_level < g1OptionList.getInt("loop"); refinement_level++)
    {
        gsMultiPatch<> multiPatch(multiPatch_init);
        multiPatch.uniformRefine_withSameRegularity(num_knots[refinement_level], g1OptionList.getInt("regularity"));
        gsInfo << "###### Level: " << refinement_level << " with " << num_knots[refinement_level] << " inner knots ###### " << "\n";

        gsMultiBasis<> mb(multiPatch);

#ifdef _OPENMP
        omp_set_num_threads( g1OptionList.getInt("threads"));
        omp_set_nested(1);
#endif

        gsG1System<real_t> g1System(multiPatch, mb, g1OptionList.getSwitch("neumann"), g1OptionList.getSwitch("twoPatch"));

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
            gsInfo << "============================================================ \n";
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
//            gsInfo << "============================================================ \n";
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

            if(g1OptionList.getSwitch("twoPatch") == true)
            {
                if (patchIndex.size() == 1)
                {
                    gsG1AuxiliaryVertexMultiplePatches singleVertex(multiPatch, patchIndex, vertIndex);
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
                        g1System.insertVertex(singleBasisFunction, patchIndex, numVer, singleVertex.get_internalDofs(), i);
                    }
                }
            }
            else
            {
                gsG1AuxiliaryVertexMultiplePatches singleVertex(multiPatch, patchIndex, vertIndex);

                singleVertex.computeG1InternalVertexBasis(g1OptionList);
                //                for (index_t i = 0; i < (g1OptionList.get("twoPatch") ? 4 : 6); i++)
                for (index_t i = 0; i < 6; i++)
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
                    g1System.insertVertex(singleBasisFunction, patchIndex, numVer, singleVertex.get_internalDofs(), i);}
            }
//            gsInfo << "============================================================ \n";
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
        gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(multiPatch, mb, bcInfo, bcInfo2, source);
        g1BiharmonicAssembler.assemble();

        gsInfo << "Computing Boundary data ... \n";

        if (!g1OptionList.getSwitch("neumann"))
            g1BiharmonicAssembler.computeDirichletDofsL2Proj(g1System); // Compute boundary values with laplace
        else
            g1BiharmonicAssembler.computeDirichletAndNeumannDofsL2Proj(g1System); // Compute boundary values with neumann

        g1System.finalize(multiPatch,mb,g1BiharmonicAssembler.get_bValue());

//        gsInfo << "Dir val: " << g1BiharmonicAssembler.get_bValue() << "\n";

        gsInfo << "Solving the system ... \n";
        gsMatrix<> solVector = g1System.solve(g1BiharmonicAssembler.matrix(), g1BiharmonicAssembler.rhs());

//        gsInfo << "Sol: " << solVector << "\n";

        gsInfo << "Plotting the solution ... \n";
        if (g1OptionList.getSwitch("plot"))
        {
            gsField<> exactField(multiPatch,solVal);
            gsWriteParaview(exactField,"G1Biharmonic_exact",2000);

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
        gsInfo << "Construct G1 solution ... \n";
        gsSparseMatrix<real_t> Sol_sparse;
        g1System.constructSparseG1Solution(solVector,Sol_sparse);

        if(sol_vector.size() == 0 || sol_vector.size() == 1)
        {
            sol_vector.push_back(Sol_sparse);
            sol_vec_basis.push_back(mb);
            sol_vec_sys.push_back(g1System);
        }
        else
        {
            sol_vector.at(0) = sol_vector.at(1);
            sol_vector.at(1) = Sol_sparse;

            sol_vec_basis.at(0) = sol_vec_basis.at(1);
            sol_vec_basis.at(1) = mb;

            sol_vec_sys.at(0) = sol_vec_sys.at(1);
            sol_vec_sys.at(1) = g1System;
        }

#ifdef _OPENMP
        omp_set_num_threads(g1OptionList.getInt("threads"));
        omp_set_nested(1);
#endif
        gsInfo << "Computing the error ... \n";
#pragma omp parallel for
        for (index_t e = 0; e < 4; ++e)
        {
            if (e == 0)
            {
                if(g1OptionList.getSwitch("residual") == true && sol_vector.size() == 2)
                {
                    gsG1ASResidualNormL2<real_t> errorL2(multiPatch, sol_vector, sol_vec_basis);
                    errorL2.compute(sol_vec_sys);
                    l2Error_vec[refinement_level] = errorL2.value();
                }
                else
                {
                    gsNormL2<real_t> errorL2(multiPatch, Sol_sparse, solVal);
                    errorL2.compute(g1System.get_numBasisFunctions());
                    l2Error_vec[refinement_level] = errorL2.value();
                }
            }

            else if (e == 1)
            {
                if(g1OptionList.getSwitch("residual") == true && sol_vector.size() == 2)
                {
                    gsG1ASResidualSeminormH1<real_t> errorSemiH1(multiPatch, sol_vector, sol_vec_basis);
                    errorSemiH1.compute(sol_vec_sys);
                    h1SemiError_vec[refinement_level] = errorSemiH1.value();
                }
                else
                {
                    gsSeminormH1<real_t> errorSemiH1(multiPatch, Sol_sparse, sol1der);
                    errorSemiH1.compute(g1System.get_numBasisFunctions());
                    h1SemiError_vec[refinement_level] = errorSemiH1.value();
                }
            }
            else if (e == 2)
            {
                if(g1OptionList.getSwitch("residual") == true && sol_vector.size() == 2)
                {
                    gsG1ASResidualSeminormH2<real_t> errorSemiH2(multiPatch, sol_vector, sol_vec_basis);
                    errorSemiH2.compute(sol_vec_sys);
                    h2SemiError_vec[refinement_level] = errorSemiH2.value();
                }
                else
                {
                    gsSeminormH2<real_t> errorSemiH2(multiPatch, Sol_sparse, sol2der);
                    errorSemiH2.compute(g1System.get_numBasisFunctions());
                    h2SemiError_vec[refinement_level] = errorSemiH2.value();
                }
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

    gsInfo << "=====================================================================\n";
    if (g1OptionList.getInt("loop") > 1)
    {
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
        gsInfo << "H1 Error: " << h1SemiError_vec[0] << "\n";
        gsInfo << "H2 Error: " << h2SemiError_vec[0] << "\n";
    }
    gsInfo << "=====================================================================\n";


} // main