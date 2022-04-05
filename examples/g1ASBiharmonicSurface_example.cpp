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



    // ======= 2D Solution =========

//    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
////    gsFunctionExpr<> source  ("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
//
//    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",2);
//    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",2);
//    gsFunctionExpr<>sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
//                             "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",2);
//    gsFunctionExpr<>sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
//                             "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
//                             " 16*pi^2*sin(4*pi*x)*sin(4*pi*y)", 2);

//    gsFunctionExpr<> source  ("(cos(2*pi*x) - pi/2) * sin(2*pi*y)",2);
//    gsFunctionExpr<> laplace ("-4*pi^2*sin(2*pi*y)*cos(2*pi*x) + 2 * pi^2 * (pi - 2 * cos(2*pi*x)) * sin(2*pi*y)",2);
//    gsFunctionExpr<> solVal("(cos(2*pi*x) - pi/2) * sin(2*pi*y)",2);
//    gsFunctionExpr<>sol1der ("-2*pi*sin(2*pi*y)*sin(2*pi*x)",
//                             "-pi * (pi - 2 * cos(2*pi*x)) * cos(2*pi*y)",2);
//    gsFunctionExpr<>sol2der ("-4*pi^2*sin(2*pi*y)*cos(2*pi*x)",
//                             "2 * pi^2 * (pi - 2 * cos(2*pi*x)) * sin(2*pi*y)",
//                             "-4*pi^2*sin(2*pi*x)*cos(2*pi*y)", 2);

//    gsFunctionExpr<> source  ("pi*pi*pi*pi*(4*cos(pi*x/2)*cos(pi*y/2) - cos(pi*x/2) - cos(pi*y/2))/16",2);
////    gsFunctionExpr<> source  ("(cos(pi*x/2) - 1) * (cos(pi*y/2) - 1)",2);
//
//    gsFunctionExpr<> laplace ("-pi*pi*(2*cos(pi*x/2)*cos(pi*y/2) - cos(pi*x/2) - cos(pi*y/2))/4",2);
//    gsFunctionExpr<> solVal("(cos(pi*x/2) - 1) * (cos(pi*y/2) - 1)",2);
//    gsFunctionExpr<>sol1der ("-pi*(cos(pi*y/2) - 1)*sin(pi*x/2)/2",
//                             "-pi*(cos(pi*x/2) - 1)*sin(pi*y/2)/2",2);
//    gsFunctionExpr<>sol2der ("-pi*pi*(cos(pi*y/2) - 1)*cos(pi*x/2)/4",
//                             "-pi*pi*(cos(pi*x/2) - 1)*cos(pi*y/2)/4",
//                             "pi*pi*sin(pi*x/2)*sin(pi*y/2)/4", 2);


//    gsFunctionExpr<> source  ("1/256 * pi^4 (-cos(4 - (pi * y)/4) + cos(4 - (pi * x)/4) * (-1 + 4 * cos(4 - (pi * y)/4)))",2);
////    gsFunctionExpr<> source  ("(cos(pi*x/4 -4) - 1) * (cos(pi*y/4 -4) - 1)",2);
//
//    gsFunctionExpr<> laplace ("1/16 * pi^2 *(cos(4 - (pi * x)/4) * (1 - 2 * cos(4 - (pi * y)/4)) + cos(4 - (pi * y)/4))",2);
//    gsFunctionExpr<> solVal("(cos(pi*x/4 -4) - 1) * (cos(pi*y/4 -4) - 1)",2);
//    gsFunctionExpr<>sol1der ("1/4 *pi* (-1 + cos(4 - (pi * y)/4)) * sin(4 - (pi * x)/4)",
//                             "1/4 *pi* (-1 + cos(4 - (pi * x)/4)) * sin(4 - (pi * y)/4)",2);
//    gsFunctionExpr<>sol2der ("-1/16 * pi^2 * cos(4 - (pi * x)/4) * (-1 + cos(4 - (pi * y)/4))",
//                             "-1/16 * pi^2 * (-1 + cos(4 - (pi * x)/4)) * cos(4 - (pi * y)/4)",
//                             "1/16 * pi^2 * sin(4 - (pi * x)/4) * sin(4 - (pi * y)/4)", 2);

//    gsFunctionExpr<> source  ("1/256 * pi^4 (-cos( - (pi * y)/4) + cos( - (pi * x)/4) * (-1 + 4 * cos( - (pi * y)/4)))",2);
////    gsFunctionExpr<> source  ("(cos(pi*x/4) - 1) * (cos(pi*y/4 ) - 1)",2);
//
//    gsFunctionExpr<> laplace ("1/16 * pi^2 *(cos( - (pi * x)/4) * (1 - 2 * cos( - (pi * y)/4)) + cos( - (pi * y)/4))",2);
//    gsFunctionExpr<> solVal("(cos(pi*x/4 ) - 1) * (cos(pi*y/4 ) - 1)",2);
//    gsFunctionExpr<>sol1der ("1/4 *pi* (-1 + cos( - (pi * y)/4)) * sin( - (pi * x)/4)",
//                             "1/4 *pi* (-1 + cos( - (pi * x)/4)) * sin( - (pi * y)/4)",2);
//    gsFunctionExpr<>sol2der ("-1/16 * pi^2 * cos( - (pi * x)/4) * (-1 + cos( - (pi * y)/4))",
//                             "-1/16 * pi^2 * (-1 + cos( - (pi * x)/4)) * cos( - (pi * y)/4)",
//                             "1/16 * pi^2 * sin( - (pi * x)/4) * sin( - (pi * y)/4)", 2);


//    gsFunctionExpr<> source("5 * cos(pi * x / 27 ) * cos(pi * y / 27 )",2);
//
//    gsFunctionExpr<> laplace ("-5/729 * pi^2 * cos((pi * x)/27) * cos((pi * y)/27) -5/729 * pi^2 * cos((pi * x)/27) * cos((pi * y)/27)",2);
//    gsFunctionExpr<> solVal("5 * cos(pi * x / 27 ) * cos(pi * y/27 )",2);
//    gsFunctionExpr<>sol1der (" -5/27 * pi * sin((pi * x)/27) * cos((pi * y)/27)",
//                             " -5/27 * pi * cos((pi * x)/27) * sin((pi * y)/27)",2);
//    gsFunctionExpr<>sol2der (" -5/729 * pi^2 * cos((pi * x)/27) * cos((pi * y)/27)",
//                             " -5/729 * pi^2 * cos((pi * x)/27) * cos((pi * y)/27)",
//                             " 5/729 * pi^2 * sin((pi * x)/27) * sin((pi * y)/27)", 2);

//    gsFunctionExpr<> source  ("1/256 * pi^4 (-cos(4 - (pi * y)/4) + cos(4 - (pi * x)/4) * (-1 + 4 * cos(4 - (pi * y)/4)))",2);
////    gsFunctionExpr<> source  ("(cos(pi*x/4 -3.5) - 1) * (cos(pi*y/4 -3.5) - 1)",2);
//
//    gsFunctionExpr<> laplace ("-1/16  * pi^2 ( -cos(3.5 - (pi * y) / 4) + cos(3.5 - pi * x / 4) * (-1 + 2 * cos(3.5 - (pi * y) / 4)))",2);
//    gsFunctionExpr<> solVal("(cos( pi*x/4 -3.5 ) - 1) * (cos( pi*y/4 -3.5 ) - 1)",2);
//    gsFunctionExpr<>sol1der ("1/4 *pi* (-1 + cos( pi * y /4 -3.5 )) * sin(3.5 - pi * x/4)",
//                             "1/4 *pi* (-1 + cos(3.5 - pi * x/4)) * sin(3.5 - pi * y/4)",2);
//    gsFunctionExpr<>sol2der ("-1/16 * pi^2 * cos(3.5 - (pi * x)/4) * (-1 + cos((pi * y)/4 -3.5))",
//                             "-1/16 * pi^2 * (-1 + cos(3.5 - (pi * x)/4)) * cos( 3.5 - pi * y /4 )",
//                             "1/16 * pi^2 * sin(3.5 - pi * x /4) * sin( 3.5 - pi * y /4 )", 2);

//    gsFunctionExpr<> source  ("0.05",2);
//    gsFunctionExpr<> laplace ("0",2);
//    gsFunctionExpr<> solVal("0 ",2);
//    gsFunctionExpr<>sol1der ("0",
//                             "0",2);
//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0", 2);

//    gsFunctionExpr<> source  ("0",2);
//    gsFunctionExpr<> laplace ("0",2);
//    gsFunctionExpr<> solVal("x",2);
//    gsFunctionExpr<>sol1der ("1",
//                             "0",2);
//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0", 2);


//    gsFunctionExpr<> source  ("(6 - x^2 - y^2)^0.5",3);
//    gsFunctionExpr<> laplace ("(y^2 - 6) / (6 - x^2 - y^2)^1.5 + (x^2 - 6) / (6 - x^2 - y^2)^1.5",3);
//    gsFunctionExpr<> solVal("(6 - x^2 - y^2)^(1/2)",3);
//    gsFunctionExpr<>sol1der ("-x / (6 - x^2 - y^2)^0.5",
//                             "- y / (6 - x^2 - y^2)^0.5",
//                            "0",3);
//    gsFunctionExpr<>sol2der ("(y^2 - 6) / (6 - x^2 - y^2)^1.5",
//                             "(x^2 - 6) / (6 - x^2 - y^2)^1.5",
//                             "0",
//                             "- x * y / (6 - x^2 - y^2)^1.5",
//                            "0",
//                            "0",3);

//    gsFunctionExpr<> source  ("(6 - x^2 - y^2)^0.5",2);
//    gsFunctionExpr<> laplace ("(y^2 - 6) / (6 - x^2 - y^2)^1.5 + (x^2 - 6) / (6 - x^2 - y^2)^1.5",2);
//    gsFunctionExpr<> solVal("(6 - x^2 - y^2)^(1/2)",2);
//    gsFunctionExpr<>sol1der ("-x / (6 - x^2 - y^2)^0.5",
//                             "- y / (6 - x^2 - y^2)^0.5",2);
//    gsFunctionExpr<>sol2der ("(y^2 - 6) / (6 - x^2 - y^2)^1.5",
//                             "(x^2 - 6) / (6 - x^2 - y^2)^1.5",
//                             "- x * y / (6 - x^2 - y^2)^1.5", 2);

//    gsFunctionExpr<> source  ("16 * (1 - x) * x * (1 - y) * y",2);
//
//    gsFunctionExpr<> laplace ("- 2 * 16 * ( (1 - x) * x + (1 - y) * y )",2);
//
//    gsFunctionExpr<> solVal("16 * (1 - x) * x * (1 - y) * y",2);
//
//    gsFunctionExpr<>sol1der ("16 * ((1 - x) * (1 - y) * y - x * (1 - y) * y)",
//                             "16 * ((1 - y) * (1 - x) * x - y * (1 - x) * x)",2);
//
//    gsFunctionExpr<>sol2der ("- 2 * 16 * (1 - y) * y",
//                             "- 2 * 16 * (1 - x) * x",
//                             "16 * ((1 - x) * (1 - y) - (1 - x) * y + x * y - x * (1 - y))", 2);

//    gsFunctionExpr<> source  ("0.00005 * (x - 1)^2 * (y - 1)^2 * x^2 * y^2",2);
//    gsFunctionExpr<> laplace ("0.00005 * 2 * (6 * x^2 - 6 * x + 1) * (y - 1)^2 * y^2 + 2 * (6 * y^2 - 6 * y + 1) * (x - 1)^2 * x^2",2);
//    gsFunctionExpr<> solVal("0.00005 *(x - 1)^2 * (y - 1)^2 * x^2 * y^2",2);
//    gsFunctionExpr<>sol1der ("0.00005 *2 * (x - 1) * x *(2 * x - 1) * (y - 1)^2 * y^2",
//                             "0.00005 * 2 * (x - 1)^2 * x^2 * (y - 1) * y * (2 * y - 1)",2);
//    gsFunctionExpr<>sol2der ("0.00005 * 2 * (6 * x^2 - 6 * x + 1) * (y - 1)^2 * y^2",
//                             "0.00005 * 2 * (6 * y^2 - 6 * y + 1) * (x - 1)^2 * x^2",
//                             "0.00005 * 4 * (x - 1) * x * (2 * x - 1) * (y - 1) * y * (2 * y - 1)", 2);


//    gsFunctionExpr<> source  ("0.00000001 * 4/9 *(290409916 + 125 *x^6 - 70748808 *y - 33996207 *y^2 + 5199403 *y^3 + \n"
//                              "   448683 *y^4 - 46191 *y^5 + 636 *y^6 + 3 *x^5 *(-9226 + 431 *y) - \n"
//                              "   3 *x^4 *(-261466 + 32525 *y + 4810 *y^2) + \n"
//                              "   x^3 *(1789314 - 926955 *y + 114060 *y^2 + 22270 *y^3) - \n"
//                              "   3 *x^2 *(13649915 - 1431401 *y - 778947 *y^2 + 19250 *y^3 + 6525 *y^4) + \n"
//                              "   3 *x *(-6476056 + 7330526 *y - 243552 *y^2 - 344275 *y^3 - 6290 *y^4 + \n"
//                              "      2195 *y^5))",2);

//    gsFunctionExpr<> source  ("0.00000001 * (16/3 - x/3 - y) * (16 - 3*x - y) * (5 - x) * (x - 7 - y) * (x/3 - 5 - y) * (-5 - x/2 - y) * (-19 - 4*x - y) * (11 + 2*x - y) * (7 + x - y) * (5 - y)",2);
//
//    gsFunctionExpr<> laplace ("0",2);
//
//    gsFunctionExpr<> solVal("0.00000001 * (16/3 - x/3 - y) * (16 - 3*x - y) * (5 - x) * (x - 7 - y) * (x/3 - 5 - y) * (-5 - x/2 - y) * (-19 - 4*x - y) * (11 + 2*x - y) * (7 + x - y) * (5 - y)",2);
//
////    gsFunctionExpr<> solVal("0",3);
//
//    gsFunctionExpr<>sol1der ("0.00000001 *  1/18 *(y - 5) *(216 *x^8 + 16 *x^7 *(y - 253) - 7* x^6 *(293 *y^2 + 370 *y + 4089) + 6 *x^5 *(23 *y^3 + 440 *y^2 + 11251 *y + 109390) \n"
//                             "+ 5 *x^4 *(712 *y^4 + 1789 *y^3 - 696 *y^2 - 27219 *y + 242206) \n"
//                             "- 4 *x^3 *(376 *y^5 + 2136 *y^4 - 28431 *y^3 - 162166 *y^2 + 971107 *y + 7561994) - 3 *x^2 *(173 *y^6 - 1476 *y^5 + 19462 *y^4 + 137317 *y^3 - 479192 *y^2 - 2096665 *y + 5849981) \n"
//                             "+ 2 *x *(63 *y^7 - 228 *y^6 + 10450 *y^5 + 155112 *y^4 - 1106231 *y^3 - 10292682 *y^2 + 25363606 *y + 207550710) + 18 *y^8 - 375 *y^7 - 2665 *y^6 - 5646 *y^5 + 360108 *y^4 \n"
//                             "+ 1678929 *y^3 - 8300949 *y^2 - 26433420 *y + 47824000)",
//                             "0.00000001 *  1/18 *(x - 5) *(24 *x^8 + 4 *x^7 *(y - 99) + x^6 *(-879 *y^2 + 2210 *y - 4219) + 4 *x^5 *(23 *y^3 - 855 *y^2 + 7288 *y + 8010) \n"
//                             "+ x^4 *(3560 *y^4 - 6624 *y^3 - 46023 *y^2 + 98282 *y + 538501) \n"
//                             "- 2 *x^3 *(1128 *y^5 - 8260 *y^4 - 61662 *y^3 + 85041 *y^2 + 1536232 *y + 6977) + x^2 *(-1211 *y^6 + 2766 *y^5 - 51610 *y^4 + 456592 *y^3 + 2646921 *y^2 - 15960910 *y - 16403076) \n"
//                             "+ 2 *x *(252 *y^7 - 4928 *y^6 + 41685 *y^5 + 128130 *y^4 - 2622102 *y^3 - 524988 *y^2 + 36924741 *y - 641350) + 6 *(27 *y^8 - 200 *y^7 - 9135 *y^6 + 77154 *y^5 + 537165 *y^4 \n"
//                             "- 4451244 *y^3 - 9222777 *y^2 + 66565010 *y + 28929600))",2);
//
////    gsFunctionExpr<>sol1der ("0",
////                             "0",
////                             "0",3);
//
//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0", 2);


    // ======= 3D Solution =========

//    gsFunctionExpr<> source  ("cos(pi * x / 4)",3);
//    gsFunctionExpr<> laplace ("-pi^2 * cos(pi * x / 4) / 16",3);
//    gsFunctionExpr<> solVal("cos(pi * x / 4)",3);
//    gsFunctionExpr<>sol1der ("-pi * sin(pi * x / 4) / 4",
//                             "0",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("-pi^2 * cos(pi * x / 4) / 16",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);

////    gsFunctionExpr<> source  ("81 * pi^4 * cos(3 * pi * x / 4) / 256",3);
//    gsFunctionExpr<> source  ("cos(3 * pi * x / 4)",3);
//    gsFunctionExpr<> laplace ("-9 * pi^2 * cos(3 * pi * x / 4) / 16",3);
//    gsFunctionExpr<> solVal("cos(3 * pi * x / 4)",3);
//    gsFunctionExpr<>sol1der ("-3 * pi * sin(3 * pi * x / 4) / 4",
//                             "0",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("-9 * pi^2 * cos(3 * pi * x / 4) / 16",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);

//    gsFunctionExpr<> source  ("sin(pi * x / 4)",3);
//    gsFunctionExpr<> laplace ("-pi^2 * sin(pi * x / 4) / 16",3);
//    gsFunctionExpr<> solVal("sin(pi * x / 4)",3);
//    gsFunctionExpr<>sol1der ("pi * cos(pi * x / 4) / 4",
//                             "0",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("-pi^2 * sin(pi * x / 4) / 16",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);

//    gsFunctionExpr<> source  ("sin(3 * pi * x / 4)",3);
//    gsFunctionExpr<> laplace ("-9 * pi^2 * sin(3 * pi * x / 4) / 16",3);
//    gsFunctionExpr<> solVal("sin(3 * pi * x / 4)",3);
//    gsFunctionExpr<>sol1der ("3 * pi * cos(3 * pi * x / 4) / 4",
//                             "0",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("-9 * pi^2 * sin(3 * pi * x / 4) / 16",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);

//    gsFunctionExpr<> source  ("y",3);
//    gsFunctionExpr<> laplace ("0",3);
//    gsFunctionExpr<> solVal("y",3);
//    gsFunctionExpr<>sol1der ("0",
//                             "1",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);

//    gsFunctionExpr<> source  ("256*pi*pi*pi*pi*(4*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",3);
////    gsFunctionExpr<> source  ("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",3); // L2 approximation RHS
//
//    gsFunctionExpr<> laplace ("-16*pi*pi*(2*cos(4*pi*x)*cos(4*pi*y) - cos(4*pi*x) - cos(4*pi*y))",3);
//    gsFunctionExpr<> solVal("(cos(4*pi*x) - 1) * (cos(4*pi*y) - 1)",3);
//    gsFunctionExpr<>sol1der ("-4*pi*(cos(4*pi*y) - 1)*sin(4*pi*x)",
//                             "-4*pi*(cos(4*pi*x) - 1)*sin(4*pi*y)",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("-16*pi^2*(cos(4*pi*y) - 1)*cos(4*pi*x)",
//                             "-16*pi^2*(cos(4*pi*x) - 1)*cos(4*pi*y)",
//                             "0",
//                             "16*pi^2*sin(4*pi*x)*sin(4*pi*y)",
//                             "0",
//                             "0", 3);



//    gsFunctionExpr<> source  ("pi*pi*pi*pi*(4*cos(pi*x/2)*cos(pi*y/2) - cos(pi*x/2) - cos(pi*y/2))/16",3); // RHS
////    gsFunctionExpr<> source  ("(cos(pi*x/2) - 1) * (cos(pi*y/2) - 1)",3); // L2 approximation RHS
//
//    gsFunctionExpr<> laplace ("-pi*pi*(2*cos(pi*x/2)*cos(pi*y/2) - cos(pi*x/2) - cos(pi*y/2))/4",3);
//    gsFunctionExpr<> solVal("(cos(pi*x/2) - 1) * (cos(pi*y/2) - 1)",3);
//    gsFunctionExpr<>sol1der ("-pi*(cos(pi*y/2) - 1)*sin(pi*x/2)/2",
//                             "-pi*(cos(pi*x/2) - 1)*sin(pi*y/2)/2",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("-pi*pi*(cos(pi*y/2) - 1)*cos(pi*x/2)/4",
//                             "-pi*pi*(cos(pi*x/2) - 1)*cos(pi*y/2)/4",
//                             "0",
//                             "pi*pi*sin(pi*x/2)*sin(pi*y/2)/4",
//                             "0",
//                             "0", 3);





// ==============================================================================================================================================
/*
 *  QUADRATIC SURFACE (u, v, u * v)
 */


//    gsFunctionExpr<> source  ("(8 * (3 * x^10 + 3 * (1 + y^2)^3 + 12 * x^2 * (1 + y^2)^2 * (1 + 5 y^2) + \n"
//                              "   x^8 * (14 + 15 * y^2) + x^6 * (25 + 36 * y^2 + 5 * y^4) + \n"
//                              "   x^4 * (23 + 96 * y^2 + 178 * y^4 + 105 * y^6))) / (1 + x^2 + y^2)^5 ", 3);
////    gsFunctionExpr<> source  ("x^4", 3);
//    gsFunctionExpr<> laplace ("0 * (12 * (x^2 + x^4) ) / (1 + x^2 + y^2)",3);
//    gsFunctionExpr<> solVal("x^4",3);
//    gsFunctionExpr<>sol1der ("(4 * x^3 * (1 + x^2)) / (1 + x^2 + y^2)",
//                             "-((4 * x^4 * y) / (1 + x^2 + y^2))",
//                             "4 * x^3 * ( -( (x^2 * y) / (1 + x^2 + y^2) ) + ( (1 + x^2) * y ) / (1 + x^2 + y^2) )", 3);
//
////    gsFunctionExpr<>sol2der ("(4 * (1 + x^2) * (3 * x^6 + 3 * x^2 * (1 + y^2) + x^4 * (6 + 7 * y^2) ) ) / (1 + x^2 + y^2)^3",
////                             "(4 * x^4 * (-1 + 4 * y^2 + 5 * y^4 + x^2 * (-1 + y^2) ) ) / (1 + x^2 + y^2)^3",
////                             "(4 * (x^4 + x^6 + 3 * x^2 * (y^2 + y^4) ) ) / (1 + x^2 + y^2)^3",
////                             "-( (4 * y * (x^7 + 4 * x^3 * (1 + y^2) + 5 * x^5 * (1 + y^2) ) ) / (1 + x^2 + y^2)^3)",
////                             "(4 * y * (3 * x^2 * (1 + y^2) + x^4 * (3 + 4 * y^2) ) ) / (1 + x^2 + y^2)^3",
////                             "-( (4 * y * (3 * x^7 + 5 * x^3 * (1 + y^2) + x^5 * (8 + 7 * y^2) ) ) / (1 + x^2 + y^2)^3)",
////                             "(4 * x^3 * (1 + x^2 - 3 * y^2 - 4 * y^4) ) / (1 + x^2 + y^2)^3",
////                             "(4 * y * (x^6 + 3 * x^2 * (1 + y^2) + x^4 * (4 + 5 * y^2) ) ) / (1 + x^2 + y^2)^3",
////                             "-( (4 * (x^7 + x^5 * (1 + y^2) + 4 * x^3 * (y^2 + y^4) ) ) / (1 + x^2 + y^2)^3)", 3);
//
//    gsFunctionExpr<>sol2der ("(12 * (x + x^3)^2 ) / (1 + x^2 + y^2)^2",
//                             "(12 * x^4 * y^2) / (1 + x^2 + y^2)^2",
//                             "(12 * x^2 * y^2) / (1 + x^2 + y^2)^2",
//                             "-((12 * x^3 * (1 + x^2) * y) / (1 + x^2 + y^2)^2)",
//                             "(12 * x^2 * (1 + x^2) * y) / (1 + x^2 + y^2)^2",
//                             "-((12 * x^3 * y^2) / (1 + x^2 + y^2)^2)", 3);


//    gsFunctionExpr<> source  ("0.00000001 * 4/9 *(290409916 + 125 *x^6 - 70748808 *y - 33996207 *y^2 + 5199403 *y^3 + \n"
//                              "   448683 *y^4 - 46191 *y^5 + 636 *y^6 + 3 *x^5 *(-9226 + 431 *y) - \n"
//                              "   3 *x^4 *(-261466 + 32525 *y + 4810 *y^2) + \n"
//                              "   x^3 *(1789314 - 926955 *y + 114060 *y^2 + 22270 *y^3) - \n"
//                              "   3 *x^2 *(13649915 - 1431401 *y - 778947 *y^2 + 19250 *y^3 + 6525 *y^4) + \n"
//                              "   3 *x *(-6476056 + 7330526 *y - 243552 *y^2 - 344275 *y^3 - 6290 *y^4 + \n"
//                              "      2195 *y^5))",3);
//
//    gsFunctionExpr<> laplace ("0",3);
//
//    gsFunctionExpr<> solVal("0.00000001 * (16/3 - x/3 - y) * (16 - 3*x - y) * (5 - x) * (x - 7 - y) * (x/3 - 5 - y) * (-5 - x/2 - y) * (-19 - 4*x - y) * (11 + 2*x - y) * (7 + x - y) * (5 - y)",3);
//
////        gsFunctionExpr<> solVal("0",3);
//
//    gsFunctionExpr<>sol1der ("0.00000001 *  1/18 *(y - 5) *(216 *x^8 + 16 *x^7 *(y - 253) - 7* x^6 *(293 *y^2 + 370 *y + 4089) + 6 *x^5 *(23 *y^3 + 440 *y^2 + 11251 *y + 109390) \n"
//                             "+ 5 *x^4 *(712 *y^4 + 1789 *y^3 - 696 *y^2 - 27219 *y + 242206) \n"
//                             "- 4 *x^3 *(376 *y^5 + 2136 *y^4 - 28431 *y^3 - 162166 *y^2 + 971107 *y + 7561994) - 3 *x^2 *(173 *y^6 - 1476 *y^5 + 19462 *y^4 + 137317 *y^3 - 479192 *y^2 - 2096665 *y + 5849981) \n"
//                             "+ 2 *x *(63 *y^7 - 228 *y^6 + 10450 *y^5 + 155112 *y^4 - 1106231 *y^3 - 10292682 *y^2 + 25363606 *y + 207550710) + 18 *y^8 - 375 *y^7 - 2665 *y^6 - 5646 *y^5 + 360108 *y^4 \n"
//                             "+ 1678929 *y^3 - 8300949 *y^2 - 26433420 *y + 47824000)",
//                             "0.00000001 *  1/18 *(x - 5) *(24 *x^8 + 4 *x^7 *(y - 99) + x^6 *(-879 *y^2 + 2210 *y - 4219) + 4 *x^5 *(23 *y^3 - 855 *y^2 + 7288 *y + 8010) \n"
//                             "+ x^4 *(3560 *y^4 - 6624 *y^3 - 46023 *y^2 + 98282 *y + 538501) \n"
//                             "- 2 *x^3 *(1128 *y^5 - 8260 *y^4 - 61662 *y^3 + 85041 *y^2 + 1536232 *y + 6977) + x^2 *(-1211 *y^6 + 2766 *y^5 - 51610 *y^4 + 456592 *y^3 + 2646921 *y^2 - 15960910 *y - 16403076) \n"
//                             "+ 2 *x *(252 *y^7 - 4928 *y^6 + 41685 *y^5 + 128130 *y^4 - 2622102 *y^3 - 524988 *y^2 + 36924741 *y - 641350) + 6 *(27 *y^8 - 200 *y^7 - 9135 *y^6 + 77154 *y^5 + 537165 *y^4 \n"
//                             "- 4451244 *y^3 - 9222777 *y^2 + 66565010 *y + 28929600))",
//                             "0",3);

//    gsFunctionExpr<>sol1der ("0",
//                             "0",
//                             "0",3);

//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);


// Residual

//    gsFunctionExpr<> source  ("5",3);
//    gsFunctionExpr<> laplace ("0",3);
//    gsFunctionExpr<> solVal("0",3);
//    gsFunctionExpr<>sol1der ("0",
//                             "0",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("0",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);


//    gsFunctionExpr<> source  ("z",3);
//    gsFunctionExpr<> laplace ("6 * x",3);
//    gsFunctionExpr<> solVal("z",3);
//    gsFunctionExpr<>sol1der ("3 * x^2",
//                             "0",
//                             "0",3);
//    gsFunctionExpr<>sol2der ("6 * x",
//                             "0",
//                             "0",
//                             "0",
//                             "0",
//                             "0", 3);


// ==============================================================================================================================================


//    gsFunctionExpr<> source  ("(x - 1)^2 * (y - 1)^2 * x^2 * y^2",3);
////    gsFunctionExpr<> source  ("8 *(3 *x^4 - 6 *x^3 + 9 *x^2 *(1 - 2 *y)^2 - 6 *x *(6 *y^2 - 6 *y + 1) + 3 *y^4 - 6 *y^3 + 9 *y^2 - 6 *y + 1)",3);
//    gsFunctionExpr<> laplace ("2 * (6 * x^2 - 6 * x + 1) * (y - 1)^2 * y^2 + 2 * (6 * y^2 - 6 * y + 1) * (x - 1)^2 * x^2",3);
//
//    gsFunctionExpr<> solVal("(x - 1)^2 * (y - 1)^2 * x^2 * y^2",3);
//
//    gsFunctionExpr<>sol1der ("2 * (x - 1) * x *(2 * x - 1) * (y - 1)^2 * y^2",
//                             "2 * (x - 1)^2 * x^2 * (y - 1) * y * (2 * y - 1)",
//                             "0",3);
//
//    gsFunctionExpr<>sol2der ("2 * (6 * x^2 - 6 * x + 1) * (y - 1)^2 * y^2",
//                             "2 * (6 * y^2 - 6 * y + 1) * (x - 1)^2 * x^2",
//                             "0",
//                             "4 * (x - 1) * x * (2 * x - 1) * (y - 1) * y * (2 * y - 1)",
//                             "0",
//                             "4 * (x - 1) * x * (2 * x - 1) * (y - 1) * y * (2 * y - 1)",
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


//    gsFunctionExpr<> source  ("(64 * x * (7 + 16 * x^6 + 2 * y^2 + 12 * y^4 + 16 * y^6 + 12 * x^4 * (1 + 4 * y^2) + x^2 * (2 + 24 * y^2 + 48 * y^4) ) ) / (1 + 4 * x^2 + 4 * y^2)^5",3);
////    gsFunctionExpr<> source  ("x",3);
//
//
//    gsFunctionExpr<> laplace ("-( (2 * x * y) / (1 + x^2 + y^2))",3);
//
//    gsFunctionExpr<> solVal("x",3);
//
//    gsFunctionExpr<>sol1der ("(1 + 4 * y^2) / (1 + 4 * x^2 + 4 * y^2)",
//                             "-((4 * x * y) / (1 + 4 * x^2 + 4 * y^2))",
//                             "-((2 * x) / (1 + 4 * x^2 + 4 * y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);



//    gsFunctionExpr<> source  ("(8 * (5120 * x^10 - 4608 * x^8 * (-1 + 2 * y^2) + 3 * (1 + 4 * y^2)^5 - \n"
//                              "   48 * x^2 * (1 + 4 * y^2)^3 * (2 + 5 * y^2) - \n"
//                              "   64 * x^6 * (-19 + 252 * y^2 + 864 * y^4) - \n"
//                              "   16 * x^4 * (7 + 512 * y^2 + 2640 * y^4 + 3712 * y^6)) ) / (1 + 4 * x^2 + 4 * y^2)^5",3);
//    gsFunctionExpr<> source  ("x^4",3);
//
//
//    gsFunctionExpr<> laplace ("-( (2 * x * y) / (1 + x^2 + y^2))",3);
//
//    gsFunctionExpr<> solVal("x^4",3);
//
//    gsFunctionExpr<>sol1der ("(4 * x^3 * (1 + 4 * y^2) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "-((16 * x^4 * y) / (1 + 4 * x^2 + 4 * y^2))",
//                             "-((8 * x^4) / (1 + 4 * x^2 + 4 * y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


//    gsFunctionExpr<> source  ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^5 ) ) * 8 * (3072 * x^14 + 3 * y^4 * (1 + 4 * y^2)^5 - 768 * x^12 * (-5 + 132 * y^2) - \n"
//                              "   12 * x^2 * y^2 * (1 + 4 * y^2)^3 * (-3 + 28 * y^2 + 132 * y^4) - 128 * x^10 * (-15 + 762 * y^2 + 200 *y^4) + \n"
//                              "   96 * x^8 * (5 - 342 * y^2 + 568 * y^4 + 4704 * y^6) + 12 * x^6 * (5 - 324 * y^2 + 3472 * y^4 + 25984 * y^6 + 37632 * y^8) + \n"
//                              "   x^4 * (3 + 96 * y^2 + 7392 * y^4 + 41664 * y^6 + 54528 * y^8 - 25600 * y^10))",3);
////    gsFunctionExpr<> source  ("x^4 * y^4",3);
//
//    gsFunctionExpr<> laplace ("(4 * x^2 * y^2 * (48 * x^6 - 8 * x^4 * (-3 + 14 * y^2) + 3 * (y + 4 * y^3)^2 + x^2 * (3 - 24 * y^2 - 112 * y^4) ) ) / (1 + 4 * x^2 + 4 * y^2)^2",3);
//
//    gsFunctionExpr<> solVal("x^4 * y^4",3);
//
//    gsFunctionExpr<>sol1der ("(4 * x^3 * y^4 * (1 - 4 * x^2 + 4 * y^2) ) / ( 1 + 4 * x^2 + 4 * y^2)",
//                             "(4 * x^4 * y^3 * (1 + 4 * x^2 - 4 * y^2) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "-( ( 16 * x^4 * y^4) / (1 + 4 * x^2 + 4 * y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


//    gsFunctionExpr<> source  ("(8 * y * (30848 * x^10 - 192 * x^2 * (1 + 3 * y^2) * (1 + 4 * y^2)^3 + \n"
//                              "   3 * (1 + 4 * y^2)^5 + 96 * x^8 * (289 + 484 * y^2) - \n"
//                              "   48 * x^6 * (-159 - 164 * y^2 + 888 * y^4) - \n"
//                              "   8 * x^4 * (-15 + 1302 * y^2 + 7956 * y^4 + 12272 * y^6) ) ) / (1 + 4 * x^2 + 4 * y^2)^5",3);
////    gsFunctionExpr<> source  ("x^4 * y^4",3);
//
//    gsFunctionExpr<> laplace ("-((4 * x^2 * y * (52 * x^4 - 3 * (1 + 4 * y^2)^2 + x^2 * (6 + 4 * y^2) ) ) / (1 + 4 * x^2 + 4 * y^2)^2)",3);
//
//    gsFunctionExpr<> solVal("x^4 * y",3);
//
//    gsFunctionExpr<>sol1der ("-( (4 * x^3 * y * (-1 + x^2 - 4 * y^2) ) / (1 + 4 * x^2 + 4 * y^2))",
//                             "(x^4 * (1 + 4 * x^2 - 16 * y^2) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "-((10 * x^4 * y) / (1 + 4 * x^2 + 4 * y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


//    gsFunctionExpr<> source  ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * 16 * (3 + 4096 * x^10 - 18 * y^2 + 184 * y^4 + 1568 * y^6 + \n"
//                              "   4224 * y^8 + 4096 * y^10 - 384 * x^8 * (-11 + 32 * y^2) - 32 * x^6 * (-49 + 528 * y^2 + 1792 * y^4) - \n"
//                              "   8 * x^4 * (-23 + 980 * y^2 + 5280 * y^4 + 7168 * y^6) - 2 * x^2 * (9 + 696 * y^2 + 3920 * y^4 + 8448 * y^6 + 6144 * y^8))",3);
////    gsFunctionExpr<> source  ("x^4 * y^4",3);
//
//    gsFunctionExpr<> laplace ("(4 * (-16 * x^6 + 3 * y^2 + 4 * y^4 - 16 * y^6 + x^4 * (4 + 80 * y^2) + x^2 * (3 + 48 * y^2 + 80 * y^4) ) ) / (1 + 4 * x^2 + 4 * y^2)^2",3);
//
//    gsFunctionExpr<> solVal("x^4 + y^4",3);
//
//    gsFunctionExpr<>sol1der ("(4 * (-4 * x * y^4 + x^3 * (1 + 4 * y^2) ) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(4 * y * (-4 * x^4 + y^2 + 4 * x^2 * y^2) ) / ( 1 + 4 * x^2 + 4 * y^2)",
//                             "-((8 * (x^4 + y^4) ) / (1 + 4 * x^2 + 4 * y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);



//    gsFunctionExpr<> source  ("(1000 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * 40 * x * y * (3072 * x^14 + 3 * y^4 * (1 + 4 * y^2)^5 - "
//                              "256 * x^12 * (-15 + 188 * y^2) - 384 * x^10 * (-5 + 118 * y^2 + 72 * y^4) - "
//                              "4 * x^2 * y^2 * (1 + 4 * y^2)^3 * (-5 + 36 * y^2 + 188 * y^4) + 32 * x^8 * (15 - 458 * y^2 + 336 * y^4 + 5600 * y^6) + "
//                              "4 * x^6 * (15 - 380 * y^2 + 3608 * y^4 + 29952 * y^6 + 44800 * y^8) + "
//                              "x^4 * (3 + 96 * y^2 + 2952 * y^4 + 14432 * y^6 + 10752 * y^8 - 27648 * y^10))",3);
////    gsFunctionExpr<> source  ("x^4 * y^4",3);
//
//    gsFunctionExpr<> laplace ("(1000 * 20 * x^3 * y^3 * (16 * x^6 + x^4 * (8 - 32 * y^2) + (y + 4 * y^3)^2 + x^2 * (1 - 6 * y^2 - 32 * y^4) ) ) / (1 + 4 * x^2 + 4 * y^2)^2",3);
//
//    gsFunctionExpr<> solVal("1000 * x^5 * y^5",3);
//
//    gsFunctionExpr<>sol1der ("(1000 * 5 * x^4 * y^5 * (1 - 4 * x^2 + 4 * y^2) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(1000 * 5 * x^5 * y^4 * (1 + 4 * x^2 - 4 * y^2) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "-1000 * ( (20 * x^5 * y^5) / (1 + 4 * x^2 + 4 * y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


//    gsFunctionExpr<> source  ("(40 * (1664 * x^11 + 3 * x * (1 + 4 * y^2)^5 - 32 * x^9 * (-39 + 196 * y^2) - \n"
//                              "   8 * (3 + 4 * y^2) * (x + 4 * x * y^2)^3 - 48 * x^7 * (-3 + 172 * y^2 + 536 * y^4) - \n"
//                              "   8 * x^5 * (13 + 414 * y^2 + 2052 * y^4 + 2864 * y^6) ) ) / (1 + 4 * x^2 + 4 * y^2)^5",3);
////    gsFunctionExpr<> source  ("x^4 * y^4",3);
//
//    gsFunctionExpr<> laplace ("(20 * x^3 * (-4 * x^4 + (1 + 4 * y^2)^2 + 2 * x^2 * (1 + 6 * y^2) ) ) / (1 + 4 * x^2 + 4 * y^2)^2",3);
//
//    gsFunctionExpr<> solVal("x^5",3);
//
//    gsFunctionExpr<>sol1der ("(5 * x^4 * (1 + 4 * y^2) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "-( (20 * x^5 * y) / (1 + 4 * x^2 + 4 * y^2))",
//                             "-( (10 * x^5) / (1 + 4 * x^2 + 4 * y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


// HYPERBOLOID

//    gsFunctionExpr<> source  ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * 4 * ( (1 + 4 * x^2 + 4 * y^2) * (357 + 5184 * x^8 + 4308 * y^2 + \n"
//                              "      24724 * y^4 + 51872 * y^6 + 40000 * y^8 + 96 * x^6 * (103 + 1008 * y^2) + 4 * x^4 * (2821 + 31704 * y^2 + 54496 * y^4) + \n"
//                              "      4 * x^2 * (453 + 4206 * y^2 + 22744 * y^4 + 41600 * y^6) ) * cos(5*x) * cos(1 - 3* y) - \n"
//                              "   80 * x * (5 - 864 * x^8 + 141 * y^2 + 782 * y^4 + 1808 * y^6 + 672 * y^8 - 48 * x^6 * (23 + 32 * y^2) + 6 * x^4 * (-51 + 168 * y^2 + 32 * y^4) + \n"
//                              "      x^2 * (-29 + 4 * y^2 * (89 + 84 * y^2 + 384 * y^4) ) ) * cos(1 - 3 * y) * sin(5 * x) - 24 * y * (2 * (3 + 3936 * x^8 + 125 * y^2 + 786 * y^4 + 2384 * y^6 + \n"
//                              "         2400 * y^8 + 48 * x^6 * (21 + 32 * y^2) - 2 * x^4 * (55 + 72 * y^2 * (-3 + 44 * y^2) ) - \n"
//                              "         x^2 * (45 + 4 * y^2 * (-7 + 444 * y^2 + 384 * y^4) ) ) * cos(5 * x) - 5 * x * (1 + 4 * x^2 + 4 * y^2) * (31 + 288 * x^6 + 242 * y^2 + 976 * y^4 + \n"
//                              "         800 * y^6 + 16 * x^4 * (45 + 86 * y^2) + 2 * x^2 * (105 + 240 * y^2 + 944 * y^4) ) * sin(5 * x) ) * sin(1 - 3 * y))",3);
////    gsFunctionExpr<> source  ("cos(5 * x) * cos(1 - 3 * y)",3);
//
//    gsFunctionExpr<> laplace ("-(1 / ( (1 + 4 * x^2 + 4 * y^2)^2 ) ) * 2 * ( cos(5 * x) * ( (17 + 72 * x^4 + 118 * y^2 + 200 * y^4 + \n"
//                              "          x^2 * (86 + 272 * y^2) ) * cos(1 - 3 * y) + 24 * y * (x^2 - y^2) * sin(1 - 3 * y) ) + \n"
//                              "          20 * x * sin(5 * x) * (2 * (x^2 - y^2) * cos(1 - 3 * y) + 3 * y * (1 + 4 * x^2 + 4 * y^2) * sin(1 - 3 * y) ) )",3);
//
//    gsFunctionExpr<> solVal("cos(5 * x) * cos(1 - 3 * y)",3);
//
//    gsFunctionExpr<>sol1der ("(-5 * (1 + 4 * y^2) * cos(1 - 3 * y) * sin(5 * x) + 12 * x * y * cos(5 * x) * sin(1 - 3 * y) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(-20 * x * y * cos(1 - 3 * y) * sin(5 * x) + 3 * (1 + 4 * x^2) * cos(5 * x) * sin(1 - 3 * y) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(10 * x * (1 + 8 * y^2) * cos(1 - 3 * y) * sin(5 * x) - 6 * (1 + 8 * x^2) * y * cos(5 * x) * sin(1 - 3 * y) ) / (1 + 4 * x^2 + 4 * y^2)",3);
////
////
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


    gsFunctionExpr<> source  ("(64 * (x - y) * (x + y) * (13 + 64 * x^4 + 12 * y^2 + 64 * y^4 + x^2 * (12 - 768 * y^2)))/(1 + 4 * x^2 + 4 * y^2)^5",3);
//    gsFunctionExpr<> source  ("z",3);

    gsFunctionExpr<> laplace ("-((8 * (x^2 - y^2))/(1 + 4 * x^2 + 4 * y^2)^2)",3);

    gsFunctionExpr<> solVal("z",3);

    gsFunctionExpr<>sol1der ("(2 * x)/(1 + 4 * x^2 + 4 * y^2)",
                             "-((2 * y)/(1 + 4 * x^2 + 4 * y^2))",
                             "((4 * (x^2 + y^2))/(1 + 4 * x^2 + 4 * y^2))",3);
//
//
    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);




//    gsFunctionExpr<> source  ("0.001 * (1 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * 24 * (3072 * x^16 * (-4 + 5 * y^2) - \n"
//                              "   256 * x^14 * (124 - 1611 * y^2 + 620 * y^4) - 128 * x^12 * (-4564 + 9877 * y^2 - 2823 * y^4 + 968 * y^6) + \n"
//                              "   32 * x^10 * (-20556 - 80909 * y^2 + 31023 * y^4 - 53064 * y^6 + 16544 * y^8) + 4 * x^8 * (-273884 + 399675 * y^2 - 2031058 * y^4 + 2325128 * y^6 - \n"
//                              "      817344 * y^8 + 132352 * y^10) + x^6 * (-470828 + 4289791 * y^2 + 10546230 * y^4 - 9879856 * y^6 + \n"
//                              "      9300512 * y^8 - 1698048 * y^10 - 123904 * y^12) + x^4 * (-203904 + 2214468 * y^2 + 11041851 * y^4 + 10546230 * y^6 - \n"
//                              "      8124232 * y^8 + 992736 * y^10 + 361344 * y^12 - 158720 * y^14) - 4 * (-3648 + 26656 * y^2 + 50976 * y^4 + 117707 * y^6 + 273884 * y^8 + \n"
//                              "      164448 * y^10 - 146048 * y^12 + 7936 * y^14 + 3072 * y^16) + x^2 * (-106624 + 124512 * y^2 + 2214468 * y^4 + 4289791 * y^6 + \n"
//                              "      1598700 * y^8 - 2589088 * y^10 - 1264256 * y^12 + 412416 * y^14 + 15360 * y^16))",3);
//    gsFunctionExpr<> source  ("0.001 * (2 + x)^3 * (2 - x)^3 * (2 + y)^3 * (2 - y)^3",3);
//
//    gsFunctionExpr<> laplace ("0.001 * (1 / ( (1 + 4 * x^2 + 4 * y^2)^2) ) * 6 * (-4 + x^2) * (-4 + y^2) * (16 * x^8 * (-4 + 5 * y^2) - \n"
//                              "   8 * x^6 * (-28 - 37 * y^2 + 18 * y^4) + x^4 * (444 - 1483 * y^2 + 616 * y^4 - 144 * y^6) - \n"
//                              "   4 * (32 + 36 * y^2 - 111 * y^4 - 56 * y^6 + 16 * y^8) + x^2 * (-144 + 176 * y^2 - 1483 * y^4 + 296 * y^6 + 80 * y^8))",3);
//
//    gsFunctionExpr<> solVal("0.001 * (2 + x)^3 * (2 - x)^3 * (2 + y)^3 * (2 - y)^3",3);
//
//    gsFunctionExpr<>sol1der ("-0.001 * ( (6 * x * (-4 + x^2)^2 * (-4 + y^2)^2 * (4 + (-1 + 4 * x^2) * y^2 - 4 * y^4) ) / (1 + 4 * x^2 + 4 * y^2))",
//                             "0.001 * (6 * (-4 + x^2)^2 * y * (-4 + y^2)^2 * (-4 + 4 * x^4 + x^2 * (1 - 4 * y^2) ) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "-0.001 * ( (24 * (-4 + x^2)^2 * (-4 + y^2)^2 * (-2 * y^2 + x^2 * (-2 + y^2) ) ) / (1 + 4 * x^2 + 4 * y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);

//    gsFunctionExpr<> source  ("0.001 * (1 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * 8 * (1696 + 3072 * x^14 - 12296 * y^2 - 22189 * y^4 - 54020 * y^6 - \n"
//                              "   136992 * y^8 - 109696 * y^10 + 36608 * y^12 + 3072 * y^14 - 2816 * x^12 * (-13 + 36 * y^2) - 128 * x^10 * (857 + 1274 * y^2 + 200 * y^4) + \n"
//                              "   32 * x^8 * (-4281 + 9438 * y^2 - 15704 * y^4 + 14112 * y^6) + 4 * x^6 * (-13505 + 135540 * y^2 + 401072 * y^4 - 151424 * y^6 + \n"
//                              "      112896 * y^8) + x^4 * (-22189 + 265696 * y^2 + 1374432 * y^4 + 1604288 * y^6 - \n"
//                              "      502528 * y^8 - 25600 * y^10) - 4 * x^2 * (3074 - 4233 * y^2 - 66424 * y^4 - 135540 * y^6 - 75504 * y^8 + \n"
//                              "      40768 * y^10 + 25344 * y^12))",3);
////    gsFunctionExpr<> source  ("x^4 * y^4",3);
//
//    gsFunctionExpr<> laplace ("0.001 * (1 / ( (1 + 4 * x^2 + 4 * y^2)^2) ) * 4 * (16 * x^8 * (-4 + 3 * y^2) - 56 * x^6 * (-4 - 5 * y^2 + 2 * y^4) + \n"
//                              "   x^4 * (316 - 1357 * y^2 + 616 * y^4 - 112 * y^6) - 4 * (32 + 44 * y^2 - 79 * y^4 - 56 * y^6 + 16 * y^8) + \n"
//                              "   x^2 * (-176 - 48 * y^2 - 1357 * y^4 + 280 * y^6 + 48 * y^8))",3);
//
//    gsFunctionExpr<> solVal("0.001 * (2 + x)^2 * (2 - x)^2 * (2 + y)^2 * (2 - y)^2",3);
//
//    gsFunctionExpr<>sol1der ("-0.001 * ( (4 * x * (-4 + y^2) * (4 * x^4 * y^2 + x^2 * (4 - 17 * y^2 - 4 * y^4) + 4 * (-4 + y^2 + 4 * y^4) ) ) / (1 + 4 * x^2 + 4 * y^2))",
//                             "0.001 * (4 * (-4 + x^2) * y * (-4 + y^2) * (-4 + 4 * x^4 + x^2 * (1 - 4 * y^2) ) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "-0.001 * ( (16 * (-4 + x^2) * (-4 + y^2) * (-2 * y^2 + x^2 * (-2 + y^2) ) ) / (1 + 4 * x^2 + 4 * y^2))",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


//    gsFunctionExpr<> source  ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * (64 * x * sin(x) * (-4 * (4 + 7 * x^2 + 33 * x^4 + 64 * x^6 + 48 * x^8 + \n"
//                              "        4 * (1 + 3 * x^2 * (3 + 8 * (x^2 + x^4) ) ) * y^2 + 3 * y^4 - 32 * (1 + 3 * x^2) * y^6 - 48 * y^8) * cos(y) + \n"
//                              "     y * (1 + 4 * x^2 + 4 * y^2) * (39 + 114 * x^2 + 144 * x^4 + 32 * x^6 + 6 * (19 + 16 * x^2 * (3 + x^2) ) * y^2 + 48 * (3 + 2 * x^2) y^4 + \n"
//                              "        32 * y^6) * sin(y) ) + 8 * cos(x) * ( (1 + 4 x^2 + 4 y^2) * (13 + 64 * x^8 + 32 * x^6 * (5 + 16 * y^2) + \n"
//                              "        4 * y^2 * (1 + y^2) * (-3 + 24 * y^2 + 16 * y^4) + x^4 * (84 + 608 * y^2 + 896 * y^4) + \n"
//                              "        4 * x^2 * (-3 + 46 * y^2 + 152 * y^4 + 128 * y^6) ) * cos(y) - 32 * y * (4 + 3 * x^4 - 48 * x^8 + 7 * y^2 + 33 * y^4 + 64 * y^6 + 48 * y^8 - \n"
//                              "        32 * x^6 * (1 + 3 * y^2) + 4 * x^2 * (1 + 3 * y^2 * (3 + 8 * (y^2 + y^4) ) ) ) * sin(y)))",3);
//gsFunctionExpr<> source  ("cos(pi * x / 2) * cos(pi * y / 2) * cos(pi * z / 2)",3);
//
//    gsFunctionExpr<> laplace ("-(1 / ( (1 + 4 * x^2 + 4 * y^2)^2) ) * 4 * ( (1 + 2 * x^2 + 2 * y^2) * cos(x) * ( (1 + "
//                              "4 * x^2 + 4 * y^2) * cos(y) - 4 * y * sin(y)) + 4 * x * sin(x) * (-(1 + 2 * x^2 + 2 * y^2) * cos(y) + \n"
//                              "       y * (1 + 4 * x^2 + 4 * y^2) * sin(y)))",3);
//
//    gsFunctionExpr<> solVal("cos(x) * cos(y)",3);
//
//    gsFunctionExpr<>sol1der ("(-2 * (1 + 4 * y^2) * cos(y) * sin(x) + 8 * x * y * cos(x) * sin(y)) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(8 * x * y * cos(y) * sin(x) - 2 * (1 + 4 * x^2) * cos(x) * sin(y)) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(4 * (x * cos(y) * sin(x) + y * cos(x) * sin(y)) ) / (1 + 4 * x^2 + 4 * y^2)",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


//    gsFunctionExpr<> source  ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * 40 * (8 * x * sin(5 * x) * (-4 * (16 + 151 * x^2 + 681 * x^4 + 1408 * x^6 + 1200 * x^8 + \n"
//                              "         4 * (19 + 153 * x^2 + 456 * x^4 + 600 * x^6) * y^2 - 3 * (23 + 192 * x^2) * y^4 - 32 * (31 + 75 * x^2) * y^6 - 1200 * y^8) * cos(5 * y) + "
//                              "         5 * y * (1 + 4 * x^2 + 4 * y^2) * (63 + 354 * x^2 + 912 * x^4 + 800 * x^6 + 6 * (59 + 304 * x^2 + 400 * x^4) * y^2 + \n"
//                              "         48 * (19 + 50 * x^2) * y^4 + 800 * y^6) * sin(5 * y) ) + cos(5 * x) * (5 * (1 + 4 * x^2 + 4 * y^2) * (37 + 276 * y^2 + \n"
//                              "         4 * (400 * x^8 + 8 * x^6 * (77 + 400 * y^2) + y^4 * (333 + 616 * y^2 + 400 * y^4) + \n"
//                              "         x^4 * (333 + 2648 * y^2 + 5600 * y^4) + x^2 * (69 + 766 * y^2 + 2648 * y^4 + 3200 * y^6) ) ) * cos(5 * y) - \n"
//                              "         32 * y * (16 - 1200 * x^8 + 151 * y^2 + 681 * y^4 - 32 * x^6 * (31 + 75 * y^2) + 16 * y^6 * (88 + 75 * y^2) - \n"
//                              "         3 * x^4 * (23 + 192 * y^2) + 4 * x^2 * (19 + 153 * y^2 + 456 * y^4 + 600 * y^6) ) * sin(5 * y)))",3);
////    gsFunctionExpr<> source  ("2 * cos(5 * x) * cos(5 * y)",3);
//
//    gsFunctionExpr<> laplace ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^2) ) * (-20 * (1 + 2 * x^2 + 2 * y^2) * cos(5 * x) * (5 * (1 + 4 * x^2 + 4 * y^2) * cos(5 * y) - "
//                              "               4 * y * sin(5 * y) ) + 80 * x * sin(5 * x) * ( (1 + 2 * x^2 + 2 * y^2) * cos(5 * y) - 5 * y * (1 + 4 * x^2 + 4 * y^2) * sin(5 * y) ) )",3);
//
//    gsFunctionExpr<> solVal("2 * cos(5 * x) * cos(5 * y)",3);
//
//    gsFunctionExpr<>sol1der ("(-10 * (1 + 4 * y^2) * cos(5 * y) * sin(5 * x) + 40 * x * y * cos(5 * x) * sin(5 * y) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(40 * x * y * cos(5 * y) * sin(5 * x) - 10 * (1 + 4 * x^2) * cos(5 * x) * sin(5 * y) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(20 * (x * cos(5 * y) * sin(5 * x) + y * cos(5 * x) * sin(5 * y) ) ) / (1 + 4 * x^2 + 4 * y^2)",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


//    gsFunctionExpr<> source  ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * 1024 * (cos(8 * x) * (2 (1 + 4 * x^2 + 4 * y^2) * (19 + 186 * y^2 + \n"
//                              "             2 * (512 * x^8 + 420 * y^4 + 776 * y^6 + 512 * y^8 + \n"
//                              "            8 * x^6 * (97 + 512 * y^2) + 4 * x^4 * (105 + 838 * y^2 + 1792 * y^4) + x^2 * (93 + 968 * y^2 + 3352 * y^4 + 4096 * y^6) ) ) * cos(8 * y) - \n"
//                              "      y * (71 - 6144 * x^8 + 770 * y^2 + 3468 * y^4 + 7184 * y^6 + 6144 * y^8 - 12 * x^4 * (31 + 252 * y^2) - 16 * x^6 * (319 + 768 * y^2) + \n"
//                              "         2 * x^2 * (193 + 1548 * y^2 + 4632 * y^4 + 6144 * y^6) ) * sin(8 * y) ) + x * sin(8 * x) * (-(71 + 6144 * x^8 + 386 * y^2 - \n"
//                              "          4 * y^4 * (3 + 4 * y^2) * (31 + 384 * y^2) + 16 * x^6 * (449 + 768 * y^2) + 12 * x^4 * (289 + 772 * y^2) + \n"
//                              "          x^2 * (770 - 24 * y^2 * (-129 + 126 * y^2 + 512 * y^4) ) ) * cos(8 * y) + 8 * y * (1 + 4 * x^2 + 4 * y^2) * (51 + 372 * y^2 + \n"
//                              "         4 * (256 * x^6 + 270 * y^4 + 256 * y^6 + 6 * x^4 * (45 + 128 * y^2) + x^2 * (93 + 540 * y^2 + 768 * y^4) ) ) * sin(8 * y)))",3);
////    gsFunctionExpr<> source  ("2 * cos(8 * x) * cos(8 * y)",3);
//
//    gsFunctionExpr<> laplace ("-(1 / ( (1 + 4 * x^2 + 4 * y^2)^2 ) ) * 128 * ( (1 + 2 * x^2 + 2 * y^2) * cos(8 * x) * ( (2 + 8 * x^2 + \n"
//                              "               8 * y^2) * cos(8 * y) - y * sin(8 * y)) + x * sin(8 * x) * (-(1 + 2 * x^2 + 2 * y^2) * cos(8 * y) + \n"
//                              "               8 * y * (1 + 4 * x^2 + 4 * y^2) * sin(8 * y)))",3);
//
//    gsFunctionExpr<> solVal("2 * cos(8 * x) * cos(8 * y)",3);
//
//    gsFunctionExpr<>sol1der ("(1 / (1 + 4 * x^2 + 4 * y^2) ) * (-16 * (1 + 4 * y^2) * cos(8 * y) * sin(8 * x) + 64 * x * y * cos(8 * x) * sin(8 * y))",
//                             "(1 / (1 + 4 * x^2 + 4 * y^2) ) * (64 * x * y * cos(8 * y) * sin(8 * x) - 16 * (1 + 4 * x^2) * cos(8 * x) * sin(8 * y))",
//                             "(32 * (x * cos(8 * y) * sin(8 * x) + y * cos(8 * x) * sin(8 * y)) ) / (1 + 4 * x^2 + 4 * y^2)",3);
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


//    gsFunctionExpr<> source  ("(1/((1 + 4 * x^2 + 4 * y^2)^5) ) * 32 * (-32 * x * sin(8 * x) * ( (57 + 428 * y^2 + \n"
//                              "         4 * (864 * x^8 + 12 * x^6 * (89 + 200 * y^2) + 3 * y^4 * (81 + 4 * y^2 - 64 * y^4) + \n"
//                              "            3 * x^4 * (191 + 716 * y^2 + 672 * y^4) + 4 * x^2 * (35 + 204 * y^2 + 273 * y^4 + 72 * y^6) ) ) * cos(1 - 6 * y)+ \n"
//                              "      24 * y * (1 + 4 * x^2 + 4 * y^2) * (11 + 144 * x^6 + 79 * y^2 + 242 * y^4 + 256 * y^6 + 2 * x^4 * (93 + 272 * y^2) + \n"
//                              "         x^2 * (72 + 428 * y^2 + 656 * y^4) ) * sin(1 - 6 * y) ) + cos(8 * x) * ( (1 + 4 * x^2 + 4 * y^2) * (775 + 8432 * y^2 + \n"
//                              "         16 * (1296 * x^8 + 96 * x^6 * (25 + 171 * y^2) + y^4 * (2597 + 5396 * y^2 + 4096 * y^4) + \n"
//                              "            4 * x^4 * (385 + 3603 * y^2 + 8260 * y^4) + x^2 * (373 + 8 * y^2 * (583 + 64 * y^2 * (34 + 43 * y^2) ) ) ) ) * cos(1 - 6 * y) + \n"
//                              "      24 * y * (57 - 8832 * x^8 + 644 * y^2 + 3132 * y^4 + 6960 * y^6 + 6144 * y^8 - 48 * x^6 * (167 + 424 * y^2) + \n"
//                              "         16 * x^2 * (11 + 99 * y^2 + 369 * y^4 + 600 * y^6) - 36 * x^4 * (43 + 28 * y^2 * (9 + 8 * y^2) ) ) * sin(1 - 6 * y)))",3);
////    gsFunctionExpr<> source  ("2 * cos(8 * x) * cos(1 - 6 * y)",3);
//
//    gsFunctionExpr<> laplace ("(1/((1 + 4 * x^2 + 4 * y^2)^2) ) * (-8 * cos(8 * x) * ( (25 + 144 * x^4 + 164 * y^2 + 256 * y^4 + \n"
//                              "        8 * x^2 * (17 + 50 * y^2) ) * cos(1 - 6 * y) + 12 * y * (1 + 2 * x^2 + 2 * y^2) * sin(1 - 6 * y)) + \n"
//                              "  128 * x * sin(8 * x) * ( (1 + 2 * x^2 + 2 * y^2) * cos(1 - 6 * y) + 6 * y * (1 + 4 * x^2 + 4 * y^2) * sin(1 - 6 * y)))",3);
//
//    gsFunctionExpr<> solVal("2 * cos(8 * x) * cos(1 - 6 * y)",3);
//
//    gsFunctionExpr<>sol1der ("-((16 * ( (1 + 4 * y^2) * cos(1 - 6 * y) * sin(8 * x) + 3 * x * y * cos(8 * x) * sin(1 - 6 * y) ) ) / (1 + 4 * x^2 + 4 * y^2))",
//                             "(64 * x * y * cos(1 - 6 * y) * sin(8 * x) + 12 * (1 + 4 * x^2) * cos(8 * x) * sin(1 - 6 * y) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(8 * (4 * x * cos(1 - 6 * y) * sin(8 * x) - 3 * y * cos(8 * x) * sin(1 - 6 * y) ) ) / (1 + 4 * x^2 + 4 * y^2)",3);
////
////
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


//    gsFunctionExpr<> source  ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * 8 * ( (1 + 4 * x^2 + 4 * y^2) * (5005 + 55732 * y^2 + \n"
//                              "      4 * (38416 * x^8 + 67765 * y^4 + 24 * y^6 * (5815 + 4374 * y^2) + 8 * x^6 * (8677 + 57232 * y^2) + \n"
//                              "         x^4 * (43637 + 397208 * y^2 + 905440 * y^4) + x^2 * (10637 + 126254 * y^2 + 467352 * y^4 + 590976 * y^6) ) ) * cos(9 * x) * cos(1 - 7 * y) - \n"
//                              "   288 * x * (36 + 367 * x^2 + 1521 * x^4 + 2880 * x^6 + 2352 * x^8 + 4 * (67 + 513 * x^2 + 1368 * x^4 + 1560 * x^6) * y^2 + \n"
//                              "      9 * (59 + 256 * (x^2 + 2 * x^4) ) * y^4 - 96 * (3 + x^2) * y^6 - 816 * y^8) * cos(1 - 7 * y) * sin(9 * x) + \n"
//                              "   56 * y * (4 * (36 - 5424 * x^8 + 415 * y^2 + 2001 * y^4 + 4416 * y^6 + 3888 * y^8 - 288 * x^6 * (17 + 43 * y^2) + \n"
//                              "         4 * x^2 * (31 + 273 * y^2 + 984 * y^4 + 1560 * y^6) - 3 * x^4 * (303 + 256 * y^2 * (7 + 6 * y^2) ) ) * cos(9 * x) - \n"
//                              "      9 * x * (1 + 4 * x^2 + 4 * y^2) * (103 + 1568 * x^6 + 786 * y^2 + 2448 * y^4 + 2592 * y^6 + 16 * x^4 * (121 + 358 * y^2) + \n"
//                              "         x^2 * (722 + 4384 * y^2 + 6752 * y^4) ) * sin(9 * x) ) * sin(1 - 7 * y) )",3);
////    gsFunctionExpr<> source  ("2 * cos(9 * x) * cos(1 - 7 * y)",3);
//
//    gsFunctionExpr<> laplace ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^2) ) * (-4 * cos(9 * x) * ( (65 + 392 * x^4 + 422 * y^2 + 648 * y^4 + \n"
//                              "        2 * x^2 * (179 + 520 * y^2) ) * cos(1 - 7 * y) + 28 * y * (1 + 2 * x^2 + 2 * y^2) * sin(1 - 7 * y) ) + \n"
//                              "  144 * x * sin(9 * x) * ( (1 + 2 * x^2 + 2 * y^2) * cos(1 - 7 * y) + 7 * y * (1 + 4 * x^2 + 4 * y^2) * sin(1 - 7 * y) ) )",3);
//
//    gsFunctionExpr<> solVal("2 * cos(9 * x) * cos(1 - 7 * y)",3);
//
//    gsFunctionExpr<>sol1der ("-( (2 * (9 * (1 + 4 * y^2) * cos(1 - 7 * y) * sin(9 * x) + 28 * x * y * cos(9 * x) * sin(1 - 7 * y) ) ) / (1 + 4 * x^2 + 4 * y^2))",
//                             "(72 * x * y * cos(1 - 7 * y) * sin(9 * x) + 14 * (1 + 4 * x^2) * cos(9 * x) * sin(1 - 7 * y) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(4 * (9 * x * cos(1 - 7 * y) * sin(9 * x) - 7 * y * cos(9 * x) * sin(1 - 7 * y) ) ) / (1 + 4 * x^2 + 4 * y^2)",3);
//
//
//    gsFunctionExpr<>sol2der ("-( ( 2 * x * (1 + x^2) * y ) / (1 + x^2 + y^2)^2)",
//                             "-( ( 2 * x * y * (1 + y^2) ) / (1 + x^2 + y^2)^2)",
//                             "(2 * x * y ) / (1 + x^2 + y^2)^2",
//                             "(1 + y^2 + x^2 * (1 + 2 * y^2) ) / (1 + x^2 + y^2)^2",
//                             "(x * (1 + x^2 - y^2) ) / (1 + x^2 + y^2)^2",
//                             "(y * (1 - x^2 + y^2) ) / (1 + x^2 + y^2)^2", 3);


//    gsFunctionExpr<> source  ("(1 / ( (1 + 4 * x^2 + 4 * y^2)^5) ) * 1280 * (cos(10 * x) * (5 * (1 + 4 * x^2 + 4 * y^2) * (14 + 147 * x^2 + 654 * x^4 + \n"
//                              "         1208 * x^6 + 800 * x^8 + (147 + 1508 * x^2 + 5224 * x^4 + 6400 * x^6) * y^2 + \n"
//                              "         2 * (327 + 2612 * x^2 + 5600 * x^4) * y^4 + 8 * (151 + 800 * x^2) * y^6 + \n"
//                              "         800 * y^8) * cos(10 * y) - y * (107 + 602 * x^2 - 4 * x^4 * (3 + 4 * x^2) * (49 + 600 * x^2) + \n"
//                              "         1202 * y^2 - 24 * x^2 * (-201 + 198 * x^2 + 800 * x^4) * y^2 + \n"
//                              "         12 * (451 + 1204 * x^2) * y^4 + 16 * (701 + 1200 * x^2) * y^6 + \n"
//                              "         9600 * y^8) * sin(10 * y)) + x * sin(10 * x) * (-(107 + 9600 * x^8 + 602 * y^2 - \n"
//                              "          4 * y^4 * (3 + 4 * y^2) * (49 + 600 * y^2) + 16 * x^6 * (701 + 1200 * y^2) + 12 * x^4 * (451 + 1204 * y^2) + \n"
//                              "          2 * x^2 * (601 + 2412 * y^2 - 2376 * y^4 - 9600 * y^6) ) * cos(10 * y) + \n"
//                              "      10 * y * (1 + 4 * x^2 + 4 * y^2) * (69 + 1600 * x^6 + 24 * x^4 * (69 + 200 * y^2) + 24 * x^2 * (23 + 138 * y^2 + 200 * y^4) + \n"
//                              "         8 * y^2 * (69 + 207 * y^2 + 200 * y^4) ) * sin(10 * y)))",3);
////    gsFunctionExpr<> source  ("2 * cos(10 * x) * cos(10 * y)",3);
//
//    gsFunctionExpr<> laplace ("-(1 / ( (1 + 4 * x^2 + 4 * y^2)^2 ) ) * 80 * ( (1 + 2 * x^2 + \n"
//                              "             2 * y^2) * cos(10 * x) * (5 * (1 + 4 * x^2 + 4 * y^2) * cos(10 * y) - 2 * y * sin(10 * y) ) - \n"
//                              "    2 * x * sin(10 * x) * ( (1 + 2 * x^2 + 2 * y^2) * cos(10 * y) - 10 * y * (1 + 4 * x^2 + 4 * y^2) * sin(10 * y) ) )",3);
//
//    gsFunctionExpr<> solVal("2 * cos(10 * x) * cos(10 * y)",3);
//
//    gsFunctionExpr<>sol1der ("(-20 * (1 + 4 * y^2) * cos(10 * y) * sin(10 * x) + 80 * x * y * cos(10 * x) * sin(10 * y) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(80 * x * y * cos(10 * y) * sin(10 * x) - 20 * (1 + 4 * x^2) * cos(10 * x) * sin(10 * y) ) / (1 + 4 * x^2 + 4 * y^2)",
//                             "(40 * (x * cos(10 * y) * sin(10 * x) + y * cos(10 * x) * sin(10 * y)) ) / (1 + 4 * x^2 + 4 * y^2)",3);
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
            numDegree = 2; // 2 == degree 3
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
        case 13:
            string_geo = "KirchhoffLoveGeo/square_singlePatch.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 14:
            string_geo = "KirchhoffLoveGeo/planar_squareHole.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 15:
            string_geo = "KirchhoffLoveGeo/planar_untrimmedGeo.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 16:
            string_geo = "KirchhoffLoveGeo/planar_G1AScubicParam.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 17:
            string_geo = "KirchhoffLoveGeo/surfaceFrom2D.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 18:
            string_geo = "KirchhoffLoveGeo/surfaceFrom2D_cubic.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 19:
            string_geo = "KirchhoffLoveGeo/surfaceFrom2DTriangular.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 20:
            string_geo = "KirchhoffLoveGeo/surfaceFromCubic2DCentralBump.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 21:
            string_geo = "KirchhoffLoveGeo/surfaceFrom2DFourPatchesEasy.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 22:
            string_geo = "KirchhoffLoveGeo/surfaceFrom2DFourPatchesCubicQuad.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 23:
            string_geo = "KirchhoffLoveGeo/geo_fivePatch.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 24:
            string_geo = "KirchhoffLoveGeo/surfaceFrom2DFivePatchesCosCos.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 25:
            string_geo = "KirchhoffLoveGeo/geo_threePatches.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 26:
            string_geo = "KirchhoffLoveGeo/geo_bumpSurf.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 27:
            string_geo = "KirchhoffLoveGeo/fivePatch_genGeo.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 28:
            string_geo = "KirchhoffLoveGeo/geoFivePatch_twoInnerKnots.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 29:
            string_geo = "KirchhoffLoveGeo/geoFivePatch_bilinearDom.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 30:
            string_geo = "KirchhoffLoveGeo/geoFivePatch_cubicDom.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 31:
            string_geo = "KirchhoffLoveGeo/surfHalfSphere.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 32:
            string_geo = "KirchhoffLoveGeo/surfHalfSphereSquare.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 33:
            string_geo = "KirchhoffLoveGeo/surfHalfSphereSquareFivePatch.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 34:
            string_geo = "KirchhoffLoveGeo/square_sixPatch.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 35:
            string_geo = "KirchhoffLoveGeo/surfHalfSphereSixPatch.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 36:
            string_geo = "KirchhoffLoveGeo/testPortionSphere.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 37:
            string_geo = "KirchhoffLoveGeo/triangularThreePatchDisturbed.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 38:
            string_geo = "KirchhoffLoveGeo/triangularThreePatchMiddShift.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 39:
            string_geo = "KirchhoffLoveGeo/regularPentagon.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 40:
            string_geo = "KirchhoffLoveGeo/squareHoleSurf.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 41:
            string_geo = "KirchhoffLoveGeo/squareSquareHoleSurf.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 42:
            string_geo = "KirchhoffLoveGeo/squareCircHoleCircSurf.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 43:
            string_geo = "KirchhoffLoveGeo/parabolicSurfCircleHole.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 44:
            string_geo = "KirchhoffLoveGeo/quadrticQuarticSurf.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 45:
            string_geo = "KirchhoffLoveGeo/squareCircHolePolinomial.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 46:
            string_geo = "KirchhoffLoveGeo/exampleCompleteSphere.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 47:
            string_geo = "KirchhoffLoveGeo/triangularSurfNotSymRegBDY.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 48:
            string_geo = "KirchhoffLoveGeo/quadrilateralSurfNotSymRegBDY.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 49:
            string_geo = "KirchhoffLoveGeo/quadrilateralFourPatchDisturbed.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 50:
            string_geo = "KirchhoffLoveGeo/pentagonalFivePatchDisturbed.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 51:
            string_geo = "KirchhoffLoveGeo/12p_3holes_planarDom.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 52:
            string_geo = "KirchhoffLoveGeo/cilinderOneHole.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 53:
            string_geo = "KirchhoffLoveGeo/planarThreeSquareHoles.xml";
            numDegree = 2; // 2 == degree 3
            break;
        case 54:
            string_geo = "KirchhoffLoveGeo/cilinderThreeHolesSurf.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 55:
            string_geo = "KirchhoffLoveGeo/cilinderThreeHolesSurf1.xml";
            numDegree = 0; // 2 == degree 3
            break;
        case 56:
            string_geo = "KirchhoffLoveGeo/4p_hyperboloid_geom.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 57:
            string_geo = "KirchhoffLoveGeo/1p_hyperboloid_geom.xml";
            numDegree = 1; // 2 == degree 3
            break;
        case 58:
            string_geo = "KirchhoffLoveGeo/4p_hyperboloid_disturbed.xml";
            numDegree = 1; // 2 == degree 3
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

    multiPatch_init.degreeElevate(g1OptionList.getInt("degree") + g1OptionList.getInt("q_tilde"));

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

    gsVector<> numBF;
    multiPatch_init.embed(3);


    //Embedding planar geometries into R3 by adding the third component =0

    for (index_t refinement_level = 0; refinement_level < g1OptionList.getInt("loop"); refinement_level++)
    {
        gsMultiPatch<> multiPatchSurf(multiPatch_init);

//        gsInfo << "KV: " << multiPatch.patch(0).basis().basis(0) << "\n";

        multiPatchSurf.uniformRefine_withSameRegularity(num_knots[refinement_level], g1OptionList.getInt("regularity"));

        gsInfo << "KV: " << multiPatchSurf.patch(0).basis() << "\n";

        gsInfo << "###### Level: " << refinement_level << " with " << num_knots[refinement_level] << " inner knots ###### " << "\n";

        gsMultiBasis<> mb(multiPatchSurf);

#ifdef _OPENMP
        omp_set_num_threads( g1OptionList.getInt("threads"));
        omp_set_nested(1);
#endif
//        gsFileData<> interfacesGD;
//        gsFileData<> boundaryGD;
//        gsFileData<> verticesGD;


        gsG1System<real_t> g1System(multiPatchSurf, mb, g1OptionList.getSwitch("neumann"), g1OptionList.getSwitch("twoPatch"), g1OptionList);

        // ########### EDGE FUNCTIONS ###########
        // Interface loop
        for (size_t numInt = 0; numInt < multiPatchSurf.interfaces().size(); numInt++ )
        {
            const boundaryInterface & item = multiPatchSurf.interfaces()[numInt];

            std::string fileName;
            std::string basename = "InterfaceBasisFunctions" + util::to_string(numInt);
            gsParaviewCollection collection(basename);

            gsG1AuxiliaryEdgeMultiplePatches singleInt(multiPatchSurf, item.first().patch, item.second().patch);
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
                    gsField<> temp_field(multiPatchSurf.patch(item.first().patch),edgeSingleBF.patch(0));
                    gsWriteParaview(temp_field,fileName,5000);
                    collection.addTimestep(fileName,i,"0.vts");
                    // Second Interface Side
                    fileName = basename + "_1_" + util::to_string(i);
                    gsField<> temp_field_1(multiPatchSurf.patch(item.second().patch),edgeSingleBF.patch(1));
                    gsWriteParaview(temp_field_1,fileName,5000);
                    collection.addTimestep(fileName,i,"0.vts");
                }
            }
            gsInfo << "============================================================ \n";
            collection.save();
        }

        // Boundaries loop
        for (size_t numBdy = 0; numBdy < multiPatchSurf.boundaries().size(); numBdy++ )
        {
            const patchSide & bit = multiPatchSurf.boundaries()[numBdy];

            std::string fileName;
            std::string basename = "BoundaryBasisFunctions" + util::to_string(numBdy);
            gsParaviewCollection collection(basename);

            gsG1AuxiliaryEdgeMultiplePatches singleBdy(multiPatchSurf, bit.patch);
             singleBdy.computeG1BoundaryBasis(g1OptionList, bit.m_index);





            for (size_t i = 0; i < singleBdy.getSinglePatch(0).getG1Basis().nPatches(); i++)
            {
                gsMultiPatch<> edgeSingleBF;

                edgeSingleBF.addPatch(singleBdy.getSinglePatch(0).getG1Basis().patch(i));

                g1System.insertBoundaryEdge(edgeSingleBF,bit,numBdy,i);

                if (g1OptionList.getSwitch("plot"))
                {
                    fileName = basename + "_0_" + util::to_string(i);
                    gsField<> temp_field(multiPatchSurf.patch(bit.patch),edgeSingleBF.patch(0));
                    gsWriteParaview(temp_field,fileName,5000);
                    collection.addTimestep(fileName,i,"0.vts");
                }
            }
//            gsInfo << "============================================================ \n";
            collection.save();
        }

        // Vertices
        for(size_t numVer=0; numVer < multiPatchSurf.vertices().size(); numVer++)
        {
            std::string fileName;
            std::string basename = "VerticesBasisFunctions" + util::to_string(numVer);
            gsParaviewCollection collection(basename);

            std::vector<patchCorner> allcornerLists = multiPatchSurf.vertices()[numVer];
            std::vector<size_t> patchIndex;
            std::vector<size_t> vertIndex;
            for(auto & allcornerList : allcornerLists)
            {
                patchIndex.push_back(allcornerList.patch);
                vertIndex.push_back(allcornerList.m_index);
            }

                gsG1AuxiliaryVertexMultiplePatches singleVertex(multiPatchSurf, patchIndex, vertIndex);
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
                            gsField<> temp_field(multiPatchSurf.patch(patchIndex[np]), singleBasisFunction.patch(np));
                            gsWriteParaview(temp_field, fileName, 5000);
                            collection.addTimestep(fileName, i, "0.vts");
                        }
                    }
                    g1System.insertVertex(singleBasisFunction, patchIndex, numVer, singleVertex.get_internalDofs(), i);}

//            gsInfo << "============================================================ \n";
            collection.save();
        }


        gsBoundaryConditions<> bcInfo, bcInfo2;
        for (gsMultiPatch<>::const_biterator bit = multiPatchSurf.bBegin(); bit != multiPatchSurf.bEnd(); ++bit)
        {
            bcInfo.addCondition( *bit, condition_type::dirichlet, &solVal );
            if (!g1OptionList.getSwitch("neumann"))
                bcInfo2.addCondition( *bit, condition_type::laplace, &laplace);
            else
                bcInfo2.addCondition(*bit, condition_type::neumann, &sol1der );
        }

        gsInfo << "Assembling the system ... \n";
        // BiharmonicAssembler
        gsG1BiharmonicAssembler<real_t> g1BiharmonicAssembler(multiPatchSurf, mb, bcInfo, bcInfo2, source, g1OptionList);
        g1BiharmonicAssembler.assemble();

        gsInfo << "Computing Boundary data ... \n";

        if (!g1OptionList.getSwitch("neumann"))
            g1BiharmonicAssembler.computeDirichletDofsL2Proj(g1System); // Compute boundary values with laplace
        else
            g1BiharmonicAssembler.computeDirichletAndNeumannDofsL2Proj(g1System); // Compute boundary values with neumann

        g1System.finalize(multiPatchSurf,mb,g1BiharmonicAssembler.get_bValue());

        gsInfo << "Solving the system ... \n";
        gsMatrix<> solVector = g1System.solve(g1BiharmonicAssembler.matrix(), g1BiharmonicAssembler.rhs());

        gsInfo << "Plotting the solution ... \n";
        if (g1OptionList.getSwitch("plot"))
        {
            gsField<> exactField(multiPatchSurf,solVal);
            gsWriteParaview(exactField,"G1Biharmonic_exact",15000);

            // construct solution: INTERIOR
            gsMultiPatch<> mpsol;
            g1BiharmonicAssembler.constructSolution(solVector.bottomRows(g1BiharmonicAssembler.matrix().dim().first),mpsol);
            gsField<> solField(multiPatchSurf, mpsol);
            // construct solution for plotting
            std::vector<gsMultiPatch<>> g1Basis;
            g1System.constructG1Solution(solVector,g1Basis, multiPatchSurf);
            g1BiharmonicAssembler.plotParaview(solField, g1Basis);

            // Pascal
//            if(refinement_level == 1)
//            {
//                gsMultiPatch<> mp_letsee;
//                gsFileData<> fd("/home/afarahat/Desktop/gismo/filedata/KirchhoffLoveGeo/cilinderThreeHolesSurf1.xml");
//                gsInfo << "Loaded file "<< fd.lastPath() <<"\n";
//
//                gsMultiPatch<> fitMP;
//                fd.getId(0, fitMP); // id=0: Multipatch domain
//                fitMP.computeTopology();
//
//                for (size_t numP = 0; numP < multiPatchSurf.nPatches(); numP++)
//                {
//                    gsMatrix<> coefsnew(multiPatchSurf.patch(numP).coefs().dim().first, 1);
//                    coefsnew = mpsol.patch(numP).coefs();
//                    for (size_t letsee = 0; letsee < g1Basis[numP].nPatches(); letsee++)
//                    {
//                        coefsnew += g1Basis[numP].patch(letsee).coefs();
//                    }
//
//                    gsMatrix<> newcontrolpoints(multiPatchSurf.patch(numP).coefs().dim().first, 3);
//                    newcontrolpoints.leftCols(2) = fitMP.patch(numP).coefs();
//                    newcontrolpoints.col(2) = coefsnew;
//                    mp_letsee.addPatch(multiPatchSurf.patch(numP));
//                    mp_letsee.patch(numP).setCoefs(newcontrolpoints);
//                }
//                gsWriteParaview(mp_letsee, "mp_surfFitting", 15000);
//
//                gsFileData<> xml;
//                xml << mp_letsee;
//                xml.save("/home/afarahat/Desktop/gismo/filedata/KirchhoffLoveGeo/cilinderThreeHolesSurf1");
//            }
//             End
        }

        // construct solution: G1 Basis
        gsInfo << "Construct G1 solution ... \n";
        gsSparseMatrix<real_t> Sol_sparse;
        g1System.constructSparseG1Solution(solVector,Sol_sparse);

        if( g1OptionList.getSwitch("residual") == true)
        {
            if(sol_vector.empty() == true || sol_vector.size() == 1)
            {
                sol_vector.push_back(Sol_sparse);
                sol_vec_basis.push_back(mb);
                sol_vec_sys.push_back(g1System);
                gsInfo << "SolVec size " << sol_vector.size() << "\n";
            }
            else
            {

                sol_vector.at(0) = sol_vector.at(1);
                sol_vector.at(1) = Sol_sparse;

                sol_vec_basis.at(0) = sol_vec_basis.at(1);
                sol_vec_basis.at(1) = mb;

                sol_vec_sys.at(0) = sol_vec_sys.at(1);
                sol_vec_sys.at(1) = g1System;
                gsInfo << "SolVec size " << sol_vector.size() << "\n";
            }
        }

#ifdef _OPENMP
        omp_set_num_threads(g1OptionList.getInt("threads"));
        omp_set_nested(1);
#endif
        gsInfo << "Computing the error ... \n";

        if(refinement_level == 0)
            numBF = g1System.get_numBasisFunctions();

        gsSparseMatrix<real_t> Sol_sparseZero = Sol_sparse;
        Sol_sparseZero.setZero();
        gsNormL2<real_t> exactErrorL2(multiPatchSurf, Sol_sparseZero, solVal);
        exactErrorL2.compute(numBF);

        gsSparseMatrix<real_t> Sol_sparseZero1 = Sol_sparse;
        Sol_sparseZero1.setZero();
        gsSeminormH1<real_t> exactErrorH1(multiPatchSurf, Sol_sparseZero1, sol1der);
        exactErrorH1.compute(numBF);

        gsSparseMatrix<real_t> Sol_sparseZero2 = Sol_sparse;
        Sol_sparseZero2.setZero();
        gsSeminormH2<real_t> exactErrorH2(multiPatchSurf, Sol_sparseZero2, laplace);
        exactErrorH2.compute(numBF);
//
//        gsInfo << "NumBF: " << numBF << "\n";

#pragma omp parallel for
        for (index_t e = 0; e < 4; ++e)
        {
            if (e == 0)
            {
                if( g1OptionList.getSwitch("residual") == true )
                {
                    if(sol_vector.size() == 2)
                    {
                        gsG1ASResidualNormL2<real_t> errorL2(multiPatchSurf, sol_vector, sol_vec_basis);
                        errorL2.compute(sol_vec_sys);
                        l2Error_vec[refinement_level] = errorL2.value();
                        gsInfo << "L2 error vec: " << l2Error_vec << "\n";
                    }
                }
                else
                {
                    gsNormL2<real_t> errorL2(multiPatchSurf, Sol_sparse, solVal);
                    errorL2.compute(g1System.get_numBasisFunctions());

//                    l2Error_vec[refinement_level] = errorL2.value();
                    l2Error_vec[refinement_level] = errorL2.value() / exactErrorL2.value();

                }
            }

            else if (e == 1)
            {
                if( g1OptionList.getSwitch("residual") == true)
                {
                    if(sol_vector.size() == 2)
                    {
                        gsG1ASResidualSeminormH1<real_t> errorSemiH1(multiPatchSurf, sol_vector, sol_vec_basis);
                        errorSemiH1.compute(sol_vec_sys);
                        h1SemiError_vec[refinement_level] = errorSemiH1.value();

                        gsInfo << "Ref lev: " << refinement_level << "\n";
                    }

                }
                else
                {
                    gsSeminormH1<real_t> errorSemiH1(multiPatchSurf, Sol_sparse, sol1der);
                    errorSemiH1.compute(g1System.get_numBasisFunctions());

//                    h1SemiError_vec[refinement_level] = errorSemiH1.value() ;
                    h1SemiError_vec[refinement_level] = errorSemiH1.value() / exactErrorH1.value();
                }
            }
            else if (e == 2)
            {
                if( g1OptionList.getSwitch("residual") == true )
                {
                    if(sol_vector.size() == 2)
                    {
                        gsG1ASResidualSeminormH2<real_t> errorSemiH2(multiPatchSurf, sol_vector, sol_vec_basis);
                        errorSemiH2.compute(sol_vec_sys);
                        h2SemiError_vec[refinement_level] = errorSemiH2.value();
                    }

                }
                else
                {
                    gsSeminormH2<real_t> errorSemiH2(multiPatchSurf, Sol_sparse, laplace);
                    errorSemiH2.compute(g1System.get_numBasisFunctions());

//                    h2SemiError_vec[refinement_level] = errorSemiH2.value() ;
                    h2SemiError_vec[refinement_level] = errorSemiH2.value() / exactErrorH2.value();
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
//    if (g1OptionList.getInt("loop") > 1)
//    {
//        gsMatrix<> rate(g1OptionList.getInt("loop") + 1,3);
//        rate.setZero();
//        printf("|%-5s|%-14s|%-5s|%-14s|%-5s|%-14s|%-5s\n", "k","L2-error", "Rate", "H1-error",
//               "Rate", "H2-error", "Rate");
//        printf("|%-5s|%-14s|%-5s|%-14s|%-5s|%-14s|%-5s\n", "-----", "--------------", "-----", "--------------",
//               "-----", "--------------", "-----");
//        printf("|%-5d|%-14.6e|%-5.2f|%-14.6e|%-5.2f|%-14.6e|%-5.2f\n", num_knots[0], l2Error_vec[0],
//               rate(0,0),h1SemiError_vec[0], rate(0,1),h2SemiError_vec[0], rate(0,2));
//        for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
//        {
//            rate(i,0) = log2(l2Error_vec[i-1] / l2Error_vec[i]);
//            rate(i,1) = log2(h1SemiError_vec[i-1] / h1SemiError_vec[i]);
//            rate(i,2) = log2(h2SemiError_vec[i-1] / h2SemiError_vec[i]);
//            printf("|%-5d|%-14.6e|%-5.2f|%-14.6e|%-5.2f|%-14.6e|%-5.2f\n", num_knots[i], l2Error_vec[i],
//                   rate(i,0),h1SemiError_vec[i], rate(i,1),h2SemiError_vec[i], rate(i,2));
//        }
//        if (g1OptionList.getSwitch("latex"))
//        {
//            printf("%-5d & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\ \n", num_knots[0],
//                   l2Error_vec[0], rate(0,0),h1SemiError_vec[0], rate(0,1), h2SemiError_vec[0], rate(0,2));
//            for (index_t i = 1; i < g1OptionList.getInt("loop"); i++)
//            {
//                rate(i,0) = log2(l2Error_vec[i-1] / l2Error_vec[i]);
//                rate(i,1) = log2(h1SemiError_vec[i-1] / h1SemiError_vec[i]);
//                rate(i,2) = log2(h2SemiError_vec[i-1] / h2SemiError_vec[i]);
//                printf("%-5d & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\ \n", num_knots[i],
//                       l2Error_vec[i], rate(i,0),h1SemiError_vec[i], rate(i,1),h2SemiError_vec[i], rate(i,2));
//            }
//        }

    const char* var1 = "|%-5d|%-14.6e|%-5.2f|%-14.6e|%-5.2f|%-14.6e|%-5.2f\n";
    const char* var2 = "%-5d & %-14.6e & %-5.2f & %-14.6e & %-5.2f & %-14.6e & %-5.2f \\\\ \n";

    if (std::string (typeid(real_t).name())== "e") // long double
    {
        var1 = "|%-5d|%-14.6Le|%-5.2Lf|%-14.6Le|%-5.2Lf|%-14.6Le|%-5.2Lf\n";
        var2 = "%-5d & %-14.6Le & %-5.2Lf & %-14.6Le & %-5.2Lf & %-14.6Le & %-5.2Lf \\\\ \n";
        gsInfo << "Long double \n";
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


    }
    else
    {
        gsInfo << "L2 Error: " << l2Error_vec[0] << "\n";
        gsInfo << "H1 Error: " << h1SemiError_vec[0] << "\n";
        gsInfo << "H2 Error: " << h2SemiError_vec[0] << "\n";
    }
    gsInfo << "=====================================================================\n";


} // main