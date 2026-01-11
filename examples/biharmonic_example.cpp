/** @file gsMakeMultiPatch.cpp

    @brief Computes correctly the boundaries and interfaces of
    a multipatch structure

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#include <gismo.h>

#include <gsIO/gsIOUtils.h>

//#include "/home/turing/theydarov/gismo/extensions/motor/jku/gsQualityMeasure2.h"

//#include "/home/turing/theydarov/gismo/extensions/motor/jku/gsMotorUtils.h"

#define PI 3.14159265358979323846
using namespace gismo;
using namespace std;
#include <fstream>
ofstream outfile;

void printTheMatrix(const gsMatrix<real_t>& matrix, const std::string& matrixName) {
    outfile << matrixName << "\n";
    for (size_t row = 0; row < matrix.rows(); ++row) {
        for (size_t col = 0; col < matrix.cols(); ++col) {
            outfile << matrix(row, col) << " ";
        }
        outfile << "\n";
    }
}

int PatchesIntersection(gsGeometry < >& geom1, gsGeometry < >& geom2,
    //                        gsVector<real_t> referencePoints, gsVector<index_t> & isToBeMoved,
    //                         gsVector<index_t>  sortedBoundaryFunctions,
    index_t numPoints = 100, real_t tolerance = 1e-6) {
    int doPatchesIntersect = 0;
    boxSide const s;
    //    gsVector<index_t> boundaryFunctions = geom1->basis().boundary(s);
    gsVector<index_t> boundaryFunctions = geom1.basis().boundary(s);
    //    gsVector<index_t> boundaryFunctions2 = geom2->basis().boundary(s);

    //    //gsInfo << boundaryFunctions;
    //    return  0;
    gsVector<index_t>  sortedBoundaryFunctions(boundaryFunctions.size());
    //    gsVector<>  auxiliaryVector(boundaryFunctions.size());
    gsVector<real_t>  referencePoints(boundaryFunctions.size());
    sortedBoundaryFunctions = boundaryFunctions;
    for (int i = 0; i < boundaryFunctions.size(); ++i) {
        //        auxiliaryVector( i) = geom1->coef(boundaryFunctions(i))(1);
        referencePoints(i) = geom1.coef(boundaryFunctions(i))(1);
    }
    for (int i = 0; i < referencePoints.size(); ++i) {
        //        sortedBoundaryFunctions(i) = auxiliaryVector.maxCoeff();
        for (int j = 0; j < referencePoints.size(); ++j) {
            if (referencePoints(i) < referencePoints(j)) {
                index_t a = sortedBoundaryFunctions(i);
                real_t  b = referencePoints(i);
                sortedBoundaryFunctions(i) = sortedBoundaryFunctions(j);
                referencePoints(i) = referencePoints(j);
                sortedBoundaryFunctions(j) = a;
                referencePoints(j) = b;
            }
        }
    }
    gsVector<index_t> isToBeMoved(sortedBoundaryFunctions.size());
    gsMatrix<> xy1;
    gsMatrix<> xy2;
    //    gsMatrix<> uv1(2, (numPoints + 1)*(referencePoints.size() - 1));
    gsMatrix<> uv1(2, (numPoints + 1) * (1));
    //    gsMatrix<> uv2(2, (numPoints + 1)*(referencePoints.size() - 1)););
    gsMatrix<> uv2(2, (numPoints + 1) * (1));
    //    gsMatrix<> pointsOfIntersection;
    for (int i = 0; i < isToBeMoved.size(); ++i) {
        isToBeMoved(i) = 0;
    }
    bool thirdIsAbove, fourthIsAbove;
    for (int i = 0; i < referencePoints.size() - 1; ++i) {
        //    for (int i = 0; i < 1; ++i) {
        for (int j = 0; j <= numPoints; ++j) {
            const real_t val = referencePoints(i) + (referencePoints(i + 1) - referencePoints(i)) * (double)j / numPoints;
            const real_t zeroWTF = 0.0;
            const real_t  oneWTF = 1.0;
            uv1(1, j) = uv2(1, j) = std::min(std::max(val, zeroWTF), oneWTF);
            //            //gsInfo << referencePoints(i) << " + (" << referencePoints(i + 1) << " - " << referencePoints(i)
            //                   << ") *" << j << "/" << numPoints << "=" << uv1(1,j) << "\n";
            uv1(0, j) = 1.0;
            uv2(0, j) = 0.0;
        }

        //        //gsInfo << "uv1 " << uv1 << "\n";
        //        //gsInfo << "uv1 " << uv1 << "\n";
        //        //gsInfo << uv2 << "\n";
        xy1 = geom1.eval(uv1);
        xy2 = geom2.eval(uv2);
        //        //gsDebugVar(xy1);
        //        //gsDebugVar(xy2);

        for (int j = 0; j < numPoints; ++j) {
            real_t x1 = xy1(0, j);
            real_t y1 = xy1(1, j);
            real_t x2 = xy1(0, j + 1);
            real_t y2 = xy1(1, j + 1);
            real_t x3 = xy2(0, j);
            real_t y3 = xy2(1, j);
            real_t x4 = xy2(0, j + 1);
            real_t y4 = xy2(1, j + 1);
            real_t val1 = (y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1);
            real_t val2 = (y4 - y1) * (x2 - x1) - (y2 - y1) * (x4 - x1);
            thirdIsAbove = val1 > -tolerance;
            fourthIsAbove = val2 > -tolerance;

            //            //gsInfo << "i = "<< i << " " <<  thirdIsAbove <<" "<< fourthIsAbove << "\n"
            //                   << x1 << " " << y1 << "\n"
            //                   << x2 << " " << y2 << "\n"
            //                   << x3 << " " << y3 << "\n"
            //                   << x4 << " " << y4 << "\n"
            //                   //            ;
            //                   << (y3 - y1)*(x2 - x1) << " " << (y2 - y1)*(x3 - x1) << " " << (y3 - y1)*(x2 - x1) - (y2 - y1)*(x3 - x1) << "\n"
            //                   << (y4 - y1)*(x2 - x1) << " " << (y2 - y1)*(x4 - x1) << " " << (y4 - y1)*(x2 - x1) - (y2 - y1)*(x4 - x1) << "\n"
            //                    ;
            if (!((thirdIsAbove || !fourthIsAbove) && (!thirdIsAbove || fourthIsAbove))) {
                //gsInfo << "i = " << i << " " << thirdIsAbove << " " << fourthIsAbove << "\n"
                //<< x1 << " " << y1 << "\n"
                //<< x2 << " " << y2 << "\n"
                //<< x3 << " " << y3 << "\n"
                //<< x4 << " " << y4 << "\n"
                ////            ;
                //<< (y3 - y1) * (x2 - x1) << " " << (y2 - y1) * (x3 - x1) << " " << (y3 - y1) * (x2 - x1) - (y2 - y1) * (x3 - x1) << "\n"
                //<< (y4 - y1) * (x2 - x1) << " " << (y2 - y1) * (x4 - x1) << " " << (y4 - y1) * (x2 - x1) - (y2 - y1) * (x4 - x1) << "\n"
                //;

                isToBeMoved(i) = isToBeMoved(i + 1) = 1;
                doPatchesIntersect = 1;

            }

        }

    }
    return doPatchesIntersect;
}


int main(int argc, char* argv[])
{
    //    std::string filename("planar/lshape2d_3patches_tens.xml");
    //std::string filename("planar/two_squares.xml");
    //std::string filename("planar/two_squares_THB.xml");
    //    std::string filename("planar/square_and_puzzle.xml");
    //    std::string filename("/home/turing/theydarov/gismoMP/gismo/filedata/domain2d/yeti_mp2.xml");
    //    std::string filename("/home/turing/theydarov/gismoMP/gismo/filedata/domain2d/yeti_mp2THB.xml");
    //std::string filename("hexagon_3p_3l_wiggly_interface.xml");
    //std::string filename("ellipse_hole.xml");
    //std::string filename("joystick.xml");
    //std::string filename("hexagonTest1.xml");
    //std::string filename("triangle2d_u3p.xml");
    //std::string filename("hexagon_3p_2l.xml");
    //std::string filename("hexagon_3p_wigglyint_01.xml");
    //std::string filename("hexagon_3p_2l_wigglyint_01_12.xml");
    //std::string filename("hexagon_3p_2l_wigglyint_01_12_test.xml");
    //std::string filename("hexagon_3p_2l_coons.xml");
    //std::string filename("hexagon_3p_2l_test.xml");
    //std::string filename("hexagon_3p_4l_coons1.xml");
    std::string filename("hexagon_3p_4l_coons1_test.xml");
    //std::string filename("hexagonTest1.xml");
    //std::string filename("tv_190325.xml");
    //std::string filename("hexagonBeispiel.xml");
    //std::string filename("tv.xml");
    std::string fileLoc = "biharmonic_example";
    outfile.open(fileLoc + ".txt");
    index_t degree;
    //    std::string filename("/home/turing/theydarov/geometries/footTHB.xml");
    //    std::string filename("/home/turing/theydarov/gismo/filedata/domain2d/lakeTHB.xml");
    //    std::string left("/home/turing/theydarov/geometries/jigsawpuzzleTP.xml");
    real_t tol = 1e-8;
    real_t gtol = 1e-8;
    bool reparam = false, gaps = true;
    index_t method = 0, nknots = 3;

    gsCmdLine cmd("Computes the topology of a set of patches, identifing interfaces and boundaries.");
    cmd.addPlainString("filename", "File containing multipatch input (.xml).", filename);
    cmd.addReal("t", "tolerance", "Tolerance for identifing patch interfaces", tol);
    cmd.addReal("g", "gap-tolerance", "Tolerance for closing gaps", gtol);
    cmd.addInt("d", "degree", "Degree of B-splines for reparameterization", degree);
    cmd.addInt("k", "knots", "Number of interior knots for reparameterization", nknots);
    cmd.addSwitch("reparam", "Reparameterize all patches using a fixed degree and number of knots", reparam);
    cmd.addSwitch("nogaps", "Close any gaps along interfaces upto tolerance \'gtol\'", gaps);

    cmd.addInt("m", "method", "method", method);

    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }

    //-------300725
    gsKnotVector <> ku(0, 1, 0, 3);
    gsTensorBSplineBasis<2, real_t> tens(ku, ku);
    gsTHBSplineBasis<2, real_t> thb1(tens);
    std::vector<index_t> box;
    box.push_back(3);
    //box.push_back(7);
    box.push_back(0);
    box.push_back(7);
    box.push_back(8);
    box.push_back(8);
    thb1.refineElements(box);
    //hb1.refineElements(box);
    //gsInfo << "HB size: " << hb1.size() << "\n";
    gsInfo << "THB size: " << thb1.size() << "\n";
    box.clear();
    box.push_back(4);
    box.push_back(0);
    box.push_back(0);
    box.push_back(16);
    box.push_back(14);
    //        box.push_back(16);
    thb1.refineElements(box);
    //hb1.refineElements(box);
    //gsInfo << "HB size: " << hb1.size() << "\n";
    gsInfo << "THB size: " << thb1.size() << "\n";
    //box.clear();
    //box.push_back(4);
    //box.push_back(0);
    //box.push_back(14);
    //box.push_back(14);
    //box.push_back(16);
    //gsInfo << "Logging started\n";
    //thb1.refineElements(box);
    //hb1.refineElements(box);
    //gsInfo << "HB size: " << hb1.size() << "\n";
    gsInfo << "THB size: " << thb1.size() << "\n";
    box.clear();
    //return 250925;
    //-------300725

    

    gsMultiPatch<>::uPtr mp = gsReadFile<>(filename);
    gsVector<> lowc1(2);
    gsVector<> uppc1(2);
    //              gsInfo << bb->support()(0,0) << "\n" << bb->support()(0,1) << "\n\n";

    lowc1(0) = 0;//bb->support()(0,0);
    lowc1(1) = 0;//bb->support()(1, 0);
    uppc1(0) = 1;//bb->support()(0,1);
    uppc1(1) = 1;//bb->support()(1,1);
//        gsDebugVar(bb->support());
    const gsVector<real_t> lowc2 = lowc1;
    gsDebugVar(lowc2);
    const gsVector<real_t> uppc2 = uppc1;
    gsDebugVar(uppc2);
    auto numpoints = 101;
    gsMatrix<> uv2 = uniformPointGrid(lowc2, uppc2, numpoints * numpoints);
    for (size_t i = 0; i < mp->nPatches(); i++)
    {
        gsMatrix<> xy2 = mp->patch(i).eval(uv2);
        for (size_t j = 0; j < uv2.cols(); j++)
        {
            auto jacdet = mp->patch(i).jacobian(uv2.col(j)).determinant();
            if (jacdet <= 0)
            {
                gsInfo << "Problem in input! i = " << i << ",j = " << j << ",\n point\n" << uv2.col(j) << ",\n jacobian = " << jacdet << "\n";
            }
        }
    }
    //    gsGeometry<>::uPtr geom = gsReadFile<>(left);

    //SOME PLAY WITH mp
    /*{
        gsVector<> lowc1(2);
        gsVector<> uppc1(2);
//        //        gsInfo << bb->support()(0,0) << "\n" << bb->support()(0,1) << "\n\n";
//
        lowc1(0) = 0;
        lowc1(1) = 0;
        uppc1(0) = 1;
        uppc1(1) = 1;
//        gsDebugVar(bb->support());
        const gsVector<real_t> lowc2 = lowc1;
//        gsDebugVar(lowc2);
        const gsVector<real_t> uppc2 = uppc1;
//        gsDebugVar(uppc2);
        gsMatrix<> uv1 = uniformPointGrid(lowc2, uppc2, 10 * 10);
        //    gsDebugVar(uv1);

        for (int patch = 0; patch < mp->nPatches(); ++patch) {
            gsDebugVar(patch);
            gsMatrix<> xy1(1, uv1.cols());
            for (int i = 0; i < uv1.cols(); ++i) {
                gsMatrix<> punto = uv1.col(i);
                //gsDebugVar(mp->patch(0).jacobian(uv1.col(i)).determinant());
                xy1(0, i) = mp->patch(patch).jacobian(uv1.col(i)).determinant();
            }
            gsDebugVar(xy1);

        }
    }*/













    ///////////////////////////////////////////////////////////
    //    gsInfo << geom->coefs()(0,0) << "Rows, " << geom->coefs()(0,1) << "\n";
    /* gsMatrix<> coefsR = mp->patch(0).coefs();
     gsDebugVar(coefsR);
     gsMatrix<> coefsL(coefsR.rows(),coefsR.cols());
     for (int i = 1; i < mp->patch(0).coefs().rows() - 1; ++i) {
         coefsL(i,0) = coefsR((i/34) * 34 + 33 - i % 34,0);
         coefsL(i,1) = coefsR(i,1);
     }
     gsTHBSplineBasis<2> * bb =
             dynamic_cast<gsTHBSplineBasis<2> *>(&mp->patch(0).basis().source());
     gsVector<> lowc1(2);
     gsVector<> uppc1(2);
     //        gsInfo << bb->support()(0,0) << "\n" << bb->support()(0,1) << "\n\n";
     real_t au = bb->support()(0,0);
     real_t bu = bb->support()(1,1);
     index_t interior = 0;//just to be uniform
     degree = bb->degree(0);
     index_t multEnd = degree + 1;
     gsKnotVector <> ku(au, bu, interior, multEnd);
     gsTensorBSplineBasis<2, real_t> tens(ku, ku);
     gsTHBSplineBasis<2, real_t> THB(tens);
     std::vector<index_t> box;
     box.push_back(2);
     box.push_back(0);
     box.push_back(0);
     box.push_back(4);
     box.push_back(4);
     THB.refineElements(box);
     box.clear();
     gsGeometry<>::uPtr geomR = THB.makeGeometry(coefsL);
     gsWrite(*geomR, "puzzleRIght");
     return 0;}*/
     ////    index_t lastMinInd = -1;
     ////    gsInfo << "MAX" << geom->coefs().col(0).maxCoeff() << "\n";
     ////    for (int i = 1; i < geom->coefs().rows() - 1; ++i) {
     //////        gsInfo << i << "\n";
     ////        if(geom->coefs()(i,0) > geom->coefs()(i - 1,0) && (geom->coefs()(i,0) - geom->coefs()(i + 1,0)) > geom->coefs().col(0).maxCoeff()/4.0){// for now
     //////            gsInfo << i << " is the locmax\n";
     ////            for (int j = lastMinInd + 1; j < i; ++j) {
     ////                coefsR(j,0) = 2*geom->coefs()(i,0) - geom->coefs()(j,0);
     ////            }
     ////            coefsR(i,0) = geom->coefs()(i,0);
     ////            lastMinInd = i;
     ////        }
     ////    }
     ////    if(geom->coefs()(geom->coefs().rows() - 2,0) < geom->coefs()(geom->coefs().rows() - 1,0) ){
     ////        for (int j = lastMinInd; j < geom->coefs().rows(); ++j) {
     ////            coefsR(j,0) = 2*geom->coefs()(geom->coefs().rows() - 1,0) - geom->coefs()(j,0);
     ////            gsInfo << geom->coefs()(geom->coefs().rows() - 1,0) << ",\t" << geom->coefs()(j,0) << "\n";
     ////        }
     ////    }
     ////    coefsR(geom->coefs().rows() - 1,0) = geom->coefs()(geom->coefs().rows() - 1,0);
     ////    gsInfo <<" Got"<< *mp <<" \n" ;
     //    gsMatrix<index_t> lowCorners;
     //    gsMatrix<index_t> upCorners;
     //    gsVector<index_t> Level;
     //    size_t d = mp->parDim();
    gsMultiPatch<> newmp;
    gsMultiPatch<> mp1;
    gsMatrix<> punto(2, 1);
    punto(0, 0) = punto(1, 0) = 0.90141560817564836;
    auto mateval = thb1.function(0).eval(punto);
    for (size_t i = 0; i < mp->nPatches(); ++i)
    {
        gsTensorBSplineBasis<2>* bb =
            dynamic_cast<gsTensorBSplineBasis<2> *>(&mp->patch(i).basis().source());
        //gsMatrix<> coefsForTHB = mp->patch(i).coefs();
        //gsWrite(coefsForTHB, "coefsForTHB");
        /*    gsTHBSplineBasis<2>* bb =
                dynamic_cast<gsTHBSplineBasis<2> *>(&mp->patch(i).basis().source());*/

        auto numpoints = 101;
        gsMatrix<> uv1 = uniformPointGrid(lowc2, uppc2, numpoints * numpoints);
        gsMatrix<> xy1 = mp->patch(i).eval(uv1);
        for (size_t j = 0; j < uv1.cols(); j++)
        {
            if (auto jacdet = mp->patch(i).jacobian(uv1.col(j)).determinant() <= 0)
            {
                gsInfo << "Problem! i = " << i << ",j = " << j << ", jacobian = " << jacdet << "\n";
            }
        }
        auto tangentSlope1 = -0.5 / 0.8660254037844386;
        //if (i == 0)
        //{
        //    for (size_t i = 0; i < xy1.cols(); i++)
        //    {

        //        xy1(0, i) *= 0.8660254037844386;// std::sin(tangentSlope1);
        //        xy1(1, i) += tangentSlope1 * xy1(0, i);// std::sin(tangentSlope1);
        //    }
        //}
        real_t au = 0;//bb->support()(0,0);
        real_t bu = 1;//bb->support()(1,1);
        index_t interior = 0;//just to be uniform
        index_t degree = 2;// bb->degree(0);
        index_t multEnd = degree + 1;
        gsKnotVector <> ku(au, bu, interior, multEnd);
        gsTensorBSplineBasis<2, real_t> tens(ku, ku);
        gsTHBSplineBasis<2, real_t> THB(tens);
        //gsTHBSplineBasis<2> tt(*bb);
        //(tt.tree()).getBoxes(lowCorners, upCorners, Level);
        //(bb->tree()).getBoxes(lowCorners, upCorners, Level);
        /*std::vector<index_t> box;
        box.push_back(4);
        box.push_back(0);
        box.push_back(0);
        box.push_back(16);
        box.push_back(16);
        THB.refineElements(box);
        box.clear();*/
        if (i == 0) {
            std::vector<index_t> box;
            box.push_back(3);
            box.push_back(7);
            box.push_back(7);
            box.push_back(8);
            box.push_back(8);
            THB.refineElements(box);
            box.clear();
            box.push_back(4);
            box.push_back(0);
            box.push_back(0);
            box.push_back(16);
            box.push_back(14);
            THB.refineElements(box);
            box.clear();
            box.push_back(4);
            box.push_back(0);
            box.push_back(14);
            box.push_back(14);
            box.push_back(16);
            THB.refineElements(box);
            box.clear();
        }
        if (i > 0) {
            std::vector<index_t> box;
            box.push_back(4);
            box.push_back(0);
            box.push_back(0);
            box.push_back(16);
            box.push_back(16);
            THB.refineElements(box);
            box.clear();
        }
        /* if (i < 1) {
         std::vector<index_t> box;
         box.push_back(3);
         box.push_back(4);
         box.push_back(0);
         box.push_back(8);
         box.push_back(4);
         THB.refineElements(box);
         box.clear();
         }*/


        gsFitting<> fit(uv1, xy1, THB);
        fit.compute();
        fit.computeErrors();
        gsInfo << fit.getBasis() << "\n";
        gsSparseMatrix<real_t>matA(THB.size(), THB.size());
        gsMatrix<real_t>matAfull(THB.size(), THB.size());
        gsMatrix<real_t>vectB(THB.size(), 2);
        matA.setZero();
        matAfull.setZero();
        vectB.setZero();
        fit.assembleSystem(matA, vectB);
        matAfull = gsEigen::MatrixXd(matA);
        printTheMatrix(matA, "matA");
        gsFileManager::open("biharmonic_example.txt");
        return 251203;
        gsMatrix<real_t>vectSol(THB.size(), 2);
        vectSol = matAfull.partialPivLu().solve(vectB);
        auto epsilon = 3e-2;// 1e-4;// 6;
        auto epsilonInterface = 1e-4;// 1e-4;// 6;
        auto oscillationFrequency = 40;
        auto oscillationFrequencyV = 60;// 40 * std::sqrt(3);
        auto pointShift = 0.00;
        auto tangentSlope = (vectSol(1, 1) - vectSol(0, 1)) / (vectSol(1, 0) - vectSol(0, 0));
        gsInfo << tangentSlope << "\n";
        auto shiftFactor = 0.00 * tangentSlope;//s0.06;
        auto innershiftFactor = 0.00 * tangentSlope;//s0.06;
        int boundarySize = (int)std::sqrt(vectSol.rows());
        auto amplitudeCoefficient = 0.05;
        //Create a wiggly boundary
        //raw,1d solution
        gsInfo << "boundarySize = " << boundarySize << "\n";
        if (i == 0)
        {
            for (int j = 1; j < boundarySize - 1; ++j)
            {
                //gsInfo << "j = " << j << "\n";
                ////EAST
                //vectSol((j + 1) * boundarySize - 1, 0) += pow(-1, j + 1) * shiftFactor;
                ////NORTH
                //vectSol((boundarySize - 1) * boundarySize + j, 0) += pow(-1, j + 1) * shiftFactor;
                ////WEST
                //vectSol(j * boundarySize, 0) += pow(-1, j + 1) * shiftFactor;
                //for (int k = 1; k < boundarySize - 1; ++k)
                //{
                //    // Extract values from vectSol with explicit logging
                //    double v1 = (double)(j) / (boundarySize - 1); //vectSol(j * boundarySize + k, 0);
                //    double u1 = (double)(k) / (boundarySize - 1); // vectSol(k * boundarySize + j, 0);
                //    double v2 = (double)(j) / (boundarySize - 1); //vectSol(j * boundarySize + k, 1);
                //    double u2 = (double)(k) / (boundarySize - 1); // vectSol(k * boundarySize + j, 1);

                //    gsInfo << "[Iteration " << k << "] j = " << j << ", k = " << k << "\n";
                //    gsInfo << "[Updating vectSol(" <<  j * boundarySize + k << ")\n";
                //    //gsInfo << "  Extracting values from vectSol:\n";
                //    //gsInfo << "    v1 (vectSol(" << j * boundarySize + k << ", 0)) = " << v1 << "\n";
                //    //gsInfo << "    u1 (vectSol(" << k * boundarySize + j << ", 0)) = " << u1 << "\n";
                //    //gsInfo << "    v2 (vectSol(" << j * boundarySize + k << ", 1)) = " << v2 << "\n";
                //    //gsInfo << "    u2 (vectSol(" << k * boundarySize + j << ", 1)) = " << u2 << "\n";

                //    // Compute interpolated values for x-coordinates
                //    auto S11 = (1 - v1) * vectSol(k, 0) + v1 * vectSol(boundarySize * (boundarySize - 1) + k, 0);
                //    auto S21 = (1 - u1) * vectSol(j * boundarySize, 0) + u1 * vectSol((j + 1) * boundarySize - 1, 0);
                //    auto S31 = (1 - u1)* (1 - v1)* vectSol(0, 0) + u1 * (1 - v1) * vectSol(boundarySize * (boundarySize - 1), 0) +
                //        v1 * (1 - u1) * vectSol(boundarySize - 1, 0) + u1 * v1 * vectSol(boundarySize * boundarySize - 1, 0);

                //    /*gsInfo << "  Computing x-coordinates:\n";
                //    gsInfo << "    S11 = (1 - " << v1 << ") * vectSol(" << k << ", 0) + " << v1 << " * vectSol(" << boundarySize * (boundarySize - 1) + k << ", 0) = " << (1 - v1) << "*" << vectSol(k, 0) << "+" << v1 << "*" << vectSol(boundarySize * (boundarySize - 1) + k, 0) << " = " << S11 << "\n";
                //    gsInfo << "    S21 = (1 - " << u1 << ") * vectSol(" << j * boundarySize << ", 0) + " << u1 << " * vectSol(" << (j + 1) * boundarySize - 1 << ", 0) = " << S21 << "\n";
                //    gsInfo << "    S31 = (1 - " << u1 << ") * (1 - " << v1 << ") * vectSol(0, 0) + ... = " << S31 << "\n";*/

                //    // Compute interpolated values for y-coordinates
                //    auto S12 = (1 - v2) * vectSol(k, 0) + v2 * vectSol(boundarySize * (boundarySize - 1) + k, 0);
                //    auto S22 = (1 - u2) * vectSol(j * boundarySize, 0) + u2 * vectSol((j + 1) * boundarySize - 1, 0);
                //    auto S32 = 0;// (1 - u2)* (1 - v2)* vectSol(0, 1) + u2 * (1 - v2) * vectSol(boundarySize * (boundarySize - 1), 1) +
                //        //v2 * (1 - u2) * vectSol(boundarySize - 1, 1) + u2 * v2 * vectSol(boundarySize * boundarySize - 1, 1);

                //    /*gsInfo << "  Computing y-coordinates:\n";
                //    gsInfo << "    S12 = (1 - " << v2 << ") * vectSol(" << k << ", 1) + " << v2 << " * vectSol(" << boundarySize * (boundarySize - 1) + k << ", 1) = " << S12 << "\n";
                //    gsInfo << "    S22 = (1 - " << u2 << ") * vectSol(" << j * boundarySize << ", 1) + " << u2 << " * vectSol(" << (j + 1) * boundarySize - 1 << ", 1) = " << S22 << "\n";
                //    gsInfo << "    S32 = (1 - " << u2 << ") * (1 - " << v2 << ") * vectSol(0, 1) + ... = " << S32 << "\n";*/

                //    gsInfo << "  Computing x-coordinates:\n";
                //    gsInfo << "    S11 = (1 - " << v1 << ") * vectSol(" << k << ", 0) + " << v1
                //        << " * vectSol(" << boundarySize * (boundarySize - 1) + k << ", 0) = "
                //        << (1 - v1) << " * " << vectSol(k, 0) << " + " << v1 << " * "
                //        << vectSol(boundarySize * (boundarySize - 1) + k, 0) << " = " << S11 << "\n";

                //    gsInfo << "    S21 = (1 - " << u1 << ") * vectSol(" << j * boundarySize << ", 0) + "
                //        << u1 << " * vectSol(" << (j + 1) * boundarySize - 1 << ", 0) = "
                //        << (1 - u1) << " * " << vectSol(j * boundarySize, 0) << " + " << u1
                //        << " * " << vectSol((j + 1) * boundarySize - 1, 0) << " = " << S21 << "\n";

                //    gsInfo << "  Computing y-coordinates:\n";
                //    gsInfo << "    S12 = (1 - " << v2 << ") * vectSol(" << k << ", 1) + " << v2
                //        << " * vectSol(" << boundarySize * (boundarySize - 1) + k << ", 1) = "
                //        << (1 - v2) << " * " << vectSol(k, 1) << " + " << v2 << " * "
                //        << vectSol(boundarySize * (boundarySize - 1) + k, 1) << " = " << S12 << "\n";

                //    gsInfo << "    S22 = (1 - " << u2 << ") * vectSol(" << j * boundarySize << ", 1) + "
                //        << u2 << " * vectSol(" << (j + 1) * boundarySize - 1 << ", 1) = "
                //        << (1 - u2) << " * " << vectSol(j * boundarySize, 1) << " + " << u2
                //        << " * " << vectSol((j + 1) * boundarySize - 1, 1) << " = " << S22 << "\n";

                //    // Compute final values
                //    auto newX = S11 + S21 - S31;
                //    auto newY = S12 + S22 - S32;

                //    gsInfo << "  Final computed values:\n";
                //    gsInfo << "    newX = " << S11 << " + " << S21 << " - " << S31 << " = " << newX << "\n";
                //    gsInfo << "    newY = " << S12 << " + " << S22 << " - " << S32 << " = " << newY << "\n";

                //    gsInfo << "  Old vectSol(" << k * boundarySize + j << ", 0) = " << vectSol(k * boundarySize + j, 0) << "\n";
                //    gsInfo << "  Old vectSol(" << k * boundarySize + j << ", 1) = " << vectSol(k * boundarySize + j, 1) << "\n";

                //    // Apply changes
                //    vectSol(k * boundarySize + j, 0) = newX;
                //    vectSol(k * boundarySize + j, 1) = newY;

                //    gsInfo << "  Updated vectSol(" << k * boundarySize + j << ", 0) = " << vectSol(k * boundarySize + j, 0) << "\n";
                //    gsInfo << "  Updated vectSol(" << k * boundarySize + j << ", 1) = " << vectSol(k * boundarySize + j, 1) << "\n";
                //    gsInfo << "-------------------------------------------\n";
                //}


            }

            //            vectSol(7, 0) += innershiftFactor;
            //            vectSol(13, 0) -= innershiftFactor;
            //            vectSol(19, 0) += innershiftFactor;
            //            vectSol(25, 0) -= innershiftFactor;
            //
            //            vectSol(8, 0) += innershiftFactor;
            //            vectSol(14, 0) -= innershiftFactor;
            //            vectSol(20, 0) += innershiftFactor;
            //            vectSol(26, 0) -= innershiftFactor;
            //
            //            vectSol(9, 0) += innershiftFactor;
            //            vectSol(15, 0) -= innershiftFactor;
            //            vectSol(21, 0) += innershiftFactor;
            //            vectSol(27, 0) -= innershiftFactor;
            //
            //            vectSol(10, 0) += innershiftFactor;
            //            vectSol(16, 0) -= innershiftFactor;
            //            vectSol(22, 0) += innershiftFactor;
            //            vectSol(28, 0) -= innershiftFactor;
            //
            /*vectSol(7, 1) += innershiftFactor;
            vectSol(8, 1) -= innershiftFactor;
            vectSol(9, 1) += innershiftFactor;
            vectSol(10, 1) -= innershiftFactor;

            vectSol(13, 1) += innershiftFactor;
            vectSol(14, 1) -= innershiftFactor;
            vectSol(15, 1) += innershiftFactor;
            vectSol(16, 1) -= innershiftFactor;

            vectSol(19, 1) += innershiftFactor;
            vectSol(20, 1) -= innershiftFactor;
            vectSol(21, 1) += innershiftFactor;
            vectSol(22, 1) -= innershiftFactor;

            vectSol(25, 1) += innershiftFactor;
            vectSol(26, 1) -= innershiftFactor;
            vectSol(27, 1) += innershiftFactor;
            vectSol(28, 1) -= innershiftFactor;*/


            /*         vectSol(1, 0) += shiftFactor;
                     vectSol(2, 0) -= shiftFactor;
                     vectSol(3, 0) += shiftFactor;
                     vectSol(4, 0) -= shiftFactor;*/
        }
        if (i == 1)
        {
            //gsInfo << "alteration, I =" << i << "\n";
            vectSol(1, 1) -= shiftFactor;
            vectSol(2, 1) += shiftFactor;
            vectSol(3, 1) -= shiftFactor;
            vectSol(4, 1) += shiftFactor;

            vectSol(6, 0) += shiftFactor;
            vectSol(12, 0) -= shiftFactor;
            vectSol(18, 0) += shiftFactor;
            vectSol(24, 0) -= shiftFactor;
            for (int j = 1; j < boundarySize - 1; ++j)
            {
                //gsInfo << "j = " << j << "\n";
                ////EAST
                //vectSol((j + 1) * boundarySize - 1, 0) += pow(-1, j + 1) * shiftFactor;
                ////NORTH
                //vectSol((boundarySize - 1) * boundarySize + j, 0) += pow(-1, j + 1) * shiftFactor;
                ////WEST
                //vectSol(j * boundarySize, 0) += pow(-1, j + 1) * shiftFactor;
                for (int k = 1; k < boundarySize - 1; ++k)
                {
                    vectSol(k * boundarySize + j, 0) += pow(-1, k) * innershiftFactor;
                    //gsInfo << k << " " << " " << j << " " << k * boundarySize + j << "\n";
                    vectSol(j * boundarySize + k, 1) += pow(-1, k) * innershiftFactor;
                    //gsInfo << j << " "  << k << " " << j * boundarySize + k << "\n";
                }
            }
        }
        if (i == 2)
        {
            vectSol(6, 0) += shiftFactor;
            vectSol(12, 0) -= shiftFactor;
            vectSol(18, 0) += shiftFactor;
            vectSol(24, 0) -= shiftFactor;

            vectSol(1, 0) += shiftFactor;
            vectSol(2, 0) -= shiftFactor;
            vectSol(3, 0) += shiftFactor;
            vectSol(4, 0) -= shiftFactor;

            for (int j = 1; j < boundarySize - 1; ++j)
            {
                //gsInfo << "j = " << j << "\n";
                ////EAST
                //vectSol((j + 1) * boundarySize - 1, 0) += pow(-1, j + 1) * shiftFactor;
                ////NORTH
                //vectSol((boundarySize - 1) * boundarySize + j, 0) += pow(-1, j + 1) * shiftFactor;
                ////WEST
                //vectSol(j * boundarySize, 0) += pow(-1, j + 1) * shiftFactor;
                for (int k = 1; k < boundarySize - 1; ++k)
                {
                    vectSol(k * boundarySize + j, 0) += pow(-1, k) * innershiftFactor;
                    //gsInfo << k << " " << " " << j << " " << k * boundarySize + j << "\n";
                    vectSol(j * boundarySize + k, 1) += pow(-1, k) * innershiftFactor;
                    //gsInfo << j << " "  << k << " " << j * boundarySize + k << "\n";
                }
            }
        }
        /*
        if (i == 0)
         {
             vectSol(6, 0) += shiftFactor;
             vectSol(12, 0) -= shiftFactor;
             vectSol(18, 0) += shiftFactor;
             vectSol(24, 0) -= shiftFactor;
         }
         if (i == 1)
         {
             vectSol(5, 0) = newmp.patch(0).coef(0, 0);
             vectSol(11, 0) = newmp.patch(0).coef(6, 0);
             vectSol(17, 0) = newmp.patch(0).coef(12, 0);
             vectSol(23, 0) = newmp.patch(0).coef(18, 0);
             vectSol(29, 0) = newmp.patch(0).coef(24, 0);
             vectSol(35, 0) = newmp.patch(0).coef(30, 0);
             vectSol(5, 1) = newmp.patch(0).coef(0, 1);
             vectSol(11,1)   = newmp.patch(0).coef(6, 1);
             vectSol(17, 1) = newmp.patch(0).coef(12, 1);
             vectSol(23, 1) = newmp.patch(0).coef(18, 1);
             vectSol(29, 1) = newmp.patch(0).coef(24, 1);
             vectSol(35, 1) = newmp.patch(0).coef(30, 1);
         }*/
         /*  if (i == 2)
           {
               vectSol(5,0) = newmp.patch(0).coef(5, 0);
               vectSol(11, 0) = newmp.patch(0).coef(4, 0);
               vectSol(17, 0) = newmp.patch(0).coef(3, 0);
               vectSol(23, 0) = newmp.patch(0).coef(2, 0);
               vectSol(29, 0) = newmp.patch(0).coef(1, 0);
               vectSol(35, 0) = newmp.patch(0).coef(0, 0);
               vectSol(5, 1) = newmp.patch(0).coef(5, 1);
               vectSol(11, 1) = newmp.patch(0).coef(4, 1);
               vectSol(17, 1) = newmp.patch(0).coef(3, 1);
               vectSol(23, 1) = newmp.patch(0).coef(2, 1);
               vectSol(29, 1) = newmp.patch(0).coef(1, 1);
               vectSol(35, 1) = newmp.patch(0).coef(0, 1);
           }*/
           //{


           //if (i == 0)
           //{
           //   /* for (auto a = 0; a < boundarySize - 1; a++)
           //    {
           //        vectSol((a + 1) * boundarySize, 0) += pow(-1, a) * 0.017*3;
           //    }*/
           //    auto a = 0;
           //    vectSol(8, 0) += pow(-1, a) * 0.017 * 3;
           //}
        if (i == 1)
        {
            for (auto a = 0; a < boundarySize; a++)
            {
                vectSol((a + 1) * boundarySize - 1, 0) = newmp.patch(0).coef(a * boundarySize, 0);
                vectSol((a + 1) * boundarySize - 1, 1) = newmp.patch(0).coef(a * boundarySize, 1);
                gsInfo << "Equalizing " << (a + 1) * boundarySize - 1 << " from the 1st and " << a * boundarySize << " from the 0th\n";
            }
        }
        if (i == 2)
        {
            for (auto a = 0; a < boundarySize - 1; a++)
            {
                vectSol((boundarySize - 1) * boundarySize + a, 1) = newmp.patch(1).coef(a, 1);
                vectSol((boundarySize - 1) * boundarySize + a, 0) = newmp.patch(1).coef(a, 0);
                gsInfo << "Equalizing " << (boundarySize - 1) * boundarySize + a << " from the 2nd and " << a << " from the 1st\n";
            }
            gsInfo << "\n";
            for (auto a = 0; a < boundarySize; a++)
            {
                vectSol((boundarySize - a) * boundarySize - 1, 1) = newmp.patch(0).coef(a, 1);
                vectSol((boundarySize - a) * boundarySize - 1, 0) = newmp.patch(0).coef(a, 0);
                gsInfo << "Equalizing " << (boundarySize - a) * boundarySize - 1 << " from the 2nd and " << a << " from the 0th\n";
            }
            //vectSol(31, 1) = newmp.patch(1).coef(1, 1);
            //vectSol(32, 1) = newmp.patch(1).coef(2, 1);
            //vectSol(33, 1) = newmp.patch(1).coef(3, 1);
            //vectSol(34, 1) = newmp.patch(1).coef(4, 1);


            //newmp.patch(0).coef(4, 0) = vectSol(11, 0);
            //newmp.patch(0).coef(4, 1) = vectSol(11, 1);
            /*  vectSol(29, 0) = newmp.patch(0).coef(1, 0);
              vectSol(23, 0) = newmp.patch(0).coef(2, 0);
              vectSol(17, 0) = newmp.patch(0).coef(3, 0);
              vectSol(11, 0) = newmp.patch(0).coef(4, 0);
              vectSol(35, 0) = newmp.patch(0).coef(0, 0);
              vectSol(29, 1) = newmp.patch(0).coef(1, 1);
              vectSol(23, 1) = newmp.patch(0).coef(2, 1);
              vectSol(17, 1) = newmp.patch(0).coef(3, 1);
              vectSol(11, 1) = newmp.patch(0).coef(4, 1);
              vectSol(35, 1) = newmp.patch(0).coef(0, 1);*/

        }

        //for (size_t coeffOfGeom = boundarySize; coeffOfGeom < vectSol.rows() - boundarySize; coeffOfGeom++)
        //{
        //    if ((coeffOfGeom % boundarySize != 0) && (coeffOfGeom % boundarySize != boundarySize - 1))
        //    {
        //        auto coeff = 5;
        //        if (i == 2) coeff *= 0.6;
        //        vectSol(coeffOfGeom,0) += pow(-1, coeffOfGeom) * 0.017 * coeff;
        //        vectSol(coeffOfGeom,1) += pow(-1, coeffOfGeom) * 0.017 * 5;
        //    }
        //}


        /*       if (i == 0)
               {
                   vectSol(4, 0) += 0.00;*/
                   //        //gsInfo << "ZURNA! " << std::sin(oscillationFrequencyV * PI * vectSol(11, 0)) << "\n";
                   //        /*vectSol(11, 0) += epsilon* std::sin(oscillationFrequencyV * PI * vectSol(11, 0));
                   //        vectSol(17, 0) += epsilon* std::sin(oscillationFrequencyV * PI * vectSol(17, 0));
                   //        vectSol(23, 0) += epsilon* std::sin(oscillationFrequencyV * PI * vectSol(23, 0));
                   //        vectSol(29, 0) += epsilon* std::sin(oscillationFrequencyV * PI * vectSol(29, 0));*/

                   //        for (size_t coeffOfGeom = 1; coeffOfGeom < vectSol.rows() - 1; coeffOfGeom++)
                   //        {
                   //            if (
                   //                //((coeffOfGeom % boundarySize) == 0) ||
                   //                 //((coeffOfGeom % boundarySize) == boundarySize - 1)
                   //                //||
                   //                //coeffOfGeom / boundarySize == 0 ||
                   //                coeffOfGeom / boundarySize == boundarySize - 1 && (coeffOfGeom % boundarySize) != 0//north
                   //                )
                   //            {
                   //                //vectSol.row(coeffOfGeom) *= 1 + (std::rand() % 1000 -1000) * epsilon;
                   //                vectSol(coeffOfGeom, 1) += std::sin(oscillationFrequency * PI * vectSol(coeffOfGeom, 1)) * epsilon;
                   //                continue;
                   //            }
                   //            if (
                   //                //((coeffOfGeom % boundarySize) == 0)
                   //                //||
                   //                ((coeffOfGeom % boundarySize) == boundarySize - 1) && (coeffOfGeom / boundarySize) != 0//east
                   //                )
                   //            {
                   //                vectSol(coeffOfGeom, 0) += epsilon * std::sin(oscillationFrequencyV * PI * vectSol(coeffOfGeom, 1));
                   //            }
                   //            //if(coeffOfGeom / boundarySize == 0 || coeffOfGeom / boundarySize == boundarySize - 1)
                   //            //{
                   //            //    //vectSol.row(coeffOfGeom) *= 1 + (std::rand() % 1000 -1000) * epsilon;
                   //            //    vectSol(coeffOfGeom,0) +=  std::sin(oscillationFrequency*PI*vectSol(coeffOfGeom,0)) * epsilon;
                   //            //}
                   //            /*if (coeffOfGeom == 1)
                   //            {
                   //                vectSol(coeffOfGeom, 1) += pointShift;
                   //            }*/
                   //        }
                   /*   }
                      if (i == 1)
                      {
                          vectSol(3, 0) += 0.00;*/
                          //        for (size_t coeffOfGeom = 1; cfoeffOfGeom < vectSol.rows() - 1; coeffOfGeom++)
                          //        {
                          //            if (
                          //                //((coeffOfGeom % boundarySize) == 0) ||
                          //                 //((coeffOfGeom % boundarySize) == boundarySize - 1)
                          //                //||
                          //                //coeffOfGeom / boundarySize == 0 ||
                          //                coeffOfGeom / boundarySize == boundarySize - 1 && (coeffOfGeom % boundarySize) != 0//north
                          //                )
                          //            {
                          //                //vectSol.row(coeffOfGeom) *= 1 + (std::rand() % 1000 -1000) * epsilon;
                          //                vectSol(coeffOfGeom, 1) += std::sin(oscillationFrequency * PI * vectSol(coeffOfGeom, 1)) * epsilon;
                          //                continue;
                          //            }
                          //            if (
                          //                //((coeffOfGeom % boundarySize) == 0)
                          //                //||
                          //                ((coeffOfGeom % boundarySize) == 0) && (coeffOfGeom / boundarySize) != boundarySize - 1//west
                          //                )
                          //            {
                          //                vectSol(coeffOfGeom, 0) += epsilon * std::sin(oscillationFrequencyV * PI * vectSol(coeffOfGeom, 1));
                          //            }
                          //            // if(coeffOfGeom / boundarySize == 0 || coeffOfGeom / boundarySize == boundarySize - 1)
                          //            //{
                          //            //    //vectSol.row(coeffOfGeom) *= 1 + (std::rand() % 1000 -1000) * epsilon;
                          //            //    vectSol(coeffOfGeom,0) +=  std::sin(oscillationFrequency*PI*vectSol(coeffOfGeom,0)) * epsilon;
                          //            //}
                          //        }
                          //}
                          //    if (i == 2)
                          //    {
                          //        for (size_t coeffOfGeom = 1; coeffOfGeom < vectSol.rows() - 1; coeffOfGeom++)
                          //        {
                          //            if (
                          //                //((coeffOfGeom % boundarySize) == 0) ||
                          //                 //((coeffOfGeom % boundarySize) == boundarySize - 1)
                          //                //||
                          //                //coeffOfGeom / boundarySize == 0 ||
                          //                coeffOfGeom / boundarySize == 0 && (coeffOfGeom % boundarySize) != boundarySize - 1//south
                          //                )
                          //            {
                          //                //vectSol.row(coeffOfGeom) *= 1 + (std::rand() % 1000 -1000) * epsilon;
                          //                vectSol(coeffOfGeom, 1) += std::sin(oscillationFrequency * PI * vectSol(coeffOfGeom, 1)) * epsilon;
                          //                continue;
                          //            }
                          //            if (
                          //                //((coeffOfGeom % boundarySize) == 0)
                          //                //||
                          //                ((coeffOfGeom % boundarySize) == 0) && (coeffOfGeom / boundarySize) != boundarySize - 1//west
                          //                )
                          //            {
                          //                vectSol(coeffOfGeom, 0) += epsilon * std::sin(oscillationFrequencyV * PI * vectSol(coeffOfGeom, 1));
                          //            }
                          //            // if(coeffOfGeom / boundarySize == 0 || coeffOfGeom / boundarySize == boundarySize - 1)
                          //            //{
                          //            //    //vectSol.row(coeffOfGeom) *= 1 + (std::rand() % 1000 -1000) * epsilon;
                          //            //    vectSol(coeffOfGeom,0) +=  std::sin(oscillationFrequency*PI*vectSol(coeffOfGeom,0)) * epsilon;
                          //            //}
                          //        }
                          //        //vectSol(29, 1) = newmp.patch(0).coef(1, 1);
                          //    }
                          //}
                          //if (i < 2)
                          //if (true)
                          //{
                          //    for (size_t coeffOfGeom = 1; coeffOfGeom < vectSol.rows() - 1; coeffOfGeom++)
                          //    {
                          //        if (((coeffOfGeom % boundarySize) == 0) || ((coeffOfGeom % boundarySize) == boundarySize - 1))
                          //        {
                          //            gsInfo << "TEST, " << coeffOfGeom << "," << vectSol(coeffOfGeom, 1) << "\n";
                          //           /* if (coeffOfGeom == 6 || coeffOfGeom == 11)
                          //            {*/
                          //                gsInfo << "YUHU, i = " << i << ", " << std::sin(oscillationFrequencyV * PI * vectSol(coeffOfGeom, 1)) << "\n";
                          //                vectSol(coeffOfGeom, 0) += epsilon * std::sin(oscillationFrequencyV * PI * vectSol(coeffOfGeom, 1));// std::sin(oscillationFrequencyV * PI * vectSol(coeffOfGeom, 1))* epsilon;
                          //                gsInfo << "ADDED" << std::sin(oscillationFrequencyV * PI * vectSol(coeffOfGeom, 1)) * epsilon << "\n";
                          //            //}
                          //        }
                          //    }
                          //}
                          /*if (i == 2)
                          {
                              for (size_t j = 0; j <= boundarySize - 1; j++)
                              {
                                  vectSol.row(vectSol.rows() - boundarySize*j - 1) = newmp.patch(0).coefs().row(j);
                                  gsInfo << "Gluing: " << vectSol.row(vectSol.rows() - boundarySize * j - 1) << "\n";
                              }
                          }*/

                          //    bb->uniformRefine(1);
                          ////        gsInfo << *bb << std::endl;
                          ////        gsInfo << "i = " << i << std::endl;
                          ////        gsInfo << "FOUND " << lowCorners.rows() << " BOXES\n";
                          ////        gsInfo << lowCorners << "\n\n" << upCorners << "\n\n" << Level << "\n\n";
                          ////        gsFitting<> fitting(coefsForTHB, )
                          ////        gsInfo << "i = " << i << "\n" << mp->patch(i).coefs();
                          ////        GISMO_ENSURE(nullptr!=bb, "Conversion error.");
                          ////        for ( size_t k = 0; k!=d; ++k )
                          ////            bb->component(k).knots().transform(0,1);
                          //


        newmp.addPatch(THB.makeGeometry(vectSol));
        //for (size_t j = 0; j < uv1.cols(); j++)
        //{
        //    auto jacdet = newmp.patch(i).jacobian(uv1.col(j)).determinant();
        //    if (jacdet <= 0)
        //    {
        //        gsInfo << "Problem! i = " << i << ",j = " << j << ",\n point\n" << uv1.col(j) << ",\n jacobian = " << jacdet << "\n";
        //    }
        //    //            gsInfo << newmp.patch(i).jacobian(uv1.col(j)) << "\n";

        //}
        gsInfo << "next\n";
    }
    gsInfo << PatchesIntersection(newmp.patch(0), newmp.patch(1)) << "\t";
    gsInfo << PatchesIntersection(newmp.patch(0), newmp.patch(2)) << "\t";
    gsInfo << PatchesIntersection(newmp.patch(1), newmp.patch(2)) << "\n";
    //
    //    if (reparam) // Reparameterize?
    //    {
    //        gsKnotVector<> kv(0,1,nknots,degree+1,1);
    //        gsTensorBSplineBasis<2> bs(kv,kv);
    //        gsMultiPatch<> newmp;
    //        // Using interpolation
    //        if (0==method)
    //        {
    //            gsMatrix<> values, anchors = bs.anchors();
    //            for ( size_t i = 0; i< mp->nPatches(); ++i )
    //            {
    //                values = mp->patch(i).eval(anchors);
    //                newmp.addPatch( bs.interpolateData(values,anchors) );
    //            }
    //        }
    //
    //        // Using Fitting
    //
    //        if (1==method)
    //        {
    //            gsVector<index_t,2> sz; //str
    //            gsVector<unsigned> sz2(2);
    //
    //            gsVector<> a(2); a.setZero();
    //            gsVector<> b(2); b.setOnes();
    //            gsMatrix<> param;
    //
    //            for ( size_t i = 0; i< mp->nPatches(); ++i )
    //            {
    //                gsTensorBSplineBasis<2> * bb =
    //                    dynamic_cast<gsTensorBSplineBasis<2> *>(&mp->patch(i).basis().source());
    //                bb->size_cwise(sz);
    //                sz2[0]= sz[0];
    //                sz2[1]= sz[1];
    //                param = gsPointGrid(a,b,sz2);
    //                gsFitting<>  fitting(param, mp->patch(i).coefs().transpose(), bs);
    //                fitting.compute(0.0000001);
    //                newmp.addPatch( *fitting.result() );
    //
    //                fitting.computeMaxNormErrors();
    //                gsInfo<< "Min error: "<< fitting.minPointError() <<"\n";
    //                gsInfo<< "Max error: "<< fitting.maxPointError() <<"\n";
    //            }
    //        }
    //
    //        if (1==method)
    //        {
    //            gsVector<index_t,2> sz; //str
    //            gsVector<unsigned> sz2(2);
    //
    //            gsVector<> a(2); a.setZero();
    //            gsVector<> b(2); b.setOnes();
    //            gsMatrix<> param;
    //
    //            for ( size_t i = 0; i< mp->nPatches(); ++i )
    //            {
    //                gsTensorBSplineBasis<2> * bb =
    //                    dynamic_cast<gsTensorBSplineBasis<2> *>(&mp->patch(i).basis().source());
    //                bb->size_cwise(sz);
    //                sz2[0]= sz[0];
    //                sz2[1]= sz[1];
    //                param = gsPointGrid(a,b,sz2);
    //                gsFitting<>  fitting(param, mp->patch(i).coefs().transpose(), bs);
    //                fitting.compute(0.0000001);
    //                newmp.addPatch( *fitting.result() );
    //
    //                fitting.computeMaxNormErrors();
    //                gsInfo<< "Min error: "<< fitting.minPointError() <<"\n";
    //                gsInfo<< "Max error: "<< fitting.maxPointError() <<"\n";
    //            }
    //        }
    //
    ////        mp->swap(newmp);
    //    }
    //
    //    gsInfo <<"Computing"<< (gaps?" corner-based ":" ")<<"topology with tolerance = "<<tol<<"... \n" ;
    //    mp->computeTopology(tol, gaps);
    //    gsInfo << * mp <<"\n" ;
    //    //gsInfo << mp->detail() <<"\n";
    //    if (gaps) //Close gaps?
    //    {
    //        gsInfo <<"Closing gaps with tolerance = "<<gtol<<"... \n" ;
    //        mp->closeGaps(gtol);
    //        gsInfo <<"Computing topology with tolerance = "<<tol<<"... \n" ;
    //        mp->computeTopology(tol, false);
    //        mp->closeGaps(gtol);//close again for any previously missed interfaces
    //        gsInfo << * mp <<"\n" ;
    //        //gsInfo << mp->detail() <<"\n";
    //    }
    ////    geom->setCoefs(coefsR);
    gsFileData<> fdG;
    ////    fdG << *geom;
    newmp.computeTopology(tol, false);
    fdG << newmp;
    //    gsFileData<> fd;
    ////    fdG.dump("GeomR");
    //fdG.dump("foot");
    //    fd<< *mp ;
    fdG.dump("makeMultipatch_output");

    //gsMatrix<> ;
    //    gsInfo <<"Resulting file is written out to makeMultipatch_output.xml\n" ;
    //    gsWriteParaview(*mp, "two_squares", 1000, true, false );
    gsWriteParaview(newmp, "makeMultipatchFile", 1000, true, false);
    gsFileManager::open("makeMultipatchFile.pvd");
    gsFileManager::open("makeMultipatch_output.xml");
    ////    gsWriteParaview(*geom, "GeomR");
    //    return EXIT_SUCCESS;
    return 0;
}