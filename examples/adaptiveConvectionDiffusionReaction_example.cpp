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


/**
 * @brief Performs THB fitting on each patch and builds a multipatch geometry.
 *
 * @param[in] mp Input multipatch geometry
 * @param[in] lowc2 Lower corner of domain
 * @param[in] uppc2 Upper corner of domain
 * @param[in] thb1 THB basis for refinement reference
 * @param[out] newmp Output multipatch with fitted geometry
 */
void performFitting(
    const gsMultiPatch<>::uPtr& mp,
    const gsVector<real_t>& lowc2,
    const gsVector<real_t>& uppc2,
    const gsTHBSplineBasis<2, real_t>& thb1,
    gsMultiPatch<>& newmp)
{
    gsMultiPatch<> mp1;
    gsMatrix<> punto(2, 1);
    punto(0, 0) = punto(1, 0) = 0.90141560817564836;
    auto mateval = thb1.function(0).eval(punto);

    for (size_t i = 0; i < mp->nPatches(); ++i)
    {
        gsTensorBSplineBasis<2>* bb =
            dynamic_cast<gsTensorBSplineBasis<2> *>(&mp->patch(i).basis().source());

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

        real_t au = 0;
        real_t bu = 1;
        index_t interior = 0;
        index_t degree = 2;
        index_t multEnd = degree + 1;
        gsKnotVector <> ku(au, bu, interior, multEnd);
        gsTensorBSplineBasis<2, real_t> tens(ku, ku);
        gsTHBSplineBasis<2, real_t> THB(tens);

        if (i == 0) {
            std::vector<index_t> box;
            box.push_back(4);
            box.push_back(0);
            box.push_back(0);
            box.push_back(16);
            box.push_back(16);
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

        gsFitting<> fit(uv1, xy1, THB);
        fit.compute();
        fit.computeErrors();
        gsInfo << fit.getBasis() << "\n";
        
        gsSparseMatrix<real_t> matA(THB.size(), THB.size());
        gsMatrix<real_t> matAfull(THB.size(), THB.size());
        gsMatrix<real_t> vectB(THB.size(), 2);
        matA.setZero();
        matAfull.setZero();
        vectB.setZero();
        fit.assembleSystem(matA, vectB);
        matAfull = gsEigen::MatrixXd(matA);
        printTheMatrix(matA, "matA");
        gsFileManager::open("biharmonic_example.txt");

        gsMatrix<real_t> vectSol(THB.size(), 2);
        vectSol = matAfull.partialPivLu().solve(vectB);
        
        auto shiftFactor = 0.00;
        int boundarySize = (int)std::sqrt(vectSol.rows());
        auto innershiftFactor = 0.00;

        gsInfo << "boundarySize = " << boundarySize << "\n";
        
        if (i == 0)
        {
            for (int j = 1; j < boundarySize - 1; ++j) {}
        }
        
        if (i == 1)
        {
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
                for (int k = 1; k < boundarySize - 1; ++k)
                {
                    vectSol(k * boundarySize + j, 0) += pow(-1, k) * innershiftFactor;
                    vectSol(j * boundarySize + k, 1) += pow(-1, k) * innershiftFactor;
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
                for (int k = 1; k < boundarySize - 1; ++k)
                {
                    vectSol(k * boundarySize + j, 0) += pow(-1, k) * innershiftFactor;
                    vectSol(j * boundarySize + k, 1) += pow(-1, k) * innershiftFactor;
                }
            }
        }

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
        }

        newmp.addPatch(THB.makeGeometry(vectSol));
        gsInfo << "next\n";
    }
}

/**
 * @brief Performs automatic THB fitting by extracting boxes separately for each patch.
 *
 * For each patch, extracts refinement boxes from that patch's THB basis
 * and automatically applies them to a new THB basis.
 *
 * @param[in] mp Input multipatch geometry (must contain THB bases)
 * @param[out] newmp Output multipatch with fitted geometry
 */
void performAutomaticFitting(
    const gsMultiPatch<>::uPtr& mp,
    gsMultiPatch<>& newmp)
{
    for (size_t patch = 0; patch < mp->nPatches(); ++patch)
    {
        // Extract domain bounds from first patch
        gsVector<real_t> lowc = mp->patch(patch).basis().support().col(0);
        gsVector<real_t> uppc = mp->patch(patch).basis().support().col(1);

        // Extract reference THB basis for this patch
        const gsTHBSplineBasis<2, real_t>* refTHB = 
            dynamic_cast<const gsTHBSplineBasis<2, real_t>*>(&mp->patch(patch).basis().source());
        
        if (!refTHB)
        {
            gsInfo << "Error: Patch " << patch << " does not contain a THB basis\n";
            continue;
        }
        
        // Extract boxes from this patch's THB basis
        gsMatrix<index_t> lowCorners, upCorners;
        gsVector<index_t> levels;
        refTHB->tree().getBoxes(lowCorners, upCorners, levels);
        
        gsInfo << "\nPatch " << patch << ": Extracted " << lowCorners.rows() << " boxes\n";

        auto numpoints = 101;
        gsMatrix<> uv1 = uniformPointGrid(lowc, uppc, numpoints * numpoints);
        gsMatrix<> xy1 = mp->patch(patch).eval(uv1);

        // Validate input Jacobian
        for (size_t j = 0; j < uv1.cols(); j++)
        {
            if (auto jacdet = mp->patch(patch).jacobian(uv1.col(j)).determinant() <= 0)
            {
                gsInfo << "Problem! patch = " << patch << ", j = " << j 
                       << ", jacobian = " << jacdet << "\n";
            }
        }

        // Create base THB basis
        real_t au = 0;
        real_t bu = 1;
        index_t interior = 0;
        index_t degree = 2;
        index_t multEnd = degree + 1;
        gsKnotVector <> ku(au, bu, interior, multEnd);
        gsTensorBSplineBasis<2, real_t> tens(ku, ku);
        gsTHBSplineBasis<2, real_t> THB(tens);

        // Automatically refine using extracted boxes for this patch
        gsInfo << "Patch " << patch << " refinement boxes:\n";
        for (index_t boxIdx = 0; boxIdx < lowCorners.rows(); ++boxIdx)
        {
            std::vector<index_t> box;
            box.push_back(levels(boxIdx));           // level
            box.push_back(lowCorners(boxIdx, 0));    // x_min
            box.push_back(lowCorners(boxIdx, 1));    // y_min
            box.push_back(upCorners(boxIdx, 0));     // x_max
            box.push_back(upCorners(boxIdx, 1));     // y_max
            
            gsInfo << "  Box " << boxIdx << ": level=" << levels(boxIdx) 
                   << " x=[" << lowCorners(boxIdx, 0) << ", " << upCorners(boxIdx, 0) << "]"
                   << " y=[" << lowCorners(boxIdx, 1) << ", " << upCorners(boxIdx, 1) << "]\n";
            
            THB.refineElements(box);
        }
        gsInfo << "Applied " << lowCorners.rows() << " refinement boxes to patch " << patch << "\n";

        // Perform fitting
        gsFitting<> fit(uv1, xy1, THB);
        fit.compute();
        fit.computeErrors();
        gsInfo << fit.getBasis() << "\n";

        gsSparseMatrix<real_t> matA(THB.size(), THB.size());
        gsMatrix<real_t> matAfull(THB.size(), THB.size());
        gsMatrix<real_t> vectB(THB.size(), 2);
        matA.setZero();
        matAfull.setZero();
        vectB.setZero();
        fit.assembleSystem(matA, vectB);
        matAfull = gsEigen::MatrixXd(matA);

        gsMatrix<real_t> vectSol(THB.size(), 2);
        vectSol = matAfull.partialPivLu().solve(vectB);

        gsInfo << "Fitted patch " << patch << " with " << THB.size() << " basis functions\n";

        newmp.addPatch(THB.makeGeometry(vectSol));
    }

    gsInfo << "Automatic fitting completed for all patches\n";
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
    //std::string filename("hexagon_3p_4l_coons1_test.xml");
    //std::string filename("savedMultipatch.xml");
    std::string filename("savedMultipatch_1.xml");
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

    gsKnotVector <> ku(0, 1, 0, 3);
    gsTensorBSplineBasis<2, real_t> tens(ku, ku);
    gsTHBSplineBasis<2, real_t> thb1(tens);
    std::vector<index_t> box;
    
    box.push_back(3);
    box.push_back(0);
    box.push_back(7);
    box.push_back(8);
    box.push_back(8);
    thb1.refineElements(box);
    gsInfo << "THB size: " << thb1.size() << "\n";
    
    box.clear();
    box.push_back(4);
    box.push_back(0);
    box.push_back(0);
    box.push_back(16);
    box.push_back(14);
    thb1.refineElements(box);
    gsInfo << "THB size: " << thb1.size() << "\n";
    box.clear();
    //-------300725



    gsMultiPatch<>::uPtr mp = gsReadFile<>(filename);
    gsVector<> lowc1(2);
    gsVector<> uppc1(2);
    lowc1(0) = 0;
    lowc1(1) = 0;
    uppc1(0) = 1;
    uppc1(1) = 1;
    const gsVector<real_t> lowc2 = lowc1;
    const gsVector<real_t> uppc2 = uppc1;
    
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
                gsInfo << "Problem in input! i = " << i << ",j = " << j << ", jacobian = " << jacdet << "\n";
            }
        }
    }

    gsMultiPatch<> newmp;
    
    // Menu for choosing fitting type
    gsInfo << "\n=== Fitting Type Selection ===\n";
    gsInfo << "1. Manual Fitting (performFitting)\n";
    gsInfo << "2. Automatic Fitting (performAutomaticFitting)\n";
    gsInfo << "Enter your choice (1 or 2): ";
    
    int choice;
    std::cin >> choice;
    
    if (choice == 2)
    {
        gsInfo << "Selected: Automatic Fitting\n";
        performAutomaticFitting(mp, newmp);
    }
    else
    {
        gsInfo << "Selected: Manual Fitting (default)\n";
        performFitting(mp, lowc2, uppc2, thb1, newmp);
    }

    gsInfo << PatchesIntersection(newmp.patch(0), newmp.patch(1)) << "\t";
    gsInfo << PatchesIntersection(newmp.patch(0), newmp.patch(2)) << "\t";
    gsInfo << PatchesIntersection(newmp.patch(1), newmp.patch(2)) << "\n";

    gsFileData<> fdG;
    newmp.computeTopology(tol, false);
    fdG << newmp;
    fdG.dump("makeMultipatch_output");

    gismo::gsWriteParaview(newmp, "makeMultipatchFile", 1000, true, false);
    gsFileManager::open("makeMultipatchFile.pvd");
    gsFileManager::open("makeMultipatch_output.xml");

    return 0;
}
