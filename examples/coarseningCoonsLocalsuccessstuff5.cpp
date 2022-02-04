#include <iostream>

#include <iostream>

#include <gismo.h>

#include "math.h"

#include <vector>

#include <gsIO/gsIOUtils.h>

#include "/home/turing/theydarov/gismo/extensions/motor/jku/gsQualityMeasure2.h"

#include "/home/turing/theydarov/gismo/extensions/motor/jku/gsMotorUtils.h"



using namespace gismo;
using namespace std;

int numRed, fineLevel = 3, coarseLevel = 2, theLev;
int numCells = 4 * (int) pow(2, coarseLevel) * 4 * (int) pow(2, coarseLevel);
int redBox[4096][5]; //LAZY SOLUTION
ofstream outfile;





void boxToDomain(int mybox[][5], double coords[][5]) {
    for (int i = 0; i < numRed; i++) {
        coords[i][0] = (double) mybox[i][0];
        coords[i][1] = (double) mybox[i][1] / ((double)(4 * pow(2, mybox[i][0])));
        coords[i][2] = (double) mybox[i][2] / ((double)(4 * pow(2, mybox[i][0])));
        coords[i][3] = (double) mybox[i][3] / ((double)(4 * pow(2, mybox[i][0])));
        coords[i][4] = (double) mybox[i][4] / ((double)(4 * pow(2, mybox[i][0])));
    }
}

void boxToDomain(int mybox[5], double coords[5]) {
    coords[0] = (double) mybox[0];
    coords[1] = (double) mybox[1] / ((double)(4 * pow(2, mybox[0])));
    coords[2] = (double) mybox[2] / ((double)(4 * pow(2, mybox[0])));
    coords[3] = (double) mybox[3] / ((double)(4 * pow(2, mybox[0])));
    coords[4] = (double) mybox[4] / ((double)(4 * pow(2, mybox[0])));

}

void saveData(const gsGeometry < > & geom,
              const std::string & output,
              const int number) {
    gsFileData < > fileData;
    fileData << geom;

    std::string out = output + "_Map_" + util::to_string(number) + ".xml";
    fileData.dump(out);
    //    gsInfo
    outfile << "geometry saved to " << out << "\n";
    gsMesh < > mesh;
    makeMesh(geom.basis(), mesh, 1600);
    geom.evaluateMesh(mesh);
    out = output + "_Mesh_" + util::to_string(number);
    gsWriteParaview(mesh, out);
    gsInfo << "Saved to " << out << "\n";
    outfile << "Saved to " << out << "\n";
}

void saveData(const gsGeometry < > & geom,
              const std::string & output
              ){
    gsFileData < > fileData;
    fileData << geom;

    std::string out = output + ".xml";
    fileData.dump(out);
    //    gsInfo
    outfile << "geometry saved to " << out << "\n";
    gsMesh < > mesh;
    makeMesh(geom.basis(), mesh, 1600);
    geom.evaluateMesh(mesh);
    out = output;
    gsWriteParaview(mesh, out);
    gsInfo << "Saved to " << out << "\n";
    outfile << "Saved to " << out << "\n";
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void boundingBox(gsTHBSplineBasis < 2 > & THBSplineBasis, gsMatrix < int > & maxLevels, gsMatrix < int > & minLevels, gsMatrix < real_t > & matrixOfTPLC, gsMatrix < real_t > & matrixOfTPUC, gsMatrix<int> basisInd) {

    //Compute the bounding box considering only the support(optimize code!!)
    gsVector < real_t > centerElement;
    gsMatrix < real_t > centerElementTP(2, 1);
    gsMatrix < index_t > resActives;
    gsMatrix < index_t > resActivesInner;
    gsMatrix < > mySupport;
    gsMatrix < index_t > resActivesTP;
    index_t totalNumberBasisFunct = THBSplineBasis.size();
    matrixOfTPUC.setZero(2, totalNumberBasisFunct);
    matrixOfTPLC.setZero(2, totalNumberBasisFunct);

    maxLevels.setZero(1, totalNumberBasisFunct);
    minLevels.setZero(1, totalNumberBasisFunct);
    vector < double > centPointsx;
    vector < double > centPointsy;

    real_t minTPLCx;
    real_t minTPLCy;
    real_t maxTPUCx;
    real_t maxTPUCy;
    //loop over basis funct
    vector < double > functionSupportTPLCx;
    vector < double > functionSupportTPLCy;
    vector < double > functionSupportTPUCx;
    vector < double > functionSupportTPUCy;
    vector < int > levels;
    int functionLevel;
    int functionLevelOfRes;
    gsMatrix < real_t > minTPLC(2, 1);
    gsMatrix < real_t > maxTPUC(2, 1);
    index_t numberActiveFunct;
    gsTensorBSplineBasis < 2 > tensorBasis;
    for (index_t ss = 0; ss != THBSplineBasis.size(); ss++)
    { // HIER. ID
        if (basisInd(ss)){
            //ss indicates the hierarchical index of the function
            functionLevel = THBSplineBasis.levelOf(ss);

            gsHDomainIterator<real_t, 2> domIter(THBSplineBasis);
            // loop over all elements
            for (; domIter.good(); domIter.next())
            {

                //Find active function in the element
                centerElement = domIter.centerPoint();

                THBSplineBasis.active_into(centerElement, resActives);
                numberActiveFunct = resActives.size();

                //For each basis funct determine if the selected basis funct is active in this element
                // loop over active functions in the element
                for (index_t ii = 0; ii != numberActiveFunct; ii++)
                {

                    //Check if myactive function is active in the element.
                    //If yes determine the support of each function active in the element
                    if (resActives(ii, 0) == ss)
                    {
                        //Find the center points of the element which forms the support of the selected basisfunct
                        centPointsx.push_back(centerElement(0, 0));
                        centPointsy.push_back(centerElement(1, 0));

                        for (index_t kk = 0; kk != numberActiveFunct; kk++)
                        {
                            functionLevelOfRes = THBSplineBasis.levelOf(resActives(kk, 0));
                            levels.push_back(functionLevelOfRes);

                        }

                    }

                }

            }

            maxLevels(0, ss) = *std::max_element(levels.begin(), levels.end());
            minLevels(0, ss) = *std::min_element(levels.begin(), levels.end());

            //Now I have the minimum level so I can initialize the corresponding TP basis from the thb basis if needed

            tensorBasis = *THBSplineBasis.getBases()[minLevels(0, ss)];

            matrixOfTPLC(0, ss) = pow(2, THBSplineBasis.maxLevel() - THBSplineBasis.levelOf(ss))
                                  * (*THBSplineBasis.getBases()[THBSplineBasis.levelOf(ss)])
                                          .elementSupport(THBSplineBasis.flatTensorIndexOf(ss))(0, 0);
            matrixOfTPLC(1, ss) = pow(2, THBSplineBasis.maxLevel() - THBSplineBasis.levelOf(ss))
                                  * (*THBSplineBasis.getBases()[THBSplineBasis.levelOf(ss)])
                                          .elementSupport(THBSplineBasis.flatTensorIndexOf(ss))(1, 0);
            matrixOfTPUC(0, ss) = pow(2, THBSplineBasis.maxLevel() - THBSplineBasis.levelOf(ss))
                                  * (*THBSplineBasis.getBases()[THBSplineBasis.levelOf(ss)])
                                          .elementSupport(THBSplineBasis.flatTensorIndexOf(ss))(0, 1);
            matrixOfTPUC(1, ss) = pow(2, THBSplineBasis.maxLevel() - THBSplineBasis.levelOf(ss))
                                  * (*THBSplineBasis.getBases()[THBSplineBasis.levelOf(ss)])
                                          .elementSupport(THBSplineBasis.flatTensorIndexOf(ss))(1, 1);
        }
    }

}

bool doOverlap(gsMatrix < real_t > box_a, gsMatrix < real_t > box_b) {
    // If one rectangle is on left side of other
    //gsInfo << box_a(0, 0) << ' ' << box_a(1, 0) << "\n" << box_a(0, 1) << ' ' << box_a(1, 1) << "\n";
    //gsInfo << "\n";
    if (box_a(0, 0) >= box_b(0, 1) || box_b(0, 0) >= box_a(0, 1))
        return false;

    // If one rectangle is above other
    if (box_a(1, 1) <= box_b(1, 0) || box_b(1, 1) <= box_a(1, 0))
        return false;
    //gsInfo << "problems with grading" << "\n";
    return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void updateCriterion(int levNow, int x1U, int y1U, int x2U, int y2U, int cellExtension, gsMatrix<int> & basisInd, gsTHBSplineBasis<2, double> THB){
    real_t cellCoords[5];
    int cellIndices[5];
    cellIndices[0] = levNow;
    cellIndices[1] = std::min(4 * (int)pow(2,levNow), std::max(x1U, 0));
    cellIndices[2] = std::min(4 * (int)pow(2,levNow),std::max(y1U, 0));
    cellIndices[3] = std::max(0,std::min(x2U, 4 * (int)pow(2,levNow)));
    cellIndices[4] = std::max(0,std::min(y2U, 4 * (int)pow(2,levNow)));
//    gsInfo << cellIndices[0] << "\t" << cellIndices[1] << "\t" << cellIndices[2] << "\t" << cellIndices[3] << "\t" << cellIndices[4] << "\n";
    boxToDomain(cellIndices, cellCoords);
//    gsInfo << cellCoords[0] << "\t" << cellCoords[1] << "\t" << cellCoords[2] << "\t" << cellCoords[3] << "\t" << cellCoords[4] << "\n";

    gsMatrix<index_t> resultNE, resultNW, resultSE, resultSW;
    gsMatrix<double> noqte(2, 1);
    //NE
    noqte(0,0) = cellCoords[1] + 0.75*(cellCoords[3] - cellCoords[1]);
    noqte(1,0) =  cellCoords[2] + 0.75*(cellCoords[4] - cellCoords[2]);
    THB.active_into(noqte, resultNE);
    gsMatrix<double> qiymet;
    for (int i = 0; i < resultNE.rows(); ++i)
    {
        THB.evalSingle_into(resultNE(i,0), noqte, qiymet);
        if(qiymet(0,0)!=0.0)
        {
//                        gsInfo << resultNE(i, 0) << "is nonzero\n";
            basisInd(resultNE(i,0)) = 1;
        }
    }
    //NW
    noqte(0,0) = cellCoords[1] + 0.25*(cellCoords[3] - cellCoords[1]);
    THB.active_into(noqte, resultNW);
//                gsInfo << resultNW << "\n";
    for (int i = 0; i < resultNW.rows(); ++i)
    {
        THB.evalSingle_into(resultNW(i,0), noqte, qiymet);
        if(qiymet(0,0)!=0.0)
        {
//                        gsInfo << resultNW(i, 0) << "is nonzero\n";
            basisInd(resultNW(i,0)) = 1;
        }
    }
    //SW
    noqte(1,0) = cellCoords[2] + 0.25*(cellCoords[4] - cellCoords[2]);
    THB.active_into(noqte, resultSW);
//                gsInfo << resultSW << "\n";
    for (int i = 0; i < resultSW.rows(); ++i)
    {
        THB.evalSingle_into(resultSW(i,0), noqte, qiymet);
        if(qiymet(0,0)!=0.0)
        {
//                        gsInfo << resultSW(i, 0) << "is nonzero\n";
            basisInd(resultSW(i,0)) = 1;
        }
    }
    //SE
    noqte(0,0) = cellCoords[1] + 0.25*(cellCoords[3] - cellCoords[1]);
    THB.active_into(noqte, resultSE);
//                gsInfo << resultSW << "\n";
    for (int i = 0; i < resultSE.rows(); ++i)
    {
        THB.evalSingle_into(resultSE(i,0), noqte, qiymet);
        if(qiymet(0,0)!=0.0)
        {
//                        gsInfo << resultSE(i, 0) << "is nonzero\n";
            basisInd(resultSE(i,0)) = 1;
        }
    }
//    gsInfo << basisInd << "\n\n";

}

gsMatrix < > computeGaussPoints(gsTHBSplineBasis < 2 > & THBSplineBasis) {
    gsHDomainIterator < real_t, 2 > domIter(THBSplineBasis); //Create domain iterator
    int numberOfElements = THBSplineBasis.numElements(); //Find number of elements
    gsGaussRule < > gr(THBSplineBasis, 2, 2);
    gsVector < > v;
    gsMatrix < > pts, vectorGaussPoints(2, gr.numNodes() * numberOfElements); //Vector that contains the coordinate of the total Gauss points of the entire mesh

    int counter = 0;
    for (; domIter.good(); domIter.next()) // loop over all elements
    {
        //Find coordinate of each Gauss point
        gr.mapTo(domIter.lowerCorner(), domIter.upperCorner(), pts, v);

        //Put coord in a vector
        for (index_t i = 0; i != pts.cols(); ++i)
            vectorGaussPoints.col(counter++) = pts.col(i);
    }

    return vectorGaussPoints;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct criterionErrorAndPosAtGauss {
    criterionErrorAndPosAtGauss(const gsTHBSplineBasis < 2 > & THBSplineBasis,
                                const gsFunction < > & fApproximate, gsGeometry < > & geom, real_t threshold,
                                const gsGaussRule < > & gr):
            myGeo(geom), f_approximate(fApproximate), thb_basis(THBSplineBasis), gauss_rule(gr), error_threshold(threshold)

    {

    }

    bool operator()(const gsDomainIterator < > & cell) const {
        gsVector < > v;
        gsMatrix < > pts, vectorGaussPoints(2, gauss_rule.numNodes());
        gauss_rule.mapTo(cell.lowerCorner(), cell.upperCorner(), pts, v);
        gsMatrix < > valueGaussPointsFitting(2, vectorGaussPoints.cols()); //Matrix that contains the value of the fitted function in the Gauss points of A SINGLE ELEMENT
        gsMatrix < > valueGaussPoints(2, vectorGaussPoints.cols()); //Matrix that contains the value of the known function in the Gauss points of A SINGLE ELEMENT
        gsMatrix < > valueDetJFitting(1, vectorGaussPoints.cols()); //Matrix that contains the value of the fitted function in the Gauss points of A SINGLE ELEMENT
        gsMatrix < > errorsVector(1, vectorGaussPoints.cols());
        gsMatrix < > errorsVectorScaled(1, vectorGaussPoints.cols());
        for (index_t i = 0; i != pts.cols(); ++i) {
            vectorGaussPoints.col(i) = pts.col(i);
        }

        //Evaluate my fitted function in Gauss nodes of the current element
        f_approximate.eval_into(vectorGaussPoints, valueGaussPointsFitting);
        //Evaluate my known function in Gauss nodes of the current element
        gsMatrix < real_t > punto(2, 1);
        gsMatrix < real_t > qiymet(2, 1);
        gsMatrix < real_t > valore1(1, 1);
        gsMatrix < real_t > valore2(1, 1);
        for (int i = 0; i != vectorGaussPoints.cols(); i++) {
            punto(0, 0) = vectorGaussPoints(0, i);
            punto(1, 0) = vectorGaussPoints(1, i);
            myGeo.eval_into(punto, qiymet);
            valueGaussPoints(0, i) = qiymet(0, 0);
            valueGaussPoints(1, i) = qiymet(1, 0);
            valueDetJFitting(0, i) = (f_approximate.jacobian(punto)).determinant();

        }

        //Evaluate errors by component
        errorsVector = (valueGaussPoints - valueGaussPointsFitting).cwiseAbs();

        //If you want to drive the refinement by the error
        //        return ( (errorsVector.array() > error_threshold ).any() );

        //If you want to drive the refinement by the positivity of the determinant of the jacobian
        return ((valueDetJFitting.array() < 0).any());

    }

    // referenced data
    const gsGeometry < > & myGeo;
    const gsFunction < > & f_approximate;
    const gsTHBSplineBasis < 2 > & thb_basis;
    const gsGaussRule < > & gauss_rule;
    real_t error_threshold;

};

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//WRITE NEW TO OVERLOAD!
vector < bool > buildMarkVector(const gsTHBSplineBasis < 2 > & THBSplineBasis,
                                const gsFunction < > & fApproximate, gsGeometry < > & geom, real_t threshold) {

    gsHDomainIterator < real_t, 2 > domIter(THBSplineBasis); //Create domain iterator
    int numberOfElements = THBSplineBasis.numElements(); //Find number of elements
    gsGaussRule < > gr(THBSplineBasis, 2, 2); //This is the gauss rule
    vector < bool > markedElements;
    markedElements.reserve(numberOfElements);
    criterionErrorAndPosAtGauss myCriterion(THBSplineBasis, fApproximate, geom, threshold, gr);

    for (; domIter.good(); domIter.next()) // loop over all elements
    {

        markedElements.push_back(myCriterion(domIter)); //marking

    }

    return markedElements;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void doRefinement(gsTHBSplineBasis < 2 > & THB, vector < bool > markedElements, int extension) {
    //Do the refinement
    const int dim = 2; //equal to 2 in this case
    // numMarked: Number of marked cells on current patch (1 patch in this case), also currently marked cell
    // poffset  : offset index for the first element on a patch
    // globalCount: counter for the current global element index
    int numMarked, poffset = 0, globalCount = 0;

    // refBoxes: contains marked boxes on a given patch
    gsMatrix < double_t > refBoxes;

    // Get number of elements to be refined on this patch
    const int numEl = THB.numElements();
    numMarked = std::count_if(markedElements.begin() + poffset,
                              markedElements.begin() + poffset + numEl,
                              std::bind2nd(std::equal_to < bool > (), true));

    refBoxes.resize(dim, 2 * numMarked);

    numMarked = 0; // counting current patch element to be refined

    // for all elements in my patch
    typename gsBasis < double_t > ::domainIter domainIter = THB.makeDomainIterator();
    for (; domainIter -> good(); domainIter -> next()) {
        if (markedElements[globalCount++]) // refine this element ?
        {

            // Construct degenerate box by setting both
            // corners equal to the center
            refBoxes.col(2 * numMarked) = domainIter -> centerPoint();
            refBoxes.col(2 * numMarked + 1) = domainIter -> centerPoint();

            // Advance marked cells counter
            numMarked++;

        }
    }

    // Refine all of the found refBoxes in this patch
    THB.refine(refBoxes, extension);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initializeLocalBasis(gsTHBSplineBasis < 2 > & thbLocal, gsTHBSplineBasis < 2 > & THB, gsMatrix < int > & maxLevels, gsMatrix < int > & minLevels, gsMatrix < real_t > & matrixOfLC, gsMatrix < real_t > & matrixOfUC, int ii) {
    //It initializes the local basis for the selected basis function in the corresponding bounding box
    //ii is the basis function number

    int myMax = maxLevels.maxCoeff(); //max level of the global thb basis
    int myMin = THB.levelOf(ii);; //min level which overlap my selected basis function

    gsTensorBSplineBasis < 2, real_t > tensorBasis = * THB.getBases()[THB.levelOf(ii)];
    //knotVector: first,last,number of interior knots, multiplicity end knots
    gsKnotVector < > u_knots(tensorBasis.knot(0, THB.maxDegree() + matrixOfLC(0, ii) / (pow(2, myMax - myMin))), tensorBasis.knot(0, THB.maxDegree() + matrixOfUC(0, ii) / (pow(2, myMax - myMin))), matrixOfUC(0, ii) / (pow(2, myMax - myMin)) - matrixOfLC(0, ii) / (pow(2, myMax - myMin)) - 1, THB.maxDegree() + 1);
    gsKnotVector < > v_knots(tensorBasis.knot(1, THB.maxDegree() + matrixOfLC(1, ii) / (pow(2, myMax - myMin))), tensorBasis.knot(1, THB.maxDegree() + matrixOfUC(1, ii) / (pow(2, myMax - myMin))), matrixOfUC(1, ii) / (pow(2, myMax - myMin)) - matrixOfLC(1, ii) / (pow(2, myMax - myMin)) - 1, THB.maxDegree() + 1);
    // Create a basis and apply initial uniform refinement
    gsTensorBSplineBasis < 2 > T_tbasis(u_knots, v_knots);

    //I have initialized my sub basis for the selected basis function. At the moment it's a tensor product basis, its level is always zero.
    //Now I have to map the indices from global to local:
    //Get knot vectors from tensorBasis of level minLevels(0,ii) and compute vectors of indices in that level
    gsKnotVector < > uglobal = tensorBasis.knots(0);
    gsKnotVector < > vglobal = tensorBasis.knots(1);
    int counter = 0;
    vector < int > uind;
    vector < int > vind;
    for (int j = THB.maxDegree(); j != uglobal.size() - THB.maxDegree(); j++) {
        uind.push_back(counter);
        vind.push_back(counter);
        counter = counter + 1;
    }

    counter = 0;
    vector < index_t > uindmapfinal;
    vector < index_t > vindmapfinal;
    //In uindmapfinal I express all the indices of the local basis in terms of the indices of the finest level of the total THB basis
    for (int j = 0; j != matrixOfUC(0, ii) - matrixOfLC(0, ii) + 1; j++) {
        uindmapfinal.push_back(matrixOfLC(0, ii) + counter);
        counter = counter + 1;
    }
    counter = 0;
    for (int j = 0; j != matrixOfUC(1, ii) - matrixOfLC(1, ii) + 1; j++) {
        vindmapfinal.push_back(matrixOfLC(1, ii) + counter);
        counter = counter + 1;
    }

    //get the boxes
    gsMatrix < index_t > bLC;
    gsMatrix < index_t > bUC;
    gsVector < index_t > myLevel;
    (THB.tree().getBoxes(bLC, bUC, myLevel));

    gsMatrix < int > box2Insert(2, 2);
    gsMatrix < int > box2InsertMap(2, 2);
    index_t level2Insert;
    bool intersection;
    std::vector < index_t > boxIn;
    for (int i = 0; i != bLC.rows(); i++) {

        //check if the box and the local basis intersect
        if (uindmapfinal[0] < bUC(i, 0) && uindmapfinal[uindmapfinal.size() - 1] > bLC(i, 0) &&
            vindmapfinal[0] < bUC(i, 1) && vindmapfinal[vindmapfinal.size() - 1] > bLC(i, 1)) {
            intersection = true;
        } else {
            intersection = false;
        }

        //if the selected box and my local basis intersect, insert the box.
        //box2Insert is expressed in terms of the indices of the global basis
        if (intersection == true) {

            level2Insert = myLevel(i) - myMin;
            box2Insert(0, 0) = std::max(uindmapfinal[0], bLC(i, 0)) / pow(2, myMax - myLevel(i, 0));
            box2Insert(0, 1) = std::max(vindmapfinal[0], bLC(i, 1)) / pow(2, myMax - myLevel(i, 0));
            box2Insert(1, 0) = std::min(uindmapfinal[uindmapfinal.size() - 1], bUC(i, 0)) /
                               pow(2, myMax - myLevel(i, 0));
            box2Insert(1, 1) = std::min(vindmapfinal[vindmapfinal.size() - 1], bUC(i, 1)) /
                               pow(2, myMax - myLevel(i, 0));

            //map indices global local
            counter = 0;
            //DON'T TOUCH THESE TWO VECTORS!!!!
            vector < int > uindmaplevel;
            vector < int > vindmaplevel;
            for (int j = 0; j != matrixOfUC(0, ii) / pow(2, myMax - myLevel(i, 0)) - matrixOfLC(0, ii) / pow(2, myMax - myLevel(i, 0)) + 1; j++) {
                uindmaplevel.push_back(matrixOfLC(0, ii) / pow(2, myMax - myLevel(i, 0)) + counter);
                counter = counter + 1;
            }
            counter = 0;
            for (int j = 0; j != matrixOfUC(1, ii) / pow(2, myMax - myLevel(i, 0)) - matrixOfLC(1, ii) / pow(2, myMax - myLevel(i, 0)) + 1; j++) {
                vindmaplevel.push_back(matrixOfLC(1, ii) / pow(2, myMax - myLevel(i, 0)) + counter);
                counter = counter + 1;
            }

            //box2InsertMap is the final box that I have to insert expressed with the indices of the local basis
            for (int r = 0; r != uindmaplevel.size(); r++) {
                if (uindmaplevel[r] == box2Insert(0, 0)) {
                    box2InsertMap(0, 0) = r;
                }
                if (uindmaplevel[r] == box2Insert(1, 0)) {
                    box2InsertMap(1, 0) = r;
                }

            }
            for (int r = 0; r != vindmaplevel.size(); r++) {
                if (vindmaplevel[r] == box2Insert(0, 1)) {
                    box2InsertMap(0, 1) = r;
                }
                if (vindmaplevel[r] == box2Insert(1, 1)) {
                    box2InsertMap(1, 1) = r;
                }

            }

            boxIn.push_back(level2Insert);
            boxIn.push_back(box2InsertMap(0, 0));
            boxIn.push_back(box2InsertMap(0, 1));
            boxIn.push_back(box2InsertMap(1, 0));
            boxIn.push_back(box2InsertMap(1, 1));
            thbLocal.refineElements(boxIn);
            boxIn.clear();

        }

    }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

gsGeometry < > ::uPtr myLocalApproximate(gsTHBSplineBasis < 2 > & THB,
                                         const gsFunctionExpr < > & f1,
                                         const gsFunctionExpr < > & f2, real_t & error, gsMatrix < real_t > & matrixOfLC, gsMatrix < real_t > & matrixOfUC, gsMatrix < int > & maxLevels, gsMatrix < int > & minLevels,
                                         real_t & au, real_t & bu, real_t & av, real_t & bv, real_t & boundaryerror, gsMatrix < int > basisInd, gsMatrix < > THBcoefs) {
    int numberOfBoundaryPoints = 100;
    real_t lowerror, uperror, lefterror, righterror;
    gsMatrix < > vectorGaussPoints = computeGaussPoints(THB); //It contains the coordinates of ALL the Gauss points
    gsMatrix < > lowboundaryPoints(2, numberOfBoundaryPoints + 1);
    gsMatrix < > upboundaryPoints(2, numberOfBoundaryPoints + 1);
    gsMatrix < > leftboundaryPoints(2, numberOfBoundaryPoints + 1);
    gsMatrix < > rightboundaryPoints(2, numberOfBoundaryPoints + 1);
    gsMatrix < > lowVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < > upVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < > leftVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < > rightVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < real_t > lowfittedVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < real_t > upfittedVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < real_t > leftfittedVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < real_t > rightfittedVal(2, numberOfBoundaryPoints + 1);
    for (int i = 0; i < numberOfBoundaryPoints + 1; i++) {
        lowboundaryPoints(0, i) = au + (bu - au) * ((double) i / numberOfBoundaryPoints);
        lowboundaryPoints(1, i) = av;
        upboundaryPoints(0, i) = au + (bu - au) * ((double) i / numberOfBoundaryPoints);
        upboundaryPoints(1, i) = bv;
        leftboundaryPoints(0, i) = au;
        leftboundaryPoints(1, i) = av + (bv - av) * ((double) i / numberOfBoundaryPoints);
        rightboundaryPoints(0, i) = bu;
        rightboundaryPoints(1, i) = av + (bv - av) * ((double) i / numberOfBoundaryPoints);
    }

    gsHDomainIterator < real_t, 2 > domIter1(THB); //Create domain iterator
    int numberOfElements = THB.numElements(); //Find number of elements

    gsGaussRule < > gr(THB, 2, 2);

    gsMatrix < > valueGaussPoints(2, gr.numNodes() * numberOfElements); //Vector for the values of my know function in Gauss points
    gsMatrix < real_t > valueGaussPointsFitting(2, gr.numNodes() * numberOfElements); //Vector for the value of the fitted function in Gauss points

    gsVector < > v; //For the gauss nodes
    vector < double > vTot; //For the global gauss nodes

    gsMatrix < > pts(2, gr.numNodes() * numberOfElements); //Number of points

    gsMatrix < real_t > unknowns(THB.size(), 2); //it will contain the solution
    unknowns.setZero(THB.size(), 2);

    gsMatrix < real_t > punto(2, 1);
    gsMatrix < real_t > valore1(1, 1);
    gsMatrix < real_t > valore2(1, 1);
    //Evaluate my function in every Gauss point and put the values in a vector (useless loop)
    for (int i = 0; i != vectorGaussPoints.cols(); i++) {
        punto(0, 0) = vectorGaussPoints(0, i);
        punto(1, 0) = vectorGaussPoints(1, i);
        (f1).eval_into(punto, valore1);
        (f2).eval_into(punto, valore2);
        valueGaussPoints(0, i) = valore1(0, 0);
        valueGaussPoints(1, i) = valore2(0, 0);

    }

    int myMax, myMin, numberOfElementsLocal;
    gsTensorBSplineBasis < 2, real_t > tensorBasis;
    gsMatrix < > LC(2, 1);
    gsMatrix < > UC(2, 1);
    gsMatrix < > vectorGaussPointsLocal;
    real_t lambda = 0; //1e-07;//lambda parameter for fitting if needed
    vector < real_t > vectorGlobalErrorsLocal;
    int functionLevelOfii, counter, functionLevelOft;
    gsTensorBSplineBasis < 2 > tensorBasisii;
    gsMatrix < int > functionSupport(2, 2);
    gsMatrix < real_t > anchorsGlobal;
    gsMatrix < real_t > anchorsLocal;
    gsTensorBSplineBasis < 2 > tensorBasisLocal;
    gsMatrix < int > functionSupportLocal(2, 2);
    int xCoeff;
    //loop over each basis function
    for (int ii = 0; ii != matrixOfLC.cols(); ii++) {
        //UPDATE IS HERE

        myMax = maxLevels.maxCoeff(); //max level of the global thb basis
        myMin = THB.levelOf(ii);

        tensorBasis = * THB.getBases()[myMin];
        //knotVector: first,last,number of interior knots, multiplicity end knots
        gsKnotVector < > u_knots(tensorBasis.knot(0, THB.maxDegree() + matrixOfLC(0, ii) / (pow(2, myMax - myMin))), tensorBasis.knot(0, THB.maxDegree() + matrixOfUC(0, ii) / (pow(2, myMax - myMin))), matrixOfUC(0, ii) / (pow(2, myMax - myMin)) - matrixOfLC(0, ii) / (pow(2, myMax - myMin)) - 1, THB.maxDegree() + 1);
        gsKnotVector < > v_knots(tensorBasis.knot(1, THB.maxDegree() + matrixOfLC(1, ii) / (pow(2, myMax - myMin))), tensorBasis.knot(1, THB.maxDegree() + matrixOfUC(1, ii) / (pow(2, myMax - myMin))), matrixOfUC(1, ii) / (pow(2, myMax - myMin)) - matrixOfLC(1, ii) / (pow(2, myMax - myMin)) - 1, THB.maxDegree() + 1);

        // Create a basis and apply initial uniform refinement
        gsTensorBSplineBasis < 2 > T_tbasis(u_knots, v_knots);
        // Create Initial local hierarchical basis
        gsTHBSplineBasis < 2 > thbLocal(T_tbasis);

        initializeLocalBasis(thbLocal, THB, maxLevels, minLevels, matrixOfLC, matrixOfUC, ii);

        //        std::string exitBasis = "thbLocal";
        //        gsWriteParaview(thbLocal, exitBasis);
        //        gsFileManager::open("thbLocal.pvd");

        //At this point I have myLocalTHB basis.
        //I have to perform local fitting. I am performing local fitting for function ii. I have to associate
        //function ii (in hierarchical index) of the global basis with the corresponding index in the local basis.
        LC(0, 0) = matrixOfLC(0, ii);
        LC(1, 0) = matrixOfLC(1, ii);
        UC(0, 0) = matrixOfUC(0, ii);
        UC(1, 0) = matrixOfUC(1, ii);

        vectorGaussPointsLocal = computeGaussPoints(thbLocal); //It contains the coordinates of ALL the Gauss points for the local basis

        gsHDomainIterator < real_t, 2 > domIterLocal(thbLocal); //Create domain iterator
        numberOfElementsLocal = thbLocal.numElements(); //Find number of elements of the local basis
        gsGaussRule < > grLocal(thbLocal, 2, 2);

        gsMatrix < > valueGaussPointsLocal(2, grLocal.numNodes() * numberOfElementsLocal); //Vector for the values of my know function in Gauss points in bbox
        gsMatrix < real_t > valueGaussPointsFitting(2, grLocal.numNodes() * numberOfElementsLocal); //Vector for the value of the fitted function in Gauss points in bbox

        gsMatrix < > pts(2, grLocal.numNodes() * numberOfElementsLocal); //Number of gasuss points in bbox

        //Evaluate my function in every Gauss point and put the values in a vector
        for (int i = 0; i != vectorGaussPointsLocal.cols(); i++) {

            punto(0, 0) = vectorGaussPointsLocal(0, i);
            punto(1, 0) = vectorGaussPointsLocal(1, i);

            (f1).eval_into(punto, valore1);
            (f2).eval_into(punto, valore2);
            valueGaussPointsLocal(0, i) = valore1(0, 0);
            valueGaussPointsLocal(1, i) = valore2(0, 0);

        }

        //        gsDebugVar(vectorGaussPointsLocal);

        //I have the points and my known function. Now  I have to approximate my function performing fitting
        gsFitting < > fittingObject(vectorGaussPointsLocal, valueGaussPointsLocal, thbLocal);

        //To use LU and therefore direct solver
        gsSparseMatrix < real_t > matA(thbLocal.size(), thbLocal.size());
        gsMatrix < real_t > matAfull(thbLocal.size(), thbLocal.size());
        gsMatrix < real_t > vectB(thbLocal.size(), 2);
        matA.setZero();
        matAfull.setZero();
        vectB.setZero();
        fittingObject.assembleSystem(matA, vectB);
        matAfull = Eigen::MatrixXd(matA);
        gsMatrix < real_t > vectSol(thbLocal.size(), 2);
        vectSol = matAfull.partialPivLu().solve(vectB);

        //Determine the coerrect coefficient for the selected basis function
        functionLevelOfii = THB.levelOf(ii);
        tensorBasisii = * THB.getBases()[functionLevelOfii];

        //put in functionsupport the support of the function ii expressed in indices of the total thb basis with rispect to the level of the selected basis function
        functionSupport(0, 0) = (THB.tensorLevel(functionLevelOfii).elementSupport(THB.flatTensorIndexOf(ii)))(0, 0);
        functionSupport(0, 1) = (THB.tensorLevel(functionLevelOfii).elementSupport(THB.flatTensorIndexOf(ii)))(1, 0);
        functionSupport(1, 0) = (THB.tensorLevel(functionLevelOfii).elementSupport(THB.flatTensorIndexOf(ii)))(0, 1);
        functionSupport(1, 1) = (THB.tensorLevel(functionLevelOfii).elementSupport(THB.flatTensorIndexOf(ii)))(1, 1);

        counter = 0;
        vector < int > uindmaplevel;
        vector < int > vindmaplevel;
        for (int j = 0; j != matrixOfUC(0, ii) / pow(2, myMax - functionLevelOfii) - matrixOfLC(0, ii) / pow(2, myMax - functionLevelOfii) + 1; j++) {
            uindmaplevel.push_back(matrixOfLC(0, ii) / pow(2, myMax - functionLevelOfii) + counter);
            counter = counter + 1;
        }
        counter = 0;
        for (int j = 0; j != matrixOfUC(1, ii) / pow(2, myMax - functionLevelOfii) - matrixOfLC(1, ii) / pow(2, myMax - functionLevelOfii) + 1; j++) {
            vindmaplevel.push_back(matrixOfLC(1, ii) / pow(2, myMax - functionLevelOfii) + counter);
            counter = counter + 1;
        }

        anchorsGlobal = THB.anchors();
        anchorsLocal = thbLocal.anchors();
        //        gsDebugVar(anchorsGlobal);
        //        gsDebugVar(anchorsLocal);
        //if (ii==28){
        //    gsInfo<<"ciao"<<std::endl;
        //
        //}

        for (int t = 0; t != thbLocal.size(); t++) {
            functionLevelOft = thbLocal.levelOf(t);

            tensorBasisLocal = * thbLocal.getBases()[functionLevelOft];

            //put in functionsupportlocal the support of the function t expressed in indices of the local thb basis with rispect to the level of the selected basis function
            functionSupportLocal(0, 0) = (thbLocal.tensorLevel(functionLevelOft).elementSupport(thbLocal.flatTensorIndexOf(t)))(0, 0);
            functionSupportLocal(0, 1) = (thbLocal.tensorLevel(functionLevelOft).elementSupport(thbLocal.flatTensorIndexOf(t)))(1, 0);
            functionSupportLocal(1, 0) = (thbLocal.tensorLevel(functionLevelOft).elementSupport(thbLocal.flatTensorIndexOf(t)))(0, 1);
            functionSupportLocal(1, 1) = (thbLocal.tensorLevel(functionLevelOft).elementSupport(thbLocal.flatTensorIndexOf(t)))(1, 1);

            //Comparison between 2 double might introduce a bug
            //            if(anchorsGlobal(0,ii)==anchorsLocal(0,t) && anchorsGlobal(1,ii)==anchorsLocal(1,t) && functionSupport(0,0)==uindmaplevel[functionSupportLocal(0,0)] && functionSupport(1,0)==uindmaplevel[functionSupportLocal(1,0)]&& functionSupport(0,1)==vindmaplevel[functionSupportLocal(0,1)] && functionSupport(1,1)==vindmaplevel[functionSupportLocal(1,1)]){
            // Comparison between 2 double with a threshold solves the problem (but the code is not robust). This is the current implemented solution
            // Correct solution would be to compare only integers (i.e compare the upper and lower corner in the local and global basis and compare the levels)
            if (abs(anchorsGlobal(0, ii) - anchorsLocal(0, t)) < 10e-12 && abs(anchorsGlobal(1, ii) - anchorsLocal(1, t)) < 10e-12 && functionSupport(0, 0) == uindmaplevel[functionSupportLocal(0, 0)] && functionSupport(1, 0) == uindmaplevel[functionSupportLocal(1, 0)] && functionSupport(0, 1) == vindmaplevel[functionSupportLocal(0, 1)] && functionSupport(1, 1) == vindmaplevel[functionSupportLocal(1, 1)]) {

                //then t is the selected basis function
                xCoeff = t;
                break;

            }

        }

        //If you use directsolver
        unknowns(ii, 0) = vectSol(xCoeff, 0);
        unknowns(ii, 1) = vectSol(xCoeff, 1);

    }
    for (int i = 0; i < numberOfBoundaryPoints + 1; i++) {
        punto(0, 0) = lowboundaryPoints(0, i);
        punto(1, 0) = lowboundaryPoints(1, i);
        (f1).eval_into(punto, valore1);
        (f2).eval_into(punto, valore2);
        lowVal(0, i) = valore1(0, 0);
        lowVal(1, i) = valore2(0, 0);
        //        gsInfo << "i = " << i << "\n";
    }

    for (int i = 0; i < numberOfBoundaryPoints + 1; i++) {
        punto(0, 0) = upboundaryPoints(0, i);
        punto(1, 0) = upboundaryPoints(1, i);
        (f1).eval_into(punto, valore1);
        (f2).eval_into(punto, valore2);
        upVal(0, i) = valore1(0, 0);
        upVal(1, i) = valore2(0, 0);
    }
    for (int i = 0; i < numberOfBoundaryPoints + 1; i++) {
        punto(0, 0) = leftboundaryPoints(0, i);
        punto(1, 0) = leftboundaryPoints(1, i);
        (f1).eval_into(punto, valore1);
        (f2).eval_into(punto, valore2);
        leftVal(0, i) = valore1(0, 0);
        leftVal(1, i) = valore2(0, 0);
    }
    for (int i = 0; i < numberOfBoundaryPoints + 1; i++) {
        punto(0, 0) = rightboundaryPoints(0, i);
        punto(1, 0) = rightboundaryPoints(1, i);
        (f1).eval_into(punto, valore1);
        (f2).eval_into(punto, valore2);
        rightVal(0, i) = valore1(0, 0);
        rightVal(1, i) = valore2(0, 0);
    }
    //Now I can create the geometry and evaluate the error
    gsGeometry < > ::uPtr fitResultGeo = THB.makeGeometry(unknowns);
    fitResultGeo -> eval_into(vectorGaussPoints, valueGaussPointsFitting);
    fitResultGeo -> eval_into(lowboundaryPoints, lowfittedVal);
    fitResultGeo -> eval_into(upboundaryPoints, upfittedVal);
    fitResultGeo -> eval_into(leftboundaryPoints, leftfittedVal);
    fitResultGeo -> eval_into(rightboundaryPoints, rightfittedVal);

    gsMatrix < > errorsVector;
    gsMatrix < > lowerrorsVector;
    gsMatrix < > uperrorsVector;
    gsMatrix < > lefterrorsVector;
    gsMatrix < > righterrorsVector;
    errorsVector = (valueGaussPoints - valueGaussPointsFitting).cwiseAbs();
    lowerrorsVector = (lowVal - lowfittedVal).cwiseAbs();
    uperrorsVector = (upVal - upfittedVal).cwiseAbs();
    lefterrorsVector = (leftVal - leftfittedVal).cwiseAbs();
    righterrorsVector = (rightVal - rightfittedVal).cwiseAbs();
    error = errorsVector.maxCoeff();
    lowerror = lowerrorsVector.maxCoeff();
    uperror = uperrorsVector.maxCoeff();
    lefterror = lefterrorsVector.maxCoeff();
    righterror = righterrorsVector.maxCoeff();
    //    for(int i = 0; i < numberOfBoundaryPoints; i++){
    ////        gsInfo << lowboundaryPoints(0,i) << "; " << lowboundaryPoints(1,i) << "\n";
    //    gsInfo << "i = " << i << " " << upboundaryPoints(0,i) << " " << upboundaryPoints(1,i) << " " << uperrorsVector(0,i) << " "  << upVal(0,i) << " " << upfittedVal(0,i) << " "  << uperrorsVector(1,i)<< " " << upVal(1,i) << " " << upfittedVal(1,i) << "\n";
    //    }
    boundaryerror = std::max({
                                     lowerror,
                                     uperror,
                                     lefterror,
                                     righterror
                             });
    ////Max componentwise error
    //cout<<"The max error is "<<error<<endl;

    return fitResultGeo -> clone();

}
gsGeometry < > ::uPtr myLocalApproximate(gsTHBSplineBasis < 2 > & THB, gsGeometry < > & geom, real_t & error, gsMatrix < real_t > & matrixOfLC, gsMatrix < real_t > & matrixOfUC, gsMatrix < int > & maxLevels, gsMatrix < int > & minLevels,
                                         real_t & au, real_t & bu, real_t & av, real_t & bv, real_t & boundaryerror, gsMatrix < int > basisInd, gsMatrix < > THBcoefs) {
    int numberOfBoundaryPoints = 100;
    real_t lowerror, uperror, lefterror, righterror;
    gsMatrix < real_t > unknowns(THB.size(), 2); //it will contain the solution
    gsMatrix < > lowboundaryPoints(2, numberOfBoundaryPoints + 1);
    gsMatrix < > upboundaryPoints(2, numberOfBoundaryPoints + 1);
    gsMatrix < > leftboundaryPoints(2, numberOfBoundaryPoints + 1);
    gsMatrix < > rightboundaryPoints(2, numberOfBoundaryPoints + 1);
    gsMatrix < > lowVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < > upVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < > leftVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < > rightVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < real_t > lowfittedVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < real_t > upfittedVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < real_t > leftfittedVal(2, numberOfBoundaryPoints + 1);
    gsMatrix < real_t > rightfittedVal(2, numberOfBoundaryPoints + 1);
    for (int i = 0; i < numberOfBoundaryPoints + 1; i++) {
        lowboundaryPoints(0, i) = au + (bu - au) * ((double) i / numberOfBoundaryPoints);
        lowboundaryPoints(1, i) = av;
        //        gsInfo << lowboundaryPoints(0,i) << "," << lowboundaryPoints(1,i) << "\n";
        upboundaryPoints(0, i) = au + (bu - au) * ((double) i / numberOfBoundaryPoints);
        upboundaryPoints(1, i) = bv;
        leftboundaryPoints(0, i) = au;
        leftboundaryPoints(1, i) = av + (bv - av) * ((double) i / numberOfBoundaryPoints);
        rightboundaryPoints(0, i) = bu;
        rightboundaryPoints(1, i) = av + (bv - av) * ((double) i / numberOfBoundaryPoints);
    }
    gsMatrix < > vectorGaussPoints = computeGaussPoints(THB); //It contains the coordinates of ALL the Gauss points

    gsHDomainIterator < real_t, 2 > domIter1(THB); //Create domain iterator
    int numberOfElements = THB.numElements(); //Find number of elements

    gsGaussRule < > gr(THB, 2, 2);

    gsMatrix < > valueGaussPoints(2, gr.numNodes() * numberOfElements); //Vector for the values of my know function in Gauss points
    gsMatrix < real_t > valueGaussPointsFitting(2, gr.numNodes() * numberOfElements); //Vector for the value of the fitted function in Gauss points

    gsVector < > v; //For the gauss nodes
    vector < double > vTot; //For the global gauss nodes

    gsMatrix < > pts(2, gr.numNodes() * numberOfElements); //Number of points


    unknowns.setZero(THB.size(), 2);

    gsMatrix < real_t > punto(2, 1);
    gsMatrix < real_t > qiymet(2, 1);
    gsMatrix < real_t > valore1(1, 1);
    gsMatrix < real_t > valore2(1, 1);
    //Evaluate my function in every Gauss point and put the values in a vector (useless loop)
    for (int i = 0; i != vectorGaussPoints.cols(); i++) {
        punto(0, 0) = vectorGaussPoints(0, i);
        punto(1, 0) = vectorGaussPoints(1, i);
        geom.eval_into(punto, qiymet);
        //        gsInfo << qiymet << "   qiymet!\n";
        //        (f2).eval_into(punto,valore2);
        valueGaussPoints(0, i) = qiymet(0, 0);
        valueGaussPoints(1, i) = qiymet(1, 0);

    }

    int myMax, myMin, numberOfElementsLocal;
    gsTensorBSplineBasis < 2, real_t > tensorBasis;
    gsMatrix < > LC(2, 1);
    gsMatrix < > UC(2, 1);
    gsMatrix < > vectorGaussPointsLocal;
    real_t lambda = 0; //1e-07;//lambda parameter for fitting if needed
    vector < real_t > vectorGlobalErrorsLocal;
    int functionLevelOfii, counter, functionLevelOft;
    gsTensorBSplineBasis < 2 > tensorBasisii;
    gsMatrix < int > functionSupport(2, 2);
    gsMatrix < real_t > anchorsGlobal;
    gsMatrix < real_t > anchorsLocal;
    gsTensorBSplineBasis < 2 > tensorBasisLocal;
    gsMatrix < int > functionSupportLocal(2, 2);
    int xCoeff;

    //loop over each basis function
    for (int ii = 0; ii != matrixOfLC.cols(); ii++) {
//        gsInfo << ii << "\n";
        if (basisInd(ii)) {
            //            if(ii <= 8 || ii >= 0){
            //                gsInfo << "i = " << ii <<" ",\n" << THB.support(ii) << "\n";
            //        }
            //if(1){
            //        {gsInfo << "Iteration "<< ii <<", I skip\n";
            //            continue;}
            //        gsInfo << "Iteration "<< ii <<", hi\n";

            myMax = maxLevels.maxCoeff(); //max level of the global thb basis
            myMin = THB.levelOf(ii);

            tensorBasis = * THB.getBases()[myMin];
            //knotVector: first,last,number of interior knots, multiplicity end knots
            gsKnotVector < > u_knots(tensorBasis.knot(0, THB.maxDegree() + matrixOfLC(0, ii) / (pow(2, myMax - myMin))),
                                     tensorBasis.knot(0, THB.maxDegree() + matrixOfUC(0, ii) / (pow(2, myMax - myMin))),
                                     matrixOfUC(0, ii) / (pow(2, myMax - myMin)) -
                                     matrixOfLC(0, ii) / (pow(2, myMax - myMin)) - 1, THB.maxDegree() + 1);
            gsKnotVector < > v_knots(tensorBasis.knot(1, THB.maxDegree() + matrixOfLC(1, ii) / (pow(2, myMax - myMin))),
                                     tensorBasis.knot(1, THB.maxDegree() + matrixOfUC(1, ii) / (pow(2, myMax - myMin))),
                                     matrixOfUC(1, ii) / (pow(2, myMax - myMin)) -
                                     matrixOfLC(1, ii) / (pow(2, myMax - myMin)) - 1, THB.maxDegree() + 1);

            // Create a basis and apply initial uniform refinement
            gsTensorBSplineBasis < 2 > T_tbasis(u_knots, v_knots);
            // Create Initial local hierarchical basis
            gsTHBSplineBasis < 2 > thbLocal(T_tbasis);

            initializeLocalBasis(thbLocal, THB, maxLevels, minLevels, matrixOfLC, matrixOfUC, ii);

            std::string exitBasis = "thbLocal";
            gsWriteParaview(thbLocal, exitBasis);
            //        gsFileManager::open("thbLocal.pvd");

            //At this point I have myLocalTHB basis.
            //I have to perform local fitting. I am performing local fitting for function ii. I have to associate
            //function ii (in hierarchical index) of the global basis with the corresponding index in the local basis.
            LC(0, 0) = matrixOfLC(0, ii);
            LC(1, 0) = matrixOfLC(1, ii);
            UC(0, 0) = matrixOfUC(0, ii);
            UC(1, 0) = matrixOfUC(1, ii);

            vectorGaussPointsLocal = computeGaussPoints(
                    thbLocal); //It contains the coordinates of ALL the Gauss points for the local basis

            gsHDomainIterator < real_t, 2 > domIterLocal(thbLocal); //Create domain iterator
            numberOfElementsLocal = thbLocal.numElements(); //Find number of elements of the local basis
            gsGaussRule < > grLocal(thbLocal, 2, 2);

            gsMatrix < > valueGaussPointsLocal(2, grLocal.numNodes() *
                                                  numberOfElementsLocal); //Vector for the values of my know function in Gauss points in bbox
            gsMatrix < real_t > valueGaussPointsFitting(2, grLocal.numNodes() *
                                                           numberOfElementsLocal); //Vector for the value of the fitted function in Gauss points in bbox

            gsMatrix < > pts(2, grLocal.numNodes() * numberOfElementsLocal); //Number of gasuss points in bbox

            //Evaluate my function in every Gauss point and put the values in a vector
            for (int i = 0; i != vectorGaussPointsLocal.cols(); i++) {

                punto(0, 0) = vectorGaussPointsLocal(0, i);
                punto(1, 0) = vectorGaussPointsLocal(1, i);

                geom.eval_into(punto, qiymet);
                valueGaussPointsLocal(0, i) = qiymet(0, 0);
                valueGaussPointsLocal(1, i) = qiymet(1, 0);

            }

            //        gsDebugVar(vectorGaussPointsLocal);

            //I have the points and my known function. Now  I have to approximate my function performing fitting
            gsFitting < > fittingObject(vectorGaussPointsLocal, valueGaussPointsLocal, thbLocal);

            //To use LU and therefore direct solver
            gsSparseMatrix < real_t > matA(thbLocal.size(), thbLocal.size());
            gsMatrix < real_t > matAfull(thbLocal.size(), thbLocal.size());
            gsMatrix < real_t > vectB(thbLocal.size(), 2);
            matA.setZero();
            matAfull.setZero();
            vectB.setZero();
            fittingObject.assembleSystem(matA, vectB);
            matAfull = Eigen::MatrixXd(matA);
            gsMatrix < real_t > vectSol(thbLocal.size(), 2);
            vectSol = matAfull.partialPivLu().solve(vectB);

            //Determine the coerrect coefficient for the selected basis function
            functionLevelOfii = THB.levelOf(ii);
            tensorBasisii = * THB.getBases()[functionLevelOfii];

            //put in functionsupport the support of the function ii expressed in indices of the total thb basis with rispect to the level of the selected basis function
            functionSupport(0, 0) = (THB.tensorLevel(functionLevelOfii).elementSupport(THB.flatTensorIndexOf(ii)))(0, 0);
            functionSupport(0, 1) = (THB.tensorLevel(functionLevelOfii).elementSupport(THB.flatTensorIndexOf(ii)))(1, 0);
            functionSupport(1, 0) = (THB.tensorLevel(functionLevelOfii).elementSupport(THB.flatTensorIndexOf(ii)))(0, 1);
            functionSupport(1, 1) = (THB.tensorLevel(functionLevelOfii).elementSupport(THB.flatTensorIndexOf(ii)))(1, 1);

            counter = 0;
            vector < int > uindmaplevel;
            vector < int > vindmaplevel;
            for (int j = 0; j != matrixOfUC(0, ii) / pow(2, myMax - functionLevelOfii) -
                                 matrixOfLC(0, ii) / pow(2, myMax - functionLevelOfii) + 1; j++) {
                uindmaplevel.push_back(matrixOfLC(0, ii) / pow(2, myMax - functionLevelOfii) + counter);
                counter = counter + 1;
            }
            counter = 0;
            for (int j = 0; j != matrixOfUC(1, ii) / pow(2, myMax - functionLevelOfii) -
                                 matrixOfLC(1, ii) / pow(2, myMax - functionLevelOfii) + 1; j++) {
                vindmaplevel.push_back(matrixOfLC(1, ii) / pow(2, myMax - functionLevelOfii) + counter);
                counter = counter + 1;
            }

            anchorsGlobal = THB.anchors();
            anchorsLocal = thbLocal.anchors();
            //        gsDebugVar(anchorsGlobal);
            //        gsDebugVar(anchorsLocal);
            //if (ii==28){
            //    gsInfo<<"ciao"<<std::endl;
            //
            //}

            for (int t = 0; t != thbLocal.size(); t++) {
                functionLevelOft = thbLocal.levelOf(t);

                tensorBasisLocal = * thbLocal.getBases()[functionLevelOft];

                //put in functionsupportlocal the support of the function t expressed in indices of the local thb basis with rispect to the level of the selected basis function
                functionSupportLocal(0, 0) = (thbLocal.tensorLevel(functionLevelOft).elementSupport(
                        thbLocal.flatTensorIndexOf(t)))(0, 0);
                functionSupportLocal(0, 1) = (thbLocal.tensorLevel(functionLevelOft).elementSupport(
                        thbLocal.flatTensorIndexOf(t)))(1, 0);
                functionSupportLocal(1, 0) = (thbLocal.tensorLevel(functionLevelOft).elementSupport(
                        thbLocal.flatTensorIndexOf(t)))(0, 1);
                functionSupportLocal(1, 1) = (thbLocal.tensorLevel(functionLevelOft).elementSupport(
                        thbLocal.flatTensorIndexOf(t)))(1, 1);

                //Comparison between 2 double might introduce a bug
                //            if(anchorsGlobal(0,ii)==anchorsLocal(0,t) && anchorsGlobal(1,ii)==anchorsLocal(1,t) && functionSupport(0,0)==uindmaplevel[functionSupportLocal(0,0)] && functionSupport(1,0)==uindmaplevel[functionSupportLocal(1,0)]&& functionSupport(0,1)==vindmaplevel[functionSupportLocal(0,1)] && functionSupport(1,1)==vindmaplevel[functionSupportLocal(1,1)]){
                // Comparison between 2 double with a threshold solves the problem (but the code is not robust). This is the current implemented solution
                // Correct solution would be to compare only integers (i.e compare the upper and lower corner in the local and global basis and compare the levels)
                if (abs(anchorsGlobal(0, ii) - anchorsLocal(0, t)) < 10e-12 &&
                    abs(anchorsGlobal(1, ii) - anchorsLocal(1, t)) < 10e-12 &&
                    functionSupport(0, 0) == uindmaplevel[functionSupportLocal(0, 0)] &&
                    functionSupport(1, 0) == uindmaplevel[functionSupportLocal(1, 0)] &&
                    functionSupport(0, 1) == vindmaplevel[functionSupportLocal(0, 1)] &&
                    functionSupport(1, 1) == vindmaplevel[functionSupportLocal(1, 1)]) {

                    //then t is the selected basis function
                    xCoeff = t;
                    break;

                }

            }

            //If you use directsolver
            unknowns(ii, 0) = vectSol(xCoeff, 0);
            unknowns(ii, 1) = vectSol(xCoeff, 1);
        }
        else {
            unknowns(ii, 0) = THBcoefs(ii, 0);
            unknowns(ii, 1) = THBcoefs(ii, 1);
        }
    }
//    gsInfo << "THE unknowns vector " << unknowns << "\n";
    for (int i = 0; i < numberOfBoundaryPoints + 1; i++) {
        punto(0, 0) = lowboundaryPoints(0, i);
        punto(1, 0) = lowboundaryPoints(1, i);
        (geom).eval_into(punto, qiymet);
        lowVal(0, i) = qiymet(0, 0);
        lowVal(1, i) = qiymet(1, 0);
        //        gsInfo << "i = " << i << " " <<  lowVal(0,i) << " " << lowVal(1,i) << "\n";
    }

    for (int i = 0; i < numberOfBoundaryPoints + 1; i++) {
        punto(0, 0) = upboundaryPoints(0, i);
        punto(1, 0) = upboundaryPoints(1, i);
        (geom).eval_into(punto, qiymet);
        //        (geom).eval_into(punto,valore2);
        upVal(0, i) = qiymet(0, 0);
        upVal(1, i) = qiymet(1, 0);
    }
    for (int i = 0; i < numberOfBoundaryPoints + 1; i++) {
        punto(0, 0) = leftboundaryPoints(0, i);
        punto(1, 0) = leftboundaryPoints(1, i);
        (geom).eval_into(punto, qiymet);
        //        (geom).eval_into(punto,valore2);
        leftVal(0, i) = qiymet(0, 0);
        leftVal(1, i) = qiymet(1, 0);
    }
    for (int i = 0; i < numberOfBoundaryPoints + 1; i++) {
        punto(0, 0) = rightboundaryPoints(0, i);
        punto(1, 0) = rightboundaryPoints(1, i);
        (geom).eval_into(punto, qiymet);
        //        (geom).eval_into(punto,valore2);
        rightVal(0, i) = qiymet(0, 0);
        rightVal(1, i) = qiymet(1, 0);
    }
    //Now I can create the geometry and evaluate the error
    gsGeometry < > ::uPtr fitResultGeo = THB.makeGeometry(unknowns);
    fitResultGeo -> eval_into(vectorGaussPoints, valueGaussPointsFitting);
    fitResultGeo -> eval_into(lowboundaryPoints, lowfittedVal);
    fitResultGeo -> eval_into(upboundaryPoints, upfittedVal);
    fitResultGeo -> eval_into(leftboundaryPoints, leftfittedVal);
    fitResultGeo -> eval_into(rightboundaryPoints, rightfittedVal);
    gsMatrix < > errorsVector;
    gsMatrix < > lowerrorsVector;
    gsMatrix < > uperrorsVector;
    gsMatrix < > lefterrorsVector;
    gsMatrix < > righterrorsVector;
    errorsVector = (valueGaussPoints - valueGaussPointsFitting).cwiseAbs();
    lowerrorsVector = (lowVal - lowfittedVal).cwiseAbs();
    uperrorsVector = (upVal - upfittedVal).cwiseAbs();
    lefterrorsVector = (leftVal - leftfittedVal).cwiseAbs();
    righterrorsVector = (rightVal - rightfittedVal).cwiseAbs();
    error = errorsVector.maxCoeff();
    lowerror = lowerrorsVector.maxCoeff();
    uperror = uperrorsVector.maxCoeff();
    lefterror = lefterrorsVector.maxCoeff();
    righterror = righterrorsVector.maxCoeff();
    //    for(int i=0;i < numberOfBoundaryPoints+1;i++){
    //        gsInfo << "i = " << i << " " << lowVal(1,i) << " " << lowfittedVal(1,i) << "\n";
    //    }
    ////Max componentwise error
    //cout
    gsInfo << "The max error is "<<error<<endl;
    gsInfo << "The boundary error is "<<boundaryerror<<endl;
    boundaryerror = std::max({
                                     lowerror,
                                     uperror,
                                     lefterror,
                                     righterror
                             });

    return fitResultGeo -> clone();

}

void checkAccuracySufficiency(gsGeometry < > & geom, gsGeometry < > & geomForCheck, real_t & error, real_t & boundaryerror, real_t & au, real_t & bu, real_t & av, real_t & bv) {
    int numberOfBoundaryPoints = 100, numCols = 19, hsize = 0;
    double coords[4 * (int) pow(2, coarseLevel) * 4 * (int) pow(2, coarseLevel)][5];
    boxToDomain(redBox, coords);

    int numberOfInnerPoints = numCols * numCols;
    gsMatrix < > innerPoints(2, numberOfInnerPoints);
    gsMatrix < > innerVal(2, numberOfInnerPoints);
    gsMatrix < > innerFittedVal(2, numberOfInnerPoints);

    gsMatrix < > punto(2, 1);
    gsMatrix < > qiymet(2, 1);

    for (int i = 0; i < numCols; i++) {
        for (int j = 0; j < numCols; j++) {
            innerPoints(0, i * numCols + j) = au + (bu - au) * ((double)(i + 1) / (numCols + 1));
            innerPoints(1, i * numCols + j) = av + (bu - au) * ((double)(j + 1) / (numCols + 1));
        }
    }
    //    gsMatrix<real_t> hpoint(2, 4*numRed*(numberOfBoundaryPoints));

    for (int i = 0; i < numRed; i++) {
        if (redBox[i][1] == 0) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                //                    hpoint(0, hsize) = au;
                //                    hpoint(1, hsize) = coords[i][ 2] + (coords[i][ 4] - coords[i][2])*((double)j/numberOfBoundaryPoints);
                hsize++;
            }
        }
        if (redBox[i][2] == 0) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                //                    hpoint(0, hsize) = coords[i][1] + (coords[i][3] - coords[i][1]) * ((double) j / numberOfBoundaryPoints);
                //                    hpoint(1, hsize) = av;
                hsize++;
            }
        }
        if (redBox[i][3] == 4 * pow(2, coarseLevel + 1)) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                //                    hpoint(0, hsize) = bu;
                //                    hpoint(1, hsize) = coords[i][2] + (coords[i][4] - coords[i][2]) * ((double) j / numberOfBoundaryPoints);
                hsize++;
            }
        }
        //            else{gsInfo << redBox[i][3] << " neq " << 4*pow(2,fineLevel);}
        if (redBox[i][4] == 4 * pow(2, coarseLevel + 1)) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                //                    hpoint(0, hsize) = coords[i][1] + (coords[i][3] - coords[i][1]) * ((double) j / numberOfBoundaryPoints);
                //                    hpoint(1, hsize) = bv;
                hsize++;
            }
        }
        //            else{gsInfo << redBox[i][4] << " neq " << 4*pow(2,fineLevel) << "\n";}
    }
    gsMatrix < real_t > hpoint(2, hsize);
    int k = 0;
    for (int i = 0; i < numRed; i++) {
        if (redBox[i][1] == 0) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                hpoint(0, k) = au;
                hpoint(1, k) = coords[i][2] + (coords[i][4] - coords[i][2]) * ((double) j / numberOfBoundaryPoints);
                k++;
            }
        }
        if (redBox[i][2] == 0) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                hpoint(0, k) = coords[i][1] + (coords[i][3] - coords[i][1]) * ((double) j / numberOfBoundaryPoints);
                hpoint(1, k) = av;
                k++;
            }
        }
        if (redBox[i][3] == 4 * pow(2, coarseLevel + 1)) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                hpoint(0, k) = bu;
                hpoint(1, k) = coords[i][2] + (coords[i][4] - coords[i][2]) * ((double) j / numberOfBoundaryPoints);
                k++;
            }
        }
        //            else { gsInfo << redBox[i][3] << " neq " << 4 * pow(2, fineLevel); }
        if (redBox[i][4] == 4 * pow(2, coarseLevel + 1)) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                //                gsInfo << coords[0][1] << ", " << coords[0][3] << "\n";
                hpoint(0, k) = coords[i][1] + (coords[i][3] - coords[i][1]) * ((double) j / numberOfBoundaryPoints);
                hpoint(1, k) = bv;
                k++;
            }
        }
    }

    //    gsInfo << hsize << "\n";
    gsMatrix < > hVal(2, hsize);
    gsMatrix < > hfittedVal(2, hsize);

    //Compute geometries
//    gsInfo << "hpoint:\n" << hpoint << "\n";
    //    gsInfo << "hsize = " << hsize << "\n";
    for (int i = 0; i < hsize; i++) {
        punto(0, 0) = hpoint(0, i);
        punto(1, 0) = hpoint(1, i);
        (geomForCheck).eval_into(punto, qiymet);
        //        (geom).eval_into(punto,valore2);
        hVal(0, i) = qiymet(0, 0);
        hVal(1, i) = qiymet(1, 0);
        (geom).eval_into(punto, qiymet);
        hfittedVal(0, i) = qiymet(0, 0);
        hfittedVal(1, i) = qiymet(1, 0);
    }

//    for (int i = 0; i < hsize; i++) {
//        outfile << "hpoint " << i << ":\t";
//        outfile << hpoint(0, i) << "\t" << hpoint(1, i) << ",\t the init.value: \t" << hVal(0, i) << "\t" << hVal(1, i) <<
//                ",\t fitted value: \t" << hfittedVal(0, i) << "\t" << hfittedVal(1, i) << "\n";
//
//    }

    for (int i = 0; i < numberOfInnerPoints; i++) {
        punto(0, 0) = innerPoints(0, i);
        punto(1, 0) = innerPoints(1, i);
        (geomForCheck).eval_into(punto, qiymet);
        //        (geom).eval_into(punto,valore2);
        innerVal(0, i) = qiymet(0, 0);
        innerVal(1, i) = qiymet(1, 0);
        (geom).eval_into(punto, qiymet);
        //        gsInfo << "i = " << i << "\n";
        innerFittedVal(0, i) = qiymet(0, 0);
        innerFittedVal(1, i) = qiymet(1, 0);
    }
    gsMatrix < > errorsVector;
    gsMatrix < > herrorsVector;
    errorsVector = (innerVal - innerFittedVal).cwiseAbs();
    herrorsVector = (hVal - hfittedVal).cwiseAbs();
//    gsInfo << "hval:\n" << hVal << "\n";
//    gsInfo << "hfittedval:\n" << hfittedVal << "\n";
//    gsInfo << herrorsVector << "\n";

    gsMatrix < > errorsAbs(1,errorsVector.cols());
    gsMatrix < > herrorsAbs(1,herrorsVector.cols());
    for (int i = 0; i < errorsVector.cols(); ++i) {
        errorsAbs(0,i) = math::sqrt(errorsVector(0,i)*errorsVector(0,i) + errorsVector(1,i)*errorsVector(1,i));
    }
    for (int i = 0; i < herrorsVector.cols(); ++i) {
        herrorsAbs(0,i) = math::sqrt(herrorsVector(0,i)*herrorsVector(0,i) + herrorsVector(1,i)*herrorsVector(1,i));
    }
    boundaryerror = herrorsAbs.maxCoeff();
    error = std::max(errorsAbs.maxCoeff(), boundaryerror);
    //    for(int i=0;i < numberOfBoundaryPoints+1;i++){
    //        gsInfo << "i = " << i << " " << lowVal(1,i) << " " << lowfittedVal(1,i) << "\n";
    //    }
    ////Max componentwise error
    cout << "The max error is " <<error << "\n";
    cout << "Epsilon f is" << boundaryerror << "\n";

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void findSelectedBasis(gsTHBSplineBasis < 2 > & THB,
                       const gsFunction < > & fittedFunction, double tolerance, std::vector < int > & selectedBasisFunctions) {
    gsHDomainIterator < real_t, 2 > domIter(THB);
    gsGaussRule < > gaussRule(THB, 2, 2);
    gsMatrix < > mappedGP(2, gaussRule.numNodes() * THB.numElements());
    gsVector < > mappedW;
    gsMatrix < > activeFunctions;
    gsMatrix < > myPoint(2, 1);
    bool isSelected;
    for (; domIter.good(); domIter.next()) {
        gaussRule.mapTo(domIter.lowerCorner(), domIter.upperCorner(), mappedGP, mappedW);

        //mark the cells
        for (int n = 0; n < mappedGP.cols(); n++) {
            myPoint(0, 0) = mappedGP(0, n);
            myPoint(1, 0) = mappedGP(1, n);
            if (fittedFunction.jacobian(myPoint).determinant() < tolerance) {
                //Then the cell is marked and we have to find the basis functions that overlap
                //selectedBasisFunctions contains the basis functions which overlap the marked cells.
                //We do not include the boundary basis functions
                for (int l = 0; l < THB.active(myPoint).rows(); l++) {
                    isSelected = true;

                    for (int h = 0; h < THB.allBoundary().rows(); h++) {
                        if (THB.active(myPoint)(l, 0) == THB.allBoundary()(h, 0)) {
                            isSelected = false;
                            break;
                        }
                    }

                    //Here I push back if the flag allows me
                    if (isSelected == true) {
                        selectedBasisFunctions.push_back(THB.active(myPoint)(l, 0));
                    }

                }

                break;

            }

        }
        //Get rid of useless duplications
        sort(selectedBasisFunctions.begin(), selectedBasisFunctions.end());
        selectedBasisFunctions.erase(unique(selectedBasisFunctions.begin(), selectedBasisFunctions.end()),
                                     selectedBasisFunctions.end());

    }

}

void local_opt(gsTHBSplineBasis < 2 > & THB, gsGeometry < > & fittedFunction, std::vector < int > & selectedBasisFunctions, int optIterations) {

    gsVector < > mappedW;
    std::vector < double > xCorr;
    std::vector < double > yCorr;
    gsMatrix < > myPoint(2, 1);
    double Csx, Csy, correctionNx, correctionNy, correctionDx, correctionDy;
    for (int iterations = 0; iterations < optIterations; iterations++) {
        gsGaussRule < > gaussRule(THB, 2, 2);
        gsMatrix < > mappedGP(2, gaussRule.numNodes() * THB.numElements());
        gsHDomainIterator < real_t, 2 > domIter1(THB);

        //Dummy vectors to save the corrections

        //loop thorough selected basis functions
        for (auto & s: selectedBasisFunctions) {

            //Save control points(I will apply the correction later)
            Csx = fittedFunction.coefs()(s, 0);
            Csy = fittedFunction.coefs()(s, 1);

            correctionNx = 0;
            correctionNy = 0;
            correctionDx = 0;
            correctionDy = 0;

            //Set the selected control points to zero(only for the loops)
            fittedFunction.coefs()(s, 0) = 0;
            fittedFunction.coefs()(s, 1) = 0;

            //Loop through the elements
            for (; domIter1.good(); domIter1.next()) {

                //If the central point of the element is inside the support of the selected bas fun perform integration
                if (domIter1.centerPoint()(0, 0) > THB.support(s)(0, 0) &&
                    domIter1.centerPoint()(0, 0) < THB.support(s)(0, 1) &&
                    domIter1.centerPoint()(1, 0) > THB.support(s)(1, 0) &&
                    domIter1.centerPoint()(1, 0) < THB.support(s)(1, 1)) {

                    gaussRule.mapTo(domIter1.lowerCorner(), domIter1.upperCorner(), mappedGP, mappedW);
                    //loop thorough gauss nodes
                    for (int n = 0; n < mappedGP.cols(); n++) {
                        myPoint(0, 0) = mappedGP(0, n);
                        myPoint(1, 0) = mappedGP(1, n);

                        correctionNx = ((fittedFunction.deriv(myPoint)(0, 0) * THB.derivSingle(s, myPoint)(0, 0) +
                                         fittedFunction.deriv(myPoint)(1, 0) * THB.derivSingle(s, myPoint)(1, 0))) *
                                       mappedW(n, 0) + correctionNx;
                        correctionNy = ((fittedFunction.deriv(myPoint)(2, 0) * THB.derivSingle(s, myPoint)(0, 0) +
                                         fittedFunction.deriv(myPoint)(3, 0) * THB.derivSingle(s, myPoint)(1, 0))) *
                                       mappedW(n, 0) + correctionNy;
                        correctionDx = (((THB.derivSingle(s, myPoint)(0, 0) * THB.derivSingle(s, myPoint)(0, 0)) +
                                         (THB.derivSingle(s, myPoint)(1, 0) * THB.derivSingle(s, myPoint)(1, 0)))) *
                                       mappedW(n, 0) + correctionDx;
                        correctionDy = (((THB.derivSingle(s, myPoint)(0, 0) * THB.derivSingle(s, myPoint)(0, 0)) +
                                         (THB.derivSingle(s, myPoint)(1, 0) * THB.derivSingle(s, myPoint)(1, 0)))) *
                                       mappedW(n, 0) + correctionDy;

                    }

                }

            }

            //Reset domain iterator
            domIter1.next();
            domIter1.reset();

            //Save corrections
            xCorr.push_back(correctionNx / correctionDx);
            yCorr.push_back(correctionNy / correctionDy);

            //Restore correct control points
            fittedFunction.coefs()(s, 0) = Csx;
            fittedFunction.coefs()(s, 1) = Csy;

            gsInfo << "Opt_CHECK" << std::endl;
        }

        //Apply the corrections
        for (int r = 0; r < selectedBasisFunctions.size(); r++) {
            fittedFunction.coefs()(selectedBasisFunctions[r], 0) = -(xCorr[r]);
            fittedFunction.coefs()(selectedBasisFunctions[r], 1) = -(yCorr[r]);
        }
    }
}

int pickCell(gsVector < int > nonCheckedCells, int & currArrayIndex, int levNow, int & x1U, int & y1U, int & x2U, int & y2U) {
    currArrayIndex = rand() % nonCheckedCells.size();
    int currCellIndex = nonCheckedCells(currArrayIndex);
    x1U = (4 * (int) pow(2, levNow) - 1) - (currCellIndex % (int)(4 * (int) pow(2, levNow)));
    x2U = x1U + 1;
    y1U = (4 * (int) pow(2, levNow) * 4 * (int) pow(2, levNow) - 1 - currCellIndex) / (4 * (int) pow(2, levNow));
    y2U = y1U + 1;
    return currCellIndex;
}

int pickCell(int attempt, int levNow, int & x1U, int & y1U, int & x2U, int & y2U) {
    int currCellIndex = attempt;
    x1U = (4 * (int) pow(2, levNow) - 1) - (currCellIndex % (int)(4 * (int) pow(2, levNow)));
    x2U = x1U + 1;
    y1U = (4 * (int) pow(2, levNow) * 4 * (int) pow(2, levNow) - 1 - currCellIndex) / (4 * (int) pow(2, levNow));
    y2U = y1U + 1;
    return currCellIndex;
}



void restoreTheHierarchy(int & createdBoxNum, int & lastNonzeroRow, int  boxMat[][5], int levNow, int centerInd, int ourBox[], int successfullAttempts) {
    for (int l = 0; l < createdBoxNum; l++) {
        boxMat[lastNonzeroRow - l][0] = levNow; //Preparation for multilevel meshes
        //                            gsInfo
        outfile << "updated coordinates of " << lastNonzeroRow - l << "box: " <<
                boxMat[lastNonzeroRow - l][0];
        for (int m = 1; m < 5; m++) {
            boxMat[lastNonzeroRow - l][m] = 0;
            //                                gsInfo
            outfile << "\t" << boxMat[lastNonzeroRow - l][m];
        }
        //                            gsInfo
        outfile << "\n";
    }
    outfile << "current " << centerInd << "box is now:";
    for (int k = 0; k < 5; k++) {
        boxMat[centerInd][k] = ourBox[k];
        //                            gsInfo
        outfile << "\t" << boxMat[centerInd][k];
    }
    outfile << "\n";
    if (successfullAttempts == 0) {
        outfile << "0th position IS RETURNED TO\n";
        for (int k = 0; k < 5; k++) {
            boxMat[0][k] = ourBox[k];
            //                                gsInfo
            outfile << "\t" << boxMat[0][k];
        }
        outfile << "\n";
    }
    lastNonzeroRow = lastNonzeroRow - createdBoxNum;
    createdBoxNum = 0;
}

int rebuildTheHierarchy(int boxMat[][5], int row, int x1U, int x1Bi, int x2U, int x2Bi, int y1U, int y1Bi, int y2U, int y2Bi, int levNow, int & lastNonzeroRow,
                        int & createdBoxNum, int & centerInd, int ourBox[], int & needToEscape) {

    int RTH = 0;
    if (boxMat[row][1] <= x1U * pow(2, boxMat[row][0] - levNow) &&
        x1U * pow(2, boxMat[row][0] - levNow) <= boxMat[row][3] &&
        boxMat[row][2] <= y1U * pow(2, boxMat[row][0] - levNow) &&
        y1U * pow(2, boxMat[row][0] - levNow) <= boxMat[row][4] &&
        boxMat[row][1] <= x2U * pow(2, boxMat[row][0] - levNow) &&
        x2U * pow(2, boxMat[row][0] - levNow) <= boxMat[row][3] &&
        boxMat[row][2] <= y2U * pow(2, boxMat[row][0] - levNow) &&
        y2U * pow(2, boxMat[row][0] - levNow) <= boxMat[row][4]
            ) {
        //        gsInfo
        outfile << "ALERT: OUR CELL OF LEVEL " << levNow << "\n" << x1U << " " << y1U << " " << x2U <<
                " " << y2U <<
                " INTERSECTS WITH [" << row << "]th box of level " <<
                boxMat[row][0] << ":\n" << boxMat[row][1] << " " << boxMat[row][2] << " " <<
                boxMat[row][3] << " " << boxMat[row][4] << "\n";
        //        wasRebuilt = 0;
        //                for(int checkRedBox = 0; checkRedBox < numRed; checkRedBox++){
        //                    if(boxMat[row][0] == redBox[checkRedBox][0] &&
        //                        redBox[checkRedBox][1] <= x1U * pow(2, redBox[checkRedBox][0] - coarseLevel)
        //                        && x1U * pow(2, redBox[checkRedBox][0] - coarseLevel) <= redBox[checkRedBox][3]
        //                        && redBox[checkRedBox][2] <= y1U * pow(2, redBox[checkRedBox][0] - coarseLevel)
        //                        && y1U * pow(2, redBox[checkRedBox][0] - coarseLevel) <= redBox[checkRedBox][4]
        //                        &&  redBox[checkRedBox][1] <= x2U * pow(2, redBox[checkRedBox][0] - coarseLevel)
        //                        && x2U * pow(2, redBox[checkRedBox][0] - coarseLevel) <= redBox[checkRedBox][3]
        //                        && redBox[checkRedBox][2] <= y2U * pow(2, redBox[checkRedBox][0] - coarseLevel)
        //                        && y2U * pow(2, redBox[checkRedBox][0] - coarseLevel) <= redBox[checkRedBox][4]
        //                        ){
        ////                                gsInfo  << x1U * pow(2, redBox[checkRedBox][0] - 3) << " " <<  y1U * pow(2, redBox[checkRedBox][0] - 3) << " " <<  x2U * pow(2, redBox[checkRedBox][0] - 3) << " " <<  y2U * pow(2, redBox[checkRedBox][0] - 3) <<
        ////                                                                " DOES NOT İNTERSECT WITH " << row << "th box of level " << redBox[checkRedBox][0] << " "<< redBox[checkRedBox][1] << " " << redBox[checkRedBox][2] << " "
        ////                                                                                                       << redBox[checkRedBox][3] << " " << redBox[checkRedBox][4] << "\n";
        //                        beingInRed = 1;
        //                        gsInfo  << "beingInRed:" << beingInRed << "\n";
        //                        break;
        //                    }
        ////                    else{
        ////                        gsInfo  << x1U * pow(2, redBox[checkRedBox][0] - coarseLevel) << " " <<  y1U * pow(2, redBox[checkRedBox][0] - coarseLevel) << " " <<  x2U * pow(2, redBox[checkRedBox][0] - coarseLevel) << " " <<  y2U * pow(2, redBox[checkRedBox][0] - coarseLevel) <<
        ////                                " DOES NOT İNTERSECT WITH RED" << row << "th box of level " << redBox[checkRedBox][0] << " "<< redBox[checkRedBox][1] << " " << redBox[checkRedBox][2] << " "
        ////                                << redBox[checkRedBox][3] << " " << redBox[checkRedBox][4] << "\t RedCheck finished" << "\n";
        ////                        beingInRed = 0;
        ////                    }
        //                }
        if (boxMat[row][0] <= levNow) {
            //            gsInfo
            outfile << "NO NEED TO REBUILD \n";
        }
            //HERE I START TO WRITE CODE IN MORE GENERAL WAY
        else if (boxMat[row][0] > levNow + 1) {
            //            gsInfo
            outfile << "I DON'T TOUCH THESE BOX AT THIS ITERATION!\n";
        }
            //                        else
        else if (boxMat[row][0] == levNow + 1) {
            x1Bi = boxMat[row][1];
            y1Bi = boxMat[row][2];
            x2Bi = boxMat[row][3];
            y2Bi = boxMat[row][4];
            //            gsInfo
            outfile << "REBUILD\n";
            RTH = 1;
            //                            needToEscape = 1;
            //SW W NW
            if (x1U * pow(2, boxMat[row][0] - levNow) - x1Bi > 0) {
                if (y1U * pow(2, boxMat[row][0] - levNow) - y1Bi > 0) {
                    //                    gsInfo
                    outfile << "CREATE SW BOX:\n";
                    boxMat[lastNonzeroRow + 1][0] = levNow + 1;
                    boxMat[lastNonzeroRow + 1][1] = x1Bi; //*pow(2, boxMat[row][0] - 3);
                    boxMat[lastNonzeroRow + 1][2] = y1Bi; //*pow(2, boxMat[row][0] - 3);
                    boxMat[lastNonzeroRow + 1][3] = x1U * pow(2, boxMat[row][0] - (levNow));
                    boxMat[lastNonzeroRow + 1][4] = y1U * pow(2, boxMat[row][0] - (levNow));
                    //                    gsInfo
                    outfile << "created the following box of level " << levNow + 1 << ": " <<
                            boxMat[lastNonzeroRow + 1][1] << " " << boxMat[lastNonzeroRow + 1][2] << " " <<
                            boxMat[lastNonzeroRow + 1][3] << " " << boxMat[lastNonzeroRow + 1][4] << "\n";
                    lastNonzeroRow++;
                    createdBoxNum++;
                }
                //                gsInfo
                outfile << "CREATE W BOX:\n";
                boxMat[lastNonzeroRow + 1][0] = levNow + 1;
                boxMat[lastNonzeroRow + 1][1] = x1Bi; //*pow(2, boxMat[row][0] - 3);
                boxMat[lastNonzeroRow + 1][2] = y1U * pow(2, boxMat[row][0] - levNow);
                boxMat[lastNonzeroRow + 1][3] = x1U * pow(2, boxMat[row][0] - levNow);
                boxMat[lastNonzeroRow + 1][4] = y2U * pow(2, boxMat[row][0] - levNow);
                //                gsInfo
                outfile << "created the following box of level " << levNow + 1 << ": " <<
                        boxMat[lastNonzeroRow + 1][1] << " " << boxMat[lastNonzeroRow + 1][2] << " " <<
                        boxMat[lastNonzeroRow + 1][3] << " " << boxMat[lastNonzeroRow + 1][4] << "\n";
                //                                gsInfo << "row is " << lastNonzeroRow + 1 << "\n";
                lastNonzeroRow++;
                createdBoxNum++;
                if (y2Bi - y2U * pow(2, boxMat[row][0] - levNow) > 0) {
                    //                    gsInfo
                    outfile << y2Bi << "\t" << y2U << "\t" << boxMat[row][0] << "\t" <<
                            y2Bi - y2U * pow(2, boxMat[row][0] - levNow) << "\n";
                    //                    gsInfo
                    outfile << "CREATE NW BOX:\n";
                    boxMat[lastNonzeroRow + 1][0] = levNow + 1;
                    boxMat[lastNonzeroRow + 1][1] = x1Bi; //*pow(2, boxMat[row][0] - 3);
                    boxMat[lastNonzeroRow + 1][2] = y2U * pow(2, boxMat[row][0] - levNow);
                    boxMat[lastNonzeroRow + 1][3] = x1U * pow(2, boxMat[row][0] - levNow);
                    boxMat[lastNonzeroRow + 1][4] = y2Bi;
                    //                    gsInfo
                    outfile << "created the following box of level " << levNow + 1 << ": " <<
                            boxMat[lastNonzeroRow + 1][1] << " " << boxMat[lastNonzeroRow + 1][2] << " " <<
                            boxMat[lastNonzeroRow + 1][3] << " " << boxMat[lastNonzeroRow + 1][4] << "\n";
                    lastNonzeroRow++;
                    createdBoxNum++;
                }
            }
            //S N
            if (y1U * pow(2, boxMat[row][0] - levNow) - y1Bi > 0) {
                //                gsInfo
                outfile << "CREATE S BOX:\n";
                boxMat[lastNonzeroRow + 1][0] = levNow + 1;
                boxMat[lastNonzeroRow + 1][1] = x1U * pow(2, boxMat[row][0] - levNow);
                boxMat[lastNonzeroRow + 1][2] = y1Bi; //*pow(2, boxMat[row][0] - 3);
                boxMat[lastNonzeroRow + 1][3] = x2U * pow(2, boxMat[row][0] - levNow);
                boxMat[lastNonzeroRow + 1][4] = y1U * pow(2, boxMat[row][0] - levNow);
                //                gsInfo
                outfile << "created the following box of level " << levNow + 1 << ": " <<
                        boxMat[lastNonzeroRow + 1][1] << " " << boxMat[lastNonzeroRow + 1][2] << " " <<
                        boxMat[lastNonzeroRow + 1][3] << " " << boxMat[lastNonzeroRow + 1][4] << "\n";
                lastNonzeroRow++;
                createdBoxNum++;
            }
            if (y2Bi - y2U * pow(2, boxMat[row][0] - levNow) > 0) {
                //                gsInfo
                outfile << "CREATE N BOX:\n";
                boxMat[lastNonzeroRow + 1][0] = levNow + 1;
                boxMat[lastNonzeroRow + 1][1] = x1U * pow(2, boxMat[row][0] - levNow);
                boxMat[lastNonzeroRow + 1][2] = y2U * pow(2, boxMat[row][0] - levNow);
                boxMat[lastNonzeroRow + 1][3] = x2U * pow(2, boxMat[row][0] - levNow);
                boxMat[lastNonzeroRow + 1][4] = y2Bi; //y2U * pow(2, boxMat[row][0] - 3);
                //                gsInfo
                outfile << "created the following box of level " << levNow + 1 << ": " <<
                        boxMat[lastNonzeroRow + 1][1] << " " << boxMat[lastNonzeroRow + 1][2] << " " <<
                        boxMat[lastNonzeroRow + 1][3] << " " << boxMat[lastNonzeroRow + 1][4] << "\n";
                lastNonzeroRow++;
                createdBoxNum++;
            }
            //SE E NE
            if (x2Bi - x2U * pow(2, boxMat[row][0] - levNow) > 0) {
                //                gsInfo
                outfile << x2Bi << ", x2U*2 = " << x2U * pow(2, boxMat[row][0] - levNow) << "\n";
                if (y1U * pow(2, boxMat[row][0] - levNow) - y1Bi > 0) {
                    //                    gsInfo
                    outfile << "CREATE SE BOX:\n";
                    boxMat[lastNonzeroRow + 1][0] = levNow + 1;
                    boxMat[lastNonzeroRow + 1][1] = x2U * pow(2, boxMat[row][0] - levNow);
                    boxMat[lastNonzeroRow + 1][2] = y1Bi; //*pow(2, boxMat[row][0] - 3);
                    boxMat[lastNonzeroRow + 1][3] = x2Bi; //*pow(2, boxMat[row][0] - 3);
                    boxMat[lastNonzeroRow + 1][4] = y1U * pow(2, boxMat[row][0] - levNow);
                    //                    gsInfo
                    outfile << "created the following box of level" << levNow + 1 << ": " <<
                            boxMat[lastNonzeroRow + 1][1] << " " << boxMat[lastNonzeroRow + 1][2] << " " <<
                            boxMat[lastNonzeroRow + 1][3] << " " << boxMat[lastNonzeroRow + 1][4] << "\n";
                    lastNonzeroRow++;
                    createdBoxNum++;
                }
                //                gsInfo
                outfile << "CREATE E BOX:\n";
                boxMat[lastNonzeroRow + 1][0] = levNow + 1;
                boxMat[lastNonzeroRow + 1][1] = x2U * pow(2, boxMat[row][0] - levNow);
                boxMat[lastNonzeroRow + 1][2] = y1U * pow(2, boxMat[row][0] - levNow);
                boxMat[lastNonzeroRow + 1][3] = x2Bi; //*pow(2, boxMat[row][0] - 3);
                boxMat[lastNonzeroRow + 1][4] = y2U * pow(2, boxMat[row][0] - levNow);
                //                gsInfo
                outfile << "created the following box of level" << levNow + 1 << ": " <<
                        boxMat[lastNonzeroRow + 1][1] << " " << boxMat[lastNonzeroRow + 1][2] << " " <<
                        boxMat[lastNonzeroRow + 1][3] << " " << boxMat[lastNonzeroRow + 1][4] << "\n";
                lastNonzeroRow++;
                createdBoxNum++;
                if (y2Bi - y2U * pow(2, boxMat[row][0] - levNow) > 0) {
                    //                    gsInfo
                    outfile << "CREATE NE BOX:\n";
                    boxMat[lastNonzeroRow + 1][0] = levNow + 1;
                    boxMat[lastNonzeroRow + 1][1] = x2U * pow(2, boxMat[row][0] - levNow);
                    boxMat[lastNonzeroRow + 1][2] = y2U * pow(2, boxMat[row][0] - levNow);
                    boxMat[lastNonzeroRow + 1][3] = x2Bi; //*pow(2, boxMat[row][0] - 3);
                    boxMat[lastNonzeroRow + 1][4] = y2Bi; //*pow(2, boxMat[row][0] - 3);
                    //                    gsInfo
                    outfile << "created the following box of level" << levNow + 1 << ": " <<
                            boxMat[lastNonzeroRow + 1][1] << " " << boxMat[lastNonzeroRow + 1][2] << " " <<
                            boxMat[lastNonzeroRow + 1][3] << " " << boxMat[lastNonzeroRow + 1][4] << "\n";
                    lastNonzeroRow++;
                    createdBoxNum++;
                }
            }
            //            gsInfo
            outfile << "CREATE C BOX:\n";
            outfile << "ourBOX is ";
            for (int k = 0; k < 5; k++) {
                ourBox[k] = boxMat[row][k];
                outfile << ourBox[k] << " ";
                centerInd = row;
            }
            outfile << "\n The centerInd is " << centerInd << "\n";
            boxMat[row][0] = levNow;
            boxMat[row][1] = x1U;
            boxMat[row][2] = y1U;
            boxMat[row][3] = x2U;
            boxMat[row][4] = y2U;
            //            gsInfo
            outfile << "created the following box of level" << levNow << ": " << x1U << " " << y1U <<
                    " " << x2U << " " << y2U << "\n";
            //                       boxMat[lastNonzeroRow + 1]

        }
    } else {
        outfile << "OUR CELL OF LEVEL " << levNow << "\n" << x1U << " " << y1U << " " << x2U <<
                " " << y2U <<
                " DOES NOT INTERSECT WITH [" << row << "]th box of level " <<
                boxMat[row][0] << ":\n" << boxMat[row][1] << " " << boxMat[row][2] << " " <<
                boxMat[row][3] << " " << boxMat[row][4] << "\n";
    }
    return RTH;
}



void setSupports2(gsTHBSplineBasis < 2, real_t > ::uPtr THB, gsMatrix<> supps){
    for (int i = 0; i < THB->size(); ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                supps(j,2*i + k) = THB->support(i)(j,k);
            }
        }
    }
}

void setSupports(gsTHBSplineBasis<2, real_t> THB, gsMatrix<> supps){
    for (int i = 0; i < THB.size(); i++)
    {
        for (int j = 0; j < THB.support(i).rows(); j++)
        {
            for (int k = 0; k < THB.support(i).cols(); k++)
            {
                supps(j,2*i + k) = THB.support(i)(j,k);
            }
        }
    }
}

real_t getSupports(int i, int j, int k, gsMatrix<> supps){
    return supps(j, 2*i + k);
}



int main(int argc, char * argv[]) {
    gsStopwatch clock;
    int row, acceptedsize, attempt = 0;
//    std::string xmlFile = "/home/turing/theydarov/xmlFiles/resultLocalsuccessstuff";
    std::string badFile = "/home/turing/theydarov/xmlFiles/badgeoLocalsuccessstuff";
    std::string pvdFile = "/home/turing/theydarov/pvdFiles/resultLocalsuccessstuff";
    int los = 0, nlos = 0, proj = 0;
    //DON'T FORGET ABOUT GLOBAL VARIABLES!
    int numCells = 4 * (int) pow(2, coarseLevel) * 4 * (int) pow(2, coarseLevel);
    std::string givenGeo;
    int gradingExtent;
    //    real_t epsilon_g = 0.003, epsilon_f = 0.001;
    //    real_t epsilon_g = 0.004, epsilon_f = 0.004;
    //    real_t epsilon_g = 0.0114, epsilon_f = 0.0114195;
    //    real_t epsilon_g = 1e+6, epsilon_f = 0.003;
    real_t epsilon_g = 1e+6, epsilon_f = 5e-5;
    real_t lcx,lcy,ucx,ucy;
    //    real_t epsilon_g = 1e+6, epsilon_f = 1e+6;
    //    ofstream myfile;
    //    myfile.open ("/home/theydarov/october/example.txt",  fstream::app);
    std::string acCond = to_string(epsilon_g) + "and" + to_string(epsilon_f);
    int bads = 0;

    std::string inputInitial("/home/turing/theydarov/geometries/jigsawpuzzleTP.xml");
//    std::string inputInitial("/home/turing/theydarov/geometries/indiana_thb.xml");

    givenGeo = "jigsawpuzzleTP";
//    givenGeo = "jigsawpuzzleTPEB";
    //    givenGeo = "wigglysurface";
//        givenGeo = "Indiana";
//    givenGeo+=  "WigglyCoonsIntDist";
    //    givenGeo = "Indiana";
    //    givenGeo+=  "WigglySquareIntDist";
    //    givenGeo+=  "WigglyCoonsPatch3lev";
    //    givenGeo+= "lexicog/**/raphic";
    //    givenGeo+="upperhalf";
    //    givenGeo+="degree2";
    givenGeo += "L2";
    givenGeo += "LO";
    givenGeo += "NLO";
    //    givenGeo+="[0,1]cross[0,0dot25]";
    //    givenGeo = "thb_map";
    std::string fileLoc = "/home/turing/theydarov/txtFiles/resultLocalsuccessstuff";
    std::string input("/home/theydarov/june21/300621/" + givenGeo + acCond + ".xml" );
    outfile.open(fileLoc + ".txt");
    std::string inputPtsParams("/home/turing/theydarov/gismo/extensions/motor/filedata/jku/thb_parameters_and_points.xml");
    //    gsInfo << "I REQUEST " << numCells << " cells\n";
    outfile << "I REQUEST " << numCells << " cells\n";

    //Iterations
    int iter = -1;
    int successOfOpt = 0;
    int successfullAttempts = 0, totalAttempts = 0;
    //Error threshold
    real_t threshold = 1e-2;

    //    int boxMat[4*(int)pow(2,coarseLevel)*4*(int)pow(2,coarseLevel)][5];
    gsInfo << "GIVE A NUMBER OF RED BOXES\n";
    outfile << "GIVE A NUMBER OF RED BOXES\n";
    gsInfo << "I HOPE YOU HAVE /home/theydarov/gismo_unsupported/gismo/extensions/motor/filedata/" + givenGeo + acCond + ".xml\n";
    outfile << "I HOPE YOU HAVE /home/theydarov/gismo_unsupported/gismo/extensions/motor/filedata/" + givenGeo + acCond + ".xml\n";
    int beingInRed;
    cin >> numRed;
    if (numRed == 0) {
        redBox[0][0] = fineLevel;
        redBox[0][1] = 0;
        redBox[0][2] = 0;
        redBox[0][3] = 4 * (int) pow(2, fineLevel);
        redBox[0][4] = 4 * (int) pow(2, fineLevel);
        //        gsInfo << "I ASSUME THAT REDBOX HAS THE COORDINATES:\n"
        outfile << "I ASSUME THAT REDBOX HAS THE COORDINATES:\n"

                <<
                redBox[0][0] << " " << redBox[0][1] << " " << redBox[0][2] << " " << redBox[0][3] << " " << redBox[0][4] << "\n";

    }
    real_t boundaryerror = 0;

    //    gsMatrix<> redBox(numRed, 5);
    //    int redBox[numRed][5];
    for (int i = 0; i < numRed; i++) {
        gsInfo << "GIVE RED BOX\n";
        cin >> redBox[i][0] >> redBox[i][1] >> redBox[i][2] >> redBox[i][3] >> redBox[i][4];
    }
    gsInfo << "Give the grading extent:\n";
    cin >> gradingExtent;
    //    clock.restart();
    numRed = 1;

    int RTH = 0, centerInd = 0;
    int lastNonzeroRow = 0; //71;//92;
    int minusnumber, createdBoxNum = 0, currCellIndex, wasRebuilt;
    int x1U, x2U, y1U, y2U; // = rand()%32;
    int x1Bi, x2Bi, y1Bi, y2Bi;
    int needToEscape = 0;
    int createSpline, failed;
    int ourBox[5], currArrayIndex;
    gsMatrix < index_t > lowCorners;
    gsMatrix < index_t > upCorners;
    gsVector < index_t > myLevel;
    currArrayIndex = 0;
    //    std::string  inputInitial;
    //    std::string
    //        theName = inputInitial =   "/home/theydarov/gismo_unsupported/gismo/extensions/motor/filedata/" + givenGeo + acCond + ".xml";
    gsFileData < > data0(inputInitial);

    std::string output("results_optimized");
    //    if (levNow != coarseLevel || wasRebuilt)//if not the first iteration
    //    {
    //        data0.clear();
    //        data0.read(theName);
    //        gsInfo << "Changed the file\n";
    //        outfile << "Changed the file\n";
    //    }
    gsGeometry < > ::uPtr geom;
    gsGeometry < > ::uPtr geomhat; //accepted geometry
    gsGeometry < > ::uPtr geomtmp;
    gsGeometry<>::uPtr geomForCheck;
    gsMatrix < > coefs;
    gsMatrix < > coefshat;
    //                gsGeometry<>::uPtr geom;
    if (data0.has < gsGeometry < > > ()) {
        geom = data0.getFirst < gsGeometry < > > ();
        geomForCheck = data0.getFirst<gsGeometry<> >();
        coefs = geom -> coefs();
        //        geomhat->setCoefs(coefs);//        geomhat = geom;

        //        geomForCheck = data0.getFirst<gsGeometry<> >();
        //        if(!wasRebuilt)
        //        {
        //            gsInfo << inputInitial << " is chosen\n";
        //            outfile << inputInitial << " is chosen\n";
        //        }
        //        if (levNow != coarseLevel || wasRebuilt)//if not the first iteration
        //        {
        //            gsInfo << theName << " is chosen\n";
        //            outfile << theName << " is chosen\n";
        //        }
    }
    coefshat = coefs;

    gsTHBSplineBasis < 2, real_t > ::uPtr THBFromGeo;

//    gsTHBSplineBasis < 2, real_t > ::uPtr THBAccepted;

    THBFromGeo = data0.getAnyFirst < gsTHBSplineBasis < 2 >> ();
    acceptedsize = THBFromGeo->size();
    gsMatrix<> anmat2;
    gsMatrix<> anmat ;
    THBFromGeo->anchors_into(anmat2);
//    gsInfo << anmat2.rows() << " " << anmat2.cols() <<  "anmat:" << anmat2 <<"\n";
//    gsInfo << THBFromGeo->support().cols() << "a\n";
    gsMatrix<> supps(THBFromGeo->support().rows(),THBFromGeo->support().cols()+ 2*THBFromGeo->size());
//    for (int k = 0; k < THBFromGeo->size(); ++k)
//                {
////
//                    gsInfo << k <<"-th support: " << THBFromGeo->support(k) << "\n";
//                }

    //setSupports2(THBFromGeo, supps);
    for (int i = 0; i < THBFromGeo->size(); ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                supps(j,2*i + k) = THBFromGeo->support(i)(j,k);
            }
        }
    }
    (THBFromGeo -> tree().getBoxes(lowCorners, upCorners, myLevel));
    gsInfo << "FOUND " << lowCorners.rows() << " BOXES\n";
    outfile << "FOUND " << lowCorners.rows() << " BOXES\n";
    int boxMat[2*numCells][5];//LAZY SOLUTION FOR NOW
    for (int i = 0; i < lowCorners.rows(); i++) {
        boxMat[i][0] = myLevel(i);
        boxMat[i][1] = (real_t) lowCorners(i, 0) / pow(2, THBFromGeo -> maxLevel() - myLevel(i));
        boxMat[i][2] = (real_t) lowCorners(i, 1) / pow(2, THBFromGeo -> maxLevel() - myLevel(i));
        boxMat[i][3] = (real_t) upCorners(i, 0) / pow(2, THBFromGeo -> maxLevel() - myLevel(i));
        boxMat[i][4] = (real_t) upCorners(i, 1) / pow(2, THBFromGeo -> maxLevel() - myLevel(i));
        gsInfo << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
               boxMat[i][4] << "\n";
        outfile << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                boxMat[i][4] << "\n";
    }
    lastNonzeroRow = lowCorners.rows() - 1;
    real_t au = 0; // starting knot
    real_t bu = 1; // ending knot
    real_t av = 0; //-1;//0;                   // starting knot
    real_t bv = 1; //0.5;//1;                   // ending knot

    unsigned interior = 3; // number of interior// knots for geometry by Ondine
    real_t x1start = 0.6, x2start = 0.6, y1start = 0.6, y2start = 0.6;
    int degree = 2;
    unsigned multEnd = degree + 1; // multiplicity at the two end knots
    gsKnotVector <> ku(au, bu, interior, multEnd);
    gsKnotVector <> kv(av, bv, interior, multEnd);
    gsTensorBSplineBasis < 2, real_t > tens(ku, kv);
    gsTHBSplineBasis < 2, real_t > THBAccepted(tens);
    gsTHBSplineBasis < 2, real_t > THBTemporary(tens);
    //    gsTHBSplineBasis<2, real_t> THBAccepted(tens);
    wasRebuilt = 0;
    real_t maxError = 0;
    THBAccepted = *THBFromGeo;
    //    *geomhat = *geom;
//    gsInfo << "TEST";
    gsGeometry<>::uPtr fittedFunction;
    int iteration = 0;
    std::string xmlFile = "/home/turing/theydarov/output/" + givenGeo + acCond;
    //    = geomhat->basis();
    for (int levNow = coarseLevel; levNow >= 0; levNow--)
//                for(int levNow = coarseLevel; levNow >= coarseLevel; levNow--)
    {
        iteration = 0;
        theLev = levNow;
//        for (int i = 0; i < ; ++i) {
//
//        }

        if (levNow < coarseLevel && wasRebuilt)
        {
            std::string input(xmlFile + ".xml" );
            gsFileData < > data1(input);
            THBAccepted = *data1.getAnyFirst<gsTHBSplineBasis<2 >>();
//            gsMatrix<> anmat2 = THBFromGeo->anchors();

            (THBAccepted.tree().getBoxes(lowCorners, upCorners, myLevel));
            gsInfo << "FOUND " << lowCorners.rows() << " BOXES\n";
            outfile << "FOUND " << lowCorners.rows() << " BOXES\n";
            for (int i = 0; i < lowCorners.rows(); i++)
            {
                boxMat[i][0] = myLevel(i);
                boxMat[i][1] = (real_t) lowCorners(i, 0) / pow(2, std::min((int)THBFromGeo->maxLevel(), (int)THBAccepted.maxLevel()) - myLevel(i));
                boxMat[i][2] = (real_t) lowCorners(i, 1) / pow(2, std::min((int)THBFromGeo->maxLevel(), (int)THBAccepted.maxLevel()) - myLevel(i));
                boxMat[i][3] = (real_t) upCorners(i, 0) / pow(2, std::min((int)THBFromGeo->maxLevel(), (int)THBAccepted.maxLevel()) - myLevel(i));
                boxMat[i][4] = (real_t) upCorners(i, 1) / pow(2, std::min((int)THBFromGeo->maxLevel(), (int)THBAccepted.maxLevel()) - myLevel(i));
                gsInfo << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                       boxMat[i][4] << "\n";
                outfile << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                        boxMat[i][4] << "\n";
            }
            lastNonzeroRow = lowCorners.rows() - 1;
        }

        //        gsInfo << "working with level " << levNow << " as a coarseLevel\n";
        outfile << "working with level " << levNow << " as a coarseLevel\n";
        //                gsFunctionExpr<> f1("sin(20*x)/10 - cos((3*x)/2) - cos(3/2) + sin(20)/10 - (cos(3/2)*sin(20))/10 - (cos((3*x)/2)*sin(20*x))/10 + 2", 2);//Coons patch based on Bezier 3rd degree curve
        //                gsFunctionExpr<>f1("sin(20*x)/50 - cos(3*x) - cos(3) + sin(20)/50 - (cos(3)*sin(20))/50 - (cos(3*x)*sin(20*x))/50 + 2", 2);//EXAMPLE ON WHICH WE RECEIVE GLOBAL REFINED 1-LEVEL MESH
        //                gsFunctionExpr<> f1("(31*sin(10*x))/1250 - cos(3*x) - cos(3) + (31*sin(10))/1250 - (31*cos(3)*sin(10))/1250 - (31*cos(3*x)*sin(10*x))/1250 + 2",2);
        //                gsFunctionExpr<> f1("(73*sin(30*x))/5000 - (73*sin(27*x))/10000 - cos(3*x) - (73*sin(33*x))/10000 - cos(3) - (73*sin(27))/10000 + (73*sin(30))/5000 - (73*sin(33))/10000 + 2",2);
        //                gsFunctionExpr<> f1("2*x^3 - 2*x^2 + 4*x*y^2 - 4*x*y + 3*x - y^2 + y - 1",2);
        gsFunctionExpr<> f1("x", 2);

        //                gsFunctionExpr<> f2("y + sin((3*x)/2) + sin(20*x)/10 + sin(3/2) + sin(20)/10 + (y*sin(20*x))/10 + (sin(3/2)*sin(20))/10 + (x*sin(20))/10 + (sin((3*x)/2)*sin(20*x))/10 + 1", 2);//Coons patch based on Bezier 3rd degree curve
        //                gsFunctionExpr<>f2("y + sin(3*x) + sin(20*x)/50 + sin(3) + sin(20)/50 + (y*sin(20*x))/50 + (sin(3)*sin(20))/50 + (x*sin(20))/50 + (sin(3*x)*sin(20*x))/50 + 1", 2);//EXAMPLE ON WHICH WE RECEIVE GLOBAL REFINED 1-LEVEL MESH
        //                gsFunctionExpr<> f2("y + sin(3*x) + (31*sin(10*x))/1250 + sin(3) + (31*sin(10))/1250 + (31*y*sin(10*x))/1250 + (31*sin(3)*sin(10))/1250 + (31*x*sin(10))/1250 + (31*sin(3*x)*sin(10*x))/1250 + 1",2);
        //gsFunctionExpr<> f2("y", 2);
        //                gsFunctionExpr<> f2("y + (73*cos(27*x))/10000 - (73*cos(33*x))/10000 + sin(3*x) + (73*sin(30*x))/5000 + (73*cos(27))/10000 - (73*cos(33))/10000 + sin(3) + (73*sin(30))/5000 + (73*y*sin(30*x))/5000 + (73*x*sin(30))/5000 + 1", 2);
        //                gsFunctionExpr<> f2("- (11*x^2*y)/10 + (81*x^2)/10 - 2*x*y^3 + 4*x*y^2 - (9*x*y)/10 - (81*x)/10 + 4*y^3 - 2*y^2 + y + 1",2);
        //                gsFunctionExpr<> f2("if(y<0, (y+1)/(1+0.5), ((y+1)/(1+0.5))+((y+1)/(1+0.5))^2 *sin(22*pi*x)/20)", 2);
        gsFunctionExpr<> f2("if(y<1/(1+0.5), y, y+y^2 *sin(22*pi*x)/20)", 2);
        int extension = 2;

        gsVector<int> nonCheckedCells(4 * pow(2, levNow) * 4 * pow(2, levNow));
        for (int i = 0; i < nonCheckedCells.size(); i++)
        {
            nonCheckedCells(i) = i;
        }
        int readjustmentNumber = 0;
        for (int i = 0; i < 4 * pow(2, levNow) * 4 * pow(2, levNow); ++i) {
            currCellIndex = pickCell(i, levNow, x1U, y1U, x2U, y2U);
            if(x1U <= x1start && x1start <= x2U && y1U <= y1start && y1start <= y2U)
            {
                readjustmentNumber = i;
                currCellIndex = 0;
                break;
            }
        }

        int success = 1;

//        for (int attempt = 0; attempt < 4 * pow(2, levNow) * 4 * pow(2, levNow); attempt++)
//            for (int attempt = 0; attempt < 2; attempt++)
        attempt = 0;
        while(nonCheckedCells.size() != 0 && success)
        {
            iteration++;
            success = 0;
            gsVector<int> vectorS = nonCheckedCells;

            while(vectorS.size() != 0) {
                lcx = 1.0;
                lcy = 1.0;
                ucx = 0.0;
                ucy = 0.0;
                failed = 0;
                //        currCellIndex = rand() % nonCheckedCells.size();
                outfile << "\n";
                outfile << "\n";
                outfile << "\n";
                outfile << "The boxes\n";
                for (int i = 0; i <= lastNonzeroRow; i++)
                {
                    for (int j = 0; j < 5; j++)
                    {

                        outfile << boxMat[i][j] << "\t";
                    }
                    outfile << "\n";
                }
                createdBoxNum = 0;
//            currCellIndex = pickCell(nonCheckedCells, currArrayIndex, levNow, x1U, y1U, x2U, y2U);
                currCellIndex = pickCell(attempt, levNow, x1U, y1U, x2U, y2U);
                int jopa = 4*(int)pow(2, levNow)*4*(int)pow(2, levNow);

//            currCellIndex = pickCell((attempt + readjustmentNumber)%(4 * (int)pow(2, levNow) * 4 * (int)pow(2, levNow)), levNow, x1U, y1U, x2U, y2U);
                //WE
                gsVector<int> cellstack;
                outfile << "Number of nonchecked cells: " << vectorS.size() << "\n";
                //                gsInfo << "Initial basis:" << *THBFromGeo << std::endl;
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //Creates manually your mesh(put this stuff in a function)

                // ...a 2D-tensor-B-spline basis with this knot vector...

                // ...and a 2D-THB-spline basis out of the tensor-B-spline basis.

                //! [constBasis]

                gsFieldCreator<real_t> JakaField;
                //    for(int attempt = 0; attempt < 256; attempt++){

                gsInfo << "attempt " << attempt << ", CURRENT INDEX IS " << vectorS(currArrayIndex) <<
                       ",\t the coordinates of the box are\n" << x1U << " " << y1U << " " << x2U << " " << y2U << "\n";
                outfile << "attempt " << attempt << ", CURRENT INDEX IS " << vectorS(currArrayIndex) <<
                        ",\t the coordinates of the box are\n" << x1U << " " << y1U << " " << x2U << " " << y2U << "\n";
                //        gsInfo << "attempt " << attempt << "\n" << x1U << " " << y1U << " " << x2U << " " << y2U << "\n";
                //                gsInfo << "need is" << needToEscape << "\n";
                //                THB.getLevelUniqueSpanAtPoints
                //                needToEscape = 1;

                vectorS.removeElement(currArrayIndex);
                attempt = (attempt + 1)%(jopa);
                createSpline = 0;
//            if(attempt!=0)
//            {

//            }
//            gsInfo << "supps before: \n";
//            for (int n = 0; n < acceptedsize; ++n)
//            {
//                gsInfo << supps(0,n) << "; ";
//            }
//            gsInfo << "\n";
//            for (int n = 0; n < acceptedsize; ++n)
//            {
//                gsInfo << supps(1,n) << "; ";
//            }
                for (int currentrow = 0; currentrow < lastNonzeroRow + 1; currentrow++)
                {
                gsInfo << "TRYING TO REBUILD\n";
                    //                minusnumber = 1;
                    row = currentrow;

                    RTH =
                            rebuildTheHierarchy(boxMat, row, x1U, x1Bi, x2U, x2Bi, y1U, y1Bi, y2U, y2Bi, levNow, lastNonzeroRow,
                                                createdBoxNum, centerInd, ourBox, needToEscape);

                    if (RTH) createSpline = 1;
                }

                for (int i = 0; i <= lastNonzeroRow; ++i) {
                    gsInfo << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                                boxMat[i][4] << "\n";
                }
                //            for(int i = 0; i < lastNonzeroRow; i++)
                //            {
                //                gsInfo << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " "
                //                       << boxMat[i][4] << "\n";
                //            }
                if (createSpline == 1)
                {
                    outfile << "creating the spline\n";
                    gsInfo  << "creating the spline\n";
                    gsTHBSplineBasis<2, real_t> THB(tens);
                    gsInfo << THB << "\n";

                    for (int therow = 0; therow < lastNonzeroRow + 1; therow++)
                    {
                        //                gsInfo
                        outfile << therow << "; ";
                        std::vector<index_t> box;
                        for (int column = 0; column < 5; column++)
                        {
                            box.push_back(boxMat[therow][column]);
                            //                    gsInfo
                            outfile << boxMat[therow][column] << "; ";
                        }
                        //                                THBFromGeo->refineElements(box);
                        THB.refineElements(box);
                        box.clear();
                        //                gsInfo
                        outfile << "\n";
                    }

                    int isAdequate = 0;
                    for (int i = 0; i < THB.size(); ++i) {
                        int levOf = THB.levelOf(i);
                        if(levOf >= levNow + 1){
                            isAdequate = 1;
                        }
                    }
                    if(!isAdequate){
                        {
                            outfile << "The basis is inadequate, rewriting\n";
                            std::string input(xmlFile+".xml" );
                            gsFileData < > data2(input);
                            THBTemporary = *data2.getAnyFirst<gsTHBSplineBasis<2 >>();
//            gsMatrix<> anmat2 = THBFromGeo->anchors();

                            (THBTemporary.tree().getBoxes(lowCorners, upCorners, myLevel));
                            gsInfo << "FOUND " << lowCorners.rows() << " BOXES\n";
                            outfile << "FOUND " << lowCorners.rows() << " BOXES\n";
                            for (int i = 0; i < lowCorners.rows(); i++)
                            {
                                boxMat[i][0] = myLevel(i);
                                boxMat[i][1] = (real_t) lowCorners(i, 0) / pow(2, std::min((int)THBFromGeo->maxLevel(), (int)THBTemporary.maxLevel()) - myLevel(i));
                                boxMat[i][2] = (real_t) lowCorners(i, 1) / pow(2, std::min((int)THBFromGeo->maxLevel(), (int)THBTemporary.maxLevel()) - myLevel(i));
                                boxMat[i][3] = (real_t) upCorners(i, 0) / pow(2, std::min((int)THBFromGeo->maxLevel(), (int)THBTemporary.maxLevel()) - myLevel(i));
                                boxMat[i][4] = (real_t) upCorners(i, 1) / pow(2, std::min((int)THBFromGeo->maxLevel(), (int)THBTemporary.maxLevel()) - myLevel(i));
                                gsInfo << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                                       boxMat[i][4] << "\n";
                                outfile << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                                        boxMat[i][4] << "\n";
                            }
                            for (int i = lowCorners.rows(); i < 2*numCells; ++i) {
                                boxMat[i][0] = boxMat[i][1] = boxMat[i][2] = boxMat[i][3] = boxMat[i][5] = 0;
                            }
                            lastNonzeroRow = lowCorners.rows() - 1;
                        }
                        attempt--;
                        continue;
                    }

                    //                            gsInfo << "after refinement" << *THBFromGeo << "\n";
                    //                            std::string exitBasis1 = "THBFromGeo:";
                    //                            gsWriteParaview(*THBFromGeo, exitBasis1);

                    //    gsFileManager::open("THB1p.pvd");
                    //For the mesh
                    //                            std::string exitMesh1 = "outputMesh1";
                    //                            gsWrite(*THBFromGeo, exitMesh1);

                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    //                if(wasRebuilt == 0)
                    //                {
                    //                    continue;
                    //                }
                    //                gsInfo
                    outfile << THB << "\n";
                    //                                cout << "Before refinement:" << endl;
                    gsInfo << THB << "\n";

                    THB.anchors_into(anmat);

                    //                THB2 = THB;

                    gsMatrix<int> basisInd(1, THB.size()); //1 if function is not presented
                    for (int l = 0; l < THB.size(); ++l)
                    {
                        basisInd(0, l) = 0;
                    }
//                gsMatrix<> anmat = THB.anchors();
//                gsInfo << "supports:\n";
//                for (int k = 0; k < THB.size(); ++k)
//                {
//
//                    gsInfo << THBAccepted.support(k) << "\t";
//                }
//                if (attempt != 0)
//                {
//                    gsMatrix<> anmat2 = THBAccepted.anchors();
//                }
                    gsMatrix<> THBcoefs(THB.size(), 2);
                    for (int k = 0; k < THB.size(); ++k)
                    {
                        THBcoefs(k,0) = THBcoefs(k,1) = 0.0;
                    }
                    for (int j = 0; j < THB.size(); j++)
                    {
                        for (int i = 0; i < anmat2.cols(); i++)
                        {
                            //                        gsInfo  << *THBAccepted  << "\n";
                            if (
                                    THB.support(j)(0, 0) == getSupports(i, 0, 0, supps)
                                    && THB.support(j)(1, 0) == getSupports(i, 1, 0, supps)
                                    && THB.support(j)(0, 1) == getSupports(i, 0, 1, supps)
                                    && THB.support(j)(1, 1) == getSupports(i, 1, 1, supps)
//                                &&
//                                anmat2(0,i) == anmat(0,j) && anmat2(1,i) == anmat(1,j)
                                    )
                            {
                                //                            gsInfo << suppcounter <<", i = " << i << ", j = " << j << "\n";
                                //                            suppcounter++;
                                //                if(coarseanc(j) == fineanc(i))  gsInfo << "And their anchors are identical\n";
                                THBcoefs(j, 0) = coefshat(i, 0);
                                THBcoefs(j, 1) = coefshat(i, 1);

                                basisInd(j) = 0;
                                //                            fineIndIsUsed(i) = 1;
                                break;
                            }
                            //                        else if(coarseanc(j) == fineanc(i) && !fineIndIsUsed(i)){
                            //                            gsInfo << suppcounter <<", i = " << i << ", j = " << j << "BUT ONLY ANCHORS\n";
                            //                        }
                        }
                    }
//                gsInfo << "THBFROmGeo supp:\n";
//                for (int n = 0; n < acceptedsize; ++n)
//                {
//                    gsInfo <<  getSupports(n, 0, 0, supps)  << ";" << getSupports(n, 1, 0, supps)
//                           << ";" << getSupports(n, 0, 1, supps)  << ";" << getSupports(n, 1, 1, supps)
//                           << ";" << anmat2(0,n) << ";" << anmat2(1,n) << "\n";
//                }
//                gsInfo << "THB:\n";
//                for (int n = 0; n < THB.size(); ++n)
//                {
//                    gsInfo << THB.support(n)(0,0) << ";" << THB.support(n)(1,0)
//                           << ";" << THB.support(n)(0,1)  << ";" << THB.support(n)(1,1)
//                           << ";" << anmat(0,n)<< ";" << anmat(1,n)<< "\n";
//                }
                    for (int m = 0; m < THB.size(); m++)
                    {
                        if(THBcoefs(m,0) == 0.0 || THBcoefs(m,1) == 0.0)
                        {
//                        gsInfo << m << "IS PROBLEMATIC\n";
                            basisInd(0,m) = 1;
                        }
                    }

                    for (int i = -gradingExtent; i < gradingExtent; ++i) {
                        for (int j = -gradingExtent; j < gradingExtent; ++j) {
                            updateCriterion(levNow, x1U + i, y1U + j, x1U + i + 1, y1U + j + 1, 1, basisInd, THB);
                        }
                    }




//                for (int j = 0; j < THB.size(); j++)
//                {
//                    if(!basisInd(j)){
//
//                        //                        gsInfo  << *THBAccepted  << "\n";
//
//                        //                            gsInfo << suppcounter <<", i = " << i << ", j = " << j << "\n";
//                        //                            suppcounter++;
//                        //                if(coarseanc(j) == fineanc(i))  gsInfo << "And their anchors are identical\n";
//                        THBcoefs(j, 0) = coefshat(j, 0);
//                        THBcoefs(j, 1) = coefshat(j, 1);
//
////                            basisInd(j) = 0;
//                        //                            fineIndIsUsed(i) = 1;
////                            break;
//
//                        //                        else if(coarseanc(j) == fineanc(i) && !fineIndIsUsed(i)){
//                        //                            gsInfo << suppcounter <<", i = " << i << ", j = " << j << "BUT ONLY ANCHORS\n";
//                        //                        }
//                    }
//                }

                    gsInfo << basisInd << "\n";
                    outfile << basisInd << "\n";
                    for (int i1 = 0; i1 < basisInd.cols(); ++i1)
                    {
                        if(basisInd(i1)){
                            if(THB.support(i1)(0,0) < lcx)  lcx = THB.support(i1)(0,0);
                            if(THB.support(i1)(1,0) < lcy)  lcy = THB.support(i1)(1,0);
                            if(THB.support(i1)(0,1) > ucx)  ucx = THB.support(i1)(0,1);
                            if(THB.support(i1)(1,1) > ucy)  ucy = THB.support(i1)(1,1);
                        }
                    }
                    //                                if(THBFromGeo->sizeByLevel(std::cout,fineLevel)==0){
//                if (THB.sizeByLevel(std::cout, levNow + 1) == 0 && THB.maxLevel() == levNow + 1)
//                {
//                    gsInfo << "INADEQUATE BASIS.\n";
//                    outfile << "INADEQUATE BASIS.\n";
//                    restoreTheHierarchy(createdBoxNum,
//                                        lastNonzeroRow,
//                                        boxMat,
//                                        levNow,
//                                        centerInd,
//                                        ourBox,
//                                        successfullAttempts);
//                    continue;
//                }
//                    gsInfo << "reached here\n";
                    index_t totalNumberBasisFunct = THB.size();
                    gsMatrix<real_t> matrixOfTPLC(2, totalNumberBasisFunct);
                    gsMatrix<real_t> matrixOfTPUC(2, totalNumberBasisFunct);
                    gsMatrix<int> maxLevels(1, totalNumberBasisFunct);
                    gsMatrix<int> minLevels(1, totalNumberBasisFunct);
                    boundingBox(THB, maxLevels, minLevels, matrixOfTPLC, matrixOfTPUC, basisInd);
                    fittedFunction = myLocalApproximate(THB,
                                                        *geom,
                                                        maxError,
                                                        matrixOfTPLC,
                                                        matrixOfTPUC,
                                                        maxLevels,
                                                        minLevels,
                                                        au,
                                                        bu,
                                                        av,
                                                        bv,
                                                        boundaryerror,
                                                        basisInd,
                                                        THBcoefs);


                    int minusnumber = checkJacobianDeterminant(*fittedFunction, 10000, true, "points", 0);//
                    checkAccuracySufficiency(*fittedFunction, *geomForCheck, maxError, boundaryerror, au, bu, av, bv);
                    gsInfo << "the errors: " << maxError << " " << boundaryerror << "\n";
                    outfile << "the errors: " << maxError << " " << boundaryerror << "\n";
                    if (minusnumber == 0 && !(epsilon_g < maxError || epsilon_f < boundaryerror))
                    {
                        outfile << "Success! iteration = " << iteration << ", coarselevel = " << coarseLevel << "\n";
                        success = 1;

                        successfullAttempts++;
                        totalAttempts++;
                        gsInfo << "Success! iteration = " << iteration << ", coarselevel = " << coarseLevel << "\n";
                        gsInfo << "SIZE: " << THB.size() << "\n";
                        outfile  << "SIZE: " << THB.size() << "\n";
                        saveData(*fittedFunction, xmlFile, 0);
                        acceptedsize = THB.size();
                        for (int l = 0; l <= lastNonzeroRow; l++) {

                            //                            gsInfo
                            for (int m = 1; m < 5; m++) {
                                outfile << boxMat[l][m];
                                //                                gsInfo
                                outfile << "\t";
                                gsInfo << boxMat[l][m];
                                //                                gsInfo
                                gsInfo << "\t";
                            }
                            //                            gsInfo
                            outfile << "\n";
                            gsInfo << "\n";
                        }
                        proj++;
                        gsWrite(*fittedFunction,
                                xmlFile);//MAYBE, THERE IS A WAY OF IMPROVING
                        setSupports(THB, supps);
//                    gsInfo << "supports after:\n";
//                    for (int n = 0; n < THB.size(); ++n)
//                    {
//                        gsInfo << supps(0,n) << "; ";
//                    }
//                    gsInfo << "\n";
                        for (int n = 0; n < THB.size(); ++n)
                        {
                            gsInfo << supps(1,n) << "; ";
                        }
                        gsInfo << "\n";
                        anmat2 = THB.anchors();
                        wasRebuilt = 1;
                        nonCheckedCells.removeElement(currArrayIndex);
                        continue;
//                    gsWrite(*fittedFunction, "/home/theydarov/june21/300621/badgeoLocalGrading.xml");
//                    return 0;
                    }
                    else
                    {
                        if(minusnumber != 0){    gsInfo << "NONREGULAR:" << minusnumber << "\n";
                            outfile << "NONREGULAR:" << minusnumber << "\n";}

                        else{
                            gsInfo << "WITHDRAW, CAUSE THE PARAMETRIZATION IS REGULAR\n";
                            outfile << "WITHDRAW, CAUSE THE PARAMETRIZATION IS REGULAR\n";
                            totalAttempts++;
                            gsInfo << maxError << "\t" << boundaryerror << "\n";
                            outfile << maxError << "\t" << boundaryerror << "\n";
//                        if(minusnumber!=0){
//                            gsWrite(*fittedFunction, "/home/theydarov/june21/300621/badgeoLocalGrading.xml");
//                        }


//                        for (int l = 0; l < lastNonzeroRow; l++) {
//
//                            //                            gsInfo
//                            for (int m = 1; m < 5; m++) {
//                                outfile << boxMat[l][m];
//                                //                                gsInfo
//                                outfile << "\t";
//                                gsInfo << boxMat[l][m];
//                                //                                gsInfo
//                                gsInfo << "\t";
//                            }
//                            //                            gsInfo
//                            outfile << "\n";
//                            gsInfo << "\n";
//                        }

                            restoreTheHierarchy(createdBoxNum,
                                                lastNonzeroRow,
                                                boxMat,
                                                levNow,
                                                centerInd,
                                                ourBox,
                                                successfullAttempts);


//                    return 0;
                            for (int l = 0; l <= lastNonzeroRow; l++) {

                                //                            gsInfo
                                for (int m = 0; m < 5; m++) {
                                    outfile << boxMat[l][m];
                                    //                                gsInfo
                                    outfile << "\t";
                                    gsInfo << boxMat[l][m];
                                    //                                gsInfo
                                    gsInfo << "\t";
                                }
                                //                            gsInfo
                                outfile << "\n";
                                gsInfo << "\n";
                            }


                            continue;
                        }
                        gsInfo << "ALERT! TRYİNG TO USE LİNEAR OPTİMİZATİON!!!\n";
                        outfile << "ALERT! TRYİNG TO USE LİNEAR OPTİMİZATİON!!!\n";
                        gsWrite(*fittedFunction, badFile);
                        int iterations = 1;
                        real_t fitting = 0;//0.99999;
                        real_t orthogonality = 0;//1e-5;//1e-12;//1e+18;//1e-12;//1e+12;//1e-10;//1e-10;//1e-13;//1e-8;
                        real_t tim = 0;
                        real_t skewness = 0;//1e-1;//1e+18;//0;//1e-10;//1e+10;//1e+13;//8;//18;//1e-18//5e+17;//1e-7;
                        real_t eccentricity = 0;//1;//1e+18;//1e-18;//1e-18;//8;//0;//1e-15;//1e-8;
                        real_t intersection = 0;//1;//1e-10;//1e-9;
                        real_t uniformity = 1;//1e+12;//e-12;//1e-2;//1e-3;//1e-17;//17;//1e-17;//5e+13;//1e-8;
                        real_t area = 0;//1e-8;//11;//1e-18;//1e+15;//1e+16;//1e-8;//0.00001;
                        real_t length = 0;//1e+6;//1e+3;//1e-18;//1e+18;//8;//1e+13;//1e+17;//1e+18;//0.000002;
                        real_t epsilon = 1e-7;
                        int jacPts = 101 * 101;
                        bool dumped = false;
                        gsFileData<> fd_in(inputPtsParams);
                        gsVector<> lowc1(2);
                        lowc1(0) = lcx;//0 + 1e-6;
                        lowc1(1) = lcy;//0 + 1e-6;
                        const gsVector<real_t> lowc2 = lowc1;
                        gsVector<> uppc1(2);
                        uppc1(0) = ucx;//1 - 1e-6;
                        uppc1(1) = ucy;//1 - 1e-6;
                        const gsVector<real_t> uppc2 = uppc1;
                        gsMatrix<> uv1 = uniformPointGrid(lowc2, uppc2, 400);
                        gsMatrix<> xy1 = fittedFunction->eval(uv1);
                        gsQualityMeasure2<real_t> optimization(*fittedFunction, uv1, xy1, basisInd, false);
                        for (int it = 0; it != iterations; it++)
                        {
//                                               gsInfo
                            outfile << "iteration: " << it << " / " << iterations - 1 << "\n";
                            optimization.optimize(fitting, orthogonality, tim, skewness,
                                                  eccentricity, uniformity,
                                                  length, area,
                                                  intersection, epsilon, dumped);
                            //gsInfo
                            outfile << "it = " << it << " Value of functional: "
                                    << optimization.functional(fitting, orthogonality, tim, skewness,
                                                               eccentricity, uniformity,
                                                               length, area,
                                                               intersection, epsilon)
                                    << "\n";
//                        saveData(*fittedFunction, output, it + 1);
                            minusnumber = checkJacobianDeterminant(*fittedFunction, jacPts, true, output, it + 1);
                            // Giving different names for files with (possibly) final results
                            saveData(*fittedFunction, output + "_final", 0);
                            if (minusnumber == 0) break;
                        }
                        gsInfo << "Number of negative points " << minusnumber << "\n";
                        outfile << minusnumber << "\n";
                        checkAccuracySufficiency(*fittedFunction, *geomForCheck, maxError, boundaryerror, au, bu, av, bv);
//                        if(epsilon_g < maxError || epsilon_f < boundaryerror)
//                        {
//                            gsInfo << "NO NEED\n";
//                            gsInfo << "WITHDRAW\n";
//                            totalAttempts++;
//                            gsInfo << "the errors after LO: " << maxError << " " << boundaryerror << "\n";
//                            outfile << "the errors after LO: " << maxError << " " << boundaryerror << "\n";
//                            outfile << "WITHDRAW\n";
//                            if(minusnumber!=0){
//                                gsWrite(*fittedFunction, "/home/theydarov/june21/300621/badgeoLocalGrading.xml");
//                            }
//                            restoreTheHierarchy(createdBoxNum,
//                                                lastNonzeroRow,
//                                                boxMat,
//                                                levNow,
//                                                centerInd,
//                                                ourBox,
//                                                successfullAttempts);
//                            continue;
//
//                        }
                        if (minusnumber == 0 && !(epsilon_g < maxError || epsilon_f < boundaryerror))
                        {

//                        for (int therow = 0; therow < lastNonzeroRow + 1; therow++)
//                        {
//                            //                gsInfo
//                            outfile << therow << "; CREATING THE TMP SPLINE ";
//                            std::vector<unsigned> box;
//                            for (int column = 0; column < 5; column++)
//                            {
//                                box.push_back(boxMat[therow][column]);
//                                //                    gsInfo
//                                outfile << boxMat[therow][column] << "; ";
//                            }
//                            //                                THBFromGeo->refineElements(box);
//                            THBAccepted.refineElements(box);
//                            box.clear();
//                            //                gsInfo
//                            outfile << "\n";
//                        }
                            successfullAttempts++;
                            success = 1;
                            totalAttempts++;
                            gsWrite(*fittedFunction,
                                    xmlFile);//MAYBE, THERE IS A WAY OF IMPROVING
                            gsInfo << "LO WORKED!\n";
                            outfile << "LO WORKED!\n";
                            los++;
                            gsInfo << "SIZE: " << THB.size() << "\n";
                            outfile  << "SIZE: " << THB.size() << "\n";
                            acceptedsize = THB.size();
                            setSupports(THB, supps);
                            anmat2 = THB.anchors();
                            wasRebuilt = 1;
                            nonCheckedCells.removeElement(currArrayIndex);
                            continue;
//                    gsWrite(*fittedFunction, "/home/theydarov/june21/300621/badgeoLocalGrading.xml");
//                    return 0;
                        }
                        else
                        {
                            std::string input = badFile + ".xml";
                            gsFileData<> data(input);
                            gsGeometry<>::uPtr geom;
                            if (data.has< gsGeometry<> >())
                            {
                                geom = data.getFirst< gsGeometry<> >();
//        gsInfo << "GOT!\n";
                            }

                            if (geom == NULL)
                            {
                                gsInfo << "Didn't get the input geometry. Aborting..." << "\n";
                            }
                            gsInfo << "ALERT AGAIN! TRYİNG TO USE NON-LİNEAR OPTİMİZATİON!!!\n";
                            outfile << "ALERT AGAIN! TRYİNG TO USE NON-LİNEAR OPTİMİZATİON!!!\n";
                            iterations = 10;
                            fitting = 0;//0.99999;
                            orthogonality = 0;//1e-5;//1e-12;//1e+18;//1e-12;//1e+12;//1e-10;//1e-10;//1e-13;//1e-8;
                            tim = 0;
                            skewness = 1;//1e-1;//1e+18;//0;//1e-10;//1e+10;//1e+13;//8;//18;//1e-18//5e+17;//1e-7;
                            eccentricity = 1;//1;//1e+18;//1e-18;//1e-18;//8;//0;//1e-15;//1e-8;
                            intersection = 0;//1;//1e-10;//1e-9;
                            uniformity = 0;//1e+12;//e-12;//1e-2;//1e-3;//1e-17;//17;//1e-17;//5e+13;//1e-8;
                            area = 1;//1e-8;//11;//1e-18;//1e+15;//1e+16;//1e-8;//0.00001;
                            length = 0;//1e+6;//1e+3;//1e-18;//1e+18;//8;//1e+13;//1e+17;//1e+18;//0.000002;
                            epsilon = 1e-7;
                            jacPts = 101 * 101;
                            dumped = false;
                            gsMatrix<> xy1 = geom->eval(uv1);
                            gsQualityMeasure2<real_t> optimization(*geom, uv1, xy1, basisInd, false);
                            for (int it = 0; it != iterations; it++)
                            {
//                                               gsInfo
                                outfile << "iteration: " << it << " / " << iterations - 1 << "\n";
                                optimization.optimize(fitting, orthogonality, tim, skewness,
                                                      eccentricity, uniformity,
                                                      length, area,
                                                      intersection, epsilon, dumped);
                                //gsInfo
                                outfile << "it = " << it << " Value of functional: "
                                        << optimization.functional(fitting, orthogonality, tim, skewness,
                                                                   eccentricity, uniformity,
                                                                   length, area,
                                                                   intersection, epsilon)
                                        << "\n";
                                saveData(*geom, output, it + 1);
                                minusnumber = checkJacobianDeterminant(*geom, jacPts, true, output, it + 1);
                                // Giving different names for files with (possibly) final results
                                outfile << minusnumber << "\n";
                                saveData(*geom, output + "_final", 0);
                                if (minusnumber == 0)
                                {
                                    gsInfo << "the parametrization become regular\n";
                                    outfile << "the parametrization become regular\n";
                                    checkAccuracySufficiency(*geom,
                                                             *geomForCheck,
                                                             maxError,
                                                             boundaryerror,
                                                             au,
                                                             bu,
                                                             av,
                                                             bv);
                                    if (epsilon_g > maxError && epsilon_f > boundaryerror)
                                    {
//                                    for (int therow = 0; therow < lastNonzeroRow + 1; therow++)
//                                    {
//                                        //                gsInfo
//                                        outfile << therow << "; CREATING THE TMP SPLINE ";
//                                        std::vector<unsigned> box;
//                                        for (int column = 0; column < 5; column++)
//                                        {
//                                            box.push_back(boxMat[therow][column]);
//                                            //                    gsInfo
//                                            outfile << boxMat[therow][column] << "; ";
//                                        }
//                                        //                                THBFromGeo->refineElements(box);
//                                        THBAccepted.refineElements(box);
//                                        box.clear();
//                                        //                gsInfo
//                                        outfile << "\n";
//                                    }
                                        gsInfo << "NLO WORKED!\n";
                                        outfile << "NLO WORKED!\n";
                                        success = 1;
                                        successfullAttempts++;
                                        totalAttempts++;
                                        nlos++;
                                        setSupports(THB, supps);
                                        anmat2 = THB.anchors();
                                        gsInfo << "SIZE: " << THB.size() << "\n";
                                        outfile << "SIZE: " << THB.size() << "\n";
                                        acceptedsize = THB.size();
                                        gsWrite(*fittedFunction,
                                                xmlFile);//MAYBE, THERE IS A WAY OF IMPROVING
                                        wasRebuilt = 1;
                                        nonCheckedCells.removeElement(currArrayIndex);
                                        break;
                                    }
                                }
                            }
                            if (minusnumber != 0 || epsilon_g < maxError || epsilon_f < boundaryerror)
                            {
                                gsInfo << "the errors after NLO: " << maxError << " " << boundaryerror << "\n";
                                outfile << "the errors after NLO: " << maxError << " " << boundaryerror << "\n";
                                gsInfo << "WITHDRAW\n";
                                outfile << "WITHDRAW\n";
                                totalAttempts++;
                                if(minusnumber!=0){
                                    gsWrite(*geom, badFile);
                                }
                                restoreTheHierarchy(createdBoxNum,
                                                    lastNonzeroRow,
                                                    boxMat,
                                                    levNow,
                                                    centerInd,
                                                    ourBox,
                                                    successfullAttempts);
                                continue;
//                    return 0;
                            }

                        }


                        ////                        gsInfo
                        //                        outfile << output << "_final\n";
                    }


                    nonCheckedCells.removeElement(currArrayIndex); //REWRITE MAKING THIS AUTOMATIC DEPENDING ON THE METHOD

                }
            }

        }

    }
    gsInfo << "SIZE: " << acceptedsize << "\n";
    outfile << "SIZE: " << acceptedsize << "\n";
    gsInfo << "proj: " << proj  << ", linear: " << los << ", nonlinear: " << nlos << "\n";
    outfile << "proj: " << proj  << ", linear: " << los << ", nonlinear: " << nlos << "\n";
    std::string geomParaviewFile = givenGeo + acCond + "_Mesh_0";
    gsWriteParaview(*fittedFunction, geomParaviewFile);

    //    gsFileManager::open("results_optimized_final_Mesh_0.pvd");
    outfile << "totalAttempts: " << totalAttempts << "\n";
    outfile << "successfullAttempts: " << successfullAttempts << "\n";
    outfile << "SUCCESS RATE: " << 100 * successfullAttempts / totalAttempts << "\n";
    outfile << "You can find the geometry in " + givenGeo + acCond + "_Mesh_0.pvd\n";
    gsInfo << "You can find the geometry in " + givenGeo + acCond + "_Mesh_0.pvd ; DON'T FORGET .vtp file!\n";
    //    gsFileManager::open("accepted_Mesh_0.pvd");
    //    gsFileManager::open(givenGeo + acCond + "_Mesh_0.pvd");
    outfile << "Total execution time:" << clock.stop() << "\n";
    gsInfo  << "Total execution time:" << clock.stop() << "\n";
    saveData(*fittedFunction, pvdFile, 0);
    gsInfo << gradingExtent;
    outfile << gradingExtent;
    return 0;
}
