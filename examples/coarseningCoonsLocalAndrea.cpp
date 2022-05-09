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

int numRed, theLev;
int redBox[4096][5]; //LAZY SOLUTION
ofstream outfile;

void boxToDomain(int mybox[][5], double coords[][5], index_t interior) {
    for (int i = 0; i < numRed; i++) {
        coords[i][0] = (double) mybox[i][0];
        coords[i][1] = (double) mybox[i][1] / ((double)((interior + 1)* pow(2, mybox[i][0])));
        coords[i][2] = (double) mybox[i][2] / ((double)((interior + 1)* pow(2, mybox[i][0])));
        coords[i][3] = (double) mybox[i][3] / ((double)((interior + 1)* pow(2, mybox[i][0])));
        coords[i][4] = (double) mybox[i][4] / ((double)((interior + 1)* pow(2, mybox[i][0])));
    }
}

void boxToDomain(int mybox[5], double coords[5], index_t interior) {
    coords[0] = (double) mybox[0];
    coords[1] = (double) mybox[1] / ((double)((interior + 1)* pow(2, mybox[0])));
    coords[2] = (double) mybox[2] / ((double)((interior + 1)* pow(2, mybox[0])));
    coords[3] = (double) mybox[3] / ((double)((interior + 1)* pow(2, mybox[0])));
    coords[4] = (double) mybox[4] / ((double)((interior + 1)* pow(2, mybox[0])));

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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void updateCriterion(int levNow, int x1U, int y1U, int x2U, int y2U, int cellExtension, gsMatrix<int> & basisInd, gsTHBSplineBasis<2, double> THB, index_t interior){
    real_t cellCoords[5];
    int cellIndices[5];
    cellIndices[0] = levNow;
    cellIndices[1] = std::min((interior + 1)* (int)pow(2,levNow), std::max(x1U, 0));
    cellIndices[2] = std::min((interior + 1)* (int)pow(2,levNow),std::max(y1U, 0));
    cellIndices[3] = std::max(0,std::min(x2U, (interior + 1)* (int)pow(2,levNow)));
    cellIndices[4] = std::max(0,std::min(y2U, (interior + 1)* (int)pow(2,levNow)));
//    gsInfo << cellIndices[0] << "\t" << cellIndices[1] << "\t" << cellIndices[2] << "\t" << cellIndices[3] << "\t" << cellIndices[4] << "\n";
    boxToDomain(cellIndices, cellCoords, interior);
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
    for (int i = 0; i < vectorGaussPoints.cols(); ++i) {
//        gsInfo << vectorGaussPoints(0,i) << ";" << vectorGaussPoints(1,i) << "\n";
    }

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
//            gsWriteParaview(thbLocal, exitBasis);
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void checkAccuracySufficiency(gsGeometry < > & geom, gsGeometry < > & geomForCheck, real_t & error, real_t & boundaryerror, real_t & au, real_t & bu, real_t & av, real_t & bv, int fineLevel, index_t interior) {
    int numberOfBoundaryPoints = 100, numCols = 19, hsize = 0;
    int coarseLevel = fineLevel - 1;
    double coords[(interior + 1)* (int) pow(2, coarseLevel) * (interior + 1)* (int) pow(2, coarseLevel)][5];
    boxToDomain(redBox, coords, interior);

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
        if (redBox[i][4] == (interior + 1)* pow(2, coarseLevel + 1)) {
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
        if (redBox[i][3] == (interior + 1)* pow(2, coarseLevel + 1)) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                hpoint(0, k) = bu;
                hpoint(1, k) = coords[i][2] + (coords[i][4] - coords[i][2]) * ((double) j / numberOfBoundaryPoints);
                k++;
            }
        }
        //            else { gsInfo << redBox[i][3] << " neq " << (interior + 1)* pow(2, fineLevel); }
        if (redBox[i][4] == (interior + 1)* pow(2, coarseLevel + 1)) {
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int pickCell(gsVector < int > vectorS, int & currArrayIndex, int levNow, int & x1U, int & y1U, int & x2U, int & y2U, index_t interior) {
    currArrayIndex = rand() % vectorS.size();
    gsInfo << "nonCheckedCells.size() " << vectorS.size() << "\n";
    int currCellIndex = vectorS(currArrayIndex);
    x1U = ((interior + 1)* (int) pow(2, levNow) - 1) - (currCellIndex % (int)((interior + 1)* (int) pow(2, levNow)));
    x2U = x1U + 1;
    y1U = ((interior + 1)* (int) pow(2, levNow) * (interior + 1)* (int) pow(2, levNow) - 1 - currCellIndex) / ((interior + 1)* (int) pow(2, levNow));
    y2U = y1U + 1;
    return currCellIndex;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int pickCell(gsVector < int > vectorS, int attempt, int levNow, int & x1U, int & y1U, int & x2U, int & y2U, int lexicographic, index_t interior) {
    index_t currArrayIndex = 0;//attempt % nonCheckedCells.size();
    gsInfo << "nonCheckedCells.size() " << vectorS.size() << "\n";
    int currCellIndex = vectorS(currArrayIndex);
    x1U = ((interior + 1)* (int) pow(2, levNow) - 1) - (currCellIndex % (int)((interior + 1)* (int) pow(2, levNow)));
    x2U = x1U + 1;
    y1U = ((interior + 1)* (int) pow(2, levNow) * (interior + 1)* (int) pow(2, levNow) - 1 - currCellIndex) / ((interior + 1)* (int) pow(2, levNow));
    y2U = y1U + 1;
    return currCellIndex;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_t getSupports(int i, int j, int k, gsMatrix<> supps){
    return supps(j, 2*i + k);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char * argv[]) {
    gsStopwatch clock;
    const int dim =2;
    int row, acceptedsize, attempt = 0, fineLevel, coarseLevel;
//    std::string xmlFile = "/home/turing/theydarov/xmlFiles/resultLocalAndrea";
    std::string badFile = "/home/turing/theydarov/xmlFiles/badgeoLocalAndrea";
    std::string pvdFile = "/home/turing/theydarov/pvdFiles/resultLocalAndrea";
    int los = 0, nlos = 0, proj = 0;
    std::string givenGeo;
    int nonValid;
    int gradingExtent;
    real_t epsilon_g = 5e-3, epsilon_f = 5e-3;
    real_t lcx,lcy,ucx,ucy;
    std::string acCond = to_string(epsilon_g) + "and" + to_string(epsilon_f);
//    std::string inputInitial("/home/turing/theydarov/geometries/WigglyCoonsPatchIntDist.xml");
    std::string inputInitial("/home/turing/theydarov/geometries/jigsawpuzzleTP.xml");
//    std::string inputInitial("/home/turing/theydarov/output/jigsawpuzzleTPL2LONLO0.005000and0.005000.xml");
//    std::string inputInitial("/home/turing/theydarov/geometries/jigsawpuzzleTPEB.xml");
//    std::string inputInitial("/home/turing/theydarov/geometries/WigglySquareIntDist.xml");
//    std::string inputInitial("/home/turing/theydarov/geometries/WigglyCoonsPatch3lev.xml");
//    std::string inputInitial("/home/turing/theydarov/geometries/indiana_thb.xml");
    givenGeo = "jigsawpuzzleTP";
//    givenGeo = "jigsawpuzzleTPEB";
    //    givenGeo = "wigglysurface";
//        givenGeo = "Indiana";
//    givenGeo+=  "WigglyCoonsIntDist";
    //    givenGeo = "Indiana";
//    givenGeo+=  "WigglySquareIntDist";
//    givenGeo+=  "WigglyCoonsPatch3lev";
//    givenGeo+= "lexicographic";
//    givenGeo+= "semi-random";
//    givenGeo+="upperhalf";
    givenGeo += "L2";
    givenGeo += "LO";
    givenGeo += "NLO";
    std::string fileLoc = "/home/turing/theydarov/txtFiles/resultLocalAndrea";
    outfile.open(fileLoc + ".txt");
    std::string inputPtsParams("/home/turing/theydarov/gismo/extensions/motor/filedata/jku/thb_parameters_and_points.xml");
    gsFileData < > data0(inputInitial);
    int iter = -1;
    int successfullAttempts = 0, totalAttempts = 0;
    gsInfo << "GIVE A NUMBER OF RED BOXES\n";
    outfile << "GIVE A NUMBER OF RED BOXES\n";
    gsInfo << "I HOPE YOU HAVE /home/theydarov/gismo_unsupported/gismo/extensions/motor/filedata/" + givenGeo + acCond + ".xml\n";
    outfile << "I HOPE YOU HAVE /home/theydarov/gismo_unsupported/gismo/extensions/motor/filedata/" + givenGeo + acCond + ".xml\n";
    gsGeometry<>::uPtr geomForCheck;
    gsTHBSplineBasis < 2, real_t > ::uPtr THBFromGeo;
    THBFromGeo = data0.getAnyFirst < gsTHBSplineBasis < 2 >> ();
    unsigned interior = THBFromGeo->numKnots(0,1) - 2*(THBFromGeo->degree(1) + 1);
    acceptedsize = THBFromGeo->size();
    gsMatrix < index_t > lowCorners;
    gsMatrix < index_t > upCorners;
    gsVector < index_t > myLevel;
    cin >> numRed;
    gsInfo << "SELECT PICKCELL METHOD:\n";
    char method;
    cin >> method;
    gsInfo  <<  "METHOD: " << method << "\n";
    outfile <<  "METHOD: " << method << "\n";
    gsGeometry < > ::uPtr geom;
    gsMatrix < > coefs;
    if (data0.has < gsGeometry < > > ()) {
        geom = data0.getFirst < gsGeometry < > > ();
        geomForCheck = data0.getFirst<gsGeometry<> >();
        coefs = geom -> coefs();
    }
    (THBFromGeo -> tree().getBoxes(lowCorners, upCorners, myLevel));
    fineLevel = myLevel.maxCoeff();
    coarseLevel = fineLevel - 1;
    int numCells = (interior + 1)* (int) pow(2, coarseLevel) * (interior + 1)* (int) pow(2, coarseLevel);
    outfile << "I REQUEST " << numCells << " cells\n";
    if (numRed == 0) {
        redBox[0][0] = fineLevel;
        redBox[0][1] = 0;
        redBox[0][2] = 0;
        redBox[0][3] = (interior + 1)* (int) pow(2, fineLevel);
        redBox[0][4] = (interior + 1)* (int) pow(2, fineLevel);
        outfile << "I ASSUME THAT REDBOX HAS THE COORDINATES:\n"
                <<
                redBox[0][0] << " " << redBox[0][1] << " " << redBox[0][2] << " " << redBox[0][3] << " " << redBox[0][4] << "\n";
    }
    real_t boundaryerror = 0;
    for (int i = 0; i < numRed; i++) {
        gsInfo << "GIVE RED BOX\n";
        cin >> redBox[i][0] >> redBox[i][1] >> redBox[i][2] >> redBox[i][3] >> redBox[i][4];
    }
    gsInfo << "Give the grading extent:\n";
    cin >> gradingExtent;
    numRed = 1;
    int RTH = 0, centerInd = 0;
    int lastNonzeroRow = 0; //71;//92;
    int minusnumber, createdBoxNum = 0, currCellIndex, wasRebuilt;
    int x1U, x2U, y1U, y2U; // = rand()%32;
    int x1Bi, x2Bi, y1Bi, y2Bi;
    int needToEscape = 0;
    int createSpline, failed;
    int ourBox[5], currArrayIndex;
    currArrayIndex = 0;
    std::string output("results_optimized");
    gsGeometry < > ::uPtr geomhat; //accepted geometry
    gsMatrix < > coefshat;
    coefshat = coefs;
    gsMatrix<> anmat2;
    gsMatrix<> anmat ;
    THBFromGeo->anchors_into(anmat2);
    gsMatrix<> supps(THBFromGeo->support().rows(),THBFromGeo->support().cols()+ 2*THBFromGeo->size());
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
    int degree = THBFromGeo->maxDegree();
    unsigned multEnd = degree + 1; // multiplicity at the two end knots
    gsKnotVector <> ku(au, bu, interior, multEnd);
    gsKnotVector <> kv(av, bv, interior, multEnd);
    gsTensorBSplineBasis < 2, real_t > tens(ku, kv);
    gsTHBSplineBasis < 2, real_t > THBAccepted(tens);
    gsTHBSplineBasis < 2, real_t > THBTemporary(tens);
    wasRebuilt = 0;
    real_t maxError = 0;
    THBAccepted = *THBFromGeo;
    gsGeometry<>::uPtr fittedFunction;
    int iteration = 0;
    std::string xmlFile = "/home/turing/theydarov/output/" + givenGeo + acCond;
    int fullMonty;
//    for (int levNow = coarseLevel; levNow >= 0; levNow--)
    for (int levNow = coarseLevel; levNow >= coarseLevel - 1; levNow--)
    {
        iteration = 0;
        theLev = levNow;
        if (levNow < coarseLevel && wasRebuilt)
        {
            std::string input(xmlFile + ".xml" );
            gsFileData < > data1(input);
            THBAccepted = *data1.getAnyFirst<gsTHBSplineBasis<2 >>();
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
        outfile << "working with level " << levNow << " as a coarseLevel\n";
        int pickedOne = -1;
        gsVector<int> nonCheckedCells((interior + 1)* pow(2, levNow) * (interior + 1)* pow(2, levNow));
        gsMatrix<int> pickedCells(1, (interior + 1)* pow(2, levNow) * (interior + 1)* pow(2, levNow));
        for(int i = 0; i < nonCheckedCells.size(); i++){
            nonCheckedCells(i) = i;
            pickedCells(0, i) = 0;
        }
        int success = 1;
        attempt = 0;
        while(nonCheckedCells.size() != 0 && success)
        {
            iteration++;
            success = 0;
            gsVector<index_t> vectorS = nonCheckedCells;
            fullMonty = 0;
            while(vectorS.size() != 0) {
                gsInfo << "vectorS.size(): " << vectorS.size() << "\n";
                lcx = 1.0;
                lcy = 1.0;
                ucx = 0.0;
                ucy = 0.0;
                failed = 0;
                outfile << "\n";
                outfile << "\n";
                outfile << "\n";
                outfile << "The boxes\n";
                for (int i = 0; i <= lastNonzeroRow; i++) {
                    for (int j = 0; j < 5; j++) {
                        outfile << boxMat[i][j] << "\t";
                    }
                    outfile << "\n";
                }
                createdBoxNum = 0;
                if(method == 'r'){
                    currCellIndex = pickCell(vectorS, currArrayIndex, levNow, x1U, y1U, x2U, y2U, interior);
                }
                else if(method == 'l'){
                    currCellIndex = pickCell(vectorS, attempt, levNow, x1U, y1U, x2U, y2U, 1, interior);
                }
                else if(method == 's'){
                    if(success && !pickedCells(0, pickedOne + 1)) {
                        pickCell(vectorS, attempt, levNow, x1U, y1U, x2U, y2U, 1, interior);
                    }
                    else{
                        currCellIndex = pickCell(vectorS, currArrayIndex, levNow, x1U, y1U, x2U, y2U, interior);
                    }
                    pickedOne = currCellIndex;
                    pickedCells(0, currCellIndex) = 1;
                }
                else{
                    gsInfo << "UNKNOWN METHOD.\n";
                    return 1;
                }
                int jopa = 4*(int)pow(2, levNow)*4*(int)pow(2, levNow);
                outfile << "Number of nonchecked cells: " << vectorS.size() << "\n";
                gsInfo << "attempt " << attempt << ", CURRENT INDEX IS " << vectorS(currArrayIndex) <<
                       ",\t the coordinates of the box are\n" << x1U << " " << y1U << " " << x2U << " " << y2U << "\n";
                outfile << "attempt " << attempt << ", CURRENT INDEX IS " << vectorS(currArrayIndex) <<
                        ",\t the coordinates of the box are\n" << x1U << " " << y1U << " " << x2U << " " << y2U << "\n";
                gsInfo << currArrayIndex << "\n";
                vectorS.removeElement(currArrayIndex);
                gsInfo << "DELETED\n";
                attempt = (attempt + 1)%(jopa);
                createSpline = 0;
                for (int currentrow = 0; currentrow < lastNonzeroRow + 1; currentrow++)
                {
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
                if (createSpline == 1)
                {
                    outfile << "creating the spline\n";
                    gsInfo  << "creating the spline\n";
                    gsTHBSplineBasis<2, real_t> THB(tens);
                    gsInfo << THB << "\n";
                    for (int therow = 0; therow < lastNonzeroRow + 1; therow++)
                    {
                        outfile << therow << "; ";
                        std::vector<index_t> box;
                        for (int column = 0; column < 5; column++)
                        {
                            box.push_back(boxMat[therow][column]);
                            outfile << boxMat[therow][column] << "; ";
                        }
                        THB.refineElements(box);
                        box.clear();
                        outfile << "\n";
                    }
                    (THB.tree().getBoxes(lowCorners, upCorners, myLevel));
                    for (int i = 0; i < lowCorners.rows(); i++) {
                        boxMat[i][0] = myLevel(i);
                        boxMat[i][1] = (real_t) lowCorners(i, 0) / pow(2, THB . maxLevel() - myLevel(i));
                        boxMat[i][2] = (real_t) lowCorners(i, 1) / pow(2, THB . maxLevel() - myLevel(i));
                        boxMat[i][3] = (real_t) upCorners(i, 0) / pow(2, THB . maxLevel() - myLevel(i));
                        boxMat[i][4] = (real_t) upCorners(i, 1) / pow(2, THB . maxLevel() - myLevel(i));
                        gsInfo << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                               boxMat[i][4] << "\n";
                        outfile << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                                boxMat[i][4] << "\n";
                    }
                    lastNonzeroRow = lowCorners.rows() - 1;
                    for (int therow = 0; therow < lastNonzeroRow + 1; therow++)
                    {
                        std::vector<index_t> box;
                        for (int column = 0; column < 5; column++)
                        {
                            box.push_back(boxMat[therow][column]);
//                            outfile << boxMat[therow][column] << "; ";
                        }
                        THB.refineElements(box);
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
                            gsInfo << "The basis is inadequate, rewriting\n";
                            std::string input(xmlFile+".xml" );
                            gsFileData < > data2(input);
                            THBTemporary = *data2.getAnyFirst<gsTHBSplineBasis<2 >>();
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
                        attempt = std::max(attempt -1 , 0);
                        continue;
                    }
                    outfile << THB << "\n";
                    gsInfo << THB << "\n";
                    THB.anchors_into(anmat);
                    gsMatrix<int> basisInd(1, THB.size()); //1 if function is not presented
                    for (int l = 0; l < THB.size(); ++l)
                    {
                        basisInd(0, l) = 0;
                    }
                    gsMatrix<> THBcoefs(THB.size(), 2);
                    for (int k = 0; k < THB.size(); ++k)
                    {
                        THBcoefs(k,0) = THBcoefs(k,1) = 0.0;
                    }
                    for (int j = 0; j < THB.size(); j++)
                    {
                        for (int i = 0; i < anmat2.cols(); i++)
                        {
                            if (
                                    THB.support(j)(0, 0) == getSupports(i, 0, 0, supps)
                                    && THB.support(j)(1, 0) == getSupports(i, 1, 0, supps)
                                    && THB.support(j)(0, 1) == getSupports(i, 0, 1, supps)
                                    && THB.support(j)(1, 1) == getSupports(i, 1, 1, supps)
                                    )
                            {

                                THBcoefs(j, 0) = coefshat(i, 0);
                                THBcoefs(j, 1) = coefshat(i, 1);
                                basisInd(j) = 0;

                                break;
                            }
                        }
                    }
                    for (int m = 0; m < THB.size(); m++)
                    {
                        if(THBcoefs(m,0) == 0.0 || THBcoefs(m,1) == 0.0)
                        {
                            basisInd(0,m) = 1;
                        }
                    }
                    for (int i = -gradingExtent; i < gradingExtent; ++i) {
                        for (int j = -gradingExtent; j < gradingExtent; ++j) {
                            updateCriterion(levNow, x1U + i, y1U + j, x1U + i + 1, y1U + j + 1, 1, basisInd, THB, interior);
                        }
                    }
//                    gsInfo << basisInd << "\n";
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
                    checkAccuracySufficiency(*fittedFunction, *geomForCheck, maxError, boundaryerror, au, bu, av, bv, fineLevel, interior);
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
                        acceptedsize = THB.size();
                        for (int l = 0; l <= lastNonzeroRow; l++) {
                            for (int m = 1; m < 5; m++) {
                                outfile << boxMat[l][m];
                                outfile << "\t";
                                gsInfo << boxMat[l][m];
                                gsInfo << "\t";
                            }
                            outfile << "\n";
                            gsInfo << "\n";
                        }
                        proj++;
                        gsWrite(*fittedFunction,
                                xmlFile);//MAYBE, THERE IS A WAY OF IMPROVING
                        setSupports(THB, supps);
                        for (int n = 0; n < THB.size(); ++n)
                        {
                            gsInfo << supps(1,n) << "; ";
                        }
                        gsInfo << "\n";
                        anmat2 = THB.anchors();
                        wasRebuilt = 1;
                        gsInfo << "nonCheckedCells.size():" << nonCheckedCells.size() << "\n";
                        nonCheckedCells.removeElement(currArrayIndex);
                        gsInfo << "nonCheckedCells.size():" << nonCheckedCells.size() << "\n";
                        continue;
                    }
                    else
                    {
                        if(minusnumber != 0){    gsInfo << "NONREGULAR:" << minusnumber << "\n";
                            outfile << "NONREGULAR:" << minusnumber << "\n";}

                        else{
                            gsInfo << "WITHDRAW, CAUSE THE PARAMETRIZATION IS REGULAR\n";
                            outfile << "WITHDRAW, CAUSE THE PARAMETRIZATION IS REGULAR\n";
                            nonCheckedCells.removeElement(currArrayIndex); //REWRITE MAKING THIS AUTOMATIC DEPENDING ON THE METHOD
                            totalAttempts++;
                            gsInfo << maxError << "\t" << boundaryerror << "\n";
                            outfile << maxError << "\t" << boundaryerror << "\n";
//                            restoreTheHierarchy(createdBoxNum,
//                                                lastNonzeroRow,
//                                                boxMat,
//                                                levNow,
//                                                centerInd,
//                                                ourBox,
//                                                successfullAttempts);
                            std::vector<index_t> box;
                            box.push_back(ourBox[0]);
                            box.push_back(ourBox[1]);
                            box.push_back(ourBox[2]);
                            box.push_back(ourBox[3]);
                            box.push_back(ourBox[4]);
                            THB.refineElements(box);

                            (THB.tree().getBoxes(lowCorners, upCorners, myLevel));
                            for (int i = 0; i < lowCorners.rows(); i++) {
                                boxMat[i][0] = myLevel(i);
                                boxMat[i][1] = (real_t) lowCorners(i, 0) / pow(2, THB . maxLevel() - myLevel(i));
                                boxMat[i][2] = (real_t) lowCorners(i, 1) / pow(2, THB . maxLevel() - myLevel(i));
                                boxMat[i][3] = (real_t) upCorners(i, 0) / pow(2, THB . maxLevel() - myLevel(i));
                                boxMat[i][4] = (real_t) upCorners(i, 1) / pow(2, THB . maxLevel() - myLevel(i));
                                gsInfo << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                                       boxMat[i][4] << "\n";
                                outfile << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                                        boxMat[i][4] << "\n";
                            }
                            lastNonzeroRow = lowCorners.rows() - 1;
                            for (int l = 0; l <= lastNonzeroRow; l++) {
                                for (int m = 0; m < 5; m++) {
                                    outfile << boxMat[l][m];
                                    outfile << "\t";
                                    gsInfo << boxMat[l][m];
                                    gsInfo << "\t";
                                }
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
                        gsQualityMeasure2<real_t> optimization(*fittedFunction, uv1, xy1, basisInd, true);
                        for (int it = 0; it != iterations; it++)
                        {
                            outfile << "iteration: " << it << " / " << iterations - 1 << "\n";
                            optimization.optimize(fitting, orthogonality, tim, skewness,
                                                  eccentricity, uniformity,
                                                  length, area,
                                                  intersection, epsilon, dumped);
                            outfile << "it = " << it << " Value of functional: "
                                    << optimization.functional(fitting, orthogonality, tim, skewness,
                                                               eccentricity, uniformity,
                                                               length, area,
                                                               intersection, epsilon)
                                    << "\n";
                            minusnumber = checkJacobianDeterminant(*fittedFunction, jacPts, true, output, it + 1);
                            if (minusnumber == 0) break;
                        }
                        gsInfo << "Number of negative points " << minusnumber << "\n";
                        outfile << minusnumber << "\n";
                        checkAccuracySufficiency(*fittedFunction, *geomForCheck, maxError, boundaryerror, au, bu, av, bv, fineLevel, interior);
                        if(epsilon_g < maxError || epsilon_f < boundaryerror)
                        {
                            gsInfo << "NO NEED\n";
                            gsInfo << "WITHDRAW\n";
                            totalAttempts++;
                            gsInfo << "the errors after LO: " << maxError << " " << boundaryerror << "\n";
                            outfile << "the errors after LO: " << maxError << " " << boundaryerror << "\n";
                            outfile << "WITHDRAW\n";
                            if(minusnumber!=0){
                                gsWrite(*fittedFunction, "/home/turing/theydarov/auxiliary/badgeoLocalGrading.xml");

                            }
//                            restoreTheHierarchy(createdBoxNum,
//                                                lastNonzeroRow,
//                                                boxMat,
//                                                levNow,
//                                                centerInd,
//                                                ourBox,
//                                                successfullAttempts);
                            std::vector<index_t> box;
                            box.push_back(ourBox[0]);
                            box.push_back(ourBox[1]);
                            box.push_back(ourBox[2]);
                            box.push_back(ourBox[3]);
                            box.push_back(ourBox[4]);
                            THB.refineElements(box);
                            (THB.tree().getBoxes(lowCorners, upCorners, myLevel));
                            for (int i = 0; i < lowCorners.rows(); i++) {
                                boxMat[i][0] = myLevel(i);
                                boxMat[i][1] = (real_t) lowCorners(i, 0) / pow(2, THB . maxLevel() - myLevel(i));
                                boxMat[i][2] = (real_t) lowCorners(i, 1) / pow(2, THB . maxLevel() - myLevel(i));
                                boxMat[i][3] = (real_t) upCorners(i, 0) / pow(2, THB . maxLevel() - myLevel(i));
                                boxMat[i][4] = (real_t) upCorners(i, 1) / pow(2, THB . maxLevel() - myLevel(i));
                                gsInfo << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                                       boxMat[i][4] << "\n";
                                outfile << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                                        boxMat[i][4] << "\n";
                            }
                            lastNonzeroRow = lowCorners.rows() - 1;
                            continue;
                        }
                        if (minusnumber == 0 && !(epsilon_g < maxError || epsilon_f < boundaryerror))
                        {
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
                            gsInfo << "nonCheckedCells.size():" << nonCheckedCells.size() << "\n";
                            nonCheckedCells.removeElement(currArrayIndex);
                            gsInfo << "nonCheckedCells.size():" << nonCheckedCells.size() << "\n";
                            continue;
                        }
                        else
                        {
                            gsWrite(*fittedFunction, badFile);
                            std::string input = badFile + ".xml";
                            gsFileData<> data(input);
                            gsGeometry<>::uPtr geom;
                            if (data.has< gsGeometry<> >())
                            {
                                geom = data.getFirst< gsGeometry<> >();
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
                            gsQualityMeasure2<real_t> optimization(*fittedFunction, uv1, xy1, basisInd, true);
                            for (int it = 0; it != iterations; it++)
                            {
                                outfile << "iteration: " << it << " / " << iterations - 1 << "\n";
                                optimization.optimize(fitting, orthogonality, tim, skewness,
                                                      eccentricity, uniformity,
                                                      length, area,
                                                      intersection, epsilon, dumped);
                                outfile << "it = " << it << " Value of functional: "
                                        << optimization.functional(fitting, orthogonality, tim, skewness,
                                                                   eccentricity, uniformity,
                                                                   length, area,
                                                                   intersection, epsilon)
                                        << "\n";
                                minusnumber = checkJacobianDeterminant(*fittedFunction, jacPts, true, output, it + 1);
                                outfile << minusnumber << "\n";
                                if (minusnumber == 0)
                                {
                                    gsInfo << "the parametrization become regular\n";
                                    outfile << "the parametrization become regular\n";
                                    checkAccuracySufficiency(*fittedFunction,
                                                             *geomForCheck,
                                                             maxError,
                                                             boundaryerror,
                                                             au,
                                                             bu,
                                                             av,
                                                             bv, fineLevel, interior);

                                    if (epsilon_g > maxError && epsilon_f > boundaryerror)
                                    {
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
                                        gsInfo << "nonCheckedCells.size():" << nonCheckedCells.size() << "\n";
                                        nonCheckedCells.removeElement(currArrayIndex);
                                        gsInfo << "nonCheckedCells.size():" << nonCheckedCells.size() << "\n";
                                        fullMonty = 1;
                                        break;
                                    }
                                    gsInfo << "BREAKPOINT1\n";
                                }
                                gsInfo << "BREAKPOINT2\n";
                            }
                            if(fullMonty)   continue;
                            if (minusnumber != 0 || epsilon_g < maxError || epsilon_f < boundaryerror)
                            {
                                gsInfo << "the errors after NLO: " << maxError << " " << boundaryerror << "\n";
                                outfile << "the errors after NLO: " << maxError << " " << boundaryerror << "\n";
                                gsInfo << "WITHDRAW\n";
                                outfile << "WITHDRAW\n";
                                nonCheckedCells.removeElement(currArrayIndex); //REWRITE MAKING THIS AUTOMATIC DEPENDING ON THE METHOD
                                totalAttempts++;
                                if(minusnumber!=0){
                                    gsWrite(*geom, badFile);

                                }
//                                restoreTheHierarchy(createdBoxNum,
//                                                    lastNonzeroRow,
//                                                    boxMat,
//                                                    levNow,
//                                                    centerInd,
//                                                    ourBox,
//                                                    successfullAttempts);
                                std::vector<index_t> box;
                                box.push_back(ourBox[0]);
                                box.push_back(ourBox[1]);
                                box.push_back(ourBox[2]);
                                box.push_back(ourBox[3]);
                                box.push_back(ourBox[4]);
                                THB.refineElements(box);
                                (THB.tree().getBoxes(lowCorners, upCorners, myLevel));
                                for (int i = 0; i < lowCorners.rows(); i++) {
                                    boxMat[i][0] = myLevel(i);
                                    boxMat[i][1] = (real_t) lowCorners(i, 0) / pow(2, THB . maxLevel() - myLevel(i));
                                    boxMat[i][2] = (real_t) lowCorners(i, 1) / pow(2, THB . maxLevel() - myLevel(i));
                                    boxMat[i][3] = (real_t) upCorners(i, 0) / pow(2, THB . maxLevel() - myLevel(i));
                                    boxMat[i][4] = (real_t) upCorners(i, 1) / pow(2, THB . maxLevel() - myLevel(i));
                                    gsInfo << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                                           boxMat[i][4] << "\n";
                                    outfile << boxMat[i][0] << " " << boxMat[i][1] << " " << boxMat[i][2] << " " << boxMat[i][3] << " " <<
                                            boxMat[i][4] << "\n";
                                }
                                lastNonzeroRow = lowCorners.rows() - 1;
                                continue;
                            }
                        }
                    }
                    gsInfo << "BREAKPOINT1\n";
                    gsInfo << "nonCheckedCells.size():" << nonCheckedCells.size() << "\n";
                    nonCheckedCells.removeElement(currArrayIndex); //REWRITE MAKING THIS AUTOMATIC DEPENDING ON THE METHOD
                    gsInfo << "nonCheckedCells.size():" << nonCheckedCells.size() << "\n";
                }
            }
        }
    }
    gsInfo << "SIZE: " << acceptedsize << "\n";
    outfile << "SIZE: " << acceptedsize << "\n";
    gsInfo << "proj: " << proj  << ", linear: " << los << ", nonlinear: " << nlos << "\n";
    outfile << "proj: " << proj  << ", linear: " << los << ", nonlinear: " << nlos << "\n";
    std::string geomParaviewFile = givenGeo + acCond + "_Mesh_0";
    outfile << "totalAttempts: " << totalAttempts << "\n";
    outfile << "successfullAttempts: " << successfullAttempts << "\n";
    outfile << "SUCCESS RATE: " << 100 * successfullAttempts / totalAttempts << "\n";
    outfile << "You can find the geometry in " + givenGeo + acCond + "_Mesh_0.pvd\n";
    gsInfo << "You can find the geometry in " + givenGeo + acCond + "_Mesh_0.pvd ; DON'T FORGET .vtp file!\n";
    outfile << "Total execution time:" << clock.stop() << "\n";
    gsInfo  << "Total execution time:" << clock.stop() << "\n";
    saveData(*fittedFunction, pvdFile, 0);
    gsInfo << gradingExtent << "\n";
    outfile << gradingExtent;
    time_t tm;
    time(&tm);
    printf("Current Date/Time = %s", ctime(&tm));
    return 0;
}