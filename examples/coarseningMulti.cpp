#include <iostream>

#include <iostream>

#include <gismo.h>

#include "math.h"

#include <vector>

#include <gsIO/gsIOUtils.h>

#include "/home/turing/theydarov/gismo/extensions/motor/jku/gsQualityMeasure2.h"

#include "/home/turing/theydarov/gismo/extensions/motor/jku/gsMotorUtils.h"

#include <gsLocalFitting/gsLocalFitting.h>
//#include <gsLocalFitting/gsLocalHBFitting.h>
#include <gsLocalFitting/gsLocalTHBFitting.h>

using namespace gismo;
using namespace std;

int numRed, fineLevel = 2, coarseLevel = 1, theLev;
//int numCells = 4 * (int) pow(2, coarseLevel) * 4 * (int) pow(2, coarseLevel);
int redBox[4096][5]; //LAZY SOLUTION
ofstream outfile;





void boxToDomain(int mybox[][5], double coords[][5]) {
    for (int i = 0; i < numRed; i++) {
        coords[i][0] = (double) mybox[i][0];
        coords[i][1] = (double) mybox[i][1] / ((double)(4 * pow(2, mybox[i][0])));
        coords[i][2] = (double) mybox[i][2] / ((double)(4 * pow(2, mybox[i][0])));
        coords[i][3] = (double) mybox[i][3] / ((double)(4 * pow(2, mybox[i][0])));
        coords[i][4] = (double) mybox[i][4] / ((double)(4 * pow(2, mybox[i][0])));
//        coords[i][1] = (double) mybox[i][1] / ((double)(1 * pow(2, mybox[i][0])));
//        coords[i][2] = (double) mybox[i][2] / ((double)(1 * pow(2, mybox[i][0])));
//        coords[i][3] = (double) mybox[i][3] / ((double)(1 * pow(2, mybox[i][0])));
//        coords[i][4] = (double) mybox[i][4] / ((double)(1 * pow(2, mybox[i][0])));
    }
}

void boxToDomain(int mybox[5], double coords[5]) {
    coords[0] = (double) mybox[0];
    coords[1] = (double) mybox[1] / ((double)(4 * pow(2, mybox[0])));
    coords[2] = (double) mybox[2] / ((double)(4 * pow(2, mybox[0])));
    coords[3] = (double) mybox[3] / ((double)(4 * pow(2, mybox[0])));
    coords[4] = (double) mybox[4] / ((double)(4 * pow(2, mybox[0])));
//    coords[1] = (double) mybox[1] / ((double)(1 * pow(2, mybox[0])));
//    coords[2] = (double) mybox[2] / ((double)(1 * pow(2, mybox[0])));
//    coords[3] = (double) mybox[3] / ((double)(1 * pow(2, mybox[0])));
//    coords[4] = (double) mybox[4] / ((double)(1 * pow(2, mybox[0])));

}

void saveData(const gsGeometry < > & geom,
              const std::string & output
) {
    gsFileData < > fileData;
    fileData << geom;

    std::string out = output + "_Map_" + ".xml";
    fileData.dump(out);
    //    gsInfo
    outfile << "geometry saved to " << out << "\n";
    gsMesh < > mesh;
    makeMesh(geom.basis(), mesh, 1600);
    geom.evaluateMesh(mesh);
    out = output + "_Mesh_" ;
    gsWriteParaview(mesh, out);
    gsInfo << "Saved to " << out << "\n";
    outfile << "Saved to " << out << "\n";
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void updateCriterion(int levNow, int x1U, int y1U, int x2U, int y2U, int cellExtension, gsMatrix<real_t> & basisInd, gsTHBSplineBasis<2, double> THB){
    real_t cellCoords[5];
    int cellIndices[5];
    cellIndices[0] = levNow;
    cellIndices[1] = std::min(4 * (int)pow(2,levNow), std::max(x1U, 0));
    cellIndices[2] = std::min(4 * (int)pow(2,levNow),std::max(y1U, 0));
    cellIndices[3] = std::max(0,std::min(x2U, 4 * (int)pow(2,levNow)));
    cellIndices[4] = std::max(0,std::min(y2U, 4 * (int)pow(2,levNow)));
//    cellIndices[1] = std::min(1 * (int)pow(2,levNow), std::max(x1U, 0));
//    cellIndices[2] = std::min(1 * (int)pow(2,levNow),std::max(y1U, 0));
//    cellIndices[3] = std::max(0,std::min(x2U, 1 * (int)pow(2,levNow)));
//    cellIndices[4] = std::max(0,std::min(y2U, 1 * (int)pow(2,levNow)));
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


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void checkAccuracySufficiency(gsGeometry < > & geom, gsGeometry < > & geomForCheck, real_t & error, real_t & boundaryerror, real_t & au, real_t & bu, real_t & av, real_t & bv) {
    int numberOfBoundaryPoints = 400, numCols = 38, hsize = 0;
    double coords[4 * (int) pow(2, coarseLevel) * 4 * (int) pow(2, coarseLevel)][5];
    boxToDomain(redBox, coords);

    int numberOfInnerPoints = numCols * numCols;
    gsMatrix <real_t> innerPoints(2, numberOfInnerPoints);
    gsMatrix <real_t > innerVal(2, numberOfInnerPoints);
    gsMatrix <real_t > innerFittedVal(2, numberOfInnerPoints);

    gsMatrix <real_t > punto(2, 1);
    gsMatrix <real_t > qiymet(2, 1);

    for (int i = 0; i < numCols; i++) {
        for (int j = 0; j < numCols; j++) {
            innerPoints(0, i * numCols + j) = au + (bu - au) * ((double)(i + 1) / (numCols + 1));
            innerPoints(1, i * numCols + j) = av + (bu - au) * ((double)(i + 1) / (numCols + 1));
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
//        if (redBox[i][3] == 4 * pow(2, coarseLevel + 1)) {
        if (redBox[i][3] == 1 * pow(2, coarseLevel + 1)) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                //                    hpoint(0, hsize) = bu;
                //                    hpoint(1, hsize) = coords[i][2] + (coords[i][4] - coords[i][2]) * ((double) j / numberOfBoundaryPoints);
                hsize++;
            }
        }
        //            else{gsInfo << redBox[i][3] << " neq " << 4*pow(2,fineLevel);}
//        if (redBox[i][4] == 4 * pow(2, coarseLevel + 1)) {
        if (redBox[i][4] == 1 * pow(2, coarseLevel + 1)) {
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
//        if (redBox[i][3] == 4 * pow(2, coarseLevel + 1)) {
        if (redBox[i][3] == 1 * pow(2, coarseLevel + 1)) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                hpoint(0, k) = bu;
                hpoint(1, k) = coords[i][2] + (coords[i][4] - coords[i][2]) * ((double) j / numberOfBoundaryPoints);
                k++;
            }
        }
        //            else { gsInfo << redBox[i][3] << " neq " << 4 * pow(2, fineLevel); }
//        if (redBox[i][4] == 4 * pow(2, coarseLevel + 1)) {
        if (redBox[i][4] == 1 * pow(2, coarseLevel + 1)) {
            for (int j = 0; j < numberOfBoundaryPoints + 1; j++) {
                //                gsInfo << coords[0][1] << ", " << coords[0][3] << "\n";
                hpoint(0, k) = coords[i][1] + (coords[i][3] - coords[i][1]) * ((double) j / numberOfBoundaryPoints);
                hpoint(1, k) = bv;
                k++;
            }
        }
    }

    //    gsInfo << hsize << "\n";
    gsMatrix < real_t> hVal(2, hsize);
    gsMatrix < real_t> hfittedVal(2, hsize);

    //Compute geometries

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
    gsMatrix < real_t> errorsVector;
    gsMatrix < real_t> herrorsVector;
    errorsVector = (innerVal - innerFittedVal).cwiseAbs();
    herrorsVector = (hVal - hfittedVal).cwiseAbs();
    gsInfo << herrorsVector << "\n";
    error = errorsVector.maxCoeff();
    boundaryerror = herrorsVector.maxCoeff();
    //    for(int i=0;i < numberOfBoundaryPoints+1;i++){
    //        gsInfo << "i = " << i << " " << lowVal(1,i) << " " << lowfittedVal(1,i) << "\n";
    //    }
    ////Max componentwise error
    cout << "The max error is " <<error << "\n";
    cout << "Epsilon f is" << boundaryerror << "\n";

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



int pickCell(gsVector < int > nonCheckedCells, int & currArrayIndex, int levNow, int & x1U, int & y1U, int & x2U, int & y2U) {
    currArrayIndex = rand() % nonCheckedCells.size();
    int currCellIndex = nonCheckedCells(currArrayIndex);
    x1U = (4 * (int) pow(2, levNow) - 1) - (currCellIndex % (int)(4 * (int) pow(2, levNow)));
//    x1U = (1 * (int) pow(2, levNow) - 1) - (currCellIndex % (int)(1 * (int) pow(2, levNow)));
    x2U = x1U + 1;
    y1U = (4 * (int) pow(2, levNow) * 4 * (int) pow(2, levNow) - 1 - currCellIndex) / (4 * (int) pow(2, levNow));
//    y1U = (1 * (int) pow(2, levNow) * 1 * (int) pow(2, levNow) - 1 - currCellIndex) / (1 * (int) pow(2, levNow));
    y2U = y1U + 1;
    return currCellIndex;
}

int pickCell(int attempt, int levNow, int & x1U, int & y1U, int & x2U, int & y2U) {
    int currCellIndex = attempt;
    x1U = (4 * (int) pow(2, levNow) - 1) - (currCellIndex % (int)(4 * (int) pow(2, levNow)));
//    x1U = (1 * (int) pow(2, levNow) - 1) - (currCellIndex % (int)(1 * (int) pow(2, levNow)));
    x2U = x1U + 1;
    y1U = (4 * (int) pow(2, levNow) * 4 * (int) pow(2, levNow) - 1 - currCellIndex) / (4 * (int) pow(2, levNow));
//    y1U = (1 * (int) pow(2, levNow) * 1 * (int) pow(2, levNow) - 1 - currCellIndex) / (1 * (int) pow(2, levNow));
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
                gsInfo<< "ALERT: OUR CELL OF LEVEL " << levNow << "\n" << x1U << " " << y1U << " " << x2U <<
                        " " << y2U <<
                        " INTERSECTS WITH [" << row << "]th box of level " <<
                        boxMat[row][0] << ":\n" << boxMat[row][1] << " " << boxMat[row][2] << " " <<
                        boxMat[row][3] << " " << boxMat[row][4] << "\n";
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
                        gsInfo  << "x1Bi: " << x1Bi << "\t" << y1Bi << "\t"<< x2Bi << "\t"<< y2Bi << "\n";
            outfile  << "x1Bi: " << x1Bi << "\t" << y1Bi << "\t"<< x2Bi << "\t"<< y2Bi << "\t";
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
                    gsInfo << boxMat[row][0] << std::endl;
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
                                gsInfo << "CREATE E BOX:\n";
                outfile << "CREATE E BOX:\n";
                boxMat[lastNonzeroRow + 1][0] = levNow + 1;
                int a = x2U * pow(2, boxMat[row][0] - levNow);//boxMat[lastNonzeroRow + 1][1];
                gsInfo << x2U << " * pow(2, " << boxMat[row][0] << " - " << levNow << ") = " << a << "\n";
                boxMat[lastNonzeroRow + 1][1] = a;//x2U * pow(2, boxMat[row][0] - levNow);
                int b = y1U * pow(2, boxMat[row][0] - levNow);
                gsInfo << y1U << " * pow(2, " << boxMat[row][0] << " - " << levNow << ") = " << b << "\n";
                boxMat[lastNonzeroRow + 1][2] = b;//y1U * pow(2, boxMat[row][0] - levNow);
                int c = x2Bi;
                //boxMat[lastNonzeroRow + 1][3];
                gsInfo << "c: " << c << "\n";
                gsInfo << x2Bi  << " = " << c << "\n";
                boxMat[lastNonzeroRow + 1][3] = c;//x2Bi; //*pow(2, boxMat[row][0] - 3);
                gsInfo << "boxMat[lastNonzeroRow + 1][3]: " << boxMat[lastNonzeroRow + 1][3] << "\n";
                int d = y2U * pow(2, boxMat[row][0] - levNow);//boxMat[lastNonzeroRow + 1][4];
                boxMat[lastNonzeroRow + 1][4] = d;//y2U * pow(2, boxMat[row][0] - levNow);
                gsInfo << "a: " << a << "\t"  << b << "\t"  << c << "\t" << d << "\n";
                outfile << x2U << "\t" << y1U << "\t" << x2Bi << "\t" << y2U << "\t" << boxMat[row][0] << "\t" << levNow << "\n";
                                gsInfo<< "created the following box of level" << levNow + 1 << ": " <<
                                        boxMat[lastNonzeroRow + 1][1] << " " << boxMat[lastNonzeroRow + 1][2] << " " <<
                                        boxMat[lastNonzeroRow + 1][3] << " " << boxMat[lastNonzeroRow + 1][4] << "\n";
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



int main(int argc, char * argv[]){
    gsStopwatch clock;
    int row, acceptedsize;
    std::string xmlFile = "/home/turing/theydarov/xmlFiles/resultLocalMulti";
    std::string badFile = "/home/turing/theydarov/xmlFiles/badgeoLocalMulti";
    std::string pvdFile = "/home/turing/theydarov/pvdFiles/resultLocalMulti";
    int los = 0, nlos = 0, proj = 0;
    //DON'T FORGET ABOUT GLOBAL VARIABLES!
    int numCells = 4 * (int) pow(2, coarseLevel) * 4 * (int) pow(2, coarseLevel);
//    int numCells = 1 * (int) pow(2, coarseLevel) * 1 * (int) pow(2, coarseLevel);
    std::string givenGeo;
    int gradingExtent;
    //    real_t epsilon_g = 0.003, epsilon_f = 0.001;
    //    real_t epsilon_g = 0.004, epsilon_f = 0.004;
    //    real_t epsilon_g = 0.0114, epsilon_f = 0.0114195;
    //    real_t epsilon_g = 1e+6, epsilon_f = 0.003;
    real_t epsilon_g = 1e+6, epsilon_f = 5e-3;
    real_t lcx,lcy,ucx,ucy;
    //    real_t epsilon_g = 1e+6, epsilon_f = 1e+6;
    //    ofstream myfile;
    //    myfile.open ("/home/theydarov/october/example.txt",  fstream::app);
    std::string acCond = to_string(epsilon_g) + "and" + to_string(epsilon_f);
    int bads = 0;
    const int dim=2;

    //std::string inputInitial("/home/theydarov/gismo_unsupported/gismo/extensions/unsupported/filedata/surfaces/thb_small_jigsaw_puzzle_1.xml");
    //        std::string inputInitial("/home/theydarov/gismo_unsupported/gismo/extensions/motor/filedata/CoonsPatchleveltwo.xml");
//    std::string inputInitial("/home/turing/theydarov/geometries/indiana_thb.xml");
//    std::string inputInitial("/home/theydarov/september21/060921/geom.xml");
    //    std::string inputInitial("/home/theydarov/march21/190321/indianaDistortedMax.xml");
    //    std::string inputInitial("/home/theydarov/gismo_unsupported/gismo/extensions/motor/filedata/wigglysurface.xml");
    //    std::string inputInitial("/home/theydarov/february/100221/coonspatch2.xml");
    //    std::string inputInitial("/home/theydarov/february/120221/coonspatch3.xml");
    //    std::string inputInitial("/home/theydarov/march21/100321/coonspatch.xml");
    //        std::string inputInitial("/home/theydarov/gismo_unsupported/gismo/extensions/motor/filedata/BeziercurveCoonsPatch.xml");
    //    std::string inputInitial("/home/theydarov/gismo_unsupported/gismo/extensions/unsupported/filedata/planar/indiana_thb.xml");
//        std::string inputInitial("/home/turing/theydarov/geometries/WigglySquareIntDist.xml");
//    std::string inputInitial("/home/turing/theydarov/geometries/jigsawpuzzleTPEB.xml");
//    std::string inputInitial("/home/turing/theydarov/geometries/jigsawpuzzleTP.xml");
//    std::string inputInitial("/home/turing/theydarov/geometries/WigglyCoonsPatchIntDist.xml");
//    std::string inputInitial("/home/turing/theydarov/gismoMP/gismo/filedata/planar/square_and_puzzle.xml");
//    std::string inputInitial("/home/turing/theydarov/gismoMP/gismo/filedata/planar/two_squares_THB.xml");
//    std::string inputInitial("/home/turing/theydarov/gismoMP/gismo/filedata/planar/square_and_puzzle.xml");
//    std::string inputInitial("/home/turing/theydarov/gismoMP/gismo/filedata/domain2d/yeti_mp2THB.xml");
    std::string inputInitial("/home/turing/theydarov/geometries/footTHBLevel2.xml");
//    std::string inputInitial("/home/turing/theydarov/geometries/indiana_thb.xml");
//    std::string inputInitial("/home/turing/theydarov/auxiliary/resultLocalMulti.xml");
//    std::string inputInitial("/home/turing/theydarov/auxiliary/resultLocalMulti_Map_0.xml");
    //    std::string inputInitial("/home/theydarov/gismo_unsupported/gismo/extensions/motor/filedata/jku/thb_map.xml");
    //        givenGeo = "indianamultilevel";
    //        std::string inputInitial("/home/theydarov/gismo_unsupported/gismo/extensions/unsupported/filedata/planar/wireCutter_thb.xml");
    //    std::string inputInitial("/home/theydarov/gismo_unsupported/gismo/extensions/unsupported/filedata/planar/indiana_thb.xml");
    givenGeo = "jigsawpuzzleTPMulti";
//    givenGeo = "jigsawpuzzleTPEB";
    //    givenGeo = "wigglysurface";
//        givenGeo = "Indiana";
//    givenGeo+=  "WigglyCoonsIntDist";
    //    givenGeo = "Indiana";
//        givenGeo+=  "WigglySquareIntDist";
//        givenGeo+=  "WigglyCoonsPatch3lev";
    //    givenGeo+= "lexicog/**/raphic";
//    givenGeo+="upperhalf";
    //    givenGeo+="degree2";
    givenGeo += "L2";
    givenGeo += "LO";
    givenGeo += "NLO";
    //    givenGeo+="[0,1]cross[0,0dot25]";
    //    givenGeo = "thb_map";
    std::string fileLoc = "/home/turing/theydarov/txtFiles/resultLocalMulti";
    std::string input("/home/theydarov/june21/300621/" + givenGeo + acCond + ".xml" );
    outfile.open(fileLoc + ".txt");
    std::string inputPtsParams("/home/turing/theydarov/gismo/extensions/motor/filedata/jku/thb_parameters_and_points.xml");
    //    gsInfo << "I REQUEST " << numCells << " cells\n";
    outfile << "I REQUEST " << numCells << " cells\n";

    //Iterations
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
    gsFileData < > data0(inputInitial);
//    gsGeometry < > ::uPtr geom;
    gsGeometry < > ::uPtr geomhat; //accepted geometry
    gsGeometry < > ::uPtr geomtmp;
    gsGeometry<> * geomForCheck;
    gsGeometry<> *aPatch;
    gsMatrix < > coefs;
    gsMatrix < > coefshat;
    gsMultiPatch<>::uPtr mp = gsReadFile<>(inputInitial);
    gsMultiPatch<> newmp;
    gsInfo << "Give the grading extent:\n";
    for (int i = 0; i < numRed; i++) {
        gsInfo << "GIVE RED BOX\n";
        cin >> redBox[i][0] >> redBox[i][1] >> redBox[i][2] >> redBox[i][3] >> redBox[i][4];
    }
    cin >> gradingExtent;
    for ( size_t i = 0; i< mp->nPatches(); ++i )
    {
        gsInfo << "Patch number: " << i << "\n";
        geomForCheck = & mp->patch(i);
        gsTHBSplineBasis < 2, real_t > * THBFromGeo =
                dynamic_cast<gsTHBSplineBasis<2> *>(&mp->patch(i).basis().source());
        real_t boundaryerror = 0;


        //REDBOX IS TO BE RETHOUGHT
//        if (numRed == 0) {
//            redBox[0][0] = fineLevel;
//            redBox[0][1] = 0;
//            redBox[0][2] = 0;
//            redBox[0][3] = (THBFromGeo->numKnots(0,0) - 2*(THBFromGeo->degree(0) + 1) + 1)  * (int) pow(2, fineLevel);
//            redBox[0][4] = (THBFromGeo->numKnots(0,1) - 2*(THBFromGeo->degree(1) + 1) + 1)  * (int) pow(2, fineLevel);
//            gsInfo << "I ASSUME THAT REDBOX HAS THE COORDINATES:\n"
//                   <<
//                   redBox[0][0] << " " << redBox[0][1] << " " << redBox[0][2] << " " << redBox[0][3] << " " << redBox[0][4] << "\n";
//            outfile << "I ASSUME THAT REDBOX HAS THE COORDINATES:\n"
//
//                    <<
//                    redBox[0][0] << " " << redBox[0][1] << " " << redBox[0][2] << " " << redBox[0][3] << " " << redBox[0][4] << "\n";
//
//        }
        numRed = 1;
        int lastNonzeroRow = 0; //71;//92;
        int minusnumber, createdBoxNum = 0, currCellIndex, wasRebuilt, RTH = 0, centerInd = 0, needToEscape = 0;
        int x1U, x2U, y1U, y2U; // = rand()%32;
        int x1Bi, x2Bi, y1Bi, y2Bi;
        int createSpline, failed, ourBox[5], currArrayIndex;
        gsMatrix < index_t > lowCorners;
        gsMatrix < index_t > upCorners;
        gsVector < index_t > myLevel;
        currArrayIndex = 0;
        std::string output("results_optimized");
        coefs = mp->patch(i).coefs();
        coefshat = coefs;
        int inadequacyCount = 0;
        acceptedsize = THBFromGeo->size();
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
        real_t au = THBFromGeo->knot(0,0,0);// 0; // starting knot
        real_t bu = THBFromGeo->knot(THBFromGeo->maxLevel(),0,THBFromGeo->numKnots(THBFromGeo->maxLevel(),0)-1);//1; // ending knot
        real_t av = THBFromGeo->knot(0,0,0);//0; //-1;//0;                   // starting knot
        real_t bv = THBFromGeo->knot(THBFromGeo->maxLevel(),1,THBFromGeo->numKnots(THBFromGeo->maxLevel(),1) - 1);//1; //0.5;//1;                   // ending knot
        gsInfo << au << "\t" << av << "\t" << bu << "\t" << bv << "\n";
        unsigned interior = THBFromGeo->numKnots(0,1) - 2*(THBFromGeo->degree(1) + 1); // number of interior// knots for geometry by Ondine
        int degree = THBFromGeo->degree(0); //2;
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
        for(int levNow = coarseLevel; levNow >= 0; levNow--)
        {
            theLev = levNow;
            if (levNow < coarseLevel && wasRebuilt)
            {
                std::string input("/home/turing/theydarov/xmlFiles/resultLocalMulti.xml" );
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
            gsFunctionExpr<> f1("x", 2);
            gsFunctionExpr<> f2("if(y<1/(1+0.5), y, y+y^2 *sin(22*pi*x)/20)", 2);
            int extension = 2;
            gsVector<int> nonCheckedCells(4 * pow(2, levNow) * 4 * pow(2, levNow));
//            gsVector<int> nonCheckedCells(1 * pow(2, levNow) * 1 * pow(2, levNow));
            for (int i = 0; i < nonCheckedCells.size(); i++)
            {
                nonCheckedCells(i) = i;
            }
            for (int attempt = 0; attempt < 4 * pow(2, levNow) * 4 * pow(2, levNow); attempt++)
//            for (int attempt = 0; attempt < 1 * pow(2, levNow) * 1 * pow(2, levNow); attempt++)
            {
//                if(attempt > 5)     return 0;
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
                currCellIndex = pickCell(nonCheckedCells, currArrayIndex, levNow, x1U, y1U, x2U, y2U);
//                currCellIndex = pickCell(attempt, levNow, x1U, y1U, x2U, y2U);
                //WE
                gsVector<int> cellstack;
                outfile << "Number of nonchecked cells: " << nonCheckedCells.size() << "\n";
                //                gsInfo << "Initial basis:" << *THBFromGeo << std::endl;
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //Creates manually your mesh(put this stuff in a function)
                // ...a 2D-tensor-B-spline basis with this knot vector...
                // ...and a 2D-THB-spline basis out of the tensor-B-spline basis.
                //! [constBasis]
                gsFieldCreator<real_t> JakaField;
                //    for(int attempt = 0; attempt < 256; attempt++){
                gsInfo << "attempt " << attempt << ", CURRENT INDEX IS " << nonCheckedCells(currArrayIndex) <<
                       ",\t the coordinates of the box are\n" << x1U << " " << y1U << " " << x2U << " " << y2U << "\n";
                outfile << "attempt " << attempt << ", CURRENT INDEX IS " << nonCheckedCells(currArrayIndex) <<
                        ",\t the coordinates of the box are\n" << x1U << " " << y1U << " " << x2U << " " << y2U << "\n";
                //        gsInfo << "attempt " << attempt << "\n" << x1U << " " << y1U << " " << x2U << " " << y2U << "\n";
                //                gsInfo << "need is" << needToEscape << "\n";
                //                THB.getLevelUniqueSpanAtPoints
                //                needToEscape = 1;
                createSpline = 0;
                for (int currentrow = 0; currentrow < lastNonzeroRow + 1; currentrow++)
                {
                    row = currentrow;
                    RTH =
                            rebuildTheHierarchy(boxMat, row, x1U, x1Bi, x2U, x2Bi, y1U, y1Bi, y2U, y2Bi, levNow, lastNonzeroRow,
                                                createdBoxNum, centerInd, ourBox, needToEscape);
                    if (RTH) createSpline = 1;
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
                            break;
                        }
                    }
                    if(!isAdequate){
                        {
                            gsInfo << "The basis is inadequate, rewriting\n";
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
                        inadequacyCount++;
                        if(inadequacyCount <= 1) {
                            attempt--;
                            continue;
                        }
                        else{
                            inadequacyCount = 0;
                            break;
                        }
                    }
                    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                    outfile << THB << "\n";
                    //                                cout << "Before refinement:" << endl;
                    gsInfo << THB << "\n";
                    THB.anchors_into(anmat);
                    gsInfo << "REACHED HERE\n";
                    gsMatrix<real_t> basisInd(1, THB.size()); //1 if function is not presented
                    gsInfo << "REACHED HERE\n";
                    for (int l = 0; l < THB.size(); ++l)
                    {
                        basisInd(0, l) = 0.0;
                    }
                    gsInfo << "REACHED HERE\n";
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

                                basisInd(j) = 0.0;
                                //                            fineIndIsUsed(i) = 1;
                                break;
                            }
                        }
                    }
                    for (int m = 0; m < THB.size(); m++)
                    {
                        if(THBcoefs(m,0) == 0.0 || THBcoefs(m,1) == 0.0)
                        {
//                        gsInfo << m << "IS PROBLEMATIC\n";
                            basisInd(0,m) = 1.0;
                        }
                    }
                    gsInfo << "REACHED HERE\n";
                    for (int i = -gradingExtent; i < gradingExtent; ++i) {
                        for (int j = -gradingExtent; j < gradingExtent; ++j) {
                            updateCriterion(levNow, x1U + i, y1U + j, x1U + i + 1, y1U + j + 1, 1, basisInd, THB);
                        }
                    }

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
                    index_t totalNumberBasisFunct = THB.size();
                    gsLocalTHBFitting<dim,real_t>locTHBF(THB,mp->patch(i), basisInd);
                    locTHBF.doLocalFittingWithGeometry();
                    gsMatrix <real_t > punto(2, 1);
                    punto(0,0) = 1;punto(1,0) = 1;
                    int minusnumber = checkJacobianDeterminant(*locTHBF.result(), 10000, true, "points", 0);//
                    minusnumber = 0;
                    checkAccuracySufficiency(*locTHBF.result(), *geomForCheck, maxError, boundaryerror, au, bu, av, bv);
                    gsInfo << "the errors: " << maxError << " " << boundaryerror << "\n";
                    outfile << "the errors: " << maxError << " " << boundaryerror << "\n";
                    saveData(*locTHBF.result(), "ZURNA");
                    if (minusnumber == 0 && !(epsilon_g < maxError || epsilon_f < boundaryerror))
                    {
                        outfile << "Success!\n";
                        nonCheckedCells.removeElement(currArrayIndex); //REWRITE MAKING THIS AUTOMATIC DEPENDING ON THE METHOD
                        successfullAttempts++;
                        totalAttempts++;
                        gsInfo << "Success!\n";
                        saveData(*locTHBF.result(), output );
                        //fittedFunction = THBAccepted.makeGeometry(supps);
                        gsInfo << "SIZE: " << THB.size() << "\n";
                        outfile  << "SIZE: " << THB.size() << "\n";
                        acceptedsize = THB.size();
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
                        proj++;
//                    gsWrite(*fittedFunction, xmlFile);//MAYBE, THERE IS A WAY OF IMPROVING
                        gsWrite(*locTHBF.result(), xmlFile);//MAYBE, THERE IS A WAY OF IMPROVING
                        aPatch = locTHBF.result();

//                        gsInfo << "GEOMSIZE: " << aPatch->size();
                        //*fittedFunction = *locTHBF.result();
                        setSupports(THB, supps);
                        anmat2 = THB.anchors();
                        wasRebuilt = 1;
                        continue;

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
                            restoreTheHierarchy(createdBoxNum,
                                                lastNonzeroRow,
                                                boxMat,
                                                levNow,
                                                centerInd,
                                                ourBox,
                                                successfullAttempts);

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
                        int iterations = 1;
                        real_t fitting = 1;//0.99999;
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
//                    gsMatrix<> xy1 = fittedFunction->eval(uv1);
                        gsMatrix<> xy1 = locTHBF.result()->eval(uv1);
//                    gsQualityMeasure2<real_t> optimization(*fittedFunction, uv1, xy1, basisInd, true);
                        gsQualityMeasure2<real_t> optimization(*locTHBF.result(), uv1, xy1, basisInd, true);
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
//                        minusnumber = checkJacobianDeterminant(*fittedFunction, jacPts, true, output, it + 1);
                            minusnumber = checkJacobianDeterminant(*locTHBF.result(), jacPts, true, output, it + 1);
                            // Giving different names for files with (possibly) final results
//                        saveData(*fittedFunction, output + "_final", 0);
                            //savedata(*locTHBF.result(), output + "_final", 0);
                            //fittedFunction = THBAccepted.makeGeometry(supps);
                            if (minusnumber == 0) break;
                        }
                        gsInfo << "Number of negative points " << minusnumber << "\n";
                        outfile << minusnumber << "\n";
//                    checkAccuracySufficiency(*fittedFunction, *geomForCheck, maxError, boundaryerror, au, bu, av, bv);
                        checkAccuracySufficiency(*locTHBF.result(), *geomForCheck, maxError, boundaryerror, au, bu, av, bv);
                        if(epsilon_g < maxError || epsilon_f < boundaryerror)
                        {
                            gsInfo << "NO NEED\n";
                            gsInfo << "WITHDRAW\n";
                            totalAttempts++;
                            gsInfo << "the errors after LO: " << maxError << " " << boundaryerror << "\n";
                            outfile << "the errors after LO: " << maxError << " " << boundaryerror << "\n";
                            outfile << "WITHDRAW\n";
                            if(minusnumber!=0){
//                            gsWrite(*fittedFunction, "/home/turing/theydarov/xmlFiles/badgeoLocalMulti.xml");
                                gsWrite(*locTHBF.result(), "/home/turing/theydarov/xmlFiles/badgeoLocalMulti.xml");
                            }
                            restoreTheHierarchy(createdBoxNum,
                                                lastNonzeroRow,
                                                boxMat,
                                                levNow,
                                                centerInd,
                                                ourBox,
                                                successfullAttempts);
                            continue;

                        }
                        if (minusnumber == 0 && !(epsilon_g < maxError || epsilon_f < boundaryerror))
                        {

                            successfullAttempts++;
                            totalAttempts++;
//                        gsWrite(*fittedFunction,xmlFile);//MAYBE, THERE IS A WAY OF IMPROVING
                            gsWrite(*locTHBF.result(),xmlFile);//MAYBE, THERE IS A WAY OF IMPROVING
                            aPatch = locTHBF.result();
//                            gsInfo << "GEOMSIZE: " << aPatch->size();
//                        *fittedFunction = *locTHBF.result()->clone();
                            gsInfo << "LO WORKED!\n";
                            outfile << "LO WORKED!\n";
                            los++;
                            gsInfo << "SIZE: " << THB.size() << "\n";
                            outfile  << "SIZE: " << THB.size() << "\n";
                            acceptedsize = THB.size();
                            setSupports(THB, supps);
                            anmat2 = THB.anchors();
                            wasRebuilt = 1;
                            saveData(*locTHBF.result(), output );
                            continue;
//                    gsWrite(*fittedFunction, "/home/theydarov/june21/300621/badgeoLocalMulti.xml");
//                    return 0;
                        }
                        else
                        {
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
//                        gsMatrix<> xy1 = fittedFunction->eval(uv1);
                            gsMatrix<> xy1 = locTHBF.result()->eval(uv1);
//                        gsQualityMeasure2<real_t> optimization(*fittedFunction, uv1, xy1, basisInd, false);
                            gsQualityMeasure2<real_t> optimization(*locTHBF.result(), uv1, xy1, basisInd, false);
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
//                            saveData(*fittedFunction, output, it + 1);
                                //savedata(*locTHBF.result(), output, it + 1);
                                //fittedFunction = THBAccepted.makeGeometry(supps);
//                            minusnumber = checkJacobianDeterminant(*fittedFunction, jacPts, true, output, it + 1);
                                minusnumber = checkJacobianDeterminant(*locTHBF.result(), jacPts, true, output, it + 1);
                                // Giving different names for files with (possibly) final results
                                outfile << minusnumber << "\n";
//                            saveData(*fittedFunction, output + "_final", 0);
                                //savedata(*locTHBF.result(), output + "_final", 0);
                                //fittedFunction = THBAccepted.makeGeometry(supps);
                                if (minusnumber == 0)
                                {
                                    gsInfo << "the parametrization become regular\n";
                                    outfile << "the parametrization become regular\n";
//                                checkAccuracySufficiency(*fittedFunction,
//                                                         *geomForCheck,
//                                                         maxError,
//                                                         boundaryerror,
//                                                         au,
//                                                         bu,
//                                                         av,
//                                                         bv);
                                    checkAccuracySufficiency(*locTHBF.result(),
                                                             *geomForCheck,
                                                             maxError,
                                                             boundaryerror,
                                                             au,
                                                             bu,
                                                             av,
                                                             bv);
                                    if (epsilon_g > maxError && epsilon_f > boundaryerror)
                                    {

                                        gsInfo << "NLO WORKED!\n";
                                        outfile << "NLO WORKED!\n";
                                        successfullAttempts++;
                                        totalAttempts++;
                                        nlos++;
                                        setSupports(THB, supps);
                                        anmat2 = THB.anchors();
                                        gsInfo << "SIZE: " << THB.size() << "\n";
                                        outfile << "SIZE: " << THB.size() << "\n";
                                        acceptedsize = THB.size();
//                                    gsWrite(*fittedFunction,                                            xmlFile);//MAYBE, THERE IS A WAY OF IMPROVING
                                        gsWrite(*locTHBF.result(),                                            xmlFile);//MAYBE, THERE IS A WAY OF IMPROVING
                                        aPatch = locTHBF.result();
//                                        gsInfo << "GEOMSIZE: " << aPatch->size();
//                                    *fittedFunction = *locTHBF.result()->clone();
                                        wasRebuilt = 1;
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
//                                gsWrite(*fittedFunction, badFile);
                                    gsWrite(*locTHBF.result(), badFile);
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
        const gsGeometry<>::uPtr thePatch = gsReadFile<>(xmlFile + ".xml");
//        gsInfo <<  aPatch->size();
        newmp.addPatch(*thePatch);
    }
    gsFileData<> fd;
    fd<< newmp ;
    fd.dump("makeMultipatch_output1");
    gsWriteParaview(newmp, "two_puzzles",1000, true, false );
    gsInfo << "SIZE: " << acceptedsize << "\n";
    outfile << "SIZE: " << acceptedsize << "\n";
    gsInfo << "proj: " << proj  << ", linear: " << los << ", nonlinear: " << nlos << "\n";
    outfile << "proj: " << proj  << ", linear: " << los << ", nonlinear: " << nlos << "\n";
    std::string geomParaviewFile = givenGeo + acCond + "_Mesh_0";
//    gsWriteParaview(*fittedFunction, geomParaviewFile);
//    gsWriteParaview(*locTHBF->result(), geomParaviewFile);

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
//    saveData(*fittedFunction, pvdFile);
    gsInfo << gradingExtent;
    outfile << gradingExtent;
    return 0;
}