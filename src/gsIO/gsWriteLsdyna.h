/** @file gsWriteLsdyna2.h

    @brief  Provides interface between G+Smo and LS-DYNA.

    Should be moved in devel...

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): J. Speh
*/

#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <ios>


#include <gsCore/gsLinearAlgebra.h>

//#include <algorithm>

namespace gismo
{

template <typename T>
class gsWriteLsdyna
{
public:
    gsWriteLsdyna(const gsGeometry<T>& surf)
        :    
        mSurf(surf),
        mBasis(surf.basis()),
        mCrack(false),
        mStartPartId1(1001),
        mStartFormulationId1(1001),
        mStartPartId2(8001),
        mStartFormulationId2(98),
        mStartNodeId1(1000001),
        mStartNodeId2(8002),
        mStartElementId(1),
        mStartElementShellId(200001),
        pGeoEval(NULL),
        pDomIt(NULL),
        mFormId(0),
        mElementId(0),
        mNodeId2(0),
        mElementShellId(0),
        mInterpolNodes(),
        mNodesMap(),
        mThickBasis(NULL)
    {
        if (surf.coefDim() != 3)
        {
            GISMO_ERROR("Input must be a surface in 3D space.\n"
                        "To do: fix this error.\n");
        }
         
        if (surf.parDim() != 2)
        {
            GISMO_ERROR("Input must be a surface.\n"
                        "We are doing thin shell analysis.\n");
        }
            
    }

    // solid -- if we model a shell as thin solid 
    // degree -- degree of BSplines in thickness direction (only relevant if 
    //           solid is true)
    void produceLsdynaInputFile(std::ostream& out,
                                const std::string name = "out",
                                const bool solid = false,
                                const int degree = 2,
                                const bool crack = false)
    {

        if (solid)
        {
            produceLsdynaSolidInputFile(out, name, degree, crack);
        }
        else 
        { 
            produceLsdynaShellInputFile(out, name);
        }
    }


    void produceLsdynaShellInputFile(std::ostream& out, 
                                     const std::string& name)
    {
        out << std::setprecision(16);
        out << std::scientific;
        mCrack = false;

        // writing 
        
        writeHeader(out, name);
        
        writePartAndSection(out);
        
        writeNodes(out);
        
        // I suspect that we can leave this part out
        writeCurve(out);

        writeShellAndInterpolationElements(out);
        
        out << "*END\n";

    }


    void produceLsdynaSolidInputFile(std::ostream& out, 
                                     const std::string& name, 
                                     const int degree, 
                                     const bool crack)
    {
        out << std::setprecision(16);
        out << std::scientific;
        mCrack = crack;

        createThroughThicknessBasis(degree);
        gsMatrix<T> newCoefs;
        createNewCoefficients(0.5, newCoefs);

        writeHeader(out, name);

        writePartAndSection(out, true);
        
        writeNodes(out, newCoefs);

        // I suspect that we can leave this part out
        writeCurve(out);

        writeSolidAndInterpolationElements(out, newCoefs);
        
        out << "*END\n";
        
        delete mThickBasis;
    }


// ----------------------------------------------------------------------
// private member functions
// ----------------------------------------------------------------------

private:

    void writeSolidAndInterpolationElements(std::ostream& out, 
                                            const gsMatrix<T>& newCoefs)
    {
        pGeoEval = mSurf.evaluator(NEED_VALUE | NEED_MEASURE);
        pDomIt = mBasis.makeDomainIterator().release();

        gsVector<int> cwiseIntPoints(2);
        cwiseIntPoints[0] = mBasis.degree(0) + 1;
        cwiseIntPoints[1] = mBasis.degree(1) + 1;

        pDomIt->computeQuadratureRule( cwiseIntPoints );

        mFormId = mStartFormulationId1;
        mElementId = mStartElementId;
        mNodeId2 = mStartNodeId2;
        mElementShellId = mStartElementShellId;
        
        for(; pDomIt->good(); pDomIt->next())
        {
            // integration points through thickness
            gsDomainIterator<T>* pThickDomIt = 
                mThickBasis->makeDomainIterator().release();
            gsVector<int> thickIntPoints(1); 
            thickIntPoints[0] = mThickBasis->degree() + 1;
            pThickDomIt->computeQuadratureRule(thickIntPoints);
            pThickDomIt->evaluateBasis(1); // value and first partial derivs

            // normal integration
            pDomIt->computeActiveFunctions();
            pDomIt->evaluateBasis(1); 
            pGeoEval->evaluateAt(pDomIt->quNodes);
            
            defineElementGeneralizedSolid(out, newCoefs, pThickDomIt);

            elementGeneralizedSolid(out, pThickDomIt);
            
            int saveNodeId2 = mNodeId2;
            constrainedNodeInterpolation(out, newCoefs, pThickDomIt);
            
            elementInterpolationSolid(out, newCoefs, saveNodeId2, pThickDomIt);
            
            delete pThickDomIt;
        }

        delete pDomIt;
        delete pGeoEval;
        mNodesMap.clear();
        mInterpolNodes.clear();
    }

    void elementInterpolationSolid(std::ostream& out,
                                   const gsMatrix<T>& newCoefs,
                                   const int saveNodeId2,
                                   gsDomainIterator<T>* pThickDomIt)
    {
        const gsMatrix<unsigned>& actFunctions = pDomIt->computeActiveFunctions();
        const index_t numAct = actFunctions.rows();

        const index_t thNumAct = pThickDomIt->computeActiveFunctions().rows();
        
        gsMatrix<T> newNodes, params, thParams, values, thValues;
        getInterpolationNodes(newCoefs, pThickDomIt, params, thParams, 
                              values, thValues, newNodes); 

        const index_t numParams = params.cols();
        gsMatrix<T> nodeVals(numAct * thNumAct, numParams * thParams.cols());
        for (index_t act = 0; act != numAct; ++act)
        {
            for (index_t thAct = 0; thAct != thNumAct; thAct++)
            {
                for (index_t thCol = 0; thCol != thParams.cols(); thCol++)
                {
                    const index_t row = act * thNumAct + thAct;
                    nodeVals.block(row, thCol * numParams, 1, numParams) 
                        = values.row(act) * thValues(thAct, thCol);
                }
            }
        }

        const index_t active = numAct * thNumAct; // global active functions
        const index_t glNumIntPts = pDomIt->numQuNodes() * pThickDomIt->numQuNodes();
        gsMatrix<T> M, G;
        M.setZero(active, active);
        G.setZero(active, glNumIntPts);

        gsMatrix<T> glGrads; // global gradients

        // returns gradients evaluatade on quadrature nodes of the cube
        getGrads(glGrads, newCoefs, pThickDomIt); 
        
        for (index_t k = 0; k != pDomIt->numQuNodes(); ++k)
        {
            for (index_t thK = 0; thK != pThickDomIt->numQuNodes(); ++thK)
            {
                const T qWeight = pDomIt->quWeights[k] * pThickDomIt->quWeights[thK]
                    * glGrads.block(0, 3 * (3 * k + thK), 3, 3).determinant();
                
                const gsVector<T>& vals = pDomIt->basisValues().col(k);
                const gsVector<T>& thVals = pThickDomIt->basisValues().col(thK);
                gsVector<T> values(vals.rows() * thVals.rows());
                for (index_t r = 0; r != thVals.rows(); r++)
                {
                    values.block(r * vals.rows(), 0, vals.rows(), 1) 
                        = vals * thVals(r);
                }

                M.noalias() += qWeight * (values * values.transpose());
                G.col(k * pThickDomIt->numQuNodes() + thK) = qWeight * values;
            }
        }

        gsVector<T> diag = M.rowwise().sum();
        M.setZero();
        for (index_t col = 0; col != M.cols(); ++col)
        {
            M(col, col) = diag(col);
        }

        G = M.colPivHouseholderQr().solve(G);
        G = nodeVals.transpose() * G; // row = nodalPoint, col = intPoint


        gsVector<int> cwiseIntPoints(2);
        cwiseIntPoints[0] = mBasis.degree(0) + 1;
        cwiseIntPoints[1] = mBasis.degree(1) + 1;
        gsVector<int> cwiseIntEls = cwiseIntPoints - gsVector<int>::Ones(2);

        gsMatrix<T> nodeWeights(1, glNumIntPts);

        gsVector<int> cE;
        cE.setZero(2);
        int str = cwiseIntEls[0] + 1;
        int r = 0;

        do
        {
            const int id = cE[1] * str + cE[0];
            const index_t N = thParams.cols() - 1;
            
            for (index_t n = 0; n != N; n++)
            {
                nodeWeights.row(0) =
                    // 4 lower cube points
                    0.125 * (G.row(N * id + n) + 
                             G.row(N * (id + 1) + n) + 
                             G.row(N * (id + str) + n) + 
                             G.row(N * (id + str + 1) + n))
                    // 4 upper cube points
                    + 0.125 * (G.row(N * id + 1 + n) + 
                               G.row(N * (id + 1) + 1 + n) + 
                               G.row(N * (id + str) + 1 + n) + 
                               G.row(N * (id + str + 1) + 1 + n));

                out << "*ELEMENT_SOLID\n"
                    << std::setw(8) << mElementShellId 
                    << std::setw(8) << mStartPartId2 << "\n" 
                    << std::setw(8) << mNodesMap[saveNodeId2 + N * id + n] 
                    << std::setw(8) << mNodesMap[saveNodeId2 + N * (id + 1) + n] 
                    << std::setw(8) << mNodesMap[saveNodeId2 + N * (id + str + 1) + n] 
                    << std::setw(8) << mNodesMap[saveNodeId2 + N * (id + str) + n] 
                    << std::setw(8) << mNodesMap[saveNodeId2 + N * id + 1 + n] 
                    << std::setw(8) << mNodesMap[saveNodeId2 + N * (id + 1) + 1 + n] 
                    << std::setw(8) 
                    << mNodesMap[saveNodeId2 + N * (id + str + 1) + 1 + n] 
                    << std::setw(8) << mNodesMap[saveNodeId2 + N * (id + str) + 1 + n] 
                    << "\n"
                    << "*ELEMENT_INTERPOLATION_SOLID\n"
                    << std::setw(10) << mElementShellId
                    << std::setw(10) << mElementId 
                    << std::setw(10) << glNumIntPts << "\n";
                
                int counter = 1;
                for (index_t w = 0; w != glNumIntPts; ++w)
                {
                    out << std::setw(10) << w + 1 
                        << std::setw(10) << nodeWeights(0, w);

                    // we can have just 8 numbers in a row (4 * 2 = 8)
                    if (counter == 4)
                    {
                        out << "\n";
                        counter = 0;
                    }
                    counter++;
                }

                if (counter != 1)
                    out << "\n";
                
                mElementShellId++;
            }
            r++;
        } while (nextLexicographic(cE, cwiseIntEls));

        mElementId++;
        mFormId++;
    }

    void getGrads(gsMatrix<T>& glGrads,
                  const gsMatrix<T>& newCoefs,
                  gsDomainIterator<T>* pThickDomIt)
    {
        const gsMatrix<unsigned> actFunctions = pDomIt->computeActiveFunctions();
        const index_t numAct = actFunctions.rows();
        const index_t thNumAct = pThickDomIt->computeActiveFunctions().rows();


        glGrads.setZero(3, pDomIt->numQuNodes() * pThickDomIt->numQuNodes() * 3); 
        glGrads.setZero();

        gsMatrix<T> grad , thGrad, value, thValue;

        for (index_t par = 0; par != pDomIt->numQuNodes(); par++)
        {
            for (index_t act = 0; act != numAct; act++)
            {
                unsigned index = actFunctions(act, 0);
                mBasis.evalSingle_into(index, pDomIt->quNodes.col(par), value);
                mBasis.derivSingle_into(index, pDomIt->quNodes.col(par), grad);

                for (index_t thAct = 0; thAct != thNumAct; thAct++)
                {
                    for (index_t thPar = 0; 
                         thPar != pThickDomIt->numQuNodes(); 
                         thPar++)
                    {
                        const gsMatrix<T>& p = pThickDomIt->quNodes.col(thPar);
                        mThickBasis->evalSingle_into(thAct, p, thValue);
                        mThickBasis->derivSingle_into(thAct, p, thGrad);
                        
                        const index_t parameter = 
                            par * pThickDomIt->numQuNodes() + thPar;
                        const index_t coefInd = act * thNumAct + thAct;
                        const gsVector<T>& coefficient = 
                            newCoefs.row(coefInd).transpose();
                        
                        glGrads.col(3 * parameter) +=
                            grad(0, 0) * thValue(0, 0) * coefficient;

                        glGrads.col(3 * parameter + 1) +=
                            grad(1, 0) * thValue(0, 0) * coefficient;
                        
                        glGrads.col(3 * parameter + 2) +=
                            value(0, 0) * thGrad(0, 0) * coefficient;
                    }
                }
            }
        }

        
    }

    void constrainedNodeInterpolation(std::ostream& out, 
                                      const gsMatrix<T>& newCoefs,
                                      gsDomainIterator<T>* pThickDomIt)
    {
        gsMatrix<unsigned> actFunctions = pDomIt->computeActiveFunctions();
        gsMatrix<T> newNodes, params, thParams, values, thValues;
        getInterpolationNodes(newCoefs, pThickDomIt, params, thParams, 
                              values, thValues, newNodes); 
        
        const index_t thSize = mThickBasis->size();
        for (index_t col = 0; col != params.cols(); col++)
        {
            // first we must check if this node already exists
            std::pair<T, T> p(params(0, col), params(1, col));
            typename std::map<std::pair<T, T>, int>::iterator it = 
                mInterpolNodes.find(p);

            out << std::setprecision(9);
            
            if (it != mInterpolNodes.end()) // node found
            {
                out << "$\n";
                for (index_t thCol = 0; thCol != thParams.cols(); thCol++)
                {
                    index_t index = thParams.cols() * col + thCol;
                    out << "$ node " << mNodeId2 <<  " duplicates node " 
                        << mInterpolNodes[p] + thCol
                        << "$ coordinate: " 
                        << newNodes.col(index).transpose() << "\n"
                        << "$\n";

                    mNodesMap[mNodeId2] = mInterpolNodes[p] + thCol;
                    mNodeId2++;
                }
                out << "$\n";
                continue;
            }
            
            // else 
            // we have a new node
            
            for (index_t thCol = 0; thCol != thParams.cols(); ++thCol)
            {
                mNodesMap[mNodeId2 + thCol] = mNodeId2 + thCol;
            }
            mInterpolNodes[p] = mNodeId2;
            
            for (index_t thCol = 0; thCol != thParams.cols(); ++thCol)
            {
                out << std::setprecision(9);
                out << "*NODE\n"
                    << std::setw(8) << mNodeId2;
                
                for (index_t  row = 0; row < newNodes.rows(); ++row)
                {
                    out << std::setw(16) 
                        << newNodes(row, thParams.cols() * col + thCol);
                }
                
                out << std::fixed << std::setprecision(0) << std::setw(7) 
                    << 7.0 << "." << std::setw(7) << 7.0 << ".\n";

                out << std::scientific << std::setprecision(9);

                unsigned nmbActWeights = 0; // number of active functions at weights
                gsMatrix<T> weights(actFunctions.rows() * thSize, 1);

                for (index_t act = 0; act != actFunctions.rows(); ++act)
                {
                    for (index_t thAct = 0; thAct != thSize; thAct++)
                    {
                        T w = values(act, col) * thValues(thAct, thCol);
                        weights(thSize * act + thAct, 0) = w;
                        if (w != 0.0)
                        {
                            nmbActWeights++;
                        }
                    }
                }

                out << "*CONSTRAINED_NODE_INTERPOLATION\n"
                    << std::setw(10) << mNodeId2 
                    << std::setw(10) << nmbActWeights << "\n";
                
                int counter = 1;
                out << std::setprecision(3);

                for (index_t act = 0; act != actFunctions.rows(); act++)
                {
                    for (index_t thAct = 0; thAct != thSize; thAct++)
                    {
                        const index_t index = thSize * act + thAct;
                    
                        if (weights(index, 0) != 0.0)
                        {
                            const int t = mStartNodeId1 + 
                                thSize * actFunctions(col, 0) + thAct;
                            
                            out << std::setw(10) << t 
                                << std::setw(10) << weights(index, 0);
                        
                            // we can have just 8 numbers in a row (4 * 2 = 8)
                            if (counter == 4)
                            {
                                out << "\n";
                                counter = 0;
                            }
                            counter++;
                        }
                    }
                }
                
                if (counter != 1)
                    out << "\n";
                
                mNodeId2++;
            }
        }
    }
    

    void getInterpolationNodes(const gsMatrix<T>& newCoefs,
                               gsDomainIterator<T>* pThickDomIt,
                               gsMatrix<T>& params,
                               gsMatrix<T>& thParams,
                               gsMatrix<T>& values, // numAct x params
                               gsMatrix<T>& thValues, // thNumAct x thParams
                               gsMatrix<T>& newNodes)
    {
        // ======================================================================
        // surface evaluation
        // ======================================================================
        
        gsVector<unsigned> uCwiseIntPoints(2);
        uCwiseIntPoints(0) = static_cast<unsigned>(mBasis.degree(0) + 1);
        uCwiseIntPoints(1) = static_cast<unsigned>(mBasis.degree(1) + 1);
        
        gsVector<T> low = pDomIt->lowerCorner();
        gsVector<T> upp = pDomIt->upperCorner();
        params = gsPointGrid<T>(low, upp, uCwiseIntPoints);
        
        gsMatrix<unsigned> actFunctions = pDomIt->computeActiveFunctions();
        
        values.resize(actFunctions.rows(), params.cols());
        values.setZero();
        gsMatrix<T> tmp;

        for (index_t act = 0; act != actFunctions.rows(); act++)
        {
            for (index_t col = 0; col != params.cols(); col++)
            {
                mBasis.evalSingle_into(actFunctions(act, 0), params.col(col), tmp);
                values(act, col) = tmp(0, 0);
            }
        }

        // ======================================================================
        // thickness 
        // ======================================================================

        gsVector<T> thLow = pThickDomIt->lowerCorner();
        gsVector<T> thUpp = pThickDomIt->upperCorner();
        
        gsVector<unsigned> thIntPoints(1);
        thIntPoints(0) = mThickBasis->degree() + 1;
        
        thParams = gsPointGrid<T>(thLow, thUpp, thIntPoints);
        thValues.resize(mThickBasis->size(), thParams.cols());
        thValues.setZero();

        for (index_t act = 0; act != mThickBasis->size(); act++)
        {
            for (index_t col = 0; col != thParams.cols(); col++)
            {
                mThickBasis->evalSingle_into(act, thParams.col(col), tmp);
                thValues(act, col) = tmp(0, 0);
            }
        }


        // ======================================================================
        // thick shell nodes
        // ======================================================================

        const index_t numNewNodes = params.cols() * thParams.cols();
        newNodes.resize(3, numNewNodes);
        newNodes.setZero();
        gsVector<T> val(3);
        gsVector<T> tmp1;
    

        for (index_t col = 0; col != params.cols(); col++)
        {
            for (index_t thCol = 0; thCol != thParams.cols(); thCol++)
            {
                val.setZero();
                for (index_t act = 0; act != actFunctions.rows(); act++)
                {
                    for (index_t thAct = 0; thAct != mThickBasis->size(); thAct++)
                    {
                        const unsigned index = actFunctions(act, 0);
                        const index_t coefIndex= mThickBasis->size() * index + thAct;
                        val += values(act, col) * thValues(thAct, thCol) * 
                            newCoefs.row(coefIndex).transpose();
                    }
                }
                newNodes.col(thParams.cols() * col + thCol) = val;
            }
        }

    }

    void elementGeneralizedSolid(std::ostream& out,
                                 gsDomainIterator<T>* pThickDomIt)
    {
        // change this function if you would like to have more 
        // elements over a thickness

        gsMatrix<unsigned> actFunctions = pDomIt->computeActiveFunctions();

        const index_t numAct = actFunctions.rows();
        const index_t numThickAct = pThickDomIt->computeActiveFunctions().rows();
        
        out << "*ELEMENT_GENERALIZED_SOLID\n"
            << std::setw(10) << mElementId 
            << std::setw(10) << mFormId 
            << std::setw(10) << numAct * numThickAct << "\n";

        int counter = 1;
        for (int index = 0; index != numAct; index++)
        {
            for (int thickI = 0; thickI != numThickAct; thickI++)
            {
                const int nodeId = mStartNodeId1 + 
                    numThickAct * actFunctions(index) + thickI;

                out << std::setw(10) << nodeId;

                if (counter == 8) // we can write only 8 ids on one line
                {
                    counter = 0;
                    out << "\n";
                }
                counter++;
            }
        }
        
        if (counter != 1)
            out << "\n";
    }

    void defineElementGeneralizedSolid(std::ostream& out,
                                       const gsMatrix<T>& newCoefs,
                                       gsDomainIterator<T>* pThickDomIt)
    {
        const index_t numAct = pDomIt->computeActiveFunctions().rows();
        const index_t numThickAct = pThickDomIt->computeActiveFunctions().rows();

        out << "*DEFINE_ELEMENT_GENERALIZED_SOLID\n"
            << std::setw(10) << mFormId 
            << std::setw(10) << pDomIt->numQuNodes() * pThickDomIt->numQuNodes()
            << std::setw(10) << numAct  * numThickAct 
            << std::setw(10) << 0 << "\n";
        
        for (index_t k = 0; k < pDomIt->numQuNodes(); ++k)
        {
            gsMatrix<T> derivs = pDomIt->basisDerivs(1).col(k);
            derivs.resize(2, numAct);
            
            const gsVector<T>& values = pDomIt->basisValues().col(k);

            for (index_t thickK = 0; thickK < pThickDomIt->numQuNodes(); ++thickK)
            {
                gsVector<T> thickDerivs = pThickDomIt->basisDerivs(1).col(thickK);
                const gsVector<T>& thickValues = 
                    pThickDomIt->basisValues().col(thickK);
                
                index_t gaussPointNumber = 
                    pThickDomIt->numQuNodes() * k + thickK + 1;

                out << "$ gauss point " << gaussPointNumber << "\n";            
                
                out << std::setprecision(13);
                out << std::setw(20) 
                    << pDomIt->quWeights[k] * pThickDomIt->quWeights[thickK]
                    << "\n";
                
                for (index_t i = 0; i != numAct; ++i)
                {
                    for (index_t thickI = 0; thickI != numThickAct; ++thickI)
                    {
                        out << std::setw(20) << values[i] * thickValues[thickI]
                            << std::setw(20) << derivs(0, i) * thickValues[thickI]
                            << std::setw(20) << derivs(1, i) * thickValues[thickI]
                            << std::setw(20) << values[i] * thickDerivs[thickI]
                            << "\n";
                    }
                }
            }
        }
    }

    void writeNodes(std::ostream& out,
                    const gsMatrix<T>& newCoefs)
    {
        gsMatrix<unsigned>* pBoundary = mBasis.boundary();
        
        int nodeId = mStartNodeId1;

        // width for writing id
        const int idw = 8;
        // width for writing coefficients
        const int coefsw = 16;

        int counter = 0; // counter for boundary values

        out << std::setprecision(9);
        out << std::setiosflags(std::ios::uppercase);

        out << "*NODES\n";

        // for loop over coefficients
        const unsigned nmbOfSurfaceCoefs = static_cast<unsigned>(mBasis.size());
        const int nmbOfThickCoefs = mThickBasis->size();
        for (unsigned i = 0; i < nmbOfSurfaceCoefs; ++i)
        {
            // here I heavily depend on fact that boundary matrix is
            // ordered
            T translConst = 0.0;
            if (counter < (*pBoundary).rows() && i == (*pBoundary)(counter))
            {
                counter++;
                if (!mCrack)
                {
                    translConst = 7.0;
                }
            }

            for (int j = 0; j != nmbOfThickCoefs; ++j)
            {
                out << std::setw(idw) << nodeId;
                nodeId++;

                // for loop over x, y, z
                for (int dim = 0; dim < newCoefs.cols(); ++dim)
                {
                    out << std::setw(coefsw) 
                        << newCoefs(i * nmbOfThickCoefs + j, dim);
                }

                out << std::fixed << std::setprecision(0) << std::setw(idw - 1) 
                    << translConst << "." << std::setw(idw - 1) << 0. << ".\n";
            
                out << std::scientific << std::setprecision(9);
            }
        }            
        delete pBoundary;
    }


    void createNewCoefficients(const T h, // thickness
                               gsMatrix<T>& newCoefs)
    {
        // compute the normals at greville points
        // new control points will be the same as old one, just moved in 
        // the direction of normals
        gsMatrix<> grevile = mBasis.anchors();
    
        gsMatrix<> deriv;
        mSurf.deriv_into(grevile, deriv);
        
        const index_t n = grevile.cols();
        gsMatrix<> normals(3, n);
        normals.setZero();

        for (index_t i = 0; i != n; i++)
        {
            const gsVector<real_t, 3>& first = deriv.col(2 * i);
            const gsVector<real_t, 3>& second = deriv.col(2 * i + 1);
            const gsVector<> vec = first.cross(second);
            normals.col(i) = vec / vec.norm();
        }
        
        // compute the scalars for points
        const int numPtsThick = mThickBasis->size(); // nmb points through thickness
        gsVector<T> scalar(numPtsThick);
        const T diff = 2. / (numPtsThick - 1);
        for (index_t i = 0; i != numPtsThick; ++i)
        {
            scalar(i) = -1 + i * diff;
        }
        
        const gsMatrix<T> coefs = mSurf.coefs();
        newCoefs.setZero(coefs.rows() * numPtsThick, coefs.cols());



        const double halfH = h / 2;
        for (index_t i = 0; i != coefs.rows(); ++i)
        {
            const gsVector<T>& node = coefs.row(i);
            const gsVector<T>& normal = normals.col(i).transpose();
            
            for (int j = 0; j != numPtsThick; j++)
            {
                newCoefs.row(numPtsThick * i + j) = 
                    node + (scalar(j)  * halfH * normal);
            }
        }
    }
    
    void createThroughThicknessBasis(int degree)
    {
        gsKnotVector<T> kv(0, 1, 0, degree + 1); 
        mThickBasis = new gsBSplineBasis<T>(kv);

    }

    void writeShellAndInterpolationElements(std::ostream& out)
    {
        pGeoEval = mSurf.evaluator(NEED_JACOBIAN | NEED_VALUE | NEED_MEASURE 
                                   | NEED_GRAD_TRANSFORM);
        pDomIt = mBasis.makeDomainIterator().release();
        
        gsVector<int> cwiseIntPoints(2);
        cwiseIntPoints[0] = mBasis.degree(0) + 1;
        cwiseIntPoints[1] = mBasis.degree(1) + 1;

        pDomIt->computeQuadratureRule( cwiseIntPoints );

        mFormId = mStartFormulationId1;
        mElementId = mStartElementId;
        mNodeId2 = mStartNodeId2;
        mElementShellId = mStartElementShellId;
        
        for(; pDomIt->good(); pDomIt->next())
        {
            pDomIt->computeActiveFunctions();
            pDomIt->evaluateBasis(1); // value and first partial derivatives
            
            pGeoEval->evaluateAt(pDomIt->quNodes);
            
            defineElementGeneralizedShell(out);
            
            elementGeneralizedShell(out);
            
            int saveNodeId2 = mNodeId2;
            constrainedNodeInterpolation(out);
            
            elementInterpolationShell(out, saveNodeId2);
        }
        
        delete pDomIt;
        delete pGeoEval;
        mNodesMap.clear();
        mInterpolNodes.clear();
    }        

    void elementInterpolationShell(std::ostream& out, 
                                   const int saveNodeId2)
    {
        const gsMatrix<unsigned>& actFunctions = pDomIt->computeActiveFunctions();
        const index_t numAct = actFunctions.rows();
        
        gsVector<unsigned> cwiseIntPoints(2);
        cwiseIntPoints[0] = static_cast<unsigned>(mBasis.degree(0) + 1);
        cwiseIntPoints[1] = static_cast<unsigned>(mBasis.degree(1) + 1);

        gsVector<T> low = pDomIt->lowerCorner();
        gsVector<T> upp = pDomIt->upperCorner();
        gsMatrix<T> params = gsPointGrid<T>(low, upp, cwiseIntPoints);

        gsMatrix<T> nodeVals(numAct, params.cols());

        gsMatrix<T> tmp;        
        for (index_t row = 0; row != numAct; row++)
        {
            mBasis.evalSingle_into(actFunctions(row, 0), params, tmp);
            nodeVals.row(row) = tmp.row(0);
        }
        
        gsMatrix<T> M, G(numAct, pDomIt->numQuNodes());
        M.setZero(numAct, numAct);
        
        for (index_t k = 0; k < pDomIt->numQuNodes(); ++k) // loop over quadrature nodes
        {
            // compute the "jacobi" at the gauss point: jacobi of the surface geometry,
            // augmented by the scaled (by h/2) normal vector as the last column

            // const T qweight1 = pDomIt->quWeights[k] * pGeoEval->measure(k);
            
            
            const T rthick = 0.05;
/*            const gsMatrix<T>& coefs = mSurf.coefs();
            gsMatrix<T> actCoefs(coefs.cols(),numAct);

            for (index_t j = 0; j < numAct; ++j)
            {
                actCoefs.col(j) = coefs.row(actFunctions(j,0)).transpose();
            }

            gsMatrix<T> derivs = pDomIt->basisDerivs(1).col(k);
            derivs.resize(2,numAct);

            const gsMatrix<T> Jg1 = actCoefs * derivs.transpose(); */
            // upper code is equaivalent to pGeoEval->jacobian(k)
            const gsMatrix<T>& Jg1 = pGeoEval->jacobian(k);
            // std::cout << "Jg2: \n " << Jg2 << "\n"
            //           << "Jg1: \n " << Jg1 << "\n\n" << std::endl;
            const gsVector<real_t,3>& t1 = Jg1.col(0);
            const gsVector<real_t,3>& t2 = Jg1.col(1);
            const gsVector<T> vec = t1.cross(t2);

            gsMatrix<T> Jg(Jg1.rows(),Jg1.cols() + 1);
            Jg.block(0,0,Jg1.rows(),Jg1.cols()) = Jg1;
            Jg.col(Jg1.cols()) = (rthick / 2) * vec / vec.norm();

            const T qweight = pDomIt->quWeights[k] * Jg.determinant();
            // std::cout << "qweight1: " << qweight1 << "\n"
            //           << "qweight: " << qweight << std::endl;
            

            M.noalias() += qweight * (pDomIt->basisValues().col(k) *
                                       pDomIt->basisValues().col(k).transpose());

            G.col(k) = qweight * pDomIt->basisValues().col(k) ;
        }

        // solve it as diagonal matrix
        gsVector<T> diag = M.rowwise().sum();

        M.setZero();
        for (int col = 0; col != M.cols(); ++col)
        {
            M(col, col) = diag(col);
        }

        G = M.colPivHouseholderQr().solve(G); 
        G = nodeVals.transpose() * G;

        // loop over all squares inside this element
        gsVector<int> cwiseIntEls(2);
        cwiseIntEls(0) = mBasis.degree(0);
        cwiseIntEls(1) = mBasis.degree(1);
        
        gsMatrix<T> nodeWeights(1, pDomIt->numQuNodes());
        
        gsVector<int> cE;
        cE.setZero(2);

        int str = cwiseIntEls[0] + 1; // stride
        int r = 0;

        do
        {
            int id = cE[1] * str + cE[0];

            nodeWeights.row(0) =
                0.25 * (G.row(id)         +
                        G.row(id + 1)     +
                        G.row(id + str)   +
                        G.row(id + str + 1));
            
            out << "*ELEMENT_SHELL\n"
                 << std::setw(8) << mElementShellId 
                 << std::setw(8) << mStartPartId2 <<
                std::setw(8) << mNodesMap[saveNodeId2 + id] <<
                std::setw(8) << mNodesMap[saveNodeId2 + id + 1] <<
                std::setw(8) << mNodesMap[saveNodeId2 + id + str + 1] <<
                std::setw(8) << mNodesMap[saveNodeId2 + id + str] << "\n"
                 << "*ELEMENT_INTERPOLATION_SHELL\n"
                 << std::setw(10) << mElementShellId << 
                std::setw(10) << mElementId <<
                std::setw(10) << pDomIt->numQuNodes() << "\n";
            
            int counter = 1;
            for (index_t w = 0; w < pDomIt->numQuNodes(); ++w)
            {
                out << std::setw(10) << w + 1 <<
                    std::setw(10) << nodeWeights(0, w);

                // we can have just 8 numbers in a row (4 * 2 = 8)
                if (counter == 4)
                {
                    out << "\n";
                    counter = 0;
                }
                counter++;
            }
            if (counter != 1)
            {
                out << "\n";
            }
            
            mElementShellId++;
            r++;
        } while (nextLexicographic(cE, cwiseIntEls));

        mElementId++;
        mFormId++;
    }


            
    void constrainedNodeInterpolation(std::ostream& out)
    {
        gsVector<unsigned> cwiseIntPoints(2);
        cwiseIntPoints[0] = static_cast<unsigned>(mBasis.degree(0) + 1);
        cwiseIntPoints[1] = static_cast<unsigned>(mBasis.degree(1) + 1);

        gsVector<T> low = pDomIt->lowerCorner();
        gsVector<T> upp = pDomIt->upperCorner();
        gsMatrix<T> params = gsPointGrid<T>(low, upp, cwiseIntPoints);
        
        gsMatrix<T> newNodes;
        mSurf.eval_into(params, newNodes);
        
        const gsMatrix<unsigned>& actFunctions = pDomIt->computeActiveFunctions();
        
        for (index_t col = 0; col != newNodes.cols(); ++col)
        {
            std::pair<T, T> p(params(0, col), params(1, col));
            typename std::map<std::pair<T, T>, int>::iterator it =
                mInterpolNodes.find(p);

            out << std::setprecision(9);
            
            if (it != mInterpolNodes.end()) // node found
            {
                out << "$\n"
                     << "$ node " << mNodeId2 << " duplicates node " <<
                    mInterpolNodes[p] << "\n"
                     << "$ coordinates: " << newNodes.col(col).transpose() << "\n"
                     << "$\n";

                mNodesMap[mNodeId2] = mInterpolNodes[p];
                mNodeId2++;
                continue;
            }
            
            // we have a new node
            
            mNodesMap[mNodeId2] = mNodeId2;
            mInterpolNodes[p] = mNodeId2;
            
            // --------------------------------------------------
            // printing information about the node
            // --------------------------------------------------

            out << "*NODE\n"
                 << std::setw(8) << mNodeId2;
            for (int row = 0; row < newNodes.rows(); ++row)
            {
                out << std::setw(16) << newNodes(row, col);
            }

            out << std::fixed << std::setprecision(0) << std::setw(7) 
                 << 7.0 << "." << std::setw(7) << 7.0 << ".\n";

            out << std::scientific << std::setprecision(9);
            
            // --------------------------------------------------
            // printing information about the interpolation
            // --------------------------------------------------
            
            unsigned nmbActWeights = 0; // number of active functions at weights
            gsVector<T> weights(actFunctions.rows());
            for (index_t row = 0; row != actFunctions.rows(); ++row)
            {
                gsMatrix<T> w;
                mBasis.evalSingle_into(actFunctions(row, 0), params.col(col), w);
                weights(row) = w(0, 0);
                if (w(0, 0) != 0.0)
                {
                    nmbActWeights++;
                }
            }
            
            out << "*CONSTRAINED_NODE_INTERPOLATION\n"
                 << std::setw(10) << mNodeId2 << std::setw(10) 
                 << nmbActWeights << "\n";
            
            int counter = 1;
            out << std::setprecision(3);

            for (index_t index = 0; index != weights.rows(); index++)
            {
                if (weights(index) != 0.0)
                {
                    out << std::setw(10) << actFunctions(index, 0) + mStartNodeId1
                         << std::setw(10) << weights(index);

                    // we can have just 8 numbers in a row (4 * 2 = 8)
                    if (counter == 4)
                    {
                        out << "\n";
                        counter = 0;
                    }
                    counter++;
                }
            }

            if (counter != 1)
            {
                out << "\n";
            }

            mNodeId2++;
        }
    }

  
    void elementGeneralizedShell(std::ostream& out)
    {

        const gsMatrix<unsigned>& actFunctions = pDomIt->computeActiveFunctions();
        const index_t numAct = actFunctions.rows();

        out << "*ELEMENT_GENERALIZED_SHELL\n"
             << std::setw(10) << mElementId << std::setw(10) << mFormId <<
            std::setw(10) << numAct << "\n";
        
        int counter = 1;
        for (index_t index = 0; index != actFunctions.size(); index++)
        {
            out << std::setw(10) << mStartNodeId1 + actFunctions(index);

            if (counter == 8) // we can write only 8 ids on one line
            {
                counter = 0;
                out << "\n";
            }
            counter++;
        }
        if (counter != 1)
            out << "\n";
    }

    
    void defineElementGeneralizedShell(std::ostream& out)
    {
        const index_t numAct = pDomIt->computeActiveFunctions().rows();
        
        out << "*DEFINE_ELEMENT_GENERALIZED_SHELL\n"
             << std::setw(10) << mFormId << std::setw(10) << pDomIt->numQuNodes()
             << std::setw(10) << numAct << std::setw(10) << 0 << ",&IFORM\n";

        // ---------------------------------------------------------------------
        // Quadrature weight and values
        // ---------------------------------------------------------------------

        for (index_t k = 0; k != pDomIt->numQuNodes(); ++k)
        {
            gsMatrix<T> derivs = pDomIt->basisDerivs(1).col(k);
            derivs.resize(2, numAct);
            
            const gsVector<T>& values = pDomIt->basisValues().col(k);

            out << "$ gauss point " << k + 1 << "\n";
            out << std::setprecision(13);
            out << std::setw(20) << pDomIt->quWeights[k]
                 << "\n";

            for (index_t i = 0; i != numAct; ++i)
            {
                out << std::setw(20) << values(i);
                for (index_t u = 0; u != 2; ++u) // parametric directions
                {
                    out << std::setw(20) << derivs(u, i);
                }
                out << "\n";
            }
            
            checkValuesAndDerivatives(k, values, derivs);
        }
        
        // ---------------------------------------------------------------------
        // Nodal differential quadrature values
        // ---------------------------------------------------------------------


        writeNodalDifferentialValues2(out);
    }
                                                          
    
    void writeNodalDifferentialValues(std::ostream& out)
    {
        const index_t numAct = pDomIt->computeActiveFunctions().rows();
        
        const gsMatrix<T>& A = pDomIt->basisValues().transpose();
        // A -- rows = quadrature points
        //   -- cols = active functions
        
        const gsMatrix<T>& B = pDomIt->basisDerivs(1).transpose();
        // B -- rows = quadrature points
        //   -- cols = 2 * active functions

        gsMatrix<T> nodalValues = A.fullPivHouseholderQr().solve(B);
        
        bool ok = true;
        for (index_t k = 0; k != numAct; ++k)
        {
            bool zeroRow = true;

            gsMatrix<T> derivs = nodalValues.row(k);
            derivs.resize(2, numAct);

            for (index_t row = 0; row != derivs.rows(); row++)
            {
                for (index_t col = 0; col != derivs.cols(); col++)
                {
                    if (1e-9 < std::abs(derivs(row, col)))
                    {
                        zeroRow = false;
                    }
                }
                if (zeroRow)
                {
                    ok = false;
                    break;
                }
            }
        }
            
        if (!ok)
        {
            const index_t cols = A.cols();
            const index_t rows = A.rows();
            const index_t size = rows + cols;
            gsMatrix<double> mat(size, size);
            mat.block(   0, 0   , rows, cols) = A;
            mat.block(rows, 0   , cols, cols) = 
                2 * Eigen::MatrixXd::Identity(cols, cols);

            mat.block(rows, cols, cols, rows) = A.transpose();
            mat.block(   0, cols, rows, rows).setZero();

            gsMatrix<double> rhs(rows + cols, B.cols());
            rhs.setZero();
            rhs.block(0, 0, rows, B.cols()) = B;

            nodalValues = mat.fullPivHouseholderQr().solve(rhs);
        }

        for (index_t k = 0; k < numAct; ++k)
        {
            gsMatrix<T> derivs = nodalValues.row(k);
            derivs.resize(2, numAct);

            if (!ok)
            {
                bool zeroRow = true;
                for (index_t row = 0; row != derivs.rows(); row++)
                {
                    for (index_t col = 0; col != derivs.cols(); col++)
                    {
                        if (1e-9 < std::abs(derivs(row, col)))
                        {
                            zeroRow = false;
                        }
                    }

                    if (zeroRow)
                    {
                        GISMO_ERROR("Zero row");
                    }
                }
            }


            for (index_t row = 0; row != derivs.rows(); ++row)
            {
                if (1e-9 < derivs.row(row).sum())
                {
                    gsWarn << "$ ** Warning \n"
                           << "$    at gauss point: " << k + 1 << "\n"
                           << "$    sum of derivatives of basis functions is not "
                        "equal to 0 in direction " << row << "\n"
                           << "$   sum is: " << derivs.row(row).sum() << "\n";
                }
            }

            out << "$ node: " << k + 1 <<"\n";
            for (index_t col = 0; col != derivs.cols(); ++col)
            {
                for (index_t row = 0; row != derivs.rows(); ++row)
                {
                    out << std::setw(20) << derivs(row, col);
                }
                out << "\n";
            }
        }
    }

    void writeNodalDifferentialValues2(std::ostream& out)
    {
        const index_t numAct = pDomIt->computeActiveFunctions().rows();

        const gsMatrix<T>& A = pDomIt->basisValues().transpose();
        // A -- rows = quadrature points
        //   -- cols = active functions

        const gsMatrix<T>& B = pDomIt->basisDerivs(1).transpose();
        // B -- rows = quadrature points
        //   -- cols = 2 * active functions

        gsMatrix<T> nodalValues(numAct, 2 * numAct);

        const int dimension = (mBasis.degree(0) + 1) * (mBasis.degree(1) + 1);
        if (numAct > dimension) // too many functions active, do the extrapolation
        {
            // 1. construct the polynomial (find the coefficients)
            gsMatrix<T> F(pDomIt->numQuNodes(),dimension); // evaluation of the monomials in the Gauss points (rows = Gauss points, cols = monomials)

            for (index_t k = 0; k != pDomIt->numQuNodes(); ++k)
            {
                const gsVector<T> pt = pDomIt->quNodes.col(k);
                int col = 0;

                for (index_t i = 0; i < mBasis.degree(0) + 1; ++i)
                {
                    for (index_t j = 0; j < mBasis.degree(1) + 1; ++j)
                    {
                        F( k,col ) = pow(pt[0],i) * pow(pt[1],j);
                        col++;
                    }
                }
            }

            const gsMatrix<T> polyCoefs = F.fullPivHouseholderQr().solve(B);
            // polyCoeffs -- rows = monomials
            //            -- cols = 2 * active functions


            // 2. fit the polynomial with THB splines

            // get missing uniformly distributed points in the domain

            const gsMatrix<T> supp = mBasis.support();
            // supp.col(0) = left bottom corner
            // supp.col(1) = right upper corner

            const int size = mBasis.size();
            const int weneed = size - dimension;
            const gsVector<T>& lower = supp.col(0);
            const gsVector<T>& upper = supp.col(1);
            const gsMatrix<T> points = uniformPointGrid(lower, upper, weneed);

            // evaluate the polynomial and the THB basis in the new points

            gsMatrix<T> value;

            gsMatrix<T> FFup(dimension,size);
            gsMatrix<T> FFlo(weneed,size);

            FFup.setZero();

            const gsVector<unsigned>& actFunctions = pDomIt->computeActiveFunctions();

            for (index_t i = 0; i != numAct; ++i)
            {
                FFup.col(actFunctions[i]) = A.col(i);
            }

            for (index_t i = 0; i != size;++i)
            {
                for (index_t pt = 0; pt != weneed; ++pt)
                {
                    // evaluate each basis function in the uniform points
                    mBasis.evalSingle_into( i, points.col(pt) ,value );
                    FFlo(pt,i) = value(0,0);
                }
            }

            // right hand-side: upper part is B, lower part is evaluations of monomials in the new points
            gsMatrix<T> Qlo(weneed,B.cols());
            gsMatrix<T> Frhs(weneed, dimension);

            for (index_t k = 0; k < weneed; ++k)
            {
                double x = points(0,k);
                double y = points(1,k);

                int col = 0;
                for (index_t i = 0; i < mBasis.degree(0) + 1; ++i)
                {
                    for (index_t j = 0; j < mBasis.degree(1) + 1; ++j)
                    {
                        Frhs( k, col ) = pow(x,i) * pow(y,j);
                        col++;
                    }
                }
            }
            Qlo = Frhs * polyCoefs;


            // set the system with the old (Gauss) and the new points

            gsMatrix<T> FF(size,size);
            FF.block(0,0,dimension, size) = FFup;
            FF.block(dimension,0,weneed, size) = FFlo;

            gsMatrix<T> Q(size, B.cols());
            Q.block(0,0,dimension,B.cols()) = B;
            Q.block(dimension,0,weneed,B.cols()) = Qlo;

            const gsMatrix<T> thbCoefs = FF.fullPivHouseholderQr().solve(Q);

            // extract the coefficients of the active THB splines

            for (index_t i = 0; i != numAct; ++i)
            {
                nodalValues.row(i) = thbCoefs.row(actFunctions[i]);
            }

        }
        else
        {
            nodalValues = A.fullPivHouseholderQr().solve(B);
        }

        // 3. write coefficients for node k (they correspond to N_k in the expression of the derivatives dN_l)
        for (index_t k = 0; k < numAct; ++k)
        {
            gsMatrix<T> derivs = nodalValues.row(k);
            derivs.resize(2, numAct);

            for (index_t row = 0; row != derivs.rows(); ++row)
            {
                if (1e-9 < derivs.row(row).sum())
                {
                    gsWarn << "$ ** Warning \n"
                           << "$    at gauss point: " << k + 1 << "\n"
                           << "$    sum of derivatives of basis functions is not "
                        "equal to 0 in direction " << row << "\n"
                           << "$   sum is: " << derivs.row(row).sum() << "\n";
                }
            }

            out << "$ node: " << k + 1 <<"\n";
            for (index_t col = 0; col != derivs.cols(); ++col)
            {
                for (index_t row = 0; row != derivs.rows(); ++row)
                {
                    out << std::setw(20) << derivs(row, col);
                }
                out << "\n";
            }
        }
    }


    void writeCurve(std::ostream& out)
    {
        out << std::setprecision(13);

        out << "*DEFINE_CURVE\n"
             << std::setw(10) << 123 << "\n"
             << std::setw(20) << 0.0 << std::setw(20) << 0.0 << "\n"
             << std::setw(20) << 5e-3 << std::setw(20) << 1.0 << "\n";
    }

    void writeNodes(std::ostream& out)
    {
        const gsMatrix<T>& coefs = mSurf.coefs();
        gsMatrix<unsigned>* pBoundary = mBasis.boundary();
        
        int nodeId = mStartNodeId1;
        
        // width for writing id
        const int idw = 8;
        // width for writing coefficients
        const int coefsw = 16;

        int counter = 0; // counter for boundary values

        out << std::setprecision(9);
        out << std::setiosflags(std::ios::uppercase);


        out << "*NODES\n";

        // for loop over coefficients
        for (unsigned i = 0; i < static_cast<unsigned>(coefs.rows()); ++i)
        {
            out << std::setw(idw) << nodeId;
            nodeId++;

            // for loop over x, y, z
            for (int dim = 0; dim < coefs.cols(); ++dim)
            {
                out << std::setw(coefsw) << coefs(i, dim);
            }

            // here I heavily depend on fact that boundary matrix is
            // ordered
            T translConst = 0.0;
            if (counter < (*pBoundary).rows() && i == (*pBoundary)(counter))
            {
                counter++;
                translConst = 7.0;
            }

            out << std::fixed << std::setprecision(0) << std::setw(idw - 1) 
                 << translConst << "." << std::setw(idw - 1) << 0. << ".\n";
            
            out << std::scientific << std::setprecision(9);
        }            
        
        delete pBoundary;
    }


    void writePartAndSection(std::ostream& out, 
                             const bool solid = false)
    {
        const int numElements = mBasis.numElements();
        
        // width of the output
        const int w = 10;
        

        // =====================================================================
        // define all variables
        // =====================================================================
        
        // ******************
        // in case of solid == true, this variables are not used
        // ******************

        // shell tickness at all nodes
        T rthick = 0.05;

        // number of through hickness integration points
        unsigned inip = 2;

        // shear correction factor
        T rshrf = 1.0;

        // rotational dof
        unsigned iiform = 1; // whitout rotational degree of freedom

        // --------------------------------------------------------------------
        // data for first mass matrix
        // --------------------------------------------------------------------

        // material identification -first
        unsigned materialId1 = 1;

        // mass density
        T massDen1 = 1.0;

        // youngs modulus
        T youngsMod1 = 10000000;

        // poisson ratio
        T possRatio1 = 0.3;


        // --------------------------------------------------------------------
        // data for second mass matrix
        // --------------------------------------------------------------------

        // material identification -first
        unsigned materialId2 = 2;

        // mass density
        T massDen2 = 1.0e-10;

        // youngs modulus
        T youngsMod2 = 1.0e-10;

        // poisson ratio
        T possRatio2 = 1.0e-10;


        // ---------------------------------------------------------------------
        // data identification and formulation id - first
        // ---------------------------------------------------------------------

        // part identification
        int partId = mStartPartId1;

        // element formulation id
        int formulationId = mStartFormulationId1;

        // --------------------------------------------------------------------
        // data identification and formulation id - second
        // --------------------------------------------------------------------

        // part identification
        int partId2 = mStartPartId2;

        // element formulation id
        int formulationId2 = mStartFormulationId2;

        out << "$\n"
             << "$ section parameters\n"
             << "$\n"
             << "*PARAMETER\n"
             << "RTHICK, "
             << std::setprecision(3) << std::scientific << rthick << "\n"
             << "*PARAMETER\n"
             << "INIP,  " << inip << "\n"
             << "*PARAMETER\n"
             << "RSHRF, " << rshrf << "\n"
             << "*PARAMETER\n"
             << "IIFORM, " << iiform << "\n"
             << "$\n"
             << "*MAT_ELASTIC\n"
             << std::setw(w) << materialId1 <<
            std::setw(w) << massDen1 <<
            std::setw(w) << youngsMod1 <<
            std::setw(w) << possRatio1 << "\n";

        for (int el = 0; el < numElements; ++el)
        {
            
            out << "*PART\n\n"; // here we must leave a blak line
                                 // part wants to have the title

            // section id is the same as part id
            out << std::setw(w) << partId 
                << std::setw(w) << partId 
                << std::setw(w) << materialId1 << "\n";
            
            if (solid)
            {
                out << "*SECTION_SOLID\n"
                    << std::setw(w) << partId 
                    << std::setw(w) << formulationId 
                    << std::setw(w) << 0 << "\n";
            }
            else 
            {
                out << "*SECTION_SHELL\n"
                    << std::setw(w) << partId 
                    << std::setw(w) << formulationId << ",&SHRF,&NIP\n"
                    << "&THICK\n";
            }
            partId++;
            formulationId++;
        }

        out << "*MAT_ELASTIC\n"
            << std::setw(w) << materialId2 
            << std::setw(w) << massDen2 
            << std::setw(w) << youngsMod2 
            << std::setw(w) << possRatio2 << "\n"
            << "*PART\n\n"
            << std::setw(w) << partId2 
            << std::setw(w) << partId2 
            << std::setw(w) << materialId2 << "\n";
        
        if (solid)
        {
            out << "*SECTION_SOLID\n"
                << std::setw(w) << partId2 
                << std::setw(w) << formulationId2 
                << std::setw(w) << 0 << "\n";
        }
        else 
        {
            out << "*SECTION_SHELL\n"
                << std::setw(w) << partId2 
                << std::setw(w) << formulationId2 << ",&SHRF,&NIP\n"
                << "&THICK\n";
        }
    }


    void writeHeader(std::ostream& out, 
                     const std::string& name)
    {
        // keyword for allocating the memory (default was "100M")
        std::string keyword = "1000M";

        // number of eigenvalues
        unsigned nEigenvalues = 200;

        // use the consistent mass matrix
        unsigned consistentMass = 1;

        int imflag = 1;

        // Initial time step size for implicit analysis
        double DT0 = (mCrack) ? 1.0 : 0.001;

        // termination time
        double termination = (mCrack) ? 1.0 : 0.005;

        // time interval between outputs ???
        double timeInterval = (mCrack) ? 1.0 : 0.001;

        // something about saving
        double saving = (mCrack) ? 1.0 : 0.001;

        // =========================================================================
        // write
        // =========================================================================

        out << "*KEYWORD " << keyword << "\n";
        out << "*TITLE\n";
        out << name << " Mesh\n";
        if (!mCrack)
        {
            out << "*CONTROL_IMPLICIT_EIGENVALUE\n";
            out << "       " << nEigenvalues << "\n";
            out << "$ add the following line for eigenvector output:\n";
            out << "$ but it can produce large files...\n";
            out << "$0,0,0,0,0,1\n";
            out << "*CONTROL_IMPLICIT_CONSISTENT_MASS\n";
            out << consistentMass << "\n";
        }
        out << "*CONTROL_IMPLICIT_GENERAL\n";
        out << "        " << imflag << " ";
        out << std::setprecision(3) << std::scientific;
        out << DT0 << "\n";
        out << "*CONTROL_TERMINATION\n";
        out << termination << "\n";
        out << "*DATABASE_BINARY_D3PLOT\n";
        out << timeInterval << "\n";
        out << "*DATABASE_GLSTAT\n";
        out << saving << "\n";
    }           

// ----------------------------------------------------------------------
// some helper functions
// ----------------------------------------------------------------------
private:
    static void checkValuesAndDerivatives(const index_t k,
                                          const gsVector<T>& values,
                                          const gsMatrix<T>& derivs)
    {

        if (1e-9 < std::abs(values.sum() - 1))
        {
            gsWarn << "$ ** Warning \n"
                   << "$   at gauss point: " << k + 1 << "\n"
                   << "$   sum of values of basis functions is not equal to 1\n"
                   << "$   sum is: " << values.sum() << "\n";
        }

        for (int u = 0; u != derivs.rows(); ++u)
        {
            if (1e-9 < derivs.row(u).sum())
            {
                gsWarn << "$** Warning \n"
                       << "$   at gauss point: " << k + 1 << "\n"
                       << "$   sum of derivatives of basis functions is not "
                    "equal to 0 in direction " << u << "\n"
                       << "$   sum is: " << derivs.row(u).sum() << "\n";
            }
        }

    }


// ----------------------------------------------------------------------
// private data memebers
// ----------------------------------------------------------------------
private:
    
    // surface geometry -- the ls-dyna file will present this geometry
    const gsGeometry<T>& mSurf;

    // ls-dyna output file is written with this geometry
    // std::ostream out;
    
    // basis of mSurf 
    const gsBasis<T>& mBasis;
    
    // are we modeling the crack problem
    bool mCrack;

    // starting numbers for different keywords
    int mStartPartId1;
    int mStartFormulationId1;
    int mStartPartId2;
    int mStartFormulationId2;
    int mStartNodeId1;
    int mStartNodeId2;
    int mStartElementId;
    int mStartElementShellId;

    // pointer to evaluators
    gsGeometryEvaluator<T>* pGeoEval;
    gsDomainIterator<T>* pDomIt;

    // some counters
    int mFormId;
    int mElementId;
    int mNodeId2;
    int mElementShellId;

    // mappers
    // maps the parameter value to the id of the interpolation node
    std::map< std::pair<T, T>, int> mInterpolNodes;
    
    // maps interpolation node ids to inderpolation node ids
    // this map is identity, exept it maps nodes with the same
    // coordinate to one specific node
    std::map<int, int> mNodesMap;

    // through tickness basis
    gsBSplineBasis<T>* mThickBasis;

};

} // namespace gismo
