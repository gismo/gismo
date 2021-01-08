/** @file gsRemapInterface.h

    @brief Provides a mapping from the patch side of geometry one to the corresponding patch side of geometry two

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Seiler, R. Schneckenleitner
    Created on: 2018-06-12
*/


#pragma once

#include <gsIO/gsOptionList.h>

#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsAffineFunction.h>
#include <gsUtils/gsSortedVector.h>
#include <gsAssembler/gsQuadRule.h>

#include <gsNurbs/gsBSpline.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsNurbs/gsTensorNurbs.h>
#include <gsNurbs/gsTensorNurbsBasis.h>

#include <gsModeling/gsCurveFitting.h>
#include <gsModeling/gsFitting.h>

namespace gismo {


template <class T>
gsRemapInterface<T>::gsRemapInterface(const gsMultiPatch<T>   & mp,
                                      const gsMultiBasis<T>   & basis,
                                      const boundaryInterface & bi)
                                      : m_g1(mp[bi.first().patch]),
                                        m_g2(mp[bi.second().patch]),
                                        m_b1(basis[bi.first().patch]),
                                        m_b2(basis[bi.second().patch])// in most cases they are just the other way around
{
    //gsInfo << "patches: " << bi.first().patch << " and " << bi.second().patch << "\n";
    m_flipSide2 = false;
    m_isMatching = checkIfMatching();

    if (domainDim() > 2)
    {
        GISMO_ENSURE(m_isMatching == true, "Can handle non-matching interfaces only for 2 dimensions.");
    }

    //gsInfo << "equal corners: " << sameCorners << "\n";

    if (!m_isMatching)
    {
        //gsInfo << "the follwing patches do not match: " << m_g1.id() << " and " << m_g2.id() << "\n";
        findInterface(bi);
        changeDir(bi);
    }
    else
    {

        m_side1 = bi.first();
        m_side2 = bi.second();

        std::vector<boxCorner> corners;
        gsMatrix<T> inversCorners;

        m_parameterbounds.first.resize(domainDim(), 2);
        m_parameterbounds.second.resize(domainDim(), 2);

        //gsInfo << "the follwing patches match: " << m_g1.id() << " and " << m_g2.id() << "\n";

        switch (m_side1.index()) {
            case 1:
                m_parameterbounds.first.row(1) = mp[bi.first().patch].parameterRange().row(1);
                m_parameterbounds.first(0, 0) = //mp[bi.first().patch].parameterRange()(0, 0);
                m_parameterbounds.first(0, 1) = mp[bi.first().patch].parameterRange()(0, 0);

                if(domainDim() == 3)
                    m_parameterbounds.first.row(2) = mp[bi.first().patch].parameterRange().row(2);

                break;
            case 2:
                m_parameterbounds.first.row(1) = mp[bi.first().patch].parameterRange().row(1);
                m_parameterbounds.first(0, 0) = //mp[bi.first().patch].parameterRange()(0, 1);
                m_parameterbounds.first(0, 1) = mp[bi.first().patch].parameterRange()(0, 1);

                if(domainDim() == 3)
                    m_parameterbounds.first.row(2) = mp[bi.first().patch].parameterRange().row(2);

                break;
            case 3:
                m_parameterbounds.first.row(0) = mp[bi.first().patch].parameterRange().row(0);
                m_parameterbounds.first(1, 0) = //mp[bi.first().patch].parameterRange()(1, 0);
                m_parameterbounds.first(1, 1) = mp[bi.first().patch].parameterRange()(1, 0);

                if(domainDim() == 3)
                    m_parameterbounds.first.row(2) = mp[bi.first().patch].parameterRange().row(2);

                break;
            case 4:
                m_parameterbounds.first.row(0) = mp[bi.first().patch].parameterRange().row(0);
                m_parameterbounds.first(1, 0) = //mp[bi.first().patch].parameterRange()(1, 1);
                m_parameterbounds.first(1, 1) = mp[bi.first().patch].parameterRange()(1, 1);

                if(domainDim() == 3)
                    m_parameterbounds.first.row(2) = mp[bi.first().patch].parameterRange().row(2);

                break;
            case 5: // for 3D case
                m_parameterbounds.first.row(0) = mp[bi.first().patch].parameterRange().row(0);
                m_parameterbounds.first.row(1) = mp[bi.first().patch].parameterRange().row(1);
                m_parameterbounds.first(2, 0) = //mp[bi.first().patch].parameterRange()(2, 0);
                m_parameterbounds.first(2, 1) = mp[bi.first().patch].parameterRange()(2, 0);
                break;
            case 6: // for 3D case
                m_parameterbounds.first.row(0) = mp[bi.first().patch].parameterRange().row(0);
                m_parameterbounds.first.row(1) = mp[bi.first().patch].parameterRange().row(1);
                m_parameterbounds.first(2, 0) = //mp[bi.first().patch].parameterRange()(2, 1);
                m_parameterbounds.first(2, 1) = mp[bi.first().patch].parameterRange()(2, 1);
                break;
        }

        switch (m_side2.index()) {

            case 1:
                m_parameterbounds.second.row(1) = mp[bi.second().patch].parameterRange().row(1);
                m_parameterbounds.second(0, 0) = //mp[bi.second().patch].parameterRange()(0, 0);
                m_parameterbounds.second(0, 1) = mp[bi.second().patch].parameterRange()(0, 0);

                if(domainDim() == 3)
                    m_parameterbounds.second.row(2) = mp[bi.second().patch].parameterRange().row(2);

                break;
            case 2:
                m_parameterbounds.second.row(1) = mp[bi.second().patch].parameterRange().row(1);
                m_parameterbounds.second(0, 0) = //mp[bi.second().patch].parameterRange()(0, 1);
                m_parameterbounds.second(0, 1) = mp[bi.second().patch].parameterRange()(0, 1);

                if(domainDim() == 3)
                    m_parameterbounds.second.row(2) = mp[bi.second().patch].parameterRange().row(2);

                break;
            case 3:
                m_parameterbounds.second.row(0) = mp[bi.second().patch].parameterRange().row(0);
                m_parameterbounds.second(1, 0) = //mp[bi.second().patch].parameterRange()(1, 0);
                m_parameterbounds.second(1, 1) = mp[bi.second().patch].parameterRange()(1, 0);

                if(domainDim() == 3)
                    m_parameterbounds.second.row(2) = mp[bi.second().patch].parameterRange().row(2);

                break;
            case 4:
                m_parameterbounds.second.row(0) = mp[bi.second().patch].parameterRange().row(0);
                m_parameterbounds.second(1, 0) = //mp[bi.second().patch].parameterRange()(1, 1);
                m_parameterbounds.second(1, 1) = mp[bi.second().patch].parameterRange()(1, 1);

                if(domainDim() == 3)
                    m_parameterbounds.second.row(2) = mp[bi.second().patch].parameterRange().row(2);

                break;
            case 5: // for 3D case
                m_parameterbounds.second.row(0) = mp[bi.second().patch].parameterRange().row(0);
                m_parameterbounds.second.row(1) = mp[bi.second().patch].parameterRange().row(1);
                m_parameterbounds.second(2, 0) = //mp[bi.second().patch].parameterRange()(2, 0);
                m_parameterbounds.second(2, 1) = mp[bi.second().patch].parameterRange()(2, 0);
                break;
            case 6: // for 3D case
                m_parameterbounds.second.row(0) = mp[bi.second().patch].parameterRange().row(0);
                m_parameterbounds.second.row(1) = mp[bi.second().patch].parameterRange().row(1);
                m_parameterbounds.second(2, 0) = //mp[bi.second().patch].parameterRange()(2, 1);
                m_parameterbounds.second(2, 1) = mp[bi.second().patch].parameterRange()(2, 1);
                break;
        }
        //gsInfo << "parameterRange: \n" << mp[bi.second().patch].parameterRange() << "\n";
        //gsInfo << "parameter bound one : \n" << m_parameterbounds.first << "\n";
        //gsInfo << "parameter bound two: \n" << m_parameterbounds.second << "\n";
        //gsInfo << "side: \n" << m_side2.index() << "\n";

        changeDir(bi);
    }

    constructReparam();
    if(!m_isMatching)
        constructBreaks();
}


template<class T>
const std::vector<T> gsRemapInterface<T>::getPointsOnInterface() const
{
    const short_t fixedDir = m_side1.direction();
    std::vector<T> vec(m_breakpoints.cols());

    Eigen::Map<Eigen::Matrix<T, 1, -1> >(&vec[0], m_breakpoints.cols()) = m_breakpoints.row(fixedDir == 1 ? 0 : 1);

    return vec;
}

template<class T>
void gsRemapInterface<T>::constructBreaks() {
    // computes break points per element

    //const gsBasis <T> &B1 = m_g1.basis();
    //const gsBasis <T> &B2 = m_g2.basis();

    // Get the interface of the patches
    // TODO: can be moved to the constructor?!
    gsMultiPatch<T> firstPatch(m_g1);
    gsMultiPatch<T> secondPatch(m_g2);
    firstPatch.computeTopology(); secondPatch.computeTopology();
    std::vector<patchSide > boundariesPatch1 = firstPatch.boundaries();
    std::vector<patchSide > boundariesPatch2 = secondPatch.boundaries();
    patchSide patchSide1, patchSide2;

    typename gsBasis<T>::domainIter domIt1, domIt2;
    gsMatrix<T> startPatch1, startPatch2;

    // Check which sides contain a part of the common interface
    for(size_t i = 0; i < boundariesPatch1.size(); i++)
        if(boundariesPatch1[i].index() == m_side1.index())
        {
            domIt1 = m_b1.makeDomainIterator( boundariesPatch1[i] );
            patchSide1 = boundariesPatch1[i];
            startPatch1 = m_parameterbounds.first.col(0);
        }

    for(size_t i = 0; i < boundariesPatch2.size(); i++)
        if(boundariesPatch2[i].index() == m_side2.index())
        {
            domIt2 = m_b2.makeDomainIterator(boundariesPatch2[i]);
            patchSide2 = boundariesPatch2[i];
            if(m_flipSide2)
                startPatch2 = m_parameterbounds.second.col(1);
            else
                startPatch2 = m_parameterbounds.second.col(0);
        }


    // Compute interface knots in physical domain by evaluating left and right geometry maps at the knot values
    size_t numelP1 = domIt1->numElements();
    size_t numelP2 = domIt2->numElements();
    gsMatrix <T> physicalKnotsP1(m_g1.geoDim(), numelP1 + 1), physicalKnotsP2(m_g2.geoDim(), numelP2 + 1), dummy;

    domIt1->reset();
    domIt2->reset();
    index_t numBreaksPatch1 = 1, numBreaksPatch2 = 1; // vars to count the entries in the physical breakpoints

    // evaluate the first point of the interface
    m_g1.eval_into(startPatch1, dummy);
    physicalKnotsP1.col(0) = dummy;

    // loop over all elements of the boundary with interface part, but evaluate only element corners on the real interface
    for (; domIt1->good(); domIt1->next())
    {
        if(m_side1.index() == 3 || m_side1.index() == 4) // v is fix
        {
            if(domIt1->lowerCorner()(0,0) > startPatch1(0,0) && domIt1->lowerCorner()(0,0) <= m_parameterbounds.first(0,1))
            {
                // An ansatz for 3D??
                //if(domainDim() == 2 || (domainDim() == 3 &&
                //                        domIt1->lowerCorner()(2,0) > startPatch1(2,0) &&
                //                        domIt1->lowerCorner()(2,0) <= m_parameterbounds.first(2,1)))
                //{
                m_g1.eval_into(domIt1->lowerCorner(), dummy);
                physicalKnotsP1.col(numBreaksPatch1) = dummy;
                numBreaksPatch1++;
                //}
            }
        }
        else
        {
            if(m_side1.index() == 1 || m_side1.index() == 2) // u is fix
            {
                if(domIt1->lowerCorner()(1,0) > startPatch1(1,0) && domIt1->lowerCorner()(1,0) <= m_parameterbounds.first(1,1))
                {
                    m_g1.eval_into(domIt1->lowerCorner(), dummy);
                    physicalKnotsP1.col(numBreaksPatch1) = dummy;
                    numBreaksPatch1++;
                }
            }
            // for the 3D case??
            //else // w is fix
            //{
            //    if((domIt1->lowerCorner()(0,0) > startPatch1(0,0) && domIt1->lowerCorner()(0,0) <= m_parameterbounds.first(0,1))
            //            && (domIt1->lowerCorner()(1,0) > startPatch1(1,0) && domIt1->lowerCorner()(1,0) <= m_parameterbounds.first(1,1)))
            //    {
            //        m_g1.eval_into(domIt1->lowerCorner(), dummy);
            //        physicalKnotsP1.col(numBreaksPatch1) = dummy;
            //        numBreaksPatch1++;
            //    }
            //}
        }
        //domIt1->next();
    }

    // evaluate the last point of the interface, i.e., this last point must also be within the parameter bound
    if(m_side1.index() == 3 || m_side1.index() == 4)
    {
        if(domIt1->upperCorner()(0,0) <= m_parameterbounds.first(0,1))
        {
            m_g1.eval_into(domIt1->upperCorner(), dummy);
            physicalKnotsP1.col(numBreaksPatch1) = dummy;
            numBreaksPatch1++;
        }
    }
    else
    {
        if(m_side1.index() == 1 || m_side1.index() == 2)
        {
            if(domIt1->upperCorner()(1,0) <= m_parameterbounds.first(1,1))
            {
                m_g1.eval_into(domIt1->upperCorner(), dummy);
                physicalKnotsP1.col(numBreaksPatch1) = dummy;
                numBreaksPatch1++;
            }
        }
    }
    //gsInfo << "physical knots 1: \n" << physicalKnotsP1 << "\n";



    /// do the same for patch 2 as above
    m_g2.eval_into(startPatch2, dummy);
    physicalKnotsP2.col(0) = dummy;

    for (; domIt2->good(); domIt2->next()) // for (index_t i = 0; i < numelP2; i++)
    {
        if(m_side2.index() == 3 || m_side2.index() == 4)
        {
            if (domIt2->lowerCorner()(0,0) > startPatch2(0,0) && domIt2->lowerCorner()(0,0) < std::max(m_parameterbounds.second(0,1), m_parameterbounds.second(0,0)))
            {
                m_g2.eval_into(domIt2->lowerCorner(), dummy);
                physicalKnotsP2.col(numBreaksPatch2) = dummy;
                numBreaksPatch2++;
            }
        }
        else
        {
            if(m_side2.index() == 1 || m_side2.index() == 2)
            {
                if (domIt2->lowerCorner()(1,0) > startPatch2(1,0) && domIt2->lowerCorner()(1,0) < std::max(m_parameterbounds.second(1,1), m_parameterbounds.second(1,0)))
                {
                    m_g2.eval_into(domIt2->lowerCorner(), dummy);
                    physicalKnotsP2.col(numBreaksPatch2) = dummy;
                    numBreaksPatch2++;
                }
            }
        }
        //domIt2->next();
    }

    // add only the breakpoints within the parameter bounds
    if(m_side2.index() == 3 || m_side2.index() == 4)
    {
        if(domIt2->upperCorner()(0,0) <= std::max(m_parameterbounds.second(0,1), m_parameterbounds.second(0,0)))
        {
            m_g2.eval_into(domIt2->upperCorner(), dummy);
            physicalKnotsP2.col(numBreaksPatch2) = dummy;
            numBreaksPatch2++;// to get the number of entries
        }
    }
    else
    {
        if(m_side2.index() == 1 || m_side2.index() == 2)
        {
            if(domIt2->upperCorner()(1,0) <= std::max(m_parameterbounds.second(1,1), m_parameterbounds.second(1,0)))
            {
                m_g2.eval_into(domIt2->upperCorner(), dummy);
                physicalKnotsP2.col(numBreaksPatch2) = dummy;
                numBreaksPatch2++;
            }
        }
    }
    //gsInfo << "physical knots 2: \n" << physicalKnotsP2 << "\n";

    /// store all the physical points in one vector
    gsMatrix<T> physicalBreaks(domainDim(), numBreaksPatch1+numBreaksPatch2); // Assume m_g1.geoDim() == m_g2.geoDim()

    for(index_t c = 0; c < numBreaksPatch1; c++)
        physicalBreaks.col(c) = physicalKnotsP1.col(c);

    for(index_t c = 0; c < numBreaksPatch2; c++)
        physicalBreaks.col(numBreaksPatch1+c) = physicalKnotsP2.col(c);

    // compute the corresponding parameter values in one patch, here of patch1
    if(domainDim() == 2)
    {
        gsSortedVector<T> parameterBreaks;

        // Determine fixed coordinate of patch2 -> Use here patch2 because we compute the Interfacemap of patch1!!!
        // fixedDir ==  0 corresponds to fixed u and 1 corresponds to a fixed v
        index_t fixedDir = patchSide1.direction();

        gsMatrix<T> G2_parametric_LC;
        for (index_t i = 0; i < physicalBreaks.cols(); i++) {
            // computes the preimages of the breakpoints for each of the two patches
            m_g1.invertPoints(physicalBreaks.col(i), G2_parametric_LC); // not exact, we have rounding errors
            // taking care of the rounding errors by iterating over the vector and checking the absolute value between the current
            // preimage and the already available ones
            if (fixedDir == 1) {
                if (parameterBreaks.size() == 0)
                    parameterBreaks.push_sorted_unique(G2_parametric_LC(0, 0)); // sort w.r.t. u direction
                else {
                    index_t j = 0;
                    gsVector<bool> roundingError = gsVector<bool>::Constant(parameterBreaks.size(), true);

                    for (typename gsSortedVector<T>::iterator it = parameterBreaks.begin();
                         it != parameterBreaks.end(); it++) {

                        if (math::abs(G2_parametric_LC(0, 0) - *it) > 1.e-4) {
                            roundingError(j) = false;
                        }
                        j++;
                    }
                    if (( false == roundingError.array() ).all())
                        parameterBreaks.push_sorted_unique(G2_parametric_LC(0, 0)); // sort w.r.t. u direction


                }
            } else {
                if (parameterBreaks.size() == 0)
                    parameterBreaks.push_sorted_unique(G2_parametric_LC(1, 0)); // sort w.r.t. v direction
                else {
                    index_t j = 0;
                    gsVector<bool> roundingError = gsVector<bool>::Constant(parameterBreaks.size(), true);

                    for (typename gsSortedVector<T>::iterator it = parameterBreaks.begin();
                         it != parameterBreaks.end(); it++) {
                        if (math::abs(G2_parametric_LC(1, 0) - *it) > 1.e-4)
                            roundingError(j) = false;

                        j++;
                    }
                    if (( false == roundingError.array() ).all())
                        parameterBreaks.push_sorted_unique(G2_parametric_LC(1, 0)); // sort w.r.t. v direction
                }
            }

        }

        m_breakpoints = gsMatrix<T>(domainDim(), parameterBreaks.size()); // Assume m_g1.geoDim() == m_g2.geoDim()
        for (size_t i = 0; i < parameterBreaks.size(); i++) {
            if (fixedDir)
                m_breakpoints.col(i) << parameterBreaks[i], G2_parametric_LC(1, 0);
            else
                m_breakpoints.col(i) << G2_parametric_LC(0, 0), parameterBreaks[i];
        }

        // only for tests
        //gsMatrix<T> result;
        //this->eval_into(m_breakpoints, result);

        //gsInfo << "Mapped: \n" << result << "\n";
    }
    else
    {
        // ?
    }

}


template<class T>
void gsRemapInterface<T>::constructReparam()
{
    const index_t numIntervals = 11; // ?
    const index_t numGeometries = 2;


    if(m_isMatching)
    {
        boundaryInterface iFace(m_side1, m_side2, domainDim());

        gsAffineFunction<T> interfaceMap(iFace.dirMap(m_side1), iFace.dirOrientation(m_side1), m_parameterbounds.first, m_parameterbounds.second);

        m_fittedInterface = interfaceMap.clone();
    }
    else // if the patches are not matching then do the fitting process
    {

        // Assume tensor structure
        // now create samples for both patches
        // the knot intervals can be different, e.g.,
        //----------------------------
        //-                          -
        //----------------------------
        //    --------
        //    -      -
        //    --------
        gsMatrix<T> t_vals = gsMatrix<T>::Zero(numGeometries, numIntervals);
        T firstKnot, lastKnot;
        gsVector<T> upper(1), lower(1);
        gsVector<unsigned> numPoints(1);
        numPoints << numIntervals;

        //gsInfo << "parameterbounds: \n" << m_parameterbounds.second << "\n";
        //gsInfo << "patch: \n" << m_g2.id() << "\n";
        for (index_t np = 0; np < numGeometries; np++) {
            if (np == 0) {
                if (m_side1.index() == 3 || m_side1.index() == 4) // v is fixed
                {
                    firstKnot = m_parameterbounds.first(0, 0);
                    lastKnot = m_parameterbounds.first(0, 1);
                } else // u is fixed
                {
                    firstKnot = m_parameterbounds.first(1, 0);
                    lastKnot = m_parameterbounds.first(1, 1);
                }
            } else {
                if (m_side2.index() == 3 || m_side2.index() == 4) // v is fixed
                {
                    firstKnot = m_parameterbounds.second(0, 0);
                    lastKnot = m_parameterbounds.second(0, 1);
                } else // u is fixed
                {
                    firstKnot = m_parameterbounds.second(1, 0);
                    lastKnot = m_parameterbounds.second(1, 1);
                }
            }

            lower(0) = firstKnot;
            upper(0) = lastKnot;

            //gsInfo << "lower:\n" << firstKnot << "\n upper:\n" << lastKnot << std::endl;

            //t_vals.row(np) = uniformPointGrid(lower, upper, numIntervals); // uniformly distributed samples between the overlapping part of the interface
            t_vals.row(np) = gsPointGrid(lower, upper, numPoints);

        }

        gsMatrix<T> samples_left, samples_right;
        gsMatrix<T> find_start_value;

        // Get the corresponding edges
        //Edge 1, {(u,v) : u = 0}
        //Edge 2, {(u,v) : u = 1}
        //Edge 3, {(u,v) : v = 0}
        //Edge 4, {(u,v) : v = 1}

        //gsInfo << "left boundary: " << m_interfacePatch1 << "\n";
        //gsInfo << "right boundary: " << m_interfacePatch2 << "\n";

        //gsMatrix<T> vals2dPatch1(t_vals.rows()+1, t_vals.cols()), vals2dPatch2(t_vals.rows()+1, t_vals.cols());
        // TODO: use already available information
        gsMatrix<T> vals2dPatch1, vals2dPatch2;
        enrichToVector(m_side1.index(), m_g1, t_vals.row(0), vals2dPatch1);
        enrichToVector(m_side2.index(), m_g2, t_vals.row(1), vals2dPatch2);

        m_g1.eval_into(vals2dPatch1, samples_left);
        m_g2.eval_into(vals2dPatch2, samples_right);

        //gsInfo << "vals2dPatch1:\n" << GEO_L_ref.coefs() << "\n vals2dPatch2:\n" << GEO_R_ref.coefs() << std::endl;
        //gsInfo << "vals2dPatch1:\n" << vals2dPatch1 << "\n vals2dPatch2:\n" << vals2dPatch2 << std::endl;
        //std::cout << "samples left:\n" << samples_left << "\n samples right:\n" << samples_right << std::endl;

        gsMatrix<T> B(numIntervals, m_g1.geoDim());

        for (index_t i = 0; i < t_vals.cols(); i++) {
            // find a suitable start value for the Newton iteration
            find_start_value = (samples_right.colwise()) - samples_left.col(i);

            size_t row, col;

            find_start_value.colwise().squaredNorm().minCoeff(&row, &col);

            gsVector<T> b_null = samples_right.col(col);

            // Pass on g2 if one wants to find a mapping from interface1 to interface2
            //gsMatrix<T> b = closestPoint(b_null, g2, samples_left.col(i));

            // this gives the same result as above
            m_g2.newtonRaphson(samples_left.col(i), b_null, true, 10e-6, 100);
            //gsInfo << "newton: " << b_null << "\n";

            // TODO: Check if the order of the coefficients has an impact on the mapping regarding assembling aso.
            B.row(i) = b_null.transpose(); // to be in the correct order

        }

        // the coefficients to fit
        //std::cout << "B:\n" << B << std::endl;

        // check the error
        // assume that the right map is the identity
        gsMatrix<T> eval_orig, eval_fit, B2, id;

        gsKnotVector<T> KV(t_vals(0, 0), t_vals(0, numIntervals - 1), 5, 4);

        gsCurveFitting<T> fit(t_vals.row(0).transpose(), B, KV);

        fit.compute();
        m_fittedInterface = fit.curve().clone();
        std::cout << "Hi, I'm the resulting curve: \n" << *m_fittedInterface << std::endl;

        unsigned errorInterval = 10;
        gsVector<unsigned > errorSamples(1);
        errorSamples << errorInterval;

        gsMatrix<T> eval_points;// = gsMatrix<T>::Zero(numGeometries, errorInterval);

        for (index_t np = 0; np < numGeometries; np++) {
            gsVector<T> lowerVal(1), upperVal(1);
            lowerVal << t_vals(np, 0);
            upperVal << t_vals(np, numIntervals - 1);
            //gsMatrix<T> grid = uniformPointGrid(lowerVal, upperVal, errorInterval);
            gsMatrix<T> grid = gsPointGrid(lowerVal, upperVal, errorSamples);
            eval_points.conservativeResize(np + 1, grid.cols()); // to check the error
            eval_points.row(np) = grid;
        }

        m_fittedInterface->eval_into(eval_points.row(0), eval_fit);
        //eval_fit(0,0) -= 0.0001; // do a nasty slight correction since the first entry is out of the domain of definition due to rounding errors

        // TODO: also here use already available information
        enrichToVector(m_side2.index(), m_g2, eval_points.row(1), id);

        m_g2.eval_into(eval_fit, eval_orig);
        m_g2.eval_into(id, B2);
        //gsInfo << "b2: \n" << id.transpose() << " and eval_orig: \n" << eval_fit.transpose() << "\n";

        // do test
        /*
        for(index_t c = 0; c < eval_fit.cols(); c++)
        {
            if(std::isnan(eval_orig.col(c).squaredNorm()))
            {
                switch (m_side2.index())
                {
                    case 1 :
                        // u = value of the first knot in u direction
                        eval_fit(0, c) += 10e-4;
                        eval_fit(1, c) += 10e-4;
                        break;
                    case 2 :
                        //u = value of the last knot in u direction
                        eval_fit(0, c) -= 10e-4;
                        eval_fit(1, c) += 10e-4;
                        break;
                    case 3 :
                        //v = value of the first knot in v direction;
                        eval_fit(0, c) += 10e-4;
                        eval_fit(1, c) += 10e-4;
                        break;
                    case 4 :
                        //v = value of the last knot in v direction
                        eval_fit(0, c) += 10e-4;
                        eval_fit(1, c) -= 10e-4;
                        break;
                }
                m_g2.eval_into(eval_fit, eval_orig);
            }
        }
*/
        //end test

        T error = 0;

        for (index_t i = 0; i < eval_points.cols(); i++)
            error += (id.col(i) - eval_fit.col(i)).squaredNorm();
            //error += (eval_orig.col(i) - B2.col(i)).squaredNorm();

        error = math::sqrt(error);

        //if(error > 0.5)
        //    gsInfo << "patch 1: \n" << eval_orig << " and patch 2: \n" << B2 << "\n";

        std::cout << "Error: " << error << std::endl;
    }

}

template <typename T>
void gsRemapInterface<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT(u.rows() == domainDim(), "The rows of the evaluation points must be equal to the dimension of the domain!");
    //result.reshape(u.rows(), u.cols());

    index_t fixedDir = m_side1.direction();

    if(m_isMatching)
    {
        m_fittedInterface->eval_into(u, result);
    }
    else
    {
        // v is fixed => loop over u values
        //gsInfo << "the patches: " << m_g1.id() << " and " << m_g2.id() << "\n";
        m_fittedInterface->eval_into(checkIfInBound(u.row(!fixedDir)), result); /// ????

        const short_t m_side2dir = m_side2.direction();

        //gsInfo << "result before: \n" << result << "\n";
        // need here the second basis since result store points in the second geometry
        if(const gsTensorBSplineBasis<2, T> *tb = dynamic_cast<const gsTensorBSplineBasis<2, T> * >(&m_g2.basis()))
        {
            /*if(m_side2.direction() == 1 && m_side2.parameter() == 0) // v = 0
                result.row(1).setConstant(tb->knots(1).first());

            if(m_side2.direction() == 1 && m_side2.parameter() == 1) // v = 1
                result.row(1).setConstant(tb->knots(1).last());

            if(m_side2.direction() == 0 && m_side2.parameter() == 0) // u = 0
                result.row(0).setConstant(tb->knots(0).first());

            if(m_side2.direction() == 0 && m_side2.parameter() == 1) // u = 1
                result.row(0).setConstant(tb->knots(0).last());*/
            if(m_side2dir == 0 || m_side2dir == 1)
            {
                result.row(m_side2dir).setConstant(m_side2.parameter() ?
                                                   tb->knots(m_side2dir).last() :
                                                   tb->knots(m_side2dir).first());
            }

            return;
        }

        if(const gsTensorNurbsBasis<2, T> * ntb = dynamic_cast<const gsTensorNurbsBasis<2, T> * >(&(m_g2.basis())))
        {
            /*if(m_side2.direction() == 1 && m_side2.parameter() == 0) // v = 0
                result.row(1).setConstant(ntb->source().knots(1).first());

            if(m_side2.direction() == 1 && m_side2.parameter() == 1) // v = 1
                result.row(1).setConstant(ntb->source().knots(1).last());

            if(m_side2.direction() == 0 && m_side2.parameter() == 0) // u = 0
                result.row(0).setConstant(ntb->source().knots(0).first());

            if(m_side2.direction() == 0 && m_side2.parameter() == 1) // u = 1
                result.row(0).setConstant(ntb->source().knots(0).last());*/
            if(m_side2dir == 0 || m_side2dir == 1)
            {
                result.row(m_side2dir).setConstant(m_side2.parameter() ?
                                                   ntb->source().knots(m_side2dir).last() :
                                                   ntb->source().knots(m_side2dir).first());
            }

            return;
        }

        GISMO_ERROR("Unfitted interface not supported");
    }

}

template<class T>
memory::unique_ptr< gsDomainIterator<T> > gsRemapInterface<T>::makeDomainIterator() const
{
    if (m_isMatching) return m_b1.makeDomainIterator(m_side1);

    gsTensorDomainBoundaryIterator<T> * tdi = new gsTensorDomainBoundaryIterator<T> (m_b1, m_side1);

    std::vector<T> newBreaks = getPointsOnInterface();
    gsInfo << "newBreaks: \n";
    for(index_t i = 0; i < m_breakpoints.cols(); i++)
        gsInfo << newBreaks[i] << "\t";

    gsInfo << "\n";

    // the input must be the direction which is moving
    //tdi->setBreaks(newBreaks, m_side1.direction()); -> gives the fixed direction
    // workaround: only works for 2 dimensions

    if (m_side1.direction() == 1)
        tdi->setBreaks(newBreaks, 0);
    else //m_side1.direction() == 0
        tdi->setBreaks(newBreaks, 1);


    return domainIterUPtr(tdi);
}

// Function to enhance a sequence of 1D points in an interval to 2D points in the parameter domain
// A matrix pts with the already correct dimensions is expected
// boundarySide is the boundary side according to gismo's boundary::side
// intervals is the number of samples to enrich to more dimensions
// pts is the matrix populated with the corresponding points for more dimensions
// only works for 2d at the moment!!
// TODO: Generalize for arbitrary dimensions
template<class T>
void gsRemapInterface<T>::enrichToVector(const short_t         boundarySide,
                                         const gsGeometry<T> & geo,
                                         const gsMatrix<T>   & intervals,
                                               gsMatrix <T>  & pts)
{
    pts.resize(geo.geoDim(), intervals.cols());

    //const gsTensorBSplineBasis<2, T> *tb = dynamic_cast<const gsTensorBSplineBasis<2, T> * >(&geo.basis());

    //if(tb == NULL)
    //const gsTensorNurbsBasis<2, T> * tb = dynamic_cast<const gsTensorNurbsBasis<2, T> * >(&(geo.basis()));


    if(const gsTensorBSplineBasis<2, T> *tb = dynamic_cast<const gsTensorBSplineBasis<2, T> * >(&geo.basis()))
    {
        switch (boundarySide)
        {
            case 1 :
                // u = value of the first knot in u direction
                pts.row(0) = gsMatrix<T>::Constant(1, intervals.cols(), tb->knots(0).first());
                pts.row(1) = intervals;
                break;
            case 2 :
                //u = value of the last knot in u direction
                pts.row(0) = gsMatrix<T>::Constant(1, intervals.cols(), tb->knots(0).last());
                pts.row(1) = intervals;
                break;
            case 3 :
                //v = value of the first knot in v direction;
                pts.row(0) = intervals;
                pts.row(1) = gsMatrix<T>::Constant(1, intervals.cols(), tb->knots(1).first());
                break;
            case 4 :
                //v = value of the last knot in v direction
                pts.row(0) = intervals;
                pts.row(1) = gsMatrix<T>::Constant(1, intervals.cols(), tb->knots(1).last());
                break;
        }
    }
    else
    {
        const gsTensorNurbsBasis<2, T> * ntb = dynamic_cast<const gsTensorNurbsBasis<2, T> * >(&(geo.basis()));

        switch (boundarySide)
        {
            case 1 :
                // u = value of the first knot in u direction
                pts.row(0) = gsMatrix<T>::Constant(1, intervals.cols(), ntb->source().knots(0).first());
                pts.row(1) = intervals;
                break;
            case 2 :
                //u = value of the last knot in u direction
                pts.row(0) = gsMatrix<T>::Constant(1, intervals.cols(), ntb->source().knots(0).last());
                pts.row(1) = intervals;
                break;
            case 3 :
                //v = value of the first knot in v direction;
                pts.row(0) = intervals;
                pts.row(1) = gsMatrix<T>::Constant(1, intervals.cols(), ntb->source().knots(1).first());
                break;
            case 4 :
                //v = value of the last knot in v direction
                pts.row(0) = intervals;
                pts.row(1) = gsMatrix<T>::Constant(1, intervals.cols(), ntb->source().knots(1).last());
                break;
        }
    }
}

// Member to find the interface of the 2 input geometries
template<class T>
void gsRemapInterface<T>::findInterface()
{
    //const gsGeometry<T> & patch1 = m_domain.patch(0); // m_g1
    //const gsGeometry<T> & patch2 = m_domain.patch(1); // m_g2

    gsMultiPatch<T> geo1(m_g1); geo1.computeTopology();
    gsMultiPatch<T> geo2(m_g2); geo2.computeTopology();

    std::vector<patchSide> sidesPatch1 = geo1.boundaries();
    std::vector<patchSide> sidesPatch2 = geo2.boundaries();

    // vectors to store the rows indices of the coefficient matrices
    gsVector<index_t > boundaryPatch1, boundaryPatch2;

    // first find the sides of the patches which belong to the interface
    const index_t nCorners = 1<<m_g1.geoDim();
    gsMatrix<T> inversMaps(m_g1.geoDim(), 4); // matrix to store the preimages of the corresponding verices on the interface
    gsVector<index_t> corners(2); // vector to store the index of the corner which lies on the other patch, 0-th entry is the index of the corner for the first patch and vice versa
    corners.setZero();

    bool completeOnPatch2 = false, completeOnPatch1 = false; // check if one side of the patches is completely contained in the other side

    // two columns for the lower and the upper bound
    m_parameterbounds.first.resize(m_g1.geoDim(), 2);
    m_parameterbounds.second.resize(m_g1.geoDim(), 2);

    // matrix to store the coefficients of the corner values
    gsMatrix<T> c = gsMatrix<T>::Zero(1, m_g1.geoDim());
    gsMatrix<T> preIm;

    // find the side for the second patch, and the corner(s) of the first patch which lies on the interface
    gsVector<bool> onGeo;

    // this is the case if there are no matching corners
    for(index_t i = 1; i <= nCorners; i++)
    {
        c = m_g1.coefAtCorner(i).transpose();

        // Be aware of two matching vertices-> due to rounding errors on of the vertices can be on two sides maybe!
        gsMatrix<T> parLoc;

        std::vector<boxSide> side = m_g2.locateOn(c, onGeo, parLoc, true);

        //for(size_t i = 0; i < side.size(); i++)
        //    gsInfo << "boundary patch: " << side[i] << "\n";

        if(onGeo(0) == true)
        {
            //inversMaps.col(1) = parLoc;

            if(!corners(0))
            {
                inversMaps.col(1) = parLoc;
                corners(0) = i;
            }
            else
            {
                // in case 2 corners lie on the interface of patch 2
                inversMaps.col(3) = parLoc;
                completeOnPatch2 = true; // side of patch1 is completely contained in the corresponding side of patch2 or the patches are fully matching
            }

            m_side2.m_index = side[0].index();
            m_side2.patch = m_g2.id();
        }

    }

    // do the same for the first patch, i.e., finding the side for patch1 and the corner(s) for patch2
    for(index_t i = 1; i <= nCorners; i++)
    {
        c = m_g2.coefAtCorner(i).transpose();

        gsMatrix<T> parLoc;
        std::vector<boxSide> side = m_g1.locateOn(c, onGeo, parLoc, true);
        //for(size_t i = 0; i < side.size(); i++)
        //    gsInfo << " side: \n" << side[i] << "\n";

        if(onGeo(0) == 1)
        {
            //inversMaps.col(0) = parLoc;

            if(!corners(1))
            {
                inversMaps.col(0) = parLoc;
                corners(1) = i;
            }
            else
            {   // in case 2 corners lie on the interface of patch 1
                inversMaps.col(2) = parLoc;
                completeOnPatch1 = true;
            }

            m_side1.m_index = side[0].index();
            m_side1.patch = m_g1.id();
        }

    }



    // store the parametric bounds of the overlap for each patch, order the points w.r.t. the non fixed boundary side
    // so far the code is not very nice -> room for improvement!!!
    // Maybe considering the exact values instead of approximations would be better
    gsMatrix<T> parIm;
    if(completeOnPatch2 == false && corners(0) != 0) // if the interface does not overlap entirely or if the side of patch2 is not a proper subset of the corresponding side of patch1
    {
        m_g1.invertPoints(m_g1.coefAtCorner(corners(0)).transpose(), parIm);

        index_t bound = 0;
        if(parIm.isApprox(inversMaps.col(0), 1e-6))
            bound = 2;

        if(m_side1.index() == 3 || m_side1.index() == 4) //sort w.r.t u
        {
            if(parIm(0,0) < inversMaps.col(bound)(0,0))
            {
                m_parameterbounds.first.col(0) = parIm;
                m_parameterbounds.first.col(1) = inversMaps.col(bound);
            }
            else
            {
                m_parameterbounds.first.col(1) = parIm;
                m_parameterbounds.first.col(0) = inversMaps.col(bound);
            }
        }
        else
        {
            if(m_side1.index() == 1 || m_side1.index() == 2) // sort w.r.t v
            {
                if(parIm(1,0) < inversMaps.col(bound)(1,0))
                {
                    m_parameterbounds.first.col(0) = parIm;
                    m_parameterbounds.first.col(1) = inversMaps.col(bound);
                }
                else
                {
                    m_parameterbounds.first.col(1) = parIm;
                    m_parameterbounds.first.col(0) = inversMaps.col(bound);
                }
            }
        }
    }
    else // one side of patch1 is completely contained in patch2
    {

        if(m_side1.index() == 3 || m_side1.index() == 4) //sort w.r.t u
        {
            m_parameterbounds.first.row(0) = m_g1.parameterRange().row(0);
            m_parameterbounds.first(1, 0) = inversMaps.col(0)(1);
            m_parameterbounds.first(1, 1) = inversMaps.col(0)(1);
        }
        else
        {
            if(m_side1.index() == 1 || m_side1.index() == 2) // sort w.r.t v
            {
                m_parameterbounds.first(0,0) = inversMaps.col(0)(0);
                m_parameterbounds.first(0,1) = inversMaps.col(0)(0);
                m_parameterbounds.first.row(1) = m_g1.parameterRange().row(1);
            }
        }
    }

    // do the same for patch2 as before for patch1
    if(completeOnPatch1 == false && corners(1) != 0)
    {
        m_g2.invertPoints(m_g2.coefAtCorner(corners(1)).transpose(), parIm);

        index_t bound = 1;
        if(parIm.isApprox(inversMaps.col(1), 1e-6))
            bound = 3;

        if(m_side2.index() == 3 || m_side2.index() == 4) //sort w.r.t u
        {
            if(parIm(0,0) < inversMaps.col(bound)(0,0))
            {
                m_parameterbounds.second.col(0) = parIm;
                m_parameterbounds.second.col(1) = inversMaps.col(bound);
            }
            else
            {
                m_parameterbounds.second.col(1) = parIm;
                m_parameterbounds.second.col(0) = inversMaps.col(bound);
            }
        }
        else
        {
            if(m_side2.index() == 1 || m_side2.index() == 2) // sort w.r.t v
            {
                if(parIm(1,0) < inversMaps.col(1)(1,0))
                {
                    m_parameterbounds.second.col(0) = parIm;
                    m_parameterbounds.second.col(1) = inversMaps.col(bound);
                }
                else
                {
                    m_parameterbounds.second.col(1) = parIm;
                    m_parameterbounds.second.col(0) = inversMaps.col(bound);
                }
            }
        }
    }
    else
    {
        if(m_side2.index() == 3 || m_side2.index() == 4) //sort w.r.t u
        {
            m_parameterbounds.second.row(0) = m_g2.parameterRange().row(0);
            m_parameterbounds.second(1, 0) = inversMaps.col(1)(1);
            m_parameterbounds.second(1, 1) = inversMaps.col(1)(1);
        }
        else
        {
            if(m_side2.index() == 1 || m_side2.index() == 2) // sort w.r.t v
            {
                m_parameterbounds.second(0,0) = inversMaps.col(1)(0);
                m_parameterbounds.second(0,1) = inversMaps.col(1)(0);
                m_parameterbounds.second.row(1) = m_g2.parameterRange().row(1);
            }
        }
    }

    //gsInfo << "Init parameterbounds: \n" << m_parameterbounds.first << " \n and \n " << m_parameterbounds.second << "\n";
}


template<class T>
void gsRemapInterface<T>::findInterface(const boundaryInterface& bi)
{
    //const gsGeometry<T> & patch1 = m_domain.patch(0); // m_g1
    //const gsGeometry<T> & patch2 = m_domain.patch(1); // m_g2

    m_side1 = bi.first();
    m_side2 = bi.second();

    //gsInfo << "first: " << m_side1.index() << " and second " << m_side2.index() << "\n";
    //gsInfo << "the corresponding patches " << m_g1.id() << " and " << m_g2.id() << "\n";

    gsMultiPatch<T> geo1(m_g1);
    geo1.computeTopology();
    gsMultiPatch<T> geo2(m_g2);
    geo2.computeTopology();

    std::vector<patchSide> sidesPatch1 = geo1.boundaries();
    std::vector<patchSide> sidesPatch2 = geo2.boundaries();

    // vectors to store the rows indices of the coefficient matrices
    gsVector<index_t> boundaryPatch1, boundaryPatch2;

    // first find the sides of the patches which belong to the interface
    const index_t nCorners = 1 << m_g1.geoDim();
    gsMatrix<T> inversMaps = gsMatrix<T>::Zero(m_g1.geoDim(), 4); // matrix to store the preimages of the corresponding verices on the interface
    gsVector<index_t> corners(2); // vector to store the index of the corner which lies on the other patch, 0-th entry is the index of the corner for the first patch and vice versa
    corners.setZero();

    bool completeOnPatch2 = false, completeOnPatch1 = false; // check if one side of the patches is completely contained in the other side

    // two columns for the lower and the upper bound
    m_parameterbounds.first.resize(m_g1.geoDim(), 2);
    m_parameterbounds.second.resize(m_g1.geoDim(), 2);

    // matrix to store the coefficients of the corner values
    gsMatrix<T> c = gsMatrix<T>::Zero(1, m_g1.geoDim());
    gsMatrix<T> preIm;

    // find the side for the second patch, and the corner(s) of the first patch which lies on the interface
    gsVector<bool> onGeo;

    // this is the case if there are no matching corners
    for (index_t i = 1; i <= nCorners; i++) {
        c = m_g1.coefAtCorner(i).transpose();

        // Be aware of two matching vertices-> due to rounding errors on of the vertices can be on two sides maybe!
        gsMatrix<T> parLoc;

        std::vector<boxSide> side = m_g2.locateOn(c, onGeo, parLoc, false); // gives an empty side

        //for(size_t i = 0; i < side.size(); i++)
        //    gsInfo << "boundary patch: " << side[i] << "\n";

        if (onGeo(0) == true) {
            //inversMaps.col(1) = parLoc;

            if (corners(0) == 0) {
                inversMaps.col(1) = parLoc;
                corners(0) = i;
                //gsInfo << "corner: " << i << " and parametric loc.: \n" << parLoc << "\n";
            } else {
                // in case 2 corners lie on the interface of patch 2
                inversMaps.col(3) = parLoc;
                completeOnPatch2 = true; // side of patch1 is completely contained in the corresponding side of patch2 or the patches are fully matching
                //gsInfo << "corner: " << i << " and parametric loc.: \n" << parLoc << "\n";
            }
        }

    }

    // do the same for the first patch, i.e., finding the side for patch1 and the corner(s) for patch2
    for (index_t i = 1; i <= nCorners; i++) {
        c = m_g2.coefAtCorner(i).transpose();

        gsMatrix<T> parLoc;
        std::vector<boxSide> side = m_g1.locateOn(c, onGeo, parLoc, false);
        //gsInfo << " onGeo: \n" << onGeo(0) << "\n";
        //for(size_t i = 0; i < side.size(); i++)
        //    gsInfo << " side: \n" << side[i] << "\n";

        if (onGeo(0) == 1) {

            if (corners(1) == 0) {
                inversMaps.col(0) = parLoc;
                corners(1) = i;
                //gsInfo << "corner: " << i << " and parametric loc.: \n" << parLoc << "\n";
            } else {   // in case 2 corners lie on the interface of patch 1
                inversMaps.col(2) = parLoc;
                completeOnPatch1 = true;
            }
        }

    }



    // store the parametric bounds of the overlap for each patch, order the points w.r.t. the non fixed boundary side
    // so far the code is not very nice -> room for improvement!!!
    // Maybe considering the exact values instead of approximations would be better
    gsMatrix<T> parIm;
    if (completeOnPatch2 == false
                   && corners(0) != 0) // if the interface does not overlap entirely or if the side of patch2 is not a proper subset of the corresponding side of patch1
    {
        m_g1.invertPoints(m_g1.coefAtCorner(corners(0)).transpose(), parIm);

        index_t bound = 0;
        if (parIm.isApprox(inversMaps.col(0), 1e-6))
            bound = 2;

        //gsInfo << "bound: \n" << inversMaps.col(bound) << " and parameter Image: \n" << parIm << "\n";
        if (m_side1.index() == 3 || m_side1.index() == 4) //sort w.r.t u
        {
            if (parIm(0, 0) < inversMaps.col(bound)(0, 0)) {
                m_parameterbounds.first.col(0) = parIm;
                m_parameterbounds.first.col(1) = inversMaps.col(bound);
            } else {
                m_parameterbounds.first.col(1) = parIm;
                m_parameterbounds.first.col(0) = inversMaps.col(bound);
            }
        } else {
            if (m_side1.index() == 1 || m_side1.index() == 2) // sort w.r.t v
            {
                if (parIm(1, 0) < inversMaps.col(bound)(1, 0)) {
                    m_parameterbounds.first.col(0) = parIm;
                    m_parameterbounds.first.col(1) = inversMaps.col(bound);
                } else {
                    m_parameterbounds.first.col(1) = parIm;
                    m_parameterbounds.first.col(0) = inversMaps.col(bound);
                }
            }
        }
    } else // one side of patch1 is completely contained in patch2
    {

        if(completeOnPatch2 == false && corners(0) == 0)
        {
            if (m_side1.index() == 3 || m_side1.index() == 4) //sort w.r.t u
            {
                //m_parameterbounds.first.row(0) = m_g1.parameterRange().row(0);
                m_parameterbounds.first(0, 0) = inversMaps.col(0)(0);
                m_parameterbounds.first(0, 1) = inversMaps.col(2)(0);
                m_parameterbounds.first(1, 0) = inversMaps.col(0)(1);
                //m_parameterbounds.first(1, 1) = inversMaps.col(0)(1);
                m_parameterbounds.first(1, 1) = inversMaps.col(2)(1);
            } else {
                if (m_side1.index() == 1 || m_side1.index() == 2) // sort w.r.t v
                {
                    m_parameterbounds.first(0, 0) = inversMaps.col(0)(0);
                    m_parameterbounds.first(0, 1) = inversMaps.col(0)(0);
                    m_parameterbounds.first(1, 0) = inversMaps.col(0)(1) < inversMaps.col(2)(1) ? inversMaps.col(0)(1) : inversMaps.col(2)(1);
                    m_parameterbounds.first(1, 1) = inversMaps.col(0)(1) < inversMaps.col(2)(1) ? inversMaps.col(2)(1) : inversMaps.col(0)(1);
                }
            }
        } else
        {
            if (m_side1.index() == 3 || m_side1.index() == 4) //sort w.r.t u
            {
                m_parameterbounds.first.row(0) = m_g1.parameterRange().row(0);
                m_parameterbounds.first(1, 0) = m_side1.index() == 3 ? m_g1.parameterRange()(1,0) : m_g1.parameterRange()(1,1); //inversMaps.col(0)(1);
                m_parameterbounds.first(1, 1) = m_parameterbounds.first(1,0); //inversMaps.col(0)(1);
            } else {
                if (m_side1.index() == 1 || m_side1.index() == 2) // sort w.r.t v
                {
                    m_parameterbounds.first(0, 0) = m_side1.index() == 1 ? m_g1.parameterRange()(0,0) : m_g1.parameterRange()(0,1); //inversMaps.col(0)(0);
                    m_parameterbounds.first(0, 1) = m_parameterbounds.first(0, 0);//inversMaps.col(0)(0);
                    m_parameterbounds.first.row(1) = m_g1.parameterRange().row(1);
                }
            }
        }
    }

    // do the same for patch2 as before for patch1
    if (completeOnPatch1 == false && corners(1) != 0)
    {
        m_g2.invertPoints(m_g2.coefAtCorner(corners(1)).transpose(), parIm);

        index_t bound = 1;
        if (parIm.isApprox(inversMaps.col(1), 1e-6))
            bound = 3;

        if (m_side2.index() == 3 || m_side2.index() == 4) //sort w.r.t u
        {
            if (parIm(0, 0) < inversMaps.col(bound)(0, 0)) {
                m_parameterbounds.second.col(0) = parIm;
                m_parameterbounds.second.col(1) = inversMaps.col(bound);
            } else {
                m_parameterbounds.second.col(1) = parIm;
                m_parameterbounds.second.col(0) = inversMaps.col(bound);
            }
        } else {
            if (m_side2.index() == 1 || m_side2.index() == 2) // sort w.r.t v
            {
                if (parIm(1, 0) < inversMaps.col(bound)(1, 0)) {
                    m_parameterbounds.second.col(0) = parIm;
                    m_parameterbounds.second.col(1) = inversMaps.col(bound);
                } else
                {

                    m_parameterbounds.second.col(1) = parIm;
                    m_parameterbounds.second.col(0) = inversMaps.col(bound);
                }
            }
        }
    } else
    {
        if(completeOnPatch1 == false && corners(1) == 0)
        {
            if (m_side2.index() == 3 || m_side2.index() == 4) //sort w.r.t u
            {
                m_parameterbounds.second(0, 0) = inversMaps.col(1)(0);
                m_parameterbounds.second(0, 1) = inversMaps.col(3)(0);
                //m_parameterbounds.second.row(0) = m_g2.parameterRange().row(0);
                m_parameterbounds.second(1, 0) = inversMaps.col(1)(1);
                m_parameterbounds.second(1, 1) = inversMaps.col(3)(1);
            } else {
                if (m_side2.index() == 1 || m_side2.index() == 2) // sort w.r.t v
                {
                    m_parameterbounds.second(0, 0) = inversMaps.col(1)(0);
                    m_parameterbounds.second(0, 1) = inversMaps.col(1)(0);
                    m_parameterbounds.second(1, 0) = inversMaps.col(1)(1) < inversMaps.col(3)(1) ? inversMaps.col(1)(1) : inversMaps.col(3)(1);
                    m_parameterbounds.second(1, 1) = inversMaps.col(1)(1) < inversMaps.col(3)(1) ? inversMaps.col(3)(1) : inversMaps.col(1)(1);
                }
            }
        }
        else
        {
            if (m_side2.index() == 3 || m_side2.index() == 4) //sort w.r.t u
            {
                m_parameterbounds.second.row(0) = m_g2.parameterRange().row(0);
                m_parameterbounds.second(1, 0) = m_side2.index() == 3 ? m_g2.parameterRange()(1, 0) : m_g2.parameterRange()(1, 1);//inversMaps.col(1)(1);
                m_parameterbounds.second(1, 1) = m_parameterbounds.second(1, 0); //inversMaps.col(1)(1);
            } else {
                if (m_side2.index() == 1 || m_side2.index() == 2) // sort w.r.t v
                {
                    m_parameterbounds.second(0, 0) = m_side2.index() == 1 ? m_g2.parameterRange()(0, 0) : m_g2.parameterRange()(0, 1);//inversMaps.col(1)(0);
                    m_parameterbounds.second(0, 1) = m_parameterbounds.second(0, 0); //inversMaps.col(1)(0);
                    m_parameterbounds.second.row(1) = m_g2.parameterRange().row(1);
                }
            }
        }
    }

    //gsInfo << "Init parameterbounds: \n" << m_parameterbounds.first << " \n and \n " << m_parameterbounds.second << "\n";
}

template<class T>
bool gsRemapInterface<T>::checkIfMatching()
{
    short_t sameCorners = 0;

    for(short_t i = 1; i <= 1<<domainDim(); i++)
    {
        gsMatrix<T> c1 = m_g1.coefAtCorner(i).transpose();

        for(short_t j = 1; j <= 1<<domainDim(); j++)
        {
            gsMatrix<T> c2 = m_g2.coefAtCorner(j).transpose();
            if((c1-c2).squaredNorm() < 1.e-6)
                sameCorners++;
        }
    }

    if(sameCorners == 1 << (domainDim()-1))
        return true;
    else
        return false;
}

template <class T>
gsMatrix<T> gsRemapInterface<T>::checkIfInBound(const gsMatrix<T> & u) const
{
    // Here u contains only the coordinates in one direction
    gsMatrix<T> evalpts = u;

    T begin = m_parameterbounds.first(!m_side1.direction(), 0);
    T end = m_parameterbounds.first(!m_side1.direction(), 1);

    for(index_t c = 0; c < u.cols(); c++)
        if(u(0,c) - begin < 0)
            evalpts(0,c) += (begin - u(0,c));
        else
            break;


    for(index_t c = u.cols()-1; c > -1; c--)
        if(u(0, c) - end > 0)
            evalpts(0, c) -= (u(0, c) - end);
        else
            break;

    /*
    if(u(0,0) - begin < 0)
        evalpts(0,0) += (begin - u(0,0));


    if(u(0, u.cols()-1) - end > 0)
        evalpts(0, u.cols()-1) -= (u(0, u.cols()-1) - end);
    */

    return evalpts;
}

template<typename T>
void gsRemapInterface<T>::changeDir(const boundaryInterface & bi)
{
    T tmp = 0;
    index_t row = -1;

    if(m_side1.index() == 3 || m_side1.index() == 4)
    {
        if(bi.dirOrientation()(0) == 0)
        {
            if(m_side2.index() == 1 || m_side2.index() == 2) //change v parameters
                row = 1;
            else
                if(m_side2.index() == 3 || m_side2.index() == 4) // change u parameters
                    row = 0;
        }
    }

    if(m_side1.index() == 1 || m_side1.index() == 2)
    {
        if(bi.dirOrientation()(1) == 0)
        {
            if(m_side2.index() == 1 || m_side2.index() == 2) //change v parameters
                row = 1;
            else
                if(m_side2.index() == 3 || m_side2.index() == 4) // change u parameters
                    row = 0;
        }
    }

    if(row == -1) // everything is fine
        return;
    else // change the directions
    {
        tmp = m_parameterbounds.second(row,0);
        m_parameterbounds.second(row,0) =  m_parameterbounds.second(row,1);
        m_parameterbounds.second(row,1) =  tmp;
        m_flipSide2 = true;
    }
}

} // End namespace gismo
