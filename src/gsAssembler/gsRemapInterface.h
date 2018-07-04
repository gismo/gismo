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
#include <gsAssembler/gsAssembler.h>
#include <gsNurbs/gsTensorBSplineBasis.h>
#include <gsCore/gsAffineFunction.h>
#include <gsUtils/gsSortedVector.h>
#include <gsAssembler/gsQuadRule.h>
#include <gsNurbs/gsBSpline.h>



namespace gismo {


template <class T>
class gsRemapInterface : public gsFunction<T>
{
public:
    /// Shared pointer for gsAffineFunction
    typedef memory::shared_ptr< gsRemapInterface > Ptr;

    /// Unique pointer for gsAffineFunction
    typedef memory::unique_ptr< gsRemapInterface > uPtr;

    // Default constructor
    // TODO: include basis
    gsRemapInterface(): m_g1(NULL), m_g2(NULL), m_b1(NULL), m_b2(NULL) {}

    // Constructor which takes a multipatch and a boundary interface, useful if the interface is fully matching
    gsRemapInterface(const gsMultiPatch<T> & mp, const gsMultiBasis<T> & basis, const boundaryInterface & bi) : m_g1(mp[bi.first().patch]),
                                                                                                                m_g2(mp[bi.second().patch]),
                                                                                                                m_b1(basis[bi.first().patch]),
                                                                                                                m_b2(basis[bi.second().patch])// in most cases they are just the other way around
    {
        m_side1 = bi.first();
        m_side2 = bi.second();

        std::vector<boxCorner> corners;
        gsMatrix<T> inversCorners;

        m_parameterbounds.first.resize(m_g1.geoDim(), 2);
        m_parameterbounds.second.resize(m_g1.geoDim(), 2);


        switch (m_side1.index()) {
            case 1:
                m_parameterbounds.first.row(1) = mp[bi.first().patch].parameterRange().row(1);
                m_parameterbounds.first(0, 0) = mp[bi.first().patch].parameterRange()(0, 0);
                m_parameterbounds.first(0, 1) = mp[bi.first().patch].parameterRange()(0, 0);
                break;
            case 2:
                m_parameterbounds.first.row(1) = mp[bi.first().patch].parameterRange().row(1);
                m_parameterbounds.first(0, 0) = mp[bi.first().patch].parameterRange()(0, 1);
                m_parameterbounds.first(0, 1) = mp[bi.first().patch].parameterRange()(0, 1);
                break;
            case 3:
                m_parameterbounds.first.row(0) = mp[bi.first().patch].parameterRange().row(0);
                m_parameterbounds.first(1, 0) = mp[bi.first().patch].parameterRange()(1, 0);
                m_parameterbounds.first(1, 1) = mp[bi.first().patch].parameterRange()(1, 0);
                break;
            case 4:
                m_parameterbounds.first.row(0) = mp[bi.first().patch].parameterRange().row(0);
                m_parameterbounds.first(1, 0) = mp[bi.first().patch].parameterRange()(1, 1);
                m_parameterbounds.first(1, 1) = mp[bi.first().patch].parameterRange()(1, 1);
                break;
        }

        switch (m_side2.index()) {
            case 1:
                m_parameterbounds.second.row(1) = mp[bi.second().patch].parameterRange().row(1);
                m_parameterbounds.second(0, 0) = mp[bi.second().patch].parameterRange()(0, 0);
                m_parameterbounds.second(0, 1) = mp[bi.second().patch].parameterRange()(0, 0);
                break;
            case 2:
                m_parameterbounds.second.row(1) = mp[bi.second().patch].parameterRange().row(1);
                m_parameterbounds.second(0, 0) = mp[bi.second().patch].parameterRange()(0, 1);
                m_parameterbounds.second(0, 1) = mp[bi.second().patch].parameterRange()(0, 1);
                break;
            case 3:
                m_parameterbounds.second.row(0) = mp[bi.second().patch].parameterRange().row(0);
                m_parameterbounds.second(1, 0) = mp[bi.second().patch].parameterRange()(1, 0);
                m_parameterbounds.second(1, 1) = mp[bi.second().patch].parameterRange()(1, 0);
                break;
            case 4:
                m_parameterbounds.second.row(0) = mp[bi.second().patch].parameterRange().row(0);
                m_parameterbounds.second(1, 0) = mp[bi.second().patch].parameterRange()(1, 1);
                m_parameterbounds.second(1, 1) = mp[bi.second().patch].parameterRange()(1, 1);
                break;
        }

    }

    // Constructor for the class which takes two geometries
    gsRemapInterface(const gsGeometry<T> & g1, const gsGeometry<T> & g2, const gsBasis<T>* basis1, const gsBasis<T>* basis2) : m_g1(g1), m_g2(g2), m_b1(*basis1), m_b2(*basis2)
    {
        GISMO_ASSERT(m_g1.geoDim() == m_g2.geoDim(), "The two given geometries do not have the same dimension!");

        findInterface();

        // if the two patches match, maybe for testing purposes
        //m_domain.computeTopology();
        //gsInfo << "detail: " << m_domain.detail() << "\n";
    }

    // Destructor for the class
    ~gsRemapInterface() {};

    // Member for constructing the reparametrization
    // Fills in m_reparamInterfaceMap
    void constructReparam();

    // Memeber function for constructing the breakpoints
    // Fills in m_breakpoints
    void constructBreaks();

    // Helper to compute the closest point to lti on the other patch via Newton's method
    gsMatrix<T> closestPoint(const gsMatrix<T> b_null, const gsGeometry<T> & R, const gsMatrix<T> & lti);

    // rename: getPointsOnInterface() --> check eval_into, then evaluate both g1 and g2 and check equality

    const std::vector<T> getPointsOnInterface() const
    {
        index_t fixedDir = m_side1.direction();
        std::vector<T> vec(m_breakpoints.cols());

        if(fixedDir == 1)
        {
            for(index_t i = 0; i < m_breakpoints.cols(); i++)
                vec[i] = m_breakpoints.row(0)[i];
        }
        else
        {
            for(index_t i = 0; i < m_breakpoints.cols(); i++)
                vec[i] = m_breakpoints.row(1)[i];

        }

        return vec;
    }
    const gsBSpline<T> & giveInterfaceMap() const { return m_fittedCurve; }
    //const boxSide & giveSide() const { return m_side; }

    // Member to compute parametric values on the boundary from patch1 to the corresponding boundary values on patch2
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    // Overloaded member to get the parameter dimension of the domain
    virtual int domainDim() const { return m_g1.geoDim(); }

    typedef memory::unique_ptr< gsDomainIterator<T> > domainIter;
    domainIter makeDomainIterator() const
    {
        //GISMO_ASSERT(m_breakpoints.cols() == 0, "reparamInterface probably not called!, call before making a domain iterator");

        // ! For now it is the geometric basis
        //const gsTensorBSplineBasis<2, T> * tb = dynamic_cast<const gsTensorBSplineBasis<2, T> * >(&m_g1.basis());
        //GISMO_ASSERT(nullptr!=tb, "Bad pointer");

        gsTensorDomainBoundaryIterator<T> * tdi = new gsTensorDomainBoundaryIterator<T> (m_b1, m_side1);

        std::vector<T> newBreaks = getPointsOnInterface();
        gsInfo << "newBreaks: \n";
        for(size_t i = 0; i < m_breakpoints.cols(); i++)
            gsInfo << newBreaks[i] << "\t";

        gsInfo << "\n";

        // the input must be the direction which is moving
        //tdi->setBreaks(newBreaks, m_side1.direction()); -> gives the fixed direction
        // workaround: only works for 2 dimensions
        if(domainDim() == 2)
        {
            if (m_side1.direction() == 1)
                tdi->setBreaks(newBreaks, 0);
            else //m_side1.direction() == 0
                tdi->setBreaks(newBreaks, 1);
        }

        return domainIter(tdi);
    }

private:
    // The geometries to consider
    const gsGeometry<T> & m_g1;
    const gsGeometry<T> & m_g2;

    // The basis to consider
    const gsBasis<T> & m_b1;
    const gsBasis<T> & m_b2;

    gsMatrix<T> m_breakpoints;
    gsBSpline<T> m_fittedCurve;
    //boundaryInterface m_interface; // maybe not used anymore

    // Member to store which boundary is the interface for patch 1 and patch 2, respectively
    patchSide m_side1;
    patchSide m_side2;

    // Member to store the parameter bounds for both patches
    // A single matrix has the structure [lower, upper]^T
    std::pair<gsMatrix<T>, gsMatrix<T> > m_parameterbounds;

    // Member to enrich a matrix of 1D points to a matrix of m_domain.geoDim() points
    void enrichToVector(const int boundarySide, const gsGeometry<T> & geo, const gsMatrix<T> & intervals, gsMatrix<T> & pts);

    // Member to find the interface between the two incoming patches
    void findInterface();



}; // End gsRemapInterface



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
            startPatch2 = m_parameterbounds.second.col(0);
        }


    // Compute interface knots in physical domain by evaluating left and right geometry maps at the knot values --------------------------------------
    int numelP1 = domIt1->numElements();
    int numelP2 = domIt2->numElements();
    gsMatrix <T> physicalKnotsP1(m_g1.geoDim(), numelP1 + 1), physicalKnotsP2(m_g2.geoDim(), numelP2 + 1), dummy;

    domIt1->reset();
    domIt2->reset();
    int numBreaksPatch1 = 1, numBreaksPatch2 = 1; // vars to count the entries in the physical breakpoints

    // evaluate the first point of the interface
    m_g1.eval_into(startPatch1, dummy);
    physicalKnotsP1.col(0) = dummy;

    // loop over all elements of the boundary with interface part, but evaluate only element corners on the real interface
    for (; domIt1->good(); domIt1->next())
    {
        if(m_side1.index() == 3 || m_side1.index() == 4)
        {
            if(domIt1->lowerCorner()(0,0) > startPatch1(0,0) && domIt1->lowerCorner()(0,0) < m_parameterbounds.first(0,1))
            {
                m_g1.eval_into(domIt1->lowerCorner(), dummy);
                physicalKnotsP1.col(numBreaksPatch1) = dummy;
                numBreaksPatch1++;
            }
        }
        else
        {
            if(domIt1->lowerCorner()(1,0) > startPatch1(1,0) && domIt1->lowerCorner()(1,0) < m_parameterbounds.first(1,1))
            {
                m_g1.eval_into(domIt1->lowerCorner(), dummy);
                physicalKnotsP1.col(numBreaksPatch1) = dummy;
                numBreaksPatch1++;
            }
        }
        //domIt1->next();
    }

    // evaluate the last point of the interface
    m_g1.eval_into(domIt1->upperCorner(), dummy);
    physicalKnotsP1.col(numBreaksPatch1) = dummy;
    numBreaksPatch1++;
    gsInfo << "physical knots 1: \n" << physicalKnotsP1 << "\n";

    // do the same for patch 2 as above
    m_g2.eval_into(startPatch2, dummy);
    physicalKnotsP2.col(0) = dummy;

    for (; domIt2->good(); domIt2->next()) // for (int i = 0; i < numelP2; i++)
    {
        if(m_side2.index() == 3 || m_side2.index() == 4)
        {
            if (domIt2->lowerCorner()(0,0) > startPatch2(0,0) && domIt2->lowerCorner()(0,0) < m_parameterbounds.second(0,1))
            {
                m_g2.eval_into(domIt2->lowerCorner(), dummy);
                physicalKnotsP2.col(numBreaksPatch2) = dummy;
                numBreaksPatch2++;
            }
        }
        else
        {
            if (domIt2->lowerCorner()(1,0) > startPatch2(1,0) && domIt2->lowerCorner()(1,0) < m_parameterbounds.second(1,1))
            {
                m_g2.eval_into(domIt2->lowerCorner(), dummy);
                physicalKnotsP2.col(numBreaksPatch2) = dummy;
                numBreaksPatch2++;
            }
        }
        //domIt2->next();
    }

    m_g2.eval_into(domIt2->upperCorner(), dummy);
    physicalKnotsP2.col(numBreaksPatch2) = dummy;
    numBreaksPatch2++; // to get the number of entries
    gsInfo << "physical knots 2: \n" << physicalKnotsP2 << "\n";

    // store all the physical points in one vector
    gsMatrix<T> physicalBreaks(m_g1.geoDim(), numBreaksPatch1+numBreaksPatch2); // Assume m_g1.geoDim() == m_g2.geoDim()

    for(int c = 0; c < numBreaksPatch1; c++)
        physicalBreaks.col(c) = physicalKnotsP1.col(c);

    for(int c = 0; c < numBreaksPatch2; c++)
        physicalBreaks.col(numBreaksPatch1+c) = physicalKnotsP2.col(c);

    // compute the corresponding parameter values in one patch, here of patch2
    gsSortedVector<T> parameterBreaks;

    // Determine fixed coordinate of patch2 -> Use here patch2 because we compute the Interfacemap of patch1!!!
    // fixedDir ==  0 corresponds to fixed u and 1 corresponds to a fixed v
    index_t fixedDir = patchSide1.direction();

    gsMatrix<T> G2_parametric_LC;
    for (int i = 0;  i < physicalBreaks.cols(); i++)
    {
        // computes the preimages of the breakpoints for each of the two patches
        m_g1.invertPoints(physicalBreaks.col(i),G2_parametric_LC); // not exact, we have rounding errors
        // taking care of the rounding errors by iterating over the vector and checking the absolute value between the current
        // preimage and the already available ones
        if(fixedDir)
        {
            if(parameterBreaks.size() == 0)
                parameterBreaks.push_sorted_unique(G2_parametric_LC(0,0)); // sort w.r.t. u direction
            else
            {
                int j = 0;
                gsVector<bool> roundingError = gsVector<bool>::Constant(parameterBreaks.size(), true);

                for(typename gsSortedVector<T>::iterator it = parameterBreaks.begin(); it != parameterBreaks.end(); it++)
                {

                    if (math::abs(G2_parametric_LC(0, 0) - *it) > 1.e-4) {
                        roundingError(j) = false;
                    }
                    j++;
                }
                if((roundingError.array() == false).all())
                    parameterBreaks.push_sorted_unique(G2_parametric_LC(0, 0)); // sort w.r.t. u direction


            }
        }
        else
        {
            if(parameterBreaks.size() == 0)
                parameterBreaks.push_sorted_unique(G2_parametric_LC(1,0)); // sort w.r.t. v direction
            else
            {
                int j = 0;
                gsVector<bool> roundingError = gsVector<bool>::Constant(parameterBreaks.size(), true);

                for (typename gsSortedVector<T>::iterator it = parameterBreaks.begin();
                     it != parameterBreaks.end(); it++)
                {
                    if (math::abs(G2_parametric_LC(1, 0) - *it) > 1.e-4)
                        roundingError(j) = false;

                    j++;
                }
                if((roundingError.array() == false).all())
                    parameterBreaks.push_sorted_unique(G2_parametric_LC(1, 0)); // sort w.r.t. v direction
            }
        }

    }

    m_breakpoints = gsMatrix<T> (m_g1.geoDim(), parameterBreaks.size() ); // Assume m_g1.geoDim() == m_g2.geoDim()
    for (int i = 0;  i < parameterBreaks.size(); i++)
    {
        if(fixedDir)
            m_breakpoints.col(i) << parameterBreaks[i], G2_parametric_LC(1,0);
        else
            m_breakpoints.col(i) << G2_parametric_LC(0,0), parameterBreaks[i];
    }

    // only for tests
    gsMatrix<T> result;
    this->eval_into(m_breakpoints, result);

    gsInfo << "Mapped: \n" << result << "\n";

}


template<class T>
void gsRemapInterface<T>::constructReparam()
{
    const int numIntervals = 11;
    const int numGeometries = 2;

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

    for(index_t np = 0; np < numGeometries; np++)
    {
        if(np == 0)
        {
            if (m_side1.index() == 3 || m_side1.index() == 4) // v is fixed
            {
                firstKnot = m_parameterbounds.first(0,0);
                lastKnot = m_parameterbounds.first(0,1);
            } else // u is fixed
            {
                firstKnot = m_parameterbounds.first(1,0);
                lastKnot = m_parameterbounds.first(1,1);
            }
        }
        else
        {
            if (m_side2.index() == 3 || m_side2.index() == 4) // v is fixed
            {
                firstKnot = m_parameterbounds.second(0,0);
                lastKnot = m_parameterbounds.second(0,1);
            } else // u is fixed
            {
                firstKnot = m_parameterbounds.second(1,0);
                lastKnot = m_parameterbounds.second(1,1);
            }
        }

        upper << lastKnot;
        lower << firstKnot;

        t_vals.row(np) = uniformPointGrid(lower, upper, numIntervals); // uniformly distributed samples between the overlapping part of the interface

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
    enrichToVector(m_side1.index(), m_g2, t_vals.row(1), vals2dPatch2);

    m_g1.eval_into(vals2dPatch1,samples_left);
    m_g2.eval_into(vals2dPatch2,samples_right);

    //gsInfo << "vals2dPatch1:\n" << GEO_L_ref.coefs() << "\n vals2dPatch2:\n" << GEO_R_ref.coefs() << std::endl;
    //gsInfo << "vals2dPatch1:\n" << vals2dPatch1 << "\n vals2dPatch2:\n" << vals2dPatch2 << std::endl;
    std::cout << "samples left:\n" << samples_left << "\n samples right:\n" << samples_right << std::endl;

    gsMatrix<> B(numIntervals, m_g1.geoDim());

    for (int i=0; i<t_vals.cols();i++)
    {
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
    std::cout << "B:\n" << B << std::endl;

    // check the error
    // assume that the right map is the identity
    gsMatrix<T> eval_orig, eval_fit, B2, id;

    gsKnotVector<T> KV(t_vals(0,0),t_vals(0, numIntervals-1),0,4);

    gsCurveFitting<T> fit(t_vals.row(0).transpose(),B,KV);

    fit.compute();
    m_fittedCurve = fit.curve();
    std::cout << "Hi, I'm the resulting curve: \n" << m_fittedCurve << std::endl;

    int errorInterval = 1000;
    gsMatrix<T> eval_points = gsMatrix<T>::Zero(numGeometries, errorInterval);

    for(index_t np = 0; np < numGeometries; np++)
    {
        gsVector<T> lower(1), upper(1);
        lower << t_vals(np, 0);
        upper << t_vals(np, numIntervals-1);
        eval_points.row(np) = uniformPointGrid(lower, upper, errorInterval); // to check the error
    }

    m_fittedCurve.eval_into(eval_points.row(0), eval_fit);
    //eval_fit(0,0) -= 0.0001; // do a nasty slight correction since the first entry is out of the domain of definition due to rounding errors

    // TODO: also here use already available information
    enrichToVector(m_side2.index(), m_g2, eval_points.row(1), id);

    m_g1.eval_into(eval_fit, eval_orig);
    m_g2.eval_into(id,B2);
    //gsInfo << "b2: \n" << id.transpose() << " and eval_orig: \n" << eval_fit.transpose() << "\n";

    double error = 0;

    for (int i = 0; i < eval_points.cols(); i++)
        error += (id.col(i) - eval_fit.col(i)).squaredNorm();

    error = std::sqrt(error);

    std::cout << "Error: " << error << std::endl;


}


template<class T>
gsMatrix<T> gsRemapInterface<T>::closestPoint(const gsMatrix<T> b_null, const gsGeometry<T> & intMap, const gsMatrix<T> & lti)
{
    double stopping_criterion = std::numeric_limits<T>::max();
    double eps = 0.0001;
    gsMatrix<T> b_neu = b_null;

    gsMatrix<T> evalGeomap, diffPatch1Patch2, jacPatch2, deriv2Geomap, b_alt;

    // iterate until stopping criterion is satisfied
    while (stopping_criterion > eps)
    {
        b_alt = b_neu;

        // compute value of function to minimize
        intMap.eval_into(b_alt,evalGeomap);
        intMap.jacobian_into(b_alt,jacPatch2);
        intMap.deriv2_into(b_alt, deriv2Geomap);

        diffPatch1Patch2 = (lti - evalGeomap);

        jacPatch2.transposeInPlace();

        //gsInfo << "derivative: " << jacPatch2 << "\n";

        // generate the objective g
        gsMatrix<T> g(m_g1.geoDim(), 1);
        g = jacPatch2 * diffPatch1Patch2;

        // compute the jacobian of g
        // be carefull, is very error-prone
        gsMatrix<T> g_deriv = gsMatrix<T>::Zero(m_g1.geoDim(), m_g1.geoDim());
        g_deriv(0,0) = deriv2Geomap(0,0) * diffPatch1Patch2(0,0) - jacPatch2(0,0) * jacPatch2(0,0) + deriv2Geomap(3,0) * diffPatch1Patch2(1,0) - jacPatch2(0,1)*jacPatch2(0,1);
        g_deriv(0,1) = deriv2Geomap(2,0) * diffPatch1Patch2(0,0) - jacPatch2(0,0) * jacPatch2(1,0) + deriv2Geomap(5,0) * diffPatch1Patch2(1,0) - jacPatch2(0,1)*jacPatch2(1,1);
        g_deriv(1,0) = deriv2Geomap(2,0) * diffPatch1Patch2(0,0) - jacPatch2(0,1) * jacPatch2(0,0) + deriv2Geomap(5,0) * diffPatch1Patch2(1,0) - jacPatch2(1,1)*jacPatch2(0,1);
        g_deriv(1,1) = deriv2Geomap(1,0) * diffPatch1Patch2(0,0) - jacPatch2(0,1) * jacPatch2(0,1) + deriv2Geomap(4,0) * diffPatch1Patch2(1,0) - jacPatch2(1,1)*jacPatch2(1,1);

        // compute the next Newton step
        b_neu = b_alt + g_deriv.fullPivHouseholderQr().solve(-g);

        stopping_criterion = (b_alt - b_neu).squaredNorm();
    }

    gsMatrix<T> par_b_neu;
    intMap.invertPoints(b_neu, par_b_neu); // invert the obtained point to get the best fitting point in the parameter domain, we are looking for the argument -> see Paper
    return par_b_neu; //b_neu

}

template <typename T>
void gsRemapInterface<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT(u.rows() == m_g1.geoDim(), "The rows of the evaluation points must be equal to the dimension of the domain!");

    index_t fixedDir = m_side1.direction();
    gsMatrix<T> evalPoints;

    if(fixedDir == 1) // v is fixed => loop over u values
        evalPoints = u.row(0);
    else
        evalPoints = u.row(1);

    m_fittedCurve.eval_into(evalPoints, result);
}

// Function to enhance a sequence of 1D points in an interval to 2D points in the parameter domain
// A matrix pts with the already correct dimensions is expected
// boundarySide is the boundary side according to gismo's boundary::side
// intervals is the number of samples to enrich to more dimensions
// pts is the matrix populated with the corresponding points for more dimensions
// only works for 2d at the moment!!
// TODO: Generalize for arbitrary dimensions
template<class T>
void gsRemapInterface<T>::enrichToVector(const int boundarySide, const gsGeometry<T> & geo, const gsMatrix<T> & intervals, gsMatrix <T> &pts)
{
    pts.resize(geo.geoDim(), intervals.cols());

    const gsTensorBSplineBasis<2, T> *tb = dynamic_cast<const gsTensorBSplineBasis<2, T> * >(&geo.basis());

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

// Member to find the interface of the 2 input geometries
template<class T>
void gsRemapInterface<T>::findInterface()
{
    //const gsGeometry<T> & patch1 = m_domain.patch(0); // m_g1
    //const gsGeometry<T> & patch2 = m_domain.patch(1); // m_g2
    //m_side1.patch =

    gsMultiPatch<T> geo1(m_g1); geo1.computeTopology();
    gsMultiPatch<T> geo2(m_g2); geo2.computeTopology();

    std::vector<patchSide> sidesPatch1 = geo1.boundaries();
    std::vector<patchSide> sidesPatch2 = geo2.boundaries();

    //gsMultiPatch<T> domain;
    //domain.addPatch(m_g1); domain.addPatch(m_g2);

    // vectors to store the rows indices of the coefficient matrices
    gsVector<index_t > boundaryPatch1, boundaryPatch2;

    // first find the sides of the patches which belong to the interface
    const int nCorners = 1<<m_g1.geoDim();
    gsMatrix<T> inversMaps(m_g1.geoDim(), 2); // matrix to store the preimages of the corresponding verices on the interface
    gsVector<T> corners(2); // vector to store the index of the corner which lies on the other patch, 0-th entry is the index of the corner for the first patch and vice versa
    corners.setZero();

    bool completeOnPatch2 = false, completeOnPatch1 = false; // check if one side of the patches is completely contained in the other side

    // two columns for the lower and the upper bound
    m_parameterbounds.first.resize(m_g1.geoDim(), 2);
    m_parameterbounds.second.resize(m_g1.geoDim(), 2);

    // matrix to store the coefficients of the corner values
    gsMatrix<T> c = gsMatrix<T>::Zero(1, m_g1.geoDim());
    gsMatrix<T> preIm;

    // find the side for the second patch, and the corner(s) of the first patch which lies on the interface
    const gsTensorBSplineBasis<2, T> *tb = dynamic_cast<const gsTensorBSplineBasis<2, T> * >(&m_g2.basis());
    gsVector<bool> onGeo;

    for(int i = 1; i <= nCorners; i++)
    {
        c = m_g1.coefAtCorner(i).transpose();

        // Be aware of two matching vertices-> due to rounding errors on of the vertices can be on two sides maybe!
        gsMatrix<T> parLoc;
        std::vector<boxSide> side = m_g2.locateOn(c, onGeo, parLoc, true);

        //for(size_t i = 0; i < side.size(); i++)
        //    gsInfo << "boundary patch: " << side[i] << "\n";

        if(onGeo(0) == true)
        {
            inversMaps.col(1) = parLoc;

            if(!corners(0))
                corners(0) = i;
            else
                completeOnPatch2 = true; // side of patch1 is completely contained in the corresponding side of patch2 or the patches are fully matching


            m_side2.m_index = side[0].index();
            m_side2.patch = m_g2.id();
        }

    }

    // do the same for the first patch, i.e., finding the side for patch1 and the corner(s) for patch2
    tb = dynamic_cast<const gsTensorBSplineBasis<2, T> * >(&m_g1.basis());
    for(int i = 1; i <= nCorners; i++)
    {
        c = m_g2.coefAtCorner(i).transpose();

        //domain.locatePoints(c, 1, boundaryPatch2, preIm);

        /// new approach??
        gsMatrix<T> parLoc;
        std::vector<boxSide> side = m_g1.locateOn(c, onGeo, parLoc, true);
        //for(size_t i = 0; i < side.size(); i++)
        //    gsInfo << " side: \n" << side[i] << "\n";

        if(onGeo(0) == 1)
        {
            inversMaps.col(0) = parLoc;

            if(!corners(1))
                corners(1) = i;
            else
                completeOnPatch1 = true;

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
        if(m_side1.index() == 3 || m_side1.index() == 4) //sort w.r.t u
        {
            if(parIm(0,0) < inversMaps.col(0)(0,0))
            {
                m_parameterbounds.first.col(0) = parIm;
                m_parameterbounds.first.col(1) = inversMaps.col(0);
            }
            else
            {
                m_parameterbounds.first.col(1) = parIm;
                m_parameterbounds.first.col(0) = inversMaps.col(0);
            }
        }
        else
        {
            if(m_side1.index() == 1 || m_side1.index() == 2) // sort w.r.t v
            {
                if(parIm(1,0) < inversMaps.col(0)(1,0))
                {
                    m_parameterbounds.first.col(0) = parIm;
                    m_parameterbounds.first.col(1) = inversMaps.col(0);
                }
                else
                {
                    m_parameterbounds.first.col(1) = parIm;
                    m_parameterbounds.first.col(0) = inversMaps.col(0);
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
        if(m_side2.index() == 3 || m_side2.index() == 4) //sort w.r.t u
        {
            if(parIm(0,0) < inversMaps.col(1)(0,0))
            {
                m_parameterbounds.second.col(0) = parIm;
                m_parameterbounds.second.col(1) = inversMaps.col(1);
            }
            else
            {
                m_parameterbounds.second.col(1) = parIm;
                m_parameterbounds.second.col(0) = inversMaps.col(1);
            }
        }
        else
        {
            if(m_side2.index() == 1 || m_side2.index() == 2) // sort w.r.t v
            {
                if(parIm(1,0) < inversMaps.col(1)(1,0))
                {
                    m_parameterbounds.second.col(0) = parIm;
                    m_parameterbounds.second.col(1) = inversMaps.col(1);
                }
                else
                {
                    m_parameterbounds.second.col(1) = parIm;
                    m_parameterbounds.second.col(0) = inversMaps.col(1);
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

} // End namespace::gismo
