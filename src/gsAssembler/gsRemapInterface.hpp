/** @file gsRemapInterface.h

    @brief Provides a mapping between the corresponding sides of two patches sharing an interface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Seiler, R. Schneckenleitner, S. Takacs
    Created on: 2018-06-12
*/

#pragma once

#include <gsCore/gsMultiPatch.h>
#include <gsUtils/gsSortedVector.h>
#include <gsNurbs/gsTensorNurbsBasis.h>
#include <gsModeling/gsCurveFitting.h>

namespace gismo {

template <class T>
gsRemapInterface<T>::gsRemapInterface(const gsMultiPatch<T>   & mp,
                                      const gsMultiBasis<T>   & basis,
                                      const boundaryInterface & bi,
                                      index_t checkAffine)
    : m_g1(&(mp[bi.first().patch])), m_g2(&(mp[bi.second().patch])),
      m_b1(&(basis[bi.first().patch])), m_b2(&(basis[bi.second().patch])),
      m_side1(bi.first()), m_side2(bi.second()),
      m_isMatching(true), m_isAffine(true), m_flipSide2(false)
{
    GISMO_ASSERT( m_g1->geoDim()==m_g2->geoDim(), "gsRemapInterface: Dimensions do not agree." );

    // First we construct the affine mapping
    computeBoundingBox();

    // Setup the affine mapping
    m_intfMap = gsAffineFunction<T>::make(bi.dirMap(m_side1), bi.dirOrientation(m_side1), m_parameterBounds1, m_parameterBounds2);

    // Next, we check (if so desired by called) if the affine mapping coincides with the real mapping
    GISMO_ASSERT( checkAffine >= notAffine, "gsRemapInterface: Parameter checkAffine has invalid value." );
    if (checkAffine==alwaysAffine)
        m_isAffine = false;
    else if (checkAffine > 0)
        m_isAffine = checkIfAffine(checkAffine);

    if (m_isAffine)
    {
        // For computing the breaks, we need the reverse mapping as well
        // TODO: do not need this to be a member
        m_intfMap_inverse = gsAffineFunction<T>::make(bi.dirMap(m_side2), bi.dirOrientation(m_side2), m_parameterBounds2, m_parameterBounds1);
        constructBreaksAffine();
    }
    else
    {
        m_isMatching = false;
        GISMO_ENSURE(m_isAffine || domainDim() <= 2, "gsRemapInterface: Can handle non-matching interfaces only for 2 dimensions.");
        findInterface(bi);
        changeDir(bi);
        constructReparam();
        constructBreaksNotAffine();
    }

}

// Used for the affine case

template <class T>
gsMatrix<T> gsRemapInterface<T>::parameterBounds(const gsGeometry<T> & geo, boxSide s, index_t dim)
{
    gsMatrix<T> result(dim, 2);
    for (index_t i = 0; i<dim; ++i)
    {
        gsMatrix<T> pr = geo.parameterRange();
        if (s.direction()==i)
            result(i,0) = result(i,1) = pr( i, s.parameter() == false ? 0 : 1 );
        else
            result.row(i) = pr.row(i);
    }
    return result;
}


template <class T>
bool gsRemapInterface<T>::checkIfAffine(index_t steps)
{
    gsVector<T> lower = m_parameterBounds1.col(0);
    gsVector<T> upper = m_parameterBounds1.col(1);
    gsVector<unsigned> numberGridPoints = gsVector<unsigned>::Constant(domainDim(),2+steps);
    numberGridPoints[m_side1.direction()] = 1;
    gsMatrix<T> points = gsPointGrid(lower,upper,numberGridPoints);

    return  (
                m_g1->eval(points)
                -
                m_g2->eval(m_intfMap->eval(points))
            ).norm() < (T)(1.e-6);
}

template <class T>
void gsRemapInterface<T>::computeBoundingBox()
{
    const T tolerance = 0; // TODO: Consider interface only as non-matching if tolerance is exceeded
    // TODO: Simplify. It is enough to call newtonRaphson only twice (never cll m_g2->newtonRaphson)

    m_parameterBounds1 = parameterBounds( *m_g1, m_side1, domainDim() );
    m_parameterBounds2 = parameterBounds( *m_g2, m_side2, domainDim() );

    gsMatrix<T> phys2 = m_g2->eval(m_parameterBounds2);

    gsVector<T> parameterBoundsFrom2_lower = m_parameterBounds1.col(0);
    m_g1->newtonRaphson( phys2.col(0), parameterBoundsFrom2_lower, true, 10e-6, 100 );

    gsVector<T> parameterBoundsFrom2_upper = m_parameterBounds1.col(1);
    m_g1->newtonRaphson( phys2.col(1), parameterBoundsFrom2_upper, true, 10e-6, 100 );

    gsMatrix<T> phys1 = m_g1->eval(m_parameterBounds1);

    gsVector<T> parameterBoundsFrom1_lower = m_parameterBounds2.col(0);
    m_g2->newtonRaphson( phys1.col(0), parameterBoundsFrom1_lower, true, 10e-6, 100 );

    gsVector<T> parameterBoundsFrom1_upper = m_parameterBounds2.col(1);
    m_g2->newtonRaphson( phys1.col(1), parameterBoundsFrom1_upper, true, 10e-6, 100 );

    for (index_t i=0; i<domainDim(); ++i)
    {
        if (parameterBoundsFrom2_upper[i]<parameterBoundsFrom2_lower[i])
            std::swap( parameterBoundsFrom2_upper[i], parameterBoundsFrom2_lower[i] );
        if (parameterBoundsFrom2_lower[i]>m_parameterBounds1(i,0))
            { m_parameterBounds1(i,0) = parameterBoundsFrom2_lower[i]; m_isMatching = false; }
        if (parameterBoundsFrom2_upper[i]<m_parameterBounds1(i,1))
            { m_parameterBounds1(i,1) = parameterBoundsFrom2_upper[i]; m_isMatching = false; }

        if (parameterBoundsFrom1_upper[i]<parameterBoundsFrom1_lower[i])
            std::swap( parameterBoundsFrom1_upper[i], parameterBoundsFrom1_lower[i] );
        if (parameterBoundsFrom1_lower[i]>m_parameterBounds2(i,0))
            { m_parameterBounds2(i,0) = parameterBoundsFrom1_lower[i]; m_isMatching = false; }
        if (parameterBoundsFrom1_upper[i]<m_parameterBounds2(i,1))
            { m_parameterBounds2(i,1) = parameterBoundsFrom1_upper[i]; m_isMatching = false; }

    }

}

namespace {
template <class T, class Vector>
inline void addBreaks( std::vector< std::vector<T> >& breaks, const gsMatrix<T>& parameterBounds, const Vector& point )
{
    const T tolerance = 1.e-5;
    const index_t dim = point.rows();
    for (index_t d=0; d<dim; ++d)
    {
        const T t = point(d,0);
        if ( parameterBounds(d,0) <= t && t <= parameterBounds(d,1) )
        {
            // As in gsSortedVector::push_sorted_unique
            typename std::vector<T>::iterator pos = std::lower_bound(breaks[d].begin(), breaks[d].end(), t-tolerance );
            if ( pos == breaks[d].end() || *pos > t+tolerance ) // If not found
                breaks[d].insert(pos, t);
        }
    }
}
}

template <class T>
void gsRemapInterface<T>::constructBreaksAffine()
{
    m_breakpoints.resize(domainDim());

    const typename gsBasis<T>::domainIter domIt1 = m_b1->makeDomainIterator(m_side1);
    addBreaks(m_breakpoints, m_parameterBounds1, m_parameterBounds1.col(0));
    for (; domIt1->good(); domIt1->next())
        addBreaks(m_breakpoints, m_parameterBounds1, domIt1->upperCorner());
    addBreaks(m_breakpoints, m_parameterBounds1, m_parameterBounds1.col(1));

    const typename gsBasis<T>::domainIter domIt2 = m_b2->makeDomainIterator(m_side2);
    for (; domIt2->good(); domIt2->next())
        addBreaks(m_breakpoints, m_parameterBounds1, m_intfMap_inverse->eval(domIt2->upperCorner()));

}


// Used for the non-affine case


template <class T>
void gsRemapInterface<T>::constructBreaksNotAffine() {
    GISMO_ENSURE(domainDim()==2, "Not implemented for d!=2.");


    // computes break points per element

    const typename gsBasis<T>::domainIter domIt1 = m_b1->makeDomainIterator(static_cast<boxSide>(m_side1));
    const typename gsBasis<T>::domainIter domIt2 = m_b2->makeDomainIterator(static_cast<boxSide>(m_side2));

    const gsMatrix<T> startPatch1 = m_parameterBounds1.col(0);
    const gsMatrix<T> startPatch2 = m_parameterBounds2.col(m_flipSide2 ? 1 : 0);

    // Compute interface knots in physical domain by evaluating left and right geometry maps at the knot values
    const size_t numelP1 = domIt1->numElements();
    const size_t numelP2 = domIt2->numElements();
    gsMatrix <T> physicalKnotsP1(m_g1->geoDim(), numelP1 + 1), physicalKnotsP2(m_g2->geoDim(), numelP2 + 1), dummy;

    domIt1->reset();
    domIt2->reset();
    index_t numBreaksPatch1 = 1, numBreaksPatch2 = 1; // vars to count the entries in the physical breakpoints

    // evaluate the first point of the interface
    m_g1->eval_into(startPatch1, dummy);
    physicalKnotsP1.col(0) = dummy;

    // loop over all elements of the boundary with interface part, but evaluate only element corners on the real interface
    for (; domIt1->good(); domIt1->next())
    {
        if (m_side1.index() == 3 || m_side1.index() == 4) // v is fix
        {
            if (domIt1->lowerCorner()(0,0) > startPatch1(0,0) && domIt1->lowerCorner()(0,0) <= m_parameterBounds1(0,1))
            {
                // An ansatz for 3D??
                //if(domainDim() == 2 || (domainDim() == 3 &&
                //                        domIt1->lowerCorner()(2,0) > startPatch1(2,0) &&
                //                        domIt1->lowerCorner()(2,0) <= m_parameterBounds1(2,1)))
                //{
                m_g1->eval_into(domIt1->lowerCorner(), dummy);
                physicalKnotsP1.col(numBreaksPatch1) = dummy;
                numBreaksPatch1++;
                //}
            }
        }
        else
        {
            if (m_side1.index() == 1 || m_side1.index() == 2) // u is fix
            {
                if (domIt1->lowerCorner()(1,0) > startPatch1(1,0) && domIt1->lowerCorner()(1,0) <= m_parameterBounds1(1,1))
                {
                    m_g1->eval_into(domIt1->lowerCorner(), dummy);
                    physicalKnotsP1.col(numBreaksPatch1) = dummy;
                    numBreaksPatch1++;
                }
            }
            // for the 3D case??
            //else // w is fix
            //{
            //    if((domIt1->lowerCorner()(0,0) > startPatch1(0,0) && domIt1->lowerCorner()(0,0) <= m_parameterBounds1(0,1))
            //            && (domIt1->lowerCorner()(1,0) > startPatch1(1,0) && domIt1->lowerCorner()(1,0) <= m_parameterBounds1(1,1)))
            //    {
            //        m_g1->eval_into(domIt1->lowerCorner(), dummy);
            //        physicalKnotsP1.col(numBreaksPatch1) = dummy;
            //        numBreaksPatch1++;
            //    }
            //}
        }
        //domIt1->next();
    }

    // evaluate the last point of the interface, i.e., this last point must also be within the parameter bound
    if (m_side1.index() == 3 || m_side1.index() == 4)
    {
        if (domIt1->upperCorner()(0,0) <= m_parameterBounds1(0,1))
        {
            m_g1->eval_into(domIt1->upperCorner(), dummy);
            physicalKnotsP1.col(numBreaksPatch1) = dummy;
            numBreaksPatch1++;
        }
    }
    else
    {
        if (m_side1.index() == 1 || m_side1.index() == 2)
        {
            if (domIt1->upperCorner()(1,0) <= m_parameterBounds1(1,1))
            {
                m_g1->eval_into(domIt1->upperCorner(), dummy);
                physicalKnotsP1.col(numBreaksPatch1) = dummy;
                numBreaksPatch1++;
            }
        }
    }
    //gsInfo << "physical knots 1: \n" << physicalKnotsP1 << "\n";



    // do the same for patch 2 as above
    m_g2->eval_into(startPatch2, dummy);
    physicalKnotsP2.col(0) = dummy;

    for (; domIt2->good(); domIt2->next()) // for (index_t i = 0; i < numelP2; i++)
    {
        if (m_side2.index() == 3 || m_side2.index() == 4)
        {
            if (domIt2->lowerCorner()(0,0) > startPatch2(0,0) && domIt2->lowerCorner()(0,0) < std::max(m_parameterBounds2(0,1), m_parameterBounds2(0,0)))
            {
                m_g2->eval_into(domIt2->lowerCorner(), dummy);
                physicalKnotsP2.col(numBreaksPatch2) = dummy;
                numBreaksPatch2++;
            }
        }
        else
        {
            if(m_side2.index() == 1 || m_side2.index() == 2)
            {
                if (domIt2->lowerCorner()(1,0) > startPatch2(1,0) && domIt2->lowerCorner()(1,0) < std::max(m_parameterBounds2(1,1), m_parameterBounds2(1,0)))
                {
                    m_g2->eval_into(domIt2->lowerCorner(), dummy);
                    physicalKnotsP2.col(numBreaksPatch2) = dummy;
                    numBreaksPatch2++;
                }
            }
        }
        //domIt2->next();
    }

    // add only the breakpoints within the parameter bounds
    if (m_side2.index() == 3 || m_side2.index() == 4)
    {
        if (domIt2->upperCorner()(0,0) <= std::max(m_parameterBounds2(0,1), m_parameterBounds2(0,0)))
        {
            m_g2->eval_into(domIt2->upperCorner(), dummy);
            physicalKnotsP2.col(numBreaksPatch2) = dummy;
            numBreaksPatch2++;// to get the number of entries
        }
    }
    else
    {
        if (m_side2.index() == 1 || m_side2.index() == 2)
        {
            if(domIt2->upperCorner()(1,0) <= std::max(m_parameterBounds2(1,1), m_parameterBounds2(1,0)))
            {
                m_g2->eval_into(domIt2->upperCorner(), dummy);
                physicalKnotsP2.col(numBreaksPatch2) = dummy;
                numBreaksPatch2++;
            }
        }
    }
    //gsInfo << "physical knots 2: \n" << physicalKnotsP2 << "\n";

    // store all the physical points in one vector
    gsMatrix<T> physicalBreaks(domainDim(), numBreaksPatch1+numBreaksPatch2); // Assume m_g1->geoDim() == m_g2->geoDim()

    for (index_t c = 0; c < numBreaksPatch1; c++)
        physicalBreaks.col(c) = physicalKnotsP1.col(c);

    for (index_t c = 0; c < numBreaksPatch2; c++)
        physicalBreaks.col(numBreaksPatch1+c) = physicalKnotsP2.col(c);

    // compute the corresponding parameter values in one patch, here of patch1
        gsSortedVector<T> parameterBreaks;

        // Determine fixed coordinate of patch2 -> Use here patch2 because we compute the Interfacemap of patch1!!!
        // fixedDir ==  0 corresponds to fixed u and 1 corresponds to a fixed v
        index_t fixedDir = m_side1.direction();

        gsMatrix<T> G2_parametric_LC;
        for (index_t i = 0; i < physicalBreaks.cols(); i++) {
            // computes the preimages of the breakpoints for each of the two patches
            m_g1->invertPoints(physicalBreaks.col(i), G2_parametric_LC); // not exact, we have rounding errors
            // taking care of the rounding errors by iterating over the vector and checking the absolute value between the current
            // preimage and the already available ones
            if (fixedDir == 1)
            {
                if (parameterBreaks.size() == 0)
                    parameterBreaks.push_sorted_unique(G2_parametric_LC(0, 0)); // sort w.r.t. u direction
                else
                {
                    index_t j = 0;
                    gsVector<bool> roundingError = gsVector<bool>::Constant(parameterBreaks.size(), true);

                    for (typename gsSortedVector<T>::iterator it = parameterBreaks.begin();
                         it != parameterBreaks.end(); it++)
                    {
                        if (math::abs(G2_parametric_LC(0, 0) - *it) > 1.e-4) {
                            roundingError(j) = false;
                        }
                        j++;
                    }
                    if (( false == roundingError.array() ).all())
                        parameterBreaks.push_sorted_unique(G2_parametric_LC(0, 0)); // sort w.r.t. u direction


                }
            }
            else
            {
                if (parameterBreaks.size() == 0)
                    parameterBreaks.push_sorted_unique(G2_parametric_LC(1, 0)); // sort w.r.t. v direction
                else
                {
                    index_t j = 0;
                    gsVector<bool> roundingError = gsVector<bool>::Constant(parameterBreaks.size(), true);

                    for (typename gsSortedVector<T>::iterator it = parameterBreaks.begin();
                         it != parameterBreaks.end(); it++)
                    {
                        if (math::abs(G2_parametric_LC(1, 0) - *it) > 1.e-4)
                            roundingError(j) = false;

                        j++;
                    }
                    if (( false == roundingError.array() ).all())
                        parameterBreaks.push_sorted_unique(G2_parametric_LC(1, 0)); // sort w.r.t. v direction
                }
            }

        }

        m_breakpoints.resize(2);
        m_breakpoints[1-fixedDir] = parameterBreaks;
}


template <class T>
void gsRemapInterface<T>::constructReparam()
{
    const index_t numIntervals = 11; // ?
    const index_t numGeometries = 2;

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

        //gsInfo << "parameterbounds: \n" << m_parameterBounds2 << "\n";
        //gsInfo << "patch: \n" << m_g2->id() << "\n";
        for (index_t np = 0; np < numGeometries; np++) {
            if (np == 0)
            {
                if (m_side1.index() == 3 || m_side1.index() == 4) // v is fixed
                {
                    firstKnot = m_parameterBounds1(0, 0);
                    lastKnot = m_parameterBounds1(0, 1);
                }
                else // u is fixed
                {
                    firstKnot = m_parameterBounds1(1, 0);
                    lastKnot = m_parameterBounds1(1, 1);
                }
            } else {
                if (m_side2.index() == 3 || m_side2.index() == 4) // v is fixed
                {
                    firstKnot = m_parameterBounds2(0, 0);
                    lastKnot = m_parameterBounds2(0, 1);
                }
                else // u is fixed
                {
                    firstKnot = m_parameterBounds2(1, 0);
                    lastKnot = m_parameterBounds2(1, 1);
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
        enrichToVector(m_side1, *m_g1, t_vals.row(0), vals2dPatch1);
        enrichToVector(m_side2, *m_g2, t_vals.row(1), vals2dPatch2);

        m_g1->eval_into(vals2dPatch1, samples_left);
        m_g2->eval_into(vals2dPatch2, samples_right);

        //gsInfo << "vals2dPatch1:\n" << GEO_L_ref.coefs() << "\n vals2dPatch2:\n" << GEO_R_ref.coefs() << std::endl;
        //gsInfo << "vals2dPatch1:\n" << vals2dPatch1 << "\n vals2dPatch2:\n" << vals2dPatch2 << std::endl;
        //std::cout << "samples left:\n" << samples_left << "\n samples right:\n" << samples_right << std::endl;

        gsMatrix<T> B(numIntervals, m_g1->geoDim());

        for (index_t i = 0; i < t_vals.cols(); i++)
        {
            // find a suitable start value for the Newton iteration
            find_start_value = (samples_right.colwise()) - samples_left.col(i);

            size_t row, col;

            find_start_value.colwise().squaredNorm().minCoeff(&row, &col);

            gsVector<T> b_null = samples_right.col(col);

            // Pass on g2 if one wants to find a mapping from interface1 to interface2
            //gsMatrix<T> b = closestPoint(b_null, g2, samples_left.col(i));

            // this gives the same result as above
            m_g2->newtonRaphson(samples_left.col(i), b_null, true, 10e-6, 100);
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
        m_intfMap = fit.curve().clone();
        std::cout << "Hi, I'm the resulting curve: \n" << *m_intfMap << std::endl;

        unsigned errorInterval = 10;
        gsVector<unsigned > errorSamples(1);
        errorSamples << errorInterval;

        gsMatrix<T> eval_points;// = gsMatrix<T>::Zero(numGeometries, errorInterval);

        for (index_t np = 0; np < numGeometries; np++)
        {
            gsVector<T> lowerVal(1), upperVal(1);
            lowerVal << t_vals(np, 0);
            upperVal << t_vals(np, numIntervals - 1);
            //gsMatrix<T> grid = uniformPointGrid(lowerVal, upperVal, errorInterval);
            gsMatrix<T> grid = gsPointGrid(lowerVal, upperVal, errorSamples);
            eval_points.conservativeResize(np + 1, grid.cols()); // to check the error
            eval_points.row(np) = grid;
        }

        m_intfMap->eval_into(eval_points.row(0), eval_fit);
        //eval_fit(0,0) -= 0.0001; // do a nasty slight correction since the first entry is out of the domain of definition due to rounding errors

        // TODO: also here use already available information
        enrichToVector(m_side2, *m_g2, eval_points.row(1), id);

        m_g2->eval_into(eval_fit, eval_orig);
        m_g2->eval_into(id, B2);
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
                m_g2->eval_into(eval_fit, eval_orig);
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

template <class T>
void gsRemapInterface<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT(u.rows() == domainDim(), "gsRemapInterface<T>::eval_into: "
        "The rows of the evaluation points must be equal to the dimension of the domain.");

    if (m_isAffine)
    {
        m_intfMap->eval_into(u, result);
    }
    else
    {
        // TODO: What is going on here?

        const index_t fixedDir = m_side1.direction();
        // v is fixed => loop over u values

        m_intfMap->eval_into(checkIfInBound(u.row(!fixedDir)), result); // ????

        // need here the second basis since result store points in the second geometry
        if (const gsTensorBSplineBasis<2, T> *tb = dynamic_cast<const gsTensorBSplineBasis<2, T> * >(&m_g2->basis()))
        {
            const short_t direction = m_side2.direction();
            result.row(direction).setConstant(m_side2.parameter() ?
                                               tb->knots(direction).last() :
                                               tb->knots(direction).first());

            return;
        }

        if (const gsTensorNurbsBasis<2, T> * ntb = dynamic_cast<const gsTensorNurbsBasis<2, T> * >(&(m_g2->basis())))
        {
            const short_t direction = m_side2.direction();
            result.row(direction).setConstant(m_side2.parameter() ?
                                               ntb->source().knots(direction).last() :
                                               ntb->source().knots(direction).first());
            return;
        }

        GISMO_ERROR("Unfitted interface not supported");
    }

}

template <class T>
typename gsDomainIterator<T>::uPtr gsRemapInterface<T>::makeDomainIterator() const
{
    gsTensorDomainBoundaryIterator<T> * tdi = new gsTensorDomainBoundaryIterator<T> (*m_b1, m_side1);
    for (index_t i=0; i<domainDim(); ++i)
    {
        if (i!=m_side1.direction())
            tdi->setBreaks(m_breakpoints[i],i);
    }
    return typename gsDomainIterator<T>::uPtr(tdi);
}

// Function to enhance a sequence of 1D points in an interval to 2D points in the parameter domain
// A matrix pts with the already correct dimensions is expected
// boundarySide is the boundary side according to gismo's boundary::side
// intervals is the number of samples to enrich to more dimensions
// pts is the matrix populated with the corresponding points for more dimensions
// only works for 2d at the moment!!
// TODO: Generalize for arbitrary dimensions
template <class T>
void gsRemapInterface<T>::enrichToVector(const boxSide         boundarySide,
                                         const gsGeometry<T> & geo,
                                         const gsMatrix<T>   & intervals,
                                               gsMatrix <T>  & pts)
{
    pts.resize(geo.geoDim(), intervals.cols());

    //const gsTensorBSplineBasis<2, T> *tb = dynamic_cast<const gsTensorBSplineBasis<2, T> * >(&geo.basis());

    //if(tb == NULL)
    //const gsTensorNurbsBasis<2, T> * tb = dynamic_cast<const gsTensorNurbsBasis<2, T> * >(&(geo.basis()));

    const index_t dim = 2;

    if(const gsTensorBSplineBasis<2, T> *tb = dynamic_cast<const gsTensorBSplineBasis<2, T> * >(&geo.basis()))
    {
        for (index_t i=0; i<dim; ++i)
        {
            if (boundarySide.direction()==i)
                pts.row(i) = gsMatrix<T>::Constant(1, intervals.cols(),
                    boundarySide.parameter()==0 ? tb->knots(0).first() : tb->knots(i).last()
                );
            else
                pts.row(i) = intervals;
        }
    }
    else
    {
        const gsTensorNurbsBasis<2, T> * ntb = dynamic_cast<const gsTensorNurbsBasis<2, T> * >(&(geo.basis()));

        for (index_t i=0; i<dim; ++i)
        {
            if (boundarySide.direction()==i)
                pts.row(i) = gsMatrix<T>::Constant(1, intervals.cols(),
                    boundarySide.parameter()==0 ? ntb->source().knots(0).first() : ntb->source().knots(i).last()
                );
            else
                pts.row(i) = intervals;
        }
    }
}

template <class T>
void gsRemapInterface<T>::findInterface(const boundaryInterface& bi)
{

    GISMO_UNUSED(bi);

    // first find the sides of the patches which belong to the interface
    const index_t nCorners = 1 << m_g1->geoDim();
    gsMatrix<T> inverseMaps = gsMatrix<T>::Zero(m_g1->geoDim(), 4); // matrix to store the preimages of the corresponding verices on the interface
    gsVector<index_t> corners(2); // vector to store the index of the corner which lies on the other patch, 0-th entry is the index of the corner for the first patch and vice versa
    corners.setZero();

    bool completeOnPatch2 = false, completeOnPatch1 = false; // check if one side of the patches is completely contained in the other side

    // two columns for the lower and the upper bound
    m_parameterBounds1.resize(m_g1->geoDim(), 2);
    m_parameterBounds2.resize(m_g1->geoDim(), 2);

    // matrix to store the coefficients of the corner values
    gsMatrix<T> c = gsMatrix<T>::Zero(1, m_g1->geoDim());

    // find the side for the second patch, and the corner(s) of the first patch which lies on the interface
    gsVector<bool> onGeo;

    // this is the case if there are no matching corners
    for (index_t i = 1; i <= nCorners; i++)
    {
        c = m_g1->coefAtCorner(i).transpose();

        // Be aware of two matching vertices-> due to rounding errors on of the vertices can be on two sides maybe!
        gsMatrix<T> parLoc;

        std::vector<boxSide> side = m_g2->locateOn(c, onGeo, parLoc, false); // gives an empty side

        //for(size_t i = 0; i < side.size(); i++)
        //    gsInfo << "boundary patch: " << side[i] << "\n";

        if (onGeo(0) == true)
        {
            //inverseMaps.col(1) = parLoc;

            if (corners(0) == 0)
            {
                inverseMaps.col(1) = parLoc;
                corners(0) = i;
                //gsInfo << "corner: " << i << " and parametric loc.: \n" << parLoc << "\n";
            }
            else
            {
                // in case 2 corners lie on the interface of patch 2
                inverseMaps.col(3) = parLoc;
                completeOnPatch2 = true; // side of patch1 is completely contained in the corresponding side of patch2 or the patches are fully matching
                //gsInfo << "corner: " << i << " and parametric loc.: \n" << parLoc << "\n";
            }
        }

    }

    // do the same for the first patch, i.e., finding the side for patch1 and the corner(s) for patch2
    for (index_t i = 1; i <= nCorners; i++)
    {
        c = m_g2->coefAtCorner(i).transpose();

        gsMatrix<T> parLoc;
        std::vector<boxSide> side = m_g1->locateOn(c, onGeo, parLoc, false);
        //gsInfo << " onGeo: \n" << onGeo(0) << "\n";
        //for(size_t i = 0; i < side.size(); i++)
        //    gsInfo << " side: \n" << side[i] << "\n";

        if (onGeo(0) == 1)
        {
            if (corners(1) == 0)
            {
                inverseMaps.col(0) = parLoc;
                corners(1) = i;
                //gsInfo << "corner: " << i << " and parametric loc.: \n" << parLoc << "\n";
            } else {   // in case 2 corners lie on the interface of patch 1
                inverseMaps.col(2) = parLoc;
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
        m_g1->invertPoints(m_g1->coefAtCorner(corners(0)).transpose(), parIm);

        index_t bound = 0;
        if (parIm.isApprox(inverseMaps.col(0), 1e-6))
            bound = 2;

        //gsInfo << "bound: \n" << inverseMaps.col(bound) << " and parameter Image: \n" << parIm << "\n";
        if (m_side1.index() == 3 || m_side1.index() == 4) //sort w.r.t u
        {
            if (parIm(0, 0) < inverseMaps.col(bound)(0, 0))
            {
                m_parameterBounds1.col(0) = parIm;
                m_parameterBounds1.col(1) = inverseMaps.col(bound);
            } else {
                m_parameterBounds1.col(1) = parIm;
                m_parameterBounds1.col(0) = inverseMaps.col(bound);
            }
        }
        else if (m_side1.index() == 1 || m_side1.index() == 2) // sort w.r.t v
        {
            if (parIm(1, 0) < inverseMaps.col(bound)(1, 0))
            {
                m_parameterBounds1.col(0) = parIm;
                m_parameterBounds1.col(1) = inverseMaps.col(bound);
            } else {
                m_parameterBounds1.col(1) = parIm;
                m_parameterBounds1.col(0) = inverseMaps.col(bound);
            }
        }
    }
    else // one side of patch1 is completely contained in patch2
    {

        if(completeOnPatch2 == false && corners(0) == 0)
        {
            if (m_side1.index() == 3 || m_side1.index() == 4) //sort w.r.t u
            {
                //m_parameterBounds1.row(0) = m_g1->parameterRange().row(0);
                m_parameterBounds1(0, 0) = inverseMaps.col(0)(0);
                m_parameterBounds1(0, 1) = inverseMaps.col(2)(0);
                m_parameterBounds1(1, 0) = inverseMaps.col(0)(1);
                //m_parameterBounds1(1, 1) = inverseMaps.col(0)(1);
                m_parameterBounds1(1, 1) = inverseMaps.col(2)(1);
            }
            else if (m_side1.index() == 1 || m_side1.index() == 2) // sort w.r.t v
            {
                m_parameterBounds1(0, 0) = inverseMaps.col(0)(0);
                m_parameterBounds1(0, 1) = inverseMaps.col(0)(0);
                m_parameterBounds1(1, 0) = inverseMaps.col(0)(1) < inverseMaps.col(2)(1) ? inverseMaps.col(0)(1) : inverseMaps.col(2)(1);
                m_parameterBounds1(1, 1) = inverseMaps.col(0)(1) < inverseMaps.col(2)(1) ? inverseMaps.col(2)(1) : inverseMaps.col(0)(1);
            }
        }
        else
        {
            if (m_side1.index() == 3 || m_side1.index() == 4) //sort w.r.t u
            {
                m_parameterBounds1.row(0) = m_g1->parameterRange().row(0);
                m_parameterBounds1(1, 0) = m_side1.index() == 3 ? m_g1->parameterRange()(1,0) : m_g1->parameterRange()(1,1); //inverseMaps.col(0)(1);
                m_parameterBounds1(1, 1) = m_parameterBounds1(1,0); //inverseMaps.col(0)(1);
            }
            else if (m_side1.index() == 1 || m_side1.index() == 2) // sort w.r.t v
            {
                m_parameterBounds1(0, 0) = m_side1.index() == 1 ? m_g1->parameterRange()(0,0) : m_g1->parameterRange()(0,1); //inverseMaps.col(0)(0);
                m_parameterBounds1(0, 1) = m_parameterBounds1(0, 0);//inverseMaps.col(0)(0);
                m_parameterBounds1.row(1) = m_g1->parameterRange().row(1);
            }
        }
    }

    // do the same for patch2 as before for patch1
    if (completeOnPatch1 == false && corners(1) != 0)
    {
        m_g2->invertPoints(m_g2->coefAtCorner(corners(1)).transpose(), parIm);

        index_t bound = 1;
        if (parIm.isApprox(inverseMaps.col(1), 1e-6))
            bound = 3;

        if (m_side2.index() == 3 || m_side2.index() == 4) //sort w.r.t u
        {
            if (parIm(0, 0) < inverseMaps.col(bound)(0, 0))
            {
                m_parameterBounds2.col(0) = parIm;
                m_parameterBounds2.col(1) = inverseMaps.col(bound);
            } else {
                m_parameterBounds2.col(1) = parIm;
                m_parameterBounds2.col(0) = inverseMaps.col(bound);
            }
        }
        else if (m_side2.index() == 1 || m_side2.index() == 2) // sort w.r.t v
        {
            if (parIm(1, 0) < inverseMaps.col(bound)(1, 0))
            {
                m_parameterBounds2.col(0) = parIm;
                m_parameterBounds2.col(1) = inverseMaps.col(bound);
            }
            else
            {
                m_parameterBounds2.col(1) = parIm;
                m_parameterBounds2.col(0) = inverseMaps.col(bound);
            }
        }
    }
    else
    {
        if(completeOnPatch1 == false && corners(1) == 0)
        {
            if (m_side2.index() == 3 || m_side2.index() == 4) //sort w.r.t u
            {
                m_parameterBounds2(0, 0) = inverseMaps.col(1)(0);
                m_parameterBounds2(0, 1) = inverseMaps.col(3)(0);
                //m_parameterBounds2.row(0) = m_g2->parameterRange().row(0);
                m_parameterBounds2(1, 0) = inverseMaps.col(1)(1);
                m_parameterBounds2(1, 1) = inverseMaps.col(3)(1);
            }
            else if (m_side2.index() == 1 || m_side2.index() == 2) // sort w.r.t v
            {
                m_parameterBounds2(0, 0) = inverseMaps.col(1)(0);
                m_parameterBounds2(0, 1) = inverseMaps.col(1)(0);
                m_parameterBounds2(1, 0) = inverseMaps.col(1)(1) < inverseMaps.col(3)(1) ? inverseMaps.col(1)(1) : inverseMaps.col(3)(1);
                m_parameterBounds2(1, 1) = inverseMaps.col(1)(1) < inverseMaps.col(3)(1) ? inverseMaps.col(3)(1) : inverseMaps.col(1)(1);
            }
        }
        else
        {
            if (m_side2.index() == 3 || m_side2.index() == 4) //sort w.r.t u
            {
                m_parameterBounds2.row(0) = m_g2->parameterRange().row(0);
                m_parameterBounds2(1, 0) = m_side2.index() == 3 ? m_g2->parameterRange()(1, 0) : m_g2->parameterRange()(1, 1);//inverseMaps.col(1)(1);
                m_parameterBounds2(1, 1) = m_parameterBounds2(1, 0); //inverseMaps.col(1)(1);
            }
            else if (m_side2.index() == 1 || m_side2.index() == 2) // sort w.r.t v
            {
                m_parameterBounds2(0, 0) = m_side2.index() == 1 ? m_g2->parameterRange()(0, 0) : m_g2->parameterRange()(0, 1);//inverseMaps.col(1)(0);
                m_parameterBounds2(0, 1) = m_parameterBounds2(0, 0); //inverseMaps.col(1)(0);
                m_parameterBounds2.row(1) = m_g2->parameterRange().row(1);
            }
        }
    }

    //gsInfo << "Init parameterbounds: \n" << m_parameterBounds1 << " \n and \n " << m_parameterBounds2 << "\n";
}

template <class T>
gsMatrix<T> gsRemapInterface<T>::checkIfInBound(const gsMatrix<T> & u) const
{
    // Here u contains only the coordinates in one direction
    gsMatrix<T> evalpts = u;

    const T begin = m_parameterBounds1(!m_side1.direction(), 0);
    const T end = m_parameterBounds1(!m_side1.direction(), 1);

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

template <class T>
void gsRemapInterface<T>::changeDir(const boundaryInterface & bi)
{

    if(m_side1.direction() == 1)
    {
        if(bi.dirOrientation()(0) == 0)
        {
            if(m_side2.direction() == 0) //change v parameters
                { std::swap( m_parameterBounds2(1,0), m_parameterBounds2(1,1) ); m_flipSide2=true; }
            else if(m_side2.direction() == 1) // change u parameters
                { std::swap( m_parameterBounds2(0,0), m_parameterBounds2(0,1) ); m_flipSide2=true; }
        }
    }

    if(m_side1.direction() == 0)
    {
        if(bi.dirOrientation()(1) == 0)
        {
            if(m_side2.direction() == 0) //change v parameters
                { std::swap( m_parameterBounds2(1,0), m_parameterBounds2(1,1) ); m_flipSide2=true; }
            else if(m_side2.direction() == 1) // change u parameters
                { std::swap( m_parameterBounds2(0,0), m_parameterBounds2(0,1) ); m_flipSide2=true; }
        }
    }


}

template <class T>
std::ostream& gsRemapInterface<T>::print(std::ostream& os) const
{
    os << "gsRemapInterface:"
       << "\n    First side:      " << m_side1
       << "\n    Second side:     " << m_side2
       << "\n    Is Affine:       " << ( m_isAffine   ? "yes" : "no")
       << "\n    Matching:        " << ( m_isMatching ? "yes" : "no")
       << "\n    Flip side 2:     " << ( m_flipSide2  ? "yes" : "no")
       << "\n    Bounding box 1:\n" << m_parameterBounds1
       << "\n    Bounding box 2:\n" << m_parameterBounds2;


    for (size_t i=0; i<m_breakpoints.size(); ++i)
    {
        os << "\n    Beakpoints:    ";
        if ( m_breakpoints[i].size() <= 10 )
        {
            for (size_t j=0; j<m_breakpoints[i].size(); ++j)
                os << "  " << m_breakpoints[i][j];
        }
        else
        {
            for (size_t j=0; j<5; ++j)
                os << "  " << m_breakpoints[i][j];
            os << "  ...";
            for (size_t j=m_breakpoints[i].size()-5; j<m_breakpoints[i].size(); ++j)
                os << "  " << m_breakpoints[i][j];
        }
    }
    os << "\n";
    return os;
}


} // End namespace gismo
