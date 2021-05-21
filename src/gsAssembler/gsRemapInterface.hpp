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
                                      index_t checkAffine,
                                      T equalityTolerance,
                                      T newtonTolerance)
    : m_g1(&(mp[bi.first().patch])), m_g2(&(mp[bi.second().patch])),
      m_b1(&(basis[bi.first().patch])), m_b2(&(basis[bi.second().patch])),
      m_bi(bi),
      m_isMatching(true), m_isAffine(true),
      m_equalityTolerance(equalityTolerance), m_newtonTolerance(newtonTolerance)
{
    GISMO_ASSERT( m_g1->geoDim()==m_g2->geoDim(), "gsRemapInterface: Dimensions do not agree." );

    // First we construct the affine mapping
    constructInterfaceBox();

    // Setup the affine mapping
    m_intfMap = gsAffineFunction<T>::make(
        bi.dirMap(m_bi.first()),
        bi.dirOrientation(m_bi.first()),
        m_parameterBounds1,
        m_parameterBounds2
        );

    // Next, we check (if so desired by called) if the affine mapping coincides with the real mapping
    GISMO_ASSERT( checkAffine > 0 || checkAffine == neverAffine || checkAffine == alwaysAffine,
        "gsRemapInterface: Parameter checkAffine has invalid value:" << checkAffine );
    if (checkAffine==neverAffine)
        m_isAffine = false;
    else if (checkAffine > 0)
        m_isAffine = checkIfAffine(checkAffine);

    if (!m_isAffine)
    {
        GISMO_ENSURE( m_isAffine || domainDim() <= 2,
            "gsRemapInterface: Can handle non-affine interfaces only for 2 dimensions." );

        constructFittingCurve();
    }

    constructBreaks();

}

namespace {
template <class T>
gsMatrix<T> determineParameterBounds(const gsGeometry<T> & geo, boxSide s)
{
    gsMatrix<T> pr = geo.parameterRange();
    const index_t dim = pr.rows();
    gsMatrix<T> result(dim, 2);
    for (index_t i = 0; i<dim; ++i)
    {
        if (s.direction()==i)
            result(i,0) = result(i,1) = pr( i, s.parameter() == false ? 0 : 1 );
        else
            result.row(i) = pr.row(i);
    }
    return result;
}

template <class T>
gsMatrix<T> transferParameterBounds(const gsGeometry<T> & g1, const gsGeometry<T> & g2,
                const gsMatrix<T> & parameterBounds1, const gsMatrix<T> & parameterBounds2, T solverTolerance)
{
    gsVector<T> transfered[2];
    for (index_t j=0; j<2; ++j)
    {
        transfered[j] = parameterBounds1.col(j); // initial guess
        g1.newtonRaphson( g2.eval(parameterBounds2).col(j), transfered[j], true, solverTolerance, 100 );
    }

    gsMatrix<T> result(transfered[0].rows(), 2);
    for (index_t i=0; i<transfered[0].rows(); ++i)
    {
        result(i,0) = std::min( transfered[0][i], transfered[1][i] );
        result(i,1) = std::max( transfered[0][i], transfered[1][i] );
    }
    return result;

}

} // end anonymous namespace

template <class T>
void gsRemapInterface<T>::constructInterfaceBox()
{
    m_parameterBounds1 = determineParameterBounds(*m_g1,m_bi.first());
    m_parameterBounds2 = determineParameterBounds(*m_g2,m_bi.second());

    gsMatrix<T> parameterBounds1trafsferedTo2
        = transferParameterBounds(*m_g2,*m_g1,m_parameterBounds2,m_parameterBounds1,m_newtonTolerance);
    gsMatrix<T> parameterBounds2trafsferedTo1
        = transferParameterBounds(*m_g1,*m_g2,m_parameterBounds1,m_parameterBounds2,m_newtonTolerance);

    for (index_t i=0; i<domainDim(); ++i)
    {
        if (parameterBounds2trafsferedTo1(i,0)>m_parameterBounds1(i,0)+m_equalityTolerance)
        {
            m_parameterBounds1(i,0) = parameterBounds2trafsferedTo1(i,0);
            m_isMatching = false;
        }
        if (parameterBounds2trafsferedTo1(i,1)<m_parameterBounds1(i,1)-m_equalityTolerance)
        {
            m_parameterBounds1(i,1) = parameterBounds2trafsferedTo1(i,1);
            m_isMatching = false;
        }

        if (parameterBounds1trafsferedTo2(i,0)>m_parameterBounds2(i,0)+m_equalityTolerance)
        {
            m_parameterBounds2(i,0) = parameterBounds1trafsferedTo2(i,0);
            m_isMatching = false;
        }
        if (parameterBounds1trafsferedTo2(i,1)<m_parameterBounds2(i,1)-m_equalityTolerance)
        {
            m_parameterBounds2(i,1) = parameterBounds1trafsferedTo2(i,1);
            m_isMatching = false;
        }
    }
}

template <class T>
bool gsRemapInterface<T>::checkIfAffine(index_t steps)
{
    gsVector<T> lower = m_parameterBounds1.col(0);
    gsVector<T> upper = m_parameterBounds1.col(1);
    gsVector<unsigned> numberGridPoints = gsVector<unsigned>::Constant(domainDim(),2+steps);
    numberGridPoints[m_bi.first().direction()] = 1;
    gsMatrix<T> points = gsPointGrid(lower,upper,numberGridPoints);

    return  (
                m_g1->eval(points)
                -
                m_g2->eval(m_intfMap->eval(points))
            ).norm() < (T)(m_equalityTolerance);
}

namespace {

// Takes the univariate parameter on interface and returns the corresponding (bivariate)
// points on parameter domain
template <class T>
void enrichToVector(const gsMatrix<T>   & u,
                    const index_t         direction,
                    const gsMatrix<T>   & parameterBounds,
                    gsMatrix <T>        & result)
{
    const index_t dim = 2;
    result.resize(dim, u.cols());
    for (index_t i=0; i<dim; ++i)
    {
        if (direction==i)
            result.row(i).setConstant( parameterBounds(i,0) );
        else
            result.row(i) = u;
    }
}

} // End anonyomous namespace

template <class T>
void gsRemapInterface<T>::constructFittingCurve()
{
    // assert domainDim()==2 taken care by constructor

    const index_t numIntervals = 11; // ?
    const index_t numGeometries = 2;

    // Check if side 2 is to be flipped
    bool flipSide2 = false;

    if(m_bi.first().direction() == 1)
    {
        if(m_bi.dirOrientation()(0) == 0)
        {
            if(m_bi.second().direction() == 0) //change v parameters
                { flipSide2=true; }
            else if(m_bi.second().direction() == 1) // change u parameters
                { flipSide2=true; }
        }
    }

    if(m_bi.first().direction() == 0)
    {
        if(m_bi.dirOrientation()(1) == 0)
        {
            if(m_bi.second().direction() == 0) //change v parameters
                { flipSide2=true; }
            else if(m_bi.second().direction() == 1) // change u parameters
                { flipSide2=true; }
        }
    }

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
            if (m_bi.first().direction() == 1) // v is fixed
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
            if (m_bi.second().direction() == 1) // v is fixed
            {
                firstKnot = m_parameterBounds2(0, flipSide2 ? 1 : 0);
                lastKnot = m_parameterBounds2(0, flipSide2 ? 0 : 1);
            }
            else // u is fixed
            {
                firstKnot = m_parameterBounds2(1, flipSide2 ? 1 : 0);
                lastKnot = m_parameterBounds2(1, flipSide2 ? 0 : 1);
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
    enrichToVector<T>(t_vals.row(0), m_bi.first().direction(),  m_parameterBounds1, vals2dPatch1);
    enrichToVector<T>(t_vals.row(1), m_bi.second().direction(), m_parameterBounds2, vals2dPatch2);

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

        //gsVector<T> b_null = samples_right.col(col);
        gsVector<T> b_null = vals2dPatch2.col(col);

        // Pass on g2 if one wants to find a mapping from interface1 to interface2
        //gsMatrix<T> b = closestPoint(b_null, g2, samples_left.col(i));

        // this gives the same result as above
        m_g2->newtonRaphson(samples_left.col(i), b_null, true, m_newtonTolerance, 100);
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
    enrichToVector<T>(eval_points.row(1), m_bi.second().direction(), m_parameterBounds2, id);

    m_g2->eval_into(eval_fit, eval_orig);
    m_g2->eval_into(id, B2);
    //gsInfo << "b2: \n" << id.transpose() << " and eval_orig: \n" << eval_fit.transpose() << "\n";

    // do test
    /*
    for(index_t c = 0; c < eval_fit.cols(); c++)
    {
        if(std::isnan(eval_orig.col(c).squaredNorm()))
        {
            switch (m_bi.second().index())
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

namespace {

template <class T, class Vector>
inline void addBreaks(std::vector< std::vector<T> >& breaks, const gsMatrix<T>& parameterBounds,
    const Vector& point, const T equalityTolerance)
{
    const index_t dim = point.rows();
    for (index_t d=0; d<dim; ++d)
    {
        const T t = point(d,0);
        if ( parameterBounds(d,0) <= t && t <= parameterBounds(d,1) )
        {
            // As in gsSortedVector::push_sorted_unique
            typename std::vector<T>::iterator pos = std::lower_bound(breaks[d].begin(), breaks[d].end(), t-equalityTolerance );
            if ( pos == breaks[d].end() || *pos > t+equalityTolerance ) // If not found
                breaks[d].insert(pos, t);
        }
    }
}

} // end anonymous namespace

template <class T>
void gsRemapInterface<T>::constructBreaks()
{
    m_breakpoints.resize(domainDim());

    const typename gsBasis<T>::domainIter domIt1 = m_b1->makeDomainIterator(m_bi.first());
    addBreaks(m_breakpoints, m_parameterBounds1, m_parameterBounds1.col(0), m_equalityTolerance);
    for (; domIt1->good(); domIt1->next())
        addBreaks(m_breakpoints, m_parameterBounds1, domIt1->upperCorner(), m_equalityTolerance);
    addBreaks(m_breakpoints, m_parameterBounds1, m_parameterBounds1.col(1), m_equalityTolerance);

    const typename gsBasis<T>::domainIter domIt2 = m_b2->makeDomainIterator(m_bi.second());
    if (m_isAffine)
    {
        gsAffineFunction<T> intfMap_inverse(m_bi.dirMap(m_bi.second()), m_bi.dirOrientation(m_bi.second()),
            m_parameterBounds2, m_parameterBounds1);

        for (; domIt2->good(); domIt2->next())
            addBreaks(m_breakpoints, m_parameterBounds1, intfMap_inverse.eval(domIt2->upperCorner()), m_equalityTolerance);
    }
    else
    {
        gsVector<T> breakpoint_transfered = m_parameterBounds1.col(0); // initial guess
        for (; domIt2->good(); domIt2->next())
        {
            m_g1->newtonRaphson( m_g2->eval(domIt2->upperCorner()), breakpoint_transfered, true, m_newtonTolerance, 100 );
            addBreaks(m_breakpoints, m_parameterBounds1, breakpoint_transfered, m_equalityTolerance);
        }
    }

}

namespace {

// Check if incoming parameters are located on interface
template <class T>
bool checkIfOnInterface(const gsMatrix<T> & u, const gsMatrix<T> & parameterBounds, const T equalityTolerance)
{
    for (index_t r=0; r<u.rows(); ++r)
    {
        const T begin = parameterBounds(r, 0);
        const T end = parameterBounds(r, 1);
        for (index_t c=0; c<u.cols(); ++c)
        {
            if ( u(r,c) < begin-equalityTolerance || u(r,c) > end+equalityTolerance )
                return false;
        }
    }

    return true;
}

// Map the parameterization onto interface
template <class T>
void projectOntoInterface(gsMatrix<T> & u, const gsMatrix<T> & parameterBounds)
{
    for(index_t r = 0; r < u.rows(); r++)
    {
        const T begin = parameterBounds(r, 0);
        const T end = parameterBounds(r, 1);

        for(index_t c = 0; c < u.cols(); c++)
            u(r,c) = std::min( std::max( u(r,c), begin ), end );
    }
}


} // End anonyomous namespace

template <class T>
void gsRemapInterface<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT( u.rows() == domainDim(), "gsRemapInterface::eval_into: "
        "The rows of the evaluation points must be equal to the dimension of the domain." );

    GISMO_ASSERT( checkIfOnInterface( u, m_parameterBounds1, m_equalityTolerance ),
        "gsRemapInterface::eval_into: The incoming coefficients are not located on the interface."
        <<"\nu=\n"<<u<<"\nm_parameterBounds1=\n"<<m_parameterBounds1);

    if (m_isAffine)
    {
        m_intfMap->eval_into(u, result);
    }
    else
    {
        gsMatrix<T> uprojected = u;
        projectOntoInterface(uprojected, m_parameterBounds1);

        // Note that the fitting curve is univariate, so we only consider the direction
        // that is orthogonal to the normal vector's direction. Since we are in 2D,
        // that direction is obtained using !direction1
        const index_t direction1 = m_bi.first().direction();
        m_intfMap->eval_into(uprojected.row(!direction1),result);

        projectOntoInterface(result, m_parameterBounds2);
    }
}


template <class T>
typename gsDomainIterator<T>::uPtr gsRemapInterface<T>::makeDomainIterator() const
{
    gsTensorDomainBoundaryIterator<T> * tdi = new gsTensorDomainBoundaryIterator<T> (*m_b1, m_bi.first());
    for (index_t i=0; i<domainDim(); ++i)
    {
        if (i!=m_bi.first().direction())
            tdi->setBreaks(m_breakpoints[i],i);
    }
    return typename gsDomainIterator<T>::uPtr(tdi);
}


template <class T>
std::ostream& gsRemapInterface<T>::print(std::ostream& os) const
{
    os << "gsRemapInterface:"
       << "\n    First side:         " << m_bi.first()
       << "\n    Second side:        " << m_bi.second()
       << "\n    Is Affine:          " << ( m_isAffine   ? "yes" : "no")
       << "\n    Matching:           " << ( m_isMatching ? "yes" : "no")
       << "\n    Bounding box 1 min: " << m_parameterBounds1.transpose().row(0)
       << "\n                   max: " << m_parameterBounds1.transpose().row(1)
       << "\n    Bounding box 2 min: " << m_parameterBounds2.transpose().row(0)
       << "\n                   max: " << m_parameterBounds2.transpose().row(1);

    for (size_t i=0; i<m_breakpoints.size(); ++i)
    {
        os << "\n    Beakpoints " << (i<3?char('x'+i):' ') << ":     ";
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
