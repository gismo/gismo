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
                                      const index_t             checkAffine,
                                      const index_t             numSamplePoints,
                                      const index_t             intervalsOfFittingCurve,
                                      const index_t             degreeOfFittingCurve,
                                      const T                   equalityTolerance,
                                      const T                   newtonTolerance)
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
        m_isAffine = estimateReparamError(checkAffine) < m_equalityTolerance;

    if (!m_isAffine)
    {
        GISMO_ENSURE( m_isAffine || domainDim() <= 2,
            "gsRemapInterface: Can handle non-affine interfaces only for 2 dimensions." );

        constructFittingCurve(numSamplePoints, intervalsOfFittingCurve, degreeOfFittingCurve);
    }

    constructBreaks();

}

namespace {
template <class T>
gsMatrix<T> determineParameterBounds(const gsGeometry<T> & geo, const boxSide s)
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
gsMatrix<T> transferParameterBounds(const index_t         direction1,
                                    const gsGeometry<T> & g1,
                                    const gsGeometry<T> & g2,
                                    const gsMatrix<T>   & parameterBounds1,
                                    const gsMatrix<T>   & parameterBounds2,
                                    const T               solverTolerance)
{
    gsVector<T> transfered[2];
    for (index_t j=0; j<2; ++j)
    {
        transfered[j] = parameterBounds1.col(j); // initial guess
        g1.newtonRaphson( g2.eval(parameterBounds2).col(j), transfered[j], true, solverTolerance, 100 );
        //TODO: We might add such a statement:
        //transfered[j][direction1] = parameterBounds1(direction1,j);
    }

    gsMatrix<T> result(transfered[0].rows(), 2);
    for (index_t i=0; i<transfered[0].rows(); ++i)
    {
        result(i,0) = std::min( transfered[0][i], transfered[1][i] );
        result(i,1) = std::max( transfered[0][i], transfered[1][i] );
        // TODO: store orientation and use it for setup of affine mapping
    }
    return result;

}

template <class T>
void setToIntersection(bool & matching, T & min1, T & max1, const T min2, const T max2, const T tol)
{
    GISMO_ASSERT( min2<max1+tol && max2>min1-tol, "gsRemapInterface: Cannot find interface: "
        "(" << min1 << "," << max1 << ") and (" << min2 << "," << max2 << ") do not overlap.");
    if (min2>min1+tol)
    {
        min1 = std::min( min2, max1 );
        matching = false;
    }
    if (max2<max1-tol)
    {
        max1 = std::max( max2, min1 );
        matching = false;
    }
}

} // end anonymous namespace

template <class T>
void gsRemapInterface<T>::constructInterfaceBox()
{
    m_parameterBounds1 = determineParameterBounds(*m_g1,m_bi.first());
    m_parameterBounds2 = determineParameterBounds(*m_g2,m_bi.second());

    gsMatrix<T> parameterBounds1transferredTo2 = transferParameterBounds(m_bi.second().direction(),
        *m_g2,*m_g1,m_parameterBounds2,m_parameterBounds1,m_newtonTolerance);
    gsMatrix<T> parameterBounds2transferredTo1 = transferParameterBounds(m_bi.first().direction(),
        *m_g1,*m_g2,m_parameterBounds1,m_parameterBounds2,m_newtonTolerance);

    gsInfo << "m_parameterBounds1:\n" << m_parameterBounds1 << "\n\n";
    gsInfo << "m_parameterBounds2:\n" << m_parameterBounds2 << "\n\n";
    gsInfo << "parameterBounds1transferredTo2:\n" << parameterBounds1transferredTo2 << "\n\n";
    gsInfo << "parameterBounds2transferredTo1:\n" << parameterBounds2transferredTo1 << "\n\n";

    for (index_t i=0; i<domainDim(); ++i)
    {
        setToIntersection(
            m_isMatching, // in&out
            m_parameterBounds1(i,0), m_parameterBounds1(i,1), // in&out
            parameterBounds2transferredTo1(i,0), parameterBounds2transferredTo1(i,1), // in
            m_equalityTolerance
        );
        setToIntersection(
            m_isMatching, // in&out
            m_parameterBounds2(i,0), m_parameterBounds2(i,1), // in&out
            parameterBounds1transferredTo2(i,0), parameterBounds1transferredTo2(i,1), // in
            m_equalityTolerance
        );
    }
}

template <class T>
T gsRemapInterface<T>::estimateReparamError(const index_t steps) const
{
    gsVector<T> lower = m_parameterBounds1.col(0);
    gsVector<T> upper = m_parameterBounds1.col(1);
    gsVector<unsigned> numberGridPoints = gsVector<unsigned>::Constant(domainDim(),2+steps);
    numberGridPoints[m_bi.first().direction()] = 1;
    gsMatrix<T> points1 = gsPointGrid(lower,upper,numberGridPoints);
    gsMatrix<T> points2;
    eval_into(points1, points2);

    return  (
                m_g1->eval(points1)
                -
                m_g2->eval(points2)
            ).template lpNorm<Eigen::Infinity>();
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
void gsRemapInterface<T>::constructFittingCurve(const index_t numSamplePoints,
                                                const index_t intervalsOfFittingCurve,
                                                const index_t degreeOfFittingCurve)
{
    // assert domainDim()==2 taken care by constructor

    const short_t direction1 = m_bi.first().direction();

    // In 2d, there is one tangential direction
    const short_t tangential1 = ! direction1;

    // In 2d, the second side could have same or flipped direction
    const bool flipSide2 = (m_bi.dirOrientation()(tangential1) == 0);

    // Setup of point grid on patch1
    gsVector<unsigned> numPoints = gsVector<unsigned>::Constant(2,numSamplePoints);
    numPoints[direction1] = 1;
    gsVector<T> lower = m_parameterBounds1.col(0);
    gsVector<T> upper = m_parameterBounds1.col(1);
    gsMatrix<T> points = gsPointGrid(lower, upper, numPoints);

    // Map points to physical domain
    gsMatrix<T> physPoints = m_g1->eval(points);

    // Map the points to patch2
    gsMatrix<T> B(numSamplePoints, 2);
    // Find a suitable start value for the Newton iteration. After the first iteration,
    // we just continue where we are.
    gsVector<T> pointMapped =  m_parameterBounds2.col(flipSide2 ? 1 : 0);
    for (index_t i = 0; i < numSamplePoints; ++i)
    {
        m_g2->newtonRaphson(physPoints.col(i), pointMapped, true, m_newtonTolerance, 100);
        B.row(i) = pointMapped.transpose();
    }

    // Use equidstant knot vector for fitting curve
    gsKnotVector<T> KV(m_parameterBounds1(tangential1, 0), m_parameterBounds1(tangential1, 1),
        intervalsOfFittingCurve, degreeOfFittingCurve+1);
    gsCurveFitting<T> fit(points.row(tangential1).transpose(), B, KV);
    fit.compute();

    m_intfMap = fit.curve().clone();

}

namespace {

template <class T, class Vector>
inline void addBreaks(std::vector< std::vector<T> > & breaks,
                      const gsMatrix<T>             & parameterBounds,
                      const Vector                  & point,
                      const T                         equalityTolerance)
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
bool checkIfOnInterface(const gsMatrix<T> & u,
                        const gsMatrix<T> & parameterBounds,
                        const T             equalityTolerance)
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
void gsRemapInterface<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
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
        gsMatrix<T> uProjected = u;
        projectOntoInterface(uProjected, m_parameterBounds1);

        // Since we are in 2D, there is exactly one tangential direction
        const index_t direction1 = m_bi.first().direction();
        const index_t tangential1 = ! direction1;
        m_intfMap->eval_into(uProjected.row(tangential1),result);

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
std::ostream& gsRemapInterface<T>::print(std::ostream & os) const
{
    os << "gsRemapInterface:"
       << "\n    First side:         " << m_bi.first()
       << "\n    Second side:        " << m_bi.second()
       << "\n    Is Affine:          " << ( m_isAffine   ? "yes" : "no")
       << "\n    Is Matching:        " << ( m_isMatching ? "yes" : "no")
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
    if (!m_isAffine)
        os << "\n    Fitting curve:\n" << *m_intfMap;
    os << "\n    Error of reparam:   " << estimateReparamError(20);

    os << "\n";
    return os;
}


} // End namespace gismo
