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
                                      const gsOptionList      & opt)
    : m_g1(&(mp[bi.first().patch])), m_g2(&(mp[bi.second().patch])),
      m_b1(&(basis[bi.first().patch])), m_b2(&(basis[bi.second().patch])),
      m_bi(bi),
      m_isMatching(true), m_isAffine(true),
      m_equalityTolerance(opt.getReal("EqualityTolerance")), m_newtonTolerance(opt.getReal("NewtonTolerance"))
{
    const index_t checkAffine             = opt.getInt("CheckAffine");
    const index_t numSamplePoints         = opt.getInt("NumSamplePoints");
    const index_t intervalsOfFittingCurve = opt.getInt("IntervalsOfFittingCurve");
    const index_t degreeOfFittingCurve    = opt.getInt("DegreeOfFittingCurve");

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
gsMatrix<T> determineCorners(const gsGeometry<T> & geo, const boxSide s)
{
    gsMatrix<T> pr = geo.parameterRange();
    const index_t dim = pr.rows();
    gsVector<T> lower(dim);
    gsVector<T> upper(dim);
    gsVector<unsigned> numberGridPoints(dim);
    for (index_t i = 0; i<dim; ++i)
    {
        if (s.direction()==i)
        {
            lower[i] = upper[i] = pr( i, s.parameter() );
            numberGridPoints[i] = 1;
        }
        else
        {
            lower[i] = pr(i, 0);
            upper[i] = pr(i, 1);
            numberGridPoints[i] = 2;
        }
    }
    return gsPointGrid(lower,upper,numberGridPoints);
}

template <class T>
void transferCorners(const gsMatrix<T> &corners, const gsGeometry<T> &g1, const gsGeometry<T> &g2,
    const T newtonTolerance, const T equalityTolerance, gsMatrix<T> &cornerstransferred, gsVector<bool> &converged)
{
    cornerstransferred.resize(corners.rows(), corners.cols());
    converged.resize(corners.cols());

    gsVector<T> cornerTransferred = g2.parameterRange().col(0);
    for (index_t i=0; i<corners.cols(); ++i)
    {
        gsVector<T> cornerPhysical = g1.eval(corners.col(i));
        g2.newtonRaphson(cornerPhysical, cornerTransferred, true, newtonTolerance, 100);
        cornerstransferred.col(i) = cornerTransferred;
        converged[i] = (cornerPhysical - g2.eval(cornerTransferred)).norm() < equalityTolerance;
    }
}

template <class Vector, class T>
void widenParameterBounds(const Vector &point, gsMatrix<T> &parameterBounds)
{
    if (parameterBounds.cols()==0)
    {
        parameterBounds.resize(point.rows(),2);
        parameterBounds.col(0) = point;
        parameterBounds.col(1) = point;
    }
    else
    {
        for (index_t i=0; i<parameterBounds.rows(); ++i)
        {
            parameterBounds(i,0) = std::min( parameterBounds(i,0), point(i,0) );
            parameterBounds(i,1) = std::max( parameterBounds(i,1), point(i,0) );
        }
    }

}

} // end anonymous namespace

template <class T>
void gsRemapInterface<T>::constructInterfaceBox()
{
    gsMatrix<T> corners1 = determineCorners(*m_g1,m_bi.first());
    gsMatrix<T> corners2 = determineCorners(*m_g2,m_bi.second());

    // Transfer the corners to the other patch and determine if the Newton
    // did converge
    gsMatrix<T> corners1transferredTo2, corners2transferredTo1;
    gsVector<bool> converged1, converged2;
    transferCorners(corners1, *m_g1, *m_g2, m_newtonTolerance, m_equalityTolerance, corners1transferredTo2, converged1);
    transferCorners(corners2, *m_g2, *m_g1, m_newtonTolerance, m_equalityTolerance, corners2transferredTo1, converged2);

    // The parameter bounds are set such that all corners that are located on
    // the (common part of the) interface are in the parameterBound
    m_isMatching = true;
    for (index_t i=0; i<corners1.cols(); ++i)
    {
        if (converged1[i])
        {
            widenParameterBounds(corners1.col(i), m_parameterBounds1);
            widenParameterBounds(corners1transferredTo2.col(i), m_parameterBounds2);
        }
        else
        {
            m_isMatching = false;
        }
    }

    for (index_t i=0; i<corners2.cols(); ++i)
    {
        if (converged2[i])
        {
            widenParameterBounds(corners2.col(i), m_parameterBounds2);
            widenParameterBounds(corners2transferredTo1.col(i), m_parameterBounds1);
        }
        else
        {
            m_isMatching = false;
        }
    }

    // Make sure that the proper value is chosen in the normal direction
    const index_t d1 = m_bi.first().direction();
    m_parameterBounds1.row(d1).setConstant(m_g1->parameterRange()(d1,m_bi.first().parameter()));

    const index_t d2 = m_bi.second().direction();
    m_parameterBounds2.row(d2).setConstant(m_g2->parameterRange()(d2,m_bi.second().parameter()));

    GISMO_ASSERT ( m_parameterBounds1.cols()&&m_parameterBounds2.cols(),
        "gsRemapInterface<T>::constructInterfaceBox: Could not find an interface.");
    for (index_t j=0; j<domainDim(); ++j)
    {
        GISMO_ASSERT ( j==m_bi.first().direction() ||m_parameterBounds1(j,0) < m_parameterBounds1(j,1),
            "gsRemapInterface<T>::constructInterfaceBox: Could not find an interface.");
        GISMO_ASSERT ( j==m_bi.second().direction()||m_parameterBounds2(j,0) < m_parameterBounds2(j,1),
            "gsRemapInterface<T>::constructInterfaceBox: Could not find an interface.");
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
        GISMO_ASSERT ( (physPoints.col(i)-m_g2->eval(pointMapped)).norm() <= m_equalityTolerance,
            "gsRemapInterface<T>::constructFittingCurve: Newton did not converge." );
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

template <class T, class Vector>
inline bool inRange(const gsMatrix<T> & range,
                    const Vector      & point)
{
    for (index_t i=0; i<range.rows(); ++i)
        if ( point(i,0) < range(i,0) || point(i,0) > range(i,1) )
            return false;
    return true;
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
        gsVector<T> breakpoint_transferred = m_parameterBounds1.col(0); // initial guess
        for (; domIt2->good(); domIt2->next())
        {
            if (inRange(m_parameterBounds2, domIt2->upperCorner()))
            {
                gsMatrix<T> breakpoint_phys = m_g2->eval(domIt2->upperCorner());

                m_g1->newtonRaphson( breakpoint_phys, breakpoint_transferred, true, m_newtonTolerance, 100 );
                GISMO_ASSERT ( (breakpoint_phys-m_g1->eval(breakpoint_transferred)).norm() <= m_equalityTolerance,
                    "gsRemapInterface::constructBreaks: Newton method did not find point in parameter domain."
                    << (breakpoint_phys-m_g1->eval(breakpoint_transferred)).norm() );

                addBreaks(m_breakpoints, m_parameterBounds1, breakpoint_transferred, m_equalityTolerance);
            }
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
