/** @file

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/

#pragma once

#include <iostream>

# include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

enum gsRBFtype
{
    Constant = 0,
    Gaussian = 1,
    InverseQuadratic = 2,
    InverseMultiquadric = 3,
    // PolyHarmonic = 4,
    ThinPlate = 5,
    Bump = 6,
    Hat = 7
};

// template <class T>
// typename std::enable_if<gsRBFtype != RBF_GAUSSIAN, gsMatrix<>>::type
// gsRBF_impl(const gsMatrix<T> & r)

template <class T, gsRBFtype Type>
typename std::enable_if<Type == Constant, gsMatrix<T>>::type
gsRBF_impl(const gsMatrix<T> & r, T eps)
{
    gsMatrix<T> phi(1, r.cols());
    for (index_t k=0; k!=r.cols(); k++)
        phi(0,k) = (r(0,k) < eps) ? 1.0 : 0.0;
    return phi;
}

template <class T, gsRBFtype Type>
typename std::enable_if<Type == Gaussian, gsMatrix<T>>::type
gsRBF_impl(const gsMatrix<T> & r, T eps)
{
    return (-(eps*r).array().square()).exp().matrix();
}

template <class T, gsRBFtype Type>
typename std::enable_if<Type == InverseQuadratic, gsMatrix<T>>::type
gsRBF_impl(const gsMatrix<T> & r, T eps)
{
    return (1.0 / (1+ (eps*r).array().square())).matrix();
}

template <class T, gsRBFtype Type>
typename std::enable_if<Type == InverseMultiquadric, gsMatrix<T>>::type
gsRBF_impl(const gsMatrix<T> & r, T eps)
{
    return (1.0 / ((1+ (eps*r).array().square()).sqrt())).matrix();
}

// template <class T, gsRBFtype Type>
// typename std::enable_if<Type == PolyHarmonic, gsMatrix<T>>::type
// gsRBF_impl(const gsMatrix<T> & r, T eps)
// {

// }

template <class T, gsRBFtype Type>
typename std::enable_if<Type == ThinPlate, gsMatrix<T>>::type
gsRBF_impl(const gsMatrix<T> & r, T /* eps */)
{
    return (r.array().square() * r.array().log()).matrix();
}

template <class T, gsRBFtype Type>
typename std::enable_if<Type == Bump, gsMatrix<T>>::type
gsRBF_impl(const gsMatrix<T> & r, T eps)
{
    gsMatrix<T> phi(1, r.cols());
    for (index_t k=0; k!=r.cols(); k++)
        phi(0,k) = (r(0,k) < 1./eps) ? math::exp(-1./(1-math::pow(eps*r(0,k),2))) : 0.;

    phi = ((-1.0 / (1 - (eps*r).array().square()))).exp().matrix(); //.cwiseMax(0)?
    phi = (r.array() < 1./eps).select(phi, 0);
    return phi;
}

template <class T, gsRBFtype Type>
typename std::enable_if<Type == Hat, gsMatrix<T>>::type
gsRBF_impl(const gsMatrix<T> & r, T beta)
{
    gsMatrix<T> phi(1, r.cols());
    for (index_t k=0; k!=r.cols(); k++)
        phi(0,k) = (r(0,k) <= beta) ? 1 - r(0,k)/beta : 0.;
    return phi;
}

template <class T, gsRBFtype Type>
gsMatrix<T> gsRBF(const gsMatrix<T> & r, T eps)
{
    return gsRBF_impl<T, Type>(r, eps);
}

template <class T, gsRBFtype Type>
class gsRBFCurve : public gsFunction<T>
{
public:

    gsRBFCurve( const gsGeometry<T> &    curve,
                const T                  beta,
                const T                  eps,
                const T                  scaling = 1.0)
    :
    gsRBFCurve(gsMultiPatch<T>(curve),
               beta,
               eps,
               scaling)
    {
    }


    gsRBFCurve( const gsMultiPatch<T> &  curves,
                const T                  eps,
                const T                  beta,
                const T                  scaling = 1.0)
    :
    m_curves(curves),
    m_eps(eps),
    m_beta(beta),
    m_scaling(scaling)
    {
        // GISMO_ASSERT(m_curves.parDim() == 1, "RBF curve only works in 1D");
        m_curves.boundingBox(m_bbox);
        // tolerance for closest point search determined by the (volume)^1/geoDim()*default
        m_tol = math::pow((m_bbox.col(1) - m_bbox.col(0)).prod(), 1.0/m_curves.geoDim()) * 1e-6;
    }

    short_t domainDim() const { return m_curves.geoDim(); }
    short_t targetDim() const { return 1; }

    void setBeta(T beta) const { m_beta = beta; }
    void setEpsilon(T eps) const { m_eps = eps; }
    void setScaling(T scaling) const { m_scaling = scaling; }

    void eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
    {
        GISMO_ASSERT(u.rows() == m_curves.geoDim(), "Invalid input point.");
        result.resize(1, u.cols());
        gsVector<T> par;
        gsMatrix<T> dist(1, u.cols());
        bool outside = false; // no patch
        for (index_t k=0; k<u.cols(); ++k)
        {
            // Check if the point is inside the bounding box of the curve
            for (index_t i=0; i<m_bbox.rows(); ++i)
            {
                if (u(i,k) > m_bbox(i,0)-m_beta && u(i,k) < m_bbox(i,1) + m_beta)
                    outside = false;
                else
                {
                    outside = true;
                    break;
                }
            }
            if (outside)
                result(0,k) = 0;
            else
            {
                dist(0,k) = math::limits::max();
                for (size_t p=0; p!=m_curves.nPatches(); p++)
                    dist(0,k) = math::min(dist(0,k),m_curves.patch(p).closestPointTo(u.col(k),par,1e-3,false));
                result(0,k) = gsRBF_impl<T, Type>(dist.col(k), m_eps).value();
            }
        }

        result *= m_scaling;

        // result.row(0) = gsRBF<T, Type>(dist, m_eps);
    }

protected:
    const gsMultiPatch<T> m_curves;
    const T m_eps;
    const T m_beta;
    const T m_scaling;
          T m_tol;
    gsMatrix<T> m_bbox;
};

} // namespace gismo
