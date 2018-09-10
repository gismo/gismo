/** @file gsField.hpp

    @brief Provides implementation of the Field class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Vogl
*/

#pragma once

namespace gismo
{


/// Computes the L2-distance between the two fields, on the physical domain
template <class T>
T gsField<T>::distanceL2(gsField<T> const & field, int numEvals) const
{
    return 0; //computeL2Distance(*this, field, numEvals);
}

/// Computes the L2-distance between the field and a function \a func on the physical domain
template <class T>
T gsField<T>::distanceL2(gsFunction<T> const & func,
                         bool isFunc_param,
                         int numEvals) const
{
    if (m_parametric) // isogeometric field
        return 0; //igaFieldL2Distance(*this, func, isFunc_param);
    else
        return 0; //computeL2Distance(*this, func, isFunc_param, numEvals);
}

/// Computes the L2-distance between the field and a function \a
/// func on the physical domain, using mesh from B
template <class T>
T gsField<T>::distanceL2(gsFunction<T> const & func,
                         gsMultiBasis<T> const & B,
                         bool isFunc_param,
                         int numEvals) const
{
    if (m_parametric) // isogeometric field
        return 0;//igaFieldL2Distance(*this, func, B, isFunc_param);
    else
        return 0;//computeL2Distance(*this, func, isFunc_param, numEvals);
}

/// Computes the H1-distance between the field and a function \a
/// func on the physical domain
template <class T>
T gsField<T>::distanceH1(gsFunction<T>
                         const & func,
                         bool isFunc_param,
                         int) const
{
    if (m_parametric) // isogeometric field
        return 0; //igaFieldH1Distance(*this, func, isFunc_param);
    else
    {
        gsWarn << "H1 seminorm not implemented.\n";
        return -1;
    }
}

/// Computes the H1-distance between the field and a function \a
/// func on the physical domain, using mesh from B
template <class T>
T gsField<T>::distanceH1(gsFunction<T>
                         const & func,
                         gsMultiBasis<T> const & B,
                         bool isFunc_param,
                         int) const
{
    if (m_parametric) // isogeometric field
        return 0; //igaFieldH1Distance(*this, func, B,isFunc_param);
    else
    {
        gsWarn << "H1 seminorm not implemented.\n";
        return -1;
    }
}

/// Computes the DG-distance between the field and a function \a
/// func on the physical domain
template <class T>
T gsField<T>::distanceDG(gsFunction<T>
                         const & func,
                         bool isFunc_param,
                         int) const
{
    if (m_parametric) // isogeometric field
        return 0; //igaFieldDGDistance(*this, func, isFunc_param);
    else
    {
        gsWarn << "DG norm not implemented.\n";
        return -1;
    }
}

} // namespace gismo
