/** @file gsFieldDistance.h

    @brief Distances between a gsField and a function (L2, H1- and
    H2-seminorm), integrated with the expression evaluator.

    These used to be implemented inside gsField (gsContainers). They are
    the only reason gsContainers needed gsExpressions/gsAssembler, so they
    live here as free functions; the gsField::distance* members forward
    to them and keep working unchanged.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Vogl
*/

#pragma once

#include <gsContainers/gsField.h>

namespace gismo
{

// The first declarations of these function templates (and their default
// arguments) are in gsContainers/gsField.h, ahead of class gsField.

/// Computes the L2-distance between two fields, on the physical domain
template<class T>
T distanceL2(const gsField<T> & field, const gsField<T> & other);

/// Computes the L2-distance between \a field and a function \a func
/// on the physical domain
template<class T>
T distanceL2(const gsField<T> & field, const gsFunctionSet<T> & func,
             bool isFunc_param);

/// Computes the L2-distance between \a field and a function \a func
/// on the physical domain, using the mesh of \a B
template<class T>
T distanceL2(const gsField<T> & field, const gsFunctionSet<T> & func,
             const gsMultiBasis<T> & B, bool isFunc_param);

/// Computes the H1-seminorm of the difference between \a field and a
/// function \a func on the physical domain
template<class T>
T distanceH1(const gsField<T> & field, const gsFunctionSet<T> & func,
             bool isFunc_param);

/// Computes the H1-seminorm of the difference between \a field and a
/// function \a func on the physical domain, using the mesh of \a B
template<class T>
T distanceH1(const gsField<T> & field, const gsFunctionSet<T> & func,
             const gsMultiBasis<T> & B, bool isFunc_param);

/// Computes the H2-seminorm of the difference between \a field and a
/// function \a func on the physical domain
template<class T>
T distanceH2(const gsField<T> & field, const gsFunctionSet<T> & func,
             bool isFunc_param);

/// Computes the DG-distance between \a field and a function \a func
/// on the physical domain (not implemented)
template<class T>
T distanceDG(const gsField<T> & field, const gsFunctionSet<T> & func,
             bool isFunc_param);

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFieldDistance.hpp)
#endif
