/** @file gsAutoDiffEigen.h

    @brief Provides Eigen integration for autodiff types
    
    This file should be included AFTER gsLinearAlgebra.h to ensure
    Eigen plugins are properly set up before autodiff's eigen.hpp
    files include Eigen/Core.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

*/

#pragma once

// This header must be included AFTER gsLinearAlgebra.h
// It provides Eigen traits for autodiff types

// Ensure Eigen plugins are defined if not already
#ifndef EIGEN_MATRIXBASE_PLUGIN
#define EIGEN_MATRIXBASE_PLUGIN <gsMatrix/gsMatrixAddons.h>
#endif
#ifndef EIGEN_PLAINOBJECTBASE_PLUGIN
#define EIGEN_PLAINOBJECTBASE_PLUGIN <gsMatrix/gsPlainObjectBaseAddons.h>
#endif

#ifndef Eigen
#define Eigen gsEigen
#define GISMO_AUTODIFF_DEFINED_EIGEN
#endif

#include <Eigen/Core>
#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/reverse/var/eigen.hpp>

#include <gsAutoDiff/gsVarAdaptor.h>

namespace Eigen {

// Delegate NumTraits to the underlying autodiff::var type
template<typename T>
struct NumTraits<gismo::VarAdaptor<T>> : NumTraits<typename gismo::VarAdaptor<T>::Variable>
{
    typedef gismo::VarAdaptor<T> Real;
    typedef gismo::VarAdaptor<T> NonInteger;
    typedef gismo::VarAdaptor<T> Nested;
    typedef gismo::VarAdaptor<T> Literal;
    enum { RequireInitialization = 1 };
};

// Define binary op traits for mixed arithmetic with scalar
template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<gismo::VarAdaptor<T>, T, BinaryOp> {
  typedef gismo::VarAdaptor<T> ReturnType;
};

template<typename T, typename BinaryOp>
struct ScalarBinaryOpTraits<T, gismo::VarAdaptor<T>, BinaryOp> {
  typedef gismo::VarAdaptor<T> ReturnType;
};

} // namespace Eigen

#ifdef GISMO_AUTODIFF_DEFINED_EIGEN
#undef Eigen
#undef GISMO_AUTODIFF_DEFINED_EIGEN
#endif
