/** @file gsAutoDiffEigen.h

    @brief Provides Eigen integration for autodiff types
    
    This file provides Eigen support for autodiff types.
    The Eigen library already provides NumTraits specializations for autodiff types
    when including the appropriate autodiff headers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

*/

#pragma once

// This header must be included AFTER gsLinearAlgebra.h to ensure
// Eigen plugins are properly set up before autodiff's eigen.hpp
// files include Eigen/Core.

// Ensure Eigen plugins are defined if not already
#ifndef EIGEN_MATRIXBASE_PLUGIN
#define EIGEN_MATRIXBASE_PLUGIN <gsMatrix/gsMatrixAddons.h>
#endif
#include <Eigen/Core>

// Suppress warnings from autodiff library headers
#ifdef __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"
#endif

#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/reverse/var/eigen.hpp>

#ifdef __GNUC__
#pragma GCC diagnostic pop
#endif

