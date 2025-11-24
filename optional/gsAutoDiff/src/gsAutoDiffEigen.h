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

#define Eigen gsEigen
#include <autodiff/forward/real/eigen.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#undef Eigen
