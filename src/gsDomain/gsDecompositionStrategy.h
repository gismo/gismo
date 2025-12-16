/** @file gsDecompositionStrategy.h

    @brief Provides an enum for decomposition strategies.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
*/

#pragma once

namespace gismo
{

/// \brief Enum to control the decomposition strategy in gsCompositeDomain
/// \ingroup Domain
enum class decompositionStrategy
{
    tensor,                 //< Decompose based on tensor-product structure (default)
    optimalBalancing,       //< Decompose for optimal load balancing, allowing subdomains to span patches
    localOptimalBalancing   //< Decompose for optimal load balancing, but restricted to within patches
};

} // namespace gismo
