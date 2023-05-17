/** @file gsHNeighborhood.h

    @brief Provides a struct for the neighborhood for admissibility

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst (2019 - ..., TU Delft)
*/

#pragma once

namespace gismo
{

enum class gsHNeighborhood 
{
    None = -1,
    Automatic = 0,
    T = 1,
    H = 2
};

} // namespace gismo