/** @file gsHSplinesParaview.h

    @brief Paraview export of gsHSplines types (moved from gsIO/gsWriteParaview,
    modularization stream S3 step A8: type-specific visualization lives
    with the type's module, the base IO module stays type-blind).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsWriteParaview.h>
#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBoxContainer.h>

#define NS 1000 // default sampling, mirrors gsWriteParaview.h

namespace gismo
{

/// \brief Export gsHBox to paraview files
///
/// \param box a gsHBox
/// \param fn filename where paraview file is written
/// \param mode controls the output format: 0 (colored by level), 1 (colored by error)
template<short_t d, class T>
void gsWriteParaview(const gsHBox<d,T> & box, std::string const & fn, short_t mode = 0);

/// \param mode controls the output format: 0 (colored by level), 1 (colored by error)
template<short_t d, class T>
void gsWriteParaview(const gsHBoxContainer<d,T> & box, std::string const & fn, short_t mode = 0);

/// Export a gsHBox
template<short_t d, class T>
void writeSingleHBox(const gsHBox<d,T> & box, std::string const & fn);

} // namespace gismo

#undef NS

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHSplinesParaview.hpp)
#endif
