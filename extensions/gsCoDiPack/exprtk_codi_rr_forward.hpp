/** @file exprtk_codi_rr_forward.hpp

    @brief Provides an exprtk adaptor for the CoDiPack
    RealReverse arithmetic types of autodiff

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <string>
#include <gsCoDiPack/gsCoDiPack.h>

namespace exprtk
{

namespace details
{
namespace numeric { namespace details {

struct codi_rr_type_tag;

template <typename T> inline T const_pi_impl(codi_rr_type_tag);
template <typename T> inline T const_e_impl (codi_rr_type_tag);
} } // namespace details // namespace numeric

inline bool is_true (const codi::RealReverse& v);
inline bool is_false(const codi::RealReverse& v);

template <typename Iterator>
inline bool string_to_real(Iterator& itr_external, const Iterator end, codi::RealReverse& t, numeric::details::codi_rr_type_tag);

}

namespace helper
{
namespace details
{
inline void print_type(const std::string&, const codi::RealReverse& v, exprtk::details::numeric::details::codi_rr_type_tag);
} } // namespace details // namespace helper

using details::is_true;
}
