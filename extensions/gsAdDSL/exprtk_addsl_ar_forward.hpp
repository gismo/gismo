/** @file exprtk_addsl_ar_forward.hpp

    @brief Provides an exprtk adaptor for the adDSL
    DSLActiveReal arithmetic types of autodiff

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <string>
#include <gsAdDSL/gsAdDSL.h>

namespace exprtk
{

namespace details
{
namespace numeric { namespace details {

struct addsl_ar_type_tag;

template <typename T> inline T const_pi_impl(addsl_ar_type_tag);
template <typename T> inline T const_e_impl (addsl_ar_type_tag);
} } // namespace details // namespace numeric

inline bool is_true (const adDSL::DSLActiveReal& v);
inline bool is_false(const adDSL::DSLActiveReal& v);

template <typename Iterator>
inline bool string_to_real(Iterator& itr_external, const Iterator end, adDSL::DSLActiveReal& t, numeric::details::addsl_ar_type_tag);

}

namespace helper
{
namespace details
{
inline void print_type(const std::string&, const adDSL::DSLActiveReal& v, exprtk::details::numeric::details::addsl_ar_type_tag);
} } // namespace details // namespace helper

using details::is_true;
}
