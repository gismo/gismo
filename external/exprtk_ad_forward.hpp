/** @file exprtk_ad_forward.hpp

    @brief Provides an exprtk adaptor for the DScalar1/DScalar2
    arithmetic types of autodiff

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <string>
#include <gsAutoDiff.h>

namespace exprtk
{

#ifndef DScalar
typedef typename gismo::ad::DScalar1<real_t, -1> DScalar;
//typedef gismo::ad::DScalar1<real_t, 2> DScalar;
//typedef gismo::ad::DScalar1<real_t, 3> DScalar;

//typedef gismo::ad::DScalar2<real_t, -1> DScalar;
//typedef gismo::ad::DScalar2<real_t, 2> DScalar;
//typedef gismo::ad::DScalar2<real_t, 3> DScalar;
#endif

typedef DScalar::Scalar DScalarValue;

namespace details
{
namespace numeric { namespace details {

struct ad_type_tag;

template <typename T> inline T const_pi_impl(ad_type_tag);
template <typename T> inline T const_e_impl (ad_type_tag);
} } // namespace details // namespace numeric

inline bool is_true (const DScalar& v);
inline bool is_false(const DScalar& v);

template <typename Iterator>
inline bool string_to_real(Iterator& itr_external, const Iterator end, DScalar& t, numeric::details::ad_type_tag);

} // namespace details

namespace rtl { namespace io { namespace details {
inline void print_type(const std::string&, const DScalar& v, exprtk::details::numeric::details::ad_type_tag);
} } } // namespace details // namespace io // namespace rtl

using details::is_true;
}


