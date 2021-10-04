/** @file exprtk_unum_posit_forward.hpp

    @brief Provides an exprtk adaptor for the Unum Posit arithmetic
    type

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#define TYPE_TAG_(UNUM_TYPE) UNUM_TYPE ## _type_tag
#define TYPE_TAG(UNUM_TYPE) TYPE_TAG_(UNUM_TYPE)
#define UNUM_TYPE_TAG TYPE_TAG(UNUM_TYPE)

namespace exprtk
{

namespace details
{
namespace numeric { namespace details {

struct UNUM_TYPE_TAG;

template <typename T> inline T const_pi_impl(UNUM_TYPE_TAG);
template <typename T> inline T const_e_impl (UNUM_TYPE_TAG);
} } // namespace details // namespace numeric

inline bool is_true (const UNUM_TYPE& v);
inline bool is_false(const UNUM_TYPE& v);

template <typename Iterator>
inline bool string_to_real(Iterator& itr_external, const Iterator end, UNUM_TYPE& t, numeric::details::UNUM_TYPE_TAG);

} // namespace details

namespace rtl { namespace io { namespace details {
inline void print_type(const std::string&, const UNUM_TYPE& v, exprtk::details::numeric::details::UNUM_TYPE_TAG);
} } } // namespace details // namespace io // namespace rtl

using details::is_true;
}
