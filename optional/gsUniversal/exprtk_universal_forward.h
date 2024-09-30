/** @file exprtk_universal_forward.hpp

    @brief Provides an exprtk adaptor for the Unum Posit arithmetic
    type

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#define TYPE_TAG_(UNIVERSAL_TYPE) UNIVERSAL_TYPE ## _type_tag
#define TYPE_TAG(UNIVERSAL_TYPE) TYPE_TAG_(UNIVERSAL_TYPE)
#define UNIVERSAL_TYPE_TAG TYPE_TAG(UNIVERSAL_TYPE)

namespace exprtk
{

namespace details
{
namespace numeric { namespace details {

struct UNIVERSAL_TYPE_TAG;

template <typename T> inline T const_pi_impl(UNIVERSAL_TYPE_TAG);
template <typename T> inline T const_e_impl (UNIVERSAL_TYPE_TAG);
} } // namespace details // namespace numeric

inline bool is_true (const UNIVERSAL_TYPE& v);
inline bool is_false(const UNIVERSAL_TYPE& v);

template <typename Iterator>
inline bool string_to_real(Iterator& itr_external, const Iterator end, UNIVERSAL_TYPE& t, numeric::details::UNIVERSAL_TYPE_TAG);

} // namespace details

namespace rtl { namespace io { namespace details {
inline void print_type(const std::string&, const UNIVERSAL_TYPE& v, exprtk::details::numeric::details::UNIVERSAL_TYPE_TAG);
} } } // namespace details // namespace io // namespace rtl

using details::is_true;
}
