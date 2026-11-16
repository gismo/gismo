/** @file exprtk_autodiff_forward.h

    @brief Provides an exprtk adaptor for autodiff
    arithmetic types

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
*/

#define TYPE_TAG_(AUTODIFF_TYPE) AUTODIFF_TYPE ## _type_tag
#define TYPE_TAG(AUTODIFF_TYPE) TYPE_TAG_(AUTODIFF_TYPE)
#define AUTODIFF_TYPE_TAG TYPE_TAG(AUTODIFF_TYPE)

namespace exprtk
{
namespace details
{
namespace numeric
{
namespace details
{

struct AUTODIFF_TYPE_TAG;

template <typename T> inline T const_pi_impl(AUTODIFF_TYPE_TAG);
template <typename T> inline T const_e_impl (AUTODIFF_TYPE_TAG);

} } // namespace details // namespace numeric

inline bool is_true (const AUTODIFF_TYPE& v);
inline bool is_false(const AUTODIFF_TYPE& v);

template <typename Iterator>
inline bool string_to_real(Iterator& itr_external, const Iterator end, AUTODIFF_TYPE& t, numeric::details::AUTODIFF_TYPE_TAG);

} // namespace details

namespace helper
{
namespace details
{
inline void print_type(const std::string&, const AUTODIFF_TYPE& v, exprtk::details::numeric::details::AUTODIFF_TYPE_TAG);
} } // namespace details // namespace helper

using details::is_true;
} // namespace exprtk
