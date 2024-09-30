/** @file exprtk_codi_forward.h

    @brief Provides an exprtk adaptor for CoDiPack
    arithmetic types of autodiff

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#define TYPE_TAG_(CODI_TYPE) CODI_TYPE ## _type_tag
#define TYPE_TAG(CODI_TYPE) TYPE_TAG_(CODI_TYPE)
#define CODI_TYPE_TAG TYPE_TAG(CODI_TYPE)

namespace exprtk
{
namespace details
{
namespace numeric
{
namespace details
{

struct CODI_TYPE_TAG;

template <typename T> inline T const_pi_impl(CODI_TYPE_TAG);
template <typename T> inline T const_e_impl (CODI_TYPE_TAG);
} } // namespace details // namespace numeric

inline bool is_true (const CODI_TYPE& v);
inline bool is_false(const CODI_TYPE& v);

template <typename Iterator>
inline bool string_to_real(Iterator& itr_external, const Iterator end, CODI_TYPE& t, numeric::details::CODI_TYPE_TAG);

} // namespace details

namespace helper
{
namespace details
{
inline void print_type(const std::string&, const CODI_TYPE& v, exprtk::details::numeric::details::CODI_TYPE_TAG);
} } // namespace details // namespace helper

using details::is_true;
} // namespace exprtk
