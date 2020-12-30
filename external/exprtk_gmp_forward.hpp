/** @file exprtk_gmp_forward.hpp

    @brief Provides an exprtk adaptor for the GMP/mpq_class number type

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <string>
#include <gsMultiPrecision/gsMultiPrecision.h>

namespace exprtk
{

/*
// in case a math function is missing
#define GMP_EXTRA_STD_UNARY_FUNCTION(std_fun) template <class U> \
inline mpq_class std_fun(const __gmp_expr<mpq_t, U> & expr)      \
{return std::std_fun(mpq_class(expr).get_d());}
// GMP_EXTRA_STD_UNARY_FUNCTION(sqrt)
#undef GMP_EXTRA_STD_UNARY_FUNCTION
*/

namespace details
{
namespace numeric { namespace details {

struct mpq_type_tag;

template <typename T> inline T const_pi_impl(mpq_type_tag);
template <typename T> inline T const_e_impl (mpq_type_tag);

} }

inline bool is_true (const mpq_class& v);
inline bool is_false(const mpq_class& v);

template <typename Iterator>
inline bool string_to_real(Iterator& itr_external, const Iterator end, mpq_class& t, numeric::details::mpq_type_tag);

}

namespace rtl { namespace io {
namespace details
{
inline void print_type(const std::string&, const mpq_class& v, exprtk::details::numeric::details::mpq_type_tag);
}
} }

using details::is_true;
}
