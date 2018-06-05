/** @file exprtk_unum_posit_forward.hpp

    @brief Provides an exprtk adaptor for the Unum Posit arithmetic
    type

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): M. Moller
*/

#pragma once

#include <string>
#include <gsUnum/gsUnum.h>

typedef sw::unum::posit<32,2> posit_32_2;

namespace exprtk
{

namespace details
{
namespace numeric { namespace details { struct unum_posit_type_tag; } }

inline bool is_true (const posit_32_2& v);
inline bool is_false(const posit_32_2& v);

template <typename Iterator>
inline bool string_to_real(Iterator& itr_external, const Iterator end, posit_32_2& t, numeric::details::unum_posit_type_tag);

}

namespace helper
{
namespace details
{
inline void print_type(const std::string&, const posit_32_2& v, exprtk::details::numeric::details::unum_posit_type_tag);
}
}

using details::is_true;
}
