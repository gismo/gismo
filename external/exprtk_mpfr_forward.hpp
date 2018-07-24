/*
 **************************************************************
 *         C++ Mathematical Expression Toolkit Library        *
 *                                                            *
 * MPFR Adaptor type                                          *
 * Authors: Arash Partow and Pavel Holoborodko                *
 * URL: http://www.partow.net/programming/exprtk/index.html   *
 *                                                            *
 * Copyright notice:                                          *
 * Free use of the Mathematical Expression Toolkit Library is *
 * permitted under the guidelines and in accordance with the  *
 * most current version of the MIT License.                   *
 * http://www.opensource.org/licenses/MIT                     *
 *                                                            *
 **************************************************************
*/


#ifndef EXPRTK_MPFRREAL_FORWARD_HPP
#define EXPRTK_MPFRREAL_FORWARD_HPP


#include <string>
#include <mpreal.h>

namespace exprtk
{
   namespace details
   {
      namespace numeric { namespace details
      {
         struct mpfrreal_type_tag;

         template <typename T> inline T const_pi_impl(mpfrreal_type_tag);
         template <typename T> inline T const_e_impl (mpfrreal_type_tag);
      }}

      inline bool is_true (const mpfr::mpreal& v);
      inline bool is_false(const mpfr::mpreal& v);

      template <typename Iterator>
      inline bool string_to_real(Iterator& itr_external, const Iterator end, mpfr::mpreal& t, numeric::details::mpfrreal_type_tag);

   }

   namespace rtl { namespace io
   {
      namespace details
      {
         inline void print_type(const std::string&, const mpfr::mpreal& v, exprtk::details::numeric::details::mpfrreal_type_tag);
      }
   }}

   using details::is_true;
}

#endif
