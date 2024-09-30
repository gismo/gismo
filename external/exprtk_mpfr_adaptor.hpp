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


#ifndef EXPRTK_MPFRREAL_ADAPTOR_HPP
#define EXPRTK_MPFRREAL_ADAPTOR_HPP


#include <string>
#include <mpreal.h>

#include "exprtk_mpfr_forward.hpp"
#include "exprtk.hpp"

namespace exprtk
{
   namespace details
   {

      namespace constant_mpfr
      {
         static const mp_prec_t  mpfr_precision = 128;
         static const mp_rnd_t   mpfr_round     = mpfr::mpreal::get_default_rnd();

         static const mpfr::mpreal e       = mpfr::const_euler(mpfr_precision,mpfr_round);
         static const mpfr::mpreal pi      = mpfr::const_pi   (mpfr_precision,mpfr_round);
         static const mpfr::mpreal pi_2    = mpfr::const_pi   (mpfr_precision,mpfr_round) / mpfr::mpreal(  2.0);
         static const mpfr::mpreal pi_4    = mpfr::const_pi   (mpfr_precision,mpfr_round) / mpfr::mpreal(  4.0);
         static const mpfr::mpreal pi_180  = mpfr::const_pi   (mpfr_precision,mpfr_round) / mpfr::mpreal(180.0);
         static const mpfr::mpreal _1_pi   = mpfr::mpreal(  1.0) / mpfr::const_pi(mpfr_precision,mpfr_round);
         static const mpfr::mpreal _2_pi   = mpfr::mpreal(  2.0) / mpfr::const_pi(mpfr_precision,mpfr_round);
         static const mpfr::mpreal _180_pi = mpfr::mpreal(180.0) / mpfr::const_pi(mpfr_precision,mpfr_round);
         static const mpfr::mpreal log2    = mpfr::const_log2 (mpfr_precision,mpfr_round);
         static const mpfr::mpreal sqrt2   = mpfr::sqrt(mpfr::mpreal(2.0));
      }

      namespace numeric
      {
         namespace details
         {
            struct mpfrreal_type_tag {};

            template<> struct number_type<mpfr::mpreal> { typedef mpfrreal_type_tag type; };

            template <>
            struct epsilon_type<mpfrreal_type_tag>
            {
               static inline mpfr::mpreal value()
               {
                  static const mpfr::mpreal epsilon =
                  #ifndef exprtk_use_mpfr_epsilon
                  mpfr::mpreal(1.0) / mpfr::mpreal (1e+20);
                  #else
                  mpfr::machine_epsilon();
                  #endif
                  return epsilon;
               }
            };

            inline bool is_nan_impl(const mpfr::mpreal& v, mpfrreal_type_tag)
            {
               return mpfr::isnan(v);
            }

            template <typename T>
            inline int to_int32_impl(const T& v, mpfrreal_type_tag)
            {
               return static_cast<int>(v.toLong());
            }

            template <typename T>
            inline long long to_int64_impl(const T& v, mpfrreal_type_tag)
            {
               return static_cast<long long int>(v.toLLong());
            }

            template <typename T> inline T   abs_impl(const T& v, mpfrreal_type_tag) { return mpfr::abs  (v); }
            template <typename T> inline T  acos_impl(const T& v, mpfrreal_type_tag) { return mpfr::acos (v); }
            template <typename T> inline T acosh_impl(const T& v, mpfrreal_type_tag) { return mpfr::acosh(v); }
            template <typename T> inline T  asin_impl(const T& v, mpfrreal_type_tag) { return mpfr::asin (v); }
            template <typename T> inline T asinh_impl(const T& v, mpfrreal_type_tag) { return mpfr::asinh(v); }
            template <typename T> inline T  atan_impl(const T& v, mpfrreal_type_tag) { return mpfr::atan (v); }
            template <typename T> inline T atanh_impl(const T& v, mpfrreal_type_tag) { return mpfr::atanh(v); }
            template <typename T> inline T  ceil_impl(const T& v, mpfrreal_type_tag) { return mpfr::ceil (v); }
            template <typename T> inline T   cos_impl(const T& v, mpfrreal_type_tag) { return mpfr::cos  (v); }
            template <typename T> inline T  cosh_impl(const T& v, mpfrreal_type_tag) { return mpfr::cosh (v); }
            template <typename T> inline T   exp_impl(const T& v, mpfrreal_type_tag) { return mpfr::exp  (v); }
            template <typename T> inline T floor_impl(const T& v, mpfrreal_type_tag) { return mpfr::floor(v); }
            template <typename T> inline T   log_impl(const T& v, mpfrreal_type_tag) { return mpfr::log  (v); }
            template <typename T> inline T log10_impl(const T& v, mpfrreal_type_tag) { return mpfr::log10(v); }
            template <typename T> inline T  log2_impl(const T& v, mpfrreal_type_tag) { return mpfr::log2 (v); }
            template <typename T> inline T   neg_impl(const T& v, mpfrreal_type_tag) { return -v;             }
            template <typename T> inline T   pos_impl(const T& v, mpfrreal_type_tag) { return  v;             }
            template <typename T> inline T   sin_impl(const T& v, mpfrreal_type_tag) { return mpfr::sin  (v); }
            template <typename T> inline T  sinh_impl(const T& v, mpfrreal_type_tag) { return mpfr::sinh (v); }
            template <typename T> inline T  sqrt_impl(const T& v, mpfrreal_type_tag) { return mpfr::sqrt (v); }
            template <typename T> inline T   tan_impl(const T& v, mpfrreal_type_tag) { return mpfr::tan  (v); }
            template <typename T> inline T  tanh_impl(const T& v, mpfrreal_type_tag) { return mpfr::tanh (v); }
            template <typename T> inline T   cot_impl(const T& v, mpfrreal_type_tag) { return mpfr::cot  (v); }
            template <typename T> inline T   sec_impl(const T& v, mpfrreal_type_tag) { return mpfr::sec  (v); }
            template <typename T> inline T   csc_impl(const T& v, mpfrreal_type_tag) { return mpfr::csc  (v); }
            template <typename T> inline T   r2d_impl(const T& v, mpfrreal_type_tag) { return (v  * exprtk::details::constant_mpfr::_180_pi); }
            template <typename T> inline T   d2r_impl(const T& v, mpfrreal_type_tag) { return (v  * exprtk::details::constant_mpfr::pi_180 ); }
            template <typename T> inline T   d2g_impl(const T& v, mpfrreal_type_tag) { return (v  * mpfr::mpreal(20.0/9.0)); }
            template <typename T> inline T   g2d_impl(const T& v, mpfrreal_type_tag) { return (v  * mpfr::mpreal(9.0/20.0)); }
            template <typename T> inline T  notl_impl(const T& v, mpfrreal_type_tag) { return (v != mpfr::mpreal(0) ? mpfr::mpreal(0) : mpfr::mpreal(1)); }
            template <typename T> inline T  frac_impl(const T& v, mpfrreal_type_tag) { return mpfr::frac (v); }
            template <typename T> inline T trunc_impl(const T& v, mpfrreal_type_tag) { return mpfr::trunc(v); }

            template <typename T> inline T const_pi_impl(mpfrreal_type_tag) { return mpfr::const_pi   (1024, exprtk::details::constant_mpfr::mpfr_round); }
            template <typename T> inline T const_e_impl (mpfrreal_type_tag) { return mpfr::const_euler(1024, exprtk::details::constant_mpfr::mpfr_round); }

            inline bool is_true_impl (const mpfr::mpreal& v)
            {
               return v.toBool();
            }

            inline bool is_false_impl(const mpfr::mpreal& v)
            {
               return !is_true_impl(v);
            }

            template <typename T>
            inline T expm1_impl(const T& v, mpfrreal_type_tag)
            {
               return mpfr::expm1(v);
            }

            template <typename T>
            inline T min_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
                using std::min;
                return min(v0,v1);
            }

            template <typename T>
            inline T max_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               using std::max;
               return max(v0,v1);
            }

            template <typename T>
            inline T nequal_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               const T epsilon  = epsilon_type<mpfrreal_type_tag>::value();
               const T eps_norm = (mpfr::max(T(1),mpfr::max(mpfr::abs(v0),mpfr::abs(v1))) * epsilon);
               return (mpfr::abs(v0 - v1) > eps_norm) ? T(1) : T(0);
            }

            template <typename T>
            inline T sgn_impl(const T& v, mpfrreal_type_tag)
            {
                    if (v > T(0)) return T(+1);
               else if (v < T(0)) return T(-1);
               else               return T( 0);
            }

            template <typename T>
            inline T log1p_impl(const T& v, mpfrreal_type_tag)
            {
               return mpfr::log1p(v);
            }

            template <typename T>
            inline T erf_impl(const T& v, mpfrreal_type_tag)
            {
               return mpfr::erf(v);
            }

            template <typename T>
            inline T erfc_impl(const T& v, mpfrreal_type_tag)
            {
               return mpfr::erfc(v);
            }

            template <typename T>
            inline T ncdf_impl(const T& v, mpfrreal_type_tag)
            {
               T cnd = T(0.5) * (T(1) + erf_impl(
                                           mpfr::abs(v) /
                                           T(constant_mpfr::sqrt2),mpfrreal_type_tag()));
               return  (v < T(0)) ? (T(1) - cnd) : cnd;
            }

            template <typename T>
            inline T modulus_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return mpfr::fmod(v0,v1);
            }

            template <typename T>
            inline T pow_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return mpfr::pow(v0,v1);
            }

            template <typename T>
            inline T logn_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return mpfr::log(v0) / mpfr::log(v1);
            }

            template <typename T>
            inline T sinc_impl(const T& v, mpfrreal_type_tag)
            {
               if (mpfr::abs(v) >= epsilon_type<mpfrreal_type_tag>::value())
                   return(mpfr::sin(v) / v);
               else
                  return T(1);
            }

            template <typename T>
            inline T xor_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return (is_false_impl(v0) != is_false_impl(v1)) ? T(1) : T(0);
            }

            template <typename T>
            inline T xnor_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               const bool v0_true = is_true_impl(v0);
               const bool v1_true = is_true_impl(v1);
               if ((v0_true &&  v1_true) || (!v0_true && !v1_true))
                  return T(1);
               else
                  return T(0);
            }

            template <typename T>
            inline T equal_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               const T epsilon  = epsilon_type<mpfrreal_type_tag>::value();
               const T eps_norm = (mpfr::max(T(1),mpfr::max(mpfr::abs(v0),mpfr::abs(v1))) * epsilon);
               return (mpfr::abs(v0 - v1) <= eps_norm) ? T(1) : T(0);
            }

            template <typename T>
            inline T round_impl(const T& v, mpfrreal_type_tag)
            {
               return mpfr::round(v);
            }

            template <typename T>
            inline T roundn_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               const T p10 = mpfr::pow(T(10),mpfr::floor(v1));
               if (v0 < T(0))
                  return T(mpfr::ceil ((v0 * p10) - T(0.5)) / p10);
               else
                  return T(mpfr::floor((v0 * p10) + T(0.5)) / p10);
            }

            template <typename T>
            inline bool is_integer_impl(const T& v, mpfrreal_type_tag)
            {
               return mpfr::isint(v);
            }

            template <typename T>
            inline T root_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return mpfr::pow(v0,T(1) / v1);
            }

            template <typename T>
            inline T hypot_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
                return mpfr::hypot(v0,v1);
            }

            template <typename T>
            inline T atan2_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return mpfr::atan2(v0,v1);
            }

            template <typename T>
            inline T shr_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return v0 * (T(1) / mpfr::pow(T(2.0),v1));
            }

            template <typename T>
            inline T shl_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return v0 * mpfr::pow(T(2.0),v1);
            }

            template <typename T>
            inline T and_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return (is_true_impl(v0) && is_true_impl(v1)) ? T(1) : T(0);
            }

            template <typename T>
            inline T nand_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return (is_false_impl(v0) || is_false_impl(v1)) ? T(1) : T(0);
            }

            template <typename T>
            inline T or_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return (is_true_impl(v0) || is_true_impl(v1)) ? T(1) : T(0);
            }

            template <typename T>
            inline T nor_impl(const T& v0, const T& v1, mpfrreal_type_tag)
            {
               return (is_false_impl(v0) && is_false_impl(v1)) ? T(1) : T(0);
            }
         }
      }

      template <typename Iterator>
      inline bool string_to_real(Iterator& itr_external, const Iterator end, mpfr::mpreal& t, numeric::details::mpfrreal_type_tag)
      {
         t = mpfr::mpreal(std::string(itr_external,end));
         return true;
      }

      inline bool is_true (const mpfr::mpreal& v) { return details::numeric::details::is_true_impl (v); }
      inline bool is_false(const mpfr::mpreal& v) { return details::numeric::details::is_false_impl(v); }
   }

   namespace rtl { namespace io
   {
      namespace details
      {
         inline void print_type(const std::string&, const mpfr::mpreal& v, exprtk::details::numeric::details::mpfrreal_type_tag)
         {
            printf("%s",v.toString().c_str());
         }
      }
   }}
}

#endif
