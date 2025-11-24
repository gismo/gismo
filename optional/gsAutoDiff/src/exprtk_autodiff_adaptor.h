/** @file exprtk_autodiff_adaptor.h

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

struct AUTODIFF_TYPE_TAG {};

template<> struct number_type<AUTODIFF_TYPE> { typedef AUTODIFF_TYPE_TAG type; };

template <>
struct epsilon_type<AUTODIFF_TYPE_TAG>
{
    static inline AUTODIFF_TYPE value()
    {
        static const AUTODIFF_TYPE epsilon = AUTODIFF_TYPE(std::numeric_limits<real_t>::epsilon());
        return epsilon;
    }
};

inline bool is_nan_impl(const AUTODIFF_TYPE& v, AUTODIFF_TYPE_TAG)
{
    return std::isnan(v.val);
}

template <typename T>
inline int to_int32_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    return static_cast<int>(v.val);
}

template <typename T>
inline long long to_int64_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    return static_cast<long long int>(v.val);
}

template <typename T>
inline unsigned long long to_uint64_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    return static_cast<unsigned long long int>(v.val);
}

template <typename T> inline T   abs_impl(const T& v, AUTODIFF_TYPE_TAG) { using autodiff::detail::abs; return abs(v); }
template <typename T> inline T  acos_impl(const T& v, AUTODIFF_TYPE_TAG) { using autodiff::detail::acos; return acos(v); }
template <typename T> inline T acosh_impl(const T& v, AUTODIFF_TYPE_TAG) {
    // acosh not available for Dual, use the mathematical identity: acosh(x) = log(x + sqrt(x^2 - 1))
    using autodiff::detail::log;
    using autodiff::detail::sqrt;
    return log(v + sqrt(v * v - T(1.0)));
}
template <typename T> inline T  asin_impl(const T& v, AUTODIFF_TYPE_TAG) { using autodiff::detail::asin; return asin(v); }
template <typename T> inline T asinh_impl(const T& v, AUTODIFF_TYPE_TAG) {
    // asinh not available for Dual, use the mathematical identity: asinh(x) = log(x + sqrt(x^2 + 1))
    using autodiff::detail::log;
    using autodiff::detail::sqrt;
    return log(v + sqrt(v * v + T(1.0)));
}
template <typename T> inline T  atan_impl(const T& v, AUTODIFF_TYPE_TAG) { using autodiff::detail::atan; return atan(v); }
template <typename T> inline T atanh_impl(const T& v, AUTODIFF_TYPE_TAG) {
    // atanh not available for Dual, use the mathematical identity: atanh(x) = 0.5 * log((1+x)/(1-x))
    using autodiff::detail::log;
    return T(0.5) * log((T(1.0) + v) / (T(1.0) - v));
}
template <typename T> inline T ceil_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    typedef typename autodiff::detail::DualType<T> DualT;
    DualT evaluated(v);
    return T(DualT(std::ceil(evaluated.val)));
}
template <typename T> inline T   cos_impl(const T& v, AUTODIFF_TYPE_TAG) { return autodiff::detail::cos(v); }
template <typename T> inline T  cosh_impl(const T& v, AUTODIFF_TYPE_TAG) {
    // cosh not available for Dual, use the identity: cosh(x) = (exp(x) + exp(-x))/2
    using autodiff::detail::exp;
    return (exp(v) + exp(-v)) * T(0.5);
}
template <typename T> inline T   exp_impl(const T& v, AUTODIFF_TYPE_TAG) { return autodiff::detail::exp(v); }
template <typename T> inline T floor_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    // For any autodiff type (Dual or expression), extract scalar and return T
    typedef typename autodiff::detail::DualType<T> DualT;
    DualT evaluated(v);  // Evaluate expression if needed
    return T(DualT(std::floor(evaluated.val)));
}
template <typename T> inline T   log_impl(const T& v, AUTODIFF_TYPE_TAG) { return autodiff::detail::log(v); }
template <typename T> inline T log10_impl(const T& v, AUTODIFF_TYPE_TAG) { return autodiff::detail::log10(v); }
template <typename T> inline T  log2_impl(const T& v, AUTODIFF_TYPE_TAG) { return autodiff::detail::log(v) / autodiff::detail::log(T(2.0)); }
template <typename T> inline T   neg_impl(const T& v, AUTODIFF_TYPE_TAG) { return -v; }
template <typename T> inline T   pos_impl(const T& v, AUTODIFF_TYPE_TAG) { return v; }
template <typename T> inline T   sin_impl(const T& v, AUTODIFF_TYPE_TAG) { return autodiff::detail::sin(v); }
template <typename T> inline T  sinh_impl(const T& v, AUTODIFF_TYPE_TAG) {
    // sinh not available for Dual, use the identity: sinh(x) = (exp(x) - exp(-x))/2
    using autodiff::detail::exp;
    return (exp(v) - exp(-v)) * T(0.5);
}
template <typename T> inline T  sqrt_impl(const T& v, AUTODIFF_TYPE_TAG) { return autodiff::detail::sqrt(v); }
template <typename T> inline T   tan_impl(const T& v, AUTODIFF_TYPE_TAG) { return autodiff::detail::tan(v); }
template <typename T> inline T  tanh_impl(const T& v, AUTODIFF_TYPE_TAG) {
    // tanh not available for Dual, use the identity: tanh(x) = (exp(2x) - 1)/(exp(2x) + 1)
    using autodiff::detail::exp;
    const T e2x = exp(T(2.0) * v);
    return (e2x - T(1.0)) / (e2x + T(1.0));
}
template <typename T> inline T   cot_impl(const T& v, AUTODIFF_TYPE_TAG) { return T(1.0) / tan_impl(v, AUTODIFF_TYPE_TAG()); }
template <typename T> inline T   sec_impl(const T& v, AUTODIFF_TYPE_TAG) { return T(1.0) / cos_impl(v, AUTODIFF_TYPE_TAG()); }
template <typename T> inline T   csc_impl(const T& v, AUTODIFF_TYPE_TAG) { return T(1.0) / sin_impl(v, AUTODIFF_TYPE_TAG()); }
template <typename T> inline T   r2d_impl(const T& v, AUTODIFF_TYPE_TAG) { return (v * exprtk::details::constant_autodiff::_180_pi); }
template <typename T> inline T   d2r_impl(const T& v, AUTODIFF_TYPE_TAG) { return (v * exprtk::details::constant_autodiff::pi_180); }
template <typename T> inline T   d2g_impl(const T& v, AUTODIFF_TYPE_TAG) { return (v * T(10.0/9.0)); }
template <typename T> inline T   g2d_impl(const T& v, AUTODIFF_TYPE_TAG) { return (v * T(9.0/10.0)); }
template <typename T> inline T  notl_impl(const T& v, AUTODIFF_TYPE_TAG) {
    // notl not available for Dual, use without derivative
    return T(v.val == real_t(0) ? real_t(1) : real_t(0));
}
template <typename T> inline T  frac_impl(const T& v, AUTODIFF_TYPE_TAG) { return v - floor_impl(v, AUTODIFF_TYPE_TAG()); }
template <typename T> inline T trunc_impl(const T& v, AUTODIFF_TYPE_TAG) {
    // trunc not available for Dual, use without derivative (math namespace)
    return gismo::math::trunc(v);
}

template <typename T> inline T const_pi_impl(AUTODIFF_TYPE_TAG) { return T(exprtk::details::constant_autodiff::pi); }
template <typename T> inline T const_e_impl (AUTODIFF_TYPE_TAG) { return T(exprtk::details::constant_autodiff::e); }

template <typename T>
inline T expm1_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::exp;
    return exp(v) - T(1.0);
}

template <typename T>
inline T min_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    return v0.val < v1.val ? v0 : v1;
}

template <typename T>
inline T max_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    return v0.val > v1.val ? v0 : v1;
}

template <typename T>
inline T nequal_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    return T(v0.val != v1.val ? T(1) : T(0));
}

template <typename T>
inline T sgn_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    if (v.val > T(0)) return T(+1);
    else if (v.val < T(0)) return T(-1);
    else               return T( 0);
}

template <typename T>
inline T log1p_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::log;
    return log(T(1.0) + v);
}

template <typename T>
inline T erf_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::erf;
    return erf(v);
}

template <typename T>
inline T erfc_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::erf;
    return T(1.0) - erf(v);
}

template <typename T>
inline T ncdf_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::abs;
    using autodiff::detail::erf;
    T cnd = T(0.5) * (T(1) + erf(
                          abs(v) /
                          T(exprtk::details::constant_autodiff::sqrt2)));
    return  (v.val < T(0)) ? (T(1) - cnd) : cnd;
}

template <typename T>
inline T modulus_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    typedef typename autodiff::detail::DualType<T> DualT;
    DualT vv0(v0);
    DualT vv1(v1);
    // Evaluate the division to a Dual before accessing .val
    DualT div_result = vv0 / vv1;
    const DualT q = DualT(std::floor(div_result.val));
    return T(vv0 - vv1 * q);
}

template <typename T>
inline T pow_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::pow;
    return pow(v0, v1);
}

template <typename T>
inline T logn_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::log;
    return log(v0) / log(v1);
}

template <typename T>
inline T sinc_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::abs;
    using autodiff::detail::sin;
    if (abs(v) >= epsilon_type<AUTODIFF_TYPE_TAG>::value())
        return(sin(v) / v);
    else
        return T(1);
}

template <typename T>
inline T xor_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    return (is_false_impl(v0) != is_false_impl(v1)) ? T(1) : T(0);
}

template <typename T>
inline T xnor_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    const bool v0_true = is_true_impl(v0);
    const bool v1_true = is_true_impl(v1);
    if ((v0_true &&  v1_true) || (!v0_true && !v1_true))
        return T(1);
    else
        return T(0);
}

template <typename T>
inline T equal_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::abs;
    using autodiff::detail::max;
    const T epsilon  = epsilon_type<AUTODIFF_TYPE_TAG>::value();
    const T eps_norm = (max(T(1),max(abs(v0),abs(v1))) * epsilon);
    return (abs(T(v0-v1)) <= eps_norm) ? T(1) : T(0);
}

template <typename T>
inline T round_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    typedef typename autodiff::detail::DualType<T> DualT;
    DualT evaluated(v);
    return T(DualT(std::round(evaluated.val)));
}

template <typename T>
inline T roundn_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::pow;
    using gismo::math::floor;
    using gismo::math::ceil;
    const T p10 = pow(T(10),floor(v1));
    if (v0 < T(0))
        return T(ceil (T((v0 * p10) - T(0.5))) / p10);
    else
        return T(floor(T((v0 * p10) + T(0.5))) / p10);
}

template <typename T>
inline bool is_integer_impl(const T& v, AUTODIFF_TYPE_TAG)
{
    typedef typename autodiff::detail::DualType<T> DualT;
    DualT evaluated(v);
    return std::ceil(evaluated.val) == evaluated.val;
}

template <typename T>
inline T root_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::pow;
    return pow(v0,T(1) / v1);
}

template <typename T>
inline T hypot_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    // hypot(x,y) = sqrt(x^2 + y^2)
    using autodiff::detail::sqrt;
    return sqrt(v0 * v0 + v1 * v1);
}

template <typename T>
inline T atan2_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    // Use autodiff's atan2 if available, otherwise implement manually
    typedef typename autodiff::detail::DualType<T> DualT;
    DualT vv0(v0);
    DualT vv1(v1);

    // Compute atan2 using the identity: atan2(y,x) = atan(y/x) with proper quadrant handling
    using autodiff::detail::atan;
    using autodiff::detail::sqrt;

    // atan2(y,x) = 2*atan(y/(sqrt(x^2+y^2)+x))
    const DualT r = sqrt(vv1 * vv1 + vv0 * vv0);
    return T(DualT(2.0) * atan(vv0 / (r + vv1)));
}

template <typename T>
inline T shr_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::pow;
    return v0 * (T(1) / pow(T(2.0),v1));
}

template <typename T>
inline T shl_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    using autodiff::detail::pow;
    return v0 * pow(T(2.0),v1);
}

template <typename T>
inline T and_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    return (is_true_impl(v0) && is_true_impl(v1)) ? T(1) : T(0);
}

template <typename T>
inline T nand_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    return (is_false_impl(v0) || is_false_impl(v1)) ? T(1) : T(0);
}

template <typename T>
inline T or_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    return (is_true_impl(v0) || is_true_impl(v1)) ? T(1) : T(0);
}

template <typename T>
inline T nor_impl(const T& v0, const T& v1, AUTODIFF_TYPE_TAG)
{
    return (is_false_impl(v0) && is_false_impl(v1)) ? T(1) : T(0);
}
}
}

template <typename Iterator>
inline bool string_to_real(Iterator& itr_external,
                           const Iterator end, AUTODIFF_TYPE& t,
                           numeric::details::AUTODIFF_TYPE_TAG)
{
    const std::string num(itr_external,end);
    t = AUTODIFF_TYPE(::atof(num.c_str()));
    return true;
}

inline bool is_true (const AUTODIFF_TYPE& v)
{ return details::numeric::details::is_true_impl (v); }

inline bool is_false(const AUTODIFF_TYPE& v)
{ return details::numeric::details::is_false_impl(v); }
}

namespace helper
{
namespace details
{
inline void print_type(const std::string&, const AUTODIFF_TYPE& v,
                       exprtk::details::numeric::details::AUTODIFF_TYPE_TAG)
{
    printf("%f",v.val);
}
} // namespace details
} // namespace helper
} // namespace exprtk

#undef TYPE_TAG_
#undef TYPE_TAG
#undef AUTODIFF_TYPE_TAG
