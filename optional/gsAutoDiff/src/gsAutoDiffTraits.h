/** @file gsAutoDiffTraits.h

    @brief Provides traits for autodiff types
    
    This file provides Eigen support for autodiff types.
    The Eigen library already provides NumTraits specializations for autodiff types
    when including the appropriate autodiff headers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <ostream>
#include <istream>
#include <string>
#include <cmath>

// Suppress warnings from autodiff library
#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wparentheses"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

// Include autodiff headers (Eigen integration is handled elsewhere)
#include <autodiff/forward/dual.hpp>
#include <autodiff/reverse/var.hpp>

#if defined(__GNUC__) || defined(__clang__)
#pragma GCC diagnostic pop
#endif

// Compatibility traits for autodiff types

// ============================================================================
// Type detection trait for autodiff types (used in compile-time decisions)
// ============================================================================

namespace gismo {

// Type trait to detect if T is an autodiff type (for compile-time decisions)
// Base definition: always available and defaults to false for non-autodiff types
template<typename U>
struct is_autodiff_type : std::false_type {};

// Specializations for autodiff types (only when autodiff is enabled)
#if defined(gsAutoDiff_ENABLED)

// Specializations for autodiff::detail::Dual
template<typename Coeff, typename Grad>
struct is_autodiff_type<autodiff::detail::Dual<Coeff, Grad>> : std::true_type {};

// Specializations for autodiff::reverse::detail::Variable
template<typename Coeff>
struct is_autodiff_type<autodiff::reverse::detail::Variable<Coeff>> : std::true_type {};

#endif // gsAutoDiff_ENABLED

} // namespace gismo

// ============================================================================
// ExprTk compatibility for autodiff types
// ============================================================================

#ifdef GISMO_WITH_EXPRTK

// Include forward declarations
#include <gsAutoDiff/exprtk_autodiff_forward.h>

namespace exprtk {
namespace details {
namespace numeric {
namespace details {

// Type tag for autodiff::detail::Dual
template<typename T, typename G>
struct number_type<autodiff::detail::Dual<T, G>> {
    typedef gismo::detail::autodiff_type_tag type;
};

// Epsilon for autodiff types
template <>
struct epsilon_type<gismo::detail::autodiff_type_tag> {
    template<typename T>
    static inline T value() {
        return std::numeric_limits<T>::epsilon();
    }
};

// is_nan implementation
template<typename T, typename G>
inline bool is_nan_impl(const autodiff::detail::Dual<T, G>& v, gismo::detail::autodiff_type_tag) {
    return std::isnan(v.val);
}

// Conversion functions
template<typename T, typename G>
inline int to_int32_impl(const autodiff::detail::Dual<T, G>& v, gismo::detail::autodiff_type_tag) {
    return static_cast<int>(v.val);
}

template<typename T, typename G>
inline long long to_int64_impl(const autodiff::detail::Dual<T, G>& v, gismo::detail::autodiff_type_tag) {
    return static_cast<long long>(v.val);
}

template<typename T, typename G>
inline unsigned long long to_uint64_impl(const autodiff::detail::Dual<T, G>& v, gismo::detail::autodiff_type_tag) {
    return static_cast<unsigned long long>(v.val);
}

// Math functions for exprtk
template<typename T, typename G> inline autodiff::detail::Dual<T,G> abs_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return abs(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> acos_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return acos(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> acosh_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return acosh(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> asin_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return asin(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> asinh_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return asinh(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> atan_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return atan(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> atanh_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return atanh(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> ceil_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { using std::ceil; return autodiff::detail::Dual<T,G>(ceil(v.val)); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> cos_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return cos(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> cosh_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return cosh(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> exp_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return exp(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> expm1_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return expm1(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> floor_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { using std::floor; return autodiff::detail::Dual<T,G>(floor(v.val)); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> log_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return log(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> log10_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return log10(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> log1p_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return log1p(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> log2_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return log2(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> neg_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return -v; }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> pos_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return v; }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> round_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { using std::round; return autodiff::detail::Dual<T,G>(round(v.val)); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> sin_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return sin(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> sinh_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return sinh(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> sqrt_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return sqrt(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> tan_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return tan(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> tanh_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { return tanh(v); }
template<typename T, typename G> inline autodiff::detail::Dual<T,G> trunc_impl(const autodiff::detail::Dual<T,G>& v, gismo::detail::autodiff_type_tag) { using std::trunc; return autodiff::detail::Dual<T,G>(trunc(v.val)); }

// Binary operations that exprtk needs
template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> mod_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    using std::fmod; 
    return autodiff::detail::Dual<T,G>(fmod(v0.val, v1.val)); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> pow_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return pow(v0, v1); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> atan2_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return atan2(v0, v1); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> min_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return (v0.val < v1.val) ? v0 : v1; 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> max_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return (v0.val > v1.val) ? v0 : v1; 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> equal_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return autodiff::detail::Dual<T,G>((v0.val == v1.val) ? T(1) : T(0)); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> nequal_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return autodiff::detail::Dual<T,G>((v0.val != v1.val) ? T(1) : T(0)); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> hypot_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return sqrt(v0*v0 + v1*v1); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> shr_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return autodiff::detail::Dual<T,G>(static_cast<T>(static_cast<long long>(v0.val) >> static_cast<long long>(v1.val))); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> shl_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return autodiff::detail::Dual<T,G>(static_cast<T>(static_cast<long long>(v0.val) << static_cast<long long>(v1.val))); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> and_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return autodiff::detail::Dual<T,G>((v0.val != T(0) && v1.val != T(0)) ? T(1) : T(0)); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> nand_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return autodiff::detail::Dual<T,G>((v0.val != T(0) && v1.val != T(0)) ? T(0) : T(1)); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> or_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return autodiff::detail::Dual<T,G>((v0.val != T(0) || v1.val != T(0)) ? T(1) : T(0)); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> nor_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return autodiff::detail::Dual<T,G>((v0.val != T(0) || v1.val != T(0)) ? T(0) : T(1)); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> xor_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return autodiff::detail::Dual<T,G>(((v0.val != T(0)) != (v1.val != T(0))) ? T(1) : T(0)); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> xnor_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return autodiff::detail::Dual<T,G>(((v0.val != T(0)) == (v1.val != T(0))) ? T(1) : T(0)); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> root_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    return pow(v0, T(1)/v1.val); 
}

template<typename T, typename G> 
inline autodiff::detail::Dual<T,G> roundn_impl(const autodiff::detail::Dual<T,G>& v0, const autodiff::detail::Dual<T,G>& v1, gismo::detail::autodiff_type_tag) { 
    using std::pow; using std::round;
    T factor = pow(T(10), v1.val);
    return autodiff::detail::Dual<T,G>(round(v0.val * factor) / factor); 
}

} // namespace details
} // namespace numeric

} // namespace details
} // namespace exprtk

#endif // GISMO_WITH_EXPRTK

// ============================================================================
// std::to_string support for autodiff types
// ============================================================================

namespace std {

template<typename T, typename G>
inline std::string to_string(const autodiff::detail::Dual<T,G>& v) {
    return std::to_string(autodiff::val(v));
}

template<typename T>
inline std::string to_string(const autodiff::reverse::detail::Variable<T>& v) {
    return std::to_string(autodiff::val(v));
}

// std:: namespace overloads for is* functions
template<typename T, typename G>
inline bool isnan(const autodiff::detail::Dual<T,G>& v) {
    using autodiff::val;
    return std::isnan(val(v));
}

template<typename T, typename G>
inline bool isinf(const autodiff::detail::Dual<T,G>& v) {
    using autodiff::val;
    return std::isinf(val(v));
}

template<typename T, typename G>
inline bool isfinite(const autodiff::detail::Dual<T,G>& v) {
    using autodiff::val;
    return std::isfinite(val(v));
}

template<typename T>
inline bool isnan(const autodiff::reverse::detail::Variable<T>& v) {
    using autodiff::val;
    return std::isnan(val(v));
}

template<typename T>
inline bool isinf(const autodiff::reverse::detail::Variable<T>& v) {
    using autodiff::val;
    return std::isinf(val(v));
}

template<typename T>
inline bool isfinite(const autodiff::reverse::detail::Variable<T>& v) {
    using autodiff::val;
    return std::isfinite(val(v));
}

} // namespace std

// Reopen std namespace for generic autodiff expression overloads
namespace std {

// Generic overload for any autodiff expression that can be converted to a value
// This handles BinaryExpr and other expression templates
template<typename Expr>
inline auto isnan(const Expr& expr)
    -> decltype(std::isnan(autodiff::val(expr)))
{
    using autodiff::val;
    return std::isnan(val(expr));
}

template<typename Expr>
inline auto isinf(const Expr& expr)
    -> decltype(std::isinf(autodiff::val(expr)))
{
    using autodiff::val;
    return std::isinf(val(expr));
}

template<typename Expr>
inline auto isfinite(const Expr& expr)
    -> decltype(std::isfinite(autodiff::val(expr)))
{
    using autodiff::val;
    return std::isfinite(val(expr));
}

// operator<< for autodiff expression types
template<typename Op, typename ExprType>
inline std::ostream& operator<<(std::ostream& os, const autodiff::detail::UnaryExpr<Op, ExprType>& expr) {
    using autodiff::val;
    return os << val(expr);
}

template<typename Op, typename L, typename R>
inline std::ostream& operator<<(std::ostream& os, const autodiff::detail::BinaryExpr<Op, L, R>& expr) {
    using autodiff::val;
    return os << val(expr);
}

} // namespace std

// ============================================================================
// gismo::math::ceil support for autodiff expression types
// ============================================================================

namespace gismo {
namespace math {

// Overloads for autodiff::detail::Dual
template<typename T, typename G>
inline T trunc(const autodiff::detail::Dual<T,G>& v) {
    using std::trunc;
    return trunc(autodiff::val(v));
}

// Overloads for autodiff::reverse::detail::Variable
template<typename T>
inline T trunc(const autodiff::reverse::detail::Variable<T>& v) {
    using std::trunc;
    return trunc(autodiff::val(v));
}

// log10 for autodiff::reverse::detail::Variable
template<typename T>
inline autodiff::reverse::detail::ExprPtr<T> log10(const autodiff::reverse::detail::Variable<T>& v) {
    using std::log10;
    return log10(v.expr);
}

// log10 for autodiff forward-mode expression types (forward to autodiff's implementation)
template<typename E,
         typename std::enable_if<autodiff::detail::traits::isExpr<E>::value
                              && !autodiff::detail::traits::isDual<E>::value, int>::type = 0>
inline auto log10(const E& expr) -> decltype(autodiff::detail::log10(expr)) {
    using autodiff::detail::log10;
    return log10(expr);
}

// Note: floor, ceil, round for autodiff types are defined in gsAutoDiffUtils.h

// max/min for mixed Dual + expression types (forward to autodiff's implementations)
template<typename T, typename G, typename E,
         typename std::enable_if<autodiff::detail::traits::isExpr<E>::value
                              && !autodiff::detail::traits::isDual<E>::value, int>::type = 0>
inline autodiff::detail::Dual<T,G> max(const autodiff::detail::Dual<T,G>& a, const E& b) {
    using autodiff::detail::max;
    return max(a, b);
}

template<typename E, typename T, typename G,
         typename std::enable_if<autodiff::detail::traits::isExpr<E>::value
                              && !autodiff::detail::traits::isDual<E>::value, int>::type = 0>
inline autodiff::detail::Dual<T,G> max(const E& a, const autodiff::detail::Dual<T,G>& b) {
    using autodiff::detail::max;
    return max(a, b);
}

template<typename T, typename G, typename E,
         typename std::enable_if<autodiff::detail::traits::isExpr<E>::value
                              && !autodiff::detail::traits::isDual<E>::value, int>::type = 0>
inline autodiff::detail::Dual<T,G> min(const autodiff::detail::Dual<T,G>& a, const E& b) {
    using autodiff::detail::min;
    return min(a, b);
}

template<typename E, typename T, typename G,
         typename std::enable_if<autodiff::detail::traits::isExpr<E>::value
                              && !autodiff::detail::traits::isDual<E>::value, int>::type = 0>
inline autodiff::detail::Dual<T,G> min(const E& a, const autodiff::detail::Dual<T,G>& b) {
    using autodiff::detail::min;
    return min(a, b);
}

// max/min for Variable and ExprPtr
template<typename T>
inline autodiff::reverse::detail::Variable<T> max(const autodiff::reverse::detail::Variable<T>& a, const autodiff::reverse::detail::ExprPtr<T>& b) {
    using autodiff::val;
    return (val(a) > val(b)) ? a : autodiff::reverse::detail::Variable<T>(b);
}

template<typename T>
inline autodiff::reverse::detail::Variable<T> max(const autodiff::reverse::detail::ExprPtr<T>& a, const autodiff::reverse::detail::Variable<T>& b) {
    using autodiff::val;
    return (val(a) > val(b)) ? autodiff::reverse::detail::Variable<T>(a) : b;
}

template<typename T>
inline autodiff::reverse::detail::Variable<T> min(const autodiff::reverse::detail::Variable<T>& a, const autodiff::reverse::detail::ExprPtr<T>& b) {
    using autodiff::val;
    return (val(a) < val(b)) ? a : autodiff::reverse::detail::Variable<T>(b);
}

template<typename T>
inline autodiff::reverse::detail::Variable<T> min(const autodiff::reverse::detail::ExprPtr<T>& a, const autodiff::reverse::detail::Variable<T>& b) {
    using autodiff::val;
    return (val(a) < val(b)) ? autodiff::reverse::detail::Variable<T>(a) : b;
}

} // namespace math
} // namespace gismo

// Also add ceil in autodiff::detail namespace for ADL (Argument Dependent Lookup)
namespace autodiff {
namespace detail {

template<typename Op, typename ExprType>
inline auto ceil(const UnaryExpr<Op, ExprType>& expr)
    -> decltype(std::ceil(val(expr)))
{
    using std::ceil;
    return ceil(val(expr));
}

template<typename Op, typename L, typename R>
inline auto ceil(const BinaryExpr<Op, L, R>& expr)
    -> decltype(std::ceil(val(expr)))
{
    using std::ceil;
    return ceil(val(expr));
}

} // namespace detail
} // namespace autodiff

// ============================================================================
// Stream operators for autodiff types
// ============================================================================

namespace std {

// Input stream operator for autodiff::detail::Dual
template<typename T, typename G>
inline std::istream& operator>>(std::istream& is, autodiff::detail::Dual<T,G>& d) {
    T val;
    is >> val;
    d = autodiff::detail::Dual<T,G>(val);
    return is;
}

// Input stream operator for autodiff::reverse::detail::Variable
template<typename T>
inline std::istream& operator>>(std::istream& is, autodiff::reverse::detail::Variable<T>& v) {
    T val;
    is >> val;
    v = val;
    return is;
}

} // namespace std
