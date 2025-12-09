/** @file gsAutoDiffTraits.h

    @brief Traits and helpers for autodiff types

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <autodiff/forward/dual.hpp>
#include <autodiff/reverse/var.hpp>
#include <gsCore/gsForwardDeclarations.h>

namespace autodiff {
namespace detail {

// Stream operators for forward-mode expressions to allow printing
template<typename Op, typename L, typename R>
std::ostream& operator<<(std::ostream& os, const BinaryExpr<Op, L, R>& expr) {
    return os << val(expr);
}

template<typename Op, typename R>
std::ostream& operator<<(std::ostream& os, const UnaryExpr<Op, R>& expr) {
    return os << val(expr);
}

// Input stream operator for Dual
template<typename T, typename G>
std::istream& operator>>(std::istream& is, Dual<T, G>& x) {
    return is >> x.val;
}

// Math function overloads for expressions
template<typename T, typename G>
auto floor(const Dual<T, G>& x) {
    using std::floor;
    return floor(val(x));
}

template<typename Op, typename L, typename R>
auto floor(const BinaryExpr<Op, L, R>& expr) {
    using std::floor;
    return floor(val(expr));
}

template<typename Op, typename R>
auto floor(const UnaryExpr<Op, R>& expr) {
    using std::floor;
    return floor(val(expr));
}

template<typename T, typename G>
auto ceil(const Dual<T, G>& x) {
    using std::ceil;
    return ceil(val(x));
}

template<typename Op, typename L, typename R>
auto ceil(const BinaryExpr<Op, L, R>& expr) {
    using std::ceil;
    return ceil(val(expr));
}

template<typename Op, typename R>
auto ceil(const UnaryExpr<Op, R>& expr) {
    using std::ceil;
    return ceil(val(expr));
}

template<typename T, typename G>
auto round(const Dual<T, G>& x) {
    using std::round;
    return round(val(x));
}

template<typename Op, typename L, typename R>
auto round(const BinaryExpr<Op, L, R>& expr) {
    using std::round;
    return round(val(expr));
}

template<typename Op, typename R>
auto round(const UnaryExpr<Op, R>& expr) {
    using std::round;
    return round(val(expr));
}

// Check functions
template<typename T, typename G>
bool isnan(const Dual<T, G>& x) {
    using std::isnan;
    return isnan(val(x));
}

template<typename T, typename G>
bool isinf(const Dual<T, G>& x) {
    using std::isinf;
    return isinf(val(x));
}

template<typename T, typename G>
bool isfinite(const Dual<T, G>& x) {
    using std::isfinite;
    return isfinite(val(x));
}

template<typename Op, typename L, typename R>
bool isnan(const BinaryExpr<Op, L, R>& expr) {
    using std::isnan;
    return isnan(val(expr));
}

template<typename Op, typename L, typename R>
bool isinf(const BinaryExpr<Op, L, R>& expr) {
    using std::isinf;
    return isinf(val(expr));
}

template<typename Op, typename L, typename R>
bool isfinite(const BinaryExpr<Op, L, R>& expr) {
    using std::isfinite;
    return isfinite(val(expr));
}

template<typename Op, typename R>
bool isnan(const UnaryExpr<Op, R>& expr) {
    using std::isnan;
    return isnan(val(expr));
}

template<typename Op, typename R>
bool isinf(const UnaryExpr<Op, R>& expr) {
    using std::isinf;
    return isinf(val(expr));
}

template<typename Op, typename R>
bool isfinite(const UnaryExpr<Op, R>& expr) {
    using std::isfinite;
    return isfinite(val(expr));
}

} // namespace detail

   // Helper to extract the base scalar type from autodiff types
   template<class T>
   struct autodiff_real_helper {
      typedef T type;
   };

   // Specialization for Dual numbers - recursively extract the value type
   template<typename T, typename G>
   struct autodiff_real_helper<autodiff::detail::Dual<T, G>> {
      typedef typename autodiff_real_helper<T>::type type;
   };

   template<class T>
   using autodiff_real = typename autodiff_real_helper<T>::type;

} // namespace autodiff

namespace std {
    template<typename T, typename G>
    std::string to_string(const autodiff::detail::Dual<T, G>& x) {
        return std::to_string(val(x));
    }
}

// Inject into gismo::math to allow math::floor etc. to work
namespace gismo {
namespace math {
    using autodiff::detail::floor;
    using autodiff::detail::ceil;
    using autodiff::detail::round;
    using autodiff::detail::isnan;
    using autodiff::detail::isinf;
    using autodiff::detail::isfinite;

    // Add trunc support
    template<typename T, typename G>
    auto trunc(const autodiff::detail::Dual<T, G>& x) {
        using std::trunc;
        return trunc(val(x));
    }
    
    template<typename Op, typename L, typename R>
    auto trunc(const autodiff::detail::BinaryExpr<Op, L, R>& expr) {
        using std::trunc;
        return trunc(val(expr));
    }
    
    template<typename Op, typename R>
    auto trunc(const autodiff::detail::UnaryExpr<Op, R>& expr) {
        using std::trunc;
        return trunc(val(expr));
    }

    // Add log10 support
    template<typename T, typename G>
    auto log10(const autodiff::detail::Dual<T, G>& x) {
        using std::log10;
        return log10(val(x));
    }
    
    template<typename Op, typename L, typename R>
    auto log10(const autodiff::detail::BinaryExpr<Op, L, R>& expr) {
        using std::log10;
        return log10(val(expr));
    }
    
    template<typename Op, typename R>
    auto log10(const autodiff::detail::UnaryExpr<Op, R>& expr) {
        using std::log10;
        return log10(val(expr));
    }

    // Overloads for max/min with autodiff types to handle expressions
    // Use enable_if to avoid ambiguity when both arguments are Dual (handled by std::max/min)
    template<typename T, typename G, typename E,
             typename std::enable_if<!std::is_same<typename std::decay<E>::type, autodiff::detail::Dual<T, G>>::value, int>::type = 0>
    auto max(const autodiff::detail::Dual<T, G>& a, const E& b) {
        return std::max(a, autodiff::detail::Dual<T, G>(b));
    }
    
    template<typename T, typename G, typename E,
             typename std::enable_if<!std::is_same<typename std::decay<E>::type, autodiff::detail::Dual<T, G>>::value, int>::type = 0>
    auto max(const E& a, const autodiff::detail::Dual<T, G>& b) {
        return std::max(autodiff::detail::Dual<T, G>(a), b);
    }

    template<typename T, typename G, typename E,
             typename std::enable_if<!std::is_same<typename std::decay<E>::type, autodiff::detail::Dual<T, G>>::value, int>::type = 0>
    auto min(const autodiff::detail::Dual<T, G>& a, const E& b) {
        return std::min(a, autodiff::detail::Dual<T, G>(b));
    }
    
    template<typename T, typename G, typename E,
             typename std::enable_if<!std::is_same<typename std::decay<E>::type, autodiff::detail::Dual<T, G>>::value, int>::type = 0>
    auto min(const E& a, const autodiff::detail::Dual<T, G>& b) {
        return std::min(autodiff::detail::Dual<T, G>(a), b);
    }
}
}

// Avoid duplicate definitions on multiple includes
#ifndef GISMO_AUTODIFF_TO_STRING_OVERLOADS
#define GISMO_AUTODIFF_TO_STRING_OVERLOADS
namespace gismo {
namespace util {

template <class C, class = void>
struct has_val : std::false_type {};
template <class C>
struct has_val<C, std::void_t<decltype(std::declval<const C&>().val)>> : std::true_type {};

// Convert autodiff expression-like types to a `real_t` and stringify that.
// Disable for types that are already `real_t` to avoid recursion.
template <typename E, std::enable_if_t<has_val<E>::value && !std::is_same<std::decay_t<E>, real_t>::value, int> = 0>
inline std::string to_string(const E &value)
{
   return util::to_string(static_cast<real_t>(value.val));
}

// ExprPtr special-case (reverse-mode expression pointers)
template <typename T>
inline std::string to_string(const autodiff::reverse::detail::ExprPtr<T> &expr)
{
   if (expr) return util::to_string(static_cast<real_t>(expr->val));
   return std::string("null");
}

} // namespace util
} // namespace gismo
#endif
