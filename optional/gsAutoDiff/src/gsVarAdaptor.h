/**
   Adaptor for autodiff::var (reverse-mode automatic differentiation)
   
   This adaptor provides value semantics compatible with GISMO's template system.
   autodiff::var operations return ExprPtr which breaks GISMO's type consistency,
   so we wrap it in a type that converts results back to the wrapper type.

   Copyright (c) 2025 H.M. Verhelst

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.
*/

#pragma once

#include <autodiff/reverse/var.hpp>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace gismo {

// Import standard math functions to prevent name hiding
using std::abs;
using std::sqrt;
using std::pow;
using std::exp;
using std::log;
using std::log10;
using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan;
using std::atan2;
using std::sinh;
using std::cosh;
using std::tanh;
using std::erf;
using std::floor;
using std::ceil;
using std::round;
using std::trunc;
using std::max;
using std::min;
using std::isinf;
using std::isnan;
using std::isfinite;

/**
 * @brief Wrapper for autodiff::var that provides value semantics
 * 
 * This wrapper ensures that all operations return VarAdaptor instead of ExprPtr,
 * making it compatible with GISMO's template system which expects consistent types.
 */
template<typename T = double>
class VarAdaptor
{
public:
    using Variable = autodiff::reverse::detail::Variable<T>;
    using ExprPtr = autodiff::reverse::detail::ExprPtr<T>;
    
    Variable var;
    
    // Constructors
    VarAdaptor() : var(0.0) {}
    VarAdaptor(const VarAdaptor& other) : var(other.var) {}
    VarAdaptor(VarAdaptor&& other) : var(std::move(other.var)) {}
    
    template<typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
    VarAdaptor(const U& val) : var(val) {}
    
    VarAdaptor(const Variable& v) : var(v) {}
    VarAdaptor(Variable&& v) : var(std::move(v)) {}
    
    // Construct from ExprPtr (this is the key to fixing the type consistency)
    VarAdaptor(const ExprPtr& e) : var(e) {}
    
    ~VarAdaptor() {}

    // Assignment operators
    VarAdaptor& operator=(const VarAdaptor& other) { var = other.var; return *this; }
    VarAdaptor& operator=(VarAdaptor&& other) { var = std::move(other.var); return *this; }
    
    template<typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
    VarAdaptor& operator=(const U& val) { var = val; return *this; }
    
    VarAdaptor& operator=(const Variable& v) { var = v; return *this; }
    VarAdaptor& operator=(const ExprPtr& e) { var = e; return *this; }
    
    // Arithmetic assignment operators
    VarAdaptor& operator+=(const VarAdaptor& other) { var += other.var; return *this; }
    VarAdaptor& operator-=(const VarAdaptor& other) { var -= other.var; return *this; }
    VarAdaptor& operator*=(const VarAdaptor& other) { var *= other.var; return *this; }
    VarAdaptor& operator/=(const VarAdaptor& other) { var /= other.var; return *this; }
    
    template<typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
    VarAdaptor& operator+=(const U& val) { var += val; return *this; }
    
    template<typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
    VarAdaptor& operator-=(const U& val) { var -= val; return *this; }
    
    template<typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
    VarAdaptor& operator*=(const U& val) { var *= val; return *this; }
    
    template<typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
    VarAdaptor& operator/=(const U& val) { var /= val; return *this; }
    
    // Conversion operators
    explicit operator T() const { return static_cast<T>(var); }
    
    // Access to underlying variable
    Variable& get() { return var; }
    const Variable& get() const { return var; }
    
    // Value extraction
    T val() const { return autodiff::val(var); }
    
    // Update methods
    void update() { var.update(); }
    void update(T value) { var.update(value); }
};

// Type alias for double version
using VarAdaptorD = VarAdaptor<double>;

// ============================================================================
// Arithmetic operators
// ============================================================================

// Unary operators
template<typename T>
inline VarAdaptor<T> operator+(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(+x.var);
}

template<typename T>
inline VarAdaptor<T> operator-(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(-x.var);
}

// Binary operators: VarAdaptor op VarAdaptor
template<typename T>
inline VarAdaptor<T> operator+(const VarAdaptor<T>& l, const VarAdaptor<T>& r)
{
    return VarAdaptor<T>(l.var + r.var);
}

template<typename T>
inline VarAdaptor<T> operator-(const VarAdaptor<T>& l, const VarAdaptor<T>& r)
{
    return VarAdaptor<T>(l.var - r.var);
}

template<typename T>
inline VarAdaptor<T> operator*(const VarAdaptor<T>& l, const VarAdaptor<T>& r)
{
    return VarAdaptor<T>(l.var * r.var);
}

template<typename T>
inline VarAdaptor<T> operator/(const VarAdaptor<T>& l, const VarAdaptor<T>& r)
{
    return VarAdaptor<T>(l.var / r.var);
}

// Binary operators: VarAdaptor op scalar
template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline VarAdaptor<T> operator+(const VarAdaptor<T>& l, const U& r)
{
    return VarAdaptor<T>(l.var + r);
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline VarAdaptor<T> operator-(const VarAdaptor<T>& l, const U& r)
{
    return VarAdaptor<T>(l.var - r);
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline VarAdaptor<T> operator*(const VarAdaptor<T>& l, const U& r)
{
    return VarAdaptor<T>(l.var * r);
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline VarAdaptor<T> operator/(const VarAdaptor<T>& l, const U& r)
{
    return VarAdaptor<T>(l.var / r);
}

// Binary operators: scalar op VarAdaptor
template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline VarAdaptor<T> operator+(const U& l, const VarAdaptor<T>& r)
{
    return VarAdaptor<T>(l + r.var);
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline VarAdaptor<T> operator-(const U& l, const VarAdaptor<T>& r)
{
    return VarAdaptor<T>(l - r.var);
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline VarAdaptor<T> operator*(const U& l, const VarAdaptor<T>& r)
{
    return VarAdaptor<T>(l * r.var);
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline VarAdaptor<T> operator/(const U& l, const VarAdaptor<T>& r)
{
    return VarAdaptor<T>(l / r.var);
}

// ============================================================================
// Comparison operators
// ============================================================================

template<typename T>
inline bool operator==(const VarAdaptor<T>& l, const VarAdaptor<T>& r)
{
    return autodiff::val(l.var) == autodiff::val(r.var);
}

template<typename T>
inline bool operator!=(const VarAdaptor<T>& l, const VarAdaptor<T>& r)
{
    return autodiff::val(l.var) != autodiff::val(r.var);
}

template<typename T>
inline bool operator<(const VarAdaptor<T>& l, const VarAdaptor<T>& r)
{
    return autodiff::val(l.var) < autodiff::val(r.var);
}

template<typename T>
inline bool operator<=(const VarAdaptor<T>& l, const VarAdaptor<T>& r)
{
    return autodiff::val(l.var) <= autodiff::val(r.var);
}

template<typename T>
inline bool operator>(const VarAdaptor<T>& l, const VarAdaptor<T>& r)
{
    return autodiff::val(l.var) > autodiff::val(r.var);
}

template<typename T>
inline bool operator>=(const VarAdaptor<T>& l, const VarAdaptor<T>& r)
{
    return autodiff::val(l.var) >= autodiff::val(r.var);
}

// Comparisons with scalars
template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator==(const VarAdaptor<T>& l, const U& r)
{
    return autodiff::val(l.var) == r;
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator==(const U& l, const VarAdaptor<T>& r)
{
    return l == autodiff::val(r.var);
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator<(const VarAdaptor<T>& l, const U& r)
{
    return autodiff::val(l.var) < r;
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator<(const U& l, const VarAdaptor<T>& r)
{
    return l < autodiff::val(r.var);
}

// ============================================================================
// Mathematical functions
// ============================================================================

template<typename T>
inline VarAdaptor<T> abs(const VarAdaptor<T>& x)
{
    using std::abs;
    return VarAdaptor<T>(abs(x.var));
}

template<typename T>
inline VarAdaptor<T> sqrt(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(sqrt(x.var));
}

template<typename T>
inline VarAdaptor<T> pow(const VarAdaptor<T>& x, const VarAdaptor<T>& y)
{
    return VarAdaptor<T>(pow(x.var, y.var));
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline VarAdaptor<T> pow(const VarAdaptor<T>& x, const U& y)
{
    return VarAdaptor<T>(pow(x.var, y));
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline VarAdaptor<T> pow(const U& x, const VarAdaptor<T>& y)
{
    return VarAdaptor<T>(pow(x, y.var));
}

template<typename T>
inline VarAdaptor<T> exp(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(exp(x.var));
}

template<typename T>
inline VarAdaptor<T> log(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(log(x.var));
}

template<typename T>
inline VarAdaptor<T> log10(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(log10(x.var));
}

template<typename T>
inline VarAdaptor<T> sin(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(sin(x.var));
}

template<typename T>
inline VarAdaptor<T> cos(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(cos(x.var));
}

template<typename T>
inline VarAdaptor<T> tan(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(tan(x.var));
}

template<typename T>
inline VarAdaptor<T> asin(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(asin(x.var));
}

template<typename T>
inline VarAdaptor<T> acos(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(acos(x.var));
}

template<typename T>
inline VarAdaptor<T> atan(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(atan(x.var));
}

template<typename T>
inline VarAdaptor<T> atan2(const VarAdaptor<T>& y, const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(atan2(y.var, x.var));
}

template<typename T>
inline VarAdaptor<T> sinh(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(sinh(x.var));
}

template<typename T>
inline VarAdaptor<T> cosh(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(cosh(x.var));
}

template<typename T>
inline VarAdaptor<T> tanh(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(tanh(x.var));
}

template<typename T>
inline VarAdaptor<T> erf(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(erf(x.var));
}

// Floor, ceil, round return non-differentiable operations (just the value)
template<typename T>
inline VarAdaptor<T> floor(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(std::floor(autodiff::val(x.var)));
}

template<typename T>
inline VarAdaptor<T> ceil(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(std::ceil(autodiff::val(x.var)));
}

template<typename T>
inline VarAdaptor<T> round(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(std::round(autodiff::val(x.var)));
}

template<typename T>
inline VarAdaptor<T> trunc(const VarAdaptor<T>& x)
{
    return VarAdaptor<T>(std::trunc(autodiff::val(x.var)));
}

// Min/max
template<typename T>
inline VarAdaptor<T> max(const VarAdaptor<T>& a, const VarAdaptor<T>& b)
{
    return autodiff::val(a.var) > autodiff::val(b.var) ? a : b;
}

template<typename T>
inline VarAdaptor<T> min(const VarAdaptor<T>& a, const VarAdaptor<T>& b)
{
    return autodiff::val(a.var) < autodiff::val(b.var) ? a : b;
}

// Utility functions
template<typename T>
inline bool isinf(const VarAdaptor<T>& x)
{
    return std::isinf(autodiff::val(x.var));
}

template<typename T>
inline bool isnan(const VarAdaptor<T>& x)
{
    return std::isnan(autodiff::val(x.var));
}

template<typename T>
inline bool isfinite(const VarAdaptor<T>& x)
{
    return std::isfinite(autodiff::val(x.var));
}

// Stream output
template<typename T>
inline std::ostream& operator<<(std::ostream& os, const VarAdaptor<T>& v)
{
    return os << autodiff::val(v.var);
}

// Stream input
template<typename T>
inline std::istream& operator>>(std::istream& is, VarAdaptor<T>& v)
{
    T val;
    if (is >> val)
        v = VarAdaptor<T>(val);
    return is;
}

// Missing mixed comparison operators
template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator>(const VarAdaptor<T>& l, const U& r)
{
    return autodiff::val(l.var) > r;
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator>(const U& l, const VarAdaptor<T>& r)
{
    return l > autodiff::val(r.var);
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator<=(const VarAdaptor<T>& l, const U& r)
{
    return autodiff::val(l.var) <= r;
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator<=(const U& l, const VarAdaptor<T>& r)
{
    return l <= autodiff::val(r.var);
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator>=(const VarAdaptor<T>& l, const U& r)
{
    return autodiff::val(l.var) >= r;
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator>=(const U& l, const VarAdaptor<T>& r)
{
    return l >= autodiff::val(r.var);
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator!=(const VarAdaptor<T>& l, const U& r)
{
    return autodiff::val(l.var) != r;
}

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic<U>::value>>
inline bool operator!=(const U& l, const VarAdaptor<T>& r)
{
    return l != autodiff::val(r.var);
}

namespace math {
    using gismo::abs;
    using gismo::sqrt;
    using gismo::pow;
    using gismo::exp;
    using gismo::log;
    using gismo::log10;
    using gismo::sin;
    using gismo::cos;
    using gismo::tan;
    using gismo::asin;
    using gismo::acos;
    using gismo::atan;
    using gismo::atan2;
    using gismo::sinh;
    using gismo::cosh;
    using gismo::tanh;
    using gismo::erf;
    using gismo::floor;
    using gismo::ceil;
    using gismo::round;
    using gismo::trunc;
    using gismo::max;
    using gismo::min;
    using gismo::isinf;
    using gismo::isnan;
    using gismo::isfinite;
}

} // namespace gismo

// Add to std namespace for compatibility
namespace std {

template<typename T>
inline string to_string(const gismo::VarAdaptor<T>& v)
{
    std::ostringstream oss;
    oss << autodiff::val(v.var);
    return oss.str();
}

} // namespace std

// Add to autodiff namespace for val() extraction
namespace autodiff {

template<typename T>
inline T val(const gismo::VarAdaptor<T>& x)
{
    return val(x.var);
}

template<typename T>
inline T val(gismo::VarAdaptor<T>& x)
{
    return val(x.var);
}

} // namespace autodiff
