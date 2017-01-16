/**
   Automatic differentiation data type for C++, depends on the Eigen
   linear algebra library.

   Copyright (c) 2012 by Wenzel Jakob. Based on code by Jon Kaldor
   and Eitan Grinspun.

   Modifications for G+Smo, Angelos Mantzaflaris, 2015

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.
*/

#pragma once

//#include <cmath>
//#define EIGEN_DONT_PARALLELIZE 
//#define EIGEN_NO_DEBUG

namespace gismo
{

namespace ad
{

/**
 * \brief Automatic differentiation scalar with first-order derivatives
 *
 * This class provides an instrumented "scalar" value, which may be dependent on
 * a number of independent variables. The implementation keeps tracks of
 * first -order derivatives with respect to these variables using a set
 * of overloaded operations and implementations of special functions (sin,
 * tan, exp, ..).
 *
 * This is extremely useful for numerical zero-finding, particularly when
 * analytic derivatives from programs like Maple or Mathematica suffer from
 * excessively complicated expressions.
 *
 * The class relies on templates, which makes it possible to fix the
 * number of independent variables at compile-time so that instances can
 * be allocated on the stack. Otherwise, they will be placed on the heap.
 *
 * This is an extended C++ port of Jon Kaldor's implementation, which is
 * based on a C version by Eitan Grinspun at Caltech)
 *
 * \sa DScalar2
 * \author Wenzel Jakob
 */
template <typename _Scalar, int d = -1> 
struct DScalar1
{
public:
    typedef _Scalar                    Scalar;
    typedef Eigen::Matrix<_Scalar,d,1> Gradient_t;
    
    // ======================================================================
    /// @{ \name Constructors and accessors
    // ======================================================================

    /// Create a new constant automatic differentiation scalar
    explicit DScalar1(const Scalar _value = Scalar(0.0) ) 
    : value(_value) 
    {
        // Note: number of variables might be still unknown if d==-1
        // We will recover them during a binary operation
        grad.setZero();
    }

    /// Construct a new scalar with the specified value and one first derivative set to 1
    DScalar1(size_t index, const Scalar _value)
    : value(_value) 
    {
        // Note: number of variables is expected to be d > 0
        assert( d != -1 );
        assert( index < d && "Index must be less than the number of variables");
        grad.setZero();
        grad(index) = 1;
    }

    /// Construct a new scalar with the specified value and one first
    /// derivative set to 1, with \a numVars variables
    DScalar1(size_t index, size_t numVars, const Scalar _value)
    {
        setVariable(index,numVars,_value);
    }

    /// Construct a scalar associated with the given gradient
    DScalar1(const Scalar _value, const Gradient_t & grad)
    : value(_value), grad(grad) 
    { }

    /// Copy constructor
    DScalar1(const DScalar1 & s)
    : value(s.value), grad(s.grad) 
    { }

    inline const Scalar     & getValue()    const { return value; }

    inline const Gradient_t & getGradient() const 
    { 
        GISMO_ASSERT(0!=grad.size(), "Gradient is empty (use gradient_into), value= "<< value);
        return grad; 
    }

    template <typename Derived>
    inline void gradient_into(const Eigen::DenseBase<Derived> & res) const
    {
        if ( 0==grad.size() )
            grad.setZero( res.rows() );

        // Note: Eigen hack to write on expression
        const_cast<Eigen::DenseBase<Derived>&>(res) = grad;
    }


    inline size_t             numVars()     const { return grad.size(); }

    inline void setVariable(size_t index, size_t numVars, const Scalar & _value)
    {
        value = _value;
        assert( d == -1 || d == numVars );
        assert( index < numVars && "Index must be less than the number of variables");
        grad.setZero(numVars);
        grad(index) = 1;
    }

    template<int _Rows, int _Cols>
    static void Initialize(const Eigen::Matrix<Scalar,_Rows,_Cols> & values,
                           Eigen::Matrix<DScalar1,_Rows,_Cols> & result)
    {
        result.resize(values.rows(), values.cols());
        const int numVars = values.size();
        DScalar1 * data = result.data();
        for ( int i = 0; i!= numVars; ++i, ++data)
            data->setVariable(i, numVars, values[i]);
    }

    // ======================================================================
    /// @{ \name Addition
    // ======================================================================
    friend DScalar1 operator+(const DScalar1 &lhs, const DScalar1 &rhs) 
    {
        prepare(lhs,rhs);
        return DScalar1(lhs.value+rhs.value, lhs.grad+rhs.grad);
    }

    friend DScalar1 operator+(const DScalar1 &lhs, const Scalar &rhs) 
    {
        return DScalar1(lhs.value+rhs, lhs.grad);
    }

    friend DScalar1 operator+(const Scalar &lhs, const DScalar1 &rhs) 
    {
        return DScalar1(rhs.value+lhs, rhs.grad);
    }

    inline DScalar1& operator+=(const DScalar1 &s) 
    {
        value += s.value;
        grad += s.grad;
        return *this;
    }

    inline DScalar1& operator+=(const Scalar &v) 
    {
        value += v;
        return *this;
    }

    /// @}
    // ======================================================================

    // ======================================================================
    /// @{ \name Subtraction
    // ======================================================================

    friend DScalar1 operator-(const DScalar1 &lhs, const DScalar1 &rhs) 
    {
        prepare(lhs,rhs);
        return DScalar1(lhs.value-rhs.value, lhs.grad-rhs.grad);
    }

    friend DScalar1 operator-(const DScalar1 &lhs, const Scalar &rhs) 
    {
        return DScalar1(lhs.value-rhs, lhs.grad);
    }

    friend DScalar1 operator-(const Scalar &lhs, const DScalar1 &rhs) 
    {
        return DScalar1(lhs-rhs.value, -rhs.grad);
    }

    friend DScalar1 operator-(const DScalar1 &s) 
    {
        return DScalar1(-s.value, -s.grad);
    }

    inline DScalar1& operator-=(const DScalar1 &s) 
    {
        assert( grad.size() == s.grad.size() && 
                "You mixed diff-scalars with different number of variables." );
        value -= s.value;
        grad  -= s.grad;
        return *this;
    }

    inline DScalar1& operator-=(const Scalar &v) 
    {
        value -= v;
        return *this;
    }
    /// @}
    // ======================================================================

    // ======================================================================
    /// @{ \name Division
    // ======================================================================
    friend DScalar1 operator/(const DScalar1 &lhs, const Scalar &rhs) 
    {
        if (rhs == 0)
            throw std::runtime_error("DScalar1: Division by zero!");
        Scalar inv = 1.0f / rhs;
        return DScalar1(lhs.value*inv, lhs.grad*inv);
    }

    friend DScalar1 operator/(const Scalar &lhs, const DScalar1 &rhs) 
    {
        return lhs * inverse(rhs);
    }

    friend DScalar1 operator/(const DScalar1 &lhs, const DScalar1 &rhs) 
    {
        prepare(lhs,rhs);
        return lhs * inverse(rhs);
    }

    friend DScalar1 inverse(const DScalar1 &s) 
    {
        Scalar valueSqr = s.value*s.value,
            invValueSqr = (Scalar) 1 / valueSqr;

        // vn = 1/v, Dvn = -1/(v^2) Dv
        return DScalar1((Scalar) 1 / s.value, s.grad * -invValueSqr);
    }

    inline DScalar1& operator/=(const Scalar &v) 
    {
        value /= v;
        grad /= v;
        return *this;
    }

    inline DScalar1& operator/=(const DScalar1 &v) 
    {
        *this = (*this) / v ;
        return *this;
    }

    /// @}
    // ======================================================================

    // ======================================================================
    /// @{ \name Multiplication
    // ======================================================================
    inline friend DScalar1 operator*(const DScalar1 &lhs, const Scalar &rhs) 
    {
        return DScalar1(lhs.value*rhs, lhs.grad*rhs);
    }

    inline friend DScalar1 operator*(const Scalar &lhs, const DScalar1 &rhs) 
    {
        return DScalar1(rhs.value*lhs, rhs.grad*lhs);
    }

    inline friend DScalar1 operator*(const DScalar1 &lhs, const DScalar1 &rhs) 
    {
        prepare(lhs,rhs);

        // Product rule
        return DScalar1(lhs.value*rhs.value,
                        rhs.grad * lhs.value + lhs.grad * rhs.value);
    }

    inline DScalar1& operator*=(const Scalar &v) 
    {
        value *= v;
        grad *= v;
        return *this;
    }

    inline DScalar1& operator*=(const DScalar1 &v) 
    {
        *this = *this * v ;
        return (*this);
    }

    /// @}
    // ======================================================================

    // ======================================================================
    /// @{ \name Miscellaneous functions
    // ======================================================================

    friend bool isnan(const DScalar1 &s) {return std::isnan(s.value); }

    friend DScalar1 abs(const DScalar1 &s) 
    {
        DScalar1 result(math::abs(s.value));
        result.grad = (s.value < 0 ? -1 : 1) * s.grad;
        return result;
    }

    friend DScalar1 sqrt(const DScalar1 &s) 
    {
        Scalar sqrtVal = math::sqrt(s.value),
            temp    = (Scalar) 1 / ((Scalar) 2 * sqrtVal);

        // vn = sqrt(v)
        // Dvn = 1/(2 sqrt(v)) Dv
        return DScalar1(sqrtVal, s.grad * temp);
    }

    friend DScalar1 pow(const DScalar1 &s, const Scalar &a)
    {
        Scalar powVal = math::pow(s.value, a),
            temp   = a * math::pow(s.value, a-1);
        // vn = v ^ a, Dvn = a*v^(a-1) * Dv
        return DScalar1(powVal, s.grad * temp);
    }

    friend DScalar1 pow(const DScalar1 &s, const DScalar1 &a)
    {
        gsWarn<<"pow not implemented.\n";
        return s;
    }

    friend DScalar1 exp(const DScalar1 &s) 
    {
        Scalar expVal = math::exp(s.value);

        // vn = exp(v), Dvn = exp(v) * Dv
        return DScalar1(expVal, s.grad * expVal);
    }

    friend DScalar1 log(const DScalar1 &s) 
    {
        Scalar logVal = math::log(s.value);

        // vn = log(v), Dvn = Dv / v
        return DScalar1(logVal, s.grad / s.value);
    }

    friend DScalar1 sin(const DScalar1 &s) 
    {
        // vn = sin(v), Dvn = cos(v) * Dv
        return DScalar1(math::sin(s.value), s.grad * math::cos(s.value));
    }

    friend DScalar1 cos(const DScalar1 &s) 
    {
        // vn = cos(v), Dvn = -sin(v) * Dv
        return DScalar1(math::cos(s.value), s.grad * -math::sin(s.value));
    }

    friend DScalar1 acos(const DScalar1 &s) 
    {
        if (math::abs(s.value) >= 1)
            throw std::runtime_error("acos: Expected a value in (-1, 1)");

        Scalar temp = -math::sqrt((Scalar) 1 - s.value*s.value);

        // vn = acos(v), Dvn = -1/sqrt(1-v^2) * Dv
        return DScalar1(math::acos(s.value),
                        s.grad * ((Scalar) 1 / temp));
    }

    friend DScalar1 asin(const DScalar1 &s) 
    {
        if (math::abs(s.value) >= 1)
            throw std::runtime_error("asin: Expected a value in (-1, 1)");

        Scalar temp = math::sqrt((Scalar) 1 - s.value*s.value);

        // vn = asin(v), Dvn = 1/sqrt(1-v^2) * Dv
        return DScalar1(math::asin(s.value),
                        s.grad * ((Scalar) 1 / temp));
    }

    friend DScalar1 atan2(const DScalar1 &y, const DScalar1 &x) 
    {
        prepare(y,x);

        const Scalar denom = x.value*x.value + y.value*y.value;
        
        // vn = atan2(y, x), Dvn = (x*Dy - y*Dx) / (x^2 + y^2)
        return DScalar1(math::atan2(y.value, x.value),
                        y.grad * (x.value / denom) - x.grad * (y.value / denom));
    }

    friend DScalar1 atan(const DScalar1 &y) 
    {
        // vn = atan(y), Dvn = Dy / (1 + y^2)
        return DScalar1(math::atan(y.value), y.grad / ((Scalar)1.0 + y.value*y.value) );
    }

    friend DScalar1 tanh(const DScalar1 &y) 
    {
        const Scalar th = math::tanh(y.value);
        // vn = tanh(y)
        // Dvn = (1 - tanh(y)^2)*Dy
        return DScalar2(th, (1-th*th) * y.grad );
    }

    friend DScalar1 cosh(const DScalar1 &y) 
    {
        const Scalar ch = math::cosh(y.value);
        const Scalar sh = math::sinh(y.value);
        // vn = cosh(y)
        // Dvn = sinh(y)*Dy
        return DScalar1(ch, sh * y.grad );
    }

    friend DScalar1 sinh(const DScalar1 &y) 
    {
        const Scalar ch = math::cosh(y.value);
        const Scalar sh = math::sinh(y.value);
        // vn = cosh(y)
        // Dvn = sinh(y)*Dy
        return DScalar1(sh, ch * y.grad );
    }


    /// @}
    // ======================================================================

    // ======================================================================
    /// @{ \name Comparison and assignment
    // ======================================================================

    inline void operator= (const DScalar1& s) { value = s.value; grad = s.grad; }
    inline void operator= (const Scalar &v  ) { value = v; grad.setZero(); }
    inline bool operator< (const DScalar1& s) const { return value < s.value; }
    inline bool operator< (const Scalar& s  ) const { return value < s; }
    inline bool operator<=(const DScalar1& s) const { return value <= s.value; }
    inline bool operator<=(const Scalar& s  ) const { return value <= s; }
    inline bool operator>=(const DScalar1& s) const { return value >= s.value; }
    inline bool operator>=(const Scalar& s  ) const { return value >= s; }
    inline bool operator> (const DScalar1& s) const { return value > s.value; }
    inline bool operator> (const Scalar& s  ) const { return value > s; }
    inline bool operator==(const DScalar1& s) const { return value == s.value; }
    inline bool operator==(const Scalar& s  ) const { return value == s; }
    inline bool operator!=(const DScalar1& s) const { return value != s.value; }
    inline bool operator!=(const Scalar& s  ) const { return value != s; }

    friend bool operator<(const Scalar s1, const DScalar1  & s2) { return  s1 <  s2.value; }
    friend bool operator>(const Scalar s1, const DScalar1  & s2) { return  s1 >  s2.value; }
    friend bool operator>=(const Scalar s1, const DScalar1 & s2) { return  s1 >= s2.value; }
    friend bool operator<=(const Scalar s1, const DScalar1 & s2) { return  s1 <= s2.value; }
    friend bool operator==(const Scalar s1, const DScalar1 & s2) { return  s1 == s2.value; }
    friend bool operator!=(const Scalar s1, const DScalar1 & s2) { return  s1 != s2.value; }

    /// @}
    // ======================================================================

protected:

    // prepare two scalars for a binary operation
    static inline void prepare(const DScalar1 &lhs, const DScalar1 &rhs) 
    {
        if ( ! lhs.grad.size() )
        {
            const size_t nv = rhs.grad.size();
            lhs.grad.setZero(nv    );
            return;
        }
        if ( ! rhs.grad.size() )
        {
            const size_t nv = lhs.grad.size();
            rhs.grad.setZero(nv    );
            return;
        }
    }

    // for d>0
    // static inline void prepare<Scalar,-1>(const DScalar1 &lhs, const DScalar1 &rhs) 
    // { }

protected:

    Scalar value;

    // mutable for treating incomplete constants
    mutable Gradient_t grad;
};

template <typename Scalar, int d>
std::ostream &operator<<(std::ostream &out, const DScalar1<Scalar,d> &s) 
{
	out << "[" << s.getValue()
		<< ", grad=" << s.getGradient().format(Eigen::IOFormat(4, 1, ", ", "; ", "", "", "[", "]"))
		<< "]";
	return out;
}

/**
 * \brief Automatic differentiation scalar with first- and second-order derivatives
 *
 * This class provides an instrumented "scalar" value, which may be dependent on
 * a number of independent variables. The implementation keeps tracks of first
 * and second-order drivatives with respect to these variables using a set
 * of overloaded operations and implementations of special functions (sin,
 * tan, exp, ..).
 *
 * This is extremely useful for numerical optimization, particularly when
 * analytic derivatives from programs like Maple or Mathematica suffer from
 * excessively complicated expressions.
 *
 * The class relies on templates, which makes it possible to fix the
 * number of independent variables at compile-time so that instances can
 * be allocated on the stack. Otherwise, they will be placed on the heap.
 *
 * This is an extended C++ port of Jon Kaldor's implementation, which is
 * based on a C version by Eitan Grinspun at Caltech)
 *
 * \sa DScalar1
 * \author Wenzel Jakob
 */
template <typename _Scalar, int d = -1> // todo: nder = 1, 2
struct DScalar2
{
public:
    typedef _Scalar                    Scalar;
    typedef Eigen::Matrix<_Scalar,d,1> Gradient_t;
    typedef Eigen::Matrix<_Scalar,d,d> Hessian_t;
    
    // ======================================================================
    /// @{ \name Constructors and accessors
    // ======================================================================

    /// Create a new constant automatic differentiation scalar
    explicit DScalar2(Scalar value = (Scalar) 0) : value(value) 
    {
        // Note: number of variables might be still unknown if d==-1
        // We will recover them during a binary operation
        grad.setZero();
        hess.setZero();
    }

    DScalar2(size_t index, const Scalar & _value)
    : value(_value) 
    {
        // Note: number of variables is expected to be d > 0
        assert( d != -1 );
        assert( index < d && "Index must be less than the number of variables");
        grad.setZero();
        grad(index) = 1;
        hess.setZero();
    }

    /// Construct a new scalar with the specified value and one first
    /// derivative set to 1, with \a numVars variables
    DScalar2(size_t index, size_t numVars, const Scalar & _value)
    {
        setVariable(index,numVars,_value);
    }

    /// Construct a scalar associated with the given gradient and Hessian
    DScalar2(Scalar _value, const Gradient_t &grad, const Hessian_t &hess)
    : value(_value), grad(grad), hess(hess) { }

    /// Copy constructor
    DScalar2(const DScalar2 & s)
    : value(s.value), grad(s.grad), hess(s.hess) { }

    inline const Scalar     & getValue()    const { return value; }

    inline const Gradient_t & getGradient() const
    {
        GISMO_ASSERT(0!=grad.size(), "Gradient is empty (use gradient_into), value= "<< value);
        return grad;
    }

    template <typename Derived>
    inline void gradient_into(const Eigen::DenseBase<Derived> & res) const
    {
        if ( 0==grad.size() )
            grad.setZero( res.rows() );

        // Note: Eigen hack to write on expression
        const_cast<Eigen::DenseBase<Derived>&>(res) = grad;
    }

    inline const Hessian_t  & getHessian()  const 
    { 
        GISMO_ASSERT(0!=grad.size(), "Hessian is empty (use hessian_into), value= "<< value);
        return hess; 
    }

    template <typename Derived>
    inline void hessian_into(const Eigen::DenseBase<Derived> & res) const
    {
        if ( 0==hess.size() )
            hess.setZero( res.rows(), res.cols() );

        // Note: Eigen hack to write on expression
        const_cast<Eigen::DenseBase<Derived>&>(res) = hess;
    }

    inline size_t             numVars()     const { return grad.size(); }

    inline void setVariable(size_t index, size_t numVars, const Scalar & _value)
    {
        value = _value;
        assert( d == -1 || d == static_cast<int>(numVars) );
        assert( index < numVars && "Index must be less than the number of variables");
        grad.setZero(numVars);
        grad(index) = 1;
        hess.setZero(numVars,numVars);
    }

    template<int _Rows, int _Cols>
    static void Initialize(const Eigen::Matrix<Scalar,_Rows,_Cols> & values,
                           Eigen::Matrix<DScalar2,_Rows,_Cols> & result)
    {
        result.resize(values.rows(), values.cols());
        const int numVars = values.size();
        DScalar2 * data = result.data();
        for ( int i = 0; i!= numVars; ++i, ++data)
            data->setVariable(i, numVars, values[i]);
    }

    // ======================================================================
    /// @{ \name Addition
    // ======================================================================
    friend DScalar2 operator+(const DScalar2 &lhs, const DScalar2 &rhs) 
    {
        prepare(lhs,rhs);
        return DScalar2(lhs.value+rhs.value,
                        lhs.grad+rhs.grad, lhs.hess+rhs.hess);
    }

    friend DScalar2 operator+(const DScalar2 &lhs, const Scalar &rhs) 
    {
        return DScalar2(lhs.value+rhs, lhs.grad, lhs.hess);
    }

    friend DScalar2 operator+(const Scalar &lhs, const DScalar2 &rhs) 
    {
        return DScalar2(rhs.value+lhs, rhs.grad, rhs.hess);
    }

    inline DScalar2& operator+=(const DScalar2 &s) 
    {
        prepare(*this,s);
        value += s.value;
        grad  += s.grad;
        hess  += s.hess;
        return *this;
    }

    inline DScalar2& operator+=(const Scalar &v) 
    {
        value += v;
        return *this;
    }

    /// @}
    // ======================================================================

    // ======================================================================
    /// @{ \name Subtraction
    // ======================================================================

    friend DScalar2 operator-(const DScalar2 &lhs, const DScalar2 &rhs) 
    {
        prepare(lhs,rhs);
        return DScalar2(lhs.value-rhs.value, lhs.grad-rhs.grad, lhs.hess-rhs.hess);
    }

    friend DScalar2 operator-(const DScalar2 &lhs, const Scalar &rhs) 
    {
        return DScalar2(lhs.value-rhs, lhs.grad, lhs.hess);
    }

    friend DScalar2 operator-(const Scalar &lhs, const DScalar2 &rhs) 
    {
        return DScalar2(lhs-rhs.value, -rhs.grad, -rhs.hess);
    }

    friend DScalar2 operator-(const DScalar2 &s) 
    {
        return DScalar2(-s.value, -s.grad, -s.hess);
    }

    inline DScalar2& operator-=(const DScalar2 &s) 
    {
        value -= s.value;
        grad -= s.grad;
        hess -= s.hess;
        return *this;
    }

    inline DScalar2& operator-=(const Scalar &v) 
    {
        value -= v;
        return *this;
    }
    /// @}
    // ======================================================================

    // ======================================================================
    /// @{ \name Division
    // ======================================================================
    friend DScalar2 operator/(const DScalar2 &lhs, const Scalar &rhs) 
    {
        if (rhs == 0)
            throw std::runtime_error("DScalar2: Division by zero!");
        Scalar inv = 1.0f / rhs;
        return DScalar2(lhs.value*inv, lhs.grad*inv, lhs.hess*inv);
    }

    friend DScalar2 operator/(const Scalar &lhs, const DScalar2 &rhs) 
    {
        return lhs * inverse(rhs);
    }

    friend DScalar2 operator/(const DScalar2 &lhs, const DScalar2 &rhs) 
    {
        prepare(lhs,rhs);
        return lhs * inverse(rhs);
    }

    friend DScalar2 inverse(const DScalar2 &s) 
    {
        Scalar valueSqr = s.value*s.value,
            valueCub = valueSqr * s.value,
            invValueSqr = (Scalar) 1 / valueSqr;

        // vn = 1/v
        DScalar2 result((Scalar) 1 / s.value);

        // Dvn = -1/(v^2) Dv
        result.grad = s.grad * -invValueSqr;

        // D^2vn = -1/(v^2) D^2v + 2/(v^3) Dv Dv^T
        result.hess = s.hess * -invValueSqr;
        result.hess += s.grad * s.grad.transpose()
            * ((Scalar) 2 / valueCub);

        return result;
    }

    inline DScalar2& operator/=(const Scalar &v) 
    {
        value /= v;
        grad /= v;
        hess /= v;
        return *this;
    }

    inline DScalar2& operator/=(const DScalar2 &v) 
    {
        *this = (*this) / v ;
        return *this;
    }

    /// @}
    // ======================================================================

    // ======================================================================
    /// @{ \name Multiplication
    // ======================================================================
    friend DScalar2 operator*(const DScalar2 &lhs, const Scalar &rhs) 
    {
        return DScalar2(lhs.value*rhs, lhs.grad*rhs, lhs.hess*rhs);
    }

    friend DScalar2 operator*(const Scalar &lhs, const DScalar2 &rhs) 
    {
        return DScalar2(rhs.value*lhs, rhs.grad*lhs, rhs.hess*lhs);
    }

    friend DScalar2 operator*(const DScalar2 &lhs, const DScalar2 &rhs) 
    {
        prepare(lhs,rhs);

        DScalar2 result(lhs.value*rhs.value);

        /// Product rule
        result.grad = rhs.grad * lhs.value + lhs.grad * rhs.value;

        // (i,j) = g*F_xixj + g*G_xixj + F_xi*G_xj + F_xj*G_xi
        result.hess = rhs.hess * lhs.value;
        result.hess += lhs.hess * rhs.value;
        result.hess += lhs.grad * rhs.grad.transpose();
        result.hess += rhs.grad * lhs.grad.transpose();

        return result;
    }

    inline DScalar2& operator*=(const Scalar &v) 
    {
        value *= v;
        grad  *= v;
        hess  *= v;
        return *this;
    }

    inline DScalar2& operator*=(const DScalar2 &v) 
    {
        *this = (*this) * v ;
        return *this;
    }

    /// @}
    // ======================================================================

    // ======================================================================
    /// @{ \name Miscellaneous functions
    // ======================================================================

    friend bool isnan(const DScalar2 &s) {return std::isnan(s.value); }

    friend DScalar2 abs(const DScalar2 &s) 
    {
        DScalar2 result(math::abs(s.value));
        result.grad = (s.value < 0 ? -1 : 1) * s.grad;
        result.hess.setZero(s.hess.rows(), s.hess.cols());
        return result;
    }

    friend DScalar2 sqrt(const DScalar2 &s) 
    {
        Scalar sqrtVal = std::sqrt(s.value),
            temp    = (Scalar) 1 / ((Scalar) 2 * sqrtVal);

        // vn = sqrt(v)
        DScalar2 result(sqrtVal);

        // Dvn = 1/(2 sqrt(v)) Dv
        result.grad = s.grad * temp;

        // D^2vn = 1/(2 sqrt(v)) D^2v - 1/(4 v*sqrt(v)) Dv Dv^T
        result.hess = s.hess * temp;
        result.hess += s.grad * s.grad.transpose()
            * (-(Scalar) 1 / ((Scalar) 4 * s.value * sqrtVal));

        return result;
    }

    friend DScalar2 pow(const DScalar2 &s, const Scalar &a) 
    {
        Scalar powVal = std::pow(s.value, a),
            temp   = a * std::pow(s.value, a-1);
        // vn = v ^ a
        DScalar2 result(powVal);

        // Dvn = a*v^(a-1) * Dv
        result.grad = s.grad * temp;

        // D^2vn = a*v^(a-1) D^2v - 1/(4 v*sqrt(v)) Dv Dv^T
        result.hess = s.hess * temp;
        result.hess += s.grad * s.grad.transpose()
            * (a * (a-1) * std::pow(s.value, a-2));

        return result;
    }

    friend DScalar2 exp(const DScalar2 &s) 
    {
        Scalar expVal = std::exp(s.value);

        // vn = exp(v)
        DScalar2 result(expVal);

        // Dvn = exp(v) * Dv
        result.grad = s.grad * expVal;

        // D^2vn = exp(v) * Dv*Dv^T + exp(v) * D^2v
        result.hess = (s.grad * s.grad.transpose()
                       + s.hess) * expVal;

        return result;
    }

    friend DScalar2 log(const DScalar2 &s) 
    {
        Scalar logVal = math::log(s.value);

        // vn = log(v)
        DScalar2 result(logVal);

        // Dvn = Dv / v
        result.grad = s.grad / s.value;

        // D^2vn = (v*D^2v - Dv*Dv^T)/(v^2)
        result.hess = s.hess / s.value -
            (s.grad * s.grad.transpose() / (s.value*s.value));

        return result;
    }

    friend DScalar2 sin(const DScalar2 &s) 
    {
        Scalar sinVal = math::sin(s.value),
            cosVal = math::cos(s.value);

        // vn = sin(v)
        DScalar2 result(sinVal);

        // Dvn = cos(v) * Dv
        result.grad = s.grad * cosVal;

        // D^2vn = -sin(v) * Dv*Dv^T + cos(v) * Dv^2
        result.hess = s.hess * cosVal;
        result.hess += s.grad * s.grad.transpose() * -sinVal;

        return result;
    }

    friend DScalar2 cos(const DScalar2 &s) 
    {
        Scalar sinVal = math::sin(s.value),
            cosVal = math::cos(s.value);
        // vn = cos(v)
        DScalar2 result(cosVal);

        // Dvn = -sin(v) * Dv
        result.grad = s.grad * -sinVal;

        // D^2vn = -cos(v) * Dv*Dv^T - sin(v) * Dv^2
        result.hess = s.hess * -sinVal;
        result.hess += s.grad * s.grad.transpose() * -cosVal;

        return result;
    }

    friend DScalar2 acos(const DScalar2 &s) 
    {
        if (math::abs(s.value) >= 1.0)
            throw std::runtime_error("acos: Expected a value in (-1, 1)");

        Scalar temp = -math::sqrt((Scalar) 1 - s.value*s.value);

        // vn = acos(v)
        DScalar2 result(math::acos(s.value));

        // Dvn = -1/sqrt(1-v^2) * Dv
        result.grad = s.grad * ((Scalar) 1 / temp);

        // D^2vn = -1/sqrt(1-v^2) * D^2v - v/[(1-v^2)^(3/2)] * Dv*Dv^T
        result.hess = s.hess * ((Scalar) 1 / temp);
        result.hess += s.grad * s.grad.transpose()
            * s.value / (temp*temp*temp);

        return result;
    }

    friend DScalar2 asin(const DScalar2 &s) 
    {
        if (math::abs(s.value) >= 1.0)
            throw std::runtime_error("asin: Expected a value in (-1, 1)");

        Scalar temp = math::sqrt((Scalar) 1 - s.value*s.value);

        // vn = asin(v)
        DScalar2 result(math::asin(s.value));

        // Dvn = 1/sqrt(1-v^2) * Dv
        result.grad = s.grad * ((Scalar) 1 / temp);

        // D^2vn = 1/sqrt(1-v*v) * D^2v + v/[(1-v^2)^(3/2)] * Dv*Dv^T
        result.hess = s.hess * ((Scalar) 1 / temp);
        result.hess += s.grad * s.grad.transpose()
            * s.value / (temp*temp*temp);

        return result;
    }

    friend DScalar2 atan2(const DScalar2 &y, const DScalar2 &x) 
    {
        prepare(y,x);

        // vn = atan2(y, x)
        DScalar2 result(math::atan2(y.value, x.value));

        // Dvn = (x*Dy - y*Dx) / (x^2 + y^2)
        Scalar denom = x.value*x.value + y.value*y.value,
            denomSqr = denom*denom;
        result.grad = y.grad * (x.value / denom)
            - x.grad * (y.value / denom);

        // D^2vn = (Dy*Dx^T + xD^2y - Dx*Dy^T - yD^2x) / (x^2+y^2)
        //    - [(x*Dy - y*Dx) * (2*x*Dx + 2*y*Dy)^T] / (x^2+y^2)^2
        result.hess = (y.hess*x.value
                       + y.grad * x.grad.transpose()
                       - x.hess*y.value
                       - x.grad*y.grad.transpose()
            ) / denom;

        result.hess -=
            (y.grad*(x.value/denomSqr) - x.grad*(y.value/denomSqr)) *
            (x.grad*((Scalar) 2 * x.value) + y.grad*((Scalar) 2 * y.value)).transpose();

        return result;
    }

    friend DScalar2 atan(const DScalar2 &y) 
    {
        const Scalar denom =  y.value*y.value + 1.0;
        // vn = atan(y)
        // Dvn = Dy / (1 + y^2)
        // D^2vn = D^2y / (1+y^2) - [ 2*y * Dy*Dy^T] / (1+y^2)^2
        return DScalar2(math::atan(y.value), y.grad / denom, 
               (y.hess - 2.0 * y.value * y.grad * y.grad.transpose() / denom ) / denom );
    }

    friend DScalar2 tanh(const DScalar2 &y) 
    {
        const Scalar th = math::tanh(y.value);
        const Scalar tmp =  1 -  th*th;
        // vn = tanh(y)
        // Dvn = (1 - tanh(y)^2)*Dy
        // D^2vn = (1-tanh(y)^2) * [ D^2y - 2*tanh(y)* Dy*Dy^T ]
        return DScalar2(th, tmp * y.grad,
               tmp * (y.hess - 2.0 * th * y.grad * y.grad.transpose()) );
    }

    friend DScalar2 cosh(const DScalar2 &y) 
    {
        const Scalar ch = math::cosh(y.value);
        const Scalar sh = math::sinh(y.value);
        // vn = cosh(y)
        // Dvn = sinh(y)*Dy
        // D^2vn = cosh(y)*Dy*Dy^T + sinh(y)*D^2y
        return DScalar2(ch, sh * y.grad, ch * y.grad * y.grad.transpose() + sh * y.hess );
    }

    friend DScalar2 sinh(const DScalar2 &y) 
    {
        const Scalar ch = math::cosh(y.value);
        const Scalar sh = math::sinh(y.value);
        // vn = cosh(y)
        // Dvn = sinh(y)*Dy
        // D^2vn = cosh(y)*Dy*Dy^T + sinh(y)*D^2y
        return DScalar2(sh, ch * y.grad, sh * y.grad * y.grad.transpose() + ch * y.hess );
    }

    /// @}
    // ======================================================================

    // ======================================================================
    /// @{ \name Comparison and assignment
    // ======================================================================

    inline void operator=(const DScalar2& s ) { value = s.value; grad = s.grad; hess = s.hess; }
    inline void operator=(const Scalar &v   ) { value = v; grad.setZero(); hess.setZero(); }
    inline bool operator<(const DScalar2& s ) const { return value < s.value; }
    inline bool operator<(const Scalar& s   ) const { return value < s; }
    inline bool operator<=(const DScalar2& s) const { return value <= s.value; }
    inline bool operator<=(const Scalar& s  ) const { return value <= s; }
    inline bool operator>(const DScalar2& s ) const { return value > s.value; }
    inline bool operator>(const Scalar& s   ) const { return value > s; }
    inline bool operator>=(const DScalar2& s) const { return value >= s.value; }
    inline bool operator>=(const Scalar& s  ) const { return value >= s; }
    inline bool operator==(const DScalar2& s) const { return value == s.value; }
    inline bool operator==(const Scalar& s  ) const { return value == s; }
    inline bool operator!=(const DScalar2& s) const { return value != s.value; }
    inline bool operator!=(const Scalar& s  ) const { return value != s; }

    friend bool operator<(const Scalar s1, const DScalar2 & s2 ) { return  s1 <  s2.value; }
    friend bool operator>(const Scalar s1, const DScalar2 & s2 ) { return  s1 >  s2.value; }
    friend bool operator>=(const Scalar s1, const DScalar2 & s2) { return  s1 >= s2.value; }
    friend bool operator<=(const Scalar s1, const DScalar2 & s2) { return  s1 <= s2.value; }
    friend bool operator==(const Scalar s1, const DScalar2 & s2) { return  s1 == s2.value; }
    friend bool operator!=(const Scalar s1, const DScalar2 & s2) { return  s1 != s2.value; }

    /// @}
    // ======================================================================


protected:

    // prepare two scalars for a binary operation 
    static inline void prepare(const DScalar2 &lhs, const DScalar2 &rhs) 
    {
        if ( ! lhs.grad.size() )
        {
            const size_t nv = rhs.grad.size();
            lhs.grad.setZero(nv    );
            lhs.hess.setZero(nv, nv);
            return;
        }
        if ( ! rhs.grad.size() )
        {
            const size_t nv = lhs.grad.size();
            rhs.grad.setZero(nv    );
            rhs.hess.setZero(nv, nv);
            return;
        }
    }

protected:

    Scalar value;

    // mutable for treating incomplete constants
    mutable Gradient_t grad;

    mutable Hessian_t hess;
};


template <typename Scalar, int d>
std::ostream &operator<<(std::ostream &out, const DScalar2<Scalar,d> & s) 
{
	out << "[" << s.getValue()
		<< ", grad=" << s.getGradient().format(Eigen::IOFormat(4, 1, ", ", "; ", "", "", "[", "]"))
		<< ", hess=" << s.getHessian().format(Eigen::IOFormat(4, 0, ", ", "; ", "", "", "[", "]"))
		<< "]";
	return out;
}

}//namespace ad


}//namespace gismo
