/** @file gsFunctionExpr.h

    @brief Provides declaration of FunctionExpr class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunction.h>

namespace gismo
{

/**
    @brief Class defining a multivariate (real or vector) function
    given by a string mathematical expression.

    Numerous forms of functional and logic processing semantics are
    supported. The class is based on The C++ Mathematical Expression
    Toolkit Library (ExprTk), see

    http://www.partow.net/programming/exprtk

    and

    https://github.com/ArashPartow/exprtk/blob/master/readme.txt

    for more details.

    \ingroup function
    \ingroup Core
*/
template<typename T>
class gsFunctionExpr : public gsFunction<T>
{
public:
    typedef T Scalar_t;

    /// Shared pointer for gsFunctionExpr
    typedef memory::shared_ptr< gsFunctionExpr > Ptr;

    /// Unique pointer for gsFunctionExpr
    typedef memory::unique_ptr<gsFunctionExpr> uPtr;

public:

    /// Default empty constructor
    gsFunctionExpr();

    ///\brief Constructor taking an expression string and the domain dimension (scalar function)
    gsFunctionExpr(const std::string & expression_string, short_t ddim);

    ///\brief Constructor taking two expression strings (2D vector valued function)
    gsFunctionExpr(const std::string & expression_string1,
                   const std::string & expression_string2,
                   short_t ddim);

    ///\brief Constructor taking three expression strings (3D vector valued function)
    gsFunctionExpr(const std::string & expression_string1,
                   const std::string & expression_string2,
                   const std::string & expression_string3,
                   short_t ddim);

    ///\brief Constructor taking four expression strings (4D vector valued function) used for matrix coefficients
    gsFunctionExpr(const std::string & expression_string1,
                   const std::string & expression_string2,
                   const std::string & expression_string3,
                   const std::string & expression_string4,
                   short_t ddim);

    ///\brief Constructor taking nine expression strings (9D vector valued function) used for (3x3) matrix coefficients
    gsFunctionExpr(const std::string & expression_string1,
                   const std::string & expression_string2,
                   const std::string & expression_string3,
                   const std::string & expression_string4,
                   const std::string & expression_string5,
                   const std::string & expression_string6,
                   const std::string & expression_string7,
                   const std::string & expression_string8,
                   const std::string & expression_string9,
                   short_t ddim);

    gsFunctionExpr(const std::vector<std::string> & expression_string, short_t ddim);

    gsFunctionExpr(const gsFunctionExpr& other);

    ///\brief Make function taking an expression string and the domain dimension (scalar function)
    static uPtr make(const std::string & expression_string, short_t ddim)
    { return uPtr(new gsFunctionExpr(expression_string, ddim)); }

    ///\brief Make function taking two expression strings (2D vector valued function)
    static uPtr make(const std::string & expression_string1, const std::string & expression_string2, short_t ddim)
    { return uPtr(new gsFunctionExpr(expression_string1, expression_string2, ddim)); }

    ///\brief Make function taking two expression strings (3D vector valued function)
    static uPtr make(const std::string & expression_string1, const std::string & expression_string2,
                     const std::string & expression_string3, short_t ddim)
    { return uPtr(new gsFunctionExpr(expression_string1, expression_string2, expression_string3, ddim)); }

    ///\brief Make function taking four expression strings (4D vector valued function) used for matrix coefficients
    static uPtr make(const std::string & expression_string1, const std::string & expression_string2,
                     const std::string & expression_string3, const std::string & expression_string4, short_t ddim)
    { return uPtr(new gsFunctionExpr(expression_string1, expression_string2, expression_string3, expression_string4, ddim)); }

    ///\brief Make function taking nine expression strings (9D vector valued function) used for (3x3) matrix coefficients
    static uPtr make(const std::string & expression_string1, const std::string & expression_string2,
                     const std::string & expression_string3, const std::string & expression_string4,
                     const std::string & expression_string5, const std::string & expression_string6,
                     const std::string & expression_string7, const std::string & expression_string8,
                     const std::string & expression_string9, short_t ddim)
    { return uPtr(new gsFunctionExpr(expression_string1, expression_string2, expression_string3, expression_string4,
                expression_string5, expression_string6, expression_string7, expression_string8, expression_string9, ddim)); }

#if EIGEN_HAS_RVALUE_REFERENCES
    gsFunctionExpr(gsFunctionExpr&& other);

    gsFunctionExpr& operator=(const gsFunctionExpr& other);

    gsFunctionExpr& operator=(gsFunctionExpr&& other);
#else
    gsFunctionExpr& operator=(gsFunctionExpr other);
#endif

    ~gsFunctionExpr();

    GISMO_CLONE_FUNCTION(gsFunctionExpr)

    /// \brief Adds another component to this (vector) function
    void addComponent(const std::string & strExpression);

private:

    // initializes the symbol table
    void init(const short_t dim);

public:

    // Function expression can be used as a global function defined
    // for any real value, on any subdomain
    virtual const gsFunctionExpr & piece(const index_t) const
    {
        return *this;
    }

    // Documented in gsFunction class
    short_t domainDim() const;

    // Documented in gsFunction class
    short_t targetDim() const;

    // returns the string expression for component \a i
    const std::string & expression(int i = 0) const;

    /// Sets the symbol "x" to a value
    void set_x (T const & v) const;
    /// Sets the symbol "y" to a value
    void set_y (T const & v) const;
    /// Sets the symbol "z" to a value
    void set_z (T const & v) const;
    /// Sets the symbol "w" to a value
    void set_w (T const & v) const;
    /// Sets the symbol "u" to a value
    void set_u (T const & v) const;
    /// Sets the symbol "v" to a value
    void set_v (T const & v) const;
    /// Sets the symbol "t" to a value
    void set_t (T const & t) const;

    // see gsFunction for documentation
    virtual void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const;

    // see gsFunction for documentation
    virtual void eval_component_into(const gsMatrix<T>& u,
                                     const index_t comp,
                                     gsMatrix<T>& result) const;

    // see gsFunction for documentation
    virtual void deriv_into(const gsMatrix<T>& u,
                            gsMatrix<T>& result) const;

    // see gsFunction for documentation
    virtual void deriv2_into(const gsMatrix<T>& u,
                             gsMatrix<T>& result) const;

    // see gsFunction for documentation
    gsMatrix<T> hess(const gsMatrix<T>& u, unsigned coord = 0) const;

    // see gsFunction for documentation
    gsMatrix<T> laplacian(const gsMatrix<T>& u) const;

    ///Mixed derivative wrt variables k and j
    gsMatrix<T> * mderiv(const gsMatrix<T>& u, const index_t k, const index_t j) const;

    // see gsFunction for documentation
    std::ostream &print(std::ostream &os) const;

    void swap(gsFunctionExpr& other)
    { std::swap(this->my, other.my); }

// Data members
private:
    class gsFunctionExprPrivate;
    typedef gsFunctionExprPrivate PrivateData_t;

    PrivateData_t * my;

}; // class gsFunctionExpr


} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFunctionExpr.hpp)
#endif
