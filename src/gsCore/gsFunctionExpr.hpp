/** @file gsFunctionExpr.hpp

    @brief Provides implementation of FunctionExpr class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

// CRITICAL: Include autodiff Eigen support BEFORE gsLinearAlgebra.h
// to ensure Eigen NumTraits specializations are defined before Eigen is used
#ifdef gsAutoDiff_ENABLED
#include <gsAutoDiff/gsAutoDiffEigen.h>
#include <gsAutoDiff/gsAutoDiff.h>
#endif

#include <gsCore/gsLinearAlgebra.h>

// Include autodiff types AFTER gsLinearAlgebra.h when autodiff is enabled
// (gsAutoDiff_ENABLED is defined in gsConfigExt.h which is included via gsLinearAlgebra.h)
#ifdef gsAutoDiff_ENABLED
#include <gsAutoDiff/gsAutoDiff.h>
#endif

// Include autodiff type specializations (requires the optional module in include path)
#ifdef gsAutoDiff_ENABLED
#include <gsAutoDiff/src/gsAutoDiffTraits.h>
#endif

/* ExprTk options */

//This define will enable printing of debug information to stdout during
//the compilation process.
//#define exprtk_enable_debugging

// This define will disable the ability for expressions to have comments.
// Expressions that have comments when parsed with a build that has this
// option, will result in a compilation failure.
#define exprtk_disable_comments

// This define will disable the loop-wise 'break' and 'continue'
// capabilities. Any expression that contains those keywords will result
// in a compilation failure.
#define exprtk_disable_break_continue

// This define will disable the short-circuit '&' (and) and '|' (or)
// operators
#define exprtk_disable_sc_andor

// This define will disable all enhanced features such as strength
// reduction and special function optimisations and expression specific
// type instantiations. This feature will reduce compilation times and
// binary sizes but will also result in massive performance degradation
// of expression evaluations.
#if !defined(NDEBUG) || defined(__MINGW32__)
#define exprtk_disable_enhanced_features
#endif

// This define will disable all string processing capabilities. Any
// expression that contains a string or string related syntax will result
// in a compilation failure.
#define exprtk_disable_string_capabilities

#define exprtk_disable_rtl_io_file
#define exprtk_disable_rtl_vecops

// The order in which header files are included is essential.
//
// It is important that all forward declaration file
// "exprtk_X_forward.hpp"" are included BEFORE the header file
// "exprtk.hpp" so that specializations for is_true(), is_false(),
// etcetera are not known for ALL types that should be supported. All
// adaptor files "exprtk_X_adaptor.hpp" have to be included AFTER the
// file "exprtk.hpp".

// gsAutoDiff_ENABLED uses dual_t/var_t from the autodiff library
#if defined(gsMpfr_ENABLED)
#include <exprtk_mpfr_forward.hpp>
#endif

#if defined(gsGmp_ENABLED)
#include <exprtk_gmp_forward.hpp>
#endif

#if defined(gsCoDiPack_ENABLED)
#include <gsCoDiPack/exprtk_codi_forward.hpp>
#endif

#if defined(gsUniversal_ENABLED)
#include <gsUniversal/exprtk_universal_forward.hpp>
#endif

#if defined(gsAutoDiff_ENABLED)
#include <gsAutoDiff/exprtk_autodiff_forward.hpp>
#endif
#include <exprtk.hpp>

#if defined(gsMpfr_ENABLED)
#include <exprtk_mpfr_adaptor.hpp>
#endif

#if defined(gsGmp_ENABLED)
#include <exprtk_gmp_adaptor.hpp>
#endif

#if defined(gsCoDiPack_ENABLED)
#include <gsCoDiPack/exprtk_codi_adaptor.hpp>
#endif

#if defined(gsUniversal_ENABLED)
#include <gsUniversal/exprtk_universal_adaptor.hpp>
#endif

#if defined(gsAutoDiff_ENABLED)
#include <gsAutoDiff/exprtk_autodiff_adaptor.hpp>
#endif

#include <gsIO/gsXml.h>

namespace
{

// addition of mixed derivative for expressions
// see https://en.wikipedia.org/wiki/Finite_difference_coefficient
template <typename T>
T mixed_derivative(const exprtk::expression<T>& e,
                     T& x, T& y,
                     const double& h = 0.00001)
{
    T num = (T)(0.0), tmp;
    T x_init = x;
    T y_init = y;

    x = x_init + (T)(2.0) * h;
    y = y_init + (T)(2.0) * h;
    num += e.value();
    y = y_init - (T)(2.0) * h;
    num -= e.value();
    x = x_init - (T)(2.0) * h;
    num += e.value();
    y = y_init + (T)(2.0) * h;
    num -= e.value();

    x = x_init + h;
    y = y_init + h;
    tmp = e.value();
    y = y_init - h;
    tmp -= e.value();
    x = x_init - h;
    tmp += e.value();
    y = y_init + h;
    tmp -= e.value();
    num += (T)(64.0) * tmp;

    x = x_init + (T)(2.0) * h;
    y = y_init - h;
    tmp = e.value();
    y = y_init + h;
    tmp -= e.value();
    x = x_init - (T)(2.0) * h;
    tmp += e.value();
    y = y_init - h;
    tmp -= e.value();

    y = y_init + (T)(2.0) * h;
    x = x_init - h;
    tmp += e.value();
    x = x_init + h;
    tmp -= e.value();
    y = y_init - (T)(2.0) * h;
    tmp += e.value();
    x = x_init - h;
    tmp -= e.value();
    num += (T)(8.0) * tmp;

    x = x_init;
    y = y_init;
    return num / ( (T)(144.0)*h*h );
}

} //namespace

#define N_VARS 7

namespace gismo
{

template<typename T> class gsFunctionExpr<T>::gsFunctionExprPrivate
{
public:

// Numeric_t selection:
// - If gsAutoDiff_ENABLED: use dual2nd_t for computing first AND second derivatives via AD
// - Otherwise: use T directly
#if defined(gsAutoDiff_ENABLED)
    typedef dual2nd_t Numeric_t;
#else
    typedef T Numeric_t;
#endif

    typedef exprtk::symbol_table<Numeric_t>  SymbolTable_t;
    typedef exprtk::expression<Numeric_t>    Expression_t;
    typedef exprtk::parser<Numeric_t>        Parser_t;

public:

    gsFunctionExprPrivate(const short_t _dim)
    : vars(), dim(_dim)
    {
        GISMO_ENSURE( dim <= N_VARS, "The number of variables can be at most 7 (x,y,z,w,u,v,t)." );
        init();
    }

    gsFunctionExprPrivate(const gsFunctionExprPrivate & other)
    : vars(), dim(other.dim)
    {
        GISMO_ASSERT ( string.size() == expression.size(), "Corrupted FunctionExpr");
        init();
        //copy_n(other.vars, N_VARS+1, vars);
        string    .reserve(string.size());
        expression.reserve(string.size());
        for (size_t i = 0; i!= other.string.size(); ++i)
            addComponent(other.string[i]);
    }

    void addComponent(const std::string & strExpression)
    {
        string.push_back( strExpression );// Keep string data
        std::string & str = string.back();
        str.erase(std::remove(str.begin(), str.end(),' '), str.end() );
        gismo::util::string_replace(str, "**", "^");

        // String expression
        expression.push_back(Expression_t());
        Expression_t & expr = expression.back();
        //expr.release();
        expr.register_symbol_table(symbol_table);

        // Parser
        Parser_t parser;
        //Collect variable symbols
        //parser.dec().collect_variables() = true;
        bool success = parser.compile(str, expr);
        if ( ! success )
            gsWarn<<"gsFunctionExpr error: " <<parser.error() <<" while parsing "<<str<<"\n";
        /*
           typedef typename exprtk::parser_t::
           dependent_entity_collector::symbol_t symbol_t;

           std::deque<symbol_t> symbol_list;
           parser.dec().symbols(symbol_list);
           for (size_t i = 0; i != symbol_list.size(); ++i)
           {
           symbol_t& symbol = symbol_list[i];
           // do something
           }
        */
    }

    void init()
    {
        //symbol_table.clear();
        // Identify symbol table
        symbol_table.add_variable("x",vars[0]);
        symbol_table.add_variable("y",vars[1]);
        symbol_table.add_variable("z",vars[2]);
        symbol_table.add_variable("w",vars[3]);
        symbol_table.add_variable("u",vars[4]);
        symbol_table.add_variable("v",vars[5]);
        symbol_table.add_variable("t",vars[6]);
        //symbol_table.remove_variable("w",vars[3]);
        symbol_table.add_pi();
        //symbol_table.add_constant("C", 1);

#       if defined(gsAutoDiff_ENABLED)
        // Warn when using backward AD with expression parser
        if constexpr (std::is_same_v<T, var_t>) {
            static bool warned = false;
            if (!warned) {
                gsWarn << "Warning: gsFunctionExpr instantiated with backward AD type (var_t). "
                       << "Derivatives will be computed using forward-mode autodiff internally. "
                       << "Backward AD does not propagate through expression parsing.\n";
                warned = true;
            }
        }
#       endif
    }

public:
    mutable Numeric_t         vars[N_VARS];
    SymbolTable_t             symbol_table;
    std::vector<Expression_t> expression;
    std::vector<std::string>  string;
    short_t dim;

private:
    gsFunctionExprPrivate();
    gsFunctionExprPrivate operator= (const gsFunctionExprPrivate & other);
};

/* / /AM: under construction
template<typename T> class gsSymbolListPrivate
{
public:
    exprtk::symbol_table<T> symbol_table;
    std::vector<T> vars  ;
    std::vector<T> params;
};

template<typename T>
gsSymbolList<T>::gsSymbolList() : my(new gsSymbolListPrivate<T>) {}

template<typename T>
gsSymbolList<T>::~gsSymbolList()
{
delete my;
}

template<typename T>
void gsSymbolList<T>::setDefault()
{
    my->symbol_table.clear();

    // Identify symbol table
    my->symbol_table.add_variable("x",my->vars[0]);
    my->symbol_table.add_variable("y",my->vars[1]);
    my->symbol_table.add_variable("z",my->vars[2]);
    my->symbol_table.add_variable("w",my->vars[3]);
    my->symbol_table.add_variable("u",my->vars[4]);
    my->symbol_table.add_variable("v",my->vars[5]);
    //my->symbol_table.remove_variable("w",my->vars[3]);
    my->symbol_table.add_pi();
    //my->symbol_table.add_constant("C", 1);
}

template<typename T>
bool gsSymbolList<T>::addConstant(const std::string & constant_name, const T& value)
{
return my->symbol_table.add_constant(constant_name, value);
}

template<typename T>
bool gsSymbolList<T>::addVariable(const std::string & variable_name)
{
    my->vars.push_back(0.0);
    return my->symbol_table.add_variable(variable_name, my->vars.back() );
}

template<typename T>
bool gsSymbolList<T>::addParameter(const std::string & variable_name)
{
    my->params.push_back(0.0);
    return my->symbol_table.add_variable(variable_name, my->params.back() );
}

template<typename T>
bool gsSymbolList<T>::hasSymbol(const std::string& symbol_name)
{
return true;
}
*/

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr() : my(new PrivateData_t(0))
{ }

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::string & expression_string, short_t ddim)
: my(new PrivateData_t(ddim))
{
    my->addComponent(expression_string);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::string & expression_string1,
                                  const std::string & expression_string2,
                                  short_t ddim)
: my(new PrivateData_t(ddim))
{
    my->addComponent(expression_string1);
    my->addComponent(expression_string2);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::string & expression_string1,
                                  const std::string & expression_string2,
                                  const std::string & expression_string3,
                                  short_t ddim)
: my(new PrivateData_t(ddim))
{
    my->addComponent(expression_string1);
    my->addComponent(expression_string2);
    my->addComponent(expression_string3);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::string & expression_string1,
                   const std::string & expression_string2,
                   const std::string & expression_string3,
                   const std::string & expression_string4,
                   short_t ddim)
: my(new PrivateData_t(ddim))
{
    my->addComponent(expression_string1);
    my->addComponent(expression_string2);
    my->addComponent(expression_string3);
    my->addComponent(expression_string4);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::string & expression_string1,
                   const std::string & expression_string2,
                   const std::string & expression_string3,
                   const std::string & expression_string4,
                   const std::string & expression_string5,
                   const std::string & expression_string6,
                   const std::string & expression_string7,
                   const std::string & expression_string8,
                   const std::string & expression_string9,
                   short_t ddim)
: my(new PrivateData_t(ddim))
{
    my->addComponent(expression_string1);
    my->addComponent(expression_string2);
    my->addComponent(expression_string3);
    my->addComponent(expression_string4);
    my->addComponent(expression_string5);
    my->addComponent(expression_string6);
    my->addComponent(expression_string7);
    my->addComponent(expression_string8);
    my->addComponent(expression_string9);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::vector<std::string> & expression_string,
                                  short_t ddim)
: my(new PrivateData_t(ddim))
{
    for (size_t i = 0; i!= expression_string.size(); ++i)
        my->addComponent(expression_string[i]);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const gsFunctionExpr & other)
{
    my = new PrivateData_t(*other.my);
}
#if EIGEN_HAS_RVALUE_REFERENCES

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(gsFunctionExpr && other)
{
    my = other.my; other.my = NULL;
}

template<typename T>
gsFunctionExpr<T> & gsFunctionExpr<T>::operator=(const gsFunctionExpr& other)
{
    if (this != &other)
    {
        delete my;
        my = new PrivateData_t(*other.my);
    }
    return *this;
}

template<typename T>
gsFunctionExpr<T> & gsFunctionExpr<T>::operator=(gsFunctionExpr&& other)
{
    if (this != &other)
    {
        delete my;
        my = other.my; other.my = NULL;
    }
    return *this;
}

#else
template<typename T>
gsFunctionExpr<T> & gsFunctionExpr<T>::operator=(gsFunctionExpr other)
{
    std::swap(my,other.my);
    return *this;
}
#endif

#if defined(gsAutoDiff_ENABLED)
// Conversion helpers for backward AD support
namespace internal {
    template <typename T>
    static inline dual2nd_t to_dual2nd(const T& val) {
        if constexpr (std::is_same_v<T, dual2nd_t>) {
            return val;
        } else if constexpr (std::is_same_v<T, dual_t>) {
            // dual_t: construct dual2nd_t with value=val and gradient=0
            // Dual struct has public val and grad members
            dual2nd_t result;
            result.val = val;
            result.grad = dual_t(0.0);
            return result;
        } else if constexpr (std::is_same_v<T, var_t>) {
            // var_t is reverse-mode AD, extract value and construct dual2nd_t
            dual2nd_t result;
            result.val = dual_t(autodiff::val(val));
            result.grad = dual_t(0.0);
            return result;
        } else {
            // Scalar types: construct dual2nd_t with value only
            dual2nd_t result;
            result.val = dual_t(static_cast<double>(val));
            result.grad = dual_t(0.0);
            return result;
        }
    }

    template <typename T>
    static inline T from_dual2nd(const dual2nd_t& val) {
        if constexpr (std::is_same_v<T, dual2nd_t>) {
            return val;
        } else if constexpr (std::is_same_v<T, dual_t>) {
            // Extract value and first derivative from dual2nd_t
            // dual2nd_t.val is the value (dual_t), dual2nd_t.grad is the gradient (dual_t)
            // For first derivative: use grad.val (the value component of the gradient)
            dual_t result;
            result.val = val.val.val;
            result.grad = val.grad.val;
            return result;
        } else if constexpr (std::is_same_v<T, var_t>) {
            // Return value only; backward AD derivatives don't propagate through expression parser
            return var_t(val.val.val);
        } else {
            // Scalar types: extract just the value
            return static_cast<T>(val.val.val);
        }
    }
}
#endif

template<typename T>
gsFunctionExpr<T>::~gsFunctionExpr()
{
    delete my;
}

template<typename T>
short_t gsFunctionExpr<T>::domainDim() const
{
    return my->dim;
}

template<typename T>
short_t gsFunctionExpr<T>::targetDim() const
{
    return static_cast<short_t>(my->string.size());
}

template<typename T>
void gsFunctionExpr<T>::addComponent(const std::string & strExpression)
{
    my->addComponent(strExpression);
}

template<typename T>
const std::string & gsFunctionExpr<T>::expression(int i) const
{
    return my->string[i];
}

template<typename T>
void gsFunctionExpr<T>::set_x (T const & v) const {
#ifdef gsAutoDiff_ENABLED
    my->vars[0]= internal::to_dual2nd(v);
#else
    my->vars[0]= v;
#endif
}

template<typename T>
void gsFunctionExpr<T>::set_y (T const & v) const {
#ifdef gsAutoDiff_ENABLED
    my->vars[1]= internal::to_dual2nd(v);
#else
    my->vars[1]= v;
#endif
}

template<typename T>
void gsFunctionExpr<T>::set_z (T const & v) const {
#ifdef gsAutoDiff_ENABLED
    my->vars[2]= internal::to_dual2nd(v);
#else
    my->vars[2]= v;
#endif
}

template<typename T>
void gsFunctionExpr<T>::set_w (T const & v) const {
#ifdef gsAutoDiff_ENABLED
    my->vars[3]= internal::to_dual2nd(v);
#else
    my->vars[3]= v;
#endif
}

template<typename T>
void gsFunctionExpr<T>::set_u (T const & v) const {
#ifdef gsAutoDiff_ENABLED
    my->vars[4]= internal::to_dual2nd(v);
#else
    my->vars[4]= v;
#endif
}

template<typename T>
void gsFunctionExpr<T>::set_v (T const & v) const {
#ifdef gsAutoDiff_ENABLED
    my->vars[5]= internal::to_dual2nd(v);
#else
    my->vars[5]= v;
#endif
}

template<typename T>
void gsFunctionExpr<T>::set_t (T const & t) const {
#ifdef gsAutoDiff_ENABLED
    my->vars[6]= internal::to_dual2nd(t);
#else
    my->vars[6]= t;
#endif
}

template<typename T>
void gsFunctionExpr<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")\n"<< *this);

    const short_t n = targetDim();
    result.resize(n, u.cols());

#pragma omp critical (gsFunctionExpr_run)
    for ( index_t p = 0; p!=u.cols(); p++ ) // for all evaluation points
    {
#ifdef gsAutoDiff_ENABLED
        for (index_t k = 0; k!=my->dim; ++k)
            my->vars[k] = internal::to_dual2nd(u(k,p));
#else
        copy_n(u.col(p).data(), my->dim, my->vars);
#endif

        for (short_t c = 0; c!= n; ++c) // for all components
#       if defined(gsAutoDiff_ENABLED)
            result(c,p) = internal::from_dual2nd<T>(my->expression[c].value());
#       else
            result(c,p) = my->expression[c].value();
#       endif
    }
}

template<typename T>
void gsFunctionExpr<T>::eval_component_into(const gsMatrix<T>& u, const index_t comp, gsMatrix<T>& result) const
{
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")");

    GISMO_ASSERT (comp < targetDim(),
                  "Given component number is higher then number of components");

    result.resize(1, u.cols());
#   pragma omp critical (gsFunctionExpr_run)
    for ( index_t p = 0; p!=u.cols(); ++p )
    {
#ifdef gsAutoDiff_ENABLED
        for (index_t k = 0; k!=my->dim; ++k)
            my->vars[k] = internal::to_dual2nd(u(k,p));
#else
        copy_n(u.col(p).data(), my->dim, my->vars);
#endif

#       if defined(gsAutoDiff_ENABLED)
            result(0,p) = internal::from_dual2nd<T>(my->expression[comp].value());
#       else
            result(0,p) = my->expression[comp].value();
#       endif
    }
}

template<typename T>
void gsFunctionExpr<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    const short_t d = domainDim();
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")");

    const short_t n = targetDim();
    result.resize(d*n, u.cols());
#   pragma omp critical (gsFunctionExpr_run)
    for ( index_t p = 0; p!=u.cols(); p++ ) // for all evaluation points
    {
#       if defined(gsAutoDiff_ENABLED)
        // Using dual2nd_t: seed each variable and extract first derivative
        // dual2nd_t.val is dual_t, dual2nd_t.grad is dual_t
        // First derivative: d f / d x_j = result.grad.val when x_j is seeded
        for (short_t c = 0; c!= n; ++c) // for all components
        {
            for (short_t j = 0; j!=d; ++j) // for all variables
            {
                // Set all variables to their values (seed for first derivative)
                for (short_t k = 0; k!=d; ++k)
                {
                    // val is dual_t, grad is dual_t
                    dual2nd_t tmp = internal::to_dual2nd(u(k,p));
                    my->vars[k].val.val = tmp.val.val;
                    my->vars[k].val.grad = tmp.val.grad;
                    my->vars[k].grad.val = (k==j ? 1.0 : 0.0);  // Seed only variable j
                    my->vars[k].grad.grad = 0.0;
                }
                // Evaluate and extract first derivative
                dual2nd_t expr_val = my->expression[c].value();
                result(c*d + j, p) = expr_val.grad.val;
            }
        }
#       else
        copy_n(u.col(p).data(), my->dim, my->vars);
        for (short_t c = 0; c!= n; ++c) // for all components
            for ( short_t j = 0; j!=d; j++ ) // for all variables
                result(c*d + j, p) =
                    exprtk::derivative<T>(my->expression[c], my->vars[j], 0.00001 ) ;
#       endif
    }
}

template<typename T>
void gsFunctionExpr<T>::deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    const short_t d = domainDim();
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")");

    const short_t n = targetDim();
    const index_t stride = d + d*(d-1)/2;
    result.resize(stride*n, u.cols() );
#   pragma omp critical (gsFunctionExpr_run)
    for ( index_t p = 0; p!=u.cols(); p++ ) // for all evaluation points
    {
        for (short_t c = 0; c!= n; ++c) // for all components
        {
#           if defined(gsAutoDiff_ENABLED)
            // Using dual2nd_t for second derivatives via AD
            // For H_{k,l}, we need to seed x_k in grad and x_l in grad.grad
            short_t m = d;
            for (short_t k = 0; k!=d; ++k)
            {
                // H_{k,k} - second derivative d²f/dx_k²
                for (short_t v = 0; v!=d; ++v)
                {
                    dual2nd_t tmp = internal::to_dual2nd(u(v,p));
                    my->vars[v].val.val = tmp.val.val;
                    my->vars[v].val.grad = (v==k ? 1.0 : 0.0);  // seed for inner derivative
                    my->vars[v].grad.val = (v==k ? 1.0 : 0.0);  // seed for outer derivative  
                    my->vars[v].grad.grad = 0.0;
                }
                dual2nd_t expr_val = my->expression[c].value();
                result(c*stride + k, p) = expr_val.grad.grad;

                for (short_t l=k+1; l<d; ++l)
                {
                    // H_{k,l} - mixed derivative d²f/(dx_k dx_l)
                    for (short_t v = 0; v!=d; ++v)
                    {
                        dual2nd_t tmp = internal::to_dual2nd(u(v,p));
                        my->vars[v].val.val = tmp.val.val;
                        my->vars[v].val.grad = (v==l ? 1.0 : 0.0);  // seed x_l for inner
                        my->vars[v].grad.val = (v==k ? 1.0 : 0.0);  // seed x_k for outer
                        my->vars[v].grad.grad = 0.0;
                    }
                    expr_val = my->expression[c].value();
                    result(c*stride + m++, p) = expr_val.grad.grad;
                }
            }
#           else
            copy_n(u.col(p).data(), my->dim, my->vars);
            short_t m = d;
            for (short_t k = 0; k!=d; ++k)
            {
                // H_{k,k}
                result(c*stride + k,p) = exprtk::
                    second_derivative<T>(my->expression[c], my->vars[k], 0.00001);

                for (short_t l=k+1; l<d; ++l)
                {
                    // H_{k,l}
                    result(c*stride + m++,p) =
                        mixed_derivative<T>( my->expression[c], my->vars[k],
                                             my->vars[l], 0.00001 );
                }
            }
#           endif
        }
    }
}

template<typename T>
gsMatrix<T>
gsFunctionExpr<T>::hess(const gsMatrix<T>& u, unsigned coord) const
{
    GISMO_ENSURE(coord == 0, "Error, function is real");
    GISMO_ASSERT ( u.cols() == 1, "Need a single evaluation point." );
    const index_t d = u.rows();
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")");

    gsMatrix<T> res(d, d);

#   pragma omp critical (gsFunctionExpr_run)
{
#   if defined(gsAutoDiff_ENABLED)
    // Using dual2nd_t for Hessian via AD
    for( index_t j=0; j!=d; ++j )
    {
        // H_{j,j}
        for (index_t v = 0; v!=d; ++v)
        {
            dual2nd_t tmp = internal::to_dual2nd(u(v,0));
            my->vars[v].val.val = tmp.val.val;
            my->vars[v].val.grad = (v==j ? 1.0 : 0.0);
            my->vars[v].grad.val = (v==j ? 1.0 : 0.0);
            my->vars[v].grad.grad = 0.0;
        }
        dual2nd_t expr_val = my->expression[coord].value();
        res(j,j) = expr_val.grad.grad;

        for( index_t k = 0; k!=j; ++k )
        {
            // H_{k,j} = H_{j,k}
            for (index_t v = 0; v!=d; ++v)
            {
                dual2nd_t tmp = internal::to_dual2nd(u(v,0));
                my->vars[v].val.val = tmp.val.val;
                my->vars[v].val.grad = (v==j ? 1.0 : 0.0);
                my->vars[v].grad.val = (v==k ? 1.0 : 0.0);
                my->vars[v].grad.grad = 0.0;
            }
            expr_val = my->expression[coord].value();
            res(k,j) = res(j,k) = expr_val.grad.grad;
        }
    }
#   else
    copy_n(u.data(), my->dim, my->vars);
    for( index_t j=0; j!=d; ++j )
    {
        res(j,j) = exprtk::
            second_derivative<T>( my->expression[coord], my->vars[j], 0.00001);

        for( index_t k = 0; k!=j; ++k )
            res(k,j) = res(j,k) =
                mixed_derivative<T>( my->expression[coord], my->vars[k],
                                     my->vars[j], 0.00001 );
    }
#   endif
}
    return res;
}

template<typename T>
gsMatrix<T> * gsFunctionExpr<T>::mderiv(const gsMatrix<T> & u,
                                        const index_t k,
                                        const index_t j) const
{
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point size.");
    const short_t n = targetDim();
    gsMatrix<T> * res= new gsMatrix<T>(n,u.cols()) ;

#   pragma omp critical (gsFunctionExpr_run)
    for( index_t p=0; p!=res->cols(); ++p )
    {
        for (short_t c = 0; c!= n; ++c) // for all components
        {
#           if defined(gsAutoDiff_ENABLED)
            // Mixed derivative d²f/(dx_k dx_j) via dual2nd_t
            for (index_t v = 0; v!=my->dim; ++v)
            {
                dual2nd_t tmp = internal::to_dual2nd(u(v,p));
                my->vars[v].val.val = tmp.val.val;
                my->vars[v].val.grad = (v==j ? 1.0 : 0.0);
                my->vars[v].grad.val = (v==k ? 1.0 : 0.0);
                my->vars[v].grad.grad = 0.0;
            }
            dual2nd_t expr_val = my->expression[c].value();
            (*res)(c,p) = expr_val.grad.grad;
#           else
            copy_n(u.col(p).data(), my->dim, my->vars);
            (*res)(c,p) =
                mixed_derivative<T>( my->expression[c], my->vars[k], my->vars[j], 0.00001 ) ;
#           endif
        }
    }
    return res;
}

template<typename T>
gsMatrix<T> gsFunctionExpr<T>::laplacian(const gsMatrix<T>& u) const
{
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point size.");
    const short_t n = targetDim();
    gsMatrix<T> res(n,u.cols());
    res.setZero();

    #   pragma omp critical (gsFunctionExpr_run)
    for( index_t p = 0; p != res.cols(); ++p )
    {
        for (short_t c = 0; c!= n; ++c) // for all components
        {
#           if defined(gsAutoDiff_ENABLED)
            // Laplacian = sum of d²f/dx_j² via dual2nd_t
            T & val = res(c,p);
            for ( index_t j = 0; j!=my->dim; ++j )
            {
                for (index_t v = 0; v!=my->dim; ++v)
                {
                    dual2nd_t tmp = internal::to_dual2nd(u(v,p));
                    my->vars[v].val.val = tmp.val.val;
                    my->vars[v].val.grad = (v==j ? 1.0 : 0.0);
                    my->vars[v].grad.val = (v==j ? 1.0 : 0.0);
                    my->vars[v].grad.grad = 0.0;
                }
                dual2nd_t expr_val = my->expression[c].value();
                val += expr_val.grad.grad;
            }
#           else
            copy_n(u.col(p).data(), my->dim, my->vars);
            T & val = res(c,p);
            for ( index_t j = 0; j!=my->dim; ++j )
                val += exprtk::
                    second_derivative<T>( my->expression[c], my->vars[j], 0.00001 );
#           endif
        }
    }
    return  res;
}

template<typename T>
std::ostream & gsFunctionExpr<T>::print(std::ostream &os) const
{
    os <<"[ ";
    if( my->string.empty() )
        os << "empty";
    else
    {
        os << my->string[0];

        for (short_t k = 1; k<targetDim(); ++k)
            os <<", " << my->string[k];
    }
    os <<" ]";
    return os;
}

namespace internal
{

/// @brief Get a FunctionsExpr from XML data
template<class T>
class gsXml< gsFunctionExpr<T> >
{
private:
    gsXml() { }
    typedef gsFunctionExpr<T> Object;
public:
    GSXML_COMMON_FUNCTIONS(Object);
    GSXML_GET_INTO(Object);
    static std::string tag ()  { return "Function"; }
    static std::string type () { return "FunctionExpr"; }

    static Object * get(gsXmlNode * node)
    {
        GISMO_ASSERT( ( !strcmp( node->name(),"Function") )
                    &&  ( !strcmp(node->first_attribute("type")->value(),
                                internal::gsXml<Object>::type().c_str() ) ),
                    "Reading gsFunctionExpr XML: No Function found" );

        GISMO_ASSERT( node->first_attribute("dim"), "Reading gsFunctionExpr XML: No dim found" ) ;
        const int d = atoi( node->first_attribute("dim")->value() );

        std::vector< std::string > expr_strings;

        gsXmlNode * child = node->first_node("c");

        if (child != NULL )
        {
            for (; child; child = child->next_sibling() )
                expr_strings.push_back(  child->value() );
        }
        else
            expr_strings.push_back(  node->value() );

        return new Object(expr_strings, d);
    }

    static gsXmlNode * put (const Object & obj,
                            gsXmlTree & data )
    {
        // Add a new node
        gsXmlNode* node = internal::makeNode("Function" , data);
        node->append_attribute( makeAttribute("type",
                                            internal::gsXml< Object >::type().c_str(), data) );
        node->append_attribute(makeAttribute("dim", obj.domainDim(), data));


        const short_t tdim = obj.targetDim();

        if ( tdim == 1)
        {
            node->value( makeValue(obj.expression(), data) );
        }
        else
        {
            gsXmlNode * cnode;
            for (short_t c = 0; c!=tdim; ++c)
            {
                cnode = makeNode("c", obj.expression(c), data);
                node->append_node(cnode);
            }
        }

        return node;
    }
};

} // internal

}; // namespace gismo
