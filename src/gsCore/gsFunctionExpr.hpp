/** @file gsFunctionExpr.hpp

    @brief Provides implementation of FunctionExpr class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

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

#if defined(GISMO_WITH_ADIFF)
#define DScalar gismo::ad::DScalar2<real_t,-1>
#include <exprtk_ad_forward.hpp>
#endif

#if defined(GISMO_WITH_MPFR)
#include <exprtk_mpfr_forward.hpp>
#endif

#if defined(GISMO_WITH_GMP)
#include <exprtk_gmp_forward.hpp>
#endif

#if defined(GISMO_WITH_CODIPACK)
#include <gsCoDiPack/exprtk_codi_rf_forward.hpp>
#include <gsCoDiPack/exprtk_codi_rr_forward.hpp>
#endif

#if defined(GISMO_WITH_UNUM)
#include <gsUnum/exprtk_unum_posit_forward.hpp>
#endif

#include <exprtk.hpp>

#if defined(GISMO_WITH_ADIFF)
#include <exprtk_ad_adaptor.hpp>
#endif

#if defined(GISMO_WITH_MPFR)
#include <exprtk_mpfr_adaptor.hpp>
#endif

#if defined(GISMO_WITH_GMP)
#include <exprtk_gmp_adaptor.hpp>
#endif

#if defined(GISMO_WITH_CODIPACK)
#include <gsCoDiPack/exprtk_codi_rf_adaptor.hpp>
#include <gsCoDiPack/exprtk_codi_rr_adaptor.hpp>
#endif

#if defined(GISMO_WITH_UNUM)
#include <gsUnum/exprtk_unum_posit_adaptor.hpp>
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
    T num = T(0.0), tmp;
    T x_init = x;
    T y_init = y;

    x = x_init + T(2.0) * h;
    y = y_init + T(2.0) * h;
    num += e.value();
    y = y_init - T(2.0) * h;
    num -= e.value();
    x = x_init - T(2.0) * h;
    num += e.value();
    y = y_init + T(2.0) * h;
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
    num += 64* tmp;

    x = x_init + T(2.0) * h;
    y = y_init - h;
    tmp = e.value();
    y = y_init + h;
    tmp -= e.value();
    x = x_init - T(2.0) * h;
    tmp += e.value();
    y = y_init - h;
    tmp -= e.value();

    y = y_init + T(2.0) * h;
    x = x_init - h;
    tmp += e.value();
    x = x_init + h;
    tmp -= e.value();
    y = y_init - T(2.0) * h;
    tmp += e.value();
    x = x_init - h;
    tmp -= e.value();
    num += 8* tmp;

    x = x_init;
    y = y_init;
    return num / ( T(144.0)*h*h );
}

} //namespace

#define N_VARS 7

namespace gismo
{

template<typename T> class gsFunctionExpr<T>::gsFunctionExprPrivate
{
public:

#ifdef GISMO_WITH_ADIFF
    typedef DScalar Numeric_t;
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
        //*/
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
//*/

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
void gsFunctionExpr<T>::set_x (T const & v) const { my->vars[0]= v; }

template<typename T>
void gsFunctionExpr<T>::set_y (T const & v) const { my->vars[1]= v; }

template<typename T>
void gsFunctionExpr<T>::set_z (T const & v) const { my->vars[2]= v; }

template<typename T>
void gsFunctionExpr<T>::set_w (T const & v) const { my->vars[3]= v; }

template<typename T>
void gsFunctionExpr<T>::set_u (T const & v) const { my->vars[4]= v; }

template<typename T>
void gsFunctionExpr<T>::set_v (T const & v) const { my->vars[5]= v; }

template<typename T>
void gsFunctionExpr<T>::set_t (T const & t) const { my->vars[6]= t; }

template<typename T>
void gsFunctionExpr<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")\n"<< *this);

    const short_t n = targetDim();
    result.resize(n, u.cols());

    const PrivateData_t & expr =
#   ifdef _OPENMP
        omp_in_parallel() ? PrivateData_t(*my) :
#   endif
        *my;

    for ( index_t p = 0; p!=u.cols(); p++ ) // for all evaluation points
    {
        copy_n(u.col(p).data(), expr.dim, expr.vars);

        for (short_t c = 0; c!= n; ++c) // for all components
#           ifdef GISMO_WITH_ADIFF
            result(c,p) = expr.expression[c].value().getValue();
#           else
            result(c,p) = expr.expression[c].value();
#           endif
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
    for ( index_t p = 0; p!=u.cols(); ++p )
    {
        copy_n(u.col(p).data(), my->dim, my->vars);

#           ifdef GISMO_WITH_ADIFF
            result(0,p) = my->expression[comp].value().getValue();
#           else
            result(0,p) = my->expression[comp].value();
#           endif
    }
}

template<typename T>
void gsFunctionExpr<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    //gsDebug<< "Using finite differences (gsFunctionExpr::deriv_into) for derivatives.\n";
    const short_t d = domainDim();
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")");

    const short_t n = targetDim();
    result.resize(d*n, u.cols());

    const PrivateData_t & expr =
#   ifdef _OPENMP
        omp_in_parallel() ? PrivateData_t(*my) :
#   endif
        *my;

    for ( index_t p = 0; p!=u.cols(); p++ ) // for all evaluation points
    {
#       ifdef GISMO_WITH_ADIFF
        for (short_t k = 0; k!=d; ++k)
            expr.vars[k].setVariable(k,d,u(k,p));
        for (short_t c = 0; c!= n; ++c) // for all components
            expr.expression[c].value().gradient_into(result.block(c*d,p,d,1));
            //result.block(c*d,p,d,1) = expr.expression[c].value().getGradient(); //fails on constants
#       else
        copy_n(u.col(p).data(), expr.dim, expr.vars);
        for (short_t c = 0; c!= n; ++c) // for all components
            for ( short_t j = 0; j!=d; j++ ) // for all variables
                result(c*d + j, p) =
                    exprtk::derivative<T>(expr.expression[c], expr.vars[j], 0.00001 ) ;
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

    const PrivateData_t & expr =
#   ifdef _OPENMP
        omp_in_parallel() ? PrivateData_t(*my) :
#   endif
        *my;

    for ( index_t p = 0; p!=u.cols(); p++ ) // for all evaluation points
    {
#       ifndef GISMO_WITH_ADIFF
        copy_n(u.col(p).data(), expr.dim, expr.vars);
#       endif

        for (short_t c = 0; c!= n; ++c) // for all components
        {
#           ifdef GISMO_WITH_ADIFF
            for (index_t v = 0; v!=d; ++v)
                expr.vars[v].setVariable(v,d,u(v,p));
            const DScalar &            ads  = expr.expression[c].value();
            const DScalar::Hessian_t & Hmat = ads.getHessian(); // note: can fail

            for ( index_t k=0; k!=d; ++k)
            {
                result(k,p) = Hmat(k,k);
                index_t m = d;
                for ( index_t l=k+1; l<d; ++l)
                    result(m++,p) = Hmat(k,l);
            }
#           else
            for (short_t k = 0; k!=d; ++k)
            {
                // H_{k,k}
                result(k,p) = exprtk::
                    second_derivative<T>(expr.expression[c], expr.vars[k], 0.00001);

                short_t m = d;
                for (short_t l=k+1; l<d; ++l)
                {
                    // H_{k,l}
                    result(m++,p) =
                        mixed_derivative<T>( expr.expression[c], expr.vars[k],
                                             expr.vars[l], 0.00001 );
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
    //gsDebug<< "Using finite differences (gsFunctionExpr::hess) for Hessian.\n";
    GISMO_ENSURE(coord == 0, "Error, function is real");
    GISMO_ASSERT ( u.cols() == 1, "Need a single evaluation point." );
    const index_t d = u.rows();
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")");

    gsMatrix<T> res(d, d);

    const PrivateData_t & expr =
#   ifdef _OPENMP
        omp_in_parallel() ? PrivateData_t(*my) :
#   endif
        *my;

#   ifdef GISMO_WITH_ADIFF
    for (index_t v = 0; v!=d; ++v)
        expr.vars[v].setVariable(v, d, u(v,0) );
    expr.expression[coord].value().hessian_into(res);
#   else
    copy_n(u.data(), expr.dim, expr.vars);
    for( index_t j=0; j!=d; ++j )
    {
        res(j,j) = exprtk::
            second_derivative<T>( expr.expression[coord], expr.vars[j], 0.00001);

        for( index_t k = 0; k!=j; ++k )
            res(k,j) = res(j,k) =
                mixed_derivative<T>( expr.expression[coord], expr.vars[k],
                                     expr.vars[j], 0.00001 );
    }
#   endif

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

    const PrivateData_t & expr =
#   ifdef _OPENMP
        omp_in_parallel() ? PrivateData_t(*my) :
#   endif
        *my;

    for( index_t p=0; p!=res->cols(); ++p )
    {
#       ifndef GISMO_WITH_ADIFF
        copy_n(u.col(p).data(), expr.dim, expr.vars);
#       endif

        for (short_t c = 0; c!= n; ++c) // for all components
        {
#           ifdef GISMO_WITH_ADIFF
            for (index_t v = 0; v!=expr.dim; ++v)
                expr.vars[v].setVariable(v, expr.dim, u(v,p) );
            (*res)(c,p) = expr.expression[c].value().getHessian()(k,j); //note: can fail
#           else
            (*res)(c,p) =
                mixed_derivative<T>( expr.expression[c], expr.vars[k], expr.vars[j], 0.00001 ) ;
#           endif
        }
    }
    return res;
}

template<typename T>
gsMatrix<T> gsFunctionExpr<T>::laplacian(const gsMatrix<T>& u) const
{
    //gsDebug<< "Using finite differences (gsFunction::laplacian) for Laplacian.\n";
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point size.");
    const short_t n = targetDim();
    gsMatrix<T> res(n,u.cols());

    const PrivateData_t & expr =
#   ifdef _OPENMP
        omp_in_parallel() ? PrivateData_t(*my) :
#   endif
        *my;

    for( index_t p = 0; p != res.cols(); ++p )
    {
#       ifndef GISMO_WITH_ADIFF
        copy_n(u.col(p).data(), expr.dim, expr.vars);
#       endif

        for (short_t c = 0; c!= n; ++c) // for all components
        {
#           ifdef GISMO_WITH_ADIFF
            for (index_t v = 0; v!=expr.dim; ++v)
                expr.vars[v].setVariable(v, expr.dim, u(v,p) );
            res(c,p) = expr.expression[c].value().getHessian().trace();
#           else
            T & val = res(c,p);
            for ( index_t j = 0; j!=expr.dim; ++j )
                val += exprtk::
                    second_derivative<T>( expr.expression[c], expr.vars[j], 0.00001 );
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
    static std::string tag ()  { return "Function"; }
    static std::string type () { return "FunctionExpr"; }

    GSXML_GET_POINTER(Object);

    static void get_into (gsXmlNode * node, Object & obj)
    {
        getFunctionFromXml(node, obj);
    }

    static gsXmlNode * put (const Object & obj,
                            gsXmlTree & data )
    {
        gsXmlNode * func = makeNode("FunctionExpr", data);
        func->append_attribute(makeAttribute("dim", obj.domainDim(), data));

        const short_t tdim = obj.targetDim();

        if ( tdim == 1)
        {
            func->value( makeValue(obj.expression(), data) );
        }
        else
        {
            gsXmlNode * cnode;
            for (short_t c = 0; c!=tdim; ++c)
            {
                cnode = makeNode("c", obj.expression(c), data);
                func->append_node(cnode);
            }
        }

        return func;
    }
};

} // internal

}; // namespace gismo
