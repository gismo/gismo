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
#define exprtk_disable_enhanced_features

// This define will disable all string processing capabilities. Any
// expression that contains a string or string related syntax will result
// in a compilation failure.
#define exprtk_disable_string_capabilities

//#define GISMO_USE_AUTODIFF
#ifdef GISMO_USE_AUTODIFF
  /* Optional automatic differentiation */
  #define DScalar ad::DScalar2<real_t,-1>
  #include <exprtk_ad_adaptor.hpp> // external file
#else
  #include <exprtk.hpp>            // external file
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

// replaces appeareances of \a oldStr with \a newStr inside the string
// \a str
inline void stringReplace(std::string& str, 
                          const std::string& oldStr, 
                          const std::string& newStr)
{
    size_t pos = 0;
    while((pos = str.find(oldStr, pos)) != std::string::npos)
    {
        str.replace(pos, oldStr.length(), newStr);
        pos += newStr.length();
    }
}

} //namespace

namespace gismo
{

template<typename T> class gsFunctionExprPrivate 
{
public:

#ifdef GISMO_USE_AUTODIFF
    typedef DScalar Numeric_t;
#else
    typedef T Numeric_t;
#endif

    typedef exprtk::symbol_table<Numeric_t>  SymbolTable_t;
    typedef exprtk::expression<Numeric_t>    Expression_t;
    typedef exprtk::parser<Numeric_t>        Parser_t;

public:
    Numeric_t                 vars[6];
    SymbolTable_t             symbol_table;
    std::vector<Expression_t> expression;
    std::vector<std::string>  string; 
    index_t dim;
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
gsFunctionExpr<T>::gsFunctionExpr() : my(new PrivateData_t) 
{ }

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::string & expression_string, int ddim)
: my(new PrivateData_t)
{
    init(ddim);
    addComponent(expression_string);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::string & expression_string1, 
                                  const std::string & expression_string2,
                                  int ddim) 
: my(new PrivateData_t)
{
    init(ddim);
    addComponent(expression_string1);
    addComponent(expression_string2);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::string & expression_string1, 
                                  const std::string & expression_string2,
                                  const std::string & expression_string3,
                                  int ddim) 
: my(new PrivateData_t)
{
    init(ddim);
    addComponent(expression_string1);
    addComponent(expression_string2);
    addComponent(expression_string3);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::string & expression_string1,
                   const std::string & expression_string2,
                   const std::string & expression_string3,
                   const std::string & expression_string4,
                   int ddim)
: my(new PrivateData_t)
{
    init(ddim);
    addComponent(expression_string1);
    addComponent(expression_string2);
    addComponent(expression_string3);
    addComponent(expression_string4);
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
                   int ddim)
: my(new PrivateData_t)
{
    init(ddim);
    addComponent(expression_string1);
    addComponent(expression_string2);
    addComponent(expression_string3);
    addComponent(expression_string4);
    addComponent(expression_string5);
    addComponent(expression_string6);
    addComponent(expression_string7);
    addComponent(expression_string8);
    addComponent(expression_string9);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::vector<std::string> & expression_string, 
                                  int ddim)
: my(new PrivateData_t)
{
    init(ddim);
    for (std::size_t i = 0; i!= expression_string.size(); ++i)
        addComponent(expression_string[i]);
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const gsFunctionExpr & other)
{
    my = new PrivateData_t;
    init(other.my->dim);
    copyRange(other.my->vars, my->vars, 6);
    for (std::size_t i = 0; i!= other.my->string.size(); ++i)
        addComponent(other.my->string[i]);

    //my->symbol_table = other.my->symbol_table;
    //my->expression   = other.my->expression;
    //my->string       = other.my->string;
    //my->dim          = other.my->dim;
}

template<typename T>
gsFunctionExpr<T> & gsFunctionExpr<T>::operator=(gsFunctionExpr other)
{
    std::swap(my,other.my);
    return *this;
}

template<typename T>
gsFunctionExpr<T>::~gsFunctionExpr()
{ 
    delete my;
}


// template<typename T>
// void gsFunctionExpr<T>::init(variables,constants)


template<typename T>
void gsFunctionExpr<T>::init(const int dim)
{
    my->dim = dim;
    GISMO_ASSERT ( dim < 7, "The number of variables can be at most 6 (x,y,z,u,v,w)." );
    
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
void gsFunctionExpr<T>::addComponent(const std::string & strExpression)
{ 
    typedef typename gsFunctionExprPrivate<T>::Expression_t Expression_t;
    typedef typename gsFunctionExprPrivate<T>::Parser_t     Parser_t;

    // String
    my->string.push_back( strExpression );// Keep string data
    std::string & str = my->string.back();
    str.erase(std::remove(str.begin(), str.end(),' '), str.end() );
    stringReplace(str, "**", "^");

    // String expression
    my->expression.push_back(Expression_t());
    Expression_t & expr = my->expression.back();
    //expr.release();
    expr.register_symbol_table(my->symbol_table);

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
    for (std::size_t i = 0; i != symbol_list.size(); ++i)
    {
        symbol_t& symbol = symbol_list[i];
        // do something
    }
//*/

}


template<typename T>
int gsFunctionExpr<T>::domainDim() const
{ 
    return my->dim;
}

template<typename T>
int gsFunctionExpr<T>::targetDim() const
{ 
    return static_cast<int>(my->string.size());
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

/*
// to do: remove
template<typename T>
void gsFunctionExpr<T>::set_x_der (std::string expression_string)
{ 
    my->der_exp[0].register_symbol_table(my->symbol_table);
    exprtk::parser<T> parser;
    bool success = parser.compile(expression_string, my->der_exp[0] );        
    if ( ! success )
        std::cout<<"gsFunctionExpr set_x_der(.) error: " <<parser.error() <<std::endl;
    else
        my->der_flag[0]=true;
}

template<typename T>
void gsFunctionExpr<T>::set_y_der (std::string expression_string)
{ 
    my->der_exp[1].register_symbol_table(my->symbol_table);
    exprtk::parser<T> parser;
    bool success = parser.compile(expression_string, my->der_exp[1] );        
    if ( ! success )
        std::cout<<"gsFunctionExpr set_x_der(.) error: " <<parser.error() <<std::endl;
    else
        my->der_flag[1]=true;
}
//*/

template<typename T>
void gsFunctionExpr<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")");

    const int n = targetDim();
    result.resize(n, u.cols());

    for ( int p = 0; p!=u.cols(); p++ ) // for all evaluation points
    {
        copyRange(u.col(p).data(), my->vars, my->dim);

        for ( int j = 0; j!=my->dim; ++j )
            my->vars[j] =u(j,p);

        for (int c = 0; c!= n; ++c) // for all components
#           ifdef GISMO_USE_AUTODIFF
            result(c,p) = my->expression[c].value().getValue();
#           else
            result(c,p) = my->expression[c].value();
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
    for ( int p = 0; p!=u.cols(); ++p )
    {
        copyRange(u.col(p).data(), my->vars, my->dim);

#           ifdef GISMO_USE_AUTODIFF
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
    const index_t d = domainDim();
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")");

    const int n = targetDim();
    result.resize(d*n, u.cols());
    
    for ( int p = 0; p!=u.cols(); p++ ) // for all evaluation points
    {
#       ifdef GISMO_USE_AUTODIFF
        for (index_t k = 0; k!=d; ++k)
            my->vars[k].setVariable(k,d,u(k,p));
        for (int c = 0; c!= n; ++c) // for all components
            result.block(c*d,p,d,1) = my->expression[c].value().getGradient();
#       else
        copyRange(u.col(p).data(), my->vars, my->dim);
        for (int c = 0; c!= n; ++c) // for all components
            for ( int j = 0; j!=d; j++ ) // for all variables
                result(c*d + j, p) = 
                    exprtk::derivative<T>(my->expression[c], my->vars[j], 0.00001 ) ;
#       endif
    }
}

template<typename T>
void gsFunctionExpr<T>::deriv2_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    const index_t d = domainDim();
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")");
    
    const int n = targetDim();
    const unsigned stride = d + d*(d-1)/2;
    result.resize(stride*n, u.cols() );
    
    for ( int p = 0; p!=u.cols(); p++ ) // for all evaluation points
    {
#       ifndef GISMO_USE_AUTODIFF
        copyRange(u.col(p).data(), my->vars, my->dim);
#       endif

        for (int c = 0; c!= n; ++c) // for all components
        {
#           ifdef GISMO_USE_AUTODIFF
            for (index_t v = 0; v!=d; ++v)
                my->vars[v].setVariable(v,d,u(v,p));
            const DScalar::Hessian_t & Hmat = my->expression[c].value().getHessian();

            for ( index_t k=0; k!=d; ++k)
            {
                result(k,p) = Hmat(k,k);
                index_t m = d;
                for ( index_t l=k+1; l<d; ++l)
                    result(m++,p) = Hmat(k,l);
            }
#           else
            for (index_t k = 0; k!=d; ++k)
            {
                // H_{k,k}
                result(k,p) = exprtk::
                    second_derivative<T>(my->expression[c], my->vars[k], 0.00001);
                
                index_t m = d;
                for (index_t l=k+1; l<d; ++l)
                {
                    // H_{k,l}
                    result(m++,p) =
                        mixed_derivative<T>( my->expression[c], my->vars[k], 
                                             my->vars[l], 0.00001 );
                }
            }
#           endif
        }
    }
}

template<typename T>
typename gsFunction<T>::uMatrixPtr
gsFunctionExpr<T>::hess(const gsMatrix<T>& u, unsigned coord) const 
{ 
    //gsDebug<< "Using finite differences (gsFunctionExpr::hess) for Hessian.\n";
    GISMO_ENSURE(coord == 0, "Error, function is real");
    GISMO_ASSERT ( u.cols() == 1, "Need a single evaluation point." );
    const index_t d = u.rows();
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point dimension (expected: "
                   << my->dim <<", got "<< u.rows() <<")");
    
    gsMatrix<T> * res = new gsMatrix<T>(d,d);

#   ifdef GISMO_USE_AUTODIFF
    for (index_t v = 0; v!=d; ++v)
        my->vars[v].setVariable(v, d, u(v,0) );
    *res = my->expression[coord].value().getHessian();
#   else
    copyRange(u.data(), my->vars, d);    
    for( int j=0; j!=d; ++j )
    {
        (*res)(j,j) = exprtk::
            second_derivative<T>( my->expression[coord], my->vars[j], 0.00001);

        for( int k = 0; k!=j; ++k )
            (*res)(k,j) = (*res)(j,k) =
                mixed_derivative<T>( my->expression[coord], my->vars[k], 
                                     my->vars[j], 0.00001 );
    }
#       endif

    return typename gsFunction<T>::uMatrixPtr(res); 
}

template<typename T>
gsMatrix<T> * gsFunctionExpr<T>::mderiv(const gsMatrix<T> & u, 
                                        const index_t k, 
                                        const index_t j) const 
{    
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point size.");
    const int n = targetDim();
    gsMatrix<T> * res= new gsMatrix<T>(n,u.cols()) ;
    
    for( index_t p=0; p!=res->cols(); ++p )
    {
#       ifndef GISMO_USE_AUTODIFF
        copyRange(u.col(p).data(), my->vars, my->dim);
#       endif

        for (int c = 0; c!= n; ++c) // for all components
        {
#           ifdef GISMO_USE_AUTODIFF
            for (index_t v = 0; v!=my->dim; ++v)
                my->vars[v].setVariable(v, my->dim, u(v,p) );
            my->expression[c].value().getHessian();
            const DScalar::Hessian_t & Hmat = my->expression[c].value().getHessian();
            (*res)(c,p) = Hmat(k,j);
#           else
            (*res)(c,p) =
                mixed_derivative<T>( my->expression[c], my->vars[k], my->vars[j], 0.00001 ) ;
#           endif
        }
    }
    return res; 
}

template<typename T>
gsMatrix<T> * gsFunctionExpr<T>::laplacian(const gsMatrix<T>& u) const
{
    //gsDebug<< "Using finite differences (gsFunction::laplacian) for Laplacian.\n";
    GISMO_ASSERT ( u.rows() == my->dim, "Inconsistent point size.");
    const int n = targetDim();
    gsMatrix<T> * res= new gsMatrix<T>(n,u.cols()) ;
    
    for( index_t p = 0; p != res->cols(); ++p )
    {
#       ifndef GISMO_USE_AUTODIFF
        copyRange(u.col(p).data(), my->vars, my->dim);
#       endif

        for (int c = 0; c!= n; ++c) // for all components
        {
#           ifdef GISMO_USE_AUTODIFF
            for (index_t v = 0; v!=my->dim; ++v)
                my->vars[v].setVariable(v, my->dim, u(v,p) );
            (*res)(c,p) = my->expression[c].value().getHessian().trace();
#           else
            T & val = (*res)(c,p);
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
        os << my->string[0];

    for (int k = 1; k!= targetDim(); ++k)
        os <<", " << my->string[k];
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
    static std::string type () { return "expr"; }

    GSXML_GET_POINTER(Object);

    static void get_into (gsXmlNode * node, Object & obj)
    {
        getFunctionFromXml(node, obj);
    }
    
    static gsXmlNode * put (const Object & obj, 
                            gsXmlTree & data )
    {
        gsXmlNode * func = makeNode("Function", data);
        func->append_attribute(makeAttribute("dim", obj.domainDim(), data));

        const int tdim = obj.targetDim();
        
        if ( tdim == 1)
        {
            func->value( makeValue(obj.expression(), data) );
        }
        else
        {
            gsXmlNode * cnode;
            for (int c = 0; c!=tdim; ++c)
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
