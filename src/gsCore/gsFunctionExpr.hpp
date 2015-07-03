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

#include <exprtk.hpp>// External file

namespace gismo
{

template<typename T> class gsFunctionExprPrivate 
{
public:
    T vars[6];
    exprtk::symbol_table<T> symbol_table;
    exprtk::expression<T> expression;
    exprtk::expression<T> der_exp[6];
    bool der_flag[6];
    std::string  string; 
    int dim;
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
gsFunctionExpr<T>::gsFunctionExpr() : my(new gsFunctionExprPrivate<T>) {}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const std::string & expression_string, int ddim)
: my(new gsFunctionExprPrivate<T>)
{
    // Keep string data
    my->string = expression_string;
    my->string.erase(std::remove(my->string.begin(),my->string.end(),' '),my->string.end());
    stringReplace(my->string, "**", "^");
    my->dim = ddim;
    init();
}

template<typename T>
gsFunctionExpr<T>::gsFunctionExpr(const gsFunctionExpr& other)
{
    my = new gsFunctionExprPrivate<T>;
    my->string = other.my->string;
    my->dim = other.my->dim;
    init();
}

template<typename T>
gsFunctionExpr<T>& gsFunctionExpr<T>::operator=(const gsFunctionExpr& other)
{
    if (this != &other)
    {
        my->string = other.my->string;
        my->dim = other.my->dim;
        init();
    }
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
void gsFunctionExpr<T>::init()
{
    my->symbol_table.clear();
    my->expression.release();

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
    my->expression.register_symbol_table(my->symbol_table);

    exprtk::parser<T> parser;
    //parser.cache_symbols() = true; 
    bool success = parser.compile(my->string, my->expression);        
    if ( ! success )
        std::cout<<"gsFunctionExpr error: " <<parser.error() <<std::endl;

/*
    AM: Changed in recent versions.
    std::vector<std::string> varlist;
    parser.expression_symbols(varlist);
    varlist.erase(std::remove(varlist.begin(), varlist.end(), "pi"), varlist.end());
    my->dim = varlist.size();
*/

    my->der_flag[0]=false;
    my->der_flag[1]=false;
    my->der_flag[2]=false;
    my->der_flag[3]=false;
    my->der_flag[4]=false;
    my->der_flag[5]=false;
}


template<typename T>
int gsFunctionExpr<T>::domainDim() const
{ 
    return my->dim;
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

template<typename T>
void gsFunctionExpr<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT ( u.rows() < 7 && u.rows() > 0, "Inconsistent point size." );
    
    result.resize(1, u.cols());
    for ( int i = 0; i<u.cols(); i++ )
    {
        for ( int j = 0; j<u.rows(); j++ )
        {
            my->vars[j] =u(j,i);
        }
        result(0,i) = my->expression.value();
    }
}

template<typename T>
void gsFunctionExpr<T>::eval_component_into(const gsMatrix<T>& u, const index_t comp, gsMatrix<T>& result) const
{
    GISMO_ASSERT (comp == 0, "Given component number is too high. Only one component");
    result.resize(1, u.cols());
    for ( int i = 0; i<u.cols(); i++ )
    {
        for ( int j = 0; j<u.rows(); j++ )
        {
            my->vars[j] =u(j,i);
        }
        result(0,i) = my->expression.value();
    }
}

template<typename T>
void gsFunctionExpr<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    int n = u.rows();
    GISMO_ASSERT ( u.rows() < 7 && u.rows() > 0, "Inconsistent point size." );

    result.resize(1, n * u.cols());

    for( int i=0; i < u.cols(); ++i )
    {
        for ( int j = 0; j<n; j++ ) 
            my->vars[j] =u(j,i);
        for( int j=0; j< n; ++j )
            if ( my->der_flag[j] )
                result(n*i+j) = my->der_exp[j].value();
            else
                result(n*i+j) = exprtk::derivative<T>( my->expression, my->vars[j], 0.00001 ) ;
    }
}

template<typename T>
typename gsFunction<T>::uMatrixPtr
gsFunctionExpr<T>::hess(const gsMatrix<T>& u, unsigned coord) const 
{ 
    GISMO_ENSURE(coord == 0, "Error, function is real");
    // one column only..
    int n = u.rows();
    GISMO_ASSERT ( u.rows() < 7 && u.rows() > 0, "Inconsistent point size." );
    
    gsMatrix<T> * res= new gsMatrix<T>(n,n) ;
    
    for ( int j = 0; j<n; j++ ) 
        my->vars[j] =u(j,0);
    for( int j=0; j< n; ++j )
    {
        (*res)(j,j) = exprtk::second_derivative<T>( my->expression, my->vars[j], 0.00001 ) ;
        for( int k=0; k<j; ++k )
            (*res)(k,j) = (*res)(j,k) =
                exprtk::mixed_derivative<T>( my->expression, my->vars[k], my->vars[j], 0.00001 ) ;
    }
    return typename gsFunction<T>::uMatrixPtr(res); 
} ;
template<typename T>
gsMatrix<T> * gsFunctionExpr<T>::mderiv(const gsMatrix<T>& u, const index_t &k, const index_t &j ) const 
{    
    gsMatrix<T> * res= new gsMatrix<T>(1,u.cols()) ;
    
    for( index_t i=0; i< res->cols(); ++i )
    {
      	for ( int t = 0; t<u.rows(); t++ )
            my->vars[t] =u(t,i);
        (*res)(0,i) =
            exprtk::mixed_derivative<T>( my->expression, my->vars[k], my->vars[j], 0.00001 ) ;
    }
    return  res; 
} ;

template<typename T>
gsMatrix<T> * gsFunctionExpr<T>::laplacian(const gsMatrix<T>& u) const
{
    gsMatrix<T> * res= new gsMatrix<T>(1,u.cols()) ;
    res->setZero();
    int n = u.rows();
    
    for( index_t i=0; i< res->cols(); ++i )
        for ( int j = 0; j<n; j++ ) 
        { 
            my->vars[j] =u(j,i);
            (*res)(0,i) += 
                exprtk::second_derivative<T>( my->expression, my->vars[j], 0.00001 ) ;        
        }
    return  res; 
}

template<typename T>
T gsFunctionExpr<T>::value() const
{ 
    return my->expression.value() ; 
}

template<typename T>
std::ostream & gsFunctionExpr<T>::print(std::ostream &os) const
{ os << my->string ; return os; }



}; // namespace gismo
