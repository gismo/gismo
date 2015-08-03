#pragma once

#include <string>
#include <vector>

#include <exprtk.hpp>// External file

namespace gismo
{

template<typename T> class gsMFunctionExprPrivate 
{
public:
  T vars[6]; // const members can change this
  exprtk::symbol_table<T> symbol_table;
  exprtk::expression<T> expression[9];
  std::string  string[9];
  int domainDim;
  int targetDim;
};

template<typename T>
gsMFunctionExpr<T>::gsMFunctionExpr() : my(new gsMFunctionExprPrivate<T>) {};

template<typename T>
gsMFunctionExpr<T>::gsMFunctionExpr(const std::string & expression_string1,
                                    int ddim) :
    my(new gsMFunctionExprPrivate<T>)
{
    my->domainDim = ddim;
    my->targetDim = 1;
    // Copy string data
        my->string[0] = expression_string1;
for (int i = 0; i!= my->targetDim; ++i)
      my->string[i].erase(std::remove(my->string[i].begin(),my->string[i].end(),' '),
                                        my->string[i].end());
    init();
}

template<typename T>
gsMFunctionExpr<T>::gsMFunctionExpr(const std::string & expression_string1, 
                                    const std::string & expression_string2,
                                    int ddim) : 
    my(new gsMFunctionExprPrivate<T>)
{
    my->domainDim = ddim;
    my->targetDim = 2;
    // Copy string data
        my->string[0] = expression_string1;
        my->string[1] = expression_string2;
for (int i = 0; i!= my->targetDim; ++i)
      my->string[i].erase(std::remove(my->string[i].begin(),my->string[i].end(),' '),
                                        my->string[i].end());
    init();
}

template<typename T>
gsMFunctionExpr<T>::gsMFunctionExpr(const std::string & expression_string1, 
                                    const std::string & expression_string2,
                                    const std::string & expression_string3,
                                    int ddim) : 
    my(new gsMFunctionExprPrivate<T>)
{
    my->domainDim = ddim;
    my->targetDim = 3;
    // Copy string data
        my->string[0] = expression_string1;
        my->string[1] = expression_string2;
        my->string[2] = expression_string3;
for (int i = 0; i!= my->targetDim; ++i)
      my->string[i].erase(std::remove(my->string[i].begin(),my->string[i].end(),' '),
                                        my->string[i].end());
    init();
}

template<typename T>
gsMFunctionExpr<T>::gsMFunctionExpr(const std::string & expression_string1,
                                    const std::string & expression_string2,
                                    const std::string & expression_string3,
                                    const std::string & expression_string4,
                                    int ddim) :
    my(new gsMFunctionExprPrivate<T>)
{
    my->domainDim = ddim;
    my->targetDim = 4;
    // Copy string data
        my->string[0] = expression_string1;
        my->string[1] = expression_string2;
        my->string[2] = expression_string3;
        my->string[3] = expression_string4;
for (int i = 0; i!= my->targetDim; ++i)
      my->string[i].erase(std::remove(my->string[i].begin(),my->string[i].end(),' '),
                                        my->string[i].end());
    init();
}

template<typename T>
gsMFunctionExpr<T>::gsMFunctionExpr(const std::string & expression_string1,
                                    const std::string & expression_string2,
                                    const std::string & expression_string3,
                                    const std::string & expression_string4,
                                    const std::string & expression_string5,
                                    const std::string & expression_string6,
                                    int ddim) :
    my(new gsMFunctionExprPrivate<T>)
{
    my->domainDim = ddim;
    my->targetDim = 6;
    // Copy string data
        my->string[0] = expression_string1;
        my->string[1] = expression_string2;
        my->string[2] = expression_string3;
        my->string[3] = expression_string4;
        my->string[4] = expression_string5;
        my->string[5] = expression_string6;
for (int i = 0; i!= my->targetDim; ++i)
      my->string[i].erase(std::remove(my->string[i].begin(),my->string[i].end(),' '),
                                        my->string[i].end());
    init();
}

template<typename T>
gsMFunctionExpr<T>::gsMFunctionExpr(const std::string & expression_string1,
                                    const std::string & expression_string2,
                                    const std::string & expression_string3,
                                    const std::string & expression_string4,
                                    const std::string & expression_string5,
                                    const std::string & expression_string6,
                                    const std::string & expression_string7,
                                    const std::string & expression_string8,
                                    const std::string & expression_string9,
                                    int ddim) :
    my(new gsMFunctionExprPrivate<T>)
{
    my->domainDim = ddim;
    my->targetDim = 9;
    // Copy string data
        my->string[0] = expression_string1;
        my->string[1] = expression_string2;
        my->string[2] = expression_string3;
        my->string[3] = expression_string4;
        my->string[4] = expression_string5;
        my->string[5] = expression_string6;
        my->string[6] = expression_string7;
        my->string[7] = expression_string8;
        my->string[8] = expression_string9;
for (int i = 0; i!= my->targetDim; ++i)
      my->string[i].erase(std::remove(my->string[i].begin(),my->string[i].end(),' '),
                                        my->string[i].end());
    init();
}

template<typename T>
gsMFunctionExpr<T>::gsMFunctionExpr(const std::vector<std::string> & expression_string, 
                                    int ddim) : 
    my(new gsMFunctionExprPrivate<T>)
{
    my->domainDim = ddim;
    my->targetDim = expression_string.size();
    GISMO_ASSERT(my->targetDim<7, "Target Dimension upto 6 supported");
    // Copy string data
    for (int i = 0; i!= my->targetDim; ++i)
    {
        my->string[i] = expression_string[i];
        my->string[i].erase(std::remove(my->string[i].begin(),my->string[i].end(),' '),
                                        my->string[i].end());
    }
    init();
}

template<typename T>
gsMFunctionExpr<T>::gsMFunctionExpr(const gsMFunctionExpr& other) : my(new gsMFunctionExprPrivate<T>)
{
    *my = *other.my;
    init();
}

template<typename T>
gsMFunctionExpr<T>& gsMFunctionExpr<T>::operator=(const gsMFunctionExpr& other)
{
    if (this != &other)
    {
        *my = *other.my;
        return *this;
    }
    return *this;
}

template<typename T>
gsMFunctionExpr<T>::~gsMFunctionExpr()
	{ 
	    delete my;
	}

template<typename T>
void gsMFunctionExpr<T>::init()
{
    my->symbol_table.clear();
    // Identify symbol table
    my->symbol_table.add_variable("x",my->vars[0]);
    my->symbol_table.add_variable("y",my->vars[1]);
    my->symbol_table.add_variable("z",my->vars[2]);
    my->symbol_table.add_variable("w",my->vars[3]);
    my->symbol_table.add_variable("u",my->vars[4]);
    my->symbol_table.add_variable("v",my->vars[5]);
    my->symbol_table.add_pi();
    //my->symbol_table.add_constant("C", 1);
    
    exprtk::parser<T> parser;
    for (int i = 0; i!= my->targetDim; ++i)
    {
        my->expression[i].release();
        my->expression[i].register_symbol_table( my->symbol_table);
        
        //parser.cache_symbols() = true; 
        bool success = parser.compile(my->string[i], my->expression[i]);        
        if ( ! success )
            gsInfo<<"gsMFunctionExpr error: " <<parser.error() <<std::endl;
    }
}


template<typename T>
int gsMFunctionExpr<T>::domainDim() const
{ 
   return my->domainDim;
}

template<typename T>
int gsMFunctionExpr<T>::targetDim() const
{ 
   return my->targetDim;
}

template<typename T>
void gsMFunctionExpr<T>::eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    GISMO_ASSERT ( u.rows() == my->domainDim, "Wrong point dimension " );
    
    result.resize(my->targetDim, u.cols());
    for ( int i = 0; i<u.cols(); i++ )
    {
        for ( int j = 0; j<my->domainDim; j++ )
            my->vars[j] =u(j,i);

        for (int k = 0; k!= my->targetDim; ++k)
            result(k,i) = my->expression[k].value();
    }
}

template<typename T>
void gsMFunctionExpr<T>::eval_component_into(const gsMatrix<T>& u, const index_t comp, gsMatrix<T>& result) const
{
    GISMO_ASSERT ( u.rows() == my->domainDim, "Wrong point dimension " );
    GISMO_ASSERT (comp < my->targetDim, 
                  "Given component number is higher then number of components");
    result.resize(1, u.cols());
    for ( int i = 0; i<u.cols(); i++ )
    {
        for ( int j = 0; j<my->domainDim; j++ )
            my->vars[j] =u(j,i);

        result(0,i) = my->expression[comp].value();
    }
}

template<typename T>
void gsMFunctionExpr<T>::deriv_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
{
    //gsDebug<< "Using finite differences for derivatives.\n";
    const int d = my->domainDim;
    GISMO_ASSERT ( u.rows() == d, "Wrong point dimension" );
    
    result.resize(d*my->targetDim, u.cols());
    
    for ( int p = 0; p!=u.cols(); p++ ) // for all evaluation points
    {
        for ( int j = 0; j<d; j++ )
            my->vars[j] =u(j,p);

        for ( int j = 0; j<d; j++ ) // for all variables
            for (int c = 0; c!= my->targetDim; ++c) // for all components
                result(c*d + j, p) = 
                    exprtk::derivative<T>(my->expression[c], my->vars[j], 0.00001 ) ;
    }
}

template<typename T>
T gsMFunctionExpr<T>::value(int k) const
{ 
    return my->expression[k].value() ; 
}

template<typename T>
std::ostream & gsMFunctionExpr<T>::print(std::ostream &os) const
{ 
    os <<"( ";
    os << my->string[0] ;
    for (int k = 1; k!= my->targetDim; ++k)
        os <<", " << my->string[k] ;
    os <<" )";
    return os;
}


}; // namespace gismo
