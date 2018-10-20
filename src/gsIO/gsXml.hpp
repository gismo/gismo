/** @file gsXml.hpp

    @brief Provides implementation of XML helper functions.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <sstream>
#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsFunctionExpr.h>

namespace gismo {

namespace internal {

/*
template<class Object>
std::string gsXml<Object>::tag ()
{
    // Next line will produce compile-time error
    // if gsXml is not specialized for Object
    Object::Object_does_not_exist_ERROR;
    return "";
}
*/

template<class T>
gsXmlNode * makeNode( const std::string & name,
                      const gsMatrix<T> & value, gsXmlTree & data,
                      bool transposed)
{
    return data.allocate_node(rapidxml::node_element ,
                              data.allocate_string(name.c_str() ),
                              makeValue(value,data,transposed) );
}

template<class T>
char * makeValue(const gsMatrix<T> & value, gsXmlTree & data,
                 bool transposed)
{
    std::ostringstream oss;
    // Set precision
    oss << std::setprecision(data.getFloatPrecision());

    // Read/Write is RowMajor
    if ( transposed )
        for ( index_t j = 0; j< value.cols(); ++j)
        {
            for ( index_t i = 0; i< value.rows(); ++i)
                oss << value(i,j) <<" ";
            oss << "\n";
        }
    else
        for ( index_t i = 0; i< value.rows(); ++i)
        {
            for ( index_t j = 0; j< value.cols(); ++j)
                oss << value(i,j) <<" ";
            oss << "\n";
        }

    return data.allocate_string( oss.str().c_str() );
}

template<class T>
void getFunctionFromXml ( gsXmlNode * node, gsFunctionExpr<T> & result )
{
    //gsWarn<<"Reading "<< node->name() <<" function\n";

    GISMO_ASSERT( node->first_attribute("dim"), "getFunctionFromXml: No dim found" ) ;
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

    result = gsFunctionExpr<T>( expr_strings, d );
}

template<class T>
void getMatrixFromXml ( gsXmlNode * node, unsigned const & rows,
                        unsigned const & cols, gsMatrix<T> & result )
{
    //gsWarn<<"Reading "<< node->name() <<" matrix of size "<<rows<<"x"<<cols<<"Geometry..\n";
    std::istringstream str;
    str.str( node->value() );
    result.resize(rows,cols);

    for (unsigned i=0; i<rows; ++i) // Read is RowMajor
        for (unsigned j=0; j<cols; ++j)
            //if ( !(str >> result(i,j) ) )
              if (! gsGetValue(str,result(i,j)) )
            {
                gsWarn<<"XML Warning: Reading matrix of size "<<rows<<"x"<<cols<<" failed.\n";
                gsWarn<<"Tag: "<< node->name() <<", Matrix entry: ("<<i<<", "<<j<<").\n";
                return;
            }
}

template<class T>
gsXmlNode * putMatrixToXml ( gsMatrix<T> const & mat, gsXmlTree & data, std::string name)
{
    // Create XML tree node
    gsXmlNode* new_node = internal::makeNode(name, mat, data);
    return new_node;
}

template<class T>
gsXmlNode * putSparseMatrixToXml ( gsSparseMatrix<T> const & mat,
                                   gsXmlTree & data, std::string name)
{
    typedef typename gsSparseMatrix<T>::InnerIterator cIter;

    std::ostringstream str;
    str << std::setprecision(data.getFloatPrecision());
    const index_t nCol = mat.cols();

    for (index_t j=0; j != nCol; ++j) // for all columns
        for ( cIter it(mat,j); it; ++it ) // for all non-zeros in column
        {
            // Write the matrix entry
            str <<it.index() <<" "<<j<<" " << it.value() << "\n";
        }

    // Create XML tree node
    gsXmlNode* new_node = internal::makeNode(name, str.str(), data);
    return new_node;
}

template<class T>
void getSparseEntriesFromXml ( gsXmlNode * node,
                              gsSparseEntries<T> & result )
{
    result.clear();

    std::istringstream str;
    str.str( node->value() );
    index_t r,c;
    T val;

    //while( (str >> r) && (str >> c) && (str >> val) )
    while( (str >> r) && (str >> c) && ( gsGetValue(str,val)) )
        result.add(r,c,val);
}

}// end namespace internal

}// end namespace gismo
