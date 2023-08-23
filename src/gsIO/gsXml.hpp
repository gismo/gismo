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
#include <gsIO/gsIOUtils.h>

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
    oss << std::setprecision(data.getFloatPrecision()) << "\n";

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

template <class T>
void getMatrixFromXml(gsXmlNode* node, unsigned const& rows,
                      unsigned const& cols, gsMatrix<T>& result,
                      const std::string& base_type_flag) {
    if (base_type_flag == "ASCII") {
        std::istringstream str;
        str.str(node->value());
        result.resize(rows, cols);
        for (unsigned i = 0; i < rows; ++i)  // Read is RowMajor
            for (unsigned j = 0; j < cols; ++j)
              // if ( !(str >> result(i,j) ) )
              if (!gsGetValue(str, result(i, j))) {
                gsWarn << "XML Warning: Reading matrix of size " << rows << "x"
                       << cols << " failed.\n";
                gsWarn << "Tag: " << node->name() << ", Matrix entry: (" << i
                       << ", " << j << ").\n";
                return;
              }
    } else {
        // Perform type checks (no integral to floting point conversion)
        if (std::is_integral<T>::value ^
            (base_type_flag.find("Int") != std::string::npos)) {
            GISMO_ERROR(
                "Conversions from integral to floating type and vice-versa is "
                "not allowed!");
        }
        // Hide the transformation in a polymorphic lambda function to avoid
        // duplicate code
        result.resize(rows, cols);
        auto copy_input_to_gsMatrix = [&result, &rows,
                                       &cols](const auto& input_array) -> void {
          // Size check
          if (input_array.size() != (rows * cols)) {
            GISMO_ERROR(
                "Input array has the wrong size or could not be converted");
          }
          // Converting into gsMatrix (manipulating directly on gsMatrix is
          // more efficient if T=InputType)
          result.resize(rows, cols);
          for (unsigned i = 0; i < rows; ++i) {
            for (unsigned j = 0; j < cols; ++j) {
              result(i, j) = static_cast<T>(input_array[i * cols + j]);
            }
          }
        };
        // Get string
        std::string b64_string(node->value());
        // Perform the actual input (using the proper encoding type)
        if (base_type_flag == "B64Uint16") {  // Unsigned int 16
            copy_input_to_gsMatrix(
                Base64::Decode<uint16_t>(b64_string));
        } else if (base_type_flag == "B64Uint32") {  // Unsigned int 32
            copy_input_to_gsMatrix(
                Base64::Decode<uint32_t>(b64_string));
        } else if (base_type_flag == "B64Uint64") {  // Unsigned int 64
            copy_input_to_gsMatrix(
                Base64::Decode<uint64_t>(b64_string));
        } else if (base_type_flag == "B64Int16") {  // Int 16
            copy_input_to_gsMatrix(
                Base64::Decode<int16_t>(b64_string));
        } else if (base_type_flag == "B64Int32") {  // Int 32
            copy_input_to_gsMatrix(
                Base64::Decode<int32_t>(b64_string));
        } else if (base_type_flag == "B64Int64") {  // Int 64
            copy_input_to_gsMatrix(
                Base64::Decode<int64_t>(b64_string));
        } else if (base_type_flag == "B64Float32") {  // Float 32
            copy_input_to_gsMatrix(
                Base64::Decode<float>(b64_string));
        } else if (base_type_flag == "B64Float64") {  // Float 64
            copy_input_to_gsMatrix(
                Base64::Decode<double>(b64_string));
        } else {
            GISMO_ERROR("Reading matrix from XML found unknown type");
        }
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
