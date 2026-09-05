/** @file gsMatrixXml.hpp

    @brief XML serialization for the gsMatrix types (implementation
    header).

    Moved here from gsIO/gsXmlUtils.hpp: each module owns its gsXml
    specializations (modularization stream S3, step A6).

    ONLY include this file from gsMatrixXml_.cpp (lib mode) or through
    gsMatrixXml_.cpp in header-only mode - a second inclusion in another
    translation unit of the library would poison the symbol visibility
    of the explicit instantiations.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>
#include <gsCore/gsLinearAlgebra.h>

namespace gismo {
namespace internal {

template<class T>
class gsXml< gsMatrix<T> >
{
private:
    gsXml() { }
    typedef gsMatrix<T> Object;

public:
    GSXML_COMMON_FUNCTIONS(Object);
    static std::string tag () { return "Matrix"; }
    static std::string type() { return ""; }

    GSXML_GET_POINTER(Object);

    static void get_into (gsXmlNode * node, Object & obj)
    {
        GISMO_ASSERT( !strcmp( node->name(),"Matrix"),
                      "Something went wrong. Expected Matrix tag." );

        unsigned rows = atoi(node->first_attribute("rows")->value());
        unsigned cols = atoi(node->first_attribute("cols")->value());
        gsXmlAttribute *format = node->first_attribute("format");
        std::string format_flag = format ? format->value() : "ascii";
        getMatrixFromXml<T>(node, rows, cols, obj, format_flag);
    }

    static gsXmlNode * put (const gsMatrix<T> & obj,
                            gsXmlTree & data )
    {
        gsXmlNode * mat_data = putMatrixToXml(obj,data);
        // Record matrix dimensions
        mat_data->append_attribute(
            makeAttribute("rows", obj.rows(), data) );
        mat_data->append_attribute(
            makeAttribute("cols", obj.cols(), data) );

        return mat_data;
    }
};

template<class T, int _Options, typename _Index>
class gsXml< gsSparseMatrix<T,_Options,_Index> >
{
private:
    gsXml() { }
    typedef gsSparseMatrix<T,_Options,_Index> Object;

public:
    GSXML_COMMON_FUNCTIONS(Object);
    static std::string tag () { return "SparseMatrix"; }
    static std::string type() { return ""; }

    GSXML_GET_POINTER(Object);

    static void get_into (gsXmlNode * node, Object & obj)
    {
        GISMO_ASSERT( !strcmp( node->name(),"SparseMatrix"),
                      "Something went wrong. Expected SparseMatrix tag." );

        const index_t rows  = atoi ( node->first_attribute("rows")->value() ) ;
        const index_t cols  = atoi ( node->first_attribute("cols")->value() ) ;

        gsSparseEntries<T> entries;
        getSparseEntriesFromXml<T>(node, entries);

        obj.resize(rows,cols);
        obj.setFrom(entries);
    }

    static gsXmlNode * put (const Object & obj,
                            gsXmlTree & data )
    {
        gsXmlNode * mat_data = putSparseMatrixToXml(obj,data);

        mat_data->append_attribute(
            makeAttribute("rows", obj.rows(), data) );
        mat_data->append_attribute(
            makeAttribute("cols", obj.cols(), data) );

        return mat_data;
    }
};
} // namespace internal
} // namespace gismo
