/** @file gsJSON.h

    @brief Wrapper for the JSON library

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. M. Verhelst, J. Li
*/

#pragma once

#include <JSON/single_include/nlohmann/json.hpp>
#include <gsCore/gsGeometry.h>

namespace gismo
{

// Add the nlohmann::json namespace to the gismo namespace
using json = nlohmann::json;

/**
 * @brief Writes a gsVector to JSON
 * @param j JSON object
 * @param vec gsVector to be written
 */
template<class T>
void to_json(json& j, const gsVector<T> & vec)
{
    if (vec.rows()==0)
        j = nullptr;
    else
    {
        j["rows"] = vec.rows();
        std::vector<T> data(vec.data(), vec.data() + vec.size());
        j["data"] = data;
    }
}

/**
 * @brief Reads a gsVector from JSON
 * @param j JSON object
 * @param vec gsVector to be read
 */
template<class T>
void from_json(const json& j, gsVector<T>& vec)
{
    if (j.is_null())
        vec = gsVector<T>();
    else
    {
        GISMO_ASSERT(j.contains("rows"), "JSON object does not contain 'rows' field");
        GISMO_ASSERT(j.contains("data"), "JSON object does not contain 'data' field");
        GISMO_ASSERT(j["rows"].is_number_integer(),"rows is not an integer");
        GISMO_ASSERT(j["data"].is_array(),"data is not an array");
        vec.resize(j["rows"]);
        for (index_t I = 0; I < vec.rows(); I++)
            vec[I] = j["data"][I];
    }
}

/**
 * @brief Writes a gsMatrix to JSON
 * @param j JSON object
 * @param mat gsMatrix to be written
 */
template<class T>
void to_json(json& j, const gsMatrix<T> & mat)
{
    if (mat.rows()==0 || mat.cols()==0)
        j = nullptr;
    else
    {
        j["rows"] = mat.rows();
        j["cols"] = mat.cols();
        std::vector<T> data(mat.data(), mat.data() + mat.size());
        j["data"] = data;
    }
}

/**
 * @brief Reads a gsMatrix from JSON
 * @param j JSON object
 * @param mat gsMatrix to be read
 */
template<class T>
void from_json(const json& j, gsMatrix<T>& mat)
{
    if (j.is_null())
        mat = gsMatrix<T>();
    else
    {
        GISMO_ASSERT(j.contains("rows"), "JSON object does not contain 'rows' field");
        GISMO_ASSERT(j.contains("cols"), "JSON object does not contain 'cols' field");
        GISMO_ASSERT(j.contains("data"), "JSON object does not contain 'data' field");
        GISMO_ASSERT(j["rows"].is_number_integer(),"rows is not an integer");
        GISMO_ASSERT(j["cols"].is_number_integer(),"cols is not an integer");
        GISMO_ASSERT(j["data"].is_array(),"data is not an array");
        std::vector<T> data = j["data"].get<std::vector<T> >();
        mat = gsAsMatrix<T>(data,j["rows"], j["cols"]);
    }
}

/**
 * @brief Writes a gsSparseMatrix to JSON
 * @param j JSON object
 * @param mat gsSparseMatrix to be written
 */
template<class T>
void to_json(json &j, const gsSparseMatrix<T> & mat)
{
    if (mat.rows()==0 || mat.cols()==0)
        j = nullptr;
    else
    {
        j["rowIndices"] = json::array();
        j["colIndices"] = json::array();
        j["values"] = json::array();
        for (index_t I = 0; I < mat.outerSize(); I++)
            for (typename gsSparseMatrix<T>::InnerIterator it(mat,I); it; ++it)
            {
                j["rowIndices"].push_back(it.row());
                j["colIndices"].push_back(it.col());
                j["values"].push_back(it.value());
            }
    }
}

/**
 * @brief Reads a gsSparseMatrix from JSON
 * @param j JSON object
 * @param mat gsSparseMatrix to be read
 */
template<class T>
void from_json(const json &j, gsSparseMatrix<T> & mat)
{
    if (j.is_null())
        mat = gsSparseMatrix<T>();
    else
    {
        GISMO_ASSERT(j.contains("rowIndices"), "JSON object does not contain 'rowIndices' field");
        GISMO_ASSERT(j.contains("colIndices"), "JSON object does not contain 'colIndices' field");
        GISMO_ASSERT(j.contains("values"), "JSON object does not contain 'values' field");
        GISMO_ASSERT(j["rowIndices"].is_array(),"rowIndices is not an array");
        GISMO_ASSERT(j["colIndices"].is_array(),"colIndices is not an array");
        GISMO_ASSERT(j["values"].is_array(),"values is not an array");
        GISMO_ASSERT(j["rowIndices"].size() == j["colIndices"].size() && j["colIndices"].size() == j["values"].size(),
            "rowIndices, colIndices, and values must have the same size");
        mat.resize(j["rowIndices"].size(), j["colIndices"].size());
        for (index_t I = 0; I < j["rowIndices"].size(); I++)
            mat.insert(j["rowIndices"][I], j["colIndices"][I]) = j["values"][I];
    }
}

/**
 * @brief Writes a gsOptionList to JSON
 * @param j JSON object
 * @param opt gsOptionList to be written
 */
void to_json(json &j, const gsOptionList & opt)
{
    typedef gsOptionList::OptionListEntry Entry;
    std::vector<Entry> entries = opt.getAllEntries();
    for (typename std::vector<Entry>::const_iterator it = entries.begin(); it!=entries.end(); it++)
    {
        if (strcmp("int", it->type.c_str()) == 0)
            j[it->label] = std::stoi(it->val);
        else if (strcmp("real", it->type.c_str()) == 0)
            j[it->label] = std::stod(it->val);
        else if (strcmp("string", it->type.c_str()) == 0)
            j[it->label] = it->val;
        else if (strcmp("bool", it->type.c_str()) == 0)
            j[it->label] = (it->val == "true");
        else
            GISMO_ERROR("Type of "<<it->label<<" not recognized");
    }
}

/**
 * @brief Reads a gsOptionList from JSON
 * @param j JSON object
 * @param opt gsOptionList to be read
 */
void from_json(const json &j, gsOptionList & opt)
{
    GISMO_ASSERT(j.is_object(), "JSON object is not an object");
    // If the JSON object is null, create an empty gsOptionList
    if      (j.is_null())
        opt = gsOptionList();
    else
    {
        for (json::const_iterator it = j.begin(); it != j.end(); ++it)
        {
            if      (it->is_number_integer() ||
                    it->is_number_unsigned() )
                opt.addInt(it.key(), "", it.value().get<int>());
            else if (it->is_number_float())
                opt.addReal(it.key(), "", it.value().get<double>());
            else if (it->is_string())
                opt.addString(it.key(), "", it.value().get<std::string>());
            else if (it->is_boolean())
                opt.addSwitch(it.key(), "", it.value().get<bool>());
            else if (it->is_array())
            {
                index_t k=0;
                for (auto & val : it.value())
                {
                    if      (val.is_number_integer() ||
                                val.is_number_unsigned() )
                        opt.addInt(it.key() + "[" + std::to_string(k) + "]", "", val.get<int>());
                    else if (val.is_number_float())
                        opt.addReal(it.key() + "[" + std::to_string(k) + "]", "", val.get<double>());
                    else if (val.is_string())
                        opt.addString(it.key() + "[" + std::to_string(k) + "]", "", val.get<std::string>());
                    else if (val.is_boolean())
                        opt.addSwitch(it.key() + "[" + std::to_string(k) + "]", "", val.get<bool>());
                    else
                        GISMO_ERROR("Type of "<<it.key()<<" not recognized");
                    k++;
                }
            }
            else if (j.is_object())
            {
                gsOptionList sublist, sublist2;
                from_json(it.value(), sublist);
                sublist2 = sublist.wrapIntoGroup(it.key());
                opt.update(sublist2, gsOptionList::addIfUnknown);
            }
        }
    }
}

/**
 * @brief Writes a gsKnotVector to JSON
 * @param j JSON object
 * @param kv gsKnotVector to be written
 */
template<class T>
void to_json(json &j, const gsKnotVector<T> & kv)
{
    j["degree"] = kv.degree();

    std::vector<T> data(kv.data(), kv.data() + kv.size());
    j["knots"] = data;
}

/**
 * @brief Reads a gsKnotVector from JSON
 * @param j JSON object
 * @param kv gsKnotVector to be read
 */
template<class T>
void from_json(const json &j, gsKnotVector<T> & kv)
{
    GISMO_ASSERT(j.contains("knots"), "JSON object does not contain 'knots' field");
    GISMO_ASSERT(j.contains("degree"), "JSON object does not contain 'degree' field");
    GISMO_ASSERT(j["knots"].is_array(), "Field 'knots' is not an array");
    GISMO_ASSERT(j["degree"].is_number_integer(), "Field 'degree' is not an integer");
    std::vector<T> knots = j["knots"].get<std::vector<T> >();
    kv = gsKnotVector<T>(j["degree"], knots.begin(), knots.end());
}

/**
 * @brief Writes a gsBSplineBasis to JSON
 * @param j JSON object
 * @param basis gsBSplineBasis to be written
 */
template<class T>
void to_json(json &j, const gsBSplineBasis<T> & basis)
{
    j["type"] = "BSplineBasis";
    j["knots"] = basis.knots();
}

/**
 * @brief Reads a gsBSplineBasis from JSON
 * @param j JSON object
 * @param basis gsBSplineBasis to be read
 */
template<class T>
void from_json(const json &j, gsBSplineBasis<T> & basis)
{
    GISMO_ASSERT(j.contains("type"), "JSON object does not contain 'type' field");
    GISMO_ASSERT(j.contains("knots"), "JSON object does not contain 'knots' field");
    GISMO_ASSERT(j["type"].is_string(), "Field 'type' is not a string");
    GISMO_ASSERT(j["type"]=="BSplineBasis","Type of basis is not BSplineBasis");
    gsKnotVector<T> kv = j["knots"].get<gsKnotVector<T> >();
    basis = gsBSplineBasis<T>(kv);
}

/**
 * @brief Writes a gsTensorBSplineBasis to JSON
 * @param j JSON object
 * @param basis gsTensorBSplineBasis to be written
 */
template<short_t d, class T>
void to_json(json &j, const gsTensorBSplineBasis<d,T> & basis)
{
    j["type"] = "TensorBSplineBasis"+util::to_string(d);
    for (unsigned D=0; D!=d; D++)
        j["component"+util::to_string(D)] = basis.component(D);
}

/**
 * @brief Reads a gsTensorBSplineBasis from JSON
 * @param j JSON object
 * @param basis gsTensorBSplineBasis to be read
 */
template<short_t d, class T>
void from_json(const json &j, gsTensorBSplineBasis<d,T> & basis)
{
    GISMO_ASSERT(j.contains("type"), "JSON object does not contain 'type' field");
    GISMO_ASSERT(j.contains("component0"), "JSON object does not contain 'component0' field");
    GISMO_ASSERT(j["type"].is_string(), "Field 'type' is not a string");
    GISMO_ASSERT(j["type"]=="TensorBSplineBasis"+util::to_string(d),"Type of basis is not TensorBSplineBasis"+util::to_string(d));

    std::vector<gsKnotVector<T> > KVs(d);
    gsBSplineBasis<T> componentBasis;
    for (unsigned D=0; D!=d; D++)
    {
        GISMO_ASSERT(j.contains("component"+util::to_string(D)), "JSON object does not contain 'component"+util::to_string(D)+"' field");
        GISMO_ASSERT(j["component"+util::to_string(D)].is_object(), "Field 'component"+util::to_string(D)+"' is not an object");
        GISMO_ASSERT(j["component"+util::to_string(D)+""]["type"].is_string(), "Field 'component"+util::to_string(D)+"' does not contain 'type' string");
        GISMO_ASSERT(j["component"+util::to_string(D)+""]["type"]=="BSplineBasis","Type of component "+util::to_string(D)+" is not BSplineBasis");
        from_json(j["component"+util::to_string(D)], componentBasis);
        KVs[D] = componentBasis.knots();
    }

    basis = gsTensorBSplineBasis<d,T>(KVs);
}

/**
 * @brief Writes a gsTensorBSpline to JSON
 * @param j JSON object
 * @param geo gsTensorBSpline to be written
 */
template<class T>
void to_json(json &j, const gsBasis<T> & basis)
{
    if      ( const gsBSplineBasis<T> * b = dynamic_cast<const gsBSplineBasis<T> *>( &basis ) )
        to_json(j, *b);
    else if ( const gsTensorBSplineBasis<2,T> * b = dynamic_cast<const gsTensorBSplineBasis<2,T> *>( &basis ) )
        to_json(j, *b);
    else if ( const gsTensorBSplineBasis<3,T> * b = dynamic_cast<const gsTensorBSplineBasis<3,T> *>( &basis ) )
        to_json(j, *b);
    else if ( const gsTensorBSplineBasis<4,T> * b = dynamic_cast<const gsTensorBSplineBasis<4,T> *>( &basis ) )
        to_json(j, *b);
    else
        GISMO_ERROR("No known basis type");
}

/**
 * @brief Writes a gsGeometry to JSON
 * @param j JSON object
 * @param geo gsGeometry to be written
 */
template<class T>
void to_json(json &j, const gsGeometry<T> & geo)
{
    j["basis"] = geo.basis();
    j["coefs"] = geo.coefs();
}

/**
 * @brief Reads a gsTensorBSpline from JSON
 * @param j JSON object
 * @param geo gsTensorBSpline to be read
 */
template<class T>
void from_json(const json &j, gsBSpline<T> & geo)
{
    GISMO_ASSERT(j.contains("basis"), "JSON object does not contain 'basis' field\nj="<<j);
    GISMO_ASSERT(j.contains("coefs"), "JSON object does not contain 'coefs' field\nj="<<j);
    GISMO_ASSERT(j["basis"].is_object(), "Field 'basis' is not an object.\nj="<<j);
    gsBSplineBasis<T> basis = j["basis"].get<gsBSplineBasis<T> >();
    gsMatrix<T> coefs = j["coefs"].get<gsMatrix<T> >();
    geo = gsBSpline<T>(basis, coefs);
}

/**
 * @brief Reads a gsTensorBSpline from JSON
 * @param j JSON object
 * @param geo gsTensorBSpline to be read
 */
template<short_t d, class T>
void from_json(const json &j, gsTensorBSpline<d,T> & geo)
{
    gsTensorBSplineBasis<d,T> basis = j["basis"].get<gsTensorBSplineBasis<d,T> >();
    gsMatrix<T> coefs = j["coefs"].get<gsMatrix<T> >();
    geo = gsTensorBSpline<d,T>(basis, coefs);
}

/**
 * @brief Reads a gsGeometry from JSON
 * @param j JSON object
 * @return Unique pointer to gsGeometry
 */
template <class T>
typename gsGeometry<T>::uPtr get_Geometry(const json &j)
{
    if (j.contains("basis"))
    {
        if (j["basis"]["type"] == "BSplineBasis")
        {
            gsBSpline<T> bspline;
            from_json(j, bspline);
            return bspline.clone();
        }
        else if (j["basis"]["type"] == "TensorBSplineBasis2")
        {
            gsTensorBSpline<2,T> tensorBspline;
            from_json(j, tensorBspline);
            return tensorBspline.clone();
        }
        else if (j["basis"]["type"] == "TensorBSplineBasis3")
        {
            gsTensorBSpline<3,T> tensorBspline;
            from_json(j, tensorBspline);
            return tensorBspline.clone();
        }
        else if (j["basis"]["type"] == "TensorBSplineBasis4")
        {
            gsTensorBSpline<4,T> tensorBspline;
            from_json(j, tensorBspline);
            return tensorBspline.clone();
        }
        else
        {
            GISMO_ERROR("Unsupported basis type: " + j["basis"]["type"].get<std::string>());
        }
    }
    else if (j.contains("path"))
    {
        GISMO_ERROR("Path not implemented for geometry via JSON.");
        // std::string path = j["path"].get<std::string>();
        // GISMO_ASSERT(gsFileManager::fileExists(path),"File does not exist: " + path);
        // GISMO_ASSERT(gsFileManager::getExtension(path) == "xml","File extension is not .xml: " + path);
        // gsReadFile<T>(path, geo);
    }
    else if (j.contains("create"))
    {
        GISMO_ASSERT(j["create"].contains("type"), "Create JSON must contain 'type' field");
        if (j["create"]["type"] == "BSplineSquare")
        {
            T L = (j["create"].contains("length") ? j["create"]["length"].get<T>() : 1.0);
            T x = (j["create"].contains("x") ? j["create"]["x"].get<T>() : 0.0);
            T y = (j["create"].contains("y") ? j["create"]["y"].get<T>() : 0.0);
            return gsNurbsCreator<T>::BSplineSquare(L, x, y);
        }
        if (j["create"]["type"] == "BSplineCube")
        {
            T L = (j["create"].contains("length") ? j["create"]["length"].get<T>() : 1.0);
            T x = (j["create"].contains("x") ? j["create"]["x"].get<T>() : 0.0);
            T y = (j["create"].contains("y") ? j["create"]["y"].get<T>() : 0.0);
            T z = (j["create"].contains("z") ? j["create"]["z"].get<T>() : 0.0);
            return gsNurbsCreator<T>::BSplineCube(L, x, y, z);
        }
        if (j["create"]["type"] == "BSplineRectangle")
        {
            T xlow = (j["create"].contains("xlow") ? j["create"]["xlow"].get<T>() : 0.0);
            T ylow = (j["create"].contains("ylow") ? j["create"]["ylow"].get<T>() : 0.0);
            T xhigh= (j["create"].contains("xhigh") ? j["create"]["xhigh"].get<T>() : 1.0);
            T yhigh= (j["create"].contains("yhigh") ? j["create"]["yhigh"].get<T>() : 1.0);
            return gsNurbsCreator<T>::BSplineRectangle(xlow,ylow,xhigh,yhigh);
        }
    }
    else
    {
        GISMO_ERROR("JSON object must contain 'basis' or 'path' field");
    }
    return nullptr; // This should never be reached, but added to avoid compiler warnings
}

/**
 * @brief Reads a gsGeometry from JSON
 * @param j JSON object
 * @param geo gsGeometry to be read
 */
template <class T>
void from_json(const json &j, gsGeometry<T> & geo)
{
    geo = *get_Geometry<T>(j);
}
/**
 * @brief Convert boundary string to boundary::side enum
 * @param str String representation of boundary
 * @return boundary::side enum value
 */
inline boundary::side getBoundaryFromString(const std::string& str) {
    static const std::unordered_map<std::string, boundary::side> boundaryMap = {
        {"north", boundary::north},
        {"south", boundary::south},
        {"east", boundary::east},
        {"west", boundary::west},
        {"front", boundary::front},
        {"back", boundary::back},
        {"left", boundary::west},
        {"right", boundary::east},
        {"top", boundary::north},
        {"bottom", boundary::south}
    };

    auto it = boundaryMap.find(str);
    if (it != boundaryMap.end()) {
        return it->second;
    }

    gsWarn << "Unknown boundary: " << str << ", defaulting to south\n";
    return boundary::south;
}

/**
 * @brief Convert boundary::side enum to string
 * @param side boundary::side enum value
 * @return String representation of boundary
 */
inline std::string getBoundaryString(boundary::side side) {
    switch (side) {
        case boundary::north: return "north";
        case boundary::south: return "south";
        case boundary::east: return "east";
        case boundary::west: return "west";
        case boundary::front: return "front";
        case boundary::back: return "back";
        default: return "unknown";
    }
}

/**
 * @brief JSON serialization for boundary::side
 * @param j JSON object
 * @param side boundary::side to be written
 */
inline void to_json(json& j, const boundary::side& side) {
    j = getBoundaryString(side);
}

/**
 * @brief JSON deserialization for boundary::side
 * @param j JSON object
 * @param side boundary::side to be read
 */
inline void from_json(const json& j, boundary::side& side) {
    side = getBoundaryFromString(j.get<std::string>());
}

/**
 * @brief Puts a gsMultiPatch into JSON
 *
 * This function converts a gsMultiPatch object into a JSON structure with the following components:
 * - "patches": An array containing all patches of the multipatch
 * - "boundaries": An array of boundary objects, each containing patch index and side information
 * - "interfaces": An array of interface objects, each containing patch indices, sides, direction mappings,
 *   and orientation information for the connection between two patches
 *
 * @tparam T The numeric type used in the gsMultiPatch (e.g., real_t, double)
 * @param[out] j The JSON object to store the serialized multipatch data
 * @param[in] mp The gsMultiPatch object to serialize
 */
template<class T>
void to_json(json &j, const gsMultiPatch<T> & mp)
{
    j["patches"] = json::array();
    for (size_t p = 0; p!= mp.nPatches(); ++p)
        j["patches"].push_back(mp.patch(p));

    j["boundaries"] = json::array();
    for (typename gsMultiPatch<T>::const_biterator it = mp.bBegin(); it != mp.bEnd(); ++it)
    {
        json b;
        b["patch"] = it->patch();
        b["side"] = it->side();
        j["boundaries"].push_back(b);
    }
    j["interfaces"] = json::array();
    for (typename gsMultiPatch<T>::const_iiterator it = mp.iBegin(); it != mp.iEnd(); ++it)
    {
        json b;
        b["patch1"] = it->first().patch;
        b["patch2"] = it->second().patch;
        b["side1"] = it->first().side;
        b["side2"] = it->second().side;
        b["directionMap"] = it->dirMap().transpose();
        b["directionOrientation"] = it->dirOrientation().transpose();
        j["interfaces"].push_back(b);
    }
}

/**
 * @brief Reads a gsMultiPatch from JSON
 * @param j JSON object
 * @param mp gsMultiPatch to be read
 *
 * Supports three input formats:
 * 1. "patches" array - basis types:
 *    - "BSplineBasis" == gsBSpline<T>
 *    - "TensorBSplineBasis2" == gsTensorBSpline<2,T>
 *    - "TensorBSplineBasis3" == gsTensorBSpline<3,T>
 *    - "TensorBSplineBasis4" == gsTensorBSpline<4,T>
 * 2. "path" - XML file path, can read a gsMultiPatch from an XML file
 * 3. "create" - supports nurbs creator with (BSplineSquare, BSplineCube, BSplineRectangle)
 */

template<class T>
void from_json(const json &j, gsMultiPatch<T> & mp)
{
    mp.clear();
    if (j.contains("patches"))
    {
        for (const auto & patch : j["patches"])
            mp.addPatch(get_Geometry<T>(patch));
    }
    else if (j.contains("path"))
    {
        std::string path = j["path"].get<std::string>();
        GISMO_ASSERT(gsFileManager::fileExists(path),"File does not exist: " + path);
        GISMO_ASSERT(gsFileManager::getExtension(path) == "xml","File extension is not .xml: " + path);
        gsReadFile<T>(path, mp);
    }
    else if (j.contains("create"))
    {
        GISMO_ASSERT(j["create"].contains("type"),"A creator string must contain a type");
        if (j["create"]["type"] == "BSplineStar")
        {
            T N  = (j["create"].contains("N") ? j["create"]["N"].get<T>() : 1.0);
            T R0 = (j["create"].contains("R0") ? j["create"]["R0"].get<T>() : 0.0);
            T R1 = (j["create"].contains("R1") ? j["create"]["R1"].get<T>() : 0.0);
            mp = gsNurbsCreator<>::BSplineStar(N, R0, R1);
            return;
        }
    }
    else
    {
        GISMO_ERROR("To read a gsMultiPatch from JSON, it must either contain a 'patches' array, a 'path' to an XML file, or a 'type' with parameters for a specific patch type.");
    }
}
/**
 * @brief Writes a gsConstantFunction to JSON
 * @param j JSON object
 * @param bc function to be written
 */
template<class T>
void to_json(json &j, const gsConstantFunction<T> & fun)
{
    j["type"] = "ConstantFunction";
    j["dim"]  = fun.domainDim();

    const short_t tdim = fun.targetDim();

    if (tdim==1)
        j["value"] = fun.value(0);
    else
    {
        j["value"] = json::array();
        for (short_t i=0; i<tdim; ++i)
            j["value"].push_back(fun.value(i));
    }
}

/**
 * @brief Reads a gsConstantFunction from JSON
 * @param j JSON object
 * @param bc gsConstantFunction to be read
 */
template<class T>
void from_json(const json &j, gsConstantFunction<T> & fun)
{
    GISMO_ASSERT(j["type"]=="ConstantFunction","Type of function is not ConstantFunction");
    GISMO_ASSERT(j["dim"].is_number_integer(),"dim is not an integer");

    short_t dim = j["dim"].get<short_t>();
    gsVector<T> value;

    if (j["value"].is_array())
    {
        value.resize(j["value"].size());
        for (size_t i=0; i<j["value"].size(); ++i)
            value[i] = j["value"][i].get<T>();
    }
    else
    {
        value.resize(1);
        value[0] = j["value"].get<T>();
    }

    fun = gsConstantFunction<T>(value, dim);
}

/**
 * @brief Writes a gsFunctionExpr to JSON
 * @param j JSON object
 * @param bc function to be written
 */
template<class T>
void to_json(json &j, const gsFunctionExpr<T> & fun)
{
    j["type"] = "FunctionExpr";
    j["dim"]  = fun.domainDim();

    const short_t tdim = fun.targetDim();

    if (tdim==1)
        j["value"] = fun.expression();
    else
    {
        // Make array of strings
        j["value"] = json::array();
        for (short_t c=0; c<tdim; ++c)
            j["value"].push_back(fun.expression(c));
    }

}


/**
 * @brief Reads a gsFunctionExpr from JSON
 * @param j JSON object
 * @param bc gsFunctionExpr to be read
 */
template<class T>
void from_json(const json &j, gsFunctionExpr<T> & fun)
{
    GISMO_ASSERT(j["type"]=="FunctionExpr","Type of function is not FunctionExpr");
    GISMO_ASSERT(j["dim"].is_number_integer(),"dim is not an integer");

    short_t dim = j["dim"].get<short_t>();
    std::vector<std::string> expressions;

    if (j["value"].is_array())
    {
        expressions.resize(j["value"].size());
        for (index_t i=0; i<j["value"].size(); ++i)
            expressions[i] = j["value"][i].get<std::string>();
    }
    else
        expressions.push_back(j["value"].get<std::string>());

    fun = gsFunctionExpr<T>(expressions, dim);
}

/**
 * @brief Writes a gsFunction to JSON
 * @param j JSON object
 * @param bc function to be written
 */
template<class T>
void to_json(json &j, const gsFunctionSet<T> & fun)
{
    const gsFunctionSet<T> * ptr = & fun;

    if      (const gsConstantFunction<T> * f =
             dynamic_cast<const gsConstantFunction<T> *>(ptr))
        to_json(j, *f);

    else if (const gsFunctionExpr<T> * f =
             dynamic_cast<const gsFunctionExpr<T> *>(ptr))
        to_json(j, *f);
    else
        GISMO_ERROR("No known function type");
}

/**
 * @brief Reads a gsFunction from JSON
 * @param j JSON object
 * @param bc gsFunction to be read
 */
template<class T>
void from_json(const json &j, gsFunctionSet<T> & fun)
{
    GISMO_ASSERT(j["type"].is_string(),"Type of function is not a string");

    if      (j["type"] == "ConstantFunction")
        from_json(j, static_cast<gsConstantFunction<T> &>(fun));
    else if (j["type"] == "FunctionExpr")
        from_json(j, static_cast<gsFunctionExpr<T> &>(fun));
    else
        GISMO_ERROR("No known function type");
}

/**
 * @brief Writes gsBoundaryConditions to JSON
 * @param j JSON object
 * @param bc boundary conditions to be written
 */
template<class T>
void to_json(json &j, const gsBoundaryConditions<T> & bc)
{
    // Add a multi-patch index place holder
    j["multipatch"] = 0;

    // inventory of functions
    typedef typename gsBoundaryConditions<T>::const_bciterator bctype_it;

    std::vector<typename gsFunctionSet<T>::Ptr> fun;
    typedef typename std::vector<const boundary_condition<T>*> bctype_vec;
    typedef typename std::map<int, bctype_vec> bctype_map;
    std::map<std::string, bctype_map> fi;

    for (bctype_it it = bc.beginAll(); it != bc.endAll(); ++it)
    {
        std::string label = it->first;
        bctype_map map;
        for (typename gsBoundaryConditions<T>::const_iterator bc = it->second.begin(); bc != it->second.end(); ++bc)
        {
            typename gsFunctionSet<T>::Ptr ptr = bc->function();
            bool contains = std::find(fun.begin(), fun.end(), ptr)
                    != fun.end();
            if (!contains)
            {
                fun.push_back(ptr);
            }

            int index = std::find(fun.begin(), fun.end(), ptr)
                    - fun.begin();
//                fun.push_sorted_unique(ptr);
//                int index = fun.getIndex(ptr);
            std::vector<const boundary_condition<T>*> vec = map[index];
            const boundary_condition<T>* b = &(*bc);
            vec.push_back(b);
            map[index] = vec;
        }
        std::pair<std::string, bctype_map> pair(label, map);
        fi.insert(pair);
    }

    // Create a function array
    j["functions"] = json::array();
    // Add the functions
    int count = 0;
    typedef typename std::vector<typename gsFunctionSet<T>::Ptr>::const_iterator fun_it;
    for (fun_it fit = fun.begin(); fit != fun.end(); ++fit)
    {
        GISMO_ASSERT((dynamic_cast<gsFunctionSet<T> *>(fit->get())),"Function is not of type gsFunction<T>");
        to_json(j["functions"][count], *(fit->get()));
        ++count;
    }

    // for all bcs, append bc, cv
    typedef typename std::map<std::string, bctype_map>::const_iterator bctype_map_it;
    typedef typename std::map<int, bctype_vec>::const_iterator bctype_iv_it;
    typedef typename bctype_vec::const_iterator bctype_vec_it;

    // Create a BCs array
    j["sides"] = json::array();
    // Add the boundary conditions
    count = 0;
    for (bctype_map_it it = fi.begin(); it != fi.end(); ++it)
    {
        std::string label = it->first;
        bctype_map map = it->second;

        for (bctype_iv_it bcV = map.begin(); bcV != map.end(); ++bcV)
        {
            int index = bcV->first;
            bctype_vec vec = bcV->second;

            j["sides"][count]["type"] = label;
            j["sides"][count]["function"] = index;
            bool first = true;
            std::ostringstream oss;
            for (bctype_vec_it bc = vec.begin(); bc != vec.end(); ++bc)
            {
                const boundary_condition<T> b = (**bc);
                if (first)
                {
                    j["sides"][count]["unknown"] = b.m_unknown;
                    j["sides"][count]["component"] = b.unkComponent();
                    j["sides"][count]["sides"] = json::array();
                    first = false;
                }
                j["sides"][count]["sides"].push_back(std::make_pair(b.ps.patch, b.ps.m_index));
            }
            ++count;
        }
    }
    // Create a BCs array
    j["corners"] = json::array();
    // Add the boundary conditions
    count = 0;
    for (typename gsBoundaryConditions<T>::const_citerator ci = bc.cornerValues().begin(); ci != bc.cornerValues().end(); ci++)
    {
        corner_value<T> c = *ci;
        j["corners"][count]["unknown"] = c.unknown;
        j["corners"][count]["component"] = c.component;
        j["corners"][count]["patch"] = c.patch;
        j["corners"][count]["corner"] = c.corner.m_index;
        j["corners"][count]["value"] = c.value;
        count++;
    }
}

/**
 * @brief Reads gsBoundaryConditions from JSON
 * @param j JSON object
 * @param bc gsBoundaryConditions to be read
 */
template<class T>
void from_json(const json &j, gsBoundaryConditions<T> & bc)
{
    std::istringstream str;
    std::map<int, int> ids;

    // Check if any of the BCs is defined on a boundary set name
    // const int mp_index = j["multipatch"].get<int>();
    // gsXmlNode* toplevel = node->parent();
    std::vector< patchSide > allboundaries;
    for (json::const_iterator it = j["sides"].begin(); it != j["sides"].end(); ++it)
    {
        if (it->contains("name"))
            GISMO_ERROR("Boundary condition with name not supported in JSON format. "
                        "Please use the multipatch index instead.");
    }

    // Read function inventory using a map, in case indices are not
    // consecutive
    index_t index;
    std::map<int, typename gsFunctionExpr<T>::Ptr> function_map{};
    for (json::const_iterator it = j["functions"].begin();
         it != j["functions"].end(); ++it)
    {
        GISMO_ASSERT(it->contains("type"), "Function type not specified");
        std::string type = it->at("type").get<std::string>();
        index = std::distance(j["functions"].begin(), it);
        if (type == "ConstantFunction")
        {
            gsConstantFunction<T> fun;
            from_json(*it, fun);
            // Make a gsFunctionExpr
            std::vector<std::string> expressions;
            for (short_t i = 0; i < fun.targetDim(); ++i)
                expressions.push_back(util::to_string(fun.value(i)));

            gsFunctionExpr<T> funExpr(expressions, fun.domainDim());
            function_map[index] = memory::make_shared(funExpr.clone().release());
        }
        else if (type == "FunctionExpr")
        {
            gsFunctionExpr<T> fun;
            from_json(*it, fun);
            function_map[index] = memory::make_shared(fun.clone().release());
        }
        else
        {
            GISMO_ERROR("Unknown function type: " << type);
        }
    }

    // Read boundary conditions
    std::vector<patchSide> boundaries;
    for (json::const_iterator it = j["sides"].begin(); it != j["sides"].end(); ++it)
    {
        const int uIndex = it->at("unknown").get<int>();
        const int fIndex = it->at("function").get<int>();
        int cIndex = -1;
        if (it->contains("component"))  cIndex = it->at("component").get<int>();

        bool ispar = false;
        if (it->contains("parametric")) ispar = it->at("parametric").get<bool>();

        if (it->contains("name"))
            gsWarn<<"Boundary condition with name not supported in JSON format. "
                  "Please use the multipatch index instead.";

        // Get the boundary sides
        GISMO_ASSERT(it->contains("sides"), "No sides specified for boundary condition");
        for (const auto &side : it->at("sides"))
        {
            index_t patchIndex = side[0].get<index_t>();
            index_t sideIndex;

            // Handle both numeric and string boundary specifications
            if (side[1].is_number())
            {
                sideIndex = side[1].get<index_t>();
            }
            else if (side[1].is_string())
            {
                // Convert string boundary name to index
                boundary::side s = getBoundaryFromString(side[1].get<std::string>());
                sideIndex = static_cast<index_t>(s);
            }
            else
            {
                GISMO_ERROR("Boundary side must be either a number or a string");
            }

            boundaries.push_back(patchSide(patchIndex, sideIndex));
        }

        if (boundaries.size() == 0)
        {
            gsWarn << "Boundary condition without boundary to apply to. The"
                      " following bc will be unused\n"
                   << *it << std::endl;
        }

        GISMO_ASSERT(it->contains("type"), "No type provided for boundary condition");
        std::string bctype = it->at("type").get<std::string>();

        // Make BC type case-insensitive by capitalizing first letter
        if (!bctype.empty())
        {
            bctype[0] = std::toupper(bctype[0]);
            for (size_t i = 1; i < bctype.length(); ++i)
                bctype[i] = std::tolower(bctype[i]);
        }

        // Add the boundary conditions
        for (std::vector<patchSide>::const_iterator itb = boundaries.begin(); itb != boundaries.end(); ++itb)
            bc.add(itb->patch, itb->side(), bctype, function_map[fIndex], uIndex, cIndex, ispar);
    }

    // Read corner values
    for (json::const_iterator it = j["corners"].begin(); it != j["corners"].end(); ++it)
    {
        GISMO_ASSERT(it->contains("corner"), "No corner specified for corner value");
        GISMO_ASSERT(it->contains("patch"), "No patch specified for corner value");
        const int cornIndex = it->at("corner").get<int>();
        const int pIndex = it->at("patch").get<int>();

        int uIndex = 0;
        if (it->contains("unknown"))    uIndex = it->at("unknown").get<int>();

        int cIndex = -1;
        if (it->contains("component"))  cIndex = it->at("component").get<int>();

        T value = 0.0;
        if (it->contains("value"))      value = it->at("value").get<T>();

        bc.addCornerValue(cornIndex, value, pIndex, uIndex, cIndex);
    }
}

/**
 * @brief JSON class
 */
class gsJSON
{
public:
    typedef json::iterator iterator;
    typedef json::const_iterator const_iterator;


public:

    /**
     * @brief Default constructor
     * @param opt Options
     */
    gsJSON()
    {
        m_data = json::object();
    }

    /**
     * @brief Constructor from a \ref gsOptionList
     * @param opt Options
     */
    gsJSON(const gsOptionList & opt)
    :
    m_data(opt)
    {
    }

    /**
     * @brief Constructor from a file
     * @param filename File name
     */
    gsJSON(const std::string & filename)
    {
        // Open the file
        std::ifstream file(filename);
        // if (!file.is_open()) {
        //     throw std::runtime_error("Could not open file: " + filename);
        // }

        GISMO_ENSURE(file.is_open(), "Could not open file: " + filename);

        // Read the entire file content
        std::string content((std::istreambuf_iterator<char>(file)),
                                std::istreambuf_iterator<char>());
        file.close();

        // Parse the file content
        m_data = json::parse(content,
                            /* callback */ nullptr,
                            /* allow exceptions */ true,
                            /* ignore_comments */ true);
    }

    /// add an object
    json::reference operator[](const std::string & key)
    {
        return m_data[key];
    }

    /// add objects from an initializer list


    /// get an object
    json::const_reference operator[](const std::string & key) const
    {
        return m_data[key];
    }

    /// add an object
    template <class U>
    void add(const std::string & key, const U & value)
    {
        m_data[key] = value;
    }

    /// get an object
    template <class U>
    U get(const std::string & key) const
    {
        return m_data[key];
    }

    /// get an object
    template <class U>
    U get() const
    {
        return m_data.get<U>();
    }

    /// get an object
    template <class U>
    void get_to(U & obj) const
    {
        m_data.get_to(obj);
    }

    /// get the size
    size_t size() const
    {
        return m_data.size();
    }

    /// check if the data is empty
    bool empty() const
    {
        return m_data.empty();
    }

    /// clear the data
    void clear()
    {
        m_data.clear();
    }

    /// check if the data contains
    bool contains(const std::string & key) const
    {
        return m_data.contains(key);
    }

    /// find an item
    iterator find(const std::string & key)
    {
        return m_data.find(key);
    }

    /// count the number of items
    size_t count(const std::string & key)
    {
        return m_data.count(key);
    }

    /// erase an item
    void erase(const std::string & key)
    {
        m_data.erase(key);
    }

    std::ostream& print( std::ostream& os ) const
    {
        os<<m_data.dump(4);
        return os;
    }

    // ITERATORS
    /// begin iterator
    iterator begin()
    {
        return m_data.begin();
    }

    /// end iterator
    iterator end()
    {
        return m_data.end();
    }

    /// const begin iterator
    const_iterator begin() const
    {
        return m_data.begin();
    }

    /// const end iterator
    const_iterator end() const
    {
        return m_data.end();
    }

    void save(const std::string & filename) const
    {
        std::ofstream file(filename);
        file << m_data.dump(4);
        file.close();
    }

    /// Get the raw nlohmann::json object
    const json& getRaw() const
    {
        return m_data;
    }

    // template<class Object>
    // Object * getObject(const std::string & key)
    // {
    //     if (m_data[key].is_null())
    //         return nullptr;
    //     else
    //     {
    //         Object * obj = new Object();
    //         m_data[key].get_to(*obj);
    //         return obj;
    //     }
    // }

protected:
    json m_data;

};

std::ostream& operator<<(std::ostream& os, const gsJSON& data)
{
    return data.print(os);
}

} //namespace gismo


