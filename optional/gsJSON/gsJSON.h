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
    j["degree"] = basis.degree();
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
    GISMO_ASSERT(j["type"]=="TensorBSplineBasis"+util::to_string(d),"Type of basis is not TensorBSplineBasis"+util::to_string(d));

    std::vector<gsKnotVector<T> > KVs(d);
    gsBSplineBasis<T> componentBasis;
    for (unsigned D=0; D!=d; D++)
    {
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
 * @brief Reads a gsTensorBSpline from JSON
 * @param j JSON object
 * @param geo gsTensorBSpline to be read
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


