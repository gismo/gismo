/** @file gsJSON.h

    @brief Header file for using Spectra extension

    https://spectralib.org/doc

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <JSON/include/nlohmann/json.hpp>
using json = nlohmann::json;

namespace gismo
{

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
        mat.resize(j["rows"], j["cols"]);
        for (index_t I = 0; I < mat.rows(); I++)
            for (index_t J = 0; J < mat.cols(); J++)
                mat(I,J) = J["data"][I*mat.cols() + J];
    }
}

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

template<class T>
void to_json(json &j, const gsKnotVector<T> & kv)
{
    j["degree"] = kv.degree();
    std::vector<T> data(kv.data(), kv.data() + kv.size());
    j["knots"] = data;
}

template<class T>
void from_json(const json &j, gsKnotVector<T> & kv)
{
    kv = gsKnotVector<T>(j["degree"], j["knots"]);
}

template<class T>
void to_json(json &j, const gsBSplineBasis<T> & basis)
{
    j["type"] = "BSplineBasis";
    j["degree"] = basis.degree();
    j["knots"] = basis.knots();
}

template<class T>
void from_json(const json &j, gsBSplineBasis<T> & basis)
{
    gsKnotVector<T> kv;
    j.get_to(kv);
    basis = gsBSplineBasis<T>(kv);
}

template<short_t d, class T>
void to_json(json &j, const gsTensorBasis<d,T> & basis)
{
    j["type"] = "TensorBasis";
    for (unsigned D=0; D!=d; D++)
        j[util::to_string(D)] = basis.component(D);
}

template<class T>
void to_json(json &j, const gsBasis<T> & basis)
{
    if      ( const gsBSplineBasis<T> * b = dynamic_cast<const gsBSplineBasis<T> *>( &basis ) )
        j = json(*b);
    else if ( const gsTensorBSplineBasis<2,T> * b = dynamic_cast<const gsTensorBSplineBasis<2,T> *>( &basis ) )
        j = json(*b);
    else if ( const gsTensorBSplineBasis<3,T> * b = dynamic_cast<const gsTensorBSplineBasis<3,T> *>( &basis ) )
        j = json(*b);
    else if ( const gsTensorBSplineBasis<4,T> * b = dynamic_cast<const gsTensorBSplineBasis<4,T> *>( &basis ) )
        j = json(*b);
    else
        GISMO_ERROR("No known basis type");
}

template<class T>
void from_json(const json &j, gsBasis<T> & basis)
{
    GISMO_ERROR("HOW TO DO THIS?");
}

template<class T>
void to_json(json &j, const gsGeometry<T> & geo)
{
    j["basis"] = geo.basis();
    j["coefs"] = geo.coefs();
}

template<class T>
void from_json(const json &j, gsGeometry<T> & geo)
{
    GISMO_ERROR("HOW TO DO THIS?");
}

class gsJSON
{
    public:
        typedef json::iterator iterator;
        typedef json::const_iterator const_iterator;


    public:

        gsJSON()
        {
            m_data = json::object();
        }

        gsJSON(const gsOptionList & opt)
        :
        m_data(opt)
        {
        }

        gsJSON(const std::string & filename)
        {
            m_data = json::parse(filename,
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

    protected:
        json m_data;

};

std::ostream& operator<<(std::ostream& os, const gsJSON& data)
{
    return data.print(os);
}

} //namespace gismo


