/** @file gsConstantFunction.hpp

    @brief Provides implementation of FunctionExpr class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
               H.M. Verhelst
*/


#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsUtils/gsCombinatorics.h>
#include <gsCore/gsConstantBasis.h>
#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>


namespace gismo
{

template <class T>
gsConstantFunction<T>::gsConstantFunction(const gsVector<T> & val, short_t domainDim)
:
m_domainDim(domainDim)
{
    m_coefs = val.transpose();
}

template <class T>
gsConstantFunction<T>::gsConstantFunction(T x, T y, short_t domainDim)
:
m_domainDim(domainDim)
{
    m_coefs.resize(1,2);
    m_coefs(0,0) = x;
    m_coefs(0,1) = y;
}

template <class T>
gsConstantFunction<T>::gsConstantFunction(T x, T y, T z, short_t domainDim)
:
m_domainDim(domainDim)
{
    m_coefs.resize(1,3);
    m_coefs(0,0) = x;
    m_coefs(0,1) = y;
    m_coefs(0,2) = z;
}

template <class T>
gsConstantFunction<T>::gsConstantFunction(T x, T y, T z, T w, short_t domainDim)
:
m_domainDim(domainDim)
{
    m_coefs.resize(1,4);
    m_coefs(0,0) = x;
    m_coefs(0,1) = y;
    m_coefs(0,2) = z;
    m_coefs(0,3) = w;
}

template <class T>
gsConstantFunction<T>::gsConstantFunction(const gsConstantBasis<T> & cb, const gsMatrix<T> & coef)
:
m_domainDim(1)
{
    m_coefs = cb.value()*coef;
}

template <class T>
gsConstantFunction<T>::gsConstantFunction(const gsConstantFunction<T> & o)
:
gsGeometry<T>(),
m_domainDim(o.m_domainDim)
{
    m_coefs = o.m_coefs;
}

template <class T>
gsConstantFunction<T>::gsConstantFunction(gsConstantFunction<T> && o)
:
m_domainDim(o.m_domainDim)
{
    m_coefs.swap(o.m_coefs);
}

template <class T>
gsConstantFunction<T> & gsConstantFunction<T>::operator=(const gsConstantFunction<T> & o)
{
    if (this != &o)
    {
        m_domainDim = o.m_domainDim;
        m_coefs = o.m_coefs;
    }
    return *this;
}

template <class T>
gsConstantFunction<T> & gsConstantFunction<T>::operator=(gsConstantFunction<T> && o)
{
    if (this != &o)
    {
        m_domainDim = o.m_domainDim;
        m_coefs.swap(o.m_coefs);
    }
    return *this;
}

template <class T>
void gsConstantFunction<T>::eval_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                  << ", expected "<< m_domainDim);
    result = m_coefs.transpose().rowwise().replicate( u.cols() );
}

template <class T>
void gsConstantFunction<T>::deriv_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                  << ", expected "<< m_domainDim);
    result = gsMatrix<T>::Zero(this->targetDim()*m_domainDim, u.cols());
}

template <class T>
void gsConstantFunction<T>::deriv2_into(const gsMatrix<T> & u, gsMatrix<T> & result) const
{
    GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                  << ", expected "<< m_domainDim);
    result = gsMatrix<T>::Zero(this->targetDim()*(m_domainDim*(m_domainDim+1))/2, u.cols());
}

template <class T>
void gsConstantFunction<T>::evalAllDers_into(const gsMatrix<T> & u, int n,
                                             std::vector<gsMatrix<T> > & result,
                                             bool sameElement) const
{
    GISMO_UNUSED(sameElement);
    GISMO_ASSERT(u.rows() == m_domainDim, "Wrong domain dimension "<< u.rows()
                 << ", expected "<< m_domainDim);

    result.resize(n+1,gsMatrix<T>());
    eval_into(u,result.front());
    for (int i = 1; i<=n; ++i)
    {
        result[i] = gsMatrix<T>::Zero(this->targetDim()*m_domainDim, u.cols());
    }
}

namespace internal
{

/// @brief Get a gsConstantFunction from XML data
template<class T>
class gsXml< gsConstantFunction<T> >
{
private:
    gsXml() { }
    typedef gsConstantFunction<T> Object;
public:
    GSXML_COMMON_FUNCTIONS(Object);
    GSXML_GET_INTO(Object);
    static std::string tag () { return "Function"; }
    static std::string type () { return "ConstantFunction"; }

    static Object * get (gsXmlNode * node)
    {
        GISMO_ASSERT( ( !strcmp( node->name(),"Function") )
                    &&  ( !strcmp(node->first_attribute("type")->value(),
                                internal::gsXml<Object>::type().c_str() ) ),
                    "Reading gsConstantFunction XML: No Function found" );

        GISMO_ASSERT( node->first_attribute("dim"), "Reading gsConstantFunction XML: No dim found" ) ;
        const int d = atoi( node->first_attribute("dim")->value() );

        gsVector<T> val;
        gsXmlNode * child = node->first_node("c");
        if (child != NULL )
        {
            val.resize(countByTag("c",node));
            index_t idx = 0;
            std::istringstream str;
            str.str(child->value());
            for (; child; child = child->next_sibling(), idx++ )
                gsGetValue(str, val.at(idx));
        }
        else
        {
            val.resize(1);
            std::istringstream str;
            str.str(node->value());
            gsGetValue(str, val.at(0));
        }

        return new Object(val, d);
    }

    static gsXmlNode * put (const Object & obj,
                            gsXmlTree & data )
    {
        // Add a new node
        gsXmlNode* node = internal::makeNode("Function" , data);
        node->append_attribute( makeAttribute("type",
                                            internal::gsXml< Object >::type().c_str(), data) );
        node->append_attribute(makeAttribute("dim", obj.domainDim(), data));

        const short_t tdim = obj.targetDim();

        if ( tdim == 1)
        {
            node->value( makeValue(std::to_string(obj.value(0)), data) );
        }
        else
        {
            gsXmlNode * cnode;
            for (short_t c = 0; c!=tdim; ++c)
            {
                cnode = makeNode("c", std::to_string(obj.value(c)), data);
                node->append_node(cnode);
            }
        }

        return node;
    }
};

} // namespace internal


} // namespace gismo