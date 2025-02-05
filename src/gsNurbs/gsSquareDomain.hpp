/** @file gsSquareDomain.hpp

    @brief Implementation of the gsSquareDomain class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDofMapper.h>
#include <gsIO/gsOptionList.h>
#include <gsCore/gsBoxTopology.h>
// #include <gsIO/gsXml.h>
// #include <gsIO/gsXmlGenericUtils.hpp>

namespace gismo
{

template <class T>
gsSquareDomain<T>::gsSquareDomain(const gsGeometry<T> & domain)
{
    m_domain = memory::make_unique(domain.clone().release());
    this->_initMapper(*m_domain,m_mapper);
    this->_initIndices(*m_domain,m_mapper,m_indices);
}

template <class T>
gsSquareDomain<T>::gsSquareDomain(const gsBasis<T> & basis)
{
    gsMatrix<T> coefs = basis.anchors();
    coefs.transposeInPlace();
    m_domain = basis.makeGeometry(give(coefs));
    this->_initMapper(*m_domain,m_mapper);
    this->_initIndices(*m_domain,m_mapper,m_indices);
}

// template <class T>
// gsSquareDomain<T>::gsSquareDomain(short_t DIM, index_t numElevation, index_t numRefine)
// {
//     // Create a linear square domain
//     gsKnotVector<T> KV({0,0,1,1},1);
//     std::vector<gsKnotVector<T>> KVs(DIM,KV);
//     if (DIM==2)
//         gsTensorBSplineBasis<DIM,T> basis(KVs);
//     gsMatrix<T> coefs = basis.anchors();
//     coefs.transposeInPlace();
//     m_domain = basis.makeGeometry(give(coefs));

//     m_domain->degreeElevate(numElevation);
//     index_t numKts = pow(2, numRefine) - 1;
//     m_domain->uniformRefine(numKts);

//     this->_initMapper(*m_domain,m_mapper);
//     this->_initIndices(*m_domain,m_mapper,m_indices);
// }

// Copy constructor
template <class T>
gsSquareDomain<T>::gsSquareDomain(const gsSquareDomain<T> & other)
{
    m_domain = memory::make_unique(other.domain().clone().release());
    m_mapper = other.mapper();
    m_indices = other.m_indices;
}

// Copy assignment
template <class T>
gsSquareDomain<T> & gsSquareDomain<T>::operator=(const gsSquareDomain<T> & other)
{
    if (this != &other)
    {
        m_domain = memory::make_unique(other.domain().clone().release());
        m_mapper = other.mapper();
        m_indices = other.m_indices;
    }
    return *this;
}

template <class T>
void gsSquareDomain<T>::addInterface(const boundaryInterface & bi)
{
    m_interfaces.push_back(bi);
}

template <class T>
gsOptionList & gsSquareDomain<T>::options()
{
    return m_options;
}

template <class T>
void gsSquareDomain<T>::applyOptions()
{
    this->_initMapper(*m_domain,m_mapper);
    this->_initIndices(*m_domain,m_mapper,m_indices);
}

template <class T>
const gsGeometry<T> & gsSquareDomain<T>::domain() const
{
    return *m_domain;
}

template <class T>
const gsDofMapper & gsSquareDomain<T>::mapper() const
{
    return m_mapper;
}

template <class T>
gsMatrix<T> gsSquareDomain<T>::support() const
{
    return m_domain->support();
}

template <class T>
short_t gsSquareDomain<T>::maxDegree() const
{
    return m_domain->basis().maxDegree();
}

template <class T>
short_t gsSquareDomain<T>::domainDim() const
{
    return m_domain->domainDim();
}

template <class T>
short_t gsSquareDomain<T>::targetDim() const
{
    return m_domain->targetDim();
}

template <class T>
void gsSquareDomain<T>::eval_into(const gsMatrix<T> & param, gsMatrix<T> & result) const
{
    m_domain->eval_into(param,result);
}

template <class T>
void gsSquareDomain<T>::deriv_into(const gsMatrix<T> & param, gsMatrix<T> & result) const
{
    m_domain->deriv_into(param,result);
}

template <class T>
void gsSquareDomain<T>::setControls(const gsVector<T> & controls)
{
    GISMO_ASSERT((size_t)controls.rows()==m_indices.size(),"Wrong size of controls vector");
    index_t ii;
    for (index_t d = 0; d!=m_domain->domainDim(); d++)
        for (index_t k = 0; k!=m_domain->coefs().rows(); k++)
        {
            ii = m_mapper.index(k,0,d);
            if (m_mapper.is_free_index(ii))
                m_domain->coefs()(k,d) = controls[ii];
        }
}

template <class T>
gsVector<T> gsSquareDomain<T>::getControls() const
{
    gsVector<T> controls(m_indices.size());
    index_t ii;
    for (index_t d = 0; d!=m_domain->domainDim(); d++)
        for (index_t k = 0; k!=m_domain->coefs().rows(); k++)
        {
            ii = m_mapper.index(k,0,d);
            if (m_mapper.is_free_index(ii))
                controls[ii] = m_domain->coefs()(k,d);
        }
    return controls;
}

// template <class T>
// const T & gsSquareDomain<T>::control(index_t i) const
// {
//     return m_domain->coefs()(m_indices[i].first,m_indices[i].second);
// }

// template <class T>
// T & gsSquareDomain<T>::control(index_t i)
// {
//     return m_domain->coefs()(m_indices[i].first,m_indices[i].second);
// }

template <class T>
size_t gsSquareDomain<T>::nControls() const
{
    return m_indices.size();
}

template <class T>
void gsSquareDomain<T>::control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const
{
    gsMatrix<T> tmp;

    result.resize(targetDim()*nControls(), points.cols());
    result.setZero();
    for (index_t p = 0; p!=points.cols(); p++)
    {
        gsAsMatrix<T> res = result.reshapeCol(p,nControls(),targetDim());
        for (index_t k = 0; k!=m_domain->coefs().rows(); k++)
            for (index_t d = 0; d!=m_domain->targetDim(); d++)
                if (m_mapper.is_free(k,0,d))
                {
                    m_domain->basis().evalSingle_into(k,points.col(p),tmp); // evaluate basis function k
                    res(m_mapper.index(k,0,d),d) = tmp(0,0); // tmp is a single value (1 point, 1 basis function)
                }
    }
}

template <class T>
void gsSquareDomain<T>::perturb(T factor)
{
    gsVector<T> controls = getControls();
    controls += factor * gsVector<T>::Random(m_indices.size());
    setControls(controls);
}

template <class T>
void gsSquareDomain<T>::_initMapper(const gsGeometry<T> & domain, gsDofMapper & mapper) const
{
    mapper = gsDofMapper(domain.basis(),domain.domainDim());

    gsBoxTopology topology(domain.domainDim(),1);
    topology.addAutoBoundaries();
    // gsMatrix<index_t> boundary = domain.basis().allBoundary();
    for (typename gsBoxTopology::biterator it = topology.bBegin(); it != topology.bEnd(); ++it)
    {
        gsMatrix<index_t> boundary = domain.basis().boundary(*it);
        if (m_options.askSwitch("Slide",false))
            mapper.markBoundary(0,boundary,it->direction());
        else
            for (index_t d = 0; d!=domain.domainDim(); d++)
                mapper.markBoundary(0,boundary,d);

    }

    gsMatrix<index_t> int10, int20, int11, int21;
    for (typename std::vector<boundaryInterface>::const_iterator ifc = m_interfaces.cbegin();
                                                                 ifc!= m_interfaces.cend(); ++ifc)
    {
        GISMO_ENSURE(ifc->first().patch == 0 && ifc->second().patch == 0,"Interfaces are only supported for single patch domains");
        GISMO_ENSURE(ifc->first().direction() == ifc->second().direction(),"Interfaces are only supported for same direction");

        int10 = domain.basis().boundaryOffset(ifc->first().side(), 0);
        int20 = domain.basis().boundaryOffset(ifc->second().side(), 0);
        int11 = domain.basis().boundaryOffset(ifc->first().side(), 1);
        int21 = domain.basis().boundaryOffset(ifc->second().side(), 1);
        GISMO_ASSERT(int10.rows() == int20.rows() && int10.rows() == int11.rows() && int10.rows() == int21.rows(),"Boundary offsets have different sizes");
        index_t sz = int10.rows();
        for (index_t d = 0; d!=domain.domainDim(); d++)
        {
            if (d!=ifc->first().direction())
            {
                for ( index_t k=0; k<sz; ++k)
                {
                    mapper.matchDof( ifc->first().patch, int10(k,0), ifc->second().patch, int20(k,0), d );
                    mapper.matchDof( ifc->first().patch, int10(k,0), ifc->second().patch, int11(k,0), d );
                    mapper.matchDof( ifc->first().patch, int20(k,0), ifc->second().patch, int21(k,0), d );
                }
            }
        }
    }

    mapper.finalize();
}

template <class T>
void gsSquareDomain<T>::_initIndices(const gsGeometry<T> & domain, const gsDofMapper & mapper, std::vector<std::pair<index_t,index_t>> & indices) const
{
    indices.resize(mapper.freeSize());
    // std::vector<index_t> i(mapper.freeSize());
    // std::vector<index_t> j(mapper.freeSize());
    index_t ii;
    for (index_t k = 0; k!=domain.coefs().rows(); k++)
        for (index_t d = 0; d!=domain.domainDim(); d++)
        {
            ii = mapper.index(k,0,d);
            if (mapper.is_free_index(ii))
                indices[mapper.index(k,0,d)] = std::make_pair(k,d);
        }
}


namespace internal
{

/// @brief Get a gsSquareDomain from XML data
template<class T>
class gsXml< gsSquareDomain<T> >
{
private:
    gsXml() { }
public:
    // GSXML_COMMON_FUNCTIONS(gsSquareDomain<TMPLA2(T)>);
    static std::string tag () { return gsXml<gsGeometry<T>>::tag(); }
    static std::string type () { return gsXml<gsGeometry<T>>::type(); }

    // static gsGeometry<T> * get (gsXmlNode * node)
    // {
    //     return getGeometryFromXml< gsGeometry<T> >( node );
    // }

    static gsXmlNode * put (const gsSquareDomain<T> & obj,
                            gsXmlTree & data )
    {
        return putGeometryToXml(obj.domain(),data);
    }
};

} // internal

} // namespace gismo

