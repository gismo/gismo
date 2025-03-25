/** @file gsComposedGeometry.hpp

    @brief Implementation of the gsComposedGeometry class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst
        S. Imperatore
*/

#pragma once

#include <gsIO/gsXml.h>
#include <gsIO/gsXmlGenericUtils.hpp>
#include <gsUtils/gsMesh/gsMesh.h>

namespace gismo
{

template <class T>
gsComposedGeometry<T>::gsComposedGeometry()
:
m_domainDim(0)
{
}

template <class T>
gsComposedGeometry<T>::gsComposedGeometry(const gsComposedBasis<T> & basis, const gsMatrix<T> & coefs)
:
Base(basis, coefs ),
m_domainDim(basis.composition().domainDim())
{
}

template <class T>
gsComposedGeometry<T>::gsComposedGeometry(const gsFunction<T> & composition, const gsGeometry<T> & geom)
:
Base(gsComposedBasis<T>(composition,geom.basis()), geom.coefs() ),
m_domainDim(geom.domainDim())
{
    GISMO_ASSERT(geom.domainDim()==composition.targetDim(),"Domain dimension of the geometry does not correspond with the target dimension of the composition!");
}

// template <class T>
// gsComposedGeometry<T>::gsComposedGeometry(const gsComposedGeometry<T> & other)
// :
// Base(other),
// m_composition(memory::make_shared(other.composition().clone().release())),
// m_geom(give(other.m_geom)),
// m_domainDim(other.domainDim())
// {
// }

// template <class T>
// gsComposedGeometry<T>::gsComposedGeometry( gsComposedGeometry&& other )
// :
// Base(give(other)),
// m_composition(other.m_composition),
// m_geom(give(other.m_geom)),
// m_domainDim(other.m_domainDim)
// { }


// template <class T>
// gsComposedGeometry<T> & gsComposedGeometry<T>::operator=(const gsComposedGeometry<T> & other)
// {
//     if (this != &other)
//     {
//         Base::operator=(other);
//         m_composition = memory::make_shared(other.composition().clone().release());
//         m_geom = give(other.m_geom);
//         m_domainDim = other.domainDim();
//     }
//     return *this;
// }

// template <class T>
// gsComposedGeometry<T> & gsComposedGeometry<T>::operator=(gsComposedGeometry<T> && other)
// {
//     if (this != &other)
//     {
//         m_composition = other.m_composition;
//         m_geom = give(other.m_geom);
//         m_domainDim = other.m_domainDim;
//         Base::operator=(other);
//     }
//     return *this;
// }

template <class T>
const typename gsComposedGeometry<T>::CompositionT & gsComposedGeometry<T>::composition() const
{
    return static_cast<gsComposedBasis<T> *>(m_basis)->composition();
}

template <class T>
typename gsComposedGeometry<T>::CompositionT & gsComposedGeometry<T>::composition()
{
    return static_cast<gsComposedBasis<T> *>(m_basis)->composition();
}

template <class T>
short_t gsComposedGeometry<T>::domainDim() const
{
    return m_domainDim;
}

template <class T>
void gsComposedGeometry<T>::compute(const gsMatrix<T> & in, gsFuncData<T> & out) const
{
    unsigned flags = NEED_ACTIVE;
    if (out.flags & NEED_VALUE)  flags |= NEED_VALUE;
    if (out.flags & NEED_DERIV)  flags |= NEED_DERIV;
    if (out.flags & NEED_DERIV2) flags |= NEED_DERIV2;

    gsFuncData<T> tmp(flags);
    Base::compute(in,tmp);

    out.dim = tmp.dim;
    out.values = tmp.values;
}

/* @hverhelst: I am not sure if this is needed. If so, we need to enable gsFunction::nControls and gsFunction::control_deriv_into

template <class T>
void gsComposedGeometry<T>::control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const
{
    const CompositionT & composition = static_cast<gsComposedBasis<T> *>(m_basis)->composition();
    // The number of rows is the target dimension times the number of controls
    // The number of cols is the number of points
    result.resize(targetDim()*composition.nControls(),points.cols());

    // Pre-compute the coordinates of the composition, the derivatives of G and the derivatives of the composition
    gsMatrix<T> c, dc, dG;
    composition.eval_into(points,c);
    composition.control_deriv_into(points,dc);   // This is dc/dpi (pi is a control of c)
    m_geom->deriv_into(c,dG);                       // This is dG/dc evaluated on c

    // Store some sizes
    index_t nControls = composition.nControls();
    index_t dd = m_geom->domainDim();
    index_t td = m_geom->targetDim();

    // Loop over the points
    for (index_t k=0; k!=points.cols(); k++)
    {
        // We need to compute dG/dpi = dG/dc * dc/dpi
        gsAsMatrix<T> DG = result.reshapeCol(k,nControls,td);
        DG = dc.reshapeCol(k,nControls,dd) * dG.reshapeCol(k,dd,td);
    }
}
*/


template <class T>
void gsComposedGeometry<T>::evaluateMesh(gsMesh<T>& mesh) const
{
    const int pDim = this->parDim();
    const int gDim = this->geoDim();

    gsMatrix<T> tmp;

    // For all vertices of the mesh, push forward the value by the
    // geometry mapping
    if (1==gDim && 3>pDim) // Plot a graph
        for (size_t i = 0; i!= mesh.numVertices(); ++i)
        {
            // m_composition->invertPoints(tmp,tmp);
            this->eval_into( mesh.vertex(i).topRows(pDim), tmp );
            mesh.vertex(i).middleRows(pDim, gDim) = tmp;
        }
    else // Plot mesh on a mapping
        for (size_t i = 0; i!= mesh.numVertices(); ++i)
        {
            this->eval_into( mesh.vertex(i).topRows(pDim), tmp );
            // m_composition->invertPoints(tmp,tmp);
            const index_t gd = math::min(3,gDim);
            mesh.vertex(i).topRows(gd) = tmp.topRows(gd);
        }
}

namespace internal
{

/// @brief Get a Tensor BSpline from XML data
///
/// \ingroup Nurbs
template<class T>
class gsXml< gsComposedGeometry<T> >
{
private:
    gsXml() { }
public:
    GSXML_COMMON_FUNCTIONS(gsComposedGeometry<T>);
    GSXML_GET_INTO(gsComposedGeometry<T>);
    static std::string tag ()  { return "Geometry"; }
    static std::string type () { return "ComposedGeometry"; }

    static gsComposedGeometry<T> * get (gsXmlNode * node)
    {
        return getGeometryFromXml< gsComposedGeometry<T> >( node );
    }

    static gsXmlNode * put (const gsComposedGeometry<T> & obj,
                            gsXmlTree & data)
    {
        return putGeometryToXml(obj,data);
    }
};

}// namespace internal

};// namespace gismo
