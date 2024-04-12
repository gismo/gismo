/** @file gsComposedGeometry.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
        H.M. Verhelst   (2019-..., TU Delft)
        A. Mantzaflaris (2019-..., Inria)
*/

#pragma once

#include <gsCore/gsComposedBasis.h>

namespace gismo
{

template<class T>
class gsComposedGeometry : public gsGeometry<T>
{

    using Base = gsGeometry<T>;

    typedef gsComposedBasis<T> Basis;
    typedef typename gsComposedBasis<T>::CompositionT CompositionT;

    typedef memory::shared_ptr< gsComposedGeometry > Ptr;
    typedef memory::unique_ptr< gsComposedGeometry > uPtr;

    GISMO_CLONE_FUNCTION(gsComposedGeometry)

public:
    gsComposedGeometry()
    :
    m_geom(nullptr)
    {

    }

    /**
     * @brief      XXXX
     *
     * @param[in]  basis  The basis
     * @param[in]  coefs  The coefs
     */
    gsComposedGeometry(const gsComposedBasis<T> & basis, const gsMatrix<T> & coefs)
    :
    Base(basis, give(coefs) ),
    m_composition(basis.composition()),
    m_geom(basis.basis()->makeGeometry(coefs).release()),
    m_domainDim(basis.domainDim())
    { }

    /**
     * @brief      XXXX
     *
     * @param[in]  composition  The composition
     * @param[in]  geom         The geometry
     */
    gsComposedGeometry(const gsFunction<T> & composition, const gsGeometry<T> & geom)
    :
    Base(gsComposedBasis<T>(composition,geom.basis()), give(geom.coefs()) ),
    m_composition(&composition),
    m_geom(geom.clone().release()),
    m_domainDim(geom.domainDim())
    {
        GISMO_ASSERT(geom.domainDim()==composition.targetDim(),"Domain dimension of the geometry does not correspond with the target dimension of the composition!");
    }

    ~gsComposedGeometry()
    {
        if (m_geom!=nullptr)
            delete m_geom;
    }

    short_t domainDim() const override { return m_domainDim; }

    /**
     * @brief      Gives the control point derivatives of \a *this. See gsFunction for more details
     *
     * @param[in]  points  The points in the parameter domain (of the composition)
     * @param      result  The control point derivatives
     */
    void control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const override
    {
        // The number of rows is the target dimension times the number of controls
        // The number of cols is the number of points
        result.resize(targetDim()*m_composition->nControls(),points.cols());

        // Pre-compute the coordinates of the composition, the derivatives of G and the derivatives of the composition
        gsMatrix<T> c, dc, dG;
        m_composition->eval_into(points,c);
        m_composition->control_deriv_into(points,dc);   // This is dc/dpi (pi is a control of c)
        m_geom->deriv_into(c,dG);                       // This is dG/dc evaluated on c

        // Store some sizes
        index_t nControls = m_composition->nControls();
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

    using Base::targetDim;

    GISMO_BASIS_ACCESSORS;

protected:
    // Map from parametric domain to geometry
    const CompositionT * m_composition;
    using Base::m_basis;

    // Map from composition to geometry
    gsGeometry<T> * m_geom;

    short_t m_domainDim;
};
}
