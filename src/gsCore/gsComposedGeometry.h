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

#include <gsCore/gsGeometry.h>
#include <gsCore/gsComposedBasis.h>
#include <gsCore/gsFuncData.h>

namespace gismo
{

template<class T>
class gsComposedGeometry : public gsGeometry<T>
{

    using Base = gsGeometry<T>;

    GISMO_OVERRIDE_CLONE_FUNCTION(gsComposedGeometry)

public:

    typedef memory::shared_ptr< gsComposedGeometry > Ptr;
    typedef memory::unique_ptr< gsComposedGeometry > uPtr;

    typedef gsComposedBasis<T> Basis;
    typedef gsComposedBasis<T> BasisT;
    typedef typename BasisT::CompositionT CompositionT;

public:


    gsComposedGeometry();

    /**
     * @brief      XXXX
     *
     * @param[in]  basis  The basis
     * @param[in]  coefs  The coefs
     */
    gsComposedGeometry( const gsComposedBasis<T> & basis,
                        const gsMatrix<T> & coefs);

    /**
     * @brief      XXXX
     *
     * @param[in]  composition  The composition
     * @param[in]  geom         The geometry
     */
    gsComposedGeometry( const gsFunction<T> & composition,
                        const gsGeometry<T> & geom);

    // /// Copy constructor (makes deep copy)
    // gsComposedGeometry(const gsComposedGeometry& other);

    // /// Move constructor
    // gsComposedGeometry( gsComposedGeometry&& other );

    // /// Assignment operator
    // gsComposedGeometry& operator= ( const gsComposedGeometry& other );

    // /// Move assignment operator
    // gsComposedGeometry& operator= ( gsComposedGeometry&& other );

    /// See \ref gsGeometry::domainDim for more details
    short_t domainDim() const override;

    /// See \ref gsGeometry::targetDim for more details
    using Base::targetDim;

    /// See \ref gsGeometry::eval_into for more details
    void compute(const gsMatrix<T> & in, gsFuncData<T> & out) const override;

    /**
     * @brief      Gives the control point derivatives of \a *this. See gsFunction for more details
     *
     * @param[in]  points  The points in the parameter domain (of the composition)
     * @param      result  The control point derivatives
     */
    void control_deriv_into(const gsMatrix<T> & points, gsMatrix<T> & result) const override;

    /// Evaluates the mesh
    void evaluateMesh(gsMesh<T>& mesh) const override;

    GISMO_OVERRIDE_BASIS_ACCESSORS;

    const CompositionT & composition() const;
          CompositionT & composition()      ;

protected:
    // Map from parametric domain to geometry
    typename CompositionT::Ptr m_composition;
    using Base::m_basis;

    // for compute();
    using Base::m_coefs;

    // Map from composition to geometry
    typename gsGeometry<T>::Ptr m_geom;

    short_t m_domainDim;
};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsComposedGeometry.hpp)
#endif
