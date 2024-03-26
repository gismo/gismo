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

    typedef memory::shared_ptr< gsComposedGeometry > Ptr;
    typedef memory::unique_ptr< gsComposedGeometry > uPtr;

    GISMO_CLONE_FUNCTION(gsComposedGeometry)

public:
    gsComposedGeometry(const gsComposedBasis<T> & basis, const gsMatrix<T> & coefs)
    :
    Base(basis, give(coefs) ),
    m_domainDim(basis.domainDim())
    { }

    gsComposedGeometry(const gsFunction<T> & composition, const gsGeometry<T> & geom)
    :
    Base(gsComposedBasis<T>(composition,geom.basis()), give(geom.coefs()) ),
    m_domainDim(geom.domainDim())
    { }

    short_t domainDim() const override { return m_domainDim; }

    GISMO_BASIS_ACCESSORS;

protected:
    using Base::m_basis;
    short_t m_domainDim;
};
}
