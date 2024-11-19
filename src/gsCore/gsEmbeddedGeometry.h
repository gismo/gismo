/** @file gsGeometry.h

    @brief Provides declaration of Geometry abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsCore/gsBoundary.h>



namespace gismo
{

template<class T>
class gsEmbeddedGeometry : public gsFunction<T>
{

public:
    /// Shared pointer for gsGeometry
    typedef memory::shared_ptr< gsGeometry > Ptr;

    /// Unique pointer for gsGeometry
    typedef memory::unique_ptr< gsGeometry > uPtr;

    typedef typename gsMultiPatch<T>::PatchContainer PatchContainer;


public:

    gsEmbeddedGeometry()
    {}

    gsEmbeddedGeometry(const gsGeometry<T> & superdomain, const gsMultiPatch<T> & embedded)
    :
    m_superdomain(superdomain),
    m_embedded(embedded)
    {}

    gsEmbeddedGeometry(const gsGeometry<T> & superdomain, const PatchContainer & embedded)
    :
    m_superdomain(superdomain),
    m_embedded(embedded)
    {

    }

    gsEmbeddedGeometry(const gsGeometry<T> & superdomain, const gsGeometry<T> & embedded)
    :
    m_superdomain(superdomain),
    m_embedded(embedded)
    {

    }

    gsEmbeddedGeometry(const gsGeometry<T> & superdomain)
    :
    m_superdomain(superdomain)
    {

    }

    void addEmbedded(const gsGeometry<T> & embedded)
    {
        m_embedded.addPatch(embedded);
    }

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const override
    {
        // todo
    }

    short_t domainDim () const override
    {
        return m_superdomain.domainDim();
    }

    const gsGeometry<T> & superdomain() const { return m_superdomain; }
    const gsMultiPatch<T> & embedded() const { return m_embedded; }

    gsMatrix<T> quadraturePoints() const
    {
        //
    }


protected:

    const gsGeometry<T> & m_superdomain;
    gsMultiPatch<T> m_embedded;

};

#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsGeometry
   */
  void pybind11_init_gsGeometry(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace gismo


// #ifndef GISMO_BUILD_LIB
// #include GISMO_HPP_HEADER(gsEmbeddedGeometry.hpp)
// #endif
