/** @file gsPolyField.h

    @brief Provides declaration of the PolyField class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Shamanskiy
*/

#pragma once

#include <gsCore/gsFunctionSet.h>
#include <gsCore/gsMultiPatch.h>
#include <gsCore/gsGeometry.h>

namespace gismo
{

/**
 * \brief An object containing several scalar or vector fields defined on the same geometry.
 *
 * The class extends the idea behind gsField to a case of several fields. The primary purpose
 * of the class is to serve as a convinient container which can be used to plot several
 * physical fields defined on the same geometry to one Paraview file which can be useful
 * for problems like incompressible Navier-Stokes (velocity + pressure), elasticity (deformation + stresses),
 * thermal expansion (temperature + deformation), etc. Can be passed as an input argument to gsWriteParaview.
 *
 * \ingroup Core
 */
template <class T>
class gsPolyField
{

public:

    gsPolyField(const gsGeometry<T> & geo)
        : m_patches(memory::make_shared_not_owned(&geo))
    {}

    gsPolyField(const gsMultiPatch<T> & patches)
        : m_patches(memory::make_shared_not_owned(&patches))
    {}

    //=========================================================================//

    void addField(const std::string & name, const gsMultiPatch<T> & f)
    {
        GISMO_ASSERT(std::find(m_names.begin(),m_names.end(),name) == m_names.end(),
                     "Field with the name *" + name + "* already exists.\n");
        GISMO_ASSERT(nPatches() == f.nPatches(),
                     "The geometry has " + util::to_string(nPatches()) + " patch(es) and field *"
                     + name + "* has " + util::to_string(f.nPatches()) + " patch(es).\n");
        m_fields.push_back(memory::make_shared_not_owned(&f));
        m_names.push_back(name);
        m_isparams.push_back(true);
    }

    void addField(const std::string & name, const gsFunctionSet<T> & f, const bool isparam = false)
    {
        GISMO_ASSERT(std::find(m_names.begin(),m_names.end(),name) == m_names.end(),
                     "Field with the name *" + name + "* already exists.\n");
        GISMO_ASSERT(nPatches() == f.nPieces(),
                     "The geometry has " + util::to_string(nPatches()) + " patch(es) and field *"
                     + name + "* has " + util::to_string(f.nPieces()) + " patch(es).\n");
        m_fields.push_back(memory::make_shared_not_owned(&f));
        m_names.push_back(name);
        m_isparams.push_back(isparam);
    }

    /// \brief Number of patches in the geometry
    index_t nPatches() const { return m_patches->nPieces(); }

    /// \brief Number of fields in the container
    size_t nFields() const { return m_fields.size(); }

    /// \brief Returns the geometry patch number \a p.
    const gsGeometry<T> & patch(index_t p) const
    {
        GISMO_ASSERT(p >= 0 && p < nPatches(),
                     "Invalid patch index " + util::to_string(p) +
                     ". Number of patches is " + util::to_string(nPatches()) +".\n");
        return static_cast<const gsGeometry<T> &>(m_patches->piece(p));
    }

    /// \brief Returns the function corresponding to the geometry patch \a p from the field \a i.
    const gsFunction<T> & function(index_t p, size_t i) const
    {
        GISMO_ASSERT(p >= 0 && p < nPatches(),
                     "Invalid patch index " + util::to_string(p) +
                     ". Number of patches is " + util::to_string(nPatches()) +".\n");
        GISMO_ASSERT(i >= 0 && i < nFields(),
                     "Invalid field index " + util::to_string(i) +
                     ". Number of fields is " + util::to_string(nFields()) +".\n");
        return static_cast<const gsFunction<T> &>(m_fields[i]->piece(p));
    }

    const gsGeometry<T> & igaFunction(index_t p, size_t i) const
    {
        GISMO_ASSERT(p >= 0 && p < nPatches(),
                     "Invalid patch index " + util::to_string(p) +
                     ". Number of patches is " + util::to_string(nPatches()) +".\n");
        GISMO_ASSERT(i >= 0 && i < nFields(),
                     "Invalid field index " + util::to_string(i) +
                     ". Number of fields is " + util::to_string(nFields()) +".\n");
        GISMO_ASSERT(m_isparams[i],"Field " + util::to_string(i) + " is not parametric.\n");
        return static_cast<const gsGeometry<T> &>(m_fields[i]->piece(p));
    }

    /// \brief Returns the name of the field number \a i.
    const std::string & name(size_t i) const
    {
        GISMO_ASSERT(i >= 0 && i < nFields(),
                     "Invalid field index " + util::to_string(i) +
                     ". Number of fields is " + util::to_string(nFields()) +".\n");
        return m_names[i];
    }

    /// \brief Returns true if the field number \a i is parametric.
    bool isParametric(size_t i) const
    {
        GISMO_ASSERT(i >= 0 && i < nFields(),
                     "Invalid field index " + util::to_string(i) +
                     ". Number of fields is " + util::to_string(nFields()) +".\n");
        return m_isparams[i];
    }

    /// \brief Returns true if the field number \a i is a BSpline function.
    bool isParametrized(size_t i) const // bad function.
    {
        GISMO_ASSERT(i >= 0 && i < nFields(),
                     "Invalid field index " + util::to_string(i) +
                     ". Number of fields is " + util::to_string(nFields()) +".\n");
         { return m_isparams[i] && dynamic_cast<const gsGeometry<T>*>(&m_fields[i]->piece(0));}
    }

    /// \brief Returns the function corresponding to the geometry patch \a p from the field named \a name.
    const gsFunction<T> & function(index_t p, const std::string & name) const
    {
        GISMO_ASSERT(p >= 0 && p < nPatches(),
                     "Invalid patch index " + util::to_string(p) +
                     ". Number of patches is " + util::to_string(nPatches()) +".\n");
        std::vector<std::string>::const_iterator it = std::find(m_names.begin(),m_names.end(),name);
        GISMO_ASSERT(it != m_names.end(), "Invalid field name *" + name + "*.\n");
        return function(p,it - m_names.begin());
    }

    /// \brief Returns the function corresponding to the geometry patch \a p from the field named \a name.
    const gsGeometry<T> & igaFunction(index_t p, const std::string & name) const
    {
        GISMO_ASSERT(p >= 0 && p < nPatches(),
                     "Invalid patch index " + util::to_string(p) +
                     ". Number of patches is " + util::to_string(nPatches()) +".\n");
        std::vector<std::string>::const_iterator it = std::find(m_names.begin(),m_names.end(),name);
        GISMO_ASSERT(it != m_names.end(), "Invalid field name *" + name + "*.\n");
        return igaFunction(p,it - m_names.begin());
    }

    /// \brief Returns true if the field named \a name is parametric.
    bool isParametric(const std::string & name) const
    {
        std::vector<std::string>::const_iterator it = std::find(m_names.begin(),m_names.end(),name);
        GISMO_ASSERT(it != m_names.end(), "Invalid field name *" + name + "*.\n");
        return isParametric(it - m_names.begin());
    }


private:

    /// Geometry object
    typename gsFunctionSet<T>::Ptr m_patches;

    /// Fields. The numbering of elements is consistent across the following three vectors (m_fields, m_names, m_isparams),
    /// meaning that n_names[i] contains the name of the field in m_fields[i] and m_isparams[i] indicates its parametric type.
    std::vector<typename gsFunctionSet<T>::Ptr> m_fields;
    /// Names of the fields, e.g. velocity or pressure. Used to name fields in Paraview.
    /// Can also be used to access a specific field (if its name is known).
    std::vector<std::string> m_names;
    /// Flags indicating whether a corresponding field is parametric or not.
    /// If yes, the field has to be evaluated on the parametric domain.
    /// If no, the field has to be evaluated on the physical domain.
    std::vector<bool> m_isparams;

};

} // namespace gismo ends

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPolyField.hpp)
#endif
