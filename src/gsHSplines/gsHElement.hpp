/** @file gsHElement.hpp

    @brief Provides a class representing an element in a hierarchical mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsHSplines/gsHElement.h>

namespace gismo
{

    template <short_t d, class T>
    gsHElement<d,T>::gsHElement()
    :
    Base(),
    m_level(-1)
    {}

    template <short_t d, class T>
    gsHElement<d,T>::gsHElement(const point_t & low, const point_t & upp, level_t level, size_t pid)
    :
    Base(low, upp, pid),
    m_level(level)
    {
        GISMO_ASSERT(level >= 0, "Level must be non-negative.");
    }

    template <short_t d, class T>
    gsHElement<d,T>::gsHElement(const box_t & box, level_t level, index_t pid)
    :
    gsHElement(box.col(0), box.col(1), level, pid)
    {
    }

    template <short_t d, class T>
    gsHElement<d,T>::gsHElement(const Base & element, level_t level)
    :
    gsHElement(element.lowerCorner(), element.upperCorner(), level, element.patch())
    {
    }

    template <short_t d, class T>
    std::ostream& gsHElement<d,T>::print(std::ostream& os) const
    {
        os << "Patch: " << m_patch << ", Level: " << m_level << ", Box: ["<< "(" << this->lowerCorner().transpose() << "), (" << this->upperCorner().transpose() << ")]";
        return os;
    }

}; //namespace gismo
