/** @file gsElement.hpp

    @brief Provides a class representing an element in a mesh

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsHSplines/gsHElement.h>
#include <gsUtils/gsCombinatorics.h>

namespace gismo
{

    template <short_t d, class T>
    gsElement<d,T>::gsElement()
    :
    m_low(point_t::Constant(-1)),
    m_upp(point_t::Constant(-1)),
    m_patch(-1)
    {
    }

    template <short_t d, class T>
    gsElement<d,T>::gsElement(const point_t & low, const point_t & upp, size_t pid)
    :
    m_low(low),
    m_upp(upp),
    m_patch(pid)
    {
        GISMO_ASSERT((low.array()<= upp.array()).all(), "Lower corner must be less than or equal to upper corner.");
    }

    template <short_t d, class T>
    gsElement<d,T>::gsElement(const box_t & box, size_t pid)
    :
    gsElement(box.col(0), box.col(1), pid)
    {
    }

    template <short_t d, class T>
    gsElement<d,T>::gsElement(const gsHElement<d,T> & element)
    :
    gsElement(element.lowerCorner(), element.upperCorner(), element.patch())
    {
    }

    template <short_t d, class T>
    std::ostream& gsElement<d,T>::print(std::ostream& os) const
    {
        os << "Patch: " << m_patch << ", Box: ["<< "(" << this->lowerCorner().transpose() << "), (" << this->upperCorner().transpose() << ")]";
        return os;
    }

}; //namespace gismo
