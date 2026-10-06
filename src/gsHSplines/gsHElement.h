/** @file gsHElement.h

    @brief Provides a class representing an element in a hierarchical mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsHSplines/gsElement.h>

namespace gismo
{

template <short_t d, class T>
class gsHElement : public gsElement<d,T>
{

    using Base = gsElement<d,T>;

// Comparison and equality operators
public:

    struct Compare
    {
        bool operator()(const gsHElement<d,T> & a, const gsHElement<d,T> & b) const;
    };

    struct Equal
    {
        bool operator()(const gsHElement<d,T> & a, const gsHElement<d,T> & b) const;
    };

public:

    typedef typename Base::box_t box_t;
    typedef typename Base::point_t point_t;
    typedef size_t level_t;

public:

    /// @brief Default constructor
    gsHElement();

    /// @brief Constructor with lower and upper corner and level
    /// @param low  Lower corner of the box in element indices
    /// @param upp  Upper corner of the box in element indices
    /// @param level  Level of the element in the hierarchy
    /// @param pid  Patch ID of the element, if applicable
    gsHElement(const point_t & low, const point_t & upp, level_t level, size_t pid = 0);

    /// @brief Constructor with box and level
    /// @param box  Box of the element, defined by its lower and upper corner in element indices
    /// @param level  Level of the element in the hierarchy
    /// @param pid  Patch ID of the element, if applicable
    gsHElement(const box_t & box, level_t level, index_t pid = 0);

    /// @brief Constructor from a gsElement
    /// @param element  The gsElement to construct from
    /// @param level  Level of the element in the hierarchy
    gsHElement(const Base & element, level_t level);

    /// @brief Copy constructor
    gsHElement(const gsHElement<d,T> & element) = default;

    /// @brief Move constructor
    gsHElement(gsHElement<d,T> && element) noexcept = default;

    /// @brief Copy assignment operator
    gsHElement<d,T>& operator=(const gsHElement<d,T> & element) = default;

    /// @brief Move assignment operator
    gsHElement<d,T>& operator=(gsHElement<d,T> && element) noexcept = default;

public:

    /// Get the level of the element
    level_t level() const
    {
        return m_level;
    }

public:

    /// Print the element information
    virtual std::ostream& print( std::ostream& os ) const;

protected:
    level_t m_level; ///< The level of the element.
    using Base::m_patch; ///< The patch ID of the element, inherited from gsElement.
};

template<short_t d, class T>
std::ostream& operator<<( std::ostream& os, const gsHElement<d,T>& b )
{
    return b.print( os );
}

template <short_t d, class T>
bool operator==(const gsHElement<d,T> & a, const gsHElement<d,T> & b)
{
    typename gsHElement<d,T>::Equal equal;
    return equal(a, b);
}

template <short_t d, class T>
bool operator!=(const gsHElement<d,T> & a, const gsHElement<d,T> & b)
{
    typename gsHElement<d,T>::Equal equal;
    return !equal(a, b);
}

template <short_t d, class T>
bool operator<(const gsHElement<d,T> & a, const gsHElement<d,T> & b)
{
    typename gsHElement<d,T>::Compare compare;
    return compare(a, b);
}

}; //namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsHElement.hpp)
#endif