/** @file gsElement.h

    @brief Provides a class representing an element in a mesh

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

template <short_t d, class T>
class gsElement
{
// Comparison and equality operators
public:

    struct Compare
    {
        bool operator()(const gsElement<d,T> & a, const gsElement<d,T> & b) const
        {
            return
             (a.patch() < b.patch())
            ||
            ((a.patch() == b.patch()) &&
             std::lexicographical_compare(  a.lowerCorner().begin(), a.lowerCorner().end(),
                                            b.lowerCorner().begin(), b.lowerCorner().end())   )
            ||
            ((a.patch() == b.patch()) &&
             (a.lowerCorner() == b.lowerCorner()) &&
             std::lexicographical_compare(  a.upperCorner().begin(), a.upperCorner().end(),
                                            b.upperCorner().begin(), b.upperCorner().end())    );
        }
    };

    struct Equal
    {
        bool operator()(const gsElement<d,T> & a, const gsElement<d,T> & b) const
        {
            return (a.patch() == b.patch()) && (a.lowerCorner() == b.lowerCorner()) && (a.upperCorner() == b.upperCorner());
        }
    };

public:
    typedef gsMatrix<index_t,d,2> box_t;
    typedef gsVector<index_t,d> point_t;

public:

    /// @brief Default constructor
    gsElement();

    /// @brief Constructor with lower and upper corner and patch ID
    /// @param low  Lower corner of the box in element indices
    /// @param upp  Upper corner of the box in element indices
    /// @param pid  Patch ID of the element, if applicable
    gsElement(const point_t & low, const point_t & upp, size_t pid = 0);

    /// @brief Constructor with box and patch ID
    /// @param box  Box of the element, defined by its lower and upper corner in element indices
    /// @param pid  Patch ID of the element, if applicable
    gsElement(const box_t & box, size_t pid = 0);

    /// @brief Constructor from a gsHElement
    /// @param element  The gsHElement to construct from
    gsElement(const gsHElement<d,T> & element);

    /// @brief Copy constructor
    gsElement(const gsElement<d,T> & element) = default;

    /// @brief Move constructor
    gsElement(gsElement<d,T> && element) noexcept = default;

    /// @brief Copy assignment operator
    gsElement<d,T>& operator=(const gsElement<d,T> & element) = default;

    /// @brief Move assignment operator
    gsElement<d,T>& operator=(gsElement<d,T> && element) noexcept = default;

public:

    /// @brief Get the patch ID of the element
    index_t patch() const
    {
        return m_patch;
    }

    /// @brief Get the lower corner of the box by reference
    const point_t & lowerCorner() const
    {
        return m_low;
    }

    /// @brief Get the upper corner of the box by reference
    const point_t & upperCorner() const
    {
        return m_upp;
    }

public:

    /// Print the element information
    virtual std::ostream& print( std::ostream& os ) const;

protected:
    // const gsMatrix<T,d,2> m_box; ///< The box of the element, defined by its lower and upper corner.
    point_t m_low; ///< The lower corner of the box.
    point_t m_upp; ///< The upper corner of the box.
    index_t m_patch; ///< The patch ID of the element, if applicable.
};

template<short_t d, class T>
std::ostream& operator<<( std::ostream& os, const gsElement<d,T>& b )
{
    return b.print( os );
}

template<short_t d, class T>
bool operator==(const gsElement<d,T>& a, const gsElement<d,T>& b)
{
    typename gsElement<d,T>::Equal equal;
    return equal(a, b);
}

template<short_t d, class T>
bool operator!=(const gsElement<d,T>& a, const gsElement<d,T>& b)
{
    typename gsElement<d,T>::Equal equal;
    return !equal(a, b);
}

template<short_t d, class T>
bool operator<(const gsElement<d,T>& a, const gsElement<d,T>& b)
{
    typename gsElement<d,T>::Compare compare;
    return compare(a, b);
}

}; //namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElement.hpp)
#endif