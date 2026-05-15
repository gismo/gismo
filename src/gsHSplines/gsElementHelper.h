/** @file gsElementHelper.h

    @brief Provides a class for creating and managing elements in a mesh.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M.Verhelst
*/

#pragma once

#include <gsHSplines/gsElement.h>
#include <gsNurbs/gsTensorBSplineBasis.h>

namespace gismo
{

template <short_t d, class T>
class gsElementHelper
{
public:
    typedef          gsElement<d,T>             element_t;
    typedef typename gsElement<d,T>::box_t      box_t;
    typedef typename gsElement<d,T>::point_t    point_t;
    typedef typename std::set<element_t,typename element_t::Compare>     ElementContainer;

public:
    /// @brief Constructor
    /// @param basis The basis to use for element creation
    explicit gsElementHelper(const gsTensorBSplineBasis<d,T> & basis);

public:

    /// @brief  Convert a box defined by parametric lower and upper corners to an element
    /// @param low  Lower corner of the box in parameter coordinates
    /// @param upp  Upper corner of the box in parameter coordinates
    /// @param patch  Patch ID of the element, if applicable
    /// @return An element representing the box in the mesh
    element_t toElement(const gsVector<T,d> & low, const gsVector<T,d> & upp, index_t patch = 0) const;

    /// @brief  Convert the basis domain to a set of elements
    /// @param patch  Patch ID of the elements, if applicable
    /// @return A set of elements representing the basis domain
    ElementContainer toElements(index_t patch = 0) const;

    /// @brief Convert an element to a box defined by its lower and upper corners in parameter coordinates
    /// @param element The element to convert
    /// @return A box representing the element in parameter coordinates
    gsMatrix<T> toBox(const element_t & element) const;

    /// @brief Convert a set of elements to a set of boxes
    /// @param elements The set of elements to convert
    /// @return A matrix where each two columns represent a box defined by its lower and upper corners
    ///         in parameter coordinates
    gsMatrix<T> toBoxes(const ElementContainer & elements) const;

    /**
     * @brief Computes the support extension of a given element
     *
     * This function calculates the support extension box for a given element.
     * The support extension is the smallest box that contains all basis function
     * supports that have non-zero contribution to the element.
     *
     * The algorithm works by:
     * 1. Finding the knot indices that limit the element box in each parametric direction
     * 2. Computing the minimum and maximum support boundaries across these limit indices
     *
     * @param element The input element for which to compute the support extension
     * @return box_t A box representing the support extension of the element
     * @note This function assumes that the functions supported in the element have a support of degree + 1,
     *      hence the 'band' of elements around the element is of size degree.
     */
    box_t getSupportExtension(const element_t & element) const;

    /// @brief Explode a box defined by its lower and upper corners into a set of elements
    /// @param low  Lower corner of the box in element indices
    /// @param upp  Upper corner of the box in element indices
    /// @param patch  Patch ID of the elements, if applicable
    /// @return A set of elements that cover the box
    ElementContainer explode(const point_t & low, const point_t & upp, index_t patch = 0) const;

    /// @brief Explode a box defined by its lower and upper corners into a set of elements
    /// @param box  Box of the element, defined by its lower and upper corner in element indices
    /// @param patch  Patch ID of the elements, if applicable
    /// @return A set of elements that cover the box
    ElementContainer explode(const box_t & box, index_t patch = 0) const;

protected:
    const gsTensorBSplineBasis<d,T> & m_basis; ///< The basis of the elements.

};

}; //namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsElementHelper.hpp)
#endif