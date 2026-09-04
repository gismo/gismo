/** @file gsDomainFaceIterator.h

    @brief Abstract base for iteration over the codimension-1 faces
    (element-to-element interfaces) of a domain.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomainIterator.h>

namespace gismo
{

/**
 * @brief Abstract iterator over the interior faces of a domain: the
 * codimension-1 boxes shared by two neighbouring elements.
 *
 * A concrete implementation (e.g. gsTensorDomainFaceIterator, a future
 * gsHDomainFaceIterator) must respect the following contract:
 *
 * - iterates all directions, dir-major: all faces normal to direction 0
 *   first, then direction 1, and so on;
 * - lowerCorner()/upperCorner() have zero extent in direction();
 *   consequently volume() and getMinCellLength() (gsDomainIterator.h)
 *   degenerate to 0 on a face -- use getPerpendicularCellSize() (left
 *   neighbour thickness h_left) and getPerpendicularCellSizeRight()
 *   (h_right) instead;
 * - side() is boxSide(direction(), true), the upper side of the left
 *   (lower-index) element, so the outward normal of the left element,
 *   nv(G.left()), equals +e_direction;
 * - id() counts faces; leftElementId()/rightElementId() are flat element
 *   indices in the parent domain's numbering -- the two indices are never
 *   to be conflated with id();
 * - leftSign()/rightSign() are -1 interior, 0 cut, +1 exterior; domains
 *   without a trimming concept report -1 for every element;
 * - the count used for the end sentinel must equal the number of faces
 *   next() actually yields, and the zero-face case (no two elements share
 *   a face in any direction) must be safe: constructing an iterator and
 *   comparing it to the end sentinel must neither throw nor loop.
 *
 * \ingroup Domain
 */
template<class T>
class gsDomainFaceIterator : public gsDomainIterator<T>
{
public:
    explicit gsDomainFaceIterator(index_t _id = 0) : gsDomainIterator<T>(_id) { }

    /// Normal direction of the current face
    short_t direction() const { return this->side().direction(); }

    virtual const T getPerpendicularCellSize()      const override = 0;
    virtual const T getPerpendicularCellSizeRight() const override = 0;
    virtual size_t  leftElementId()  const override = 0;
    virtual size_t  rightElementId() const override = 0;
    virtual short_t leftSign()  const override { return -1; }
    virtual short_t rightSign() const override { return -1; }
};

} // namespace gismo
