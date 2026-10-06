/** @file gsElementHelper.hpp

    @brief Provides a class for creating and managing elements in a mesh.

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
    gsElementHelper<d,T>::gsElementHelper(const gsTensorBSplineBasis<d,T> & basis)
    :
    m_basis(basis)
    {
    }

    template <short_t d, class T>
    typename gsElementHelper<d,T>::element_t gsElementHelper<d,T>::toElement(const gsVector<T,d> & low, const gsVector<T,d> & upp, index_t patch) const
    {
        GISMO_ASSERT((low.array() <= upp.array()).all(), "Lower corner must be less than or equal to upper corner.");
        gsVector<index_t, d> lowIdx, uppIdx;
        for(index_t j = 0; j < d; j++)
        {
            // Convert the parameter coordinates to (unique) knot indices
            const gsKnotVector<T> & kv = m_basis.knots(j);
            lowIdx(j) = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd(),low[j] ) - 1).uIndex();
            uppIdx(j) = (std::upper_bound(kv.domainUBegin(), kv.domainUEnd()+1,upp[j] ) - 1).uIndex();
        }
        return element_t(lowIdx, uppIdx, patch);
    }

    template <short_t d, class T>
    typename gsElementHelper<d,T>::ElementContainer gsElementHelper<d,T>::toElements(index_t patch) const
    {
        ElementContainer elements;
        for (const auto & elem : m_basis.domain()->allElements())
        {
            // Convert knot coordinates to indices
            elements.insert(this->toElement(elem.lowerCorner(), elem.upperCorner(), patch));
        }
        return elements;
    }

    template <short_t d, class T>
    gsMatrix<T> gsElementHelper<d,T>::toBox(const typename  gsElementHelper<d,T>::element_t & element) const
    {
        gsMatrix<T> box(d, 2);
        const point_t & low = element.lowerCorner();
        const point_t & upp = element.upperCorner();
        for(index_t j = 0; j < d; j++)
        {
            // Convert the knot indices back to parameter coordinates
            const gsKnotVector<T> & kv = m_basis.knots(j);
            box(j,0) = kv.uValue(low(j));
            box(j,1) = kv.uValue(upp(j));
        }
        return box;
    }

    template <short_t d, class T>
    gsMatrix<T> gsElementHelper<d,T>::toBoxes(const ElementContainer & elements) const
    {
        gsMatrix<T> boxes(d, 2*elements.size());
        index_t i = 0;
        for (const auto & elem : elements)
        {
            // Convert the knot indices back to coordinates
            boxes.middleCols(i,2) = this->toBox(elem);
            i += 2;
        }
        return boxes;
    }

    template <short_t d, class T>
    typename gsElementHelper<d,T>::box_t gsElementHelper<d,T>::getSupportExtension(const typename  gsElementHelper<d,T>::element_t & element) const
    {
        box_t result;
        const point_t & low = element.lowerCorner();
        const point_t & upp = element.upperCorner();

        short_t degree;
        // short_t support;
        // short_t halfsupport;
        for( int i = 0; i < d; ++i )
        {
            const gsKnotVector<T> & knots = m_basis.knots(i);
            degree = m_basis.degree(i);
            // support = degree + 1; // support is degree + 1
            // halfsupport = (support + 2 - 1) / 2; // support / 2, rounded up
            // gsDebug<<"knots = "<<knots.asMatrix()<<", low: "<<low(i)<<", upp: "<<upp(i)
            //        <<", degree = "<<degree<<", support = "<<support<<", halfsupport = "<<halfsupport<<"\n";

            result(i,0) = math::max(0,low(i) - degree);
            result(i,1) = math::min((index_t)knots.uSize()-1, upp(i) + degree);
        }
        return result;
    }

    template <short_t d, class T>
    typename gsElementHelper<d,T>::ElementContainer gsElementHelper<d,T>::explode(const point_t & low, const point_t & upp, index_t patch) const
    {
        GISMO_ASSERT((low.array() < upp.array()).all(), "Lower corner must be less than or equal to upper corner.");
        ElementContainer result;

        // Split the box into all combinations of knot indices
        point_t curr, currpp, start, end;
        start = curr = low;
        end = upp;
        do
        {
            currpp.array() = curr.array() + 1; // upper corner of the current box
            // Create the element and add it to the result
            result.emplace(curr, currpp, patch);

        } while (nextLexicographic(curr, start, end));
        return result;
    }

    template <short_t d, class T>
    typename gsElementHelper<d,T>::ElementContainer gsElementHelper<d,T>::explode(const box_t & box, index_t patch) const
    {
        return explode(box.col(0), box.col(1), patch);
    }

}; //namespace gismo
