/** @file gsPointDomain.h

    @brief Iterator over the elements of a tensor-structured grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDomain.h>

namespace gismo
{

template<class T, int D>
class gsPointDomain : public gsDomain<T>
{
private:
    // typedef typename std::vector<T>::const_iterator  uiter;
    typename gsTensorBasis<d,T>::domainIter domainIter;

public:

    gsPointDomain(const gsMatrix<T> & points)
    :
    m_points(points),
    d( points.rows() )
    {
    }

    gsDomainIterator<T>::uPtr iterator(index_t i, const boxSid s = boundary::none) override
    {
        GISMO_ENSURE(s==boundary::none, "Point domain does not have boundary elements.");
        return domainIter(new gsPointDomainIterator<T,d>(*this));
    }

    // Look at gsBasis class for a description
    size_t numElements(boxSide const & s = boundary::none) const override
    {
        GISMO_ENSURE(s==boundary::none, "Point domain does not have boundary elements.");
        return m_points.cols();
    }

    short_t dim() const override { return d; }

    // NOTE: Do we need this one? Does it make sense?
    // gsMatrix<T> boundingBox() const override
    // {
    //     gsMatrix<T> result(2,d);
    //     for (short_t i = 0; i < d; ++i)
    //     {
    //         result(0,i) = *std::min_element(m_points.row(i).data(), m_points.row(i).data()+m_points.cols());
    //         result(1,i) = *std::max_element(m_points.row(i).data(), m_points.row(i).data()+m_points.cols());
    //     }
    //     return result;
    // }

    // virtual gsMesh<T> mesh() const override
    // {
    //     // gsMesh<T> mesh;
    //     // mesh.setDimension(d);
    //     // mesh.setBasis(m_basis);
    //     // return mesh;
    // }

// Specific for gsPointDomain
const gsMatrix<T> & points() const { return m_points; }

protected:
    const gsMatrix<T> & m_points;

};

} // namespace gismo
