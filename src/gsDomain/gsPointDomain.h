/** @file gsPointDomain.h

    @brief @hverhelst UPDATE

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomain.h>

namespace gismo
{

template<class T>
class gsPointDomain : public gsDomain<T>
{
private:
    typedef gsDomainIteratorWrapper<T> domainIter;

public:

    gsPointDomain(const gsMatrix<T> & points)
    :
    m_points(points)
    {
    }

    /// Documentation in gsDomain.h
    domainIter beginAll() const override
    {
        return domainIter(new gsPointDomainIterator<T>(*this));
    }

    /// Documentation in gsDomain.h
    domainIter beginBdr(const boxSide bs) const override
    {
        GISMO_ENSURE(bs==boundary::none, "Point domain does not have boundary elements.");
        return domainIter(new gsPointDomainIterator<T>(*this));
    }

    /// Documentation in gsDomain.h
    size_t numElements() const override
    {
        return m_points.cols();
    }

    /// Documentation in gsDomain.h
    size_t numElementsBdr(boxSide const & s = boundary::none) const override
    {
        GISMO_ENSURE(s==boundary::none, "Point domain does not have boundary elements.");
        return 0;
    }

    short_t dim() const override { return m_points.rows(); }

public:
    // Specific for gsPointDomain
    const gsMatrix<T> & points() const { return m_points; }

protected:
    const gsMatrix<T> & m_points;
};

} // namespace gismo
