/** @file gsDomainContainer.h

    @brief Container for multiple domains.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsDomain.h>

#include <gsUtils/gsCombinatorics.h>

#include <gsAssembler/gsGaussRule.h>

namespace gismo
{
// Documentation in gsDomain.h
template<class T>
class gsDomainContainer : public gsDomain<T>
{
private:

    typedef typename std::vector<gsDomain<T>::Ptr> domainContainer;

public:

    /** @brief Constructor from a \ref gsFunctionSet.
     *
     * @param multiBasis The multiBasis from which the domains are extracted.
     */
    gsDomainContainer(const gsFunctionSet<T> & multiBasis)
    :
    m_domains(multiBasis.size())
    {
        for (index_t i = 0; i < multiBasis.nPieces(); ++i)
        {
            m_domains[i] = multiBasis.basis(i).domain();
        }
    }

    /// See \ref gsDomain.h for documentation.
    gsDomainIterator<T>::uPtr iterator(index_t i, const boxSid s = boundary::none) override
    {
        return m_domains[i]->iterator(s);
    }

    /// See \ref gsDomain.h for documentation.
    size_t numElements(boxSide const & s = boundary::none) const override
    {
        size_t size = 0;
        for (index_t i = 0; i < m_domains.size(); ++i)
            size += m_domains[i]->numElements(s);
        return size;
    }

    /// See \ref gsDomain.h for documentation.
    short_t dim() const override { return m_domains.front().dim(); }

    /// See \ref gsDomain.h for documentation.
    gsMatrix<T> boundingBox() const override
    {
        GISMO_NO_IMPLEMENTATION;
    }

    /// See \ref gsDomain.h for documentation.
    virtual gsMesh<T> mesh() const override
    {
        // gsMesh<T> mesh;
        // mesh.setDimension(d);
        // mesh.setBasis(m_basis);
        // return mesh;
    }

    /// See \ref gsDomain.h for documentation.
    std::ostream &print(std::ostream &os) const override
    {
        os << "Domain container with " << m_domains.size() << " domains.";
        return os;
    }

protected:
    domainContainer m_domains;

}; // class gsDomainContainer


} // namespace gismo
