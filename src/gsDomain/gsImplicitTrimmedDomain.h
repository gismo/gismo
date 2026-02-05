/** @file gsTrimmedDomain.h

    @brief TODO

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsTrimmedDomain.h>


namespace gismo
{

/**
 * @brief TODO
 *
 * \ingroup Domain
 */

/// Class representing an implicit trimmed domain
template<short_t d, class T, class Z=size_t>
class gsImplicitTrimmedDomain : public gsTrimmedDomain<d,T,Z>
{
    static bool constexpr implicit = true;
    typedef gsTrimmedDomain<d,T,Z> Base;
private:
    typename gsFunction<T>::Ptr m_implFunction; // implicit function

public:
    gsImplicitTrimmedDomain(const gsFunction<T> & fnc,
                        const gsTensorBSplineBasis<d,T> & tbasis) :
    m_implFunction(memory::make_shared_not_owned(&fnc))
    {
        this->init(tbasis, 3);
    }

    inline gsVector<short_t> sign(const gsMatrix<T> & u)
    {
        gsVector<T> val = m_implFunction->eval(u).row(0);
        return gsVector<short_t>(val.array().sign().template cast<short_t>());
    }

    gsMatrix<T> boundingBox() const override
    {
        //always exists?
        return m_implFunction->support(); 
    }

};

} // namespace gismo
