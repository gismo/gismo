/** @file gsC1Basis.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): F. Buchegger
*/

/*
    TO DO (@Pascal)
        - Put everything from gsC1AgyrisBasis.h in here and in gsC1Basis.hpp
 */

#pragma once

#include<gsCore/gsBasis.h>
#include<gsIO/gsOptionList.h>

namespace gismo
{

template<short_t d,class T>
class gsC1Basis  : public gsBasis<T>
{

public:
    /// Shared pointer for gsC1Basis
    typedef memory::shared_ptr< gsC1Basis > Ptr;

    /// Unique pointer for gsC1Basis
    typedef memory::unique_ptr< gsC1Basis > uPtr;

    // gsC1Basis() : m_something(nullptr)
    // { }

    gsC1Basis(gsMultiPatch<T> const & mp, index_t patchID );

    // gsC1Basis( const gsC1Basis& other );

    virtual ~gsC1Basis();

    gsOptionList options() {return m_options;}
    void defaultOptions();
    void getOptions();

// implementations of gsBasis
public:
    short_t domainDim()
    {

    }

    // Data members
protected:
    /// Class options
    gsOptionList m_options;

    /// The multipatch
    const gsMultiPatch<T> m_patches;

    /// THe ID of the single basis
    index_t m_patchID;


};

}

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsC1Basis.hpp)
#endif
