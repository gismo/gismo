/** @file gsL2Projection.h

    @brief Class that performs an L2 projection

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/


#pragma once

#include <gsAssembler/gsExprAssembler.h>

namespace gismo {

/** \brief Class that performs an L2 projection

    \tparam T coefficient type
 */

template <class T>
struct gsL2Projection
{

protected:
    typedef gsExprAssembler<>::geometryMap geometryMap;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;
    typedef gsExprAssembler<>::element     element;

public:
    /**
     * @brief      Projects a \a source geometry onto \a basis and returns it in \a result
     *
     * @param[in]  basis   The basis to project on
     * @param[in]  source  The source geometry
     * @param      result  The result geometry
     */
    static void projectGeometry(    const gsMultiBasis<T> & basis,
                                    const gsMultiPatch<T> & source,
                                    gsMultiPatch<T> & result);

}; //struct

} // gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsL2Projection.hpp)
#endif
