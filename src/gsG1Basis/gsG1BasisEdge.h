/** @file gsG1Basis_mp.h

    @brief Provides assembler for a G1 Basis for multiPatch.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmüller
*/

#pragma once

#include <gsG1Basis/gsGluingData.h>


namespace gismo
{
template<class T>
class gsG1BasisEdge : public gsAssembler<T>
{
public:
    typedef gsAssembler<T> Base;

public:
    gsG1BasisEdge(gsMultiPatch<T> const & mp,
                 gsMultiBasis<T> & mb,
                 gsOptionList & optionList)
        : m_mp(mp), m_mb(mb), m_optionList(optionList)
    {
        // Computing the gluing data
        gsGluingData<T> gluingData(m_mp,m_mb,m_optionList);
        //m_gD = gluingData;

        // Computing the G1 - basis function at the edge
        // TODO

    }






protected:

    // Input
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    gsOptionList m_optionList;

    // Gluing data
    gsGluingData<T> m_gD;


}; // class gsG1BasisEdge


} // namespace gismo