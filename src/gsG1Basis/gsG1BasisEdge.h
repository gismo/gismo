/** @file gsG1BasisEdge.h

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
        m_gD = gluingData;

        // Computing the G1 - basis function at the edge
        // Spaces for computing the g1 basis
        index_t m_r = m_optionList.getInt("regularity");

        gsBSplineBasis<> basis_1 = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(0).component(1)); // v
        gsBSplineBasis<> basis_2 = dynamic_cast<gsBSplineBasis<> &>(m_mb.basis(1).component(0)); // u

        index_t m_p; // Minimum degree at the interface
        if (basis_1.degree() >= basis_2.degree())
            m_p = basis_2.maxDegree();
        else
            m_p = basis_1.maxDegree();

        // first,last,interior,mult_ends,mult_interior
        gsKnotVector<T> kv_tilde(0,1,0,m_p+1,m_p-1-m_r); // p,r+1 //-1 bc r+1
        gsBSplineBasis<> basis_plus(kv_tilde);

        if (basis_1.numElements() <= basis_2.numElements()) //
            for (size_t i = basis_1.degree()+1; i < basis_1.knots().size() - (basis_1.degree()+1); i = i+(basis_1.degree()-m_r))
                basis_plus.insertKnot(basis_1.knot(i),m_p-1-m_r);
        else
            for (size_t i = basis_2.degree()+1; i < basis_2.knots().size() - (basis_2.degree()+1); i = i+(basis_2.degree()-m_r))
                basis_plus.insertKnot(basis_2.knot(i),m_p-1-m_r);

        m_basis_plus = basis_plus;


        // TODO

    }






protected:

    // Input
    gsMultiPatch<T> m_mp;
    gsMultiBasis<T> m_mb;
    gsOptionList m_optionList;

    // Gluing data
    gsGluingData<T> m_gD;

    // Basis for getting the G1 Basis
    gsBSplineBasis<> m_basis_plus;
    gsBSplineBasis<> m_basis_minus;

    // Basis for the G1 Basis
    gsMultiBasis<> m_basis_g1;

}; // class gsG1BasisEdge


} // namespace gismo