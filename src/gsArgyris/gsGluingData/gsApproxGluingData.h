/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

# include <gsArgyris/gsGluingData/gsApproxGluingDataAssembler.h>


namespace gismo
{

template<short_t d, class T>
class gsApproxGluingData
{
private:
    typedef typename std::vector<gsC1ArgyrisAuxiliaryPatch<d,T>> ArgyrisAuxPatchContainer;

public:
    gsApproxGluingData()
    { }


    gsApproxGluingData(ArgyrisAuxPatchContainer const & auxPatchContainer,
                       gsOptionList const & optionList)
        : m_auxPatches(auxPatchContainer), m_optionList(optionList)
    {

        if (m_auxPatches.size() == 2)
        {
            setGlobalGluingData(0,m_auxPatches[0].side()); // Order is important!!!
            setGlobalGluingData(1,m_auxPatches[1].side());
        }
        else
            gsInfo << "SOMETHING WENT WRONG \n";

    }

    // Computed the gluing data globally
    void setGlobalGluingData(index_t patchID = 0,  index_t side = 1);

    gsBSpline<T> & alphaS(index_t patchID) { return alphaSContainer[patchID]; }
    gsBSpline<T> & betaS(index_t patchID) { return betaSContainer[patchID]; }

protected:

    // Spline space for the gluing data + multiPatch
    ArgyrisAuxPatchContainer m_auxPatches;

    const gsOptionList m_optionList;

    // Result
    std::vector<gsBSpline<T>> alphaSContainer, betaSContainer;

}; // class gsGluingData


template<short_t d, class T>
void gsApproxGluingData<d, T>::setGlobalGluingData(index_t patchID, index_t side)
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsBSplineBasis<T> bsp_gD = m_auxPatches[patchID].getArygrisBasisRotated().getBasisGluingData(side);
    gsInfo << "Gluing data: " << bsp_gD.knots().asMatrix() << "\n";

    index_t dir = patchID == 0 ? 1 : 0;

    gsApproxGluingDataAssembler<T> approxGluingDataAssembler(m_auxPatches[patchID].getPatch(), bsp_gD, dir, m_optionList);
    alphaSContainer.push_back(approxGluingDataAssembler.getAlphaS());
    betaSContainer.push_back(approxGluingDataAssembler.getBetaS());

    if (patchID == 0)
        gsWriteParaview(approxGluingDataAssembler.getAlphaS(), "alpha_R", 1000);
    if (patchID == 1)
        gsWriteParaview(approxGluingDataAssembler.getAlphaS(), "alpha_L", 1000);

    if (patchID == 0)
        gsWriteParaview(approxGluingDataAssembler.getBetaS(), "beta_R", 1000);
    if (patchID == 1)
        gsWriteParaview(approxGluingDataAssembler.getBetaS(), "beta_L", 1000);

} // setGlobalGluingData


} // namespace gismo

