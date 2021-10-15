/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

# include <gsUnstructuredSplines/gsApproxGluingDataAssembler.h>


namespace gismo
{

template<short_t d, class T>
class gsApproxGluingData
{
private:
    typedef typename std::vector<gsPatchReparameterized<d,T>> C1AuxPatchContainer;

public:
    gsApproxGluingData()
    { }


    gsApproxGluingData(C1AuxPatchContainer const & auxPatchContainer,
                       gsOptionList const & optionList,
                       std::vector<index_t> sidesContainer = std::vector<index_t>{},
                       std::vector<bool> isInterface = std::vector<bool>{})
        : m_auxPatches(auxPatchContainer), m_optionList(optionList)
    {
        alphaSContainer.resize(2);
        betaSContainer.resize(2);
        if (m_auxPatches.size() == 2) // Interface
        {
            setGlobalGluingData(1,m_auxPatches[1].side(), 0); // u
            setGlobalGluingData(0,m_auxPatches[0].side(), 1); // v
        }
        else if (m_auxPatches.size() == 1) // Vertex
        {
            for (size_t dir = 0; dir < sidesContainer.size(); dir++)
            {
                index_t localSide = auxPatchContainer[0].getMapIndex(sidesContainer[dir]);
                //gsInfo << "Global: " << sidesContainer[dir] << " : " << localSide << "\n";
                index_t localDir = localSide < 3 ? 1 : 0;
                if(isInterface[dir]) // West
                    setGlobalGluingData(0, sidesContainer[dir], localDir);
                else
                {
                    // empty
                }
            }
        }
        else
            gsInfo << "Something went wrong \n";

    }

    // Computed the gluing data globally
    void setGlobalGluingData(index_t patchID = 0,  index_t side = 1, index_t dir = 1);

    gsBSpline<T> & alphaS(index_t patchID) { return alphaSContainer[patchID]; }
    gsBSpline<T> & betaS(index_t patchID) { return betaSContainer[patchID]; }

protected:

    // Spline space for the gluing data + multiPatch
    C1AuxPatchContainer m_auxPatches;

    const gsOptionList m_optionList;

    // Result
    std::vector<gsBSpline<T>> alphaSContainer, betaSContainer;

}; // class gsGluingData


template<short_t d, class T>
void gsApproxGluingData<d, T>::setGlobalGluingData(index_t patchID, index_t globalSide, index_t dir)
{
    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsBSplineBasis<T> bsp_gD = dynamic_cast<gsBSplineBasis<T>&>(m_auxPatches[patchID].getC1BasisRotated().getHelperBasis(globalSide-1, 3));

    gsApproxGluingDataAssembler<T> approxGluingDataAssembler(m_auxPatches[patchID].getPatch(), bsp_gD, dir, m_optionList);
    alphaSContainer[dir] = approxGluingDataAssembler.getAlphaS();
    betaSContainer[dir] = approxGluingDataAssembler.getBetaS();
/*
    if (m_auxPatches.size() == 1)
        gsInfo << m_auxPatches[patchID].getC1BasisRotated().getPatchID() << "\n";

    gsMatrix<> points(1,6);
    points << 0, 0.2, 0.4, 0.6, 0.8, 1;
    gsInfo << "Beta S " << approxGluingDataAssembler.getBetaS().eval(points) << "\n";
    gsInfo << "Alpha " << approxGluingDataAssembler.getAlphaS().eval(points) << "\n";

    gsMatrix<> ones = points;
    ones.setOnes();
    gsInfo << "Beta " << -29403.0/20000.0 * ones + 29403.0*points/10000.0 - 29403.0/20000.0 * points.cwiseProduct(points) << "\n";

    if (patchID == 0)
        gsWriteParaview(approxGluingDataAssembler.getAlphaS(), "alpha_R", 1000);
    if (patchID == 1)
        gsWriteParaview(approxGluingDataAssembler.getAlphaS(), "alpha_L", 1000);

    if (patchID == 0)
        gsWriteParaview(approxGluingDataAssembler.getBetaS(), "beta_R", 1000);
    if (patchID == 1)
        gsWriteParaview(approxGluingDataAssembler.getBetaS(), "beta_L", 1000);
*/
} // setGlobalGluingData


} // namespace gismo

