/** @file gsGluingData.h

    @brief Compute the gluing data for one interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/gsApproxC1Utils.h>


namespace gismo
{

template<short_t d, class T>
class gsApproxGluingData
{
private:
    typedef typename std::vector<gsPatchReparameterized<d,T>> C1AuxPatchContainer;

    /// Shared pointer for gsApproxGluingData
    typedef memory::shared_ptr<gsApproxGluingData> Ptr;

    /// Unique pointer for gsApproxGluingData
    typedef memory::unique_ptr<gsApproxGluingData> uPtr;

public:
    gsApproxGluingData()
    { }


    gsApproxGluingData(C1AuxPatchContainer const & auxPatchContainer,
                       gsOptionList const & optionList,
                       std::vector<index_t> sidesContainer,
                       std::vector<bool> isInterface = std::vector<bool>{})
        : m_auxPatches(auxPatchContainer), m_optionList(optionList)
    {
        alphaSContainer.resize(2);
        betaSContainer.resize(2);
        if (m_auxPatches.size() == 2) // Interface
        {
            setGlobalGluingData(1,sidesContainer[1], 0); // u
            setGlobalGluingData(0,sidesContainer[0], 1); // v
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
    gsBSplineBasis<T> bsp_gD = dynamic_cast<gsBSplineBasis<T>&>(m_auxPatches[patchID].getBasisRotated().getHelperBasis(globalSide-1, 3));

    //! [Problem setup]
    gsSparseSolver<real_t>::LU solver;
    gsExprAssembler<> A(1,1);

    typedef gsExprAssembler<>::variable    variable;
    typedef gsExprAssembler<>::space       space;
    typedef gsExprAssembler<>::solution    solution;

    // Elements used for numerical integration
    gsMultiBasis<T> BsplineSpace(bsp_gD);
    A.setIntegrationElements(BsplineSpace);
    gsExprEvaluator<> ev(A);

    // Set the discretization space
    space u = A.getSpace(BsplineSpace);

    gsBoundaryConditions<> bc_empty;
    u.setup(bc_empty, dirichlet::homogeneous, 0);
    A.initSystem();

    gsAlpha<real_t> alpha(m_auxPatches[patchID].getPatchRotated(), dir);
    auto aa = A.getCoeff(alpha);

    A.assemble(u * u.tr(),u * aa);

    solver.compute( A.matrix() );
    gsMatrix<> solVector = solver.solve(A.rhs());

    solution u_sol = A.getSolution(u, solVector);
    gsMatrix<> sol;
    u_sol.extract(sol);

    gsGeometry<>::uPtr tilde_temp;
    tilde_temp = bsp_gD.makeGeometry(sol);
    alphaSContainer[dir] = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    gsBeta<real_t> beta(m_auxPatches[patchID].getPatchRotated(), dir);
    auto bb = A.getCoeff(beta);
    A.initSystem();

    A.assemble(u * u.tr(),u * bb);

    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());

    solution u_sol2 = A.getSolution(u, solVector);
    u_sol2.extract(sol);

    tilde_temp = bsp_gD.makeGeometry(sol);
    betaSContainer[dir] = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

} // setGlobalGluingData


} // namespace gismo