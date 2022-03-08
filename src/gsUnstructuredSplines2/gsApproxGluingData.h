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
                       std::vector<patchSide> sidesContainer,
                       std::vector<bool> isInterface = std::vector<bool>{})
        : m_auxPatches(auxPatchContainer), m_optionList(optionList)
    {
        alphaSContainer.resize(2);
        betaSContainer.resize(2);
        if (m_auxPatches.size() == 2) // Interface
        {
            setGlobalGluingData(1,0); // u
            setGlobalGluingData(0,1); // v
        }
        else if (m_auxPatches.size() == 1 && sidesContainer.size() == 2) // Vertex
        {
            for (size_t dir = 0; dir < sidesContainer.size(); dir++)
            {
                index_t localSide = auxPatchContainer[0].getMapIndex(sidesContainer[dir].index());
                //gsInfo << "Global: " << sidesContainer[dir] << " : " << localSide << "\n";
                index_t localDir = localSide < 3 ? 1 : 0;
                if(isInterface[dir]) // West
                    setGlobalGluingData(0, localDir);
                else
                {
                    // empty
                }
            }
        }
        //else
        //    gsInfo << "I am here \n";

    }

    // Computed the gluing data globally
    void setGlobalGluingData(index_t patchID = 0,  index_t dir = 1);

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
void gsApproxGluingData<d, T>::setGlobalGluingData(index_t patchID, index_t dir)
{
    // Interpolate boundary yes or no //
    bool interpolate_boundary = false;
    // Interpolate boundary yes or no //

    // ======== Space for gluing data : S^(p_tilde, r_tilde) _k ========
    gsBSplineBasis<T> bsp_gD;
    createGluingDataSpace(m_auxPatches[patchID].getPatchRotated(), m_auxPatches[patchID].getBasisRotated().piece(0),
                          dir, bsp_gD, m_optionList.getInt("gluingDataDegree"), m_optionList.getInt("gluingDataSmoothness"));

    //! [Problem setup]
    gsSparseSolver<real_t>::LU solver;
    gsExprAssembler<T> A(1,1);

    // Elements used for numerical integration
    gsMultiBasis<T> BsplineSpace(bsp_gD);
    A.setIntegrationElements(BsplineSpace);
    gsExprEvaluator<> ev(A);

    gsAlpha<real_t> alpha(m_auxPatches[patchID].getPatchRotated(), dir);
    auto aa = A.getCoeff(alpha);

    // Set the discretization space
    auto u = A.getSpace(BsplineSpace);

    // Create Mapper
    gsDofMapper map(BsplineSpace);
    gsMatrix<index_t> act(2,1);
    act(0,0) = 0;
    act(1,0) = BsplineSpace[0].size()-1; // First and last
    if (interpolate_boundary)
        map.markBoundary(0, act); // Patch 0
    map.finalize();

    u.setupMapper(map);

    gsMatrix<T> & fixedDofs = const_cast<expr::gsFeSpace<T>&>(u).fixedPart();
    fixedDofs.setZero( u.mapper().boundarySize(), 1 );

    // For the boundary
    gsMatrix<> points_bdy(1,2);
    points_bdy << 0.0, 1.0;
    if (interpolate_boundary)
        fixedDofs = alpha.eval(points_bdy).transpose();

    A.initSystem();
    A.assemble(u * u.tr(), u * aa);

    solver.compute( A.matrix() );
    gsMatrix<> solVector = solver.solve(A.rhs());

    auto u_sol = A.getSolution(u, solVector);
    gsMatrix<> sol;
    u_sol.extract(sol);

    gsGeometry<>::uPtr tilde_temp;
    tilde_temp = bsp_gD.makeGeometry(sol);
    alphaSContainer[dir] = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

    gsBeta<real_t> beta(m_auxPatches[patchID].getPatchRotated(), dir);
    auto bb = A.getCoeff(beta);

    // For the boundary
    if (interpolate_boundary)
        fixedDofs = beta.eval(points_bdy).transpose();

    A.initSystem();
    A.assemble(u * u.tr(), u * bb);

    solver.compute( A.matrix() );
    solVector = solver.solve(A.rhs());

    auto u_sol2 = A.getSolution(u, solVector);
    u_sol2.extract(sol);

    tilde_temp = bsp_gD.makeGeometry(sol);
    betaSContainer[dir] = dynamic_cast<gsBSpline<T> &> (*tilde_temp);

} // setGlobalGluingData


} // namespace gismo
