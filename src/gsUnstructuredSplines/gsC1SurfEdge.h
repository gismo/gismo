/** @file gsC1SurfEdge.h

    @brief Creates the (approx) C1 Edge space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Farahat & P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/gsG1AuxiliaryPatch.h>
#include <gsUnstructuredSplines/gsC1SurfBasisEdge.h>

#include <gsUnstructuredSplines/gsC1SurfGluingData.h>

namespace gismo
{
template<short_t d, class T>
class gsC1SurfEdge
{

private:

    /// Shared pointer for gsC1SurfEdge
    typedef memory::shared_ptr<gsC1SurfEdge> Ptr;

    /// Unique pointer for gsC1SurfEdge
    typedef memory::unique_ptr<gsC1SurfEdge> uPtr;


public:
    /// Empty constructor
    ~gsC1SurfEdge() { }

    gsC1SurfEdge(const gsMultiPatch<> & mp, const boundaryInterface & item){
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(item.first().patch), item.first().patch));
        auxGeom.push_back(gsG1AuxiliaryPatch(mp.patch(item.second().patch), item.second().patch));
    }

    gsC1SurfEdge(const gsMultiPatch<> & sp, const patchSide & item){
        auxGeom.push_back(gsG1AuxiliaryPatch(sp.patch(item.patch), item.patch));
    }

    void computeG1InterfaceBasis()
    {
        basisEdgeResult.clear();

        gsMultiPatch<> mp_init;
        mp_init.addPatch(auxGeom[0].getPatch());// Right -> 0 ====> v along the interface
        mp_init.addPatch(auxGeom[1].getPatch()); // Left -> 1 ====> u along the interface

        gsMultiPatch<> test_mp(reparametrizeInterface()); // auxGeom contains now the reparametrized geometry
        gsMultiBasis<> test_mb(test_mp);
        gsMultiPatch<> g1Basis_0, g1Basis_1;

        gsC1SurfGluingData<real_t> g1BasisEdge(test_mp, test_mb);
        gsC1SurfBasisEdge<real_t> g1BasisEdge_0(test_mp.patch(0), test_mb.basis(0), 1, false, g1BasisEdge);
        gsC1SurfBasisEdge<real_t> g1BasisEdge_1(test_mp.patch(1), test_mb.basis(1), 0, false, g1BasisEdge);
        g1BasisEdge_0.setG1BasisEdge(g1Basis_0);
        g1BasisEdge_1.setG1BasisEdge(g1Basis_1);

//      Patch 0 -> Right
        auxGeom[0].parametrizeBasisBack(g1Basis_0);
//      Patch 1 -> Left
        auxGeom[1].parametrizeBasisBack(g1Basis_1);

        basisEdgeResult.push_back(auxGeom[0].getG1Basis());
        basisEdgeResult.push_back(auxGeom[1].getG1Basis());
    }

    void computeG1BoundaryBasis(const int boundaryInd)
    {
        basisEdgeResult.clear();

        gsMultiPatch<> test_mp(reparametrizeBoundary(boundaryInd));
        gsMultiBasis<> test_mb(test_mp);
        gsMultiPatch<> g1Basis_edge;

        gsC1SurfGluingData<real_t> bdyGD; // Empty constructor creates the sol and solBeta in a suitable way to manage the GD on the boundary
        gsC1SurfBasisEdge<real_t> g1BasisEdge(test_mp, test_mb, 1, true, bdyGD);
        g1BasisEdge.setG1BasisEdge(g1Basis_edge);

        auxGeom[0].parametrizeBasisBack(g1Basis_edge);

        basisEdgeResult.push_back(auxGeom[0].getG1Basis());
    }

    gsG1AuxiliaryPatch & getSinglePatch(const unsigned i)
    {
        return auxGeom[i];
    }

    std::vector<gsMultiPatch<T>> getBasis(){return basisEdgeResult;}

protected:

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisEdgeResult;

    std::vector<gsG1AuxiliaryPatch> auxGeom;

private:

    // Compute topology
    // After computeTopology() the patches will have the same patch-index as the position-index in auxGeom
    // EXAMPLE: global patch-index-order inside auxGeom: [2, 3, 4, 1, 0]
    //          in auxTop: 2->0, 3->1, 4->2, 1->3, 0->4
    gsMultiPatch<> computeAuxTopology();

    gsMultiPatch<> reparametrizeInterface();

    gsMultiPatch<> reparametrizeBoundary(index_t side);


}; // Class gsC1SurfEdge

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsC1SurfEdge.hpp)
#endif