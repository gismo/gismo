/** @file gsApproxC1Edge.h

    @brief Creates the (approx) C1 Edge space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once

#include <gsUnstructuredSplines/gsContainerBasis.h>
#include <gsUnstructuredSplines/gsPatchReparameterized.h>

#include <gsUnstructuredSplines/gsApproxGluingData.h>

#include <gsUnstructuredSplines/gsApproxC1Utils.h>

namespace gismo
{


template<short_t d, class T>
class gsApproxC1Edge
{

private:
    typedef gsContainerBasis<d, T> Basis;
    typedef typename std::vector<Basis> BasisContainer;
    typedef typename std::vector<gsPatchReparameterized<d,T>> C1AuxPatchContainer;

    /// Shared pointer for gsApproxC1Edge
    typedef memory::shared_ptr<gsApproxC1Edge> Ptr;

    /// Unique pointer for gsApproxC1Edge
    typedef memory::unique_ptr<gsApproxC1Edge> uPtr;


public:
    /// Empty constructor
    ~gsApproxC1Edge() { }


    gsApproxC1Edge(gsMultiPatch<T> const & mp,
                   BasisContainer & bases,
                const boundaryInterface & item,
                size_t & numInt,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_optionList(optionList)
    {
        m_auxPatches.clear();
        m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(item.first().patch), m_bases[item.first().patch]));
        m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(item.second().patch), m_bases[item.second().patch]));

        std::vector<patchSide> sidesContainer(2);
        sidesContainer[0] = item.first();
        sidesContainer[1] = item.second();

        reparametrizeInterfacePatches();

        compute(sidesContainer);

        if (m_optionList.getSwitch("plot"))
        {
            std::string fileName;
            std::string basename = "InterfaceBasisFunctions" + util::to_string(numInt);
            gsParaviewCollection collection(basename);

            for (size_t i = 0; i< basisEdgeResult[0].nPatches(); i++)
            {
                // First Interface Side
                fileName = basename + "_0_" + util::to_string(i);
                gsField<> temp_field(m_mp.patch(item.first().patch), basisEdgeResult[0].patch(i));
                gsWriteParaview(temp_field, fileName, 5000);
                collection.addTimestep(fileName, i, "0.vts");
                // Second Interface Side
                fileName = basename + "_1_" + util::to_string(i);
                gsField<> temp_field_1(m_mp.patch(item.second().patch), basisEdgeResult[1].patch(i));
                gsWriteParaview(temp_field_1, fileName, 5000);
                collection.addTimestep(fileName, i, "0.vts");
            }
            collection.save();

            //gsWriteParaview(basisEdgeResult[0], "interface_basis", 20000);
        }
    }

    gsApproxC1Edge(gsMultiPatch<T> const & mp,
                   BasisContainer & bases,
                const patchSide & item,
                size_t & numBdy,
                const gsOptionList & optionList)
                : m_mp(mp), m_bases(bases), m_optionList(optionList)
    {
        m_auxPatches.clear();
        m_auxPatches.push_back(gsPatchReparameterized<d,T>(m_mp.patch(item.patch), m_bases[item.patch]));

        std::vector<patchSide> sidesContainer(1);
        sidesContainer[0] = item;

        reparametrizeSinglePatch(item.side().index());

        compute(sidesContainer);

        if (m_optionList.getSwitch("plot")) {
            std::string fileName;
            std::string basename = "BoundaryBasisFunctions" + util::to_string(numBdy);
            gsParaviewCollection collection(basename);

            for (size_t i = 0; i < basisEdgeResult[0].nPatches(); i++) {
                // First Interface Side
                fileName = basename + "_0_" + util::to_string(i);
                gsField<> temp_field(m_mp.patch(item.patch), basisEdgeResult[0].patch(i));
                gsWriteParaview(temp_field, fileName, 5000);
                collection.addTimestep(fileName, i, "0.vts");
            }
            collection.save();
        }
    }

    std::vector<gsMultiPatch<T>> getEdgeBasis() { return basisEdgeResult; };

protected:

    // Input
    gsMultiPatch<T> const & m_mp;
    BasisContainer & m_bases;

    const gsOptionList & m_optionList;

    // Need for rotation, etc.
    C1AuxPatchContainer m_auxPatches;

    // Store temp solution
    std::vector<gsMultiPatch<T>> basisEdgeResult;

private:

    // Compute topology
    // After computeTopology() the patches will have the same patch-index as the position-index in auxGeom
    // EXAMPLE: global patch-index-order inside auxGeom: [2, 3, 4, 1, 0]
    //          in auxTop: 2->0, 3->1, 4->2, 1->3, 0->4
    void computeAuxTopology();

    void reparametrizeInterfacePatches();

    void reparametrizeSinglePatch(index_t side);

    void compute(std::vector<patchSide> & sidesContainer);

}; // Class gsApproxC1Edge

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsApproxC1Edge.hpp)
#endif
