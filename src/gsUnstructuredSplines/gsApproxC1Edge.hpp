/** @file gsApproxC1Edge.hpp

    @brief Creates the approx C1 space.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller & A. Farahat
*/

#pragma once

#include <gsUnstructuredSplines/gsApproxC1Edge.h>

#include <gsUnstructuredSplines/gsPatchReparameterized.h>

#include <gsUnstructuredSplines/gsApproxGluingData.h>

namespace gismo
{


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::computeAuxTopology()
    {
        for(unsigned i = 0; i <  m_auxPatches.size(); i++)
        {
            if(m_auxPatches[i].getPatchRotated().orientation() == -1)
            {
                m_auxPatches[i].swapAxis();
                //gsInfo << "Changed axis on patch: " << i << "\n";
            }
        }
    }


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::reparametrizeInterfacePatches()
    {
        computeAuxTopology();

        gsMultiPatch<> temp_mp;
        for(unsigned i = 0; i <  m_auxPatches.size(); i++)
            temp_mp.addPatch(m_auxPatches[i].getPatchRotated());

        temp_mp.computeTopology();

        // Right patch along the interface. Patch 0 -> v coordinate. Edge west along interface
        switch (temp_mp.interfaces()[0].second().side().index())
        {
            case 1:
                //gsInfo << "Global patch: " << patch_2 << "\tLocal patch: " << temp_mp.interfaces()[0].second().patch << " not rotated\n";
                break;
            case 4: m_auxPatches[0].rotateParamClock();
                //gsInfo << "Global patch: " << patch_2 <<"\tLocal patch: " << temp_mp.interfaces()[0].second().patch << " rotated clockwise\n";
                break;
            case 3: m_auxPatches[0].rotateParamAntiClock();
                //gsInfo << "Global patch: " << patch_2 <<"\tLocal patch: " << temp_mp.interfaces()[0].second().patch << " rotated anticlockwise\n";
                break;
            case 2: m_auxPatches[0].rotateParamAntiClockTwice();
                //gsInfo << "Global patch: " << patch_2 <<"\tLocal patch: " << temp_mp.interfaces()[0].second().patch << " rotated twice anticlockwise\n";
                break;
            default:
                break;
        }

        // Left patch along the interface. Patch 1 -> u coordinate. Edge south along interface
        switch (temp_mp.interfaces()[0].first().side().index())
        {
            case 3:
                //gsInfo << "Global patch: " << patch_1 <<"\tLocal patch: " << temp_mp.interfaces()[0].first().patch << " not rotated\n";
                break;
            case 4: m_auxPatches[1].rotateParamAntiClockTwice();
                //gsInfo << "Global patch: " << patch_1 <<"\tLocal patch: " << temp_mp.interfaces()[0].first().patch << " rotated twice anticlockwise\n";
                break;
            case 2: m_auxPatches[1].rotateParamAntiClock();
                //gsInfo << "Global patch: " << patch_1 <<"\tLocal patch: " << temp_mp.interfaces()[0].first().patch << " rotated anticlockwise\n";
                break;
            case 1: m_auxPatches[1].rotateParamClock();
                //gsInfo << "Global patch: " << patch_1 <<"\tLocal patch: " << temp_mp.interfaces()[0].first().patch << " rotated clockwise\n";
                break;
            default:
                break;
        }
    } // reparametrizeInterfacePatches


    template<short_t d,class T>
    void gsApproxC1Edge<d,T>::reparametrizeSinglePatch(index_t side)
    {
        computeAuxTopology();

        if(m_auxPatches[0].getOrient())
        {
            switch (side)
            {
                case 3:
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " not rotated\n";
                    break;
                case 2:
                    m_auxPatches[0].rotateParamClock();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated clockwise\n";
                    break;
                case 4:
                    m_auxPatches[0].rotateParamAntiClockTwice();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated twice anticlockwise\n";
                    break;
                case 1:
                    m_auxPatches[0].rotateParamAntiClock();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated anticlockwise\n";
                    break;
            }
        }
        else
        {
            switch (side)
            {
                case 1:
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " not rotated\n";
                    break;
                case 4:
                    m_auxPatches[0].rotateParamClock();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated clockwise\n";
                    break;
                case 2:
                    m_auxPatches[0].rotateParamAntiClockTwice();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated twice anticlockwise\n";
                    break;
                case 3:
                    m_auxPatches[0].rotateParamAntiClock();
                    //gsInfo << "Global patch: " << patch_1 << " with side " << side << " rotated anticlockwise\n";
                    break;
            }
        }
    } // reparametrizeSinglePatch

} // namespace gismo