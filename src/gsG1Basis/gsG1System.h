/** @file gsG1System.h

    @brief Create a G1-System for a Biharmonic equation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): P. Weinmueller
*/

#pragma once


namespace gismo
{
template<class T>
class gsG1System
{
public:

    gsG1System()
    {

    }

    void plotParaview(gsMultiPatch<T> & mp,
                      std::vector<gsG1AuxiliaryPatch> & interface,
                      std::vector<gsG1AuxiliaryPatch> & boundaries,
                      std::vector<gsG1AuxiliaryPatch> & vertices,
                      std::string basename);

}; // class gsG1System

template<class T>
void gsG1System<T>::plotParaview(gsMultiPatch<T> & mp,
                                 std::vector<gsG1AuxiliaryPatch> & interface,
                                 std::vector<gsG1AuxiliaryPatch> & boundaries,
                                 std::vector<gsG1AuxiliaryPatch> & vertices,
                                 std::string basename)
{
    gsParaviewCollection collection(basename);
    std::string fileName;
    index_t iter = 0;
    for (std::vector<gsG1AuxiliaryPatch>::iterator auxGeo = interface.begin(); auxGeo != interface.end(); ++auxGeo)
    {
        for (size_t i = 2; i < auxGeo->getG1Basis().nPatches()-2; i++)
        {
            fileName = basename + "_" + util::to_string(iter);
            gsField<> temp_field(mp.patch(auxGeo->getGlobalPatchIndex()),auxGeo->getG1Basis().patch(i));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addTimestep(fileName,iter,"0.vts");
            iter ++;
        }
    }
    for (std::vector<gsG1AuxiliaryPatch>::iterator auxGeo = boundaries.begin(); auxGeo != boundaries.end(); ++auxGeo)
    {
        for (size_t i = 0; i < auxGeo->getG1Basis().nPatches(); i++)
        {
            fileName = basename + "_" + util::to_string(iter);
            gsField<> temp_field(mp.patch(auxGeo->getGlobalPatchIndex()),auxGeo->getG1Basis().patch(i));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addTimestep(fileName,iter,"0.vts");
            iter ++;
        }
    }
    for (std::vector<gsG1AuxiliaryPatch>::iterator auxGeo = vertices.begin(); auxGeo != vertices.end(); ++auxGeo)
    {
        for (size_t i = 0; i < auxGeo->getG1Basis().nPatches(); i++)
        {
            fileName = basename + "_" + util::to_string(iter);
            gsField<> temp_field(mp.patch(auxGeo->getGlobalPatchIndex()),auxGeo->getG1Basis().patch(i));
            gsWriteParaview(temp_field,fileName,5000);
            collection.addTimestep(fileName,iter,"0.vts");
            iter ++;
        }
    }

    collection.save();
}

} // namespace
