/** @file gsParaviewCollection.cpp

    @brief Provides a helper class to write Paraview collection (.pvd) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gsIO/gsParaviewCollection.h>


namespace gismo
{
    // EVERY PATCH NEEDS TO BE PUT INTO ITS OWN "PART" THUS ITS OWN <DATASET>
    // The part does not need to be specified as long as the <DataSet> appear
    // in the same order for each timestep

    // A gsParaviewDataSet is meant to be an abstraction for multiple <DataSet> tags in paraview, 
    // that all stem from the same gsGeometryMap, and refer to the same timestep.
    void gsParaviewCollection::addDataSet(gsParaviewDataSet dataSet, real_t time)
    {
        dataSet.save(); // the actual files are written to disk/finalized
        std::vector<std::string> filenames( dataSet.filenames() );

        time = time==-1 ? m_step_count : time;
        mfile << "<!-- Time = " << time << " -->\n"; 
        for (size_t i=0; i!=filenames.size(); i++)
        {
            // Keep only the filename because <DataSet> tags in the .pvd
            // file use relative paths
            addPart( gsFileManager::getBasename(filenames[i]) +"." +
                     gsFileManager::getExtension(filenames[i]), -1, time); 
        }
    }

    void gsParaviewCollection::newTimeStep(gsExprHelper<real_t>::geometryMap * geoMap, real_t time)
    {   
        if (-1 == time )
        {
            m_step_count += 1.0;
            time = m_step_count;
        }
        else { m_step_count = time; }
        std::string name = gsFileManager::getPath(mfn) + gsFileManager::getBasename(mfn);
        name += "_t" + std::to_string(time);
       
        m_dataset = new gsParaviewDataSet(name, geoMap, m_evaluator);
    }
}

