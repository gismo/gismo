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

        time = time==-1 ? m_time : time;
        mfile << "<!-- Time = " << time << " -->\n"; 
        for (size_t i=0; i!=filenames.size(); i++)
        {
            // Keep only the filename because <DataSet> tags in the .pvd
            // file use relative paths
            std::string name;
            if (filenames[i].find("_mesh") != std::string::npos)
                name="Mesh";
            else if (filenames[i].find("_cnet") != std::string::npos)
                name="Control Net";
            else if (filenames[i].find("_patch") != std::string::npos)
                name="Geometry";
            else
                name="";
            
            // This is so only the relative part of the path in filenames[] is kept
            addPart( gsFileManager::getBasename(filenames[i]) +"." +
                     gsFileManager::getExtension(filenames[i]), -1, time, name); 
        }
    }

    void gsParaviewCollection::newTimeStep(gsMultiPatch<real_t> * geometry, real_t time)
    {   
        if (-1 == time )
        {
            m_time += 1.0;
            time = m_time;
        }
        else { m_time = time; }
        std::string name = gsFileManager::getPath(m_filename) + gsFileManager::getBasename(m_filename);
        name += "_t" + std::to_string(time);
       
        m_dataset = gsParaviewDataSet(name, geometry, m_evaluator, m_options);
    }
}

