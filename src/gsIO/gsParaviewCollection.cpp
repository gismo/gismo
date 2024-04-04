/** @file gsParaviewCollection.cpp

    @brief Provides a helper class to write Paraview collection (.pvd) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, C. Karampatzakis
*/

#include <gsIO/gsParaviewCollection.h>


namespace gismo
{
    // EVERY PATCH NEEDS TO BE PUT INTO ITS OWN "PART" THUS ITS OWN <DATASET>
    // The part does not need to be specified as long as the <DataSet> appear
    // in the same order for each timestep

    // A gsParaviewDataSet is meant to be an abstraction for multiple <DataSet> tags in paraview, 
    // that all stem from the same gsGeometryMap, and refer to the same timestep.
    void gsParaviewCollection::addDataSet(gsParaviewDataSet & dataSet, real_t time)
    {
        GISMO_ENSURE(!dataSet.isEmpty(), "The gsParaviewDataSet you are trying to add is empty!");
        GISMO_ASSERT(time>=0, "Time should be a non-negative real number.");

        if (! dataSet.isSaved()) dataSet.save(); // the actual files are written to disk/finalized
        std::vector<std::string> filenames( dataSet.filenames() );

        time = time==-1 ? m_time : time;
        mfile << "<!-- Time = " << time << " -->\n"; 
        for (size_t i=0; i!=filenames.size(); i++)
        {
            // Keep only the filename because <DataSet> tags in the .pvd
            // file use relative paths
            std::string name;
            if (filenames[i].find("_mesh") != std::string::npos)
                name="Mesh"+std::to_string(i);
            else if (filenames[i].find("_cnet") != std::string::npos)
                name="Control Net"+std::to_string(i);
            else if (filenames[i].find("_patch") != std::string::npos)
                name="Geometry"+std::to_string(i);
            else
                name="";

            // This is so only the relative part of the path in filenames[] is kept
            std::string relativeFilename = gsFileManager::makeRelative(
                            gsFileManager::getCurrentPath() + gsFileManager::getPath(m_filename),
                            filenames[i]);
            
            addPart( relativeFilename, time, name); 
        }
    }

    void gsParaviewCollection::newTimeStep(gsMultiPatch<real_t> * geometry, real_t time)
    {   
        GISMO_ASSERT( m_dataset.isEmpty() || m_dataset.isSaved(), "Previous timestep has not been saved. try running saveTimeStep() before newTimeStep().");
        GISMO_ASSERT(-1==time || time>=0, "Time should be a non-negative real number.");

        if (-1 == time )
        {
            m_time += 1.0;
            time = m_time;
        }
        else { m_time = cast<real_t,int>( time ); }

        std::string name;
        if ( m_options.askSwitch("makeSubfolder",true) )
        {
            std::string subfolder = m_options.askString("subfolder");
            subfolder = ("" == subfolder ) ? gsFileManager::getBasename(m_filename)+"_pvd" : subfolder;
            gsFileManager::mkdir( gsFileManager::getPath(m_filename) + subfolder ); 
            char sep = gsFileManager::getNativePathSeparator();

            name = gsFileManager::getPath(m_filename) + subfolder + sep + gsFileManager::getBasename(m_filename);
        }
        else
        {
            name = gsFileManager::getPath(m_filename) + gsFileManager::getBasename(m_filename);
        }


        name += "_t" + std::to_string( cast<real_t, double>( time ) );
       
        m_dataset = gsParaviewDataSet(name, geometry, m_evaluator, m_options);
    }
}

