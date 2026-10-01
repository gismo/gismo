/** @file gsParaviewCollection.hpp

    @brief Provides a helper class to write Paraview collection (.pvd) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsParaviewCollection.h>
#include <gsIO/gsParaviewDataSet.h>

namespace gismo
{

template <class T>
gsParaviewCollection<T>::gsParaviewCollection(const String& fn)
    : m_filename(fn)
    , m_isSaved(false)
    , m_time(-1)
    , m_evaluator(nullptr)
    , m_options(gsParaviewDataSet<T>::defaultOptions())
    , counter(0)
{
    std::string path = gsFileManager::getPath(m_filename);

    if (!gsFileManager::isFullyQualified(path))
        path = gsFileManager::getCurrentPath() + path;

    m_filename = path + gsFileManager::getBasename(m_filename) + ".pvd";
    gsFileManager::mkdir(path);

    mfile << "<?xml version=\"1.0\"?>\n";
    mfile << "<VTKFile type=\"Collection\" version=\"0.1\">\n";
    mfile << "<Collection>\n";
}

template <class T>
gsParaviewCollection<T>::gsParaviewCollection(const String& fn, gsExprEvaluator<T>& evaluator)
    : gsParaviewCollection<T>(fn)
{
    m_evaluator = &evaluator;
}

template <class T>
gsParaviewCollection<T>::~gsParaviewCollection() = default;

template <class T>
void gsParaviewCollection<T>::saveTimeStep()
{
    GISMO_ENSURE(m_dataset && !m_dataset->isEmpty(),
        "The gsParaviewDataSet, stored internally by gsParaviewCollection, is empty! "
        "Try running newTimestep() before saveTimeStep().");
    addDataSet(*m_dataset, m_time);
}

template <class T>
void gsParaviewCollection<T>::addDataSet(gsParaviewDataSet<T>& dataSet, T time)
{
    GISMO_ENSURE(!dataSet.isEmpty(), "The gsParaviewDataSet you are trying to add is empty!");
    GISMO_ASSERT(time >= 0, "Time should be a non-negative real number.");

    if (!dataSet.isSaved())
        dataSet.save();
    std::vector<std::string> filenames(dataSet.filenames());

    time = (time == static_cast<T>(-1)) ? m_time : time;
    mfile << "<!-- Time = " << time << " -->\n";
    for (size_t i = 0; i != filenames.size(); ++i)
    {
        std::string name;
        if (filenames[i].find("_mesh") != std::string::npos)
            name = "Mesh" + std::to_string(i);
        else if (filenames[i].find("_cnet") != std::string::npos)
            name = "Control Net" + std::to_string(i);
        else if (filenames[i].find("_patch") != std::string::npos)
            name = "Geometry" + std::to_string(i);
        else
            name = "";

        std::string relativeFilename = gsFileManager::makeRelative(
            gsFileManager::getPath(m_filename), filenames[i]);
        addPart(relativeFilename, time, name);
    }
}

template <class T>
void gsParaviewCollection<T>::newTimeStep(const gsMultiPatch<T>& geometry, T time)
{
    GISMO_ASSERT(!m_dataset || m_dataset->isEmpty() || m_dataset->isSaved(),
                 "Previous timestep has not been saved. try running saveTimeStep() before newTimeStep().");
    GISMO_ASSERT(time == static_cast<T>(-1) || time >= 0, "Time should be a non-negative real number.");

    if (time == static_cast<T>(-1))
    {
        m_time += 1.0;
        time = m_time;
    }
    else
    {
        m_time = time;
    }

    std::string name;
    if (m_options.askSwitch("makeSubfolder", true))
    {
        std::string subfolder = m_options.askString("subfolder");
        subfolder = ("" == subfolder) ? gsFileManager::getBasename(m_filename) + "_pvd" : subfolder;
        gsFileManager::mkdir(gsFileManager::getPath(m_filename) + subfolder);
        char sep = gsFileManager::getNativePathSeparator();
        name = gsFileManager::getPath(m_filename) + subfolder + sep + gsFileManager::getBasename(m_filename);
    }
    else
    {
        name = gsFileManager::getPath(m_filename) + gsFileManager::getBasename(m_filename);
    }

    name += "_t" + std::to_string(cast<T, double>(time));

    if (m_evaluator)
        m_dataset = memory::make_unique(new gsParaviewDataSet<T>(name, geometry, *m_evaluator, m_options));
    else
        m_dataset = memory::make_unique(new gsParaviewDataSet<T>(name, geometry, m_options));
}

} // namespace gismo
