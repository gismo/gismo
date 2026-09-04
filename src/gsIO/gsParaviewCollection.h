/** @file gsParaviewCollection.h

    @brief Provides a helper class to write Paraview collection (.pvd) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, C. Karampatzakis
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsIO/gsOptionList.h>
#include <gsIO/gsFileManager.h>
#include <gsIO/gsParaviewDataSet.h>

#include <fstream>
#include <utility>

namespace gismo
{

template <class T>
class GISMO_EXPORT gsParaviewCollection
{
public:
    typedef memory::unique_ptr<gsParaviewCollection<T> > uPtr;
    typedef std::string String;

public:
    /// Constructor using a filename.
    explicit gsParaviewCollection(const String& fn);

    /// Constructor using a filename and evaluator (for expression fields).
    gsParaviewCollection(const String& fn, gsExprEvaluator<T>& evaluator);

    ~gsParaviewCollection();

    void addPart(const String& fn, T tStep = -1, std::string name = "", index_t part = -1)
    {
        std::string ext = "";
        if (gsFileManager::getExtension(fn) == "")
        {
            if (name == "Mesh" || name == "mesh")
                ext = ".vtp";
            else if (name == "Geometry" || name == "geometry" || name == "Solution" || name == "solution")
                ext = ".vts";
            else
                GISMO_ERROR("No extension could be found for file " << fn << ". Try to add an extension or add name 'Mesh','Solution','Geometry'");
        }

        GISMO_ASSERT(!m_isSaved, "Error: collection has been already saved.");
        mfile << "<DataSet ";
        if (part != -1) mfile << "part=\"" << part << "\" ";
        if (tStep != static_cast<T>(-1)) mfile << "timestep=\"" << tStep << "\" ";
        if (name != "") mfile << "name=\"" << name << "\" ";
        mfile << "file=\"" << fn + ext << "\"/>\n";
    }

    GISMO_DEPRECATED void addTimestep(const String& fn, double tstep, const String& ext)
    {
        addPart(fn + ext, static_cast<T>(tstep));
    }

    GISMO_DEPRECATED void addTimestep(const String& fn, int part, double tstep, const String& ext)
    {
        addPart(fn + "_" + std::to_string(part) + ext, static_cast<T>(tstep), "", part);
    }

    GISMO_DEPRECATED void addPart(const String& fn, String extension)
    {
        addPart(fn + extension);
    }

    /// @brief Adds all files relevant to a gsParaviewDataSet to this collection.
    void addDataSet(gsParaviewDataSet<T>& dataSet, T time = -1);

    /// @brief Creates a new time step where all information will be added to.
    void newTimeStep(const gsMultiPatch<T>& geometry, T time = -1);

    /// @brief All arguments are forwarded to gsParaviewDataSet::addField().
    template <typename... Rest>
    void addField(Rest && ... rest)
    {
        GISMO_ENSURE(m_dataset && !m_dataset->isEmpty(),
            "The gsParaviewDataSet, stored internally by gsParaviewCollection, is empty! "
            "Try running newTimestep() before addField().");
        m_dataset->addField(std::forward<Rest>(rest)...);
    }

    /// @brief All arguments are forwarded to gsParaviewDataSet::addFields().
    template <typename... Rest>
    void addFields(Rest && ... rest)
    {
        GISMO_ENSURE(m_dataset && !m_dataset->isEmpty(),
            "The gsParaviewDataSet, stored internally by gsParaviewCollection, is empty! "
            "Try running newTimestep() before addFields().");
        m_dataset->addFields(std::forward<Rest>(rest)...);
    }

    /// @brief The current timestep is saved and files written to disk.
    void saveTimeStep();

    /// Finalizes the collection by closing the XML tags.
    void save()
    {
        GISMO_ASSERT(!m_isSaved, "Error: gsParaviewCollection::save() already called.");
        if (!m_isSaved)
        {
            mfile << "</Collection>\n";
            mfile << "</VTKFile>\n";

            gsDebug << "Exporting to " << m_filename << "\n";
            std::ofstream f(m_filename.c_str());
            GISMO_ASSERT(f.is_open(), "Error creating " << m_filename);
            f << mfile.rdbuf();
            f.close();
            mfile.str("");
            m_isSaved = true;
            counter = -1;
        }
    }

    gsOptionList& options() { return m_options; }

private:
    std::stringstream mfile;
    std::string m_filename;
    bool m_isSaved;
    real_t m_time;
    gsExprEvaluator<T>* m_evaluator;
    typename gsParaviewDataSet<T>::uPtr m_dataset;
    gsOptionList m_options;
    index_t counter;

private:
    gsParaviewCollection();
};

/// Fast creation of a collection using base filename \a fn and extension \a ext.
inline void makeCollection(std::string const& fn, std::string const& ext, int n = 0)
{
    gsParaviewCollection<real_t> pc(fn);
    const std::string base = gsFileManager::getFilename(fn);
    if (n > 0)
    {
        for (int i = 0; i < n; ++i)
            pc.addPart(base + std::to_string(i) + ext);
    }
    else
        pc.addPart(base + ext);

    pc.save();
}

} // end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsParaviewCollection.hpp)
#endif
