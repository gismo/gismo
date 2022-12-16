/** @file gsParaviewCollection.h

    @brief Provides a helper class to write Paraview collection (.pvd) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsAssembler/gsExprHelper.h>
#include <gsIO/gsFileManager.h>
#include <gsIO/gsParaviewDataSet.h>



#include<fstream>

namespace gismo {
/**
    \brief This class is used to create a Paraview .pvd (collection)
    file.

    A collection is an XML file that contains a list of other
    files to be opened in Paraview.

    Typical usage is
    \verbatim
    gsParaviewCollection pc(fn); // Initialize collection

    // add files ("parts"), make sure that these files exist
    pc.addPart(filename1);
    pc.addPart(filename2);

    pc.save() // finalize and close the file
    \endverbatim

    The above creates a file with extension pvd. When opening this
    file with Paraview, the containts of all parts in the list are
    loaded.

    \ingroup IO
*/
class GISMO_EXPORT gsParaviewCollection
{
public:
    typedef std::string String;
public:

    /// Constructor using a filename and an (optional) evaluator.
    gsParaviewCollection(String const  &fn,
                         gsExprEvaluator<> * evaluator=nullptr)
                        : m_filename(fn),
                        counter(0),
                        m_step_count(-1),
                        m_evaluator(evaluator),
                        m_options(gsParaviewDataSet::defaultOptions())
    {
        m_filename = gsFileManager::getPath(m_filename) + gsFileManager::getBasename(m_filename) + ".pvd";
        gsInfo << m_filename << "\n";
        gsFileManager::mkdir( gsFileManager::getPath(m_filename) );
        // if ( "" != m_filename.parent_path())
        // GISMO_ENSURE( fsystem::exists( m_filename.parent_path() ), 
        //     "The specified folder " << m_filename.parent_path() << " does not exist, please create it first.");  
        mfile <<"<?xml version=\"1.0\"?>\n";
        mfile <<"<VTKFile type=\"Collection\" version=\"0.1\">\n";
        mfile <<"<Collection>\n";
    }
    
    // Full filename ( with extension )
    void addPart(String const & fn, index_t part=-1, real_t tStep=-1)
    {   
        GISMO_ASSERT( gsFileManager::getExtension(fn) != "" , "File without extension");
        GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        // GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        mfile << "<DataSet ";
        if (part != -1)   mfile << "part=\""<< part <<"\" ";
        if (tStep != -1)  mfile << "timestep=\""<< tStep <<"\" ";
        mfile << "file=\"" << fn <<"\"/>\n";
    }

    void addDataSet(gsParaviewDataSet dataSet, real_t time=-1);

    void newTimeStep(gsMultiPatch<real_t> * geometry, real_t time=-1);

 
    template <typename... Rest>
    void addField(Rest... rest)
    {
        m_dataset.addField(rest...);
    }

    template <typename... Rest>
    void addFields(Rest... rest)
    {
        m_dataset.addFields(rest...);
    }

    void saveTimeStep(){
        addDataSet(m_dataset,m_step_count);
    };

    /// Finalizes the collection by closing the XML tags, always call
    /// this function (once) when you finish adding files
    void save()
    {
        GISMO_ASSERT(counter!=-1, "Error: gsParaviewCollection::save() already called." );
        mfile <<"</Collection>\n";
        mfile <<"</VTKFile>\n";

        gsInfo << "Exporting to " << m_filename << "\n";
        std::ofstream f( m_filename.c_str() );
        GISMO_ASSERT(f.is_open(), "Error creating "<< m_filename );
        f << mfile.rdbuf();
        f.close();
        mfile.str("");
        counter = -1;
    }

    gsOptionList & options() {return m_options;}

private:
    /// Pointer to char stream
    std::stringstream mfile;

    /// File name
    std::string m_filename;

    /// Counter for the number of parts (files) added in the collection
    int counter;

    int m_step_count;

    gsExprEvaluator<> * m_evaluator;

    gsParaviewDataSet m_dataset;

    gsOptionList m_options;

private:
    // Construction without a filename is not allowed
    gsParaviewCollection();
};

//=================================================================================================



/// Fast creation of a collection using base filename \a fn, extension
/// \a ext.  The collection will contain the files fn_0.ext,
/// fn_1.ext,...,fn_{n-1}.ext In the special case of n=0, the
/// collection is just fn.pvd and contains only the file fn.ext
inline void makeCollection(std::string const & fn, std::string const & ext, int n = 0)
{
    gsParaviewCollection pc(fn);
    if ( n > 0)
    {
        for (int i=0; i<n ; i++)
            pc.addPart(fn + std::to_string(i) + ext);
    }
    else
        pc.addPart(fn + ext);

    pc.save();
}


} // end namespace gismo
