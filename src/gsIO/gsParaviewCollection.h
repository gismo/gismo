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

    Or when working with gsExpressions:
    \verbatim
    gsParaviewCollection pc(fn, &evaluator ); // Initialize collection

    // Specify plotting options
    pc.options().setInt("numPoints", 100);
    pc.options().setInt("precision",5);
    pc.options().setSwitch("plotElements",true);
    pc.options().setSwitch("plotControlNet",true);

    // In your solution loop:
    while ( ... )
    {
        // Create new file(s) for this timestep
        pc.newTimestep(&multiPatch)
        // Solve here 

        // Write solution fields ( e.g Pressure, Temperature, Stress )
        pc.addField( gsExpression_1, "label_1" );
        pc.addField( gsExpression_2, "label_2" );
        pc.saveTimeStep(); // Save file(s) for this timestep
    } // end solution loop

    pc.save() // finalize and close the .pvd file
    \endverbatim


    The above creates a file with extension pvd. When opening this
    file with Paraview, the contents of all parts in the list are
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
                        m_isSaved(false),
                        m_time(-1),
                        m_evaluator(evaluator),
                        m_options(gsParaviewDataSet::defaultOptions()),
                        counter(0)
    {
        m_filename = gsFileManager::getPath(m_filename) + gsFileManager::getBasename(m_filename) + ".pvd";
        gsFileManager::mkdir( gsFileManager::getPath(m_filename) );
        // if ( "" != m_filename.parent_path())
        // GISMO_ENSURE( fsystem::exists( m_filename.parent_path() ), 
        //     "The specified folder " << m_filename.parent_path() << " does not exist, please create it first.");  
        mfile <<"<?xml version=\"1.0\"?>\n";
        mfile <<"<VTKFile type=\"Collection\" version=\"0.1\">\n";
        mfile <<"<Collection>\n";
    }

    /// @brief Appends a file to the Paraview collection (.pvd file).
    /// @param fn Filename to be added. Can also be a path relative to the where the collection file is. 
    /// @param tStep Time step ( optional )
    /// @param name An optional name for this part. This is the name that will show up  in Paraview's MultiBlock inspector for this Part. It is not the filename.
    /// @param part Part ID ( optional )
    void addPart(String const & fn, real_t tStep=-1, std::string name="", index_t part=-1)
    {   
        std::string ext = "";
        if (gsFileManager::getExtension(fn) == "")
        {
            if (name=="Mesh" || name=="mesh")
                ext = ".vtp";
            else if (name=="Geometry" || name=="geometry" || name=="Solution" || name=="solution")
                ext = ".vts";
            else
                GISMO_ERROR("No extension could be found for file "<<fn<<". Try to add an extension or add name 'Mesh','Solution','Geometry'");
        }

        GISMO_ASSERT( !m_isSaved , "Error: collection has been already saved." );
        mfile << "<DataSet ";
        if (part != -1)   mfile << "part=\""<< part <<"\" ";
        if (tStep != -1)  mfile << "timestep=\""<< tStep <<"\" ";
        if (name != "") mfile << "name=\"" << name << "\" ";
        mfile << "file=\"" << fn+ext <<"\"/>\n";
    }
    // CAUTION! 
    // The previous 3 versions of gsParaviewCollection::addPart() have been combined into the one above
    // since they were all doing basiacally the same thing. Below you can see a 'conversion table' that
    // can help you adapt your code to the new syntax. For questions contact C. Karampatzakis (Github @ckarampa )
    // OLD SYNTAX ( Deprecated )                              |   NEW SYNTAX 
    //-----------------------------------------------------------------------------------------------------------
    // addPart(String const & fn)                             |   addPart( fn, -1, "", counter++);               |
    //                                                                                                           |
    // addPart(String const & fn, String const & ext)         |   addPart( fn+ext, -1, "", counter++);           |
    //                                                                                                           |
    // addPart(String const & fn, int i, String const & ext)  |   addPart( fn+std::to_string(i)+ext, -1, "", i); |
    // ----------------------------------------------------------------------------------------------------------

    // The following functions are all deprecated, only here for backwards compatibility!

    GISMO_DEPRECATED void addTimestep(String const & fn, double tstep, String const & ext)
    {
        // mfile << "<DataSet timestep=\""<<tstep<<"\" file=\""<<fn<<ext<<"\"/>\n";
        addPart( fn+ext, tstep);
    }

    GISMO_DEPRECATED void addTimestep(String const & fn, int part, double tstep, String const & ext)
    {
        // mfile << "<DataSet part=\""<<part<<"\" timestep=\""<<tstep<<"\" file=\""<<fn<<"_"<<part<<ext<<"\"/>\n";
        addPart( fn+"_"+std::to_string(part)+ext, tstep, "", part);
    }

    // End of deprecated functions

    GISMO_DEPRECATED void addPart(String const & fn, String extension)
    {
        addPart( fn+extension);
    }

    /// @brief Adds all the files relevant to a gsParaviewDataSet, to the collection.
    /// @param dataSet The gsParaviewDataSet to be added.
    /// @param time Time step (optional, else an internal integer counter is used)
    void addDataSet(gsParaviewDataSet & dataSet, real_t time=-1);

    /// @brief Creates a new time step where all information will be added to.
    /// @param geometry A gsMultiPatch of the geometry where the solution fields are defined.
    /// @param time Value of time for this timestep (optional, else an internal integer counter is used)
    void newTimeStep(gsMultiPatch<real_t> * geometry, real_t time=-1);

 
    /// @brief All arguments are forwarded to gsParaviewDataSet::addField().
    template <typename... Rest>
    void addField(Rest... rest)
    {
        GISMO_ENSURE( !m_dataset.isEmpty(), "The gsParaviewDataSet, stored internally by gsParaviewCollection, is empty! Try running newTimestep() before addField().");
        m_dataset.addField(rest...);
    }

    /// @brief All arguments are forwarded to gsParaviewDataSet::addFields().
    template <typename... Rest>
    void addFields(Rest... rest)
    {
        GISMO_ENSURE( !m_dataset.isEmpty(), "The gsParaviewDataSet, stored internally by gsParaviewCollection, is empty! Try running newTimestep() before addFields().");
        m_dataset.addFields(rest...);
    }

    /// @brief The current timestep is saved and files written to disk.
    void saveTimeStep(){
        GISMO_ENSURE( !m_dataset.isEmpty(), "The gsParaviewDataSet, stored internally by gsParaviewCollection, is empty! Try running newTimestep() before saveTimeStep().");
        addDataSet(m_dataset,m_time);
    };

    /// Finalizes the collection by closing the XML tags, always call
    /// this function (once) when you finish adding files
    void save()
    {
        GISMO_ASSERT(!m_isSaved, "Error: gsParaviewCollection::save() already called." );
        if (!m_isSaved)
        {
            mfile <<"</Collection>\n";
            mfile <<"</VTKFile>\n";

            gsDebug << "Exporting to " << m_filename << "\n";
            std::ofstream f( m_filename.c_str() );
            GISMO_ASSERT(f.is_open(), "Error creating "<< m_filename );
            f << mfile.rdbuf();
            f.close();
            mfile.str("");
            m_isSaved=true;
            counter = -1;
        }
    }

    /// @brief Accessor to the current options.
    gsOptionList & options() {return m_options;}

private:
    /// Pointer to char stream
    std::stringstream mfile;

    /// File name
    std::string m_filename;

    /// Flag for checking if collection is already saved.
    bool m_isSaved;

    int m_time;

    gsExprEvaluator<> * m_evaluator;

    gsParaviewDataSet m_dataset;

    gsOptionList m_options;

    index_t counter;

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
