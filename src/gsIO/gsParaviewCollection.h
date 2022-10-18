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
// #include <gsAssembler/gsExprHelper.h>
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
class gsParaviewCollection
{
public:
    typedef std::string String;
public:

    /// Constructor using a filename ( without extension ).
    gsParaviewCollection(String const  &fn)
    : mfn(fn), counter(0), m_evaluator(nullptr)
    {
        mfile <<"<?xml version=\"1.0\"?>\n";
        mfile <<"<VTKFile type=\"Collection\" version=\"0.1\">";
        mfile <<"<Collection>\n";
    }

    /// Constructor using a filename ( without extension ) and an evaluator.
    gsParaviewCollection(String const  &fn, gsExprEvaluator<> * evaluator)
    : mfn(fn), counter(0), m_evaluator(evaluator)
    {
        mfile <<"<?xml version=\"1.0\"?>\n";
        mfile <<"<VTKFile type=\"Collection\" version=\"0.1\">";
        mfile <<"<Collection>\n";
    }
    
    /// Adds a part in the collection, with complete filename (including extension) \a fn
    void addPart(String const & fn)
    {
        GISMO_ASSERT(fn.find_last_of(".") != String::npos, "File without extension");
        GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        mfile << "<DataSet part=\""<<counter++<<"\" file=\""<<fn<<"\"/>\n";
    }

    /// Adds a part in the collection, with filename \a fn with extension \a ext appended
    void addPart(String const & fn, String const & ext)
    {
        GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        mfile << "<DataSet part=\""<<counter++<<"\" file=\""<<fn<<ext<<"\"/>\n";
    }

    /// Adds a part in the collection, with filename \a fni and extension \a ext appended
    void addPart(String const & fn, int i, String const & ext)
    {
        GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        mfile << "<DataSet part=\""<<i<<"\" file=\""<<fn<<i<<ext<<"\"/>\n";
    }
    

    // Full filename ( with extension )
    void addPart(String const & fn, index_t part=-1, real_t tStep=-1)
    {   
        // GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        mfile << "<DataSet ";
        if (part != -1)   mfile << "part=\""<< part <<"\" ";
        if (tStep != -1)  mfile << "timestep=\""<< tStep <<"\" ";
        mfile << "file=\"" << fn <<"\"/>\n";
    }

    // to do: make time collections as well
	// ! i is not included in the filename, must be in included fn !
    void addTimestep(String const & fn, int tstep, String const & ext)
    {
        mfile << "<DataSet timestep=\""<<tstep<<"\" file=\""<<fn<<ext<<"\"/>\n";
    }

    // EVERY PATCH NEEDS TO BE PUT INTO ITS OWN "PART" THUS ITS OWN <DATASET>

    // The part does not need to be specified as long as the <DataSet> appear
    // in the same order for each timestep

    // A data set is meant to be an abstraction for multiple <DataSet> tags in paraview, 
    // that all stem from the same gsGeometryMap, and refer to the same timestep.
    void addDataSet(gsParaviewDataSet dataSet, real_t time)
    {
        dataSet.save(); // the actual files are written to disk/finalized
        std::vector<std::string> filenames( dataSet.filenames() );

        for (size_t i=0; i!=filenames.size(); i++)
        {
            addPart( mfn + std::to_string(i), time);
        }
    }

    void /*gsParaviewDataSet*/ newTimeStep(real_t time=-1)
    {
        // TODO:
        // This will return an empty dataSet that has a proper filename according
        // to the internal timestep numbering

        // The user will then add all the desired fields to it, and then execute
        // addDataSet() to append it to the pvd file.
        // return gsParaviewDataSet("", ,m_evaluator); 
    }

    /// Finalizes the collection by closing the XML tags, always call
    /// this function (once) when you finish adding files
    void save()
    {
        // GISMO_ASSERT(counter!=-1, "Error: gsParaviewCollection::save() already called." );
        // mfile <<"</Collection>\n";
        // mfile <<"</VTKFile>\n";

        // mfn.append(".pvd");
        // std::ofstream f( mfn.c_str() );
        // GISMO_ASSERT(f.is_open(), "Error creating "<< mfn );
        // f << mfile.rdbuf();
        // f.close();
        // mfile.str("");
        // counter = -1;
    }

private:
    /// Pointer to char stream
    std::stringstream mfile;

    /// File name
    String mfn;

    /// Counter for the number of parts (files) added in the collection
    int counter;

    int m_step_count=-1;

    gsExprEvaluator<> * m_evaluator;

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
            pc.addPart(fn, i, ext);
    }
    else
        pc.addPart(fn, ext);

    pc.save();
}


} // end namespace gismo
