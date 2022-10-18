/** @file gsParaviewDataSet<real_t>.hpp

    @brief Provides a helper class to write Paraview (.vts) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Karampatzakis, A. Mantzaflaris
*/

#include<gsIO/gsParaviewDataSet.h>



namespace gismo
{
    gsParaviewDataSet::gsParaviewDataSet(std::string basename,
                    gsExprHelper<real_t>::geometryMap * geoMap,
                    gsExprEvaluator<real_t> * eval)
                    :m_basename(basename),
                    m_geoMap(geoMap),
                    m_evaltr(eval),
                    m_numPatches(geoMap->source().nPieces())
    {
        unsigned nPts = m_evaltr->options().askInt("plot.npts", 1000);

        // QUESTION: Can I be certain that the ids are consecutive?
        std::vector<std::string> fnames = filenames();
        for ( index_t k=0; k!=m_numPatches; k++) // For every patch.
        {
            gsMatrix<real_t> activeBases = m_geoMap->source().piece(k).support();
            gsGridIterator<real_t,CUBE> pt(activeBases, nPts);

            const gsVector<index_t> & np( pt.numPointsCwise() );
            index_t np1 = (np.size()>1 ? np(1)-1 : 0);
            index_t np2 = (np.size()>2 ? np(2)-1 : 0);

            // initializes individual .vts files
            // for every patch
            std::ofstream file(fnames[k].c_str());
            file << std::fixed; // no exponents
            file << std::setprecision(5); // PLOT_PRECISION
            file <<"<?xml version=\"1.0\"?>\n";
            file <<"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
            file <<"<StructuredGrid WholeExtent=\"0 "<< np(0)-1<<" 0 "<< np1 <<" 0 "
                << np2 <<"\">\n";
            file <<"<Piece Extent=\"0 "<< np(0)-1<<" 0 "<<np1<<" 0 "
                << np2 <<"\">\n";
            file <<"<PointData>\n";
            file.close();
        }
    }

    
    std::vector<std::string> gsParaviewDataSet::filenames()
    {   std::vector<std::string> names;
        for ( index_t k=0; k!=m_numPatches; k++) // For every patch.
        {
            names.push_back( m_basename + "_patch" +std::to_string(k)+".vts" );
        }
        return names;
    }


    void gsParaviewDataSet::save()
    {
        std::vector<std::string> points = m_evaltr->geoMap2vtk(*m_geoMap);
        // QUESTION: Can I be certain that the ids are consecutive?
        for ( index_t k=0; k!=m_numPatches; k++) // For every patch.
        {
            std::string filename;
            filename = m_basename + "_patch" +std::to_string(k)+".vts";
            std::ofstream file;
            file.open(filename.c_str(), std::ios_base::app); // Append to file 
            file << std::fixed; // no exponents
            file << std::setprecision(5); // PLOT_PRECISION


            file <<"</PointData>\n";
            file << points[k];
            file << "</Piece>\n</StructuredGrid>\n</VTKFile>";
            file.close();
        }
        // output text files for each part.
    }
} // End namespace gismo