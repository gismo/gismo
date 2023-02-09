/** @file gsParaviewDataSet<real_t>.hpp

    @brief Provides a helper class to write Paraview (.vts) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Karampatzakis, A. Mantzaflaris
*/

#include<gsIO/gsParaviewDataSet.h>
#include<gsIO/gsWriteParaview.h>



namespace gismo
{
    gsParaviewDataSet::gsParaviewDataSet(std::string basename,
                    gsMultiPatch<real_t> * const geometry,
                    gsExprEvaluator<real_t> * eval,
                    gsOptionList options)
                    :m_basename(basename),
                    m_geometry(geometry),
                    m_evaltr(eval),
                    m_options(options)
    {
        unsigned nPts = m_options.askInt("numPoints",1000);

        // QUESTION: Can I be certain that the ids are consecutive?
        initFilenames();
        for ( index_t k=0; k!=m_geometry->nPieces(); k++) // For every patch.
        {
            gsMatrix<real_t> activeBases = m_geometry->piece(k).support();
            gsGridIterator<real_t,CUBE> pt(activeBases, nPts);

            const gsVector<index_t> & np( pt.numPointsCwise() );
            index_t np1 = (np.size()>1 ? np(1)-1 : 0);
            index_t np2 = (np.size()>2 ? np(2)-1 : 0);

            // initializes individual .vts files
            // for every patch
            std::ofstream file(m_filenames[k].c_str());
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

    
    const std::vector<std::string> gsParaviewDataSet::filenames()
    {
        return m_filenames;
    }


    void gsParaviewDataSet::save()
    {
        unsigned nPts = m_options.askInt("numPoints",1000);
        unsigned precision = m_options.askInt("precision",5);
        bool plotElements   = m_options.askSwitch("plotElements", false);
        bool plotControlNet = m_options.askSwitch("plotControlNet", false);

        std::vector<std::string> points = toVTK(*m_geometry,nPts,precision); //m_evaltr->geoMap2vtk(*m_geometry,nPts, precision);
        // QUESTION: Can I be certain that the ids are consecutive?
        for ( index_t k=0; k!=m_geometry->nPieces(); k++) // For every patch.
        {
            std::ofstream file;
            file.open(m_filenames[k].c_str(), std::ios_base::app); // Append to file 
            file <<"</PointData>\n\n\n<!-- GEOMETRY -->\n<Points>\n";
            file << points[k];
            file << "</Points>\n</Piece>\n</StructuredGrid>\n</VTKFile>";
            file.close();
            if (plotControlNet)
            {
                writeSingleControlNet( m_geometry->piece(k), m_basename + "_cnet" + std::to_string(k));
                m_filenames.push_back( m_basename + "_cnet" + std::to_string(k)+".vtp");
            } 
            if ( plotElements)
            {
                int numPoints = m_options.getInt("plotElements.resolution");
                if (-1 == numPoints )
                {
                    const real_t evalPtsPerElem = 16 * (1.0 / m_geometry->piece(k).basis().numElements());

                    // copied from gsWriteParaview
                    numPoints = cast<real_t,int>(
                        static_cast<real_t>(math::max( m_geometry->piece(k).basis().maxDegree()-1, (short_t)1))
                        * math::pow(evalPtsPerElem, (real_t)(1.0)/static_cast<real_t>(m_geometry->domainDim())) );
                }
                gsMesh<real_t> msh( gsMultiBasis<real_t>(*m_geometry).basis(k), numPoints);
                static_cast<const gsGeometry<real_t>&>(m_geometry->piece(k)).evaluateMesh(msh);
                gsWriteParaview(msh, m_basename + "_mesh" + std::to_string(k), false);
                m_filenames.push_back( m_basename + "_mesh" + std::to_string(k)+".vtp");
            }
        }
        // output text files for each part.
    }

    void gsParaviewDataSet::initFilenames()
    {   
        std::vector<std::string> names;
        for ( index_t k=0; k!=m_geometry->nPieces(); k++) // For every patch.
        {
            names.push_back( m_basename + "_patch" +std::to_string(k)+".vts" );
        }
        m_filenames = names;
    }
} // End namespace gismo
