/** @file gsParaviewDataSet.hpp

    @brief Provides a helper class to write Paraview (.vts) files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Karampatzakis, A. Mantzaflaris
*/

#pragma once

#include <gsIO/gsParaviewDataSet.h>

#include <gsIO/gsWriteParaview.h>
#include <gsMSplines/gsMappedBasis.h>
#include <gsCore/gsDofMapper.h>
#include <gsAssembler/gsExprEvaluator.h>

namespace gismo
{

template <class T>
unsigned gsParaviewDataSet<T>::getNumPoints(const gsOptionList& opts)
{
    const index_t canonical = opts.askInt("numPoints", 1000);
    const index_t legacy = opts.askInt("plot.npts", 1000);
    if (canonical != 1000)
        return static_cast<unsigned>(canonical);
    if (legacy != 1000)
        return static_cast<unsigned>(legacy);
    return 1000u;
}

template <class T>
bool gsParaviewDataSet<T>::getPlotElements(const gsOptionList& opts)
{
    const bool canonical = opts.askSwitch("plotElements", false);
    const bool legacy = opts.askSwitch("plot.elements", false);
    if (canonical != false)
        return canonical;
    if (legacy != false)
        return legacy;
    return false;
}

template <class T>
gsParaviewDataSet<T>::gsParaviewDataSet(std::string basename,
                                        const gsMultiPatch<T>& geometry,
                                        gsOptionList options)
    : m_basename(basename)
    , m_geometry(&geometry)
    , m_evaltr(nullptr)
    , m_options(options)
    , m_isSaved(false)
{
    const unsigned nPts = getNumPoints(m_options);

    const bool export_base64 = m_options.askSwitch("base64", true);
    const bool is_little_endian = []() -> bool {
        const int n{1};
        return *(char*)&n == 1;
    }();

    initFilenames();
    for (index_t k = 0; k != m_geometry->nPieces(); ++k)
    {
        gsMatrix<T> activeBases = m_geometry->piece(k).support();
        gsGridIterator<T, CUBE> pt(activeBases, nPts);

        const gsVector<index_t>& np(pt.numPointsCwise());
        index_t np1 = (np.size() > 1 ? np(1) - 1 : 0);
        index_t np2 = (np.size() > 2 ? np(2) - 1 : 0);

        std::ofstream file(m_filenames[k].c_str());
        file << std::fixed;
        file << std::setprecision(m_options.askInt("precision", 5));
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\"";
        if (export_base64)
        {
            file << " byte_order=\""
                 << (is_little_endian ? "LittleEndian" : "BigEndian")
                 << "\" header_type=\"UInt64\"";
        }
        file << ">\n";
        file << "<StructuredGrid WholeExtent=\"0 " << np(0) - 1 << " 0 "
             << np1 << " 0 " << np2 << "\">\n";
        file << "<Piece Extent=\"0 " << np(0) - 1 << " 0 " << np1 << " 0 "
             << np2 << "\">\n";
        file << "<PointData>\n";
        file.close();
    }
}

template <class T>
gsParaviewDataSet<T>::gsParaviewDataSet(std::string basename,
                                        const gsMultiPatch<T>& geometry,
                                        gsExprEvaluator<T>& eval,
                                        gsOptionList options)
    : gsParaviewDataSet<T>(basename, geometry, options)
{
    m_evaltr = &eval;
}

template <class T>
const std::vector<std::string> gsParaviewDataSet<T>::filenames()
{
    return m_filenames;
}

template <class T>
void gsParaviewDataSet<T>::save()
{
    GISMO_ASSERT(!m_isSaved, "gsParaviewDataSet already saved.");
    if (!m_isSaved)
    {
        m_isSaved = true;
        const gsMultiPatch<T>& geometry = *m_geometry;

        const unsigned nPts = getNumPoints(m_options);
        const unsigned precision = static_cast<unsigned>(m_options.askInt("precision", 5));
        const bool plotElements = getPlotElements(m_options);
        const bool plotControlNet = m_options.askSwitch("plotControlNet", false);
        const bool export_base64 = m_options.askSwitch("base64", false);

        const std::vector<std::string> points =
            toParaview(geometry, nPts, precision, "", export_base64);

        for (index_t k = 0; k != geometry.nPieces(); ++k)
        {
            std::ofstream file;
            file.open(m_filenames[k].c_str(), std::ios_base::app);
            file << "</PointData>\n\n\n<!-- GEOMETRY -->\n<Points>\n";
            file << points[k];
            file << "</Points>\n</Piece>\n</StructuredGrid>\n</VTKFile>";
            file.close();
            if (plotControlNet)
            {
                writeSingleControlNet(geometry.piece(k), m_basename + "_cnet" + std::to_string(k));
                m_filenames.push_back(m_basename + "_cnet" + std::to_string(k) + ".vtp");
            }
            if (plotElements)
            {
                int numPoints = m_options.getInt("plotElements.resolution");
                if (-1 == numPoints)
                {
                    const T evalPtsPerElem = 16 * (1.0 / geometry.piece(k).basis().numElements());
                    numPoints = cast<T, int>(
                        static_cast<T>(math::max(geometry.piece(k).basis().maxDegree() - 1, (short_t)1))
                        * math::pow(evalPtsPerElem, static_cast<T>(1.0) / static_cast<T>(geometry.domainDim())));
                }

                if (m_evaltr)
                {
                    gsMesh<T> msh(*m_evaltr->exprData()->domain().subdomain(k), numPoints);
                    static_cast<const gsGeometry<T>&>(geometry.piece(k)).evaluateMesh(msh);
                    gsWriteParaview(msh, m_basename + "_mesh" + std::to_string(k), false);
                }
                else
                {
                    gsMesh<T> msh(geometry.piece(k).basis(), numPoints);
                    static_cast<const gsGeometry<T>&>(geometry.piece(k)).evaluateMesh(msh);
                    gsWriteParaview(msh, m_basename + "_mesh" + std::to_string(k), false);
                }
                m_filenames.push_back(m_basename + "_mesh" + std::to_string(k) + ".vtp");
            }
        }
    }
}

template <class T>
bool gsParaviewDataSet<T>::isEmpty()
{
    return false;
}

template <class T>
bool gsParaviewDataSet<T>::isSaved()
{
    return m_isSaved;
}

template <class T>
void gsParaviewDataSet<T>::initFilenames()
{
    std::vector<std::string> names;
    for (index_t k = 0; k != m_geometry->nPieces(); ++k)
    {
        names.push_back(m_basename + "_patch" + std::to_string(k) + ".vts");
    }
    m_filenames = names;
}

} // End namespace gismo
