/** @file gsParaview.hpp

    @brief Provides implementation of gsParaview class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsIO/gsParaview.h>
#include <gsIO/gsWriteParaview.h>

#include <gsCore/gsGeometry.h>
#include <gsCore/gsGeometrySlice.h>
#include <gsCore/gsField.h>
#include <gsCore/gsDebug.h>

#include <gsCore/gsMultiPatch.h>

#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsSolid.h>

#include <gsMesh2/gsSurfMesh.h>

#include <gsHSplines/gsHBoxContainer.h>

namespace gismo
{

template<class T>
gsParaview<T>::gsParaview()
    : m_options(defaultOptions())
{ }

template<class T>
gsParaview<T>::gsParaview(const gsOptionList & options)
    : m_options(options)
{ }

template<class T>
gsParaview<T>::~gsParaview() = default;

template<class T>
gsOptionList gsParaview<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addInt   ("numPoints",      "Number of sampling points per patch", 1000);
    opt.addInt   ("precision",      "Output file floating-point precision", 5);
    opt.addString("patchDelimiter", "Delimiter between patch indices in filenames", "_");
    opt.addSwitch("plotElements",   "Plot the parameter mesh", false);
    opt.addSwitch("plotControlNet", "Plot the control net", false);
    opt.addSwitch("show",           "Open Paraview after writing", false);
    opt.addSwitch("bezier",         "Use Bezier elements export for MultiPatch", false);
    opt.addSwitch("boundary",       "Write boundaries for MultiPatch", false);
    opt.addSwitch("interfaces",     "Write interfaces for MultiPatch", false);
    opt.addSwitch("trimCurves",     "Plot trim curves (for gsTrimSurface)", false);
    opt.addSwitch("graph",          "Plot function as graph", true);
    opt.addSwitch("fullSupport",    "Plot basis over whole domain (for mapped basis)", false);
    opt.addInt   ("hboxMode",       "Mode for gsHBox/gsHBoxContainer: 0=level, 1=error, 2=projectedError", 0);
    opt.addSwitch("writePvd",       "Wrap one-shot writes in a .pvd collection file", true);
    opt.addSwitch("base64",         "Write unstructured-grid arrays in base64 binary format", false);
    return opt;
}

template<class T>
void gsParaview<T>::openIfRequested(const std::string & fn) const
{
    if (m_options.getSwitch("show"))
    {
        gsFileManager::open(fn + ".pvd");
    }
}

// ============================================================
// Write methods implementations
// ============================================================

template<class T>
void gsParaview<T>::write(const gsGeometry<T> & geo, const std::string & fn) const
{
    gsWriteParaview(geo, fn,
                    m_options.getInt("numPoints"),
                    m_options.getSwitch("plotElements"),
                    m_options.getSwitch("plotControlNet"),
                    !m_options.getSwitch("writePvd"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsMultiPatch<T> & mp, const std::string & fn) const
{
    const unsigned npts = m_options.getInt("numPoints");
    const bool mesh     = m_options.getSwitch("plotElements");
    const bool ctrlNet  = m_options.getSwitch("plotControlNet");
    const std::string pDelim = m_options.getString("patchDelimiter");
    const bool exportBase64 = m_options.getSwitch("base64");
    const bool skipPvd = !m_options.getSwitch("writePvd");

    if (m_options.getSwitch("bezier"))
    {
        gsWriteParaviewBezier(mp, fn, ctrlNet);
    }
    else
    {
        // NOTE: this gsMultiPatch convenience overload has no skipPvd parameter
        // (it forwards to the std::vector<gsGeometry*> overload, which doesn't
        // have one either); "writePvd" is not honored for this write() overload.
        gsWriteParaview(mp, fn, npts, mesh, ctrlNet, pDelim);
    }

    if (m_options.getSwitch("boundary"))
    {
        gsWriteParaviewBdr(mp, fn + "_boundary", npts, ctrlNet);
    }

    if (m_options.getSwitch("interfaces"))
    {
        gsWriteParaviewIfc(mp, fn + "_interfaces", npts, ctrlNet);
    }

    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(std::vector<gsGeometry<T>*> const & geos, const std::string & fn) const
{
    gsWriteParaview(geos, fn,
                    m_options.getInt("numPoints"),
                    m_options.getSwitch("plotElements"),
                    m_options.getSwitch("plotControlNet"),
                    m_options.getString("patchDelimiter"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsField<T> & field, const std::string & fn) const
{
    gsWriteParaview(field, fn,
                    m_options.getInt("numPoints"),
                    m_options.getSwitch("plotElements"),
                    m_options.getString("patchDelimiter"),
                    !m_options.getSwitch("writePvd"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsBasis<T> & basis, const std::string & fn) const
{
    gsWriteParaview(basis, fn,
                    m_options.getInt("numPoints"),
                    m_options.getSwitch("plotElements"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsBasis<T> & basis, const std::vector<index_t> & indices,
                          const std::string & fn) const
{
    gsWriteParaview(basis, indices, fn,
                    m_options.getInt("numPoints"),
                    m_options.getSwitch("plotElements"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsFunctionSet<T> & func, const std::string & fn) const
{
    gsWriteParaview(func, fn, m_options.getInt("numPoints"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsFunction<T> & func, const gsMatrix<T> & supp,
                          const std::string & fn) const
{
    gsWriteParaview(func, supp, fn,
                    m_options.getInt("numPoints"),
                    m_options.getSwitch("graph"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsFunctionSet<T> & geo, const gsFunctionSet<T> & func,
                          const std::string & fn) const
{
    gsWriteParaview(geo, func, fn,
                    m_options.getInt("numPoints"),
                    m_options.getString("patchDelimiter"),
                    !m_options.getSwitch("writePvd"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsGeometrySlice<T> & geo, const std::string & fn) const
{
    gsWriteParaview(geo, fn, m_options.getInt("numPoints"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsMesh<T> & mesh, const std::string & fn) const
{
    gsWriteParaview(mesh, fn, true);
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsSurfMesh & mesh, const std::string & fn,
                          std::initializer_list<std::string> props) const
{
    gsWriteParaview(mesh, fn, props);
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsMesh<T> & mesh, const gsMatrix<T> & params,
                          const std::string & fn) const
{
    gsWriteParaview(mesh, fn, params);
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const std::vector<gsMesh<T>> & meshes, const std::string & fn) const
{
    gsWriteParaview(meshes, fn);
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsTrimSurface<T> & ts, const std::string & fn) const
{
    gsWriteParaview(ts, fn,
                    m_options.getInt("numPoints"),
                    m_options.getSwitch("trimCurves"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsPlanarDomain<T> & pdomain, const std::string & fn) const
{
    gsWriteParaview(pdomain, fn, m_options.getInt("numPoints"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsCurveLoop<T> & cloop, const std::string & fn) const
{
    gsWriteParaview(cloop, fn, m_options.getInt("numPoints"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsVolumeBlock<T> & volBlock, const std::string & fn) const
{
    gsWriteParaview(volBlock, fn, m_options.getInt("numPoints"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsSolid<T> & solid, const std::string & fn) const
{
    gsWriteParaviewSolid(solid, fn, m_options.getInt("numPoints"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsMultiPatch<T> & mp, const gsMultiBasis<T> & mb,
                          const std::string & fn) const
{
    gsWriteParaview(mp, mb, fn, m_options.getInt("numPoints"),
                    !m_options.getSwitch("writePvd"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsMappedSpline<2,T> & mspline, const std::string & fn) const
{
    gsWriteParaview(mspline, fn, m_options.getInt("numPoints"),
                    !m_options.getSwitch("writePvd"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsFunctionSet<T> & geom, const gsMappedBasis<2,T> & mbasis,
                          const std::string & fn,
                          const std::vector<index_t> & indices) const
{
    gsWriteParaview(geom, mbasis, fn,
                    m_options.getInt("numPoints"),
                    m_options.getSwitch("fullSupport"),
                    indices,
                    !m_options.getSwitch("writePvd"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsMultiBasis<T> & mb, const gsMultiPatch<T> & domain,
                          const std::string & fn) const
{
    gsWriteParaview(mb, domain, fn, m_options.getInt("numPoints"));
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsMultiPatch<T> & patches,
                          const typename gsBoundaryConditions<T>::bcContainer & bcs,
                          const std::string & fn) const
{
    gsWriteParaview(patches, bcs, fn,
                    m_options.getInt("numPoints"),
                    m_options.getSwitch("plotControlNet"));
    openIfRequested(fn);
}

// ============================================================
// Write methods for boxes
// ============================================================

template<class T>
void gsParaview<T>::write(const gsMatrix<T> & boxes, const std::string & fn,
                          const std::vector<T> & values) const
{
    gsWriteParaview(boxes, fn, values);
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::writeBoxes(const gsMatrix<T> & boxes, const std::string & fn,
                               const std::vector<T> & values) const
{
    gsWriteParaview(boxes, fn, values);
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::write(const gsMatrix<T> & boxes, const gsVector<T> & values,
                          const std::string & fn) const
{
    gsWriteParaview(boxes, fn, values);
    openIfRequested(fn);
}

template<class T>
template<short_t d>
void gsParaview<T>::write(const gsHBox<d,T> & box, const std::string & fn) const
{
    gsWriteParaview(box, fn, static_cast<short_t>(m_options.getInt("hboxMode")));
    openIfRequested(fn);
}

template<class T>
template<short_t d>
void gsParaview<T>::write(const gsHBoxContainer<d,T> & boxes, const std::string & fn) const
{
    gsWriteParaview(boxes, fn, static_cast<short_t>(m_options.getInt("hboxMode")));
    openIfRequested(fn);
}

// ============================================================
// Write methods for points
// ============================================================

template<class T>
void gsParaview<T>::writePoints(const gsMatrix<T> & points, const std::string & fn) const
{
    gsWriteParaviewPoints(points, fn);
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::writePoints(const gsMatrix<T> & points, const std::string & fn,
                                const gsVector<T> & values) const
{
    gsWriteParaview(points, fn, values);
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::writePoints(const gsMatrix<T> & X, const gsMatrix<T> & Y,
                                const std::string & fn) const
{
    gsWriteParaviewPoints(X, Y, fn);
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::writePoints(const gsMatrix<T> & X, const gsMatrix<T> & Y,
                                const gsMatrix<T> & Z, const std::string & fn) const
{
    gsWriteParaviewPoints(X, Y, Z, fn);
    openIfRequested(fn);
}

template<class T>
void gsParaview<T>::writePoints(const gsMatrix<T> & X, const gsMatrix<T> & Y,
                                const gsMatrix<T> & Z, const gsMatrix<T> & V,
                                const std::string & fn) const
{
    gsWriteParaviewPoints(X, Y, Z, V, fn);
    openIfRequested(fn);
}

// ============================================================
// Write methods for single basis function
// ============================================================

template<class T>
void gsParaview<T>::writeBasisFunction(int i, const gsBasis<T> & basis,
                                       const std::string & fn) const
{
    gsWriteParaview_basisFnct(i, basis, fn, m_options.getInt("numPoints"));
    openIfRequested(fn);
}

// ============================================================
// Write methods for trimmed curves
// ============================================================

template<class T>
void gsParaview<T>::writeTrimmedCurve(const gsTrimSurface<T> & surf,
                                      unsigned idLoop, unsigned idCurve,
                                      const std::string & fn) const
{
    gsWriteParaviewTrimmedCurve(surf, idLoop, idCurve, fn,
                                m_options.getInt("numPoints"));
    openIfRequested(fn);
}

} // namespace gismo
