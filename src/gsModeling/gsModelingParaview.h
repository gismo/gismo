/** @file gsModelingParaview.h

    @brief Paraview export of gsModeling types (moved from gsIO/gsWriteParaview,
    modularization stream S3 step A8: type-specific visualization lives
    with the type's module, the base IO module stays type-blind).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsIO/gsWriteParaview.h>
#include <gsModeling/gsSolid.h>
#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsPlanarDomain.h>
#include <gsModeling/gsCurveLoop.h>

#define NS 1000 // default sampling, mirrors gsWriteParaview.h

namespace gismo
{

/// \brief Depicting edge graph of each volume of one gsSolid with a segmenting loop
///
/// \param sl a gsMesh object
/// \param fn filename where paraview file is written
/// \param numPoints_for_eachCurve number of points used for sampling each curve
/// \param vol_Num ID of face(s), that should be written
/// \param edgeThick thickness of edges
/// \param translate "translate" vector, toward the volume is translated
/// \param color_convex Color, if face is convex and not eloop.
/// \param color_nonconvex Color, if face is not convex
/// \param color_eloop Color, if is in heSet and convex
/// \param eloop     a vector of ID numbers of vertices, often for representing a segmenting loop
/// \todo please document
template <class T>
void gsWriteParaview(gsSolid<T> const& sl, std::string const & fn,
                     unsigned numPoints_for_eachCurve=50, int vol_Num=0,
                     T edgeThick=0.01, gsVector3d<T> const & translate=gsVector3d<T>(0,0,0),
                     int color_convex=0, int color_nonconvex=20, int color_eloop=10,
                     std::vector<unsigned> const & eloop=std::vector<unsigned>());

/// Export a gsSolid to Paraview file
template <class T>
void gsWriteParaviewSolid(gsSolid<T> const& sl,
                          std::string const & fn,
                          unsigned numSamples = NS);

/// \brief Visualizing a gsCurveLoop
///
/// \param cloop the curve loop
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each curve
template<class T>
void gsWriteParaview(gsCurveLoop<T> const & cloop, std::string const & fn, unsigned npts)
{
    std::vector<gsGeometry<T> *> all_curves;
    for(index_t j =0; j< cloop.numCurves() ; j++)
        all_curves.push_back( const_cast<gsCurve<T> *>(cloop.curve(j)) );

    gsWriteParaview( all_curves, fn, npts);
}

/// \brief Visualizing a gsPlanarDomain
///
/// \param pdomain the planar domain
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling
template<class T>
void gsWriteParaview(gsPlanarDomain<T> const & pdomain,
                     std::string const & fn, unsigned npts=NS);

/// Visualizing a gsTrimSurface
template<class T>
void gsWriteParaview(const gsTrimSurface<T> & ts, std::string const & fn,
                     unsigned npts=NS, bool trimCurves = false);

/// \brief Export a volumeBlock.
///
/// Currently: output file shows boundary curves of this block.
///
/// \param volBlock pointer to the volume block
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling of a curve
template<typename T>
void gsWriteParaview(const gsVolumeBlock<T>& volBlock,
                     std::string const & fn,
                     unsigned npts = NS);

/// \brief Export a boundary/hole curve in trimmed surface
///
/// \param surf trimmed surface
/// \param idLoop curve loop number of a curve (0 - boundary, > 0 - hole)
/// \param idCurve curve number in a curve loop
/// \param fn filename (output paraview file)
/// \param npts number of points used for sampling a curve
template<typename T>
void gsWriteParaviewTrimmedCurve(const gsTrimSurface<T>& surf,
                                 const unsigned idLoop,
                                 const unsigned idCurve,
                                 const std::string fn,
                                 unsigned npts = NS);

} // namespace gismo

#undef NS

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsModelingParaview.hpp)
#endif
