/** @file gsWriteParaview.h

    @brief Provides declaration of functions writing Paraview files.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsCore/gsGeometry.h>
#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsExport.h>

#include <sstream>
#include <fstream>

#pragma once

#define NS 1000


namespace gismo {


/// \brief Export a gsGeometry (without scalar information) to paraview file
///
/// \param Geo a geometry object
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each patch
/// \param mesh if true, the parameter mesh is plotted as well
/// \param ctrlNet if true, the control net is plotted as well
///
/// \ingroup IO
template<class T>
void gsWriteParaview(const gsGeometry<T> & Geo, std::string const & fn, 
                     unsigned npts=NS, bool mesh = false, bool ctrlNet = false);

/// \brief Export a mesh to paraview file
///
/// \param sl a gsMesh object
/// \param fn filename where paraview file is written
/// \param pvd if true, a .pvd file is generated (for compatibility)
template <class T>
void gsWriteParaview(gsMesh<T> const& sl, std::string const & fn, bool pvd = true);

/// \brief Export a vector of meshes, each mesh in its own file.
///
/// \param meshes vector of gsMesh objects
/// \param fn filename
template <typename T>
void gsWriteParaview(const std::vector<gsMesh<T> >& meshes, std::string const& fn);


/// \brief Write a file containing a solution field (as color on its geometry) to paraview file
///
/// \param field a field object
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each patch
/// \param mesh if true, the parameter mesh is plotted as well
template<class T>
void gsWriteParaview(const gsField<T> & field, std::string const & fn, 
                     unsigned npts=NS, bool mesh = false);

/// \brief Export a multipatch Geometry (without scalar information) to paraview file
///
/// \param Geo a multipatch object
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each patch
/// \param mesh if true, the parameter mesh is plotted as well
/// \param ctrlNet if true, the control net is plotted as well
template<class T>
void gsWriteParaview(const gsMultiPatch<T> & Geo, std::string const & fn, 
                     unsigned npts=NS, bool mesh = false, bool ctrlNet = false)
{
    gsWriteParaview( Geo.patches(), fn, npts, mesh, ctrlNet);
}

/// \brief Export a multipatch Geometry (without scalar information) to paraview file
///
/// \param Geo a vector of the geometries to be plotted
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each geometry
/// \param mesh if true, the parameter mesh is plotted as well
/// \param ctrlNet if true, the control net is plotted as well
template<class T>
void gsWriteParaview( std::vector<gsGeometry<T> *> const & Geo, 
                      std::string const & fn, unsigned npts=NS,
                      bool mesh = false, bool ctrlNet = false);

/// \brief Export a computational mesh to paraview file
template<class T>
void gsWriteParaview(const gsMultiBasis<T> & mb, const gsMultiPatch<T> & domain,
                     std::string const & fn, unsigned npts);

/// \brief Export a composite Geometry to paraview file
///
/// \param Geo a composite geometry
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each curve
/// \param mesh if true, the parameter mesh is plotted as well
/// \param ctrNet if true, the control net is plotted as well
template<class T>
void gsWriteParaview( gsCompositeGeom<2,T> const & Geo, 
                      std::string const & fn, unsigned npts=NS, bool mesh = false, bool ctrNet = false );


/// \brief Export i-th Basis function to paraview file
///
/// \param i index of a basis function
/// \param basis a basis object
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each curve
template<class T>
void gsWriteParaview_basisFnct(int i, gsBasis<T> const& basis, 
                               std::string const & fn, unsigned npts =NS);


/// \brief Export a Geometry slice to paraview file
///
/// \param Geo a gsGeometrySlice
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each curve
template<class T>
void gsWriteParaview(const gsGeometrySlice<T> & Geo,
                     std::string const & fn, unsigned npts =NS);


/// \brief Export a function plot to paraview file
///
/// \param func a function object
/// \param supp a matrix with two columns defining (lower and upper
/// corner of) a box: this is the domain over which we shall evaluate
/// the function, after sampling
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling the domain
template<class T>
void gsWriteParaview(gsFunction<T> const& func, 
                     gsMatrix<T> const& supp, 
                     std::string const & fn, 
                     unsigned npts =NS);


/// \brief Export Basis functions to paraview files
///
/// \param basis a basis object
/// \param fn filename where paraview file is written
/// \param npts number of points used for sampling each curve
/// \param mesh if true, the parameter mesh is plotted as well
template<class T>
void gsWriteParaview(gsBasis<T> const& basis, std::string const & fn, 
                     unsigned npts =NS, bool mesh = false);


/// \brief Export 2D Point set to Paraview file
///
/// \param X  1 times n matrix of values for x direction
/// \param Y  1 times n matrix of values for y direction
/// \param fn filename where paraview file is written
template<class T>
void gsWriteParaviewPoints(gsMatrix<T> const& X, 
                           gsMatrix<T> const& Y, 
                           std::string const & fn);

/// \brief Export 3D Point set to Paraview file
///
/// \param X  1 times n matrix of values for x direction
/// \param Y  1 times n matrix of values for y direction
/// \param Z  1 times n matrix of values for z-direction
/// \param fn filename where paraview file is written
template<class T>
void gsWriteParaviewPoints(gsMatrix<T> const& X,
                           gsMatrix<T> const& Y,
                           gsMatrix<T> const& Z,
                           std::string const & fn);

/// \brief Export Point set to Paraview file
///
/// \param points matrix that contain 2D or 3D points, points are columns
/// \param fn filename where paraview file is written
template<class T>
void gsWriteParaviewPoints(gsMatrix<T> const& points, std::string const & fn);

/// \brief Export tensor-structured point set with field data to Paraview file
///
/// \param points matrix that contain 2D or 3D points, points are columns
/// \param data
/// \param np
/// \param fn filename where paraview file is written
template<class T>
void gsWriteParaviewTPgrid(gsMatrix<T> const& points,
                           gsMatrix<T> const& data,
                           const gsVector<index_t> & np,
                           std::string const & fn);

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

// function to plot a field on a single patch
template<class T>
void writeSinglePatchField(const gsFunction<T> & geometry,
                           const gsFunction<T> & parField,
                           const bool isParam,
                           std::string const & fn, unsigned npts);

// Please document
template <class T>
void plot_errors(const gsMatrix<T> & orig, 
                 const gsMatrix<T> & comp,
                 std::vector<T> const & errors,
                 std::string const & fn);


} // namespace gismo


#undef NS


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsWriteParaview.hpp)
#endif
