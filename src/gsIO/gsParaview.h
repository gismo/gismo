/** @file gsParaview.h

    @brief Provides a class for exporting to Paraview file format.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsIO/gsOptionList.h>
#include <gsIO/gsFileManager.h>
#include <gsCore/gsForwardDeclarations.h>

namespace gismo
{

/**
    @brief Class for exporting geometry and data to Paraview file format.

    This class provides an object-oriented interface for Paraview export
    functionality, using gsOptionList to manage all configuration options.

    Typical usage:
    \verbatim
    gsMultiPatch<> mp = ...;

    gsParaview<real_t> pv;
    pv.options().setInt("numPoints", 100);
    pv.options().setInt("precision", 5);
    pv.options().setSwitch("plotElements", true);
    pv.options().setSwitch("plotControlNet", true);

    pv.write(mp, "output_file");
    // or with default filename:
    pv.write(mp);  // writes to "multipatch.pvd"
    \endverbatim

    \ingroup IO
*/
template<class T>
class gsParaview
{
public:
    typedef T Scalar;

public:
    /// @brief Default constructor
    gsParaview();

    /// @brief Constructor with custom options
    gsParaview(const gsOptionList & options);

    /// @brief Destructor
    ~gsParaview();

    /// @brief Returns a reference to the options list
    gsOptionList & options() { return m_options; }

    /// @brief Returns a const reference to the options list
    const gsOptionList & options() const { return m_options; }

    /// @brief Returns the default options
    static gsOptionList defaultOptions();

    // ============================================================
    // Write methods for different geometry types
    // ============================================================

    /// @brief Export a gsGeometry to Paraview file
    void write(const gsGeometry<T> & geo,
               const std::string & fn = "geometry") const;

    /// @brief Export a gsMultiPatch to Paraview file
    /// @note If bezier option is true, uses Bezier element export.
    ///       If boundary option is true, also writes boundaries.
    ///       If interfaces option is true, also writes interfaces.
    void write(const gsMultiPatch<T> & mp,
               const std::string & fn = "multipatch") const;

    /// @brief Export a gsMultiPatch with options to Paraview file
    void write(const gsMultiPatch<T> & Geo, const std::string & fn,
               const std::vector<std::string> & props) const;

    /// @brief Export a vector of geometries to Paraview file
    void write(std::vector<gsGeometry<T>*> const & geos,
               const std::string & fn = "geometries") const;

    /// @brief Export a gsField to Paraview file
    void write(const gsField<T> & field,
               const std::string & fn = "field") const;

    /// @brief Export a gsBasis to Paraview file
    void write(const gsBasis<T> & basis,
               const std::string & fn = "basis") const;

    /// @brief Export selected basis functions to Paraview file
    void write(const gsBasis<T> & basis, const std::vector<index_t> & indices,
               const std::string & fn = "basis_functions") const;

    /// @brief Export a gsFunctionSet to Paraview file
    void write(const gsFunctionSet<T> & func,
               const std::string & fn = "functionset") const;

    /// @brief Export a gsFunction with support to Paraview file
    void write(const gsFunction<T> & func, const gsMatrix<T> & supp,
               const std::string & fn = "function") const;

    /// @brief Export a solution field defined by geometry and function to Paraview file
    void write(const gsFunctionSet<T> & geo, const gsFunctionSet<T> & func,
               const std::string & fn = "solution") const;

    /// @brief Export a gsGeometrySlice to Paraview file
    void write(const gsGeometrySlice<T> & geo,
               const std::string & fn = "geometry_slice") const;

    /// @brief Export a gsMesh to Paraview file
    void write(const gsMesh<T> & mesh,
               const std::string & fn = "mesh") const;

    /// @brief Export a gsSurfMesh to Paraview file
    void write(const gsSurfMesh<T> & mesh,
               const std::string & fn = "surfmesh",
               std::initializer_list<std::string> props = {}) const;

    /// @brief Export a gsMultiPatch in Bezier format to Paraview file
    void writeBezier(const gsMultiPatch<T> & mp,
                     const std::string & fn = "multipatch_bezier") const;

    /// @brief Export a gsMesh with parameters to Paraview file
    void write(const gsMesh<T> & mesh, const gsMatrix<T> & params,
               const std::string & fn = "mesh_param") const;

    /// @brief Export a vector of meshes to Paraview files
    void write(const std::vector<gsMesh<T>> & meshes,
               const std::string & fn = "meshes") const;

    /// @brief Export a gsTrimSurface to Paraview file
    void write(const gsTrimSurface<T> & ts,
               const std::string & fn = "trim_surface") const;

    /// @brief Export a gsPlanarDomain to Paraview file
    void write(const gsPlanarDomain<T> & pdomain,
               const std::string & fn = "planar_domain") const;

    /// @brief Export a gsCurveLoop to Paraview file
    void write(const gsCurveLoop<T> & cloop,
               const std::string & fn = "curve_loop") const;

    /// @brief Export a gsVolumeBlock to Paraview file
    void write(const gsVolumeBlock<T> & volBlock,
               const std::string & fn = "volume_block") const;

    /// @brief Export a gsSolid to Paraview file
    void write(const gsSolid<T> & solid,
               const std::string & fn = "solid") const;

    /// @brief Export basis functions of a gsMultiBasis on a gsMultiPatch to Paraview file
    void write(const gsMultiPatch<T> & mp, const gsMultiBasis<T> & mb,
               const std::string & fn = "multibasis") const;

    /// @brief Export a gsMappedSpline to Paraview file
    void write(const gsMappedSpline<2,T> & mspline,
               const std::string & fn = "mapped_spline") const;

    /// @brief Export a gsMappedBasis over a geometry to Paraview file
    void write(const gsFunctionSet<T> & geom, const gsMappedBasis<2,T> & mbasis,
               const std::string & fn = "mapped_basis",
               const std::vector<index_t> & indices = std::vector<index_t>()) const;

    /// @brief Export a gsMultiBasis (computational mesh) to Paraview file
    void write(const gsMultiBasis<T> & mb, const gsMultiPatch<T> & domain,
               const std::string & fn = "computational_mesh") const;

    /// @brief Export boundary conditions on a gsMultiPatch to Paraview file
    void write(const gsMultiPatch<T> & patches,
               const typename gsBoundaryConditions<T>::bcContainer & bcs,
               const std::string & fn = "boundary_conditions") const;

    // ============================================================
    // Write methods for boxes
    // ============================================================

    /// @brief Export boxes (gsMatrix) to Paraview file
    void write(const gsMatrix<T> & boxes,
               const std::string & fn,
               const std::vector<T> & values) const;

    /// @brief Export boxes (gsMatrix) to Paraview file
    void writeBoxes(const gsMatrix<T> & boxes,
                    const std::string & fn = "boxes",
                    const std::vector<T> & values = std::vector<T>()) const;

    /// @brief Export boxes (gsMatrix) with gsVector values to Paraview file
    void write(const gsMatrix<T> & boxes,
               const gsVector<T> & values,
               const std::string & fn = "boxes") const;

    /// @brief Export a gsHBox to Paraview file
    template<short_t d>
    void write(const gsHBox<d,T> & box,
               const std::string & fn = "hbox") const;

    /// @brief Export a gsHBoxContainer to Paraview file
    template<short_t d>
    void write(const gsHBoxContainer<d,T> & boxes,
               const std::string & fn = "hbox_container") const;

    // ============================================================
    // Write methods for points
    // ============================================================

    /// @brief Export points (gsMatrix) to Paraview file
    void writePoints(const gsMatrix<T> & points,
                     const std::string & fn = "points") const;

    /// @brief Export points (gsMatrix) with values to Paraview file
    void writePoints(const gsMatrix<T> & points,
                     const std::string & fn,
                     const gsVector<T> & values) const;

    /// @brief Export 2D points to Paraview file
    void writePoints(const gsMatrix<T> & X, const gsMatrix<T> & Y,
                     const std::string & fn = "points2d") const;

    /// @brief Export 3D points to Paraview file
    void writePoints(const gsMatrix<T> & X, const gsMatrix<T> & Y,
                     const gsMatrix<T> & Z,
                     const std::string & fn = "points3d") const;

    /// @brief Export 3D points with values to Paraview file
    void writePoints(const gsMatrix<T> & X, const gsMatrix<T> & Y,
                     const gsMatrix<T> & Z, const gsMatrix<T> & V,
                     const std::string & fn = "points3d_values") const;

    // ============================================================
    // Write methods for single basis function
    // ============================================================

    /// @brief Export i-th basis function to Paraview file
    void writeBasisFunction(int i, const gsBasis<T> & basis,
                            const std::string & fn = "basis_function") const;

    // ============================================================
    // Write methods for trimmed curves
    // ============================================================

    /// @brief Export a trimmed curve to Paraview file
    void writeTrimmedCurve(const gsTrimSurface<T> & surf,
                           unsigned idLoop, unsigned idCurve,
                           const std::string & fn = "trimmed_curve") const;

private:
    /// @brief Opens the Paraview file if "show" option is enabled
    void openIfRequested(const std::string & fn) const;

private:
    gsOptionList m_options;

}; // class gsParaview

#ifdef GISMO_WITH_PYBIND11

/**
 * @brief Initializes the Python wrapper for the class: gsParaview
 */
void pybind11_init_gsParaview(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsParaview.hpp)
#endif
