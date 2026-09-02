/** @file IO_vol.h

    @brief File readers and writers for gsVolMesh: Gmsh, VTK and VTU.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsMesh2/gsVolMesh.h>

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

namespace gismo {

//=============================================================================
// Gmsh (.msh).  Reads ASCII versions 2.2 and 4.1, writes ASCII 2.2.
//=============================================================================

template <class T>
bool GISMO_EXPORT read_msh(gsVolMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT write_msh(const gsVolMesh<T>& mesh, const std::string& filename);

//=============================================================================
// Legacy VTK (.vtk), ASCII DATASET UNSTRUCTURED_GRID
//=============================================================================

template <class T>
bool GISMO_EXPORT read_vtk(gsVolMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT write_vtk(const gsVolMesh<T>& mesh, const std::string& filename);

//=============================================================================
// VTK XML unstructured grid (.vtu), inline ASCII data
//=============================================================================

template <class T>
bool GISMO_EXPORT read_vtu(gsVolMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT write_vtu(const gsVolMesh<T>& mesh, const std::string& filename);

//=============================================================================
// Dispatch on the file extension
//=============================================================================

template <class T>
bool GISMO_EXPORT read_volmesh(gsVolMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT write_volmesh(const gsVolMesh<T>& mesh, const std::string& filename);

namespace internal {

/// lower-case extension of \a filename, without the dot; empty if there is none
inline std::string volmesh_extension(const std::string& filename)
{
    const std::string::size_type dot = filename.rfind('.');
    if (std::string::npos == dot) return std::string();
    std::string ext = filename.substr(dot+1);
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return ext;
}

} // namespace internal

/** @brief add one cell of VTK type \a vtktype, whose vertices are given in the
    canonical order of that type.

    The Gmsh element types 4, 5, 6 and 7 use the same vertex order as the VTK
    types 10, 12, 13 and 14, so all three readers share this function.
    \returns an invalid cell if the type is unsupported or the vertex count
    does not match it.
*/
template <class T>
typename gsVolMesh<T>::Cell
add_typed_cell(gsVolMesh<T>& mesh, int vtktype,
               const std::vector<typename gsVolMesh<T>::Vertex>& v);

/// the VTK cell type corresponding to Gmsh element type \a gmshtype, or 0 if
/// that element is not a 3D cell (points, lines, triangles and quadrangles are
/// boundary markers in a Gmsh file, not cells)
inline int vtk_type_of_gmsh(int gmshtype)
{
    switch (gmshtype)
    {
    case 4: return 10;   // tetrahedron
    case 5: return 12;   // hexahedron
    case 6: return 13;   // prism / wedge
    case 7: return 14;   // pyramid
    default: return 0;
    }
}

/// the Gmsh element type corresponding to VTK cell type \a vtktype, or 0
inline int gmsh_type_of_vtk(int vtktype)
{
    switch (vtktype)
    {
    case 10: return 4;
    case 12: return 5;
    case 13: return 6;
    case 14: return 7;
    default: return 0;
    }
}

/// the number of vertices of VTK cell type \a vtktype, or 0 if it is not one
/// of the four standard volumetric types
inline unsigned int vtk_type_size(int vtktype)
{
    switch (vtktype)
    {
    case 10: return 4;
    case 12: return 8;
    case 13: return 6;
    case 14: return 5;
    default: return 0;
    }
}

//=============================================================================
} // namespace gismo
//=============================================================================

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(IO_vol.hpp)
#include GISMO_HPP_HEADER(IO_msh.hpp)
#include GISMO_HPP_HEADER(IO_vtk.hpp)
#include GISMO_HPP_HEADER(IO_vtu.hpp)
#endif
