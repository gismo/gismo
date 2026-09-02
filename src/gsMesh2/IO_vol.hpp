/** @file IO_vol.hpp

    @brief Shared pieces of the gsVolMesh file readers and writers.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#pragma once

#include <gsMesh2/IO_vol.h>

#include <gsCore/gsDebug.h>

#include <algorithm>
#include <cctype>
#include <string>

namespace gismo {

template <class T>
typename gsVolMesh<T>::Cell
add_typed_cell(gsVolMesh<T>& mesh, int vtktype,
               const std::vector<typename gsVolMesh<T>::Vertex>& v)
{
    typedef typename gsVolMesh<T>::Cell Cell;

    if (v.size() != vtk_type_size(vtktype))
        return Cell();

    switch (vtktype)
    {
    case 10: return mesh.add_tet    (v[0],v[1],v[2],v[3]);
    case 14: return mesh.add_pyramid(v[0],v[1],v[2],v[3],v[4]);
    case 13: return mesh.add_prism  (v[0],v[1],v[2],v[3],v[4],v[5]);
    case 12: return mesh.add_hex    (v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7]);
    default: return Cell();
    }
}

namespace internal {

/** report how many cells of \a mesh are inverted.

    Files written by meshers do contain inverted elements.  Reordering them
    silently would hide a broken input, so the readers only count them and warn.
*/
template <class T>
unsigned int count_inverted(const gsVolMesh<T>& mesh)
{
    unsigned int n = 0;
    for (auto c : mesh.cells())
        if (mesh.volume(c) <= (T)0) ++n;
    return n;
}

} // namespace internal

template <class T>
bool read_volmesh(gsVolMesh<T>& mesh, const std::string& filename)
{
    const std::string ext = internal::volmesh_extension(filename);

    if ("msh" == ext) return read_msh(mesh, filename);
    if ("vtk" == ext) return read_vtk(mesh, filename);
    if ("vtu" == ext) return read_vtu(mesh, filename);

    gsWarn << "gsVolMesh: no reader for the extension \"" << ext << "\".\n";
    return false;
}

template <class T>
bool write_volmesh(const gsVolMesh<T>& mesh, const std::string& filename)
{
    const std::string ext = internal::volmesh_extension(filename);

    if ("msh" == ext) return write_msh(mesh, filename);
    if ("vtk" == ext) return write_vtk(mesh, filename);
    if ("vtu" == ext) return write_vtu(mesh, filename);

    gsWarn << "gsVolMesh: no writer for the extension \"" << ext << "\".\n";
    return false;
}

} // namespace gismo
