#pragma once

//== INCLUDES =================================================================


#include <gsMesh2/gsSurfMesh.h>

#include <string>


//== NAMESPACE ================================================================


namespace gismo {


//=============================================================================

template <class T>
bool GISMO_EXPORT read_off_ascii(gsSurfMesh<T>& mesh, char * node);

template <class T>
bool GISMO_EXPORT read_mesh(gsSurfMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT read_off(gsSurfMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT read_obj(gsSurfMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT read_poly(gsSurfMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT read_stl(gsSurfMesh<T>& mesh, const std::string& filename);

template <class T>
bool GISMO_EXPORT write_mesh(const gsSurfMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT write_off(const gsSurfMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT write_obj(const gsSurfMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT write_poly(const gsSurfMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT write_stl(const gsSurfMesh<T>& mesh, const std::string& filename);
// template <class T>
//bool GISMO_EXPORT write_vtk(const gsSurfMesh<T>& mesh, const std::string& filename);


//=============================================================================
} // namespace gismo
//=============================================================================

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(IO.hpp)
#include GISMO_HPP_HEADER(IO_off.hpp)
#include GISMO_HPP_HEADER(IO_obj.hpp)
#include GISMO_HPP_HEADER(IO_poly.hpp)
#include GISMO_HPP_HEADER(IO_stl.hpp)
#include GISMO_HPP_HEADER(IO_vtk.hpp)
#endif