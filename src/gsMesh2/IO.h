
#pragma once

#include <gsMesh2/gsSurfMesh.h>

#include <string>

namespace gismo {


//=============================================================================

template <class T>
bool GISMO_EXPORT read_off_ascii(gsSurfMesh<T>& mesh, char * node);
template <class T>
bool GISMO_EXPORT read_off(gsSurfMesh<T>& mesh, const std::string& filename);

template <class T>
bool GISMO_EXPORT read_obj(gsSurfMesh<T>& mesh, std::istream& in);
template <class T>
bool GISMO_EXPORT read_stl(gsSurfMesh<T>& mesh, const std::string& filename);

template <class T>
bool GISMO_EXPORT write_mesh(const gsSurfMesh<T>& mesh, const std::string& filename);

template <class T>
bool GISMO_EXPORT write_off(const gsSurfMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT write_obj(const gsSurfMesh<T>& mesh, const std::string& filename);
template <class T>
bool GISMO_EXPORT write_stl(const gsSurfMesh<T>& mesh, const std::string& filename);

//=============================================================================
} // namespace gismo
//=============================================================================

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(IO.hpp)
#include GISMO_HPP_HEADER(IO_off.hpp)
#include GISMO_HPP_HEADER(IO_obj.hpp)
#include GISMO_HPP_HEADER(IO_stl.hpp)
#endif
