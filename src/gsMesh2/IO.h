#pragma once

//== INCLUDES =================================================================


#include <gsMesh2/gsSurfMesh.h>

#include <string>


//== NAMESPACE ================================================================


namespace gismo {


//=============================================================================


bool GISMO_EXPORT read_mesh(gsSurfMesh& mesh, const std::string& filename);
bool GISMO_EXPORT read_off(gsSurfMesh& mesh, const std::string& filename);
bool GISMO_EXPORT read_obj(gsSurfMesh& mesh, const std::string& filename);
bool GISMO_EXPORT read_poly(gsSurfMesh& mesh, const std::string& filename);
bool GISMO_EXPORT read_stl(gsSurfMesh& mesh, const std::string& filename);

bool GISMO_EXPORT write_mesh(const gsSurfMesh& mesh, const std::string& filename);
bool GISMO_EXPORT write_off(const gsSurfMesh& mesh, const std::string& filename);
bool GISMO_EXPORT write_obj(const gsSurfMesh& mesh, const std::string& filename);
bool GISMO_EXPORT write_poly(const gsSurfMesh& mesh, const std::string& filename);
bool GISMO_EXPORT write_stl(const gsSurfMesh& mesh, const std::string& filename);
//bool GISMO_EXPORT write_vtk(const gsSurfMesh& mesh, const std::string& filename);


//=============================================================================
} // namespace gismo
//=============================================================================
