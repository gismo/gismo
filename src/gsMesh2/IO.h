
#pragma once

#include <gsMesh2/gsSurfMesh.h>

#include <string>

namespace gismo {

bool GISMO_EXPORT read_off_ascii(gsSurfMesh& mesh, char * node);
bool GISMO_EXPORT read_off(gsSurfMesh& mesh, const std::string& filename);

bool GISMO_EXPORT read_obj(gsSurfMesh& mesh, std::istream& in);
bool GISMO_EXPORT read_stl(gsSurfMesh& mesh, const std::string& filename);
    
bool GISMO_EXPORT write_mesh(const gsSurfMesh& mesh, const std::string& filename);

bool GISMO_EXPORT write_off(const gsSurfMesh& mesh, const std::string& filename);
bool GISMO_EXPORT write_obj(const gsSurfMesh& mesh, const std::string& filename);
bool GISMO_EXPORT write_stl(const gsSurfMesh& mesh, const std::string& filename);

}
