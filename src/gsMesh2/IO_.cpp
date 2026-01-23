#include <gsCore/gsTemplateTools.h>
#include <gsMesh2/IO.h>
#include <gsMesh2/IO.hpp>
#include <gsMesh2/IO_off.hpp>
#include <gsMesh2/IO_obj.hpp>
#include <gsMesh2/IO_poly.hpp>
#include <gsMesh2/IO_stl.hpp>
#include <gsMesh2/IO_vtk.hpp>

namespace gismo 
{
    TEMPLATE_INST
    bool read_mesh(gsSurfMesh<real_t> & mesh, const std::string& filename);
    TEMPLATE_INST
    bool read_off(gsSurfMesh<real_t> & mesh, const std::string& filename);
    TEMPLATE_INST
    bool read_obj(gsSurfMesh<real_t> & mesh, const std::string& filename);
    TEMPLATE_INST
    bool read_poly(gsSurfMesh<real_t> & mesh, const std::string& filename);
    TEMPLATE_INST
    bool read_stl(gsSurfMesh<real_t> & mesh, const std::string& filename);
    TEMPLATE_INST
    bool read_vtk(gsSurfMesh<real_t> & mesh, const std::string& filename);

    TEMPLATE_INST
    bool write_mesh(const gsSurfMesh<real_t> & mesh, const std::string& filename);
    TEMPLATE_INST
    bool write_off(const gsSurfMesh<real_t> & mesh, const std::string& filename);
    TEMPLATE_INST
    bool write_obj(const gsSurfMesh<real_t> & mesh, const std::string& filename);
    TEMPLATE_INST
    bool write_poly(const gsSurfMesh<real_t> & mesh, const std::string& filename);
    TEMPLATE_INST
    bool write_stl(const gsSurfMesh<real_t> & mesh, const std::string& filename);
    TEMPLATE_INST
    bool write_vtk(const gsSurfMesh<real_t> & mesh, const std::string& filename);
}
