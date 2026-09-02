#include <gsCore/gsTemplateTools.h>

#include <gsMesh2/gsVolMesh.h>
#include <gsMesh2/gsVolMesh.hpp>

#include <gsMesh2/IO_vol.h>
#include <gsMesh2/IO_vol.hpp>
#include <gsMesh2/IO_msh.hpp>
#include <gsMesh2/IO_vtk.hpp>
#include <gsMesh2/IO_vtu.hpp>

namespace gismo
{

CLASS_TEMPLATE_INST gsVolMesh<real_t>;

template bool read_msh <real_t>(gsVolMesh<real_t>&, const std::string&);
template bool write_msh<real_t>(const gsVolMesh<real_t>&, const std::string&);
template bool read_vtk <real_t>(gsVolMesh<real_t>&, const std::string&);
template bool write_vtk<real_t>(const gsVolMesh<real_t>&, const std::string&);
template bool read_vtu <real_t>(gsVolMesh<real_t>&, const std::string&);
template bool write_vtu<real_t>(const gsVolMesh<real_t>&, const std::string&);
template bool read_volmesh <real_t>(gsVolMesh<real_t>&, const std::string&);
template bool write_volmesh<real_t>(const gsVolMesh<real_t>&, const std::string&);

} // namespace gismo
