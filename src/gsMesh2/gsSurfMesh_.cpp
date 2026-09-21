#include <gsCore/gsTemplateTools.h>

#include <gsMesh2/gsSurfMesh.h>
#include <gsMesh2/gsSurfMesh.hpp>
#include <gsMesh2/IO.hpp>
#include <gsMesh2/IO_stl.hpp>
#include <gsMesh2/IO_off.hpp>
#include <gsMesh2/IO_obj.hpp>

namespace gismo {
    // Explicit instantiation for real_t (templated specialization lives in headers)
    template class gsSurfMesh<real_t>;
    template bool write_mesh<real_t>(const gsSurfMesh<real_t>&, const std::string&);
    template bool read_stl<real_t>(gsSurfMesh<real_t>&, const std::string&);
    template bool write_stl<real_t>(const gsSurfMesh<real_t>&, const std::string&);
    template bool read_off<real_t>(gsSurfMesh<real_t>&, const std::string&);
    template bool write_off<real_t>(const gsSurfMesh<real_t>&, const std::string&);
    template bool read_obj<real_t>(gsSurfMesh<real_t>&, std::istream&);
    template bool write_obj<real_t>(const gsSurfMesh<real_t>&, const std::string&);
    CLASS_TEMPLATE_INST internal::gsXml< gsSurfMesh<real_t> >;
}



