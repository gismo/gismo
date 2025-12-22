#include <gsCore/gsTemplateTools.h>

#include <gsMesh2/gsSurfMesh.h>
#include <gsMesh2/gsSurfMesh.hpp>

namespace gismo {
    // Instantiate for real_t
    CLASS_TEMPLATE_INST gsSurfMesh<real_t>;
    CLASS_TEMPLATE_INST internal::gsXml< gsSurfMesh<real_t> >;
}
