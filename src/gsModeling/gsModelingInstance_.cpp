#include <gsCore/gsTemplateTools.h>

#include <gsModeling/gsSolid.h>
#include <gsModeling/gsSolid.hpp>

#include <gsModeling/gsTriMeshToSolid.h>
#include <gsModeling/gsTriMeshToSolid.hpp>

namespace gismo
{

// Explicit instantiation

#define T real_t
#define uZ unsigned
#define Z int


CLASS_TEMPLATE_INST gsSolid<T>;
CLASS_TEMPLATE_INST gsTriMeshToSolid<T>;


#undef T
#undef uZ
#undef Z

} // namespace gismo
