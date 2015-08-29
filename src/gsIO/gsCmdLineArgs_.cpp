#include <gsCore/gsTemplateTools.h>

#include <string>
#include <gsCore/gsForwardDeclarations.h>
#include <gsIO/gsCmdLineArgs.h>
#include <gsIO/gsCmdLineArgs.hpp>

namespace gismo
{

// Explicit instantiation

#define T real_t
#define uZ unsigned
#define Z int
#define S std::string

CLASS_TEMPLATE_INST gsArgMultiVal<T>;
CLASS_TEMPLATE_INST gsArgMultiVal<uZ>;
CLASS_TEMPLATE_INST gsArgMultiVal<Z>;
CLASS_TEMPLATE_INST gsArgMultiVal<S>;

#undef T
#undef uZ
#undef Z
#undef S

} // end namespace gismo
