#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsForwardDeclarations.h>

#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBoxContainer.h>

#include <gsHSplines/gsHBoxUtils.h>
#include <gsHSplines/gsHBoxUtils.hpp>

namespace gismo
{

STRUCT_TEMPLATE_INST gsHBoxUtils<1,real_t>;
STRUCT_TEMPLATE_INST gsHBoxUtils<2,real_t>;
STRUCT_TEMPLATE_INST gsHBoxUtils<3,real_t>;
STRUCT_TEMPLATE_INST gsHBoxUtils<4,real_t>;

STRUCT_TEMPLATE_INST gsHBoxCompare<1,real_t>;
STRUCT_TEMPLATE_INST gsHBoxCompare<2,real_t>;
STRUCT_TEMPLATE_INST gsHBoxCompare<3,real_t>;
STRUCT_TEMPLATE_INST gsHBoxCompare<4,real_t>;

STRUCT_TEMPLATE_INST gsHBoxEqual<1,real_t>;
STRUCT_TEMPLATE_INST gsHBoxEqual<2,real_t>;
STRUCT_TEMPLATE_INST gsHBoxEqual<3,real_t>;
STRUCT_TEMPLATE_INST gsHBoxEqual<4,real_t>;

STRUCT_TEMPLATE_INST gsHBoxContains<1,real_t>;
STRUCT_TEMPLATE_INST gsHBoxContains<2,real_t>;
STRUCT_TEMPLATE_INST gsHBoxContains<3,real_t>;
STRUCT_TEMPLATE_INST gsHBoxContains<4,real_t>;

STRUCT_TEMPLATE_INST gsHBoxIsContained<1,real_t>;
STRUCT_TEMPLATE_INST gsHBoxIsContained<2,real_t>;
STRUCT_TEMPLATE_INST gsHBoxIsContained<3,real_t>;
STRUCT_TEMPLATE_INST gsHBoxIsContained<4,real_t>;

}
