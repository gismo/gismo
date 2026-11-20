#include <gsCore/gsTemplateTools.h>

#include <gsHSplines/gsHTree.h>
#include <gsHSplines/gsHTree.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsHTree<1,index_t>;
    CLASS_TEMPLATE_INST gsHTree<2,index_t>;
    CLASS_TEMPLATE_INST gsHTree<3,index_t>;
    CLASS_TEMPLATE_INST gsHTree<4,index_t>;
}
