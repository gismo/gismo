#include <gsCore/gsTemplateTools.h>

#include <gsDomain/gsHTree.h>
#include <gsDomain/gsHTree.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsHTree<1,index_t>;
    CLASS_TEMPLATE_INST gsHTree<2,index_t>;
    CLASS_TEMPLATE_INST gsHTree<3,index_t>;
    CLASS_TEMPLATE_INST gsHTree<4,index_t>;

/*
    // Explicit member instansiations. Quite ugly for now..
    // Should be done inside gsHTree
    #define HDOMAIN1 gsHTree<1,index_t>
    TEMPLATE_INST
    HDOMAIN1::numLeaves_visitor::return_type
    HDOMAIN1::leafSearch<HDOMAIN1::numLeaves_visitor>() const;
    TEMPLATE_INST
    HDOMAIN1::printLeaves_visitor::return_type
    HDOMAIN1::leafSearch<HDOMAIN1::printLeaves_visitor>() const;
    TEMPLATE_INST
    HDOMAIN1::levelUp_visitor::return_type
    HDOMAIN1::leafSearch<HDOMAIN1::levelUp_visitor>() const;
    TEMPLATE_INST
    HDOMAIN1::levelDown_visitor::return_type
    HDOMAIN1::leafSearch<HDOMAIN1::levelDown_visitor>() const;
    TEMPLATE_INST
    HDOMAIN1::numNodes_visitor::return_type
    HDOMAIN1::nodeSearch<HDOMAIN1::numNodes_visitor>() const;
    TEMPLATE_INST
    HDOMAIN1::liftCoordsOneLevel_visitor::return_type
    HDOMAIN1::nodeSearch<HDOMAIN1::liftCoordsOneLevel_visitor>() const;
    #undef HDOMAIN1
*/
}
