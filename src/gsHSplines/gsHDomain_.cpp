#include <gsCore/gsTemplateTools.h>

#include <gsHSplines/gsHDomain.h>
#include <gsHSplines/gsHDomain.hpp>

namespace gismo
{

    CLASS_TEMPLATE_INST gsHDomain<1,unsigned>;
    CLASS_TEMPLATE_INST gsHDomain<2,unsigned>;
    CLASS_TEMPLATE_INST gsHDomain<3,unsigned>;
    CLASS_TEMPLATE_INST gsHDomain<4,unsigned>;


    // Explicit member instansiations. Quite ugly for now..
    // Should be done inside gsHDomain
    #define HDOMAIN1 gsHDomain<1,unsigned> 
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
    #define HDOMAIN1 gsHDomain<2,unsigned> 
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
    #define HDOMAIN1 gsHDomain<3,unsigned> 
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
    #define HDOMAIN1 gsHDomain<4,unsigned> 
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

}
