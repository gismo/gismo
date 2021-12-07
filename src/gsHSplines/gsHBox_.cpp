#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsForwardDeclarations.h>
#include <gsCore/gsLinearAlgebra.h>

// #include <gsHSplines/gsHBoxContainer.h>

#include <gsHSplines/gsHBox.h>
#include <gsHSplines/gsHBox.hpp>

#include <gsHSplines/gsHBoxContainer.h>
#include <gsHSplines/gsHBoxContainer.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsHBoxContainer<1,real_t>;
    CLASS_TEMPLATE_INST gsHBoxContainer<2,real_t>;
    CLASS_TEMPLATE_INST gsHBoxContainer<3,real_t>;
    CLASS_TEMPLATE_INST gsHBoxContainer<4,real_t>;

    CLASS_TEMPLATE_INST gsHBox<1,real_t>;
    CLASS_TEMPLATE_INST gsHBox<2,real_t>;
    CLASS_TEMPLATE_INST gsHBox<3,real_t>;
    CLASS_TEMPLATE_INST gsHBox<4,real_t>;
}
