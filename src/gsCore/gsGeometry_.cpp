
#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsGeometry.h>
#include <gsCore/gsGeometry.hpp>

#include <gsCore/gsCurve.h>
#include <gsCore/gsCurve.hpp>

#include <gsCore/gsSurface.h>
#include <gsCore/gsSurface.hpp>

#include <gsCore/gsVolume.h>
#include <gsCore/gsVolume.hpp>

#include <gsCore/gsBulk.h>
#include <gsCore/gsBulk.hpp>

namespace gismo
{
    CLASS_TEMPLATE_INST gsGeometry<real_t>;
    CLASS_TEMPLATE_INST gsCurve   <real_t> ;
    CLASS_TEMPLATE_INST gsSurface <real_t> ;
    CLASS_TEMPLATE_INST gsVolume  <real_t> ;
    CLASS_TEMPLATE_INST gsBulk    <real_t> ;
}
