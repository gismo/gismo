#include <gsCore/gsTemplateTools.h>

//Prerequisits
#include<gsAssembler/gsQuadrature.h>
#include<gsHSplines/gsHTensorBasis.h>
#include<gsCore/gsLinearAlgebra.h>
#include<gsCore/gsFunction.h>
#include<gsNurbs/gsBSpline.h>
#include<gsUtils/gsCombinatorics.h>

#include <gsHSplines/gsQuasiInterpolate.h>
#include <gsHSplines/gsQuasiInterpolate.hpp>

namespace gismo
{

STRUCT_TEMPLATE_INST gsQuasiInterpolate<real_t>;

}
