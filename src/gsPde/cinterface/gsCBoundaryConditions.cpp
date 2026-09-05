#include <gsPde/gsBoundaryConditions.h>
#include <gsCore/gsFunctionSet.h>
#include <gsCInterface/gsCTypes.h>
#include <gsMatrix/cinterface/gsCMemory.h>
#include <gsCInterface/gsMacros.h>
#include <gsCInterface/gsCError.h>

#ifdef __cplusplus
extern "C"
{
#endif

using namespace gismo;

GISMO_EXPORT gsCBoundaryConditions * gsBoundaryConditions_create()
{
    GISMO_CAPI_BEGIN
    return RICAST_CBC(new gsBoundaryConditions<double>());
    GISMO_CAPI_END(NULL)
}

GISMO_EXPORT void gsBoundaryConditions_addCondition(gsCBoundaryConditions * bc,
                                                    int patch,
                                                    int side,
                                                    int ctype,
                                                    gsCFunctionSet * fun,
                                                    int unknown,
                                                    bool parametric,
                                                    int component)
{
    GISMO_CAPI_BEGIN
    boxSide bside(side);
    gsFunctionSet<double> * f_ptr = RICAST_F(fun);
    RICAST_BC(bc)->addCondition(patch,
                                bside,
                                (condition_type::type)ctype,
                                f_ptr,
                                unknown,
                                parametric,
                                component);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBoundaryConditions_addCornerValue(gsCBoundaryConditions * bc,
                                                      int corner,
                                                      double value,
                                                      int patch,
                                                      int unknown,
                                                      int component)
{
    GISMO_CAPI_BEGIN
    boxCorner bcorner(corner);
    RICAST_BC(bc)->addCornerValue(bcorner, value, patch, unknown, component);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBoundaryConditions_setGeoMap(gsCBoundaryConditions * bc, gsCFunctionSet * gm)
{
    GISMO_CAPI_BEGIN
    RICAST_BC(bc)->setGeoMap(*RICAST_F(gm));
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBoundaryConditions_print(gsCBoundaryConditions * bc)
{
    GISMO_CAPI_BEGIN
    RICAST_BC(bc)->print(gsInfo);
    GISMO_CAPI_END_VOID
}

GISMO_EXPORT void gsBoundaryConditions_delete(gsCBoundaryConditions * bc)
{
    GISMO_CAPI_BEGIN
    delete RICAST_BC(bc);
    GISMO_CAPI_END_VOID
}

#ifdef __cplusplus
}
#endif
