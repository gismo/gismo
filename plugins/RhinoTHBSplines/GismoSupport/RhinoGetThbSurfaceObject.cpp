#include "stdafx.h"
#include "RhinoGetThbSurfaceObject.h"
#include "ThbSurfaceObject.h"

CRhinoGetThbSurfaceObject::CRhinoGetThbSurfaceObject()
{
}


CRhinoGetThbSurfaceObject::~CRhinoGetThbSurfaceObject()
{
}

bool CRhinoGetThbSurfaceObject::CustomGeometryFilter(const CRhinoObject* object, const ON_Geometry* geometry, ON_COMPONENT_INDEX component_index) const
{
    const CThbSurfaceObject* thb = CThbSurfaceObject::Cast(object);
    return thb != nullptr;
}