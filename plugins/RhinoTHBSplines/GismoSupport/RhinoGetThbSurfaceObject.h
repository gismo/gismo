#pragma once
class CRhinoGetThbSurfaceObject : public CRhinoGetObject
{
public:
    CRhinoGetThbSurfaceObject();
    ~CRhinoGetThbSurfaceObject();

    bool CustomGeometryFilter(const CRhinoObject* object, const ON_Geometry* geometry, ON_COMPONENT_INDEX component_index) const override;
};

