/** @file  ThbSurfaceObject.h
    @brief The custom object to manage a Truncated Hierarchical B-Spline surface in Rhino.
*/
#pragma once

class CThbSurfaceObject :
    public CRhinoBrepObject
{
    ON_OBJECT_DECLARE(CThbSurfaceObject);

public:
    CThbSurfaceObject();
    CThbSurfaceObject(const ON_3dmObjectAttributes&);
    CThbSurfaceObject(const CThbSurfaceObject&);
    CThbSurfaceObject& operator =(const CThbSurfaceObject&);
    ~CThbSurfaceObject(void);

    const wchar_t* ShortDescription(bool bPlural) const override;

    void Draw(CRhinoDisplayPipeline& cp) const override;

    void EnableGrips(bool bGripsOn) override;

    // after an object is added to the doc
    // the model object id is stored on the 
    // geometry user data
    void AddToDocNotification() override;

    /// set the hierarchical surface. It is stored on 
    /// a user data object, and attached to the geometry
    void SetHierarchicalSurface(gsTHBSpline2* thb);

#ifndef RHINO_V6_READY
  int GetMeshes( ON::mesh_type, ON_SimpleArray<const ON_Mesh *>& ) const override;
  void DestroyMeshes( ON::mesh_type mesh_type, bool bDeleteMeshes ) override;
  int CreateMeshes( ON::mesh_type type, const ON_MeshParameters& param, bool bIgnoreCustom ) ;
#endif

private:
    void CommonConstructor();
    void CommonCopyConstructor(const CThbSurfaceObject& src);
    class CThbSurfaceDrawer* m_drawer;
  
  ON_SimpleArray<const ON_Mesh*> m_renderMeshes;
};

