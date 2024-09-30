/** @file ThbSurfaceDrawer.h
@brief Provides drawing functionality to draw a Truncated Hierarchical B-Spline surface

*/
#pragma once
class CThbSurfaceDrawer
{
public:
    CThbSurfaceDrawer(const gsTHBSpline2& thb);
    ~CThbSurfaceDrawer();    

    /// Update the cached representation of the THB surface
    void UpdateHierarchicalSurface(const gsTHBSpline2& thb);

    /// Draw the cached surface
    void Draw(CRhinoDisplayPipeline& dp) const;

private:
  friend class CThbSurfaceObject;
    // these are the patches that are cached for drawing.
    ON_Brep m_patches;    
#ifdef RHINO_V6_READY
    CRhinoCacheHandle m_cache;
#endif
};

