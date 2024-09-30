#include "stdafx.h"
#include "ThbSurfaceDrawer.h"
#include "ThbSurfaceUtils.h"

CThbSurfaceDrawer::CThbSurfaceDrawer(const gsTHBSpline2& thb)
{
    ON_GismoUtils::BrepForm(thb, m_patches);
    UpdateHierarchicalSurface(thb);
}

CThbSurfaceDrawer::~CThbSurfaceDrawer()
{
}

#ifdef RHINO_V6_READY

// in Rhino5, when the object needs to be drawn in shaded mode, Rhino uses the meshes
// of the object itself. In Rhino6, however, drawing in shaded mode is completely covered
// here. 
void CThbSurfaceDrawer::UpdateHierarchicalSurface(const gsTHBSpline2& thb)
{
    m_cache = CRhinoCacheHandle();
    ON_GismoUtils::BrepForm(thb, m_patches);
}

void CThbSurfaceDrawer::Draw(CRhinoDisplayPipeline& dp) const
{
  // pass 1 of the drawer
    if (dp.ObjectsShouldDrawShadedMeshes())
    {
        if (dp.GetRhinoVP()->DisplayModeIsShaded())
            dp.DrawShadedBrep(&m_patches, dp.DisplayAttrs()->m_pMaterial, const_cast<CRhinoCacheHandle*>(&m_cache));
        
    }

    // pass 2 of the drawer
    if (dp.ObjectsShouldDrawWires()) 
    {
        ON_Color c = c = dp.DisplayAttrs()->m_ObjectColor;
        dp.DrawBrep(&m_patches, c, 1, false, const_cast<CRhinoCacheHandle*>(&m_cache));
    }
}
#else
void CThbSurfaceDrawer::UpdateHierarchicalSurface(const gsTHBSpline2& thb)
{
  ON_GismoUtils::BrepForm(thb, m_patches);  
}

void CThbSurfaceDrawer::Draw(CRhinoDisplayPipeline& dp) const
{
  dp.DrawBrep(m_patches, dp.DisplayAttrs()->m_ObjectColor);
}

#endif
