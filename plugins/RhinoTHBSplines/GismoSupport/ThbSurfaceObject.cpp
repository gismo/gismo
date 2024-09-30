#include "stdafx.h"
#include "ThbSurfaceUtils.h"
#include "ThbSurfaceObject.h"
#include "ThbSurfaceGrips.h"
#include "ThbSurfaceGrip.h"
#include "ThbSurfaceUserData.h"
#include "stdafx.h"
#include "ThbSurfaceDrawer.h"

ON_OBJECT_IMPLEMENT(CThbSurfaceObject, CRhinoBrepObject, "163f2f14-26bf-41fa-aa4b-311f77b3f46d")

void CThbSurfaceObject::CommonConstructor()
{
    m_drawer = nullptr; // will be initialized in first call to Draw  
  m_renderMeshes.SetCount(0);
}

void CThbSurfaceObject::CommonCopyConstructor(const CThbSurfaceObject& src)
{
    const ON_UserData* ud = src.GetGeometryUserData(CThbSurfaceUserData::Id());  
    const CThbSurfaceUserData* sud = CThbSurfaceUserData::Cast(ud);
    if (!sud) return;
    
    auto thb = sud->HierarchicalSurface();
    gsTHBSpline2* newThb = new gsTHBSpline2(*thb);
    SetHierarchicalSurface(newThb);

    m_drawer = nullptr; // will be initialized in first call to Draw
}


void CThbSurfaceObject::SetHierarchicalSurface(gsTHBSpline2* thb)
{
  ON_Brep* b = ON_Brep::New();
    if (0 != ON_GismoUtils::BrepForm(*thb, *b))
    {
    ON_wString err;
    ON_TextLog log(err);
    if (!b->IsValid(&log))
      RhinoApp().Print(L"Invalid BREP: %s\n", err);
        SetBrep(b);
    }
    else
    {
        delete b;
    }

    // add the user data here, so that a Transform operation
    // is applied to the user data when necessary
    ON_UUID id = ModelObjectId(); // note: id is most probably GUID::EMPTY here
    CThbSurfaceUserData* ud = new CThbSurfaceUserData(thb, id);

    if (0 == AttachGeometryUserData(ud))
    {
        delete ud;
    }
    if (m_drawer)
        m_drawer->UpdateHierarchicalSurface(*thb);
}

void CThbSurfaceObject::AddToDocNotification()
{
    // set the model object id _after_ the object has been added to the document
    // before, it is set to GUID::EMPTY.
    CThbSurfaceUserData* pUd = CThbSurfaceUserData::Cast(GetGeometryUserData(CThbSurfaceUserData::Id()));
    if (pUd)
    {
        pUd->SetModelObjectId( ModelObjectId() );
    }
}

CThbSurfaceObject::CThbSurfaceObject()
    : CRhinoBrepObject()
{
    CommonConstructor();    
}

CThbSurfaceObject::~CThbSurfaceObject()
{
    CThbSurfaceUserData* pUd = CThbSurfaceUserData::Cast(GetGeometryUserData(CThbSurfaceUserData::Id()));
    // I get a crash here if a file with a ThbSurface object is open and that same file is
    // read during Audit3dmfile. That is a bit of a corner case, I know, but may apply in 
    // other cases too (e.g. when  using a reference file).
    if (pUd)
    {
        delete pUd->HierarchicalSurface();
    }
    DetachUserData(pUd);

    if (m_drawer)
        delete m_drawer;
    m_drawer = nullptr;
}

CThbSurfaceObject::CThbSurfaceObject(const ON_3dmObjectAttributes& a)
    :CRhinoBrepObject(a)
{
    CommonConstructor();
}

CThbSurfaceObject::CThbSurfaceObject(const CThbSurfaceObject& src)
    : CRhinoBrepObject(src)
{
    CommonCopyConstructor(src);
}

CThbSurfaceObject& CThbSurfaceObject::operator=(const CThbSurfaceObject& src)
{
    if (&src != this)
    {
        CRhinoBrepObject::operator=(src);
        CommonCopyConstructor(src);
    }

    return *this;
}

void CThbSurfaceObject::Draw(CRhinoDisplayPipeline& dp) const
{
    if (!GripsOn()) // we get to draw ourselves if the grips are not on.
    {
        // late initialization of the drawer object, because if we 
        // initialize it here, the potential transformation of the 
        // user data has been done and the THB surface on the user
        // data is on the correct position. Any earlier and you may
        // miss the transformation and draw the surface on its original
        // position instead of the transformed position. 
        if (!m_drawer)
        {
            const ON_UserData* ud = GetGeometryUserData(CThbSurfaceUserData::Id());
            const CThbSurfaceUserData* sud = CThbSurfaceUserData::Cast(ud);
            if (!sud) return;

            auto thb = sud->HierarchicalSurface();
      if (!thb){ RhinoApp().Print("NULLPTR SUKKEL\n"); return; }
            const_cast<CThbSurfaceObject*>(this)->m_drawer = new CThbSurfaceDrawer(*thb);
        }
        m_drawer->Draw(dp);
    }
    else
    {
        // we don't get do draw ourselves, but ask the drawer on the
        // grips object to draw the surface "as it is being edited".
        const CThbSurfaceGrips* g = dynamic_cast<CThbSurfaceGrips*>(m_grips);
        if (g)
        {
            g->Drawer()->Draw(dp);
        }
    }
}

const wchar_t* CThbSurfaceObject::ShortDescription(bool bPlural) const
{
    return bPlural ? L"hierarchical surfaces" : L"hierarchical surface";
}

void CThbSurfaceObject::EnableGrips(bool bGripsOn)
{
    if (bGripsOn)
    {
        CTHBSplineGripsRegistration::GetInstance()->RegisterGrips();
    }

    if (!bGripsOn || (m_grips && 0 == CThbSurfaceGrips::THBSplineGrips(m_grips)))
    {
        // turn off wrong kind of grips
        CRhinoObject::EnableGrips(false);
    }

    if (bGripsOn && !m_grips)
    {
        ON_UserData* ud = GetGeometryUserData(CThbSurfaceUserData::Id());
        CThbSurfaceUserData* sud = CThbSurfaceUserData::Cast(ud);
        if (!sud) return;


        auto m_thb = sud->HierarchicalSurface();
        if (m_thb)
        {
            // turn on rectangle grips
            CThbSurfaceGrips* grips = new CThbSurfaceGrips(CRhinoGripObject::custom_grip, CThbSurfaceGrips::m_thbspline_grips_id);
            if (grips->CreateGrips(*m_thb))
                CRhinoObject::EnableCustomGrips(grips);
            else
                delete grips;
        }
    }
}

#ifndef RHINO_V6_READY

// in Rhino5, when the object needs to be drawn in shaded mode, Rhino uses the meshes
// of the object itself. In Rhino6, however, drawing in shaded mode is completely covered
// by the CThbSurfaceDrawer (it gets two passes).

int CThbSurfaceObject::GetMeshes( ON::mesh_type type, ON_SimpleArray<const ON_Mesh *>& meshes) const
{
  if (!GripsOn()) return CRhinoBrepObject::GetMeshes(type, meshes);
  
  const CThbSurfaceGrips* g = dynamic_cast<CThbSurfaceGrips*>(m_grips);

    if (g)
    {
    const_cast<CThbSurfaceObject*>(this)->DestroyMeshes(ON::any_mesh, true);
    const_cast<CThbSurfaceObject*>(this)->CreateMeshes(ON::any_mesh, ON_MeshParameters::FastRenderMesh, false);
    int i;
      for(i = 0; i < m_renderMeshes.Count(); ++i)
          meshes.Append(m_renderMeshes[i]);
    
      return i;
    }
  return 0;
}

void CThbSurfaceObject::DestroyMeshes( ON::mesh_type mesh_type, bool bDeleteMeshes ) 
{
  for(int i = 0; i < m_renderMeshes.Count(); ++i)
        delete m_renderMeshes[i];

    m_renderMeshes.SetCount(0);    
}

int CThbSurfaceObject::CreateMeshes( ON::mesh_type type, const ON_MeshParameters& param, bool bIgnoreCustom ) 
{
  if (!GripsOn())
  {
    return CRhinoBrepObject::CreateMeshes(type, param, bIgnoreCustom);
  }
  const CThbSurfaceGrips* g = dynamic_cast<CThbSurfaceGrips*>(m_grips);

    if (g)
  {
    ON_MeshParameters toUse(param);
      if (!bIgnoreCustom)
      {
          this->GetRenderMeshParameters(toUse);
      }
    ON_SimpleArray<ON_Mesh*> a;
    g->Drawer()->m_patches.CreateMesh(toUse, a);
    for(int i = 0; i < a.Count(); ++i)
      m_renderMeshes.Append(a[i]);
  }
  return m_renderMeshes.Count();
}

#endif