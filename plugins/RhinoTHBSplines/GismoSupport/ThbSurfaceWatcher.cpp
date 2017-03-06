#include "stdafx.h"
#include "ThbSurfaceWatcher.h"
#include "ThbSurfaceUserData.h"
#include "ThbSurfaceObject.h"


static CThbSurfaceWatcher* _instance = nullptr;
CThbSurfaceWatcher* CThbSurfaceWatcher::GetInstance()
{
    if (!_instance) _instance = new CThbSurfaceWatcher();
    return _instance;
}

void CThbSurfaceWatcher::DestroyInstance()
{
    delete _instance;
    _instance = nullptr;
}

CThbSurfaceWatcher::CThbSurfaceWatcher()
    : CRhinoIsIdle(ON_UuidFromString("2A1816A1-7798-41E7-AFE9-AA80A125DD71"))
{
    Register();
    Enable(false);
}

void CThbSurfaceWatcher::Notify(const CRhinoIsIdle::CParameters& params)
{
    Enable(false);
    CRhinoDoc* doc = RhinoApp().ActiveDoc();
    ON_SimpleArray<ON_UUID> list;
    int n = _uuids.GetUuids(list);
    
    for (int i = 0; i < n; ++i)
    {
        const ON_UUID& id = list[i];
        CRhinoObjRef objRef(id);
        
        CThbSurfaceUserData* sud = (CThbSurfaceUserData*)objRef.Object()
            ->GetGeometryUserData(CThbSurfaceUserData::Id());
        if (!sud) continue; // how did this happen

        gsTHBSpline2* thb = sud->HierarchicalSurface(); // we can now safely use this object, because the user data does not delete it.
        
        //ON_ThbSurface* srf = new ON_ThbSurface();
        //srf->GiveTHBSpline(thb);

        CThbSurfaceObject* obj = new CThbSurfaceObject();
        obj->SetHierarchicalSurface(thb);

        doc->ReplaceObject(objRef, obj);
    }
    
    doc->Regen();
    _uuids.Empty();
}