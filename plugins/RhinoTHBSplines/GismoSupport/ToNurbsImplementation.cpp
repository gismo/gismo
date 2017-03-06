#include "stdafx.h"
#include "ToNurbsImplementation.h"
#include "ThbSurfaceObject.h"
#include "ThbSurfaceUserData.h"
#include "ThbSurfaceUtils.h"

CToNurbsImplementation::CToNurbsImplementation()
{
}


CToNurbsImplementation::~CToNurbsImplementation()
{
}

CRhinoCommand::result CToNurbsImplementation::RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand)
{
    CRhinoGetObject go;
    go.SetCommandPrompt(L"Select hierarchical surface to convert to NURBS");
    go.SetGeometryFilter(CRhinoGetObject::polysrf_object);
    bool deleteInput(false);
#ifdef RHINO_V6_READY
    callingCommand.Settings().GetBool(L"thbToNurbs_DELETEINPUT", deleteInput, false);
#endif
    int nOptDeleteInput = go.AddCommandOptionToggle(RHCMDOPTNAME(L"DeleteInput"), RHCMDOPTVALUE(L"No"), RHCMDOPTVALUE(L"Yes"), deleteInput, &deleteInput);

    CRhinoGet::result gr = go.GetObjects();
    if (gr != CRhinoGet::object)
        return CRhinoCommand::success;

    CRhinoObjRef oRef = go.Object(0);
    const CRhinoObject* obj = oRef.Object();
    const CThbSurfaceObject* thb = CThbSurfaceObject::Cast(obj);
    if (thb)
    {
        CThbSurfaceUserData* sud = CThbSurfaceUserData::Cast(thb->GetGeometryUserData(CThbSurfaceUserData::Id()));
        if (sud)
        {
            ON_NurbsSurface ns;
            if (0 == ON_GismoUtils::NurbForm(*sud->HierarchicalSurface(), ns)) return CRhinoCommand::success;
            if (deleteInput)
                context.m_doc.ReplaceObject(oRef, ns);
            else
                context.m_doc.AddSurfaceObject(ns);
        }
        else
        {
            RhinoApp().Print("Unable to find hierarchical surface definition on hierarchical surface object. That should not happen.");
            return CRhinoCommand::success;
        }
    }
#ifdef RHINO_V6_READY
    callingCommand.Settings().SetBool(L"thbToNurbs_DELETEINPUT", deleteInput);
#endif
    context.m_doc.Redraw();
    return CRhinoCommand::success;
}