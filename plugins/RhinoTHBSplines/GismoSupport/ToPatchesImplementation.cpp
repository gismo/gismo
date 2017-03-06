#include "stdafx.h"
#include "ToPatchesImplementation.h"
#include "ThbSurfaceObject.h"
#include "ThbSurfaceUserData.h"
#include "ThbSurfaceUtils.h"

CToPatchesImplementation::CToPatchesImplementation()
{
}


CToPatchesImplementation::~CToPatchesImplementation()
{
}

CRhinoCommand::result CToPatchesImplementation::RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand) 
{
    CRhinoGetObject go;
    go.SetCommandPrompt(L"Select hierarchical surface to convert to NURBS patches");
    go.SetGeometryFilter(CRhinoGetObject::polysrf_object);
    bool join(true);
    bool deleteInput(false);
#ifdef RHINO_V6_READY
    callingCommand.Settings().GetBool(L"thbToPatches_JOIN", join, true);
    callingCommand.Settings().GetBool(L"thbToPatches_DELETEINPUT", deleteInput, false);
#endif
    int nOptJoin = go.AddCommandOptionToggle(RHCMDOPTNAME(L"Join"), RHCMDOPTVALUE(L"No"), RHCMDOPTVALUE(L"Yes"), join, &join);
    int nOptDeleteInput = go.AddCommandOptionToggle(RHCMDOPTNAME(L"DeleteInput"), RHCMDOPTVALUE(L"No"), RHCMDOPTVALUE(L"Yes"), deleteInput, &deleteInput);

    while (true) 
    {
        CRhinoGet::result gr = go.GetObjects();
        if (gr == CRhinoGet::object)
            break;
        if (gr == CRhinoGet::cancel)
            return CRhinoCommand::cancel;
    }

    CRhinoObjRef oRef = go.Object(0);
    const CRhinoObject* obj = oRef.Object();
    const CThbSurfaceObject* thbObj = CThbSurfaceObject::Cast(obj);
    if (thbObj)
    {
        CThbSurfaceUserData* ud = CThbSurfaceUserData::Cast(thbObj->GetGeometryUserData(CThbSurfaceUserData::Id()));
        if (ud)
        {
            
            gsTHBSpline2& thb = *ud->HierarchicalSurface();
            
            ON_Brep brep;
            if (0 == ON_GismoUtils::BrepForm(thb, brep))
                return CRhinoCommand::failure;

            if (join)
            {
                ON_SimpleArray<ON_Brep*> breps;
                breps.SetCapacity(brep.m_F.Count());
                for (int i = 0; i < brep.m_F.Count(); ++i)
                {
                    ON_Brep* b = brep.DuplicateFace(i, false);
                    breps.Append(b);
                }
                ON_Brep* joined = RhinoJoinBreps(breps, context.m_doc.AbsoluteTolerance());
                if (joined) 
                {
                    context.m_doc.AddBrepObject(*joined);
                    delete joined;
                }
                for (int i = 0; i < breps.Count(); ++i)
                    delete breps[i];
            }
            else
            {
                for (int i = 0; i < brep.m_S.Count(); ++i)
                {
                    context.m_doc.AddSurfaceObject(*brep.m_S[i]);
                }
            }
        }
        else
        {
            RhinoApp().Print("Unable to find hierarchical surface definition on hierarchical surface object. That should not happen.");
            return CRhinoCommand::success;
        }
    }

#ifdef RHINO_V6_READY
    callingCommand.Settings().SetBool(L"thbToPatches_JOIN", join);
    callingCommand.Settings().SetBool(L"thbToPatches_DELETEINPUT", deleteInput);
#endif

    context.m_doc.Redraw();
  return CRhinoCommand::success;
}