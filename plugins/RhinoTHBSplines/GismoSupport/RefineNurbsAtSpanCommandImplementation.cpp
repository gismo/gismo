#include "stdafx.h"
#include "RefineNurbsAtSpanCommandImplementation.h"
#include "ThbSurfaceObject.h"
#include "ThbSurfaceUserData.h"
#include "ThbSurfaceUtils.h"

CRefineNurbsAtSpanCommandImplementation::CRefineNurbsAtSpanCommandImplementation()
{
}


CRefineNurbsAtSpanCommandImplementation::~CRefineNurbsAtSpanCommandImplementation()
{
}

CRhinoCommand::result CRefineNurbsAtSpanCommandImplementation::RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand)
{
    CRhinoGetObject go;
    go.SetCommandPrompt(L"Select surface to refine");
    go.SetGeometryFilter(CRhinoGetObject::surface_object | CRhinoGetObject::polysrf_object);
    CRhinoGet::result gr = go.GetObjects();
    if (gr != CRhinoGet::object)
        return CRhinoCommand::success;

    int u(-1), v(-1);

    CRhinoObjRef oRef = go.Object(0);
    const CRhinoObject*pObj = oRef.Object();
    const CThbSurfaceObject* pSrf = CThbSurfaceObject::Cast(pObj);
    if (pSrf)
    {
        CThbSurfaceUserData* sud = CThbSurfaceUserData::Cast(pSrf->GetGeometryUserData(CThbSurfaceUserData::Id()));
        if (sud)
        {
            const gsTHBSpline2* orig = sud->HierarchicalSurface();
            const gsTHBSplineBasis<2>& b = orig->basis();
            // this seems to work but looks like a bit of a hack
            // b.numBreaks does not work, at least not on an unrefined surface, though.
            u = b.getBases()[0]->component(0).size();
            u -= b.degree(0);
            v = b.getBases()[0]->component(1).size();
            v -= b.degree(1);
        }
    }
    else
    {
        const ON_Surface* srf = oRef.Surface();
        u = srf->SpanCount(0);
        v = srf->SpanCount(1);
    }


    std::vector<unsigned int> boxes;
    while (true)
    {
        CRhinoGetInteger gi;
        ON_wString p;
        p.Format(L"Give refinement level %s", boxes.size() > 0 ? L"[ENTER] to continue" : L"");
        gi.SetCommandPrompt(p);
        gi.AcceptNothing(true);
        gr = gi.GetNumber();

        if (gr == CRhinoGet::nothing)
            break;

        if (gr != CRhinoGet::number)
            return CRhinoCommand::success;
        int level = (int)gi.Number();
        int uSpanCount = (1 << level) * u;
        int vSpanCount = (1 << level) * v;


        p.Format(L"Give u-start coordinate [0 - %i]", uSpanCount);
        gi.SetCommandPrompt(p);
        gr = gi.GetInteger();
        if (gr != CRhinoGet::number)
            return CRhinoCommand::success;
        int uStart = (int)gi.Number();

        p.Format(L"Give u-end coordinate <%i - %i]", uStart, uSpanCount);
        gi.SetCommandPrompt(p);
        gr = gi.GetInteger();
        if (gr != CRhinoGet::number)
            return CRhinoCommand::success;
        int uEnd = (int)gi.Number();

        p.Format(L"Give v-start coordinate [0 - %i]", vSpanCount);
        gi.SetCommandPrompt(p);
        gr = gi.GetInteger();
        if (gr != CRhinoGet::number)
            return CRhinoCommand::success;
        int vStart = gi.Number();

        p.Format(L"Give v-end coordinate <%i - %i]", vStart, vSpanCount);
        gi.SetCommandPrompt(p);
        gr = gi.GetInteger();
        if (gr != CRhinoGet::number)
            return CRhinoCommand::success;
        int vEnd = (int)gi.Number();

        boxes.push_back(level);
        boxes.push_back(uStart);
        boxes.push_back(vStart);
        boxes.push_back(uEnd);
        boxes.push_back(vEnd);
    }

    CThbSurfaceObject *obj(nullptr);
    if (pSrf)
    {
        CThbSurfaceUserData* sud = CThbSurfaceUserData::Cast(pSrf->GetGeometryUserData(CThbSurfaceUserData::Id()));
        if (sud)
        {
            const gsTHBSpline2* orig = sud->HierarchicalSurface();
            gsTHBSpline2* refined = new gsTHBSpline2(*orig);
            refined->refineElements(boxes);
            obj = new CThbSurfaceObject();
            obj->SetHierarchicalSurface(refined);
        }
    }
    else
    {
        obj = new CThbSurfaceObject();
        gsTHBSpline2* thb = new gsTHBSpline2();
        const ON_Surface* srf = oRef.Surface();
        ON_GismoUtils::FromSurface(*srf, boxes, *thb);
        obj->SetHierarchicalSurface(thb);
    }

    context.m_doc.ReplaceObject(oRef, obj);

    context.m_doc.Redraw();
    return CRhinoCommand::success;
}