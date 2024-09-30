#include "stdafx.h"
#include "RefineNurbsCommandImplementation.h"
#include "ThbSurfaceObject.h"
#include "ThbSurfaceUserData.h"
#include "ThbSurfaceUtils.h"

CRefineNurbsCommandImplementation::CRefineNurbsCommandImplementation()
{
}


CRefineNurbsCommandImplementation::~CRefineNurbsCommandImplementation()
{
}

bool SurfaceContainsKnot(const ON_NurbsSurface* s, int direction, double knot)
{
    for (int i = 0; i < s->KnotCount(direction); ++i)
    {
        double a = s->Knot(direction, i);
        if (abs(knot - a) < 1e-8)
            return true;
    }
    return false;
}

class CThbGetSurfaceRefinement : public CRhinoGetPoint
{
public:
    CThbGetSurfaceRefinement(ON_NurbsSurface* pSrf, int level)
        : m_srf(pSrf)
    {
        m_level = level;
        m_direction = 0;
        Constrain(*m_srf);
    }

    ~CThbGetSurfaceRefinement()
    {
    }

    void SetDirection(int uv)
    {
        m_direction = uv;
    }

    void AddLineAtSelectedPoint(int& span)
    {
    int dir = (m_direction + 1) % 2;

        ON_3dPoint p = Point();
        if (p.IsValid())
        {
            double s, t;
            m_srf->GetClosestPoint(p, &s, &t);
            
      if (!SurfaceContainsKnot(m_srf, dir, m_direction == 0 ? t : s))
            {
                m_srf->InsertKnot(dir, m_direction == 0 ? t : s);                
            }
      ON_Interval iv;
      // find the unique knot index
      for(int i = 0, j = 0; i < m_srf->KnotCount(dir); ++j)
      {
        double k = m_srf->Knot(dir, i);
        if (abs((m_direction == 0 ? t : s) - k) < 1e-8)
        {  
          span = j;
          break;          
        }
        i += m_srf->KnotMultiplicity(dir, i);
      }
      }
  }

#ifdef RHINO_V6_READY    
    void DynamicDraw(CRhinoDisplayPipeline& dp,    const ON_3dPoint& at) override
    {
        auto old = dp.ObjectColor();
        
        dp.SetObjectColor(ON_Color::White);
        dp.DrawSurface(*m_srf);
        dp.SetObjectColor(old);
        
    
    // draw the iso-curve at which the refinement will start or end
        double s, t;
        m_srf->GetClosestPoint(at, &s, &t);
        ON_Curve* iso(nullptr);
        if (m_direction == 0)
        {
            iso = m_srf->IsoCurve(0, t);
        }
        else
        {
            iso = m_srf->IsoCurve(1, s);
        }
        if (iso) 
        {
            dp.DrawCurve(*iso, ON_Color::SaturatedRed);
            delete iso;
        }

    // draw a hint arrow to show the user the order in which to pick the points
    // the order matters, because if the higher point were picked before the 
    // lower point, the box indices are no longer correct.
    ON_3dVector derU, derV;

    ON_3dPoint _;
    m_srf->Ev1Der(s, t, _, derU, derV);
    dp.SetObjectColor(RGB(150,75,75));
    double scale;
    dp.VP().GetWorldToScreenScale(at, &scale);
    derU /= scale/10;
    derV /= scale/10;
    
    if (m_direction == 1)
      dp.DrawDirectionArrow(at, derU);
    dp.SetObjectColor(RGB(75,150,75));

    if(m_direction == 0)
      dp.DrawDirectionArrow(at, derV);

    dp.SetObjectColor(old);
    }
#else
  void DynamicDraw(HDC hdc, CRhinoViewport& vp,    const ON_3dPoint& at) override
    {
    CRhinoDisplayPipeline* dp = vp.DisplayPipeline();
    COLORREF old = dp->DisplayAttrs()->m_ObjectColor;
    // set white
    dp->SetObjectColor(RGB(255,255,255));        
    vp.DrawSurface(*m_srf);
        
    // draw the iso-curve at which the refinement will start or end
        double s, t;
        m_srf->GetClosestPoint(at, &s, &t);
        ON_Curve* iso(nullptr);
    if (m_direction == 0)
        {
            iso = m_srf->IsoCurve(0, t);      
        }
        else
        {
            iso = m_srf->IsoCurve(1, s);
        }
        if (iso) 
        {
      dp->SetObjectColor(RGB(255,0,0));
            vp.DrawCurve(*iso);      
            delete iso;
        }
    
    // draw a hint arrow to show the user the order in which to pick the points
    // the order matters, because if the higher point were picked before the 
    // lower point, the box indices are no longer correct.
    ON_3dVector derU, derV;

    ON_3dPoint _;
    m_srf->Ev1Der(s, t, _, derU, derV);
        vp.SetDrawColor(RGB(150,75,75));
    double scale;
    vp.VP().GetWorldToScreenScale(at, &scale);
    derU /= scale/10;
    derV /= scale/10;
    
    if (m_direction == 1)
      vp.DrawDirectionArrow(at, derU);
    vp.SetDrawColor(RGB(75,150,75));

    if(m_direction == 0)
      vp.DrawDirectionArrow(at, derV);

    dp->SetObjectColor(old);
    }
#endif
private:
    int m_level;
    int m_direction;

    ON_NurbsSurface* m_srf;
};

CRhinoCommand::result CRefineNurbsCommandImplementation::RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand)
{
    CRhinoGetObject go;
    go.SetCommandPrompt(L"Select surface to refine");
    go.SetGeometryFilter(CRhinoGetObject::surface_object | CRhinoGetObject::polysrf_object);
  CRhinoGet::result gr = go.GetObjects();
    if (gr != CRhinoGet::object)
        return CRhinoCommand::success;

    CRhinoObjRef oRef = go.Object(0);

    const ON_Surface* srf = oRef.Surface();
    ON_NurbsSurface nurbs;
    srf->GetNurbForm(nurbs);


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

        int uSpanCount = (1 << level) * nurbs.SpanCount(0);
        int vSpanCount = (1 << level) * nurbs.SpanCount(1);

    int uStart(-1), uEnd(-1), vStart(-1), vEnd(1);
        
        CThbGetSurfaceRefinement gp(&nurbs, level);
        gp.SetDirection(0);
        gp.SetCommandPrompt(L"Select start [ENTER for start edge]");
        gp.AcceptNothing(true);
        auto gr = gp.GetPoint();
        if (gr == CRhinoGet::nothing)
        {
            vStart = 0;
        }
        else if (gr == CRhinoGet::point)
        {
            gp.AddLineAtSelectedPoint(vStart);
        }
        else
        {
            return gp.CommandResult();
        }
        gp.SetCommandPrompt(L"Select end [ENTER for end edge]");
        gr = gp.GetPoint();
        if (gr == CRhinoGet::nothing)
        {
            vEnd = vSpanCount;
        }
        else if (gr == CRhinoGet::point)
        {
            gp.AddLineAtSelectedPoint(vEnd);
        }
        else
        {
            return gp.CommandResult();
        }
        gp.SetDirection(1);
        gr = gp.GetPoint();
        if (gr == CRhinoGet::nothing)
        {
            uStart = 0;
        }
        else if (gr == CRhinoGet::point)
        {
            gp.AddLineAtSelectedPoint(uStart);
        }
        else
        {
            return gp.CommandResult();
        }
        gp.SetCommandPrompt(L"Select end [ENTER for start edge]");
        gr = gp.GetPoint();
        if (gr == CRhinoGet::nothing)
        {
            uEnd = uSpanCount;
        }
        else if (gr == CRhinoGet::point)
        {
            gp.AddLineAtSelectedPoint(uEnd);
        }
        else
        {
            return gp.CommandResult();
        }

        int uMin = std::min(uStart, uEnd);
        int uMax = std::max(uStart, uEnd);
        int vMin = std::min(vStart, vEnd);
        int vMax = std::max(vStart, vEnd);        

        boxes.push_back(level);
        boxes.push_back((1 << level) *(uMin));
        boxes.push_back((1 << level) *(vMin));
        boxes.push_back((1 << level) *(uMax));
        boxes.push_back((1 << level) *(vMax));
    }
    
    CThbSurfaceObject* obj = new CThbSurfaceObject();
    gsTHBSpline2* thb = new gsTHBSpline2();
    ON_GismoUtils::FromSurface(nurbs, boxes, *thb);
    obj->SetHierarchicalSurface(thb);
        
    context.m_doc.ReplaceObject(oRef, obj);

    context.m_doc.Redraw();
    return CRhinoCommand::success;
}