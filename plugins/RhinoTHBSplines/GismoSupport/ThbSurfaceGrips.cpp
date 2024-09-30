#include "stdafx.h"
#include "ThbSurfaceUtils.h"
#include "ThbSurfaceGrip.h"
#include "ThbSurfaceGrips.h"
#include "ThbSurfaceObject.h"
#include "ThbSurfaceDrawer.h"
#include "ThbGripsState.h"
#include "ThbControlNetDrawer.h"

// {59383B12-A82A-4EEC-B43C-B708F50FE0E2}
const ON_UUID CThbSurfaceGrips::m_thbspline_grips_id =
{ 0x59383b12, 0xa82a, 0x4eec,{ 0xb4, 0x3c, 0xb7, 0x8, 0xf5, 0xf, 0xe0, 0xe2 } };

CThbSurfaceGrips* CThbSurfaceGrips::THBSplineGrips(CRhinoObjectGrips* grips)
{
    if (!grips)
        return nullptr;
    if (CRhinoGripObject::custom_grip != grips->m__grips_type)
        return nullptr;
    if (ON_UuidCompare(CThbSurfaceGrips::m_thbspline_grips_id, grips->m_grips_id))
        return nullptr;
    return static_cast<CThbSurfaceGrips*>(grips);
}

CThbSurfaceGrips::CThbSurfaceGrips(CRhinoGripObject::GRIP_TYPE type, const ON_UUID& customGripsId)
    : CRhinoObjectGrips(type)
    , m_editSrf(nullptr)
    , m_origSrf(nullptr)
    , m_drawer(nullptr)
    , m_cnDrawer(nullptr)
{
    m_grips_id = customGripsId;    
}


CThbSurfaceGrips::~CThbSurfaceGrips()
{
    if (m_editSrf) 
        delete m_editSrf;
    m_editSrf = nullptr;

    if (m_drawer)
        delete m_drawer;
    m_drawer = nullptr;

    if (m_cnDrawer)
        delete m_cnDrawer;
    m_cnDrawer = nullptr;

    // note: we don't get to delete the m_origSrf
    //       pointer, this is managed by another object.
}

void CThbSurfaceGrips::Reset()
{
    // assign the original surface to reset the edit surface
    *m_editSrf = *m_origSrf;

    // set all grip base points and grip locations to original positions.
    const gsMatrix<>& coefs = m_editSrf->coefs();
    bool rational = coefs.cols() == 4;
    for (int row = 0; row < coefs.rows(); ++row)
    {
        CRhinoGripObject* pGrip = m_grip_list[row];
        if (rational) 
        {
            ON_4dPoint origPt(coefs(row, 0), coefs(row, 1), coefs(row, 2), coefs(row, 3));
            pGrip->m_base_point = ON_3dPoint(origPt);
            pGrip->SetPoint(origPt);
        }
        else
        {
            ON_3dPoint origPt(coefs(row, 0), coefs(row, 1), coefs(row, 2));
            pGrip->m_base_point = origPt;
            pGrip->SetPoint(origPt);
        }
    }

    m_drawer->UpdateHierarchicalSurface(*m_editSrf);
    m_cnDrawer->UpdateHierarchicalSurface(*m_editSrf);
}

CRhinoObject* CThbSurfaceGrips::NewObject()
{
    if (!m_editSrf){ RhinoApp().Print("NULLPTR SUKKEL\n"); return nullptr; }
    gsTHBSpline2* editedSrf = new gsTHBSpline2(*m_editSrf);
  ON_3dmObjectAttributes attr(m_owner_object->Attributes());
  CThbSurfaceObject* newObj = new CThbSurfaceObject(attr);
    newObj->SetHierarchicalSurface(editedSrf); 
    return newObj;
}
/// <summary>
/// Creates the grips.
/// </summary>
/// <param name="srf">The SRF.</param>
/// <returns></returns>
bool CThbSurfaceGrips::CreateGrips(const gsTHBSpline2& srf)
{
    if (m_editSrf)
    {
        delete m_editSrf;
        m_editSrf = nullptr;
    }

    if (!m_editSrf)
    {        
        m_editSrf = new gsTHBSpline2(srf);
    }

    if (m_drawer)
    {
        delete m_drawer;
        m_drawer = nullptr;
    }
    if (!m_drawer)
    {
        m_drawer = new CThbSurfaceDrawer(*m_editSrf);
    }
    if (!m_cnDrawer)
    {
        m_cnDrawer = new CThbControlNetDrawer();
    }

    m_origSrf = &srf;

    const gsMatrix<>& coefs = m_editSrf->coefs();
    bool rational = coefs.cols() == 4;

    m_grip_list.SetCapacity(coefs.rows());
    for (int row = 0; row < coefs.rows(); ++row)
    {
        int levelOf = m_editSrf->basis().levelOf(row);
        CThbSurfaceGrip* pGrip = new CThbSurfaceGrip(levelOf);

        double w = rational ? coefs(row, 3) : 1.0;
        for (int i = 0; i < 3; ++i)
            pGrip->m_base_point[i] = coefs(row, i) / w;

        CRhinoGripObject::GRIP_TYPE type = pGrip->GripType();
        if (!pGrip->SetWeight(w))
        {
            return false;
        }
        // important: set the weight first and then the grip_owner
        // otherwise, the SetWeight will call the grip_owner::UpdateGrips leading to a crash
        pGrip->m_GripNumber = row;
        //pGrip->m_grips = this; - don't do this, it makes the grips un-selectable in Rhino5. Use m_grip_owner instead
        pGrip->m_grip_owner = this;
        pGrip->m_grip_index = row;

        m_grip_list.Append(pGrip);
    }

  m_drawer->UpdateHierarchicalSurface(*m_editSrf);
    m_cnDrawer->UpdateHierarchicalSurface(*m_editSrf);
    return true;
}

void CThbSurfaceGrips::UpdateGrips(int gripIndex)
{
    CRhinoGripObject* pGrip = m_grip_list[gripIndex];
    ON_3dPoint newLocation = pGrip->Point();
    double w = pGrip->Weight();

    bool rational = m_editSrf->coefs().cols() == 4;
    if (!rational && w != 1.0)
    {
        // make rational by replacing coefs
        gsMatrix<> ratCoefs(m_editSrf->coefs().rows(), 4);
        for (int i = 0; i < m_editSrf->coefs().rows(); ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                ratCoefs(i, j) = m_editSrf->coef(i, j);
            }
            double wt = i == gripIndex ? w : 1.0;
            ratCoefs(i, 3) = wt;
        }
        m_editSrf->setCoefs(ratCoefs);
    }
    else
    {
        m_editSrf->coefs()(gripIndex, 0) = newLocation.x*w;
        m_editSrf->coefs()(gripIndex, 1) = newLocation.y*w;
        m_editSrf->coefs()(gripIndex, 2) = newLocation.z*w;
        m_editSrf->coefs()(gripIndex, 3) = w;
    }
    m_drawer->UpdateHierarchicalSurface(*m_editSrf);
    m_cnDrawer->UpdateHierarchicalSurface(*m_editSrf);
}

void CThbSurfaceGrips::Draw(CRhinoDrawGripsSettings &dgs)
{
  int maxLevel = m_editSrf->basis().maxLevel();
    
    int activeLevel = CThbGripsState::GetInstance()->ActiveLevel;
    int activeRow = CThbGripsState::GetInstance()->ActiveRow;
    int activeCol = CThbGripsState::GetInstance()->ActiveColumn;
    
    for (int lvl = 0; lvl <= maxLevel; ++lvl)
    {
        if (activeLevel >= 0 && lvl != activeLevel) continue;
        ON_Color c;
        switch (lvl)
        {
            case (0):
            {
        c.SetRGBA(0,255,255,0);//.SaturatedCyan;
                break;
            }
            case (1):
            {
                c.SetRGBA(255,0,0,0);
                break;
            }
            case (2):
            {
                c.SetRGBA(0,0,255,0);
                break;
            }
            default:
            {
                c.SetRGBA(0,0,0,0);
                break;
            }
        }
        dgs.m_grip_color = c;
        for (int i = 0; i < m_grip_list.Count(); ++i)
        {
            CThbSurfaceGrip* pGrip = CThbSurfaceGrip::Cast(m_grip_list[i]);

            if (pGrip->Level != lvl)
            {
                dgs.m_grip_status[i].m_bCulled = true;                
                dgs.m_grip_status[i].m_bVisible = false;                
            }
            else
            {
                dgs.m_grip_status[i].m_bCulled = false;                
                dgs.m_grip_status[i].m_bVisible = true;
            }
        }
        /* TODO: look into drawing of the gumball for manipulation of selected control points.
        if (dgs.m_bDrawDynamicStuff)
        {
            CRhinoGumball gb;
            gb.SetToDefaultGumball();
            gb.DrawDynamicGumball(dgs.m_dp, gb_mode_nothing, ON_Line::UnsetLine, )
        }
        else if (dgs.m_bDrawStaticStuff)
        {
            CRhinoGumball gb;
            gb.SetToDefaultGumball();
            gb.DrawStaticGumball(dgs.m_dp);
        }*/
        CRhinoObjectGrips::Draw(dgs);
    }

    m_cnDrawer->Draw(dgs);

}

#pragma region THBSpline Grips Enabler

CTHBSplineGripsEnabler::CTHBSplineGripsEnabler()
{
    m_grips_id = CThbSurfaceGrips::m_thbspline_grips_id;
}

void CTHBSplineGripsEnabler::TurnOnGrips(CRhinoObject* object) const
{
    CThbSurfaceObject* obj = CThbSurfaceObject::Cast(object);
    if (obj)
    {
        obj->EnableGrips(true);
    }
}

#pragma endregion

#pragma region THBSpline Grips Registration

static CTHBSplineGripsRegistration* _instance = nullptr;

CTHBSplineGripsRegistration* CTHBSplineGripsRegistration::GetInstance()
{
    if (!_instance)
        _instance = new CTHBSplineGripsRegistration();
    return _instance;
}

void CTHBSplineGripsRegistration::DestroyInstance()
{
    delete _instance;
    _instance = nullptr;
}

CTHBSplineGripsRegistration::CTHBSplineGripsRegistration()
    : m_pGripsEnabler(nullptr)
{
}

CTHBSplineGripsRegistration::~CTHBSplineGripsRegistration()
{
    if (m_pGripsEnabler)
    {
        delete m_pGripsEnabler;
        m_pGripsEnabler = nullptr;
    }
}

void CTHBSplineGripsRegistration::RegisterGrips()
{
    if (!m_pGripsEnabler)
    {
        // register once and only once
        m_pGripsEnabler = new CTHBSplineGripsEnabler();
        RhinoApp().RegisterGripsEnabler(m_pGripsEnabler);
    }
}

#pragma endregion