#include "stdafx.h"

#include "ThbSurfaceGrip.h"
#include "ThbSurfaceGrips.h"
#include "ThbGripsState.h"

ON_OBJECT_IMPLEMENT(CThbSurfaceGrip, CRhinoGripObject, "5e1cea4c-016d-4c0d-9bcf-40effe3156c1")


CThbSurfaceGrip::~CThbSurfaceGrip()
{
}

CThbSurfaceGrip::CThbSurfaceGrip()
    : _Level(-1)
    , m_weight(1.0)
{

}


CThbSurfaceGrip::CThbSurfaceGrip(int level) 
  : m_weight(1.0)
{
    _Level = level;
}

bool CThbSurfaceGrip::IsActive() const
{
    int al = CThbGripsState::GetInstance()->GetActiveLevel();
    return al < 0 || al == _Level;
}

int CThbSurfaceGrip::Pick(const CRhinoPickContext& pc, class CRhinoObjRefArray& pl) const
{
  if (!IsActive()) return 0;
    return CRhinoGripObject::Pick(pc, pl);
}

bool CThbSurfaceGrip::IsVisibleInViewport(const CRhinoViewport& vp) const
{
    return IsActive();
}

const wchar_t* CThbSurfaceGrip::ShortDescription(bool bPlural) const
{
    if (bPlural) return L"hierarchical grips";
    return L"hierarchical grip";
}

void CThbSurfaceGrip::NewLocation()
{
    if (GripBasePoint() != GripLocation())
    {
        m_grip_owner->m_bGripsMoved = true;
        m_grip_owner->m_bNewLocation = true;
        
        CThbSurfaceGrips* owner =
            dynamic_cast<CThbSurfaceGrips*>(m_grip_owner);
        if (owner)
        {
            owner->UpdateGrips(this->m_GripNumber);
        }
    }
}

bool CThbSurfaceGrip::SetWeight(double weight)
{    
    m_weight = weight; 
    CThbSurfaceGrips* owner =
        dynamic_cast<CThbSurfaceGrips*>(m_grip_owner);
    if (owner)
    {
        owner->UpdateGrips(this->m_GripNumber);
    }
    return true; 
    
}