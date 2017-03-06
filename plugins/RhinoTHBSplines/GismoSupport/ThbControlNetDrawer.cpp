#include "stdafx.h"
#include "ThbControlNetDrawer.h"
#include "ThbGripsState.h"

CThbControlNetDrawer::CThbControlNetDrawer()
{
    
}


CThbControlNetDrawer::~CThbControlNetDrawer()
{
    for (auto it = m_levelEdges.begin(); it != m_levelEdges.end(); ++it)
    {
        delete it->second;
    }
}


void CThbControlNetDrawer::UpdateHierarchicalSurface(const gsTHBSpline2& thb)
{
    for (auto it = m_levelEdges.begin(); it != m_levelEdges.end(); ++it)
    {
        auto edges = it->second;        
        edges->SetCount(0);
    }

    int activeLevel = CThbGripsState::GetInstance()->ActiveLevel;

    for (unsigned int lvl = 0; lvl <= thb.basis().maxLevel(); ++lvl)
    {
        gsMesh<> cn;

        if (thb.coefs().cols() == 4) // homogeneous coordinates; does not work
        {
            gsMatrix<> cvs(thb.coefs().rows(), 3);
            for (int i = 0; i < thb.coefs().rows(); ++i)
            {
                double w = thb.coef(i, 3);
                for (int k = 0; k < 3; ++k)
                    cvs(i, k) = thb.coef(i, k) / w;
            }
            thb.basis().connectivity(cvs, lvl, cn);
        }

        else
        {
            thb.basis().connectivity(thb.coefs(), lvl, cn);
        }

        for (auto it = cn.edge.begin(); it != cn.edge.end(); ++it)
        {
            const gsEdge<>& e = *it;
            gsVertex<>* from = e.source;
            gsVertex<>* to = e.target;

            ON_3dPoint f(from->coords.x(), from->coords.y(), from->coords.z());
            ON_3dPoint t(to->coords.x(), to->coords.y(), to->coords.z());

            std::map<int, ON_SimpleArray<ON_Line>* >::iterator at = m_levelEdges.find(lvl);
            if (at == m_levelEdges.end())
            {
                ON_SimpleArray<ON_Line>* edges = new ON_SimpleArray<ON_Line>();
                edges->Append(ON_Line(f, t));
                m_levelEdges[lvl] = edges;
            }
            else
            {
                at->second->Append(ON_Line(f, t));
            }
        }
    }
}

void CThbControlNetDrawer::Draw(CRhinoDrawGripsSettings& dgs) const
{
    int activeLevel = CThbGripsState::GetInstance()->ActiveLevel;
    if (activeLevel >= 0)
    {
        auto at = m_levelEdges.find(activeLevel);
        if (at == m_levelEdges.end()) return;

        DrawLevel(at->second, dgs);
        
    }
    else 
    {
        for (auto it = m_levelEdges.begin(); it != m_levelEdges.end(); ++it)
        {
            DrawLevel(it->second, dgs);
        }
    }
}

// let's hope the compiler inlines this function.
void CThbControlNetDrawer::DrawLevel(ON_SimpleArray<ON_Line>* lines, CRhinoDrawGripsSettings& dgs) 
{
    for (int i = 0; i < lines->Count(); ++i)
    {
        const ON_Line& l = (*lines)[i];
#ifdef RHINO_V6_READY
    dgs.m_dp.DrawLine(l.from, l.to, ON_Color::Gray126, 1, 0xFFFFFFFF);
#else
    dgs.m_dp.DrawLine(l.from, l.to, RGB(128,128,128), 1, 0xFFFFFFFF);
#endif
    }    
}