/** @file ThbControlNetDrawer.h
    @brief This draws the control nets at each level.

*/
#pragma once

class CThbControlNetDrawer
{
public:
    CThbControlNetDrawer();
    ~CThbControlNetDrawer();

    void UpdateHierarchicalSurface(const gsTHBSpline2& thb);

    void Draw(CRhinoDrawGripsSettings& dp) const;

private:
    std::map<int, ON_SimpleArray<ON_Line>* > m_levelEdges;
    static void DrawLevel(ON_SimpleArray<ON_Line>*, CRhinoDrawGripsSettings& dp);

};

