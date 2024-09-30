/** @file ThbSurfaceGrip.h
    @brief A single control point (grip) of a Truncated Hierarchical B-Spline surface.

*/
#pragma once
class CThbSurfaceGrip : public CRhinoGripObject
{
    ON_OBJECT_DECLARE(CThbSurfaceGrip);

public:

    const wchar_t* ShortDescription(bool bPlural) const override;

    void NewLocation() override;    

    CThbSurfaceGrip();

    CThbSurfaceGrip(int level);

  GetSetPropertyMacro(Level, int);

  bool IsActive() const;

    int Pick(const CRhinoPickContext& pick_context, class CRhinoObjRefArray& pick_list) const;
    bool IsVisibleInViewport(const CRhinoViewport& vp) const;
    double Weight() const override { return m_weight; }
    bool SetWeight(double weight) override;

    ~CThbSurfaceGrip();
    int m_GripNumber;
private:
    int _Level;
    double m_weight;
};

