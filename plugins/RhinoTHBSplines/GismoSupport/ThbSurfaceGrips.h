/** @file ThbSurfaceGrips.h
    @brief A collection of THB surface grips. Provides drawing functionality while control points are on.

*/
#pragma once


class CThbSurfaceGrips : public CRhinoObjectGrips
{
public:
    static
        const ON_UUID m_thbspline_grips_id;
    static
        CThbSurfaceGrips* THBSplineGrips(CRhinoObjectGrips *grips);
public:
    CThbSurfaceGrips(CRhinoGripObject::GRIP_TYPE type, const ON_UUID& customGripsId);
    ~CThbSurfaceGrips();

    /// resets the grip location, e.g. when the user presses the Esc key during CP dragging
    void Reset() override;

    /// provides a new object based on the latest control point locations
    CRhinoObject* NewObject() override;

    /// draws the control points and the control point net
    void Draw(CRhinoDrawGripsSettings &) override;

    /// creates the grip objects for the given THB surface
    bool CreateGrips(const gsTHBSpline2&);

    /// Update the grip locations after the grip with given index has changed.
    void UpdateGrips(int gripIndex);

  /// get the drawer
    const class CThbSurfaceDrawer* Drawer() const { return m_drawer; }

private:
    // the edit surface is the object that is constantly changing while
    // control point(s) are being dragged.
    gsTHBSpline2* m_editSrf;

    // the drawer provides drawing for the edit surface. Note that this
    // drawer is called by the owner of the grips (CThbSurfaceObject) using
    // the Drawer() function above
    class CThbSurfaceDrawer* m_drawer;
    class CThbControlNetDrawer* m_cnDrawer;

    // a pointer to the original surface
    const gsTHBSpline2* m_origSrf;
};

class CTHBSplineGripsEnabler : public CRhinoGripsEnabler
{
public:
    CTHBSplineGripsEnabler();
    void TurnOnGrips(CRhinoObject* object) const override;
};

class CTHBSplineGripsRegistration
{
public:
    static CTHBSplineGripsRegistration* GetInstance();
    static void DestroyInstance();
    void RegisterGrips();

private:
    CTHBSplineGripsRegistration();
    ~CTHBSplineGripsRegistration();
    CTHBSplineGripsEnabler    * m_pGripsEnabler;
};

