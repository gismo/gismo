/** @file  ThbSurfaceUserData.h
    @brief The user data attached to the CThbSurfaceObject's geometry. 
    This is used to facilitate reading and writing of the 
    THB surface definition to 3dm file.
*/
#pragma once

class CThbSurfaceUserData : public ON_UserData
{
public:
    ON_OBJECT_DECLARE(CThbSurfaceUserData);

public:
    CThbSurfaceUserData();
    CThbSurfaceUserData(gsTHBSpline2* thb, ON_UUID& objId);
    ~CThbSurfaceUserData();
    CThbSurfaceUserData& CThbSurfaceUserData::operator=(const CThbSurfaceUserData& src);

    static ON_UUID Id();

    static ON_UUID PlugInId();

    BOOL GetDescription(ON_wString& description) override;

    ON_BOOL32 Archive() const override;
    ON_BOOL32 Write(ON_BinaryArchive&) const override;
    ON_BOOL32 Read(ON_BinaryArchive&) override;

    ON_BOOL32 Transform(const ON_Xform&) override;

    gsTHBSpline2* HierarchicalSurface() const { return m_thb; }
    void SetModelObjectId(ON_UUID& id) { m_id = id; }

private:
    // TODO: really think about when this object gets deleted.
    //       currently, it gets deleted by the CThbSurfaceObject destructor
    //       but that does not seem logical (anymore).
    gsTHBSpline2* m_thb;
    ON_UUID m_id;
};

