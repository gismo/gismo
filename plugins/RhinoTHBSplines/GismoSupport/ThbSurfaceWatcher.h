/** @file  ThbSurfaceWatcher.h
    @brief This idle-event watcher is used by the CThbSurfaceUserData
           to replace an object with a THB surface with a CThbSurfaceObject
*/
#pragma once


class CThbSurfaceWatcher :
    public CRhinoIsIdle
{
public:
    static CThbSurfaceWatcher* GetInstance();
    static void DestroyInstance();

    void Notify(const class CRhinoIsIdle::CParameters& params) override;

    void AddUuid(const ON_UUID& id)
    {
        if (!_uuids.AddUuid(id))
        {
            RhinoApp().Print("%s - UUID not added %d-%d-%d-(...).\n", __FUNCTION__, id.Data1, id.Data2, id.Data3);
        }
    }

    

private:
    CThbSurfaceWatcher();

    ON_UuidList _uuids;
};

