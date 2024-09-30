#include "stdafx.h"
#include "gsStreamData.h"
#include "ThbSurfaceUserData.h"
#include "ThbSurfaceObject.h"
#include "ThbSurfaceWatcher.h"
#include "ThbSurfaceUtils.h"

ON_OBJECT_IMPLEMENT(CThbSurfaceUserData, ON_UserData, "d1190fa6-0b45-4502-8be6-e1a29f0ecadf");

CThbSurfaceUserData::CThbSurfaceUserData()
    : m_thb(nullptr)
{
    m_userdata_uuid = Id();
    m_application_uuid = PlugInId();
    m_userdata_copycount = 1;
}

CThbSurfaceUserData::CThbSurfaceUserData(gsTHBSpline2* thb, ON_UUID& objId)
: m_thb(thb)
, m_id(objId)
{
    m_userdata_uuid = Id();
    m_application_uuid = PlugInId();
    m_userdata_copycount = 1;
}

CThbSurfaceUserData& CThbSurfaceUserData::operator=(const CThbSurfaceUserData& src)
{
    if (&src != this)
    {
        m_thb = src.m_thb;
    }

    return *this;
}

CThbSurfaceUserData::~CThbSurfaceUserData()
{
    // if at some point this pointer needs to managed by this object
    // look at ThbSurfaceWatcher::Notify what happens there (it takes this pointer) and adjust as needed.
    //m_thb = nullptr;
}

#ifdef RHINO_V6_READY
ON_UUID CThbSurfaceUserData::Id()
{
    return CThbSurfaceUserData::m_CThbSurfaceUserData_class_rtti.Uuid();
}
#else
ON_UUID CThbSurfaceUserData::Id()
{
  return CThbSurfaceUserData::m_CThbSurfaceUserData_class_id.Uuid();
}

#endif

ON_UUID CThbSurfaceUserData::PlugInId()
{
    // THBSplines PlugIn Id
    return ON_UuidFromString("2A1816A1-7798-41E7-AFE9-AA80A125DD71");
}


BOOL CThbSurfaceUserData::GetDescription(ON_wString& description)
{
    description = L"hierarchical surface data";
    return true;
}

ON_BOOL32 CThbSurfaceUserData::Archive() const
{
    return true;
}

const int major = 0, minor = 0;
ON_BOOL32 CThbSurfaceUserData::Write(ON_BinaryArchive& a) const
{    
    // write version
    a.Write3dmChunkVersion(major, minor);
    // write object id
    a.WriteUuid(m_id);
    
    // write XML representation of THBSpline
    std::stringstream ss;
    gsStreamData xml;
  gsTHBSpline<2>& temp = *m_thb;
    xml << temp;
    xml.WriteXmlStream(ss);
    const std::string& str = ss.str();
    
    ON_wString wStr(str.c_str());
#ifdef RHINO_V6_READY
  bool rc = a.WriteWideString(wStr);
#else
  bool rc = a.WriteString(wStr);
#endif
    //*/
    return rc;
}



ON_BOOL32 CThbSurfaceUserData::Read(ON_BinaryArchive& a)
{
    // read version
    int maj(-1), min(-1);
    if (!a.Read3dmChunkVersion(&maj, &min)) return false;

    m_thb = nullptr;
    if (maj == 0 && min == 0)
    {
        // read object id
        if (!a.ReadUuid(m_id)) return false;
        
        // read xml representation of THBSpline
        ON_wString wStr;
#ifdef RHINO_V6_READY
        if (!a.ReadWideString(wStr)) return false;
#else
    if (!a.ReadString(wStr)) return false;
#endif

        const wchar_t* wChars = wStr.Array();

        char* chars = new char[wcslen(wChars) + 1];
        wcstombs(chars, wChars, wcslen(wChars));
        chars[wcslen(wChars)] = '\0';
        std::stringstream ss;
        ss << chars;

        gsStreamData fd;
        if (!fd.ReadXmlStream(ss)) return false;
        if (fd.has<gsTHBSpline<2> >())
        {
      auto temp = fd.getFirst<gsTHBSpline<2> >();
            m_thb = new gsTHBSpline2(*temp);
        }

        delete[] chars;
        //*/
    }

    if (m_thb)
    {
        CThbSurfaceWatcher* w = CThbSurfaceWatcher::GetInstance();
        w->AddUuid(m_id);
        w->Enable(true);
    }
    
    return m_thb != nullptr;
}

ON_BOOL32 CThbSurfaceUserData::Transform(const ON_Xform& x)
{
    if (!m_thb) return false;
    const gsMatrix<>& coefs = m_thb->coefs();
    gsMatrix<> transformed(coefs);

    ON_3dPoint cv;
    for (int c = 0; c < coefs.rows(); ++c)
    {
        cv.x = coefs(c, 0);
        cv.y = coefs(c, 1);
        cv.z = coefs(c, 2);
        cv.Transform(x);
        transformed(c, 0) = cv.x;
        transformed(c, 1) = cv.y;
        transformed(c, 2) = cv.z;
    }

    m_thb->setCoefs(transformed);
    return ON_UserData::Transform(x);
}