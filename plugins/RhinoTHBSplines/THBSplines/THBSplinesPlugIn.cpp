// THBSplinesPlugIn.cpp : defines the initialization routines for the plug-in.
//

#include "StdAfx.h"
//#include "rhinoSdkPlugInDeclare.h"
#include "THBSplinesPlugIn.h"

// The plug-in object must be constructed before any plug-in classes derived
// from CRhinoCommand. The #pragma init_seg(lib) ensures that this happens.
#pragma warning(push)
#pragma warning(disable : 4073)
#pragma init_seg(lib)
#pragma warning(pop)

// Rhino plug-in declaration
#ifdef RHINO_V6_READY
#include "C:\Program Files\Rhino 6.0 SDK\inc\rhinoSdkPlugInDeclare.h"
#else
RHINO_PLUG_IN_DECLARE
#endif

// Rhino plug-in name
// Provide a short, friendly name for this plug-in.
RHINO_PLUG_IN_NAME(L"THB Splines");

// Rhino plug-in id
// don't change this, it is used in userdata
RHINO_PLUG_IN_ID(L"2A1816A1-7798-41E7-AFE9-AA80A125DD71");

// Rhino plug-in version
// Provide a version number string for this plug-in
#include "PluginVersion.h"

// Rhino plug-in developer declarations
// TODO: fill in the following developer declarations with
// your company information. Note, all of these declarations
// must be present or your plug-in will not load.
//
RHINO_PLUG_IN_DEVELOPER_ORGANIZATION(L"MARIN");
RHINO_PLUG_IN_DEVELOPER_ADDRESS(L"PO Box 28\r\n6700AA Wageningen");
RHINO_PLUG_IN_DEVELOPER_COUNTRY(L"The Netherlands");
RHINO_PLUG_IN_DEVELOPER_PHONE(L"+31 317 493911");
RHINO_PLUG_IN_DEVELOPER_FAX(L"-");
RHINO_PLUG_IN_DEVELOPER_EMAIL(L"support@marin.nl");
RHINO_PLUG_IN_DEVELOPER_WEBSITE(L"http://www.marin.nl");
RHINO_PLUG_IN_UPDATE_URL(L"-");


// The one and only CTHBSplinesPlugIn object
static class CTHBSplinesPlugIn thePlugIn;

CTHBSplinesPlugIn* CTHBSplinesPlugIn::GetInstance()
{
    return &thePlugIn;
}

/////////////////////////////////////////////////////////////////////////////
// CTHBSplinesPlugIn definition

CTHBSplinesPlugIn& THBSplinesPlugIn()
{
  // Return a reference to the one and only CTHBSplinesPlugIn object
  return thePlugIn;
}

CTHBSplinesPlugIn::CTHBSplinesPlugIn()    
{
  // Description:
  //   CTHBSplinesPlugIn constructor. The constructor is called when the
  //   plug-in is loaded and "thePlugIn" is constructed. Once the plug-in
  //   is loaded, CTHBSplinesPlugIn::OnLoadPlugIn() is called. The
  //   constructor should be simple and solid. Do anything that might fail in
  //   CTHBSplinesPlugIn::OnLoadPlugIn().

  // TODO: Add construction code here
  m_plugin_version = RhinoPlugInVersion();
}

CTHBSplinesPlugIn::~CTHBSplinesPlugIn()
{
}

/////////////////////////////////////////////////////////////////////////////
// Required overrides

const wchar_t* CTHBSplinesPlugIn::PlugInName() const
{
  // Description:
  //   Plug-in name display string.  This name is displayed by Rhino when
  //   loading the plug-in, in the plug-in help menu, and in the Rhino
  //   interface for managing plug-ins.

  // TODO: Return a short, friendly name for the plug-in.
  return RhinoPlugInName();
}

const wchar_t* CTHBSplinesPlugIn::PlugInVersion() const
{
  // Description:
  //   Plug-in version display string. This name is displayed by Rhino
  //   when loading the plug-in and in the Rhino interface for managing
  //   plug-ins.

  // TODO: Return the version number of the plug-in.
  return m_plugin_version;
}

GUID CTHBSplinesPlugIn::PlugInID() const
{
  // Description:
  //   Plug-in unique identifier. The identifier is used by Rhino to
  //   manage the plug-ins.

  
  // {2A1816A1-7798-41E7-AFE9-AA80A125DD71}
  return ON_UuidFromString(RhinoPlugInId());
}

/////////////////////////////////////////////////////////////////////////////
// Additional overrides

BOOL CTHBSplinesPlugIn::OnLoadPlugIn()
{
  // Description:
  //   Called after the plug-in is loaded and the constructor has been
  //   run. This is a good place to perform any significant initialization,
  //   license checking, and so on.  This function must return TRUE for
  //   the plug-in to continue to load.

  // Remarks:
  //    Plug-ins are not loaded until after Rhino is started and a default document
  //    is created.  Because the default document already exists
  //    CRhinoEventWatcher::On????Document() functions are not called for the default
  //    document.  If you need to do any document initialization/synchronization then
  //    override this function and do it here.  It is not necessary to call
  //    CPlugIn::OnLoadPlugIn() from your derived class.

  // TODO: Add plug-in initialization code here.
    //if (pWatcher)
    //{
    //    ((CRhinoOnTransformObject*)pWatcher)->Register();
    //    ((CRhinoOnTransformObject*)pWatcher)->Enable(true);
    //    ((CRhinoAfterTransformObject*)pWatcher)->Register();
    //    ((CRhinoAfterTransformObject*)pWatcher)->Enable(true);
    //}
  return TRUE;
}

void CTHBSplinesPlugIn::OnUnloadPlugIn()
{
  // Description:
  //    Called one time when plug-in is about to be unloaded. By this time,
  //    Rhino's mainframe window has been destroyed, and some of the SDK
  //    managers have been deleted. There is also no active document or active
  //    view at this time. Thus, you should only be manipulating your own objects.
  //    or tools here.

  // TODO: Add plug-in cleanup code here.
}

