// GismoImport.h : main header file for the GismoImport DLL.
//

#pragma once

#ifndef __AFXWIN_H__
    #error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h" // main symbols

// CGismoImportApp
// See GismoImportApp.cpp for the implementation of this class
//

class CGismoImportApp : public CWinApp
{
public:
  // CRITICAL: DO NOT CALL RHINO SDK FUNCTIONS HERE!
  // Only standard MFC DLL instance construction belongs here. 
  // All other significant initialization should take place in
  // CGismoImportPlugIn::OnLoadPlugIn().
  CGismoImportApp() {};

  // CRITICAL: DO NOT CALL RHINO SDK FUNCTIONS HERE!
  // Only standard MFC DLL instance initialization belongs here. 
  // All other significant initialization should take place in
  // CGismoImportPlugIn::OnLoadPlugIn().
    BOOL InitInstance() override;
  
  // CRITICAL: DO NOT CALL RHINO SDK FUNCTIONS HERE!
  // Only standard MFC DLL instance clean up belongs here. 
  // All other significant cleanup should take place in either
  // CGismoImportPlugIn::OnSaveAllSettings() or
  // CGismoImportPlugIn::OnUnloadPlugIn().  
    int ExitInstance() override;
  
    DECLARE_MESSAGE_MAP()
};
