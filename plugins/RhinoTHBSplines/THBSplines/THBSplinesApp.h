// THBSplines.h : main header file for the THBSplines DLL.
//

#pragma once

#ifndef __AFXWIN_H__
    #error "include 'stdafx.h' before including this file for PCH"
#endif

#include "resource.h" // main symbols

// CTHBSplinesApp
// See THBSplinesApp.cpp for the implementation of this class
//

class CTHBSplinesApp : public CWinApp
{
public:
  // CRITICAL: DO NOT CALL RHINO SDK FUNCTIONS HERE!
  // Only standard MFC DLL instance construction belongs here. 
  // All other significant initialization should take place in
  // CTHBSplinesPlugIn::OnLoadPlugIn().
  CTHBSplinesApp() {};

  // CRITICAL: DO NOT CALL RHINO SDK FUNCTIONS HERE!
  // Only standard MFC DLL instance initialization belongs here. 
  // All other significant initialization should take place in
  // CTHBSplinesPlugIn::OnLoadPlugIn().
    BOOL InitInstance() override;
  
  // CRITICAL: DO NOT CALL RHINO SDK FUNCTIONS HERE!
  // Only standard MFC DLL instance clean up belongs here. 
  // All other significant cleanup should take place in either
  // CTHBSplinesPlugIn::OnSaveAllSettings() or
  // CTHBSplinesPlugIn::OnUnloadPlugIn().  
    int ExitInstance() override;
  
    DECLARE_MESSAGE_MAP()
};
