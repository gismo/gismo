/////////////////////////////////////////////////////////////////////////////
// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently

#pragma once

#ifndef VC_EXTRALEAN
#define VC_EXTRALEAN              // Exclude rarely-used stuff from Windows headers
#endif

#if _MSC_VER == 1600
#pragma message("Compiler is Visual Studio 2010 - assuming Rhino 5. See stdafx.h for more details.")
#elif _MSC_VER == 1900
#pragma message("Compiler is Visual Studio 2015 - assuming Rhino 6. See stdafx.h for more details.")
#define RHINO_V6_READY
#ifdef _DEBUG
#define RHINO_LIB_DIR "C:\\Program Files\\Rhino 6.0 SDK\\lib\\Debug\\"
#else
#define RHINO_LIB_DIR "C:\\Program Files\\Rhino 6.0 SDK\\lib\\Release\\"
#endif
#else 
#error("Compiler not recognized. See stdafx.h to define RHINO_V6_READY")
#endif

// Disable warning C4100: 'identifier' : unreferenced formal parameter
#pragma warning( disable:4100 )

// Rhino plug-ins must always use the release MFC that Rhino uses.
// Rhino plug-ins that require debugging information need to be built
// with RHINO_DEBUG_PLUGIN defined.
#if defined(RHINO_DEBUG_PLUGIN) && defined(_DEBUG)
// V4 PseudoDebug plug-ins should define RHINO_DEBUG_PLUGIN, 
// but not define _DEBUG in the .vcproj file.
#error Do not define _DEBUG - use RHINO_DEBUG_PLUGIN instead
#endif

// Rhino SDK Preamble
#ifdef RHINO_V6_READY
#include "C:\Program Files\Rhino 6.0 SDK\inc\RhinoSdkStdafxPreamble.h"
#else
#include "C:\Program Files (x86)\Rhino 5.0 x64 SDK\Inc\RhinoSdkStdafxPreamble.h"
#endif

#define _ATL_CSTRING_EXPLICIT_CONSTRUCTORS    // some CString constructors will be explicit

#include <afxwin.h>           // MFC core and standard components
#include <afxext.h>           // MFC extensions

#ifndef _AFX_NO_OLE_SUPPORT
#include <afxole.h>           // MFC OLE classes
#include <afxodlgs.h>         // MFC OLE dialog classes
#include <afxdisp.h>          // MFC Automation classes
#endif // _AFX_NO_OLE_SUPPORT

#ifndef _AFX_NO_DB_SUPPORT
#include <afxdb.h>                  // MFC ODBC database classes
#endif // _AFX_NO_DB_SUPPORT

#ifndef _AFX_NO_DAO_SUPPORT
#include <afxdao.h>                  // MFC DAO database classes
#endif // _AFX_NO_DAO_SUPPORT

#include <afxdtctl.h>              // MFC support for Internet Explorer 4 Common Controls
#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>                  // MFC support for Windows Common Controls
#endif // _AFX_NO_AFXCMN_SUPPORT

//
// TODO: include additional commonly used header files here
//

#if defined(_M_X64) && defined(WIN32) && defined(WIN64)
//   The afxwin.h includes afx.h, which includes afxver_.h, 
//   which unconditionally defines WIN32  This is a bug.
//   Note, all Windows builds (32 and 64 bit) define _WIN32.
//   Only 64-bit builds define _WIN64. Never define/undefine
//   _WIN32 or _WIN64.  Only define EXACTLY one of WIN32 or WIN64.
//   See the MSDN "Predefined Macros" help file for details.
#undef WIN32
#endif

// Rhino Plug-in
#ifdef RHINO_V6_READY
#include "C:\Program Files\Rhino 6.0 SDK\inc\RhinoSdk.h"
#else
#include "C:\Program Files (x86)\Rhino 5.0 x64 SDK\Inc\RhinoSdk.h"
#endif



// Render Development Kit.
#ifdef RHINO_V6_READY
#include "C:\Program Files\Rhino 6.0 SDK\inc\RhRdkHeaders.h"
#else
#include "C:\Program Files (x86)\Rhino 5.0 x64 SDK\Inc\RhRdkHeaders.h"
#endif

#if defined(RHINO_DEBUG_PLUGIN)
// Now that all the system headers are read, we can
// safely define _DEBUG so the developers can test
// for _DEBUG in their code.
#define _DEBUG
#endif

// Rhino Plug-in Linking Pragmas
#ifdef RHINO_V6_READY
#include "C:\Program Files\Rhino 6.0 SDK\inc\rhinoSdkPlugInLinkingPragmas.h"
#else
#include "C:\Program Files (x86)\Rhino 5.0 x64 SDK\Inc\rhinoSdkPlugInLinkingPragmas.h"
#endif

