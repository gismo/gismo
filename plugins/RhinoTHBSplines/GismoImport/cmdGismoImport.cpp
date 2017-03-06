// cmdGismoImport.cpp : command file
//

#include "StdAfx.h"
#include "GismoImportPlugIn.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// BEGIN GismoImport command
//

#pragma region GismoImport command

// Do NOT put the definition of class CCommandGismoImport in a header
// file. There is only ONE instance of a CCommandGismoImport class
// and that instance is the static theGismoImportCommand that appears
// immediately below the class definition.

class CCommandGismoImport : public CRhinoCommand
{
public:
  // The one and only instance of CCommandGismoImport is created below.
  // No copy constructor or operator= is required.
  // Values of member variables persist for the duration of the application.

  // CCommandGismoImport::CCommandGismoImport()
  // is called exactly once when static theGismoImportCommand is created.
  CCommandGismoImport() {};

  // CCommandGismoImport::~CCommandGismoImport()
  // is called exactly once when static theGismoImportCommand is destroyed.
  // The destructor should not make any calls to the Rhino SDK. 
  // If your command has persistent settings, then override 
  // CRhinoCommand::SaveProfile and CRhinoCommand::LoadProfile.
  ~CCommandGismoImport() {};

  // Returns a unique UUID for this command.
  // If you try to use an id that is already being used, then
  // your command will not work. Use GUIDGEN.EXE to make unique UUID.
  UUID CommandUUID() override
  {
    // {2FFFFF26-5660-414D-A114-F1CF94343A44}
    static const GUID GismoImportCommand_UUID =
    { 0x2FFFFF26, 0x5660, 0x414D, { 0xA1, 0x14, 0xF1, 0xCF, 0x94, 0x34, 0x3A, 0x44 } };
    return GismoImportCommand_UUID;
  }

  // Returns the English command name.
  // If you want to provide a localized command name, then override 
  // CRhinoCommand::LocalCommandName.
  const wchar_t* EnglishCommandName() override { return L"GismoImport"; }

  // Rhino calls RunCommand to run the command.
  CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;
};

// The one and only CCommandGismoImport object
// Do NOT create any other instance of a CCommandGismoImport class.
static class CCommandGismoImport theGismoImportCommand;

CRhinoCommand::result CCommandGismoImport::RunCommand(const CRhinoCommandContext& context)
{
  // CCommandGismoImport::RunCommand() is called when the user
  // runs the "GismoImport".

  // TODO: Add command code here.

  // Rhino command that display a dialog box interface should also support
  // a command-line, or scriptable interface.

  ON_wString str;
  str.Format(L"The \"%s\" command is under construction.\n", EnglishCommandName());
  if (context.IsInteractive())
    RhinoMessageBox(str, GismoImportPlugIn().PlugInName(), MB_OK);
  else
    RhinoApp().Print(str);

  // TODO: Return one of the following values:
  //   CRhinoCommand::success:  The command worked.
  //   CRhinoCommand::failure:  The command failed because of invalid input, inability
  //                            to compute the desired result, or some other reason
  //                            computation reason.
  //   CRhinoCommand::cancel:   The user interactively canceled the command 
  //                            (by pressing ESCAPE, clicking a CANCEL button, etc.)
  //                            in a Get operation, dialog, time consuming computation, etc.

  return CRhinoCommand::success;
}

#pragma endregion

//
// END GismoImport command
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
