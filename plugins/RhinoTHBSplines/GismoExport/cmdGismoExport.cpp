// cmdGismoExport.cpp : command file
//

#include "StdAfx.h"
#include "GismoExportPlugIn.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// BEGIN GismoExport command
//

#pragma region GismoExport command

// Do NOT put the definition of class CCommandGismoExport in a header
// file. There is only ONE instance of a CCommandGismoExport class
// and that instance is the static theGismoExportCommand that appears
// immediately below the class definition.

class CCommandGismoExport : public CRhinoCommand
{
public:
  // The one and only instance of CCommandGismoExport is created below.
  // No copy constructor or operator= is required.
  // Values of member variables persist for the duration of the application.

  // CCommandGismoExport::CCommandGismoExport()
  // is called exactly once when static theGismoExportCommand is created.
  CCommandGismoExport() {};

  // CCommandGismoExport::~CCommandGismoExport()
  // is called exactly once when static theGismoExportCommand is destroyed.
  // The destructor should not make any calls to the Rhino SDK. 
  // If your command has persistent settings, then override 
  // CRhinoCommand::SaveProfile and CRhinoCommand::LoadProfile.
  ~CCommandGismoExport() {};

  // Returns a unique UUID for this command.
  // If you try to use an id that is already being used, then
  // your command will not work. Use GUIDGEN.EXE to make unique UUID.
  UUID CommandUUID() override
  {
    // {6AE6959D-1364-41BF-B64C-3A838879DA8F}
    static const GUID GismoExportCommand_UUID =
    { 0x6AE6959D, 0x1364, 0x41BF, { 0xB6, 0x4C, 0x3A, 0x83, 0x88, 0x79, 0xDA, 0x8F } };
    return GismoExportCommand_UUID;
  }

  // Returns the English command name.
  // If you want to provide a localized command name, then override 
  // CRhinoCommand::LocalCommandName.
  const wchar_t* EnglishCommandName() override { return L"GismoExport"; }

  // Rhino calls RunCommand to run the command.
  CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;
};

// The one and only CCommandGismoExport object
// Do NOT create any other instance of a CCommandGismoExport class.
static class CCommandGismoExport theGismoExportCommand;

CRhinoCommand::result CCommandGismoExport::RunCommand(const CRhinoCommandContext& context)
{
  // CCommandGismoExport::RunCommand() is called when the user
  // runs the "GismoExport".

  // TODO: Add command code here.

  // Rhino command that display a dialog box interface should also support
  // a command-line, or scriptable interface.

  ON_wString str;
  str.Format(L"The \"%s\" command is under construction.\n", EnglishCommandName());
  if (context.IsInteractive())
    RhinoMessageBox(str, GismoExportPlugIn().PlugInName(), MB_OK);
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
// END GismoExport command
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
