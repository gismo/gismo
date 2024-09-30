// cmdTHBSplines.cpp : command file
//

#include "stdafx.h"
#include "THBSplinesPlugIn.h"
#include "TestCommandImplementation.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// BEGIN THBSplines command
//

#pragma region THBSplines command

// Do NOT put the definition of class CCommandTHBSplines in a header
// file. There is only ONE instance of a CCommandTHBSplines class
// and that instance is the static theTHBSplinesCommand that appears
// immediately below the class definition.

class CCommandTHBSplines : public CRhinoCommand
{
public:
  // The one and only instance of CCommandTHBSplines is created below.
  // No copy constructor or operator= is required.
  // Values of member variables persist for the duration of the application.

  // CCommandTHBSplines::CCommandTHBSplines()
  // is called exactly once when static theTHBSplinesCommand is created.
  CCommandTHBSplines() {};

  // CCommandTHBSplines::~CCommandTHBSplines()
  // is called exactly once when static theTHBSplinesCommand is destroyed.
  // The destructor should not make any calls to the Rhino SDK. 
  // If your command has persistent settings, then override 
  // CRhinoCommand::SaveProfile and CRhinoCommand::LoadProfile.
  ~CCommandTHBSplines() {};

  // Returns a unique UUID for this command.
  // If you try to use an id that is already being used, then
  // your command will not work. Use GUIDGEN.EXE to make unique UUID.
  UUID CommandUUID() override
  {
    // {7D8735EB-6354-403E-803D-DA508AE99073}
    static const GUID THBSplinesCommand_UUID =
    { 0x7D8735EB, 0x6354, 0x403E, { 0x80, 0x3D, 0xDA, 0x50, 0x8A, 0xE9, 0x90, 0x73 } };
    return THBSplinesCommand_UUID;
  }

  // Returns the English command name.
  // If you want to provide a localized command name, then override 
  // CRhinoCommand::LocalCommandName.
  const wchar_t* EnglishCommandName() override { return L"thbTest"; }

  // Rhino calls RunCommand to run the command.
  CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;

private:
    CTestCommandImplementation m_impl;
};

// The one and only CCommandTHBSplines object
// Do NOT create any other instance of a CCommandTHBSplines class.
static class CCommandTHBSplines theTHBSplinesCommand;

CRhinoCommand::result CCommandTHBSplines::RunCommand(const CRhinoCommandContext& context)
{
    return m_impl.RunActualCommand(context, *this);
}

#pragma endregion

//
// END THBSplines command
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
