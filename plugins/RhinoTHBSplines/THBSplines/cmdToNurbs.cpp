// cmdToNurbs.cpp : command file
//

#include "stdafx.h"
#include "THBSplinesPlugIn.h"
#include "ToNurbsImplementation.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// BEGIN ToNurbs command
//

#pragma region ToNurbs command

// Do NOT put the definition of class CCommandToNurbs in a header
// file. There is only ONE instance of a CCommandToNurbs class
// and that instance is the static theToNurbsCommand that appears
// immediately below the class definition.

class CCommandToNurbs : public CRhinoCommand
{
public:
    // The one and only instance of CCommandToNurbs is created below.
    // No copy constructor or operator= is required.
    // Values of member variables persist for the duration of the application.

    // CCommandToNurbs::CCommandToNurbs()
    // is called exactly once when static theToNurbsCommand is created.
  CCommandToNurbs() {};

    // CCommandToNurbs::~CCommandToNurbs()
    // is called exactly once when static theToNurbsCommand is destroyed.
    // The destructor should not make any calls to the Rhino SDK. 
    // If your command has persistent settings, then override 
    // CRhinoCommand::SaveProfile and CRhinoCommand::LoadProfile.
  ~CCommandToNurbs() {};

    // Returns a unique UUID for this command.
    // If you try to use an id that is already being used, then
    // your command will not work. Use GUIDGEN.EXE to make unique UUID.
    UUID CommandUUID() override
    {
        // {D02FFD69-F4DD-49A8-A84A-8FA405530450}
        static const GUID ToNurbsCommand_UUID =
        { 0xd02ffd69, 0xf4dd, 0x49a8,{ 0xa8, 0x4a, 0x8f, 0xa4, 0x5, 0x53, 0x4, 0x50 } };
        return ToNurbsCommand_UUID;
    }

    // Returns the English command name.
    // If you want to provide a localized command name, then override 
    // CRhinoCommand::LocalCommandName.
    const wchar_t* EnglishCommandName() override { return L"thbToNurbs"; }

    // Rhino calls RunCommand to run the command.
    CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;
private:
    CToNurbsImplementation m_impl;
};

// The one and only CCommandToNurbs object
// Do NOT create any other instance of a CCommandToNurbs class.
static class CCommandToNurbs theToNurbsCommand;

CRhinoCommand::result CCommandToNurbs::RunCommand(const CRhinoCommandContext& context)
{
    return m_impl.RunActualCommand(context, *this);
}

#pragma endregion

//
// END ToNurbs command
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
