// cmdToPatches.cpp : command file
//

#include "stdafx.h"
#include "THBSplinesPlugIn.h"
#include "ToPatchesImplementation.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// BEGIN ToPatches command
//

#pragma region ToPatches command

// Do NOT put the definition of class CCommandToPatches in a header
// file. There is only ONE instance of a CCommandToPatches class
// and that instance is the static theToPatchesCommand that appears
// immediately below the class definition.

class CCommandToPatches : public CRhinoCommand
{
public:
    // The one and only instance of CCommandToPatches is created below.
    // No copy constructor or operator= is required.
    // Values of member variables persist for the duration of the application.

    // CCommandToPatches::CCommandToPatches()
    // is called exactly once when static theToPatchesCommand is created.
  CCommandToPatches() {};

    // CCommandToPatches::~CCommandToPatches()
    // is called exactly once when static theToPatchesCommand is destroyed.
    // The destructor should not make any calls to the Rhino SDK. 
    // If your command has persistent settings, then override 
    // CRhinoCommand::SaveProfile and CRhinoCommand::LoadProfile.
  ~CCommandToPatches() {};

    // Returns a unique UUID for this command.
    // If you try to use an id that is already being used, then
    // your command will not work. Use GUIDGEN.EXE to make unique UUID.
    UUID CommandUUID() override
    {
        // {39CD9237-A438-4FB5-BB93-540E2A0EB723}
        static const GUID ToPatchesCommand_UUID =
        { 0x39cd9237, 0xa438, 0x4fb5,{ 0xbb, 0x93, 0x54, 0xe, 0x2a, 0xe, 0xb7, 0x23 } };
        return ToPatchesCommand_UUID;
    }

    // Returns the English command name.
    // If you want to provide a localized command name, then override 
    // CRhinoCommand::LocalCommandName.
    const wchar_t* EnglishCommandName() override { return L"thbToPatches"; }

    // Rhino calls RunCommand to run the command.
    CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;
private:
    CToPatchesImplementation m_impl;
};

// The one and only CCommandToPatches object
// Do NOT create any other instance of a CCommandToPatches class.
static class CCommandToPatches theToPatchesCommand;

CRhinoCommand::result CCommandToPatches::RunCommand(const CRhinoCommandContext& context)
{
    return m_impl.RunActualCommand(context, *this);
}

#pragma endregion

//
// END ToPatches command
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
