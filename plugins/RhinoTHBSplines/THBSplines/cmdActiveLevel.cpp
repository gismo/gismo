// cmdActiveLevel.cpp : command file
//

#include "stdafx.h"
#include "THBSplinesPlugIn.h"
#include "ActiveLevelCommandImplementation.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// BEGIN ActiveLevel command
//

#pragma region ActiveLevel command

// Do NOT put the definition of class CCommandActiveLevel in a header
// file. There is only ONE instance of a CCommandActiveLevel class
// and that instance is the static theActiveLevelCommand that appears
// immediately below the class definition.

class CCommandActiveLevel : public CRhinoCommand
{
public:
    // The one and only instance of CCommandActiveLevel is created below.
    // No copy constructor or operator= is required.
    // Values of member variables persist for the duration of the application.

    // CCommandActiveLevel::CCommandActiveLevel()
    // is called exactly once when static theActiveLevelCommand is created.
  CCommandActiveLevel() {};

    // CCommandActiveLevel::~CCommandActiveLevel()
    // is called exactly once when static theActiveLevelCommand is destroyed.
    // The destructor should not make any calls to the Rhino SDK. 
    // If your command has persistent settings, then override 
    // CRhinoCommand::SaveProfile and CRhinoCommand::LoadProfile.
  ~CCommandActiveLevel() {};

    // Returns a unique UUID for this command.
    // If you try to use an id that is already being used, then
    // your command will not work. Use GUIDGEN.EXE to make unique UUID.
    UUID CommandUUID() override
    {
        // {10249054-F753-4E1A-9AE5-69F71B1B719E}
        static const GUID ActiveLevelCommand_UUID =
        { 0x10249054, 0xf753, 0x4e1a,{ 0x9a, 0xe5, 0x69, 0xf7, 0x1b, 0x1b, 0x71, 0x9e } };
        return ActiveLevelCommand_UUID;
    }

    // Returns the English command name.
    // If you want to provide a localized command name, then override 
    // CRhinoCommand::LocalCommandName.
    const wchar_t* EnglishCommandName() override { return L"thbLevel"; }

    // Rhino calls RunCommand to run the command.
    CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;
private:
    CActiveLevelCommandImplementation m_impl;
};

// The one and only CCommandActiveLevel object
// Do NOT create any other instance of a CCommandActiveLevel class.
static class CCommandActiveLevel theActiveLevelCommand;

CRhinoCommand::result CCommandActiveLevel::RunCommand(const CRhinoCommandContext& context)
{
    return m_impl.RunActualCommand(context, *this);
}

#pragma endregion

//
// END ActiveLevel command
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
