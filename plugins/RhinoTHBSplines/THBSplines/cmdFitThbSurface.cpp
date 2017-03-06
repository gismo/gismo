// cmdThbFitSurface.cpp : command file
//

#include "stdafx.h"
#include "THBSPlinesPlugIn.h"
#include "ThbFitCommandImplementation.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// BEGIN ThbFitSurface command
//

#pragma region ThbFitSurface command

// Do NOT put the definition of class CCommandThbFitSurface in a header
// file. There is only ONE instance of a CCommandThbFitSurface class
// and that instance is the static theThbFitSurfaceCommand that appears
// immediately below the class definition.

class CCommandThbFitSurface : public CRhinoCommand
{
public:
    // The one and only instance of CCommandThbFitSurface is created below.
    // No copy constructor or operator= is required.
    // Values of member variables persist for the duration of the application.

    // CCommandThbFitSurface::CCommandThbFitSurface()
    // is called exactly once when static theThbFitSurfaceCommand is created.
  CCommandThbFitSurface() {};

    // CCommandThbFitSurface::~CCommandThbFitSurface()
    // is called exactly once when static theThbFitSurfaceCommand is destroyed.
    // The destructor should not make any calls to the Rhino SDK. 
    // If your command has persistent settings, then override 
    // CRhinoCommand::SaveProfile and CRhinoCommand::LoadProfile.
  ~CCommandThbFitSurface() {};

    // Returns a unique UUID for this command.
    // If you try to use an id that is already being used, then
    // your command will not work. Use GUIDGEN.EXE to make unique UUID.
    UUID CommandUUID() override
    {
        // {07ADE302-40CE-47E0-AAB3-3D3F79611590}
        static const GUID ThbFitSurfaceCommand_UUID =
        { 0x7ade302, 0x40ce, 0x47e0,{ 0xaa, 0xb3, 0x3d, 0x3f, 0x79, 0x61, 0x15, 0x90 } };
        return ThbFitSurfaceCommand_UUID;
    }

    // Returns the English command name.
    // If you want to provide a localized command name, then override 
    // CRhinoCommand::LocalCommandName.
    const wchar_t* EnglishCommandName() override { return L"thbFitSurface"; }

    // Rhino calls RunCommand to run the command.
    CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;

private:
    CThbFitCommandImplementation m_impl;
};

// The one and only CCommandThbFitSurface object
// Do NOT create any other instance of a CCommandThbFitSurface class.
static class CCommandThbFitSurface theThbFitSurfaceCommand;

CRhinoCommand::result CCommandThbFitSurface::RunCommand(const CRhinoCommandContext& context)
{
    return m_impl.RunActualCommand(context, *this);
}

#pragma endregion

//
// END ThbFitSurface command
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
