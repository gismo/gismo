// cmdRefineNurbs.cpp : command file
//

#include "stdafx.h"
#include "THBSplinesPlugIn.h"
#include "RefineNurbsCommandImplementation.h"

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// BEGIN RefineNurbs command
//

#pragma region RefineNurbs command

// Do NOT put the definition of class CCommandRefineNurbs in a header
// file. There is only ONE instance of a CCommandRefineNurbs class
// and that instance is the static theRefineNurbsCommand that appears
// immediately below the class definition.

class CCommandRefineNurbs : public CRhinoCommand
{
public:
    // The one and only instance of CCommandRefineNurbs is created below.
    // No copy constructor or operator= is required.
    // Values of member variables persist for the duration of the application.

    // CCommandRefineNurbs::CCommandRefineNurbs()
    // is called exactly once when static theRefineNurbsCommand is created.
  CCommandRefineNurbs() {};

    // CCommandRefineNurbs::~CCommandRefineNurbs()
    // is called exactly once when static theRefineNurbsCommand is destroyed.
    // The destructor should not make any calls to the Rhino SDK. 
    // If your command has persistent settings, then override 
    // CRhinoCommand::SaveProfile and CRhinoCommand::LoadProfile.
  ~CCommandRefineNurbs() {};

    // Returns a unique UUID for this command.
    // If you try to use an id that is already being used, then
    // your command will not work. Use GUIDGEN.EXE to make unique UUID.
    UUID CommandUUID() override
    {
        // {22CD3FA5-15E3-4351-99A9-AF87C69C02D7}
        static const GUID RefineNurbsCommand_UUID =
        { 0x22cd3fa5, 0x15e3, 0x4351,{ 0x99, 0xa9, 0xaf, 0x87, 0xc6, 0x9c, 0x2, 0xd7 } };
        return RefineNurbsCommand_UUID;
    }

    // Returns the English command name.
    // If you want to provide a localized command name, then override 
    // CRhinoCommand::LocalCommandName.
    const wchar_t* EnglishCommandName() override { return L"thbRefine"; }

    // Rhino calls RunCommand to run the command.
    CRhinoCommand::result RunCommand(const CRhinoCommandContext& context) override;
private:
    CRefineNurbsCommandImplementation m_impl;
};

// The one and only CCommandRefineNurbs object
// Do NOT create any other instance of a CCommandRefineNurbs class.
static class CCommandRefineNurbs theRefineNurbsCommand;

CRhinoCommand::result CCommandRefineNurbs::RunCommand(const CRhinoCommandContext& context)
{
    return m_impl.RunActualCommand(context, *this);
}

#pragma endregion

//
// END RefineNurbs command
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
