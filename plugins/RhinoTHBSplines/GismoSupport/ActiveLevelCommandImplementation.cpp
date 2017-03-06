#include "stdafx.h"
#include "ActiveLevelCommandImplementation.h"
#include "ThbGripsState.h"

CActiveLevelCommandImplementation::CActiveLevelCommandImplementation()
{
}


CActiveLevelCommandImplementation::~CActiveLevelCommandImplementation()
{
}

CRhinoCommand::result CActiveLevelCommandImplementation::RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand)
{
    int level = -1;
#ifdef RHINO_V6_READY
    callingCommand.Settings().GetInteger(L"thbActiveLevel", level, -1, -1, INT_MAX);
#endif
    CRhinoGetInteger gi;
    gi.SetCommandPrompt(L"Set active level");
    int nOptAll = gi.AddCommandOption(RHCMDOPTNAME(L"All"));

    CRhinoGet::result gr = gi.GetInteger();
    if (gr == CRhinoGet::number)
    {
        level = (int)gi.Number();
    }
    else if (gr == CRhinoGet::option && gi.OptionIndex() == nOptAll)
    {
        level = -1;
    }
    else
    {
        return CRhinoCommand::success;
    }

#ifdef RHINO_V6_READY
    callingCommand.Settings().SetInteger(L"thbActiveLevel", level);
#endif
    CThbGripsState::GetInstance()->ActiveLevel = level;

    context.m_doc.Regen();
    return CRhinoCommand::success;
}