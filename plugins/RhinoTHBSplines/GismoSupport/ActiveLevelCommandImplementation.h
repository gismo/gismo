#pragma once
#include "CommandImplementation.h"
class GISMO_SUPPORT_DLLEXPORT CActiveLevelCommandImplementation :
    public ICommandImplementation
{
public:
    CActiveLevelCommandImplementation();
    ~CActiveLevelCommandImplementation();

    CRhinoCommand::result RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand);
};

