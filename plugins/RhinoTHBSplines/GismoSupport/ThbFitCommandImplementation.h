#pragma once
#include "CommandImplementation.h"
class GISMO_SUPPORT_DLLEXPORT CThbFitCommandImplementation :
    public ICommandImplementation
{
public:
    CThbFitCommandImplementation();
    ~CThbFitCommandImplementation();

    CRhinoCommand::result RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand);
};

