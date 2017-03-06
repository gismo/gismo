/** @file ToPatchesImplementation.h
@brief Provides the implementation for the thbToPatches command

*/
#pragma once
#include "CommandImplementation.h"
class GISMO_SUPPORT_DLLEXPORT CToPatchesImplementation :
    public ICommandImplementation
{
public:
    CToPatchesImplementation();
    ~CToPatchesImplementation();

    CRhinoCommand::result RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand) override;
};

