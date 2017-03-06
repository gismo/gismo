
/** @file TestCommandImplementation.h
@brief Provides the implementation for the thbTest command

*/
#pragma once

#include "CommandImplementation.h"
class GISMO_SUPPORT_DLLEXPORT CTestCommandImplementation : public ICommandImplementation
{
public:
    CTestCommandImplementation();
    ~CTestCommandImplementation();
    CRhinoCommand::result RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand) override;
};

