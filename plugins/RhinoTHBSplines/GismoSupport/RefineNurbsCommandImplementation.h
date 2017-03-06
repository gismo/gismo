/** @file RefineNurbsCommandImplementation.h
@brief Provides the implemementation for the thbRefineNurbs command

*/
#pragma once

#include "CommandImplementation.h"

class GISMO_SUPPORT_DLLEXPORT CRefineNurbsCommandImplementation : public ICommandImplementation
{
public:
    CRefineNurbsCommandImplementation();
    ~CRefineNurbsCommandImplementation();

    CRhinoCommand::result RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand) override;
};

