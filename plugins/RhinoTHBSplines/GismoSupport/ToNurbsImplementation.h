/** @file ToNurbsImplementation.h
@brief Provides the implementation for the thbToNurbs command

*/
#pragma once
#include "CommandImplementation.h"

class GISMO_SUPPORT_DLLEXPORT CToNurbsImplementation : public ICommandImplementation
{
public:
    CToNurbsImplementation();
    ~CToNurbsImplementation();

    CRhinoCommand::result RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand) override;

};

 