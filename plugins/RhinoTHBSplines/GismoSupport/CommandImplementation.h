/** @file CommandImplementation.h
@brief Provides an interface for classes providing the implementation of a Rhino command.

*/
#pragma once

#ifdef DO_EXPORTS
#define GISMO_SUPPORT_DLLEXPORT __declspec(dllexport)
#else
#define GISMO_SUPPORT_DLLEXPORT __declspec(dllimport)
#endif

/** @brief interface (only pure virtual) for Rhino command implementation. 
           Don't forget to tag your class with GISMO_SUPPORT_DLLEXPORT.
*/
class GISMO_SUPPORT_DLLEXPORT ICommandImplementation
{
public:

    /// Run the actual command. Inheritors need to implement this.
    virtual CRhinoCommand::result RunActualCommand(const CRhinoCommandContext& context, CRhinoCommand& callingCommand) = 0;
};
