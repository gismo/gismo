/** @file FileImportImplementation.h
@brief Provides an interface for classes providing the implementation of file import.

*/
#pragma once

#ifdef DO_EXPORTS
#define GISMO_SUPPORT_DLLEXPORT __declspec(dllexport)
#else
#define GISMO_SUPPORT_DLLEXPORT __declspec(dllimport)
#endif

/** @brief interface (only pure virtual) for Rhino file immport implementation.
Don't forget to tag your class with GISMO_SUPPORT_DLLEXPORT.
*/
class GISMO_SUPPORT_DLLEXPORT IFileImportImplementation
{
public:

    /// Read the file
    virtual BOOL ReadFile(const wchar_t* filename, int index, CRhinoDoc& doc, const CRhinoFileReadOptions& options, CRhinoPlugIn& plugIn) = 0;
    
    /// Add a supported file type
    virtual void AddFileType(ON_ClassArray<CRhinoFileType>& extensions, const CRhinoFileReadOptions& options, CRhinoPlugIn& plugIn) = 0;
};
