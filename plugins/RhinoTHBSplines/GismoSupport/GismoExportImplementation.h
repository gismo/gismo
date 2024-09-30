/** @file GismoExportImplementation.h
@brief Provides an the implementation for writing G+SMO XML files

*/
#pragma once
#include "FileExportImplementation.h"

/**
@brief Provides an the implementation for reading G + SMO XML files

*/

class GISMO_SUPPORT_DLLEXPORT CGismoExportImplementation : IFileExportImplementation
{
public:
    CGismoExportImplementation();
    ~CGismoExportImplementation();

    /// Read the file
    virtual BOOL WriteFile(const wchar_t* filename, int index, CRhinoDoc& doc, const CRhinoFileWriteOptions& options, CRhinoPlugIn& plugIn);

    /// Add a supported file type
    virtual void AddFileType(ON_ClassArray<CRhinoFileType>& extensions, const CRhinoFileWriteOptions& options, CRhinoPlugIn& plugIn);

};

