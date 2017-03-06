/** @file GismoImportImplementation.h
@brief Provides an the implementation for reading G+SMO XML files

*/
#pragma once
#include "FileImportImplementation.h"

/**
@brief Provides an the implementation for reading G + SMO XML files

*/

class GISMO_SUPPORT_DLLEXPORT CGismoImportImplementation : IFileImportImplementation
{
public:
    CGismoImportImplementation();
    ~CGismoImportImplementation();

    BOOL ReadFile(const wchar_t* filename, int index, CRhinoDoc& doc, const CRhinoFileReadOptions& options, CRhinoPlugIn& plugIn) override;
    void AddFileType(ON_ClassArray<CRhinoFileType>& extensions, const CRhinoFileReadOptions& options, CRhinoPlugIn& plugIn) override;
};

