#include "stdafx.h"
#include "GismoExportImplementation.h"
#include "RhinoGetThbSurfaceObject.h"
#include "ThbSurfaceObject.h"
#include "ThbSurfaceUserData.h"

CGismoExportImplementation::CGismoExportImplementation()
{
}


CGismoExportImplementation::~CGismoExportImplementation()
{
}

/// Read the file
BOOL CGismoExportImplementation::WriteFile(const wchar_t* filename, int index, CRhinoDoc& doc, const CRhinoFileWriteOptions& options, CRhinoPlugIn& plugIn)
{
    BOOL rc = FALSE;

    CRhinoGetThbSurfaceObject go;
    go.SetCommandPrompt(L"Select THB spline surface(s)");
    CRhinoGet::result gr = go.GetObjects(1, 0);
    if (gr != CRhinoGet::object)
        return rc;

    std::wstring ws = std::wstring(filename);
    std::string fn(ws.begin(), ws.end());
    gsFileData<> fd;

    for (int i = 0; i < go.ObjectCount(); ++i)
    {
        CRhinoObjRef ref = go.Object(i);
        const CRhinoObject* obj = ref.Object();
        const CThbSurfaceObject* thb = CThbSurfaceObject::Cast(obj);
        if (thb == nullptr) continue; // should not happen

        CThbSurfaceUserData* ud = CThbSurfaceUserData::Cast(thb->GetGeometryUserData(CThbSurfaceUserData::Id()));
        if (ud == nullptr) continue; // shold not happen either

        gsTHBSpline<2>* s = ud->HierarchicalSurface();
        if (s == nullptr) continue;

        fd << *s;
    }
    
    fd.save(fn, index == 1);
    return TRUE;
}

/// Add a supported file type
void CGismoExportImplementation::AddFileType(ON_ClassArray<CRhinoFileType>& extensions, const CRhinoFileWriteOptions& options, CRhinoPlugIn& plugIn)
{
    CRhinoFileType ftXml(plugIn.PlugInID(), L"G+SMO Xml files (*.xml)", L"xml");
    extensions.Append(ftXml);
    //this gives double extensions somehow. 
    //CRhinoFileType ftXmlGz(plugIn.PlugInID(), L"G+SMO Compressed Xml Files (*.xml.gz)", L"xml.gz");
    //extensions.Append(ftXmlGz);
}
