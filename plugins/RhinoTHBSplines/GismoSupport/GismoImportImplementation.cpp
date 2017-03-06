#include "stdafx.h"
#include "GismoImportImplementation.h"
#include "ThbSurfaceUtils.h"
#include "ThbSurfaceObject.h"

CGismoImportImplementation::CGismoImportImplementation()
{
}


CGismoImportImplementation::~CGismoImportImplementation()
{
}


void CGismoImportImplementation::AddFileType(ON_ClassArray<CRhinoFileType>& extensions, const CRhinoFileReadOptions& options, CRhinoPlugIn& plugIn)
{
    CRhinoFileType ft(plugIn.PlugInID(), L"G+SMO Xml files (*.xml)", L"xml");
    extensions.Append(ft); 
}

BOOL CGismoImportImplementation::ReadFile(const wchar_t* filename, int index, CRhinoDoc& doc, const CRhinoFileReadOptions& options, CRhinoPlugIn& plugIn)
{
    gsFileData<> fd;

    size_t len = wcslen(filename);
    char* chars = new char[len + 1];
    wcstombs(chars, filename, len + 1);

    std::string fn(chars);

    fd.read(fn);

    std::vector<gsTHBSpline<2>*> thbs = fd.getAll<gsTHBSpline<2> >();
    for (auto it = thbs.begin(); it != thbs.end(); ++it)
    {
        gsTHBSpline<2>* thb = *it;
        auto obj = new CThbSurfaceObject();
    gsTHBSpline2* temp = new gsTHBSpline2(*thb);
        obj->SetHierarchicalSurface(temp);

        if (!doc.AddObject(obj))
        {
            delete obj;
        }
        obj = nullptr;
    }

    std::vector<gsTensorBSpline2*> nurbs = fd.getAll<gsTensorBSpline2 >();
    for (auto it = nurbs.begin(); it != nurbs.end(); ++it)
    {
        gsTensorBSpline2* nurb = *it;
        ON_NurbsSurface ns;
        ON_GismoUtils::NurbForm(*nurb, ns);
        delete nurb;
        ON_wString err;
        ON_TextLog log(err);
        if (!ns.IsValid(&log))
        {
            RhinoApp().Print(err);
            continue;
        }
        doc.AddSurfaceObject(ns);    
    }
    std::vector<gsTensorNurbs<2>*> rnurbs = fd.getAll<gsTensorNurbs<2> >();
    for (auto it = rnurbs.begin(); it != rnurbs.end(); ++it)
    {
        gsTensorNurbs<2>* rnurb = *it;
        ON_NurbsSurface ns;
        ON_GismoUtils::NurbForm(*rnurb, ns);
        delete rnurb;
        ON_wString err;
        ON_TextLog log(err);
        if (!ns.IsValid(&log))
        {
            RhinoApp().Print(err);
            continue;
        }
        doc.AddSurfaceObject(ns);
    }

    //std::vector<gsMultiPatch<>* > multi = fd.getAll<gsMultiPatch<> >();
    //for (auto it = multi.begin(); it != multi.end(); ++it)
    //{
    //    gsMultiPatch<>* p = *it;
    //    for (auto git = p->begin(); git != p->end(); ++it)
    //    {
    //        gsGeometry<>* g = *git;
    //        gsTensorNurbs2* tn = dynamic_cast<gsTensorNurbs2*>(g);
    //        if (tn != nullptr)
    //        {
    //            ON_NurbsSurface ns;
    //            ON_GismoUtils::NurbForm(*tn, ns);
    //            ON_wString err;
    //            ON_TextLog log(err);
    //            if (!ns.IsValid(&log))
    //            {
    //                RhinoApp().Print(err);
    //                continue;
    //            }
    //            doc.AddSurfaceObject(ns);
    //        }
    //        gsTensorBSpline2* tb = dynamic_cast<gsTensorBSpline2*>(g);
    //        if (tb != nullptr)
    //        {
    //            ON_NurbsSurface ns;
    //            ON_GismoUtils::NurbForm(*tb, ns);
    //            ON_wString err;
    //            ON_TextLog log(err);
    //            if (!ns.IsValid(&log))
    //            {
    //                RhinoApp().Print(err);
    //                continue;
    //            }
    //            doc.AddSurfaceObject(ns);
    //        }
    //        gsTHBSpline2* thb = dynamic_cast<gsTHBSpline2*>(g);
    //        if (thb != nullptr)
    //        {
    //            auto obj = new CThbSurfaceObject();
    //            obj->SetHierarchicalSurface(thb);
    //            if (!doc.AddObject(obj))
    //            {
    //                delete obj;
    //            }
    //        }
    //    }
    //    delete p;
    //}// end iterate multi-patch

    std::vector<gsBSpline<>*> curves = fd.getAll<gsBSpline<> >();
    for (auto it = curves.begin(); it != curves.end(); ++it)
    {
        gsBSpline<>* c = *it;
        ON_NurbsCurve crv;
        ON_GismoUtils::NurbForm(*c, crv);
        doc.AddCurveObject(crv);
        delete c;
    }
    // PlanarDomain

    //CThbSurfaceObject* obj(nullptr);
    //if (fd.has< gsTHBSpline2 >())
    //{
    //    gsTHBSpline2* thb = fd.getFirst<gsTHBSpline2 >();
    //    obj = new CThbSurfaceObject();
    //    obj->SetHierarchicalSurface(thb);

    //    if (!doc.AddObject(obj))
    //    {
    //        delete obj;            
    //    }
    //}
    // todo: more support for other geometry types
    // NurbsCurve, NurbsSurface, Mesh, Point, etc.


    delete[] chars;
    return true;
}