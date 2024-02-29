/** @file gsReadBrep.cpp

    @brief Reading OpenCascade .brep files

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsOpenCascade/gsWriteOcct.h>
#include <gsCore/gsSurface.h>
#include <gsCore/gsMultiPatch.h>
#include <gsNurbs/gsBSplineBasis.h>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>

#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <AIS_TexturedShape.hxx>
#include <BRepTools.hxx>
#include <Graphic3d_Texture2D.hxx>
#include <BRep_Tool.hxx>

#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Face.hxx>
#include <V3d_DirectionalLight.hxx>

#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <BRepLib.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <ShapeBuild_ReShape.hxx>
#include <BRepCheck_Analyzer.hxx>

#include "ShapeFix_Shape.hxx"
#include "ShapeFix_Wireframe.hxx"
#include "ShapeFix.hxx"
#include "ShapeFix_FixSmallFace.hxx"
#include "ShapeFix_Wireframe.hxx"

#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRepMesh.hxx>
#include <BRepMesh_IncrementalMesh.hxx>

#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_NurbsConvert.hxx>
#include "BRepBuilderAPI_MakeVertex.hxx"
#include "BRepBuilderAPI_MakeShell.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include <BRepBuilderAPI_MakeFace.hxx>

#include "BRepOffsetAPI_Sewing.hxx"

#include <GeomAdaptor_Curve.hxx>
#include <Geom_Line.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomConvert.hxx>

#include <Geom2dAdaptor_Curve.hxx>
#include <Geom2dConvert.hxx>
#include <Geom2d_TrimmedCurve.hxx>
#include <Geom2d_Line.hxx>

#include <IGESControl_Controller.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>
#include <IGESData_GlobalSection.hxx>
#include <IGESData_IGESModel.hxx>

#include <STEPControl_Controller.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>

#include <StlAPI_Writer.hxx>
#include <VrmlAPI_Writer.hxx>

#include <Interface_Static.hxx>

namespace gismo {

namespace extensions {



Handle(Geom_Surface) makeOcctGeom(const gsSurface<real_t> & srf)
{
    TColgp_Array2OfPnt cps(1,srf.basis().component(0).size(),
                           1,srf.basis().component(1).size() );

    const gsKnotVector<real_t> & kv1 = 
        dynamic_cast<const gsBSplineBasis<real_t>&>( srf.basis().component(0) ).knots();
    const gsKnotVector<real_t> & kv2 = 
        dynamic_cast<const gsBSplineBasis<real_t>&>( srf.basis().component(1) ).knots();
    TColStd_Array1OfReal UKnots(1,kv1.uSize()), VKnots(1,kv2.uSize());
    TColStd_Array1OfInteger UMults(1,kv1.uSize()), VMults(1,kv2.uSize());

    Standard_Integer theIndex = 1;
    for(auto  i = kv1.ubegin(); i!=kv1.uend(); ++i)
    {
        UMults.SetValue(theIndex, i.multiplicity());
        UKnots.SetValue(theIndex++, i.value());
    }
    theIndex = 1;
    for(auto  i = kv2.ubegin(); i!=kv2.uend(); ++i)
    {
        VMults.SetValue(theIndex, i.multiplicity());
        VKnots.SetValue(theIndex++, i.value());
    }
    theIndex = 0;
    const bool dd = ( 3 == srf.targetDim() );
    for(index_t  j = 0; j!=srf.basis().component(1).size(); ++j)
        for(index_t  i = 0; i!=srf.basis().component(0).size(); ++i)
        {
            cps.SetValue(i+1,j+1, gp_Pnt( srf.coef(theIndex,0),
                                          srf.coef(theIndex,1),
                                          dd ? srf.coef(theIndex,2) : 0.0));
            ++theIndex;
        }

    Handle(Geom_Surface) bsp;
    if (srf.basis().isRational() )
    {
        TColStd_Array2OfReal wgt(1,srf.basis().component(0).size(),
                                   1,srf.basis().component(1).size() );
        theIndex = 0;
        for(index_t  j = 0; j!=srf.basis().component(1).size(); ++j)
            for(index_t  i = 0; i!=srf.basis().component(0).size(); ++i)
                wgt.SetValue(i+1,j+1, srf.basis().weights().at(theIndex++) );

        bsp.reset(new Geom_BSplineSurface(cps, wgt,  UKnots, VKnots, UMults, VMults,
                                          srf.basis().degree(0), srf.basis().degree(1) ) );
    }
    else
        bsp.reset(new Geom_BSplineSurface(cps, UKnots, VKnots, UMults, VMults,
                                          srf.basis().degree(0), srf.basis().degree(1) ) );
return bsp;
}

bool writeOcctIges(const gsSurface<real_t> & srf, const std::string & name)
{
    IGESControl_Writer writer;
    // Interface_Static::SetCVal("write.iges.header.product", "G+Smo");
    // Interface_Static::SetCVal("write.iges.header.receiver", "G+SmoFriend");
    // Interface_Static::SetCVal("write.iges.header.author", "G+SmoUser");
    // Interface_Static::SetCVal("write.iges.header.company", "Inria");
    IGESData_GlobalSection header = writer.Model()->GlobalSection();
    header.SetAuthorName(new TCollection_HAsciiString("G+Smo"));
    header.SetCompanyName(new TCollection_HAsciiString("Inria"));
    header.SetSendName(new TCollection_HAsciiString("G+SmoSender"));
    header.SetSystemId(new TCollection_HAsciiString("G+Smo"));
    header.SetUnitName(new TCollection_HAsciiString("G+Smo"));
    header.SetReceiveName(new TCollection_HAsciiString("G+Smo"));
    header.SetFileName(new TCollection_HAsciiString(name.c_str()));
    writer.Model()->SetGlobalSection(header);                
    writer.AddGeom( makeOcctGeom(srf) );
    const std::string fname = name + ".igs";
    return writer.Write(fname.c_str());;
}

bool writeOcctIgesMp(const gsMultiPatch<real_t> & mp, const std::string & name)
{
    IGESControl_Writer writer;
    // Interface_Static::SetCVal("write.iges.header.product", "G+Smo");
    // Interface_Static::SetCVal("write.iges.header.receiver", "G+SmoFriend");
    // Interface_Static::SetCVal("write.iges.header.author", "G+SmoUser");
    // Interface_Static::SetCVal("write.iges.header.company", "Inria");
    IGESData_GlobalSection header = writer.Model()->GlobalSection();
    header.SetAuthorName(new TCollection_HAsciiString("G+Smo"));
    header.SetCompanyName(new TCollection_HAsciiString("Inria"));
    header.SetSendName(new TCollection_HAsciiString("G+SmoSender"));
    header.SetSystemId(new TCollection_HAsciiString("G+Smo"));
    header.SetUnitName(new TCollection_HAsciiString("G+Smo"));
    header.SetReceiveName(new TCollection_HAsciiString("G+Smo"));
    header.SetFileName(new TCollection_HAsciiString(name.c_str()));
    writer.Model()->SetGlobalSection(header);                

    for (size_t i = 0; i!= mp.nPatches(); ++i)
    {
        writer.AddGeom( makeOcctGeom( (const gsSurface<real_t>&)mp.patch(i) ) );
    }

    const std::string fname = name + ".igs";
    return writer.Write(fname.c_str());;
}

bool writeOcctStep(const gsSurface<real_t> & srf, const std::string & name)
{
    STEPControl_Writer writer;
    // Interface_Static::SetCVal("xstep.cascade.unit", unit.c_str());
    // Interface_Static::SetCVal("write.step.unit", unit.c_str());
    // Interface_Static::SetIVal("write.step.nonmanifold", 1);
    TopoDS_Face face = BRepBuilderAPI_MakeFace(makeOcctGeom(srf),1e-16);
    IFSelect_ReturnStatus status = writer.Transfer(face, STEPControl_AsIs);
    GISMO_ENSURE(status, "Error transferring shape to STEP.");
    const std::string fname = name + ".stp";
    return writer.Write(fname.c_str());
}

}// namespace extensions

}// namespace gismo

