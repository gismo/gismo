/** @file gsReadBrep.cpp

    @brief Reading OpenCascade .brep files

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsOpenCascade/gsReadBrep.h>


#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <AIS_TexturedShape.hxx>
#include <BRepTools.hxx>
#include <Graphic3d_Texture2D.hxx>
#include <BRep_Tool.hxx>
#include <TopoDS.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Face.hxx>
#include <V3d_DirectionalLight.hxx>


#include "BRepTools.hxx"
#include "BRep_Tool.hxx"
#include "BRep_Builder.hxx"
#include "BRepLib.hxx"
#include "BRepGProp.hxx"
#include "GProp_GProps.hxx"

#include "ShapeBuild_ReShape.hxx"
#include "BRepCheck_Analyzer.hxx"

#include "ShapeFix_Shape.hxx"
#include "ShapeFix_Wireframe.hxx"
#include "ShapeFix.hxx"
#include "ShapeFix_FixSmallFace.hxx"
#include "ShapeFix_Wireframe.hxx"

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <TopExp_Explorer.hxx>


#include <IGESControl_Controller.hxx>
#include <IGESControl_Reader.hxx>
#include <IGESControl_Writer.hxx>

#include <STEPControl_Controller.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>

#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <BRepMesh.hxx>
#include <BRepMesh_IncrementalMesh.hxx>


#include <StlAPI_Writer.hxx>
#include <VrmlAPI_Writer.hxx>


#include "BRepBuilderAPI_MakeVertex.hxx"
#include "BRepBuilderAPI_MakeShell.hxx"
#include "BRepBuilderAPI_MakeSolid.hxx"
#include "BRepOffsetAPI_Sewing.hxx"



#include <Geom_BSplineSurface.hxx>



#include <GeomAdaptor_Curve.hxx>
#include <Geom2dAdaptor_Curve.hxx>

#include <Geom_Line.hxx>


#include <GeomConvert.hxx>
#include <Geom2dConvert.hxx>

#include <Geom2d_TrimmedCurve.hxx>

#include <Geom2d_Line.hxx>


#include <BRepBuilderAPI_NurbsConvert.hxx>


namespace gismo {

namespace extensions {


bool gsReadBrep( const char * filename, internal::gsXmlTree & data)
{
    // Read
    BRep_Builder Builder;
    TopoDS_Shape shape;
    BRepTools::Read(shape, filename, Builder);

    return readTopoDS_Shape(shape, data);
}

bool readTopoDS_Shape( const TopoDS_Shape & inputShape, internal::gsXmlTree & data  )
{
    // Note: some edges used to be lines, this makes them all b-splines
    BRepBuilderAPI_NurbsConvert conv;
    conv.Perform(inputShape);
    const TopoDS_Shape & shape = conv.Shape();
    //healGeometry(1e-6,true,true,true,false,false);
    BRepTools::Clean(shape);


    //-------------- Phase I -- make map

    TopTools_IndexedMapOfShape solids_map, shells_map, faces_map, wires_map,  edges_map, vertices_map;
    TopExp_Explorer exp_solids, exp_shells, exp_faces, exp_wires, exp_edges, exp_vertices;
    
    for (exp_solids.Init(shape, TopAbs_SOLID); exp_solids.More(); exp_solids.Next())
    {
        TopoDS_Solid solid = TopoDS::Solid(exp_solids.Current());
        if (solids_map.FindIndex(solid) < 1)
        {
            solids_map.Add(TopoDS::Solid(exp_solids.Current()));
            for (exp_shells.Init(exp_solids.Current(), TopAbs_SHELL); exp_shells.More(); exp_shells.Next())
            {
                TopoDS_Shell shell = TopoDS::Shell(exp_shells.Current().Composed(solid.Orientation()));
                if (shells_map.FindIndex(shell) < 1)
                {
                    shells_map.Add(shell);
                    for (exp_faces.Init(shell, TopAbs_FACE); exp_faces.More(); exp_faces.Next())
                    {
                        TopoDS_Face face = TopoDS::Face(exp_faces.Current().Composed(shell.Orientation()));
                        if (faces_map.FindIndex(face) < 1)
                        {
                            faces_map.Add(face);
                            for (exp_wires.Init(exp_faces.Current(), TopAbs_WIRE); exp_wires.More(); exp_wires.Next())
                            {
                                TopoDS_Wire wire = TopoDS::Wire(exp_wires.Current().Composed(face.Orientation()));
                                if (wires_map.FindIndex(wire) < 1)
                                {
                                    wires_map.Add(wire);
                                    for (exp_edges.Init(exp_wires.Current(), TopAbs_EDGE); exp_edges.More(); exp_edges.Next())
                                    {
                                        TopoDS_Edge edge = TopoDS::Edge(exp_edges.Current().Composed(wire.Orientation()));
                                        if (edges_map.FindIndex(edge) < 1)
                                        {
                                            edges_map.Add(edge);
                                            for (exp_vertices.Init(exp_edges.Current(), TopAbs_VERTEX); exp_vertices.More(); exp_vertices.Next())
                                            {
                                                TopoDS_Vertex vertex = TopoDS::Vertex(exp_vertices.Current().Composed(edge.Orientation()));
                                                if (vertices_map.FindIndex(vertex) < 1)
                                                    vertices_map.Add(vertex);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    
    /*

    // Lone shells
    for (exp_shells.Init(shape, TopAbs_SHELL, TopAbs_SOLID); exp_shells.More(); exp_shells.Next())
    {
    TopoDS_Shape shell = exp_shells.Current().Composed(exp_solids.Current().Orientation());
    if (shells_map.FindIndex(shell) < 1)
    {
    shells_map.Add(shell);
    for (exp_faces.Init(shell, TopAbs_FACE); exp_faces.More(); exp_faces.Next())
    {
    TopoDS_Face face = TopoDS::Face(exp_faces.Current().Composed(shell.Orientation()));
    if (faces_map.FindIndex(face) < 1)
    {
    faces_map.Add(face);
    for (exp_wires.Init(exp_faces.Current(), TopAbs_WIRE); exp_wires.More(); exp_wires.Next())
    {
    TopoDS_Wire wire = TopoDS::Wire(exp_wires.Current().Composed(face.Orientation()));
    if (wires_map.FindIndex(wire) < 1)
    {
    wires_map.Add(wire);
    for (exp_edges.Init(exp_wires.Current(), TopAbs_EDGE); exp_edges.More(); exp_edges.Next())
    {
    TopoDS_Edge edge = TopoDS::Edge(exp_edges.Current().Composed(wire.Orientation()));
    if (edges_map.FindIndex(edge) < 1)
    {
    edges_map.Add(edge);
    for (exp_vertices.Init(exp_edges.Current(), TopAbs_VERTEX); exp_vertices.More(); exp_vertices.Next())
    {
    TopoDS_Vertex vertex = TopoDS::Vertex(exp_vertices.Current().Composed(edge.Orientation()));
    if (vertices_map.FindIndex(vertex) < 1)
    vertices_map.Add(vertex);
    }
    }
    }
    }
    }
    }
    }
    }
    }

    // Lone faces
    for (exp_faces.Init(shape, TopAbs_FACE, TopAbs_SHELL); exp_faces.More(); exp_faces.Next())
    {
    TopoDS_Face face = TopoDS::Face(exp_faces.Current().Composed(exp_faces.Current().Orientation()));
    if (faces_map.FindIndex(face) < 1)
    {
    faces_map.Add(face);
    for (exp_wires.Init(exp_faces.Current(), TopAbs_WIRE); exp_wires.More(); exp_wires.Next())
    {
    TopoDS_Wire wire = TopoDS::Wire(exp_wires.Current().Composed(face.Orientation()));
    if (wires_map.FindIndex(wire) < 1)
    {
    wires_map.Add(wire);
    for (exp_edges.Init(exp_wires.Current(), TopAbs_EDGE); exp_edges.More(); exp_edges.Next())
    {
    TopoDS_Edge edge = TopoDS::Edge(exp_edges.Current().Composed(wire.Orientation()));
    if (edges_map.FindIndex(edge) < 1)
    {
    edges_map.Add(edge);
    for (exp_vertices.Init(exp_edges.Current(), TopAbs_VERTEX); exp_vertices.More(); exp_vertices.Next())
    {
    TopoDS_Vertex vertex = TopoDS::Vertex(exp_vertices.Current().Composed(edge.Orientation()));
    if (vertices_map.FindIndex(vertex) < 1)
    vertices_map.Add(vertex);
    }
    }
    }
    }
    }
    }
    }

    // Free Wires
    for (exp_wires.Init(shape, TopAbs_WIRE, TopAbs_FACE); exp_wires.More(); exp_wires.Next())
    {
    TopoDS_Wire wire = TopoDS::Wire(exp_wires.Current().Composed(exp_wires.Current().Orientation()));
    if (wires_map.FindIndex(wire) < 1)
    {
    wires_map.Add(wire);
    for (exp_edges.Init(exp_wires.Current(), TopAbs_EDGE); exp_edges.More(); exp_edges.Next())
    {
    TopoDS_Edge edge = TopoDS::Edge(exp_edges.Current().Composed(wire.Orientation()));
    if (edges_map.FindIndex(edge) < 1)
    {
    edges_map.Add(edge);
    for (exp_vertices.Init(exp_edges.Current(), TopAbs_VERTEX); exp_vertices.More(); exp_vertices.Next())
    {
    TopoDS_Vertex vertex = TopoDS::Vertex(exp_vertices.Current().Composed(edge.Orientation()));
    if (vertices_map.FindIndex(vertex) < 1)
    vertices_map.Add(vertex);
    }
    }
    }
    }
    }

    // Lone edges
    for (exp_edges.Init(shape, TopAbs_EDGE, TopAbs_WIRE); exp_edges.More(); exp_edges.Next())
    {
    TopoDS_Edge edge = TopoDS::Edge(exp_edges.Current().Composed(exp_edges.Current().Orientation()));
    if (edges_map.FindIndex(edge) < 1)
    {
    edges_map.Add(edge);
    for (exp_vertices.Init(exp_edges.Current(), TopAbs_VERTEX); exp_vertices.More(); exp_vertices.Next())
    {
    TopoDS_Vertex vertex = TopoDS::Vertex(exp_vertices.Current().Composed(edge.Orientation()));
    if (vertices_map.FindIndex(vertex) < 1)
    vertices_map.Add(vertex);
    }
    }
    }

    // Lone vertices
    for (exp_vertices.Init(shape, TopAbs_VERTEX, TopAbs_EDGE); exp_vertices.More(); exp_vertices.Next())
    {
    TopoDS_Vertex vertex = TopoDS::Vertex(exp_vertices.Current());
    if (vertices_map.FindIndex(vertex) < 1)
    vertices_map.Add(vertex);
    }

    */


    //-------------- Phase II -- Read in data
    
    int nvertices = vertices_map.Extent();

    for (int i = 1; i <= nvertices; ++i)   // starts at 1 !!
    {
        TopoDS_Vertex v = TopoDS::Vertex(vertices_map(i));
        gp_Pnt pt = BRep_Tool::Pnt(v);
        std::cout<<"Vertex "<<i <<": ("<< pt.X() <<", "<< pt.Y() <<", "<< pt.Z() <<" )\n";
    }
    std::cout << nvertices << " Vertices imported" << std::endl;
    std::cout << "----------------------------" << std::endl;
  
    int nedges = edges_map.Extent();

    for (int i = 1; i <= nedges; i++)
    {
        TopoDS_Edge edgo= TopoDS::Edge(edges_map(i));
        int from = vertices_map.FindIndex(TopExp::FirstVertex(edgo));
        int to = vertices_map.FindIndex(TopExp::LastVertex(edgo));
        Standard_Real aFirst, aLast;
        Handle(Geom_Curve) aCurve3d = BRep_Tool::Curve(edgo, aFirst, aLast);
        //std::cout << fromcovn << " " << i << " " << from << " " << aFirst << std::endl;
//    ICoVertexN fromcov(*this,fromcovn,from,i,aFirst);
//    iterator_covertices fromcov_ptr= AddCoVertex(fromcov).first;
//    std::cout << fromcov_ptr->num << " " << fromcov_ptr->GetEdge()->num << " " << fromcov_ptr->GetVertex()->num << " " << fromcov_ptr->U() << std::endl;
//        int tocovn=numcov++;
        //std::cout << tocovn << " " << i << " " << to << " " << aLast << std::endl;
//    ICoVertexN tocov(*this,tocovn,to,i,aLast);
//    iterator_covertices tocov_ptr=AddCoVertex(tocov).first;
//    std::cout << tocov_ptr->num << " " << tocov_ptr->GetEdge()->num << " " << tocov_ptr->GetVertex()->num << " " << tocov_ptr->U() << std::endl;
        bool dege=BRep_Tool::Degenerated(edgo);
//    iterator_curves crv=end_curves();
//    if (!dege)
//      crv=AddCurve(ICurveOCC(ncrv,aCurve3d,BRep_Tool::Degenerated(edgo))).first;
//    IPoint pst,pend;
        std::cout << "* edge " << i << " fromv=" << from<<  " tov=" << to << " [start,end]= [" << aFirst << "," << aLast<<"]"<<std::endl;
//    edg->SetFromCoVertex(fromcov_ptr);
//    edg->SetToCoVertex(tocov_ptr);
//    if (!dege) edg->SetSupport(crv);
//    IEdge *edg=AddEdge(iedgo).first;
//    std::cout << "edge " << iedgo.num << " " << edg->num << std::endl;
    }
    std::cout << nedges << " Edges imported" << std::endl;
    std::cout << "----------------------------" << std::endl;
  
    int nfaces = faces_map.Extent();
    int nwires=0;
    int nsurfaces=0;

    for (int i = 1; i <= nfaces; i++)
    {
        std::cout << "** Face " << i <<"\n";
      
        TopoDS_Face faceo=TopoDS::Face(faces_map(i));
        Handle(Geom_Surface) aSurface = BRep_Tool::Surface(faceo);

        Handle(Geom_BSplineSurface) bspSurface = Handle(Geom_BSplineSurface)::DownCast(aSurface);
        assert( !bspSurface.IsNull() );
        std::cout <<"degree1="<< bspSurface->UDegree() <<", degree2="<< bspSurface->VDegree() <<"\n";
        std::cout <<"periodic1="<< bspSurface->IsUPeriodic() <<", periodic2="<< bspSurface->IsVPeriodic() <<"\n";
    
        const TColStd_Array1OfReal & k1 = bspSurface->UKnots();
        const TColStd_Array1OfReal & k2 = bspSurface->VKnots();
        const TColStd_Array1OfInteger & mult1 = bspSurface->UMultiplicities();
        const TColStd_Array1OfInteger & mult2 = bspSurface->VMultiplicities();    
        std::cout <<"knots1="<< k1.Size() <<" knots2="<< k2.Size() <<"\n";
        std::copy(k1.begin(), k1.end(), std::ostream_iterator<double>(std::cout, " ") );std::cout<<std::endl;
        std::copy(mult1.begin(), mult1.end(), std::ostream_iterator<double>(std::cout, " ") );std::cout<<std::endl;
        std::copy(k2.begin(), k2.end(), std::ostream_iterator<double>(std::cout, " ") );std::cout<<std::endl;
        std::copy(mult2.begin(), mult2.end(), std::ostream_iterator<double>(std::cout, " ") );std::cout<<std::endl;
    
        const TColStd_Array2OfReal & w = *bspSurface->Weights();
        const TColgp_Array2OfPnt   & cf = bspSurface->Poles();
        std::cout <<"coefs1="<< bspSurface->NbUPoles() <<" coefs2="<< bspSurface->NbVPoles()<<" (i.e. "<<cf.Size() <<" ceofs.)\n";

        // loop over coefficients
        /* 
           int aa = 0;
           for (int r = cf.LowerRow(); r <= cf.UpperRow(); r++)
           for(int c = cf.LowerCol(); c <= cf.UpperCol(); c++)
           {
           ++aa;
           const gp_Pnt & q = cf(r,c);
           }
        */
        std::cout<<std::endl;

        nsurfaces++;
        TopExp_Explorer exp_wires,exp_edges;
        for (exp_wires.Init(faceo, TopAbs_WIRE); exp_wires.More(); exp_wires.Next())
        {
            TopoDS_Wire wire = TopoDS::Wire(exp_wires.Current().Composed(faceo.Orientation()));
            nwires++;
            std::cout <<"loop\n";
            for (exp_edges.Init(wire, TopAbs_EDGE); exp_edges.More(); exp_edges.Next())
            {
                Standard_Real aFirst, aLast;
                TopoDS_Edge edgeo = TopoDS::Edge(exp_edges.Current().Composed(wire.Orientation()));
                Handle(Geom2d_Curve) aCurve2d = BRep_Tool::CurveOnSurface(edgeo,faceo, aFirst, aLast);
                assert( !aCurve2d.IsNull() );
/*          
//https://www.opencascade.com/doc/occt-6.9.1/refman/html/_geom_abs___curve_type_8hxx.html#af25c179d5cabd33fddebe5a0dc96971c
Geom2dAdaptor_Curve ad (aCurve2d);
GeomAbs_CurveType adtype = ad.GetType();
std::cout <<" "<< adtype<<"\n";
// 0: GeomAbs_Line
// 6: GeomAbs_BSplineCurve
if ( 0 == adtype )
{
Handle(Geom2d_Line) l2d = Handle(Geom2d_Line)::DownCast(aCurve2d);
assert( !l2d.IsNull() );
}
*/
    int numedg=edges_map.FindIndex(edgeo);

/*
  if (BRep_Tool::Degenerated(edgeo))
  {
  iterator_edges ptredg=GetEdge(numedg);
  IBrep::iterator_covertices from=ptredg->GetFromCoVertex();
  IBrep::iterator_covertices to=ptredg->GetToCoVertex();
  from->U() =aFirst;
  to->U() =aLast;
  }
*/
//        ICurveOCC coo(numcrv,aCurve2d,BRep_Tool::Degenerated(edgeo));
//        AddCurve(coo);
//        IPoint pst,pend;
//        coo.P(aFirst,pst);
//        coo.P(aLast,pend);
    std::cout << "   edge " << numedg << ", orif=" << faceo.Orientation() << ", oriw=" << wire.Orientation() << ", orie= " << edgeo.Orientation() <<  ", param=[" << aFirst
        //<< " " << pst[0] << " " << pst[1] << " " << pst[2] << " "
              << ", " <<  aLast <<"] "
        //<< " " << pend[0] << " " << pend[1] << " " << pend[2]
              <<std::endl;
//        coo.SetEdge();
//        ICoEdgeN coe(*this,ncoedges,numedg, (1-edgeo.Orientation() *2),i);
//        ptredg->AddCoEdge(ncoedges);
//        coe.SetSupport(numcrv);
//        IBrep::iterator_coedges coep= AddCoEdge(coe).first;
//        loop.push_back(coep);        

    // get the bspline trimming curve
    Handle(Geom2d_BSplineCurve) bsp2d =  Handle(Geom2d_BSplineCurve)::DownCast(aCurve2d);
    assert( !bsp2d.IsNull() );
    std::cout <<"degree="<< bsp2d->Degree() <<", periodic="<<bsp2d->IsPeriodic()<<", rational="<<bsp2d->IsRational()<<", numCoefs="<<bsp2d->NbPoles()<<", knots:\n";

    const TColStd_Array1OfReal & rKnots = bsp2d->KnotSequence();
    //const TColStd_Array1OfReal & uKnots = bsp2d->Knots();
    //const TColStd_Array1OfInteger & Kmult = bsp2d->Multiplicities();
    std::copy(rKnots.begin(), rKnots.end(), std::ostream_iterator<double>(std::cout, " ") );std::cout<<std::endl;

    const TColgp_Array1OfPnt2d & ccfs = bsp2d->Poles();
    const TColStd_Array1OfReal & cwgts = *bsp2d->Weights();

    for(int c = ccfs.Lower(); c <= ccfs.Upper(); c++)
    {
        const gp_Pnt2d & q = ccfs(c);
        std::cout<<" ["<<q.X() <<", "<< q.Y()<< "]";
    }
    std::cout<<std::endl;
                

            }
//      ifaceo.AddLoop(loop);
        }
//    AddFace(ifaceo).first;
    }
    std::cout << nfaces << " Faces imported" << std::endl;
    std::cout << nsurfaces << " Surfaces imported" << std::endl;
    std::cout << nwires << " Loops imported" << std::endl;
    std::cout << "----------------------------" << std::endl;

    int nvolumes = solids_map.Extent();
    int nshells=0;
    int ncofaces=0;

    for(int i = 1; i <= nvolumes; i++)
    {

        std::cout << "*** Volume " << i <<"\n";
            
        TopoDS_Solid solido=TopoDS::Solid(solids_map(i));
//    IVolumeOCC ivol(i,solido);
        TopExp_Explorer exp_shells,exp_faces;
//    iterator_volumes vol=AddVolume(ivol).first;
        for (exp_shells.Init(solido, TopAbs_SHELL); exp_shells.More(); exp_shells.Next())
        {
            TopoDS_Shell shell = TopoDS::Shell(exp_shells.Current().Composed(solido.Orientation()));
            nshells++;
//      IShellPtr ishell;
            for (exp_faces.Init(shell, TopAbs_FACE); exp_faces.More(); exp_faces.Next())
            {
                TopoDS_Face faceo = TopoDS::Face(exp_faces.Current().Composed(shell.Orientation()));

                int numface=faces_map.FindIndex(faceo);
                gsInfo<<" "<<numface;
//        iterator_faces ptrface=GetFace(numface);
//        ICoFaceN cof(*this,ncofaces,numface, (1-solido.Orientation() *2),i);
//        iterator_cofaces pcof=AddCoFace(cof).first;
//        ishell.push_back(pcof);
//        ptrface->AddCoFace(pcof);
                ncofaces++;
            }
//      vol->AddShell(ishell);
        }
        std::cout << std::endl;
    }
    std::cout << nvolumes << " Volumes imported" << std::endl;
    std::cout << nshells << " Shells imported" << std::endl;
    std::cout << ncofaces << " Cofaces imported" << std::endl;
    std::cout << "----------------------------" << std::endl;
//  UpdateAll();

    return true;
}






}// namespace extensions

}// namespace gismo

