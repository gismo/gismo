/** @file gsReadBrep.cpp

    @brief Reading OpenCascade .brep files

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsOpenCascade/gsReadBrep.h>


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

#include <STEPControl_Controller.hxx>
#include <STEPControl_Reader.hxx>
#include <STEPControl_Writer.hxx>

#include <StlAPI_Writer.hxx>
#include <VrmlAPI_Writer.hxx>

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
        gsInfo<<"Vertex "<<i <<": ("<< pt.X() <<", "<< pt.Y() <<", "<< pt.Z() <<" )\n";
    }
    gsInfo << nvertices << " Vertices imported" << gsEndl;
    gsInfo << "----------------------------" << gsEndl;
  
    int nedges = edges_map.Extent();

    for (int i = 1; i <= nedges; i++)
    {
        TopoDS_Edge edgo= TopoDS::Edge(edges_map(i));
        int from = vertices_map.FindIndex(TopExp::FirstVertex(edgo));
        int to = vertices_map.FindIndex(TopExp::LastVertex(edgo));
        Standard_Real aFirst, aLast;
        Handle(Geom_Curve) aCurve3d = BRep_Tool::Curve(edgo, aFirst, aLast);
        bool dege = BRep_Tool::Degenerated(edgo);
        gsInfo << "* edge " << i << " fromv=" << from<<  " tov=" << to << " [start,end]= ["
               << aFirst << "," << aLast<<"]" << (dege ? " !" : "")<<gsEndl;
    }
    gsInfo << nedges << " Edges imported" << gsEndl;
    gsInfo << "----------------------------" << gsEndl;
  
    int nfaces = faces_map.Extent();
    int nwires=0;
    int nsurfaces=0;

    for (int i = 1; i <= nfaces; i++)
    {
        gsInfo << "** Face " << i <<"\n";
      
        TopoDS_Face faceo=TopoDS::Face(faces_map(i));
        Handle(Geom_Surface) aSurface = BRep_Tool::Surface(faceo);
        
        readGeom_Surface(aSurface, data);

        nsurfaces++;
        TopExp_Explorer exp_wires,exp_edges;
        for (exp_wires.Init(faceo, TopAbs_WIRE); exp_wires.More(); exp_wires.Next())
        {
            TopoDS_Wire wire = TopoDS::Wire(exp_wires.Current().Composed(faceo.Orientation()));
            nwires++;
            gsInfo <<"loop\n";
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
                gsInfo <<" "<< adtype<<"\n";
                // 0: GeomAbs_Line
                // 6: GeomAbs_BSplineCurve
                if ( 0 == adtype )
                {
                Handle(Geom2d_Line) l2d = Handle(Geom2d_Line)::DownCast(aCurve2d);
                assert( !l2d.IsNull() );
                }
                */
                int numedg=edges_map.FindIndex(edgeo);

                gsInfo << "   edge " << numedg << ", orif=" << faceo.Orientation() << ", oriw=" << wire.Orientation() << ", orie= " << edgeo.Orientation() <<  ", param=[" << aFirst
                    //<< " " << pst[0] << " " << pst[1] << " " << pst[2] << " "
                       << ", " <<  aLast <<"] "
                    //<< " " << pend[0] << " " << pend[1] << " " << pend[2]
                       <<gsEndl;

                // get the bspline trimming curve
                readGeom2d_Curve(aCurve2d, data);
            }

            //###addloop
        }
        
        //###addface
        
    }
    gsInfo << nfaces << " Faces imported" << gsEndl;
    gsInfo << nsurfaces << " Surfaces imported" << gsEndl;
    gsInfo << nwires << " Loops imported" << gsEndl;
    gsInfo << "----------------------------" << gsEndl;

    int nvolumes = solids_map.Extent();
    int nshells=0;
    //int ncofaces=0;

    for(int i = 1; i <= nvolumes; i++)
    {

        gsInfo << "*** Volume " << i <<"\n";
            
        TopoDS_Solid solido=TopoDS::Solid(solids_map(i));
        TopExp_Explorer exp_shells,exp_faces;
        for (exp_shells.Init(solido, TopAbs_SHELL); exp_shells.More(); exp_shells.Next())
        {
            TopoDS_Shell shell = TopoDS::Shell(exp_shells.Current().Composed(solido.Orientation()));
            nshells++;
            for (exp_faces.Init(shell, TopAbs_FACE); exp_faces.More(); exp_faces.Next())
            {
                TopoDS_Face faceo = TopoDS::Face(exp_faces.Current().Composed(shell.Orientation()));
                int numface=faces_map.FindIndex(faceo);
                gsInfo<<" "<<numface;
                //ncofaces++;
            }
            //###addshell
        }
        gsInfo << gsEndl;
    }
    gsInfo << nvolumes << " Volumes imported" << gsEndl;
    gsInfo << nshells  << " Shells imported" << gsEndl;
    //gsInfo << ncofaces << " Cofaces imported" << gsEndl;
    gsInfo << "----------------------------" << gsEndl;

    return true;
}


bool readGeom_Surface(const Handle(Geom_Surface)& S, internal::gsXmlTree & data  )
{
    Handle(Geom_BSplineSurface) ds = Handle(Geom_BSplineSurface)::DownCast(S);
    if ( !ds.IsNull() ) return readGeom_BSplineSurface(ds, data);
    return false;
}


bool readGeom_BSplineSurface( const opencascade::handle<Geom_BSplineSurface> & S, internal::gsXmlTree & data  )
{
    internal::gsXmlNode* parent = data.first_node("xml") ;

    internal::gsXmlNode* g = internal::makeNode("Geometry", data);
    parent->append_node(g);
    internal::gsXmlNode* basis = internal::makeNode("Basis", data);
    internal::gsXmlNode* basis_rat = 0;

    if ( S->IsURational() || S->IsVRational() )
      {
        g->append_attribute( internal::makeAttribute("type", "TensorNurbs2", data) );
        basis_rat = internal::makeNode("Basis", data);
        basis_rat->append_attribute( internal::makeAttribute("type", "TensorNurbsBasis2", data) );
        basis = internal::makeNode("Basis", data);
        basis->append_attribute( internal::makeAttribute("type", "TensorBSplineBasis2", data) );
        basis_rat->append_node(basis);
        g->append_node(basis_rat);
      }
    else
      {
        g->append_attribute( internal::makeAttribute("type", "TensorBSpline2", data) );
        basis->append_attribute( internal::makeAttribute("type", "TensorBSplineBasis2", data) );
        g->append_node(basis);
      }

    if(S->IsUPeriodic())
        S->SetUNotPeriodic();
    if(S->IsVPeriodic())
        S->SetVNotPeriodic();

    int deg1 = S->UDegree();
    int deg2 = S->VDegree();
    const TColStd_Array1OfReal & knotSeq1 = S->UKnotSequence();
    const TColStd_Array1OfReal & knotSeq2 = S->VKnotSequence();
    int knot_count1 = knotSeq1.Length();
    int knot_count2 = knotSeq2.Length();

    for (int dir = 0; dir < 2; dir++ )
    {
        std::ostringstream tmp;
        tmp << std::setprecision(16);
        internal::gsXmlNode* b = internal::makeNode("Basis", data);
        b->append_attribute( internal::makeAttribute("type", "BSplineBasis", data) );
        b->append_attribute( internal::makeAttribute("index", dir, data) );
        basis->append_node(b);
        int deg;

        if(dir == 0)
        {
            deg = deg1;
            for(int k=1; k<=knot_count1; k++)
                tmp <<" "<< knotSeq1(k);

        }
        else
        {
            deg = deg2;
            for(int k=1; k<=knot_count2; k++)
                tmp <<" "<< knotSeq2(k);
        }
        internal::gsXmlNode* kn = internal::makeNode("KnotVector", tmp.str(), data);
        kn->append_attribute( internal::makeAttribute("degree", deg, data ) ) ;
        b->append_node(kn);
    }


    std::ostringstream tmp;
    tmp << std::setprecision(16);
    const int pole_countU = S->NbUPoles();
    const int pole_countV = S->NbVPoles();
    if ( pole_countU*pole_countV == 0 )
    {
        gsWarn<<"gsReadBrep: NULL Poles (control points) array\n";
    }
    else // read control points (reversed lex order wrt Gismo!)
    {
        tmp<<"\n";
        for(int i=1; i<=pole_countV; i++)
        {
            for(int j=1; j<=pole_countU; j++)
            {
                gp_Pnt p = S->Pole(j,i);
                tmp <<p.X()<<" "<<p.Y()<<" "<<p.Z()<<"\n";
            }
        }
    }

    internal::gsXmlNode* c = internal::makeNode("coefs", tmp.str(), data);
    c->append_attribute( internal::makeAttribute("geoDim", 3, data ) );
    g->append_node(c);

    if ( S->IsURational() || S->IsVRational() )
      {
        std::ostringstream weight;
        weight << std::setprecision(16);
        weight <<"\n";
        for(int i=1; i<=pole_countV; i++)
        {
            for(int j=1; j<=pole_countU; j++)
            {
                double w = S->Weight(j,i);
                weight <<w<<"\n";
            }
        }
        c = internal::makeNode("weights", weight.str(), data);
        basis_rat->append_node(c);
      }

    gsInfo <<"degree1="<< S->UDegree() <<", degree2="<< S->VDegree() <<"\n";
    gsInfo <<"periodic1="<< S->IsUPeriodic() <<", periodic2="<< S->IsVPeriodic() <<"\n";
    
    const TColStd_Array1OfReal & k1 = S->UKnots();
    const TColStd_Array1OfReal & k2 = S->VKnots();
    const TColStd_Array1OfInteger & mult1 = S->UMultiplicities();
    const TColStd_Array1OfInteger & mult2 = S->VMultiplicities();
    gsInfo <<"knots1="<< k1.Size() <<" knots2="<< k2.Size() <<"\n";
    std::copy(k1.begin(), k1.end(), std::ostream_iterator<double>(gsInfo, " ") );gsInfo<<gsEndl;
    std::copy(mult1.begin(), mult1.end(), std::ostream_iterator<double>(gsInfo, " ") );gsInfo<<gsEndl;
    std::copy(k2.begin(), k2.end(), std::ostream_iterator<double>(gsInfo, " ") );gsInfo<<gsEndl;
    std::copy(mult2.begin(), mult2.end(), std::ostream_iterator<double>(gsInfo, " ") );gsInfo<<gsEndl;


    if(S->IsURational() || S->IsVRational() )
    {
        const TColStd_Array2OfReal & w = *S->Weights();
    }

    const TColgp_Array2OfPnt   & cf = S->Poles(); // control points
    gsInfo <<"coefs1="<< S->NbUPoles() <<" coefs2="<< S->NbVPoles()<<" (i.e. "<<cf.Size() <<" ceofs.)\n";
    
    // loop over coefficients

       int aa = 0;
       for (int r = cf.LowerRow(); r <= cf.UpperRow(); r++)
       for(int c = cf.LowerCol(); c <= cf.UpperCol(); c++)
       {
       ++aa;
       const gp_Pnt & q = cf(r,c);
       gsInfo<<" ["<<q.X() <<", "<< q.Y()<< ", "<< q.Z() << "]";
       }

    gsInfo<<gsEndl;
    return true;
}

bool readGeom2d_Curve( const opencascade::handle<Geom2d_Curve> & C, internal::gsXmlTree & data  )
{
    Handle(Geom2d_BSplineCurve) ds = Handle(Geom2d_BSplineCurve)::DownCast(C);
    if ( !ds.IsNull() ) return readGeom2d_BSplineCurve(ds, data);
    return false;
}

bool readGeom2d_BSplineCurve( const opencascade::handle<Geom2d_BSplineCurve> & bsp2d, internal::gsXmlTree & data  )
{
/*

    internal::gsXmlNode* parent = data.first_node("xml") ;
    internal::gsXmlNode* g = internal::makeNode("Geometry", data);
    parent->append_node(g);
    internal::gsXmlNode* basis = internal::makeNode("Basis", data);
    internal::gsXmlNode* basis_rat = 0;

    if ( bsp2d->IsRational() )
      {
        g->append_attribute( internal::makeAttribute("type", "Nurbs", data) );
        basis_rat = internal::makeNode("Basis", data);
        basis_rat->append_attribute( internal::makeAttribute("type", "NurbsBasis", data) );
        basis = internal::makeNode("Basis", data);
        basis->append_attribute( internal::makeAttribute("type", "BSplineBasis", data) );
        basis_rat->append_node(basis);
        g->append_node(basis_rat);
      }
    else
      {
        g->append_attribute( internal::makeAttribute("type", "BSpline", data) );
        basis->append_attribute( internal::makeAttribute("type", "BSplineBasis", data) );
        g->append_node(basis);
      }

    std::ostringstream tmp;
    tmp << std::setprecision(16);
    int knot_count = bsp2d->KnotSequence().Length();
    const TColStd_Array1OfReal & knots = bsp2d->KnotSequence();
    for(int k=1; k<=knot_count; k++)
        tmp <<" "<< knots(k);

    internal::gsXmlNode* kn = internal::makeNode("KnotVector", tmp.str(), data);
    kn->append_attribute( internal::makeAttribute("degree", bsp2d->Degree(), data ) ) ;
    basis->append_node(kn);


    tmp.clear();
    tmp.str("");
    int pole_count = bsp2d->NbPoles();
    if ( pole_count == 0 )
    {
        gsWarn<<"gsReadBrep: NULL Poles (control points) array\n";
    }
    else
    {
        const TColgp_Array1OfPnt2d P = bsp2d->Poles();
        tmp <<"\n";
        for(int k=1; k<=pole_count; k++)
            tmp << P(k).X()<<" "<< P(k).Y()<<" \n";
    }

    internal::gsXmlNode* c = internal::makeNode("coefs", tmp.str(), data);
    c->append_attribute( internal::makeAttribute("geoDim", 2, data ) );
    g->append_node(c);


    if ( bsp2d->IsRational() )
      {
        std::ostringstream weight;
        weight << std::setprecision(16);
        int pole_count = bsp2d->NbPoles();
        weight <<"\n";
        for(int k=1; k<=pole_count; k++)
            weight << bsp2d->Weight(k)<<" ";
        c = internal::makeNode("weights", weight.str(), data);
        basis_rat->append_node(c);
      }

*/

//    gsInfo <<"degree="<< bsp2d->Degree() <<", periodic="<<bsp2d->IsPeriodic()<<", rational="<<bsp2d->IsRational()<<", numCoefs="<<bsp2d->NbPoles()<<", knots:\n";

//    const TColStd_Array1OfReal & rKnots = bsp2d->KnotSequence();
//    //const TColStd_Array1OfReal & uKnots = bsp2d->Knots();
//    //const TColStd_Array1OfInteger & Kmult = bsp2d->Multiplicities();
//    std::copy(rKnots.begin(), rKnots.end(), std::ostream_iterator<double>(gsInfo, " ") );gsInfo<<gsEndl;

//    const TColgp_Array1OfPnt2d & ccfs = bsp2d->Poles();

//    if (bsp2d->IsRational())
//    {
//        const TColStd_Array1OfReal & cwgts = *bsp2d->Weights();
//    }

//    for(int c = ccfs.Lower(); c <= ccfs.Upper(); c++)
//    {
//        const gp_Pnt2d & q = ccfs(c);
//        gsInfo<<" ["<<q.X() <<", "<< q.Y()<< "]";
//    }
//    gsInfo <<gsEndl;

    return true;
}

}// namespace extensions

}// namespace gismo

