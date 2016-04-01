/** @file gsReadOpenNurbs.cpp

    @brief Implementation of function for data input from the Rhinoceros 3DM file format.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#include <gsOpennurbs/gsReadOpenNurbs.h>

// Note: for building a dll use:
// #define ON_DLL_EXPORTS
// #define ON_DLL_IMPORTS
#include <onurbs/opennurbs.h>
//#include "gsOpennurbs/examples_linking_pragmas.h"

#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>


namespace gismo {

namespace extensions {

/*
Returns:
  True if .3dm file was successfully read into an ONX_Model class.
*/
  bool gsReadOpenNurbs( const char * arg, internal::gsXmlTree & data )
{
  ON_wString ws_arg = arg;
  const wchar_t* sFileName = ws_arg;


  // Opennurbs Model object
  ONX_Model model;
  //ON_TextLog dump;

  // open file containing opennurbs archive
  FILE* archive_fp = ON::OpenFile( sFileName, L"rb");
  if ( !archive_fp ) 
  {
    gsWarn<< "  Unable to open 3dm file: file not valid.\n";
    return false;
  }

  // create achive object from file pointer
  ON_BinaryFile archive( ON::read3dm, archive_fp );

  // read the contents of the file into "model"
  bool rc = model.Read( archive );// optionally ON_TextLog

  // close the file
  ON::CloseFile( archive_fp );
  
  if ( ! rc )
  {
      gsWarn<< "Error while reading 3dm file.\n";
      return false;
  } 

  ON_TextLog log;
  if ( ! model.IsValid(&log) )
  {
    gsWarn<< "OpenNurbs reported an invalid 3dm file.\n";
    return false;
  }

  for (int i = 0; i < model.m_object_table.Count(); i++ )
  {
    const ON_Geometry* pGeometry = ON_Geometry::Cast( model.m_object_table[i].m_object );
  if ( pGeometry ) 
  {
    if ( ON_Curve::Cast(pGeometry) ) 
      {
	const ON_Curve * pcurve = ON_Curve::Cast(pGeometry);
	rc = rc & readON_Curve( pcurve, data );
      }
    else if ( ON_Surface::Cast(pGeometry) ) 
      {
	const ON_Surface * psurface = ON_Surface::Cast(pGeometry);
	rc = rc & readON_Surface( psurface, data );
      }    
    else if ( ON_Brep::Cast(pGeometry) ) 
    {
      // pGeometry is derived from ON_Brep
      const ON_Brep* brep = ON_Brep::Cast(pGeometry);
      rc = rc & readON_Brep( brep, data );
    }
    else if ( ON_NurbsCage::Cast(pGeometry) ) 
    {
      const ON_NurbsCage* cage = ON_NurbsCage::Cast(pGeometry);
      rc = rc & readON_NurbsCage( cage, data );
    }
    else if ( ON_MorphControl::Cast(pGeometry) ) 
    {
      const ON_MorphControl* mc = ON_MorphControl::Cast(pGeometry);
      rc = rc & readON_MorphControl( mc, data );
    }
    else if ( ON_Mesh::Cast(pGeometry) ) 
    {
      gsInfo<< "Found an ON_Mesh Object.\n";
      //const ON_Mesh* pMesh = ON_Mesh::Cast(pGeometry);
    }
    else if ( ON_Extrusion::Cast(pGeometry) ) 
      {
	// pGeometry is derived from ON_Extrusion
	//  Note: 
	//   ON_Extrusion::BrepForm() will return a brep form
	//   if you don't want to explicitly handle extrusions. 
	gsInfo<< "In fact it is a ON_Extrusion.\n";
	//const ON_Extrusion* extrusion = ON_Extrusion::Cast(pGeometry);
      }
    else if ( pGeometry->HasBrepForm() ) 
    {
      gsInfo<< "Found an object which HasBrepForm (ON_Brep).\n";
      // pGeometry is not derived from ON_Brep but its geometry can
      // be represented by an ON_Brep.
      ON_Brep* brep = pGeometry->BrepForm();
      rc = rc & readON_Brep( brep, data );
      // you manage the ON_Brep returned by pGeometry->BrepForm();
      delete brep;
    }
    else{
      gsInfo<< "Found a geometry object (discarded).\n";
      ON_TextLog dump;
      pGeometry->Dump(dump);
    }
  }
  else{
    gsWarn<< "Discarding non-geometry object ("<<i<<") from 3dm file.\n";
  }

  }

  // OPTIONAL: Call just before your application exits to clean
  //           up opennurbs class definition information.
  //           Opennurbs will not work correctly after ON::End()
  //           is called.
  ON::End();
  
  return rc;
}


bool readON_NurbsSurface( const ON_NurbsSurface * psurface, internal::gsXmlTree & data  )
{
  //ON_TextLog dump; psurface->Dump(dump);
  internal::gsXmlNode* parent = data.first_node("xml") ;
  
  internal::gsXmlNode* g = internal::makeNode("Geometry", data);    
  parent->append_node(g);
  internal::gsXmlNode* basis = internal::makeNode("Basis", data);
  internal::gsXmlNode* basis_rat = 0;
  
  if ( psurface->m_is_rat ) 
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
  
  int dir;
  for ( dir = 0; dir < 2; dir++ ) 
    {
      std::ostringstream tmp;
      tmp << std::setprecision(16);
      internal::gsXmlNode* b = internal::makeNode("Basis", data);
      b->append_attribute( internal::makeAttribute("type", "BSplineBasis", data) );
      b->append_attribute( internal::makeAttribute("index", dir, data) );
      basis->append_node(b);
      
      int knot_count = ON_KnotCount( psurface->m_order[dir], psurface->m_cv_count[dir] );
      int k=-1;
      const double* knot = psurface->m_knot[dir];
      tmp << psurface->SuperfluousKnot(dir,0);
      while ( ++k < knot_count ) 
	{
	  tmp <<" "<< knot[k];	  
	}
      tmp <<" " << psurface->SuperfluousKnot(dir,1);
      internal::gsXmlNode* kn = internal::makeNode("KnotVector", tmp.str(), data);
      kn->append_attribute( internal::makeAttribute("degree", psurface->m_order[dir]-1, data ) ) ;
      b->append_node(kn);
    }
  
  std::ostringstream tmp;
  tmp << std::setprecision(16);
  std::ostringstream weight;
  weight << std::setprecision(16);
  if ( !psurface->m_cv ) 
    {
      gsWarn<<"gsReadOpenNurbs: NULL cv array\n";
    }
  else 
    {
      const int cvdim =  psurface->m_dim;
      // read control points (reversed lex order wrt Gismo!)
      if ( psurface->m_is_rat ) 
	{	  
	  double w;
	  for (int i = 0; i < psurface->m_cv_count[1]; i++ ) 
	    {
	      const double * P = psurface->CV(0,i);
	      for (int k = 0; k < psurface->m_cv_count[0]; k++ ) 
		{
		  w = P[cvdim];
		  weight << " "<< w ;
		  if ( w != 0.0 ) // account for points at infinity
		    w = 1.0/w;
		  for (int j = 0; j < cvdim; j++ ) 
		    tmp <<" "<< w * P[j]  ;
		  P += psurface->m_cv_stride[0];
		}
	    }
	}
      else
	{
	  for (int i = 0; i < psurface->m_cv_count[1]; i++ ) 
	    {
	      const double * P = psurface->CV(0,i);	      
	      for (int k = 0; k < psurface->m_cv_count[0]; k++ ) 
		{
		  tmp << " "<< P[0];
		  for (int j = 1; j < cvdim; j++ ) 
		    tmp <<" "<< P[j];
		  P += psurface->m_cv_stride[0];
		}
	    }
	}
    }
  internal::gsXmlNode* c = internal::makeNode("coefs", tmp.str(), data);
  c->append_attribute( internal::makeAttribute("geoDim", psurface->m_dim, data ) );
  g->append_node(c);
  if ( psurface->m_is_rat ) 
    {
      c = internal::makeNode("weights", weight.str(), data);
      basis_rat->append_node(c);
    }
  
  return true;
}



bool readON_NurbsCurve( const ON_NurbsCurve * pcurve, internal::gsXmlTree & data  )
{
  //ON_TextLog dump; pcurve->Dump(dump);

  internal::gsXmlNode* parent = data.first_node("xml") ;
  
  internal::gsXmlNode* g = internal::makeNode("Geometry", data);    
  parent->append_node(g);
  internal::gsXmlNode* basis = internal::makeNode("Basis", data);
  internal::gsXmlNode* basis_rat = 0;
  
  if ( pcurve->m_is_rat ) 
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
  int knot_count = ON_KnotCount( pcurve->m_order, pcurve->m_cv_count );
  int k=-1;
  const double* knot = pcurve->m_knot;
  tmp << pcurve->SuperfluousKnot(0);
  while ( ++k < knot_count ) 
    {
      tmp <<" "<< knot[k];	  
    }
  tmp <<" " << pcurve->SuperfluousKnot(1);
  internal::gsXmlNode* kn = internal::makeNode("KnotVector", tmp.str(), data);
  kn->append_attribute( internal::makeAttribute("degree", pcurve->m_order-1, data ) ) ;
  basis->append_node(kn);
  
  tmp.clear();
  tmp.str("");
  std::ostringstream weight;
  weight << std::setprecision(16);
  if ( !pcurve->m_cv ) 
    {
      gsWarn<<"gsReadOpenNurbs: NULL cv array\n";
    }
  else 
    {
      const int cvdim =  pcurve->m_dim;
      if ( pcurve->m_is_rat ) 
	{
	  double w;
	  const double * P = pcurve->m_cv;
	  for (int k = 0; k != pcurve->m_cv_count; k++ ) 
	    {
	      w = P[cvdim];
	      weight << " "<< w ;
	      if ( w != 0.0 ) // account for points at infinity
		w = 1.0/w;
	      for (int j = 0; j < cvdim; j++ ) 
		tmp <<" "<< w * P[j]  ;
	      P += cvdim+1;
	    }
	}
      else
	{
	  const double * P = pcurve->m_cv;
	  for (int k = 0; k < pcurve->m_cv_count; k++ ) 
	    {
	      for (int j = 0; j < cvdim; j++ ) 
		tmp <<" "<< P[k];
	      P += cvdim;
	    }
	}
    }

  internal::gsXmlNode* c = internal::makeNode("coefs", tmp.str(), data);
  c->append_attribute( internal::makeAttribute("geoDim", pcurve->m_dim, data ) );
  g->append_node(c);
  if ( pcurve->m_is_rat ) 
    {
      c = internal::makeNode("weights", weight.str(), data);
      basis_rat->append_node(c);
    }
  
  return true;
}

// NOTE;
// Surfaces in Rhino are read as ON_Breps
// Internally, Rhino stores all surfaces as some type of b-rep and the
// openNURBS toolkit reads these objects as b-reps.  To see if an entire
// ON_Brep is really just a surface, use
//    BOOL ON_Brep::IsSurface().
// If ON_Brep::IsSurface() returns TRUE, then the b-rep geometry
// is the same as the surface ON_Brep::m_S[0].
// To see of a particular face in a b-rep is really just a surface,
// use 
//    BOOL ON_Brep::FaceIsSurface( face_index ).
// If ON_Brep::FaceIsSurface( face_index ) returns TRUE, then the
// face's geometry is the same as the surface ON_Brep::m_S[face.m_si].
bool readON_Brep( const ON_Brep * pbrep, internal::gsXmlTree & data  )
{
  //ON_TextLog dump; pbrep->Dump(dump);

  int n = pbrep->m_F.Count();
  for (int  i = 0; i != n ; i++ )
    {
      if ( pbrep->FaceIsSurface(i) )
      {
	const ON_Surface* g = pbrep->m_F[i].SurfaceOf();
	if ( g )
	  readON_Surface(g, data);
      else
	  gsWarn<<"Null surface in Brep\n";
      }
      else// trimmed surf.
	{
	  gsWarn<<"Reading trimmed surf.\n";
	  // pbrep->m_F[i].LoopCount();
	  // ON_BrepLoop* lp = pbrep->m_F[i].OuterLoop()
	  // ON_BrepLoop* lp = pbrep->m_F[i].Loop(i)
	  // 1256
	  // lp->TrimCount()
	  // ON_BrepTrim* tm = lp->Trim(i) 
	  //const ON_Curve* tm->TrimCurveOf();

	  const ON_Surface* g = pbrep->m_F[i].SurfaceOf();
	  if ( g )
	    readON_Surface(g, data);
	  else
	    gsWarn<<"Null surface in Brep\n";
	}
    }


#if 0
  //Look for surfaces (geometries on faces)
  int n = pbrep->m_S.Count();
  for (int  i = 0; i != n ; i++ )
    {
      const ON_Surface* g = pbrep->m_S[i];
      if ( g )
	  readON_Surface(g, data);
      else
	  gsWarn<<"Null surface in Brep\n";
    }

  //Look for 3D curves (geometries on egdes)
  n = pbrep->m_C3.Count();
  for (int  i = 0; i != n ; i++ )
    {
      const ON_Curve* g = pbrep->m_C3[i];
      if ( g )
	readON_Curve(g, data);
      else
	gsWarn<<"Null 3D curve in Brep\n";
    }

  //Look for 2D curves (trimming curves)
  n = pbrep->m_C2.Count();
  for (int  i = 0; i != n ; i++ )
    {
      const ON_Curve* g = pbrep->m_C2[i];
      if ( g )
	readON_Curve(g, data);
      else
	gsWarn<<"Null 2D curve in Brep\n";
    }
#endif


  return true;
}

bool readON_Surface ( const ON_Surface * psurface, internal::gsXmlTree & data  )
{
  if ( ON_NurbsSurface::Cast(psurface) ) {
    const ON_NurbsSurface* p = ON_NurbsSurface::Cast(psurface);
    return readON_NurbsSurface( p, data);
  }
  else if ( psurface->HasNurbForm() ) {
    ON_NurbsSurface nf;
    psurface->GetNurbForm(nf);
    return readON_NurbsSurface( &nf, data );
  }
  else if ( ON_PlaneSurface::Cast(psurface) ) {
    //gsInfo<< "In fact it is a ON_PlaneSurface.\n";
    //const ON_PlaneSurface* p = ON_PlaneSurface::Cast(psurface);
  }
  else if ( ON_RevSurface::Cast(psurface) ) {
    //gsInfo<< "In fact it is a ON_RevSurface.\n";
    //const ON_RevSurface* p = ON_RevSurface::Cast(psurface);
  }
  else if ( ON_BrepFace::Cast(psurface) ) {
    //gsInfo<< "In fact it is a ON_BrepFace.\n";
    //const ON_BrepFace* p = ON_BrepFace::Cast(psurface);
  }
  else if ( ON_SurfaceProxy::Cast(psurface) ) {
    //gsInfo<< "In fact it is a ON_SurfaceProxy.\n";
    //const ON_SurfaceProxy* p = ON_SurfaceProxy::Cast(psurface);
  }
  else {
    ON_TextLog dump; psurface->Dump(dump);
    gsInfo<< "Unidentified ON_Surface.\n";
  }

  return false;    
}

bool readON_Curve( const ON_Curve   * pcurve  , internal::gsXmlTree & data  )
{
  if ( ON_NurbsCurve::Cast(pcurve) ) {
    //gsInfo<< "In fact it is a ON_NurbsCurve.\n";
    const ON_NurbsCurve* p = ON_NurbsCurve::Cast(pcurve);
    return readON_NurbsCurve( p, data);
  }
  else if ( pcurve->HasNurbForm() ) {
    ON_NurbsCurve nf;
    pcurve->GetNurbForm(nf);
    return readON_NurbsCurve( &nf, data );
  }
  else if ( ON_ArcCurve::Cast(pcurve) ) {
    //gsInfo<< "In fact it is a ON_ArcCurve.\n";
    //const ON_ArcCurve* p = ON_ArcCurve::Cast(pcurve);
  }
  else if ( ON_CurveOnSurface::Cast(pcurve) ) {
    //gsInfo<< "In fact it is a ON_CurveOnSurface.\n";
    //const ON_CurveOnSurface* p = ON_CurveOnSurface::Cast(pcurve);
  }
  else if ( ON_BrepEdge::Cast(pcurve) ) {
    //gsInfo<< "In fact it is a ON_BrepEdge.\n";
    //const ON_BrepEdge* p = ON_BrepEdge::Cast(pcurve);
  }
  else if ( ON_LineCurve::Cast(pcurve) ) {
    //gsInfo<< "In fact it is a ON_LineCurve.\n";
    //const ON_LineCurve* p = ON_LineCurve::Cast(pcurve);
  }
  else if ( ON_PolyCurve::Cast(pcurve) ) {
    //gsInfo<< "In fact it is a ON_PolyCurve.\n";
    //const ON_PolyCurve* p = ON_PolyCurve::Cast(pcurve);
  }
  else if ( ON_PolylineCurve::Cast(pcurve) ) {
    //gsInfo<< "In fact it is a ON_PolylineCurve.\n";
    //const ON_PolylineCurve* p = ON_PolylineCurve::Cast(pcurve);
  }
  else if ( ON_CurveProxy::Cast(pcurve) ) {
    //gsInfo<< "In fact it is a ON_CurveProxy.\n";
    //const ON_CurveProxy* p = ON_CurveProxy::Cast(pcurve);
  }
  else {
    ON_TextLog dump; pcurve->Dump(dump);
    gsInfo<< "Unidentified ON_Curve.\n";
  }

  return false;    
}


bool readON_MorphControl( const ON_MorphControl * mc  , internal::gsXmlTree & data  )
{
  switch(mc->m_varient)
  {
  case 1:
    //mc->m_nurbs_curve0;// ref. object
    readON_NurbsCurve( &mc->m_nurbs_curve, data);
    return true;
    break;
  case 2:
    //mc->m_nurbs_surface0;// ref. object
    readON_Surface( &mc->m_nurbs_surface, data);
    return true;
    break;
  case 3:
    //mc->m_nurbs_cage0;// ref object
    readON_NurbsCage( & mc->m_nurbs_cage, data);
    return true;
    break;
  default:
    return false;
    break;
  }
}


bool readON_NurbsCage( const ON_NurbsCage * pcage, internal::gsXmlTree & data  )
{
  //ON_TextLog dump; pcage->Dump(dump);
  internal::gsXmlNode* parent = data.first_node("xml") ;
  
  internal::gsXmlNode* g = internal::makeNode("Geometry", data);    
  parent->append_node(g);
  internal::gsXmlNode* basis = internal::makeNode("Basis", data);
  internal::gsXmlNode* basis_rat = 0;
  
  if ( pcage->m_is_rat ) 
    {
      g->append_attribute( internal::makeAttribute("type", "TensorNurbs3", data) );
      basis_rat = internal::makeNode("Basis", data);
      basis_rat->append_attribute( internal::makeAttribute("type", "TensorNurbsBasis3", data) );
      basis = internal::makeNode("Basis", data);
      basis->append_attribute( internal::makeAttribute("type", "TensorBSplineBasis3", data) );
      basis_rat->append_node(basis);
      g->append_node(basis_rat);
    }
  else
    {
      g->append_attribute( internal::makeAttribute("type", "TensorBSpline3", data) );
      basis->append_attribute( internal::makeAttribute("type", "TensorBSplineBasis3", data) );
      g->append_node(basis);
    }
  
  int dir;
  for ( dir = 0; dir < 3; dir++ ) 
    {
      std::ostringstream tmp;
      tmp << std::setprecision(16);
      internal::gsXmlNode* b = internal::makeNode("Basis", data);
      b->append_attribute( internal::makeAttribute("type", "BSplineBasis", data) );
      b->append_attribute( internal::makeAttribute("index", dir, data) );
      basis->append_node(b);
      
      int knot_count = ON_KnotCount( pcage->m_order[dir], pcage->m_cv_count[dir] );
      int k=-1;
      const double* knot = pcage->m_knot[dir];

      tmp << ON_SuperfluousKnot(pcage->m_order[dir],pcage->m_cv_count[dir],
				pcage->m_knot[dir], 0 );
      while ( ++k < knot_count ) 
	{
	  tmp <<" "<< knot[k];	  
	}
      tmp <<" "<< ON_SuperfluousKnot(pcage->m_order[dir],pcage->m_cv_count[dir],
				pcage->m_knot[dir], 1 );
      internal::gsXmlNode* kn = internal::makeNode("KnotVector", tmp.str(), data);
      kn->append_attribute( internal::makeAttribute("degree", pcage->m_order[dir]-1, data ) ) ;
      b->append_node(kn);
    }
  
  std::ostringstream tmp;
  tmp << std::setprecision(16);
  std::ostringstream weight;
  weight << std::setprecision(16);
  if ( !pcage->m_cv ) 
    {
      gsWarn<<"gsReadOpenNurbs: NULL cv array\n";
    }
  else 
    {
      const int cvdim =  3;//pcage->m_dim;
      // read control points (reversed lex order wrt Gismo!)
      if ( pcage->m_is_rat ) 
	{	  
	  double w;
	  for (int t = 0; t < pcage->m_cv_count[2]; t++ ) 
	    for (int i = 0; i < pcage->m_cv_count[1]; i++ ) 
	      {
		const double * P = pcage->CV(0,i,t);
		for (int k = 0; k < pcage->m_cv_count[0]; k++ ) 
		  {
		    w = P[cvdim];
		    weight << " "<< w ;
		    if ( w != 0.0 ) // account for points at infinity
		      w = 1.0/w;
		    for (int j = 0; j < cvdim; j++ ) 
		      tmp <<" "<< w * P[j]  ;
		    P += pcage->m_cv_stride[0];
		  }
	      }
	}
      else
	{
	  for (int t = 0; t < pcage->m_cv_count[2]; t++ ) 
	    for (int i = 0; i < pcage->m_cv_count[1]; i++ ) 
	      {
		const double * P = pcage->CV(0,i,t);	      
		for (int k = 0; k < pcage->m_cv_count[0]; k++ ) 
		  {
		    tmp << " "<< P[0];
		    for (int j = 1; j < cvdim; j++ ) 
		      tmp <<" "<< P[j];
		    P += pcage->m_cv_stride[0];
		  }
	      }
	}
    }
  internal::gsXmlNode* c = internal::makeNode("coefs", tmp.str(), data);
  c->append_attribute( internal::makeAttribute("geoDim", pcage->m_dim, data ) );
  g->append_node(c);
  if ( pcage->m_is_rat ) 
    {
      c = internal::makeNode("weights", weight.str(), data);
      basis_rat->append_node(c);
    }
  
  return true;
}



}// namespace extensions

}// namespace gismo
