/** @file gsWriteOpenNurbs.cpp

    @brief Implementation of function for data output to the Rhinoceros 3DM file format.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#include <gsOpennurbs/gsWriteOpenNurbs.h>

#include <gsCore/gsTemplateTools.h>

#include <gsCore/gsMultiPatch.h>
#include <gsModeling/gsPlanarDomain.h>

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsKnotVector.h>

#include <onurbs/opennurbs.h>

#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>


namespace gismo {

namespace extensions {


/// Writes a Curve to OpenNurbs file
template<class T>
bool writeON_NurbsCurve( const gsCurve<T> & curve, ONX_Model & model, const std::string & name)
{
      // write a wiggly cubic curve on the "green NURBS wiggle" layer
      ON_NurbsCurve* wiggle = new ON_NurbsCurve(
        3, // dimension
        false, // true if rational
        curve.degree()+1,     // order = degree+1
        curve.coefsSize()  // number of control vertices
        );

      for (int k = 0; k < wiggle->CVCount(); k++ ) 
      {
          ON_3dPoint pt( cast<T,double>(curve.coef(k,0)), cast<T,double>(curve.coef(k,1)), 0.0  ); // pt = some 3d point
          wiggle->SetCV( k, pt );
      }

      const gsKnotVector<T> & kv = 
          dynamic_cast<const gsBSplineBasis<T>&>( curve.basis() ).knots();

      // ON_NurbsCurve's have order+cv_count-2 knots.
      for (size_t k = 1; k < kv.size()-1; k++ ) 
      {
          wiggle->SetKnot(k-1, cast<T,double>(kv[k]) );
      }
      
      ON_TextLog log;
      if ( wiggle->IsValid(&log) ) 
      {
        ONX_Model_Object& mo = model.m_object_table.AppendNew();
        mo.m_object = wiggle;
        mo.m_bDeleteObject = true;
        mo.m_attributes.m_layer_index = 0;
        mo.m_attributes.m_name = name.c_str();
        //mo.m_attributes.m_uuid = ON_UUID();
      }
      else
          delete wiggle;

      return true;
}

/// Writes a Surface to OpenNurbs file
template<class T>
bool writeON_NurbsSurface( const gsSurface<T> & surface, 
                           ONX_Model & model, const std::string & name)
{
    ON_NurbsSurface* onsurf = new ON_NurbsSurface(
        3, // dimension
        false, // true if rational
        surface.basis().degree(0)+1,     // order u
        surface.basis().degree(1)+1,     // order v
        surface.basis().component(0).size(),  // number of control vertices in u
        surface.basis().component(1).size()  // number of control vertices in v
        );
    
    int c = 0;

    for ( int j = 0; j < onsurf->CVCount(1); j++ )
        for ( int i = 0; i < onsurf->CVCount(0); i++ )
        {
            ON_3dPoint pt(cast<T,double>(surface.coef(c,0)), cast<T,double>(surface.coef(c,1)), cast<T,double>(surface.coef(c,2)) );
            //ON_3dPoint pt( surface.coef(c,0), surface.coef(c,1), 0 );
            onsurf->SetCV( i, j, pt );//Note: j runs faster than i for CP(i,j)
            c++;
        }
    
      const gsKnotVector<T> & kv1 = 
          dynamic_cast<const gsBSplineBasis<T>&>( surface.basis().component(0) ).knots();
      const gsKnotVector<T> & kv2 = 
          dynamic_cast<const gsBSplineBasis<T>&>( surface.basis().component(1) ).knots();
      //Note: ON_NurbsSurface's have order+cv_count-2 knots per direction.
      for (size_t k = 1; k < kv1.size()-1; k++ ) 
          onsurf->SetKnot(0, k-1, cast<T,double>(kv1[k]) );

      for (size_t k = 1; k < kv2.size()-1; k++ ) 
          onsurf->SetKnot(1, k-1, cast<T,double>(kv2[k]) );

      ON_TextLog log;
      if ( onsurf->IsValid(&log) ) 
      {
        ONX_Model_Object& mo = model.m_object_table.AppendNew();
        mo.m_object = onsurf;
        mo.m_bDeleteObject = true;
        mo.m_attributes.m_layer_index = 0;
        mo.m_attributes.m_name = name.c_str();
        //mo.m_attributes.m_uuid = ON_UUID();
      }
      else
          delete onsurf;

      return true;
}

/// Writes a MultiPatch to OpenNurbs file
template<class T>
bool writeON_MultiPatch( const gsMultiPatch<T> & patches)
{
  //bool rc;
  const char* filename;

  ON::Begin();
  // If you want to learn to write b-rep models, first work through
  // this example paying close attention to write_trimmed_surface_example(),
  // then examime example_brep.cpp.

  // The OpenNURBS toolkit will write version 2 and 3 and read
  // version 1, 2 and 3 of the 3DM file format.
  //
  // version 1 is the legacy Rhino I/O tookit format and was used by Rhino 1.x.
  // version 2 is the OpenNURBS format (released 1 July 2000) and is used by Rhino 2.x
  // version 3 is the OpenNURBS format (released 1 November 2002) and is used by Rhino 3.x
  // version 4 is the OpenNURBS format (released September 2006) and is used by Rhino 4.x
  // version 5 is the OpenNURBS format (released September 2009) and is used by Rhino 5.x

  // version to write
  int version = 0; // version will be ON_BinaryArchive::CurrentArchiveVersion()

  // errors printed to stdout
  ON_TextLog error_log;

  // messages printed to stdout
  ON_TextLog message_log;

  // errors logged in text file
  //FILE* error_log_fp = ON::OpenFile("error_log.txt","w");
  //ON_TextLog error_log(error_log_fp);

  filename = "mp.3dm";
  FILE* fp = ON::OpenFile( filename, "wb" );

  // example demonstrates how to write a NURBS curve, line, and circle
  ONX_Model model;

  // set revision history information
  model.m_properties.m_RevisionHistory.NewRevision();

  // set application information
  model.m_properties.m_Application.m_application_name = "OpenNURBS write_curves_example() function";
  model.m_properties.m_Application.m_application_URL = "http://www.opennurbs.org";
  model.m_properties.m_Application.m_application_details = "Example program in OpenNURBS toolkit.";

  // some notes
  model.m_properties.m_Notes.m_notes = "This file was made with the OpenNURBS write_curves_example() function.";
  model.m_properties.m_Notes.m_bVisible = true;


  // file settings (units, tolerances, views, ...)
  model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::inches;
  model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.001;
  model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
  model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%


/*  ON_Layer layer;
  layer.SetLayerName("Default");
  layer.SetVisible(true);
  layer.SetLocked(false);
  layer.SetColor( ON_Color(0,0,0) );
  model.m_layer_table.Append(layer);
*/

  for(size_t i = 0; i < patches.nPatches(); ++i)
  {          
      //gsInfo<< "Write patch "<< i << "\n";
      std::stringstream  nm("patch");
      nm << i ;
      
      if ( const gsCurve<T> * c = dynamic_cast<const gsCurve<T>*>(&patches.patch(i) ) )
      {
          writeON_NurbsCurve(*c, model, nm.str() );
      }

      if ( const gsSurface<T> * c = dynamic_cast<const gsSurface<T>*>(&patches.patch(i) ) )
      {      
          writeON_NurbsSurface(*c, model, nm.str() );
      }
  }

  ON_BinaryFile archive( ON::write3dm, fp ); // fp = pointer from fopoen(...,"wb")
  // start section comment
  const char* sStartSectionComment = __FILE__ "write_points_example()" __DATE__;
  // Set uuid's, indices, etc.
  model.Polish();
  // writes model to archive
  bool ok = model.Write(archive, version, sStartSectionComment, &error_log );

  ON::CloseFile( fp );
  if (ok)
    message_log.Print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  ON::End();

    return true;
}

/// Writes a planar domain to OpenNurbs file
template<class T>
bool writeON_PlanarDomain( const gsPlanarDomain<T> & pd)
{
  //bool rc;
  const char* filename;

  ON::Begin();
  // If you want to learn to write b-rep models, first work through
  // this example paying close attention to write_trimmed_surface_example(),
  // then examime example_brep.cpp.

  // The OpenNURBS toolkit will write version 2 and 3 and read
  // version 1, 2 and 3 of the 3DM file format.
  //
  // version 1 is the legacy Rhino I/O tookit format and was used by Rhino 1.x.
  // version 2 is the OpenNURBS format (released 1 July 2000) and is used by Rhino 2.x
  // version 3 is the OpenNURBS format (released 1 November 2002) and is used by Rhino 3.x
  // version 4 is the OpenNURBS format (released September 2006) and is used by Rhino 4.x
  // version 5 is the OpenNURBS format (released September 2009) and is used by Rhino 5.x

  // version to write
  int version = 0; // version will be ON_BinaryArchive::CurrentArchiveVersion()

  // errors printed to stdout
  ON_TextLog error_log;

  // messages printed to stdout
  ON_TextLog message_log;

  // errors logged in text file
  //FILE* error_log_fp = ON::OpenFile("error_log.txt","w");
  //ON_TextLog error_log(error_log_fp);

  filename = "pd.3dm";
  FILE* fp = ON::OpenFile( filename, "wb" );

  // example demonstrates how to write a NURBS curve, line, and circle
  ONX_Model model;

  // set revision history information
  model.m_properties.m_RevisionHistory.NewRevision();

  // set application information
  model.m_properties.m_Application.m_application_name = "OpenNURBS write_curves_example() function";
  model.m_properties.m_Application.m_application_URL = "http://www.opennurbs.org";
  model.m_properties.m_Application.m_application_details = "Example program in OpenNURBS toolkit.";

  // some notes
  model.m_properties.m_Notes.m_notes = "This file was made with the OpenNURBS write_curves_example() function.";
  model.m_properties.m_Notes.m_bVisible = true;


  // file settings (units, tolerances, views, ...)
  model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::inches;
  model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.001;
  model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
  model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%


/*  ON_Layer layer;
  layer.SetLayerName("Default");
  layer.SetVisible(true);
  layer.SetLocked(false);
  layer.SetColor( ON_Color(0,0,0) );
  model.m_layer_table.Append(layer);
*/

  for(index_t i =0; i<pd.numLoops();i++)
      for(index_t j =0; j< pd.loop(i).numCurves() ; j++)
  {
      const gsCurve<T> & c = pd.loop(i).curve(j);
      gsInfo<< "Write loop "<< i <<", curve "<<j<<"\n";

      std::stringstream  nm("curve");
      nm << i <<"_"<<j;
      
      writeON_NurbsCurve(c, model, nm.str() );

    }

  ON_BinaryFile archive( ON::write3dm, fp ); // fp = pointer from fopoen(...,"wb")
  // start section comment
  const char* sStartSectionComment = __FILE__ "write_points_example()" __DATE__;
  // Set uuid's, indices, etc.
  model.Polish();
  // writes model to archive
  bool ok = model.Write(archive, version, sStartSectionComment, &error_log );

  ON::CloseFile( fp );
  if (ok)
    message_log.Print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  ON::End();

    return true;
}


}// namespace extensions

}// namespace gismo
