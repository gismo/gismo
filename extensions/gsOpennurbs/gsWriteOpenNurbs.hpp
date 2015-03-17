/** @file gsWriteOpenNurbs.cpp

    @brief Implementation of function for data output to the Rhinoceros 3DM file format.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


#include <gsCore/gsTemplateTools.h>

#include <gsModeling/gsPlanarDomain.h>

#include <gsNurbs/gsBSplineBasis.h>
#include <gsNurbs/gsKnotVector.h>


#include <gsOpennurbs/gsWriteOpenNurbs.h>

#include <onurbs/opennurbs.h>

#include <sstream>
#include <string>
#include <fstream>
#include <iomanip>


namespace gismo {

namespace extensions {

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

      // write a wiggly cubic curve on the "green NURBS wiggle" layer
      ON_NurbsCurve* wiggle = new ON_NurbsCurve(
        3, // dimension
        false, // true if rational
        c.degree()+1,     // order = degree+1
        c.coefsSize()  // number of control vertices
        );

      for (int k = 0; k < wiggle->CVCount(); k++ ) 
      {
          ON_3dPoint pt( c.coef(k,0), c.coef(k,1), 0.0  ); // pt = some 3d point
          wiggle->SetCV( k, pt );
      }

      const gsKnotVector<T> & kv = 
          dynamic_cast<const gsBSplineBasis<T>&>( c.basis() ).knots();
      // ON_NurbsCurve's have order+cv_count-2 knots.
      for (int k = 1; k < kv.size()-1; k++ ) 
      {
          wiggle->SetKnot(k-1, kv[k] );
      }
      
      if ( wiggle->IsValid() ) 
      {
        ONX_Model_Object& mo = model.m_object_table.AppendNew();
        mo.m_object = wiggle;
        mo.m_bDeleteObject = true;
        mo.m_attributes.m_layer_index = 0;
        std::stringstream  nm("curve");
        nm << i <<"_"<<j;
        mo.m_attributes.m_name = nm.str().c_str();
        //mo.m_attributes.m_uuid = ON_UUID();
      }
      else
      {
          gsInfo<< "Invalid !! \n";      
          delete wiggle;
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


}// namespace extensions

}// namespace gismo
