/** @file fileIo_example.cpp

    @brief Testing file reading and writing

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

#ifdef gsOpennurbs_ENABLED
#include "gsOpennurbs/gsWriteOpenNurbs.h"
#endif

#ifdef gsOpenCascade_ENABLED
#include "gsOpenCascade/gsReadOcct.h"
#include "gsOpenCascade/gsWriteOcct.h"
#endif


using namespace gismo;


template<class object>
void lookFor( const gsFileData<> & data )
{
  if ( data.has< object >() )
  {
      gsInfo<<"* There are "<< data.count< object >() <<" "
          <<data.type<object>()<<" "<< data.tag<object>() <<" in the file."<<"\n";

    memory::unique_ptr<object> o = data.getFirst< object >();
    gsInfo<< "  Read the first one in address "<< o.get() << "\n";
    gsInfo<< "  "<<*o ;


    gsInfo<< "  Write back to a new file \"dump_write.xml\"" << "\n";
    gsFileData<> newdata;
    newdata << *o ;

    // Compressed XML format xml.gz
    newdata.saveCompressed("dump_write");

    // Non-compressed XML
    newdata.save("dump_write");
  }
}


int main(int argc, char *argv[])
{
  
#ifdef gsOpennurbs_ENABLED
    gsInfo << "Using opennurbs.\n";
#endif

#ifdef gsOpenCascade_ENABLED
    gsInfo << "Using OpenCascade.\n";
#endif

  std::string fn;
  //fn = "basis_thbs_01.xml"; //default example
  
  gsCmdLine cmd("Hi, give me a file (.xml, .txt, .axl) and I will read the contents!");
    cmd.addPlainString("filename", "Filename containing readable data (.xml or third party)", fn); 
  
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  
  if (fn.empty() )
  {
      gsInfo<< cmd.getMessage();
      gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
      return 0;
  }
  
  gsFileData<>  data( fn );
  data.save();
  gsInfo<< "File contains "<< data.numTags() <<" objects, total size is "<< data.bufferSize() <<"\n";
  //gsInfo<< data.contents() <<"\n";
  gsInfo<< "Dumped all file data in dump.xml without allocating Gismo objects."<<"\n";
  
  lookFor< gsBasis<>      > (data) ;
  lookFor< gsGeometry<>   > (data) ;
  lookFor< gsBSpline<>    > (data) ;
  lookFor< gsNurbs<>      > (data) ;
  //lookFor< gsBezier<>     > (data) ;
  lookFor< gsTensorBSpline<2,real_t> > (data) ;
  lookFor< gsTensorBSpline<3,real_t> > (data) ;
  lookFor< gsTensorNurbs<  2,real_t> > (data) ;
  lookFor< gsTensorNurbs<  3,real_t> > (data) ;
  //lookFor< gsTensorBezier<2, real_t> > (data) ;
  //lookFor< gsTensorBezier<2, real_t> > (data) ;
  
  lookFor< gsMesh<> > (data) ;

  lookFor< gsSolid<> > (data) ;

  lookFor< gsKnotVector<> > (data) ;

  lookFor< gsMatrix<> > (data) ;

  lookFor< gsSparseMatrix<> > (data) ;

  lookFor< gsPoissonPde<> > (data) ;

  lookFor< gsPlanarDomain<> > (data) ;

  lookFor< gsCurveLoop<> > (data) ;

  lookFor< gsTrimSurface<> > (data) ;

  lookFor< gsMultiPatch<> > (data) ;

  lookFor< gsFunctionExpr<> > (data);

  lookFor< gsHBox<2,real_t> > (data);

  lookFor< gsHBoxContainer<2,real_t> > (data);

    
#if defined(gsOpennurbs_ENABLED) || defined(gsOpenCascade_ENABLED)
  if ( data.has< gsPlanarDomain<> >() )
  {
      gsInfo<<"* There is a "<< data.count< gsPlanarDomain<> >() <<" "
          <<data.type<gsPlanarDomain<> >()<<" "<< data.tag<gsPlanarDomain<> >() 
          <<" in the file.\n";

    gsPlanarDomain<>::uPtr o = data.getFirst< gsPlanarDomain<> >();
    gsInfo<< "  Read it ..\n";
    gsInfo<< "  "<<*o ;

#ifdef gsOpennurbs_ENABLED
    gsInfo<< "  Write back to pd.3dm\n";
    extensions::writeON_PlanarDomain(*o,"pd");
#endif
  }

  if ( data.has< gsMultiPatch<> >() )
  {
      gsInfo<<"* There is "<< data.count< gsMultiPatch<> >() <<" "
          <<data.type< gsMultiPatch<> >()<<" "<< data.tag< gsMultiPatch<> >() 
          <<" in the file.\n";

      gsMultiPatch<>::uPtr o = data.getFirst< gsMultiPatch<> >();
      gsInfo<< "  Read it ..\n";
      gsInfo<< "  "<<*o ;

#ifdef gsOpennurbs_ENABLED
      gsInfo<< "  Write back to mp.3dm\n";
      extensions::writeON_MultiPatch(*o,"mp");
#endif
#ifdef gsOpenCascade_ENABLED
      gsInfo<< "  Write back to mp.igs\n";
      extensions::writeON_MultiPatch(*o,"mp");
#endif
  }

  if ( data.has< gsGeometry<> >() )
  {
      gsInfo<<"* There is "<< data.count< gsGeometry<> >() <<" "
          <<data.type< gsGeometry<> >()<<" "<< data.tag< gsGeometry<> >() 
          <<" in the file.\n";

      gsGeometry<>::uPtr o = data.getFirst< gsGeometry<> >();
      gsInfo<< "  Read it ..\n";
      gsInfo<< "  "<<*o ;

      if ( gsSurface<> * srf = dynamic_cast<gsSurface<>*>(o.get()) )
      {      
#ifdef gsOpennurbs_ENABLED
      gsInfo<< "  Write back to geo.3dm\n";
      extensions::writeON_NurbsSurface(*srf,"geo");
#endif
#ifdef gsOpenCascade_ENABLED
      gsInfo<< "  Write back to geo.igs\n";
      extensions::writeOcctIges(*srf,"geo");
      gsInfo<< "  Write back to geo.step\n";
      extensions::writeOcctStep(*srf,"geo");
#endif
      GISMO_UNUSED(srf);
      }
  }

  if ( data.has< gsMesh<> >() )
  {
      gsInfo<<"* There is "<< data.count< gsMesh<> >() <<" "
          <<data.type< gsMesh<> >()<<" "<< data.tag< gsMesh<> >()
          <<" in the file.\n";

      gsMesh<>::uPtr o = data.getFirst< gsMesh<> >();
      gsInfo<< "  Read it ..\n";
      gsInfo<< "  "<<*o ;
#ifdef gsOpennurbs_ENABLED
      gsInfo<< "  Write back to mesh.3dm\n";
      extensions::writeON_Mesh(*o,"mesh");
#endif
  }

#endif

  return EXIT_SUCCESS;
}
