/** @file gsReadWrite_test.cpp

    @brief Testing file reading and writing

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

#ifdef GISMO_WITH_ONURBS
#include "gsOpennurbs/gsWriteOpenNurbs.h"
#endif

#ifdef GISMO_WITH_OCC
#include "gsOpenCascade/gsReadBrep.h"
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
  
#ifdef GISMO_WITH_ONURBS
    gsInfo << "Using opennurbs.\n";
#endif

#ifdef GISMO_WITH_OCC
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
    
#ifdef GISMO_WITH_ONURBS
  if ( data.has< gsPlanarDomain<> >() )
  {
      gsInfo<<"* There is a "<< data.count< gsPlanarDomain<> >() <<" "
          <<data.type<gsPlanarDomain<> >()<<" "<< data.tag<gsPlanarDomain<> >() 
          <<" in the file.\n";

    gsPlanarDomain<>::uPtr o = data.getFirst< gsPlanarDomain<> >();
    gsInfo<< "  Read it ..\n";
    gsInfo<< "  "<<*o ;

    gsInfo<< "  Write back to pd.3dm\n";
    extensions::writeON_PlanarDomain(*o);
  }

  if ( data.has< gsMultiPatch<> >() )
  {
      gsInfo<<"* There is "<< data.count< gsMultiPatch<> >() <<" "
          <<data.type< gsMultiPatch<> >()<<" "<< data.tag< gsMultiPatch<> >() 
          <<" in the file.\n";

      gsMultiPatch<>::uPtr o = data.getFirst< gsMultiPatch<> >();
      gsInfo<< "  Read it ..\n";
      gsInfo<< "  "<<*o ;
      
      gsInfo<< "  Write back to mp.3dm\n";
      extensions::writeON_MultiPatch(*o);
  }
 
#endif

#ifdef GISMO_WITH_OCC

 
#endif
  
  return 0;
}
