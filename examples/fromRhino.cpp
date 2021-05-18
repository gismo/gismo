/** @file fromRhino.cpp

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

  gsCmdLine cmd("Hi, give me a file (.3dm) and I will write it to G+Smo!!");
    cmd.addPlainString("filename", "Filename containing readable data (.xml or third party)", fn);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  if (fn.empty() )
  {
      gsInfo<< cmd.getMessage();
      gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
      return 0;
  }

  gsFileData<>  data;
  data.read(fn);
  data.dump();
  gsInfo<< "File contains "<< data.numTags() <<" objects, total size is "<< data.bufferSize() <<"\n";

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
