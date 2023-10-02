/** @file getFileinMP.cpp

    @brief Testing file reading and writing

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Imperatore
*/

#include <iostream>

#include <gismo.h>
using namespace gismo;

int main(int argc, char *argv[])
{
  std::string fn;
  fn = "";

  gsCmdLine cmd("Input: .xml file with gsGeometry, to be stored as gsMultiPatch");
  cmd.addPlainString("filename", "Filename containing readable data .xml", fn);

  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

  if (fn.empty() )
  {
      gsInfo<< cmd.getMessage();
      gsInfo<<"\nType "<< argv[0]<< " -h, to get the list of command line options.\n";
      return 0;
  }

  gsFileData<>  data( fn );
  data.save();
  gsInfo<< "File contains "<< data.numTags() <<" gsGeometry objects.\n";

  gsFileData<> fd;

  std::vector< memory::unique_ptr< gsGeometry<> > >  geoContainer = data.getAll< gsGeometry<> >();
  gsMultiPatch<> mpBSplineCurves;

  for(index_t el = 0; el < geoContainer.size(); el++)
  {
    // gsInfo << *geoContainer[el] << "\n";
    mpBSplineCurves.addPatch(*geoContainer[el]);
  }

  gsInfo << mpBSplineCurves << "\n";
  fd << mpBSplineCurves ;

  fd.dump("mp_output.xml");

  return EXIT_SUCCESS;
}
