/** @file gsParasolid_test.cpp

    @brief Test for the connection to the Siemens's Parasolid geometric kernel

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <iostream>

#include <gismo.h>

#include <gsParasolid/gsReadParasolid.h>
#include <gsParasolid/gsWriteParasolid.h>

using std::cout;
using std::endl;

using namespace gismo;
using extensions::gsWriteParasolid;


int main(int argc, char *argv[])
{
    std::string fn = "surfaces/simple.xml";

    gsCmdLine cmd("Hi, give me a file and I will read the contents to/from Parasolid.");
    cmd.addPlainString("filename", "G+SMO or Parasolid file", fn);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
    
    // Read in a surface
    cout << "Read in "<< gsFileManager::getFilename(fn) <<"\n";
    gsMultiPatch<> mp;
    gsReadFile<>(fn, mp);
    cout << mp <<"\n";
    
    // Get filename and extension
    std::string name = gsFileManager::getBasename (fn);
    std::string ext  = gsFileManager::getExtension(fn);

    if ( ext == "xml" )
    {
        // Write out in parasolid format

        gsWriteParasolid(mp, name);
        name += ".xmt_txt";
        cout << "Write out "<< name <<"\n";
    }
    else // if parasolid
    {
        // Write out in G+SMO format
        gsFileData<> fd;
        cout << "Write out "<< name <<".xml \n";
        fd<< mp;
        fd.dump(name);
    }
    
    return 0;
}
