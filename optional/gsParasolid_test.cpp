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
    std::string fn;
    try 
    {
        gsCmdLine cmd("Hi, give me a file and I will read the contents to/from Parasolid.");    
        gsArgValPlain<std::string> arg_file("filename","G+SMO or Parasolid file", 
                                  false, "", "string",cmd );
        cmd.parse(argc,argv);
        fn = arg_file.getValue();

        if (fn.empty() )
        {
            fn = GISMO_SOURCE_DIR "/filedata/surfaces/simple.xml";
        }
        
    } catch ( gsArgException& e )
    { cout << "Error: " << e.error() << " " << e.argId() << endl; }
    
    // Read in a surface
    cout << "Read in "<< gsFileData<>::getFilename(fn) <<"\n";
    gsMultiPatch<> * mp = gsReadFile<>(fn);
    cout << *mp <<"\n";
    
    // Get filename and extension
    std::string name = gsFileData<>::getBasename (fn);
    std::string ext  = gsFileData<>::getExtension(fn);

    if ( ext == "xml" )
    {
        // Write out in parasolid format

        gsWriteParasolid(*mp, name);        
        name += ".xmt_txt";
        cout << "Write out "<< name <<"\n";
    }
    else // if parasolid
    {
        // Write out in G+SMO format
        gsFileData<> fd;
        cout << "Write out "<< name <<".xml \n";
        fd<< *mp;
        fd.dump(name);
    }

    delete mp;
    
    return 0;
}
