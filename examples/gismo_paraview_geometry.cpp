#include <gismo.h>
using namespace gismo;

int main(int argc, char* argv[])
{
    std::string filename("domain2d/two_patches_from_gluing_data.xml");

    gsCmdLine cmd("Example for get gluing data.");
    cmd.addString("f", "file", "G+Smo input multi patch file.", filename);

    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    gsMultiPatch<>::uPtr mp = gsReadFile<>(filename);
    if (!mp) { gsInfo << "Failed to read file.\n"; return -1; }
    
    gsInfo << "Read " << mp->nPatches() << " patches.\n";
    
    // Extract base name from XML filename (remove directory and extension)
    std::string basename = filename;
    size_t lastSlash = basename.find_last_of("/\\");
    if (lastSlash != std::string::npos)
        basename = basename.substr(lastSlash + 1);
    size_t lastDot = basename.find_last_of(".");
    if (lastDot != std::string::npos)
        basename = basename.substr(0, lastDot);
    
    // Export to ParaView format
    gsWriteParaview(*mp, basename, 100);
    gsInfo << "Written to " << basename << ".pvd\n";
    gsInfo << "\nUsage examples:\n";
    gsInfo << "  Default:     ./bin/gismo_paraview_geometry\n";
    gsInfo << "  Custom file: ./bin/gismo_paraview_geometry -f domain2d/two_bicubic_patches.xml\n";
    gsInfo << "  View in ParaView: paraview " << basename << ".pvd\n";
    
    return 0;
}