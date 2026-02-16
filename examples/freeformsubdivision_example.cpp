/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include <gismo.h>
#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSurfMesh.h>

using namespace gismo;

int main(int argc, char** argv)
{
    // Command line
    std::string filedata("off/octtorus.off");
    std::string operations("sd");

    // Inputs
    gsCmdLine cmd("Freeform subdivision");
    cmd.addPlainString("filename", "File containing mesh.", filedata);
    cmd.addString("o", "operations", "Operations to perform on the mesh. Use d for subdivision and s for (c1) smoothing", operations);
    try
    {
        cmd.getValues(argc, argv);
    }
    catch (int errorcode)
    {
        return errorcode;
    }

    gsSurfMesh mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(filedata, mesh);

    auto subdiv = gsFreeformSubdivision<5, 3>(&mesh);

    subdiv.initialize_data();
    gsWriteParaview(subdiv.multipatch(), "results/initial_data");


    size_t i(1);
    for(char c : operations){
        switch (c) {
            case 'd':
                gsInfo << "Step " << std::string(i, 'I') << ": Subdividing.\n";
                subdiv.subdivide();
                break;
            case 's':
                gsInfo << "Step " << std::string(i, 'I') << ": Smoothing.\n";
                subdiv.smooth(1);
                break;
            default:
                break;
       }
        gsWriteParaview(subdiv.multipatch(),
                        "results/step" + std::string(i++, 'I'));
    }
}
