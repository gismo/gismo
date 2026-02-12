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
    bool no_smooth(false);
    index_t steps(1);

    // Inputs
    gsCmdLine cmd("Freeform subdivision");
    cmd.addPlainString("filename", "File containing mesh.", filedata);
    cmd.addSwitch("no_smooth", "C1 smoothing before subdivision.", no_smooth);
    cmd.addInt("s", "steps", "Number of subdivision steps.", steps);
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

    auto subdiv = gsFreeformSubdivision<5, 3>();

    subdiv.initialize_data(mesh);
    gsWriteParaview(subdiv.multipatch(mesh), "results/initial_data");

    if (!no_smooth)
    {
        subdiv.smooth(mesh, 1);
        gsWriteParaview(subdiv.multipatch(mesh), "results/c1");
    }

    for (index_t i = 0; i < steps; ++i)
    {
        subdiv.subdivide(mesh);
        subdiv.smooth(mesh, 1);

        gsWriteParaview(subdiv.multipatch(mesh),
                        "results/subdiv" + std::string(i, 'a'));
    }
}
