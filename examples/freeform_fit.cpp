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
    std::string filepath("off/octtorus.off");
    std::string patchpath("freeform/original/");
    std::string operations("sd");
    bool control_net(false);
    bool optimize_fit(false);

    // Inputs
    gsCmdLine cmd("Freeform subdivision");
    cmd.addPlainString("filepath", "File containing mesh.", filepath);
    cmd.addString("o", "operations",
                  "Operations to perform on the mesh. Use d for subdivision "
                  "and s for (c1) smoothing",
                  operations);
    cmd.addString("p", "patchpath",
                  "The path to the files containing the model patches for EV "
                  "subdivision.",
                  patchpath);
    cmd.addSwitch("cnet", "Shows the control net of the patches.", control_net);
    cmd.addSwitch(
        "opt",
        "Optimizes the fit via a functional instead of linear constraints.",
        optimize_fit);
    try
    {
        cmd.getValues(argc, argv);
    }
    catch (int errorcode)
    {
        return errorcode;
    }

    gsSurfMesh mesh = gsSurfMesh();
    auto subdiv = gsFreeformSubdivision<5, 3>(&mesh);

    subdiv.options().setString("model_patch_path", patchpath);
    subdiv.options().setSwitch("optimize_fit", optimize_fit);

    std::string xml(".xml");
    std::string off(".off");
    // Check the filetype to be loaded.
    if (std::equal(filepath.begin() + filepath.size() - xml.size(),
                   filepath.end(), xml.begin()))
    {
        gsInfo << "Loading xml\n";
        subdiv.initialize_data_xml(filepath);
    }
    else if (std::equal(filepath.begin() + filepath.size() - off.size(),
                        filepath.end(), off.begin()))
    {
        gsInfo << "Loading off\n";
        subdiv.initialize_data_off(filepath);
    }
    else
    {
        gsWarn << "Unsupported Filetype!\n";
        return 1;
    }

    // Single .pvd collection for all steps.
    gsParaviewCollection collection("results/function_fit");

    subdiv.write_paraview("results/initial_data", &collection, 0, control_net);

    size_t i(1);
    for (char c : operations)
    {
        switch (c)
        {
        case 'd':
            gsInfo << "Step " << std::to_string(i) << ": Subdividing.\n";
            subdiv.subdivide();
            break;
        case 's':
            gsInfo << "Step " << std::to_string(i) << ": Smoothing.\n";
            subdiv.smooth(1);
            break;
        default:
            break;
        }

        subdiv.write_paraview("results/step" + std::to_string(i),
                              &collection, i, control_net);
        ++i;
    }

    collection.save();

    return 0;
}
