/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include "gsCore/gsFunctionExpr.h"
#include <gismo.h>
#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSurfMesh.h>

using namespace gismo;

int main(int argc, char** argv)
{
    // Command line
    std::string filepath("off/octtorus.off");
    std::string operations("sd");
    std::string function("");
    bool control_net(false);

    // Inputs
    gsCmdLine cmd("Freeform subdivision");
    cmd.addPlainString("filepath", "File containing mesh.", filepath);
    cmd.addString("o", "operations",
                  "Operations to perform on the mesh. Use d for subdivision "
                  "and s for (c1) smoothing",
                  operations);
    cmd.addString("f", "function",
                  "A function to replace the last coordinate of your loaded "
                  "object. E.g. `x^2 + y`.",
                  function);
    cmd.addSwitch("cnet", "Shows the control net of the patches.", control_net);
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

    // Take care of function
    if(function.length() > 0){
        subdiv.replace_last_coordinate_with_function(gsFunctionExpr<real_t>(function, 2));
    }

    gsWriteParaview(subdiv.multipatch(), "results/initial_data", 1000, false,
                    control_net);

    size_t i(1);
    for (char c : operations)
    {
        switch (c)
        {
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
                        "results/step" + std::string(i++, 'I'), 1000, false,
                        control_net);
    }

    return 0;
}
