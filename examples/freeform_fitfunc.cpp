/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include "gsCore/gsFunctionExpr.h"
#include <gismo.h>
#include <gsIO/gsFileData.hpp>
#include <gsIO/gsCsv.h>
#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSurfMesh.h>

using namespace gismo;

int main(int argc, char** argv)
{
    // Command line inputs
    std::string filepath("freeform/original/Val5Flat.xml");
    std::string patchpath("freeform/original/");
    index_t steps(2);
    index_t valence(-1);
    index_t samples(10);
    std::string function("x+y");
    bool control_net(false);
    bool optimize_fit(false);
    bool paraview(false);

    // Inputs
    gsCmdLine cmd("Freeform subdivision");
    cmd.addPlainString("filepath", "File containing mesh.", filepath);
    cmd.addString("f", "function",
                  "A function to replace the last coordinate of your loaded "
                  "object. E.g. `x^2 + y`. Defaults to `x+y`",
                  function);
    cmd.addString("p", "patchpath",
                  "The path to the files containing the model patches for EV "
                  "subdivision.",
                  patchpath);
    cmd.addInt("s", "steps",
               "The number of steps (subdivide, fit, smooth) to repeat.",
               steps);
    cmd.addInt("a", "samples",
               "The number of samples on each patch (will be squared).",
               samples);
    cmd.addInt(
        "v", "valence",
        "The valence of the patch to fit. Overwrites the file path, if set.",
        valence);
    cmd.addSwitch("paraview", "Outputs the fits to paraview.", paraview);
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

    if (valence > 0)
    {
        filepath =
            "freeform/original/Val" + std::to_string(valence) + "Flat.xml";
    }

    gsSurfMesh mesh = gsSurfMesh();
    auto subdiv = gsFreeformSubdivision<5, 3>(&mesh);
    subdiv.options().setString("model_patch_path", patchpath);
    subdiv.options().setSwitch("optimize_fit", optimize_fit);

    subdiv.initialize_data_xml(filepath);

    gsFunctionExpr<real_t> func(function, 2);

    gsMatrix<real_t> errors(3, steps);

    for (index_t i = 0; i < steps; ++i)
    {
        if (i > 0)
            subdiv.subdivide();

        subdiv.fit_last_coordinate_to_function(func);
        subdiv.smooth(1);
        errors.col(i) = subdiv.error(func, samples);

        if (paraview)
        {
            gsWriteParaview(subdiv.multipatch(),
                            "results/step" + std::string(i + 1, 'I'), 1000,
                            false, control_net);
        }

        gsInfo << "Finished writing Step " << std::string(i + 1, 'I') << ".\n";
    }

    // write error matrix
    gsWriteCsv("errors.csv", errors);

    return 0;
}
