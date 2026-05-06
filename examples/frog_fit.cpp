/** @file frog_fit.cpp

    @brief Loads a quad mesh and performs a user-specified sequence of
    frog subdivision and \f$C^1\f$ smoothing steps, writing each
    intermediate result to Paraview.

    The mesh is loaded from a \c .off or \c .xml file. The operation
    sequence is given as a string of characters: \c d performs one
    subdivision step and \c s performs one \f$C^1\f$ smoothing step.
    Every step (including the initial mesh) is written to a Paraview
    \c .vts file and registered in a single \c .pvd time-series
    collection at \c results/function_fit.pvd.

    \par Command-line arguments
    - \c mesh (positional, default: \c off/octtorus.off): path to the
      input mesh file. Supported formats: \c .off (3D point data) and
      \c .xml (collection of \c gsTensorBSpline<2> patches).
    - \b -o / \b --operations (default: \c "sd"): sequence of operations
      to perform. Each character is executed in order:
      - \c d — one subdivision step (\c subdivide()).
      - \c s — one \f$C^1\f$ smoothing step (\c smooth(1)).
    - \b -p / \b --frogdir (default: \c frog/bubble/): path to the
      directory containing a set of frog spline generating functions for
      each required valence.
    - \b --cnet: if set, also writes the Bézier control net alongside the
      patch surface at every step.
    - \b --opt: if set, uses kernel-space functional optimisation
      (\c fit_ev_opt) instead of file-loaded linear constraints
      (\c fit_ev) when fitting around extraordinary vertices.

    Author(s): L. Mussmaecher
*/

#include <gismo.h>
#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsFrogSplines.h>
#include <gsMesh2/gsSurfMesh.h>

using namespace gismo;

int main(int argc, char** argv)
{
    // Command line
    std::string filepath("off/octtorus.off");
    std::string frogdir("frog/bubble/");
    std::string operations("sd");
    bool control_net(false);
    bool optimize_fit(false);

    // Inputs
    gsCmdLine cmd("frog subdivision");
    cmd.addPlainString("mesh", "Path to the file containing the mesh.",
                       filepath);
    cmd.addString("o", "operations",
                  "Operations to perform on the mesh. Use d for subdivision "
                  "and s for (c1) smoothing",
                  operations);
    cmd.addString("p", "frogdir",
                  "The path to the folder containing a set of frog spline "
                  "generating functions for each required valence.",
                  frogdir);
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
    auto subdiv = gsFrogSplines<5>(&mesh, 3);

    subdiv.options().setString("frog_dir", frogdir);
    subdiv.options().setSwitch("optimize_fit", optimize_fit);

    subdiv.initialize_data(filepath, 3);

    // Single .pvd collection for all steps.
    gsParaviewCollection collection("results/fit");
    gsParaviewCollection cnet_collection("results/cnet");

    subdiv.write_paraview("results/initial_data", &collection, &cnet_collection,
                          0, control_net);

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
            // subdiv.smooth(1);
            subdiv.smooth(1);
            break;
        default:
            break;
        }

        subdiv.write_paraview("results/step" + std::to_string(i), &collection,
                              &cnet_collection, i, control_net);
        ++i;
    }

    collection.save();
    if (control_net)
        cnet_collection.save();

    return 0;
}
