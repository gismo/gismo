/** @file frog_showbasis.cpp

    @brief Exports the frog spline generating functions for a given valence
    to Paraview.

    For the given valence \f$v\f$, loads all blending functions from
    \c Val\<v\>Fcts.xml (one \c gsMultiPatch per function) and exports each
    one as a Bézier multipatch to a Paraview time-series collection
    \c results/basis.pvd.  The function index is used as the timestep so that
    individual basis functions can be inspected by scrubbing the Paraview
    timeline.

    Optionally, the Bézier control net of each basis function is exported as
    well and registered in \c results/basis_cnet.pvd.

    \par Command-line arguments
    - \b -p / \b --frogdir (default: \c frog/bubble/): path, relative
      to \c filedata/, to the directory containing a set of frog spline
      generating functions for each required valence.
    - \b -v / \b --valence (default: \c 5): extraordinary-vertex valence
      whose basis functions are to be exported.
    - \b --cnet: if set, also exports the Bézier control net of each basis
      function and registers it in \c results/basis_cnet.pvd.

    Author(s): L. Mussmaecher
*/

#include <gismo.h>
#include <gsIO/gsFileData.hpp>
#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsFrogSplines.h>
#include <gsMesh2/gsSurfMesh.h>

using namespace gismo;

int main(int argc, char** argv)
{
    // Command-line arguments
    std::string frogdir("frog/bubble/");
    index_t valence(5);
    bool control_net(false);

    gsCmdLine cmd("Exports the frog spline generating functions for a given "
                  "valence to Paraview.");
    cmd.addString("p", "frogdir",
                  "The path to the folder containing a set of frog spline "
                  "generating functions for each required valence.",
                  frogdir);
    cmd.addInt("v", "valence",
               "The extraordinary-vertex valence whose basis functions are "
               "to be exported.",
               valence);
    cmd.addSwitch("cnet",
                  "Also exports the Bézier control net of each basis function "
                  "to `results/basis_cnet.pvd`.",
                  control_net);
    try
    {
        cmd.getValues(argc, argv);
    }
    catch (int errorcode)
    {
        return errorcode;
    }

    // Load the list of all basis function patches for the given valence.
    gsFileData<real_t> fd(frogdir + "Val" + std::to_string(valence) +
                          "Fcts.xml");
    auto fcts = fd.getAll<gsMultiPatch<real_t>>();
    gsInfo << "Found " << fcts.size() << " basis functions for valence "
           << valence << ".\n\n";

    // Set up the frog spline object.
    gsSurfMesh mesh;
    auto subdiv = gsFrogSplines<5>(&mesh, 3);
    subdiv.options().setString("frog_dir", frogdir);

    // Collections for the Paraview time-series output.
    gsParaviewCollection collection("results/basis");
    gsParaviewCollection cnet_collection("results/basis_cnet");

    for (size_t function = 0; function < fcts.size(); ++function)
    {
        gsInfo << "Exporting basis function " << function << " / "
               << (fcts.size() - 1) << "...\n";

        subdiv.initialize_data_multipatch(*fcts[function]);

        const std::string name = "results/basis_fct" + std::to_string(function);

        subdiv.write_paraview(name, &collection, &cnet_collection, function + 1,
                              control_net);
    }

    collection.save();
    if (control_net)
        cnet_collection.save();

    return 0;
}
