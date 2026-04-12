/** @file freeform_pde.cpp

    @brief Convergence study for solving the Laplace-Beltrami PDE on a
    freeform subdivision surface.

    For each refinement level \f$i = 0, 1, \ldots, \mathrm{steps}-1\f$, the
    example loads the input mesh, subdivides it \f$i\f$ times, and solves the
    Laplace-Beltrami equation
    \f[
        -\Delta_\Gamma u = f
    \f]
    on the surface \f$\Gamma\f$ using the \c laplace_beltrami method. The
    solution replaces the last coordinate of each patch.

    If the exact solution \f$u^*\f$ is known and supplied via \b --exact, the
    \f$L^\infty\f$ and \f$L^2\f$ approximation errors between the computed
    solution and \f$u^*\f$ are printed for every refinement level and
    written as a matrix to \c errors.csv.

    Optionally, each refinement level can also be exported to Paraview (\c .vts
    surface patches, \c .pvd collection, and – if an exact solution is given –
    a point-wise error field, and Bézier control net).

    \par Command-line arguments
    - \b -m / \b --mesh (default: \c freeform/flat/Val5Flat.xml): path to
      the input mesh \c .xml file (collection of \c gsTensorBSpline<2>
      patches), relative to \c filedata/.
    - \b -f / \b --function (default: \c "1"): analytic expression of the
      right-hand side function \f$f(x,y)\f$ in the Laplace-Beltrami PDE.
    - \b -e / \b --exact (default: empty): analytic expression of the exact
      (analytical) solution \f$u^*(x,y)\f$ of the PDE. If provided,
   \f$L^\infty\f$ and \f$L^2\f$ errors between the computed solution and
   \f$u^*\f$ are computed and printed at every refinement level.
    - \b -p / \b --patches (default: \c freeform/bubble/): path, relative to
      \c filedata/, to the directory containing the model patch files for
      extraordinary-vertex subdivision.
    - \b -s / \b --steps (default: \c 2): number of refinement levels to test
      (i.e. maximum number of subdivision steps).
    - \b -a / \b --samples (default: \c 10): number of sample points per
      parameter direction per patch face used when computing the error norms.
    - \b -v / \b --valence (default: \c -1, disabled): if set to a positive
      value, \b --mesh is ignored and the flat model patch for the given
      extraordinary valence is loaded directly from the patch directory.
    - \b --paraview: if set, writes Paraview output (\c results/pde_*.vts
      and \c results/pde.pvd) for every refinement level. If \b --exact is
      also set, a per-face error field is written as well.
    - \b --cnet: if set (together with \b --paraview), also writes the Bézier
      control net.
      WARNING: Currently unsuable due to a gismo bug.

    Author(s): L. Mussmaecher
*/

#include "gsCore/gsFunctionExpr.h"
#include <gismo.h>
#include <gsIO/gsCsv.h>
#include <gsIO/gsFileData.hpp>
#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSurfMesh.h>

using namespace gismo;

int main(int argc, char** argv)
{
    // Command line inputs
    std::string mesh_path("freeform/flat/Val5Flat.xml");
    std::string model_patch_path("freeform/bubble/");
    index_t steps(2);
    index_t valence(-1);
    index_t samples(10);
    std::string rhs_str("1");
    std::string exact_str("");
    bool control_net(false);
    bool paraview(false);

    gsCmdLine cmd("Freeform subdivision: Laplace-Beltrami PDE solver");
    cmd.addString("m", "mesh", "File containing the mesh.", mesh_path);
    cmd.addString("f", "function",
                  "The right-hand side function f(x,y) of the Laplace-Beltrami "
                  "PDE. E.g. `1` or `x^2 + y`. Defaults to `1`.",
                  rhs_str);
    cmd.addString(
        "e", "exact",
        "The exact (analytical) solution u*(x,y) of the PDE. If set, "
        "L-inf and L2 errors between the computed solution and u* are "
        "printed at each refinement level.",
        exact_str);
    cmd.addString("p", "patches",
                  "The path to the folder containing the model patches for EV "
                  "subdivision.",
                  model_patch_path);
    cmd.addInt("s", "steps",
               "The number of refinement levels (subdivide, solve PDE) to "
               "repeat.",
               steps);
    cmd.addInt("a", "samples",
               "The number of samples per patch direction (will be squared).",
               samples);
    cmd.addInt("v", "valence",
               "The valence of the patch. Overwrites the file path, if set.",
               valence);
    cmd.addSwitch("paraview", "Outputs the results to Paraview.", paraview);
    cmd.addSwitch("cnet", "Shows the control net of the patches.", control_net);

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
        mesh_path = "freeform/flat/Val" + std::to_string(valence) + "Flat.xml";
    }

    size_t MESH_DIM = 2;
    gsSurfMesh mesh = gsSurfMesh();
    auto subdiv = gsFreeformSubdivision<5>(&mesh, MESH_DIM);
    subdiv.options().setString("model_patch_path", model_patch_path);

    subdiv.initialize_data(mesh_path);

    gsFunctionExpr<real_t> rhs(rhs_str, 2);

    const bool has_exact = !exact_str.empty();
    // When no exact solution is given, initialise with a placeholder; the
    // object is never evaluated in that case.
    gsFunctionExpr<real_t> exact(has_exact ? exact_str : "0", 3);

    gsMatrix<real_t> errors(2, steps);

    gsParaviewCollection collection("results/pde");
    gsParaviewCollection cnet_collection("results/pde_cnet");
    gsParaviewCollection error_collection("results/pde_error");

    for (index_t i = 0; i < steps; ++i)
    {
        subdiv.initialize_data(mesh_path, MESH_DIM);
        for (index_t j = 0; j < i; ++j)
        {
            subdiv.subdivide();
        }
        // this is necessary if the loaded data is not smoothed already, e.g. in case of an OFF file.
        subdiv.smooth(1);
        subdiv.laplace_beltrami(rhs);

        if (has_exact)
        {
            errors.col(i) = subdiv.error(exact, samples);

            if (paraview)
            {
                subdiv.write_paraview_error(exact, errors(0, i),
                                            "results/pde_error" +
                                                std::to_string(i + 1),
                                            &error_collection, i + 1);
            }
        }

        if (paraview)
        {
            const std::string stepname =
                "results/pde_step" + std::to_string(i + 1);
            subdiv.write_paraview(stepname, &collection, &cnet_collection,
                                  i + 1, control_net);
        }

        gsInfo << "Finished Step " << std::to_string(i + 1) << ".\n";
    }

    if (paraview)
    {
        collection.save();
        if (has_exact)
            error_collection.save();
        if (control_net)
            cnet_collection.save();
    }

    if (has_exact)
        gsWriteCsv("results/errors.csv", errors);

    return 0;
}
