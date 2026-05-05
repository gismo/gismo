/** @file frog_pde.cpp

    @brief Convergence study for solving the Laplace-Beltrami PDE on a
    frog subdivision surface.

    For each refinement level \f$i = 0, 1, \ldots, \mathrm{steps}-1\f$, the
    example loads the input mesh, subdivides it \f$i\f$ times, and solves the
    Laplace-Beltrami equation
    \f[
        -\Delta_\Gamma u = f
    \f]
    on the surface \f$\Gamma\f$ using the \c laplace_beltrami method. The
    solution replaces the last coordinate of each patch.

    If the exact solution \f$u^*\f$ is known and supplied at ID 2 of the
    input XML file, the \f$L^\infty\f$ and \f$L^2\f$ approximation errors
    between the computed solution and \f$u^*\f$ are printed for every
    refinement level and written as a matrix to \c errors.csv.

    Optionally, each refinement level can also be exported to Paraview (\c .vts
    surface patches, \c .pvd collection, and – if an exact solution is given –
    a point-wise error field, and Bézier control net).

    \par Command-line arguments
    - \b -m / \b --mesh (default: \c frog/flat/Val5Flat.xml): path to
      the input mesh \c .xml file (collection of \c gsTensorBSpline<2>
      patches), relative to \c filedata/. Ignored if \c valence is set in
      the input options (ID 0).
    - \b -d / \b --mesh-dim (default: \c 3): spatial dimension of the mesh
      (number of coordinates per vertex). Automatically set to \c 2 when
      loading a flat mesh via the \c valence option.
    - \b -f / \b --function: right-hand side function string (e.g. \c "x+y").
      Used as the RHS when \b --input is not provided. At least one of
      \b --input or \b --function must be supplied.
    - \b -i / \b --input: path to an XML file containing PDE inputs. If not
      provided, \b --function is used as the RHS instead. Contents:
      - ID 0 (\c OptionList, optional): solver options:
        - \c valence (\c index_t): if positive, overrides \b --mesh with
          the flat \c ValNFlatStraight.xml model for the given extraordinary
          valence and sets \b --mesh-dim to \c 2.
        - \c x_scale (\c real_t): multiplicative scaling applied in
          \f$x\f$.
        - \c y_scale (\c real_t): multiplicative scaling applied in
          \f$y\f$.
      - ID 1 (\c Function, optional): the right-hand side \f$f(x,y)\f$ of
        the Laplace-Beltrami PDE. Defaults to \c "1" if absent.
      - ID 2 (\c Function, optional): the exact solution \f$u^*(x,y)\f$.
        If present, \f$L^\infty\f$ and \f$L^2\f$ errors are printed and
        written to \c errors.csv.
    - \b -p / \b --frogdir (default: \c frog/bubble/): path, relative to
      \c filedata/, to the directory containing a set of frog spline
      generating functions for each required valence.
    - \b -s / \b --steps (default: \c 2): number of refinement levels
      (subdivide-then-solve iterations).
    - \b --paraview: if set, writes Paraview output (\c results/pde_*.vts
      and \c results/pde.pvd) for every refinement level. If an exact solution
      is provided, a per-face error field is written as well.
    - \b --cnet: if set (together with \b --paraview), also writes the Bézier
      control net.
      WARNING: Currently unusable due to a gismo bug.

    Author(s): L. Mussmaecher
*/

#include "gsCore/gsDebug.h"
#include "gsCore/gsFunctionExpr.h"
#include "gsIO/gsFileData.h"
#include "gsIO/gsOptionList.h"
#include "gsMatrix/gsVector.h"
#include <gismo.h>
#include <gsIO/gsCsv.h>
#include <gsIO/gsFileData.hpp>
#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsFrogSplines.h>
#include <gsMesh2/gsSurfMesh.h>

using namespace gismo;

int main(int argc, char** argv)
{
    // Command line inputs
    std::string mesh_path("frog/flat/Val5Flat.xml");
    index_t valence(-1);
    index_t mesh_dim(3);
    std::string function_str("");
    std::string input_path("");
    std::string model_patch_path("frog/bubble/");
    index_t steps(2);
    bool control_net(false);
    bool paraview(false);

    gsCmdLine cmd("Frog subdivision: Laplace-Beltrami PDE solver");

    // Function/Mesh Input
    cmd.addString("m", "mesh", "File containing the mesh.", mesh_path);
    cmd.addInt("d", "mesh-dim", "Dimension of the mesh. 3 is the default.",
               mesh_dim);
    cmd.addString("f", "function", "Right-hand side of the Laplace PDE.",
                  function_str);

    // Alternative: Input file
    cmd.addString("i", "input",
                  "Path to an XML file containing PDE inputs (OptionList at "
                  "ID 0, RHS function at ID 1, exact solution at ID 2).",
                  input_path);
    cmd.addString("p", "frogdir",
                  "The path to the folder containing a set of frog spline "
                  "generating functions for each required valence.",
                  model_patch_path);

    // Process modification
    cmd.addInt("s", "steps",
               "The number of refinement levels (subdivide, solve PDE) to "
               "repeat.",
               steps);

    // Output modification
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

    // Load input data
    gsFunctionExpr<real_t> rhs, exact;
    bool has_exact(false);
    real_t x_scale(1.0), y_scale(1.0);
    // First, try to get the data from the file
    if (!input_path.empty())
    {
        // The target function
        gsFileData<> fd(input_path);
        if (fd.hasId(1))
        {
            fd.getId(1, rhs);
        }

        // The exact known solution
        if (fd.hasId(2))
        {
            has_exact = true;
            fd.getId(2, exact);
        }

        // The options
        if (fd.hasId(0))
        {
            gsOptionList opts;
            fd.getId(0, opts);
            valence = opts.getInt("valence");
            x_scale = opts.getReal("x_scale");
            y_scale = opts.getReal("y_scale");
        }
    }
    else
    {
        rhs = gsFunctionExpr<>(function_str, mesh_dim);
    }

    // if the valence is set, mesh_dim and mesh_path are overwritten.
    if (valence > 0)
    {
        mesh_path =
            "frog/flat/Val" + std::to_string(valence) + "FlatStraight.xml";
        mesh_dim = 2;
    }

    gsVector<> scale(3);
    scale << x_scale, y_scale, 1.0;

    gsSurfMesh mesh = gsSurfMesh();
    auto subdiv = gsFrogSplines<5>(&mesh, mesh_dim);
    subdiv.options().setString("frog_dir", model_patch_path);

    gsMatrix<real_t> errors(2, steps);

    gsParaviewCollection collection("results/pde");
    gsParaviewCollection cnet_collection("results/pde_cnet");
    gsParaviewCollection error_collection("results/pde_error");

    for (index_t i = 0; i < steps; ++i)
    {
        subdiv.initialize_data(mesh_path, mesh_dim);
        for (index_t j = 0; j < i; ++j)
        {
            subdiv.subdivide();
        }
        // this is necessary if the loaded data is not smoothed already, e.g. in
        // case of an OFF file.
        subdiv.scale(scale);
        subdiv.smooth(1);
        subdiv.laplace_beltrami(rhs);

        if (has_exact)
        {
            errors.col(i) = subdiv.error(exact, 10);

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

    // The `build/errors/` directory needs to be manually created for this to work!
    if (has_exact)
        gsWriteCsv("errors/" + input_path.substr(input_path.find_last_of("/\\") + 1) +
                       "_errors.csv",
                   errors);

    return 0;
}
