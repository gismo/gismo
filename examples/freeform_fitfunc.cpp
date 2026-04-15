/** @file freeform_fitfunc.cpp

    @brief Convergence study for fitting a scalar function to the last
    coordinate of a freeform subdivision surface.

    For each refinement level \f$i = 0, 1, \ldots, \mathrm{steps}-1\f$, the
    example loads the input mesh, subdivides it \f$i\f$ times, replaces the
    last coordinate of each patch by the best-fit Bézier approximation of a
    user-supplied function of the first two coordinates, smooths the result to
    \f$C^1\f$ using the precomputed extraordinary-vertex kernel data, and
    evaluates the \f$L^\infty\f$ and \f$L^2\f$ approximation errors. Error
    values for every refinement level are printed to \c stdout and optionally
    written as a matrix to \c errors.csv.

    Optionally, each refinement level can also be exported to Paraview (\c .vts
    surface patches, \c .pvd collection, point-wise error field, and Greville
    control points for extraordinary vertices together with their point indices
    as scalar point data).

    \par Command-line arguments
    - \b -m / \b --mesh (default: \c freeform/flat/Val5Flat.xml): path to
   the input mesh \c .xml file (collection of \c gsTensorBSpline<2> patches),
      relative to \c filedata/.
    - \b -f / \b --function (default: \c "x+y"): analytic expression of the
      target function \f$f(x,y)\f$ to be fitted to the last coordinate.
    - \b -p / \b --patches (default: \c freeform/bubble/): path, relative
      to \c filedata/, to the directory containing the model patch files for
      extraordinary-vertex subdivision together with the auxiliary
      \c Val\<v\>Kernel.xml files and, unless \b --opt is set, the
      \c Val\<v\>Constraints.xml files.
    - \b -s / \b --steps (default: \c 2): number of refinement levels to test
      (i.e. maximum number of subdivision steps).
    - \b -a / \b --samples (default: \c 10): number of sample points per
      parameter direction per patch face used when computing the error norms.
    - \b -v / \b --valence (default: \c -1, disabled): if set to a positive
      value, \b --mesh is ignored and the flat model patch for the given
      extraordinary valence is loaded directly from the patch directory.
    - \b --paraview: if set, writes Paraview output (\c results/fitfunc_*.vts
      and \c results/fitfunc.pvd) for every refinement level.
    - \b --cnet: if set (together with \b --paraview), also writes the Bézier
      control net and the EV Greville control points with their point indices.
    - \b --errors: if set, writes the \f$L^\infty\f$ and \f$L^2\f$ error
      column vectors to \c errors.csv.
    - \b --weighted: if set, weights the least-squares fit around
      extraordinary vertices using a per-sample weight vector loaded from \c
      filedata/freeform/val\<v\>_weights.xml.
    - \b --opt: if set, \c smooth() selects extraordinary-vertex coefficient
      representatives by minimising the diff functional used by
      \c fit_ev_opt(); otherwise it uses the linear functionals stored in
      \c Val\<v\>Constraints.xml.

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
    std::string function("x+y");
    bool control_net(false);
    bool optimize_fit(false);
    bool weighted_fit(false);
    bool paraview(false);
    bool write_errors(false);

    // Inputs
    gsCmdLine cmd("Freeform subdivision");
    cmd.addString("m", "mesh", "File containing the mesh.", mesh_path);
    cmd.addString("f", "function",
                  "A function to replace the last coordinate of your loaded "
                  "object. E.g. `x^2 + y`. Defaults to `x+y`",
                  function);
    cmd.addString("p", "patches",
                  "The path to the folder containing the model patches for EV "
                  "subdivision.",
                  model_patch_path);
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
    cmd.addSwitch("errors",
                  "Writes an error matrix to a file and outputs the error as a "
                  "paraview patch.",
                  write_errors);
    cmd.addSwitch(
        "opt",
        "During smoothing, selects EV coefficient representatives by "
        "minimizing the diff functional instead of using "
        "`Val<v>Constraints.xml`.",
        optimize_fit);
    cmd.addSwitch("weighted",
                  "Uses a weighting for the EV fit, loaded from a weights "
                  "vector at `filedata/freeform/val<v>_weights.xml`.",
                  weighted_fit);
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

    gsSurfMesh mesh = gsSurfMesh();
    auto subdiv = gsFreeformSubdivision<5>(&mesh, 2);
    subdiv.options().setString("model_patch_path", model_patch_path);
    subdiv.options().setSwitch("optimize_fit", optimize_fit);
    subdiv.options().setSwitch("weighted_fit", weighted_fit);

    subdiv.initialize_data(mesh_path);

    gsFunctionExpr<real_t> func(function, 2);

    gsMatrix<real_t> errors(2, steps);

    // Single .pvd collection for all steps.
    gsParaviewCollection collection("results/function_fit");
    gsParaviewCollection cnet_collection("results/function_cnet");
    gsParaviewCollection error_collection("results/function_error");

    for (index_t i = 0; i < steps; ++i)
    {
        subdiv.initialize_data(mesh_path, 2);
        for (index_t j = 0; j < i; ++j)
        {
            subdiv.subdivide();
        }
        subdiv.fit_function(func);

        if (write_errors)
            errors.col(i) = subdiv.error(func, samples);

        if (paraview)
        {

            const std::string stepname = "results/step" + std::to_string(i + 1);

            subdiv.write_paraview(stepname, &collection, &cnet_collection,
                                  i + 1, control_net);

            if (write_errors)
            {
                subdiv.write_paraview_error(
                    func, errors(0, i), "results/error" + std::to_string(i + 1),
                    &error_collection, i + 1);
            }

            // If the user wants to save the control net, we also want to save
            // the Greville control points for the first EV. Only works with set valence currently.
            if (control_net && valence > 0)
            {
                auto all_coefs = subdiv.smooth(1);
                // Here is also the spot to see more control points if wanted.
                gsMatrix<real_t> ev_coefs = all_coefs.topRows(2 * valence + 1).transpose();
                gsMatrix<real_t> indexed_ev_coefs(4, ev_coefs.cols());
                // for each EV
                gsInfo << "all_coefs: " << all_coefs.rows() << "x" << all_coefs.cols() << "\n";
                gsInfo << "ev_coefs: " << ev_coefs.rows() << "x" << ev_coefs.cols() << "\n";

                indexed_ev_coefs.topRows(3) = ev_coefs;
                for (index_t k = 0; k < ev_coefs.cols(); ++k)
                    indexed_ev_coefs(3, k) = (real_t)k;

                gsWriteParaviewPoints(indexed_ev_coefs, stepname + "_greville");
                // Register that file in the time series collection
                cnet_collection.addPart("step" + std::to_string(i + 1) +
                                            "_greville"
                                            ".vtp",
                                        i + 1, "Greville");
            }
        }

        gsInfo << "Finished Step " << std::to_string(i + 1) << ".\n";
    }

    if (paraview)
    {
        collection.save();
        if (write_errors)
            error_collection.save();
        if (control_net)
            cnet_collection.save();
    }

    // Write error matrix to errors.csv.
    if (write_errors)
        gsWriteCsv("results/errors.csv", errors);

    return 0;
}
