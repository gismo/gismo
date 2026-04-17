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

    If the exact solution \f$u^*\f$ is known and supplied at ID 1 of the
    functions XML file, the \f$L^\infty\f$ and \f$L^2\f$ approximation errors
    between the computed solution and \f$u^*\f$ are printed for every
    refinement level and written as a matrix to \c errors.csv.

    Optionally, each refinement level can also be exported to Paraview (\c .vts
    surface patches, \c .pvd collection, and – if an exact solution is given –
    a point-wise error field, and Bézier control net).

    \par Command-line arguments
    - \b -m / \b --mesh (default: \c freeform/flat/Val5Flat.xml): path to
      the input mesh \c .xml file (collection of \c gsTensorBSpline<2>
      patches), relative to \c filedata/.
    - \b -f / \b --functions (default: \c input_functions.xml): path to an
      XML file containing the functions used in the PDE:
      - ID 0: the right-hand side function \f$f(x,y)\f$ of the Laplace-Beltrami
        PDE (e.g. \c "1" or \c "x^2+y"). Defaults to \c "1" if absent.
      - ID 1 (optional): the exact (analytical) solution \f$u^*(x,y)\f$. If
        present, \f$L^\infty\f$ and \f$L^2\f$ errors between the computed
        solution and \f$u^*\f$ are computed and printed at every refinement
        level.
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
      and \c results/pde.pvd) for every refinement level. If an exact solution
      is provided in the functions file, a per-face error field is written as well.
    - \b --cnet: if set (together with \b --paraview), also writes the Bézier
      control net.
      WARNING: Currently unsuable due to a gismo bug.

    Author(s): L. Mussmaecher
*/

#include "gsCore/gsFunctionExpr.h"
#include "gsIO/gsFileData.h"
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
    index_t samples(10);
    index_t valence(-1);
    index_t mesh_dim(3);
    std::string functions_path("input_functions.xml");
    bool control_net(false);
    bool paraview(false);

    gsCmdLine cmd("Freeform subdivision: Laplace-Beltrami PDE solver");
    cmd.addString("m", "mesh", "File containing the mesh.", mesh_path);
    cmd.addString(
        "f", "functions",
        "Path to an XML file containing the functions used in the PDE. "
        "ID 0: right-hand side function f(x,y) of the Laplace-Beltrami PDE "
        "(e.g. `1` or `x^2 + y`); defaults to `1` if absent. "
        "ID 1 (optional): exact (analytical) solution u*(x,y); if present, "
        "L-inf and L2 errors between the computed solution and u* are "
        "printed at each refinement level.",
        functions_path);
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
    cmd.addInt("d", "mesh-dim",
               "Dimension of the mesh. 3 is the default.",
               mesh_dim);
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

    gsSurfMesh mesh = gsSurfMesh();
    auto subdiv = gsFreeformSubdivision<5>(&mesh, mesh_dim);
    subdiv.options().setString("model_patch_path", model_patch_path);

    subdiv.initialize_data(mesh_path);

    gsFileData<> fd(functions_path);
    gsFunctionExpr<real_t> rhs, exact;
    if(fd.hasId(0)){
        fd.getId(0, rhs);
    } else {
        rhs = gsFunctionExpr<>("1", mesh_dim);
    }

    const bool has_exact = fd.hasId(1);
    if(has_exact){
        fd.getId(1, exact);
    }
    

    // When no exact solution is given, initialise with a placeholder; the
    // object is never evaluated in that case.

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
