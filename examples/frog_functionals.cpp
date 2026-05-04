/** @file frog_functionals.cpp

    @brief Precomputes the linear functionals that select constrained
    extraordinary-vertex coefficient representatives.

    For each valence \f$v\f$ from 3 up to \c valence_max (skipping \f$v=4\f$),
    the example loads all extraordinary-vertex blending functions from the
    collated file \c Val\<v\>Fcts.xml (one \c gsMultiPatch per function),
    initialises the frog-spline mesh for each function via
    \c initialize_data_multipatch, applies one round of \f$C^1\f$ smoothing
    with functional optimisation (\c optimize_fit = \c true), reads back the
    resulting extraordinary-vertex (EV) coefficient representatives, and
    assembles the space of coefficient vectors selected by that functional.
    The kernel of this space consists of the linear functionals that vanish on
    all legal representatives and is written to \c Val\<v\>Constraints.xml in
    the model patch directory.

    The number of blending functions is determined at runtime from the number
    of \c gsMultiPatch entries in \c Val\<v\>Fcts.xml and is not assumed to
    equal \f$2v+1\f$.

    The example assumes that \c Val\<v\>Kernel.xml has already been generated,
    e.g. by running \c frog_kernels.cpp first, because
    \c gsFrogSplines::smooth() uses those kernel bases internally.

    These functional files are subsequently consumed by
    \c gsFrogSplines::smooth() when \c optimize_fit is disabled and by
    \c fit_ev() during constrained EV fitting.

    \note The example assumes it is run from the build directory so that
    \c ../filedata is reachable. The patch-path option must point to a
    subdirectory of \c filedata.

    \par Command-line arguments
    - \b -p / \b --patchpath (default: \c frog/bubble/): path, relative
      to \c filedata/, to the directory that contains the collated model patch
      file \c Val\<v\>Fcts.xml.
    - \b -v / \b --valence (default: \c 9): maximum extraordinary-vertex
      valence to process. Functional matrices are generated for all valences
      \f$v \in \{3, 5, 6, \ldots, \mathrm{valence\_max}\}\f$ (valence 4 is
      regular and is skipped).

    Author(s): L. Mussmaecher
*/

#include <gismo.h>
#include <gsIO/gsFileData.hpp>
#include <gsMesh2/gsFrogSplines.h>
#include <gsMesh2/gsSurfMesh.h>
#include <string>

using namespace gismo;

int main(int argc, char** argv)
{
    // CMD arguments
    std::string patchpath("frog/bubble/");
    gsCmdLine cmd("Frog subdivision");
    index_t valence_max(6);
    cmd.addString("p", "patchpath",
                  "The path to the files containing the model patches for EV "
                  "subdivision.",
                  patchpath);
    cmd.addInt("v", "valence", "Maximal valence to calculate.", valence_max);
    try
    {
        cmd.getValues(argc, argv);
    }
    catch (int errorcode)
    {
        return errorcode;
    }

    // Basic objects
    gsSurfMesh mesh = gsSurfMesh();
    auto subdiv = gsFrogSplines<5>(&mesh, 3);
    subdiv.options().setString("model_patch_path", patchpath);
    subdiv.options().setSwitch("optimize_fit", true);

    // Iterate all valences
    for (size_t valence = 3; valence <= size_t(valence_max); ++valence)
    {
        if (valence == 4)
            continue;

        gsInfo << "=================\n    Valence " << valence
               << "\n=================\n\n";

        gsFileData<real_t> fd(patchpath + "Val" + std::to_string(valence) +
                              "Fcts.xml");
        auto fcts = fd.getAll<gsMultiPatch<real_t>>();

        gsMatrix<real_t> coeffs(fcts.size(), fcts.size());

        // Iterate all functions for that valence
        for (size_t function = 0; function < fcts.size(); ++function) 
        {
            gsInfo << "Function " << function << "\n";
            // Load the basis function patch file.
            subdiv.initialize_data_multipatch(*fcts[function]);

            // Apply functional-optimised smoothing to obtain the preferred EV
            // coefficient representative in the precomputed kernel family.
            gsMatrix<real_t> ev_coefs = subdiv.smooth(1);
            gsInfo << "ev_coefs: " << ev_coefs.rows() << "x" << ev_coefs.cols() << "\n";
            // Collect the selected EV coefficient representative for this
            // blending function.
            coeffs.row(function) = ev_coefs.topRows(fcts.size()).transpose().row(2);
            gsInfo << ev_coefs.topRows(fcts.size()).transpose().row(2) << "\n";
        }

        // Now `coeffs` spans the representative EV coefficient space selected
        // by the diff functional.
        auto legal_pl = coeffs.fullPivLu();
        legal_pl.setThreshold(1e-4);
        gsMatrix<> K = legal_pl.kernel().transpose();
        gsWrite(K, "../filedata/" + patchpath + "Val" +
                       std::to_string(valence) + "Constraints.xml");

        gsInfo << "Written representative-selection functionals for valence "
               << valence
               << " to `"
                << (gsFileManager::findInDataDir("") + patchpath + "Val" +
                    std::to_string(valence) + "Constraints.xml")
               << "`.\n";
    }

    return 0;
}
