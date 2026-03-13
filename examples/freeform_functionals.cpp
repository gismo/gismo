/** @file freeform_functionals.cpp

    @brief Precomputes the kernel-space constraint matrices used by the
    \f$C^1\f$ extraordinary-vertex fitting step of the freeform subdivision
    scheme.

    For each valence \f$v\f$ from 3 up to \c valence_max (skipping \f$v=4\f$),
    the example loads all \f$2v+1\f$ basis-function patches from the model
    patch directory (\c Val\<v\>Fct\<f\>.xml), applies one round of \f$C^1\f$
    smoothing with functional optimisation (\c optimize_fit = \c true), reads
    the resulting row of extraordinary-vertex (EV) coefficients, assembles the
    full \f$(2v+1)\times 12\f$ coefficient matrix, and computes its kernel.
    The kernel is saved as \c Val\<v\>Constraints.xml inside the model patch
    directory.

    These constraint files are subsequently consumed by \c fit_ev when
    performing \f$C^1\f$ smooth surface fitting around extraordinary vertices.

    \note The example assumes it is run from the build directory so that
    \c ../filedata is reachable. The patch-path option must point to a
    subdirectory of \c filedata.

    \par Command-line arguments
    - \b -p / \b --patchpath (default: \c freeform/original/): path, relative
      to \c filedata/, to the directory that contains the model patch files
      \c Val\<v\>Fct\<f\>.xml.
    - \b -v / \b --valence (default: \c 9): maximum extraordinary-vertex
      valence to process. Constraint matrices are generated for all valences
      \f$v \in \{3, 5, 6, \ldots, \mathrm{valence\_max}\}\f$ (valence 4 is
      regular and is skipped).

    Author(s): L. Mussmaecher
*/

#include <gismo.h>
#include <gsIO/gsFileData.hpp>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSurfMesh.h>
#include <string>

using namespace gismo;

int main(int argc, char** argv)
{
    // CMD arguments
    std::string patchpath("freeform/original/");
    gsCmdLine cmd("Freeform subdivision");
    index_t valence_max(9);
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
    auto subdiv = gsFreeformSubdivision<5, 3>(&mesh);
    subdiv.options().setString("model_patch_path", patchpath);
    subdiv.options().setSwitch("optimize_fit", true);

    // Iterate all valences
    for (size_t valence = 3; valence < size_t(valence_max); ++valence)
    {
        if (valence == 4)
            continue;

        gsInfo << "=================\n    Valence " << valence
               << "\n=================\n\n";

        gsMatrix<real_t> coeffs(2 * valence + 1, 2 * valence + 1);
        std::vector<gsMatrix<real_t>> ev_coefs;
        std::vector<gsMatrix<real_t>> ev_coefs_outer;

        // Iterate all functions for that valence
        for (size_t function = 1; function <= 2 * valence + 1; ++function)
        {
            gsInfo << "Function " << function << "\n";
            // Load the basis function patch file.
            subdiv.initialize_data_xml(patchpath + "Val" +
                                       std::to_string(valence) + "Fct" +
                                       std::to_string(function) + ".xml");

            // Apply functional-optimised smoothing to obtain the C1-constrained EV coefficients.
            subdiv.smooth(1, ev_coefs, ev_coefs_outer);
            // Collect all other coefficients in a matrix
            coeffs.row(function - 1) = ev_coefs[0].transpose().row(2);
            for(int i = 0; i < 3; i++)
                gsInfo << ev_coefs[0].transpose().row(i) << "\n";
            gsInfo << "\n";
        }

        // Now `coeffs` contains a basis for the space of C1-compatible EV coefficient vectors.
        auto legal_pl = coeffs.fullPivLu();
        legal_pl.setThreshold(1e-4);
        gsMatrix<> K = legal_pl.kernel().transpose();
        gsWrite(K, "../filedata/" + patchpath + "Val" +
                       std::to_string(valence) + "Constraints.xml");

        gsInfo << "Written functional constraints for valence " << valence
               << " to `"
               << (gsFileManager::findInDataDir("") + patchpath + "Val" +
                   std::to_string(valence) + "Constraints.xml")
               << "`.\n";
    }

    return 0;
}
