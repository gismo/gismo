/** @file freeform_sdmatrix.cpp

    @brief Goes through all the valences in the 'patchpath' subfolder of
   filedata optimizes the fit with a target functional to create a legal
   coefficient space, then takes the kernel of that to generate the appropriate
   free form functionals and places them in the correct folder. Assumes to be
   run in the 'build'-folder and to reach 'filedata' from there via
   '../filedata'.

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
    std::string patchpath("freeform/");
    gsCmdLine cmd("Freeform subdivision");
    cmd.addString("p", "patchpath",
                  "The path to the files containing the model patches for EV "
                  "subdivision.",
                  patchpath);
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
    for (size_t valence = 3; valence < 10; ++valence)
    {
        if (valence == 4)
            continue;

        gsInfo << "=================\n    Valence " << valence
               << "\n=================\n\n";

        gsMatrix<real_t> coeffs(2 * valence + 1, 2 * valence + 1);

        // Iterate all functions for that valence
        for (size_t function = 1; function <= 2 * valence + 1; ++function)
        {
            gsInfo << "Function " << function << "\n";
            // Load the coefficients
            subdiv.initialize_data_xml(patchpath + "Val" +
                                       std::to_string(valence) + "Fct" +
                                       std::to_string(function) + ".xml");

            // Now smooth to get the desired linear combinations
            auto res = subdiv.smooth(1);
            // Collect all other coefficients in a matrix
            coeffs.row(function - 1) = res[0].transpose().row(2);
            gsInfo << "\n";
        }

        // Now coeffs contains a basis of the legal coefficient space.
        auto legal_pl = coeffs.fullPivLu();
        legal_pl.setThreshold(1e-8);
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
