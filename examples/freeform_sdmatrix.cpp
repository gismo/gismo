/** @file freeform_sdmatrix.cpp

    @brief Goes through all the valences in the 'patchpath' folder of filedata
   and creates the subdivision matrices for them. Also outputs the coefficient
   values for the outer points.

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
    index_t valence_max(9);
    gsCmdLine cmd("Freeform subdivision");
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

    // Iterate all valences
    for (size_t valence = 3; valence < size_t(valence_max); ++valence)
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

            // Subdivide once
            subdiv.subdivide();
            // Now smooth to get the desired matrices
            std::vector<gsMatrix<real_t>> ev_coefs;
            std::vector<gsMatrix<real_t>> ev_coefs_outer;
            subdiv.smooth(1, ev_coefs, ev_coefs_outer);
            // For the first one, save the outer functions
            if (function == 1)
            {
                gsWrite(ev_coefs_outer[0], "CoefficientsOuterVal" +
                                               std::to_string(valence) +
                                               ".xml");
            }
            // Collect all other coefficients in a matrix
            coeffs.row(function - 1) = ev_coefs[0].transpose().row(2);
            gsInfo << "\n";
        }

        // save all collected coefficients to file
        gsWrite(coeffs,
                "CoefficientsInnerVal" + std::to_string(valence) + ".xml");
    }

    return 0;
}
