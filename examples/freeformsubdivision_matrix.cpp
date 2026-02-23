/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include <gismo.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSurfMesh.h>
#include <fstream>
#include <string>

using namespace gismo;

int main(int argc, char** argv)
{

    gsSurfMesh mesh = gsSurfMesh();
    auto subdiv = gsFreeformSubdivision<5, 3>(&mesh);

    for (size_t valence = 3; valence < 10; ++valence)
    {
        if (valence == 4)
            continue;

        gsMatrix<real_t> coeffs(2 * valence + 1, 2 * valence + 1);

        for (size_t function = 1; function <= 2 * valence + 1; ++function)
        {
            subdiv.initialize_data_xml(
                "freeformSubdivision/fitting_functions/Val" +
                std::to_string(valence) + "Fct" + std::to_string(function) +
                ".xml");

            auto res = subdiv.smooth(1);

            coeffs.row(function - 1) = res[0].transpose().row(2);
        }

        // save all coefficients to file
        std::ofstream file("CoefficientsVal" + std::to_string(valence) + ".xml");
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<xml>\n"
             << "  <Matrix rows=\"" << coeffs.rows()
             << "\" cols=\"" << coeffs.cols() << "\" format=\"ascii\">\n";
        for (index_t r = 0; r < coeffs.rows(); ++r)
        {
            file << "    ";
            for (index_t c = 0; c < coeffs.cols(); ++c)
                file << coeffs(r, c) << (c + 1 < coeffs.cols() ? " " : "");
            file << "\n";
        }
        file << "  </Matrix>\n</xml>\n";
    }

    return 0;
}
