/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include <fstream>
#include <gismo.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSurfMesh.h>
#include <string>

using namespace gismo;

void saveMatrixToFile(const gsMatrix<real_t>& matrix, const std::string& filepath)
{
    std::ofstream file(filepath);
    file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<xml>\n"
         << "  <Matrix rows=\"" << matrix.rows() << "\" cols=\""
         << matrix.cols() << "\" format=\"ascii\">\n";
    for (index_t r = 0; r < matrix.rows(); ++r)
    {
        file << "    ";
        for (index_t c = 0; c < matrix.cols(); ++c)
            file << matrix(r, c) << (c + 1 < matrix.cols() ? " " : "");
        file << "\n";
    }
    file << "  </Matrix>\n</xml>\n";
}

int main(int argc, char** argv)
{

    gsSurfMesh mesh = gsSurfMesh();
    auto subdiv = gsFreeformSubdivision<5, 3>(&mesh);

    for (size_t valence = 3; valence < 10; ++valence)
    {
        if (valence == 4)
            continue;

        gsInfo << "=================\n    Valence " << valence
               << "\n=================\n\n";

        gsMatrix<real_t> coeffs(2 * valence + 1, 2 * valence + 1);

        for (size_t function = 1; function <= 2 * valence + 1; ++function)
        {
            gsInfo << "Function " << function << "\n";
            subdiv.initialize_data_xml(
                "freeformSubdivision/fitting_functions/Val" +
                std::to_string(valence) + "Fct" + std::to_string(function) +
                ".xml");

            subdiv.subdivide();
            auto res = subdiv.smooth(1);
            if(function==1){
            saveMatrixToFile(res[1], "OuterCoefficientsVal" + std::to_string(valence) + ".xml");
            }

            coeffs.row(function - 1) = res[0].transpose().row(2);
            gsInfo << "\n";
        }

        // save all coefficients to file
        saveMatrixToFile(coeffs, "CoefficientsVal" + std::to_string(valence) + ".xml");
    }

    return 0;
}
