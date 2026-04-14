
/** @file freeform_kernel.cpp

    @brief Precomputes the kernel matrices used by the freeform subdivision
   scheme.


    \note The example assumes it is run from the build directory so that
    \c ../filedata is reachable. The patch-path option must point to a
    subdirectory of \c filedata.

    \par Command-line arguments
    - \b -p / \b --patchpath (default: \c freeform/bubble/): path, relative
      to \c filedata/, to the directory that contains the model patch files
      \c Val\<v\>Fct\<f\>.xml.
    - \b -v / \b --valence (default: \c 9): maximum extraordinary-vertex
      valence to process. Kernels are generated for all valences
      \f$v \in \{3, 5, 6, \ldots, \mathrm{valence\_max}\}\f$ (valence 4 is
      regular and is skipped).

    Author(s): L. Mussmaecher
*/

#include <gismo.h>
#include <gsIO/gsFileData.hpp>
#include <string>

using namespace gismo;

int main(int argc, char** argv)
{
    // CMD arguments
    std::string patchpath("freeform/bubble/");
    gsCmdLine cmd("Freeform subdivision");
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

    // Iterate all valences
    for (size_t valence = 3; valence <= size_t(valence_max); ++valence)
    {
        if (valence == 4)
            continue;

        gsInfo << "=================\n    Valence " << valence
               << "\n=================\n\n";

        gsMatrix<real_t> coefs(4 * 25 * valence, 2 * valence + 1);

        // Iterate all functions for that valence
        for (size_t function = 0; function < 2 * valence + 1; ++function)
        {
            gsInfo << "Function " << function << "\n";
            // Load the basis function patch file.
            auto filepath = patchpath + "Val" + std::to_string(valence) +
                            "Fct" + std::to_string(function) + ".xml";
            std::vector<std::unique_ptr<gsTensorBSpline<2, real_t>>> patches =
                gsFileData<real_t>(filepath)
                    .getAll<gsTensorBSpline<2, real_t>>();
            for (index_t p = 0; p < (index_t)patches.size(); ++p)
            {
                coefs.block(25 * p, function, 25, 1) =
                    patches[p]->coefs().col(2);
            }
        }

        // Now `coefs` contains the z control values of all functions.
        auto legal_pl = coefs.fullPivLu();
        legal_pl.setThreshold(1e-8);
        gsInfo << "Function Rank: " << legal_pl.rank() << "\n";
        gsInfo << "Kernel size:   " << (2 * valence + 1 - legal_pl.rank()) << "\n";
        gsMatrix<> K = legal_pl.kernel();
        gsWrite(K, "../filedata/" + patchpath + "Val" +
                       std::to_string(valence) + "Kernel.xml");

        gsInfo << "Written functional constraints for valence " << valence
               << " to `"
               << (gsFileManager::findInDataDir("") + patchpath + "Val" +
                   std::to_string(valence) + "Kernel.xml")
               << "`.\n";
    }

    return 0;
}
