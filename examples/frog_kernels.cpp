
/** @file frog_kernels.cpp

    @brief Precomputes kernel bases for the extraordinary-vertex blending
    functions of the frog subdivision scheme.

    For each valence \f$v\f$ from 3 up to \c valence_max (skipping \f$v=4\f$),
    the example loads all extraordinary-vertex blending functions from the
    collated file \c Val\<v\>Fcts.xml (one \c gsMultiPatch per function),
    stacks the third control-net coordinate of every patch of every function
    into one common coefficient matrix, and computes the kernel of that matrix.
    The resulting kernel basis spans all coefficient vectors that do not change
    the represented EV function and is written to \c Val\<v\>Kernel.xml in the
    same patch directory.

    The number of blending functions is determined at runtime from the number
    of \c gsMultiPatch entries in \c Val\<v\>Fcts.xml and is not assumed to
    equal \f$2v+1\f$.

    These kernel files are consumed by \c gsFrogSplines::smooth() and
    by \c frog_functionals.cpp.

    \note The example assumes it is run from the build directory so that
    \c ../filedata is reachable. The patch-path option must point to a
    subdirectory of \c filedata.

    \par Command-line arguments
    - \b -p / \b --patchpath (default: \c frog/bubble/): path, relative
      to \c filedata/, to the directory that contains the collated model patch
      file \c Val\<v\>Fcts.xml.
    - \b -v / \b --valence (default: \c 9): maximum extraordinary-vertex
      valence to process. Kernel bases are generated for all valences
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

    const size_t D = 4;

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

        gsMatrix<real_t> coefs(4 * (D + 1) * (D + 1) * valence,
                               fcts.size()); 

        // Iterate all functions for that valence
        for (size_t function = 0; function < fcts.size();
             ++function)
        {
            gsInfo << "Function " << function << "\n";

            auto patches = fcts[function]->patches();
            for (index_t p = 0; p < (index_t)patches.size(); ++p)
            {
                coefs.block(25 * p, function, 25, 1) =
                    patches[p]->coefs().col(2);
            }
        }

        // Now `coefs` contains the z control values of all functions.
        auto legal_pl = coefs.fullPivLu();
        legal_pl.setThreshold(1e-4);
        gsInfo << "Function Rank: " << legal_pl.rank() << "\n";
        gsInfo << "Kernel size:   " << (fcts.size() - legal_pl.rank())
               << "\n";
        gsMatrix<> K = legal_pl.kernel();
        gsWrite(K, "../filedata/" + patchpath + "Val" +
                       std::to_string(valence) + "Kernel.xml");

        gsInfo << "Written kernel basis for valence " << valence << " to `"
               << (gsFileManager::findInDataDir("") + patchpath + "Val" +
                   std::to_string(valence) + "Kernel.xml")
               << "`.\n";
    }

    return 0;
}
