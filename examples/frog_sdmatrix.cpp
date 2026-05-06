/** @file frog_sdmatrix.cpp

    @brief Precomputes the inner (and outer) EV coefficient matrices for
    the frog subdivision scheme.

    For each valence \f$v\f$ from 3 up to \c valence_max (skipping \f$v=4\f$),
    the example determines the number of blending functions from the collated
    file \c Val\<v\>Fcts.xml, then for each function loads the corresponding
    per-function patch file \c Val\<v\>Fct\<f\>.xml, subdivides once, applies
    one round of \f$C^1\f$ smoothing, and reads out the row of
    extraordinary-vertex (EV) Bézier coefficients from the resulting patch.
    These rows are assembled into a coefficient matrix with one row per
    blending function.

    The number of blending functions is determined at runtime from the number
    of \c gsMultiPatch entries in \c Val\<v\>Fcts.xml and is not assumed to
    equal \f$2v+1\f$.

    The full coefficient matrix is saved as
    \c SubdivisionFullVal\<v\>.xml in the current working directory.
    These files are consumed by \c fit_ev to reconstruct the correct \f$C^1\f$
    coefficient when fitting around extraordinary vertices.

    \note The example assumes it is run from the build directory so that
    \c ../filedata is reachable. The patch-path option must point to a
    subdirectory of \c filedata.

    \par Command-line arguments
    - \b -p / \b --frogdir (default: \c frog/bubble/): path, relative
      to \c filedata/, to the directory containing a set of frog spline
      generating functions for each required valence.
    - \b -v / \b --valence (default: \c 9): maximum extraordinary-vertex
      valence to process. Coefficient matrices are generated for all valences
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
    std::string frogdir("frog/bubble/");
    index_t valence_max(6);
    gsCmdLine cmd("frog subdivision");
    cmd.addString("p", "frogdir",
                  "The path to the folder containing a set of frog spline "
                  "generating functions for each required valence.",
                  frogdir);
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
    subdiv.options().setString("frog_dir", frogdir);

    // Iterate all valences
    for (size_t valence = 3; valence <= size_t(valence_max); ++valence)
    {
        if (valence == 4)
            continue;

        gsInfo << "=================\n    Valence " << valence
               << "\n=================\n\n";

        gsFileData<real_t> fd(frogdir + "Val" + std::to_string(valence) +
                              "Fcts.xml");
        auto fcts = fd.getAll<gsMultiPatch<real_t>>();

        size_t point_count(fcts.size() + 20 * valence);
        gsMatrix<real_t> coeffs(fcts.size(), point_count);

        // Iterate all functions for that valence
        for (size_t function = 0; function < fcts.size(); ++function)
        {
            gsInfo << "Function " << function << "\n";
            // Load the basis function patch file.
            subdiv.initialize_data_multipatch(*fcts[function]);

            // Subdivide once
            subdiv.subdivide();
            // Now smooth to get the desired matrices
            gsMatrix<real_t> coefs = subdiv.smooth(1);
            // Collect all coefficients in a matrix row
            coeffs.row(function) = coefs.topRows(point_count).transpose().row(2);
            gsInfo << "\n";
        }

        // save all collected coefficients to file
        gsWrite(coeffs,
                "SubdivisionFullVal" + std::to_string(valence) + ".xml");

        gsInfo << "Written functional constraints for valence " << valence
               << " to `"
               << "SubdivisionFullVal" + std::to_string(valence) + ".xml"
               << "`.\n";
    }

    return 0;
}
