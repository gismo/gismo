/** @file freeform_sdmatrix.cpp

    @brief Precomputes the inner (and outer) EV coefficient matrices for
    the freeform subdivision scheme.

    For each valence \f$v\f$ from 3 up to \c valence_max (skipping \f$v=4\f$),
    the example loads all \f$2v+1\f$ basis-function patches from the model patch
    directory (\c Val\<v\>Fct\<f\>.xml), subdivides each once, applies one
    round of \f$C^1\f$ smoothing, and reads out the row of extraordinary-vertex
    (EV) Bézier coefficients from the resulting patch. These rows are assembled
    into an inner coefficient matrix (all \f$2v+1\f$ rows) and an outer
    coefficient vector (coefficients of basis function 1 only).

    The inner matrix is saved as \c CoefficientsInnerVal\<v\>.xml and the outer
    vector as \c CoefficientsOuterVal\<v\>.xml inside the model patch directory.
    These files are consumed by \c fit_ev to reconstruct the correct \f$C^1\f$
    coefficient when fitting around extraordinary vertices.

    \note The example assumes it is run from the build directory so that
    \c ../filedata is reachable. The patch-path option must point to a
    subdirectory of \c filedata.

    \par Command-line arguments
    - \b -p / \b --patches (default: \c freeform/bubble/): path, relative
      to \c filedata/, to the directory containing the model patch files
      \c Val\<v\>Fct\<f\>.xml.
    - \b -v / \b --valence (default: \c 9): maximum extraordinary-vertex
      valence to process. Coefficient matrices are generated for all valences
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
    std::string patches("freeform/bubble/");
    index_t valence_max(6);
    gsCmdLine cmd("Freeform subdivision");
    cmd.addString("p", "patches",
                  "The path to the files containing the model patches for EV "
                  "subdivision.",
                  patches);
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
    auto subdiv = gsFreeformSubdivision<5>(&mesh, 3);
    subdiv.options().setString("model_patch_path", patches);

    // Iterate all valences
    for (size_t valence = 3; valence <= size_t(valence_max); ++valence)
    {
        if (valence == 4)
            continue;

        gsInfo << "=================\n    Valence " << valence
               << "\n=================\n\n";

        size_t point_count(2 * valence + 1 + 20 * valence);
        gsMatrix<real_t> coeffs(2 * valence + 1, point_count);

        // Iterate all functions for that valence
        for (size_t function = 0; function < 2 * valence + 1; ++function)
        {
            gsInfo << "Function " << function << "\n";
            // Load the basis function patch file.
            subdiv.initialize_data_xml(patches + "Val" +
                                       std::to_string(valence) + "Fct" +
                                       std::to_string(function) + ".xml");

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
