/** @file quadsplit_surf

    @brief quadsplit surface

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#include <gismo.h>

using namespace gismo;



int main(int argc, char** argv)
{
    index_t r(0), s(0);
    std::string fn("off/octtorus.off");
    std::string fn_patch("");
    bool plot = false;

    //! [Parse Command line]
    gsCmdLine cmd("Hi, give me a CC mesh");
    cmd.addInt("c", "ref", "Number of refinement steps", r);
    cmd.addInt("s", "sam", "Number of sampling refinement steps", s);
    cmd.addPlainString("filename", "File containing mesh", fn);
    cmd.addString("g", "geometry", "File containing multipatch geometry to compare)", fn_patch);
    cmd.addSwitch("plot", "Plot the results", plot);
    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }

    gsSurfMesh mesh;
    //mesh.read(fn);
    gsReadFile<>(fn, mesh);

    gsInfo << "Input: " << mesh.n_vertices() << " vertices, "
        << mesh.n_edges() << " edges, " << mesh.n_faces() << " faces. \n";

    //subdivision before creating ACC3 patches
    for (index_t i = 0; i < r; ++i)
        mesh.quad_split();

    mesh.write("qsplit" + util::to_string(r) + ".off");

    return EXIT_SUCCESS;
}


//=============================================================================
