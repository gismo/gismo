/** @file volMesh_example.cpp

    @brief Builds a block of hexahedra with gsVolMesh, walks its half-face
    structure and writes out the volume mesh and its boundary surface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, H.M. Verhelst
*/

#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
    index_t n = 2;
    std::string output = "volmesh";

    gsCmdLine cmd("Half-face volume mesh (gsVolMesh) example.");
    cmd.addInt("n", "cells", "Number of cells per direction", n);
    cmd.addString("o", "output", "Base name of the output files", output);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    typedef gsVolMesh<real_t>          Mesh;
    typedef Mesh::Point                Point;
    typedef Mesh::Vertex               Vertex;

    // ---------------------------------------------------------------- build
    //
    // add_hex() takes its eight corners in VTK / CGNS order.  Faces and edges
    // shared with a cell that is already there are found and glued
    // automatically, so the block below ends up with one geometric face per
    // interior interface, not two.
    Mesh mesh;

    const index_t s = n+1;
    std::vector<Vertex> id(s*s*s);
    for (index_t i = 0; i != s; ++i)
        for (index_t j = 0; j != s; ++j)
            for (index_t k = 0; k != s; ++k)
                id[(i*s+j)*s+k] = mesh.add_vertex(Point(i,j,k));

#   define ID(a,b,c) id[(((a)*s+(b))*s+(c))]
    for (index_t i = 0; i != n; ++i)
        for (index_t j = 0; j != n; ++j)
            for (index_t k = 0; k != n; ++k)
                mesh.add_hex(ID(i  ,j  ,k  ), ID(i+1,j  ,k  ),
                             ID(i+1,j+1,k  ), ID(i  ,j+1,k  ),
                             ID(i  ,j  ,k+1), ID(i+1,j  ,k+1),
                             ID(i+1,j+1,k+1), ID(i  ,j+1,k+1));
#   undef ID

    gsInfo << mesh;
    mesh.mesh_statistics();

    std::string msg;
    gsInfo << "\nTopology is "
           << (mesh.is_valid_topology(&msg) ? "sound." : "BROKEN: " + msg) << "\n";
    gsInfo << "Total volume: " << mesh.volume() << "\n";

    // ------------------------------------------------------- navigate a bit
    //
    // Everything volumetric comes out of beta3: mate() takes a dart across a
    // face into the neighbouring cell, and composing it with the inherited
    // beta2 rotates around a geometric edge, one step per incident cell.
    if (n > 1)
    {
        const Vertex centre = id[((n/2*s + n/2)*s) + n/2];
        gsInfo << "\nAround the interior vertex " << centre << ":\n"
               << "  cells    : " << mesh.valence(centre) << "\n"
               << "  edges    : " << mesh.edges(centre).size() << "\n"
               << "  faces    : " << mesh.faces(centre).size() << "\n"
               << "  boundary : " << (mesh.is_boundary(centre) ? "yes" : "no") << "\n";

        const Mesh::Edge e = mesh.find_edge(centre, id[(((n/2+1)*s + n/2)*s) + n/2]);
        if (e.is_valid())
        {
            gsInfo << "\nRadial cycle of the edge " << e << " (one dart per cell):\n  ";
            for (auto c : mesh.cells(e)) gsInfo << c << " ";
            gsInfo << "\n  valence " << mesh.valence(e)
                   << ", length " << mesh.edge_length(e)
                   << ", on the boundary: " << (mesh.is_boundary(e) ? "yes" : "no") << "\n";
        }
    }

    // the first cell, from the inside out
    const Mesh::Cell c0 = *mesh.cells().begin();
    gsInfo << "\nCell " << c0 << " has " << mesh.valence(c0) << " faces, "
           << mesh.n_edges(c0) << " edges, " << mesh.n_vertices(c0)
           << " vertices and volume " << mesh.volume(c0) << ".\n";
    for (auto hf : mesh.halffaces(c0))
    {
        gsInfo << "  half-face of face " << mesh.face(hf) << ":";
        for (auto v : mesh.vertices(hf)) gsInfo << " " << v;
        gsInfo << (mesh.is_boundary(hf) ? "   [boundary]"
                                        : "   [glued to cell " + util::to_string(mesh.opposite_cell(hf).idx()) + "]")
               << "\n";
    }

    // ---------------------------------------------------------------- write
    const std::string vtu = output + ".vtu";
    const std::string off = output + "_boundary.off";

    if (mesh.write(vtu))
        gsInfo << "\nWrote the volume mesh to " << vtu << "\n";

    gsSurfMesh<real_t> bnd = mesh.boundary_mesh();
    if (bnd.write(off))
        gsInfo << "Wrote its boundary (" << bnd.n_vertices() << " vertices, "
               << bnd.n_faces() << " faces) to " << off << "\n";

    return EXIT_SUCCESS;
}
