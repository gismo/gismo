/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsSurfMesh.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gismo.h>

using namespace gismo;

int main(int argc, char** argv)
{
    gsSurfMesh mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(std::string("off/octtorus.off"), mesh);

    auto subdiv = gsFreeformSubdivision<9>();

    subdiv.initialize_data(mesh);
    subdiv.subdivide(mesh);
    subdiv.make_c1(mesh);

    gsWriteParaview(subdiv.multipatch(mesh), "results/beziers");

    mesh.write("results/mesh_out.off");
    gsWriteParaview(mesh, "results/mesh_out", { });

}
