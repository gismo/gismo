/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include <gsCore/gsMultiPatch.h>
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

    gsMultiPatch<> patch;

    gsProperty<gsFreeformFaceData<9>> patch_data = mesh.get_face_property<gsFreeformFaceData<9>>("bezier_points");
    for(auto face : mesh.faces()) {
        patch.addPatch(patch_data.vector()[face.idx()].patch());
    }

    gsWriteParaview(patch, "results/beziers");


    mesh.write("results/mesh_out.off");
    gsWriteParaview(mesh, "results/mesh_out", { });

}
