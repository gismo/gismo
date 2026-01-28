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
    mesh.add_face_property(std::string("bezier_points"), gsFreeformFaceData());
    gsProperty<gsFreeformFaceData> patch_data = mesh.get_face_property<gsFreeformFaceData>("bezier_points");
    for (auto f : mesh.faces()){
        patch_data.vector()[f.idx()] = gsFreeformFaceData(mesh, f);
    }

    auto subdiv = gsFreeformSubdivision();
    subdiv.make_c1(mesh);

    gsMultiPatch<> patch;

    for(auto face : mesh.faces()) {
        patch.addPatch(patch_data.vector()[face.idx()].patch());
    }

    gsWriteParaview(patch, "results/beziers");

    mesh.write("results/mesh_out.off");
    gsWriteParaview(mesh, "results/mesh_out", { });

}
