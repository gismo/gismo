/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include "gsCore/gsDebug.h"
#include "gsIO/gsWriteParaview.h"
#include "gsMesh2/gsSurfMesh.h"
#include <array>
#include <gismo.h>

using namespace gismo;

class FreeformFaceData{
    public:
    std::array<gismo::gsVector3d<double>, 9> data;
    FreeformFaceData() : data({
          gismo::gsVector3d<double>(0.0, 0.0, 0.0),
          gismo::gsVector3d<double>(0.0, 0.0, 0.0),
          gismo::gsVector3d<double>(0.0, 0.0, 0.0),
          gismo::gsVector3d<double>(0.0, 0.0, 0.0),
          gismo::gsVector3d<double>(0.0, 0.0, 0.0),
          gismo::gsVector3d<double>(0.0, 0.0, 0.0),
          gismo::gsVector3d<double>(0.0, 0.0, 0.0),
          gismo::gsVector3d<double>(0.0, 1.0, 0.0),
          gismo::gsVector3d<double>(0.0, 0.0, 0.0),
      }) {}
};

int main(int argc, char** argv)
{
    gsSurfMesh mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(std::string("off/cube.off"), mesh);
    mesh.add_face_property(std::string("bezier_points"), FreeformFaceData());

    gsProperty<FreeformFaceData> x = mesh.get_face_property<FreeformFaceData>("bezier_points");


    gsInfo << "Sizes:" << std::endl;
    gsInfo << mesh.faces_size() << std::endl;
    gsInfo << x.vector().size() << std::endl << std::endl;;

    for (auto x : x.vector()) {
        gsInfo << "Face: " << "\n";
        for(auto v : x.data) {
            gsInfo << "V: " << v << "\n\n";
        }
    }


    mesh.write("mesh_out.off");
    gsWriteParaview(mesh, "mesh_out", { });

}
