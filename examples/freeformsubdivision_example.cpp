/** @file subdivisionSurfaces_example.cpp

    @brief Tests different subdivision schemes

    Author(s): A. Mantzaflaris, M.Marsala, D.Tolis, L. Mussmaecher
*/

#include "gsIO/gsWriteParaview.h"
#include "gsMesh2/gsSurfMesh.h"
#include <gismo.h>

using namespace gismo;

int main(int argc, char** argv)
{
    gsSurfMesh mesh = gsSurfMesh();

    gismo::gsVector3d<double> vecs[9] = {
      gismo::gsVector3d<double>(0.0, 0.0, 0.0),
      gismo::gsVector3d<double>(0.0, 0.0, 0.0),
      gismo::gsVector3d<double>(0.0, 0.0, 0.0),
      gismo::gsVector3d<double>(0.0, 0.0, 0.0),
      gismo::gsVector3d<double>(0.0, 0.0, 0.0),
      gismo::gsVector3d<double>(0.0, 0.0, 0.0),
      gismo::gsVector3d<double>(0.0, 0.0, 0.0),
      gismo::gsVector3d<double>(0.0, 1.0, 0.0),
      gismo::gsVector3d<double>(0.0, 0.0, 0.0),
    };

    mesh.add_face_property(std::string("bezier_points"), vecs);
    auto _readFile = gsReadFile<>(std::string("off/cube.off"), mesh);

    mesh.write("mesh_out.off");
    gsWriteParaview(mesh, "mesh_out", { });

}
