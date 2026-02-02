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

    auto _readFile = gsReadFile<>(std::string("off/polycube.off"), mesh);

    auto subdiv = gsFreeformSubdivision<5>();

    subdiv.initialize_data(mesh);
    gsWriteParaview(subdiv.multipatch(mesh), "results/initial_data");
    subdiv.make_c1(mesh);
    gsWriteParaview(subdiv.multipatch(mesh), "results/c1");
    subdiv.subdivide(mesh);
    gsWriteParaview(subdiv.multipatch(mesh), "results/subdiv");
    subdiv.subdivide(mesh);
    gsWriteParaview(subdiv.multipatch(mesh), "results/subdiv2");

    // gsWriteParaview(subdiv.multipatch(mesh), "results/beziers");

    // mesh.write("results/mesh_out.off");
    // gsWriteParaview(mesh, "results/mesh_out", { });

}

void minimal_bug_example(){
    gsSurfMesh mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(std::string("off/octtorus.off"), mesh);

    for(gsSurfMesh::Face f : mesh.faces()){

        bool has_ev1 =
            std::any_of(mesh.vertices(f).begin(), mesh.vertices(f).end(),
                        [](gsSurfMesh::Vertex v) {return v.idx() == 0;
                        });

        bool has_ev2(false);
        for (gsSurfMesh::Vertex v : mesh.vertices(f))
        {
            has_ev2 = has_ev2 || v.idx() == 0;// !is_ordinary(mesh, v);
        }

        if(has_ev1 != has_ev2){
            gsInfo << "Discrepancy!\n";
            gsInfo << "First vertex: " << *mesh.vertices(f).begin() << "\n";
            gsInfo << "Last vertex: " << *mesh.vertices(f).end() << "\n";
        }

    }
}
