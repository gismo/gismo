/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include <gismo.h>
#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSurfMesh.h>

using namespace gismo;

int main(int argc, char** argv)
{
    // Command line
    std::string filedata("off/polycube.off");
    bool no_smooth(false);
    size_t steps(1);

    // Inputs
    gsCmdLine cmd("Freeform subdivision");
    cmd.addPlainString("filename", "File containing mesh.", filedata);
    cmd.addSwitch("no_smooth", "C1 smoothing before subdivision.", no_smooth);
    cmd.addInt("steps", "Number of subdivision steps.", steps);
    try
    {
        cmd.getValues(argc, argv);
    }
    catch (int errorcode)
    {
        return errorcode;
    }

    gsSurfMesh mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(filedata, mesh);

    auto subdiv = gsFreeformSubdivision<5, 3>();

    subdiv.initialize_data(mesh);
    gsWriteParaview(subdiv.multipatch(mesh), "results/initial_data");

    if (!no_smooth)
    {
        subdiv.make_c1(mesh);
        gsWriteParaview(subdiv.multipatch(mesh), "results/c1");
    }
    for (size_t i = 0; i < steps; ++i)
    {
        subdiv.subdivide(mesh);
        gsWriteParaview(subdiv.multipatch(mesh),
                        "results/subdiv" + std::to_string(i));
    }
}

void minimal_bug_example()
{
    gsSurfMesh mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(std::string("off/octtorus.off"), mesh);

    for (gsSurfMesh::Face f : mesh.faces())
    {

        bool has_ev1 =
            std::any_of(mesh.vertices(f).begin(), mesh.vertices(f).end(),
                        [](gsSurfMesh::Vertex v) { return v.idx() == 0; });

        bool has_ev2(false);
        for (gsSurfMesh::Vertex v : mesh.vertices(f))
        {
            has_ev2 = has_ev2 || v.idx() == 0; // !is_ordinary(mesh, v);
        }

        if (has_ev1 != has_ev2)
        {
            gsInfo << "Discrepancy!\n";
            gsInfo << "First vertex: " << *mesh.vertices(f).begin() << "\n";
            gsInfo << "Last vertex: " << *mesh.vertices(f).end() << "\n";
        }
    }
}
