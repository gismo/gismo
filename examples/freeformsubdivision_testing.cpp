/** @file freeformsubdivision_testing.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include <gismo.h>
#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSurfMesh.h>

using namespace gismo;
gismo::gsTensorBSpline<2, real_t> load_patch(int degree, int valence,
                                             std::string subtype);

int main(int argc, char** argv)
{
    auto coarse = load_patch(5, 5, "coarse");
    auto fine_1 = load_patch(5, 5, "fine_1");
    auto fine_2 = load_patch(5, 5, "fine_2");
    auto fine_3 = load_patch(5, 5, "fine_3");
    auto fine_4 = load_patch(5, 5, "fine_4");

    gsMultiPatch<> patch;
    patch.addPatch(coarse);
    patch.addPatch(fine_1);
    patch.addPatch(fine_2);
    patch.addPatch(fine_3);
    patch.addPatch(fine_4);

    for (real_t s = 0.0; s <= 1.0; s += 0.25)
    {
        for (real_t t = 0.0; t <= 1.0; t += 0.25)
        {
            gsVector<real_t, 2> param = gsVector<real_t, 2>::vec(s,t);
            gsVector<real_t, 2> point = fine_1.eval(param);

            gsVector<real_t> closest_point = gsVector<real_t, 2>::vec(s / 2.,  t / 2.);
            coarse.closestPointTo(point, closest_point, 1e-6, true);
            gsVector<real_t, 2> point2 = coarse.eval(closest_point);

            gsInfo << "Value at   (" << param.transpose() << ") is (" << point.transpose()
                   << ").\n";
            gsInfo << "Closest at (" << closest_point.transpose() << ") is (" << point2.transpose()
                   << ").\n";
        }
    }

    gsWriteParaview(patch, "results/patches");
}

gismo::gsTensorBSpline<2, real_t> load_patch(int degree, int valence,
                                             std::string subtype)
{
    // build the filepath
    auto path = "freeformSubdivision/control_net_d" + std::to_string(degree) +
                "_v" + std::to_string(valence) + "_" + subtype + ".xml";

    gsInfo << "Loading `" << path << "`.\n";

    // load the file containing the correct matrix of controlp points
    gsFileData<real_t> filedata(path);
    auto mat = filedata.getFirst<gsMatrix<real_t>>();

    gsKnotVector<> kv1(0, 1, 0, degree);
    gsKnotVector<> kv2(0, 1, 0, degree);
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);

    return gsTensorBSpline<2>(basis, *mat);
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
