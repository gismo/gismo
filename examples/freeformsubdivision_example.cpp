/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include "gsCore/gsConfig.h"
#include "gsCore/gsDebug.h"
#include "gsCore/gsMultiPatch.h"
#include "gsIO/gsWriteParaview.h"
#include "gsMesh2/gsSurfMesh.h"
#include "gsNurbs/gsTensorBSpline.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <gismo.h>
#include <iterator>

using namespace gismo;

class FreeformFaceData{
    public:
    std::array<gismo::gsVector3d<real_t>, 9> data;
    FreeformFaceData() : data() {
        data.fill(gismo::gsVector3d<real_t>(0., 0., 0.));
    }

    FreeformFaceData(const gsSurfMesh& mesh, gsSurfMesh::Face face) : data() {
        std::vector<gismo::gsVector3d<real_t>> points;

        for (auto v : mesh.vertices(face)){
            points.emplace_back(mesh.position(v));
        }
        assert(points.size() == 4);

        for(size_t i = 0; i<9; ++i){
            real_t ih(i % 3 + 1);
            real_t iv(std::floor(i / 3) + 1);
          
            data[i] =
                points[0] * (4-ih) * (4-iv)/16. +
                points[1] * ih * (4-iv)/16. +
                points[3] * (4-ih) * iv/16. +
                points[2] * ih * iv/16.;
        }
        
    }
    std::array<gismo::gsVector3d<real_t>, 3> corresponding_points(const gsSurfMesh& mesh, gsSurfMesh::Face face, gsSurfMesh::Halfedge hedge)
    {
        std::array<gismo::gsVector3d<real_t>, 3> result;

        size_t hedge_position(0);

        for(auto he : mesh.halfedges(face)){
            if(he == hedge) break;
            ++hedge_position;
        }

        switch(hedge_position)
        {
            case 1:
                result[0] = data[0];
                result[1] = data[1];
                result[2] = data[2];
                break;
            case 2:
                result[0] = data[2];
                result[1] = data[5];
                result[2] = data[8];
                break;
            case 3:
                result[0] = data[8];
                result[1] = data[7];
                result[2] = data[6];
                break;
            case 0:
                result[0] = data[6];
                result[1] = data[3];
                result[2] = data[0];
                break;
            default:
                gsInfo << "Oh no\n";
        }        

        return result;
    }

    std::array<gismo::gsVector3d<real_t>, 25> full_patch(const gsSurfMesh& mesh, gsSurfMesh::Face face)
    {
        std::array<gismo::gsVector3d<real_t>, 25> result;
        result.fill(gismo::gsVector3d<real_t>(0., 0., 0.));

        // TEMPORARY
        std::vector<gismo::gsVector3d<real_t>> points;

        for (auto v : mesh.vertices(face)){
            points.emplace_back(mesh.position(v));
        }
        assert(points.size() == 4);

        for(size_t i = 0; i<25; ++i){
            real_t ih(i % 5);
            real_t iv(std::floor(i / 5));
          
            result[i] =
                points[0] * (4-ih) * (4-iv)/16. +
                points[1] * ih * (4-iv)/16. +
                points[3] * (4-ih) * iv/16. +
                points[2] * ih * iv/16.;
        }
        // TEMPORARY END

        // inner
        for(size_t i = 0; i<9;++i){
            size_t ih((i%3)+1);
            size_t iv((i/3)+1);
            result[5 * iv + ih] = this->data[i];
        }

       
        // get halfedges
        std::vector<gsSurfMesh::Halfedge> hedges;
        for(auto he : mesh.halfedges(face)){
            hedges.push_back(he);
        }
        // get other data
        auto patch_data = mesh.get_face_property<FreeformFaceData>("bezier_points");

        // 'upper'
{        auto hedge = hedges[1];
        auto opp = mesh.opposite_halfedge(hedge);
        auto opp_face = mesh.face(opp);
        auto other_points = patch_data.vector()[opp_face.idx()].corresponding_points(mesh, opp_face, opp);

        result[1] = other_points[2] * 0.5 + this->data[0] * 0.5;
        result[2] = other_points[1] * 0.5 + this->data[1] * 0.5;
        result[3] = other_points[0] * 0.5 + this->data[2] * 0.5;
}
{        auto hedge = hedges[3];
        auto opp = mesh.opposite_halfedge(hedge);
        auto opp_face = mesh.face(opp);
        auto other_points = patch_data.vector()[opp_face.idx()].corresponding_points(mesh, opp_face, opp);

        result[23] = other_points[0] * 0.5 + this->data[8] * 0.5;
        result[22] = other_points[1] * 0.5 + this->data[7] * 0.5;
        result[21] = other_points[2] * 0.5 + this->data[6] * 0.5;
}
        return result;
    }
};

int main(int argc, char** argv)
{
    gsSurfMesh mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(std::string("off/octtorus.off"), mesh);
    mesh.add_face_property(std::string("bezier_points"), FreeformFaceData());
    gsProperty<FreeformFaceData> patch_data = mesh.get_face_property<FreeformFaceData>("bezier_points");
    for (auto f : mesh.faces()){
        patch_data.vector()[f.idx()] = FreeformFaceData(mesh, f);
    }

    // Now convert into bezier stuff
    // 1. construction of a knot vector for each direction
    gsKnotVector<> kv1(0, 1, 0, 5);
    gsKnotVector<> kv2(0, 1, 0, 5);

 
    // 2. construction of a basis
    gsTensorBSplineBasis<2, real_t> basis(kv1, kv2);

    gsMultiPatch<> patch;

    for(auto face : mesh.faces()) {

        auto ffdata(patch_data.vector()[face.idx()]);
        auto full_patch = ffdata.full_patch(mesh, face);

        gsMatrix<> coeffs(25,3);

        for(size_t i = 0; i < full_patch.size(); ++i){
            coeffs(i, 0) = full_patch[i].x();
            coeffs(i, 1) = full_patch[i].y();
            coeffs(i, 2) = full_patch[i].z();
        }

        patch.addPatch(gsTensorBSpline<2, real_t>(basis, coeffs));
    }

    gsWriteParaview(patch, "results/beziers");

    mesh.write("results/mesh_out.off");
    gsWriteParaview(mesh, "results/mesh_out", { });

}
