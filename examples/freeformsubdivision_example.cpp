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
#include <gismo.h>

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
        gsInfo << "Size: " << points.size() << "\n";
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

    gismo::gsVector3d<real_t> corresponding_point(const gsSurfMesh& mesh, gsSurfMesh::Face face, gsSurfMesh::Vertex vertex)
    {
        gismo::gsVector3d<real_t> result;

        size_t vertex_position(0);

        for(auto v : mesh.vertices(face)){
            if(v == vertex) break;
            ++vertex_position;
        }

        switch(vertex_position)
        {
            case 0:
                result = data[0];
                break;
            case 1:
                result = data[2];
                break;
            case 2:
                result = data[8];
                break;
            case 3:
                result = data[6];
                break;
            default:
                gsInfo << "Oh no\n";
        }        

        return result;
    }

    std::array<gismo::gsVector3d<real_t>, 25> full_patch(const gsSurfMesh& mesh, gsSurfMesh::Face face)
    {
        // Create an array to hold the 5x5 control points for the full patch
        std::array<gismo::gsVector3d<real_t>, 25> result;
        result.fill(gismo::gsVector3d<real_t>(0., 0., 0.));

        // Prepare data of other patches
        auto patch_data = mesh.get_face_property<FreeformFaceData>("bezier_points");


        // === INNER ===
        // These are directly saved in the data structure.
        for(size_t i = 0; i<9;++i){
            size_t ih((i%3)+1);
            size_t iv((i/3)+1);
            result[5 * iv + ih] = this->data[i];
        }

   
        // === EDGES ===

        // Prepare indices.
        std::vector<std::array<size_t, 3>> hedge_result = {
            {15,10,5},
            {1,2,3},
            {9,14,19},
            {23,22,21},
        };

        // Iterator over all half edges adjacent to this one.
        size_t hedge_counter(0);
        for(auto hedge : mesh.halfedges(face)){
            // Get the two half-edges and then faces for this half-edge.
            auto hedge1 (hedge);
            auto hedge2(mesh.opposite_halfedge(hedge));
            auto face1 (mesh.face(hedge1)); // this should be the current face.
            auto face2 (mesh.face(hedge2));

            // Get the two 3-sets of corresponding points.
            auto points1 = patch_data.vector()[face1.idx()].corresponding_points(mesh, face1, hedge1);
            auto points2 = patch_data.vector()[face2.idx()].corresponding_points(mesh, face2,hedge2);
            
            // Calculate the average of each control point with its partner on the other side and store it in the appropriate control point in the result.
            for(int i = 0; i<3; ++i){
                result[hedge_result[hedge_counter][i]] = points1[i] * 0.5 + points2[2-i] * 0.5;
            }

            ++hedge_counter;
        }

        // === CORNERS ===

        // Prepare indices.
        std::vector<size_t> vertex_result = {
            0, 4, 24, 20
        };

        // Iterate over all faces adjacent to this vertex.
        size_t vertex_counter(0);
        for(auto v : mesh.vertices(face)){
            // Sum over the inner control point neares to this corner for each face.
            gsVector3d<real_t> sum(gismo::gsVector3d<real_t>(0., 0., 0.));
            real_t count(0.0);
            for(auto extra_face : mesh.faces(v)) {
                sum += patch_data.vector()[extra_face.idx()].corresponding_point(mesh, extra_face, v);
                count += 1.0;
            }
            // Store this as the corner control point of the result.
            result[vertex_result[vertex_counter]] = sum / count;
            ++vertex_counter;
        }

        return result;
    }
};

int main(int argc, char** argv)
{
    gsSurfMesh mesh = gsSurfMesh();

    auto _readFile = gsReadFile<>(std::string("off/bitorus.off"), mesh);
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
