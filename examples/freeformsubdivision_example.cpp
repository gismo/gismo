/** @file freeformsubdivision_example.cpp

    @brief Tests the free form subdivision.

    Author(s): L. Mussmaecher
*/

#include <gismo.h>
#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsFreeformSubdivision.h>
#include <gsMesh2/gsSurfMesh.h>
#include <map>
#include <tuple>
#include <cmath>

using namespace gismo;

/// \brief Converts a vector of tensor BSplines to a gsSurfMesh
///
/// This function takes a collection of TensorBSpline patches and creates a quad mesh
/// where each face corresponds to one patch. The corners of each face are positioned
/// at the corner control points of the corresponding BSpline. Adjacent patches are
/// automatically connected by detecting shared corner control points.
///
/// Each face gets a property "bezier_points" containing a gsFreeformFaceData<N,3> object
/// with the control points from the original BSpline and a back reference to the face.
/// The control points are oriented such that the (0,0) to (0,N-1) first row aligns
/// with the first halfedge of the face.
///
/// \tparam N Degree+1 of the BSpline patches in each direction
/// \param patches Vector of unique_ptr to gsTensorBSpline<2,real_t> objects
/// \return gsSurfMesh with quad faces and bezier_points property
template<size_t N>
gsSurfMesh bsplinesToSurfMesh(const std::vector<std::unique_ptr<gsTensorBSpline<2,real_t>>>& patches)
{
    gsSurfMesh mesh;
    
    // Add face property for storing gsFreeformFaceData with control points
    auto bezier_points = mesh.add_face_property<gsFreeformFaceData<N, 3>>("bezier_points");
    
    // Map from corner positions to vertex indices for detecting shared vertices
    // We use a tolerance-based comparison for floating point coordinates
    std::map<std::tuple<real_t, real_t, real_t>, gsSurfMesh::Vertex> cornerMap;
    const real_t tol = 1e-10;
    
    auto findOrCreateVertex = [&](const gsMatrix<real_t>& point) -> gsSurfMesh::Vertex {
        // Round coordinates for map lookup
        real_t x = std::round(point(0) / tol) * tol;
        real_t y = std::round(point(1) / tol) * tol;
        real_t z = std::round(point(2) / tol) * tol;
        auto key = std::make_tuple(x, y, z);
        
        auto it = cornerMap.find(key);
        if (it != cornerMap.end()) {
            return it->second;
        } else {
            gsSurfMesh::Vertex v = mesh.add_vertex(gsSurfMesh::Point(point(0), point(1), point(2)));
            cornerMap[key] = v;
            return v;
        }
    };
    
    // Process each patch
    for (size_t patchIdx = 0; patchIdx < patches.size(); ++patchIdx)
    {
        const auto& patch = patches[patchIdx];
        
        // Get control points and dimensions
        const gsMatrix<real_t>& coefs = patch->coefs();
        
        // Extract corner control points (lexicographic indexing: i + j*n_u)
        // BSpline corners: (0,0), (N-1,0), (N-1,N-1), (0,N-1) in (u,v) coordinates
        // Map to mesh vertices based on their physical positions
        std::vector<gsSurfMesh::Vertex> corners(4);
        corners[0] = findOrCreateVertex(coefs.row(0 + 0*N));         // BSpline (0,0) → v0
        corners[1] = findOrCreateVertex(coefs.row((N-1) + 0*N));     // BSpline (N-1,0) → v1
        corners[2] = findOrCreateVertex(coefs.row((N-1) + (N-1)*N)); // BSpline (N-1,N-1) → v2
        corners[3] = findOrCreateVertex(coefs.row(0 + (N-1)*N));     // BSpline (0,N-1) → v3
        
        // Add face
        gsSurfMesh::Face f = mesh.add_quad(corners[0], corners[1], corners[2], corners[3]);
        
        // Build control points matrix for gsFreeformFaceData
        // The first halfedge goes from v3→v0 (corners[3]→corners[0])
        // Control point layout should be:
        //   faceControlPoints(0,0) near v0 = BSpline (0,0)
        //   faceControlPoints(0,N-1) near v1 = BSpline (N-1,0)
        //   faceControlPoints(N-1,0) near v3 = BSpline (0,N-1)
        //   faceControlPoints(N-1,N-1) near v2 = BSpline (N-1,N-1)
        // So the mapping is: faceControlPoints(i,j) = BSpline(j, i)
        gsMatrix<gsVector<real_t, 3>, N, N> faceControlPoints;
        
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                // Map face matrix (i,j) to BSpline (u,v) = (j,i)
                index_t linearIdx = j + i * N;
                // Extract the row as a 3D point and assign it
                auto point = coefs.row(linearIdx);
                faceControlPoints(i, j) << point(0), point(1), point(2);
            }
        }
        
        // Create gsFreeformFaceData with control points and face back reference
        bezier_points[f] = gsFreeformFaceData<N, 3>(faceControlPoints, f);
    }
    
    return mesh;
}

int main(int argc, char** argv)
{
    // Command line
    std::string filedata("off/octtorus.off");
    std::string operations("sd");
    bool control_net(false);

    // Inputs
    gsCmdLine cmd("Freeform subdivision");
    cmd.addPlainString("filename", "File containing mesh.", filedata);
    cmd.addString("o", "operations", "Operations to perform on the mesh. Use d for subdivision and s for (c1) smoothing", operations);
    cmd.addSwitch("cnet", "Shows the control net of the patches.", control_net);
    try
    {
        cmd.getValues(argc, argv);
    }
    catch (int errorcode)
    {
        return errorcode;
    }

    gsSurfMesh mesh = gsSurfMesh();
    auto subdiv = gsFreeformSubdivision<5, 3>(&mesh);


    std::string xml(".xml");
    std::string off(".off");
    // Check the filetype to be loaded.
    if(std::equal(filedata.begin() + filedata.size() - xml.size(), filedata.end(), xml.begin())){
        gsInfo << "Loading xml\n";
        // An XML is assumed to be a collection of TensorBSplines that form a mesh.
        mesh = bsplinesToSurfMesh<5>(gsFileData<real_t>(
                            filedata)
                            .getAll<gsTensorBSpline<2, real_t>>());
    } else if (std::equal(filedata.begin() + filedata.size() - off.size(), filedata.end(), off.begin())){
        gsInfo << "Loading off\n";
        // An off is directly loaded
        auto _readFile = gsReadFile<>(filedata, mesh);
        subdiv.initialize_data();
    } else {
        gsWarn << "Unsupported Filetype!\n";
        return 1;
    }


    gsWriteParaview(subdiv.multipatch(), "results/initial_data", 1000, false, control_net);

    size_t i(1);
    for(char c : operations){
        switch (c) {
            case 'd':
                gsInfo << "Step " << std::string(i, 'I') << ": Subdividing.\n";
                subdiv.subdivide();
                break;
            case 's':
                gsInfo << "Step " << std::string(i, 'I') << ": Smoothing.\n";
                subdiv.smooth(1);
                break;
            default:
                break;
       }
        gsWriteParaview(subdiv.multipatch(),
                        "results/step" + std::string(i++, 'I'), 1000, false, control_net);
    }

    return 0;
}
