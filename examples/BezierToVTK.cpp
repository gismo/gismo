#include <gismo.h>
#include <fstream>
#include <iostream>
#include <gsLsdyna/gsBextFormat.h>
#include <gsIO/gsParaviewUtils.h>

using namespace gismo;

#define VTK_BEZIER_QUADRILATERAL 77

// /// @brief ID transformation between G+Smo and vtk  control point IDs 
// /// @param nU Number of control points in u parametric direction
// /// @param nV Number of control points in u parametric direction
// /// @return 
// gsSparseMatrix<> vtkIDTransform(index_t nU, index_t nV)
// {
//     // T converts coefs from G+Smo's convetnion to ParaView's convention
//     gsSparseMatrix<> T(nU*nV, nU*nV);

//     // T( Paraview , gismo  )
//     // Corners ( always 0-3 )
//     T(0,0) = 1;
//     T(1,nU-1) = 1;
//     T(2,nU*nV-1) = 1;
//     T(3,nU*nV-nU) = 1;

//     // Edges
//     for (index_t i=1;i<nU-1;++i) // Parallel to u
//     {
//         T(3+i,i) = 1;
//         T(3+ (nU-2) + (nV-2) + i, nU*(nV-1) + i  ) = 1;


//     }
//     for (index_t j=1;j<nV-1;++j) // Parallel to v
//     {
//         T( 3 + (nU-2) + j, (j+1)*nU-1) = 1 ;
//         T( 3 + 2*(nU-2) + (nV-2) + j, nU*j) = 1;
//     }
//     // Internal
//     for (index_t i=0;i<nU-2;++i)
//     {
//         for (index_t j=0;j<nV-2;++j)
//         {
//             T(2*( nU + nV) - 4 + j*(nU-2) + i , nU*(j+1)+(i+1)) = 1;
//         }
//     }
//     return T;
// }



int main(int argc, char* argv[])
{
    // gsTensorBSpline<2> geom(*gsNurbsCreator<>::BSplineSquare(10, 0, 0).get());
    // gsTensorBSpline<2> geom(*gsNurbsCreator<>::BSplineFatQuarterAnnulus().get());
    gsFileData<> inFile("surfaces/thbs_face_3levels.xml");
    gsTHBSpline<2> geom;
    inFile.getId(0, geom);

    gsMultiPatch<> mPatch;
    mPatch.addPatch(geom);
    mPatch.addPatch(geom.clone());
    mPatch.patch(1).scale(3,1);
    // geom.degreeElevate(2);
    // geom.uniformRefine(4);

    // std::map<index_t, ElementBlock> ElBlocks = BezierOperator(geom.basis() );
    std::ofstream out;
    out.open("./bezier.vtu");

    // out << "<?xml version=\"1.0\"?>\n"
    //     << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n"
    //     << "\t<UnstructuredGrid>\n";


    // for (auto const& pair : ElBlocks)
    // {
    //     gsInfo << ΕlBlock2vtk(pair.second, geom.coefs());
    // }
    // out << "\t</UnstructuredGrid>\n";
    // out << "</VTKFile>\n";
    out << BezierVTK(mPatch);
    // gsInfo << toVTK(mPatch)[0] << "\n";)
    
    out.close();

}