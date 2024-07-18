#include <gismo.h>
#include <fstream>
#include <iostream>
#include<gsLsdyna/gsBextFormat.h>

using namespace gismo;

#define VTK_BEZIER_QUADRILATERAL 77

/// @brief ID transformation between G+Smo and vtk  control point IDs 
/// @param nU Number of control points in u parametric direction
/// @param nV Number of control points in u parametric direction
/// @return 
gsSparseMatrix<> gs2vtk(index_t nU, index_t nV)
{
    // T converts coefs from G+Smo's convetnion to ParaView's convention
    gsSparseMatrix<> T(nU*nV, nU*nV);

    // T( Paraview , gismo  )
    // Corners ( always 0-3 )
    T(0,0) = 1;
    T(1,nU-1) = 1;
    T(2,nU*nV-1) = 1;
    T(3,nU*nV-nU) = 1;

    // Edges
    for (index_t i=1;i<nU-1;++i) // Parallel to u
    {
        T(3+i,i) = 1;
        T(3+ (nU-2) + (nV-2) + i, nU*(nV-1) + i  ) = 1;


    }
    for (index_t j=1;j<nV-1;++j) // Parallel to v
    {
        T( 3 + (nU-2) + j, (j+1)*nU-1) = 1 ;
        T( 3 + 2*(nU-2) + (nV-2) + j, nU*j) = 1;
    }
    // Internal
    for (index_t i=0;i<nU-2;++i)
    {
        for (index_t j=0;j<nV-2;++j)
        {
            T(2*( nU + nV) - 4 + j*(nU-2) + i , nU*(j+1)+(i+1)) = 1;
        }
    }
    return T;
}


std::string toDataArray(gsMatrix<index_t> mat, std::map<std::string, std::string> attributes={{"",""}})
{
    std::stringstream stream;

    stream.setf(std::ios::fixed);  // write floating point values in fixed-point notation.
    // stream.precision(5);
    // stream << std::setfill('0') << std::setw(16);
    // Format as vtk xml string
    stream << "<DataArray type=\"Float32\" format=\"ascii\" ";
    for (auto const& block : attributes)
    {
        if (block.first!="")
        stream << block.first <<"=\""<< block.second <<"\" ";
    }
    stream <<">\n";
    // For every point
    for (index_t i = 0; i < mat.rows(); ++i)
    {
        for (index_t j = 0; j != mat.cols(); ++j)
            stream << mat(i, j) << " ";
        for (index_t j = mat.cols(); j < 3; ++j) stream << "0 ";
        if (mat.rows()-1==i) break;
        stream << "\n";
    }
    stream << "\n</DataArray>\n";

    return stream.str();
}

std::string toDataArray(gsMatrix<real_t> mat, std::map<std::string, std::string> attributes={{"",""}})
{
    std::stringstream stream;

    stream.setf(std::ios::fixed);  // write floating point values in fixed-point notation.
    // stream.precision(5);
    // stream << std::setfill('0') << std::setw(16);
    // Format as vtk xml string
    stream << "<DataArray type=\"Float32\" format=\"ascii\" ";
    stream << "NumberOfComponents=\""<< /* TODO mat.cols() */ 3 << "\" ";
    for (auto const& block : attributes)
    {
        if (block.first!="")
        stream << block.first <<"=\""<< block.second <<"\" ";
    }
    stream <<">\n";
    // For every point
    for (index_t i = 0; i < mat.rows(); ++i)
    {
        for (index_t j = 0; j != mat.cols(); ++j)
            stream << mat(i, j) << " ";
        for (index_t j = mat.cols(); j < 3; ++j) stream << "0 ";
        if (mat.rows()-1==i) break;
        stream << "\n";
    }
    stream << "\n</DataArray>\n";

    return stream.str();
}

std::string toDataArray(index_t num, std::map<std::string, std::string> attributes={{"",""}})
{
    std::stringstream stream;

    stream.setf(std::ios::fixed);  // write floating point values in fixed-point notation.
    // stream.precision(5);
    // stream << std::setfill('0') << std::setw(16);
    // Format as vtk xml string
    stream << "<DataArray type=\"Float32\" format=\"ascii\" ";
    for (auto const& block : attributes)
    {
        if (block.first!="")
        stream << block.first <<"=\""<< block.second <<"\" ";
    }
    stream <<">\n";
    stream << num ;
    stream << "\n</DataArray>\n";

    return stream.str();
}

//TODO Delete at some point
std::string bez2vtk( gsTensorBSpline<2,real_t> bezier)
{
    GISMO_ENSURE(bezier.parDim() == 2, "Bezier Paraview Output only implemented for bi-variate splines!");
    GISMO_ASSERT( bezier.basis().component(0).isClamped() && 
                  bezier.basis().component(0).size() == bezier.degree(0)+1 &&
                  bezier.basis().component(1).size() == bezier.degree(1)+1,
                  "The patch is not Bezier!");
    std::stringstream stream;
    gsSparseMatrix<> T = gs2vtk(bezier.degree(0)+1, bezier.degree(1)+1);
    gsMatrix<> newCoefs = T * bezier.coefs();
    gsMatrix<index_t> connectivity(1,bezier.coefsSize());

    stream << "<Piece NumberOfPoints=\""<< bezier.coefsSize()<<"\" NumberOfCells=\"1\">\n";
    
    // TODO remove
    stream << "<PointData>\n";
    stream << toDataArray( newCoefs, {{"Name","Position"}});
    stream << "</PointData>\n";

    stream << "<CellData HigherOrderDegrees=\"HigherOrderDegrees\">\n";
    gsMatrix<index_t> degrees(1,3);
    degrees << bezier.degree(0), bezier.degree(1), 0; 
    // TODO: Handle indentation in a neat way
    stream << "" << toDataArray(degrees,{{"Name","HigherOrderDegrees"},{"NumberOfComponents","3"}});
    stream << "</CellData>\n";
    stream << "<Points>\n";
    
    stream << ""<< toDataArray( newCoefs, {{"Name","Points"}});
    stream << "</Points>\n";

    stream << "<Cells>\n";
    connectivity.asVector().setLinSpaced(0,bezier.coefsSize()-1);
    stream << toDataArray(connectivity, {{"Name","connectivity"}});
    stream << toDataArray(bezier.coefsSize(), {{"Name","offsets"}});
    stream << toDataArray(77, {{"Name","types"}});
    stream << "</Cells>\n";
    stream << "</Piece>\n";

    return stream.str();
}

std::string elblock2vtk( ElementBlock ElBlock, const gsMatrix<> geomCoefs)
{
    // Number of all control points of resulting Bezier elements of this block
    index_t totalPoints = ((ElBlock.PR+1) * (ElBlock.PS+1)) * ElBlock.numElements;
    std::stringstream stream;

    // Control point ID transformation matrix
    // from G+Smo notation to Paraview notation
    gsSparseMatrix<> T = gs2vtk(ElBlock.PR+1, ElBlock.PS+1);

    // Setup matrices with Cell data
    gsMatrix<index_t> degrees(ElBlock.numElements,3);
    degrees.rowwise() = (gsVector3d<index_t>()<< ElBlock.PR, ElBlock.PS, 0).finished().transpose();

    gsMatrix<index_t> connectivity(1,totalPoints);
    connectivity.asVector().setLinSpaced(0,totalPoints-1);

    gsMatrix<index_t> offsets(1,ElBlock.numElements);
    offsets.asVector().setLinSpaced((ElBlock.PR+1) * (ElBlock.PS+1),totalPoints);

    gsMatrix<index_t> types(1,ElBlock.numElements);
    types.setOnes();
    types*= VTK_BEZIER_QUADRILATERAL;

    // Loop over all elements of the Element Block
    auto Ait = ElBlock.actives.begin();        // Actives Iterator
    auto Cit = ElBlock.coefVectors.begin();    // Coefficients Iteratos

    index_t i = 0;
    gsMatrix<real_t> newCoefs(totalPoints, geomCoefs.cols());
    for(; Ait != ElBlock.actives.end() && Cit != ElBlock.coefVectors.end(); ++Ait, ++Cit)
    {
        newCoefs.middleRows(i*(ElBlock.PR+1) * (ElBlock.PS+1),(ElBlock.PR+1) * (ElBlock.PS+1)) =
             T * *Cit * geomCoefs(Ait->asVector(),gsEigen::all);
        // ID transform * bezier extraction operator * original control points
        ++i;
    } 

    // Format to xml
    stream << "<Piece NumberOfPoints=\""<< totalPoints<<"\" NumberOfCells=\""<<ElBlock.numElements<<"\">\n";
    
    // TODO remove
    stream << "<PointData>\n";
    stream << toDataArray( newCoefs, {{"Name","Position"}});
    stream << "</PointData>\n";

    stream << "<CellData HigherOrderDegrees=\"HigherOrderDegrees\">\n";
    stream << "" << toDataArray(degrees,{{"Name","HigherOrderDegrees"},{"NumberOfComponents","3"}});
    stream << "</CellData>\n";

    stream << "<Points>\n";
    stream << ""<< toDataArray( newCoefs, {{"Name","Points"}});
    stream << "</Points>\n";

    stream << "<Cells>\n";
    stream << toDataArray(connectivity, {{"Name","connectivity"}});
    stream << toDataArray(offsets, {{"Name","offsets"}});
    stream << toDataArray(types, {{"Name","types"}});
    stream << "</Cells>\n";
    stream << "</Piece>\n";

    return stream.str(); 
     // TODO: Handle indentation in a neat way
}

int main(int argc, char* argv[])
{
    // gsTensorBSpline<2> geom(*gsNurbsCreator<>::BSplineSquare(10, 0, 0).get());
    // gsTensorBSpline<2> geom(*gsNurbsCreator<>::BSplineFatQuarterAnnulus().get());
    gsFileData<> inFile("surfaces/thbs_face_3levels.xml");
    gsTHBSpline<2> geom;
    inFile.getId(0, geom);

    // geom.degreeElevate(2);
    // geom.uniformRefine(4);

    std::map<index_t, ElementBlock> ElBlocks = BezierOperator<2>(geom.basis() );
    std::ofstream out;
    out.open("./bezier.vtu");

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n"
        << "<UnstructuredGrid>\n";


    for (auto const& block : ElBlocks)
    {
        out << elblock2vtk(block.second, geom.coefs());
    }
    out << "</UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    out.close();

}