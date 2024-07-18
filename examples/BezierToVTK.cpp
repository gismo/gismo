#include <gismo.h>
#include <fstream>
#include <iostream>
#include<gsLsdyna/gsBextFormat.h>

using namespace gismo;


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
    for (auto const& pair : attributes)
    {
        if (pair.first!="")
        stream << pair.first <<"=\""<< pair.second <<"\" ";
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
    for (auto const& pair : attributes)
    {
        if (pair.first!="")
        stream << pair.first <<"=\""<< pair.second <<"\" ";
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
    for (auto const& pair : attributes)
    {
        if (pair.first!="")
        stream << pair.first <<"=\""<< pair.second <<"\" ";
    }
    stream <<">\n";
    stream << num ;
    stream << "\n</DataArray>\n";

    return stream.str();
}

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

    stream << "\t\t<Piece NumberOfPoints=\""<< bezier.coefsSize()<<"\" NumberOfCells=\"1\">\n";
    
    // TODO remove
    stream << "<PointData>\n";
    stream << toDataArray( newCoefs, {{"Name","Position"}});
    stream << "</PointData>\n";

    stream << "\t\t\t<CellData HigherOrderDegrees=\"HigherOrderDegrees\">\n";
    gsMatrix<index_t> degrees(1,3);
    degrees << bezier.degree(0), bezier.degree(1), 0; 
    // TODO: Handle indentation in a neat way
    stream << "\t\t\t\t" << toDataArray(degrees,{{"Name","HigherOrderDegrees"},{"NumberOfComponents","3"}});
    stream << "\t\t\t</CellData>\n";
    stream << "\t\t\t<Points>\n";
    
    stream << "\t\t\t\t"<< toDataArray( newCoefs, {{"Name","Points"}});
    stream << "\t\t\t</Points>\n";

    stream << "\t\t\t<Cells>\n";
    connectivity.asVector().setLinSpaced(0,bezier.coefsSize()-1);
    stream << toDataArray(connectivity, {{"Name","connectivity"}});
    stream << toDataArray(bezier.coefsSize(), {{"Name","offsets"}});
    stream << toDataArray(77, {{"Name","types"}});
    stream << "</Cells>\n";
    stream << "</Piece>\n";

    return stream.str();
}

int main(int argc, char* argv[])
{
    // gsTensorBSpline<2> bezier(*gsNurbsCreator<>::BSplineSquare(10, 0, 0).get());
    // gsTensorBSpline<2> bezier(*gsNurbsCreator<>::BSplineFatQuarterAnnulus().get());
    gsFileData<> inFile("/Users/ck/Documents/ParaViewBezier/filedata/surfaces/thbs_face_3levels.xml");
    gsTHBSpline<2> bezier;
    inFile.getId(0, bezier);

    // bezier.degreeElevate(1);
    // bezier.uniformRefine(2);

    std::map<index_t, ElementBlock> ElBlocks = BezierOperator<2>(bezier.basis() );
    std::ofstream out;
    out.open("/Users/ck/Desktop/bezier.vtu");

    out << "<?xml version=\"1.0\"?>\n"
        << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" header_type=\"UInt32\">\n"
        << "\t<UnstructuredGrid>\n";

    for (auto const& pair : ElBlocks)
    {
        // gsInfo << "Actives size (nctrlpts): " << pair.first << "\n";
        // gsInfo << "NumElements: " << pair.second.numElements << "\n";
        // gsInfo << "Actives:\n" << pair.second.actives.back().transpose()<< "\n";
        // gsInfo << "CoefVecs:\n (" << pair.second.coefVectors.front().rows()
        //         << "," << pair.second.coefVectors.front().cols()<< ")\n";
        // gsInfo << "PR, PS, PT: " << pair.second.PR << "," << pair.second.PS << "," << pair.second.PT <<'\n';
        // gsInfo << pair.second.coefVectors.front() << "\n";

        // for each element
        gsKnotVector<> kvR(pair.second.PR);
        kvR.initClamped(pair.second.PR);
        gsKnotVector<> kvS(pair.second.PS);
        kvS.initClamped(pair.second.PS);
        gsDebugVar(kvR);
        gsDebugVar(kvS);

        gsMatrix<> dummCoeffs=gsMatrix<>::Zero( (pair.second.PR+1)*(pair.second.PS+1), bezier.geoDim());
        gsTensorBSpline<2> bez(kvR,kvS,dummCoeffs);


        auto Ait = pair.second.actives.begin(); // Actives Iterator
        auto Cit = pair.second.coefVectors.begin();    // Coefficients Iteratos

        for(; Ait != pair.second.actives.end() && Cit != pair.second.coefVectors.end(); ++Ait, ++Cit)
        {
            bez.coefs() = *Cit * bezier.coefs()(Ait->asVector(),gsEigen::all);
            out << bez2vtk(bez) << "\n";

        }
    }
    out << "</UnstructuredGrid>\n";
    out << "</VTKFile>\n";
    out.close();

/* 
    gsMatrix<> dumb(bezier.basis().size(), 1);
    dumb.asVector().setLinSpaced(1,dumb.rows() );
    gsInfo << dumb.transpose() << "\n";


    std::ofstream out;
    out.open("/Users/ck/Desktop/bezier.vtu");
    out << bez2vtk(bezier);
    out.close();
*/
}