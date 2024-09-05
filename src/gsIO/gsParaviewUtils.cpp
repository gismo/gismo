/** @file gsParaviewUtils.h

    @brief ParaView output Utilities

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Zwar, C. Karampatzakis
*/

#include <gsIO/gsParaviewUtils.h>
#include <fstream>
#include <iostream>

#define VTK_BEZIER_QUADRILATERAL 77


namespace gismo
{
    gsMatrix<real_t> vtkIDTransform(index_t nU, index_t nV)
    {
        // T converts coefs from G+Smo's convetnion to ParaView's convention
        gsMatrix<real_t> T(nU*nV, nU*nV);
        T.setZero();

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
    


    /// @brief Converts an integer to a 'DataArray' xml tag, which is returned as a string.
    /// @param num The integer to be formatted
    /// @param attributes Optional, map of strings, with attribute name mapping to attribute value.
    /// @return 
    std::string toDataArray(index_t num, std::map<std::string, std::string> attributes)
    {
        std::stringstream stream;

        stream.setf(std::ios::fixed);  // write floating point values in fixed-point notation.
        // stream.precision(5);
        // stream << std::setfill('0') << std::setw(16);
        // Format as vtk xml string
        stream <<"<DataArray type=\"Float32\" format=\"ascii\" ";
        for (auto const& block : attributes)
        {
            if (block.first!="")
            stream << block.first <<"=\""<< block.second <<"\" ";
        }
        stream <<">\n" << num << "\n</DataArray>\n";

        return stream.str();
    }

    std::string ΕlBlock2vtk( ElementBlock ElBlock, const gsMatrix<> geomCoefs)
    {
        // Number of all control points of resulting Bezier elements of this block
        index_t totalPoints = ((ElBlock.PR+1) * (ElBlock.PS+1)) * ElBlock.numElements;
        std::stringstream stream;

        // Control point ID transformation matrix
        // from G+Smo notation to Paraview notation
        gsMatrix<> T = vtkIDTransform(ElBlock.PR+1, ElBlock.PS+1);

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
        stream <<"<Piece NumberOfPoints=\""<< totalPoints<<"\" NumberOfCells=\""<<ElBlock.numElements<<"\">\n";
        // TODO remove
        stream << "<PointData>\n";
        // newCoefs.transposeInPlace();
        stream << toDataArray( newCoefs.transpose(), {{"Name","Position"}});
        stream << "</PointData>\n";

        stream << "<CellData HigherOrderDegrees=\"HigherOrderDegrees\">\n";
        degrees.transposeInPlace();
        stream << "" << toDataArray(degrees,{{"Name","HigherOrderDegrees"}});
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
    }



} // namespace gismo

#undef VTK_BEZIER_QUADRILATERAL
