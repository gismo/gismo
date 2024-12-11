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
        // Transform converts coefs from G+Smo's convention to ParaView's convention
        gsMatrix<real_t> Transform(nU*nV, nU*nV);
        Transform.setZero();

        // Transform( Paraview , gismo  )
        // Corners ( always 0-3 )
        Transform(0,0) = 1;
        Transform(1,nU-1) = 1;
        Transform(2,nU*nV-1) = 1;
        Transform(3,nU*nV-nU) = 1;

        // Edges
        for (index_t i=1;i<nU-1;++i) // Parallel to u
        {
            Transform(3+i,i) = 1;
            Transform(3+ (nU-2) + (nV-2) + i, nU*(nV-1) + i  ) = 1;


        }
        for (index_t j=1;j<nV-1;++j) // Parallel to v
        {
            Transform( 3 + (nU-2) + j, (j+1)*nU-1) = 1 ;
            Transform( 3 + 2*(nU-2) + (nV-2) + j, nU*j) = 1;
        }
        // Internal
        for (index_t i=0;i<nU-2;++i)
        {
            for (index_t j=0;j<nV-2;++j)
            {
                Transform(2*( nU + nV) - 4 + j*(nU-2) + i , nU*(j+1)+(i+1)) = 1;
            }
        }
        return Transform;
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

} // namespace gismo

#undef VTK_BEZIER_QUADRILATERAL
