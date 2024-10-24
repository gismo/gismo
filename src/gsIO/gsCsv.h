/** @file gsCsv.h

    @brief Provides functions writing .csv files.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): C. Karampatzakis
*/

#include <gsCore/gsForwardDeclarations.h>

#include <fstream>

#pragma once


namespace gismo {

/// @brief Export a \a gsMatrix to a .csv (comma separated values) file
/// @tparam T 
/// @param filename path of output file
/// @param matrix a \a gsMatrix to be written to the file
/// @param headers optionally, a vector of strings to be used as column headers 
///
/// \ingroup IO
template<class T>
void gsWriteCsv(std::string const & filename, const gsMatrix<T> & matrix, const std::vector<std::string> & headers = std::vector<std::string>() )
{
    // Define format of .csv file.
    const static gsEigen::IOFormat CSVFormat(gsEigen::FullPrecision, gsEigen::Aligned, ", ", "\n");

    std::ofstream csv_file;
    csv_file.open(filename);
    GISMO_ASSERT( (headers.empty() ||  headers.size() == matrix.cols()), "The column headers should be as many as the columns of the matrix provided." );

    // If column headers are provided, write to file
    if (! headers.empty())
    {
        for ( index_t j=0 ; j <  headers.size() ; j++)
        {
            csv_file << headers[j];
            if (headers.size()-1 == j) 
                csv_file << CSVFormat.rowSeparator;
            else
                csv_file << CSVFormat.coeffSeparator ;
        }
    }
    // write matrix entries to file 
    csv_file << matrix.format(CSVFormat);
    csv_file.close();
}


} // namespace gismo
