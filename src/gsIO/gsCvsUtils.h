// I/O as comma-separated values

#pragma once

#include <gsIO/gsFileData.h>


void writeCVS( const gsMatrix<> & mat, std::string const & filename )
{
    std::string tmp = gsFileData<>::getExtension(filename);
    if (tmp != "cvs" )
        tmp = filename + ".cvs";
    else
        tmp = filename;

    std::ofstream file(tmp.c_str());

    for ( index_t i = 0 ; i != mat.rows(); ++i )  
    {
        file << mat(i,0);
        for ( index_t j = 1 ; j != mat.cols(); ++j )  
            file << ", " << mat(i,j) ;
        file<< "\n";
    }
    file.close();
}
