/** @file gsWriteParaview.cpp

    @brief Utility for plotting error / obsolete - to be removed

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): G. Kiss
*/

#include <gsCore/gsTemplateTools.h>
#include <gsCore/gsLinearAlgebra.h>

#include <fstream>

namespace gismo
{

template <class T>
void plot_errors(const gsMatrix<T> & orig, 
                 const gsMatrix<T> & comp, const std::vector<T> & errors, 
                 std::string const & fn)
{
    std::string mfn(fn);
    mfn.append(".vtu");
    std::ofstream myfile;
    myfile.open (mfn.c_str());
    if (myfile.is_open()){
        myfile<<"<?xml version=\"1.0\"?>\n";
        myfile <<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
        myfile <<"<UnstructuredGrid>\n";
        myfile << "<Piece NumberOfPoints=\""<< 2*errors.size() <<"\" NumberOfCells=\""<< errors.size() <<"\">"<<"\n";
        myfile <<"<PointData Scalars=\"scalars\">\n";
        myfile <<"<DataArray type=\"Float32\" Name=\"SolutionField\" format=\"ascii\">\n";
        for(size_t j = 0; j < errors.size();j++){
            myfile<< errors[j] <<" ";
        }
        for(size_t j = 0; j < errors.size();j++){
            myfile<< errors[j] <<" ";
        }
        myfile <<"</DataArray>\n";
        myfile <<"</PointData>\n";

        myfile <<"<CellData Scalars=\"scalars\">\n";
        myfile <<"<DataArray type=\"Float32\" Name=\"CellData\" format=\"ascii\">\n";
        for(size_t j = 0; j < errors.size();j++){
            myfile<< 1 <<" ";
        }
        myfile <<"</DataArray>\n";
        myfile <<"</CellData>\n";

        myfile <<"<Points>\n";
        myfile <<"<DataArray type=\"Float32\" NumberOfComponents=\""<< 3 <<"\" format=\"ascii\">\n";
        for(unsigned int i =0; i < errors.size();i++){
            myfile<<orig(0,i)<<" " << orig(1,i)<<" "<<orig(2,i)<<"\n";
        }
        for(unsigned int i =0; i < errors.size();i++){
            myfile<<comp(0,i)<<" " << comp(1,i)<<" "<<comp(2,i)<<"\n";
        }
        myfile <<"</DataArray>\n";
        myfile <<"</Points>\n";

        myfile <<"<Cells>\n";
        myfile <<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
        for(unsigned int i =0; i < errors.size(); i++){
            myfile<< i<<" "<< i + errors.size()<<" ";
        }
        myfile <<"</DataArray>\n";
        myfile <<"<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
        //<DataArray type="Int32" Name="offsets" Format="ascii">
        int counter =0;
        for(unsigned int i =0; i < errors.size(); i++){
            counter += 2;
            myfile<<counter<<" ";
        }
        myfile <<"</DataArray>\n";
        myfile <<"<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
        for(unsigned int i =0; i < errors.size(); i++){
            myfile<< 3 <<" ";
        }
        myfile <<"</DataArray>\n";
        myfile <<"</Cells>\n";
        myfile <<"</Piece>\n";
        myfile <<"</UnstructuredGrid>\n";

        myfile <<"</VTKFile>\n";
        myfile.close();
    }

}


TEMPLATE_INST 
void plot_errors<real_t>(const gsMatrix<real_t>&,  
                         const gsMatrix<real_t>&, 
                         const std::vector<real_t>&,
                         std::string const&); 

}



