#include <gsIO/gsWriteParaview.h>

#include <gsIO/gsParaviewCollection.h>

#include <gsThbs/gsTHBSpline.h>

namespace gismo
{

template <class T>
void plot_errors(gsMatrix<T> orig, 
                 gsMatrix<T> comp, std::vector<T> errors, 
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
        myfile << "<Piece NumberOfPoints=\""<< 2*errors.size() <<"\" NumberOfCells=\""<< errors.size() <<"\">"<<std::endl;
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
            myfile<<orig(0,i)<<" " << orig(1,i)<<" "<<orig(2,i)<<std::endl;
        }
        for(unsigned int i =0; i < errors.size();i++){
            myfile<<comp(0,i)<<" " << comp(1,i)<<" "<<comp(2,i)<<std::endl;
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
void plot_errors<real_t>(gismo::gsMatrix<real_t>, 
                         gismo::gsMatrix<real_t>, 
                         std::vector<real_t>, std::string const&); 

}
