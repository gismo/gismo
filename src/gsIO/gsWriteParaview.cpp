/** @file gsWriteParaview.cpp

    @brief Utility for plotting error / obsolete - to be removed

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): G. Kiss, A. Mantzaflaris
*/

#include <gsIO/gsWriteParaview.h>
#include <gsMesh2/gsSurfMesh.h>
#include <gsIO/gsParaviewCollection.h>

#include <fstream>
#include <initializer_list>

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



#define PLOT_PRECISION 12

template<class Scalar>
void gsWriteParaview(gsSurfMesh<Scalar> const & sm,
                     std::string const & fn)
{
    std::vector<std::string> pname = sm.vertex_properties();
    gsWriteParaview(sm, fn, pname);
}

template<class Scalar>
void gsWriteParaview(const gsSurfMesh<Scalar> & sm,
                     std::string const & fn,
                     std::vector<std::string> props)
{
    using MeshT = gsSurfMesh<Scalar>;

    std::string mfn(fn);
    mfn.append(".vtk");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        gsWarn<<"gsWriteParaview: Problem opening file \""<<fn<<"\""<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    //https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    file << "# vtk DataFile Version 4.2\n";
    file << "vtk output\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    // Vertices
    auto vpt = sm.template get_vertex_property<typename MeshT::Point>("v:point");
    file << "POINTS " << sm.n_vertices() << " float\n";
    for (auto v : sm.vertices() )
        file << vpt[v].transpose() <<"\n";
    file << "\n";

    // Triangles or quads
    file << "POLYGONS " << sm.n_faces() << " " <<
        sm.face_valence_sum() + sm.n_faces() << "\n";
    for (auto f : sm.faces())
    {
        file << sm.valence(f) <<" "; //3: triangles, 4: quads
        for (auto v : sm.vertices(f))
            file << v.idx() << " ";
        file << "\n";
    }
    file << "\n";

    //todo: count props starting with v:, f:, e:
    if (0!=props.size())
        file << "POINT_DATA " << sm.n_vertices() << "\n";//once
    for( auto & pr : propvec )
    {
        if (pr == "v:connectivity") continue;
        if (pr == "v:deleted") continue;

        if (pr == "v:normal")
        {
            auto vn = sm.template get_vertex_property<typename MeshT::Point>(pr);
            GISMO_ASSERT(vn,"No normals found");
            file << "NORMALS "<<pr<<" float\n";
            for (auto v : sm.vertices() )
                file << vn[v].transpose() <<"\n";
            file << "\n";
            continue;
        }

        auto vp = sm.template get_vertex_property<typename MeshT::Point>(pr);
        if (vp)
        {
            file << "VECTORS "<<pr<<" float\n";
            for (auto v : sm.vertices() )
                file << vp[v].transpose() <<"\n";
            file << "\n";
            continue;
        }

        auto vs = sm.template get_vertex_property<Scalar>(pr);
        if (vs)
        {
            file << "SCALARS "<<pr<<" float\nLOOKUP_TABLE default\n";
            for (auto v : sm.vertices() )
                file << vs[v] <<" ";
            file << "\n";
            continue;
        }

        auto vi = sm.template get_vertex_property<index_t>(pr);
        if (vi)
        {
            file << "SCALARS "<<pr<<" float\nLOOKUP_TABLE default\n";
            for (auto v : sm.vertices() )
                file << vi[v] <<" ";
            file << "\n";
            continue;
        }

        auto vb = sm.template get_vertex_property<bool>(pr);
        if (vb)
        {
            file << "SCALARS "<<pr<<" float\nLOOKUP_TABLE default\n";
            for (auto v : sm.vertices() )
                file << vb[v] <<" ";
            file << "\n";
            continue;
        }

        const std::type_info & ti = sm.get_vertex_property_type(pr);
        gsWarn<< "gsWriteParaview: Property "<< pr << " ignored, "<<ti.name()<<".\n";
    }

    file.close();
    //makeCollection(fn, ".vtk"); // legacy inside pvd seems to not work
}

TEMPLATE_INST
void gsWriteParaview<real_t>(const gsSurfMesh<real_t> &,
                     std::string const &,
                     std::vector<std::string>);

}//namespace gismo
