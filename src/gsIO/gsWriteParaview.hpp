/** @file gsWriteParaview.hpp

    @brief Provides implementation of functions writing Paraview files.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/


#pragma once

#include <gsIO/gsWriteParaview.h>
#include <gsIO/gsParaviewCollection.h>
#include <gsIO/gsIOUtils.h>

#include <gsCore/gsGeometry.h>
#include <gsCore/gsGeometrySlice.h>
#include <gsCore/gsField.h>
#include <gsCore/gsDebug.h>

#include <gsModeling/gsTrimSurface.h>
#include <gsModeling/gsSolid.h>
//#include <gsUtils/gsMesh/gsHeMesh.h>


#define PLOT_PRECISION 5

namespace gismo
{

// Export a 3D parametric mesh
template<class T>
void writeSingleBasisMesh3D(const gsMesh<T> & sl,
                            std::string const & fn)
{
    const unsigned numVer = sl.numVertices();
    const unsigned numEl  = numVer / 8;
    std::string mfn(fn);
    mfn.append(".vtu");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        gsWarn<<"writeSingleBasisMesh3D: Problem opening file \""<<fn<<"\""<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file <<"<UnstructuredGrid>\n";

    // Number of vertices and number of cells
    file <<"<Piece NumberOfPoints=\""<< numVer <<"\" NumberOfCells=\""<<numEl<<"\">\n";

    // Coordinates of vertices
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (typename std::vector< gsVertex<T>* >::const_iterator it=sl.vertices().begin(); it!=sl.vertices().end(); ++it)
    {
        const gsVertex<T>& vertex = **it;
        file << vertex[0] << " " << vertex[1] << " " << vertex[2] << " \n";
    }
    file << "\n";
    file <<"</DataArray>\n";
    file <<"</Points>\n";

    // Point data
    file <<"<PointData Scalars=\"CellVolume\">\n";
    file <<"<DataArray type=\"Float32\" Name=\"CellVolume\" format=\"ascii\" NumberOfComponents=\"1\">\n";
    for (typename std::vector< gsVertex<T>* >::const_iterator it=sl.vertices().begin(); it!=sl.vertices().end(); ++it)
    {
        file << (*it)->data <<" ";
    }
    file << "\n";
    file <<"</DataArray>\n";
    file <<"</PointData>\n";

    // Cells
    file <<"<Cells>\n";

    // Connectivity
    file <<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (unsigned i = 0; i!= numVer;++i)
    {
        file << i << " ";
    }
    file << "\n";
    file <<"</DataArray>\n";

    // Offsets
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (unsigned i = 1; i<= numEl;++i)
    {
        file << 8*i << " ";
    }
    file << "\n";
    file << "</DataArray>\n";

    // Type
    file << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
    for (unsigned i = 1; i<= numEl;++i)
    {
        file << "11 ";
    }
    file << "\n";
    file << "</DataArray>\n";

    file <<"</Cells>\n";
    file << "</Piece>\n";
    file <<"</UnstructuredGrid>\n";
    file <<"</VTKFile>\n";
    file.close();

    //if( pvd ) // make a pvd file
    //    makeCollection(fn, ".vtp");
}

// Export a 2D parametric mesh -- note: duplicates code from writeSingleBasisMesh3D,
//
template<class T>
void writeSingleBasisMesh2D(const gsMesh<T> & sl,
                            std::string const & fn)
{
    const unsigned numVer = sl.numVertices();
    const unsigned numEl  = numVer / 4; //(1<<dim)
    std::string mfn(fn);
    mfn.append(".vtu");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        gsWarn<<"writeSingleBasisMesh2D: Problem opening file \""<<fn<<"\""<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file <<"<UnstructuredGrid>\n";

    // Number of vertices and number of cells
    file <<"<Piece NumberOfPoints=\""<< numVer <<"\" NumberOfCells=\""<<numEl<<"\">\n";

    // Coordinates of vertices
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (typename std::vector< gsVertex<T>* >::const_iterator it=sl.vertices().begin(); it!=sl.vertices().end(); it+=4)
    {
        // order is important!
        const gsVertex<T>& vertex0 = **it;
        const gsVertex<T>& vertex1 = **(it+1);
        const gsVertex<T>& vertex3 = **(it+3);
        const gsVertex<T>& vertex2 = **(it+2);
        file << vertex0[0] << " " << vertex0[1] << " " << vertex0[2] << " \n";
        file << vertex1[0] << " " << vertex1[1] << " " << vertex1[2] << " \n";
        file << vertex3[0] << " " << vertex3[1] << " " << vertex3[2] << " \n";
        file << vertex2[0] << " " << vertex2[1] << " " << vertex2[2] << " \n";
    }
    file << "\n";
    file <<"</DataArray>\n";
    file <<"</Points>\n";

    // Point data
    file <<"<PointData Scalars=\"CellArea\">\n";
    file <<"<DataArray type=\"Float32\" Name=\"CellVolume\" format=\"ascii\" NumberOfComponents=\"1\">\n";
    for (typename std::vector< gsVertex<T>* >::const_iterator it=sl.vertices().begin(); it!=sl.vertices().end(); ++it)
    {
        file << (*it)->data <<" ";
    }
    file << "\n";
    file <<"</DataArray>\n";
    file <<"</PointData>\n";

    // Cells
    file <<"<Cells>\n";

    // Connectivity
    file <<"<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (unsigned i = 0; i!= numVer;++i)
    {
        file << i << " ";
    }
    file << "\n";
    file <<"</DataArray>\n";

    // Offsets
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (unsigned i = 1; i<= numEl;++i)
    {
        file << 4*i << " "; //step: (1<<dim)
    }
    file << "\n";
    file << "</DataArray>\n";

    // Type
    file << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
    for (unsigned i = 1; i<= numEl;++i)
    {
        file << "9 ";// 11: 3D, 9: 2D
    }
    file << "\n";
    file << "</DataArray>\n";

    file <<"</Cells>\n";
    file << "</Piece>\n";
    file <<"</UnstructuredGrid>\n";
    file <<"</VTKFile>\n";
    file.close();

    //if( pvd ) // make a pvd file
    //    makeCollection(fn, ".vtp");
}


/// Export a parametric mesh
template<class T>
void writeSingleBasisMesh(const gsBasis<T> & basis,
                         std::string const & fn)
{
    gsMesh<T> msh(basis, 0);
    if ( basis.dim() == 3)
        writeSingleBasisMesh3D(msh,fn);
    else if ( basis.dim() == 2)
        writeSingleBasisMesh2D(msh,fn);
    else
        gsWriteParaview(msh, fn, false);
}

/// Export a computational mesh
template<class T>
void writeSingleCompMesh(const gsBasis<T> & basis, const gsGeometry<T> & Geo,
                         std::string const & fn, unsigned resolution = 8)
{
    gsMesh<T> msh(basis, resolution);
    Geo.evaluateMesh(msh);

    // if ( basis.dim() == 3)
    //     writeSingleBasisMesh3D(msh,fn);
    // else if ( basis.dim() == 2)
    //     writeSingleBasisMesh2D(msh,fn);
    // else
        gsWriteParaview(msh, fn, false);
}

/// Export a control net
template<class T>
void writeSingleControlNet(const gsGeometry<T> & Geo,
                           std::string const & fn)
{
    const int d = Geo.parDim();
    gsMesh<T> msh;
    Geo.controlNet(msh);
    const unsigned n = Geo.geoDim();
    if ( n == 1 )
    {
        gsMatrix<T> anch = Geo.basis().anchors();
        // Lift vertices at anchor positions
        for (size_t i = 0; i!= msh.numVertices(); ++i)
        {
            msh.vertex(i)[d] = msh.vertex(i)[0];
            msh.vertex(i).topRows(d) = anch.col(i);
        }
    }
    else if (n>3)
    {
        gsDebug<<"Writing 4th coordinate\n";
        const gsMatrix<T> & cp = Geo.coefs();
        gsWriteParaviewPoints<T>(cp.transpose(), fn );
        return;
    }

    gsWriteParaview(msh, fn, false);
}

template<class T>
void gsWriteParaviewTPgrid(const gsMatrix<T> & eval_geo  ,
                           const gsMatrix<T> & eval_field,
                           const gsVector<index_t> & np,
                           std::string const & fn)
{
    const int n = eval_geo.rows();
    GISMO_ASSERT(eval_geo.cols()==eval_field.cols()
                 && static_cast<index_t>(np.prod())==eval_geo.cols(),
                 "Data do not match");

    std::string mfn(fn);
    mfn.append(".vts");
    std::ofstream file(mfn.c_str());
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
    file <<"<StructuredGrid WholeExtent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "
         << (np.size()>2 ? np(2)-1 : 0) <<"\">\n";
    file <<"<Piece Extent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "
         << (np.size()>2 ? np(2)-1 : 0) <<"\">\n";
    file <<"<PointData "<< ( eval_field.rows()==1 ?"Scalars":"Vectors")<<"=\"SolutionField\">\n";
    file <<"<DataArray type=\"Float32\" Name=\"SolutionField\" format=\"ascii\" NumberOfComponents=\""<< ( eval_field.rows()==1 ? 1 : 3) <<"\">\n";
    if ( eval_field.rows()==1 )
        for ( index_t j=0; j<eval_field.cols(); ++j)
            file<< eval_field.at(j) <<" ";
    else
    {
        for ( index_t j=0; j<eval_field.cols(); ++j)
        {
            for ( index_t i=0; i!=eval_field.rows(); ++i)
                file<< eval_field(i,j) <<" ";
            for ( index_t i=eval_field.rows(); i<3; ++i)
                file<<"0 ";
        }
    }
    file <<"</DataArray>\n";
    file <<"</PointData>\n";
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
    for ( index_t j=0; j<eval_geo.cols(); ++j)
    {
        for ( index_t i=0; i!=n; ++i)
            file<< eval_geo(i,j) <<" ";
        for ( index_t i=n; i<3; ++i)
            file<<"0 ";
    }
    file <<"</DataArray>\n";
    file <<"</Points>\n";
    file <<"</Piece>\n";
    file <<"</StructuredGrid>\n";
    file <<"</VTKFile>\n";

    file.close();
}

template<class T>
void writeSinglePatchField(const gsFunction<T> & geometry,
                           const gsFunction<T> & parField,
                           const bool isParam,
                           std::string const & fn, unsigned npts)
{
    const int n = geometry.targetDim();
    const int d = geometry.domainDim();

    gsMatrix<T> ab = geometry.support();
    gsVector<T> a = ab.col(0);
    gsVector<T> b = ab.col(1);

    gsVector<unsigned> np = uniformSampleCount(a, b, npts);
    gsMatrix<T> pts = gsPointGrid(a, b, np);

    gsMatrix<T> eval_geo = geometry.eval(pts);//pts
    gsMatrix<T>  eval_field = isParam ? parField.eval(pts) : parField.eval(eval_geo);

    if ( 3 - d > 0 )
    {
        np.conservativeResize(3);
        np.bottomRows(3-d).setOnes();
    }
    else if (d > 3)
    {
        gsWarn<< "Cannot plot 4D data.\n";
        return;
    }

    if ( 3 - n > 0 )
    {
        eval_geo.conservativeResize(3,eval_geo.cols() );
        eval_geo.bottomRows(3-n).setZero();
    }
    else if (n > 3)
    {
        gsWarn<< "Data is more than 3 dimensions.\n";
    }

   if ( eval_field.rows() == 2)
    {
        eval_field.conservativeResize(3,eval_geo.cols() );
        eval_field.bottomRows(1).setZero(); // 3-field.dim()
    }

    gsWriteParaviewTPgrid(eval_geo, eval_field, np.template cast<index_t>(), fn);
}

/// Write a file containing a solution field over a single geometry
template<class T>
void writeSinglePatchField(const gsField<T> & field, int patchNr,
                           std::string const & fn, unsigned npts)
{
    writeSinglePatchField(field.patch(patchNr), field.function(patchNr), field.isParametric(), fn, npts);
/*
    const int n = field.geoDim();
    const int d = field.parDim();

    gsMatrix<T> ab = field.patches().parameterRange(patchNr);
    gsVector<T> a = ab.col(0);
    gsVector<T> b = ab.col(1);

    gsVector<unsigned> np = uniformSampleCount(a, b, npts);
    gsMatrix<T> pts = gsPointGrid(a, b, np);

    gsMatrix<T> eval_geo = field.point ( pts, patchNr );//pts

    if ( 3 - d > 0 )
    {
        np.conservativeResize(3);
        np.bottomRows(3-d).setOnes();
    }
    else if (d > 3)
    {
        gsWarn<< "Cannot plot 4D data.\n";
        return;
    }

    if ( 3 - n > 0 )
    {
        eval_geo.conservativeResize(3,eval_geo.cols() );
        eval_geo.bottomRows(3-n).setZero();
    }
    else if (d > 3)
    {
        gsWarn<< "Cannot plot 4D data.\n";
        return;
    }

    gsMatrix<T>  eval_field = field.value ( pts, patchNr );//values
    GISMO_ASSERT( eval_field.rows() == field.dim(), "Error in field dimension");
    if ( eval_field.rows() > 1 )
    {
        eval_field.conservativeResize(3,eval_geo.cols() );
        eval_field.bottomRows( 3-field.dim() ).setZero();
    }

    std::string mfn(fn);
    mfn.append(".vts");
    std::ofstream file(mfn.c_str());
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
    file <<"<StructuredGrid WholeExtent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
    file <<"<Piece Extent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
    file <<"<PointData "<< ( field.dim()==1 ?"Scalars":"Vectors")<<"=\"SolutionField\">\n";
    file <<"<DataArray type=\"Float32\" Name=\"SolutionField\" format=\"ascii\" NumberOfComponents=\""<< eval_field.rows() <<"\">\n";
    for ( index_t j=0; j<eval_field.cols(); ++j)
        for ( index_t i=0; i<eval_field.rows(); ++i)
            file<< eval_field(i,j) <<" ";
    file <<"</DataArray>\n";
    file <<"</PointData>\n";
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\""<<eval_geo.rows()<<"\">\n";
    for ( index_t j=0; j<eval_geo.cols(); ++j)
        for ( index_t i=0; i<eval_geo.rows(); ++i)
            file<< eval_geo(i,j) <<" ";
    file <<"</DataArray>\n";
    file <<"</Points>\n";
    file <<"</Piece>\n";
    file <<"</StructuredGrid>\n";
    file <<"</VTKFile>\n";

    file.close();
*/
}

/// Export a geometry represented by \a func
template<class T>
void writeSingleGeometry(gsFunction<T> const& func,
                         gsMatrix<T> const& supp,
                         std::string const & fn, unsigned npts)
{
    int n = func.targetDim();
    const int d = func.domainDim();

    gsVector<T> a = supp.col(0);
    gsVector<T> b = supp.col(1);
    gsVector<unsigned> np = uniformSampleCount(a,b, npts );
    gsMatrix<T> pts = gsPointGrid(a,b,np) ;

    gsMatrix<T>  eval_func = func.eval  ( pts ) ;//pts

    if ( 3 - d > 0 )
    {
        np.conservativeResize(3);
        np.bottomRows(3-d).setOnes();
    }

    if ( 3 - n > 0 )
    {
        eval_func.conservativeResize(3,eval_func.cols() );
        eval_func.bottomRows(3-n).setZero();

        if ( n == 1 )
        {
            if (d==3)
            {
                n = 4;
                eval_func.conservativeResize(4,eval_func.cols() );
            }

            //std::swap( eval_geo.row(d),  eval_geo.row(0) );
            eval_func.row(d) = eval_func.row(0);
            eval_func.topRows(d) = pts;
        }
    }

    std::string mfn(fn);
    mfn.append(".vts");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        gsWarn<<"writeSingleGeometry: Problem opening file \""<<fn<<"\""<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);
    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
    file <<"<StructuredGrid WholeExtent=\"0 "<<np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
    file <<"<Piece Extent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
    // Add norm of the point as data
    // file <<"<PointData Scalars =\"PointNorm\">\n";
    // file <<"<DataArray type=\"Float32\" Name=\"PointNorm\" format=\"ascii\" NumberOfComponents=\""<< 1 <<"\">\n";
    // for ( index_t j=0; j<eval_geo.cols(); ++j)
    //     file<< eval_geo.col(j).norm() <<" ";
    // file <<"</DataArray>\n";
    // file <<"</PointData>\n";
    // end norm

    //---------
    if (n > 3)
    {
        //gsWarn<< "4th dimension as scalar data.\n";
        file <<"<PointData "<< "Scalars=\"Coordinate4\">\n";
        file <<"<DataArray type=\"Float32\" Name=\"Coordinate4\" format=\"ascii\" NumberOfComponents=\"1\">\n";
        for ( index_t j=0; j!=eval_func.cols(); ++j)
            file<< eval_func(3,j) <<" ";
        file <<"</DataArray>\n";
        file <<"</PointData>\n";
    }
    //---------

    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
    for ( index_t j=0; j<eval_func.cols(); ++j)
        for ( index_t i=0; i!=3; ++i)
            file<< eval_func(i,j) <<" ";
    file <<"</DataArray>\n";
    file <<"</Points>\n";
    file <<"</Piece>\n";
    file <<"</StructuredGrid>\n";
    file <<"</VTKFile>\n";
    file.close();
}

/// Export a curve geometry represented by \a func
template<class T>
void writeSingleCurve(gsFunction<T> const& func,
                      gsMatrix<T> const& supp,
                      std::string const & fn, unsigned npts)
{
    const unsigned n = func.targetDim();
    const unsigned d = func.domainDim();
    GISMO_ASSERT( d == 1, "Not a curve");

    gsVector<T> a = supp.col(0);
    gsVector<T> b = supp.col(1);
    gsVector<unsigned> np = uniformSampleCount(a,b, npts );
    gsMatrix<T> pts = gsPointGrid(a,b,np) ;

    gsMatrix<T>  eval_func = func.eval  ( pts ) ;//pts

    np.conservativeResize(3);
    np.bottomRows(2).setOnes();

    if ( 3 - n > 0 )
    {
        eval_func.conservativeResize(3,eval_func.cols() );
        eval_func.bottomRows(3-n).setZero();

        if ( n == 1 )
        {
            //std::swap( eval_geo.row(d),  eval_geo.row(0) );
            eval_func.row(d) = eval_func.row(0);
            eval_func.topRows(d) = pts;
        }
    }

    std::string mfn(fn);
    mfn.append(".vtp");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        gsWarn<<"writeSingleCurve: Problem opening file \""<<fn<<"\""<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);
    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file <<"<PolyData>\n";
    // Accounting
    file <<"<Piece NumberOfPoints=\""<< npts
         <<"\" NumberOfVerts=\"0\" NumberOfLines=\""<< npts-1
         <<"\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\""<<eval_func.rows()<<"\">\n";
    for ( index_t j=0; j<eval_func.cols(); ++j)
        for ( index_t i=0; i<eval_func.rows(); ++i)
            file<< eval_func(i,j) <<" ";
    file <<"\n</DataArray>\n";
    file <<"</Points>\n";
    // Lines
    file <<"<Lines>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<npts-1<<"\">\n";
    for (unsigned i=0; i< npts-1; ++i )
    {
        file << i << " " << i+1 << " ";
    }
    // offsets
    file <<"\n</DataArray>\n";
    unsigned offset(0);
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<npts-1<<"\">\n";
    for (unsigned i=0; i< npts-1; ++i )
    {
        offset +=2;
        file << offset << " ";
    }
    file <<"\n</DataArray>\n";
    file <<"</Lines>\n";
    // Closing
    file <<"</Piece>\n";
    file <<"</PolyData>\n";
    file <<"</VTKFile>\n";
    file.close();
}

template<class T>
void writeSingleCurve(const gsGeometry<T> & Geo, std::string const & fn, unsigned npts)
{
    gsMatrix<T> ab = Geo.parameterRange();
    writeSingleCurve( Geo, ab, fn, npts);
}

template<class T>
void writeSingleGeometry(const gsGeometry<T> & Geo, std::string const & fn, unsigned npts)
{
    /*
      gsMesh<T> msh;
      Geo.toMesh(msh, npts);
      gsWriteParaview(msh, fn, false);
      return;
    //*/
    gsMatrix<T> ab = Geo.parameterRange();
    writeSingleGeometry( Geo, ab, fn, npts);
}

template<class T>
void writeSingleTrimSurface(const gsTrimSurface<T> & surf,
                            std::string const & fn,
                            unsigned npts)
{
    typename gsMesh<T>::uPtr msh = surf.toMesh(npts);
    gsWriteParaview( *msh, fn);
}

/// Write a file containing a solution field over a geometry
template<class T>
void gsWriteParaview(const gsField<T> & field,
                     std::string const & fn,
                     unsigned npts, bool mesh)
{
    /*
    if (mesh && (!field.isParametrized()) )
    {
        gsWarn<< "Cannot plot mesh from non-parametric field.";
        mesh = false;
    }
    */

    const unsigned n = field.nPieces();
    gsParaviewCollection collection(fn);
    std::string fileName;

    for ( unsigned i=0; i < n; ++i )
    {
        const gsBasis<T> & dom = field.isParametrized() ?
            field.igaFunction(i).basis() : field.patch(i).basis();

        fileName = fn + util::to_string(i);
        writeSinglePatchField( field, i, fileName, npts );
        collection.addPart(fileName, ".vts");
        if ( mesh )
        {
            fileName+= "_mesh";
            writeSingleCompMesh(dom, field.patch(i), fileName);

            collection.addPart(fileName, ".vtp");
        }

    }
    collection.save();
}


/// Export a Geometry without scalar information
template<class T>
void gsWriteParaview(const gsGeometry<T> & Geo, std::string const & fn,
                     unsigned npts, bool mesh, bool ctrlNet)
{
    const bool curve = ( Geo.domainDim() == 1 );

    gsParaviewCollection collection(fn);

    if ( curve )
    {
        writeSingleCurve(Geo, fn, npts);
        collection.addPart(fn, ".vtp");
    }
    else
    {
        writeSingleGeometry(Geo, fn, npts);
        collection.addPart(fn, ".vts");
    }

    if ( mesh ) // Output the underlying mesh
    {
        const std::string fileName = fn + "_mesh";

	int ptsPerEdge;

	// If not using default, compute the resolution from npts.
	if(npts!=8)
	{
	    const T evalPtsPerElem = npts * (1.0 / Geo.basis().numElements());

	    // The following complicated formula should ensure similar
	    // resolution of the mesh edges and the surface. The
	    // additional multiplication by deg - 1 ensures quadratic
	    // elements to be approximated by at least two lines etc.
	    ptsPerEdge = cast<T,int>(
            math::max(Geo.basis().maxDegree()-1, (index_t)1) * math::pow(evalPtsPerElem, T(1.0)/Geo.domainDim()) );
	}
	else
	{
	    ptsPerEdge = npts;
	}

        writeSingleCompMesh(Geo.basis(), Geo, fileName, ptsPerEdge);
        collection.addPart(fileName, ".vtp");
    }

    if ( ctrlNet ) // Output the control net
    {
        const std::string fileName = fn + "_cnet";
        writeSingleControlNet(Geo, fileName);
        collection.addPart(fileName, ".vtp");
    }

    // Write out the collection file
    collection.save();
}

// Export a multibasis mesh
template<class T>
void gsWriteParaview(const gsMultiBasis<T> & mb, const gsMultiPatch<T> & domain,
                     std::string const & fn, unsigned npts)
{
    // GISMO_ASSERT sizes

    gsParaviewCollection collection(fn);

    for (size_t i = 0; i != domain.nPatches(); ++i)
    {
        const std::string fileName = fn + util::to_string(i) + "_mesh";
        writeSingleCompMesh(mb[i], domain.patch(i), fileName, npts);
        collection.addPart(fileName, ".vtp");
    }

    // Write out the collection file
    collection.save();
}

/// Export a Geometry without scalar information
template<class T>
void gsWriteParaview(const gsGeometrySlice<T> & Geo,
                     std::string const & fn,
                     unsigned npts)
{
    const gsMatrix<T> supp = Geo.parameterRange();
    writeSingleGeometry(Geo, supp, fn, npts);
    // Write out a pvd file
    makeCollection(fn, ".vts"); // make also a pvd file
}


/// Export a multipatch Geometry without scalar information
template<class T>
void gsWriteParaview( std::vector<gsGeometry<T> *> const & Geo,
                      std::string const & fn,
                      unsigned npts, bool mesh, bool ctrlNet)
{
    const size_t n = Geo.size();

    gsParaviewCollection collection(fn);
    std::string fnBase;

    for ( size_t i=0; i<n ; i++)
    {
        fnBase = fn + "_" + util::to_string(i);

        if ( Geo.at(i)->domainDim() == 1 )
        {
            writeSingleCurve(*Geo[i], fnBase, npts);
            collection.addPart(fnBase, ".vtp");
        }
        else
        {
            writeSingleGeometry( *Geo[i], fnBase, npts ) ;
            collection.addPart(fnBase, ".vts");
        }

        if ( mesh )
        {
            const std::string fileName = fnBase + "_mesh";
            writeSingleCompMesh(Geo[i]->basis(), *Geo[i], fileName);
            collection.addPart(fileName, ".vtp");
        }

        if ( ctrlNet ) // Output the control net
        {
            const std::string fileName = fnBase + "_cnet";
            writeSingleControlNet(*Geo[i], fileName);
            collection.addPart(fileName, ".vtp");
        }
    }
    collection.save();
}

/// Export i-th Basis function
template<class T>
void gsWriteParaview_basisFnct(int i, gsBasis<T> const& basis, std::string const & fn, unsigned npts)
{
    // basis.support(i) --> returns a (tight) bounding box for the
    // supp. of i-th basis func.
    int d= basis.dim();
    int n= d+1;

    gsMatrix<T> ab = basis.support(i) ;
    gsVector<T> a = ab.col(0);
    gsVector<T> b = ab.col(1);
    gsVector<unsigned> np = uniformSampleCount(a,b, npts );
    gsMatrix<T> pts = gsPointGrid(a,b,np) ;

    gsMatrix<T>  eval_geo = basis.evalSingle ( i, pts ) ;

    if ( 3 - d > 0 )
    {
        np.conservativeResize(3);
        np.bottomRows(3-d).setOnes();
    }

    if ( 2 - d > 0 )
    {
        pts.conservativeResize(2,eval_geo.cols());
        pts.bottomRows(2-d).setZero();
    }

    if ( d > 2 )
    {
//        gsWarn<<"Info: The dimension is to big, projecting into first 2 coordinatess..\n";
        d=2;
        pts.conservativeResize(2,eval_geo.cols());
    }

    if ( 3 - n > 0 )
    {
        eval_geo.conservativeResize(3,eval_geo.cols() );
        eval_geo.bottomRows(3-n).setZero();
    }

    std::string mfn(fn);
    mfn.append(".vts");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        gsWarn<<"gsWriteParaview_basisFnct: Problem opening file \""<<fn<<"\""<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);
    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
    file <<"<StructuredGrid WholeExtent=\"0 "<<np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
    file <<"<Piece Extent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
    // Scalar information
    file <<"<PointData "<< "Scalars"<<"=\"SolutionField\">\n";
    file <<"<DataArray type=\"Float32\" Name=\"SolutionField\" format=\"ascii\" NumberOfComponents=\""<<1<<"\">\n";
    for ( index_t j=0; j<eval_geo.cols(); ++j)
            file<< eval_geo(0,j) <<" ";
    file <<"</DataArray>\n";
    file <<"</PointData>\n";
    //
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\""<<3<<"\">\n";
    for ( index_t j=0; j<eval_geo.cols(); ++j)
    {
        for ( int l=0; l!=d; ++l)
            file<< pts(l,j) <<" ";
        file<< eval_geo(0,j) <<" ";
         for ( index_t l=d; l!=pts.rows(); ++l)
             file<< pts(l,j) <<" ";
    }
    file <<"</DataArray>\n";
    file <<"</Points>\n";
    file <<"</Piece>\n";
    file <<"</StructuredGrid>\n";
    file <<"</VTKFile>\n";
    file.close();
}


/// Export a function
template<class T>
void gsWriteParaview(gsFunction<T> const& func, gsMatrix<T> const& supp, std::string const & fn, unsigned npts)
{
    int d = func.domainDim(); // tested for d==2
    //int n= d+1;

    gsVector<T> a = supp.col(0);
    gsVector<T> b = supp.col(1);
    gsVector<unsigned> np = uniformSampleCount(a,b, npts );
    gsMatrix<T> pts = gsPointGrid(a,b,np);

    gsMatrix<T> ev;
    func.eval_into(pts, ev);

    if ( 3 - d > 0 )
    {
        np.conservativeResize(3);
        np.bottomRows(3-d).setOnes();
    }

    std::string mfn(fn);
    mfn.append(".vts");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        gsWarn<<"gsWriteParaview: Problem opening file \""<<fn<<"\""<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);
    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
    file <<"<StructuredGrid WholeExtent=\"0 "<<np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
    file <<"<Piece Extent=\"0 "<< np(0)-1<<" 0 "<<np(1)-1<<" 0 "<<np(2)-1<<"\">\n";
    // Scalar information
    file <<"<PointData "<< "Scalars"<<"=\"SolutionField\">\n";
    file <<"<DataArray type=\"Float32\" Name=\"SolutionField\" format=\"ascii\" NumberOfComponents=\""<<1<<"\">\n";
    for ( index_t j=0; j<ev.cols(); ++j)
            file<< ev(0,j) <<" ";
    file <<"</DataArray>\n";
    file <<"</PointData>\n";
    //
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\""<<3<<"\">\n";
    for ( index_t j=0; j<ev.cols(); ++j)
    {
        for ( int i=0; i< d; ++i)
            file<< pts(i,j) <<" ";
        file<< ev(0,j) <<" ";
//         for ( index_t i=d; i< pts.rows(); ++i)
//             file<< pts(i,j) <<" ";
    }
    file <<"</DataArray>\n";
    file <<"</Points>\n";
    file <<"</Piece>\n";
    file <<"</StructuredGrid>\n";
    file <<"</VTKFile>\n";
    file.close();
}


/// Export Basis functions
template<class T>
void gsWriteParaview(gsBasis<T> const& basis, std::string const & fn,
                     unsigned npts, bool mesh)
{
    const index_t n = basis.size();
    gsParaviewCollection collection(fn);

    for ( index_t i=0; i< n; i++)
    {
        std::string fileName = fn + util::to_string(i);
        gsWriteParaview_basisFnct<T>(i, basis, fileName, npts ) ;
        collection.addPart(fileName, ".vts");
    }

    if ( mesh )
    {
        std::string fileName = fn + "_mesh";
        writeSingleBasisMesh(basis, fileName);
        //collection.addPart(fileName, ".vtp");
        collection.addPart(fileName, ".vtu");
    }

    collection.save();
}


/// Export Point set to Paraview
template<class T>
void gsWriteParaviewPoints(gsMatrix<T> const& X, gsMatrix<T> const& Y, std::string const & fn)
{
    assert( X.cols() == Y.cols() );
    assert( X.rows() == 1 && Y.rows() == 1 );
    index_t np = X.cols();

    std::string mfn(fn);
    mfn.append(".vtp");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        gsWarn<<"gsWriteParaviewPoints: Problem opening file \""<<fn<<"\""<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);
    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file <<"<PolyData>\n";
    file <<"<Piece NumberOfPoints=\""<<np<<"\" NumberOfVerts=\"1\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    file <<"<PointData>\n";
    file <<"</PointData>\n";
    file <<"<CellData>\n";
    file <<"</CellData>\n";
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\""<<X.minCoeff()<<"\" RangeMax=\""<<X.maxCoeff()<<"\">\n";
    for (index_t i=0; i< np; ++i )
        file << X(0,i) <<" "<<Y(0,i)<<" "<< 0.0 <<"\n";
    file <<"\n</DataArray>\n";
    file <<"</Points>\n";
    file <<"<Verts>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\""<<0<<"\" RangeMax=\""<<np-1<<"\">\n";
    for (index_t i=0; i< np; ++i )
        file << i<<" ";
    file <<"\n</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n"<<np<<"\n";
    file <<"</DataArray>\n";
    file <<"</Verts>\n";
    file <<"<Lines>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<np-1<<"\">\n";
    file <<"</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n";
    file <<"</DataArray>\n";
    file <<"</Lines>\n";
    file <<"<Strips>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<np-1<<"\">\n";
    file <<"</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n";
    file <<"</DataArray>\n";
    file <<"</Strips>\n";
    file <<"<Polys>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<np-1<<"\">\n";
    file <<"</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n";
    file <<"</DataArray>\n";
    file <<"</Polys>\n";
    file <<"</Piece>\n";
    file <<"</PolyData>\n";
    file <<"</VTKFile>\n";
    file.close();

    makeCollection(fn, ".vtp"); // make also a pvd file
}

template<class T>
void gsWriteParaviewPoints(gsMatrix<T> const& X,
                           gsMatrix<T> const& Y,
                           gsMatrix<T> const& Z,
                           std::string const & fn)
{
    GISMO_ASSERT(X.cols() == Y.cols() && X.cols() == Z.cols(),
                 "X, Y and Z must have the same size of columns!");
    GISMO_ASSERT(X.rows() == 1 && Y.rows() == 1 && Z.cols(),
                 "X, Y and Z must be row matrices!");
    index_t np = X.cols();

    std::string mfn(fn);
    mfn.append(".vtp");
    std::ofstream file(mfn.c_str());

    if (!file.is_open())
    {
        gsWarn << "Problem opening " << fn << " Aborting..." << std::endl;
        return;
    }

    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file <<"<PolyData>\n";
    file <<"<Piece NumberOfPoints=\""<<np<<"\" NumberOfVerts=\"1\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    file <<"<PointData>\n"; //empty
    file <<"</PointData>\n";
    file <<"<CellData>\n";
    file <<"</CellData>\n";
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\""<<X.minCoeff()<<"\" RangeMax=\""<<X.maxCoeff()<<"\">\n";

    for (index_t i = 0; i < np; ++i)
    {
        file << X(0, i) << " " << Y(0, i) << " " << Z(0, i) << "\n";
    }

    file <<"\n</DataArray>\n";
    file <<"</Points>\n";
    file <<"<Verts>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\""<<0<<"\" RangeMax=\""<<np-1<<"\">\n";

    for (index_t i=0; i< np; ++i )
    {
        file << i << " ";
    }

    file <<"\n</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n"<<np<<"\n";
    file <<"</DataArray>\n";
    file <<"</Verts>\n";
    file <<"<Lines>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<np-1<<"\">\n";
    file <<"</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n";
    file <<"</DataArray>\n";
    file <<"</Lines>\n";
    file <<"<Strips>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<np-1<<"\">\n";
    file <<"</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n";
    file <<"</DataArray>\n";
    file <<"</Strips>\n";
    file <<"<Polys>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<np-1<<"\">\n";
    file <<"</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n";
    file <<"</DataArray>\n";
    file <<"</Polys>\n";
    file <<"</Piece>\n";
    file <<"</PolyData>\n";
    file <<"</VTKFile>\n";
    file.close();

    makeCollection(fn, ".vtp"); // make also a pvd file
}

template<class T>
void gsWriteParaviewPoints(gsMatrix<T> const& X,
                           gsMatrix<T> const& Y,
                           gsMatrix<T> const& Z,
                           gsMatrix<T> const& V,
                           std::string const & fn)
{
    GISMO_ASSERT(X.cols() == Y.cols() && X.cols() == Z.cols(),
                 "X, Y and Z must have the same size of columns!");
    GISMO_ASSERT(X.rows() == 1 && Y.rows() == 1 && Z.cols(),
                 "X, Y and Z must be row matrices!");
    index_t np = X.cols();

    std::string mfn(fn);
    mfn.append(".vtp");
    std::ofstream file(mfn.c_str());

    if (!file.is_open())
    {
        gsWarn << "Problem opening " << fn << " Aborting..." << std::endl;
        return;
    }

    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file <<"<PolyData>\n";
    file <<"<Piece NumberOfPoints=\""<<np<<"\" NumberOfVerts=\"1\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
    //---------
    file <<"<PointData "<< "Scalars=\"PointInfo\">\n";
    file <<"<DataArray type=\"Float32\" Name=\"PointInfo\" format=\"ascii\" NumberOfComponents=\"1\">\n";
    for ( index_t j=0; j<np; ++j)
        file<< V(0,j) <<" ";
    file <<"</DataArray>\n";
    file <<"</PointData>\n";
    //---------
    file <<"<CellData>\n";
    file <<"</CellData>\n";
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\" RangeMin=\""<<X.minCoeff()<<"\" RangeMax=\""<<X.maxCoeff()<<"\">\n";

    for (index_t i = 0; i < np; ++i)
    {
        file << X(0, i) << " " << Y(0, i) << " " << Z(0, i) << "\n";
    }

    file <<"\n</DataArray>\n";
    file <<"</Points>\n";
    file <<"<Verts>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\""<<0<<"\" RangeMax=\""<<np-1<<"\">\n";

    for (index_t i=0; i< np; ++i )
    {
        file << i << " ";
    }

    file <<"\n</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n"<<np<<"\n";
    file <<"</DataArray>\n";
    file <<"</Verts>\n";
    file <<"<Lines>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<np-1<<"\">\n";
    file <<"</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n";
    file <<"</DataArray>\n";
    file <<"</Lines>\n";
    file <<"<Strips>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<np-1<<"\">\n";
    file <<"</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n";
    file <<"</DataArray>\n";
    file <<"</Strips>\n";
    file <<"<Polys>\n";
    file <<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\""<<np-1<<"\">\n";
    file <<"</DataArray>\n";
    file <<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\" RangeMin=\""<<np<<"\" RangeMax=\""<<np<<"\">\n";
    file <<"</DataArray>\n";
    file <<"</Polys>\n";
    file <<"</Piece>\n";
    file <<"</PolyData>\n";
    file <<"</VTKFile>\n";
    file.close();

    makeCollection(fn, ".vtp"); // make also a pvd file
}

template<class T>
void gsWriteParaviewPoints(gsMatrix<T> const& points, std::string const & fn)
{
    const index_t rows = points.rows();
    switch (rows)
    {
    case 1:
        gsWriteParaviewPoints<T>(points.row(0), gsMatrix<T>::Zero(1, points.cols()), fn);
        break;
    case 2:
        gsWriteParaviewPoints<T>(points.row(0), points.row(1), fn);
        break;
    case 3:
        gsWriteParaviewPoints<T>(points.row(0), points.row(1), points.row(2), fn);
        break;
    case 4:
        gsWriteParaviewPoints<T>(points.row(0), points.row(1), points.row(2), points.row(3), fn);
        break;
    default:
        GISMO_ERROR("Point plotting is implemented just for 2D and 3D (rows== 1, 2 or 3).");
    }
}

// Depicting edge graph of each volume of one gsSolid with a segmenting loop
// INPUTS:
// \param eloop: a vector of ID numbers of vertices, often for representing a segmenting loop
template <class T>
void gsWriteParaview(gsSolid<T> const& sl, std::string const & fn, unsigned numPoints_for_eachCurve, int vol_Num,
                     T edgeThick, gsVector3d<T> const & translate, int color_convex,
                     int color_nonconvex, int color_eloop, std::vector<unsigned> const & eloop)
{
    // options
    int color=color_convex;

    gsSolidHalfFace<T>* face;
    int numOfCurves;
    int numOfPoints = numPoints_for_eachCurve;

    T faceThick = edgeThick;
//    T camera1 = 1;
//    T camera2 = 1;
//    T camera3 = 1;

    std::string mfn(fn);
    mfn.append(".vtp");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        gsWarn<<"gsWriteParaview: Problem opening file \""<<fn<<"\""<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);
    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file <<"<PolyData>\n";


    // collect HEs representing the edge loop
    numOfCurves = eloop.size();
    gsSolidHalfEdge<T>* he;
    std::vector< typename gsSolid<T>::gsSolidHalfEdgeHandle > heSet;
    typename gsSolid<T>::gsSolidHeVertexHandle source,target;
    if (eloop.size()>0){
    for (int iedge=0; iedge!= numOfCurves; iedge++)
    {
        source = sl.vertex[eloop[iedge]];
        target = sl.vertex[eloop[(iedge+1)%numOfCurves]];
        he = source->getHalfEdge(target);
        heSet.push_back(he);
        face = he->face;
    }}


    //gsDebug<<"\n ------------------------------------- number of hafl faces: "<< sl.nHalfFaces();
    for (int iface=0;iface!= sl.nHalfFaces();iface++)
    {
        face = sl.getHalfFaceFromID(iface);
        //gsDebug<<"\n ------------------------------------- vol of face:"<< face->vol->getId()<< " :for face: "<< iface <<"\n";
        //gsDebug << std::flush;
        if (face->vol->getId()==vol_Num)
        {
            numOfCurves=face->nCurvesOfOneLoop(0);
            //gsDebug<<"\n -----------INSIDE-------------------- vol of face:"<< face->vol->getId()<< " :for face: "<< iface <<"\n";

            for (int iedge=0; iedge!= numOfCurves; iedge++)
            {
                he = face->getHalfEdgeFromBoundaryOrder(iedge);
                // search if he is in heSet
                bool isMember(false);
                for (size_t iheSet=0;iheSet<heSet.size();iheSet++)
                {
                    if ( he->isEquiv(heSet.at(iheSet))==true || he->mate->isEquiv(heSet.at(iheSet))==true)
                    {isMember=true;
                        break;}
                }
                gsMatrix<T> curvePoints = face->surf->sampleBoundaryCurve(iedge, numPoints_for_eachCurve);
                if (iedge==0) assert( numOfPoints == curvePoints.cols());
                color=color_convex;
                if (isMember==true) color=color_eloop;
                if (face->getHalfEdgeFromBoundaryOrder(iedge)->is_convex==false){color = color_nonconvex;}
                /// Number of vertices and number of faces
                file <<"<Piece NumberOfPoints=\""<< 2*numOfPoints <<"\" NumberOfVerts=\"0\" NumberOfLines=\""<< 0
                    <<"\" NumberOfStrips=\"0\" NumberOfPolys=\""<< numOfPoints-1 << "\">\n";

                /// Coordinates of vertices
                file <<"<Points>\n";
                file <<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
                // translate the volume towards the *translate* vector
                for (index_t iCol = 0;iCol!=curvePoints.cols();iCol++)
                {
                    file << curvePoints(0,iCol) + translate(0) << " " << curvePoints(1,iCol) + translate(1) << " " << curvePoints(2,iCol) + translate(2) << " \n";
                    // translate the vertex about along the vector (faceThick,0,0)
                    file << curvePoints(0,iCol) + faceThick + translate(0) << " " << curvePoints(1,iCol) + faceThick + translate(1)
                         << " " << curvePoints(2,iCol) +faceThick + translate(2) << " \n";
                };
                file << "\n";
                file <<"</DataArray>\n";
                file <<"</Points>\n";

                /// Scalar field attached to each degenerate face on the "edge"
                file << "<CellData Scalars=\"cell_scalars\">\n";
                file << "<DataArray type=\"Int32\" Name=\"cell_scalars\" format=\"ascii\">\n";
                /// limit: for now, assign all scalars to 0
                for (index_t iCol = 0;iCol!=curvePoints.cols()-1;iCol++)
                {
                    file << color << " ";
                }
                file << "\n";
                file << "</DataArray>\n";
                file << "</CellData>\n";

                /// Which vertices belong to which faces
                file << "<Polys>\n";
                file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
                for (index_t iCol = 0;iCol<=curvePoints.cols()-2;iCol++)
                {
                    //file << iCol << " " << iCol+1 << " "<< iCol+1 << " " << iCol << " ";
                    file << 2*iCol << " " << 2*iCol+1 << " "<< 2*iCol+3 << " " << 2*iCol+2 << " ";
                }
                //file << curvePoints.cols()-1 << " " << 0 << " "<< 0 << " "<< curvePoints.cols()-1 << " ";
                file << "\n";
                file << "</DataArray>\n";
                file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
                unsigned offsets(0);
                for (index_t iCol = 0;iCol!=curvePoints.cols()-1;iCol++)
                {
                    offsets +=4;
                    file << offsets << " ";
                }
                file << "\n";
                file << "</DataArray>\n";
                file << "</Polys>\n";

                file << "</Piece>\n";

                /// Space between edges
                file << "\n";
                file << "\n";
            }
        }
    }

    file <<"</PolyData>\n";
    file <<"</VTKFile>\n";
    file.close();

    makeCollection(fn, ".vtp"); // make also a pvd file
}

template <class T>
void gsWriteParaviewSolid(gsSolid<T> const& sl,
                          std::string const & fn,
                          unsigned numSamples)
{
    const size_t n = sl.numHalfFaces;
    gsParaviewCollection collection(fn);

    // for( typename gsSolid<T>::const_face_iterator it = sl.begin();
    //      it != sl.end(); ++it)

    for ( size_t i=0; i<n ; i++)
    {
        std::string fnBase = fn + util::to_string(i);
        writeSingleTrimSurface(*sl.face[i]->surf, fnBase, numSamples);
        collection.addPart(fnBase, ".vtp");
    }

    // Write out the collection file
    collection.save();
}


/// Visualizing a mesh
template <class T>
void gsWriteParaview(gsMesh<T> const& sl, std::string const & fn, bool pvd)
{
    std::string mfn(fn);
    mfn.append(".vtp");
    std::ofstream file(mfn.c_str());
    if ( ! file.is_open() )
        gsWarn<<"gsWriteParaview: Problem opening file \""<<fn<<"\""<<std::endl;
    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    file <<"<?xml version=\"1.0\"?>\n";
    file <<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file <<"<PolyData>\n";

    /// Number of vertices and number of faces
    file <<"<Piece NumberOfPoints=\""<< sl.numVertices() <<"\" NumberOfVerts=\"0\" NumberOfLines=\""
         << sl.numEdges()<<"\" NumberOfStrips=\"0\" NumberOfPolys=\""<< sl.numFaces() << "\">\n";

    /// Coordinates of vertices
    file <<"<Points>\n";
    file <<"<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (typename std::vector< gsVertex<T>* >::const_iterator it=sl.vertices().begin(); it!=sl.vertices().end(); ++it)
    {
        const gsVertex<T>& vertex = **it;
        file << vertex[0] << " ";
        file << vertex[1] << " ";
        file << vertex[2] << " \n";
    }

    file << "\n";
    file <<"</DataArray>\n";
    file <<"</Points>\n";

    // Scalar field attached to each face
    // file << "<PointData Scalars=\"point_scalars\">\n";
    // file << "<DataArray type=\"Int32\" Name=\"point_scalars\" format=\"ascii\">\n";
    // for (typename std::vector< gsVertex<T>* >::const_iterator it=sl.vertex.begin();
    //      it!=sl.vertex.end(); ++it)
    // {
    //     file << 0 << " ";
    // }
    // file << "\n";
    // file << "</DataArray>\n";
    // file << "</PointData>\n";

    // Write out edges
    file << "<Lines>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (typename std::vector< gsEdge<T> >::const_iterator it=sl.edges().begin();
         it!=sl.edges().end(); ++it)
    {
            file << it->source->getId() << " " << it->target->getId() << "\n";
    }
    file << "</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int count=0;
    for (typename std::vector< gsEdge<T> >::const_iterator it=sl.edges().begin();
         it!=sl.edges().end(); ++it)
    {
        count+=2;
        file << count << " ";
    }
    file << "\n";
    file << "</DataArray>\n";
    file << "</Lines>\n";

    // Scalar field attached to each face (* if edges exists, this has a problem)
    // file << "<CellData Scalars=\"cell_scalars\">\n";
    // file << "<DataArray type=\"Int32\" Name=\"cell_scalars\" format=\"ascii\">\n";
    // for (typename std::vector< gsFace<T>* >::const_iterator it=sl.face.begin();
    //      it!=sl.face.end(); ++it)
    // {
    //     file << 1 << " ";
    // }
    // file << "\n";
    // file << "</DataArray>\n";
    // file << "</CellData>\n";

    /// Which vertices belong to which faces
    file << "<Polys>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (typename std::vector< gsFace<T>* >::const_iterator it=sl.faces().begin();
         it!=sl.faces().end(); ++it)
    {
        for (typename std::vector< gsVertex<T>* >::const_iterator vit= (*it)->vertices.begin();
             vit!=(*it)->vertices.end(); ++vit)
        {
            file << (*vit)->getId() << " ";
        }
        file << "\n";
    }
    file << "</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    count=0;
    for (typename std::vector< gsFace<T>* >::const_iterator it=sl.faces().begin();
         it!=sl.faces().end(); ++it)
    {
        count += (*it)->vertices.size();
        file << count << " ";
    }
    file << "\n";
    file << "</DataArray>\n";
    file << "</Polys>\n";

    file << "</Piece>\n";
    file <<"</PolyData>\n";
    file <<"</VTKFile>\n";
    file.close();

    if( pvd ) // make also a pvd file
        makeCollection(fn, ".vtp");
}

template <typename T>
void gsWriteParaview(const std::vector<gsMesh<T> >& meshes,
                     const std::string& fn)
{
    for (unsigned index = 0; index < meshes.size(); index++)
    {
        std::string file = fn + "Level" + util::to_string(index);
        gsWriteParaview(meshes[index], file, false);
    }
}

template<class T>
void gsWriteParaview(gsPlanarDomain<T> const & pdomain, std::string const & fn, unsigned npts)
{
    std::vector<gsGeometry<T> *> all_curves;
    for(index_t i =0; i<pdomain.numLoops();i++)
        for(index_t j =0; j< pdomain.loop(i).numCurves() ; j++)
            all_curves.push_back( const_cast<gsCurve<T> *>(&pdomain.loop(i).curve(j)) );

    gsWriteParaview( all_curves, fn, npts);
}

template<class T>
void gsWriteParaview(const gsTrimSurface<T> & surf, std::string const & fn,
                     unsigned npts, bool trimCurves)
{
    gsParaviewCollection collection(fn);

    writeSingleTrimSurface(surf, fn, npts);
    collection.addPart(fn, ".vtp");

    if ( trimCurves )
    {
        gsWarn<<"trimCurves: To do.\n";
    }

    // Write out the collection file
    collection.save();
}

template<typename T>
void gsWriteParaview(const gsVolumeBlock<T>& volBlock,
                     std::string const & fn,
                     unsigned npts)
{
    using util::to_string;

    gsParaviewCollection collection(fn);

    // for each face
    for (unsigned idFace = 0; idFace != volBlock.face.size(); idFace++)
    {
        typename gsVolumeBlock<T>::HalfFace* face = volBlock.face[idFace];
        gsPlanarDomain<T>& domain = face->surf->domain();

        // for each curve loop (boundary + holes)
        unsigned numLoops = static_cast<unsigned>(domain.numLoops());
        for (unsigned idLoop = 0; idLoop < numLoops; idLoop++)
        {
            gsCurveLoop<T>& curveLoop = domain.loop(idLoop);

            unsigned clSize = static_cast<unsigned>(curveLoop.size());

            // for each curve in curve loop
            for (unsigned idCurve = 0; idCurve < clSize; idCurve++)
            {
                // file name is fn_curve_Fface_Lloop_Ccurve
                std::string fileName = fn + "_curve_F";
                fileName += to_string(idFace) + "_L" +
                            to_string(idLoop) + "_C" +
                            to_string(idCurve);

                gsWriteParaviewTrimmedCurve(*(face->surf), idLoop, idCurve,
                                            fileName, npts);

                collection.addPart(fileName,".vts");

            } // for each curve
        } // for each curve loop
    } // for each face

    collection.save();
}

template<typename T>
void gsWriteParaviewTrimmedCurve(const gsTrimSurface<T>& surf,
                                 const unsigned idLoop,
                                 const unsigned idCurve,
                                 const std::string fn,
                                 unsigned npts)
{
    // computing parameters and points

    int idL = static_cast<int>(idLoop);
    int idC = static_cast<int>(idCurve);

    gsCurve<T>& curve = surf.getCurve(idL, idC);

    gsMatrix<T> ab = curve.parameterRange() ;
    gsVector<T> a = ab.col(0);
    gsVector<T> b = ab.col(1);

    gsVector<unsigned> np = uniformSampleCount(a, b, npts);
    gsMatrix<T> param = gsPointGrid(a, b, np);

    gsMatrix<T> points;
    surf.evalCurve_into(idLoop, idCurve, param, points);

    np.conservativeResize(3);
    np.bottomRows(3 - 1).setOnes();


    // writing to the file

    std::string myFile(fn);
    myFile.append(".vts");

    std::ofstream file(myFile.c_str());
    if (!file.is_open())
    {
        gsWarn << "Problem opening " << fn << " Aborting..." << std::endl;
        return;
    }

    file << std::fixed; // no exponents
    file << std::setprecision (PLOT_PRECISION);

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">\n";
    file << "<StructuredGrid WholeExtent=\"0 "<< np(0) - 1 <<
            " 0 " << np(1) - 1 << " 0 " << np(2) - 1 << "\">\n";

    file << "<Piece Extent=\"0 " << np(0) - 1 << " 0 " << np(1) - 1 << " 0 "
         << np(2) - 1 << "\">\n";

    file << "<Points>\n";
    file << "<DataArray type=\"Float32\" NumberOfComponents=\"" << points.rows()
         << "\">\n";

    for (index_t j = 0; j < points.cols(); ++j)
    {
        for (index_t i = 0; i < points.rows(); ++i)
        {
            file << points(i, j) << " ";
        }
        file << "\n";
    }

    file << "</DataArray>\n";
    file << "</Points>\n";
    file << "</Piece>\n";
    file << "</StructuredGrid>\n";
    file << "</VTKFile>\n";
    file.close();

}

} // namespace gismo


#undef PLOT_PRECISION
