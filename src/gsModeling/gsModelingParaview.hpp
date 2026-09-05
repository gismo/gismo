/** @file gsModelingParaview.hpp

    @brief Paraview export of gsModeling types (moved from gsIO/gsWriteParaview,
    modularization stream S3 step A8: type-specific visualization lives
    with the type's module, the base IO module stays type-blind).

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#pragma once

#include <gsModeling/gsModelingParaview.h>
#include <gsIO/gsParaviewCollection.h>

#define PLOT_PRECISION 12 // mirrors gsWriteParaview.hpp

// NOTE: deliberately NOT including gsIO/gsWriteParaview.hpp - helper
// templates are declared in gsWriteParaview.h and their real_t
// instantiations are exported from the gsIO module (a second inclusion
// of that .hpp would poison the exported symbol visibility).

namespace gismo
{

template<class T>
void writeSingleTrimSurface(const gsTrimSurface<T> & surf,
                            std::string const & fn,
                            unsigned npts)
{
    typename gsMesh<T>::uPtr msh = surf.toMesh(npts);
    gsWriteParaview( *msh, fn);
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
        std::string fnBase_nopath = gsFileManager::getFilename(fnBase);
        writeSingleTrimSurface(*sl.face[i]->surf, fnBase, numSamples);
        collection.addPart(fnBase_nopath + ".vtp");
    }

    // Write out the collection file
    collection.save();
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
    std::string fn_nopath = gsFileManager::getFilename(fn);
    collection.addPart(fn_nopath + ".vtp");

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
                std::string fileName_nopath = gsFileManager::getFilename(fileName);
                fileName += to_string(idFace) + "_L" +
                            to_string(idLoop) + "_C" +
                            to_string(idCurve);

                gsWriteParaviewTrimmedCurve(*(face->surf), idLoop, idCurve,
                                            fileName, npts);

                collection.addPart(fileName_nopath + ".vts");

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
