/** @file gsIOUtils.h

    @brief Input and output Utilities.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Speh
*/

#pragma once

#include <fstream>

#include <gsUtils/gsMesh/gsMesh.h>
#include <gsUtils/gsCombinatorics.h>

#include <gsHSplines/gsHTensorBasis.h>
#include <gsHSplines/gsHDomainIterator.h>

#include <gsIO/gsFileData.h>
#include <gsIO/gsFileManager.h>

namespace gismo {


/// @brief
/// Writes body part of GoTools (.g2) format to a file.
///
/// \param bspl tensor B-Spline
/// \param out file stream
///
/// \ingroup IO
template <short_t d, typename T>
void gsWriteGoToolsBodySpline(const gsTensorBSpline<d, T>& bspl,
                              std::ofstream& out)
{
    // The format is (similar for d < 3)
    //
    // - dimension of geometry space, whether or not the volume is
    //   rational: 1=rational, 0=non-rational
    // - the number of coefficients in the first parameter direction,
    //   the polynomial order in this direction (i.e. degree+1)
    // - the knot vector in the first parameter direction, multiple
    //   knots are represented by giving the knot value several times
    // - the number of coefficients in the second parameter direction,
    //   the polynomial order in this direction (i.e. degree+1)
    // - the knot vector in the second parameter direction
    // - the number of coefficients in the third parameter direction,
    //   the polynomial order in this direction (i.e. degree+1)
    // - the knot vector in the third parameter direction
    // - the volume coefficients



    out << std::setprecision(15);

    out << bspl.geoDim() << " " << 0 << "\n";

    for (short_t dim = 0; dim < bspl.parDim(); dim++)
    {
        out << bspl.basis().size(dim) << " " <<
               bspl.basis().degree(dim) + 1 << "\n";

        const gsKnotVector<T> kv = bspl.basis().knots(dim);
        typename gsKnotVector<T>::const_iterator iter;
        for (iter = kv.begin(); iter != kv.end(); iter++)
        {
            out << *iter << " ";
        }
        out << "\n";
    }

    const gsMatrix<T>& coefs = bspl.coefs();

    for (int row = 0; row < coefs.rows(); row++)
    {
        out << coefs.row(row) << "\n";

    }

    out << std::endl;
}


/// Writes bspline in GoTools (.g2) format to a file.
///
/// \param bspl tensor B-Spline
/// \param out file stream
template <short_t d, typename T>
void gsWriteGoToolsSpline(const gsTensorBSpline<d, T>& bspl,
                          std::ofstream& out)
{

    // first we write header

    // format ID VERSION, where
    // ID : Class_SplineCurve = 100,
    //      Class_SplineSurface = 200,
    //      Class_SplineVolume = 700
    // VERSION 1 0 0

    // about the headers
    // https://github.com/SINTEF-Geometry/GoTools/blob/master/gotools-core/include/GoTools/geometry/ClassType.h
    // each element in the patch must have its own header

    if (d == 1)
    {
        out << 100;
    }
    else if (d == 2)
    {
        out << 200;
    }
    else if (d == 3)
    {
        out << 700;
    }
    else
    {
        gsWarn << "Dimension is too high: GoTools does not support dimension "
                  "higher than two."
                  "Aborting ...";
        return;
    }

    out << " 1 0 0\n";

    gsWriteGoToolsBodySpline<d, T>(bspl, out);

}


/// Writes geometry in GoTools (.g2) format to a file.
///
/// \param geom geometry
/// \param fileName output file name
template <typename T>
void gsWriteGoTools(const gsGeometry<T>& geom,
                    const std::string& fileName)
{
    std::string fn(fileName);

    // check the extension
    std::string ext = gsFileManager::getExtension(fileName);
    if (ext != "g2")
    {
         fn += ".g2";
    }

    // opening file
    std::ofstream file;
    file.open(fn.c_str());
    if (!file.is_open())
    {
        gsWarn << "Can not open file: " << fileName << "\n"
                  "Aborting ...";
        return;
    }

    gsWriteGoTools(geom, file);

    file.close();

}


/// Writes geometry in GoTools (.g2) format to a file.
///
/// \param geom geometry
/// \param out file stream
template <typename T>
void gsWriteGoTools(const gsGeometry<T>& geom,
                    std::ofstream& out)
{

    // check which object we have and call appropriate function
    const gsTensorBSpline<1, T>* bspline1 =
            dynamic_cast<const gsTensorBSpline<1, T>* >(&geom);

    if (bspline1 != NULL)
    {
        gsWriteGoToolsSpline<1, T>(*bspline1, out);
        return;
    }

    const gsTensorBSpline<2, T>* bspline2 =
            dynamic_cast<const gsTensorBSpline<2, T>* >(&geom);

    if (bspline2 != NULL)
    {
        gsWriteGoToolsSpline<2, T>(*bspline2, out);
        return;
    }

    const gsTensorBSpline<3, T>* bspline3 =
            dynamic_cast<const gsTensorBSpline<3, T>* >(&geom);

    if (bspline3 != NULL)
    {
        gsWriteGoToolsSpline<3, T>(*bspline3, out);
        return;
    }

    gsWarn << "Can not write your geometry, unsupported format.\n"
              "Currently supported formats: BSplines.\n"
              "Aborting ...\n";

    return;
}


/// Writes multi patch in GoTools (.g2) format to a file.
///
/// \param multiPatch multi patch
/// \param fileName output file name
template <typename T>
void gsWriteGoTools(const gsMultiPatch<T>& multiPatch,
                    const std::string& fileName)
{
    std::string fn(fileName);

    // check the extension
    std::string ext = gsFileManager::getExtension(fileName);
    if (ext != "g2")
    {
        fn += ".g2";
    }

    // opening file
    std::ofstream file;
    file.open(fn.c_str());
    if (!file.is_open())
    {
        gsWarn << "Can not open file: " << fileName << "\n"
                  "Aborting ...";
        return;
    }

    typename gsMultiPatch<T>::const_iterator iter;
    for (iter = multiPatch.begin(); iter != multiPatch.end(); ++iter)
    {
        gsWriteGoTools<T>(**iter, file);
    }

    file.close();

}


/// \brief Returns the computational mesh of \a basis.
///
/// \param[in] basis
/// \param mesh
/// \param[in] n number of samples per element side
template<class T>
void makeMesh(const gsBasis<T>& basis, gsMesh<T> & mesh, int n = 0)
{
    const unsigned d = basis.dim();

    typedef typename gsMesh<T>::VertexHandle Vertex;
    typename gsBasis<T>::domainIter domIter = basis.makeDomainIterator();

    // variables for iterating over a cube (element is a cube)
    const gsVector<unsigned> zeros = gsVector<unsigned>::Zero(d);
    const gsVector<unsigned> ones  = gsVector<unsigned>::Ones(d);
    gsVector<unsigned> cur;

    // maps integer representation of a vertex into pointer to the
    // vertex coordinates
    std::vector<Vertex> map(1ULL<<d);

    // neighbour[i] are integer representations of certain neighbours of
    // vertex i (i counts in lexicographics order over all vertices)
    std::vector<std::vector<unsigned> > neighbour(1ULL<<d,
                                                  std::vector<unsigned>() );

    cur.setZero(d);
    int counter = 0;
    do
    {
        // set neighbour
        for (unsigned dim = 0; dim < d; dim++)
        {
            if (cur(dim) == 0)
            {
                const unsigned tmp =  counter | (1<< dim) ;
                neighbour[counter].push_back(tmp);
            }
        }
        counter++;

    } while (nextCubePoint<gsVector<unsigned> >(cur, zeros, ones));

    gsVector<T> vertex(d);

    for (; domIter->good(); domIter->next())
    {
        const gsVector<T>& low = domIter->lowerCorner();
        const gsVector<T>& upp = domIter->upperCorner();
        const T vol = domIter->volume();

        vertex.setZero();
        cur.setZero();
        counter = 0;

        // add points to the mesh
        do
        {
            // get appropriate coordinate of a point
            for (unsigned dim = 0; dim < d; dim++)
            {
                vertex(dim) = ( cur(dim) ?  upp(dim) : low(dim) );
            }

            Vertex v = mesh.addVertex(vertex);
            v->data  = vol;
            map[counter++] = v;

        } while (nextCubePoint<gsVector<unsigned> >(cur, zeros, ones));


        // add edges to the mesh (connect points)
        for (size_t index = 0; index != neighbour.size(); index++)
        {
            const std::vector<unsigned> & v = neighbour[index];

            for (size_t ngh = 0; ngh != v.size(); ngh++)
            {
                // add more vertices (n) for better physical resolution
                mesh.addLine( map[index], map[v[ngh]], n );
                //mesh.addEdge( map[index], map[v[ngh]] );
            }
        }

        // idea: instead of edges add the faces to the mesh
        // mesh->addFace( mesh.vertices().back(),
        //                *(mesh.vertices().end()-3),
        //                *(mesh.vertices.end()-4),
        //                *(mesh.vertices.end()-2)
        //     );

    }
}

namespace internal
{

/// Look at function gismo::makeHierarchicalMesh
template <short_t d, typename T>
void makeHierarchicalMesh(const gsHTensorBasis<d, T>& basis,
                          std::vector<gsMesh<T> >& meshes,
                          int n = 0)
{
    // prepare meshes
    meshes.clear();
    for (unsigned i = 0; i < basis.maxLevel() + 1; i++)
    {
        meshes.push_back(gsMesh<T>());
    }

    // variables for iterating over a cube (element is a cube)
    const gsVector<unsigned> zeros = gsVector<unsigned>::Zero(d);
    const gsVector<unsigned> ones  = gsVector<unsigned>::Ones(d);
    gsVector<unsigned> cur;

    // neighbour[i] are integer representations of certain neighbours of
    // vertex i (i counts in lexicographics order over all vertices)
    std::vector<std::vector<unsigned> > neighbour(1 << d,
                                                  std::vector<unsigned>() );


    cur.setZero(d);
    int counter = 0;
    do
    {
        // set neighbour
        for (unsigned dim = 0; dim < d; dim++)
        {
            if (cur(dim) == 0)
            {
                const unsigned tmp =  counter | (1<< dim) ;
                neighbour[counter].push_back(tmp);
            }
        }
        counter++;

    } while (nextCubePoint<gsVector<unsigned> >(cur, zeros, ones));


    // maps integer representation of a vertex into pointer to the
    // vertex coordinates
    typedef typename gsMesh<T>::VertexHandle Vertex;
    std::vector<Vertex> map(1 << d);

    gsVector<T> vertex(d);

    gsHDomainIterator<T, d> domIter(basis);

    for (; domIter.good(); domIter.next())
    {

        int level = domIter.getLevel();

        const gsVector<T>& low = domIter.lowerCorner();
        const gsVector<T>& upp = domIter.upperCorner();


        vertex.setZero();
        cur.setZero();
        counter = 0;

        // add points to the mesh
        do
        {
            // get appropriate coordinate of a point
            for (unsigned dim = 0; dim < d; dim++)
            {
                vertex(dim) = ( cur(dim) ?  upp(dim) : low(dim) );
            }

            meshes[level].addVertex(vertex);
            map[counter++] = meshes[level].vertices().back();

        } while (nextCubePoint<gsVector<unsigned> >(cur, zeros, ones));


        // add edges to the mesh (connect points)
        for (size_t index = 0; index != neighbour.size(); index++)
        {
            const std::vector<unsigned> & v = neighbour[index];

            for (size_t ngh = 0; ngh != v.size(); ngh++)
            {
                // add more vertices (n) for better physical resolution
                meshes[level].addLine( map[index], map[v[ngh]], n );
            }
        }
    }
}

} // end namespace internal



/// Constructs a series of meshes, each mesh presents on level in hierarchical
/// basis.
///
/// \param basis hierarchocal tensor basis
/// \param meshes we return meshes via reference to std::vector
/// \param n how many vertices is used in presentation of one element
///
/// \result success of the function
template <typename T>
bool makeHierarchicalMesh(const gsBasis<T>& basis,
                          std::vector<gsMesh<T> >& meshes,
                          int n = 0)
{

    const gsHTensorBasis<1, T>* hBasis1 =
            dynamic_cast<const gsHTensorBasis<1, T>* >(&basis);

    if (hBasis1 != NULL)
    {
        internal::makeHierarchicalMesh<1, T>(*hBasis1, meshes, n);
        return true;
    }


    const gsHTensorBasis<2, T>* hBasis2 =
            dynamic_cast<const gsHTensorBasis<2, T>* >(&basis);

    if (hBasis2 != NULL)
    {
        internal::makeHierarchicalMesh<2, T>(*hBasis2, meshes, n);
        return true;
    }


    const gsHTensorBasis<3, T>* hBasis3 =
            dynamic_cast<const gsHTensorBasis<3, T>* >(&basis);

    if (hBasis3 != NULL)
    {
        internal::makeHierarchicalMesh<3, T>(*hBasis3, meshes, n);
        return true;
    }

    return false;
}

} // end namespace gismo
