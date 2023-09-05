/** @file gsIOUtils.h

    @brief Input and output Utilities

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, J. Speh
*/

#pragma once

#include <gsCore/gsDebug.h>
#include <gsHSplines/gsHDomainIterator.h>
#include <gsHSplines/gsHTensorBasis.h>
#include <gsIO/gsBase64.h>

#include <algorithm>
#include <vector>

namespace gismo {

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

}  // namespace gismo
