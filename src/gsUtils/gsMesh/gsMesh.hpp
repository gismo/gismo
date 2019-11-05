/** @file gsMesh.hpp

    @brief Provides implementation of the Mesh class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, D. Mayer
*/

#pragma once

#include <gsCore/gsBasis.h>
#include <gsUtils/gsCombinatorics.h>
#include <gsCore/gsDomainIterator.h>

namespace gismo
{

template<class T>
gsMesh<T>::~gsMesh()
{
    //gsInfo << "delete gsMesh\n";
    freeAll(m_vertex);
    freeAll(m_face);
}

template<class T>
gsMesh<T>::gsMesh(const gsBasis<T> & basis, int midPts)
: MeshElement()
{
    const unsigned d = basis.dim();

    typedef typename gsMesh<T>::VertexHandle vtx;
    typename gsBasis<T>::domainIter domIter = basis.makeDomainIterator();

    // variables for iterating over a cube (element is a cube)
    const gsVector<unsigned> zeros = gsVector<unsigned>::Zero(d);
    const gsVector<unsigned> ones  = gsVector<unsigned>::Ones(d);
    gsVector<unsigned> cur;

    // maps integer representation of a vertex into pointer to the
    // vertex coordinates
    std::vector<vtx> map(1ULL<<d);

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

    gsVector<T> vv(d);

    for (; domIter->good(); domIter->next())
    {
        const gsVector<T>& low = domIter->lowerCorner();
        const gsVector<T>& upp = domIter->upperCorner();
        const T vol = domIter->volume();

        vv.setZero();
        cur.setZero();
        counter = 0;

        // Add points to the mesh.
        do
        {
            // Get the appropriate coordinate of a point.
            for (unsigned dim = 0; dim < d; dim++)
            {
                vv(dim) = ( cur(dim) ?  upp(dim) : low(dim) );
            }

            vtx v = addVertex(vv);
            v->data  = vol;
            map[counter++] = v;

        } while (nextCubePoint<gsVector<unsigned> >(cur, zeros, ones));


        // Add edges to the mesh (connect points).
        for (size_t index = 0; index != neighbour.size(); index++)
        {
            const std::vector<unsigned> & v = neighbour[index];

            for (size_t ngh = 0; ngh != v.size(); ngh++)
            {
                // Add more vertices for better physical resolution.
                addLine( map[index], map[v[ngh]], midPts );
                //addEdge( map[index], map[v[ngh]] );
            }
        }

        // idea: instead of edges add the faces to the mesh
        // addFace( mesh.vertex.back(),
        //                *(vertex.end()-3),
        //                *(vertex.end()-4),
        //                *(vertex.end()-2)
        //     );
    }
}


template<class T>
typename gsMesh<T>::VertexHandle gsMesh<T>::addVertex(scalar_t const& x, scalar_t const& y, scalar_t const& z)
{
    VertexHandle v = this->makeVertex(x,y,z);
    v->setId(m_vertex.size());
    m_vertex.push_back(v );
    return v;
}


template<class T>
typename gsMesh<T>::VertexHandle gsMesh<T>::addVertex(gsVector<T> const & u)
{
    VertexHandle v = this->makeVertex(u);
    v->setId(m_vertex.size());
    m_vertex.push_back(v);
    return v;
}


template<class T>
void gsMesh<T>::addEdge(VertexHandle v0, VertexHandle v1)
{
    m_edge.push_back( Edge(v0,v1) );
}


template<class T>
void gsMesh<T>::addEdge(int const vind0, int const vind1)
{
    GISMO_ASSERT( (size_t)vind0 < numVertices(), "Invalid vertex index "
                  << vind0 << "(numVertices="<< numVertices() <<").");
    GISMO_ASSERT( (size_t)vind1 < numVertices(), "Invalid vertex index "
                  << vind1 << "(numVertices="<< numVertices() <<").");

    addEdge(m_vertex[vind0], m_vertex[vind1]);
}


template<class T>
void gsMesh<T>::addEdge(gsVector<T> const & u0,
             gsVector<T> const & u1 )
{
    addEdge( addVertex(u0), addVertex(u1) );
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(std::vector<VertexHandle> const & vert)
{
    FaceHandle f = this->makeFace( vert );
    f->setId(m_face.size());
    m_face.push_back(f);
    return f;
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(VertexHandle const & v0, VertexHandle const & v1,
                   VertexHandle const & v2)
{
    FaceHandle f = this->makeFace( v0, v1, v2 );
    f->setId(m_face.size());
    m_face.push_back(f);
    return f;
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(VertexHandle const & v0, VertexHandle const & v1,
                   VertexHandle const & v2,  VertexHandle const & v3)
{
    FaceHandle f = this->makeFace( v0,v1,v2,v3 );
    f->setId(m_face.size());
    m_face.push_back(f);
    return f;
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(std::vector<int> const & vert)
{
    std::vector<VertexHandle> pvert; //(vert.size() );
    pvert.reserve(vert.size());
    for ( std::vector<int>::const_iterator it = vert.begin();
          it!= vert.end(); ++it )
        pvert.push_back( m_vertex[*it] );

    FaceHandle f = this->makeFace( pvert );
    f->setId(m_face.size());
    m_face.push_back(f);
    return f;
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(const int v0, const int v1, const int v2)
{
    FaceHandle f = this->makeFace( m_vertex[v0],m_vertex[v1],m_vertex[v2] );
    f->setId(m_face.size());
    m_face.push_back(f);
    return f;
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(const int v0, const int v1, const int v2, const int v3)
{
    FaceHandle f = this->makeFace( m_vertex[v0],m_vertex[v1],m_vertex[v2],m_vertex[v3] );
    f->setId(m_face.size());
    m_face.push_back(f);
    return f;
}


template<class T>
std::ostream &gsMesh<T>::print(std::ostream &os) const
{
    os<<"gsMesh with "<<numVertices()<<" vertices, "<<numEdges()<<
        " edges and "<<numFaces()<<" faces.\n";
//             for ( typename std::vector<FaceHandle>::const_iterator
//                       it = face.begin(); it!= face.end(); ++it)
//                 os<<" "<< **it ;
    return os;
}


template <class T>
gsMesh<T>& gsMesh<T>::cleanMesh()
{
    gsDebug << "Cleaning the gsMesh\n";

    // This function looks for duplicated vertex coordinates. For each
    // vector, it chooses a unique vertex having that vector, then makes
    // sure that the source and target of every edge is one of the chosen
    // vertices. The old way was more efficient but did not work for
    // non-manifold solids.

    /*gsDebug << "std::vector<> vertex before cleanMesh\n";
    for (size_t i = 0; i < m_vertex.size(); ++i)
    {
        gsDebug << i << ": " << m_vertex[i] << " id: " << m_vertex[i]->getId() << " " << *m_vertex[i];
    }
    gsDebug << "----------------------------------------\n";*/

    // build up the unique map
    std::vector<size_t> uniquemap;
    uniquemap.reserve(m_vertex.size());
    for(size_t i = 0; i < m_vertex.size(); i++)
    {
        size_t buddy = i;
        for(size_t j = 0; j < i; j++)
        {
            if(*(m_vertex[i]) == *(m_vertex[j])) // overload compares coords
            {
                buddy = j;
                break;
            }
        }
        uniquemap.push_back(buddy);
    }

    for(size_t i = 0; i < m_face.size(); i++)
    {
        for (size_t j = 0; j < 3; j++)
        {
            m_face[i]->vertices[j] = m_vertex[uniquemap[m_face[i]->vertices[j]->getId()]];
        }
    }

    for(size_t i = 0; i < m_edge.size(); i++)
    {
        m_edge[i].source = m_vertex[uniquemap[m_edge[i].source->getId()]];
        m_edge[i].target = m_vertex[uniquemap[m_edge[i].target->getId()]];
    }

    std::set<size_t> uniqueset(uniquemap.begin(), uniquemap.end());    // O(n*log(n))
    std::vector<VertexHandle> uvertex;
    uvertex.reserve(uniqueset.size());
    for(size_t i = 0; i < uniquemap.size(); i++) {     // O(n)
        if(uniqueset.find(i) != uniqueset.end())    // O(log(m)), n >> m
        {
            // re-number vertices id by new sequence - should we not do?
            m_vertex[i]->setId(uvertex.size());
            uvertex.push_back(m_vertex[i]);
        }
        else
        {
            delete m_vertex[i];
            m_vertex[i] = nullptr;
        }
    }   // O(n*log(n)+O(n)*O(log(m)) ==> O(n*log(n))
    m_vertex.swap(uvertex);

    return *this;
}

template<class T>
gsMesh<T>& gsMesh<T>::reserve(size_t vertex, size_t face, size_t edge)
{
    m_vertex.reserve(vertex);
    m_face.reserve(face);
    m_edge.reserve(edge);
    return *this;
}

template <class T>
void gsMesh<T>::addLine(gsMatrix<T> const & points)
    {
        const index_t cols = points.cols();
        if ( cols < 2 ) return;

        const bool zzero = ( points.rows()==2 );

        VertexHandle v1,
            v0 = addVertex( points(0,0), points(1,0), zzero ? 0 : points(2,0) );

        for ( index_t i = 1; i<cols; ++i)
        {
            v1 = addVertex( points(0, i), points(1, i), zzero ? 0 : points(2,0) );
            addEdge(v0 , v1);
            v0 = v1;
        }
    }

template <class T>
void gsMesh<T>::addLine(VertexHandle v0, VertexHandle v1, int midPts)
{
    const gsVector3d<T> & start = *dynamic_cast<gsVector3d<T>* >(v0);
    const T h = (*v1 - start).norm() / (midPts + 1);
    const gsVector3d<T> step = (*v1 - start).normalized();

    VertexHandle last = v0;
    VertexHandle next;
    for ( int i = 0; i<midPts; ++i )
    {
        next = addVertex(start + i*h*step);
        addEdge(last, next);
        last = next;
    }
    addEdge(last, v1);
}


//template <class T>
//gsSolid<T> *
//gsMesh<T>::meshToSolid(int kvOuterPoints, int kvAdditionalInnerPoints, bool addInner,
//                       bool plot,std::string name,int meshPoints, T wE, T wI,int closeBoundary,
//                       std::vector<std::vector<VertexHandle> > iPoints,
//                       std::vector<std::vector<VertexHandle> > oPoints,
//                       std::vector< std::vector<std::vector<VertexHandle> > > innerBdrys,
//                       std::vector< std::vector<Vertex>  > innerBdrysMassP,
//                       std::vector<std::vector<bool> > oPointsConvexFlag)
//    {
//        gsSolid<T> * s1 = new gsSolid<T>();
//        this->toSolid(*s1,iPoints,oPoints,innerBdrys,innerBdrysMassP,oPointsConvexFlag,kvOuterPoints,kvAdditionalInnerPoints,plot,name,meshPoints,addInner,wE,wI,closeBoundary);
//        return s1;
//    }

};// namespace gismo
