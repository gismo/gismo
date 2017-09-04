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
    // Delete vertices
    typename std::vector<VertexHandle>::iterator vIter;
    for(vIter = vertex.begin(); vIter != vertex.end(); vIter++)
    {
        delete *vIter;
    }
    vertex.clear();

    // Delete Faces
    typename std::vector<FaceHandle>::iterator fIter;
    for(fIter = face.begin(); fIter != face.end(); fIter++)
    {
        delete *fIter;
    }
    face.clear();
}

template<class T>
gsMesh<T>::gsMesh(const gsBasis<T> & basis, int n)
: MeshElement(), numVertices(0), numEdges(0), numFaces(0) 
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

        // add points to the mesh
        do
        {
            // get appropriate coordinate of a point
            for (unsigned dim = 0; dim < d; dim++)
            {
                vv(dim) = ( cur(dim) ?  upp(dim) : low(dim) );
            }

            vtx v = addVertex(vv);
            v->data  = vol;
            map[counter++] = v;

        } while (nextCubePoint<gsVector<unsigned> >(cur, zeros, ones));


        // add edges to the mesh (connect points)
        for (std::size_t index = 0; index != neighbour.size(); index++)
        {
            const std::vector<unsigned> & v = neighbour[index];

            for (std::size_t ngh = 0; ngh != v.size(); ngh++)
            {
                // add more vertices (n) for better physical resolution
                addLine( map[index], map[v[ngh]], n );
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
    vertex.push_back(v );
    v->setId(numVertices++);
    return v;
}


template<class T>
typename gsMesh<T>::VertexHandle gsMesh<T>::addVertex(gsVector<T> const & u)
{
    VertexHandle v = this->makeVertex(u);
    vertex.push_back(v);
    v->setId(numVertices++);
    return v;
}


template<class T>
void gsMesh<T>::addEdge(VertexHandle v0, VertexHandle v1)
{
    edge.push_back( Edge(v0,v1) );
    numEdges++;
}


template<class T>
void gsMesh<T>::addEdge(int const & vind0, int const & vind1)
{
    GISMO_ASSERT( vind0 < numVertices, "Invalid vertex index "
                  << vind0 << "(numVertices="<< numVertices <<").");
    GISMO_ASSERT( vind1 < numVertices, "Invalid vertex index "
                  << vind1 << "(numVertices="<< numVertices <<").");

    addEdge(vertex[vind0], vertex[vind1]);
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
    face.push_back( f );
    f->setId(numFaces++);
    return f;
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(VertexHandle const & v0, VertexHandle const & v1, 
                   VertexHandle const & v2)
{
    FaceHandle f = this->makeFace( v0, v1, v2 );
    face.push_back( f );
    f->setId(numFaces++);
    return f;
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(VertexHandle const & v0, VertexHandle const & v1, 
                   VertexHandle const & v2,  VertexHandle const & v3)
{
    FaceHandle f = this->makeFace( v0,v1,v2,v3 );
    face.push_back( f );
    f->setId(numFaces++);
    return f;
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(std::vector<int> const & vert)
{
    std::vector<VertexHandle> pvert; //(vert.size() );
    for ( std::vector<int>::const_iterator it = vert.begin();
          it!= vert.end(); ++it )
        pvert.push_back( vertex[*it] );
    
    FaceHandle f = this->makeFace( pvert );
    face.push_back( f );
    f->setId(numFaces++);
    return f;
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(const int & v0, const int & v1, const int & v2)
{
    FaceHandle f = this->makeFace( vertex[v0],vertex[v1],vertex[v2] );
    face.push_back( f );
    f->setId(numFaces++);
    return f;
}


template<class T>
typename gsMesh<T>::FaceHandle gsMesh<T>::addFace(const int & v0, const int & v1, const int & v2, const int & v3)
{
    FaceHandle f = this->makeFace( vertex[v0],vertex[v1],vertex[v2],vertex[v3] );
    face.push_back( f );
    f->setId(numFaces++);
    return f;
}


template<class T>
std::ostream &gsMesh<T>::print(std::ostream &os) const
{
    os<<"gsMesh with "<<numVertices<<" vertices, "<<numEdges<<
        " edges and "<<numFaces<<" faces.\n";
//             for ( typename std::vector<FaceHandle>::const_iterator
//                       it = face.begin(); it!= face.end(); ++it)
//                 os<<" "<< **it ;
    return os;
}


template <class T>
void gsMesh<T>::cleanStlMesh()
{
    gsWarn<<"Cleaning the stl mesh..."<<"\n";
    
    // This function looks for duplicated vertex coordinates. For each
    // vector, it chooses a unique vertex having that vector, then makes
    // sure that the source and target of every edge is one of the chosen
    // vertices. The old way was more efficient but did not work for
    // non-manifold solids.
    
    // build up the unique map
    std::vector<int> uniquemap;
    for(std::size_t i = 0; i < vertex.size(); i++)
    {
        std::size_t buddy = i;
        for(std::size_t j = 0; j < i; j++)
        {
            if(*(vertex[i]) == *(vertex[j])) // overload compares coords
            {
                buddy = j;
                break;
            }
        }
        uniquemap.push_back(buddy);
    }
    
    for(std::size_t i = 0; i < face.size(); i++)
    {
        for(std::size_t j = 0; j < 3; j++)
        {
            face[i]->vertices[j] = vertex[uniquemap[face[i]->vertices[j]->getId()]];
        }
    }
    
    for(std::size_t i = 0; i < edge.size(); i++)
    {
        edge[i].source = vertex[uniquemap[edge[i].source->getId()]];
        edge[i].target = vertex[uniquemap[edge[i].target->getId()]];
    }

}


template <class T>
void gsMesh<T>::addLine(gsMatrix<T> const & points)
    { 
        const index_t cols = points.cols();
        const bool zzero = ( points.rows()==2 );

        if ( cols < 2 )
            return;
        
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
void gsMesh<T>::addLine(VertexHandle v0, VertexHandle v1, int n)
{
    const gsVector3d<T> & start = v0->coords;
    const T h = (v1->coords - start).norm() / (n+1);
    const gsVector3d<T> step = (v1->coords - start).normalized();
    
    VertexHandle last = v0;
    VertexHandle next;
    for ( int i = 0; i<n; ++i )
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

