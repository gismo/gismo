/** @file gsHalfEdgeMesh.hpp

    @brief Provides implementation of the gsHalfEdgeMesh class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl
*/

#pragma once

namespace gismo
{
//struct less_than_ptr
//{
//    bool operator()(gsMesh<>::gsVertexHandle lhs, gsMesh<>::gsVertexHandle rhs)
//    {
//        return ((*lhs) < (*rhs));
//    }
//};
//
//struct equal_ptr
//{
//    bool operator()(gsMesh<>::gsVertexHandle lhs, gsMesh<>::gsVertexHandle rhs)
//    {
//        return ((*lhs) == (*rhs));
//    }
//};

template<class T>
gsHalfEdgeMesh<T>::gsHalfEdgeMesh(const gsMesh<T> &mesh, T precision)
    : gsMesh<T>(mesh), m_precision(precision)
{
    //this->cleanMesh();
    //std::sort(this->m_vertex.begin(), this->m_vertex.end(), less_than_ptr());
    //typename std::vector<gsVertex<T> *, std::allocator<gsVertex<T> *> >::iterator
    //last = std::unique(this->m_vertex.begin(), this->m_vertex.end(), equal_ptr());

    m_halfedges.reserve(3*this->m_face.size());
    for (size_t i = 0; i < this->m_face.size(); i++)
    {
        m_halfedges.push_back(getInternHalfedge(this->m_face[i], 1));
        m_halfedges.push_back(getInternHalfedge(this->m_face[i], 2));
        m_halfedges.push_back(getInternHalfedge(this->m_face[i], 3));
    }

    m_boundary = Boundary(m_halfedges);
    m_n = this->m_vertex.size() - m_boundary.getNumberOfVertices();
    sortVertices();
}

template<class T>
size_t gsHalfEdgeMesh<T>::getNumberOfVertices() const
{
    return this->m_vertex.size();
}

template<class T>
size_t gsHalfEdgeMesh<T>::getNumberOfTriangles() const
{
    return this->m_face.size();
}

template<class T>
size_t gsHalfEdgeMesh<T>::getNumberOfInnerVertices() const
{
    return m_n;
}

template<class T>
size_t gsHalfEdgeMesh<T>::getNumberOfBoundaryVertices() const
{
    return m_boundary.getNumberOfVertices();
}

template<class T>
const typename gsMesh<T>::gsVertexHandle &gsHalfEdgeMesh<T>::getVertex(const size_t vertexIndex) const
{
    GISMO_ASSERT(vertexIndex <= this->m_vertex.size(), "Vertex with index 'vertexIndex'=" << vertexIndex
        << " does not exist. There are only " << this->m_vertex.size() << " vertices.");
    return ((this->m_vertex[m_sorting[vertexIndex - 1]]));
}

/*template<class T>
size_t gsHalfEdgeMesh<T>::getVertexIndex(const gsMesh<T>::gsVertexHandle &vertex) const
{
    return m_inverseSorting[getInternVertexIndex(vertex)];
}*/

template<class T>
size_t gsHalfEdgeMesh<T>::getGlobalVertexIndex(size_t localVertexIndex, size_t triangleIndex) const
{
    return m_inverseSorting[this->m_face[triangleIndex]->vertices[localVertexIndex-1]->getId()];
}

template<class T>
T gsHalfEdgeMesh<T>::getBoundaryLength() const
{
    return m_boundary.getLength();
}

template<class T>
bool rangeCheck(const std::vector<index_t> &corners, const size_t minimum, const size_t maximum)
{
    for (std::vector<index_t>::const_iterator it = corners.begin(); it != corners.end(); it++)
    {
        if ((size_t)*it < minimum || (size_t)*it > maximum)
        { return false; }
    }
    return true;
}

template<class T>
std::vector<T> gsHalfEdgeMesh<T>::getCornerLengths(std::vector<index_t> &corners) const
{
    GISMO_ASSERT(rangeCheck<T>(corners, 1, getNumberOfBoundaryVertices()), "The corners must be <= number of boundary vertices.");

    std::sort(corners.begin(), corners.end());
    size_t s = corners.size();
    std::vector<T> lengths;
    for (size_t i = 0; i < s; i++)
    {
        lengths.push_back(m_boundary.getDistanceBetween(corners[i], corners[(i + 1) % s]));
    }
    return lengths;
}

template<class T>
T gsHalfEdgeMesh<T>::getShortestBoundaryDistanceBetween(size_t i, size_t j) const
{
    return m_boundary.getShortestDistanceBetween(i, j, m_precision);
}

template<class T>
const std::vector<T> gsHalfEdgeMesh<T>::getBoundaryChordLengths() const
{
    return m_boundary.getHalfedgeLengths();
}

template<class T>
T gsHalfEdgeMesh<T>::getHalfedgeLength(size_t originVertexIndex, size_t endVertexIndex) const
{
    GISMO_ASSERT(originVertexIndex <= this->m_vertex.size() && endVertexIndex <= this->m_vertex.size(),
                 "One of the input vertex indices " << originVertexIndex << " or " << endVertexIndex
                 << " does not exist. There are only " << this->m_vertex.size() << " vertices.");

    return gsVector3d<T>(getVertex(originVertexIndex)->x() - getVertex(endVertexIndex)->x(),
                                     getVertex(originVertexIndex)->y() - getVertex(endVertexIndex)->y(),
                                     getVertex(originVertexIndex)->z() - getVertex(endVertexIndex)->z()).norm();
}

template<class T>
short_t gsHalfEdgeMesh<T>::isTriangleVertex(size_t vertexIndex, size_t triangleIndex) const
{
    if (vertexIndex > this->m_vertex.size())
    {
        gsWarn << "gsHalfEdgeMesh::isTriangleVertex: Vertex with vertex index " << vertexIndex
                  << " does not exist. There are only " << this->m_vertex.size() << " vertices.\n";
        return 0;
    }
    if (triangleIndex > getNumberOfTriangles())
    {
        gsWarn << "gsHalfEdgeMesh::isTriangleVertex: The " << triangleIndex
                  << "-th triangle does not exist. There are only " << getNumberOfTriangles() << " triangles.\n";
        return 0;
    }
    if (*(this->m_vertex[m_sorting[vertexIndex - 1]]) == *(this->m_face[triangleIndex]->vertices[0]))
    { return 1; }
    if (*(this->m_vertex[m_sorting[vertexIndex - 1]]) == *(this->m_face[triangleIndex]->vertices[1]))
    { return 2; }
    if (*(this->m_vertex[m_sorting[vertexIndex - 1]]) == *(this->m_face[triangleIndex]->vertices[2]))
    { return 3; }
    return 0; // not contained
}

template<class T>
const std::queue<typename gsHalfEdgeMesh<T>::Halfedge>
gsHalfEdgeMesh<T>::getOppositeHalfedges(const size_t vertexIndex, const bool innerVertex) const
{
    std::queue<Halfedge> oppositeHalfedges;
    if (vertexIndex > this->m_vertex.size())
    {
        gsInfo << "gsHalfEdgeMesh::getOppositeHalfedges: The vertex with index " << vertexIndex
                  << " does not exist. There are only " << this->m_vertex.size() << " vertices.\n";
        return oppositeHalfedges;
    }
    else if (vertexIndex > m_n && innerVertex)
    {
        gsWarn << "gsHalfEdgeMesh::getOppositeHalfedges: Inner vertex with index 'vertexIndex' = " << vertexIndex
                  << "is not an inner vertex. There are only " << m_n << " inner vertices.\n";
    }

    size_t v1, v2, v3;
    for (size_t i = 0; i < getNumberOfTriangles(); i++)
    {
        switch (isTriangleVertex(vertexIndex, i))
        {
            case 1:
                v3 = getGlobalVertexIndex(3, i);
                v2 = getGlobalVertexIndex(2, i);
                oppositeHalfedges.push(Halfedge(v3, v2, getHalfedgeLength(v3, v2)));
                break;
            case 2:
                v1 = getGlobalVertexIndex(1, i);
                v3 = getGlobalVertexIndex(3, i);
                oppositeHalfedges.push(Halfedge(v1, v3, getHalfedgeLength(v1, v3)));
                break;
            case 3:
                v2 = getGlobalVertexIndex(2, i);
                v1 = getGlobalVertexIndex(1, i);
                oppositeHalfedges.push(Halfedge(v2, v1, getHalfedgeLength(v2, v1)));
                break;
            default:
                // vertex is not supposed to show up
                break;
        }
    }
    return oppositeHalfedges;
}

template<class T>
bool gsHalfEdgeMesh<T>::isBoundaryVertex(const size_t internVertexIndex) const
{
    if (internVertexIndex > this->m_vertex.size() - 1)
    {
        gsWarn << "gsHalfEdgeMesh::isBoundaryVertex: Vertex with intern vertex index = " << internVertexIndex
                  << " does not exist. There are only " << this->m_vertex.size() << " vertices.\n";
        return false;
    }
    else
        return m_boundary.isVertexContained(internVertexIndex);
}

template<class T>
const typename gsHalfEdgeMesh<T>::Halfedge
gsHalfEdgeMesh<T>::getInternHalfedge(const typename gsMesh<T>::gsFaceHandle &triangle, size_t numberOfHalfedge) const
{
    size_t index1 = triangle->vertices[0]->getId();
    size_t index2 = triangle->vertices[1]->getId();
    size_t index3 = triangle->vertices[2]->getId();
    if (numberOfHalfedge < 1 || numberOfHalfedge > 3)
    {
        gsWarn << "gsHalfEdgeMesh::getInternHalfedge: The inputted number of the halfedge " << numberOfHalfedge
                  << "  is supposed to be 1,2 or 3. Because input was not expected, first halfedge is returned.\n";
        numberOfHalfedge = 1;
    }
    if (numberOfHalfedge == 1)
    {
        return Halfedge(index2, index1,
                                         gsVector3d<T>(
                                             triangle->vertices[1]->x() - triangle->vertices[0]->x(),
                                             triangle->vertices[1]->y() - triangle->vertices[0]->y(),
                                             triangle->vertices[1]->z() - triangle->vertices[0]->z()).norm());
    }
    if (numberOfHalfedge == 2)
    {
        return Halfedge(index3, index2,
                                         gsVector3d<T>(
                                             triangle->vertices[2]->x() - triangle->vertices[1]->x(),
                                             triangle->vertices[2]->y() - triangle->vertices[1]->y(),
                                             triangle->vertices[2]->z() - triangle->vertices[1]->z()).norm());
    }
    if (numberOfHalfedge == 3)
    {
        return Halfedge(index1, index3,
                                         gsVector3d<T>(
                                             triangle->vertices[0]->x() - triangle->vertices[2]->x(),
                                             triangle->vertices[0]->y() - triangle->vertices[2]->y(),
                                             triangle->vertices[0]->z() - triangle->vertices[2]->z()).norm());
    }
    return Halfedge();
}

template<class T>
void gsHalfEdgeMesh<T>::sortVertices()
{
    // first all interior, then all boundary

    size_t numberOfInnerVerticesFound = 0;
    m_sorting.resize(this->m_vertex.size(), 0);
    m_inverseSorting.resize(this->m_vertex.size(), 0);

    for (size_t i = 0; i != this->m_vertex.size(); ++i)
    {
        if (!isBoundaryVertex(i))
        {
            numberOfInnerVerticesFound++;
            m_sorting[numberOfInnerVerticesFound - 1] = i;
            m_inverseSorting[i] = numberOfInnerVerticesFound;
        }
    }

    std::list<size_t> boundaryVertices = m_boundary.getVertexIndices();
    for (size_t i = 0; i < getNumberOfBoundaryVertices(); i++)
    {
        m_sorting[m_n + i] = boundaryVertices.front();
        m_inverseSorting[boundaryVertices.front()] = m_n + i + 1;
        boundaryVertices.pop_front();
    }
}

// nested class Boundary
/// @cond
template<class T>
gsHalfEdgeMesh<T>::Boundary::Boundary(const std::vector<typename gsHalfEdgeMesh<T>::Halfedge> &halfedges)
{
    std::list<Halfedge> unsortedNonTwinHalfedges = findNonTwinHalfedges(halfedges);
    m_boundary.appendNextHalfedge(unsortedNonTwinHalfedges.front());
    unsortedNonTwinHalfedges.pop_front();
    std::queue<Halfedge> nonFittingHalfedges;
    while (!unsortedNonTwinHalfedges.empty())
    {
        if (m_boundary.isAppendableAsNext(unsortedNonTwinHalfedges.front()))
        {
            m_boundary.appendNextHalfedge(unsortedNonTwinHalfedges.front());
            unsortedNonTwinHalfedges.pop_front();
            while (!nonFittingHalfedges.empty())
            {
                unsortedNonTwinHalfedges.push_back(nonFittingHalfedges.front());
                nonFittingHalfedges.pop();
            }
        }
        else if (m_boundary.isAppendableAsPrev(unsortedNonTwinHalfedges.front()))
        {
            m_boundary.appendPrevHalfedge(unsortedNonTwinHalfedges.front());
            unsortedNonTwinHalfedges.pop_front();
            while (!nonFittingHalfedges.empty())
            {
                unsortedNonTwinHalfedges.push_back(nonFittingHalfedges.front());
                nonFittingHalfedges.pop();
            }
        }
        else
        {
            nonFittingHalfedges.push(unsortedNonTwinHalfedges.front());
            unsortedNonTwinHalfedges.pop_front();
        }
    }
    if (!m_boundary.isClosed())
        gsWarn << "gsHalfEdgeMesh::Boundary::Boundary: Boundary is not closed although it should be. End points are:\n"
                  << m_boundary.getFirstHalfedge().getOrigin() << "\n and "
                  << m_boundary.getLastHalfedge().getEnd() << "\n";
}
/// @endcond
// private

template<class T>
const std::list<typename gsHalfEdgeMesh<T>::Halfedge> gsHalfEdgeMesh<T>::Boundary::findNonTwinHalfedges(const std::vector<typename gsHalfEdgeMesh<T>::Halfedge> &allHalfedges)
{
    std::queue<Halfedge> queue0;
    std::queue<Halfedge> queue1;
    bool actualQueue = 0;
    std::list<Halfedge> nonTwinHalfedges;
    for (size_t i = 0; i < allHalfedges.size(); ++i)
    {
        queue0.push(allHalfedges[i]);
    }
    while (!(queue0.empty() && queue1.empty()))
    {
        if (actualQueue == 0)
        {
            nonTwinHalfedges.push_back(queue0.front());
            queue0.pop();
            while (!queue0.empty())
            {
                if (nonTwinHalfedges.back().isTwin(queue0.front()))
                {
                    queue0.pop();
                    nonTwinHalfedges.pop_back();
                    while (!queue0.empty())
                    {
                        queue1.push(queue0.front());
                        queue0.pop();
                    }
                }
                else
                {
                    queue1.push(queue0.front());
                    queue0.pop();
                }
            }
            actualQueue = 1;
        }
        else if (actualQueue == 1)
        {
            nonTwinHalfedges.push_back(queue1.front());
            queue1.pop();
            while (!queue1.empty())
            {
                if (nonTwinHalfedges.back().isTwin(queue1.front()))
                {
                    queue1.pop();
                    nonTwinHalfedges.pop_back();
                    while (!queue1.empty())
                    {
                        queue0.push(queue1.front());
                        queue1.pop();
                    }
                }
                else
                {
                    queue0.push(queue1.front());
                    queue1.pop();
                }
            }
            actualQueue = 0;
        }
    }
    return nonTwinHalfedges;
}

// nested class Chain

template<class T>
bool gsHalfEdgeMesh<T>::Chain::isClosed() const
{
    if (m_chainedHalfedges.empty())
    {
        gsWarn << "gsHalfEdgeMesh<T>::Chain::isClosed: The chain does not store any halfedges yet.\n";
        return true;
    }
    return (m_chainedHalfedges.front().getOrigin() == m_chainedHalfedges.back().getEnd());
}

template<class T>
size_t gsHalfEdgeMesh<T>::Chain::getNumberOfVertices() const
{
    if (this->isClosed())
        return m_chainedHalfedges.size();
    else
        return m_chainedHalfedges.size() + 1;
}

template<class T>
T gsHalfEdgeMesh<T>::Chain::getLength() const
{
    T length = 0;
    for (typename std::list<Halfedge>::const_iterator it = m_chainedHalfedges.begin();
         it != m_chainedHalfedges.end(); ++it)
    {
        length += it->getLength();
    }
    return length;
}

template<class T>
std::vector<T> gsHalfEdgeMesh<T>::Chain::getHalfedgeLengths() const
{
    if (this->isEmpty())
    {
        gsWarn << "gsHalfEdgeMesh::Chain::getHalfedgeLengths: The chain does not store any halfedges yet.\n";
    }
    std::vector<T> lengths;
    for (typename std::list<Halfedge>::const_iterator it = m_chainedHalfedges.begin();
         it != m_chainedHalfedges.end(); ++it)
    {
        lengths.push_back(it->getLength());
    }
    return lengths;
}

template<class T>
const typename gsHalfEdgeMesh<T>::Halfedge& gsHalfEdgeMesh<T>::Chain::getFirstHalfedge() const
{
    if (this->isEmpty())
    {
        gsInfo << "gsHalfEdgeMesh::Chain::getFirstHalfedge: The chain does not store any halfedges yet.\n";
    }
    return m_chainedHalfedges.front();
}

template<class T>
const typename gsHalfEdgeMesh<T>::Halfedge& gsHalfEdgeMesh<T>::Chain::getLastHalfedge() const
{
    if (this->isEmpty())
    {
        gsInfo << "gsHalfEdgeMesh::Chain::getLastHalfedge: The chain does not store any halfedges yet.\n";
    }
    return m_chainedHalfedges.back();
}

template<class T>
const std::list<size_t> gsHalfEdgeMesh<T>::Chain::getVertexIndices() const
{
    if (this->isEmpty())
    {
        gsWarn << "gsHalfEdgeMesh::Chain::getVertexIndices: The chain does not store any halfedges yet.\n";
    }
    std::list<size_t> vertexIndices;
    for (typename std::list<Halfedge>::const_iterator it = m_chainedHalfedges.begin();
         it != m_chainedHalfedges.end(); ++it)
    {
        vertexIndices.push_back(it->getOrigin());
    }
    if (m_chainedHalfedges.size() == 1 || !(this->isClosed()))
    {
        vertexIndices.push_back(m_chainedHalfedges.back().getEnd());
    }
    return vertexIndices;
}

template<class T>
T gsHalfEdgeMesh<T>::Chain::getShortestDistanceBetween(size_t i, size_t j, T precision) const
{
    if (this->isEmpty())
    {
        gsWarn << "gsHalfEdgeMesh::Chain::getShortestDistanceBetween: The chain does not store any halfedges yet.\n";
        return 0;
    }
    if (i > this->getNumberOfVertices() || j > this->getNumberOfVertices() || i < 1 || j < 1)
    {
        gsInfo << "gsHalfEdgeMesh::Chain::getShortestDistanceBetween: FirstIndex: " << i << " and second index: "
                  << j
                  << " must be positiv integers smaller than the number of points of the chain, which is "
                  << this->getNumberOfVertices() << ".\n";
        return 0;
    }
    if (i > j) std::swap(i, j);
    T distance = 0;
    std::vector<T> l = this->getHalfedgeLengths();
    for (size_t z = i - 1; z < j - 1; z++)
    {
        distance += l[z];
    }
    if (this->isClosed() && (this->getLength() - distance < distance - precision))
    {
        distance = this->getLength() - distance;
    }
    return distance;
}

template<class T>
T gsHalfEdgeMesh<T>::Chain::getDistanceBetween(size_t i, size_t j) const
{
    if (this->isEmpty())
    {
        gsWarn << "gsHalfEdgeMesh::Chain::getDistanceBetween: The chain does not store any halfedges yet.\n";
        return 0;
    }
    if (i > this->getNumberOfVertices() || j > this->getNumberOfVertices() || i < 1 || j < 1)
    {
        gsInfo << "gsHalfEdgeMesh::Chain::getDistanceBetween: FirstIndex: " << i << " and second index: "
                  << j
                  << " must be positiv integers smaller than the number of points of the chain, which is "
                  << this->getNumberOfVertices() << ".\n";
        return 0;
    }
    bool ordered = (i < j);
    if (!ordered) std::swap(i, j);
    T distance = 0;
    std::vector<T> l = this->getHalfedgeLengths();
    for (size_t z = i - 1; z < j - 1; z++)
    {
        distance += l[z];
    }
    if (this->isClosed() && !ordered)
    {
        distance = this->getLength() - distance;
    }
    else if (!ordered)
        gsWarn << "gsHalfEdgeMesh::Chain::getDistanceBetween: The chain is supposed to be closed in case the input is not ordered.\n";
    return distance;
}

template<class T>
bool gsHalfEdgeMesh<T>::Chain::isVertexContained(const size_t &vertexIndex) const
{
    if (this->isEmpty())
    {
        gsWarn << "gsHalfEdgeMesh::Chain::isVertexContained: The chain does not store any halfedges yet.\n";
        return false;
    }
    for (typename std::list<Halfedge>::const_iterator it = m_chainedHalfedges.begin();
         it != m_chainedHalfedges.end(); ++it)
    {
        if (it->getOrigin() == vertexIndex)
            return true;
    }
    if (!(this->isClosed()) && this->getLastHalfedge().getEnd() == vertexIndex) //here ! was added
        return true;
    return false;
}

template<class T>
bool gsHalfEdgeMesh<T>::Chain::isAppendableAsPrev(const typename gsHalfEdgeMesh<T>::Halfedge &previousHalfedge) const
{
    if (m_chainedHalfedges.empty())
        return true;
    else
        return m_chainedHalfedges.front().isNext(previousHalfedge);
}

template<class T>
bool gsHalfEdgeMesh<T>::Chain::isAppendableAsNext(const typename gsHalfEdgeMesh<T>::Halfedge &nextHalfedge) const
{
    if (m_chainedHalfedges.empty())
        return true;
    else
        return m_chainedHalfedges.back().isPrev(nextHalfedge);
}

template<class T>
void gsHalfEdgeMesh<T>::Chain::appendPrevHalfedge(const typename gsHalfEdgeMesh<T>::Halfedge &prevHalfedge)
{
    if (!this->isAppendableAsPrev(prevHalfedge))
    {
        gsWarn << "gsHalfEdgeMesh::Chain::appendPrevHalfedge: This halfedge is not appendable at the beginning.\n";
        gsInfo << "The first halfedge of the chain has origin " << this->getFirstHalfedge().getOrigin()
                  << " and prevHalfedge has the end " << prevHalfedge.getEnd() << ".\n";
    }
    else
    {
        m_chainedHalfedges.push_front(prevHalfedge);
    }
}

template<class T>
void gsHalfEdgeMesh<T>::Chain::appendNextHalfedge(const typename gsHalfEdgeMesh<T>::Halfedge &nextHalfedge)
{
    if (!isAppendableAsNext(nextHalfedge))
    {
        gsWarn << "gsHalfEdgeMesh::Chain::appendNextHalfedge: This halfedge is not appendable at the end.\n";
        gsInfo << "The last halfedge of the chain has end " << this->getLastHalfedge().getEnd()
                  << " and nextHalfedge has the origin " << nextHalfedge.getOrigin() << ".\n";
    }
    else
    {
        m_chainedHalfedges.push_back(nextHalfedge);
    }
}

} // namespace gismo
