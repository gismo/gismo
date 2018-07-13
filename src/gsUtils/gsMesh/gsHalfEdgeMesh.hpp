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
//**********************************************
//************ class gsHalfEdgeMesh ************
//**********************************************
struct less_than_ptr
{
    bool operator()(gsMesh<>::gsVertexHandle lhs, gsMesh<>::gsVertexHandle rhs)
    {
        return ((*lhs) < (*rhs));
    }
};

struct equal_ptr
{
    bool operator()(gsMesh<>::gsVertexHandle lhs, gsMesh<>::gsVertexHandle rhs)
    {
        return ((*lhs) == (*rhs));
    }
};

//********************************************************************************

template<class T>
gsHalfEdgeMesh<T>::gsHalfEdgeMesh(const gsMesh<> &mesh)
    : gsMesh<>(mesh)
{
    std::sort(this->vertex.begin(), this->vertex.end(), less_than_ptr());
    std::vector<gsVertex<double> *, std::allocator<gsVertex<double> *> >::iterator
    last = std::unique(this->vertex.begin(), this->vertex.end(), equal_ptr());
    this->vertex.erase(last, this->vertex.end());
    for (std::size_t i = 0; i < this->face.size(); i++)
    {
        m_halfedges.push_back(getInternHalfedge(this->face[i], 1));
        m_halfedges.push_back(getInternHalfedge(this->face[i], 2));
        m_halfedges.push_back(getInternHalfedge(this->face[i], 3));
    }
    m_boundary = Boundary(m_halfedges);
    m_n = this->vertex.size() - m_boundary.getNumberOfVertices();
    sortVertices();
}

template<class T>
std::size_t gsHalfEdgeMesh<T>::getNumberOfVertices() const
{
    return this->vertex.size();
}

template<class T>
std::size_t gsHalfEdgeMesh<T>::getNumberOfTriangles() const
{
    return this->face.size();
}

template<class T>
std::size_t gsHalfEdgeMesh<T>::getNumberOfInnerVertices() const
{
    return m_n;
}

template<class T>
std::size_t gsHalfEdgeMesh<T>::getNumberOfBoundaryVertices() const
{
    return m_boundary.getNumberOfVertices();
}

template<class T>
const gsMesh<>::gsVertexHandle &gsHalfEdgeMesh<T>::getVertex(const std::size_t vertexIndex) const
{
    if (vertexIndex > this->vertex.size())
    {
        std::cerr << "Error: [" << __PRETTY_FUNCTION__ << "] Vertex with index 'vertexIndex'=" << vertexIndex
                  << " does not exist. There are only " << this->vertex.size() << " vertices." << std::endl;
    }
    return ((this->vertex[m_sorting[vertexIndex - 1]]));
}

template<class T>
std::size_t gsHalfEdgeMesh<T>::getVertexIndex(const gsMesh<>::gsVertexHandle &vertex) const
{
    return m_inverseSorting[getInternVertexIndex(vertex)];
}

template<class T>
std::size_t gsHalfEdgeMesh<T>::getGlobalVertexIndex(std::size_t localVertexIndex, std::size_t triangleIndex) const
{
    if ((localVertexIndex != 1 && localVertexIndex != 2 && localVertexIndex != 3)
        || triangleIndex > getNumberOfTriangles() - 1)
        std::cerr << "Error: [" << __PRETTY_FUNCTION__ << "] The 'localVertexIndex'=" << localVertexIndex
                  << " should be 1,2 or 3 and the triangle with 'triangleIndex'=" << triangleIndex << " does not exist"
                  << std::endl;
    if (localVertexIndex == 1)
        return getVertexIndex((this->face[triangleIndex]->vertices[0]));
    if (localVertexIndex == 2)
        return getVertexIndex((this->face[triangleIndex]->vertices[1]));
    return getVertexIndex((this->face[triangleIndex]->vertices[2]));
}

template<class T>
double gsHalfEdgeMesh<T>::getBoundaryLength() const
{
    return m_boundary.getLength();
}

template<class T>
bool rangeCheck(const std::vector<int> &corners, const std::size_t minimum, const std::size_t maximum)
{
    for (std::vector<int>::const_iterator it = corners.begin(); it != corners.end(); it++)
    {
        if (*it < minimum || *it > maximum)
        { return false; }
    }
    return true;
}

template<class T>
std::vector<double> gsHalfEdgeMesh<T>::getCornerLengths(std::vector<int> &corners) const
{
    std::size_t B = getNumberOfBoundaryVertices();
    if (!rangeCheck<T>(corners, 1, B))
    {
        std::cerr << "Error: [" << __PRETTY_FUNCTION__ << "] The corners must be <= number of boundary vertices."
                  << std::endl;
        std::cerr << "One of these is >= " << getNumberOfBoundaryVertices() << std::endl;
        for (std::vector<int>::const_iterator it = corners.begin(); it != corners.end(); it++)
        {
            std::cout << *it << std::endl;
        }
    }
    std::sort(corners.begin(), corners.end());
    std::size_t s = corners.size();
    std::vector<double> lengths;
    for (std::size_t i = 0; i < s; i++)
    {
        lengths.push_back(m_boundary.getDistanceBetween(corners[i], corners[(i + 1) % s]));
    }
    return lengths;
}

template<class T>
double gsHalfEdgeMesh<T>::getShortestBoundaryDistanceBetween(std::size_t i, std::size_t j) const
{
    return m_boundary.getShortestDistanceBetween(i, j);
}

template<class T>
const std::vector<double> gsHalfEdgeMesh<T>::getBoundaryChordLengths() const
{
    return m_boundary.getHalfedgeLengths();
}

template<class T>
double gsHalfEdgeMesh<T>::getHalfedgeLength(std::size_t originVertexIndex, std::size_t endVertexIndex) const
{
    if (originVertexIndex > this->vertex.size() || endVertexIndex > this->vertex.size())
    {
        std::cerr << "Error: [" << __PRETTY_FUNCTION__ << "] One of the input vertex indices " << originVertexIndex
                  << " or " << endVertexIndex << " does not exist. There are only " << this->vertex.size()
                  << " vertices." << std::endl;
    }
    return gsVector3d<real_t>(getVertex(originVertexIndex)->x() - getVertex(endVertexIndex)->x(),
                                     getVertex(originVertexIndex)->y() - getVertex(endVertexIndex)->y(),
                                     getVertex(originVertexIndex)->z() - getVertex(endVertexIndex)->z()).norm();
}

template<class T>
triangleVertexIndex gsHalfEdgeMesh<T>::isTriangleVertex(std::size_t vertexIndex, std::size_t triangleIndex) const
{
    if (vertexIndex > this->vertex.size())
    {
        std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] Vertex with vertex index " << vertexIndex
                  << " does not exist. There are only " << this->vertex.size() << " vertices." << std::endl;
        return error;
    }
    if (triangleIndex > getNumberOfTriangles())
    {
        std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The " << triangleIndex
                  << "-th triangle does not exist. There are only " << getNumberOfTriangles() << " triangles."
                  << std::endl;
        return error;
    }
    if (*(this->vertex[m_sorting[vertexIndex - 1]]) == *(this->face[triangleIndex]->vertices[0]))
    { return first; }
    if (*(this->vertex[m_sorting[vertexIndex - 1]]) == *(this->face[triangleIndex]->vertices[1]))
    { return second; }
    if (*(this->vertex[m_sorting[vertexIndex - 1]]) == *(this->face[triangleIndex]->vertices[2]))
    { return third; }
    return error;
}

template<class T>
const std::queue<typename gsHalfEdgeMesh<T>::Halfedge>
gsHalfEdgeMesh<T>::getOppositeHalfedges(const std::size_t vertexIndex, const bool innerVertex) const
{
    std::queue<Halfedge> oppositeHalfedges;
    if (vertexIndex > this->vertex.size())
    {
        std::cerr << "Error: [" << __PRETTY_FUNCTION__ << "] The vertex with index " << vertexIndex
                  << " does not exist. There are only " << this->vertex.size() << " vertices." << std::endl;
        return oppositeHalfedges;
    }
    else if (vertexIndex > m_n && innerVertex)
    {
        std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] Inner vertex with index 'vertexIndex' = " << vertexIndex
                  << "is not an inner vertex. There are only " << m_n << " inner vertices." << std::endl;
    }

    for (std::size_t i = 0; i < getNumberOfTriangles(); i++)
    {
        switch (isTriangleVertex(vertexIndex, i))
        {
            case first:
                oppositeHalfedges.push(Halfedge(getGlobalVertexIndex(3, i),
                                                                 getGlobalVertexIndex(2, i),
                                                                 getHalfedgeLength(getGlobalVertexIndex(3, i),
                                                                                   getGlobalVertexIndex(2, i))));
                break;
            case second:
                oppositeHalfedges.push(Halfedge(getGlobalVertexIndex(1, i),
                                                                 getGlobalVertexIndex(3, i),
                                                                 getHalfedgeLength(getGlobalVertexIndex(1, i),
                                                                                   getGlobalVertexIndex(3, i))));
                break;
            case third:
                oppositeHalfedges.push(Halfedge(getGlobalVertexIndex(2, i),
                                                                 getGlobalVertexIndex(1, i),
                                                                 getHalfedgeLength(getGlobalVertexIndex(2, i),
                                                                                   getGlobalVertexIndex(1, i))));
                break;
            default:
                //not supposed to show up
                break;
        }
    }
    return oppositeHalfedges;
}

//*****************************************************************************************************
//*****************************************************************************************************
//*******************THE******INTERN******FUNCTIONS******ARE******NOW******FOLLOWING*******************
//*****************************************************************************************************
//*****************************************************************************************************
template<class T>
bool gsHalfEdgeMesh<T>::isBoundaryVertex(const std::size_t internVertexIndex) const
{
    if (internVertexIndex > this->vertex.size() - 1)
    {
        std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] Vertex with intern vertex index = " << internVertexIndex
                  << " does not exist. There are only " << this->vertex.size() << " vertices." << std::endl;
        return false;
    }
    else
        return m_boundary.isVertexContained(internVertexIndex);
}

template<class T>
std::size_t gsHalfEdgeMesh<T>::getInternVertexIndex(const gsMesh<real_t>::gsVertexHandle &vertex) const
{
    std::size_t internVertexIndex = 0;
    for (std::size_t i = 0; i < this->vertex.size(); i++)
    {
        if ((*(this->vertex[i])) == *vertex)
            return internVertexIndex;
        internVertexIndex++;
    }
    if (internVertexIndex > this->vertex.size() - 1)
    {
        std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The ESS_IO::IO_Vertex 'vertex' = (" << vertex->x()
                  << ", " << vertex->y() << ", " << vertex->z() << ") is not contained in gsHalfEdgeMesh vertices"
                  << std::endl;
        return 0;
    }
    return 0;
}

template<class T>
const typename gsHalfEdgeMesh<T>::Halfedge
gsHalfEdgeMesh<T>::getInternHalfedge(const gsMesh<real_t>::gsFaceHandle &triangle, std::size_t numberOfHalfedge) const
{
    std::size_t index1 = getInternVertexIndex((triangle->vertices[0]));
    std::size_t index2 = getInternVertexIndex((triangle->vertices[1]));
    std::size_t index3 = getInternVertexIndex((triangle->vertices[2]));
    if (numberOfHalfedge < 1 || numberOfHalfedge > 3)
    {
        std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The inputted number of the halfedge " << numberOfHalfedge
                  << "  is supposed to be 1,2 or 3. Because input was not expected, first halfedge is returned."
                  << std::endl;
        numberOfHalfedge = 1;
    }
    if (numberOfHalfedge == 1)
    {
        return Halfedge(index2, index1,
                                         gsVector3d<real_t>(
                                             triangle->vertices[1]->x() - triangle->vertices[0]->x(),
                                             triangle->vertices[1]->y() - triangle->vertices[0]->y(),
                                             triangle->vertices[1]->z() - triangle->vertices[0]->z()).norm());
    }
    if (numberOfHalfedge == 2)
    {
        return Halfedge(index3, index2,
                                         gsVector3d<real_t>(
                                             triangle->vertices[2]->x() - triangle->vertices[1]->x(),
                                             triangle->vertices[2]->y() - triangle->vertices[1]->y(),
                                             triangle->vertices[2]->z() - triangle->vertices[1]->z()).norm());
    }
    if (numberOfHalfedge == 3)
    {
        return Halfedge(index1, index3,
                                         gsVector3d<real_t>(
                                             triangle->vertices[0]->x() - triangle->vertices[2]->x(),
                                             triangle->vertices[0]->y() - triangle->vertices[2]->y(),
                                             triangle->vertices[0]->z() - triangle->vertices[2]->z()).norm());
    }
    return Halfedge();
}

template<class T>
void gsHalfEdgeMesh<T>::sortVertices()
{
    std::size_t numberOfInnerVerticesFound = 0;
    std::vector<std::size_t> sorting(this->vertex.size(), 0);
    m_sorting = sorting;
    std::vector<std::size_t> inverseSorting(this->vertex.size(), 0);
    m_inverseSorting = inverseSorting;
    std::list<std::size_t> boundaryVertices = m_boundary.getVertexIndices();
    for (std::size_t i = 0; i < this->vertex.size(); i++)
    {
        if (!isBoundaryVertex(i))
        {
            numberOfInnerVerticesFound++;
            m_sorting[numberOfInnerVerticesFound - 1] = i;
            m_inverseSorting[i] = numberOfInnerVerticesFound;
        }
    }
    for (std::size_t i = 0; i < getNumberOfBoundaryVertices(); i++)
    {
        m_sorting[m_n + i] = boundaryVertices.front();
        m_inverseSorting[boundaryVertices.front()] = m_n + i + 1;
        boundaryVertices.pop_front();
    }
}

//***********************************************
//************ nested class Boundary ************
//***********************************************

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
        std::cout << "Warning: [" << __PRETTY_FUNCTION__
                  << "] Boundary is not closed although it should be. End points are: " << std::endl
                  << m_boundary.getFirstHalfedge().getOrigin() << std::endl << " and "
                  << m_boundary.getLastHalfedge().getEnd() << std::endl;
}

//********* private ***********

template<class T>
const std::list<typename gsHalfEdgeMesh<T>::Halfedge> gsHalfEdgeMesh<T>::Boundary::findNonTwinHalfedges(const std::vector<typename gsHalfEdgeMesh<T>::Halfedge> &allHalfedges)
{
    std::queue<Halfedge> queue0;
    std::queue<Halfedge> queue1;
    bool actualQueue = 0;
    std::list<Halfedge> nonTwinHalfedges;
    for (std::size_t i = 0; i < allHalfedges.size(); ++i)
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

//***********************************************
//************ nested class Chain ************
//***********************************************

template<class T>
bool gsHalfEdgeMesh<T>::Chain::isClosed() const
{
    if (m_chainedHalfedges.empty())
    {
        std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                  << std::endl;
        return true;
    }
    return (m_chainedHalfedges.front().getOrigin() == m_chainedHalfedges.back().getEnd());
}

template<class T>
std::size_t gsHalfEdgeMesh<T>::Chain::getNumberOfVertices() const
{
    if (this->isClosed())
        return m_chainedHalfedges.size();
    else
        return m_chainedHalfedges.size() + 1;
}

template<class T>
double gsHalfEdgeMesh<T>::Chain::getLength() const
{
    double length = 0;
    for (typename std::list<Halfedge>::const_iterator it = m_chainedHalfedges.begin();
         it != m_chainedHalfedges.end(); ++it)
    {
        length += it->getLength();
    }
    return length;
}

template<class T>
const std::vector<double> gsHalfEdgeMesh<T>::Chain::getHalfedgeLengths() const
{
    if (this->isEmpty())
    {
        std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                  << std::endl;
    }
    std::vector<double> lengths;
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
        std::cerr << "Error: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                  << std::endl;
    }
    return m_chainedHalfedges.front();
}

template<class T>
const typename gsHalfEdgeMesh<T>::Halfedge& gsHalfEdgeMesh<T>::Chain::getLastHalfedge() const
{
    if (this->isEmpty())
    {
        std::cerr << "Error: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                  << std::endl;
    }
    return m_chainedHalfedges.back();
}

template<class T>
const std::list<std::size_t> gsHalfEdgeMesh<T>::Chain::getVertexIndices() const
{
    if (this->isEmpty())
    {
        std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                  << std::endl;
    }
    std::list<std::size_t> vertexIndices;
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
double gsHalfEdgeMesh<T>::Chain::getShortestDistanceBetween(std::size_t i, std::size_t j) const
{
    if (this->isEmpty())
    {
        std::cout << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                  << std::endl;
        return 0;
    }
    if (i > this->getNumberOfVertices() || j > this->getNumberOfVertices() || i < 1 || j < 1)
    {
        std::cout << "Error: [" << __PRETTY_FUNCTION__ << "] FirstIndex: " << i << " and second index: "
                  << j
                  << " must be positiv integers smaller than the number of points of the chain, which is "
                  << this->getNumberOfVertices() << "." << std::endl;
        return 0;
    }
    if (i > j)
        std::swap(i, j);//myFunctions::orderIntegers(i, j);
    double distance = 0;
    std::vector<double> l = this->getHalfedgeLengths();
    for (std::size_t z = i - 1; z < j - 1; z++)
    {
        distance += l[z];
    }
    if (this->isClosed() && (this->getLength() - distance < distance - 1e-8))
    {
        distance = this->getLength() - distance;
    }
    return distance;
}

template<class T>
double gsHalfEdgeMesh<T>::Chain::getDistanceBetween(std::size_t i, std::size_t j) const
{
    if (this->isEmpty())
    {
        std::cout << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                  << std::endl;
        return 0;
    }
    if (i > this->getNumberOfVertices() || j > this->getNumberOfVertices() || i < 1 || j < 1)
    {
        std::cout << "Error: [" << __PRETTY_FUNCTION__ << "] FirstIndex: " << i << " and second index: "
                  << j
                  << " must be positiv integers smaller than the number of points of the chain, which is "
                  << this->getNumberOfVertices() << "." << std::endl;
        return 0;
    }
    bool ordered = (i < j);
    if (!ordered)
        std::swap(i, j);//myFunctions::orderIntegers(i, j);
    double distance = 0;
    std::vector<double> l = this->getHalfedgeLengths();
    for (std::size_t z = i - 1; z < j - 1; z++)
    {
        distance += l[z];
    }
    if (this->isClosed() && !ordered)
    {
        distance = this->getLength() - distance;
    }
    else if (!ordered)
        std::cout << "Warning: [" << __PRETTY_FUNCTION__
                  << "] The chain is supposed to be closed in case the input is not ordered." << std::endl;
    return distance;
}

template<class T>
bool gsHalfEdgeMesh<T>::Chain::isVertexContained(const std::size_t &vertexIndex) const
{
    if (this->isEmpty())
    {
        std::cerr << "Warning: [" << __PRETTY_FUNCTION__ << "] The chain does not store any halfedges yet."
                  << std::endl;
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
        std::cout << "Warning: [" << __PRETTY_FUNCTION__
                  << "] This halfedge is not appendable at the beginning." << std::endl;
        std::cout << "The first halfedge of the chain has origin " << this->getFirstHalfedge().getOrigin()
                  << " and prevHalfedge has the end " << prevHalfedge.getEnd() << "." << std::endl;
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
        std::cout << "Warning: [" << __PRETTY_FUNCTION__ << "] This halfedge is not appendable at the end."
                  << std::endl;
        std::cout << "The last halfedge of the chain has end " << this->getLastHalfedge().getEnd()
                  << " and nextHalfedge has the origin " << nextHalfedge.getOrigin() << "." << std::endl;
    }
    else
    {
        m_chainedHalfedges.push_back(nextHalfedge);
    }
}

// TODO: refactor rest of inner classes

} // namespace gismo