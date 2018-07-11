/** @file gsHalfEdgeMesh.hpp

    @brief Provides implementation of the gsHalfEdgeMesh class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl
*/

#pragma once

//#include <gsUtils/gsMesh/gsMesh.h>
//#include "gsHalfEdgeMesh.h"

namespace gismo
{
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

gsHalfEdgeMesh::gsHalfEdgeMesh(const gsMesh<> &mesh)
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

std::size_t gsHalfEdgeMesh::getNumberOfVertices() const
{
    return this->vertex.size();
}

std::size_t gsHalfEdgeMesh::getNumberOfTriangles() const
{
    return this->face.size();
}

std::size_t gsHalfEdgeMesh::getNumberOfInnerVertices() const
{
    return m_n;
}

std::size_t gsHalfEdgeMesh::getNumberOfBoundaryVertices() const
{
    return m_boundary.getNumberOfVertices();
}

const gsMesh<>::gsVertexHandle &gsHalfEdgeMesh::getVertex(const std::size_t vertexIndex) const
{
    if (vertexIndex > this->vertex.size())
    {
        std::cerr << "Error: [" << __PRETTY_FUNCTION__ << "] Vertex with index 'vertexIndex'=" << vertexIndex
                  << " does not exist. There are only " << this->vertex.size() << " vertices." << std::endl;
    }
    return ((this->vertex[m_sorting[vertexIndex - 1]]));
}

std::size_t gsHalfEdgeMesh::getVertexIndex(const gsMesh<>::gsVertexHandle &vertex) const
{
    return m_inverseSorting[getInternVertexIndex(vertex)];
}

std::size_t gsHalfEdgeMesh::getGlobalVertexIndex(std::size_t localVertexIndex, std::size_t triangleIndex) const
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

double gsHalfEdgeMesh::getBoundaryLength() const
{
    return m_boundary.getLength();
}

bool rangeCheck(const std::vector<int> &corners, const std::size_t minimum, const std::size_t maximum)
{
    for (std::vector<int>::const_iterator it = corners.begin(); it != corners.end(); it++)
    {
        if (*it < minimum || *it > maximum)
        { return false; }
    }
    return true;
}

std::vector<double> gsHalfEdgeMesh::getCornerLengths(std::vector<int> &corners) const
{
    std::size_t B = getNumberOfBoundaryVertices();
    if (!rangeCheck(corners, 1, B))
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

double gsHalfEdgeMesh::getShortestBoundaryDistanceBetween(std::size_t i, std::size_t j) const
{
    return m_boundary.getShortestDistanceBetween(i, j);
}

const std::vector<double> gsHalfEdgeMesh::getBoundaryChordLengths() const
{
    return m_boundary.getHalfedgeLengths();
}

double gsHalfEdgeMesh::getHalfedgeLength(std::size_t originVertexIndex, std::size_t endVertexIndex) const
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

triangleVertexIndex gsHalfEdgeMesh::isTriangleVertex(std::size_t vertexIndex, std::size_t triangleIndex) const
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

const std::queue<gsHalfEdgeMesh::Boundary::Chain::Halfedge>
gsHalfEdgeMesh::getOppositeHalfedges(const std::size_t vertexIndex, const bool innerVertex) const
{
    std::queue<Boundary::Chain::Halfedge> oppositeHalfedges;
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
                oppositeHalfedges.push(Boundary::Chain::Halfedge(getGlobalVertexIndex(3, i),
                                                                 getGlobalVertexIndex(2, i),
                                                                 getHalfedgeLength(getGlobalVertexIndex(3, i),
                                                                                   getGlobalVertexIndex(2, i))));
                break;
            case second:
                oppositeHalfedges.push(Boundary::Chain::Halfedge(getGlobalVertexIndex(1, i),
                                                                 getGlobalVertexIndex(3, i),
                                                                 getHalfedgeLength(getGlobalVertexIndex(1, i),
                                                                                   getGlobalVertexIndex(3, i))));
                break;
            case third:
                oppositeHalfedges.push(Boundary::Chain::Halfedge(getGlobalVertexIndex(2, i),
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

bool gsHalfEdgeMesh::isBoundaryVertex(const std::size_t internVertexIndex) const
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

std::size_t gsHalfEdgeMesh::getInternVertexIndex(const gsMesh<real_t>::gsVertexHandle &vertex) const
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

const gsHalfEdgeMesh::Boundary::Chain::Halfedge
gsHalfEdgeMesh::getInternHalfedge(const gsMesh<real_t>::gsFaceHandle &triangle, std::size_t numberOfHalfedge) const
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
        return Boundary::Chain::Halfedge(index2, index1,
                                         gsVector3d<real_t>(
                                             triangle->vertices[1]->x() - triangle->vertices[0]->x(),
                                             triangle->vertices[1]->y() - triangle->vertices[0]->y(),
                                             triangle->vertices[1]->z() - triangle->vertices[0]->z()).norm());
    }
    if (numberOfHalfedge == 2)
    {
        return Boundary::Chain::Halfedge(index3, index2,
                                         gsVector3d<real_t>(
                                             triangle->vertices[2]->x() - triangle->vertices[1]->x(),
                                             triangle->vertices[2]->y() - triangle->vertices[1]->y(),
                                             triangle->vertices[2]->z() - triangle->vertices[1]->z()).norm());
    }
    if (numberOfHalfedge == 3)
    {
        return Boundary::Chain::Halfedge(index1, index3,
                                         gsVector3d<real_t>(
                                             triangle->vertices[0]->x() - triangle->vertices[2]->x(),
                                             triangle->vertices[0]->y() - triangle->vertices[2]->y(),
                                             triangle->vertices[0]->z() - triangle->vertices[2]->z()).norm());
    }
    return Boundary::Chain::Halfedge();
}

void gsHalfEdgeMesh::sortVertices()
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
} // namespace gismo