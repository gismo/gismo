/** @file gsPeriodicOverlap.hpp

    @brief Provides implementation of the gsPeriodicOverlap class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris

*/

#include <gismo.h>

namespace gismo
{

template <class T>
void gsPeriodicOverlap<T>::compute()
{
    calculate(this->m_options.getInt("parametrizationMethod"));
}

template<class T>
void gsPeriodicOverlap<T>::calculate(const size_t paraMethod)
{
    size_t n = this->m_mesh.getNumberOfInnerVertices();
    size_t N = this->m_mesh.getNumberOfVertices();

    Neighbourhood neighbourhood(this->m_mesh, paraMethod);

    this->initParameterPoints();

    // Construct the twins.
    constructTwins();

    // Solve.
    constructAndSolveEquationSystem(neighbourhood, n, N);
}




template <class T>
void gsPeriodicOverlap<T>::constructTwinsBetween(size_t& currentNrAllVertices,
                                                 std::list<size_t> vertexIndices,
                                                 size_t from, size_t to,
                                                 bool rightHandSide)
{
    // TODO: The whiles do not check if the sought member is indeed in
    // the list (danger of an infinite loop).

    // Rotate the vertexIndices so as to start from from.
    while(vertexIndices.front() != from)
    {
        vertexIndices.push_back(vertexIndices.front());
        vertexIndices.pop_front();
    }

    // Push the corresponding pairs to the twin vector.
    // Note that if we would do std::prev on begin(), we should not dereference.
    // The error is visible only with -DGISMO_EXTRA_DEBUG=ON.
    for(std::list<size_t>::const_iterator it=vertexIndices.begin();
        it == vertexIndices.begin() || *std::prev(it) != to;
        ++it)
    {
        size_t twin = findTwin(*it);
        if(rightHandSide)
            m_twins.push_back(std::pair<size_t, size_t>(twin, ++currentNrAllVertices));
        else
            m_twins.push_back(std::pair<size_t, size_t>(++currentNrAllVertices, twin));
    }

}

template<class T>
void gsPeriodicOverlap<T>::constructTwins()
{
    // vertex with parameter v = 0 and lowest u value
    gsVertexHandle uMinv0 = this->m_mesh.getVertex(this->m_indicesV0.front());

    // vertex with parameter v = 0 and highest u value
    gsVertexHandle uMaxv0 = this->m_mesh.getVertex(this->m_indicesV0.back());

    // vertex with parameter v = 1 and lowest u value
    gsVertexHandle uMinv1 = this->m_mesh.getVertex(this->m_indicesV1.front());

    // vertex with parameter v = 1 and highest u value
    gsVertexHandle uMaxv1 = this->m_mesh.getVertex(this->m_indicesV1.back());

    const std::list<size_t> vertexIndices = m_overlapHEM.getBoundaryVertexIndices();
    size_t currentNrAllVertices = this->m_mesh.getNumberOfVertices();

    // Construct the twins on the right boundary of the overlap.
    constructTwinsBetween(currentNrAllVertices, vertexIndices, uMinv1, uMinv0, true);

    // Construct the twins on the left boundary of the overlap.
    constructTwinsBetween(currentNrAllVertices, vertexIndices, uMaxv0, uMaxv1, false);
}

template <class T>
void gsPeriodicOverlap<T>::constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
                                                           const size_t n,
                                                           const size_t N)
{
    size_t numTwins = m_twins.size();
    gsMatrix<T> LHS(N + numTwins, N + numTwins);
    gsMatrix<T> RHS(N + numTwins, 2);
    // prevent Valgrind warnings
    LHS.setZero();
    RHS.setZero();
    std::vector<T> lambdas;

    // interior points
    for (size_t i = 0; i < n; i++)
    {
        lambdas = neighbourhood.getLambdas(i);
        updateLambdasWithTwins(lambdas, i+1);

        for (size_t j = 0; j < N + numTwins; j++)
            LHS(i, j) = ( i==j ? T(1) : -lambdas[j] );
    }

    // points on the lower and upper boundary
    for (size_t i=n; i<N; i++)
    {
        LHS(i, i)  = T(1);
        RHS.row(i) = this->m_parameterPoints[i];
    }

    // points on the overlap
    for (size_t i=N; i<N+numTwins; i++)
    {
        size_t first   = m_twins[i-N].first-1;
        size_t second  = m_twins[i-N].second-1;

        LHS(i, first)  = T( 1);
        LHS(i, second) = T(-1);

        RHS(i, 0)      = T(-1);
        RHS(i, 1)      = T( 0);
    }

    Eigen::PartialPivLU<typename gsMatrix<T>::Base> LU = LHS.partialPivLu();
    gsMatrix<T> sol = LU.solve(RHS);
    for (size_t i = 0; i < N; i++)
    {
        this->m_parameterPoints[i] << sol(i, 0), sol(i, 1);
    }
}

template <class T>
void gsPeriodicOverlap<T>::updateLambdasWithTwins(std::vector<T>& lambdas,
                                                  size_t vertexId) const
{
    lambdas.reserve(lambdas.size() + m_twins.size());
    for(size_t i=0; i<m_twins.size(); i++)
        lambdas.push_back(0);

    // Determine, whether vertexId is on the left or right side of the overlap.
    bool isLeft  = false;
    bool isRight = false;

    for(auto it=m_twins.begin(); it!=m_twins.end(); ++it)
    {
        if(it->first == vertexId)
        {
            isRight = true;
            break;
        }
        else if(it->second == vertexId)
        {
            isLeft = true;
            break;
        }
    }


    for(size_t i=0; i<m_twins.size(); i++)
    {
        size_t first=m_twins[i].first-1;
        size_t second=m_twins[i].second-1;

        // Left vertex swaps all its right neighbours.
        if(isRight && first > second && lambdas[second] != 0)
        {
            lambdas[first] = lambdas[second];
            lambdas[second] = 0;
        }
        // Right vertex swaps all its left neighbours
        else if(isLeft && first < second && lambdas[first] != 0)
        {
            lambdas[second] = lambdas[first];
            lambdas[first] = 0;
        }
        // Nothing happens to vertices that are not on the overlap.         
    }
}

template<class T>
gsMesh<T> gsPeriodicOverlap<T>::createFlatMesh() const
{
    // Remember the vertices on the overlap boundaries.
    std::vector<size_t> left, right;
    for(auto it=m_twins.begin(); it!=m_twins.end(); ++it)
    {
        if(it->first < it->second)
            left.push_back(it->first);
        else
            right.push_back(it->second);
    }

    typename gsPeriodicParametrization<T>::FlatMesh display(createExtendedFlatMesh(left, right));
    return display.createRestrictedFlatMesh();
}

template<class T>
gsMesh<T> gsPeriodicOverlap<T>::createExtendedFlatMesh(const std::vector<size_t>& right,
                                                       const std::vector<size_t>& left) const
{
    typedef typename gsParametrization<T>::Point2D Point2D;

    gsMesh<T> midMesh;
    midMesh.reserve(3 * this->m_mesh.getNumberOfTriangles(), this->m_mesh.getNumberOfTriangles(), 0);

    for (size_t i = 0; i < this->m_mesh.getNumberOfTriangles(); i++)
    {
        size_t vInd[3];

        // the indices of the current triangle vertices that are among left or right, respectively
        std::vector<size_t> lVert, rVert;
        for (size_t j = 1; j <= 3; ++j)
        {
            vInd[j-1] = this->m_mesh.getGlobalVertexIndex(j, i);
            if(std::find(right.begin(), right.end(), vInd[j-1]) != right.end())
                rVert.push_back(j-1);
            if(std::find(left.begin(),  left.end(),  vInd[j-1]) != left.end())
                lVert.push_back(j-1);
        }

        // Is the triangle inside overlap?
        if(lVert.size() > 0 && rVert.size() > 0 && lVert.size() + rVert.size() == 3)
        {
            // Make two shifted copies of the triangle.
            typename gsMesh<T>::VertexHandle mvLft[3], mvRgt[3];
            
            for (size_t j=0; j<3; ++j)
            {
                const Point2D vertex = gsParametrization<T>::getParameterPoint(vInd[j]);
                if(std::find(rVert.begin(), rVert.end(), j) != rVert.end())
                {
                    mvLft[j] = midMesh.addVertex(vertex[0],        vertex[1]);
                    mvRgt[j] = midMesh.addVertex(vertex[0]+(T)(1), vertex[1]);
                }
                else
                {
                    mvLft[j] = midMesh.addVertex(vertex[0]-(T)(1), vertex[1]);
                    mvRgt[j] = midMesh.addVertex(vertex[0],        vertex[1]);
                }
            }
            midMesh.addFace(mvLft[0], mvLft[1], mvLft[2]);
            midMesh.addFace(mvRgt[0], mvRgt[1], mvRgt[2]);
        }
        else
        {
            // Make just one triangle.
            typename gsMesh<T>::VertexHandle mv[3];
            for (size_t j=0; j<3; ++j)
            {
                const Point2D vertex = gsParametrization<T>::getParameterPoint(vInd[j]);
                mv[j] = midMesh.addVertex(vertex[0], vertex[1]);
            }
            midMesh.addFace(mv[0], mv[1], mv[2]);
        }
    }
    return midMesh.cleanMesh();
}

} // namespace gismo
