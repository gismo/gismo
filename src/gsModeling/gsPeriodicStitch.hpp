/** @file gsPeriodicStitch.hpp

    @brief Provides implementation of the gsPeriodicStitch class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris

*/

#include <gsModeling/gsPeriodicStitch.h>

namespace gismo
{

/* Nested class Neighbourhood */

template<class T>
std::vector<size_t> gsPeriodicStitch<T>::Neighbourhood::computeCorrections(const std::vector<size_t>& stitchIndices,
                                                                           const LocalNeighbourhood& localNeighbourhood) const
{
    auto indexIt = std::find(stitchIndices.begin(), stitchIndices.end(), localNeighbourhood.getVertexIndex());

    if(indexIt == stitchIndices.end()) // Not on the stitch, nothing to do.
    {
        return std::vector<size_t>();
    }

    std::list<size_t> result;
    std::list<size_t> neighbours = localNeighbourhood.getVertexIndicesOfNeighbours();

    if(indexIt == stitchIndices.begin()) // In the beginning of the stitch.
    {
        auto nextOnStitch = std::find(neighbours.begin(), neighbours.end(), *std::next(indexIt));
        // (Assuming that the stitch has at least two vertices.)
        for(auto it=nextOnStitch; it!=neighbours.end(); ++it)
            result.push_back(*it);
    }
    else if(std::next(indexIt) == stitchIndices.end()) // In the end of the stitch.
    {
        auto prevOnStitch = std::find(neighbours.begin(), neighbours.end(), *std::prev(indexIt));
        // (Again assuming the stitch to have at least two vertices.)
        for(auto it=neighbours.begin(); it!=prevOnStitch; ++it)
            result.push_back(*it);
    }
    else // In the middle of the stitch.
    {
        while(neighbours.front() != *std::next(indexIt))
        {
            neighbours.push_back(neighbours.front());
            neighbours.pop_front();
        }

        auto prevOnStitch = std::find(neighbours.begin(), neighbours.end(), *std::prev(indexIt));
        for(auto it=neighbours.begin(); it!=prevOnStitch; ++it)
            result.push_back(*it);
    }

    // Other stitch vertices can still be present in the neighbourhood.
    for(auto it=stitchIndices.begin(); it!=stitchIndices.end(); ++it)
        result.remove(*it);

    std::vector<size_t> finalResult;
    finalResult.reserve(result.size());
    for(auto it=result.begin(); it!=result.end(); ++it)
        finalResult.push_back(*it);

    return finalResult;
}

template<class T>
gsPeriodicStitch<T>::Neighbourhood::Neighbourhood(const gsHalfEdgeMesh<T> & meshInfo,
                                                  const std::vector<size_t>& stitchIndices,
                                                  gsSparseMatrix<int>& corrections,
                                                  const size_t parametrizationMethod)
    : gsParametrization<T>::Neighbourhood(meshInfo, parametrizationMethod)
{
    // We re-do a little of the work done already in the constructor of the parent class.
    // Alternatively, we could provide a constructor of the parent class setting m_basicInfos
    // and do everything else here, much as we did when this was a part of the parent class.
    // Cf., e.g., fcacc860ee28edd608e841af5aeb74dacc90e006 for a reference.

    index_t N = meshInfo.getNumberOfVertices();
    corrections.resize(N, N);
    corrections.setZero(); // important, otherwise you might end up with artifacts

    for(size_t i=1; i <= meshInfo.getNumberOfVertices(); i++)
    {
        LocalNeighbourhood localNeighbourhood = (i <= meshInfo.getNumberOfInnerVertices()) ?
            LocalNeighbourhood(meshInfo, i) :
            LocalNeighbourhood(meshInfo, i, 0);

        std::vector<size_t> corr = computeCorrections(stitchIndices, localNeighbourhood);
        
        for(auto it=corr.begin(); it!=corr.end(); ++it)
        {
            corrections(i-1, *it-1) = 1;
            corrections(*it-1, i-1) = -1;
        }
    }
}

template <class T>
void gsPeriodicStitch<T>::compute()
{
    calculate(this->m_options.getInt("parametrizationMethod"));
}

template<class T>
void gsPeriodicStitch<T>::calculate(const size_t paraMethod)
{
    size_t n = this->m_mesh.getNumberOfInnerVertices();
    size_t N = this->m_mesh.getNumberOfVertices();

    Neighbourhood neighbourhood(this->m_mesh, m_stitchIndices, m_corrections, paraMethod);

    this->initParameterPoints();

    /// Solve.
    constructAndSolveEquationSystem(neighbourhood, n, N);
}

template <class T>
void gsPeriodicStitch<T>::constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
                                                                         const size_t n,
                                                                         const size_t N)
{
    std::vector<T> lambdas;
    gsMatrix<T> LHS(N, N);
    gsMatrix<T> RHS(N, 2);

    // prevent Valgrind warnings
    LHS.setZero();
    RHS.setZero();

    // interior points
    for (size_t i = 0; i < n; i++)
    {
        lambdas = neighbourhood.getLambdas(i);
        for (size_t j = 0; j < N; j++)
        {
            LHS(i, j) = ( i==j ? (T)(1) : -lambdas[j] );

            // If your neighbour is across the stitch, its contributions appear
            // on the right hand-side multiplied by +1 or -1. Write the equations
            // down if it is unclear. (-;
            if(m_corrections(i, j) == 1)
                RHS(i, 0) -= lambdas[j];
            else if(m_corrections(i, j) == -1)
                RHS(i, 0) += lambdas[j];
        }
    }

    // points on the lower and upper boundary
    for (size_t i=n; i<N; i++)
    {
        LHS(i, i) = (T)(1);
        RHS.row(i) = this->m_parameterPoints[i];
    }

    // Solve the system and save the results.
    gsEigen::PartialPivLU<typename gsMatrix<T>::Base> LU = LHS.partialPivLu();
    gsMatrix<T> sol = LU.solve(RHS);
    for (size_t i = 0; i < n; i++)
    {
        this->m_parameterPoints[i] << sol(i, 0), sol(i, 1);
    }
}

template<class T>
gsMesh<T> gsPeriodicStitch<T>::createUnfoldedFlatMesh() const
{
    typedef typename gsMesh<T>::VertexHandle       VertexHandle;
    typedef typename gsParametrization<T>::Point2D Point2D;

    gsMesh<T> result;
    for(size_t i=0; i<this->m_mesh.getNumberOfTriangles(); i++)
    {
        std::vector<size_t> vertices;
        for(size_t j=1; j<=3; ++j)
        {
            vertices.push_back(this->m_mesh.getGlobalVertexIndex(j, i));
        }
        bool nearStitchTriangle = (edgeIsInCorrections(vertices[0]-1, vertices[1]-1) ||
                                   edgeIsInCorrections(vertices[1]-1, vertices[2]-1) ||
                                   edgeIsInCorrections(vertices[2]-1, vertices[0]-1));

        VertexHandle v[3];
        for (size_t j = 1; j <= 3; ++j)
        {               
            const Point2D& point = gsParametrization<T>::getParameterPoint(vertices[j-1]);
            // The near-stitch triangles get their stitch vertices shifted by 1 to the left.
            if( nearStitchTriangle && isOnStitch(vertices[j-1]) )
                v[j - 1] = result.addVertex(point[0] + (T)(1), point[1]);
            else
                v[j - 1] = result.addVertex(point[0],     point[1]);
        }
        result.addFace( v[0],  v[1],  v[2]);
    }

    return result.cleanMesh();
}

} // namespace gismo
