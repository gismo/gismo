/** @file gsPeriodicParametrizationStitch.hpp

    @brief Provides implementation of the gsPeriodicParametrizationStitch class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris

*/

#include <gismo.h>
#include "gsModeling/gsPeriodicParametrizationStitch.h"

namespace gismo
{

/* Nested class Neighbourhood */



template<class T>
std::vector<size_t> gsPeriodicParametrizationStitch<T>::Neighbourhood::computeCorrections(const std::vector<size_t>& stitchIndices,
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
gsPeriodicParametrizationStitch<T>::Neighbourhood::Neighbourhood(const gsHalfEdgeMesh<T> & meshInfo,
						   const std::vector<size_t>& stitchIndices,
						   std::vector<std::vector<size_t> >& posCorrections,
						   std::vector<std::vector<size_t> >& negCorrections,
						   const size_t parametrizationMethod)
    : gsParametrization<T>::Neighbourhood(meshInfo, parametrizationMethod)
{
    // We re-do a little of the work done already in the constructor of the parent class.
    // Alternatively, we could provide a constructor of the parent class setting m_basicInfos
    // and do everything else here, much as we did when this was a part of the parent class.
    // Cf., e.g., fcacc860ee28edd608e841af5aeb74dacc90e006 for a reference.
    posCorrections.resize(meshInfo.getNumberOfVertices());
    negCorrections.resize(meshInfo.getNumberOfVertices());

    for(size_t i=1; i <= meshInfo.getNumberOfInnerVertices(); i++)
    {
	LocalNeighbourhood localNeighbourhood(meshInfo, i);

	posCorrections[i-1] = computeCorrections(stitchIndices, localNeighbourhood);
	for(size_t j=0; j < posCorrections[i-1].size(); j++)
	{
	    negCorrections[posCorrections[i-1][j]-1].push_back(i);
	}
    }

    for(size_t i=meshInfo.getNumberOfInnerVertices()+1; i<= meshInfo.getNumberOfVertices(); i++)
    {
	LocalNeighbourhood localNeighbourhood(meshInfo, i, 0);

	posCorrections[i-1] = computeCorrections(stitchIndices, localNeighbourhood);
	for(size_t j=0; j < posCorrections[i-1].size(); j++)
	{
	    negCorrections[posCorrections[i-1][j]-1].push_back(i);
	}
    }
}

template <class T>
gsPeriodicParametrizationStitch<T>& gsPeriodicParametrizationStitch<T>::compute_periodic_stitch(std::string bottomFile,
												std::string topFile,
												std::string stitchFile)
{
    // Read the indices and u-coordinates of the points with v=0.
    std::vector<size_t> indicesV0;
    std::vector<T> valuesV0;
    gsParametrization<T>::readIndicesAndValues(bottomFile, indicesV0, valuesV0);

    // Read the indices and u-coordinates of the points with v=1.
    std::vector<size_t> indicesV1;
    std::vector<T> valuesV1;
    gsParametrization<T>::readIndicesAndValues(topFile, indicesV1, valuesV1);

    // Read the indices of the points on the stitch.
    std::vector<size_t> stitchIndices = gsParametrization<T>::readIndices(stitchFile);

    // Calculation itself.
    calculate_periodic_stitch(this->m_options.getInt("parametrizationMethod"),
			      indicesV0, valuesV0, indicesV1, valuesV1, stitchIndices);
    return *this;
}

template<class T>
void gsPeriodicParametrizationStitch<T>::calculate_periodic_stitch(const size_t paraMethod,
								   const std::vector<size_t>& indicesV0,
								   const std::vector<T>& valuesV0,
								   const std::vector<size_t>& indicesV1,
								   const std::vector<T>& valuesV1,
								   const std::vector<size_t>& stitchIndices)
{
    typedef typename gsParametrization<T>::Point2D       Point2D ;

    size_t n = this->m_mesh.getNumberOfInnerVertices();
    size_t N = this->m_mesh.getNumberOfVertices();

    Neighbourhood neighbourhood(this->m_mesh, stitchIndices, m_posCorrections, m_negCorrections, paraMethod);

    this->m_parameterPoints.reserve(N);
    for (size_t i = 1; i <= n; i++)
    {
	this->m_parameterPoints.push_back(Point2D(0, 0, i));
    }

    // Add the parameters of the boundary points.
    GISMO_ASSERT(indicesV0.size() == valuesV0.size(), "Different sizes of u0.");
    GISMO_ASSERT(indicesV1.size() == valuesV1.size(), "Different sizes of u1.");
    GISMO_ASSERT(indicesV0.size() + indicesV1.size() == this->m_mesh.getNumberOfBoundaryVertices(),
		 "Not prescribing all boundary points.");

    size_t numPtsSoFar = n;
    this->m_parameterPoints.resize(n + indicesV0.size() + indicesV1.size());

    for(size_t i=0; i<indicesV0.size(); i++)
    	this->m_parameterPoints[indicesV0[i]-1] = Point2D(valuesV0[i], 0, numPtsSoFar++);

    for(size_t i=0; i<indicesV1.size(); i++)
	this->m_parameterPoints[indicesV1[i]-1] = Point2D(valuesV1[i], 1, numPtsSoFar++);

    /// Solve.
    constructAndSolveEquationSystem(neighbourhood, n, N);
}

template <class T>
void gsPeriodicParametrizationStitch<T>::constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
									 const size_t n,
									 const size_t N)
{
    std::vector<T> lambdas;
    gsMatrix<T> LHS(N, N);
    gsMatrix<T> RHS(N, 2);

    // interior points
    for (size_t i = 0; i < n; i++)
    {
        lambdas = neighbourhood.getLambdas(i);
        for (size_t j = 0; j < N; j++)
        {
            LHS(i, j) = ( i==j ? T(1) : -lambdas[j] );
	}

	for (size_t j = 0; j < m_posCorrections[i].size(); j++)
	{
	    RHS(i, 0) -= lambdas[ m_posCorrections[i][j] - 1 ];
	}
	// TODO: Actually, one can work directly with posCorrections.
	for (size_t j = 0; j < m_negCorrections[i].size(); j++)
	{
	    RHS(i, 0) += lambdas[ m_negCorrections[i][j] - 1];
	}
    }

    // points on the lower and upper boundary
    for (size_t i=n; i<N; i++)
    {
	LHS(i, i) = T(1);
	RHS.row(i) = this->m_parameterPoints[i];
    }

    // Solve the system and save the results.
    Eigen::PartialPivLU<typename gsMatrix<T>::Base> LU = LHS.partialPivLu();
    gsMatrix<T> sol = LU.solve(RHS);
    for (size_t i = 0; i < n; i++)
    {
    	this->m_parameterPoints[i] << sol(i, 0), sol(i, 1);
    }
}

template<class T>
gsMesh<T> gsPeriodicParametrizationStitch<T>::createFlatMesh(bool restrict) const
{
    // We unfold the mesh first.
    gsMesh<T> mesh;
    for (size_t i = 0; i < this->m_mesh.getNumberOfTriangles(); i++)
    {
	bool normalTriangle = false;

	bool correctors[3];
	correctors[0] = false;
	correctors[1] = false;
	correctors[2] = false;

        for (size_t j = 1; j <= 3; ++j)
        {
	    // Note: globIndex is numbered from 1 and so are the contents of m_posCorrections.
	    size_t globIndex = this->m_mesh.getGlobalVertexIndex(j, i);

	    // If we find a vertex that is neither on a stitch nor in
	    // a correction, we mark the triangle as normal.
	    if(!normalTriangle)
	    {
		// If it's not on the stitch, we look whether it is in the corrections.
		if(m_posCorrections[globIndex-1].size() == 0)
		{
		    bool isNearStitch = false;
		    for(size_t k = 0; k < m_posCorrections.size(); k++)
		    {
			for(size_t l = 0; l < m_posCorrections[k].size(); l++)
			{
			    if(m_posCorrections[k][l] == globIndex)
			    {
				isNearStitch = true;
				correctors[j-1] = true;
			    }
			}
		    }
		    if(!isNearStitch)
			normalTriangle = true;
		}
	    }
	}

	// Normal triangle means it does not intersect the domain boundary.
	if(normalTriangle)
	{
	    typename gsMesh<T>::VertexHandle v[3];
	    for (size_t j = 1; j <= 3; ++j)
	    {
		v[j - 1] = mesh.addVertex(gsParametrization<T>::getParameterPoint(this->m_mesh.getGlobalVertexIndex(j, i))[0],
					  gsParametrization<T>::getParameterPoint(this->m_mesh.getGlobalVertexIndex(j, i))[1]);
	    }
	    mesh.addFace( v[0],  v[1],  v[2]);
	}
	else
	{
	    typename gsMesh<T>::VertexHandle v[3];
	    for (size_t j = 1; j <= 3; ++j)
	    {
		if(correctors[j-1])
		{
		    v[j - 1] = mesh.addVertex(gsParametrization<T>::getParameterPoint(this->m_mesh.getGlobalVertexIndex(j, i))[0]-1,
					      gsParametrization<T>::getParameterPoint(this->m_mesh.getGlobalVertexIndex(j, i))[1]);
		}
		else
		{
		    v[j - 1] = mesh.addVertex(gsParametrization<T>::getParameterPoint(this->m_mesh.getGlobalVertexIndex(j, i))[0],
					      gsParametrization<T>::getParameterPoint(this->m_mesh.getGlobalVertexIndex(j, i))[1]);
		}
	    }
	    mesh.addFace( v[0],  v[1],  v[2]);
	}	    
    }

    gsMesh<T> unfolded = mesh.cleanMesh();
    if(restrict)
    {
	typename gsParametrization<T>::FlatMesh display(unfolded);
	return display.createRestrictedFlatMesh();
    }
    else
	return unfolded;
}

} // namespace gismo
