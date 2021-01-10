/** @file gsPeriodicParametrizationStitch.hpp

    @brief Provides implementation of the gsPeriodicParametrizationStitch class.

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
gsPeriodicParametrizationStitch<T>& gsPeriodicParametrizationStitch<T>::compute_periodic_stitch(std::string bottomFile,
												std::string topFile,
												std::string stitchFile,
												std::vector<std::vector<size_t> >& posCorrections)
{
    gsFileData<> fd_v0(bottomFile);
    gsMatrix<> pars, pts;
    fd_v0.getId<gsMatrix<> >(0, pars);
    fd_v0.getId<gsMatrix<> >(1, pts);

    GISMO_ASSERT(pars.cols() == pts.cols(), "The numbers of parameters and points of v0 differ.");
    
    std::vector<size_t> indicesV0;
    std::vector<T> valuesV0;

    for(index_t c=0; c<pts.cols(); c++)
    {
	indicesV0.push_back(this->m_mesh.findVertex(pts(0, c), pts(1, c), pts(2, c), true));
	valuesV0.push_back(pars(0, c)); // pars.cols() - c - 1 used to be here.
    }

    gsFileData<> fd_v1(topFile);
    fd_v1.getId<gsMatrix<> >(0, pars);
    fd_v1.getId<gsMatrix<> >(1, pts);

    GISMO_ASSERT(pars.cols() == pts.cols(), "The numbers of parameters and points of v1 differ.");

    std::vector<size_t> indicesV1;
    std::vector<T> valuesV1;

    for(index_t c=0; c<pts.cols(); c++)
    {
	indicesV1.push_back(this->m_mesh.findVertex(pts(0, c), pts(1, c), pts(2, c), true));
	valuesV1.push_back(pars(0, c));
    }

    // Read the stitch indices.
    gsFileData<> fd_overlap(stitchFile);
    gsMatrix<> stitchPoints;
    fd_overlap.getId<gsMatrix<> >(0, stitchPoints);
    std::vector<size_t> stitchIndices;

    for(index_t c=0; c<stitchPoints.cols(); c++)
	stitchIndices.push_back(this->m_mesh.findVertex(stitchPoints(0, c), stitchPoints(1, c), stitchPoints(2, c), true));

    // Calculation itself.
    calculate_periodic_stitch(this->m_options.getInt("parametrizationMethod"),
			      indicesV0, valuesV0, indicesV1, valuesV1, stitchIndices,
			      posCorrections);
    return *this;
}

template<class T>
void gsPeriodicParametrizationStitch<T>::calculate_periodic_stitch(const size_t paraMethod,
								   const std::vector<size_t>& indicesV0,
								   const std::vector<T>& valuesV0,
								   const std::vector<size_t>& indicesV1,
								   const std::vector<T>& valuesV1,
								   const std::vector<size_t>& stitchIndices,
								   std::vector<std::vector<size_t> >& posCorrections)
{
    typedef typename gsParametrization<T>::Point2D       Point2D ;

    size_t n = this->m_mesh.getNumberOfInnerVertices();
    size_t N = this->m_mesh.getNumberOfVertices();

    //std::vector<std::vector<size_t> > posCorrections(N);
    posCorrections = std::vector<std::vector<size_t> >(N);
    std::vector<std::vector<size_t> > negCorrections(N);
    Neighbourhood neighbourhood(this->m_mesh, stitchIndices, posCorrections, negCorrections, paraMethod);

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
    constructAndSolveEquationSystem(neighbourhood, n, N, posCorrections, negCorrections);
}

template <class T>
void gsPeriodicParametrizationStitch<T>::constructAndSolveEquationSystem(
								   const Neighbourhood &neighbourhood,
								   const size_t n,
								   const size_t N,
								   const std::vector<std::vector<size_t> >& posCorrections,
								   const std::vector<std::vector<size_t> >& negCorrections)
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

	for (size_t j = 0; j < posCorrections[i].size(); j++)
	{
	    RHS(i, 0) += lambdas[ posCorrections[i][j] - 1 ];
	}
	// TODO: Actually, one can work directly with posCorrections.
	for (size_t j = 0; j < negCorrections[i].size(); j++)
	{
	    RHS(i, 0) -= lambdas[ negCorrections[i][j] - 1];
	}
    }

    // points on the lower and upper boundary
    for (size_t i=n; i<N; i++)
    {
	LHS(i, i) = T(1);
	RHS.row(i) = this->m_parameterPoints[i];
    }

    Eigen::PartialPivLU<typename gsMatrix<T>::Base> LU = LHS.partialPivLu();
    gsMatrix<T> sol = LU.solve(RHS);
    for (size_t i = 0; i < n; i++)
    {
    	this->m_parameterPoints[i] << sol(i, 0), sol(i, 1);
    }
}

template<class T>
gsMesh<T> gsPeriodicParametrizationStitch<T>::createFlatMesh(const std::vector<std::vector<size_t> >& posCorrections,
					       bool restrict) const
{
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
	    // Note: globIndex is numbered from 1 and so are the contents of posCorrections.
	    size_t globIndex = this->m_mesh.getGlobalVertexIndex(j, i);

	    // If we find a vertex that is neither on a stitch nor in
	    // a correction, we mark the triangle as normal.
	    if(!normalTriangle)
	    {
		// If it's not on the stitch, we look whether it is in the corrections.
		if(posCorrections[globIndex-1].size() == 0)
		{
		    bool isNearStitch = false;
		    for(size_t k = 0; k < posCorrections.size(); k++)
		    {
			for(size_t l = 0; l < posCorrections[k].size(); l++)
			{
			    if(posCorrections[k][l] == globIndex)
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
		    v[j - 1] = mesh.addVertex(gsParametrization<T>::getParameterPoint(this->m_mesh.getGlobalVertexIndex(j, i))[0]+1,
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
	gsHalfEdgeMesh<T> unfoldedHEM(unfolded);
	return gsParametrization<T>::createRestrictedFlatMesh(unfoldedHEM);
    }
    else
	return unfolded;
}

} // namespace gismo
