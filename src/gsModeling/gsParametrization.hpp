/** @file gsParametrization.hpp

    @brief Provides implementation gsParametrization class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): L. Groiss, J. Vogl, D. Mokris

*/

#include <gsIO/gsOptionList.h>
#include <gsModeling/gsLineSegment.h>
#include <gismo.h>

namespace gismo
{

template<class T>
bool gsParametrization<T>::rangeCheck(const std::vector<int> &corners, const size_t minimum, const size_t maximum)
{
    for (std::vector<int>::const_iterator it = corners.begin(); it != corners.end(); it++)
    {
        if ((size_t)*it < minimum || (size_t)*it > maximum)
        { return false; }
    }
    return true;
}

template<class T>
gsOptionList gsParametrization<T>::defaultOptions()
{
    gsOptionList opt;
    opt.addInt("boundaryMethod", "boundary methodes: {1:chords, 2:corners, 3:smallest, 4:restrict, 5:opposite, 6:distributed}", 4);
    opt.addInt("parametrizationMethod", "parametrization methods: {1:shape, 2:uniform, 3:distance}", 1);
    std::vector<int> corners;
    opt.addMultiInt("corners", "vector for corners", corners);
    opt.addReal("range", "in case of restrict or opposite", 0.1);
    opt.addInt("number", "number of corners, in case of corners", 4);
    opt.addReal("precision", "precision to calculate", 1E-8);
    return opt;
}

template<class T>
gsParametrization<T>::gsParametrization(gsMesh<T> &mesh,
					const gsOptionList & list,
					bool periodic) : m_mesh(mesh, 1E-12, periodic)
{
    m_options.update(list, gsOptionList::addIfUnknown);
}

template<class T>
void gsParametrization<T>::calculate(const size_t boundaryMethod,
                                     const size_t paraMethod,
                                     const std::vector<int>& cornersInput,
                                     const T rangeInput,
                                     const size_t numberInput)
{
    GISMO_ASSERT(boundaryMethod >= 1 && boundaryMethod <= 6,
                 "The boundary method " << boundaryMethod << " is not valid.");
    GISMO_ASSERT(paraMethod >= 1 && paraMethod <= 3, "The parametrization method " << paraMethod << " is not valid.");
    size_t n = m_mesh.getNumberOfInnerVertices();
    size_t N = m_mesh.getNumberOfVertices();
    size_t B = m_mesh.getNumberOfBoundaryVertices();
    Neighbourhood neighbourhood(m_mesh, paraMethod);

    T w = 0;
    std::vector<T> halfedgeLengths = m_mesh.getBoundaryChordLengths();
    std::vector<int> corners;
    std::vector<T> lengths;

    switch (boundaryMethod)
    {
        case 1:
            m_parameterPoints.reserve(n + B);
            for (size_t i = 1; i <= n + 1; i++)
            {
                m_parameterPoints.push_back(Point2D(0, 0, i));
            }
            for (size_t i = 0; i < B - 1; i++)
            {
                w += halfedgeLengths[i] * (1. / m_mesh.getBoundaryLength()) * 4;
                m_parameterPoints.push_back(Neighbourhood::findPointOnBoundary(w, n + i + 2));
            }
            break;
        case 2:
            corners = cornersInput;
        case 3:
        case 4:
        case 5:
        case 6: // N
            if (boundaryMethod != 2)
                corners = neighbourhood.getBoundaryCorners(boundaryMethod, rangeInput, numberInput);

            m_parameterPoints.reserve(N);
            for (size_t i = 1; i <= N; i++)
            {
                m_parameterPoints.push_back(Point2D(0, 0, i));
            }

            lengths = m_mesh.getCornerLengths(corners);
            m_parameterPoints[n + corners[0] - 1] = Point2D(0, 0, n + corners[0]);

            for (size_t i = corners[0] + 1; i < corners[0] + B; i++)
            {
                w += halfedgeLengths[(i - 2) % B]
                    / findLengthOfPositionPart(i > B ? i - B : i, B, corners, lengths);
                m_parameterPoints[(n + i - 1) > N - 1 ? n + i - 1 - B : n + i - 1] =
                    Neighbourhood::findPointOnBoundary(w, n + i > N ? n + i - B : n + i);
            }
            break;
        default:
            GISMO_ERROR("boundaryMethod not valid: " << boundaryMethod);
    }

    constructAndSolveEquationSystem_2(neighbourhood, n, N);
}


template<class T>
void gsParametrization<T>::constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
                                                           const size_t n,
                                                           const size_t N)
{
    gsMatrix<T> A;
    A.resize(n, n);
    std::vector<T> lambdas;
    gsVector<T> b1(n), b2(n);
    b1.setZero(); b2.setZero();

    for (size_t i = 0; i < n; i++)
    {
        lambdas = neighbourhood.getLambdas(i);
        for (size_t j = 0; j < n; j++)
        {
            A(i, j) = ( i==j ? T(1) : -lambdas[j] );
        }

        for (size_t j = n; j < N; j++)
        {
            b1(i) += (lambdas[j]) * (m_parameterPoints[j][0]);
            b2(i) += (lambdas[j]) * (m_parameterPoints[j][1]);
        }
    }

    gsVector<T> u(n), v(n);
    Eigen::PartialPivLU<typename gsMatrix<T>::Base> LU = A.partialPivLu();
    u = LU.solve(b1);
    v = LU.solve(b2);

    for (size_t i = 0; i < n; i++)
        m_parameterPoints[i] << u(i), v(i);
}

template <class T>
void gsParametrization<T>::constructAndSolveEquationSystem_2(const Neighbourhood &neighbourhood,
							     const size_t n,
							     const size_t N)
{
    gsMatrix<T> LHS(N,N);
    gsMatrix<T> RHS(N,2);
    std::vector<T> lambdas;

    for (size_t i = 0; i < n; i++)
    {
        lambdas = neighbourhood.getLambdas(i);
        for (size_t j = 0; j < N; j++)
        {
	    // Standard way:
            // LHS(i, j) = ( i==j ? T(1) : -lambdas[j] );
	    LHS(i, j) = lambdas[j];
	    // Initial guess:
	    RHS(i, 0) = 0.5;
	    RHS(i, 1) = 0.5;
        }
    }

    for (size_t i=n; i<N; i++)
    {
	LHS(i,i) = T(1);
	RHS.row(i) = m_parameterPoints[i];
    }

    gsMatrix<T> sol;
    // Eigen::PartialPivLU<typename gsMatrix<T>::Base> LU = LHS.partialPivLu();
    // sol = LU.solve(RHS);

    for(size_t k=0; k<=100; k++)
    {
	sol = LHS * RHS;
	RHS = sol;

	for (size_t i = 0; i < n; i++)
	    m_parameterPoints[i] << sol(i, 0), sol(i, 1);

	if(k%5 == 0)
	{
	    const gsMesh<T> mesh = createFlatMesh();
	    gsWriteParaview(mesh, "mesh" + std::to_string(k));
	}
    }
}

template<class T>
const typename gsParametrization<T>::Point2D &gsParametrization<T>::getParameterPoint(size_t vertexIndex) const
{
    return m_parameterPoints[vertexIndex - 1];
}

template<class T>
gsMatrix<T> gsParametrization<T>::createUVmatrix()
{
    gsMatrix<T> m(2, m_mesh.getNumberOfVertices());
    for (size_t i = 1; i <= m_mesh.getNumberOfVertices(); i++)
    {
        m.col(i - 1) << this->getParameterPoint(i)[0], this->getParameterPoint(i)[1];
    }
    return m;
}

template<class T>
gsMatrix<T> gsParametrization<T>::createXYZmatrix()
{
    gsMatrix<T> m(3, m_mesh.getNumberOfVertices());
    for (size_t i = 1; i <= m_mesh.getNumberOfVertices(); i++)
    {
        m.col(i - 1) << m_mesh.getVertex(i)->x(), m_mesh.getVertex(i)->y(), m_mesh.getVertex(i)->z();
    }
    return m;
}

template<class T>
gsMesh<T> gsParametrization<T>::createFlatMesh() const
{
    gsMesh<T> mesh;
    mesh.reserve(3 * m_mesh.getNumberOfTriangles(), m_mesh.getNumberOfTriangles(), 0);
    for (size_t i = 0; i < m_mesh.getNumberOfTriangles(); i++)
    {
        typename gsMesh<T>::VertexHandle v[3];
        for (size_t j = 1; j <= 3; ++j)
        {
            v[j - 1] = mesh.addVertex(getParameterPoint(m_mesh.getGlobalVertexIndex(j, i))[0],
				      getParameterPoint(m_mesh.getGlobalVertexIndex(j, i))[1]);
        }
	 mesh.addFace(v[0], v[1], v[2]);
    }
    return mesh.cleanMesh();
}

template<class T>
gsMesh<T> gsParametrization<T>::createFlatMesh(const std::vector<std::vector<size_t> >& posCorrections,
					       bool restrict) const
{
    gsMesh<T> mesh;
    for (size_t i = 0; i < m_mesh.getNumberOfTriangles(); i++)
    {
	bool normalTriangle = false;

	bool correctors[3];
	correctors[0] = false;
	correctors[1] = false;
	correctors[2] = false;

        for (size_t j = 1; j <= 3; ++j)
        {
	    // Note: globIndex is numbered from 1 and so are the contents of posCorrections.
	    size_t globIndex = m_mesh.getGlobalVertexIndex(j, i);

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
		v[j - 1] = mesh.addVertex(getParameterPoint(m_mesh.getGlobalVertexIndex(j, i))[0],
					  getParameterPoint(m_mesh.getGlobalVertexIndex(j, i))[1]);
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
		    v[j - 1] = mesh.addVertex(getParameterPoint(m_mesh.getGlobalVertexIndex(j, i))[0]+1,
					      getParameterPoint(m_mesh.getGlobalVertexIndex(j, i))[1]);
		}
		else
		{
		    v[j - 1] = mesh.addVertex(getParameterPoint(m_mesh.getGlobalVertexIndex(j, i))[0],
					      getParameterPoint(m_mesh.getGlobalVertexIndex(j, i))[1]);
		}
	    }
	    mesh.addFace( v[0],  v[1],  v[2]);
	}	    
    }

    gsMesh<T> unfolded = mesh.cleanMesh();
    if(restrict)
    {
	gsHalfEdgeMesh<T> unfoldedHEM(unfolded);
	return createRestrictedFlatMesh(unfoldedHEM);
    }
    else
	return unfolded;
}

template<class T>
gsMesh<T> gsParametrization<T>::createShiftedCopy(const gsMesh<T>& original, real_t uShift) const
{
    gsHalfEdgeMesh<T> originalHEM(original);
    gsMesh<T> result;
    for(size_t i=0; i<originalHEM.getNumberOfTriangles(); i++)
    {
	typename gsMesh<T>::VertexHandle vh[3];
	for(size_t j=1; j<=3; j++)
	{
	    vh[j-1] = originalHEM.getVertex(originalHEM.getGlobalVertexIndex(j ,i));
	}
	result.addFace(
	    result.addVertex(vh[0]->x() + uShift, vh[0]->y()),
	    result.addVertex(vh[1]->x() + uShift, vh[1]->y()),
	    result.addVertex(vh[2]->x() + uShift, vh[2]->y()));
    }
    return result;
}

template<class T>
real_t gsParametrization<T>::correspondingV(const typename gsMesh<T>::VertexHandle& h0,
					    const typename gsMesh<T>::VertexHandle& h1,
					    real_t u) const
{
    real_t u0 = (*h0)[0];
    real_t u1 = (*h1)[0];
    real_t v0 = (*h0)[1];
    real_t v1 = (*h1)[1];

    real_t t = (u - u0) / (u1 - u0);

    return (1 - t) * v0 + t * v1;
}

// v1 is outside the domain, v0 and v2 inside.
template<class T>
void gsParametrization<T>::addThreeFlatTrianglesOneOut(gsMesh<T>& mesh,
						       const typename gsMesh<T>::VertexHandle& v0,
						       const typename gsMesh<T>::VertexHandle& v1,
						       const typename gsMesh<T>::VertexHandle& v2) const
{
    // Note: v are in the input mesh, w in the output.

    typename gsMesh<T>::VertexHandle w0 = mesh.addVertex(v0->x(), v0->y());
    typename gsMesh<T>::VertexHandle w2 = mesh.addVertex(v2->x(), v2->y());

    if(v1->x() < 0)
    {
	// Two triangles on the left.
	typename gsMesh<T>::VertexHandle w01 = mesh.addVertex(0, correspondingV(v0, v1, 0));
	typename gsMesh<T>::VertexHandle w12 = mesh.addVertex(0, correspondingV(v1, v2, 0));

	mesh.addFace(w0, w01, w12);
	mesh.addFace(w0, w12, w2);

	// One triangle on the right.
	typename gsMesh<T>::VertexHandle vvv01 = mesh.addVertex(1, correspondingV(v0, v1, 0));
	typename gsMesh<T>::VertexHandle vvv12 = mesh.addVertex(1, correspondingV(v1, v2, 0));
	typename gsMesh<T>::VertexHandle v1copy = mesh.addVertex(v1->x() + 1, v1->y());
	mesh.addFace(vvv01, v1copy, vvv12);	
    }
    else if(v1->x() > 1)
    {
	// Two triangles on the left.
	typename gsMesh<T>::VertexHandle w01 = mesh.addVertex(1, correspondingV(v0, v1, 1));
	typename gsMesh<T>::VertexHandle w12 = mesh.addVertex(1, correspondingV(v1, v2, 1));

	mesh.addFace(w0, w01, w12);
	mesh.addFace(w0, w12, w2);

	// One triangle on the right.
	typename gsMesh<T>::VertexHandle vvv01 = mesh.addVertex(0, correspondingV(v0, v1, 1));
	typename gsMesh<T>::VertexHandle vvv12 = mesh.addVertex(0, correspondingV(v1, v2, 1));
	typename gsMesh<T>::VertexHandle v1copy = mesh.addVertex(v1->x() - 1, v1->y());
	mesh.addFace(vvv01, v1copy, vvv12);
    }
    else
	gsWarn << "This situation of addThreeFlatTriangles should not happen, v1->x() = "
	       << v1->x() << "." << std::endl;
}

// v1 is inside the domain, v0 and v2 outside.
template<class T>
void gsParametrization<T>::addThreeFlatTrianglesTwoOut(gsMesh<T>& mesh,
						       const typename gsMesh<T>::VertexHandle& v0,
						       const typename gsMesh<T>::VertexHandle& v1,
						       const typename gsMesh<T>::VertexHandle& v2) const
{
    if(v0->x() < 0 && v2->x() < 0)
    {
	typename gsMesh<T>::VertexHandle w0 = mesh.addVertex(v0->x() + 1, v0->y());
	typename gsMesh<T>::VertexHandle w1 = mesh.addVertex(v1->x() + 1, v1->y());
	typename gsMesh<T>::VertexHandle w2 = mesh.addVertex(v2->x() + 1, v2->y());
	addThreeFlatTrianglesOneOut(mesh, w0, w1, w2);
    }
    else if(v0->x() > 1 && v2->x() > 1)
    {
	typename gsMesh<T>::VertexHandle w0 = mesh.addVertex(v0->x() - 1, v0->y());
	typename gsMesh<T>::VertexHandle w1 = mesh.addVertex(v1->x() - 1, v1->y());
	typename gsMesh<T>::VertexHandle w2 = mesh.addVertex(v2->x() - 1, v2->y());
	addThreeFlatTrianglesOneOut(mesh, w0, w1, w2);
    }
    else
	gsWarn << "This situation of addThreeFlatTrianglesTwoOut should not happen, v1->x()="
	       << v1->x() << "." << std::endl;
}

template<class T>
void gsParametrization<T>::addOneFlatTriangleNotIntersectingBoundary(gsMesh<T>& mesh,
								     typename gsMesh<T>::VertexHandle& v0,
								     typename gsMesh<T>::VertexHandle& v1,
								     typename gsMesh<T>::VertexHandle& v2) const
{
    // Note: I wanted to solve this by modifying the x-coordinates of
    // the vertex handles and recursion. However, this creates mess,
    // as the vertex handles are shared among several triangles.
    real_t v0x = v0->x();
    real_t v1x = v1->x();
    real_t v2x = v2->x();

    while(v0x > 1 && v1x > 1 && v2x > 1)
    {
	v0x -= 1;
	v1x -= 1;
	v2x -= 1;
    }

    while(v0x < 0 && v1x < 0 && v2x < 0)
    {
	v0x += 1;
	v1x += 1;
	v2x += 1;
    }

    if(v0x >= 0 && v0x <= 1 &&
       v1x >= 0 && v1x <= 1 &&
       v2x >= 0 && v2x <= 1)
    {
	mesh.addFace(
	    mesh.addVertex(v0x, v0->y()),
	    mesh.addVertex(v1x, v1->y()),
	    mesh.addVertex(v2x, v2->y()));
    }
    else
    {
	gsWarn << "This triangle does intersect the boundary.";
	gsWarn << "v0: " << v0x << ", " << v0->y() << std::endl;
	gsWarn << "v1: " << v1x << ", " << v1->y() << std::endl;
	gsWarn << "v2: " << v2x << ", " << v2->y() << std::endl;
    }
}

template<class T>
gsMesh<T> gsParametrization<T>::createMidMesh(const std::vector<size_t>& right,
					      const std::vector<size_t>& left) const
{
    gsMesh<T> midMesh;
    midMesh.reserve(3 * m_mesh.getNumberOfTriangles(), m_mesh.getNumberOfTriangles(), 0);

    for (size_t i = 0; i < m_mesh.getNumberOfTriangles(); i++)
    {
	size_t vInd[3];
	std::vector<size_t> lVert, rVert;
        for (size_t j = 1; j <= 3; ++j)
	{
	    vInd[j-1] = m_mesh.getGlobalVertexIndex(j, i);
	    if(std::find(right.begin(), right.end(), vInd[j-1]) != right.end())
		rVert.push_back(j-1);
	    if(std::find(left.begin(),  left.end(),  vInd[j-1]) != left.end())
		lVert.push_back(j-1);
	}

	if(lVert.size() > 0 && rVert.size() > 0 && lVert.size() + rVert.size() == 3)
	{
	    // Make two copies
	    typename gsMesh<T>::VertexHandle mvLft[3], mvRgt[3];
	    
	    for (size_t j=0; j<3; ++j)
	    {
		if(std::find(rVert.begin(), rVert.end(), j) != rVert.end())
		{
		    mvLft[j] = midMesh.addVertex(getParameterPoint(vInd[j])[0],
		    				 getParameterPoint(vInd[j])[1]);
		    mvRgt[j] = midMesh.addVertex(getParameterPoint(vInd[j])[0]+1,
		    				 getParameterPoint(vInd[j])[1]);
		}
		else
		{
		    mvLft[j] = midMesh.addVertex(getParameterPoint(vInd[j])[0]-1,
						getParameterPoint(vInd[j])[1]);
		    mvRgt[j] = midMesh.addVertex(getParameterPoint(vInd[j])[0],
						 getParameterPoint(vInd[j])[1]);
		}
	    }
	    midMesh.addFace(mvLft[0], mvLft[1], mvLft[2]);
	    midMesh.addFace(mvRgt[0], mvRgt[1], mvRgt[2]);
	}
	else
	{
	    // Make just one triangle
	    typename gsMesh<T>::VertexHandle mv[3];
	    for (size_t j=0; j<3; ++j)
	    {
		mv[j] = midMesh.addVertex(getParameterPoint(vInd[j])[0],
					  getParameterPoint(vInd[j])[1]);
	    }
	    midMesh.addFace(mv[0], mv[1], mv[2]);
	}
    }
    return midMesh.cleanMesh();
}

template<class T>
gsMesh<T> gsParametrization<T>::createFlatMesh(const std::vector<size_t>& left,
					       const std::vector<size_t>& right,
					       bool restrict) const
{
    if(!restrict)
	return createFlatMesh();

    gsHalfEdgeMesh<T> unfolded(createMidMesh(left, right));
    return createRestrictedFlatMesh(unfolded);
}

template<class T>
gsMesh<T> gsParametrization<T>::createRestrictedFlatMesh(const gsHalfEdgeMesh<T>& unfolded) const
{
    gsMesh<T> result;

    for(size_t i=0; i<unfolded.getNumberOfTriangles(); i++)
    {
	// Remember the corners and which of them are inside the domain.
	bool out[3];
	typename gsMesh<T>::VertexHandle vh[3];
	for(size_t j=1; j<=3; ++j)
	{
	    vh[j-1] = unfolded.getVertex(unfolded.getGlobalVertexIndex(j, i));
	    real_t u = vh[j-1]->x();

	    if(u < 0 || u > 1)
		out[j-1] = true;
	    else
		out[j-1] = false;
	}
	if( !out[0] && !out[1] && !out[2] )
	    addOneFlatTriangleNotIntersectingBoundary(result, vh[0], vh[1], vh[2]);

	else if( out[0] && !out[1] && out[2] )
	    addThreeFlatTrianglesTwoOut(result, vh[0], vh[1], vh[2]);
	else if( out[0] && out[1] && !out[2] )
	    addThreeFlatTrianglesTwoOut(result, vh[1], vh[2], vh[0]);
	else if( !out[0] && out[1] && out[2] )
	    addThreeFlatTrianglesTwoOut(result, vh[2], vh[0], vh[1]);

	else if( !out[0] && !out[1] && out[2] )
	    addThreeFlatTrianglesOneOut(result, vh[1], vh[2], vh[0]);
	else if( !out[0] && out[1] && !out[2] )
	    addThreeFlatTrianglesOneOut(result, vh[0], vh[1], vh[2]);
	else if( out[0] && !out[1] && !out[2] )
	    addThreeFlatTrianglesOneOut(result, vh[2], vh[0], vh[1]);

	else
	    addOneFlatTriangleNotIntersectingBoundary(result, vh[0], vh[1], vh[2]);
    }
    return result.cleanMesh();
}

template <class T>
void gsParametrization<T>::writeTexturedMesh(std::string filename) const
{
    gsMatrix<T> params(m_mesh.numVertices(), 2);

    for(size_t i=0; i<m_mesh.numVertices(); i++)
    {
	size_t index = m_mesh.unsorted(i);
	params.row(i) = getParameterPoint(index);
    }
    gsWriteParaview(m_mesh, "mesh", params);
}

template <class T>
void gsParametrization<T>::writeSTL(const gsMesh<T>& mesh, std::string filename) const
{
    std::string mfn(filename);
    mfn.append(".stl");
    std::ofstream file(mfn.c_str());

    gsHalfEdgeMesh<T> hMesh(mesh);

    if(!file.is_open())
	gsWarn << "Opening file " << mfn << " for writing failed." << std::endl;

    file << std::fixed;
    file << std::setprecision(12);

    file << "solid created by G+Smo" << std::endl;
    for(size_t t=0; t<hMesh.getNumberOfTriangles(); t++)
    {
	file << " facet normal 0 0 -1" << std::endl;
	file << "  outer loop" << std::endl;
	for(size_t v=0; v<3; v++)
	{
	    typename gsMesh<T>::VertexHandle handle = hMesh.getVertex(hMesh.getGlobalVertexIndex(v+1, t));
	    file << "   vertex " << handle->y() << " " << handle->x() << " " << handle->z() << std::endl;
	}
	file << "  endloop" << std::endl;
	file << " endfacet" << std::endl;
    }
    file << "endsolid" << std::endl;
}

template<class T>
gsParametrization<T>& gsParametrization<T>::setOptions(const gsOptionList& list)
{
    m_options.update(list, gsOptionList::addIfUnknown);
    return *this;
}

template<class T>
gsParametrization<T>& gsParametrization<T>::compute()
{
    calculate(m_options.getInt("boundaryMethod"),
              m_options.getInt("parametrizationMethod"),
              m_options.getMultiInt("corners"),
              m_options.getReal("range"),
              m_options.getInt("number"));

    return *this;
}

template<class T>
gsParametrization<T>& gsParametrization<T>::compute_free_boundary()
{
    calculate_free_boundary(m_options.getInt("parametrizationMethod"),
			    m_options.getString("fileCorners"));

    return *this;
}

template<class T>
T gsParametrization<T>::findLengthOfPositionPart(const size_t position,
                                                      const size_t numberOfPositions,
                                                      const std::vector<int> &bounds,
                                                      const std::vector<T> &lengths)
{
    GISMO_UNUSED(numberOfPositions);
    GISMO_ASSERT(1 <= position && position <= numberOfPositions, "The position " << position
                 << " is not a valid input. There are only " << numberOfPositions << " possible positions.");
    GISMO_ASSERT(rangeCheck(bounds, 1, numberOfPositions), "The bounds are not a valid input. They have to be out of the possible positions, which only are "
                 << numberOfPositions << ". ");
    size_t numberOfBounds = bounds.size();
    size_t s = lengths.size();
    if (position > (size_t)bounds[numberOfBounds - 1] || position <= (size_t)bounds[0])
        return lengths[s - 1];
    for (size_t i = 0; i < numberOfBounds; i++)
    {
        if (position - (size_t)bounds[0] + 1 > (size_t)bounds[i] - (size_t)bounds[0] + 1
            && position - (size_t)bounds[0] + 1 <= (size_t)bounds[(i + 1) % numberOfBounds] - (size_t)bounds[0] + 1)
            return lengths[i];
    }
    return 0;
}


//******************************************************************************************
//******************************* nested class Neighbourhood *******************************
//******************************************************************************************

template<class T>
gsParametrization<T>::Neighbourhood::Neighbourhood(const gsHalfEdgeMesh<T> & meshInfo,
						   const size_t parametrizationMethod)
    : m_basicInfos(meshInfo)
{
    m_localParametrizations.reserve(meshInfo.getNumberOfInnerVertices());
    for(size_t i=1; i <= meshInfo.getNumberOfInnerVertices(); i++)
    {
        m_localParametrizations.push_back(LocalParametrization(meshInfo, LocalNeighbourhood(meshInfo, i),
							       parametrizationMethod));
    }

    m_localBoundaryNeighbourhoods.reserve(meshInfo.getNumberOfVertices() - meshInfo.getNumberOfInnerVertices());
    for(size_t i=meshInfo.getNumberOfInnerVertices()+1; i<= meshInfo.getNumberOfVertices(); i++)
    {

        m_localBoundaryNeighbourhoods.push_back(LocalNeighbourhood(meshInfo, i, 0));
    }
}

template<class T>
std::vector<size_t> gsParametrization<T>::Neighbourhood::computeCorrections(
    const std::vector<size_t>& stitchIndices,
    const LocalNeighbourhood& localNeighbourhood) const
{
    //gsInfo << "vertex index: " << localNeighbourhood.getVertexIndex() << std::endl;
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
gsParametrization<T>::Neighbourhood::Neighbourhood(const gsHalfEdgeMesh<T> & meshInfo,
						   const std::vector<size_t>& stitchIndices,
						   std::vector<std::vector<size_t> >& posCorrections,
						   std::vector<std::vector<size_t> >& negCorrections,
						   const size_t parametrizationMethod)
    : m_basicInfos(meshInfo)
{
    // TODO: Now we have posCorrections for those on the stitch interacting with those to the left.
    // We also need the posCorrections for those to the left interacting with the stitch.
    GISMO_ASSERT(posCorrections.size() == meshInfo.getNumberOfVertices(), "posCorrections not properly initialized.");
    m_localParametrizations.reserve(meshInfo.getNumberOfInnerVertices());
    gsInfo << "Positive correction\n";
    for(size_t i=1; i <= meshInfo.getNumberOfInnerVertices(); i++)
    {
	LocalNeighbourhood localNeighbourhood(meshInfo, i);

        m_localParametrizations.push_back(LocalParametrization(meshInfo, localNeighbourhood,
							       parametrizationMethod));

	posCorrections[i-1] = computeCorrections(stitchIndices, localNeighbourhood);
	for(size_t j=0; j < posCorrections[i-1].size(); j++)
	{
	    negCorrections[posCorrections[i-1][j]-1].push_back(i);
	}
    }

    m_localBoundaryNeighbourhoods.reserve(meshInfo.getNumberOfVertices() - meshInfo.getNumberOfInnerVertices());
    for(size_t i=meshInfo.getNumberOfInnerVertices()+1; i<= meshInfo.getNumberOfVertices(); i++)
    {
	LocalNeighbourhood localNeighbourhood(meshInfo, i, 0);

        m_localBoundaryNeighbourhoods.push_back(localNeighbourhood);
	posCorrections[i-1] = computeCorrections(stitchIndices, localNeighbourhood);
	for(size_t j=0; j < posCorrections[i-1].size(); j++)
	{
	    negCorrections[posCorrections[i-1][j]-1].push_back(i);
	}
    }
}

template<class T>
const std::vector<T>& gsParametrization<T>::Neighbourhood::getLambdas(const size_t i) const
{
    return m_localParametrizations[i].getLambdas();
}

template<class T>
const std::vector<int> gsParametrization<T>::Neighbourhood::getBoundaryCorners(const size_t method, const T range, const size_t number) const
{
    std::vector<std::pair<T , size_t> > angles;
    std::vector<int> corners;
    angles.reserve(m_localBoundaryNeighbourhoods.size());
    for(typename std::vector<LocalNeighbourhood>::const_iterator it=m_localBoundaryNeighbourhoods.begin(); it!=m_localBoundaryNeighbourhoods.end(); it++)
    {
        angles.push_back(std::pair<T , size_t>(it->getInnerAngle(), it->getVertexIndex() - m_basicInfos.getNumberOfInnerVertices()));
    }
    std::sort(angles.begin(), angles.end());
    if(method == 3)
    {
        this->takeCornersWithSmallestAngles(4, angles, corners);
        std::sort(corners.begin(), corners.end());
        gsDebug << "According to the method 'smallest inner angles' the following corners were chosen:\n";
        for(std::vector<int>::iterator it=corners.begin(); it!=corners.end(); it++)
        {
            gsDebug << (*it) << "\n";
        }
    }
    else if(method == 5)
    {
        searchAreas(range, angles, corners);
        gsDebug << "According to the method 'nearly opposite corners' the following corners were chosen:\n";
        for(size_t i=0; i<corners.size(); i++)
        {
            gsDebug << corners[i] << "\n";
        }
    }
    else if(method == 4)
    {
        bool flag = true;
        corners.reserve(4);
        corners.push_back(angles.front().second);
        angles.erase(angles.begin());
        while(corners.size() < 4)
        {
            flag = true;
            for(std::vector<int>::iterator it=corners.begin(); it!=corners.end(); it++)
            {
                if(m_basicInfos.getShortestBoundaryDistanceBetween(angles.front().second, *it) < range*m_basicInfos.getBoundaryLength())
                    flag = false;
            }
            if(flag)
                corners.push_back(angles.front().second);
            angles.erase(angles.begin());
        }
        std::sort(corners.begin(), corners.end());
        for(size_t i=0; i<corners.size(); i++)
        {
            gsDebug << corners[i] << "\n";
        }
    }
    else if(method == 6)
    {
        T oldDifference = 0;
        T newDifference = 0;
        std::vector<int> newCorners;
        std::vector<T> lengths;
        angles.erase(angles.begin()+number, angles.end());
        gsDebug << "Angles:\n";
        for(size_t i=0; i<angles.size(); i++)
        {
            gsDebug << angles[i].first << ", " << angles[i].second << "\n";
        }
        newCorners.reserve((angles.size()*(angles.size()-1)*(angles.size()-2)*(angles.size()-3))/6);
        corners.reserve((angles.size()*(angles.size()-1)*(angles.size()-2)*(angles.size()-3))/6);
        for(size_t i=0; i<angles.size(); i++)   // n
        {
            for(size_t j=i+1; j<angles.size(); j++) // * (n-1)/2
            {
                for(size_t k=j+1; k<angles.size(); k++) // * (n-2)/3
                {
                    for(size_t l=k+1; l<angles.size(); l++) // * (n-3)/4
                    {
                        newCorners.push_back(angles[i].second);
                        newCorners.push_back(angles[j].second);
                        newCorners.push_back(angles[k].second);
                        newCorners.push_back(angles[l].second);
                        std::sort(newCorners.begin(), newCorners.end());
                        lengths = m_basicInfos.getCornerLengths(newCorners);
                        std::sort(lengths.begin(), lengths.end());
                        newDifference = math::abs(lengths[0] - lengths[3]);
                        if(oldDifference == 0 || newDifference < oldDifference)
                        {
                            corners.erase(corners.begin(), corners.end());
                            corners.push_back(angles[i].second);
                            corners.push_back(angles[j].second);
                            corners.push_back(angles[k].second);
                            corners.push_back(angles[l].second);
                            std::sort(corners.begin(), corners.end());
                        }
                        newCorners.erase(newCorners.begin(), newCorners.end());
                        oldDifference = newDifference;
                    }
                }
            }
        }
        gsDebug << "According to the method 'evenly distributed corners' the following corners were chosen:\n";
        for(size_t i=0; i<corners.size(); i++)
        {
            gsDebug << corners[i] << "\n";
        }
    }
    return corners;
}

template<class T>
const typename gsParametrization<T>::Point2D gsParametrization<T>::Neighbourhood::findPointOnBoundary(const T w, size_t vertexIndex)
{
    GISMO_ASSERT(0 <= w && w <= 4, "Wrong value for w.");
    if(0 <= w && w <=1)
        return Point2D(w,0, vertexIndex);
    else if(1<w && w<=2)
        return Point2D(1,w-1, vertexIndex);
    else if(2<w && w<=3)
        return Point2D(1-w+2,1, vertexIndex);
    else if(3<w && w<=4)
        return Point2D(0,1-w+3, vertexIndex);
    return Point2D();
}

//*****************************************************************************************************
//*****************************************************************************************************
//*******************THE******INTERN******FUNCTIONS******ARE******NOW******FOLLOWING*******************
//*****************************************************************************************************
//*****************************************************************************************************

template<class T>
void gsParametrization<T>::Neighbourhood::takeCornersWithSmallestAngles(size_t number, std::vector<std::pair<T , size_t> >& sortedAngles, std::vector<int>& corners) const
{
    sortedAngles.erase(sortedAngles.begin()+number, sortedAngles.end());

    corners.clear();
    corners.reserve(sortedAngles.size());
    for(typename std::vector<std::pair<T, size_t> >::iterator it=sortedAngles.begin(); it!=sortedAngles.end(); it++)
        corners.push_back(it->second);
}

template<class T>
std::vector<T> gsParametrization<T>::Neighbourhood::midpoints(const size_t numberOfCorners, const T length) const
{
    std::vector<T> midpoints;
    midpoints.reserve(numberOfCorners-1);
    T n = 1./numberOfCorners;
    for(size_t i=1; i<numberOfCorners; i++)
    {
        midpoints.push_back(i*length*n);
    }
    return midpoints;
}

template<class T>
void gsParametrization<T>::Neighbourhood::searchAreas(const T range, std::vector<std::pair<T, size_t> >& sortedAngles, std::vector<int>& corners) const
{
    T l = m_basicInfos.getBoundaryLength();
    std::vector<T> h = m_basicInfos.getBoundaryChordLengths();
    this->takeCornersWithSmallestAngles(1,sortedAngles, corners);
    std::vector<std::vector<std::pair<T , size_t> > > areas;
    areas.reserve(3);
    for(size_t i=0; i<3; i++)
    {
        areas.push_back(std::vector<std::pair<T , size_t> >());
    }
    std::vector<T> midpoints = this->midpoints(4, l);

    T walkAlong = 0;
    for(size_t i=0; i<h.size(); i++)
    {
        walkAlong += h[(corners[0]+i-1)%h.size()];
        for(int j = 2; j>=0; j--)
        {
            if(math::abs(walkAlong-midpoints[j]) <= l*range)
            {
                areas[j].push_back(std::pair<T , size_t>(m_localBoundaryNeighbourhoods[(corners[0]+i)%(h.size())].getInnerAngle(), (corners[0]+i)%h.size() + 1));
                break;
            }
        }
    }
    std::sort(areas[0].begin(), areas[0].end());
    std::sort(areas[1].begin(), areas[1].end());
    std::sort(areas[2].begin(), areas[2].end());
    bool smaller = false;
    //corners.reserve(3);
    for(size_t i=0; i<areas[0].size(); i++)
    {
        if(areas[0][i].second > (size_t)corners[0] || areas[0][i].second < (size_t)corners[0])
        {
            corners.push_back(areas[0][i].second);
            if(areas[0][i].second < (size_t)corners[0])
            {
                smaller = true;
            }
            break;
        }
    }
    for(size_t i=0; i<areas[1].size(); i++)
    {
        if(smaller)
        {
            if(areas[1][i].second > (size_t)corners[1] && areas[1][i].second < (size_t)corners[0])
            {
                corners.push_back(areas[1][i].second);
                break;
            }
        }
        else if(areas[1][i].second > (size_t)corners[1] || areas[1][i].second < (size_t)corners[0])
        {
            corners.push_back(areas[1][i].second);
            if(areas[1][i].second < (size_t)corners[0])
            {
                smaller = true;
            }
            break;
        }
    }
    for(size_t i=0; i<areas[2].size(); i++)
    {
        if(smaller)
        {
            if(areas[2][i].second > (size_t)corners[2] && areas[2][i].second < (size_t)corners[0])
            {
                corners.push_back(areas[2][i].second);
                break;
            }
        }
        else if(areas[2][i].second > (size_t)corners[2] || areas[2][i].second < (size_t)corners[0])
        {
            corners.push_back(areas[2][i].second);
            break;
        }
    }
}

//*******************************************************************************************
//**************************** nested class LocalParametrization ****************************
//*******************************************************************************************

template<class T>
gsParametrization<T>::LocalParametrization::LocalParametrization(const gsHalfEdgeMesh<T>& meshInfo,
								 const LocalNeighbourhood& localNeighbourhood,
								 const size_t parametrizationMethod)
{
    m_vertexIndex = localNeighbourhood.getVertexIndex();
    std::list<size_t> indices = localNeighbourhood.getVertexIndicesOfNeighbours();
    size_t d = localNeighbourhood.getNumberOfNeighbours();
    switch (parametrizationMethod)
    {
        case 1:
        {
            std::list<T> angles = localNeighbourhood.getAngles();
            VectorType points;
            T theta = 0;
            T nextAngle = 0;
            for(typename std::list<T>::iterator it = angles.begin(); it!=angles.end(); ++it)
            {
                theta += *it;
            }
            Point2D p(0, 0, 0);
            T length = (*meshInfo.getVertex(indices.front()) - *meshInfo.getVertex(m_vertexIndex)).norm();
            Point2D nextPoint(length, 0, indices.front());
            points.reserve(indices.size());
            points.push_back(nextPoint);
            gsVector<T> actualVector = nextPoint - p;
            indices.pop_front();
            T thetaInv = 1./theta;
            gsVector<T> nextVector;
            while(!indices.empty())
            {
                length = (*meshInfo.getVertex(indices.front()) - *meshInfo.getVertex(m_vertexIndex)).norm();
                //length =  (meshInfo.getVertex(indices.front()) - meshInfo.getVertex(m_vertexIndex) ).norm();
                nextAngle = angles.front()*thetaInv * 2 * EIGEN_PI;
                nextVector = (Eigen::Rotation2D<T>(nextAngle).operator*(actualVector).normalized()*length) + p;
                nextPoint = Point2D(nextVector[0], nextVector[1], indices.front());
                points.push_back(nextPoint);
                actualVector = nextPoint - p;
                angles.pop_front();
                indices.pop_front();
            }
            calculateLambdas(meshInfo.getNumberOfVertices(), points);
        }
            break;
        case 2:
            m_lambdas.reserve(meshInfo.getNumberOfVertices());
            for(size_t j=1; j <= meshInfo.getNumberOfVertices(); j++)
            {
                m_lambdas.push_back(0); // Lambda(m_vertexIndex, j, 0)
            }
            while(!indices.empty())
            {
                m_lambdas[indices.front()-1] += (1./d);
                indices.pop_front();
            }
            break;
        case 3:
        {
            std::list<T> neighbourDistances = localNeighbourhood.getNeighbourDistances();
            T sumOfDistances = 0;
            for(typename std::list<T>::iterator it = neighbourDistances.begin(); it != neighbourDistances.end(); it++)
            {
                sumOfDistances += *it;
            }
            T sumOfDistancesInv = 1./sumOfDistances;
            m_lambdas.reserve(meshInfo.getNumberOfVertices());
            for(size_t j=1; j <= meshInfo.getNumberOfVertices(); j++)
            {
                m_lambdas.push_back(0); //Lambda(m_vertexIndex, j, 0)
            }
            for(typename std::list<T>::iterator it = neighbourDistances.begin(); it != neighbourDistances.end(); it++)
            {
                m_lambdas[indices.front()-1] += ((*it)*sumOfDistancesInv);
                indices.pop_front();
            }
        }
            break;
        default:
            GISMO_ERROR("parametrizationMethod not valid: " << parametrizationMethod);
    }
}

template<class T>
const std::vector<T>& gsParametrization<T>::LocalParametrization::getLambdas() const
{
    return m_lambdas;
}

//*****************************************************************************************************
//*****************************************************************************************************
//*******************THE******INTERN******FUNCTIONS******ARE******NOW******FOLLOWING*******************
//*****************************************************************************************************
//*****************************************************************************************************

template<class T>
void gsParametrization<T>::LocalParametrization::calculateLambdas(const size_t N, VectorType& points)
{
    m_lambdas.reserve(N);
    for(size_t j=1; j <= N; j++)
    {
        m_lambdas.push_back(0); //Lambda(m_vertexIndex, j, 0)
    }
    Point2D p(0, 0, 0);
    size_t d = points.size();
    std::vector<T> my(d, 0);
    size_t l=1;
    size_t steps = 0;
    //size_t checkOption = 0;
    for(typename VectorType::const_iterator it=points.begin(); it != points.end(); it++)
    {
        gsLineSegment<2,T> actualLine(p, *it);
        for(size_t i=1; i < d-1; i++)
        {
            if(l+i == d)
                steps = d - 1;
            else
                steps = (l+i)%d - 1;
            //checkoption is set to another number, in case mu is negativ
            if(actualLine.intersectSegment(*(points.begin()+steps), *(points.begin()+(steps+1)%d)/*, checkOption*/))
            {
                //BarycentricCoordinates b(p, *it, *(points.begin()+steps), *(points.begin()+(steps+1)%d));
                // calculating Barycentric Coordinates
                gsMatrix<T, 3, 3> matrix;
                matrix.topRows(2).col(0) = *it;
                matrix.topRows(2).col(1) = *(points.begin()+steps);
                matrix.topRows(2).col(2) = *(points.begin()+(steps+1)%d);
                matrix.row(2).setOnes();

                gsVector3d<T> vector3d;
                vector3d << p, 1;
                gsVector3d<T> delta = matrix.partialPivLu().solve(vector3d);
                my[l-1] = delta(0);
                my[steps] = delta(1);
                my[(steps + 1)%d] = delta(2);
                break;
            }
        }
        for(size_t k = 1; k <= d; k++)
        {
            m_lambdas[points[k-1].getVertexIndex()-1] += (my[k-1]);
        }
        std::fill(my.begin(), my.end(), 0);
        l++;
    }
    for(typename std::vector<T>::iterator it=m_lambdas.begin(); it != m_lambdas.end(); it++)
    {
        *it /= d;
    }
    for(typename std::vector<T>::iterator it=m_lambdas.begin(); it != m_lambdas.end(); it++)
    {
        if(*it < 0)
            gsInfo << *it << "\n";
    }
}

//*******************************************************************************************
//***************************** nested class LocalNeighbourhood *****************************
//*******************************************************************************************

template<class T>
gsParametrization<T>::LocalNeighbourhood::LocalNeighbourhood(const gsHalfEdgeMesh<T>& meshInfo,
							     const size_t vertexIndex,
							     const bool innerVertex)
{
    GISMO_ASSERT(!((innerVertex && vertexIndex > meshInfo.getNumberOfInnerVertices()) || vertexIndex < 1),
                 "Vertex with index " << vertexIndex << " does either not exist (< 1) or is not an inner vertex (> "
                 << meshInfo.getNumberOfInnerVertices() << ").");

    m_vertexIndex = vertexIndex;
    std::queue<typename gsHalfEdgeMesh<T>::Halfedge>
        allHalfedges = meshInfo.getOppositeHalfedges(m_vertexIndex, innerVertex);
    std::queue<typename gsHalfEdgeMesh<T>::Halfedge> nonFittingHalfedges;
    m_neighbours.appendNextHalfedge(allHalfedges.front());
    m_angles.push_back((*meshInfo.getVertex(allHalfedges.front().getOrigin()) - *meshInfo.getVertex(m_vertexIndex))
                       .angle((*meshInfo.getVertex(allHalfedges.front().getEnd())
                               - *meshInfo.getVertex(vertexIndex))));
    m_neighbourDistances.push_back(allHalfedges.front().getLength());
    allHalfedges.pop();
    while (!allHalfedges.empty())
    {
        if (m_neighbours.isAppendableAsNext(allHalfedges.front()))
        {
            m_neighbours.appendNextHalfedge(allHalfedges.front());
            m_angles
                .push_back((*meshInfo.getVertex(allHalfedges.front().getOrigin()) - *meshInfo.getVertex(m_vertexIndex))
                           .angle((*meshInfo.getVertex(allHalfedges.front().getEnd())
                                   - *meshInfo.getVertex(m_vertexIndex))));
            m_neighbourDistances.push_back(allHalfedges.front().getLength());
            allHalfedges.pop();
            while (!nonFittingHalfedges.empty())
            {
                allHalfedges.push(nonFittingHalfedges.front());
                nonFittingHalfedges.pop();
            }
        }
        else if (m_neighbours.isAppendableAsPrev(allHalfedges.front()))
        {
            m_neighbours.appendPrevHalfedge(allHalfedges.front());
            m_angles
                .push_front((*meshInfo.getVertex(allHalfedges.front().getOrigin()) - *meshInfo.getVertex(m_vertexIndex))
                            .angle((*meshInfo.getVertex(allHalfedges.front().getEnd())
                                    - *meshInfo.getVertex(m_vertexIndex))));
            m_neighbourDistances.push_back(allHalfedges.front().getLength());
            allHalfedges.pop();
            while (!nonFittingHalfedges.empty())
            {
                allHalfedges.push(nonFittingHalfedges.front());
                nonFittingHalfedges.pop();
            }
        }
        else
        {
            nonFittingHalfedges.push(allHalfedges.front());
            allHalfedges.pop();
        }
    }
}

template<class T>
size_t gsParametrization<T>::LocalNeighbourhood::getVertexIndex() const
{
    return m_vertexIndex;
}

template<class T>
size_t gsParametrization<T>::LocalNeighbourhood::getNumberOfNeighbours() const
{
    return m_neighbours.getNumberOfVertices();
}

template<class T>
const std::list<size_t> gsParametrization<T>::LocalNeighbourhood::getVertexIndicesOfNeighbours() const
{
    return m_neighbours.getVertexIndices();
}

template<class T>
const std::list<T>& gsParametrization<T>::LocalNeighbourhood::getAngles() const
{
    return m_angles;
}

template<class T>
T gsParametrization<T>::LocalNeighbourhood::getInnerAngle() const
{
    T angle = 0;
    for(typename std::list<T>::const_iterator it=m_angles.begin(); it!=m_angles.end(); it++)
    {
        angle += (*it);
    }
    return angle;
}

template<class T>
std::list<T> gsParametrization<T>::LocalNeighbourhood::getNeighbourDistances() const
{
    return m_neighbourDistances;
}

template<class T>
void gsParametrization<T>::constructTwins(std::vector<std::pair<size_t, size_t> >& twins,
					  const gsMesh<T>& overlapMesh,
					  typename gsMesh<T>::gsVertexHandle uMinv0,
					  typename gsMesh<T>::gsVertexHandle uMaxv0,
					  typename gsMesh<T>::gsVertexHandle uMinv1,
					  typename gsMesh<T>::gsVertexHandle uMaxv1)
{
    // TODO: Remove the debugging messages instead of just commenting them out.
    size_t currentNrAllVertices = m_mesh.getNumberOfVertices();

    gsHalfEdgeMesh<T> overlapHEM(overlapMesh);
    std::list<size_t> vertexIndices = overlapHEM.getBoundaryVertexIndices();

    // gsInfo << "Boundary indices\n";
    // for(auto it=vertexIndices.begin(); it!=vertexIndices.end(); ++it)
    // 	gsInfo << *it << ": " << *overlapHEM.getVertexUnsorted(*it) << "; ";
    // gsInfo << std::endl;

    // Rotate vertexIndices so that they form the right boundary.
    size_t from = overlapHEM.findVertex(uMinv1);
    // gsInfo << "From: " << from << ", uMinv1: " << *uMinv1 << std::endl;
    size_t i=0;
    while(vertexIndices.front() != from)
    {
	i++;
	vertexIndices.push_back(vertexIndices.front());
	vertexIndices.pop_front();
    }

    // Push the corresponding pairs to the twin vector.
    size_t to = overlapHEM.findVertex(uMinv0);
    // gsInfo << "To: " << to << ", uMinv0: " << *uMinv0 << std::endl;
    for(std::list<size_t>::const_iterator it=vertexIndices.begin(); *std::prev(it) != to; ++it)
    {
	size_t twin = m_mesh.getVertexIndex(overlapHEM.getVertexUnsorted(*it));
	// gsInfo << "Twin input: " << twin << "\n" << *overlapHEM.getVertexUnsorted(*it)
	//         << *m_mesh.getVertex(twin);
	twins.push_back(std::pair<size_t, size_t>(twin, ++currentNrAllVertices));
    }

    // Rotate vertexIndices so that they form the left boundary.
    from = overlapHEM.findVertex(uMaxv0);
    i=0;
    while(vertexIndices.front() != from)
    {
	i++;
	vertexIndices.push_back(vertexIndices.front());
	vertexIndices.pop_front();
    }

    // Push the corresponding pairs to the twin vector.
    to = overlapHEM.findVertex(uMaxv1);
    for(std::list<size_t>::const_iterator it=vertexIndices.begin(); *std::prev(it) != to; ++it)
    {
	size_t twin = m_mesh.getVertexIndex(overlapHEM.getVertexUnsorted(*it));
	// gsInfo << "Twin input: " << twin << "\n" << *overlapHEM.getVertexUnsorted(*it)
	//        << *m_mesh.getVertex(twin);
	twins.push_back(std::pair<size_t, size_t>(++currentNrAllVertices, twin));
    }

    // gsInfo << "Twins:\n";
    // for(auto it=twins.begin(); it!=twins.end(); ++it)
    // 	gsInfo << "(" << it->first << ", " << it->second << ")\n";
}

template <class T>
gsParametrization<T>& gsParametrization<T>::compute_periodic_overlap(std::string bottomFile,
								     std::string topFile,
								     std::string overlapFile,
								     std::vector<size_t>& left,
								     std::vector<size_t>& right)
{
    // TODO: Reading of pars and pts can be extracted into an inline function.
    gsFileData<> fd_v0(bottomFile);
    gsMatrix<> pars, pts;
    fd_v0.getId<gsMatrix<> >(0, pars);
    fd_v0.getId<gsMatrix<> >(1, pts);

    GISMO_ASSERT(pars.cols() == pts.cols(), "The numbers of parameters and points of v0 differ.");
    
    std::vector<size_t> indicesV0;
    std::vector<T> valuesV0;

    for(index_t c=0; c<pts.cols(); c++)
    {
	indicesV0.push_back(m_mesh.findVertex(pts(0, c), pts(1, c), pts(2, c), true));
	valuesV0.push_back(pars(0, c));
    }

    gsFileData<> fd_v1(topFile);
    fd_v1.getId<gsMatrix<> >(0, pars);
    fd_v1.getId<gsMatrix<> >(1, pts);

    GISMO_ASSERT(pars.cols() == pts.cols(), "The numbers of parameters and points of v1 differ.");

    std::vector<size_t> indicesV1;
    std::vector<T> valuesV1;

    for(index_t c=0; c<pts.cols(); c++)
    {
	indicesV1.push_back(m_mesh.findVertex(pts(0, c), pts(1, c), pts(2, c), true));
	valuesV1.push_back(pars(0, c));
    }

    gsFileData<> fd_overlap(overlapFile);
    gsMesh<real_t>::uPtr overlapMesh = fd_overlap.getFirst<gsMesh<real_t> >();

    calculate_periodic_overlap(m_options.getInt("parametrizationMethod"),
			       indicesV0, valuesV0, indicesV1, valuesV1, *overlapMesh, left, right);
    return *this;
}

template <class T>
gsParametrization<T>& gsParametrization<T>::compute_periodic_stitch(std::string bottomFile,
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
	indicesV0.push_back(m_mesh.findVertex(pts(0, c), pts(1, c), pts(2, c), true));
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
	indicesV1.push_back(m_mesh.findVertex(pts(0, c), pts(1, c), pts(2, c), true));
	valuesV1.push_back(pars(0, c));
    }

    // Read the stitch indices.
    gsFileData<> fd_overlap(stitchFile);
    gsMatrix<> stitchPoints;
    fd_overlap.getId<gsMatrix<> >(0, stitchPoints);
    std::vector<size_t> stitchIndices;

    for(index_t c=0; c<stitchPoints.cols(); c++)
	stitchIndices.push_back(m_mesh.findVertex(stitchPoints(0, c), stitchPoints(1, c), stitchPoints(2, c), true));

    // Calculation itself.
    calculate_periodic_stitch(m_options.getInt("parametrizationMethod"),
			      indicesV0, valuesV0, indicesV1, valuesV1, stitchIndices,
			      posCorrections);
    return *this;
}

template<class T>
void gsParametrization<T>::calculate_periodic_overlap(const size_t paraMethod,
						      const std::vector<size_t>& indicesV0,
						      const std::vector<T>& valuesV0,
						      const std::vector<size_t>& indicesV1,
						      const std::vector<T>& valuesV1,
						      const gsMesh<T>& overlapMesh,
						      std::vector<size_t>& left,
						      std::vector<size_t>& right)
{
    size_t n = m_mesh.getNumberOfInnerVertices();
    size_t N = m_mesh.getNumberOfVertices();

    Neighbourhood neighbourhood(m_mesh, paraMethod);

    m_parameterPoints.reserve(N);
    for (size_t i = 1; i <= n; i++)
    {
	m_parameterPoints.push_back(Point2D(0, 0, i));
    }

    // Add the parameters of the boundary points.
    GISMO_ASSERT(indicesV0.size() == valuesV0.size(), "Different sizes of u0.");
    GISMO_ASSERT(indicesV1.size() == valuesV1.size(), "Different sizes of u1.");
    GISMO_ASSERT(indicesV0.size() + indicesV1.size() == m_mesh.getNumberOfBoundaryVertices(),
		 "Not prescribing all boundary points.");

    size_t numPtsSoFar = n;
    m_parameterPoints.resize(n + indicesV0.size() + indicesV1.size());

    for(size_t i=0; i<indicesV0.size(); i++)
    	m_parameterPoints[indicesV0[i]-1] = Point2D(valuesV0[i], 0, numPtsSoFar++);

    for(size_t i=0; i<indicesV1.size(); i++)
	m_parameterPoints[indicesV1[i]-1] = Point2D(valuesV1[i], 1, numPtsSoFar++);

    // Construct the twins.
    std::vector<std::pair<size_t, size_t> > twins;
    constructTwins(twins, overlapMesh,
		   m_mesh.getVertex(indicesV0.front()),
		   m_mesh.getVertex(indicesV0.back()),
		   m_mesh.getVertex(indicesV1.front()),
		   m_mesh.getVertex(indicesV1.back()));

    // Solve.
    constructAndSolveEquationSystem(neighbourhood, n, N, twins);

    // Remember the boundaries.
    for(auto it=twins.begin(); it!=twins.end(); ++it)
    {
	if(it->first < it->second)
	    left.push_back(it->first);
	else
	    right.push_back(it->second);
    }
}

template<class T>
void gsParametrization<T>::calculate_periodic_stitch(const size_t paraMethod,
						     const std::vector<size_t>& indicesV0,
						     const std::vector<T>& valuesV0,
						     const std::vector<size_t>& indicesV1,
						     const std::vector<T>& valuesV1,
						     const std::vector<size_t>& stitchIndices,
						     std::vector<std::vector<size_t> >& posCorrections)
{
    size_t n = m_mesh.getNumberOfInnerVertices();
    size_t N = m_mesh.getNumberOfVertices();

    //std::vector<std::vector<size_t> > posCorrections(N);
    posCorrections = std::vector<std::vector<size_t> >(N);
    std::vector<std::vector<size_t> > negCorrections(N);
    Neighbourhood neighbourhood(m_mesh, stitchIndices, posCorrections, negCorrections, paraMethod);

    m_parameterPoints.reserve(N);
    for (size_t i = 1; i <= n; i++)
    {
	m_parameterPoints.push_back(Point2D(0, 0, i));
    }

    // Add the parameters of the boundary points.
    GISMO_ASSERT(indicesV0.size() == valuesV0.size(), "Different sizes of u0.");
    GISMO_ASSERT(indicesV1.size() == valuesV1.size(), "Different sizes of u1.");
    GISMO_ASSERT(indicesV0.size() + indicesV1.size() == m_mesh.getNumberOfBoundaryVertices(),
		 "Not prescribing all boundary points.");

    size_t numPtsSoFar = n;
    m_parameterPoints.resize(n + indicesV0.size() + indicesV1.size());

    for(size_t i=0; i<indicesV0.size(); i++)
    	m_parameterPoints[indicesV0[i]-1] = Point2D(valuesV0[i], 0, numPtsSoFar++);

    for(size_t i=0; i<indicesV1.size(); i++)
	m_parameterPoints[indicesV1[i]-1] = Point2D(valuesV1[i], 1, numPtsSoFar++);

    /// Solve.
    constructAndSolveEquationSystem(neighbourhood, n, N, posCorrections, negCorrections);
}

template<class T>
std::vector<size_t> gsParametrization<T>::readIndices(const std::string& filename) const
{
    gsMatrix<> pts;
    gsFileData<> fd(filename);
    fd.getId<gsMatrix<> >(0, pts);

    std::vector<size_t> result;
    for(index_t c=0; c<pts.cols(); c++)
	result.push_back(m_mesh.findVertex(pts(0, c), pts(1, c), pts(2, c), true));

    return result;
}

template <class T>
std::vector<size_t> gsParametrization<T>::getSide(const std::list<size_t>& boundary, size_t beg, size_t end) const
{
    gsInfo << "Boundary (" << boundary.size() << " elements):\n";
    for(auto it = boundary.begin(); it != boundary.end(); ++it)
	gsInfo << *it << " ";
    gsInfo << std::endl;

    auto itBeg = std::find(boundary.begin(), boundary.end(), beg);
    auto itEnd = std::find(boundary.begin(), boundary.end(), end);

    GISMO_ASSERT(itBeg != boundary.end(), "Beg corner not found.");
    GISMO_ASSERT(itEnd != boundary.end(), "End corner not found.");

    std::vector<size_t> result;
    for(auto it = std::next(itBeg); it != itEnd; ++it)
    {
	if(it == boundary.end())
	    it = boundary.begin();
	result.push_back(*it);
    }
    return result;
}

template<class T>
void gsParametrization<T>::calculate_free_boundary(const size_t paraMethod,
						   const std::string& fileCorners)
{
    size_t n = m_mesh.getNumberOfInnerVertices();
    size_t N = m_mesh.getNumberOfVertices();

    Neighbourhood neighbourhood(m_mesh, paraMethod);

    std::vector<size_t> corners = readIndices(fileCorners);
    GISMO_ASSERT(corners.size() == 4, "Wrong number of corners.");
    for(auto it=corners.begin(); it!=corners.end(); ++it)
	gsInfo << "Corner " << *it << std::endl;

    // TODO: Actually, this conversion from unsorted to sorted
    // numbering can be done much easier, as n+1, ..., N.
    std::list<size_t> unsortedBoundary = m_mesh.getBoundaryVertexIndices();
    std::list<size_t> sortedBoundary;
    for(auto it=unsortedBoundary.begin(); it!=unsortedBoundary.end(); ++it)
	sortedBoundary.push_back(m_mesh.getVertexIndex(m_mesh.getVertexUnsorted(*it)));

    std::vector<size_t> v0 = getSide(sortedBoundary, corners[0], corners[1]);
    std::vector<size_t> u1 = getSide(sortedBoundary, corners[1], corners[2]);
    std::vector<size_t> v1 = getSide(sortedBoundary, corners[2], corners[3]);
    std::vector<size_t> u0 = getSide(sortedBoundary, corners[3], corners[0]);

    /// Solve.
    constructAndSolveEquationSystem(neighbourhood, n, N, corners, v0, u1, v1, u0);
}

template <class T>
void gsParametrization<T>::updateLambdasWithTwins(std::vector<T>& lambdas,
						  const std::vector<std::pair<size_t, size_t> >& twins,
						  size_t vertexId) const
{
    lambdas.reserve(lambdas.size() + twins.size());
    for(size_t i=0; i<twins.size(); i++)
	lambdas.push_back(0);

    // Determine, whether vertexId is on the left or right side of the overlap.
    bool isLeft  = false;
    bool isRight = false;

    for(auto it=twins.begin(); it!=twins.end(); ++it)
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


    for(size_t i=0; i<twins.size(); i++)
    {
	size_t first=twins[i].first-1;
	size_t second=twins[i].second-1;

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

// Every twin pair is a pair of vertices with (almost) the same coordinates.
// The first in the pair has the u-coordinate smaller by one than the second.
template <class T>
void gsParametrization<T>::constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
							   const size_t n,
							   const size_t N,
							   const std::vector<std::pair<size_t, size_t> >& twins)
{
    size_t numTwins = twins.size();
    gsMatrix<T> LHS(N + numTwins, N + numTwins);
    gsMatrix<T> RHS(N + numTwins, 2);
    std::vector<T> lambdas;


    // interior points
    for (size_t i = 0; i < n; i++)
    {
        lambdas = neighbourhood.getLambdas(i);
	updateLambdasWithTwins(lambdas, twins, i+1);

        for (size_t j = 0; j < N + numTwins; j++)
        {
	    // Standard way:
            LHS(i, j) = ( i==j ? T(1) : -lambdas[j] );
	    // LHS(i, j) = lambdas[j];
	    // RHS(i, 0) = 0.5;
	    // RHS(i, 1) = 0.5;
        }
    }

    // points on the lower and upper boundary
    for (size_t i=n; i<N; i++)
    {
	LHS(i, i) = T(1);
	RHS.row(i) = m_parameterPoints[i];
    }

    // points on the overlap
    for (size_t i=N; i<N+numTwins; i++)
    {
	size_t first   = twins[i-N].first-1;
	size_t second  = twins[i-N].second-1;

	LHS(i, first)  = T(1);
	LHS(i, second) = T(-1);

	RHS(i, 0)      = T(-1);
	RHS(i, 1)      = T(0);
    }

    gsMatrix<T> sol;
    Eigen::PartialPivLU<typename gsMatrix<T>::Base> LU = LHS.partialPivLu();
    sol = LU.solve(RHS);
    for (size_t i = 0; i < N; i++)
    {
    	m_parameterPoints[i] << sol(i, 0), sol(i, 1);
    }

    // for(size_t i=0; i<15; i++)
    // {
    // 	sol = LHS * RHS;
    // 	RHS = sol;

    // gsMatrix<T> uv(2, N), aux(2, twins.size());

    // for (size_t i = 0; i < N; i++)
    // {
    // 	m_parameterPoints[i] << sol(i, 0), sol(i, 1);
    // 	uv(0, i) = sol(i, 0);
    // 	uv(1, i) = sol(i, 1);
    // }

    // for(size_t i=0; i<numTwins; i++)
    // {
    // 	aux(0, i) = sol(N+i, 0);
    // 	aux(1, i) = sol(N+i, 1);
    // }

    // // Practical, keep:
    // //gsWriteParaviewPoints(uv, "uv");
    // //gsWriteParaviewPoints(aux, "aux");

    // std::vector<gsMesh<> > meshes = createThreeFlatMeshes(twins);
    // gsWriteParaview(meshes[0], "left_"  + std::to_string(i));
    // gsWriteParaview(meshes[0], "mid_"   + std::to_string(i));
    // gsWriteParaview(meshes[0], "right_" + std::to_string(i));
    // }
}

template <class T>
void gsParametrization<T>::constructAndSolveEquationSystem(
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
	RHS.row(i) = m_parameterPoints[i];
    }

    Eigen::PartialPivLU<typename gsMatrix<T>::Base> LU = LHS.partialPivLu();
    gsMatrix<T> sol = LU.solve(RHS);
    for (size_t i = 0; i < n; i++)
    {
    	m_parameterPoints[i] << sol(i, 0), sol(i, 1);
    }
}

template <class T>
void gsParametrization<T>::constructAndSolveEquationSystem(
    const Neighbourhood &neighbourhood,
    const size_t n,
    const size_t N,
    const std::vector<size_t>& corners,
    const std::vector<size_t>& botBoundary,
    const std::vector<size_t>& rgtBoundary,
    const std::vector<size_t>& topBoundary,
    const std::vector<size_t>& lftBoundary)
{
    GISMO_ASSERT(corners.size() + botBoundary.size() + rgtBoundary.size()
		 + topBoundary.size() + lftBoundary.size() == N-n,
		 "Wrong number of boundary points.");

    std::vector<T> lambdas;
    gsMatrix<T> LHSx(N, N), LHSy(N, N);
    gsMatrix<T> RHSx(N, 1), RHSy(N, 1);

    // interior points
    for (size_t i = 0; i < n; i++)
    {
        lambdas = neighbourhood.getLambdas(i);
        for (size_t j = 0; j < N; j++)
        {
            LHSx(i, j) = ( i==j ? T(1) : -lambdas[j] );
	    LHSy(i, j) = LHSx(i, j);
	}
    }

    // corners
    for (size_t i=n; i<4; i++)
    {
	size_t j = corners[i];
	LHSx(j, j) = T(1);
	LHSy(i, j) = T(1);
	switch(i) {
	case 0:
	    RHSx(j, 0) = 0;
	    RHSy(j, 0) = 0;
	    break;
	case 1:
	    RHSx(j, 0) = 0;
	    RHSy(j, 0) = 1;
	    break;
	case 2:
	    RHSx(j, 0) = 1;
	    RHSy(j, 0) = 1;
	    break;
	case 3:
	    RHSx(j, 0) = 1;
	    RHSy(j, 0) = 0;
	    break;
	default:
	    gsWarn << "This switch case should never happen." << std::endl;
	}
    }

    // Bottom boundary:
    for(auto it = botBoundary.begin(); it != botBoundary.end(); ++it)
    {
	size_t i = *it;
	lambdas = neighbourhood.getLambdas(i);
        for (size_t j = 0; j < N; j++)
	{
	    LHSx(i, j) = ( i==j ? T(1) : -lambdas[j] );
	    LHSy(i, j) = ( -lambdas[j] );
	}
    }

    // Top boundary:
    for(auto it = topBoundary.begin(); it != topBoundary.end(); ++it)
    {
	size_t i = *it;
	// TODO next time: this crashes, as lambdas don't seem to be
	// pre-computed for the boundary vertices.
	lambdas = neighbourhood.getLambdas(i);
        for (size_t j = 0; j < N; j++)
	{
	    LHSx(i, j) = ( i==j ? T(1) : -lambdas[j] );
	    LHSy(i, j) = ( -lambdas[j] );
	    RHSy(i, 0) = T(1);
	}
    }


    // Left boundary:
    for(auto it = lftBoundary.begin(); it != lftBoundary.end(); ++it)
    {
	size_t i = *it;
	lambdas = neighbourhood.getLambdas(i);
        for (size_t j = 0; j < N; j++)
	{
	    LHSx(i, j) = ( -lambdas[j] );
	    LHSy(i, j) = ( i==j ? T(1) : -lambdas[j] );
	}
    }

    // Right boundary:
    for(auto it = rgtBoundary.begin(); it != rgtBoundary.end(); ++it)
    {
	size_t i = *it;
	gsInfo << "Assembling row " << i << ".\n";
	lambdas = neighbourhood.getLambdas(i);
	// for(auto jt = lambdas.begin(); jt != lambdas.end(); ++jt)
	//     gsInfo << *jt << " ";
	// gsInfo << std::endl << std::endl;
        for (size_t j = 0; j < N; j++)
	{
	    LHSx(i, j) = ( -lambdas[j] );
	    LHSy(i, j) = ( i==j ? T(1) : -lambdas[j] );
	    RHSx(i, 0) = T(1);
	}
    }
    
    Eigen::PartialPivLU<typename gsMatrix<T>::Base> LUx = LHSx.partialPivLu();
    gsInfo << "det(x): " << LUx.determinant() << std::endl;
    gsMatrix<T> solx = LUx.solve(RHSx);

    Eigen::PartialPivLU<typename gsMatrix<T>::Base> LUy = LHSy.partialPivLu();
    gsInfo << "det(y): " << LUy.determinant() << std::endl;
    gsMatrix<T> soly = LUy.solve(RHSy);

    for (size_t i = 0; i < n; i++)
    {
    	m_parameterPoints[i] << solx(i, 0), soly(i, 0);
    }
}

} // namespace gismo
