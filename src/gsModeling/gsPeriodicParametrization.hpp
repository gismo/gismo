/** @file gsPeriodicParametrization.hpp

    @brief Provides implementation gsPeriodicParametrization class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris

*/

#include <gsModeling/gsPeriodicParametrization.h>

namespace gismo
{

/* Nested class FlatMesh */

template<class T>
real_t gsPeriodicParametrization<T>::FlatMesh::correspondingV(const VertexHandle& h0,
							      const VertexHandle& h1,
							      real_t u) const
{
    real_t u0 = (*h0)[0];
    real_t u1 = (*h1)[0];
    real_t v0 = (*h0)[1];
    real_t v1 = (*h1)[1];

    real_t t = (u - u0) / (u1 - u0);

    return (1 - t) * v0 + t * v1;
}

template<class T>
void gsPeriodicParametrization<T>::FlatMesh::addThreeFlatTrianglesOneOut(gsMesh<T>& mesh,
									 const VertexHandle& v0,
									 const VertexHandle& v1,
									 const VertexHandle& v2) const
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

template<class T>
void gsPeriodicParametrization<T>::FlatMesh::addThreeFlatTrianglesTwoOut(gsMesh<T>& mesh,
                                                                         const VertexHandle& v0,
                                                                         const VertexHandle& v1,
                                                                         const VertexHandle& v2) const
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
void gsPeriodicParametrization<T>::FlatMesh::addOneFlatTriangleNotIntersectingBoundary(gsMesh<T>& mesh,
                                                                                       const typename gsMesh<T>::VertexHandle& v0,
                                                                                       const typename gsMesh<T>::VertexHandle& v1,
                                                                                       const typename gsMesh<T>::VertexHandle& v2) const
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
gsMesh<T> gsPeriodicParametrization<T>::FlatMesh::createRestrictedFlatMesh() const
{
    gsMesh<T> result;

    for(size_t i=0; i<m_unfolded.getNumberOfTriangles(); i++)
    {
        // Remember the corners and which of them are inside the domain.
        bool out[3];
        typename gsMesh<T>::VertexHandle vh[3];
        for(size_t j=1; j<=3; ++j)
        {
            vh[j-1] = m_unfolded.getVertex(m_unfolded.getGlobalVertexIndex(j, i));
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

// back to implementing the main class

template <class T>
void gsPeriodicParametrization<T>::restrictMatrices(gsMatrix<T>& uv, const gsMatrix<T>& xyz,
						    real_t uMin, real_t uMax) const
{
    real_t uLength = uMax - uMin;
    for(index_t j=0; j<uv.cols(); j++)
    {
        real_t u = uv(0, j);

        if(u < uMin)
            uv(0, j) += uLength;
        else if(u > uMax)
            uv(0 ,j) -= uLength;
    }
}

template <class T>
void gsPeriodicParametrization<T>::initParameterPoints()
{
    typedef typename gsParametrization<T>::Point2D Point2D;

    size_t n = this->m_mesh.getNumberOfInnerVertices();
    size_t N = this->m_mesh.getNumberOfVertices();

    this->m_parameterPoints.reserve(N);
    for (size_t i = 1; i <= n; i++)
        this->m_parameterPoints.push_back(Point2D(0, 0, i));

    // Save the sizes as size_t to compare without warnings.
    size_t v0size = m_paramsV0.cols();
    size_t v1size = m_paramsV1.cols();
    GISMO_ASSERT(m_indicesV0.size() == v0size, "Different sizes of u0.");
    GISMO_ASSERT(m_indicesV1.size() == v1size, "Different sizes of u1.");
    GISMO_ASSERT(v0size + v1size == this->m_mesh.getNumberOfBoundaryVertices(),
                 "Not prescribing all boundary points.");

    size_t numPtsSoFar = n;
    this->m_parameterPoints.resize(n + v0size + v1size);

    // Set the parameter values on the v=0 boundary.
    for(size_t i=0; i<v0size; i++)
        this->m_parameterPoints[m_indicesV0[i]-1] = Point2D(m_paramsV0(0, i), 0, numPtsSoFar++);

    // Set the parameter values on the v=1 boundary.
    for(size_t i=0; i<v1size; i++)
        this->m_parameterPoints[m_indicesV1[i]-1] = Point2D(m_paramsV1(0, i), 1, numPtsSoFar++);
}

} // namespace gismo
