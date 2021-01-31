/** @file gsPeriodicParametrization.h

    @brief Abstract class with the functionality common to
    gsPeriodicParametrizationStitch and
    gsPeriodicParametrizationOverlap.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
*/

#pragma once

#include "gsModeling/gsParametrization.h"

namespace gismo
{

template <class T>
class gsPeriodicParametrization : public gsParametrization<T>
{

public:

    typedef memory::shared_ptr<gsPeriodicParametrization<T> > uPtr;

    /// Nested class for plotting flat meshes restricted to [0, 1]^2.
    class FlatMesh
    {
	typedef typename gsMesh<T>::VertexHandle VertexHandle;
    public:

	/// Constructor.
	/// @a unfolded Flat mesh possibly intersecting the domain boundaries.
	FlatMesh(const gsMesh<T>& unfolded)
	    : m_unfolded(unfolded)
	{}

	/// Trims the mesh to [0, 1]^2.
	gsMesh<T> createRestrictedFlatMesh() const;

    protected:

	/// Finds the v-coordinate of the point on the line segment
	/// (@a v0, @a v1) that has the u-coordinate @a u.
	real_t correspondingV(const VertexHandle& v0,
			      const VertexHandle& v1,
			      real_t u) const;

	/// Adds three flat triangles in the situation where one of
	/// the vertices of the original triangle is outside the
	/// domain.
	/// @a v1 is outside the domain, @a v0 and @a v2 inside.
	void addThreeFlatTrianglesOneOut(gsMesh<T>& mesh,
					 const VertexHandle& v0,
					 const VertexHandle& v1,
					 const VertexHandle& v2) const;

	/// Adds three flat triangles in the situation where two of
	/// the vertices of the original triangle are outside the
	/// domain.
	/// @a v1 is inside the domain, @a v0 and @a v2 outside.
	void addThreeFlatTrianglesTwoOut(gsMesh<T>& mesh,
					 const VertexHandle& v0,
					 const VertexHandle& v1,
					 const VertexHandle& v2) const;

	/// Adds a flat triangle and shifts it inside the domain if necessary.
	void addOneFlatTriangleNotIntersectingBoundary(gsMesh<T>& mesh,
						       const VertexHandle& v0,
						       const VertexHandle& v1,
						       const VertexHandle& v2) const;

    protected: // members
	gsHalfEdgeMesh<T> m_unfolded;///< flat mesh possibly intersecting the domain boundaries
    };

public:

    /** Constructor
     * @param mesh the surface mesh to be parametrized
     * @param verticesV0 vertices on the bottom (i.e., v = 0) boundary
     * @param paramsV0 their prescribed parameters
     * @param verticessV1 vertices on the upper (i.e., v = 1) boundary
     * @param paramsV1 their prescribed parameters
     * @param list list of the method options
     */
    gsPeriodicParametrization(const gsMesh<T>& mesh,
			      const gsMatrix<T>& verticesV0,
			      const gsMatrix<T>& paramsV0,
			      const gsMatrix<T>& verticesV1,
			      const gsMatrix<T>& paramsV1,
			      const gsOptionList &list = gsParametrization<T>::defaultOptions())
	: gsParametrization<T>(mesh, list),
	  m_paramsV0(paramsV0), m_paramsV1(paramsV1),
	  m_indicesV0(this->indices(verticesV0)),
	  m_indicesV1(this->indices(verticesV1))
	{
	    GISMO_ASSERT(this->m_paramsV0.rows() == 1, "one row expected in paramsV0");
	    GISMO_ASSERT(this->m_paramsV1.rows() == 1, "one row expected in paramsV1");

	    GISMO_ASSERT(this->m_valuesV0.rows() == 3, "three rows expected in valuesV0");
	    GISMO_ASSERT(this->m_valuesV1.rows() == 3, "three rows expected in valuesV1");

	    GISMO_ASSERT(this->m_paramsV0.cols() == this->m_valuesV0.cols(),
			 "paramsV0 and valuesV0 are required to have the same number of cols");
	    GISMO_ASSERT(this->m_paramsV1.cols() == this->m_valuesV1.cols(),
			 "paramsV1 and valuesV1 are required to have the same number of cols");
	}

    /**
     * Moves the u-coordinates of parameters outside the
     * interval [@a uMin, @a uMax] to inside the interval.
     * Note: it modifies uv!
     * @param uv Matrix of the parameters, one column per point
     * @param xyz Matrix of the coordinates, one column per point
     * @param uMin minimal desired u
     * @param uMax maximal desired u
     */
    void restrictMatrices(gsMatrix<T>& uv, const gsMatrix<T>& xyz,
			  real_t uMin = 0, real_t uMax = 1) const;

    /// Computes the periodic parametrization.
    virtual void compute() = 0;

protected:

    /**
     * Prepares m_parameterPoints for a computation of a periodic parametrization.
     * Important: Neighbourhood has to be constructed before calling this function.
     */
    void initParameterPoints();

protected: // members

    const gsMatrix<T>& m_paramsV0;         ///< u-parameters of the vertices with v=0
    const gsMatrix<T>& m_paramsV1;         ///< u-parameters of the vertices with v=1
    const std::vector<size_t> m_indicesV0; ///< indices of the vertices with v=0
    const std::vector<size_t> m_indicesV1; ///< indicesV1 indices of the vertices with v=1

};

} // namespace gismo
