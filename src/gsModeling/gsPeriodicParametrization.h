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
    class FlatMesh
    {
    public:
	FlatMesh(const gsMesh<T>& unfolded)
	    : m_unfolded(unfolded)
	{}

	gsMesh<T> createRestrictedFlatMesh() const;

    private:
	real_t correspondingV(const typename gsMesh<T>::VertexHandle& v0,
			      const typename gsMesh<T>::VertexHandle& v1,
			      real_t u) const;

	void addThreeFlatTrianglesOneOut(gsMesh<T>& mesh,
					 const typename gsMesh<T>::VertexHandle& v0,
					 const typename gsMesh<T>::VertexHandle& v1,
					 const typename gsMesh<T>::VertexHandle& v2) const;

	void addThreeFlatTrianglesTwoOut(gsMesh<T>& mesh,
					 const typename gsMesh<T>::VertexHandle& v0,
					 const typename gsMesh<T>::VertexHandle& v1,
					 const typename gsMesh<T>::VertexHandle& v2) const;

	void addOneFlatTriangleNotIntersectingBoundary(gsMesh<T>& mesh,
						       const typename gsMesh<T>::VertexHandle& v0,
						       const typename gsMesh<T>::VertexHandle& v1,
						       const typename gsMesh<T>::VertexHandle& v2) const;
    private: // members
	gsHalfEdgeMesh<T> m_unfolded;
    };

public:
    gsPeriodicParametrization(gsMesh<T>& mesh,
			      const gsOptionList &list = gsParametrization<T>::defaultOptions())
	: gsParametrization<T>(mesh, list)
	{}

    using gsParametrization<T>::defaultOptions;
};

} // namespace gismo
