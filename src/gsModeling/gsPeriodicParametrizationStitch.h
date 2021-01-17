/** @file gsPeriodicParametrizationStitch.h

    @brief Implementation of periodic Floater parametrization using a
    stitch. This is an alternative to gsPeriodicParametrizationOverlap
    class.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
*/

#pragma once

#include "gsModeling/gsPeriodicParametrizationOverlap.h"
#include "gsIO/gsOptionList.h"

namespace gismo
{

template <class T>
class GISMO_EXPORT gsPeriodicParametrizationStitch : public gsParametrization<T>
{
    class Neighbourhood : public gsParametrization<T>::Neighbourhood
    {
    public:
	typedef typename gsParametrization<T>::LocalNeighbourhood LocalNeighbourhood;

	explicit Neighbourhood(const gsHalfEdgeMesh<T> &meshInfo,
			       const std::vector<size_t>& stitchIndices,
			       std::vector<std::vector<size_t> >& posCorrections,
			       std::vector<std::vector<size_t> >& negCorrections,
			       const size_t parametrizationMethod = 2);

    private:
	std::vector<size_t> computeCorrections(const std::vector<size_t>& stitchIndices,
					       const LocalNeighbourhood& localNeighbourhood) const;
    };

    //typedef typename gsParametrization<T>::Neighbourhood Neighbourhood;

public:
    explicit gsPeriodicParametrizationStitch(gsMesh<T> &mesh,
					     const gsOptionList &list = gsParametrization<T>::defaultOptions())
	: gsParametrization<T>(mesh, list)
	{}

    /// Periodic parametrization using Pierre's trick.
    gsPeriodicParametrizationStitch<T>& compute_periodic_stitch(std::string bottomFile,
								std::string topFile,
								std::string stitchFile);

protected:
    void calculate_periodic_stitch(const size_t paraMethod,
				   const std::vector<size_t>& indicesV0,
				   const std::vector<T>& valuesV0,
				   const std::vector<size_t>& indicesV1,
				   const std::vector<T>& valuesV1,
				   const std::vector<size_t>& stitchIndices);

    /** Similar to @a constructAndSolveEquationSystem but works for periodic meshes using
     * the corrections.
     */ // TODO: Explain the parameters.
    void constructAndSolveEquationSystem(const Neighbourhood &neighbourhood,
					 const size_t n,
					 const size_t N);

    // From here on the visualisation methods
public:

    using gsParametrization<T>::createFlatMesh;
    // TODO: remove
    //using gsPeriodicParametrizationOverlap<T>::createFlatMesh;

    /**
     * Creates a flat mesh out of a periodic parametrization created by a the stitch method.
     * @param posCorrections Positive corrections from the stitch algorithm.
     * @param restrict If set to true, the mesh is restricted to [0, 1]^2.
     */
    gsMesh<T> createFlatMesh(bool restrict) const;

protected:
    std::vector<std::vector<size_t> > m_posCorrections;
    std::vector<std::vector<size_t> > m_negCorrections;
};

} // namespace gismo
