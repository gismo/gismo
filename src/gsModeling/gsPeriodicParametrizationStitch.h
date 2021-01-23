/** @file gsPeriodicParametrizationStitch.h

    @brief Implementation of periodic Floater parametrization using a
    stitch. The idea is adapted from

    Tong, Y., Alliez, P., Cohen-Steiner, D., Desbrun, M.: Designing
    quadrangulations with discrete harmonic forms, in: Sheffer, A.,
    Polthier, K. (Eds.), Symposium on Geometry Processing,
    Eurographics. pp. 201–210, 2006,

    where it was used constructing discrete harmonic mappings on
    arbitrary topology.

    This class is an alternative to gsPeriodicParametrizationOverlap.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
*/

#pragma once

#include "gsModeling/gsPeriodicParametrization.h"
#include "gsIO/gsOptionList.h"

namespace gismo
{

template <class T>
class GISMO_EXPORT gsPeriodicParametrizationStitch : public gsPeriodicParametrization<T>
{
    class Neighbourhood : public gsParametrization<T>::Neighbourhood
    {
    public:
	typedef typename gsParametrization<T>::LocalNeighbourhood LocalNeighbourhood;

	explicit Neighbourhood(const gsHalfEdgeMesh<T> &meshInfo,
			       const std::vector<size_t>& stitchIndices,
			       gsMatrix<int>& corrections,
			       const size_t parametrizationMethod = 2);

    private:
	std::vector<size_t> computeCorrections(const std::vector<size_t>& stitchIndices,
					       const LocalNeighbourhood& localNeighbourhood) const;
    };

    //typedef typename gsParametrization<T>::Neighbourhood Neighbourhood;

public:
    explicit gsPeriodicParametrizationStitch(gsMesh<T> &mesh,
					     const gsOptionList &list = gsParametrization<T>::defaultOptions())
	: gsPeriodicParametrization<T>(mesh, list)
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

    // TODO: remove
    using gsParametrization<T>::createFlatMesh;

    /**
     * Creates a flat mesh out of a periodic parametrization created by a the stitch method.
     * @param posCorrections Positive corrections from the stitch algorithm.
     * @param restrict If set to true, the mesh is restricted to [0, 1]^2.
     */
    gsMesh<T> createFlatMesh(bool restrict) const
    {
    	gsMesh<T> unfolded = createUnfoldedFlatMesh();
    	if(restrict)
    	{
    		typename gsPeriodicParametrization<T>::FlatMesh display(unfolded);
    		return display.createRestrictedFlatMesh();
    	}
    	else
    	    return unfolded;
    }

protected:
    gsMesh<T> createUnfoldedFlatMesh() const;

    // Alternatively, one can save the stitch vertices into a vector beforehand.
    // vertexIndex is numbered from 1
    bool isOnStitch(size_t vertexIndex) const
    {
	for(index_t c=0; c<m_corrections.cols(); c++)
	    if(m_corrections(vertexIndex-1, c) == 1)
		return true;
	return false;
    }

    bool edgeIsInCorrections(index_t beg, index_t end) const
    {
	return ((m_corrections(beg, end) ==  1) ||
		(m_corrections(beg, end) == -1) ||
		(m_corrections(end, beg) ==  1) ||
		(m_corrections(end, beg) == -1));
	// Actually, the first two conditions should be enough.		
    }

protected:
    gsMatrix<int> m_corrections;
};

} // namespace gismo
