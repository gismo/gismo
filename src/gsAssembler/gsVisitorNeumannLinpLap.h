/** @file gsVisitorNeumann.h

@brief Neumann conditions visitor for elliptic problems.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsFuncData.h>

namespace gismo
{
	/** @brief
	Implementation of a Neumann BC for elliptic Assembler.

	It sets up an assembler and adds the following term to the linear term.
	\f[ \alpha(\nabla u_{i-1})\nabla u \cdot \mathbf{n} = g_N  \f]
	*/

	template <class T>
	class gsVisitorNeumannLinpLap
	{
	public:

		gsVisitorNeumannLinpLap(const gsPde<T> &pde, const boundary_condition<T> & s, index_t subdiv_ = 1)
			: neudata_ptr(s.function().get()), side(s.side()), subdiv(subdiv_)
		{ 
			const gsLinpLapPde<T> * pde_ptr = static_cast<const gsLinpLapPde<T>*>(&pde);

			eps = pde_ptr->eps;
			p = pde_ptr->p;
			w = pde_ptr->w;
		}

		/** @brief
		* Constructor of the assembler object
		*
		\param[in] neudata is the Neumann boundary data.
		\param[in] s are the sides of the geometry where neumann BC are prescribed.

		\f[ \nabla u \cdot \mathbf{n} = g_N  \f]
		*/
		gsVisitorNeumannLinpLap(const gsFunction<T> & neudata, boxSide s) :
			neudata_ptr(&neudata), side(s)
		{ }

		void initialize(const gsBasis<T> & basis,
			gsQuadRule<T>    & rule)
		{
			const int dir = side.direction();
			gsVector<int> numQuadNodes(basis.dim());
			for (int i = 0; i < basis.dim(); ++i)
				numQuadNodes[i] = basis.degree(i) + 1;
			numQuadNodes[dir] = 1;

			// Setup Quadrature
			rule = gsSubdividedRule<T, gsQuadRule<T> >(gsGaussRule<T>(numQuadNodes), subdiv);// harmless slicing occurs here

												// Set Geometry evaluation flags
			md.flags = NEED_VALUE | NEED_JACOBIAN;
		}

		void initialize(const gsBasis<T>   & basis,
			const index_t,
			const gsOptionList & options,
			gsQuadRule<T>      & rule)
		{
			// Setup Quadrature (harmless slicing occurs)
			rule = gsQuadrature::get(basis, options, side.direction());

			// Set Geometry evaluation flags
			md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
		}

		// Evaluate on element.
		inline void evaluate(const gsBasis<T>       & basis, // to do: more unknowns
			const gsGeometry<T>    & geo,
			// todo: add element here for efficiency
			const gsMatrix<T>      & quNodes)
		{
			md.points = quNodes;
			// Compute the active basis functions
			// Assumes actives are the same for all quadrature points on the current element
			basis.active_into(md.points.col(0), actives);
			const index_t numActive = actives.rows();

			// Evaluate basis values and derivatives on element
			basis.evalAllDers_into(md.points, 1, basisData);

			// Compute geometry related values
			geo.computeMap(md);

			// Evaluate the Neumann data
			neudata_ptr->eval_into(md.values[0], neuData);

			// Initialize local matrix/rhs
			localRhs.setZero(numActive, neudata_ptr->targetDim());
		}

		inline void assemble(gsDomainIterator<T>    &,
			const gsVector<T>      & quWeights)
		{
			gsMatrix<T> & bGrads = basisData[1];
			const index_t numActive = actives.rows();

			gsMatrix<T> w_(numActive, 1);
			for (index_t i = 0; i < numActive; i++)
			{
				w_(i, 0) = w(actives(i), 0);
			}

			for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
			{
				// Compute physical gradients at k as a Dim x NumActive matrix
				transformGradients(md, k, bGrads, pGrads);

				//Compute physical gradient of w
				gsMatrix<T> wGrad = pGrads * w_;

				const T a = pow(eps * eps + (wGrad.transpose() * wGrad).value(), (p - 2) / 2);
				
				// Compute the outer normal vector on the side
				outerNormal(md, k, side, unormal);

				// Multiply quadrature weight by the measure of normal
				const T weight = a * quWeights[k] * unormal.norm();
				
				localRhs.noalias() += weight * basisData[0].col(k) * neuData.col(k).transpose();
			}
		}

		inline void localToGlobal(const index_t patchIndex,
			const std::vector<gsMatrix<T> > &,
			gsSparseSystem<T>               & system)
		{
			// Map patch-local DoFs to global DoFs
			system.mapColIndices(actives, patchIndex, actives);

			// Add contributions to the system matrix and right-hand side
			system.pushToRhs(localRhs, actives, 0);
		}

		void localToGlobal(const gsDofMapper & mapper,
			const gsMatrix<T> & eliminatedDofs,
			const index_t       patchIndex,
			gsSparseMatrix<T> & sysMatrix,
			gsMatrix<T>       & rhsMatrix)
		{
			// Local DoFs to global DoFs
			mapper.localToGlobal(actives, patchIndex, actives);
			const index_t numActive = actives.rows();

			// Push element contribution to the global load vector
			for (index_t j = 0; j != numActive; ++j)
			{
				// convert local DoF index to global DoF index
				const unsigned jj = actives(j);
				if (mapper.is_free_index(jj))
					rhsMatrix.row(jj) += localRhs.row(j);
			}
		}

	public:

		gsMatrix<T> w;
		T eps;
		T p;

	protected:

		index_t subdiv;

		// Neumann function
		const gsFunction<T> * neudata_ptr;
		boxSide side;

		// Basis values
		std::vector<gsMatrix<T> > basisData;
		gsMatrix<T>      pGrads;
		gsMatrix<index_t> actives;

		// Normal and Neumann values
		gsVector<T> unormal;
		gsMatrix<T> neuData;

		// Local matrix and rhs
		gsMatrix<T> localMat;
		gsMatrix<T> localRhs;

		gsMapData<T> md;
	};


} // namespace gismo

