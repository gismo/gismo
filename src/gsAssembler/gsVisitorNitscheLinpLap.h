/** @file gsVisitorNitsche.h

@brief Weak (Nitsche-type) BC imposition visitor for elliptic problems.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): A. Mantzaflaris , S. Moore
*/

#pragma once

namespace gismo
{
	/** \brief Visitor for the weak imposition of the dirichlet boundary condition.
	*
	* Adds this term to the bilinear terms
	* \f[ (\nabla u, v)_{\partial \Omega} + (u, \nabla v )_{\partial \Omega}
	*                                     + (\mu*u, v)_{\partial \Omega} \f]
	*
	* The following term is also added to the linear form
	* \f[ (g_D, \mu*v + \nabla v)_{\partial \Omega} \f],
	* where the dirichlet term is given as \f[ g_D \f].
	*/

	template <class T>
	class gsVisitorNitscheLinpLap
	{
	public:

		gsVisitorNitscheLinpLap(const gsPde<T> & pde, const boundary_condition<T> & s)
			: dirdata_ptr(s.function().get()), side(s.side())
		{
			const gsLinpLapPde<T> * pde_ptr = static_cast<const gsLinpLapPde<T>*>(&pde);

			w = pde_ptr->w;
			eps = pde_ptr->eps;
			p = pde_ptr->p;
		}

		/** @brief
		Constructor of the assembler object.

		\param[in] dirdata  is a gsBoundaryConditions object that holds boundary conditions of the form:
		\f[ \text{Dirichlet: } u = g_D \text{ on } \Gamma.\f]
		\f$ v \f$ is the test function and \f$ \Gamma\f$ is the boundary side.
		\param[in] _penalty for inputing the penalty choice
		\param[in] s
		*/
		gsVisitorNitscheLinpLap(const gsFunction<T> & dirdata, T _penalty, boxSide s) :
			dirdata_ptr(&dirdata), penalty(_penalty), side(s)
		{ }

		void initialize(const gsBasis<T> & basis,
			gsQuadRule<T> & rule)
		{
			const int dir = side.direction();
			gsVector<int> numQuadNodes(basis.dim());
			for (int i = 0; i < basis.dim(); ++i)
				numQuadNodes[i] = 4* basis.degree(i) + 1;
			numQuadNodes[dir] = 1;

			// Setup Quadrature
			rule = gsGaussRule<T>(numQuadNodes);// harmless slicing occurs here

												// Set Geometry evaluation flags
			md.flags = NEED_VALUE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM;
		}

		void initialize(const gsBasis<T> & basis,
			const index_t,
			const gsOptionList & options,
			gsQuadRule<T>    & rule)
		{
			// Setup Quadrature (harmless slicing occurs)
			rule = gsQuadrature::get(basis, options, side.direction());

			// Set Geometry evaluation flags
			md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;

			// Compute penalty parameter
			const int deg = basis.maxDegree();
			penalty = (deg + basis.dim()) * (deg + 1) * T(2.5);
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

			// Evaluate the Dirichlet data
			dirdata_ptr->eval_into(md.values[0], dirData);

			// Initialize local matrix/rhs
			localMat.setZero(numActive, numActive);
			localRhs.setZero(numActive, dirdata_ptr->targetDim());
		}

		inline void assemble(gsDomainIterator<T>    & element,
			const gsVector<T>      & quWeights)
		{
			gsMatrix<T> & bGrads = basisData[1];
			const index_t numActive = actives.rows();

			for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
			{

				const typename gsMatrix<T>::Block bVals =
					basisData[0].block(0, k, numActive, 1);

				// Compute the outer normal vector on the side
				outerNormal(md, k, side, unormal);

				// Multiply quadrature weight by the geometry measure
				const T weight = quWeights[k] * unormal.norm();

				// Compute the unit normal vector 
				unormal.normalize();

				// Compute physical gradients at k as a Dim x NumActive matrix
				transformGradients(md, k, bGrads, pGrads);

				//Access the coefficients of w which correspond to the active basis functions
				gsMatrix<T> w_(numActive, 1);
				for (index_t i = 0; i < numActive; i++)
				{
					w_(i, 0) = w(actives(i), 0);
				}

				//Compute the Gradient of the approximative function w by multiplying the coefficients of w with the physical gradients
				//Compute the point evaluation of w by multiplying the coefficients with the basis function evaluations.
				gsMatrix<T> wGrad = pGrads * w_;
				gsMatrix<T> wVal = bVals.transpose() * w_;

				// Get penalty parameter
				const T h = element.getCellSize();
				const T a_pen = pow(eps * eps + ((dirData.col(k) - wVal).transpose() * (dirData.col(k) - wVal)).value() / (h * h), (p - 2) / 2);
				const T mu = a_pen * penalty / h;
				const T a = pow(eps * eps + (wGrad.transpose() * wGrad).value(), (p - 2) / 2);

				//gsInfo << a_pen*mu<< "\n";
				
				//gsInfo << dirData.col(k).rows() << " x "<< dirData.col(k).cols() << "\n";

				// Sum up quadrature point evaluations
				localRhs.noalias() -= weight * ((-mu * a_pen * bVals)
												* dirData.col(k).transpose());

				localMat.noalias() -= weight * (a * bVals * unormal.transpose() * pGrads
								    - mu * a_pen * bVals * bVals.transpose());
			}
		}

		inline void localToGlobal(const index_t                     patchIndex,
			const std::vector<gsMatrix<T> > &,
			gsSparseSystem<T>               & system)
		{
			// Map patch-local DoFs to global DoFs
			system.mapColIndices(actives, patchIndex, actives);

			// Add contributions to the system matrix and right-hand side
			system.pushAllFree(localMat, localRhs, actives, 0);
		}

		void localToGlobal(const gsDofMapper  & mapper,
			const gsMatrix<T>  & eliminatedDofs,
			const index_t        patchIndex,
			gsSparseMatrix<T>  & sysMatrix,
			gsMatrix<T>        & rhsMatrix)
		{
			// Local DoFs to global DoFs
			mapper.localToGlobal(actives, patchIndex, actives);
			const index_t numActive = actives.rows();

			// Push element contribution to the global matrix and load vector
			for (index_t j = 0; j != numActive; ++j)
			{
				const unsigned jj = actives(j);
				rhsMatrix.row(jj) += localRhs.row(j);
				for (index_t i = 0; i != numActive; ++i)
				{
					const unsigned ii = actives(i);
					//                if ( jj <= ii ) // assuming symmetric problem
					sysMatrix(ii, jj) += localMat(i, j);
				}
			}

		}
	public:
		gsMatrix<T> w;
		T eps;
		T p;

	private:
		// Dirichlet function
		const gsFunction<T> * dirdata_ptr;

		// Penalty constant
		T penalty;

		// Side
		boxSide side;

	private:
		// Basis values
		std::vector<gsMatrix<T> > basisData;
		gsMatrix<T>      pGrads;
		gsMatrix<index_t> actives;

		// Normal and Neumann values
		gsVector<T> unormal;
		gsMatrix<T> dirData;

		// Local  matrix and rhs
		gsMatrix<T> localMat;
		gsMatrix<T> localRhs;

		gsMapData<T> md;
	};


} // namespace gismo

