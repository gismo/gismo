/** @file gsVisitorLinpLap.h
@brief LinpLap equation element visitor.
This file is part of the G+Smo library.
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
Author(s):
*/

#pragma once

#include <gsAssembler/gsQuadrature.h>
#include <gsAssembler/gsSubdividedRule.h>

namespace gismo
{

	/** \brief Visitor for the linearized p-Laplace equation.
	*
	* Assembles the bilinear terms
	* \f[ (\nabla u,\nabla v)_\Omega \text{ and } (f,v)_\Omega \f]
	* For \f[ u = g \quad on \quad \partial \Omega \f],
	*
	*/

	template <class T, bool paramCoef = false>
	class gsVisitorLinpLap
	{
	public:

		/** \brief Constructor for gsVisitorLinpLap.
		*/
		gsVisitorLinpLap(const gsPde<T> & pde,T eps_R, index_t subdiv_ = 1, bool prec_ = 0)
		{
			pde_ptr = static_cast<const gsLinpLapPde<T>*>(&pde);
			subdiv = subdiv_;
			prec = prec_;
			eps_ = eps_R;
		}

		void initialize(const gsBasis<T> & basis,
			const index_t patchIndex,
			const gsOptionList & options,
			gsQuadRule<T>    & rule)
		{
			// Grab right-hand side for current patch
			rhs_ptr = &pde_ptr->rhs()->piece(patchIndex);

			// Setup Quadrature
			//rule = gsQuadrature::get(basis, options); // harmless slicing occurs here
			//TODO: make this configurable
			rule = gsSubdividedRule<T, gsQuadRule<T> >(gsQuadrature::get(basis, options), subdiv);

			// Set Geometry evaluation flags
			md.flags = NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM;
		}

		// Evaluate on element.
		inline void evaluate(const gsBasis<T>       & basis,
			const gsGeometry<T>    & geo,
			const gsMatrix<T>      & quNodes)
		{
			md.points = quNodes;
			// Compute the active basis functions
			// Assumes actives are the same for all quadrature points on the elements
			basis.active_into(md.points.col(0), actives);
			numActive = actives.rows();

			// Evaluate basis functions on element
			basis.evalAllDers_into(md.points, 1, basisData);

			// Compute image of Gauss nodes under geometry mapping as well as Jacobians
			geo.computeMap(md);

			// Evaluate right-hand side at the geometry points paramCoef
			// specifies whether the right hand side function should be
			// evaluated in parametric(true) or physical (false)
			rhs_ptr->eval_into((paramCoef ? md.points : md.values[0]), rhsVals);

			// Initialize local matrix/rhs
			localMat.setZero(numActive, numActive);
			localRhs.setZero(numActive, rhsVals.rows());//multiple right-hand sides
			localJ = 0;
		}

		inline void assemble(gsDomainIterator<T>    &element,
			gsVector<T> const      & quWeights)
		{

			gsMatrix<T> & bVals = basisData[0];
			gsMatrix<T> & bGrads = basisData[1];

			//Access the coefficients of w which correspond to the active basis functions
			gsMatrix<T> w_(numActive, 1);

			const T h = element.getCellSize();

			for (index_t i = 0; i < numActive; i++)
			{
				w_(i, 0) = pde_ptr->w(actives(i), 0);
			}

			real_t eps = pde_ptr->eps; //simplicity

			if (prec)      //elementwise regularization
			{
				real_t sum = 0;
				for (index_t k = 0; k < quWeights.rows(); ++k)
				{
					const T weight = quWeights[k] * md.measure(k);
					transformGradients(md, k, bGrads, physGrad);
					gsMatrix<T> wGrad = physGrad * w_;
					sum += weight * (wGrad.transpose() * wGrad).value();
				}
				if (sum < 0.1*h*h) //good value for the gradient?
				{
					eps = math::max(eps_, eps);
					//gsInfo<<"yes \n";
				}
			}

			for (index_t k = 0; k < quWeights.rows(); ++k) // loop over quadrature nodes
			{
				// Multiply weight by the geometry measure
				const T weight = quWeights[k] * md.measure(k);

				// Compute physical gradients at k as a Dim x NumActive matrix
				transformGradients(md, k, bGrads, physGrad);

				//Compute the Gradient of the approximative function w by multiplying the coefficients of w with the physical gradients
				gsMatrix<T> wGrad = physGrad * w_;
				gsMatrix<T> wVal = bVals.col(k).transpose() * w_;

				eps = pde_ptr->eps;

				T a = pow(eps * eps + (wGrad.transpose() * wGrad).value(), (pde_ptr->p - 2) / 2);

				localJ += weight * (pow(eps * eps + (wGrad.transpose() * wGrad).value(), (pde_ptr->p) / 2) / (pde_ptr->p) - (rhsVals.col(k).transpose()*wVal).value());
				localRhs.noalias() += weight * (bVals.col(k) * rhsVals.col(k).transpose());
				localMat.noalias() += weight * (a * (physGrad.transpose() * physGrad) + pde_ptr->lambda * pow(wVal.value()*wVal.value(),pde_ptr->alpha/2) * bVals.col(k) * bVals.col(k).transpose());
			}
		}

		inline void localToGlobal(const index_t                     patchIndex,
			const std::vector<gsMatrix<T> > & eliminatedDofs,
			gsSparseSystem<T>               & system,
			T								& J)
		{
			// Map patch-local DoFs to global DoFs
			system.mapColIndices(actives, patchIndex, actives);

			// Add contributions to the system matrix and right-hand side
			system.push(localMat, localRhs, actives, eliminatedDofs.front(), 0, 0);

			//Add contributions to the evaluation of the energy functional in w
			J = J + localJ;
		}

		T eps_;

	protected:
		//number of subdivisions for quadrature
		index_t subdiv;
		bool prec;

	protected:
		// Pointer to the pde data
		const gsLinpLapPde<T> * pde_ptr;

	protected:
		// Basis values
		std::vector<gsMatrix<T> > basisData;
		gsMatrix<T>        physGrad;
		gsMatrix<index_t> actives;
		index_t numActive;

	protected:
		// Right hand side ptr for current patch
		const gsFunction<T> * rhs_ptr;

		// Local values of the right hand side
		gsMatrix<T> rhsVals;

	protected:
		// Local matrices
		gsMatrix<T> localMat;
		gsMatrix<T> localRhs;
		T localJ;

		gsMapData<T> md;
	};


} // namespace gismo