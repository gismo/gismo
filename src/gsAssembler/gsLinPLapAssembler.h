/** @file gsLinpLapAssembler.h
@brief Provides assembler for the linearized p-Laplace equation.
This file is part of the G+Smo library.
This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.
Author(s):
*/

#pragma once

#include <gsAssembler/gsAssembler.h>
#include <gsPde/gsLinpLapPde.h>

#include <gsAssembler/gsVisitorLinpLap.h> // Stiffness volume integrals
#include <gsAssembler/gsVisitorNeumannLinpLap.h> 
#include <gsAssembler/gsVisitorNitscheLinpLap.h> // Nitsche boundary integrals
#include <gsAssembler/gsVisitorDg.h>      // DG interface integrals


namespace gismo
{

	/** @brief
	Implementation of an (multiple right-hand side) linearized p-Laplace assembler.
	The LinpLap equation: \f$-\Delta\mathbf{u}=\mathbf{f} \f$
	It sets up an assembler and assembles the system patch wise and combines
	the patch-local stiffness matrices into a global system by various methods
	(see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
	conditions in various ways (see gismo::gsDirichletStrategy).
	\ingroup Assembler
	*/
	template <class T>
	class gsLinpLapAssembler : public gsAssembler<T>
	{
	public:
		typedef gsAssembler<T> Base;

		typedef typename gsBoundaryConditions<T>::bcContainer bcContainer;

	public:

		gsLinpLapAssembler()
		{ }

		/** @brief Main Constructor of the assembler object.
		\param[in] pde A boundary value Poisson problem
		\param[in] bases a multi-basis that contains patch-wise bases
		*/
		gsLinpLapAssembler(const gsLinpLapPde<T>          & pde,
			const gsMultiBasis<T>          & bases,
			T eps_R,
			index_t						   subdiv_ = 1,
			bool						   prec_ = 0)
		{
			Base::initialize(pde, bases, m_options);
			J = 0;
			subdiv = subdiv_;
			prec = prec_;
			eps_=eps_R;
		}


		/** @brief Main Constructor of the assembler object.
		\param[in] pde A boundary value LinpLap problem
		\param[in] bases a multi-basis that contains patch-wise bases
		\param[in] dirStrategy option for the treatment of Dirichlet boundary
		\param[in] intStrategy option for the treatment of patch interfaces
		*/
		gsLinpLapAssembler(const gsLinpLapPde<T>          & pde,
			T eps_R,
			const gsMultiBasis<T>          & bases,
			dirichlet::strategy           dirStrategy,
			iFace::strategy               intStrategy = iFace::glue,
			index_t						  subdiv_ = 1,
			bool						  prec_ = 0)
		{
			m_options.setInt("DirichletStrategy", dirStrategy);
			m_options.setInt("InterfaceStrategy", intStrategy);

			Base::initialize(pde, bases, m_options);
			J = 0;
			subdiv = subdiv_;
			prec = prec_;
			eps_ = eps_R;
		}

		/** @brief
		Constructor of the assembler object.
		\param[in] patches is a gsMultiPatch object describing the geometry.
		\param[in] basis a multi-basis that contains patch-wise bases
		\param[in] bconditions is a gsBoundaryConditions object that holds all boundary conditions.
		\param[in] rhs is the right-hand side of the linearized p-Laplace equation, \f$\mathbf{f}\f$.
		\param[in] dirStrategy option for the treatment of Dirichlet boundary
		\param[in] intStrategy option for the treatment of patch interfaces
		\ingroup Assembler
		*/
		gsLinpLapAssembler(gsMultiPatch<T> const         & patches,
			gsMultiBasis<T> const         & basis,
			gsBoundaryConditions<T> const & bconditions,
			const gsFunction<T>           & rhs,
			const T &eps,
			const T &p,
			T eps_R,
			const gsMatrix<T> &w,
			dirichlet::strategy           dirStrategy = dirichlet::elimination,
			iFace::strategy               intStrategy = iFace::glue,
			index_t						  subdiv_ = 1,
			bool						  prec_ = 0)
		{
			m_options.setInt("DirichletStrategy", dirStrategy);
			m_options.setInt("InterfaceStrategy", intStrategy);

			typename gsPde<T>::Ptr pde(new gsLinpLapPde<T>(patches, bconditions, rhs, eps, p, w));
			Base::initialize(pde, basis, m_options);
			J = 0;
			eps_ = eps_R;
			subdiv = subdiv_;
			prec = prec_;
		}

		void initialize(const gsPde<T>           & pde,
			T eps_R,
			const gsMultiBasis<T>    & bases,
			const gsOptionList & opt = Base::defaultOptions(),
			index_t subdiv_ = 1,
			bool prec_ = 0)
		{
			typename gsPde<T>::Ptr _pde = memory::make_shared_not_owned(&pde);
			Base::initialize(_pde, bases, opt);
			subdiv = subdiv_;
			eps_ = eps_R,
			prec = prec_;
		}

		virtual gsAssembler<T>* clone() const
		{
			return new gsLinpLapAssembler<T>(*this);
		}

		virtual gsAssembler<T>* create() const
		{
			return new gsLinpLapAssembler<T>();
		}

		// Refresh routine
		virtual void refresh();

		// Main assembly routine
		virtual void assemble();


	protected:

		template<class ElementVisitor>
		void push1(const ElementVisitor & visitor)
		{
			for (size_t np = 0; np < m_pde_ptr->domain().nPatches(); ++np)
			{
				ElementVisitor curVisitor = visitor;
				//Assemble (fill m_matrix and m_rhs) on patch np
				apply1(curVisitor, np);
			}
		}

		template<class BElementVisitor>
		void push1(const bcContainer & BCs)
		{
			for (typename bcContainer::const_iterator it = BCs.begin(); it != BCs.end(); ++it)
			{
				BElementVisitor visitor(*m_pde_ptr, *it, subdiv);
				//Assemble (fill m_matrix and m_rhs) contribution from this BC
				Base::apply(visitor, it->patch(), it->side());
			}
		}


		template<class ElementVisitor>
		void apply1(ElementVisitor & visitor,
			size_t patchIndex,
			boxSide side = boundary::none)
		{
			//gsDebug<< "Apply to patch "<< patchIndex <<"("<< side <<")\n";

			const gsBasisRefs<T> bases(m_bases, patchIndex);

			gsQuadRule<T> quRule; // Quadrature rule
			gsMatrix<T> quNodes; // Temp variable for mapped nodes
			gsVector<T> quWeights; // Temp variable for mapped weights

								   // Initialize reference quadrature rule and visitor data
			visitor.initialize(bases, patchIndex, m_options, quRule);

			const gsGeometry<T> & patch = m_pde_ptr->patches()[patchIndex];

			// Initialize domain element iterator -- using unknown 0
			typename gsBasis<T>::domainIter domIt = bases[0].makeDomainIterator(side);

			// Start iteration over elements
			//for (domIt->next(tid); domIt->good(); domIt->next(nt))

			for (; domIt->good(); domIt->next())
			{
				// Map the Quadrature rule to the element
				quRule.mapTo(domIt->lowerCorner(), domIt->upperCorner(), quNodes, quWeights);

				// Perform required evaluations on the quadrature nodes
				visitor.evaluate(bases, patch, quNodes);

				// Assemble on element
				visitor.assemble(*domIt, quWeights);

				// Push to global matrix and right-hand side vector
				visitor.localToGlobal(patchIndex, m_ddof, m_system, J); // omp_locks inside
			}
		}

		/// Returns an expression of the "full" assembled sparse
		/// matrix. Note that matrix() might return a lower diagonal
		/// matrix, if we exploit possible symmetry during assembly
		/// (check: m_matrix.symmetry() == true )
		Eigen::SparseSelfAdjointView< typename gsSparseMatrix<T>::Base, Lower> fullMatrix()
		{
			return m_system.matrix().template selfadjointView<Lower>();
		}

	public:

		T energy() { return J; }

	protected:

		T J;
		index_t subdiv;
		bool prec;
		T eps_;

		// Members from gsAssembler
		using Base::m_pde_ptr;
		using Base::m_bases;
		using Base::m_ddof;
		using Base::m_options;
		using Base::m_system;
	};

	template<class T>
	void gsLinpLapAssembler<T>::refresh()
	{
		// We use predefined helper which initializes the system matrix
		// rows and columns using the same test and trial space
		Base::scalarProblemGalerkinRefresh();
	}

	template<class T>
	void gsLinpLapAssembler<T>::assemble()
	{
		GISMO_ASSERT(m_system.initialized(),
			"Sparse system is not initialized, call initialize() or refresh()");

		// Reserve sparse system
		m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

		// Compute the Dirichlet Degrees of freedom (if needed by m_options)
		Base::computeDirichletDofs();

		// Clean the sparse system
		// m_system.setZero(); //<< this call leads to a quite significant performance degrade!
		J = 0;
		// Assemble volume integrals
		push1<gsVisitorLinpLap<T> >(gsVisitorLinpLap<T>(*m_pde_ptr, eps_, subdiv, prec));

		// Enforce Neumann boundary conditions
		push1<gsVisitorNeumannLinpLap<T>>(m_pde_ptr->bc().neumannSides());

		switch (m_options.getInt("DirichletStrategy"))
		{
		case dirichlet::penalize:
			Base::penalizeDirichletDofs();
			break;
		case dirichlet::nitsche:
			push1<gsVisitorNitscheLinpLap<T> >(m_pde_ptr->bc().dirichletSides());
			break;
		default:
			break;
		}

		if (m_options.getInt("InterfaceStrategy") == iFace::dg)
			Base::template pushInterface<gsVisitorDg<T> >();

		// Assembly is done, compress the matrix
		Base::finalize();
	}
} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsPoissonAssembler.hpp)
#endif