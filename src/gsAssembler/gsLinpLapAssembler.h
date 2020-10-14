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
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
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

	public:

		gsLinpLapAssembler()
		{ }

		/** @brief Main Constructor of the assembler object.

		\param[in] pde A boundary value Poisson problem
		\param[in] bases a multi-basis that contains patch-wise bases
		*/
		gsLinpLapAssembler(const gsLinpLapPde<T>          & pde,
			const gsMultiBasis<T>          & bases)
		{
			Base::initialize(pde, bases, m_options);
		}


		/** @brief Main Constructor of the assembler object.

		\param[in] pde A boundary value LinpLap problem
		\param[in] bases a multi-basis that contains patch-wise bases
		\param[in] dirStrategy option for the treatment of Dirichlet boundary
		\param[in] intStrategy option for the treatment of patch interfaces
		*/
		gsLinpLapAssembler(const gsLinpLapPde<T>          & pde,
			const gsMultiBasis<T>          & bases,
			dirichlet::strategy           dirStrategy,
			iFace::strategy               intStrategy = iFace::glue)
		{
			m_options.setInt("DirichletStrategy", dirStrategy);
			m_options.setInt("InterfaceStrategy", intStrategy);

			Base::initialize(pde, bases, m_options);
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
			const gsMatrix<T> &w,
			dirichlet::strategy           dirStrategy = dirichlet::elimination,
			iFace::strategy               intStrategy = iFace::glue)
		{
			m_options.setInt("DirichletStrategy", dirStrategy);
			m_options.setInt("InterfaceStrategy", intStrategy);

			typename gsPde<T>::Ptr pde(new gsLinpLapPde<T>(patches, bconditions, rhs, eps, p, w));
			Base::initialize(pde, basis, m_options);
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

		/// Returns an expression of the "full" assembled sparse
		/// matrix. Note that matrix() might return a lower diagonal
		/// matrix, if we exploit possible symmetry during assembly
		/// (check: m_matrix.symmetry() == true )
		Eigen::SparseSelfAdjointView< typename gsSparseMatrix<T>::Base, Lower> fullMatrix()
		{
			return m_system.matrix().template selfadjointView<Lower>();
		}

	protected:

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

		// Assemble volume integrals
		Base::template push<gsVisitorLinpLap<T> >();

		// Enforce Neumann boundary conditions
		Base::template push<gsVisitorNeumann<T> >(m_pde_ptr->bc().neumannSides());

		switch (m_options.getInt("DirichletStrategy"))
		{
		case dirichlet::penalize:
			Base::penalizeDirichletDofs();
			break;
		case dirichlet::nitsche:
			Base::template push<gsVisitorNitscheLinpLap<T> >(m_pde_ptr->bc().dirichletSides());
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
