/** @file gsRMShellAssembler.h

   @brief Provides assembler for the RM shell equation.

   This file is part of the G+Smo library.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): Y. Xia, HS. Wang

   Date:   2020-12-23
*/
#pragma once

#include <gsAssembler/gsAssembler.h>
#include "gsRMShellPDE.h"
//#include "gsRMShellBase.h"
#include "gsRMShellBoundary.hpp"
#include"gsMyBase/gsMyBase.h"
namespace gismo
{
/** @brief
	Implementation of Reissner-Mindlin shell assembler.

	It sets up an assembler and assembles the system patch wise and combines
	the patch-local stiffness matrices into a global system by various methods
	(see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
	conditions in various ways (see gismo::gsDirichletStrategy).

	\ingroup Assembler
*/
template < class T>
class gsRMShellAssembler : public gsAssembler < T>
{
public:
	typedef gsAssembler<T> Base;

public:
	// 默认构造函数
	gsRMShellAssembler()
	{}

	/** @brief Main Constructor of the assembler object.

	\param[in] pde A boundary value RM shell problem
	\param[in] bases a multi-basis that contains patch-wise bases
	*/
	// 构造函数1
	gsRMShellAssembler(const gsRMShellPde<T>& pde,
					   const gsMultiBasis<T>& bases)
	{
		Base::initialize(pde, bases, m_options);
	}

	/** @brief Main Constructor of the assembler object.

   \param[in] pde A boundary value RM shell problem
   \param[in] bases a multi-basis that contains patch-wise bases
   \param[in] dirStrategy option for the treatment of Dirichlet boundary
   \param[in] intStrategy option for the treatment of patch interfaces
   */
	// 构造函数2
	gsRMShellAssembler(const gsRMShellPde<T>& pde,
					   const gsMultiBasis<T>& bases,
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
	\param[in] rhs is the right-hand side of the RM shell equation.
	\param[in] dirStrategy option for the treatment of Dirichlet boundary
	\param[in] intStrategy option for the treatment of patch interfaces

	\ingroup Assembler
	*/
	// 构造函数3
	gsRMShellAssembler(gsMultiPatch<T> const& patches,
	                   gsMultiBasis<T> const& basis,
	                   gsBoundaryConditions<T> const& bconditions,
	                   const gsFunction<T> & rhs,
	                   dirichlet::strategy           dirStrategy = dirichlet::none,
	                   iFace::strategy               intStrategy = iFace::none)
	{
	     m_options.setInt("DirichletStrategy", dirStrategy);
	     m_options.setInt("InterfaceStrategy", intStrategy);
	
	     typename gsPde<T>::Ptr pde(new gsRMShellPde<T>(patches, bconditions, rhs));
	     Base::initialize(pde, basis, m_options);
	}

	/** @brief 
	Constructor of the assembler object.

	\param[in] patches is a gsMultiPatch object describing the geometry.
	\param[in] basis a multi-basis that contains patch-wise bases
	\param[in] bconditions is a gsBoundaryConditions object that holds all boundary conditions.
	\param[in] rhs is the right-hand side of the RM shell equation.
	\param[in] dirStrategy option for the treatment of Dirichlet boundary
	\param[in] intStrategy option for the treatment of patch interfaces

	\ingroup Assembler
	*/
	// 构造函数4
	gsRMShellAssembler(gsMultiPatch<T> const& patches,
		gsMultiBasis<T> const& basis,
		Material &m_material,
		gsRMShellBoundary<T>  & BoundaryCondition,
		ifstream& bc_stream,
		const gsFunction<T>& rhs,
		gsBoundaryConditions<T> const& bconditions,
		dirichlet::strategy     dirStrategy = dirichlet::none,
		iFace::strategy         intStrategy = iFace::none,
		index_t dof = 6,
		index_t intept = 0): m_BoundaryCondition(BoundaryCondition)
	{
		m_patches		= patches;
		m_dim			= patches.parDim(); // 参数域的维度
		m_basis			= basis;
		m_sum_nodes		= basis.size(); // 总控制点数
		material		= m_material;
		m_dof_per_node	= dof;
		m_inte			= intept; // 增加的积分点数
		
		m_options.setInt("DirichletStrategy", dirStrategy);
		m_options.setInt("InterfaceStrategy", intStrategy);
	
		typename gsPde<T>::Ptr pde(new gsRMShellPde<T>(m_patches, bconditions, rhs));
		Base::initialize(pde, basis, m_options);	

		m_BoundaryCondition.setBoundary(bc_stream);
		//setBoundary(bc_stream);
	}
	/*
	virtual gsAssembler<T> * clone() const
	{
	     return new gsRMShellAssembler<T>(*this);
	}
	virtual gsAssembler<T> * create() const
	{
	     return new gsRMShellAssembler<T>();
	}
	*/
	// Refresh routine
	virtual void refresh();
	
	// Main assembly routine
	virtual void assemble();
	
	virtual void constructSolution(gsMatrix<T>& solVector,
		gsMultiPatch<T>& result);

	Eigen::SparseSelfAdjointView< typename gsSparseMatrix<T>::Base, Lower> fullMatrix()
	{
	     return m_system.matrix().template selfadjointView<Lower>();
	}

	

public:
	Material material;
	index_t m_sum_nodes;
	index_t m_dof_per_node;
	index_t m_inte;
	index_t m_dim;
	ifstream m_bcfile;
	gsMultiPatch<T> m_patches;
	gsMultiBasis<T> m_basis;
	gsNURBSinfo<T>  m_nurbs_info;
	gsRMShellBoundary<T> m_BoundaryCondition;

protected:

	// Members from gsAssembler
	using Base::m_pde_ptr;
	using Base::m_bases;
	using Base::m_ddof;
	using Base::m_options;
	using Base::m_system;

 }; // class gsRMShellAssembler 

} // namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsRMShellAssembler.hpp)
#endif