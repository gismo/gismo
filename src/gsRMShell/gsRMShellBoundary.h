/** @file gsRMShellBoundary.h

   @brief Set Boundary condition for the RM shell equation.

   This file is part of the G+Smo library.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): Y. Xia, HS. Wang

   Date:   2020-12-23
*/
#pragma once

#include "gsRMShellBase.h"
#include "gsNURBSinfo.hpp"
#include"gsMyBase/gsMyBase.h"
namespace gismo
{
	/** @brief
		Boundary condition of Reissner-Mindlin shell .

		It sets up an assembler and assembles the system patch wise and combines
		the patch-local stiffness matrices into a global system by various methods
		(see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
		conditions in various ways (see gismo::gsDirichletStrategy).

		\ingroup Assembler
	*/
	template < class T>
	class gsRMShellBoundary 
	{
	
	public:
		// 默认构造函数
		gsRMShellBoundary()	{}
		typedef memory::shared_ptr<gsRMShellBoundary> Ptr;
		typedef memory::unique_ptr<gsRMShellBoundary> uPtr;
		
		// 构造函数
		gsRMShellBoundary(ifstream& bc_stream,
			gsMultiBasis<T> const & m_basis,
			gsNURBSinfo<T> & nurbs_info,
			index_t dof_per_node)
		{
			m_basis_ = m_basis;
			m_nurbs_info = nurbs_info;
			//bc_file = bc_stream;
			sum_nodes = nurbs_info.nurbs_node; //节点总数
			dof_node = dof_per_node;
			//setBoundary(bc_stream);
		}

 		gsRMShellBoundary(gsRMShellBoundary<T>& C)
 		{
 			m_basis_		= C.m_basis_;
			m_nurbs_info	= C.m_nurbs_info;
 			sum_nodes		= C.sum_nodes;
 			dof_node		= C.dof_node;
 // 			bc_file			= C.bc_file;
 // 			m_BCs			= C.m_BCs;
 // 			m_pLoads		= C.m_pLoads;
 // 			m_pPressure		= C.m_pPressure;
			
 		}
		

		// 析构函数
		~gsRMShellBoundary()
		{
			//gsInfo << "You have called Destructor of gsRMShellBoundary!\n";
		}

		// 0. 需要设置全局函数，避免变量被析构
		void setBoundary(ifstream& bc_stream);

		// 1. 读取边界条件
		void _bc_read(gsMultiBasis<T> const& m_basis,
			gsNURBSinfo<T> const& nurbs_info,
			ifstream& bc_file,
			vector<str_2d_SPC1>& vec_SPC1,
			vector<str_2d_FORCE>& vec_FORCE,
			vector<str_2d_PRESSURE>& vec_PRESSURE);

		// 2. 位移约束定义
		void _setBc_spc(vector<str_2d_SPC1>& vec_SPC1,
			vector<gsShellBoundaryCondition>& spcs);

		// 3. 集中载荷定义
		void _setBc_pload(vector<str_2d_FORCE>& vec_FORCE,
			gsPointLoads<real_t>& pLoads);

		// 4. 均布载荷定义
		void _setBc_pressure(vector<str_2d_PRESSURE>& vec_Pressure,
			gsDistriLoads<T>& pressure);

		// 5. 设置集中载荷
		void setpLoad(gsSparseSystem<T>& m_system);

		// 6. 设置分布载荷
		void setPressure(gsSparseSystem<T>& m_system);

		// 7. 设置位移约束
		void setDisp_constra(gsSparseSystem<T>& m_system);

	public:

		gsMultiBasis<T> m_basis_;
		gsNURBSinfo<T>	m_nurbs_info;
		index_t sum_nodes;
		index_t dof_node;
		ifstream bc_file;
		gsPointLoads<T> m_pLoads;
		gsDistriLoads<T> m_pPressure;
		vector<gsShellBoundaryCondition> m_BCs;
	};
}
