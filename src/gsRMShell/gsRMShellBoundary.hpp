#pragma once
/** @file gsRMShellBoundary.h

   @brief Set Boundary condition for the RM shell equation.

   This file is part of the G+Smo library.

   This Source Code Form is subject to the terms of the Mozilla Public
   License, v. 2.0. If a copy of the MPL was not distributed with this
   file, You can obtain one at http://mozilla.org/MPL/2.0/.

   Author(s): Y. Xia, HS. Wang

   Date:   2020-12-23
*/

#include "gsRMShellBoundary.h"

namespace gismo
{
	// >>> --------------------------------------------------------------------
	template < class T>
	void gsRMShellBoundary<T>::setBoundary(ifstream& bc_stream)
	{
		vector<str_2d_SPC1> vec_SPC1;
		vector<str_2d_FORCE> vec_FORCE;
		vector<str_2d_PRESSURE> vec_PRESSURE;

		//1.read the bc file
		_bc_read(m_basis_, m_nurbs_info, bc_stream,
			vec_SPC1, vec_FORCE, vec_PRESSURE, m_pdDomain);

		//2.set boundary conditions
		_setBc_spc(vec_SPC1, m_BCs);

		//3.set point loads
		_setBc_pload(vec_FORCE, m_pLoads);

		//4.set pressure
		_setBc_pressure(vec_PRESSURE, m_pPressure);
	}

	// >>> --------------------------------------------------------------------
	// 1. 读取边界条件
	template < class T>
	void gsRMShellBoundary<T>::_bc_read(gsMultiBasis<T> const& m_basis,
		gsNURBSinfo<T> const&	nurbsinfo,
		ifstream&				bc_file,
		vector<str_2d_SPC1>&	vec_SPC1,
		vector<str_2d_FORCE>&	vec_FORCE,
		vector<str_2d_PRESSURE>& vec_PRESSURE,
		vector<real_t>&			pdDomain)
	{
		string line;
		stringstream ss_buff;
		string buff_ss_;
		vector<real_t> tem_disp;
		real_t tem_double;
		while (!bc_file.eof())
		{
			getline(bc_file, line);

			// 如果发现了 "!" 就执行
			if (line.find("!") != std::string::npos)
			{
				//cout<<line<<endl;
			}
			// 如果发现了"SPC1" 就执行 
			else if (line.find("SPC1") != std::string::npos)
			{
				str_2d_SPC1 temp_str_SPC1;

				ss_buff.str("");
				ss_buff.str(line);

				ss_buff >> buff_ss_;
				ss_buff >> temp_str_SPC1.no_patch;
				ss_buff >> temp_str_SPC1.position[0];
				ss_buff >> temp_str_SPC1.position[1];
				ss_buff >> temp_str_SPC1.dof;
				ss_buff >> temp_str_SPC1.value;
				ss_buff >> temp_str_SPC1.is_parametric;
				ss_buff.clear();
				vec_SPC1.push_back(temp_str_SPC1);
			}
			else if (line.find("SPCSIDE") != std::string::npos)
			{
				// 使用 str("") 初始化 ss_buff，将 空字符 复制给 ss_buff 的操作
				ss_buff.str("");
				ss_buff.str(line);
				ss_buff >> buff_ss_;		// "SPCSIDE"
				int no_patch;
				ss_buff >> no_patch;		// 0
				double u1[2], u2[2];
				ss_buff >> buff_ss_;		// U1
				ss_buff >> u1[0] >> u1[1];	// 0 0
				ss_buff >> buff_ss_;		// U2
				ss_buff >> u2[0] >> u2[1];	// 0 1
				string dof;					// UX, UY, UZ, RX...
				double va;					// 0.0 value
				bool par;					// parameter = 1
				ss_buff >> dof >> va >> par;

				// 如果约束是在 U1 这条边上
				if (u1[0] != u1[1])
				{
					index_t cp_No = m_nurbs_info.nurbsdata[no_patch].knot_vector[0].knotvec_size;
					for (int i = 0; i < cp_No; ++i)
					{
						str_2d_SPC1 temp_str_SPC1;
						temp_str_SPC1.no_patch = no_patch;
						// 将参数空间分成（每边控制点数 - 1）份，使SPCSIDE 转化为 SPC1
						temp_str_SPC1.position[0] = u1[0] + i * (u1[1] - u1[0]) / (cp_No - 1);
						temp_str_SPC1.position[1] = u2[0];
						temp_str_SPC1.dof = dof;
						temp_str_SPC1.value = va;
						temp_str_SPC1.is_parametric = par;
						vec_SPC1.push_back(temp_str_SPC1);
					}
				}
				else if (u2[0] != u2[1])
				{
					index_t cp_No = m_nurbs_info.nurbsdata[no_patch].knot_vector[1].knotvec_size;
					for (int i = 0; i < cp_No; ++i)
					{
						str_2d_SPC1 temp_str_SPC1;
						temp_str_SPC1.no_patch = no_patch;
						temp_str_SPC1.position[0] = u1[0];
						temp_str_SPC1.position[1] = u2[0] + i * (u2[1] - u2[0]) / (cp_No - 1);
						temp_str_SPC1.dof = dof;
						temp_str_SPC1.value = va;
						temp_str_SPC1.is_parametric = par;
						vec_SPC1.push_back(temp_str_SPC1);
					}
				}
				ss_buff.clear();
			}
			else if (line.find("FORCE") != std::string::npos)
			{
				str_2d_FORCE temp_str_FORCE;
				ss_buff.str("");
				ss_buff.str(line);
				ss_buff >> buff_ss_;						// FORCE
				ss_buff >> temp_str_FORCE.no_patch;			// 0
				ss_buff >> temp_str_FORCE.position[0];		// 1
				ss_buff >> temp_str_FORCE.position[1];		// 1
				ss_buff >> temp_str_FORCE.dof;				// FX,FY,FZ,MX,MY,MZ
				ss_buff >> temp_str_FORCE.value;			// -100
				ss_buff >> temp_str_FORCE.is_parametric;	// 1
				ss_buff.clear();
				vec_FORCE.push_back(temp_str_FORCE);
			}
			else if (line.find("PRESSURE") != std::string::npos)
			{
				str_2d_PRESSURE temp_str_PRESS;
				ss_buff.str("");
				ss_buff.str(line);
				ss_buff >> buff_ss_;						// PRESSURE
				ss_buff >> temp_str_PRESS.no_patch;			// 0
				ss_buff >> temp_str_PRESS.dof;				// PX,PY,PZ,MX,MY,MZ
				ss_buff >> temp_str_PRESS.value;			// -90
				ss_buff >> temp_str_PRESS.section;			// SUR,UP,DOWN,LEFT,RIGHT
				ss_buff >> temp_str_PRESS.is_parametric;	// 1
				ss_buff.clear();
				vec_PRESSURE.push_back(temp_str_PRESS);
			}
			else if (line.find("PDIGADOMAIN") != std::string::npos)
			{
				ss_buff.str("");
				ss_buff.str(line);
				tem_disp.clear();
				ss_buff >> buff_ss_;				// PDIGADOMAIN 
				ss_buff >> buff_ss_;				// Xi 
				ss_buff >> tem_double;				// 0.4
				m_pdDomain.push_back(tem_double);
				ss_buff >> tem_double;				// 0.6 
				m_pdDomain.push_back(tem_double);
				ss_buff >> buff_ss_;				// Eta
				ss_buff >> tem_double;				// 0.4
				m_pdDomain.push_back(tem_double);
				ss_buff >> tem_double;				// 0.6
				m_pdDomain.push_back(tem_double);
				ss_buff.clear();
			}
		}
	}

	// 2. 位移约束
	template < class T>
	void gsRMShellBoundary<T>::_setBc_spc(vector<str_2d_SPC1>& vec_SPC1,
		vector<gsShellBoundaryCondition>& spcs)
	{
		for (unsigned int i = 0; i < vec_SPC1.size(); ++i)
		{
			spcs.push_back(gsShellBoundaryCondition(vec_SPC1[i].no_patch,
				vec_SPC1[i].position[0],
				vec_SPC1[i].position[1]));

			if (vec_SPC1[i].dof == "UX") { spcs[i].value = vec_SPC1[i].value; spcs[i].direction = 0; }
			else if (vec_SPC1[i].dof == "UY") { spcs[i].value = vec_SPC1[i].value; spcs[i].direction = 1; }
			else if (vec_SPC1[i].dof == "UZ") { spcs[i].value = vec_SPC1[i].value; spcs[i].direction = 2; }
			else if (vec_SPC1[i].dof == "RX") { spcs[i].value = vec_SPC1[i].value; spcs[i].direction = 3; }
			else if (vec_SPC1[i].dof == "RY") { spcs[i].value = vec_SPC1[i].value; spcs[i].direction = 4; }
			else if (vec_SPC1[i].dof == "RZ") { spcs[i].value = vec_SPC1[i].value; spcs[i].direction = 5; }
			else if (vec_SPC1[i].dof == "ALL") { spcs[i].value = vec_SPC1[i].value; spcs[i].direction = 6; }
		}
	}

	// 3. 集中载荷的设定
	template < class T>
	void gsRMShellBoundary<T>::_setBc_pload(vector<str_2d_FORCE>& vec_FORCE,
		gsPointLoads<real_t>& pLoads)
	{
		gsInfo << "Sum of Forces: " << vec_FORCE.size() << endl;
		for (unsigned int i = 0; i < vec_FORCE.size(); ++i)
		{
			gsVector<> point(2);
			point[0] = vec_FORCE[i].position[0];
			point[1] = vec_FORCE[i].position[1];

			gsVector<> load(6);
			load << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
			if (vec_FORCE[i].dof == "FX") load(0) = vec_FORCE[i].value;
			else if (vec_FORCE[i].dof == "FY") load(1) = vec_FORCE[i].value;
			else if (vec_FORCE[i].dof == "FZ") load(2) = vec_FORCE[i].value;
			else if (vec_FORCE[i].dof == "MX") load(3) = vec_FORCE[i].value;
			else if (vec_FORCE[i].dof == "MY") load(4) = vec_FORCE[i].value;
			else if (vec_FORCE[i].dof == "MZ") load(5) = vec_FORCE[i].value;
			pLoads.addLoad(point, load, vec_FORCE[i].no_patch, vec_FORCE[i].is_parametric);
		}
		
	}

	// 4. 均布载荷的设定
	template < class T>
	void gsRMShellBoundary<T>::_setBc_pressure(vector<str_2d_PRESSURE>& vec_Pressure,
		gsDistriLoads<T>& pressure)
	{
		gsInfo << "Sum of distributed loads: " << vec_Pressure.size() << endl;
		//Distributed loads
		for (unsigned int i = 0; i < vec_Pressure.size(); ++i)
		{
			gsVector<> load(6); load << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
			
			if (vec_Pressure[i].dof == "PX") load(0) = vec_Pressure[i].value;
			else if (vec_Pressure[i].dof == "PY") load(1) = vec_Pressure[i].value;
			else if (vec_Pressure[i].dof == "PZ") load(2) = vec_Pressure[i].value;
			else if (vec_Pressure[i].dof == "MX") load(3) = vec_Pressure[i].value;
			else if (vec_Pressure[i].dof == "MY") load(4) = vec_Pressure[i].value;
			else if (vec_Pressure[i].dof == "MZ") load(5) = vec_Pressure[i].value;
			pressure.addLoad(load, vec_Pressure[i].no_patch, vec_Pressure[i].section, vec_Pressure[i].is_parametric);
		}
	}

	// >>> --------------------------------------------------------------------
	// 5. 设置均布载荷
	template < class T>
	void gsRMShellBoundary<T>::setPressure(gsSparseSystem<T>& m_system)
	{
		index_t patch_no = 0;
		index_t nodes_row = m_nurbs_info.nurbsdata[patch_no].knot_vector[0].knotvec_size;	// 每行控制点数
		index_t nodes_col = m_nurbs_info.nurbsdata[patch_no].knot_vector[1].knotvec_size;	// 每列控制点数

		// 继续向k文件中添加载荷等数据
		// 显示指定app模式，防止已有数据被丢弃
		ofstream append("../../TestResult/RMShell.k", ofstream::app);
		append << "*LOAD_NODE_POINT\n";

		for (size_t i = 0; i < m_pPressure.numLoads(); ++i)
		{
			if (m_pPressure[i].section == "SUR")
			{
				for (index_t k = 0; k < sum_nodes; ++k)
				{
					for (index_t m = 0; m < dof_node; ++m)
					{
						m_system.rhs()(k * dof_node + m, i) *= m_pPressure[i].value[m];
					}					
				}
			}
			else if (m_pPressure[i].section == "DOWN")
			{
				for (index_t k = 0; k < sum_nodes; ++k)
				{
					if (k < nodes_row)
					{
						for (index_t m = 0; m < dof_node; ++m)
						{
							m_system.rhs()(k * dof_node + m, i) *= m_pPressure[i].value[m];
						}
					}
					else
					{
						for (index_t m = 0; m < dof_node; ++m)
						{
							m_system.rhs()(k * dof_node + m, i) = 0;
						}
					}
				}
			}
			else if (m_pPressure[i].section == "UP")
			{
				for (index_t k = 0; k < sum_nodes; ++k)
				{
					if ((nodes_col-1)*nodes_row < k)
					{
						for (index_t m = 0; m < dof_node; ++m)
						{
							m_system.rhs()(k * dof_node + m, i) *= m_pPressure[i].value[m];
						}
					}
					else
					{
						for (index_t m = 0; m < dof_node; ++m)
						{
							m_system.rhs()(k * dof_node + m, i) = 0;
						}
					}
				}
			}
			else if (m_pPressure[i].section == "LEFT")
			{
				for (index_t k = 0; k < sum_nodes; ++k)
				{
					if ( k % nodes_row == 0)
					{
						for (index_t m = 0; m < dof_node; ++m)
						{
							m_system.rhs()(k * dof_node + m, i) *= m_pPressure[i].value[m];
						}
					}
					else
					{
						for (index_t m = 0; m < dof_node; ++m)
						{
							m_system.rhs()(k * dof_node + m, i) = 0;
						}
					}
				}
			}
			else if (m_pPressure[i].section == "RIGHT")
			{
				for (index_t k = 0; k < sum_nodes; ++k)
				{
					if ((k+1) % nodes_row == 0 )
					{
						for (index_t m = 0; m < dof_node; ++m)
						{
							m_system.rhs()(k * dof_node + m, i) *= m_pPressure[i].value[m];
							
						}
					}
					else
					{
						for (index_t m = 0; m < dof_node; ++m)
						{
							m_system.rhs()(k * dof_node + m, i) = 0;
						}
					}

				}
			}
		}

		// 如果存在多个均布载荷就合并为1个
		if (1 < m_pPressure.numLoads())
		{
			for (size_t i = 1; i < m_pPressure.numLoads(); ++i)
			{
				m_system.rhs().col(0) += m_system.rhs().col(i);
			}
		}

		// 输出k文件中的均布载荷
		for (index_t k = 0; k < sum_nodes; ++k)
		{
			for (index_t m = 0; m < dof_node; ++m)
			{
				if (m_system.rhs()(k * dof_node + m, 0) != 0)
				{
					append << setw(10) << k + 1
						<< setw(10) << m + 1
						<< setw(10) << 0
						<< setw(10) << m_system.rhs()(k * dof_node + m, 0)
						<< setw(10) << 0 << "\n";
				}	
			}
		}
		append.close();

		if (m_pPressure.numLoads()==0)
		{ 
			m_system.rhs().setZero();
		}

	}

	// 6. 设置集中载荷
	template < class T>
	void gsRMShellBoundary<T>::setpLoad(gsSparseSystem<T>& m_system)
	{
		gsMatrix<T>			bVals;
		gsMatrix<index_t>	acts;

		// 继续向k文件中添加载荷等数据
		// 显示指定app模式，防止已有数据被丢弃
		ofstream append("../../TestResult/RMShell.k", ofstream::app);
		append << "*LOAD_NODE_POINT\n";

		for (size_t i = 0; i < m_pLoads.numLoads(); ++i)
		{
			//m_patches
			if (m_pLoads[i].parametric)
			{
				m_basis_.basis(m_pLoads[i].patch).active_into(m_pLoads[i].point, acts);
				m_basis_.basis(m_pLoads[i].patch).eval_into(m_pLoads[i].point, bVals);
			}
			else
			{
				gsWarn << "Point loads parametric for now.\n";
			}

			// translate patch-local indices to global dof indices
			for (index_t k = 0; k < acts.rows(); ++k)
			{
				for (index_t m = 0; m < dof_node; ++m)
				{
					m_system.rhs()(acts(k, 0) * dof_node + m, 0) +=
						bVals(k, 0) * m_pLoads[i].value[m];
					
					// k 文件 载荷数据
					if (bVals(k, 0)!=0 && m_pLoads[i].value[m]!=0)
					{
						append << setw(10) << acts(k, 0) + 1
							<< setw(10) << m + 1
							<< setw(10) << 0
							<< setw(10) << bVals(k, 0) * m_pLoads[i].value[m]
							<< setw(10) << 0 << "\n";
					}
					
				}
			}
		}
		append.close();

	}
	
	// 7. 设置位移约束
	template < class T>
	void gsRMShellBoundary<T>::setDisp_constra(gsSparseSystem<T>& m_system)
	{
		gsMatrix<T>			bVals;
		gsMatrix<index_t>	acts;
		const double		LARGE_NUMBER_MINE = 6.0e15;
		
		vector<set<index_t> > vec_boundedCp(7, set<index_t>());
		set<index_t> boundedCp;

		// 继续向k文件中添加位移约束等数据
		// 显示指定app模式，防止已有数据被丢弃
		ofstream append("../../TestResult/RMShell.k", ofstream::app);
		append << "*BOUNDARY_SPC_NODE\n";

		for (size_t i = 0; i < m_BCs.size(); ++i)
		{
			if (m_BCs[i].is_parametric)
			{
				// m_bases.front() 会返回对 m_bases 中第一个元素的引用
				// .basis() {return *m_bases[i]};
				// 所以这个地方其实是要根据 m_BCs.point 找出与之对应的控制点和基函数值吧，
				m_basis_.front().basis(m_BCs[i].patch).active_into(m_BCs[i].point, acts);
				m_basis_.front().basis(m_BCs[i].patch).eval_into(m_BCs[i].point, bVals);
			}
			else
			{
				gsWarn << "Boundary conditions parametric for now.\n";
			}

			// translate patch-local indices to global dof indices
			index_t tem = 0;
			//multiply large number method
			for (index_t k = 0; k < acts.rows(); ++k)
			{
				// (vec_SPC1[i].dof == "ALL") {spcs[i].direction = 6;}
				if (m_BCs[i].direction == 6 && bVals(k, 0) != 0
					&& (boundedCp.find(acts(k, 0)) == boundedCp.end()))
				{
					for (int m = 0; m < 6; ++m)
					{
						// 约束点约束方向对应的全局K阵对角线坐标
						tem = acts(k, 0) * dof_node + m; 
						// 将P_j 用 alpha * K_jj * a_j 代替，
						// alpha是大数，K_jj是刚度阵对角元，a_j是约束位移值
						m_system.rhs()(tem, 0) = m_system.matrix()(tem, tem)
											* LARGE_NUMBER_MINE * m_BCs[i].value;
						// 对角元素乘以大数
						m_system.matrix()(tem, tem) *= LARGE_NUMBER_MINE;
					}
					// 防止一个点重复乘大数
					boundedCp.insert(acts(k, 0));

					// k 文件 位移约束
					append << setw(10) << acts(k, 0) + 1
						<< setw(10) << " "
						<< setw(10) << 1 << setw(10) << 1 << setw(10) << 1
						<< setw(10) << 1 << setw(10) << 1 << setw(10) << 1 << "\n";
				}
				else if (m_BCs[i].direction != 6 && bVals(k, 0) != 0
					&& (vec_boundedCp[m_BCs[i].direction].find(acts(k, 0))
						== vec_boundedCp[m_BCs[i].direction].end()))
				{
					tem = acts(k, 0) * dof_node + m_BCs[i].direction;
					
					m_system.rhs()(tem, 0) = m_system.matrix()(tem, tem) 
											* LARGE_NUMBER_MINE * m_BCs[i].value;
					
					m_system.matrix()(tem, tem) *= LARGE_NUMBER_MINE;
					// 防止重复对同一个乘成大数
					vec_boundedCp[m_BCs[i].direction].insert(acts(k, 0));

					// k 文件 位移约束
					append << setw(10) << acts(k, 0) + 1
						<< setw(10) << " ";
						for (index_t d = 0; d < dof_node; ++d)
						{
							if(d == m_BCs[i].direction)
								append << setw(10) << 1;
							else
								append << setw(10) << 0;
						}
					append << "\n";
				}
			}
		}

		append << "*END";
		append.close();
	}
	// -------------------------------------------------------------------- <<<
}