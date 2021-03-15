#pragma once
/** @file gsRMShellAssembler.hpp

	@brief Provides assembler implementation for the RM shell equation.

	This file is part of the G+Smo library.

	This Source Code Form is subject to the terms of the Mozilla Public
	License, v. 2.0. If a copy of the MPL was not distributed with this
	file, You can obtain one at http://mozilla.org/MPL/2.0/.

	Author(s): Y. Xia, HS. Wang

	Date:   2020-12-23
*/

#include <gsAssembler/gsVisitorPoisson.h> // Stiffness volume integrals
#include <gsAssembler/gsVisitorNeumann.h> // Neumann boundary integrals
#include <gsAssembler/gsVisitorNitsche.h> // Nitsche boundary integrals
#include <gsAssembler/gsVisitorDg.h>      // DG interface integrals
#include "gsRMShellAssembler.h"
#include "gsHyperMeshOut.h"

namespace gismo
{
	// >>> --------------------------------------------------------------------
	template<class T>
	void gsRMShellAssembler<T>::refresh()
	{
		// We use predefined helper which initializes the system matrix
		// rows and columns using the same test and trial space
		Base::scalarProblemGalerkinRefresh();
		// Base::setFixedDofs(m_coefMatrix); 
	}

	template<class T>
	void gsRMShellAssembler<T>::assemble()
	{
		GISMO_ASSERT(m_system.initialized(),
			"Sparse system is not initialized, call initialize() or refresh()");

		// Reserve sparse system
		m_system.reserve(m_bases[0], m_options, this->pde().numRhs());

		/*1. 计算均布载荷的时候，计算表面积要用到雅可比矩阵，为了避免重复计算，
		在计算单元刚度阵的时候先把雅可比矩阵计算出来，放到rhs中，最后再统一乘载荷值*/
		index_t m_cols = 1;
		if (m_BoundaryCondition.m_pPressure.numLoads() != 0)
		{
			// 均布载荷可能没有，也可能不止一个
			m_cols = m_BoundaryCondition.m_pPressure.numLoads();
		}

		// 2. 因为没有找到m_system中的自由度是怎么定义的，所以需要这样初始化m_system的维度
		index_t alldof = m_dofpn * m_basis.size();
		m_system.matrix().resize(alldof, alldof);
		m_system.rhs().resize(alldof, m_cols);
		m_system.setZero();

		// 3. >>> 刚度阵 K 的核心计算过程（多片在 push 中遍历）
		// Assemble over all elements of the domain and applies
		gsRMShellVisitor<T> visitor(m_patches, m_dofmapper, m_material,
			m_BoundaryCondition.m_pPressure, m_SSdata, m_dofpn, m_inte);
		Base::template push<gsRMShellVisitor<T>>(visitor);
		// 如果上一句报错的话就用下面的for语句
		/*for (index_t j = 0; j < m_patches.nPatches(); ++j)
		{
			this->apply(visitor, j);
		}*/

		// 4. 设置边界条件(顺序不可颠倒)
		/*m_BoundaryCondition.setPressure(m_system);
		m_BoundaryCondition.setpLoad(m_system);
		m_BoundaryCondition.setDisp_constra(m_system);*/

		// Assembly is done, compress the matrix
		Base::finalize();

		// 输出等几何刚度阵
		/*
		ofstream igak("../../TestResult/RM_Shell_IGAk.txt");
		for (index_t i=0; i<m_system.matrix().rows(); ++i)
		{
			for (index_t j = 0; j < m_system.matrix().cols(); ++j)
			{
				igak << setw(16) << m_system.matrix().coeffRef(i, j);
			}
			igak << endl;
		}
		igak.close();*/
	}

	// >>> --------------------------------------------------------------------
	// 将位移解转换成矩阵形式 [控制点数 x 自由度]
	template<class T>
	void gsRMShellAssembler<T>::constructSolution(gsMatrix<T>& solVector,
		gsMultiPatch<T>& result)
	{

		// The final solution is the deformed shell, therefore we add the
		// solVector to the undeformed coefficients
		result = m_patches;

		const index_t dim = m_basis.dim() + 1;

		for (index_t p = 0; p < m_patches.nPatches(); ++p)
		{
			// Reconstruct solution coefficients on patch p
			const index_t sz = m_basis[p].size();
			gsMatrix<T>& coeffs = result.patch(p).coefs();

			for (index_t i = 0; i < sz; ++i)
			{
				for (index_t j = 0; j < 3; ++j)
				{
					coeffs(i, j) = solVector(i * m_dofpn + j);
				}
			}
		}
	}

	// >>> --------------------------------------------------------------------
	// 计算应力应变
	template<class T>
	void gsRMShellAssembler<T>::StressStrain(const gsMatrix<T>& solMatrix,
		gsMatrix<T>& StrainMat, gsMatrix<T>& StressMat)
	{
		// 记录节点重复次数，节点应力应变除重复计算次数才是实际的应力应变
		gsMatrix<index_t> nodeRepeats;
		nodeRepeats.setZero(m_nurbsinfo.nurbs_sumcp, 1);

		index_t sumele = 0;        // 全局单元计数
		index_t sumcps = 0;        // 全局控制点计数

		// SSdata 是逐片 push_bach 进vector中的，所以只能分片计算
		for (index_t np = 0; np < m_nurbsinfo.nurbs_sumpatch; ++np)
		{
			// 单元控制点数（考虑不同片上的基函数阶次可能不同）
			index_t eleNodes = (m_nurbsinfo.nurbsdata[np].knot_vector[0].knotvec_degree + 1)
				* (m_nurbsinfo.nurbsdata[np].knot_vector[1].knotvec_degree + 1); // 9
			// 单元的积分点数（所有片肯定都是一样的）
			index_t sumqus = m_SSdata.m_Ng[0].rows();  // 4

			// 单元对应控制点的位移
			gsVector<real_t> EleNodeDisp;
			EleNodeDisp.setZero(eleNodes * m_dofpn); // [9*6 x 1]

			// 记录的是积分点的应变，因为B是在积分点计算的
			gsMatrix<real_t> m_strain;
			gsMatrix<real_t> m_stress;
			m_strain.setZero(sumqus, m_dofpn); // [4 x 6]
			m_stress.setZero(sumqus, m_dofpn);

			// 单元对应的控制点编号的全局编号
			gsMatrix<index_t> globalpNo;
			globalpNo.setZero(eleNodes, 1); // [9 x 1]

			// 逐单元计算应力应变
			for (index_t i = 0; i < m_nurbsinfo.nurbsdata[np].patch_sumele; ++i)
			{
				// 获取单元对应的 当前片的控制点编号 的全局编号
				m_dofmapper.localToGlobal(m_SSdata.m_Pg[sumele + i], np, globalpNo);

				// 单元节点对应的位移向量
				for (index_t j = 0; j < eleNodes; ++j)
				{
					// [9*6 x 1]
					EleNodeDisp.segment(j * m_dofpn, m_dofpn) =
						solMatrix.block(globalpNo(j, 0) * m_dofpn, 0, m_dofpn, 1);
				}

				// 积分点上的应力应变
				for (index_t k = 0; k < sumqus; ++k)
				{
					// [4 x 6] = [4*6 x 9*6] x [9*6 x 1]
					m_strain.row(k) = (m_SSdata.m_Bg[sumele + i].block(k * m_dofpn, 0, m_dofpn, m_dofpn * eleNodes)
						* EleNodeDisp).transpose();
					// [4 x 6] = [4*6 x 6] x [4 x 6]
					m_stress.row(k) = m_SSdata.m_Dg[sumele + i].block(k * m_dofpn, 0, m_dofpn, m_dofpn)
						* m_strain.row(k).transpose();
				}

				// 积分点的基函数值
				gsMatrix<T> NGauss = m_SSdata.m_Ng[sumele + i];
				Eigen::MatrixXf MatrixXf_mat;
				Eigen::MatrixXf MatrixXf_bstress, MatrixXf_xstress;
				Eigen::MatrixXf MatrixXf_bstrain, MatrixXf_xstrain;
				gsMatrix2Matrixxf(NGauss, MatrixXf_mat);
				gsMatrix2Matrixxf(m_stress, MatrixXf_bstress);
				gsMatrix2Matrixxf(m_strain, MatrixXf_bstrain);

				// 求控制点的应力应变
				// [9 x 6]
				MatrixXf_xstress = MatrixXf_mat.colPivHouseholderQr().solve(MatrixXf_bstress);
				MatrixXf_xstrain = MatrixXf_mat.colPivHouseholderQr().solve(MatrixXf_bstrain);

				for (index_t m = 0; m < eleNodes; ++m)
				{
					for (index_t n = 0; n < m_dofpn; ++n)
					{
						StressMat(globalpNo(m, 0), n) += MatrixXf_xstress(m, n);
						StrainMat(globalpNo(m, 0), n) += MatrixXf_xstrain(m, n);
					}
					nodeRepeats(globalpNo(m, 0), 0)++;
				}
			}
			sumele += m_nurbsinfo.nurbsdata[np].patch_sumele;       // 单元总数
			sumcps += m_nurbsinfo.nurbsdata[np].patch_sumcp;        // 控制点总数
		}

		// 处理重复计算的点上的应力应变
		for (index_t i = 0; i < sumcps; ++i)
		{
			StressMat.row(i) /= nodeRepeats(i, 0);
			StrainMat.row(i) /= nodeRepeats(i, 0);
		}
	}

	// >>> --------------------------------------------------------------------
	// 不同矩阵形式的互相转换
	template<class T>
	void gsRMShellAssembler<T>::vector2Matrix(const gsMatrix<T>& solVector,
		gsMatrix<T>& result)
	{
		index_t sz = m_nurbsinfo.nurbs_sumcp;
		for (index_t i = 0; i < sz; ++i)
		{
			for (index_t j = 0; j < m_dofpn; ++j)
			{
				index_t tempr = i * m_dofpn + j;
				result(i, j) = solVector(tempr, 0);
			}
		}
	}

	template<class T>
	void gsRMShellAssembler<T>::VectorgoMatrix(const gsVector<T>& in_vector,
		gsMatrix<T>& result)
	{
		index_t size = in_vector.rows() / m_dofpn;
		for (index_t i = 0; i < size; ++i)
		{
			for (index_t j = 0; j < m_dofpn; ++j)
			{
				index_t tempr = i * m_dofpn + j;
				result(i, j) = in_vector(tempr, 0);
			}
		}
	}

	template <class T>
	void gsRMShellAssembler<T>::gsMatrix2Matrixxf(const gsMatrix<T>& a,
		Eigen::MatrixXf& aCopy)
	{
		int row = a.rows();
		int col = a.cols();

		aCopy.resize(row, col);
		for (int i = 0; i < row; ++i)
		{
			for (int j = 0; j < col; ++j)
			{
				aCopy(i, j) = a(i, j);
			}
		}
	}

	// 将刚度阵从gsSparseMatrix转换为Eigen::SparseMatrix,umf求解器要用
	template <class T>
	void gsRMShellAssembler<T>::gsSparseToSparse(Eigen::SparseMatrix<T>& eiK)
	{
		// 行循环
		for (index_t i = 0; i<m_system.matrix().outerSize(); ++i)
		{
			typename gsSparseMatrix<T, RowMajor,index_t>::iterator it1(m_system.matrix(), i);
			for (; it1; ++it1)
			{
				eiK.coeffRef(i, it1.col()) = it1.value();
			}
		}
		
	}

	template <class T>
	void gsRMShellAssembler<T>::gsMatrixtoVector(vector<T>& rhs_v)
	{
		// 行循环
		for (index_t i = 0; i < m_system.rhs().rows(); ++i)
		{
			rhs_v.push_back(m_system.rhs()(i,0));
		}

	}
	
	// >>> --------------------------------------------------------------------=
	template <class T>
	void patch_test(gsRMShellAssembler<T>& ass)
	{
		//const gsDofMapper& mapper = ass.m_dofmapper.front();
		//dof sequences of all patches
		//vector<index_t> dof = mapper.m_dofs;
		
		//displacements
		gsMatrix<T> dis(ass.m_sum_nodes* ass.m_dofpn, 1);
		//int no_in_patch = dof.size() / ass.m_patches.nPatches();
		for (unsigned int i = 0; i < ass.m_patches.nPatches(); i++)
		{
			gsMatrix<T> coor = (ass.m_patches).patch(i).coefs();
			for (int j = 0; j < coor.rows(); j++)
			{
				//node sequence
				//int nod = dof[i * no_in_patch + j];
				int nod = ass.m_dofmapper.index(j,0,i);
				// disp in x direction u=0.002x
				dis(nod * 6, 0) = coor(j, 0) * 0.002;
				// disp in y direction v=-0.0006y
				dis(nod * 6 + 1, 0) = coor(j, 1) * (-0.0006);

				// shear model
				//dis(nod * 6, 0) = coor(j, 1) * 0.005;
				//dis(nod * 6 + 1, 0) = coor(j, 0) * 0.005;

				dis(nod * 6 + 2, 0) = 0;
				dis(nod * 6 + 3, 0) = 0;
				dis(nod * 6 + 4, 0) = 0;
				dis(nod * 6 + 5, 0) = 0;
			}
		}
		//patch test A stiffness matrix * displacements = F
		gsMatrix<T> res = ass.matrix() * dis;
		
		ofstream out;
		out.open("../../TestResult/RM_Shell_patchtest_F.txt");
		//gsOutputSparseMatrix(ass.m_matrix, "stiff.txt");
		for (int i = 0; i < ass.m_sum_nodes; ++i)
		{
			out << setw(3) << i << setw(20) << res(6 * i, 0);
			out << setw(20) << res(6 * i + 1, 0) << "\n";
		}
		double res_xsum = 0, res_ysum = 0;
		for (int i = 0; i < ass.m_sum_nodes; ++i)
		{
			res_xsum += res(i * 6, 0);
			res_ysum += res(i * 6 + 1, 0);
		}
		out << "Sum: \t" << res_xsum << "\t" << res_ysum << endl;
		out.close();

		// calculate the strains
		gsMatrix<> strainMat, stressMat;
		strainMat.setZero(ass.m_sum_nodes, ass.m_dofpn);
		stressMat.setZero(ass.m_sum_nodes, ass.m_dofpn);
		ass.StressStrain(dis, strainMat, stressMat);

		gsMatrix<> solMat; // [n x 6]
		solMat.setZero(ass.m_sum_nodes, ass.m_dofpn);
		ass.vector2Matrix(dis, solMat);

		gsOutputHMRes("../../TestResult/RM_Shell_patchtest.ascii",
			solMat, stressMat, strainMat);

		ofstream outf("../../TestResult/RM_Shell_patchtest_F.ascii");
		outf << "ALTAIR ASCII FILE\n"
			<< "$TITLE	 = Patch test\n"
			<< "$SUBCASE	 = 1	Subcase	1\n"
			<< "$binding	 = NODE\n"
			<< "$COLUMN_INFO	 = ENTITY_ID\n"
			<< "$RESULT_TYPE	 = Force(v)\n";
		{
			outf << "$TIME = " << setw(16) << 1 << "sec\n";
			for (int j = 0; j < ass.m_sum_nodes; ++j)
			{
				outf << setw(8) << j + 1;
				//输出位移
				{
					outf << setw(16) << res(j * 6, 0);
					outf << setw(16) << res(j * 6 + 1, 0);
					outf << setw(16) << 0 << endl;
				}
			}
		}
		outf.close();
	}

	
	// -------------------------------------------------------------------- <<<
}// namespace gismo
