#pragma once
/** @file gsNURBSinfo.h

	@brief RM shell's NURBS information.

	This file is not part of the G+Smo library.

	This Source Code Form is subject to the terms of the Mozilla Public
	License, v. 2.0. If a copy of the MPL was not distributed with this
	file, You can obtain one at http://mozilla.org/MPL/2.0/.

	Author(s): HS. Wang

	Date:   2020-12-30
*/

#include <gismo.h>
#include <cctype>
#include "gsGeometryEvaluator.hpp"

namespace gismo
{
	/** @brief
		NURBS information of Reissner-Mindlin shell .

		It sets up an assembler and assembles the system patch wise and combines
		the patch-local stiffness matrices into a global system by various methods
		(see gismo::gsInterfaceStrategy). It can also enforce Dirichlet boundary
		conditions in various ways (see gismo::gsDirichletStrategy).

		\ingroup Assembler
	*/
	template < class T>
	class gsNURBSinfo
	{

	public:
		// 默认构造函数
		gsNURBSinfo()
		{}
		//typedef memory::shared_ptr<gsNURBSinfo> Ptr;
		//typedef memory::unique_ptr<gsNURBSinfo> uPtr;

		// 构造函数1
		gsNURBSinfo(gsMultiPatch<T> const& m_patch,
			gsMultiBasis<T> const& m_basis) :multi_patch(m_patch)
		{
			// 方案1 通过m_basis.basis(0)输出txt文件，然后再读入，效果不是很好，效率极低
			// 先输出
			std::ofstream out_info("NURBSinfo.txt");
			out_info << m_basis.basis(0);
			out_info.close();
			// 再读入
			getKnotVectorInfo();
		}

		// 构造函数2
		gsNURBSinfo(gsMultiPatch<T> const& m_patch,
			gsMultiBasis<T> const& m_basis,
			const gsDofMapper& dofmapr)
		{
			// 方案2 使用内部函数获取NURBS相关信息
			/* nurbsdata[patch_info[kontvec_info]]
			* 所有片的信息[ 单片的信息[ 节点向量信息]] */
			// 2.1 初始化
			knotvec_info.knotvec_initialize();			// 单个节点向量信息
			patch_info.igapatch_initialize();			// 单片的NURBS信息
			vector<IGApatch_Info>().swap(nurbsdata);	// 所有片的NURBS信息
			nurbs_sumcp = 0;							// 控制点总数
			nurbs_sumele = 0;							// 单元总数
			nurbs_sumpatch = 0;							// 片总数

			// 2.2 循环所有patch
			nurbs_sumpatch = m_patch.nPatches();		// 片总数
			for (index_t i = 0; i < nurbs_sumpatch; ++i)
			{
				// 2.3 单片的NURBS信息
				patch_info.patch_dim = m_basis[i].dim();
				patch_info.patch_sumcp = m_basis.size(i);
				patch_info.patch_sumele = m_basis[i].numElements();
				// 2.4 总控制点数 总单元数
				nurbs_sumcp += patch_info.patch_sumcp;
				nurbs_sumele += patch_info.patch_sumele;
				// 2.5 单个节点向量信息
				for (index_t j = 0; j < patch_info.patch_dim; ++j)
				{
					knotvec_info.knotvec_ele = m_basis[i].component(j).numElements(); // 节点向量对应单元数 n-p
					knotvec_info.knotvec_size = m_basis[i].component(j).size();		// 节点向量对应控制点数 n
					knotvec_info.knotvec_degree = m_basis[i].component(j).degree(0);	// 节点向量阶次 p
					real_t knot_span = 1.0 / knotvec_info.knotvec_ele;
					knotvec_info.knotvec_span = knot_span;

					// 2.6 生成节点向量（仅适用于均匀节点）
					// [0 0 
					for (index_t m = 0; m < knotvec_info.knotvec_degree; ++m)
					{
						knotvec_info.knotvector.push_back(0.0);
					}
					// 0 to 1
					for (index_t k = 0; k <= knotvec_info.knotvec_ele; ++k)
					{
						real_t knot_temp = real_t(k * knot_span);
						knotvec_info.knotvector.push_back(knot_temp);
					}
					//  1 1]
					for (index_t n = 0; n < knotvec_info.knotvec_degree; ++n)
					{
						knotvec_info.knotvector.push_back(1.0);
					}
					patch_info.knot_vector.push_back(knotvec_info);
				}

				// 2.7 控制点坐标& 编号
				patch_info.phys_coor.setZero(patch_info.patch_sumcp, 3);
				gsMatrix<T>  coeffs = m_patch.patch(i).coefs();
				for (index_t no = 0; no < patch_info.patch_sumcp; ++no)
				{
					index_t tempno = dofmapr.index(no, 0, i);
					patch_info.phys_coor.row(tempno) = coeffs.row(no);
					//patch_info.phys_coor.row(dofmapr.index(no, i)) = coeffs.row(no);
					patch_info.cp_no.push_back(tempno);
				}

				// 2.8 参数坐标
				index_t sum_para = (patch_info.knot_vector[0].knotvec_ele + 1)
					* (patch_info.knot_vector[1].knotvec_ele + 1);
				patch_info.para_coor.resize(sum_para, 2);
				sum_para = 0;
				for (index_t m = patch_info.knot_vector[1].knotvec_degree;
					m <= patch_info.knot_vector[1].knotvec_ele; ++m)
				{
					for (index_t n = patch_info.knot_vector[0].knotvec_degree;
						n <= patch_info.knot_vector[0].knotvec_ele; ++n)
					{
						patch_info.para_coor(sum_para, 0) = patch_info.knot_vector[0].knotvector[n];
						patch_info.para_coor(sum_para, 1) = patch_info.knot_vector[1].knotvector[m];
						++sum_para;
					}
				}

				// 2.9 Greville 横标点
				m_patch.basis(i).anchors_into(patch_info.anch_coor);

				nurbsdata.push_back(patch_info);
			}
		}

		// 析构函数
		~gsNURBSinfo()
		{}

		// 获取节点向量信息的函数模板
		void getKnotVectorInfo();

		// 输出K文件
		// 方案1 利用Geometry信息
		void OutputKfile(const gsGeometry<T>& Geo, const std::string& name);
		// 方案2 利用NURBS信息
		void ExportKfile(const std::string& name);

		// 字符串转换为数字的函数模板
		T trans_str_to_any(std::string s_str);
		//void trans_str_to_int(std::string s_str,T & res);

	public:
		// 方案1的数据结构
		class IGAinfo
		{
		public:
			IGAinfo()
			{
				initialize();
			}
			~IGAinfo() {}

			void initialize()
			{
				vector<real_t>().swap(knot_vector);	// 清空vector
				vector<real_t>().swap(para_coor);
				vector<real_t>().swap(phys_coor);
				knot_direction = 0;				// 参数方向 xi or eta
				knot_degree = 0;				// 节点向量阶次
				knot_size = 0;				// 节点向量数
				knot_minSpan = 0.0;				// 最小节点距
				knot_maxSpan = 0.0;				// 最大节点距
			}
		public:

			vector<real_t> knot_vector;	//	节点向量
			vector<real_t> para_coor;	//	参数坐标
			vector<real_t> phys_coor;	//	物理坐标
			index_t	knot_direction;
			index_t	knot_degree;
			index_t	knot_size;
			real_t	knot_minSpan;
			real_t  knot_maxSpan;
		private:

		};

		// 方案2的数据结构
		// 单个节点向量信息
		class KnotVec_Info
		{
		public:
			KnotVec_Info()
			{
				knotvec_initialize();
			}
			~KnotVec_Info()
			{}

			void knotvec_initialize()
			{
				knotvec_ele = 0;
				knotvec_size = 0;
				knotvec_degree = 0;
				knotvec_span = 0;
				vector<real_t>().swap(knotvector);
			}

		public:
			index_t knotvec_ele;	// 单元数 n-p
			index_t knotvec_size;	// 控制点数 n
			index_t knotvec_degree;	// 基函数阶次 p
			real_t  knotvec_span;	// 基函数节点距，假设是在[0,1]上等距分布
			vector<real_t>	knotvector; // 节点向量 [0 to 1]
		};// end of class KnotVec_Info

		// 单片的NURBS信息
		class IGApatch_Info
		{
		public:
			IGApatch_Info()
			{
				igapatch_initialize();
			}
			~IGApatch_Info() {}

			void igapatch_initialize()
			{
				vector<KnotVec_Info>().swap(knot_vector); // 清空vector
				para_coor.setZero();
				phys_coor.setZero();
				anch_coor.setZero();
				patch_dim = 0;
				patch_sumcp = 0;
				patch_sumele = 0;
			}

		public:
			vector<KnotVec_Info> knot_vector;	//	单片所有维度的节点向量
			gsMatrix<real_t> para_coor;	// 参数坐标
			gsMatrix<real_t> phys_coor;	// 物理坐标
			gsMatrix<real_t> anch_coor; // Greville 横标
			vector<index_t>	 cp_no;		// 控制点编号
			index_t patch_dim;			// 单片的维度，一般是2
			index_t	patch_sumcp;		// 单片的控制点总数 n*m
			index_t patch_sumele;		// 单片的单元总数 （n-p）*（m-q）
		};// end of class IGApatch_Info

	public:
		// 方案1
		std::string		basis_type;	// 基函数的类型
		index_t			nurbs_dim;	// 基函数的维度
		index_t			nurbs_node;	// 结点总数
		real_t			tran_number; // 用于数据类型转换的局部变量
		IGAinfo			info_iga;	// 单个节点向量的信息
		vector<IGAinfo> igainfo;	// 容纳所有的节点向量的容器
		gsMultiPatch<T> multi_patch;

	public:
		// 方案2
		KnotVec_Info		knotvec_info;	// 单个节点向量信息
		IGApatch_Info		patch_info;		// 单片的NURBS信息
		vector<IGApatch_Info> nurbsdata;	// 所有片的NURBS信息
		index_t				nurbs_sumcp;	// 控制点总数
		index_t				nurbs_sumele;	// 单元总数
		index_t				nurbs_sumpatch;	// 片总数
	};
}
