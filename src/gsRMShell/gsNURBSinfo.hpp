#pragma once
/** @file gsNURBSinfo.hpp

	@brief RM shell's NURBS information.

	This file is not part of the G+Smo library.

	This Source Code Form is subject to the terms of the Mozilla Public
	License, v. 2.0. If a copy of the MPL was not distributed with this
	file, You can obtain one at http://mozilla.org/MPL/2.0/.

	Author(s): HS. Wang

	Date:   2020-12-30
*/

#include "gsNURBSinfo.h"

namespace gismo
{
	/* 方案1 是利用multibasis.basis(0)输出的结果以数据流的形式读入，
	* 该方法的产生完全是受边界条件数据txt文件读入的启发，
	* 该方法主要是针对字符串类型的数据进行处理，
	* 比如要提取“maxSpan=0.25)”中的“0.25”，还是蛮难操作的，
	* 因为该函数对字符串数据的处理具有极高的参考价值，
	* 也是本人花费很大精力设计出来的方法，
	* 尽管有了方案2这种对于整个程序而言非常好的方案
	* 仍希望对方案1进行保留，这里面的思想总有用得上的时候。
	* 王洪帅 2021-01-30 */

	// 获取节点向量信息的函数模板
	template < class T>
	void gsNURBSinfo<T>::getKnotVectorInfo()
	{
		std::ifstream info_in("NURBSinfo.txt");

		std::string line;
		std::stringstream ss_buff;
		std::string temps;

		index_t knot_dir = 0;

		while (!info_in.eof())
		{
			getline(info_in, line);

			if (line.find("Direction") != std::string::npos)
			{
				ss_buff.str("");
				ss_buff.str(line);

				ss_buff >> temps;	// Direction
				ss_buff >> temps;	//	0:
				info_iga.knot_direction = knot_dir;

				// 节点向量
				ss_buff >> temps;	// [
				ss_buff >> temps;	// 0
				std::string judge_symbol("]");
				while (temps != judge_symbol)
				{
					tran_number = trans_str_to_any(temps);
					info_iga.knot_vector.push_back(tran_number);	//	节点向量
					ss_buff >> temps;	// 0 0 1 1 1 ]
				}
				// =============================================
				ss_buff >> temps; // (deg=2,
				std::string dig_str;
				for (auto ca : temps)
				{
					if (isdigit(ca))
					{
						dig_str += ca;
					}
				}
				info_iga.knot_degree = trans_str_to_any(dig_str);
				// ============================================
				ss_buff >> temps; // size=9,
				std::string siz_str;
				for (auto cb : temps)
				{
					if (isdigit(cb))
					{
						siz_str += cb;
					}
				}
				info_iga.knot_size = trans_str_to_any(siz_str);
				// ============================================
				ss_buff >> temps; // minSpan=0.25
				std::string mins_str;
				//char m_dot('.');	
				for (char cc : temps)
				{
					if (isdigit(cc) || cc == '.')
					{
						mins_str += cc;
					}
				}
				info_iga.knot_minSpan = trans_str_to_any(mins_str);
				// ===========================================
				ss_buff >> temps; // maxSpan=0.25)
				std::string maxs_str;
				for (char cd : temps)
				{
					if (isdigit(cd) || cd == '.')
					{
						maxs_str += cd;
					}
				}
				info_iga.knot_maxSpan = trans_str_to_any(maxs_str);

				ss_buff.clear();

				igainfo.push_back(info_iga);
				info_iga.initialize();
				//igainfo[knot_dir]=(info_iga);
				++knot_dir;
			}

			else if (!line.empty())
			{ // TensorBSplineBasis: dim=2, size=36.
				ss_buff.str("");
				ss_buff.str(line);

				ss_buff >> basis_type;	// TensorBSplineBasis:
				// ===========================================
				ss_buff >> temps; // dim=2,
				// 获取数值
				std::string dim_str;
				for (auto ce : temps)
				{
					if (isdigit(ce)) // 如果是数字
					{
						dim_str += ce;
					}
				}
				nurbs_dim = trans_str_to_any(dim_str);
				// ===========================================
				ss_buff >> temps; // size=4.
				std::string size_str;
				for (auto cf : temps)
				{
					if (isdigit(cf))
					{
						size_str += cf;
					}
				}
				nurbs_node = trans_str_to_any(size_str);

				ss_buff.clear();
			}

		}

		gsInfo << basis_type << " ";
		gsInfo << "dim=" << nurbs_dim << ", ";
		gsInfo << "size=" << nurbs_node << ".\n";

		for (index_t j = 0; j < igainfo.size(); ++j)
		{
			gsInfo << "  Direction " << j << ": [ ";
			for (index_t k = 0; k < igainfo[j].knot_vector.size(); ++k)
			{
				gsInfo << igainfo[j].knot_vector[k] << " ";
			}
			gsInfo << "] (";
			gsInfo << "deg=" << igainfo[j].knot_degree
				<< ", size=" << igainfo[j].knot_size
				<< ", minSpan=" << igainfo[j].knot_minSpan
				<< ", maxSpan=" << igainfo[j].knot_maxSpan;
			gsInfo << ")\n";
		}

	}

	// 字符串转换为数字的函数模板
	template<typename T>
	T gsNURBSinfo<T>::trans_str_to_any(std::string s_str)
	{
		T res;
		std::stringstream Repeater;
		Repeater << s_str;
		Repeater >> res;
		return res;
	}

	// >>> --------------------------------------------------------------------
	// 输出K文件
	// 方案1
	template < class T>
	void gsNURBSinfo<T>::OutputKfile(const gsGeometry<T>& Geo, const std::string& name)
	{
		gsMapData<T> map_data;
		map_data.flags = NEED_VALUE | NEED_JACOBIAN | NEED_GRAD_TRANSFORM;
		gsGenericGeometryEvaluator<T, 2, 1> geoEval(Geo, map_data.flags);
		gsMatrix<real_t> globalCp;
		globalCp = geoEval.geometry().coefs();

		ofstream outK(name.c_str());
		outK << "*KEYWORD\n*NODE\n";
		//循环每一个片 输出节点/控制点信息
		for (index_t i = 0; i < nurbs_node; ++i)
		{
			outK << setw(8) << i + 1;
			outK << setw(16) << globalCp(i, 0);
			outK << setw(16) << globalCp(i, 1);
			outK << setw(16) << globalCp(i, 2);
			outK << "\n";
		}

		outK << "*PART\n";
		outK << "*ELEMENT_SHELL\n";
		unsigned int jj = 1;
		for (unsigned int i = 0; i < multi_patch.nPatches(); ++i)
		{
			index_t ele_row = igainfo[0].knot_size - 1;
			index_t ele_col = igainfo[1].knot_size - 1;

			for (unsigned int j = 0; j < ele_col; ++j)
			{
				for (unsigned int k = 0; k < ele_row; ++k)
				{
					index_t st = k + (ele_row + 1) * j;
					outK << setw(8) << jj; // 控制网格单元编号
					outK << setw(8) << 1;
					outK << setw(8) << (st + 1) << setw(8) << (st + 2); // 控制网格节点编号
					outK << setw(8) << (st + 2 + ele_row)
						<< setw(8) << (st + 1 + ele_row) << endl;
					jj++;
				}
			}
		}
		outK << "*END";
		outK.close();
		gsInfo << "Export HyperMesh k files: " << name.c_str() << "\n";

	}
	// 方案2
	template < class T>
	void gsNURBSinfo<T>::ExportKfile(const std::string& name)
	{
		ofstream outK(name.c_str());
		//循环每一个片 输出节点/控制点信息
		outK << "*KEYWORD\n*NODE\n";
		index_t tempi = 1;
		for (index_t i = 0; i < nurbs_sumpatch; ++i)
		{
			for (index_t j = 0; j < nurbsdata[i].patch_sumcp; ++j)
			{
				outK << setw(8) << nurbsdata[i].cp_no[j] + 1;
				outK << setw(16) << nurbsdata[i].phys_coor(j, 0);
				outK << setw(16) << nurbsdata[i].phys_coor(j, 1);
				outK << setw(16) << nurbsdata[i].phys_coor(j, 2);
				outK << "\n";
				++tempi;
			}
		}
		// 输出单元信息
		//outK << "*PART\n";
		outK << "*ELEMENT_SHELL\n";
		unsigned int jj = 1;
		for (unsigned int i = 0; i < nurbs_sumpatch; ++i)
		{
			index_t cp_row = nurbsdata[i].knot_vector[0].knotvec_size;
			index_t cp_col = nurbsdata[i].knot_vector[1].knotvec_size;
			index_t ele_row = cp_row - 1;
			index_t ele_col = cp_col - 1;
			vector<int> seq = nurbsdata[i].cp_no;// 控制点编号
			for (unsigned int j = 0; j < ele_row; ++j)
			{
				for (unsigned int k = 0; k < ele_col; ++k)
				{
					index_t st = k + cp_col * j;
					outK << setw(8) << jj; // 控制网格单元编号
					outK << setw(8) << i + 1;
					// 控制网格节点编号
					outK << setw(8) << (seq[st] + 1)
						<< setw(8) << (seq[st + 1] + 1)
						<< setw(8) << (seq[st + cp_col + 1] + 1)
						<< setw(8) << (seq[st + cp_col] + 1) << endl;
					++jj;
				}
			}
		}
		//outK << "*END";
		outK.close();
		gsInfo << "Export HyperMesh k files: " << name.c_str() << "\n";
	}
	// -------------------------------------------------------------------- <<<
}




