/** @file gsNURBSinfo.hpp

	@brief RM shell's NURBS information.

	This file is not part of the G+Smo library.

	This Source Code Form is subject to the terms of the Mozilla Public
	License, v. 2.0. If a copy of the MPL was not distributed with this
	file, You can obtain one at http://mozilla.org/MPL/2.0/.

	Author(s): HS. Wang

	Date:   2020-12-30
*/
#pragma once 
#include "gsNURBSinfo.h"

namespace gismo
{

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
				while(temps != judge_symbol)
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

}


