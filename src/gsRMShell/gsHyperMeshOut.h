#pragma once
/** @file gsHyperMeshOut.h

    @brief RM shell's deformation result.

    This file is not part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): HS. Wang

    Date:   2021-01-05
*/

#include <iostream>
#include <stdio.h>
#include <gismo.h>
#include <vector>

namespace gismo
{
	/*计算结果的ascii文件*/
	template<class T>
	void gsOutputHMRes(const std::string& name,
		gsMatrix<T>& disp,
		gsMatrix<T>& stress,
		gsMatrix<T>& strain)
	{
		ofstream out(name.c_str());
		out << "ALTAIR ASCII FILE\n"
			<< "$TITLE	 = Transient analysis\n"
			<< "$SUBCASE	 = 1	Subcase	1\n"
			<< "$binding	 = NODE\n"
			<< "$COLUMN_INFO	 = ENTITY_ID\n"
			<< "$RESULT_TYPE	 = Displacement(v), Stress(t), Strain(t)\n";
		// 时间
		out << "$TIME = " << setw(16) << 1.0 << "sec\n";
		for (int j = 0; j < disp.rows(); ++j) // 遍历所有控制点
		{
			out << setw(8) << j + 1;
			// 1. 输出位移解
			for (int k = 0; k < disp.cols(); ++k) // 遍历所有自由度 6个
			{
				out << setw(16) << disp(j, k);
			}
			// 2. 输出应力应变
			if (stress.rows() == 0)
			{
				for (int k2 = 0; k2 < 12; ++k2)
				{
					out << setw(16) << 0;
				}
			}
			else
			{
				for (int k2 = 0; k2 < stress.cols(); k2++)
				{
					out << setw(16) << stress(j, k2);
				}
				for (int k2 = 0; k2 < strain.cols(); k2++)
				{
					out << setw(16) << strain(j, k2);
				}
			}
			out << "\n";
		}// end of for-loop
		out.close();
	}

}

