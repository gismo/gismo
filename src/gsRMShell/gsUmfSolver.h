#pragma once
/**************************************************************************************
*	说明：此程序为仿照郑老师的umfSolve.h写的求解器，旨在调用Umf求解器进行快速运算，
*		郑老师的umfSolver不能满足当前程序的需求，数据结构不一致
*	作者：王洪帅
*	日期：2021-02-03
***************************************************************************************/
#include <iostream>
#include <vector>
#include <map>
#include <umf/umfSolver.h>
#include "Eigen/Dense"
#include "gismo.h"

// #pragma comment(lib, "external/umf/lib/umf_debug.lib")

using namespace Eigen;
using namespace std;
using namespace gismo;

void gsumf_solver(const gsSparseMatrix<double>& A, gsMatrix<double>& mx, const gsMatrix<double>& mb)
{
	vector<int> Ap;
	vector<int> Ai;
	vector<double> Ax;
	vector<double> x;
	vector<double> b;
	x.resize(mx.rows());

	int n = (int)(A.cols());
	Ap.resize(n + 1);
	Ap[0] = 0;

	int num = (int)(A.nonZeros());
	Ai.resize(num);
	Ax.resize(num);

	int k = 0;
	for (int i = 0; i < A.outerSize(); i++)
	{
		Ap[i + 1] = Ap[i];
		for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it)
		{
			Ax[k] = it.value();

			Ai[k] = (int)(it.row());
			k++;
			Ap[i + 1]++;
		}
	}

	for (int i =0 ; i< mb.rows(); ++i)
	{
		b.push_back(mb(i, 0));
	}

	double* null = (double*)NULL;
	void* Symbolic, * Numeric;

	umfpack_di_symbolic(n, n, &Ap[0], &Ai[0], &Ax[0], &Symbolic, null, null);
	umfpack_di_numeric(&Ap[0], &Ai[0], &Ax[0], Symbolic, &Numeric, null, null);
	umfpack_di_free_symbolic(&Symbolic);
	umfpack_di_solve(UMFPACK_A, &Ap[0], &Ai[0], &Ax[0], &x[0], &b[0], Numeric, null, null);
	umfpack_di_free_numeric(&Numeric);

	for (int j =0; j<x.size(); ++j)
	{
		mx(j, 0) = x[j];
	}
}
