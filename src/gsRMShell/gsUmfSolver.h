#pragma once
/**************************************************************************************
*	˵�����˳���Ϊ����֣��ʦ��umfSolve.hд���������ּ�ڵ���Umf��������п������㣬
*		֣��ʦ��umfSolver�������㵱ǰ������������ݽṹ��һ��
*	���ߣ�����˧
*	���ڣ�2021-02-03
***************************************************************************************/
#include <iostream>
#include <vector>
#include <map>
#include "umf/umfpack.h"
#include "Eigen/Dense"
#include "gismo.h"

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