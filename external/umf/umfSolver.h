#ifndef DLUT_SAE_PERIDYNAMIC_UMF_SOLVER_HXX_20190903
#define DLUT_SAE_PERIDYNAMIC_UMF_SOLVER_HXX_20190903

#include <iostream>
#include <vector>
#include <map>
#include "umf/umfSolver.h"
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

void umf_solver(const SparseMatrix<double>& A, vector<double>& x, const vector<double>& b)
{
	vector<int> Ap;
	vector<int> Ai;
	vector<double> Ax;
	
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

	double *null = (double *)NULL;
	void *Symbolic, *Numeric;

	umfpack_di_symbolic(n, n, &Ap[0], &Ai[0], &Ax[0], &Symbolic, null, null);
	umfpack_di_numeric(&Ap[0], &Ai[0], &Ax[0], Symbolic, &Numeric, null, null);
	umfpack_di_free_symbolic(&Symbolic);
	umfpack_di_solve(UMFPACK_A, &Ap[0], &Ai[0], &Ax[0], &x[0], &b[0], Numeric, null, null);
	umfpack_di_free_numeric(&Numeric);
}

#endif