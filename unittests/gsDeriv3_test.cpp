/** @file gsDeriv3_test.cpp

    @brief Tests evaluation of third derivatives.

	@Author(s): D. Mokris

**/

#include "gismo_unittest.h"

using namespace gismo;

SUITE(deriv3_into)
{
	// Idea: Univariate cubic Bernstein basis has constant third derivatives.
	TEST(deriv3_KnownValues_Bernstein)
	{
		typedef real_t T;

		// Cubic Bernstein basis
		gsKnotVector<T> kv(0.0, 1.0, 0, 4); // clamped cubic, no internal knots
		gsBSplineBasis<T> basis(kv);

		// Sanity: should have 4 basis functions
		CHECK(basis.size() == 4);

		// Evaluation point (any point in (0,1) works)
		gsMatrix<T> u(1,1);
		u << 0.37; // arbitrary

		gsMatrix<T> result;
		basis.deriv3_into(u, result);

		// result(i,0) = third derivative of basis function i

		CHECK_CLOSE(result(0,0), -6.0,  1e-12);
		CHECK_CLOSE(result(1,0),  18.0, 1e-12);
		CHECK_CLOSE(result(2,0), -18.0, 1e-12);
		CHECK_CLOSE(result(3,0),  6.0,  1e-12);
	}

	// Idea: Third derivatives of a quadratic basis need to sum to zero in any point.
	TEST(deriv3_PartitionOfUnity)
	{
		typedef real_t T;

		gsKnotVector<T> kv(0.0, 3.0, 1, 3);
		kv.uniformRefine(2);

		gsTensorBSplineBasis<2,T> basis(kv, kv);

		gsMatrix<T> result;

		for (T uu = 0.2; uu < 2.8; uu += 0.4)
		{
			for (T vv = 0.2; vv < 2.8; vv += 0.4)
			{
				gsMatrix<T> u(2,1);
				u << uu, vv;

				basis.deriv3_into(u, result);

				// result: (#basis functions) × (#derivative components)
				// or transposed depending on G+Smo version

				// Sum over all active basis functions
				T sum_duuu = 0.0, sum_duuv = 0.0, sum_duvv = 0.0, sum_dvvv = 0.0;

				// There are four derivatives for each active function.
				for (index_t i = 0; i < result.rows(); i += 4)
				{
					sum_duuu += result(i,     0);
					sum_dvvv += result(i + 1, 0);
					sum_duuv += result(i + 2, 0);
					sum_duvv += result(i + 3, 0);
				}

				CHECK_CLOSE(sum_duuu, 0.0, 1e-10);
				CHECK_CLOSE(sum_dvvv, 0.0, 1e-10);
				CHECK_CLOSE(sum_duuv, 0.0, 1e-10);
				CHECK_CLOSE(sum_duvv, 0.0, 1e-10);
			}
		}
	}

	// Idea: Derivatives of the tensor basis must be products of univariate derivatives.
	TEST(deriv3_TensorStructure)
	{
		typedef real_t T;

		gsKnotVector<T> kv(0.0, 3.0, 1, 3);
		kv.uniformRefine(1);

		gsBSplineBasis<T> uBasis(kv);
		gsBSplineBasis<T> vBasis(kv);

		gsTensorBSplineBasis<2,T> tBasis(kv, kv);

		gsMatrix<T> u(1,1), v(1,1);
		u(0,0) = 1.3;
		v(0,0) = 1.7;

		gsMatrix<T> Nu, dNu, d2Nu, d3Nu;
		gsMatrix<T> Nv, dNv, d2Nv, d3Nv;

		uBasis.eval_into(u, Nu);
		uBasis.deriv_into(u, dNu);
		uBasis.deriv2_into(u, d2Nu);
		uBasis.deriv3_into(u, d3Nu);

		vBasis.eval_into(v, Nv);
		vBasis.deriv_into(v, dNv);
		vBasis.deriv2_into(v, d2Nv);
		vBasis.deriv3_into(v, d3Nv);

		gsMatrix<T> uv(2,1);
		uv << u(0,0), v(0,0);

		gsMatrix<T> result;
		tBasis.deriv3_into(uv, result);

		index_t nd = 4; // number of derivatives
		index_t tid = 0; // tensor id
		for(index_t j=0; j < Nv.rows(); j++)
			for(index_t i=0; i < Nu.rows(); i++)
			{
				CHECK_CLOSE(result(tid * nd,     0), d3Nu(i, 0) *   Nv(j, 0), 1e-10);
				CHECK_CLOSE(result(tid * nd + 1, 0),   Nu(i, 0) * d3Nv(j, 0), 1e-10);
				CHECK_CLOSE(result(tid * nd + 2, 0), d2Nu(i, 0) *  dNv(j, 0), 1e-10);
				CHECK_CLOSE(result(tid * nd + 3, 0),  dNu(i, 0) * d2Nv(j, 0), 1e-10);

				// Move to the next function.
				tid++;
			}
	}
}
