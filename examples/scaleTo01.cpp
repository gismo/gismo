/** @file scaleTo01.cpp

    @brief Pre-processing step for gsHLBFGS experiments: scale the
    data so that they fit to [0, 1]^3.

	TODO: Make the inverted option available as well.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): D. Mokris
*/

#include <gismo.h>

using namespace gismo;

real_t scaleTo01(real_t tMin, real_t t, real_t tMax)
{
	return (t - tMin) / (tMax - tMin);
}

template <class T>
void scaleTo01(real_t tMin, gsMatrix<T>& mT, real_t tMax)
{
	for(index_t i=0; i<mT.rows(); i++)
		for(index_t j=0; j<mT.cols(); j++)
			mT(i, j) = scaleTo01(tMin, mT(i, j), tMax);
}

real_t scaleFrom01(real_t tMin, real_t tMax, real_t t)
{
	return (tMax - tMin) * t + tMin;
}

template <class T>
void scaleFrom01(real_t tMin, gsMatrix<T>& mT, real_t tMax)
{
	for(index_t i=0; i<mT.rows(); i++)
		for(index_t j=0; j<mT.cols(); j++)
			mT(i, j) = scaleFrom01(tMin, mT(i, j), tMax);
}

int main(int argc, char *argv[])
{
	// Parsing the cmd arguments.
	std::string fin("");
	std::string fout("");

	index_t uvIdIn = 0;
	index_t uvIdOut = 1;
	index_t xyzIdIn = 1;
	index_t xyzIdOut = 0;

	gsCmdLine cmd("Scaling data to/from [0, 1]^3");
	cmd.addString("i", "fin", "filename of the input", fin);
	cmd.addString("o", "fout", "filename of the output", fout);

	cmd.addInt("u", "uvIdIn",   "id of the input matrix with the uv coordinates",   uvIdIn);
	cmd.addInt("v", "uvIdOut",  "id of the output matrix with the uv coordinates",  uvIdOut);
	cmd.addInt("x", "xyzIdIn",  "id of the input matrix with the xyz coordinates",  xyzIdIn);
	cmd.addInt("y", "xyzIdOut", "id of the output matrix with the xyz coordinates", xyzIdOut);

	try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

	// Reading the inputs.
	gsFileData<> fdIn(fin);
	gsMatrix<> xyz, uv;
	fdIn.getId<gsMatrix<>>(uvIdIn,  uv);
	fdIn.getId<gsMatrix<>>(xyzIdIn, xyz);

	// Scaling the data.
	real_t xyzMin = xyz.minCoeff();
	real_t xyzMax = xyz.maxCoeff();
	scaleTo01(xyzMin, xyz, xyzMax);
	gsInfo << "The scale was " << real_t(1.0) / (xyzMax - xyzMin) << std::endl;

	// Writing the outputs
	gsFileData<> fdOut;

	if(uvIdOut < xyzIdOut)
	{
		fdOut << uv;
		fdOut << xyz;
	}
	else
	{
		fdOut << xyz;
		fdOut << uv;
	}
	fdOut.dump(fout);

	return 0;
}
