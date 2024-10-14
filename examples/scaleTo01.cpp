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

template <class T>
void scaleTo01(gsMatrix<T>& xyz, bool verbose)
{
	real_t xyzMin = xyz.minCoeff();
	real_t xyzMax = xyz.maxCoeff();
	scaleTo01(xyzMin, xyz, xyzMax);

	if(verbose)
		gsInfo << std::setprecision(15)
			   << "Scaling TO [0, 1]^3 as follows.\n"
			   << "min:   " << xyzMin << "\n"
			   << "max:   " << xyzMax << "\n"
			   << "scale: " << real_t(1.0) / (xyzMax - xyzMin) << std::endl;
}

real_t scaleFrom01(real_t tMin, real_t t, real_t tMax)
{
	return (tMax - tMin) * t + tMin;
}

template <class T>
void scaleFrom01(real_t tMin, gsMatrix<T>& mT, real_t tMax, bool verbose)
{
	for(index_t i=0; i<mT.rows(); i++)
		for(index_t j=0; j<mT.cols(); j++)
			mT(i, j) = scaleFrom01(tMin, mT(i, j), tMax);

	if(verbose)
		gsInfo << "Scaling FROM [0, 1]^3.\n"
			   << "inverted scale: " << real_t(1.0) / (tMax - tMin) << std::scientific << std::endl;
}

template <class T>
void scaleFrom01(real_t tMin, gsGeometry<T>& geo, real_t tMax, bool verbose)
{
	gsMatrix<T> coefs = geo.coefs();
	scaleFrom01(tMin, coefs, tMax, verbose);
	geo.coefs() = coefs;
}

void scaleGeo(const std::string& fin,
			  const std::string& fout,
			  real_t tMin,
			  real_t tMax,
			  bool verbose)
{
	gsFileData<> fdIn(fin);
	gsGeometry<>::uPtr geo = fdIn.getFirst<gsGeometry<>>();
	scaleFrom01(tMin, *geo.get(), tMax, verbose);
	gsFileData<> fdOut;
	fdOut << *geo.get();
	fdOut.dump(fout);
}

void scalePts(const std::string& fin,
			  const std::string& fout,
			  index_t uvIdIn,
			  index_t uvIdOut,
			  index_t xyzIdIn,
			  index_t xyzIdOut,
			  real_t tMin,
			  real_t tMax,
			  bool verbose)
{
	// Reading the inputs.
	gsFileData<> fdIn(fin);
	gsMatrix<> xyz, uv;
	fdIn.getId<gsMatrix<>>(uvIdIn,  uv);
	fdIn.getId<gsMatrix<>>(xyzIdIn, xyz);

	// Scaling the data.
	if(tMax < tMin) // defaults, i.e., scaling down
		scaleTo01(xyz, verbose);
	else
		scaleFrom01(tMin, xyz, tMax, verbose);

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

	// Intentionally negative interval to see if options are set.
	real_t tMin = 1;
	real_t tMax = 0;

	bool geo = false;
	bool verbose = false;

	gsCmdLine cmd("Scaling data to/from [0, 1]^3.\nWhen setting m and n, the data are scaled FROM [0, 1]^d, otherwise they are scaled TO [0, 1]^d.");
	cmd.addString("i", "fin", "filename of the input", fin);
	cmd.addString("o", "fout", "filename of the output", fout);

	cmd.addInt("u", "uvIdIn",   "id of the input matrix with the uv coordinates",   uvIdIn);
	cmd.addInt("v", "uvIdOut",  "id of the output matrix with the uv coordinates",  uvIdOut);
	cmd.addInt("x", "xyzIdIn",  "id of the input matrix with the xyz coordinates",  xyzIdIn);
	cmd.addInt("y", "xyzIdOut", "id of the output matrix with the xyz coordinates", xyzIdOut);

	cmd.addReal("m", "tMin", "required minimum coordinate", tMin);
	cmd.addReal("n", "tMax", "required maximum coordinate", tMax);

	cmd.addSwitch("g", "geometry", "if set to true, a geometry is scaled, otherwise a point cloud", geo);
	cmd.addSwitch("w", "verbose", "print information", verbose);

	try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

	if(geo)
		scaleGeo(fin, fout, tMin, tMax, verbose);
	else
		scalePts(fin, fout, uvIdIn, uvIdOut, xyzIdIn, xyzIdOut, tMin, tMax, verbose);

	return 0;
}
