/** @file surfaceReparametrization_example.cpp

@brief This example demonstrates how to use the gsSurfaceReparameterization class to improve the parameterization of a given surface.
@details The parameterization of a surface is improved by minimizing a functional that measures the quality of the parameterization.
 The optimization is performed using a gradient-based optimization method.
 The resulting surface has a more uniform distribution of isoparametric lines.

This file is part of the G+Smo library.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

Author(s): Y. Ji
*/

#include <gismo.h>
#include <gsModeling/gsSurfaceReparameterization.h>
#include <gsOptimizer/gsGradientDescent.h>

using namespace gismo;

/**
 * @brief This example demonstrates how to use the gsSurfaceReparameterization class to improve the parameterization of a given surface.
 * 
 * The parameterization of a surface is improved by minimizing a functional that measures the quality of the parameterization.
 * The optimization is performed using a gradient-based optimization method.
 * The resulting surface has a more uniform distribution of isoparametric lines.
 * 
 * @tparam real_t The floating point type, e.g. double.
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 * @return 0 if successful, 1 otherwise.
 */
int main(int argc, char *argv[])
{

  //! [Parse command line]
  index_t solver  = 0;
  std::string INPUT_FILE = "surfaces/crazySurf.xml";

  gsCmdLine cmd("Demonstrates the use of optimizers.");
  cmd.addInt( "s", "solver", "Solver used. 0:gsGradientDescent, 1:gsHLBFGS, 2:IpOpt", solver);
  cmd.addPlainString("input", "Input file", INPUT_FILE);
  try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }
  //! [Parse command line]


  // Load XML file - multi-patch domain
  if (!gsFileManager::fileExists(INPUT_FILE))
  {
    gsWarn << "The file "<<INPUT_FILE<<" cannot be found!\n";
    return EXIT_FAILURE;
  }

  //! [Read geometry]
  // MultiPatch reader
  gsInfo << "Reading file: " << INPUT_FILE << "\n";
  gsMultiPatch<real_t>::uPtr mp = gsReadFile<>(INPUT_FILE);
  gsInfo << "Loaded geometry: " << *mp << "\n";
  //! [Read geometry]

  //! [Write input surface]
  gsWriteParaview(*mp, "input_surface", 1000);
  //! [Write input surface]

#ifdef gsHLBFGS_ENABLED
  //! [Set up optimizer]
  // Set up the optimizer
  gsHLBFGS<real_t> optimizer;
  optimizer.options().setReal("MinGradLen", 1e-6);
  optimizer.options().setReal("MinStepLen", 1e-6);
  optimizer.options().setInt("MaxIterations", 200);
  optimizer.options().setInt("Verbose", 0);
#else
  gsGradientDescent<real_t> optimizer;
  optimizer.options().setReal("MinGradientLength", 1e-6);
  optimizer.options().setReal("MinStepLength", 1e-6);
  optimizer.options().setInt("MaxIterations", 200);
  optimizer.options().setInt("Verbose", 0);
#endif
  //! [Set up optimizer]

  //! [Create reparametrization object]
  // Create the surface reparametrization object
  SurfaceReparameterization<real_t> reparam(*mp,optimizer);
  //! [Create reparametrization object]

  //! [Solve for optimal geometry]
  // Generate the final optimized geometry as a B-Spline surface
  gsMultiPatch<real_t> optSurface = reparam.solve();
  //! [Solve for optimal geometry]

  //! [Write optimized surface]
  // Output the resulting geometry to a Paraview file
  gsWriteParaview(optSurface, "optimized_surface", 1000);
  //! [Write optimized surface]

  return EXIT_SUCCESS;
}
