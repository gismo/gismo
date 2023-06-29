/** @file as_parameterization_example.cpp

    @brief Analysis-suitable parameterization example.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji
*/

#include <gismo.h>
#include "gsModeling/gsBarrierPatch.h"

using namespace gismo;

template<typename T>
void makeMapper(const gsMultiPatch<T> &mp,
                gsDofMapper &mapper,
                std::function<gsMatrix<index_t>(const gsBasis<T> &,
                                                const typename gsMultiPatch<>::const_biterator &)> getBoundaryIndices) {

  gsMultiBasis<T> mb(mp);
  mapper.init(mb, mp.targetDim());

  for (gsBoxTopology::const_iiterator it = mb.topology().iBegin();
       it != mb.topology().iEnd(); ++it) {
    mb.matchInterface(*it, mapper);
  }

  for (gsMultiPatch<>::const_biterator bit = mp.bBegin(); bit != mp.bEnd();
       ++bit) {
    gsMatrix<index_t> idx = getBoundaryIndices(mb.basis(bit->patch), bit);
    for (index_t idim = 0; idim != mp.targetDim(); ++idim)
      mapper.markBoundary(bit->patch, idx, idim);
  }

  mapper.finalize();

  gsInfo << "#Numb of free  variables is " << mapper.freeSize() << "\n";
  gsInfo << "#Numb of fixed variables is " << mapper.boundarySize() << "\n";
  gsInfo << "#Numb of total variables is " << mapper.size() << "\n\n\n";
}

template<typename T>
void makeMapperFixedInterface(const gsMultiPatch<T> &mp, gsDofMapper &mapper) {
  makeMapper(mp,
             mapper,
             [](const gsBasis<T> &basis,
                const typename gsMultiPatch<>::const_biterator &bit) { return basis.allBoundary(); });
}

template<typename T>
void makeMapperFreeInterface(const gsMultiPatch<T> &mp, gsDofMapper &mapper) {
  makeMapper(mp,
             mapper,
             [](const gsBasis<T> &basis,
                const typename gsMultiPatch<>::const_biterator &bit) {
               return basis.boundary(bit->index());
             });
}

template<typename T>
void outputResult(const gsMultiPatch<T> &mp, const std::string &filename) {
  gsWrite(mp, filename);

  gsExprEvaluator<T> ev;
  gsMultiBasis<T> mb(mp);
  ev.setIntegrationElements(mb);
  gsExprAssembler<>::geometryMap G = ev.getMap(mp);

  auto mdim = mp.parDim();

  gsOptionList options;
  options.addReal("quA", "Number of quadrature points: quA*deg + quB", 1.0);
  options.addInt("quB", "Number of quadrature points: quA*deg + quB", 1);
  options.addInt("quRule",
                 "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                 1);
  options.addInt("overInt", "Apply over-integration or not?", 0);
  options.addSwitch("plot.mesh",
                    "If true, plot the element mesh of the parameterization",
                    false);
  options.addSwitch("plot.net",
                    "If true, plot the control net of the parameterization",
                    false);
  if (mdim == 2)
    options.addInt("plot.npts",
                   "Number of sampling points for plotting",
                   1000);
  else
    options.addInt("plot.npts",
                   "Number of sampling points for plotting",
                   10000);
  ev.options() = options;

  gsInfo << "\n Parameterization quality info (only values for quadrature "
            "points). \n More information please refer to paraview files! \n";
//   Scaled Jacobian metric
  T minScaledJacobian = 0.0, maxScaledJacobian = 0.0,
      integralScaledJacobian = 0.0;
  if (mdim == 2) {
    auto metric_ScaledJacobian =
        jac(G).det() / (jac(G)[0].norm() * jac(G)[1].norm());
    ev.writeParaview(metric_ScaledJacobian, G,
                     filename + "_ScaledJacobian");
    minScaledJacobian = ev.template min(metric_ScaledJacobian);
    maxScaledJacobian = ev.template max(metric_ScaledJacobian);
    integralScaledJacobian = ev.template integral(metric_ScaledJacobian);
  } else if (mdim == 3) {
    auto metric_ScaledJacobian = jac(G).det() /
        (jac(G)[0].norm() * jac(G)[1].norm() *
            jac(G)[2].norm());
    ev.writeParaview(metric_ScaledJacobian, G, filename + "_ScaledJacobian");
    minScaledJacobian = ev.template min(metric_ScaledJacobian);
    maxScaledJacobian = ev.template max(metric_ScaledJacobian);
    integralScaledJacobian = ev.template integral(metric_ScaledJacobian);
  }

  gsInfo << " Scaled Jacobian:    min.: " << minScaledJacobian
         << "   max: " << maxScaledJacobian
         << "     integral: " << integralScaledJacobian / mp.nPatches()
         << "\n";

  // Uniformity metric
  T areaTotal = ev.integral(jac(G).det());
  gsConstantFunction<T> areaConstFunc(areaTotal, mdim);
  auto area = ev.getVariable(areaConstFunc);
  auto metric_Uniformity = pow((jac(G).det() - area.val()) / area.val(), 2);
  ev.writeParaview(metric_Uniformity, G, filename + "_Uniformity");

  T minUniformity = ev.template min(metric_Uniformity);
  T maxUniformity = ev.template max(metric_Uniformity);
  T integralUniformity = ev.template integral(metric_Uniformity);
  gsInfo << " Uniformity:    min.: " << minUniformity
         << "    max: " << maxUniformity
         << "    integral: " << integralUniformity / mp.nPatches() << "\n";
}

int main(int argc, char *argv[]) {
#ifdef _OPENMP
  gsInfo << "Available threads: " << omp_get_max_threads() << "\n";
#endif

  //////////////////// STEP 1: read a multi-patch file /////////////////
  // Variables initialization with more descriptive names
  index_t paramMethod = 0; // Parameterization method
  index_t verboseMode = 0; // Verbose mode: 0 - no print, 1 - summary, 2 - full print
  index_t numRefine = 0; // Number of Uniform h-refinement
  index_t numElevate = 0; // Number of degree elevation steps to perform
  index_t AAPreconditionType = 0; // Preconditioner type for Anderson Acceleration
  bool isInterfaceFree = false; // Make interfaces free
  bool isBRep = false; // Input data is B-Rep

  // Constants for input file and description
  std::string INPUT_FILE = "breps/other/TUDflame.xml";
  const std::string DESCRIPTION = "Hi, give me a file (eg: .xml) containing multi-patch and I will try to parameterize it!";

  // Read input from command line arguments
  gsCmdLine cmd(DESCRIPTION);
  cmd.addSwitch("bRep", "Input data is B-Rep", isBRep);
  cmd.addPlainString("input", "Name of the input file containing multi-patch", INPUT_FILE);
  cmd.addInt("e", "degreeElevation", "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)", numElevate);
  cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement", numRefine);
  cmd.addInt("p", "paramMethod", "Parameterization Method: 0 Barrier patch (default), 1 Penalty patch, 2 Penalty patch (2), 3: PDE patch", paramMethod);
  cmd.addInt("v", "Verbose", "Print output 0: no print, 1: summary, 2: full print", verboseMode);
  cmd.addInt("a", "preconditionType", "Preconsitioner type for Anderson Acceleration", AAPreconditionType);
  cmd.addSwitch("free", "Make interfaces free", isInterfaceFree);
  try { cmd.getValues(argc, argv); }
  catch (int rv) { return rv; }
  //! [Parse command line]

  // Load XML file - multi-patch computational domain
  //! [Read geometry]
  // Check if the input file exists
  if (!gsFileManager::fileExists(INPUT_FILE)) {
    gsWarn << "The file cannot be found!\n";
    return EXIT_FAILURE;
  }

  // MultiPatch reader
  gsInfo << "Read file \"" << INPUT_FILE << "\"\n";
  gsMultiPatch<real_t>::uPtr mp = gsReadFile<>(INPUT_FILE);
  gsInfo << " Got" << *mp << " \n";
  //! [Read geometry]

  //! [construct a initial parameterization if the input data is B-Rep]
  short_t pardim = mp->parDim();
  short_t geodim = mp->targetDim();
  // If input data is B-Rep or dimension does not match, construct an initial parameterization
  if (isBRep || pardim != geodim) {
    gsSpringPatch<real_t> springInitializer(*mp);
    gsInfo << "Created a " << springInitializer.compute() << "\n";
    mp->clear();
    *mp = springInitializer.result();
  }
  //! [construct a initial parameterization if the input data is B-Rep]

  //! [align the orientations of a multi-patch parameterization]
  mp->fixOrientation();
  //! [align the orientations of a multi-patch parameterization]
  
  //! [refine the geometry for better result]
  // Elevate and p-refine the multi-patch geometry to order p + numElevate
  if (numElevate != 0) {
    mp->degreeElevate(numElevate);
    gsInfo << numElevate << " degree elevation(s) performed!\n";
  }

  // h-refine each patch
  if (numRefine != 0) {
    for (auto r = 0; r < numRefine; ++r)
      mp->uniformRefine();
    gsInfo << numRefine << " uniform h-refinement(s) performed!\n";
  }

  if (numRefine != 0 || numElevate != 0)
    gsInfo << "Now, the multi-patch geometry is: " << *mp << "\n";
  //! [refine the geometry for better result]

  //! [perform analysis-suitable parameterization construstion]
  gsMultiPatch<> result;
  if (geodim == 2) {
    gsBarrierPatch<2, real_t> opt(*mp, !isInterfaceFree);
    opt.options().setInt("Verbose", verboseMode);
    opt.options().setInt("ParamMethod", paramMethod);
    opt.options().setInt("AAPreconditionType", AAPreconditionType);
    opt.compute();
    result = opt.result();
  } else if (geodim == 3) {
    gsBarrierPatch<3, real_t> opt(*mp, !isInterfaceFree);
    opt.options().setInt("Verbose", verboseMode);
    opt.options().setInt("ParamMethod", paramMethod);
    opt.compute();
    result = opt.result();
  } else {
    gsInfo << "current version only support pardim = geodim = 2 or 3.\n";
    return EXIT_FAILURE;
  }
  //! [perform analysis-suitable parameterization construstion]

  //! [output the parameterization result]
  std::string outputFilename = INPUT_FILE;
  outputFilename = outputFilename.substr(outputFilename.find_last_of("/\\") +1); // file name without a path
  outputFilename = outputFilename.substr(0, outputFilename.find_last_of(".\\")); // file name without an extension
  outputFilename += "_result";
  outputFilename += "R"+std::to_string(numRefine)+"E"+std::to_string(numElevate);
  if (!isInterfaceFree) outputFilename += "Fixed";
  else outputFilename += "Free";

  outputResult(result, outputFilename);
  gsInfo << "The results have been written into " << outputFilename << "\n";
  //! [output the parameterization result]

  return EXIT_SUCCESS;
}
