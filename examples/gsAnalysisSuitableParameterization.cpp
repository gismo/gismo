/** @file gsAnalysisSuitableParameterization_example.cpp

    @brief Analysis-suitable parameterization example.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Ye Ji
*/

#include <gismo.h>
#include <gsHLBFGS/gsHLBFGS.h>
#include "gsModeling/gsBarrierPatch.h"
#include "gsModeling/gsBarrierCore.h"
//#include "gsModeling/gsBarrierPatchGenerator.h"

using namespace gismo;

template<typename T>
void makeMapperFixedInterface(const gsMultiPatch<T> &mp, gsDofMapper &mapper)
{
    // control point on interfaces are fixed

    ////! [Make mapper for the design DoFs]
    // Now, we set all the inner control points of each patch as optimization variables
    // It is also possible to set only a part of them as optimization variables
    gsMultiBasis<T> mb(mp);
    mapper.init(mb, mp.targetDim());

    for (gsBoxTopology::const_iiterator it = mb.topology().iBegin();
         it != mb.topology().iEnd(); ++it) // C^0 at the interface
    {
        mb.matchInterface(*it, mapper);
    }

    for (gsMultiPatch<>::const_biterator bit = mp.bBegin();
         bit != mp.bEnd(); ++bit)
    {
        gsMatrix<index_t> idx = mb.basis(bit->patch).allBoundary();
        for (index_t idim = 0; idim != mp.targetDim(); ++idim)
            mapper.markBoundary(bit->patch, idx, idim);
    }

    mapper.finalize();

    gsInfo << "#Numb of free  variables is " << mapper.freeSize() << "\n";
    gsInfo << "#Numb of fixed variables is " << mapper.boundarySize()
           << "\n";
    gsInfo << "#Numb of total variables is " << mapper.size() << "\n\n\n";
}

template<typename T>
void makeMapperFreeInterface(const gsMultiPatch<T> &mp, gsDofMapper &mapper)
{
    // control points on the interfaces are flexible

    ////! [Make mapper for the design DoFs]
    // Now, we set all the inner control points as optimization variables, including
    // the control points on the interfaces
    // It is also possible to set only a part of them as optimization variables
    gsMultiBasis<T> mb(mp);
    mapper.init(mb, mp.targetDim());

    for (gsBoxTopology::const_iiterator it = mb.topology().iBegin();
         it != mb.topology().iEnd(); ++it) // C^0 at the interface
    {
        mb.matchInterface(*it, mapper);
    }

    for (gsMultiPatch<>::const_biterator bit = mp.bBegin();
         bit != mp.bEnd(); ++bit)
    {
        gsMatrix<index_t> idx = mb.basis(bit->patch).boundary(bit->index());
        for (index_t idim = 0; idim != mp.targetDim(); ++idim)
            mapper.markBoundary(bit->patch, idx, idim);
    }

    mapper.finalize();

    gsInfo << "#Numb of free  variables is " << mapper.freeSize() << "\n";
    gsInfo << "#Numb of fixed variables is " << mapper.boundarySize()
           << "\n";
    gsInfo << "#Numb of total variables is " << mapper.size() << "\n\n\n";
}

template<typename T>
void outputResult(const gsMultiPatch<T> &mp, const std::string &filename)
{
    // write multi-patch geometry to file
    gsWrite(mp, filename);

    // compute parameterization quality metrics and write .pvd file
    gsExprEvaluator<T> ev;
    gsMultiBasis<T> mb(mp);
    ev.setIntegrationElements(mb);
    gsExprAssembler<>::geometryMap G = ev.getMap(mp);
//    GISMO_ENSURE(mp.nPatches() == 1, "Does not yet work for multi-patch, "
//                                       "but multi-patch has "
//            << mp.nPatches() << " patches");

    // resulting mesh
    gsOptionList options;
    options.addReal("quA",
                    "Number of quadrature points: quA*deg + quB; For patchRule: Order of the target space",
                    1.0);
    options.addInt("quB",
                   "Number of quadrature points: quA*deg + quB; For patchRule: Regularity of the target space",
                   1);
    options.addInt("quRule",
                   "Quadrature rule [1:GaussLegendre,2:GaussLobatto,3:PatchRule]",
                   1);
    options.addInt("overInt", "Apply over-integration or not?", 0);
    options.addInt("plot.npts", "Number of sampling points for plotting",
                   10000);
    options.addSwitch("plot.mesh",
                      "If true, plot the element mesh of the parameterization",
                      false);
    options.addSwitch("plot.net",
                      "If true, plot the control net of the parameterization",
                      false);

    ev.options() = options;
    auto mdim = mp.targetDim();
    if (mdim == 2) ev.options().setInt("plot.npts", 1000);
    else if (mdim == 3) ev.options().setInt("plot.npts", 10000);

    // TODO: colormap in paraview has a name of SolutionField, need to change it to scaledJacobian and uniformity
    // TODO: seem the values of the metrics are only on quadrature points
    // Scaled Jacobian metric
    T minScaledJacobian = 0.0, maxScaledJacobian = 0.0, integralScaledJacobian = 0.0;
    if (mdim == 2)
    {
        auto metric_ScaledJacobian =
                jac(G).det() / (jac(G)[0].norm() * jac(G)[1].norm());
        ev.writeParaview(metric_ScaledJacobian, G,
                         filename + "_ScaledJacobian");
        minScaledJacobian = ev.template min(metric_ScaledJacobian);
        maxScaledJacobian = ev.template max(metric_ScaledJacobian);
        integralScaledJacobian = ev.template integral(
                metric_ScaledJacobian);
    } else if (mdim == 3)
    {
        auto metric_ScaledJacobian = jac(G).det() /
                                     (jac(G)[0].norm() * jac(G)[1].norm() *
                                      jac(G)[2].norm());
        ev.writeParaview(metric_ScaledJacobian, G,
                         filename + "_ScaledJacobian");
        minScaledJacobian = ev.template min(metric_ScaledJacobian);
        maxScaledJacobian = ev.template max(metric_ScaledJacobian);
        integralScaledJacobian = ev.template integral(
                metric_ScaledJacobian);
    }
    gsInfo << "Parameterization quality info (reference, not strictly):\n";
    gsInfo << " Scaled Jacobian:    min.: " << minScaledJacobian
           << "   max: " << maxScaledJacobian
           << "     integral: " << integralScaledJacobian / mp.nPatches()
           << "\n";

    // Uniformity metric
    T areaTotal = 0.0;
    if (mdim == 2)
        areaTotal = gsBarrierCore<2, T>::computeArea(mp);
    else if (mdim == 3)
        areaTotal = gsBarrierCore<3, T>::computeArea(mp);

    gsConstantFunction<T> areaConstFunc(areaTotal, mdim);
    auto area = ev.getVariable(areaConstFunc);
    auto metric_Uniformity = pow((jac(G).det() - area.val()) / area.val(), 2);
    ev.writeParaview(metric_Uniformity, G,
                     filename + "_Uniformity");

    T minUniformity = ev.template min(metric_Uniformity);
    T maxUniformity = ev.template max(metric_Uniformity);
    T integralUniformity = ev.template integral(metric_Uniformity);
    gsInfo << " Uniformity:    min.: " << minUniformity
           << "    max: " << maxUniformity
           << "    integral: " << integralUniformity / mp.nPatches() << "\n";

    // Frobenius condition number metric
//    auto metric_FrobCondNum = jac(G).norm() * jac(G).inv().norm();
//    ev.writeParaview(metric_FrobCondNum, G, filename + "_FrobCondNum");
//
//    T minFrobCondNum = ev.template min(metric_FrobCondNum);
//    T maxFrobCondNum = ev.template max(metric_FrobCondNum);
//    T integralFrobCondNum = ev.template integral(metric_FrobCondNum);
//    gsInfo << " Frobenius condition number:    min.: " << minFrobCondNum << "   max: " << maxFrobCondNum
//           << "     integral: " << integralFrobCondNum << "\n";
}

int main(int argc, char *argv[])
{
#ifdef _OPENMP
    gsInfo << "Available threads: " << omp_get_max_threads() << "\n";
#endif

    //////////////////// STEP 1: read a multi-patch file /////////////////
    index_t initialMethod = 0;
    index_t paramMethod = 1;
    index_t verbose = 0;
    index_t numRefine = 0;
    index_t numElevate = 0;
    index_t AAPreconditionType = 0;
    bool free = false, isBRep = false;

    // Load XML file containing a multi-patch
    std::string filename_input("planar/multipatch_tunnel.xml");

    // Read input from command line arguments
    //! [Parse command line]
    gsCmdLine cmd(
            "Hi, give me a file (eg: .xml) containing multi-patch"
            "and I will try to parameterize it!");

    cmd.addSwitch("bRep", "Input data is B-Rep", isBRep);
    cmd.addPlainString("input",
                       "Name of the input file containing multi-patch",
                       filename_input);
    cmd.addInt("i", "initial",
               "Initialization Method: 0 Coons' patch (default), 1 Spring patch, 2: Cross-Ap. patch",
               initialMethod);
    cmd.addInt("e", "degreeElevation",
               "Number of degree elevation steps to perform before solving (0: equalize degree in all directions)",
               numElevate);
    cmd.addInt("r", "uniformRefine", "Number of Uniform h-refinement",
               numRefine);
    cmd.addInt("p", "paramMethod",
               "Parameterization Method: 0 Barrier patch (default), 1 Penalty patch, 2: PDE patch",
               paramMethod);
    cmd.addInt("v", "Verbose",
               "Print output 0: no print, 1: summary, 2: full print",
               verbose);
    cmd.addInt("a", "preconditionType",
               "Preconsitioner type for Anderson Acceleration", AAPreconditionType);
    cmd.addSwitch("free", "Make interfaces free", free);
    try { cmd.getValues(argc, argv); }
    catch (int rv) { return rv; }
    //! [Parse command line]

    /////////////////////////////////////////////////////////
//    std::string fn_inputBase("zhanzheng_xml/");
//    std::string fn_output("times");
//
//    fn_output += std::to_string(paramMethod) + ".txt";
//
//    std::ofstream outputFile(fn_output);
//    if (outputFile.is_open()) {
//      outputFile << "ifile" << " "
//                 << "time" << std::endl;
//      //    outputFile.close();
//    } else {
//      std::cerr << "Unable to open" << fn_output << "file." << std::endl;
//      return 1;
//    }
//    for (auto ifile=0; ifile!=977; ++ifile) {
//
//      filename_input = fn_inputBase + std::to_string(ifile) + ".xml";
      ////////////////////////////////////////////////////////////

      // Load XML file - multi-patch computational domain
      //! [Read geometry]
      if (!gsFileManager::fileExists(filename_input)) {
        gsWarn << "The file cannot be found!\n";
        return EXIT_FAILURE;
      }

      gsInfo << "Read file \"" << filename_input << "\"\n";
      gsMultiPatch<real_t>::uPtr mp = gsReadFile<>(filename_input);
      gsInfo << " Got" << *mp << " \n";
      //! [Read geometry]

      //! [construct a initial parameterization if the input data is B-Rep]
      short_t pardim = mp->domainDim();
      short_t geodim = mp->targetDim();
      if (isBRep || pardim != geodim) {
        gsSpringPatch<real_t> springInitializer(*mp);
        gsInfo << "Created a " << springInitializer.compute() << "\n";
        mp->clear();
        *mp = springInitializer.result();

//        if (geodim == 2) {
//          gsBarrierPatch<2, real_t> patchGenerator(*mp);
//          patchGenerator.options().setInt("InitialMethod", initialMethod);
//          patchGenerator.initialization();
//          *mp = patchGenerator.result();
//        } else if (geodim == 3) {
//          gsBarrierPatch<3, real_t> patchGenerator(*mp);
//          patchGenerator.options().setInt("InitialMethod", initialMethod);
//          patchGenerator.initialization();
//          *mp = patchGenerator.result();
//        }
      }
      //! [construct a initial parameterization if the input data is B-Rep]

      //! [align the orientations of a multi-patch parameterization]
      // TODO: TAKE CARE! This line may change the interface information!
      // For multi-patch parameterization, some patches can have inverted orientations.
      // However, consistent orientations are mandatory in our method.
      // The next line is to align the orientations of a multi-patch parameterization.
      // Note, the following line will mark all the interfaces and boundaries automatically.
      mp->fixOrientation();
      mp->computeTopology();
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

      // TODO: another usage of gsBarrierPatch is to give it a use-defined mapper
//    gsDofMapper mapper;
//    makeMapperFixedInterface(*mp, mapper);
//
//    gsDofMapper mapper2;
//    makeMapperFreeInterface(*mp, mapper2);

      gsStopwatch timer;
      //! [perform analysis-suitable parameterization construstion]
      gsMultiPatch<> result;
      if (geodim == 2) {
        gsBarrierPatch<2, real_t> opt(*mp, !free);
        opt.options().setInt("Verbose", verbose);
        opt.options().setInt("ParamMethod", paramMethod);
        opt.options().setInt("AAPreconditionType", AAPreconditionType);
        opt.compute();
        result = opt.result();
      } else if (geodim == 3) {
        gsBarrierPatch<3, real_t> opt(*mp, !free);
        opt.options().setInt("Verbose", verbose);
        opt.options().setInt("ParamMethod", paramMethod);
        opt.compute();
        result = opt.result();
      } else {
        gsInfo << "current version only support pardim = geodim = 2 or 3.\n";
        return EXIT_FAILURE;
      }
      //! [perform analysis-suitable parameterization construstion]

      //! [output the parameterization result]
      std::string filename_output = filename_input;
      filename_output = filename_output.substr(
          filename_output.find_last_of("/\\") +
              1); // file name without a path
      filename_output = filename_output.substr(0, filename_output.find_last_of(
          ".\\")); // file name without an extension
      filename_output += "_result";
      filename_output +=
          "R" + std::to_string(numRefine) + "E" + std::to_string(numElevate);
      if (!free)
        filename_output += "Fixed";
      else
        filename_output += "Free";

      outputResult(result, filename_output);
      gsInfo << "The results have been written into " << filename_output
             << "\n";
    //! [output the parameterization result]

    return EXIT_SUCCESS;
}
