/** @file linearSolveAutoDiff.cpp

    @brief Example demonstrating automatic differentiation through linear solvers

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst

*/

//! [Include namespace]
#include <gismo.h>
#include <gsAutoDiff/gsAutoDiff2.h>
#include <gsAutoDiff/gsAutoDiffUtils.h>
// Include reverse-mode AD with full GISMO matrix support
#include <gsAutoDiff/gsAutoDiffVar.h>

#include <gsElasticity/src/gsElasticityAssembler.h>
#include <gsElasticity/src/gsIterative.h>
#include <gsElasticity/src/gsWriteParaviewMultiPhysics.h>
#include <gsElasticity/src/gsMaterialBase.h>
#include <gsElasticity/src/gsLinearMaterial.h>
#include <gsElasticity/src/gsNeoHookeLogMaterial.h>


using namespace gismo;
using namespace autodiff;

using namespace gismo;
//! [Include namespace]

using namespace autodiff;

int main()
{   
    gsInfo << "This is Cook's membrane benchmark with nonlinear elasticity solver.\n";

    //=====================================//
                // Input //
    //=====================================//

    std::string filename = "cooks.xml";
    dual youngsModulus = 240.565e6;
    dual poissonsRatio = 0.4;
    index_t numUniRef = 4;
    index_t numDegElev = 1;
    index_t numPlotPoints = 10000;

    // // minimalistic user interface for terminal
    // gsCmdLine cmd("This is Cook's membrane benchmark with nonlinear elasticity solver.");
    // cmd.addReal("p","poisson","Poisson's ratio used in the material law",poissonsRatio);
    // cmd.addInt("r","refine","Number of uniform refinement application",numUniRef);
    // cmd.addInt("d","degelev","Number of degree elevation application",numDegElev);
    // cmd.addInt("s","point","Number of points to plot to Paraview",numPlotPoints);
    // try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    //=============================================//
        // Scanning geometry and creating bases //
    //=============================================//

    autodiff::detail::seed<1>(poissonsRatio, 1.0);  // Seed poissonsRatio for differentiation

    // scanning geometry
    gsMultiPatch<dual> geometry;
    gsReadFile<dual>(filename, geometry);
    // creating bases
    gsMultiBasis<dual> basisDisplacement(geometry);
    for (index_t i = 0; i < numDegElev; ++i)
        basisDisplacement.degreeElevate();
    for (index_t i = 0; i < numUniRef; ++i)
        basisDisplacement.uniformRefine();

    //=============================================//
        // Setting loads and boundary conditions //
    //=============================================//

    // neumann BC
    gsConstantFunction<dual> f(0.,625e4,2);

    // boundary conditions
    gsBoundaryConditions<dual> bcInfo;
    for (index_t d = 0; d < 2; ++d)
        bcInfo.addCondition(0,boundary::west,condition_type::dirichlet,nullptr,d);
    bcInfo.addCondition(0,boundary::east,condition_type::neumann,&f);

    // source function, rhs
    gsConstantFunction<dual> g(0.,0.,2);

    //=============================================//
                  // Solving //
    //=============================================//

    gsNeoHookeLogMaterial<dual> materialMat(youngsModulus,poissonsRatio,2);
    // gsLinearMaterial<dual> materialMat(youngsModulus,poissonsRatio);

    // creating assembler
    gsElasticityAssembler<dual> assembler(geometry,basisDisplacement,bcInfo,g,&materialMat);
    gsInfo << "Initialized system with " << assembler.numDofs() << " dofs.\n";

    // setting Newton's method
    gsIterative<dual> solver(assembler);
    solver.options().setInt("Verbosity",solver_verbosity::all);
    solver.options().setInt("Solver",linear_solver::LDLT);

    gsInfo << "Solving...\n";
    gsStopwatch clock;
    clock.restart();
    solver.solve();
    gsInfo << "Solved the system in " << clock.stop() <<"s.\n";

    //=============================================//
                  // Output //
    //=============================================//

    // solution to the nonlinear problem as an isogeometric displacement field
    gsMultiPatch<dual> displacement;
    assembler.constructSolution(solver.solution(),solver.allFixedDofs(),displacement);
    // gsPiecewiseFunction<dual> stresses;
    // assembler.constructCauchyStresses(displacement,stresses,stress_components::von_mises);

//     if (numPlotPoints > 0)
//     {
//         // constructing an IGA field (geometry + solution)
//         gsField<dual> displacementField(assembler.patches(),displacement);
//         gsField<dual> stressField(assembler.patches(),stresses,true);
//         // creating a container to plot all fields to one Paraview file
//         std::map<std::string,const gsField<dual> *> fields;
//         fields["Displacement"] = &displacementField;
//         fields["von Mises"] = &stressField;
//         gsWriteParaviewMultiPhysics(fields,"cooks",numPlotPoints);
//         gsInfo << "Open \"cooks.pvd\" in Paraview for visualization.\n";
//     }

    // // validation
    // gsMatrix<dual> A(2,1);
    // A << 1.,1.;
    // A = displacement.patch(0).eval(A);
    // gsInfo << "X-displacement of the top-right corner: " << A.at(0) << std::endl;
    // gsInfo << "Y-displacement of the top-right corner: " << A.at(1) << std::endl;

    // Derivative of solution w.r.t. Poisson's ratio
    gsInfo << "\nComputing derivative of displacement w.r.t. Poisson's ratio:\n";
    gsMatrix<double> ddisplacement_dnu(displacement.patch(0).coefs().rows(), displacement.patch(0).coefs().cols());
    for (index_t i = 0; i < displacement.patch(0).coefs().rows(); ++i)
    {
        for (index_t j = 0; j < displacement.patch(0).coefs().cols(); ++j)
        {
            ddisplacement_dnu(i, j) = autodiff::detail::derivative<1>(displacement.patch(0).coefs()(i, j));
        }
    }
    gsInfo << "Derivative of displacement coefficients w.r.t. Poisson's ratio:\n";
    gsInfo << ddisplacement_dnu << "\n";

    // gsInfo<<rhs.transpose().cast<double>()<<"\n";
    // gsInfo<<"dRHS/dtheta (forward-mode):\n";
    // for (index_t i=0; i<rhs.rows(); ++i)
    // {
    //     gsInfo<<autodiff::detail::derivative<1>(rhs(i,0))<<" ";
    // }
    // gsInfo<<"\n\n";

    return 0;
}
