/** @file r_allen_cahn_example.cpp
    @brief Allen-Cahn equation on a moving spline mesh using an ALE formulation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. BAHARI
*/

#include <gismo.h>
#include <gsAssembler/gsAdaptiveMultiPatchBuilder.h>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false, errorsave = false, export_b64 = false, manufactured = false;
    index_t numRefine = 6, numElevate = 1, numSteps = 30, maxIter = 25, legacySteps = -1;
    index_t mappingUpdates = 3, adaptRefCrit = 1, initialId = 1995, plotEvery = 0;
    real_t intensity = 12., adaptRefParam = .5, dt = 1e-5;
    real_t epsilon = 0.007389228264793657, meshRelax = .25, meshCfl = .25, meshTimeScale = 1e-4;
    std::string filename, initialExpression;

    gsCmdLine cmd("Allen-Cahn: u_t = Delta u - (u^3-u)/epsilon^2, with homogeneous Neumann boundaries.");
    cmd.addInt("u", "uniformRefine", "Uniform refinements of the solution and moving-mesh bases", numRefine);
    cmd.addInt("e", "degreeElevation", "Degree elevation before refinement", numElevate);
    cmd.addInt("l", "numLRefine", "Number of physical time steps (no h-refinement in time)", legacySteps);
    cmd.addInt("", "steps", "Number of physical time steps", numSteps);
    cmd.addInt("i", "iter", "Initial MMPDE iteration budget", maxIter);
    cmd.addInt("m", "mappingUpdates", "Warm-started MMPDE updates per time step (0: freeze the mesh)", mappingUpdates);
    cmd.addInt("r", "adaptRefCrit", "Marking criterion [1:GARU,2:PUCA,3:BULK,4:PBULK]", adaptRefCrit);
    cmd.addReal("a", "adaptRefParam", "Marking parameter for the gradient indicator", adaptRefParam);
    cmd.addReal("f", "IntensityMAE", "Monitor intensity (1: fixed uniform mesh)", intensity);
    cmd.addReal("t", "dt", "Physical time step", dt);
    cmd.addReal("", "epsilon", "Diffuse-interface width", epsilon);
    cmd.addReal("", "meshRelax", "Maximum fraction of the proposed mesh update", meshRelax);
    cmd.addReal("", "meshTimeScale", "Physical relaxation time used to scale mesh motion with dt", meshTimeScale);
    cmd.addReal("", "meshCfl", "Maximum mesh displacement relative to the smallest sampled element altitude", meshCfl);
    cmd.addString("d", "file", "Optional XML file containing the initial condition", filename);
    cmd.addInt("", "initialId", "Initial-condition function ID in the optional XML file", initialId);
    cmd.addString("", "initial", "Initial condition as a function of x,y; default: a circular phase field", initialExpression);
    cmd.addInt("", "plotEvery", "Output interval in steps (0: approximately every 1e-4 physical time units)", plotEvery);
    cmd.addSwitch("manufactured", "Stationary cosine solution with matching forcing, for verification", manufactured);
    cmd.addSwitch("plot", "Write a ParaView time series", plot);
    cmd.addSwitch("export_b64", "Use base64 encoding for ParaView", export_b64);
    cmd.addSwitch("errorsave", "Save energy, mesh and nonlinear-solver diagnostics to CSV", errorsave);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    if (legacySteps >= 0) numSteps = legacySteps;

    GISMO_ENSURE(std::isfinite(dt) && dt > 0. && std::isfinite(epsilon) && epsilon > 0.,
                 "dt and epsilon must be finite and positive.");
    GISMO_ENSURE(numRefine >= 0 && numElevate >= 0 && numSteps >= 0 && maxIter >= 0
                 && mappingUpdates >= 0 && plotEvery >= 0 && legacySteps >= -1, "Iteration and refinement counts must be nonnegative.");
    GISMO_ENSURE(std::isfinite(intensity) && intensity >= 1. && meshRelax >= 0. && meshRelax <= 1.
                 && meshCfl > 0. && std::isfinite(meshCfl) && meshTimeScale > 0. && std::isfinite(meshTimeScale),
                 "Use intensity >= 1, meshRelax in [0,1], and finite meshCfl, meshTimeScale > 0.");
    if (plotEvery == 0)
        plotEvery = std::max<index_t>(1, static_cast<index_t>(std::ceil(1e-4 / dt)));

    std::ostringstream width;
    width << std::setprecision(17) << epsilon;
    gsFunctionExpr<> initial, forcing("0", 2);
    if (manufactured)
    {
        initialExpression = "0.25*cos(pi*x)*cos(pi*y)";
        const std::string u = "(" + initialExpression + ")";
        forcing = gsFunctionExpr<>("2*pi*pi*" + u + "+(" + u + "*" + u + "*" + u + "-" + u
                                    + ")/(" + width.str() + "*" + width.str() + ")", 2);
    }
    if (!initialExpression.empty())
        initial = gsFunctionExpr<>(initialExpression, 2);
    else if (!filename.empty())
        gsFileData<>(filename).getId(initialId, initial);
    else
        // initial = gsFunctionExpr<>("tanh((0.2-sqrt((x-0.5)^2+(y-0.5)^2))/(sqrt(2)*" + width.str() + "))", 2);
        initial = gsFunctionExpr<>("-1.+ tanh( ( 0.4 - sqrt((x-0.5)**2+ (y-0.5)**2))/(sqrt(2)*" + width.str() + ")) - tanh( ( 0.3 - sqrt((x-0.5)**2+ (y-0.5)**2))/(sqrt(2)*" + width.str() + "))", 2);
    GISMO_ENSURE(initial.domainDim() == 2 && initial.targetDim() == 1, "The initial condition must be a scalar function of x,y.");

    // The physical domain is the unit square; r-refinement moves only its interior.
    gsMultiPatch<> reference = gsNurbsCreator<>::BSplineSquareGrid(1, 1);
    reference.degreeElevate(numElevate);
    reference.computeTopology();
    gsMultiPatch<> Psi = reference;
    for (index_t r = 0; r < numRefine; ++r)
        Psi.uniformRefine();
    GISMO_ENSURE(Psi.basis(0).component(0).size() > 2 && Psi.basis(0).component(1).size() > 2,
                 "Increase uniformRefine or degreeElevation to provide interior mesh coefficients.");
    gsMultiBasis<> basis(Psi, true);
    // /***********************************************************************************************  **/
    // The moving-mesh basis is a copy of the solution basis, but with the boundary control points fixed.
    gsAdaptiveMultiPatchBuilder MAE(reference, numRefine, maxIter, intensity, 0, 0, true);

    gsExprAssembler<> A(1,1);
    A.options().setSwitch("SameElement", false);
    A.options().setReal("quA", 2.);
    A.setIntegrationElements(basis);
    gsExprEvaluator<> ev(A);
    auto v = A.getSpace(basis);
    gsBoundaryConditions<> bc; // Natural zero Neumann data require no boundary load.
    v.setup(bc, dirichlet::homogeneous, 0);
    auto G = A.getMap(Psi);
    auto initialValue = A.getCoeff(initial, G);
    auto source = A.getCoeff(forcing, G);
    gsMatrix<> coefficients;
    auto u = A.getSolution(v, coefficients);

    auto projectInitial = [&]()
    {
        A.initSystem();
        A.assemble(v * v.tr() * meas(G), v * initialValue * meas(G));
        gsSparseSolver<>::SimplicialLDLT massSolver;
        massSolver.compute(A.matrix());
        GISMO_ENSURE(massSolver.info() == gsEigen::Success, "Initial projection factorization failed.");
        coefficients = massSolver.solve(A.rhs());
        GISMO_ENSURE(massSolver.info() == gsEigen::Success && coefficients.allFinite(), "Initial projection failed.");
    };

    auto updateMonitor = [&]()
    {
        ev.integralElWise(igrad(u, G).sqNorm() * meas(G));
        const std::vector<real_t> indicators = ev.elementwise();
        // A spatially constant phase should not move the mesh due to roundoff.
        if (std::accumulate(indicators.begin(), indicators.end(), real_t(0))
            <= 1e-24 * std::max(real_t(1), coefficients.squaredNorm()))
            return false;
        std::vector<bool> marked(indicators.size());
        gsMarkElementsForRef(indicators, adaptRefCrit, adaptRefParam, marked);
        MAE.buildMultiPatch(MAE.buildDensity(basis, marked));
        return true;
    };

    // Sample each knot span at endpoints and midpoints. This is a quality check,
    // not an algebraic certificate of bijectivity.
    gsVector<index_t> counts(2);
    counts << 2*basis.basis(0).component(0).numElements()+1,
              2*basis.basis(0).component(1).numElements()+1;
    const gsMatrix<> meshGrid = gsGridIterator<real_t,CUBE>(basis.basis(0).support(), counts).toMatrix();
    auto meshQuality = [&](const gsMultiPatch<> & mesh)
    {
        const gsMatrix<> J = mesh.patch(0).jacobian(meshGrid);
        real_t minJac = std::numeric_limits<real_t>::infinity();
        real_t minSize = minJac;
        for (index_t q = 0; q < meshGrid.cols(); ++q)
        {
            const gsMatrix<> block = J.middleCols(2*q,2);
            const real_t det = block.determinant();
            if (!block.allFinite() || det <= 0.)
                return std::make_pair(real_t(-1), real_t(0));
            minJac = std::min(minJac, det);
            minSize = std::min(minSize, det / block.col(1).norm()
                              / basis.basis(0).component(0).numElements());
            minSize = std::min(minSize, det / block.col(0).norm()
                              / basis.basis(0).component(1).numElements());
        }
        return std::make_pair(minJac, minSize);
    };

    auto acceptMesh = [&](const gsMultiPatch<> & candidate, real_t fraction, bool limitDisplacement)
    {
        const gsMatrix<> previous = Psi.patch(0).coefs();
        GISMO_ENSURE(candidate.patch(0).coefs().rows() == previous.rows(), "Moving-mesh basis changed during time integration.");
        const gsMatrix<> increment = candidate.patch(0).coefs() - previous;
        const real_t displacement = (candidate.patch(0).eval(meshGrid)-Psi.patch(0).eval(meshGrid)).colwise().norm().maxCoeff();
        const auto quality = meshQuality(Psi);
        GISMO_ENSURE(quality.first > 0., "The current mesh is invalid.");
        if (limitDisplacement && displacement > 0.)
            fraction = std::min(fraction, meshCfl * quality.second / displacement);
        bool accepted = false;
        for (index_t k = 0; k < 24; ++k)
        {
            Psi.patch(0).coefs() = previous + fraction * increment;
            if (meshQuality(Psi).first >= .05 * quality.first)
            {
                accepted = true;
                break;
            }
            fraction *= .5;
        }
        if (!accepted)
        {
            Psi.patch(0).coefs() = previous;
            fraction = 0.;
        }
        MAE.MAmapping = Psi; // The next MMPDE update starts from the accepted mesh.
        return fraction;
    };

    projectInitial();
    if (intensity > 1. && meshRelax > 0.)
    {
        if (updateMonitor())
        {
            acceptMesh(MAE.MAmapping, 1., false);
            projectInitial(); // Initial data are known analytically, so reproject once.
        }
    }
    MAE.MAmapping = Psi;
    MAE.m_maxIter = mappingUpdates;

    // These objects stay alive and nonempty: expression handles reference them.
    gsMultiPatch<> previousMesh = Psi;
    gsMatrix<> previousCoefficients = coefficients;
    auto Gold = A.getMap(previousMesh);
    auto uold = A.getSolution(v, previousCoefficients);
    const real_t reactionScale = dt / (epsilon*epsilon);
    // Backward Euler in moving coordinates:
    // U_t - w.grad(u) = Delta(u) - (u^3-u)/epsilon^2 + f,
    // with dt*w = X_new-X_old. U_old uses the OLD coefficients at the SAME xi.
    auto residual = (v * (u-uold) + dt * igrad(v,G) * igrad(u,G).tr()
                    + reactionScale * v * (u*u*u-u).val()
                    - v * (igrad(u,G) * (G-Gold)).tr() - dt * v * source) * meas(G);
    auto jacobian = (v * v.tr() + dt * igrad(v,G) * igrad(v,G).tr()
                    + reactionScale * v * v.tr() * (3.*(u*u).val()-1.)
                    - v * (igrad(v,G) * (G-Gold)).tr()) * meas(G);

    gsParaviewCollection<>::uPtr collection;
    if (plot)
    {
        collection.reset(new gsParaviewCollection<>("ParaviewOutput/time_solution", ev));
        collection->options().setSwitch("elements", true);
        collection->options().setSwitch("base64", export_b64);
        collection->options().setInt("elementResolution", 8);
        collection->options().setInt("numPoints", 10000);
    }
    auto writeSnapshot = [&](real_t time)
    {
        if (!plot) return;
        collection->newTimeStep(Psi, time);
        collection->addField(u, "solution");
        collection->addField(jac(G).det(), "Jacobian");
        collection->saveTimeStep();
    };

    std::ofstream diagnostics;
    if (errorsave)
    {
        diagnostics.open("allen_cahn_diagnostics.csv");
        GISMO_ENSURE(diagnostics.is_open(), "Cannot open allen_cahn_diagnostics.csv.");
        diagnostics << "time,energy,min_u,max_u,min_jacobian,newton_iterations,residual,mesh_fraction";
        if (manufactured) diagnostics << ",L2_error";
        diagnostics << "\n" << std::setprecision(17);
    }
    auto freeEnergy = [&]()
    {
        return ev.integral((.5*igrad(u,G).sqNorm()
                                + pow((u*u).val()-1., 2.)/(4.*epsilon*epsilon)) * meas(G));
    };
    auto report = [&](real_t time, index_t iterations, real_t norm, real_t meshFraction)
    {
        const real_t energy = freeEnergy();
        const real_t lo = ev.min(u), hi = ev.max(u), minJ = meshQuality(Psi).first;
        GISMO_ENSURE(std::isfinite(energy) && minJ > 0., "Invalid solution or mesh.");
        gsInfo << std::scientific << "t=" << time << " energy=" << energy << " u=[" << lo << "," << hi
               << "] minJ=" << minJ << " Newton=" << iterations << " residual=" << norm
               << " mesh_fraction=" << meshFraction;
        if (errorsave)
            diagnostics << time << "," << energy << "," << lo << "," << hi << "," << minJ
                        << "," << iterations << "," << norm << "," << meshFraction;
        if (manufactured)
        {
            const real_t error = math::sqrt(ev.integral((u-initialValue).sqNorm()*meas(G)));
            gsInfo << " L2_error=" << error;
            if (errorsave) diagnostics << "," << error;
        }
        gsInfo << "\n";
        if (errorsave) diagnostics << "\n";
    };

    report(0., 0, 0., 0.);
    writeSnapshot(0.);
    for (index_t step = 1; step <= numSteps; ++step)
    {
        const real_t oldEnergy = freeEnergy();
        previousMesh = Psi;
        previousCoefficients = coefficients;
        real_t meshFraction = 0.;
        if (intensity > 1. && mappingUpdates > 0 && meshRelax > 0.)
        {
            // Evaluate the old solution on its old physical mesh.
            if (updateMonitor())
            {
                // Scaling with dt avoids an unbounded mesh velocity as dt tends to zero.
                const real_t fraction = meshRelax * std::min(real_t(1), dt/meshTimeScale);
                meshFraction = acceptMesh(MAE.MAmapping, fraction, true);
            }
        }

        const gsMultiPatch<> proposedMesh = Psi;
        const real_t proposedFraction = meshFraction;
        real_t norm = 0., tolerance = 0.;
        index_t iteration = 0;
        bool stepAccepted = false;
        for (index_t retry = 0; retry <= 10; ++retry)
        {
            if (retry > 0)
            {
                // Retry the entire ALE step from the unchanged old state.
                coefficients = previousCoefficients;
                Psi = previousMesh;
                const real_t blend = retry < 10 ? std::pow(.5, retry) : 0.;
                meshFraction = proposedFraction * acceptMesh(proposedMesh, blend, false);
            }
            iteration = 0;
            for (; iteration < 30; ++iteration)
            {
                A.initSystem();
                A.assemble(jacobian, residual);
                norm = A.rhs().norm();
                if (iteration == 0) tolerance = 1e-10 + 1e-8 * norm;
                if (norm <= tolerance) break;

                // Mesh advection makes the ALE matrix nonsymmetric.
                gsSparseSolver<>::LU newtonSolver;
                newtonSolver.compute(A.matrix());
                GISMO_ENSURE(newtonSolver.info() == gsEigen::Success, "Newton factorization failed.");
                const gsMatrix<> correction = newtonSolver.solve(A.rhs());
                GISMO_ENSURE(newtonSolver.info() == gsEigen::Success && correction.allFinite(), "Newton solve failed.");
                const gsMatrix<> iterate = coefficients;
                real_t damping = 1.;
                bool accepted = false;
                for (index_t backtrack = 0; backtrack < 20; ++backtrack)
                {
                    coefficients = iterate - damping * correction;
                    A.clearRhs();
                    A.assemble(residual);
                    const real_t trialNorm = A.rhs().norm();
                    if (std::isfinite(trialNorm) && trialNorm <= (1.-1e-4*damping)*norm)
                    {
                        accepted = true;
                        norm = trialNorm;
                        break;
                    }
                    damping *= .5;
                }
                GISMO_ENSURE(accepted, "Newton line search failed; reduce dt or meshRelax.");
            }
            GISMO_ENSURE(norm <= tolerance, "Newton did not converge; reduce dt or meshRelax.");
            // The unforced Neumann problem dissipates energy. Reduce mesh motion
            // if its first-order transport approximation introduces energy gain.
            if (manufactured || freeEnergy() <= oldEnergy + 1e-10*std::max(real_t(1), oldEnergy))
            {
                stepAccepted = true;
                break;
            }
            if (meshFraction == 0.) break;
        }
        GISMO_ENSURE(stepAccepted, "Energy increased even on a fixed mesh; reduce dt.");
        const real_t time = step * dt;
        report(time, iteration, norm, meshFraction);
        if (step % plotEvery == 0 || step == numSteps)
            writeSnapshot(time);
    }
    if (plot)
    {
        collection->save();
        gsInfo << "ParaView time series: ParaviewOutput/time_solution.pvd\n";
    }
    return EXIT_SUCCESS;
}
