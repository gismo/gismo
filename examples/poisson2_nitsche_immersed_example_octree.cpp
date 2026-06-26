/** @file poisson2_quarter_circle_immersed_example.cpp

    @brief Poisson on an immersed quarter-circle domain inside [0,1]^2.

    Embedding box: [0,1]x[0,1]
    Physical domain: Omega = {(x,y): x>0, y>0, x^2+y^2<1}
    ./build/bin/poisson2_nitsche_immersed_example_octree -r 2 -e 1 -L 3 --plot -g 1e8 -o output_immersed_octree -q 13

    This file is part of the G+Smo library.
*/

#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <fstream>
#include <string>

using namespace gismo;

int main(int argc, char *argv[])
{
    bool plot = false;
    index_t numRefine  = 5;
    index_t numElevate = 0;
    index_t octLevels  = 2;
    real_t penalty = 1e3;
    real_t gamma = 1e3;
    real_t alpha = 1e-10;
    index_t quRule = gsQuadrature::AlgoimRule;
    std::string outFolder = "output_immersed";

    gsCmdLine cmd("Poisson immersed quarter-circle in [0,1]^2.");
    cmd.addInt("e", "degreeElevation", "Number of degree elevation steps", numElevate);
    cmd.addInt("r", "uniformRefine", "Number of uniform refinement steps (before solve)", numRefine);
    cmd.addInt("L", "octLevels", "Number of octree subdivision levels", octLevels);
    cmd.addReal("p", "penalty", "Weak Dirichlet penalty", penalty);
    cmd.addSwitch("plot", "Create a ParaView file", plot);
    cmd.addReal("g","gamma","penalty parameter", gamma);
    cmd.addReal("a","alpha","FCM stabilization factor (alpha * full-domain stiffness)", alpha);
    cmd.addInt("q", "quRule", "Rule used for plotting only (all three are always compared): 11=CutCell, 12=Algoim, 13=Octree", quRule);
    cmd.addString("o", "output", "Output folder", outFolder);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    if (quRule != gsQuadrature::AlgoimRule && quRule != gsQuadrature::CutCellRule && quRule != gsQuadrature::OctreeRule)
    {
        gsWarn << "Unsupported quRule=" << quRule
               << ". Falling back to AlgoimRule (12).\n";
        quRule = gsQuadrature::AlgoimRule;
    }

    // Embedding box [0,1]^2
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineSquare());

    gsFunctionExpr<> u_exact("sin(pi*x)*sin(pi*y)", 2);
    gsFunctionExpr<> f_rhs("2*pi^2*sin(pi*x)*sin(pi*y)", 2);
    // gsFunctionExpr<> u_exact("if(-0.16 + (x-0.5)^2 + (y-0.5)^2<0,sin(pi*x)*sin(pi*y), 0)",2);


    // Quarter-circle level-set: inside when max(-x,-y,x^2+y^2-1) < 0
    // gsFunctionExpr<> impl_fun("max(max(-x,-y),x^2+y^2-1)", 2);
    // gsFunctionExpr<> impl_fun("1 - x^2 - y^2", 2);
    gsFunctionExpr<> impl_fun("-0.16 + (x-0.5)^2 + (y-0.5)^2", 2);
    // Square
    gsFunctionExpr<> impl_square("x*(x-1)*y*(y-1)",2);
    // gsFunctionExpr<> normal_immersed("2*(x-0.5)","2*(y-0.5)",2);
    gsFunctionExpr<> normal_immersed_unit(
        "2*(x-0.5)/sqrt(4*(x-0.5)^2 + 4*(y-0.5)^2)",
        "2*(y-0.5)/sqrt(4*(x-0.5)^2 + 4*(y-0.5)^2)",
        2
    );


    gsBoundaryConditions<> bc;
    // bc.addCondition(0, boundary::west, condition_type::weak_dirichlet, &u_exact);
    // bc.addCondition(0, boundary::south, condition_type::weak_dirichlet, &u_exact);
    // bc.addCondition(0, boundary::none,  condition_type::weak_dirichlet, &u_exact); // immersed interface
    bc.setGeoMap(mp);

    gsMultiBasis<> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);
    
    const index_t deg = dbasis.maxCwiseDegree();
    // const real_t gamma = (penalty < 0)? real_t(2.5) * (deg + 2) * (deg + 1): penalty;
    // const real_t gamma = 1e+8;


    #ifdef GISMO_WITH_PARDISO
        gsSparseSolver<>::PardisoLDLT solver;
    #else
        gsSparseSolver<>::LU solver;
    #endif

    // Compare all three immersed quadrature rules in a single run
    struct RuleInfo { std::string name; index_t id; };
    const std::vector<RuleInfo> rules = {
        { "CutCell", gsQuadrature::CutCellRule },
        { "Algoim",  gsQuadrature::AlgoimRule  },
        { "Octree",  gsQuadrature::OctreeRule  }
    };

    // Rows: rule index, Cols: refinement level
    gsMatrix<> l2err(rules.size(), numRefine + 1), h1err(rules.size(), numRefine + 1);
    gsInfo<< "Degree of the basis: " << dbasis.maxCwiseDegree() << "\n";

    std::string outPath = outFolder;
    if (outPath.empty())
        outPath = "output_immersed";
    const bool isAbsolutePath = (!outPath.empty() && outPath[0] == '/');
    if (!isAbsolutePath)
        outPath = gsFileManager::getCurrentPath() + "/" + outPath;

    std::string out = gsFileManager::getCanonicRepresentation(outPath);
    gsFileManager::mkdir(out);
    const std::string debugDir = out + "/trimmed_domain_debug";
    const std::string insideDir = debugDir + "/inside";
    const std::string boundaryDir = debugDir + "/boundarysign";
    if (plot)
    {
        gsFileManager::mkdir(debugDir);
        gsFileManager::mkdir(insideDir);
        gsFileManager::mkdir(boundaryDir);
    }

    auto makeRootPvd = [](const std::string & sourcePvd,
                          const std::string & targetPvd,
                          const std::string & prefix)
    {
        std::ifstream in(sourcePvd.c_str());
        if (!in.is_open())
            return;

        std::ofstream outFile(targetPvd.c_str());
        if (!outFile.is_open())
            return;

        std::string line;
        const std::string tag = "file=\"";
        while (std::getline(in, line))
        {
            size_t pos = line.find(tag);
            while (pos != std::string::npos)
            {
                pos += tag.size();
                line.insert(pos, prefix);
                pos = line.find(tag, pos + prefix.size());
            }
            outFile << line << "\n";
        }
    };

    for (size_t k = 0; k != rules.size(); ++k)
    {
        // Reset the basis so each rule starts from the same coarse mesh
        dbasis = gsMultiBasis<>(mp, true);
        dbasis.setDegree(dbasis.maxCwiseDegree() + numElevate);

        gsInfo << "\n=== Quadrature rule: " << rules[k].name
               << " (id " << rules[k].id << ") ===\n"
               << "(dot1=assembled, dot2=solved, dot3=error)\nDoFs: ";

    for(int r = 0; r <= numRefine; ++r)
    {        
        dbasis.uniformRefine();

        // Determine maximum mesh size
        real_t hmax = 0;
        for (size_t p=0; p!=dbasis.nBases(); p++)
        {
            hmax = math::max(hmax, dbasis.basis(p).getMaxCellLength());
        }

        gsTensorBSplineBasis<2,real_t> * tbsPtr =
            dynamic_cast<gsTensorBSplineBasis<2,real_t>*>(&dbasis.basis(0));
        GISMO_ENSURE(tbsPtr, "Basis is not a tensor B-spline basis");

        // gsImplicitTrimmedDomain<2,real_t> tr_domain(impl_fun, *tbsPtr);  // DANGLING POINTER - COMMENTED OUT!

        typedef gsExprAssembler<>::geometryMap geometryMap;
        typedef gsExprAssembler<>::space       space;
        typedef gsExprAssembler<>::solution    solution;

        gsExprAssembler<> A(1, 1);
        A.options().setInt("quRule", rules[k].id);
        A.options().addInt("octLevels", "Number of octree subdivision levels for OctreeRule", octLevels);
        // Trimmed (implicit) integration domain. It is applied AFTER the FCM
        // stabilization term below, so that the stabilization is integrated over
        // the full background mesh while the physical terms use the trimmed domain.
        memory::shared_ptr<gsImplicitTrimmedDomain<2,real_t> > tr_domain =
            memory::make_shared(new gsImplicitTrimmedDomain<2,real_t>(impl_fun, *tbsPtr));
        gsExprEvaluator<> ev(A);

        std::vector<patchSide> bdr_immersed(1); // vector of size 1
        bdr_immersed[0]= patchSide(0,boundary::none); // assign the first and only element

        geometryMap G = A.getMap(mp);
        space u = A.getSpace(dbasis);

        auto ff    = A.getCoeff(f_rhs, G);
        auto n_imm = A.getCoeff(normal_immersed_unit, G); // grad(phi) — normalized?
        auto u_ex = ev.getVariable(u_exact, G);
        auto impl = A.getCoeff(impl_fun,G);
        gsMatrix<> solVector;
        solution u_sol = A.getSolution(u, solVector);

        u.setup(bc, dirichlet::none, 0);
        // Start on the FULL background mesh so the FCM stabilization below is
        // integrated over all cells (including exterior/cut ones).
        A.setIntegrationElements(dbasis);
        A.initSystem();
        A.computePattern(igrad(u) * igrad(u).tr());
        // A.getCoeff()

        // --- FCM stabilization ---------------------------------------------
        // Trimmed cut-cell/octree/Algoim integration leaves basis functions on
        // cut/exterior cells with (near-)zero stiffness rows, making the system
        // singular/ill-conditioned (small cut-cell instability). Assemble a tiny
        // multiple alpha of the FULL background-mesh stiffness FIRST (standard
        // Gauss, before restricting to the trimmed domain) so those DoFs keep a
        // well-conditioned coupling to their neighbours. Perturbs the physical
        // solution only by O(alpha).
        A.options().setInt("quRule", gsQuadrature::GaussLegendre);
        A.assemble( alpha * igrad(u, G) * igrad(u, G).tr() * meas(G) );

        // Restrict integration to the trimmed (implicit) physical domain for the
        // remaining (physical) terms, using the selected immersed quadrature rule.
        A.setIntegrationDomain(tr_domain);
        A.options().setInt("quRule", rules[k].id);

        A.assemble(
            igrad(u, G) * igrad(u, G).tr() * meas(G),
            u * ff * meas(G)
        );

        auto g_D = A.getCoeff(u_exact, G); // prescribe solution at the immersed boundary!

        // Symmetric Nitsche terms (matrix part)
        // A.assembleBdr(
        //     bdr_immersed,
        //     - (igrad(u, G) * igrad(impl,G)) * u.tr()
        //     - u * (igrad(u, G) * igrad(impl,G)).tr()
        //     + gamma / hmax * u * u.tr() * meas(G)
        // );

        A.assembleBdr(
            bdr_immersed,
            - (igrad(u, G) * n_imm) * u.tr()
            - u * (igrad(u, G) * n_imm).tr()
            + gamma / hmax * u * u.tr() * meas(G)
        );

        // Symmetric Nitsche terms (rhs part)
        // unit outward normal: n_imm.normalized()
        // A.assembleBdr(
        //     bdr_immersed,
        //     - (igrad(u, G) * igrad(impl,G)) * g_D * meas(G)
        //     + gamma / hmax * u * g_D * meas(G)
        // );

        // if i use igrad(impl,G), I get some issues in the computation! 

        A.assembleBdr(
            bdr_immersed,
            - (igrad(u, G) * n_imm) * g_D * meas(G)
            + gamma / hmax * u * g_D * meas(G)
        );

        gsInfo << A.numDofs() << "." << std::flush;

        solver.compute(A.matrix());
        solVector = solver.solve(A.rhs());
        gsInfo << "." << std::flush;

        l2err(k, r) = math::sqrt(ev.integral((u_ex - u_sol).sqNorm() * meas(G)));
        h1err(k, r) = l2err(k, r) + math::sqrt(ev.integral((igrad(u_ex) - igrad(u_sol, G)).sqNorm() * meas(G)));
        gsInfo << ". " << std::flush;
    
        if (plot && r == numRefine && rules[k].id == quRule)
        {
            ev.options().setSwitch("plot.elements", true);
            ev.options().setInt("plot.npts", 100000);
            ev.writeParaview(u_sol, G, out+"/poisson2_immersed_circle");
            ev.writeParaview(u_ex, G, out+"/poisson2_exact_solution");
            // ===== QUADRATURE POINT VISUALIZATION =====
            {
                gsMatrix<real_t> quad_points(4, 0);
                gsMatrix<real_t> pts;
                gsVector<real_t> wts;
                const gsImplicitTrimmedDomain<2,real_t> * tr_domain =
                    dynamic_cast<const gsImplicitTrimmedDomain<2,real_t>*>(&A.domain());
                GISMO_ENSURE(tr_domain, "Expected gsImplicitTrimmedDomain for quadrature export.");

                auto l2_error_point = (u_ex - u_sol).norm();
                auto quRuleVis = gsQuadrature::getPtr(*tr_domain, A.options());

                for (auto elem_it = tr_domain->beginAll(); elem_it != tr_domain->endAll(); ++elem_it)
                {
                    quRuleVis->mapTo(elem_it.lowerCorner(), elem_it.upperCorner(), pts, wts);
                    if (pts.cols() == 0) continue;

                    gsMatrix<real_t> phys;
                    mp.patch(0).eval_into(pts, phys);

                    const index_t c = quad_points.cols();
                    quad_points.conservativeResize(4, c + pts.cols());
                    quad_points.block(0, c, 2, pts.cols()) = phys;
                    quad_points.row(2).segment(c, pts.cols()).setZero();
                    for (index_t j = 0; j < pts.cols(); ++j)
                        quad_points(3, c + j) = ev.eval(l2_error_point, pts.col(j), 0)(0,0);
                }

                if (quad_points.cols() > 0)
                {
                    const std::string base = out + "/quadrature_points_tr_domain";
                    gsInfo << "Exporting " << quad_points.cols()
                           << " trimmed-domain quadrature points to " << base << ".vtp\n";
                    gsWriteParaviewPoints(quad_points, base);
                }
            }
            // ==========================================

// L2 error
            auto l2_error_field = (u_ex - u_sol).norm();
            ev.writeParaview(l2_error_field, G, out+"/poisson2_immersed_L2_error");
            
            // Compute H1 seminorm error
            auto h1_error_field = (igrad(u_ex, G) - igrad(u_sol, G)).sqNorm();
            ev.writeParaview(h1_error_field, G, out+"/poisson2_immersed_H1_error");

            // Plot background mesh
            gsMesh<> mesh(dbasis.basis(0));
            gsWriteParaview(mesh,out+"/background_mesh");
        }

        if (plot && r == numRefine && rules[k].id == quRule)
        {
            // Debug export of trimmed-domain element classification per refinement level:
            // sign < 0 : interior cell, sign == 0 : cut/boundary cell.
            gsMatrix<> boxesInterior(2,0), boxesBoundary(2,0);

            auto appendBox = [](gsMatrix<> & boxes,
                                const gsVector<> & lo,
                                const gsVector<> & hi)
            {
                const index_t c = boxes.cols();
                boxes.conservativeResize(2, c + 2);
                boxes.col(c)     = lo;
                boxes.col(c + 1) = hi;
            };

            const gsImplicitTrimmedDomain<2,real_t> * tr_domain =
                dynamic_cast<const gsImplicitTrimmedDomain<2,real_t>*>(&A.domain());
            GISMO_ENSURE(tr_domain, "Expected gsImplicitTrimmedDomain for trimmed-domain export.");

            // Interior elements only (sign < 0)
            for (auto it = tr_domain->beginInterior(); it != tr_domain->end<InteriorSign>(); ++it)
                appendBox(boxesInterior, it.lowerCorner(), it.upperCorner());

            // Cut elements on immersed boundary (sign == 0)
            for (auto it = tr_domain->beginBdr(boundary::none); it != tr_domain->endBdr(boundary::none); ++it)
                appendBox(boxesBoundary, it.lowerCorner(), it.upperCorner());

            if (boxesInterior.cols() > 0)
            {
                const std::string insideBase = insideDir + "/mesh_inside";
                gsWriteParaview(boxesInterior, insideBase, real_t(1));
                makeRootPvd(insideBase + ".pvd", debugDir + "/inside_collection.pvd", "inside/");
            }

            if (boxesBoundary.cols() > 0)
            {
                const std::string boundaryBase = boundaryDir + "/mesh_cut";
                gsWriteParaview(boxesBoundary, boundaryBase, real_t(2));
                makeRootPvd(boundaryBase + ".pvd", debugDir + "/boundarysign_collection.pvd", "boundarysign/");
            }


        }

        gsInfo << "Penalty: " << gamma/hmax << "\n";

    } // refinement loop
    } // quadrature-rule loop

    // ---- Side-by-side comparison of the three quadrature rules ----
    gsInfo << "\n\n===== Convergence comparison =====\n";
    for (size_t k = 0; k != rules.size(); ++k)
    {
        gsInfo << "\n--- " << rules[k].name << " ---\n";
        gsInfo << "L2 error: " << std::scientific << std::setprecision(3)
               << l2err.row(k) << "\n";
        gsInfo << "H1 error: " << std::scientific << h1err.row(k) << "\n";

        if (numRefine > 0)
        {
            gsInfo << "EoC (L2): " << std::fixed << std::setprecision(2)
                   << (l2err.row(k).head(numRefine).array() / l2err.row(k).tail(numRefine).array()).log().matrix() / std::log(2.0)
                   << "\n";
            gsInfo << "EoC (H1): " << std::fixed << std::setprecision(2)
                   << (h1err.row(k).head(numRefine).array() / h1err.row(k).tail(numRefine).array()).log().matrix() / std::log(2.0)
                   << "\n";
        }
    }

    return EXIT_SUCCESS;
}
