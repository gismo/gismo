/** @file periodicBSpline_example.cpp

    @brief L2-projection convergence study on periodic and non-periodic
    B-spline spaces.

    A known function is L2-projected onto uniform B-spline spaces of degree p
    with N elements on [0,1] (d = 1) or [0,1]^2 (d = 2), and the L2 error and
    the estimated order of convergence (EoC) between successive levels are
    printed, one table row per (variant, p, N).

    Variants:
      - base : clamped, non-periodic basis,
      - C1   : periodic basis built by the periodic constructor,
      - C2   : periodic basis built from the uniform stretched knot vector
               t_k = k*h - p*h, k = 0,...,N+2p, followed by setPeriodic().
      - base-x, C1-x, C2-x : the same spaces, applied to f(x) = x (1D only).
    The C1 variants run by default; --withC1 skips them.

    Source functions:
      - 1D: f(x) = sin(2 pi x) + cos(4 pi x)/2. It is smooth and periodic,
        so a correct periodic space converges with order p+1 in L2.
      - 1D "-x" variants: f(x) = x. Its periodic extension jumps at the seam,
        so a periodic space converges only with order 1/2. A non-periodic
        basis reproduces x exactly, so its error is at round-off level.
      - 2D: f(x,y) = sin(2 pi x)(1+y^2), periodic in x (direction 0) and
        clamped in y; it is run with --with2d or an explicit --case ...,2.

    The projection is assembled without any geometry map: the parameter
    domain is the physical domain, so the measure of the map is 1. No
    geometry is built from the periodic bases; the projection is evaluated
    only through its coefficients. The integration domain, for the assembly
    and for the error, is the domain [0,1] (resp. [0,1]^2) of a separate
    non-periodic reference basis of degree p, and not the domain reported by
    the basis under test. A periodic basis whose own domain differs from
    [0,1] therefore shows up as a large error and a mismatching 'domain'
    column rather than as a silently different problem.

    Output columns: variant, p, N, domain (support of the basis under test),
    isPeriodic, size, m (numCrossingFunctions(); of component 0 in 2D),
    L2error, EoC, and a status field (solve status, timings, or the error
    message of a caught exception). EoC = log2(e_{l-1}/e_l) since N doubles
    per level; it is '-' at the first level, after a failed or skipped level,
    and when one of the errors is at round-off (< 1e-12).

    Expected EoC for a correct periodic space and the smooth function: p+1.
    Expected EoC for the -x variants of the periodic spaces: about 0.5.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include <gismo.h>

#include <memory>
#include <sstream>
#include <iomanip>

using namespace gismo;

namespace
{

/// Uniform stretched knot vector on [0,1] with N elements and degree p,
/// t_k = (k - p)/N for k = 0,...,N+2p (N+2p+1 knots, domain [t_p, t_{N+p}]).
gsKnotVector<real_t> stretchedKnots(const index_t N, const short_t p)
{
    std::vector<real_t> knots(N + 2 * p + 1);
    const real_t h = 1.0 / N;
    for (index_t k = 0; k < static_cast<index_t>(knots.size()); ++k)
        knots[k] = (k - p) * h;
    return gsKnotVector<real_t>(knots, p);
}

/// Clamped uniform knot vector on [0,1] with N elements and degree p.
gsKnotVector<real_t> clampedKnots(const index_t N, const short_t p)
{
    return gsKnotVector<real_t>(0.0, 1.0, N - 1, p + 1, 1, p);
}

std::string domainString(const gsMatrix<real_t> & s)
{
    std::ostringstream os;
    os << std::setprecision(6);
    for (index_t d = 0; d < s.rows(); ++d)
        os << (d ? "x[" : "[") << s(d, 0) << "," << s(d, 1) << "]";
    return os.str();
}

std::string sci(const real_t v)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(6) << v;
    return os.str();
}

/// Projects \a fun onto \a basis in L2 and returns the L2 error, integrated
/// over the domain of \a ref. \a ref must have degree p and outlive the call.
/// Cost: one sparse LU of the mass matrix of size basis.size().
real_t projectionError(const gsBasis<real_t> & basis,
                       const gsBasis<real_t> & ref,
                       const gsFunctionExpr<real_t> & fun,
                       bool & solveOk,
                       std::string & timings)
{
    typedef gsExprAssembler<real_t>::space    space;
    typedef gsExprAssembler<real_t>::solution solution;

    gsStopwatch sw;

    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationDomain(ref.domain());
    space u = A.getSpace(basis, 1);
    auto f = A.getCoeff(fun);

    sw.restart();
    A.initSystem();
    A.assemble(u * u.tr());
    A.assemble(u * f);
    const double tAsm = sw.stop();

    sw.restart();
    gsSparseSolver<real_t>::LU solver;
    solver.compute(A.matrix());
    gsMatrix<real_t> coefs = solver.solve(A.rhs());
    solveOk = solver.succeed();
    const double tSolve = sw.stop();

    sw.restart();
    solution sol = A.getSolution(u, coefs);
    gsExprEvaluator<real_t> ev(A);
    ev.options().setInt("quB", 3);
    const real_t err = math::sqrt(ev.integral((sol - f).sqNorm()));
    const double tErr = sw.stop();

    std::ostringstream os;
    os << "tAsm=" << tAsm << "s tSolve=" << tSolve << "s tErr=" << tErr << "s";
    timings = os.str();
    return err;
}

struct Case
{
    std::string variant; // base, C1, C2 (without the -x suffix)
    bool        linear;  // use f(x) = x
    index_t     dim;
    short_t     p;
    std::string label() const
    {
        return variant + (linear ? "-x" : "") + (dim == 2 ? "/2d" : "");
    }
};

void runCase(const Case & c, const index_t N0, const index_t lmax)
{
    const short_t p = c.p;
    gsFunctionExpr<real_t> fun;
    if (c.dim == 1)
        fun = gsFunctionExpr<real_t>(c.linear ? "x" : "sin(2*pi*x)+0.5*cos(4*pi*x)", 1);
    else
        fun = gsFunctionExpr<real_t>("sin(2*pi*x)*(1+y^2)", 2);

    real_t prevErr = -1;
    bool prevValid = false;

    for (index_t l = 0; l <= lmax; ++l)
    {
        const index_t N = N0 << l;
        std::ostringstream row;
        row << std::left << std::setw(10) << c.label() << " "
            << std::setw(2) << p << " " << std::setw(4) << N << " ";

        if (N < p + 1)
        {
            gsInfo << row.str() << "skipped (n<p+1)\n" << std::flush;
            prevValid = false;
            continue;
        }

        std::string dom = "?", per = "?", sz = "?", m = "?";
        std::string errStr = "ERROR", eocStr = "-", status;
        bool valid = false;
        real_t err = 0;
        try
        {
            std::unique_ptr<gsBasis<real_t> > basis;
            std::unique_ptr<gsBasis<real_t> > ref;

            if (c.dim == 1)
            {
                ref.reset(new gsBSplineBasis<real_t>(0.0, 1.0, N - 1, p));
                gsBSplineBasis<real_t> * b = nullptr;
                if (c.variant == "base")
                    b = new gsBSplineBasis<real_t>(0.0, 1.0, N - 1, p);
                else if (c.variant == "C1")
                    b = new gsBSplineBasis<real_t>(0.0, 1.0, N - 1, p, 1, true);
                else
                {
                    b = new gsBSplineBasis<real_t>(stretchedKnots(N, p));
                    b->setPeriodic();
                }
                basis.reset(b);
                dom = domainString(b->support());
                per = util::to_string(b->isPeriodic());
                sz  = util::to_string(b->size());
                m   = util::to_string(b->numCrossingFunctions());
            }
            else
            {
                typedef gsTensorBSplineBasis<2, real_t> Tb;
                ref.reset(new Tb(clampedKnots(N, p), clampedKnots(N, p)));
                Tb * b = nullptr;
                if (c.variant == "base")
                    b = new Tb(clampedKnots(N, p), clampedKnots(N, p));
                else if (c.variant == "C1")
                    b = new Tb(new gsBSplineBasis<real_t>(0.0, 1.0, N - 1, p, 1, true),
                               new gsBSplineBasis<real_t>(0.0, 1.0, N - 1, p));
                else
                {
                    b = new Tb(stretchedKnots(N, p), clampedKnots(N, p));
                    b->setPeriodic(0);
                }
                basis.reset(b);
                dom = domainString(b->support());
                per = util::to_string(b->isPeriodic());
                sz  = util::to_string(b->size());
                m   = util::to_string(b->component(0).numCrossingFunctions());
            }

            bool solveOk = false;
            std::string timings;
            err = projectionError(*basis, *ref, fun, solveOk, timings);
            errStr = sci(err);
            valid = solveOk;
            status = std::string("solve=") + (solveOk ? "ok" : "FAIL") + " " + timings;
            if (prevValid && valid && prevErr > 1e-12 && err > 1e-12)
            {
                std::ostringstream os;
                os << std::fixed << std::setprecision(3) << std::log(prevErr / err) / std::log(2.0);
                eocStr = os.str();
            }
        }
        catch (std::exception & e)
        {
            valid = false;
            status = std::string("exception: ") + e.what();
        }

        row << std::setw(22) << dom << " " << std::setw(10) << per << " "
            << std::setw(5) << sz << " " << std::setw(2) << m << " "
            << std::setw(13) << errStr << " " << std::setw(7) << eocStr << " "
            << status;
        gsInfo << row.str() << "\n" << std::flush;

        prevValid = valid;
        prevErr = err;
    }
}

bool parseCase(const std::string & spec, Case & c)
{
    std::vector<std::string> tok;
    std::stringstream ss(spec);
    std::string t;
    while (std::getline(ss, t, ','))
        tok.push_back(t);
    if (tok.size() != 3)
        return false;

    std::string v = tok[0];
    c.linear = false;
    if (v.size() > 2 && v.compare(v.size() - 2, 2, "-x") == 0)
    {
        c.linear = true;
        v.resize(v.size() - 2);
    }
    if (v != "base" && v != "C1" && v != "C2")
        return false;
    c.variant = v;

    try
    {
        c.p   = static_cast<short_t>(std::stoi(tok[1]));
        c.dim = std::stoi(tok[2]);
    }
    catch (std::exception &)
    {
        return false;
    }
    if (c.p < 1 || (c.dim != 1 && c.dim != 2) || (c.dim == 2 && c.linear))
        return false;
    return true;
}

} // namespace

int main(int argc, char * argv[])
{
    index_t pmin = 1, pmax = 3, lmax = 3, N0 = 4;
    bool with2d = false, withC1 = true;
    std::string caseStr;

    gsCmdLine cmd("L2-projection convergence on periodic and non-periodic B-spline spaces.");
    cmd.addInt("", "pmin", "Smallest degree", pmin);
    cmd.addInt("", "pmax", "Largest degree", pmax);
    cmd.addInt("", "lmax", "Levels 0..lmax; N = N0*2^level", lmax);
    cmd.addInt("", "N0", "Number of elements at level 0", N0);
    cmd.addSwitch("with2d", "Also run the 2D cases", with2d);
    cmd.addSwitch("withC1", "Skip the C1 and C1-x variants (periodic constructor), which run by default", withC1);
    cmd.addString("", "case", "Single case '<variant>,<p>,<dim>', e.g. C2,2,1 "
                  "(variant: base, C1, C2, base-x, C1-x, C2-x)", caseStr);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    gsInfo << "variant    p  N    domain                 isPeriodic size  m  L2error       EoC     status\n"
           << std::flush;

    if (!caseStr.empty())
    {
        Case c;
        if (!parseCase(caseStr, c))
        {
            gsInfo << "Invalid --case '" << caseStr
                   << "'; expected <variant>,<p>,<dim> with variant in "
                      "{base,C1,C2,base-x,C1-x,C2-x} (-x only for dim 1), p >= 1, dim in {1,2}\n";
            return EXIT_FAILURE;
        }
        runCase(c, N0, lmax);
        return EXIT_SUCCESS;
    }

    const char * variants[] = {"base", "C1", "C2"};
    for (index_t dim = 1; dim <= (with2d ? 2 : 1); ++dim)
        for (int lin = 0; lin <= (dim == 1 ? 1 : 0); ++lin)
            for (const char * v : variants)
            {
                if (!withC1 && std::string(v) == "C1")
                    continue;
                for (index_t p = pmin; p <= pmax; ++p)
                {
                    Case c;
                    c.variant = v;
                    c.linear  = (lin == 1);
                    c.dim     = dim;
                    c.p       = static_cast<short_t>(p);
                    runCase(c, N0, lmax);
                }
                gsInfo << "\n";
            }

    return EXIT_SUCCESS;
}
