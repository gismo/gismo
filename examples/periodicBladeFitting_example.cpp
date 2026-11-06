/** @file periodicBladeFitting_example.cpp

    @brief Least-squares fit of a turbine-blade point cloud with a tensor
    B-spline surface that is periodic around the blade sections.

    The cloud is sampled from NACA four-digit sections, stacked along the span
    with linearly varying chord and twist (or read from a file with -f). It is
    fitted with gsFitting twice per control-net size:

      (P) a tensor B-spline basis that is periodic in u (around the section) and
          clamped in v (along the span). The seam, at the trailing edge (TE), has
          multiplicity mu: the space is C^{p-mu} across the seam, so the default
          mu = p gives a C^0 seam, which matches the corner of a closed NACA TE.
      (C) a tensor B-spline basis of the same size with a clamped u knot vector
          (no seam, the basis is interpolatory at the ends).

    Periodic space. With n = mu + k functions in u (k simple interior knots per
    period [0,1]), the knot vector is the stretched periodic one
        t_j = s_{(j-m) mod n} + L floor((j-m)/n),   m = p+1-mu,  L = 1,
    with s = [0 (mu times), 1/(k+1), ..., k/(k+1)]. Basis functions crossing the
    seam are identified, so the space has n = size() degrees of freedom, and the
    clamped basis with n_u functions has the same count.

    Synthetic blade. The section is the NACA four-digit section with maximum
    camber M, its position P and thickness t (default NACA 2412): thickness
        y_t(x) = 5t [0.2969 sqrt(x) - 0.1260 x - 0.3516 x^2 + 0.2843 x^3 - c4 x^4],
    camber line y_c(x) (parabolic arcs, joined at x = P), and surfaces
        upper: x - y_t sin(theta), y_c + y_t cos(theta),
        lower: x + y_t sin(theta), y_c - y_t cos(theta),   theta = atan(y_c').
    The standard coefficient is c4 = 0.1015, which leaves a TE of finite
    thickness. Here c4 = 0.1036 is used: the five coefficients then sum to 0,
    y_t(1) = 0 and the section closes at the TE (1,0), where the seam sits.
    Reference: E. N. Jacobs, K. E. Ward, R. M. Pinkerton, "The characteristics of
    78 related airfoil sections from tests in the variable-density wind tunnel",
    NACA Report No. 460 (1933); I. H. Abbott, A. E. von Doenhoff, "Theory of
    Wing Sections", Dover (1959), chap. 6.
    Points are cosine-spaced, ordered from the TE along the upper surface to the
    leading edge, then along the lower surface back towards the TE. The section
    is scaled by the chord, rotated by the twist about the quarter-chord axis,
    and placed at z = span * v. The u parameters are the normalised chord
    lengths of the closed polygon, starting at the TE, so u in [0,1).

    Measures per fit (printed in one table, sizes refined together):
      max / rms : max and root-mean-square Euclidean distance of the cloud to the
                  surface at the (optionally corrected) parameters,
      gap       : max over v of |S(a,v) - S(b,v)|, the position jump over the seam,
      kink      : max over v of the angle (deg) between dS/du at a and at b-,
                  the tangent jump over the seam.
    The periodic space has gap = 0 up to round-off by construction. A C^0 seam
    (mu = p) can represent the TE corner, so its kink is close to the TE turning
    angle of the section (printed for the synthetic cloud). A C^{p-mu >= 1} seam
    has kink = 0 and larger errors near the TE. Both are expected. Max errors
    of both fits should decrease under the joint refinement of n_u and n_v; no
    fixed order is asserted.

    Input file format (-f): a G+Smo XML file with
      matrix id 0 : 2 x N (u,v) parameters, u in [0,1] (period [0,1], seam at 0=1),
      matrix id 1 : 3 x N points (x,y,z).
    The v domain is [min v, max v] of the data.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H. Verhelst
*/

#include <gismo.h>

#include <iomanip>
#include <sstream>

using namespace gismo;

namespace {

/// Periodic knot vector of degree p on [0,1] with seam multiplicity mu and k
/// simple interior knots per period; the resulting basis has mu + k functions.
gsKnotVector<> periodicKnots(index_t k, short_t p, index_t mu)
{
    const index_t n = mu + k, m = p + 1 - mu;
    std::vector<real_t> s(n, 0);
    for (index_t i = 1; i <= k; ++i)
        s[mu + i - 1] = real_t(i) / real_t(k + 1);
    std::vector<real_t> t(n + m + p + 1);
    for (index_t j = 0; j < static_cast<index_t>(t.size()); ++j)
    {
        index_t q = (j - m) / n, r = (j - m) - q * n; // floored division
        if (r < 0) { r += n; --q; }
        t[j] = (1 == q && r < mu) ? real_t(1) : s[r] + real_t(q);
    }
    gsKnotVector<> kv(t, p);
    kv.setPeriodic(true); // validates the knots, throws std::runtime_error if invalid
    return kv;
}

/// Point of the closed NACA four-digit section (chord 1) at abscissa x on the
/// upper (upper = true) or lower surface.
void nacaPoint(real_t x, bool upper, real_t M, real_t P, real_t t,
               real_t & X, real_t & Y)
{
    const real_t yt = 5 * t * (0.2969 * math::sqrt(x) - 0.1260 * x
                               - 0.3516 * x * x + 0.2843 * x * x * x
                               - 0.1036 * x * x * x * x);
    real_t yc = 0, dyc = 0;
    if (M != 0 && P > 0 && P < 1)
    {
        if (x < P)
        {
            yc  = M / (P * P) * (2 * P * x - x * x);
            dyc = 2 * M / (P * P) * (P - x);
        }
        else
        {
            yc  = M / ((1 - P) * (1 - P)) * ((1 - 2 * P) + 2 * P * x - x * x);
            dyc = 2 * M / ((1 - P) * (1 - P)) * (P - x);
        }
    }
    const real_t th = math::atan(dyc);
    const real_t sg = upper ? 1 : -1;
    X = x - sg * yt * math::sin(th);
    Y = yc + sg * yt * math::cos(th);
}

/// Angle in degrees between vectors a and b of size 3, stable for small angles.
real_t angleDeg(const real_t * a, const real_t * b)
{
    const real_t cx = a[1] * b[2] - a[2] * b[1],
                 cy = a[2] * b[0] - a[0] * b[2],
                 cz = a[0] * b[1] - a[1] * b[0];
    const real_t cr  = math::sqrt(cx * cx + cy * cy + cz * cz);
    const real_t dot = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    return math::atan2(cr, dot) * 180 / EIGEN_PI;
}

/// Max position jump and max tangent angle (deg) across the u-seam of a surface.
/// At u = b the left limit is evaluated, at u = a the right limit.
void seamMeasures(const gsGeometry<> & geo, real_t & gap, real_t & kink)
{
    const gsMatrix<> supp = geo.support();
    const real_t a = supp(0, 0), b = supp(0, 1);
    const index_t Ns = 101;
    gsMatrix<> A(2, Ns), B(2, Ns);
    for (index_t s = 0; s < Ns; ++s)
    {
        const real_t w = supp(1, 0) + (supp(1, 1) - supp(1, 0)) * s / (Ns - 1);
        A(0, s) = a; A(1, s) = w;
        B(0, s) = b; B(1, s) = w;
    }
    gsMatrix<> SA, SB;
    geo.eval_into(A, SA);
    geo.eval_into(B, SB);
    const gsMatrix<> JA = geo.jacobian(A), JB = geo.jacobian(B);
    gap = 0; kink = 0;
    for (index_t s = 0; s < Ns; ++s)
    {
        gap  = std::max(gap, (SA.col(s) - SB.col(s)).norm());
        kink = std::max(kink, angleDeg(JA.col(2 * s).data(), JB.col(2 * s).data()));
    }
}

/// Error statistics of one fit.
struct FitStats
{
    real_t maxErr, rmsErr, gap, kink, seconds;
};

/// Fits (uv,xyz) in the given basis. A fresh fitter is used for every call.
/// Returns false if gsFitting produced no result.
bool fitOnce(const gsMatrix<> & uv, const gsMatrix<> & xyz, gsBasis<> & basis,
             real_t lambda, index_t numPc, FitStats & st, const std::string & plotName)
{
    gsStopwatch sw;
    gsFitting<real_t> fit(uv, xyz, basis);
    fit.compute(lambda);
    if (fit.result() == NULL)
        return false;
    if (numPc > 0)
        fit.parameterCorrection(1e-8, numPc, 1e-6);
    fit.computeErrors();
    st.seconds = sw.stop();

    st.maxErr = fit.maxPointError();
    const std::vector<real_t> & e = fit.pointWiseErrors();
    real_t sq = 0;
    for (size_t i = 0; i < e.size(); ++i)
        sq += e[i] * e[i];
    st.rmsErr = math::sqrt(sq / e.size());

    seamMeasures(*fit.result(), st.gap, st.kink);

    if (!plotName.empty())
        gsWriteParaview(*fit.result(), plotName, 20000);
    return true;
}

} // namespace

int main(int argc, char * argv[])
{
    index_t degU = 3, degV = 3, seamMult = -1, numU = 8, numV = 4, levels = 4;
    index_t numPc = 0, numPtsU = 256, numSec = 65;
    real_t lambda = 0, thick = 0.12, camber = 0.02, camberPos = 0.4;
    real_t chordRoot = 1.0, chordTip = 0.5, twistRoot = 0, twistTip = 20, span = 3.0;
    std::string fn;
    bool plot = false;

    gsCmdLine cmd("Fit a turbine-blade point cloud with a tensor B-spline surface that is "
                  "periodic around the sections (seam at the trailing edge) and compare it "
                  "with a clamped fit of equal size. A user cloud is an XML file with matrix id 0 "
                  "(2 x N, (u,v) with u in [0,1]) and id 1 (3 x N points).");
    cmd.addInt("p", "degU", "degree in u (around the section)", degU);
    cmd.addInt("q", "degV", "degree in v (along the span)", degV);
    cmd.addInt("m", "seamMult", "seam multiplicity 1..p (-1: p, a C^0 seam)", seamMult);
    cmd.addInt("n", "numCoefsU", "number of coefficients in u at level 0", numU);
    cmd.addInt("k", "numCoefsV", "number of coefficients in v at level 0", numV);
    cmd.addInt("r", "levels", "number of control-net sizes (doubling)", levels);
    cmd.addReal("s", "lambda", "smoothing weight", lambda);
    cmd.addInt("c", "parcor", "steps of parameter correction", numPc);
    cmd.addInt("a", "numPtsU", "points per section (even)", numPtsU);
    cmd.addInt("b", "numSections", "number of sections along the span", numSec);
    cmd.addReal("t", "thickness", "NACA maximum thickness t", thick);
    cmd.addReal("", "camber", "NACA maximum camber M", camber);
    cmd.addReal("", "camberPos", "NACA position of maximum camber P", camberPos);
    cmd.addReal("", "chordRoot", "chord at the root", chordRoot);
    cmd.addReal("", "chordTip", "chord at the tip", chordTip);
    cmd.addReal("", "twistRoot", "twist at the root (deg)", twistRoot);
    cmd.addReal("", "twistTip", "twist at the tip (deg)", twistTip);
    cmd.addReal("", "span", "span length", span);
    cmd.addString("f", "file", "user point cloud (XML, see description)", fn);
    cmd.addSwitch("plot", "write ParaView files and the synthetic cloud", plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // Option validation
    if (seamMult < 0) seamMult = degU;
    if (degU < 1 || degV < 1)
    { gsInfo << "Degrees must be at least 1.\n"; return EXIT_FAILURE; }
    if (seamMult < 1 || seamMult > degU)
    { gsInfo << "The seam multiplicity must satisfy 1 <= mu <= p = " << degU
             << " (got " << seamMult << ").\n"; return EXIT_FAILURE; }
    if (numU < degU + 1 || numV < degV + 1)
    { gsInfo << "Need numCoefsU >= p+1 and numCoefsV >= q+1.\n"; return EXIT_FAILURE; }
    if (levels < 1)
    { gsInfo << "Need at least one level.\n"; return EXIT_FAILURE; }
    if (numPc < 0 || lambda < 0)
    { gsInfo << "parcor and lambda must be non-negative.\n"; return EXIT_FAILURE; }
    if (fn.empty() && (numPtsU < 8 || numPtsU % 2 != 0 || numSec < 2))
    { gsInfo << "Need an even number of points per section >= 8 and >= 2 sections.\n";
      return EXIT_FAILURE; }

    // Point cloud
    gsStopwatch sw;
    gsMatrix<> uv, xyz;
    real_t v0 = 0, v1 = 1;
    std::string source;
    if (fn.empty())
    {
        const index_t Nu = numPtsU, Nv = numSec;
        uv.resize(2, Nu * Nv);
        xyz.resize(3, Nu * Nv);
        std::vector<real_t> Qx(Nu), Qy(Nu);
        for (index_t j = 0; j < Nv; ++j)
        {
            const real_t v = real_t(j) / real_t(Nv - 1);
            const real_t chord = chordRoot + (chordTip - chordRoot) * v;
            const real_t phi = (twistRoot + (twistTip - twistRoot) * v) * EIGEN_PI / 180;
            for (index_t i = 0; i < Nu; ++i)
            {
                const real_t beta = 2 * EIGEN_PI * i / Nu;
                const real_t x = (1 + math::cos(beta)) / 2;
                real_t X, Y;
                nacaPoint(x, beta <= EIGEN_PI, camber, camberPos, thick, X, Y);
                X = chord * (X - 0.25);
                Y = chord * Y;
                Qx[i] = X * math::cos(phi) - Y * math::sin(phi);
                Qy[i] = X * math::sin(phi) + Y * math::cos(phi);
            }
            // Chord-length parameters of the closed polygon, starting at the TE
            std::vector<real_t> d(Nu);
            real_t D = 0;
            for (index_t i = 0; i < Nu; ++i)
            {
                const index_t i1 = (i + 1) % Nu;
                d[i] = math::sqrt((Qx[i1] - Qx[i]) * (Qx[i1] - Qx[i])
                                  + (Qy[i1] - Qy[i]) * (Qy[i1] - Qy[i]));
                D += d[i];
            }
            real_t acc = 0;
            for (index_t i = 0; i < Nu; ++i)
            {
                const index_t c = j * Nu + i;
                uv(0, c) = acc / D;
                uv(1, c) = v;
                xyz(0, c) = Qx[i];
                xyz(1, c) = Qy[i];
                xyz(2, c) = span * v;
                acc += d[i];
            }
        }
        std::ostringstream os;
        os << "synthetic NACA " << std::fixed << std::setprecision(0)
           << 100 * camber << 10 * camberPos << std::setw(2) << std::setfill('0')
           << 100 * thick;
        source = os.str();
    }
    else
    {
        gsFileData<> fd_in(fn);
        memory::unique_ptr<gsMatrix<> > uvPtr  = fd_in.getId<gsMatrix<> >(0);
        memory::unique_ptr<gsMatrix<> > xyzPtr = fd_in.getId<gsMatrix<> >(1);
        if (!uvPtr || !xyzPtr)
        { gsInfo << "Wrong input: the file must contain matrix id 0 (2 x N) and id 1 (3 x N).\n";
          return EXIT_FAILURE; }
        uv  = give(*uvPtr);
        xyz = give(*xyzPtr);
        if (uv.rows() != 2 || xyz.rows() != 3 || uv.cols() != xyz.cols() || uv.cols() == 0)
        { gsInfo << "Wrong input: expected a 2 x N and a 3 x N matrix.\n"; return EXIT_FAILURE; }
        if (uv.row(0).minCoeff() < 0 || uv.row(0).maxCoeff() > 1)
        { gsInfo << "Wrong input: u must lie in [0,1].\n"; return EXIT_FAILURE; }
        v0 = uv.row(1).minCoeff();
        v1 = uv.row(1).maxCoeff();
        if (v0 == v1)
        { gsInfo << "Wrong input: v is constant.\n"; return EXIT_FAILURE; }
        source = fn;
    }
    const double tCloud = sw.stop();

    gsInfo << "Periodic blade fitting\n"
           << "  degrees (u,v)          : " << degU << ", " << degV << "\n"
           << "  seam multiplicity      : " << seamMult << "\n"
           << "  number of points       : " << uv.cols() << "\n"
           << "  lambda                 : " << lambda << "\n"
           << "  parameter corrections  : " << numPc << "\n"
           << "  data                   : " << source << "\n";

    if (fn.empty())
        gsInfo << "  points x sections      : " << numPtsU << " x " << numSec << "\n";

    if (fn.empty())
    {
        // Turning angle at the TE of the analytic section: one-sided directions in x
        const real_t h = 1e-6;
        real_t Ux, Uy, U0x, U0y, Lx, Ly, L0x, L0y;
        nacaPoint(1 - h, true, camber, camberPos, thick, Ux, Uy);
        nacaPoint(1, true, camber, camberPos, thick, U0x, U0y);
        nacaPoint(1, false, camber, camberPos, thick, L0x, L0y);
        nacaPoint(1 - h, false, camber, camberPos, thick, Lx, Ly);
        const real_t t1[3] = { Ux - U0x, Uy - U0y, 0 }, t2[3] = { L0x - Lx, L0y - Ly, 0 };
        gsInfo << "Reference TE turning angle (deg): " << angleDeg(t1, t2) << "\n";
    }

    gsInfo << "#   nu   nv      per_max      per_rms      per_gap     per_kink"
              "      cla_max      cla_rms      cla_gap     cla_kink\n";

    std::vector<FitStats> sP(levels), sC(levels);
    for (index_t l = 0; l < levels; ++l)
    {
        const index_t nu = numU << l, nv = numV << l;
        gsKnotVector<> kvP = periodicKnots(nu - seamMult, static_cast<short_t>(degU), seamMult);
        gsKnotVector<> kvC(0.0, 1.0, nu - degU - 1, degU + 1, 1, static_cast<short_t>(degU));
        gsKnotVector<> kvV(v0, v1, nv - degV - 1, degV + 1, 1, static_cast<short_t>(degV));
        gsTensorBSplineBasis<2, real_t> bP(kvP, kvV), bC(kvC, kvV);
        GISMO_ENSURE(bP.isPeriodic(0) && !bC.isPeriodic(),
                     "Expected a periodic and a clamped u-basis.");
        GISMO_ENSURE(bP.size(0) == nu && bC.size(0) == nu,
                     "The two u-bases must have the same size.");

        const bool last = (l == levels - 1);
        if (!fitOnce(uv, xyz, bP, lambda, numPc, sP[l], last && plot ? "periodicBlade_periodic" : ""))
        { gsInfo << "gsFitting returned no result for the periodic basis.\n"; return EXIT_FAILURE; }
        if (!fitOnce(uv, xyz, bC, lambda, numPc, sC[l], last && plot ? "periodicBlade_clamped" : ""))
        { gsInfo << "gsFitting returned no result for the clamped basis.\n"; return EXIT_FAILURE; }

        std::ostringstream row;
        row << std::setw(4) << nu << " " << std::setw(4) << nv
            << std::scientific << std::setprecision(4)
            << std::setw(13) << sP[l].maxErr << std::setw(13) << sP[l].rmsErr
            << std::setw(13) << sP[l].gap    << std::setw(13) << sP[l].kink
            << std::setw(13) << sC[l].maxErr << std::setw(13) << sC[l].rmsErr
            << std::setw(13) << sC[l].gap    << std::setw(13) << sC[l].kink;
        gsInfo << "  " << row.str() << "\n";
    }

    gsInfo << "Timings [s]: cloud " << tCloud << "\n";
    for (index_t l = 0; l < levels; ++l)
        gsInfo << "  level " << l << ": periodic fit " << sP[l].seconds
               << ", clamped fit " << sC[l].seconds << "\n";

    if (plot)
    {
        gsWriteParaviewPoints(xyz, "periodicBlade_points");
        gsInfo << "Written: periodicBlade_points.vtp/.pvd, periodicBlade_periodic.pvd/.vts, "
                  "periodicBlade_clamped.pvd/.vts";
        if (fn.empty())
        {
            gsFileData<> fd;
            fd.setFloatPrecision(17);
            fd.add(uv, 0);
            fd.add(xyz, 1);
            fd.save("periodicBlade_cloud");
            gsInfo << ", periodicBlade_cloud.xml";
        }
        gsInfo << "\n";
    }
    return EXIT_SUCCESS;
}
