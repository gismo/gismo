/** @file fischerkolmogrov_prostate_example_immersed.cpp

    @brief Fisher-Kolmogorov tumour growth on an IMMERSED prostate geometry,
           integrated with the finite-cell (Algoim + octree) quadrature.

    This is fischerkolmogrov_prostate_example.cpp with exactly one thing
    changed: the way the weak form is INTEGRATED. The generalized-alpha time
    stepping, the Newton loop, the hyperbolic-tangent tumour seed, the
    biomarkers, the csv output and the xml/command line interface are the ones
    of the body-fitted example.

    Instead of a body-fitted spline patch read from the xml file, the domain is
    given implicitly as Omega = {phi < 0}, where phi is the signed distance
    function of a closed triangle mesh (a Wavefront .obj file) normalized into
    the background box [0,1]^3 -- built exactly as in
    immersed_poisson2_nitsche_geometry_octree_example.cpp and its adaptive
    variant. The background box carries a plain tensor-product B-spline basis;
    its cells are classified by gsImplicitTrimmedDomain<3> and the quadrature
    points of Omega are

      - interior cells : a plain tensor Gauss rule,
      - cut cells      : the adaptive Algoim octree rule (gsAlgoimAdaptiveRule,
                         options --quadDepth / --indicator / --indicatorTol),
      - outside cells  : dropped.

    Those points and weights ARE the quadrature rule of this example: every
    integral -- the constant matrix, the Newton residual, the Newton tangent
    and all biomarkers -- is a sum over them. Since the background map is the
    identity on [0,1]^3, meas(G) == 1 and the parametric derivatives of the
    basis are already the physical ones, so the Fisher-Kolmogorov weak form
    reads, per quadrature point,

        w * dN  +  D * grad(w).grad(N)  -  rho * w * N * (1-N).

    Boundary conditions. The immersed boundary {phi == 0} carries the natural
    zero-flux condition D grad(N).n = 0, which is what the model wants (no
    tumour cells leave the prostate) and which needs no surface term at all --
    so, unlike the immersed Poisson examples, NO Nitsche coupling is assembled.
    The boundary conditions of the xml file are read (for interface
    compatibility) but not used. Degrees of freedom whose support never meets
    {phi < 0} are kept nonsingular by the same tiny (1e-10) full-background
    regularization as the Poisson examples.

    mkdir -p 00_output_FK_prostate
    ./build/bin/fischerkolmogrov_prostate_example_immersed \
        -f pde/fisher_kolmogorov_cube_bvp.xml -m prostata_lorenzo.obj \
        -r 4 -N 10 -t 0.5 -D 1e-3 -p 1e-3 -R 0.005 -s 0.075 --Nth 0.15 \
        --quadDepth 2 -o 00_output_FK_prostate --plot --plotRes 100\
        --xc 0.6 --yc 0.3 --zc 0.6

    Author(s): L. Venta-Viñuela
*/
#include <gismo.h>
#include <gsAlgoim/gsAlgoimRule.h>
#include <gsAlgoim/gsAlgoimAdaptiveRule.h>
#include <gsDomain/gsMeshLevelSet.h>

#include <iomanip>
#include <map>
#include <vector>

using namespace gismo;


// =============================================================================
// FISHER-KOLMOGOROV MODEL PARAMETERS (user-settable, for literature calibration)
// =============================================================================
// The values below are only the fall-back defaults, used when the label is not
// found in the xml file (see main()). They are chosen for a domain of unit size;
// D and rho carry the units of the geometry, so for e.g. a prostate measured in
// mm and time in days they have to be re-calibrated ([mm^2/day] and [1/day]).
// NOTE: the immersed geometry is normalized into [0,1]^3, so these defaults are
// the right order of magnitude here.
template <class T>
struct fkParams
{
    T D     = 1.0e-3; // tumor cell diffusivity          [length^2/time]
    T rho   = 0.5;    // net proliferation rate          [1/time]
    T R     = 0.15;   // radius of the initial tumor     [length]
    T sigma = 0.05;   // interface width of the IC       [length]
    T Nth   = 0.15;   // threshold defining Omega_T      [-]
    // Centre of the initial tumor (physical coordinates).
    // NaN means: use the centre of the bounding box of the domain.
    T xc = std::numeric_limits<T>::quiet_NaN();
    T yc = std::numeric_limits<T>::quiet_NaN();
    T zc = std::numeric_limits<T>::quiet_NaN();
};

// =============================================================================
// FINITE-CELL (IMMERSED) DISCRETIZATION SETTINGS
// =============================================================================
// Everything the immersed quadrature needs; the values come from the xml file
// (assembler option list, id 5) and/or the command line, same precedence as the
// model parameters above.
struct fcmParams
{
    std::string mesh         = "prostata_lorenzo.obj"; // closed surface mesh
    real_t      fill         = 0.9;   // fraction of [0,1]^3 spanned by the mesh
    index_t     degElevate   = 1;     // background basis degree = 1 + this
    index_t     quadDepth    = 2;     // octree depth of the cut-cell quadrature
    std::string indicator    = "integralChange"; // uniform|fallback|integralChange
    real_t      indicatorTol = 1e-2;  // integralChange acceptance tolerance
    index_t     plotRes      = 64;    // sample points per direction of the N field
    bool        plotPoints   = false; // also dump N at the quadrature points
    bool        plotMesh     = true;  // also write the background grid
};

// Smooth spherical tumor seed:  N(x,0) = 0.5*(1 - tanh((||x - x_tumor|| - R)/sigma))
template <short_t dim, class T>
gsFunctionExpr<T> tumorSeed(const fkParams<T> & par, const gsMultiPatch<T> & mp)
{
    // Default seed location: centre of the bounding box of the domain
    // (the background box [0,1]^3, whose centre is also the centre of the
    //  normalized immersed geometry -- gsNormalizeToUnitBox centers the mesh)
    gsMatrix<T> bbox;
    mp.boundingBox(bbox);
    gsVector<T,3> c;
    c.setZero();
    for (short_t d = 0; d != dim; d++)
        c[d] = 0.5 * (bbox(d,0) + bbox(d,1));

    if (!math::isnan(par.xc))            c[0] = par.xc;
    if (!math::isnan(par.yc))            c[1] = par.yc;
    if (dim == 3 && !math::isnan(par.zc)) c[2] = par.zc;

    std::ostringstream expr;
    expr.precision(16);
    expr << "0.5*(1-tanh((sqrt((x-(" << c[0] << "))^2+(y-(" << c[1] << "))^2";
    if (dim == 3)
        expr << "+(z-(" << c[2] << "))^2";
    expr << ")-(" << par.R << "))/(" << par.sigma << ")))";

    gsInfo << "Initial condition: N(x,0) = " << expr.str() << "\n";
    return gsFunctionExpr<T>(expr.str(), dim);
}

// =============================================================================
// IMMERSED QUADRATURE BLOCK
// =============================================================================
// This replaces gsExprAssembler's integration elements. One block holds the
// quadrature of one background cell that meets Omega = {phi < 0}, together with
// everything the assembly needs at those points. The background grid is static
// during the whole simulation, so the blocks -- basis values and derivatives
// included -- are built once and reused by every Newton iteration of every time
// step; only the coefficient vectors change.
template <class T>
struct fkQuadBlock
{
    gsMatrix<T>          pts;   // dim x nq       nodes (parametric == physical)
    gsVector<T>          wts;   // nq             weights
    gsMatrix<index_t>    act;   // na x nq        global index of active function i at point q
    gsMatrix<index_t>    loc;   // na x nq        slot of act(i,q) inside \a glob
    std::vector<index_t> glob;  // nu             local -> global map of the block
    gsMatrix<T>          bv;    // na x nq        basis values
    gsMatrix<T>          bd;    // (na*dim) x nq  basis derivatives
};

/// Evaluate the basis at (\a pts, \a wts) and cache it in \a b, together with a
/// contiguous local numbering of the functions active anywhere in the block.
/// The local numbering is needed because a quadrature point sitting exactly on
/// a cell face can be attributed to the neighbouring knot span, so the active
/// set is not guaranteed to be identical for every point of the block.
template <class T>
void fkBuildBlock(const gsBasis<T> & basis,
                  const gsMatrix<T> & pts, const gsVector<T> & wts,
                  fkQuadBlock<T> & b)
{
    b.pts = pts;
    b.wts = wts;
    basis.active_into(pts, b.act);
    basis.eval_into  (pts, b.bv);
    basis.deriv_into (pts, b.bd);

    b.glob.clear();
    b.loc.resize(b.act.rows(), b.act.cols());
    std::map<index_t,index_t> slot;
    for (index_t q = 0; q != b.act.cols(); ++q)
        for (index_t i = 0; i != b.act.rows(); ++i)
        {
            const std::pair<std::map<index_t,index_t>::iterator,bool> ins =
                slot.insert(std::make_pair(b.act(i,q), (index_t)b.glob.size()));
            if (ins.second) b.glob.push_back(b.act(i,q));
            b.loc(i,q) = ins.first->second;
        }
}

// =============================================================================
// VOLUMETRIC FIELD OUTPUT
// =============================================================================
// Writes one VTK XML structured grid (.vts) holding BOTH the cell density N and
// the level set phi, sampled on a uniform grid of the background box. Carrying
// phi in the same file is what makes the output usable for an immersed result:
// a single ParaView "Clip -> Scalars: phi -> Value 0 -> Invert" shows N on the
// prostate only, instead of the whole background box. (gsWriteParaview has no
// overload that puts two fields in one file, hence the hand-written .vts.)
//
// N is MASKED to 0 outside {phi<0}. Outside the immersed domain the spline is
// not a solution of anything: those coefficients are only held in place by the
// tiny background regularization, so they drift to large arbitrary values
// (O(10) after a few steps) which would otherwise dominate the ParaView colour
// range and hide the tumour. Masking costs nothing -- the equation is only
// solved, and the biomarkers only integrated, on {phi<0}.
template <class T>
void fkWriteField(const std::string & fn,
                  const gsVector<index_t> & res,  // sample points per direction
                  const gsMatrix<T> & gridPts,    // 3 x np, x fastest then y then z
                  const gsMatrix<T> & Nvals,      // 1 x np
                  const gsMatrix<T> & phiVals)    // 1 x np
{
    const index_t np = gridPts.cols();
    std::ofstream f(fn.c_str());
    f << std::setprecision(7);
    f << "<?xml version=\"1.0\"?>\n"
      << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
      << "  <StructuredGrid WholeExtent=\"0 " << res[0]-1 << " 0 " << res[1]-1
      << " 0 " << res[2]-1 << "\">\n"
      << "    <Piece Extent=\"0 " << res[0]-1 << " 0 " << res[1]-1
      << " 0 " << res[2]-1 << "\">\n"
      << "      <PointData Scalars=\"N\">\n"
      << "        <DataArray type=\"Float32\" Name=\"N\" format=\"ascii\">\n";
    for (index_t i = 0; i != np; ++i)
        f << Nvals(0,i) << (0 == (i+1) % 12 ? "\n" : " ");
    f << "\n        </DataArray>\n"
      << "        <DataArray type=\"Float32\" Name=\"phi\" format=\"ascii\">\n";
    for (index_t i = 0; i != np; ++i)
        f << phiVals(0,i) << (0 == (i+1) % 12 ? "\n" : " ");
    f << "\n        </DataArray>\n"
      << "      </PointData>\n"
      << "      <Points>\n"
      << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (index_t i = 0; i != np; ++i)
        f << gridPts(0,i) << " " << gridPts(1,i) << " " << gridPts(2,i) << "\n";
    f << "        </DataArray>\n"
      << "      </Points>\n"
      << "    </Piece>\n"
      << "  </StructuredGrid>\n"
      << "</VTKFile>\n";
    f.close();
}

// =============================================================================
// SOLVER ROUTINE
// =============================================================================
template <short_t dim, class T>
void solve(gsMultiPatch<T> & mp,
           const gsFunction<T> & phi,
           gsBoundaryConditions<T> & bc,
           gsOptionList & CHopt,
           gsOptionList & TIMEopt,
           gsOptionList & Aopt,
           real_t & dt,
           index_t & maxSteps,
           index_t & numRefine,
           const fkParams<T> & par,
           const fcmParams & fcm,
           std::string out,
           bool plot)
{
    // The immersed boundary {phi==0} carries the natural zero-flux condition,
    // which contributes no term to the weak form: the xml boundary conditions
    // (which describe the boundary of the background box) are not used.
    GISMO_UNUSED(bc);

    // -------------------------------------------------------------------------
    // 1. Basis & Problem Setup (Uniform Tensor Product Basis, No Adaptivity)
    //    The patch is now the BACKGROUND box [0,1]^dim, not a body-fitted patch.
    // -------------------------------------------------------------------------
    gsMultiBasis<T> dbasis(mp, true);
    dbasis.setDegree(dbasis.maxCwiseDegree() + fcm.degElevate);
    for (index_t r = 0; r != numRefine; ++r)
        for (size_t p = 0; p != dbasis.nBases(); ++p)
        {
            dbasis.basis(p).uniformRefine();
        }

    gsTensorBSplineBasis<dim,T> * tbsPtr =
        dynamic_cast<gsTensorBSplineBasis<dim,T> *>(&dbasis.basis(0));
    GISMO_ENSURE(tbsPtr, "Expected a tensor B-spline background basis.");
    const index_t deg   = dbasis.maxCwiseDegree();
    const index_t nDofs = tbsPtr->size();

    real_t hmax    = dbasis.basis(0).getMaxCellLength();
    real_t lambda  = CHopt.askReal("lambda", 1.0 / (32.0 * pow(EIGEN_PI, 2)));
    real_t penalty = 1e4 * lambda;
    GISMO_UNUSED(hmax); GISMO_UNUSED(penalty);

    // -------------------------------------------------------------------------
    // 1b. IMMERSED QUADRATURE (this is the only part that differs from the
    //     body-fitted example). Classify the background cells and collect the
    //     quadrature of Omega = {phi < 0} once: interior cells get a tensor
    //     Gauss rule, cut cells the adaptive Algoim octree rule.
    // -------------------------------------------------------------------------
    const index_t nFallback = 4; // uniform samples/direction when Algoim returns 0 nodes

    gsImplicitTrimmedDomain<dim,T> tr_domain(phi, *tbsPtr);

    gsGaussRule<T> gauss(gsVector<index_t>::Constant(dim, deg + 1));
    gsOptionList qopt = gsAlgoimAdaptiveRule<T>::defaultOptions();
    qopt.setInt   ("maxDepth",     fcm.quadDepth);
    qopt.setInt   ("nFallback",    nFallback);
    qopt.setString("indicator",    fcm.indicator);
    qopt.setReal  ("indicatorTol", fcm.indicatorTol);
    gsAlgoimAdaptiveRule<T> volumeRule(phi, *tbsPtr, qopt);

    std::vector<fkQuadBlock<T> > blocks; // quadrature of Omega, one entry per cell
    gsMatrix<T> qpts;
    gsVector<T> qwts;
    index_t nQuadInterior = 0, nQuadCut = 0, nInteriorCells = 0, nCutCells = 0;

    gsStopwatch swQuad;
    for (auto it = tr_domain.beginInterior(); it != tr_domain.template end<InteriorSign>(); ++it)
    {
        gauss.mapTo(it.lowerCorner(), it.upperCorner(), qpts, qwts);
        ++nInteriorCells;
        if (qpts.cols() == 0) continue;
        blocks.push_back(fkQuadBlock<T>());
        fkBuildBlock(*tbsPtr, qpts, qwts, blocks.back());
        nQuadInterior += qpts.cols();
    }
    for (auto it = tr_domain.beginBdr(boundary::none); it != tr_domain.endBdr(boundary::none); ++it)
    {
        volumeRule.mapTo(it.lowerCorner(), it.upperCorner(), qpts, qwts);
        ++nCutCells;
        if (qpts.cols() == 0) continue;
        blocks.push_back(fkQuadBlock<T>());
        fkBuildBlock(*tbsPtr, qpts, qwts, blocks.back());
        nQuadCut += qpts.cols();
    }
    const double tQuad = swQuad.stop();

    // Measure of the immersed domain: the plain sum of the weights. This is the
    // volume the biomarkers below are normalized with, so it is worth printing.
    T volOmega = 0;
    for (size_t e = 0; e != blocks.size(); ++e)
        volOmega += blocks[e].wts.sum();

    // The blocks cache the basis values/derivatives at every quadrature point,
    // which is what makes the repeated Newton assemblies cheap -- but it is also
    // what the memory of a run is spent on, so report it.
    size_t cacheBytes = 0;
    for (size_t e = 0; e != blocks.size(); ++e)
        cacheBytes += sizeof(T)*(blocks[e].bv.size() + blocks[e].bd.size()
                                 + blocks[e].pts.size() + blocks[e].wts.size())
                    + sizeof(index_t)*(blocks[e].act.size() + blocks[e].loc.size());

    gsInfo << "Immersed quadrature: " << nInteriorCells << " interior cells ("
           << nQuadInterior << " pts) + " << nCutCells << " cut cells ("
           << nQuadCut << " pts), total " << (nQuadInterior+nQuadCut)
           << " points, built in " << tQuad << " s ("
           << (cacheBytes >> 20) << " MB cached).\n"
           << "  quadDepth = " << fcm.quadDepth << ", indicator = " << fcm.indicator
           << ", |Omega| = " << volOmega << " (of the unit box)\n";

    // -------------------------------------------------------------------------
    // 2. Time Integration Setup (Generalized-Alpha Method)
    // -------------------------------------------------------------------------
    real_t rho_inf = TIMEopt.askReal("rho_inf", 0.5);
    real_t alpha_m = 0.5 * (3.0 - rho_inf) / (1.0 + rho_inf);
    real_t alpha_f = 1.0 / (1.0 + rho_inf);
    real_t gamma   = 0.5 + alpha_m - alpha_f;

    // Parameters (fixed model constants, set from the command line)
    real_t D   = par.D;   // tumor cell diffusivity
    real_t rho = par.rho; // net proliferation rate
    gsInfo << "Fisher-Kolmogorov parameters: D = " << D << ", rho = " << rho
           << ", R = " << par.R << ", sigma = " << par.sigma << ", N_th = " << par.Nth << "\n";
    gsInfo << "Discretization: numRefine = " << numRefine << ", degree = " << deg
           << ", DoFs = " << nDofs << ", dt = " << dt << ", steps = " << maxSteps << "\n";

    // Initialization biomarkers
    real_t V_T, V_N, V_P, N_dash, N_T, A_p;
    real_t int_N, int_H; // \int_Omega N dOmega  and  \int_Omega H(N-N_th) dOmega

    gsFileManager::mkdir(out);
    std::ofstream csvFile;
    csvFile.open(out+"/output_data.csv");
    csvFile << "TimeStep,Time,NumDOFs,int_N,int_H,V_P,N_dash,N_T,A_p\n";

    real_t N_hat_th = par.Nth;  // Cell density threshold defining Omega_T (tumor)
    real_t k_steep = 200.0;  // Steepness factor for Heaviside step function

    // -------------------------------------------------------------------------
    // 3. Assembly on the immersed quadrature points
    //    (this replaces the gsExprAssembler expressions of the body-fitted
    //     example; the background map is the identity on [0,1]^dim, so
    //     meas(G) == 1 and the parametric derivatives are the physical ones)
    // -------------------------------------------------------------------------

    // Scatter a local dense block (indexed by fkQuadBlock::glob) into triplets.
    auto scatter = [](const fkQuadBlock<T> & b, const gsMatrix<T> & Ke,
                      gsSparseEntries<T> & trip)
    {
        const index_t nu = (index_t)b.glob.size();
        for (index_t i = 0; i != nu; ++i)
            for (index_t j = 0; j != nu; ++j)
                if ((T)0 != Ke(i,j))
                    trip.add(b.glob[i], b.glob[j], Ke(i,j));
    };

    // Constant part of the tangent, scaled by \a coeff:
    //    alpha_m * w*w^T  +  alpha_f*gamma*dt * D * grad(w).grad(w)^T
    // \a coeff is 1 on Omega and the tiny regularization weight on the
    // background cells.
    auto addConstant = [&](const fkQuadBlock<T> & b, T coeff, gsSparseEntries<T> & trip)
    {
        const index_t na = b.act.rows(), nu = (index_t)b.glob.size();
        gsMatrix<T> Ke(nu,nu);
        Ke.setZero();
        for (index_t q = 0; q != b.pts.cols(); ++q)
        {
            const T w = b.wts[q] * coeff;
            for (index_t i = 0; i != na; ++i)
            {
                const index_t li = b.loc(i,q);
                for (index_t j = 0; j != na; ++j)
                {
                    T gg = 0;
                    for (short_t d = 0; d != dim; ++d)
                        gg += b.bd(i*dim+d,q) * b.bd(j*dim+d,q);
                    Ke(li, b.loc(j,q)) += w * ( alpha_m * b.bv(i,q) * b.bv(j,q)
                                              + (alpha_f*gamma*dt) * D * gg );
                }
            }
        }
        scatter(b, Ke, trip);
    };

    // Domain residual of  dN/dt = div(D grad N) + rho*N*(1-N):
    //   \int w*dN + \int D*grad(w).grad(N) - \int rho*w*N*(1-N) = 0,
    // integrated over the immersed quadrature points of Omega. When \a trip is
    // non-null, the reaction part of the tangent,
    //   (-alpha_f*gamma*dt) * rho * w * (1 - 2N) * w^T,
    // is accumulated in the same sweep (it needs the very same N at the very
    // same points).
    auto assembleFK = [&](const gsMatrix<T> & Na, const gsMatrix<T> & dNa,
                          gsMatrix<T> & Q, gsSparseEntries<T> * trip)
    {
        Q.setZero(nDofs, 1);
        gsVector<T> gN(dim);
        for (size_t e = 0; e != blocks.size(); ++e)
        {
            const fkQuadBlock<T> & b = blocks[e];
            const index_t na = b.act.rows(), nu = (index_t)b.glob.size();
            gsVector<T> Fe(nu);
            Fe.setZero();
            gsMatrix<T> Ke;
            if (trip) Ke.setZero(nu,nu);

            for (index_t q = 0; q != b.pts.cols(); ++q)
            {
                // N, grad(N) and dN/dt at the quadrature point
                T Nv = 0, dNv = 0;
                gN.setZero();
                for (index_t i = 0; i != na; ++i)
                {
                    const index_t gi = b.act(i,q);
                    Nv  += b.bv(i,q) * Na (gi,0);
                    dNv += b.bv(i,q) * dNa(gi,0);
                    for (short_t d = 0; d != dim; ++d)
                        gN[d] += b.bd(i*dim+d,q) * Na(gi,0);
                }

                const T w = b.wts[q];
                for (index_t i = 0; i != na; ++i)
                {
                    T gg = 0;
                    for (short_t d = 0; d != dim; ++d)
                        gg += b.bd(i*dim+d,q) * gN[d];
                    Fe[b.loc(i,q)] += w * ( b.bv(i,q) * dNv + D * gg
                                          - rho * b.bv(i,q) * Nv * (1.0 - Nv) );
                }

                if (trip)
                {
                    const T cf = w * (-alpha_f*gamma*dt) * rho * (1.0 - 2.0*Nv);
                    for (index_t i = 0; i != na; ++i)
                    {
                        const index_t li = b.loc(i,q);
                        for (index_t j = 0; j != na; ++j)
                            Ke(li, b.loc(j,q)) += cf * b.bv(i,q) * b.bv(j,q);
                    }
                }
            }

            for (index_t i = 0; i != nu; ++i)
                Q(b.glob[i],0) += Fe[i];
            if (trip) scatter(b, Ke, *trip);
        }
    };

    // Quantities of interest, all of them integrals over Omega of a function of
    // N alone -- this is what computeInt() did with the expression evaluator on
    // the body-fitted patch, now evaluated on the immersed quadrature points.
    auto biomarkers = [&](const gsMatrix<T> & coefs,
                          real_t & o_int_N, real_t & o_int_H, real_t & o_V_N,
                          real_t & o_V_P, real_t & o_A_p)
    {
        o_int_N = o_int_H = o_V_N = o_V_P = o_A_p = 0;
        for (size_t e = 0; e != blocks.size(); ++e)
        {
            const fkQuadBlock<T> & b = blocks[e];
            const index_t na = b.act.rows();
            for (index_t q = 0; q != b.pts.cols(); ++q)
            {
                T Nv = 0;
                for (index_t i = 0; i != na; ++i)
                    Nv += b.bv(i,q) * coefs(b.act(i,q),0);
                // Indicator of the tumor domain Omega_T (where N >= N_th)
                const T H = 1.0 / (1.0 + math::exp(-k_steep*(Nv - N_hat_th)));
                const T w = b.wts[q];
                o_int_N += w * Nv;              // \int_Omega N
                o_int_H += w * H;               // \int_Omega H(N-N_th) = V_T
                o_V_N   += w * Nv * H;          // total tumor cell volume
                o_V_P   += w;                   // prostate volume
                o_A_p   += w * rho * Nv * (1.0 - Nv) * H; // proliferation activity
            }
        }
    };

    // -------------------------------------------------------------------------
    // 4. Initial Condition Setup (smooth hyperbolic-tangent tumor seed)
    // -------------------------------------------------------------------------
    gsFunctionExpr<T> N0 = tumorSeed<dim,T>(par, mp);

    // L2-projection of the (physical-space) seed onto the background basis.
    // There are no eliminated degrees of freedom here (the immersed boundary
    // condition is natural), so the coefficients ARE the unknowns.
    gsMatrix<T> N0coefs;
    gsInfo << "IC L2-projection error: "
           << gsL2Projection<T>::project(dbasis, mp, N0, N0coefs) << "\n";

    gsMatrix<T> Nnew, dNnew, Nold, dNold;
    gsMatrix<T> Nalpha, dNalpha, dNupdate, Q;

    Nold = N0coefs;
    dNold.setZero(nDofs, 1);

    typename gsSparseSolver<T>::uPtr solver = gsSparseSolver<T>::get("SimplicialLDLT");
    gsSparseMatrix<T> K, K_const;

    // -------------------------------------------------------------------------
    // 5. Pre-assemble Constant Terms Matrix
    // -------------------------------------------------------------------------
    {
        gsSparseEntries<T> trip;
        for (size_t e = 0; e != blocks.size(); ++e)
            addConstant(blocks[e], (T)1, trip);

        // Tiny full-background regularization: keeps the tangent nonsingular for
        // the degrees of freedom whose support never touches {phi<0}. Exactly
        // as in the immersed Poisson examples, it is added to the MATRIX only,
        // so it contributes nothing to the residual and those dofs simply stay
        // at their initial value.
        fkQuadBlock<T> bg;
        for (auto it = tbsPtr->domain()->beginAll(); it != tbsPtr->domain()->endAll(); ++it)
        {
            gauss.mapTo(it.lowerCorner(), it.upperCorner(), qpts, qwts);
            fkBuildBlock(*tbsPtr, qpts, qwts, bg);
            addConstant(bg, (T)1e-10, trip);
        }

        K_const.resize(nDofs, nDofs);
        K_const.setFrom(trip);
    }

    // -------------------------------------------------------------------------
    // 6. ParaView output of the cell density N (one time step per frame).
    //    Primary output: N as a volumetric FIELD on a uniform sample grid of the
    //    background box, written together with the level set phi in the same
    //    .vts file, so that one ParaView "Clip -> Scalars: phi -> Value 0"
    //    restricts it to the prostate. Optionally (--plotPoints) the same N
    //    sampled at the immersed quadrature points, which is the honest picture
    //    of where the weak form is actually integrated.
    // -------------------------------------------------------------------------
    std::unique_ptr<gsParaviewCollection> colField, colPts, colMesh;
    gsMatrix<T> cloud, gridPts, phiGrid, Nvals;
    gsVector<index_t> res(3);

    if (plot)
    {
        gsFileManager::mkdir(out + "/N");
        colField.reset(new gsParaviewCollection(out + "/N/N"));

        // Uniform sample grid of [0,1]^3, x fastest then y then z (VTK order).
        res.setConstant(math::max(fcm.plotRes, (index_t)2));
        gridPts.resize(3, res.prod());
        index_t c = 0;
        for (index_t k = 0; k != res[2]; ++k)
        for (index_t j = 0; j != res[1]; ++j)
        for (index_t i = 0; i != res[0]; ++i, ++c)
        {
            gridPts(0,c) = (T)i / (T)(res[0]-1);
            gridPts(1,c) = (T)j / (T)(res[1]-1);
            gridPts(2,c) = (T)k / (T)(res[2]-1);
        }
        // phi is static, so it is evaluated once and reused by every frame
        // (each evaluation is a BVH query against the triangle mesh).
        gsStopwatch swPhi;
        phi.eval_into(gridPts, phiGrid);
        gsInfo << "Plot grid: " << res[0] << "^3 = " << res.prod()
               << " points, level set evaluated in " << swPhi.stop() << " s.\n";
    }
    if (plot && fcm.plotMesh)
    {
        // The background grid the immersed quadrature lives on. The Poisson
        // examples rewrite it at every refinement level; here it is static for
        // the whole simulation, so it is written once as the single part of
        // its own collection.
        gsFileManager::mkdir(out + "/background");
        colMesh.reset(new gsParaviewCollection(out + "/background/background"));
        gsMesh<T> bgMesh(dbasis.basis(0));
        gsWriteParaview(bgMesh, out + "/background/grid");
        colMesh->addPart("grid.vtp", 0);
    }
    if (plot && fcm.plotPoints)
    {
        cloud.resize(dim+1, nQuadInterior + nQuadCut);
        gsFileManager::mkdir(out + "/Npts");
        colPts.reset(new gsParaviewCollection(out + "/Npts/Npts"));
    }

    // N sampled on the background grid -> one .vts (N and phi) per time step.
    auto writeField = [&](const gsMatrix<T> & coefs, index_t frame, real_t time)
    {
        typename gsGeometry<T>::uPtr geo = dbasis.basis(0).makeGeometry(coefs);
        geo->eval_into(gridPts, Nvals);
        for (index_t i = 0; i != Nvals.cols(); ++i) // see fkWriteField
            if (phiGrid(0,i) >= 0) Nvals(0,i) = 0;
        const std::string fn = "N_" + util::to_string(frame);
        fkWriteField<T>(out + "/N/" + fn + ".vts", res, gridPts, Nvals, phiGrid);
        colField->addPart(fn + ".vts", time);
    };

    // N sampled at the immersed quadrature points (4 x N: xyz + value).
    auto writeCloud = [&](const gsMatrix<T> & coefs, index_t frame, real_t time)
    {
        index_t c = 0;
        for (size_t e = 0; e != blocks.size(); ++e)
        {
            const fkQuadBlock<T> & b = blocks[e];
            const index_t na = b.act.rows();
            for (index_t q = 0; q != b.pts.cols(); ++q, ++c)
            {
                T Nv = 0;
                for (index_t i = 0; i != na; ++i)
                    Nv += b.bv(i,q) * coefs(b.act(i,q),0);
                cloud.col(c).head(dim) = b.pts.col(q);
                cloud(dim,c) = Nv;
            }
        }
        const std::string fn = "Npts_" + util::to_string(frame);
        gsWriteParaviewPoints(cloud, out + "/Npts/" + fn);
        colPts->addPart(fn + ".vtp", time);
    };

    // frame i of the collection = row i-1 of the csv file
    auto writeOut = [&](const gsMatrix<T> & coefs, index_t frame, real_t time)
    {
        if (!plot) return;
        writeField(coefs, frame, time);
        if (fcm.plotPoints) writeCloud(coefs, frame, time);
    };

    writeOut(Nold, 0, 0.0); // initial condition

    // -------------------------------------------------------------------------
    // 7. Time Integration & Newton-Raphson Solving Loop
    // -------------------------------------------------------------------------
    index_t maxIt = 50;
    real_t tol    = TIMEopt.askReal("tol", 1e-4);

    for (index_t step = 0; step < maxSteps; ++step)
    {
        gsInfo << "Time Step " << step + 1 << " / " << maxSteps << "\n";

        Nnew  = Nold;
        dNnew = (gamma - 1.0) / gamma * dNold;

        for (index_t it = 0; it < maxIt; ++it)
        {
            Nalpha  = Nold  + alpha_f * (Nnew  - Nold);
            dNalpha = dNold + alpha_m * (dNnew - dNold);

            // Domain Residual + reaction Tangent Assembly, on the immersed
            // quadrature points (both come out of the same sweep because they
            // need N at the same points).
            gsSparseEntries<T> trip;
            assembleFK(Nalpha, dNalpha, Q, &trip);

            if (it > 0 && Q.norm() < tol)
            {
                gsInfo << "  Converged in " << it << " iterations.\n";
                break;
            }

            K.resize(nDofs, nDofs);
            K.setFrom(trip);
            K += K_const;

            // Solve and Update
            solver->compute(K);
            dNupdate = solver->solve(-Q);

            dNnew += dNupdate;
            Nnew.noalias() += (gamma * dt) * dNupdate;
        }

        // Advance step
        Nold  = Nnew;
        dNold = dNnew;

        writeOut(Nold, step+1, (step+1)*dt);

        // ==============================================================
        // Quantities of interest of the forward simulation, integrated over the
        // immersed domain with the very same quadrature as the weak form.
        biomarkers(Nold, int_N, int_H, V_N, V_P, A_p);
        gsInfo << "  int_Omega N dOmega          = " << int_N << "\n"
               << "  int_Omega H(N-N_th) dOmega  = " << int_H << "\n";

        V_T    = int_H;         // Tumor volume
        N_dash = V_N / V_T;     // mean normalized tumor cell density
        N_T    = N_dash * V_T / V_P; // total tumor index
        // ==============================================================

        // Write in a csv file
        csvFile << step << "," << (step+1)*dt << "," << nDofs << ","
                << int_N << "," << int_H << ","
                << V_P << "," << N_dash << "," << N_T << "," << A_p << "\n";
        csvFile.flush();
    }
    csvFile.close();
    if (plot)
    {
        colField->save();
        if (fcm.plotPoints) colPts->save();
        if (fcm.plotMesh) colMesh->save();
    }

}

// =============================================================================
// MAIN FUNCTION
// =============================================================================
int main(int argc, char *argv[])
{
    std::string fn("pde/fisher_kolmogorov_cube_bvp.xml");
    std::string out("output");

    // Pre-scan the command line for the setup file, so that the option lists it
    // contains can serve as the defaults of all the other arguments below.
    // Precedence: built-in default  <  xml file  <  command line.
    for (int i = 1; i+1 != argc; ++i)
    {
        const std::string arg(argv[i]);
        if (arg == "-f" || arg == "--file") fn = argv[i+1];
    }

    // Read input problem settings
    gsFileData<> fd(fn);
    gsInfo << "Loaded configuration file: " << fd.lastPath() << "\n";

    gsBoundaryConditions<> bc;
    fd.getId(2, bc); // Boundary conditions (not used on the immersed domain)

    gsOptionList CHopt, TIMEopt, Aopt;
    fd.getId(3, CHopt);   // Physics options
    fd.getId(4, TIMEopt); // Integrator options
    fd.getId(5, Aopt);    // Assembler options

    // Model parameters (xml id 3), the labels match the long command line names
    fkParams<real_t> par;
    par.D     = CHopt.askReal("diffusivity", par.D    );
    par.rho   = CHopt.askReal("rho",         par.rho  );
    par.R     = CHopt.askReal("radius",      par.R    );
    par.sigma = CHopt.askReal("sigma",       par.sigma);
    par.Nth   = CHopt.askReal("Nth",         par.Nth  );
    par.xc    = CHopt.askReal("xc",          par.xc   );
    par.yc    = CHopt.askReal("yc",          par.yc   );
    par.zc    = CHopt.askReal("zc",          par.zc   );

    // Time integration (xml id 4) and discretization (xml id 5)
    real_t  dt        = TIMEopt.askReal  ("dt",        0.5);
    index_t maxSteps  = TIMEopt.askInt   ("Nsteps",    10 );
    index_t numRefine = Aopt.askInt      ("numRefine", 4  );
    bool    plot      = Aopt.askSwitch   ("plot",   false );

    // Immersed (finite-cell) settings, same precedence as everything else
    fcmParams fcm;
    fcm.mesh         = Aopt.askString("mesh",         fcm.mesh        );
    fcm.fill         = Aopt.askReal  ("fill",         fcm.fill        );
    fcm.degElevate   = Aopt.askInt   ("degreeElevation", fcm.degElevate);
    fcm.quadDepth    = Aopt.askInt   ("quadDepth",    fcm.quadDepth   );
    fcm.indicator    = Aopt.askString("indicator",    fcm.indicator   );
    fcm.indicatorTol = Aopt.askReal  ("indicatorTol", fcm.indicatorTol);
    fcm.plotRes      = Aopt.askInt   ("plotRes",      fcm.plotRes     );
    fcm.plotPoints   = Aopt.askSwitch("plotPoints",   fcm.plotPoints  );
    fcm.plotMesh     = Aopt.askSwitch("plotMesh",     fcm.plotMesh    );
    bool noPlotMesh  = !fcm.plotMesh;

    // Command-line interface
    gsCmdLine cmd("Fisher-Kolmogorov Generalized-Alpha Solver on an IMMERSED domain (forward problem)");
    cmd.addReal("t", "dt", "Time step size", dt);
    cmd.addInt("N", "Nsteps", "Number of time steps", maxSteps);
    cmd.addInt("r", "numRefine", "Number of uniform refinement steps", numRefine);
    cmd.addString("f", "file", "Input XML setup file", fn);
    cmd.addString( "o", "output", "Output directory", out);
    cmd.addSwitch("plot", "Write ParaView files of N for every time step", plot);
    cmd.addReal("D", "diffusivity", "Tumor cell diffusivity D [length^2/time]", par.D);
    cmd.addReal("p", "rho", "Net proliferation rate rho [1/time]", par.rho);
    cmd.addReal("R", "radius", "Radius R of the initial tumor [length]", par.R);
    cmd.addReal("s", "sigma", "Interface width sigma of the initial tumor [length]", par.sigma);
    cmd.addReal("", "Nth", "Cell density threshold N_th defining Omega_T", par.Nth);
    cmd.addReal("", "xc", "x-coordinate of the tumor centre (default: domain centre)", par.xc);
    cmd.addReal("", "yc", "y-coordinate of the tumor centre (default: domain centre)", par.yc);
    cmd.addReal("", "zc", "z-coordinate of the tumor centre (default: domain centre)", par.zc);
    // Immersed / finite-cell quadrature
    cmd.addString("m", "mesh", "Closed surface mesh (.obj) defining the immersed domain", fcm.mesh);
    cmd.addReal("", "fill", "Fill fraction of the background box [0,1]^3 (0..1)", fcm.fill);
    cmd.addInt("e", "degreeElevation", "Degree elevation of the background basis", fcm.degElevate);
    cmd.addInt("", "quadDepth", "Octree subdivision depth of the cut-cell quadrature", fcm.quadDepth);
    cmd.addString("", "indicator", "Adaptive indicator: uniform, fallback or integralChange", fcm.indicator);
    cmd.addReal("", "indicatorTol", "integralChange acceptance tolerance", fcm.indicatorTol);
    cmd.addInt("", "plotRes", "Sample points per direction of the plotted N field", fcm.plotRes);
    cmd.addSwitch("plotPoints", "Additionally dump N at the immersed quadrature points", fcm.plotPoints);
    cmd.addSwitch("noPlotMesh", "Do not write the background grid", noPlotMesh);

    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }
    fcm.plotMesh = !noPlotMesh;

    // -------------------------------------------------------------------------
    // Immersed geometry: load the closed triangle mesh, normalize it into the
    // background box [0,1]^3 and build its signed distance level set -- exactly
    // as in the immersed Poisson examples.
    // -------------------------------------------------------------------------
    gsSurfMesh mesh;
    if (!gsReadSurfMesh(fcm.mesh, mesh))
    {
        gsWarn << "Failed to read a triangle mesh from: " << fcm.mesh << "\n";
        return EXIT_FAILURE;
    }
    gsInfo << "Loaded mesh '" << fcm.mesh << "': " << mesh.n_vertices()
           << " vertices, " << mesh.n_faces() << " triangles.\n";

    const real_t scale = gsNormalizeToUnitBox(mesh, fcm.fill);
    GISMO_UNUSED(scale);

    gsMatrix<real_t> bbox(3, 2);
    bbox.col(0).setZero();
    bbox.col(1).setOnes();
    gsMeshSignedDist<real_t> impl_fun(mesh, bbox);

    // Background box [0,1]^3 (identity map: parameter space == physical space).
    // It replaces the body-fitted patch of the xml file, which is ignored.
    if (fd.hasId(0))
        gsInfo << "Note: the geometry of the xml file is ignored; the domain is "
                  "the immersed mesh inside the background box [0,1]^3.\n";
    gsMultiPatch<> mp(*gsNurbsCreator<>::BSplineCube((real_t)1));
    bc.setGeoMap(mp);

    solve<3, real_t>(mp, impl_fun, bc, CHopt, TIMEopt, Aopt,
                     dt, maxSteps, numRefine, par, fcm, out, plot);

    if (plot)
        gsWriteParaview(mesh, out + "/" + gsFileManager::getBasename(fcm.mesh) + "_geometry");

    return EXIT_SUCCESS;
}
