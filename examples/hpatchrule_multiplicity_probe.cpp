/** @file hpatchrule_multiplicity_probe.cpp

    @brief Throwaway diagnostic probe: does a THB fixture with a repeated
    interior knot actually change the gsHPatchRule target-knot-vector
    construction?

    This is NOT a test and NOT a deliverable feature. It reproduces, side by
    side, the CURRENT target-knot-vector formula used by gsHPatchRule::_build
    (uniform end/interior multiplicity from the unique breakpoints) and a NEW
    formula that instead reads the real multiplicities of the level-l knot
    vector, and reports whether any candidate fixture makes the two formulas
    disagree on a leaf that is actually reachable by the hierarchical domain.
    The two formulas are reproduced locally because the members they touch
    (_build, _ensureEven1D, CellData, ...) are private to gsHPatchRule.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>
#include <gsAssembler/gsHPatchRule.h>
#include <gsDomain/gsHDomain.h>
#include <gsNurbs/gsBSplineBasis.h>

using namespace gismo;

// Mirrors the private gsHPatchRule<T>::_ensureEven1D (gsHPatchRule.hpp:506-535),
// copied verbatim because the probe is an external example and has no access
// to it. Kept exactly as written -- including the diff array built over the
// *repeated* knots -- because the probe's job is to predict what the real
// code does, not what it should do.
template<class T>
gsKnotVector<T> ensureEven1D(const gsKnotVector<T> & kv)
{
    gsBSplineBasis<T> basis(kv);
    const index_t N = basis.size();
    if (N % 2 == 0) return kv;

    gsKnotVector<T> result = kv;

    typedef typename gsKnotVector<T>::const_iterator kIt;
    std::vector<T> diff;
    {
        kIt prev = result.begin();
        for (kIt it = result.begin() + 1; it != result.end(); ++it)
        {
            diff.push_back(*it - *prev);
            prev = it;
        }
    }

    const T maxSpan = *std::max_element(diff.begin(), diff.end());
    std::vector<index_t> maxIdx;
    for (index_t k = 0; k < (index_t)diff.size(); ++k)
        if (math::abs(diff[k] - maxSpan) / maxSpan < T(1e-15))
            maxIdx.push_back(k);

    const index_t i = maxIdx.at(
        static_cast<index_t>(std::ceil(maxIdx.size() / 2.0)) - 1);
    result.insert((result.at(i) + result.at(i + 1)) / T(2));
    return result;
}

// Reproduces unittests/gsQuadratureRules_test.cpp:341-353 verbatim.
gsTHBSplineBasis<2,real_t> makeRefinedTHB(index_t p, real_t pos, index_t mu,
                                           const std::vector<index_t> & boxes)
{
    gsKnotVector<real_t> kv0(0, 1, 0, p + 1, 1);
    kv0.uniformRefine();
    kv0.uniformRefine();          // breakpoints 0, .25, .5, .75, 1

    gsKnotVector<real_t> kv1 = kv0;
    if (mu > 1) kv0.insert(pos, mu - 1);   // direction 0 carries the repeated knot

    gsTensorBSplineBasis<2,real_t> tb(kv0, kv1);
    gsTHBSplineBasis<2,real_t> thb(tb);
    thb.uniformRefine();                                  // level-1 everywhere
    std::vector<index_t> b = boxes;
    thb.refineElements(b);
    return thb;
}

// Reproduces unittests/gsQuadratureRules_test.cpp:458-482 verbatim (thb by value).
gsMatrix<real_t> hpatchAssemble(gsTHBSplineBasis<2,real_t> thb, int form,
                                index_t quRule, index_t quA, index_t quB)
{
    gsTensorBSplineBasis<2,real_t> & l0 = thb.tensorLevel(0);
    gsMatrix<real_t> coefs = l0.anchors().transpose();
    gsTensorBSpline<2,real_t> geom(l0, coefs);
    gsMultiPatch<real_t> mp; mp.addPatch(geom);
    gsMultiBasis<real_t> mb; mb.addBasis(thb.clone().release()); mb.setTopology(mp);

    gsExprAssembler<real_t> A(1,1);
    gsOptionList opt = gsExprAssembler<real_t>::defaultOptions();
    opt.setInt("quRule", quRule);
    opt.setReal("quA", static_cast<real_t>(quA));
    opt.setInt("quB", quB);
    A.setOptions(opt);
    A.setIntegrationElements(mb);
    auto G = A.getMap(mp);
    auto u = A.getSpace(mb);
    A.initSystem();
    if (form == 0) A.assemble( u * u.tr() * meas(G) );
    else           A.assemble( igrad(u,G) * igrad(u,G).tr() * meas(G) );
    return A.matrix().toDense();
}

struct Fixture
{
    std::string name;
    index_t     mu;
    real_t      pos;
    std::vector<index_t> boxes;
};

// Per-leaf, per-direction data needed both for the table and for the verdict.
struct LeafDirData
{
    index_t leafIdx;
    index_t level;
    short_t dir;
    index_t lowerIdx, upperIdx;
    real_t  lowerParam, upperParam;
    std::vector<std::pair<real_t,index_t> > interiorMult; // (value, multiplicity)
    std::vector<real_t> repCur, repNew;
    index_t sizeCur, sizeNew;
    bool    seqDiffer;
    index_t Ncur, Nnew, n1dCur, n1dNew;
    bool    capClamped;
};

int main(int argc, char* argv[])
{
    index_t p = 3;
    bool verbose = false;

    gsCmdLine cmd("Probe: does a THB fixture with a repeated interior knot change "
                  "the gsHPatchRule target-knot-vector construction?");
    cmd.addInt("P","deg","level-0 degree",p);
    cmd.addSwitch("verbose","Print extra diagnostics",verbose);
    try { cmd.getValues(argc,argv); } catch (int rv) { return rv; }

    const index_t m_deg = 2 * p;

    std::vector<Fixture> fixtures;
    fixtures.push_back(Fixture{"A_C0_mu3_at_0.5",        3, 0.5,  {2,0,0,4,4}});
    fixtures.push_back(Fixture{"B_double_mu2_at_0.5",    2, 0.5,  {2,0,0,4,4}});
    fixtures.push_back(Fixture{"C_C0_mu3_at_0.75",       3, 0.75, {2,0,0,4,4}});
    fixtures.push_back(Fixture{"D_C0_mu3_at_0.5_coarse", 3, 0.5,  {2,0,0,2,2}});

    std::vector<index_t> mregs; mregs.push_back(p - 1); mregs.push_back(p - 2); // mass, stiffness

    std::string recipeBlock;
    bool haveWinner = false;

    for (const Fixture & fx : fixtures)
    {
        gsInfo << "\n================ Fixture " << fx.name
               << "  (p=" << p << ", mu=" << fx.mu << ", pos=" << fx.pos
               << ", boxes={" << fx.boxes[0] << "," << fx.boxes[1] << ","
               << fx.boxes[2] << "," << fx.boxes[3] << "," << fx.boxes[4]
               << "}) ================\n";

        gsTHBSplineBasis<2,real_t> thb = makeRefinedTHB(p, fx.pos, fx.mu, fx.boxes);

        for (index_t mregIdx = 0; mregIdx < (index_t)mregs.size(); ++mregIdx)
        {
            const index_t m_reg = mregs[mregIdx];
            const bool isMass = (mregIdx == 0);
            gsInfo << "\n--- m_reg = " << m_reg << " (" << (isMass ? "mass" : "stiffness") << ") ---\n";

            memory::shared_ptr<gsDomain<real_t> > dom = thb.domain();
            const gsHDomain<2,real_t,index_t> * hd =
                dynamic_cast<const gsHDomain<2,real_t,index_t>*>(dom.get());
            GISMO_ENSURE(hd != nullptr, "not a hierarchical domain");

            std::vector<LeafDirData> rows;
            index_t leafCounter = 0;
            index_t predictedFallbackCur = 0;
            index_t nCandTensorCurLeaf, nCandTensorNewLeaf, nGaussTotalLeaf;
            std::vector<index_t> leafNCandCur, leafNCandNew, leafNGauss;
            std::vector<bool> leafFallbackCur, leafFallbackNew;

            for (auto leaf = hd->tree().beginLeafIterator(); leaf.good(); leaf.next(), ++leafCounter)
            {
                const index_t l = leaf.level();
                const gsTensorBSplineBasis<2,real_t> & lvlBasis = thb.tensorLevel(l);

                nCandTensorCurLeaf = 1; nCandTensorNewLeaf = 1; nGaussTotalLeaf = 1;

                for (short_t dd = 0; dd < 2; ++dd)
                {
                    const gsKnotVector<real_t> & kvLev = lvlBasis.component(dd).knots();

                    const index_t lowerIdx = leaf.lowerCorner()[dd];
                    const index_t upperIdx = leaf.upperCorner()[dd];
                    const real_t lowerParam = kvLev(lowerIdx);
                    const real_t upperParam = kvLev(upperIdx);

                    const index_t nElemDD = upperIdx - lowerIdx;
                    nGaussTotalLeaf *= nElemDD * (thb.degree(dd) + 1);

                    LeafDirData row;
                    row.leafIdx = leafCounter; row.level = l; row.dir = dd;
                    row.lowerIdx = lowerIdx; row.upperIdx = upperIdx;
                    row.lowerParam = lowerParam; row.upperParam = upperParam;
                    row.capClamped = false;

                    // ---- CURRENT formula: unique breakpoints, uniform end/interior mult ----
                    std::vector<real_t> uKnots;
                    for (index_t i = lowerIdx; i <= upperIdx; ++i)
                        uKnots.push_back(kvLev(i));
                    gsKnotVector<real_t> targetKvCur(uKnots, (int)m_deg, (int)m_reg);
                    row.sizeCur = gsBSplineBasis<real_t>(targetKvCur).size();
                    for (auto it = targetKvCur.begin(); it != targetKvCur.end(); ++it)
                        row.repCur.push_back(*it);

                    // ---- NEW formula: real multiplicities of kvLev ----
                    std::vector<real_t> repKnots;
                    for (index_t i = lowerIdx; i <= upperIdx; ++i)
                    {
                        const real_t z = kvLev(i);
                        const bool endK = (i == lowerIdx || i == upperIdx);
                        index_t mult;
                        if (endK)
                            mult = m_deg + 1;
                        else
                        {
                            const index_t muj = kvLev.multiplicity(z);
                            row.interiorMult.push_back(std::make_pair(z, muj));
                            const index_t raw = m_deg - m_reg - 1 + muj;
                            if (raw >= m_deg + 1) row.capClamped = true;
                            mult = std::min(raw, m_deg + 1);
                        }
                        repKnots.insert(repKnots.end(), mult, z);
                    }
                    gsKnotVector<real_t> targetKvNew(repKnots, (short_t)m_deg);
                    row.sizeNew = gsBSplineBasis<real_t>(targetKvNew).size();
                    row.repNew = repKnots;

                    row.seqDiffer = (row.repCur.size() != row.repNew.size());
                    if (!row.seqDiffer)
                        for (size_t k = 0; k < row.repCur.size(); ++k)
                            if (math::abs(row.repCur[k] - row.repNew[k]) > 1e-14)
                            { row.seqDiffer = true; break; }

                    gsKnotVector<real_t> evenCur = ensureEven1D(targetKvCur);
                    gsKnotVector<real_t> evenNew = ensureEven1D(targetKvNew);
                    row.Ncur = gsBSplineBasis<real_t>(evenCur).size();
                    row.Nnew = gsBSplineBasis<real_t>(evenNew).size();

                    // V2/V5 evidence (task 03): degree/size identity of the 2-arg
                    // constructor, and whether _ensureEven1D can lift an
                    // already-capped interior multiplicity above m_deg + 1.
                    gsInfo << "    V2: targetKvNew.degree()=" << targetKvNew.degree()
                           << " m_deg=" << m_deg
                           << " sizeNew=" << row.sizeNew
                           << " repKnots.size()-m_deg-1=" << (index_t)repKnots.size() - m_deg - 1 << "\n";
                    gsInfo << "    V5: evenNew.maxInteriorMultiplicity()="
                           << evenNew.maxInteriorMultiplicity()
                           << " m_deg+1=" << (m_deg + 1) << "\n";
                    row.n1dCur = row.Ncur / 2;
                    row.n1dNew = row.Nnew / 2;

                    nCandTensorCurLeaf *= row.n1dCur;
                    nCandTensorNewLeaf *= row.n1dNew;

                    rows.push_back(row);

                    gsInfo << "leaf " << leafCounter << " lvl " << l << " dir " << dd
                           << " [" << lowerIdx << "," << upperIdx << "] span="
                           << (upperIdx - lowerIdx) << " param=[" << lowerParam << ","
                           << upperParam << "]\n";
                    gsInfo << "    interior breakpoints (value, mult@levelL):";
                    if (row.interiorMult.empty()) gsInfo << " (none)";
                    for (auto & pr : row.interiorMult)
                        gsInfo << " (" << pr.first << "," << pr.second << ")";
                    gsInfo << "\n";
                    gsInfo << "    repCUR:"; for (real_t v : row.repCur) gsInfo << " " << v; gsInfo << "\n";
                    gsInfo << "    repNEW:"; for (real_t v : row.repNew) gsInfo << " " << v; gsInfo << "\n";
                    gsInfo << "    size_cur=" << row.sizeCur << " size_new=" << row.sizeNew
                           << " sequences differ: " << (row.seqDiffer ? "yes" : "no") << "\n";
                    gsInfo << "    N_cur=" << row.Ncur << " N_new=" << row.Nnew
                           << " n1d_cur=" << row.n1dCur << " n1d_new=" << row.n1dNew
                           << " cap clamped: " << (row.capClamped ? "yes" : "no") << "\n";
                }

                const bool useFallbackCur = (nCandTensorCurLeaf >= nGaussTotalLeaf);
                const bool useFallbackNew = (nCandTensorNewLeaf >= nGaussTotalLeaf);
                if (useFallbackCur) ++predictedFallbackCur;
                leafNCandCur.push_back(nCandTensorCurLeaf);
                leafNCandNew.push_back(nCandTensorNewLeaf);
                leafNGauss.push_back(nGaussTotalLeaf);
                leafFallbackCur.push_back(useFallbackCur);
                leafFallbackNew.push_back(useFallbackNew);

                gsInfo << "  SUMMARY leaf " << leafCounter
                       << " nCandidateTensor_cur=" << nCandTensorCurLeaf
                       << " nCandidateTensor_new=" << nCandTensorNewLeaf
                       << " nGaussTotal=" << nGaussTotalLeaf
                       << " useFallback_cur=" << (useFallbackCur ? "true" : "false")
                       << " useFallback_new=" << (useFallbackNew ? "true" : "false") << "\n";
            }

            // ---- VERDICT ----
            // Among all direction-0 rows satisfying (a), prefer the first leaf that
            // satisfies (a),(b),(c) fully; otherwise report the first (a)-satisfying
            // leaf and the note explaining which of (b)/(c) it failed.
            index_t winLeaf = -1;
            bool condA = false, condB = false, condC = false;
            bool weakNote = false;
            std::string extraNote;
            for (const LeafDirData & row : rows)
            {
                if (row.dir != 0) continue;
                const index_t span = row.upperIdx - row.lowerIdx;
                bool hasMultInterior = false;
                for (auto & pr : row.interiorMult) if (pr.second >= 2) hasMultInterior = true;
                if (!(span >= 2 && hasMultInterior)) continue;

                const index_t lf = row.leafIdx;
                const bool b = (row.sizeCur != row.sizeNew);
                const bool c = !leafFallbackNew[lf];

                if (winLeaf < 0 || (b && c && !(condB && condC)))
                {
                    winLeaf = lf; condA = true; condB = b; condC = c;
                    extraNote.clear();
                    if (!b && row.seqDiffer)
                        extraNote = "NOTE: sizes coincide but target sequences differ on leaf " +
                                     std::to_string(lf);
                    if (b && !c)
                        extraNote = "NOTE: fallback predicted on leaf " + std::to_string(lf) +
                                    " -- Task 02 must assert the fallback count, not numFallbackCells()==0";
                }
                if (b && c) break; // fully passing leaf found, stop searching
            }

            const bool verdictPass = condA && condB && condC;
            gsInfo << "VERDICT " << fx.name << " m_reg=" << m_reg << " "
                   << (verdictPass ? "PASS" : "FAIL")
                   << " leaf=" << winLeaf << "\n";
            if (!extraNote.empty())
                gsInfo << extraNote << "\n";
            if (condA && condC)
                gsInfo << "useFallback_cur on winning leaf = "
                       << (winLeaf >= 0 ? (leafFallbackCur[winLeaf] ? "true" : "false") : "n/a") << "\n";

            // ---- SELFCHECK ----
            gsHPatchRule<real_t> rule(*thb.domain(), m_deg, m_reg, false);
            const index_t realFallback = rule.numFallbackCells();
            const index_t realPoints   = rule.numPointsTotal();
            const bool selfOk = (predictedFallbackCur <= realFallback);
            gsInfo << "SELFCHECK: predicted_fallback_cur=" << predictedFallbackCur
                   << " rule.numFallbackCells()=" << realFallback
                   << " rule.numPointsTotal()=" << realPoints
                   << " -> " << (selfOk ? "SELFCHECK: OK" : "SELFCHECK: MISMATCH") << "\n";

            // ---- MAXABS ----
            if (condA && condB)
            {
                const int form = isMass ? 0 : 1;
                gsMatrix<real_t> Mg = hpatchAssemble(thb, form, gsQuadrature::GaussLegendre, m_deg, m_reg);
                gsMatrix<real_t> Mh = hpatchAssemble(thb, form, gsQuadrature::HPatchRule,    m_deg, m_reg);
                const real_t maxabs = (Mg - Mh).cwiseAbs().maxCoeff();
                weakNote = (maxabs < 1e-10);
                gsInfo << "MAXABS " << fx.name << " form=" << form << " m_reg=" << m_reg
                       << " maxabs=" << maxabs << (weakNote ? " WEAK" : "") << "\n";

                if (verdictPass && !weakNote)
                {
                    haveWinner = true;
                    recipeBlock += "Fixture " + fx.name + " m_reg=" + std::to_string(m_reg) +
                        " form=" + std::to_string(form) + " maxabs=" + std::to_string(maxabs) + "\n";
                }
            }
        }
    }

    gsInfo << "\n================ FIXTURE RECIPE ================\n";
    if (haveWinner)
    {
        gsInfo << "Winning fixture(s) (paste into unittests/gsQuadratureRules_test.cpp):\n"
               << recipeBlock;
        gsInfo << "p=" << p << "; degree m_deg=2*p; boxes and mu/pos taken from the fixture table above.\n";
    }
    else
    {
        gsInfo << "No candidate fixture satisfied (a),(b),(c) with a non-WEAK maxabs. "
                  "See VERDICT / MAXABS lines above for the failing condition per fixture.\n";
    }

    return 0;
}
