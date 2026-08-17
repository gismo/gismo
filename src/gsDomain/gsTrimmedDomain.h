/** @file gsTrimmedDomain.h

    @brief A basis-independent kd-tree partition of a box into interior,
           exterior and cut (trimmed) elements, driven by the sign of a
           level set.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomain.h>
#include <gsDomain/gsCutTreeData.h>
#include <gsAssembler/gsQuadrature.h>
#include <gsHSplines/gsHTensorBasis.h>

namespace gismo
{

/*
  ===========================================================================
   DESIGN OF gsTrimmedDomain
  ===========================================================================

  WHAT IS STORED
  --------------
  Two things, and deliberately nothing else:

    m_tree    a kd-tree (gsKdTree) whose LEAVES are the elements of the
              domain. Each leaf carries an integer box [lowerCorner,
              upperCorner), a dyadic level L, and a sign in {-1, 0, +1}.
    m_breaks  m_breaks[L][j] = the break vector in direction j at level L.
              Level 0 is the background grid; level L+1 inserts the midpoint
              of every level-L interval (see _ensureLevel()).

  A leaf's corners are indices into m_breaks[L][j] AT ITS OWN LEVEL L, not
  into the level-0 grid. That is the key representational choice: it makes
  the domain SELF-CONTAINED. Once constructed it holds no reference to the
  basis (if any) it was built from, so the basis may be refined, moved or
  destroyed without invalidating the partition, and numElements() is a plain
  span product with no 2^(L*d) correction (see numElements<SignOp>()).

  THE SIGN CONVENTION
  -------------------
    -1  interior  (level set negative at every sample point)
    +1  exterior  (level set positive at every sample point)
     0  cut       (samples disagree => the surface passes through the cell)
  Leaves are selected by the predicate functors below (InteriorSign,
  BoundarySign, ...), which every iterator and count is templated on. Note
  AllSign is {interior + cut}, i.e. the ACTIVE part of an immersed domain.

  A leaf's sign is a SAMPLED verdict, not an exact one: it is decided from
  samples^d Lobatto points inside the leaf (_classifyLeaf()). A feature
  thinner than the sample spacing can therefore be missed. Lobatto rather
  than Gauss is used on purpose -- it includes the endpoints, so cells whose
  intersection with the geometry is a corner sliver are still detected.
  This sampling error is the reason Phase A below refines unconditionally.

  ---------------------------------------------------------------------------
   CLASSIFICATION PIPELINE
  ---------------------------------------------------------------------------
  There are two classifiers, and they differ in one respect only: whether
  the split decision consumes the leaf's sign.

  (1) _classifyTree(samples)  --  uniform background grid
      Used by the bbox and the tensor-B-spline init(). The target partition
      is known in advance: one leaf per level-0 element. The sign plays NO
      part in deciding where to split, which is what lets the work split
      into two phases with no interleaving:

        PHASE A  (serial, ZERO level-set evaluations)
                 Subdivide every leaf that still spans more than one element
                 until each holds exactly one. Unconditional on sign -- see
                 "why unconditional" below.
        PHASE B  (parallel, ALL the level-set evaluations)
                 Classify the terminal leaves from Phase A. Every iteration
                 writes only its own leaf's sign and reads nothing mutable,
                 so this is one flat omp for over independent work: no
                 rounds, no barriers, no reduction.

      Why the phases exist at all: the two used to be interleaved, so every
      MULTI-element node was classified and then immediately split -- but a
      split node is no longer a leaf, and signs are only ever read through
      leaf iterators. Every one of those classifications was discarded. For
      a 16^3 grid that is 4095 of 8191 classifications thrown away, at
      samples^d level-set evaluations each: exactly half the work.

      Why Phase A refines unconditionally: an exterior verdict on a large
      block may be a sampling accident (all samples happen to miss the
      surface). Freezing such a block at block resolution would delete real
      geometry. Splitting regardless of sign costs nothing here, because
      Phase A evaluates the level set zero times.

  (2) _classifyTreeAdaptive(samples, shouldStop, maxLevel)  --  adaptive
      Used by the size-based init() and by the HTB init() (there with a
      shouldStop that always returns true, i.e. classify-only: the HTB tree
      topology was already mirrored by _mirrorHTBTree()).

      Here the split decision genuinely consumes the sign, so Phase A/Phase B
      cannot be separated. Instead the traversal runs in LEVEL-SYNCHRONOUS
      ROUNDS, each round being a small A/B pair:

        per round:  serial   _ensureLevel() up to the round's max level
                             (it grows m_breaks, so it must not run inside
                             the parallel region)
                    PARALLEL classify every pending leaf -- the only step
                             that touches the level set, hence the only one
                             worth parallelising
                    serial   apply shouldStop(), promote survivors to level
                             L+1 (double the integer corners, ++level) and
                             split them; their children are the next round

      Every mutation of the tree and of m_breaks stays on one thread. The
      result is identical to the original depth-first version: each leaf is
      classified exactly once, shouldStop() is called exactly once per leaf
      (it is an arbitrary functor -- assume neither purity nor cheapness),
      and a leaf's fate depends only on its own sign and corners, so
      breadth-first ordering yields the same tree.

      maxLevel bounds the recursion. Leaves that hit it are counted and
      reported once at the end rather than warned about individually.

  THREAD SAFETY
  -------------
  The parallel steps call the virtual sign() concurrently on distinct
  leaves. sign() must therefore be thread-safe. It is declared non-const for
  historical reasons only; that is not licence to make it stateful. Both
  overrides in the tree are fine (gsTrimmedDomain's is unimplemented,
  gsImplicitTrimmedDomain's is a const eval() of the level set). Each thread
  builds its own gsLobattoRule, since mapTo() takes it by non-const ref.
*/

// Cells with sign 0
struct BoundarySign { bool operator()(short_t s) const { return s==0; } };
// Cells with sign -1
struct InteriorSign { bool operator()(short_t s) const { return s< 0; } };
// Cells with sign  1
struct ExteriorSign { bool operator()(short_t s) const { return s> 0; } };
// Cells with signs 0 or -1
struct AllSign      { bool operator()(short_t s) const { return s<=0; } };
// Cells with any sign
struct AnySign      { bool operator()(short_t s) const { return true; } };

/**
 * @brief TODO
 *
 * \ingroup Domain
 */

template<short_t d, class T, class Z=size_t>
class gsTrimmedDomain : public gsDomain<T>
{
public:
    typedef gsCutTreeData<d,Z> TData_t;
    typedef gsKdTree<d,Z,TData_t > Tree_t;
    typedef gsDomainIteratorWrapper<T> domainIter;
    typedef typename Tree_t::const_literator leafIterator;
    typedef typename Tree_t::point_t point_t;

    template <typename SignOp, short_t _d, class _T, class _Z>
    friend class gsTrimmedDomainIterator;

protected:
    Tree_t m_tree; ///< The KdTree partition characterizing in/out/cut elements
    /// m_breaks[level][dim] = break vector at that level for that direction.
    /// Level 0 is the coarsest (background) grid; higher levels are dyadic refinements.
    /// Leaf coordinates (lowerCorner, upperCorner) are stored as indices into m_breaks[L][j]
    /// for their respective level L, making the domain self-contained and independent of any basis.
    std::vector< std::vector< std::vector<T> > > m_breaks;

    /// The polynomial degree used to size quadrature rules over this domain.
    ///
    /// \b IMPORTANT \b -- \b THIS \b MEMBER \b SHOULD \b NOT \b EXIST. It is
    /// the sharpest instance of the design defect documented at
    /// gsDomain::degree(), and removing it is a high-priority item for the next
    /// release. Read that note first; the summary is that gsQuadrature asks the
    /// DOMAIN for a degree that only the BASIS knows, so every domain has to
    /// carry one whether or not it has any business knowing it.
    ///
    /// A trimmed domain is a kd-tree over a box. Two of its four construction
    /// paths never see a basis at all -- the bounding-box+cell-count init() and
    /// the element-size init() -- so on those paths the degree cannot be
    /// derived and must be passed in by the caller. Consequences to be aware of
    /// while the defect stands:
    ///
    ///   - The default is 1. A caller who builds a basis-free trimmed domain and
    ///     forgets \a deg gets linear-order quadrature silently, and every
    ///     integral assembled on that domain is under-integrated with no error
    ///     and no warning at the point of the mistake.
    ///   - Nothing here checks \a deg against the space actually assembled on
    ///     this domain. They can disagree without complaint.
    ///   - One domain cannot serve two spaces of different degree.
    ///
    /// When gsQuadrature::numNodes() is changed to take the order as an
    /// argument, this member, its \a deg parameters, the degree() accessor and
    /// the invariant below all delete together.
    ///
    /// INVARIANT (until then): *every* init() overload must assign m_deg. The
    /// in-class initialiser guarantees a well-defined value meanwhile. The four
    /// sites are init(T,T,index_t,short_t) [element-size, from \a deg],
    /// init(gsMatrix,gsVector,index_t,short_t) [bounding box, from \a deg],
    /// init(gsTensorBSplineBasis,index_t) [from tbasis.maxDegree()] and
    /// init(gsHTensorBasis,index_t) [from htbasis.maxDegree()].
    short_t m_deg = 1;

public: // virtual interface

    virtual gsMatrix<T> boundingBox() const override
    {
        GISMO_NO_IMPLEMENTATION
    }

    /// Computes the sign at given points \a u in the domain
    virtual gsVector<short_t> sign(const gsMatrix<T> & u)
    {
        GISMO_NO_IMPLEMENTATION
    }

public: // non-virtual interface

    inline std::vector<char> inDomain(const gsMatrix<T> & u)
    {
        std::vector<char> res(u.cols());
        gsVector<short_t> val = this->sign(u);
        char * r = res.data();
        for (short_t * a = val.data(); a != val.data()+val.size(); ++a)
            *(r++) = (*a<0 ? false : true);
        return res;
    }

    inline std::vector<char> onBoundary(const gsMatrix<T> & u)
    {
        std::vector<char> res(u.cols());
        gsVector<short_t> val = this->sign(u);
        char * r = res.data();
        for (short_t * a = val.data(); a != val.data()+val.size(); ++a)
            *(r++) = (*a!=0 ? false : true);
        return res;
    }

    domainIter beginAll() const override
    {
        return begin<AllSign>();
    }

    domainIter beginBdr(const boxSide bs) const override
    {
        return begin<BoundarySign>();
    }

    domainIter beginActive() const
    {
        return beginAll();
    }

    domainIter beginExterior() const
    {
        return begin<ExteriorSign>();
    }

    domainIter beginInterior() const
    {
        return begin<InteriorSign>();
    }

    domainIter beginAny() const
    {
        return begin<AnySign>();
    }

    template<typename SignOp>
    domainIter begin() const
    {
        return domainIter(new gsTrimmedDomainIterator<SignOp,d,T,Z>(*this));
    }

    template<typename SignOp>
    domainIter end() const
    {
        return domainIter(new gsDomainIteratorEnd<T>(this->template numElements<SignOp>()));
    }

    /// Returns the number of elements whose sign satisfies SignOp.
    /// Each leaf at level L spanning [lc,uc) in level-0 coords contains
    /// (uc-lc).prod() * 2^(L*d) sub-elements (one per level-L break cell).
    template<typename SignOp>
    size_t numElements() const
    {
        leafIterator it = m_tree.beginLeafIterator();
        size_t nel(0);
        while (it.good())
        {
            if (SignOp()(it.data().sign()))
            {
                const point_t & uc = it.data().upperCorner();
                const point_t & lc = it.data().lowerCorner();
                // With level-L leaf coordinates, element count is simply the span.
                // No need to multiply by 2^L since coordinates index the level-L break vector.
                size_t leafNel = 1;
                for (short_t j = 0; j < d; ++j)
                    leafNel *= static_cast<size_t>(uc[j] - lc[j]);
                nel += leafNel;
            }
            it.next();
        }
        return nel;
    }

    size_t numElements() const override
    {
        return numElements<AllSign>();
    }

    size_t numElementsBdr(boxSide const & s = boundary::none) const override
    {
        return numElements<BoundarySign>();
    }

    short_t dim() const override { return d; }

    short_t degree(short_t i = 0) const
    {GISMO_UNUSED(i); return m_deg;}

    const Tree_t & tree() const { return m_tree; }

    /// Returns break vector at given \a level for direction \a j.
    const std::vector<T> & breaks(unsigned level, short_t j) const
    {
        GISMO_ASSERT(level < m_breaks.size(), "Level " << level << " out of range");
        return m_breaks[level][j];
    }

    unsigned numLevels() const { return static_cast<unsigned>(m_breaks.size()); }

protected:

    /// Dyadically subdivides breaks up to and including \a level.
    /// m_breaks[0] must already be populated; each subsequent level
    /// inserts midpoints between every consecutive pair of parent breaks.
    void _ensureLevel(unsigned level)
    {
        while (m_breaks.size() <= static_cast<size_t>(level))
        {
            const unsigned parentL = static_cast<unsigned>(m_breaks.size()) - 1;
            m_breaks.push_back(std::vector<std::vector<T>>(d));
            for (short_t j = 0; j < d; ++j)
            {
                const std::vector<T> & parentBk = m_breaks[parentL][j];
                std::vector<T> & childBk = m_breaks.back()[j];
                childBk.reserve(parentBk.size() * 2 - 1);
                for (size_t i = 0; i + 1 < parentBk.size(); ++i)
                {
                    childBk.push_back(parentBk[i]);
                    childBk.push_back((parentBk[i] + parentBk[i+1]) / T(2));
                }
                childBk.push_back(parentBk.back());
            }
        }
    }

    /// Returns the parametric lower and upper corners of \a node from m_breaks.
    void _getParamCorners(const Tree_t * node,
                          gsVector<T,d> & u1, gsVector<T,d> & u2) const
    {
        const point_t & lc = node->nodeData().lowerCorner();
        const point_t & uc = node->nodeData().upperCorner();
        const index_t   L  = node->nodeData().level();
        for (short_t j = 0; j < d; ++j)
        {
            u1[j] = m_breaks[L][j][static_cast<index_t>(lc[j])];
            u2[j] = m_breaks[L][j][static_cast<index_t>(uc[j])];
        }
    }

    /// Evaluates sign at quadrature points within \a node's parametric box and
    /// assigns the result to the node's sign field. Uses \a QuRule for quadrature.
    ///
    /// \note Called concurrently from the classification loops below, on
    /// distinct \a node pointers. It touches no member state other than
    /// reading m_breaks (via _getParamCorners) and calling the virtual
    /// sign(), so a derived class whose sign() override mutates its own
    /// members must synchronise them itself. The two overrides in the tree --
    /// gsTrimmedDomain's (unimplemented) and gsImplicitTrimmedDomain's (a
    /// const eval() of the level set) -- are both fine. Note that sign() is
    /// declared non-const for historical reasons; that is not licence to make
    /// it stateful.
    void _classifyLeaf(Tree_t * node, gsLobattoRule<T> & QuRule)
    {
        gsVector<T,d>  u1, u2;
        gsVector<T>    wts;
        gsMatrix<T>    pts;
        _getParamCorners(node, u1, u2);
        QuRule.mapTo(u1, u2, pts, wts);
        gsVector<short_t> sgn = this->sign(pts);
        node->nodeData().sign() =
            (sgn.array() ==  1).all() ?  1 :
            (sgn.array() == -1).all() ? -1 : 0;
    }

    /// Classifies all leaves of m_tree using m_breaks[0] for index->parametric
    /// conversion. Mixed leaves spanning more than one level-0 element are split
    /// recursively; single-element mixed leaves are marked as boundary (sign=0).
    void _classifyTree(index_t samples)
    {
        // Two phases, no interleaving -- see "CLASSIFICATION PIPELINE (1)" in
        // the design block at the top of this file for why the subdivision and
        // the classification must NOT be merged back together.

        // --- Phase A (serial, no level-set evaluation): subdivide until every
        // leaf holds a single element, regardless of sign.
        std::vector<Tree_t*> leaves;
        std::vector<Tree_t*> stack;
        stack.reserve(std::pow(4, d));
        stack.push_back(&m_tree);

        while (!stack.empty())
        {
            Tree_t * cur = stack.back();
            stack.pop_back();

            if (!cur->isLeaf())
            {
                stack.push_back(cur->left);
                stack.push_back(cur->right);
                continue;
            }

            const point_t & lc = cur->nodeData().lowerCorner();
            const point_t & uc = cur->nodeData().upperCorner();
            if ((uc - lc).prod() > 1)
            {
                cur->anyMidSplit(1);
                stack.push_back(cur->left);
                stack.push_back(cur->right);
                continue;
            }
            leaves.push_back(cur);   // terminal: classified below, never split
        }

        // --- Phase B (parallel): classify the terminal leaves. Data-race free
        // given a thread-safe sign(); see "THREAD SAFETY" in the design block.
        gsVector<index_t> numNodes(d);
        numNodes.setConstant(samples);
        const index_t nLeaves = static_cast<index_t>(leaves.size());

#       pragma omp parallel
        {
            gsLobattoRule<T> QuRule(numNodes);
#           pragma omp for schedule(static)
            for (index_t i = 0; i < nLeaves; ++i)
                _classifyLeaf(leaves[i], QuRule);
        }
    }

    /// Adaptive refinement of the cut-tree driven by a stopping predicate.
    ///
    /// Traverses all leaves. For each leaf, evaluates the sign at \a samples
    /// quadrature points per direction. Assigns sign (-1/0/1). If the predicate
    /// \a shouldStop returns false, splits the leaf and promotes it to the next
    /// dyadic level (doubling the index coords and incrementing level by 1).
    ///
    /// \a shouldStop is called as shouldStop(sign, u1, u2) and should return
    /// true when no further refinement is needed for this leaf.
    /// \a maxLevel caps the recursion depth (refinement halves the leaf volume
    /// per level, so the default of 20 admits a ~1e-6 reduction of the root box
    /// while bounding the leaf count at ~1e6).
    template<typename StopPred>
    void _classifyTreeAdaptive(index_t samples, StopPred shouldStop, unsigned maxLevel = 20)
    {
        // LEVEL-SYNCHRONOUS ROUNDS: classify a round's leaves in parallel, then
        // serially apply the predicate and split. See "CLASSIFICATION PIPELINE
        // (2)" in the design block at the top of this file for why the phases
        // cannot be separated here as they are in _classifyTree(), and for the
        // equivalence with the original depth-first traversal.
        gsVector<index_t> numNodes(d);
        numNodes.setConstant(samples);

        std::vector<Tree_t*> pending;
        std::vector<Tree_t*> stack;
        stack.reserve(std::pow(4, d));
        stack.push_back(&m_tree);
        while (!stack.empty())
        {
            Tree_t * cur = stack.back();
            stack.pop_back();
            if (cur->isLeaf())
                pending.push_back(cur);
            else
            {
                stack.push_back(cur->left);
                stack.push_back(cur->right);
            }
        }

        index_t nCapped = 0;

        while (!pending.empty())
        {
            // Serial: materialise every break level this round will read.
            // _ensureLevel() grows m_breaks, so it must not run inside the
            // parallel region -- hoisting it to the round's maximum level is
            // equivalent, since it fills in all levels up to its argument.
            index_t maxL = 0;
            for (size_t i = 0; i < pending.size(); ++i)
                maxL = std::max(maxL, pending[i]->nodeData().level());
            _ensureLevel(static_cast<unsigned>(maxL));

            // Parallel: the only part that evaluates the level set.
            const index_t nPending = static_cast<index_t>(pending.size());
#           pragma omp parallel
            {
                gsLobattoRule<T> QuRule(numNodes);
#               pragma omp for schedule(static)
                for (index_t i = 0; i < nPending; ++i)
                    _classifyLeaf(pending[i], QuRule);
            }

            // Serial: predicate, promotion and split.
            std::vector<Tree_t*> next;
            for (size_t i = 0; i < pending.size(); ++i)
            {
                Tree_t * cur = pending[i];

                gsVector<T,d> u1, u2;
                _getParamCorners(cur, u1, u2);
                if (shouldStop(cur->nodeData().sign(), u1, u2))
                    continue;

                const index_t L = cur->nodeData().level();
                if (static_cast<unsigned>(L) >= maxLevel)
                {
                    ++nCapped;
                    continue;
                }

                const index_t newL = L + 1;
                _ensureLevel(static_cast<unsigned>(newL));

                point_t & lc = cur->nodeData().lowerCorner();
                point_t & uc = cur->nodeData().upperCorner();
                for (short_t j = 0; j < d; ++j)
                {
                    lc[j] = static_cast<Z>(lc[j]) * 2;
                    uc[j] = static_cast<Z>(uc[j]) * 2;
                }
                cur->nodeData().level() = newL;

                short_t splitAxis = 0;
                Z maxSpan = 0;
                for (short_t j = 0; j < d; ++j)
                {
                    Z span = uc[j] - lc[j];
                    if (span > maxSpan) { maxSpan = span; splitAxis = j; }
                }
                Z splitPos = lc[splitAxis] + (uc[splitAxis] - lc[splitAxis]) / 2;
                cur->split(splitAxis, splitPos);
                next.push_back(cur->left);
                next.push_back(cur->right);
            }
            pending.swap(next);
        }

        if (nCapped)
            gsWarn << "gsTrimmedDomain: " << nCapped << " leaf/leaves reached the "
                      "maximum refinement level " << maxLevel << "; refinement stopped "
                      "there. Increase maxLevel or relax the element-size criterion.\n";
    }

    /// Size-based adaptive init. Refines until each leaf has volume <= maxElementSize,
    /// or is a cut cell with volume <= minElementSize. Leaves are tracked at their
    /// current dyadic level; m_breaks is extended on demand via _ensureLevel().
    /// Since there is no basis on this path, the polynomial degree used to size
    /// quadrature rules (see degree()) must be supplied explicitly by \a deg.
    /// \note This overload expects a \c d x 2 matrix from boundingBox() /
    /// support(): row j = direction j, column 0 = lower corner, column 1 = upper
    /// corner. The bbox overload below expects the transpose (2 x d, caller-supplied).
    void init(T maxElementSize = T(1), T minElementSize = T(0.1),
              index_t samples  = 10,   short_t deg      = 1)
    {
        GISMO_ASSERT(deg >= 0, "Polynomial degree must be non-negative, got "
                     << deg << ".");
        m_deg = deg;
        gsMatrix<T> bb = this->boundingBox();

        GISMO_ENSURE(bb.size() != 0,
                     "gsTrimmedDomain: the element-size init() needs a finite bounding "
                     "box, but boundingBox() returned an empty matrix. For an implicit "
                     "domain this means the level-set function has no support() (e.g. "
                     "gsFunctionExpr, whose domain of definition is all of R^d). Use the "
                     "gsImplicitTrimmedDomain(fnc, bbox, numCells, samples, deg) "
                     "constructor and supply the box explicitly (bbox is 2 x d: row 0 = "
                     "lower corners, row 1 = upper corners).");
        GISMO_ENSURE(bb.rows() == d && bb.cols() == 2,
                     "gsTrimmedDomain: boundingBox() must be d x 2 (row = direction, "
                     "col 0/1 = lower/upper corner), got " << bb.rows() << " x "
                     << bb.cols() << " for d = " << d << ".");

        m_breaks.clear();
        m_breaks.resize(1);
        m_breaks[0].resize(d);
        for (short_t j = 0; j < d; ++j)
        {
            m_breaks[0][j].push_back(static_cast<T>(bb(j, 0)));
            m_breaks[0][j].push_back(static_cast<T>(bb(j, 1)));
        }

        point_t rootUpper;
        for (short_t j = 0; j < d; ++j)
            rootUpper[j] = static_cast<Z>(1);

        m_tree = Tree_t(TData_t(point_t::Zero(), rootUpper, 0));

        _classifyTreeAdaptive(samples, [maxElementSize, minElementSize]
            (short_t leafSign, const gsVector<T,d> & u1, const gsVector<T,d> & u2) -> bool
        {
            const T volume = (u2 - u1).prod();
            if (leafSign != 0)
                return volume <= maxElementSize;
            else
                return volume <= minElementSize;
        });
    }

    /// Initializes from a bounding box + uniform cell count per direction.
    /// Since there is no basis on this path, the polynomial degree used to size
    /// quadrature rules (see degree()) must be supplied explicitly by \a deg.
    /// \note This overload expects a \c 2 x d matrix (\a bbox): row 0 = lower
    /// corners, row 1 = upper corners. The element-size init above expects the
    /// transpose (d x 2, from support()).
    void init(const gsMatrix<T>         & bbox,
              const gsVector<index_t,d> & numCells,
              index_t samples = 5,
              short_t deg     = 1)
    {
        GISMO_ASSERT(deg >= 0, "Polynomial degree must be non-negative, got "
                     << deg << ".");
        m_deg = deg;
        m_breaks.clear();
        m_breaks.resize(1);
        m_breaks[0].resize(d);
        for (short_t j = 0; j < d; ++j)
        {
            const T lo = bbox(0, j);
            const T hi = bbox(1, j);
            gsKnotVector<T> kv(lo, hi, numCells[j] - 1, 2);
            m_breaks[0][j] = kv.breaks();
        }

        point_t upperIndex;
        for (short_t j = 0; j < d; ++j)
            upperIndex[j] = static_cast<Z>(numCells[j]);

        m_tree = Tree_t(TData_t(point_t::Zero(), upperIndex, 0));
        _classifyTree(samples);
    }

    /// Initializes from a tensor B-spline basis (uses its knot vectors as the level-0 grid).
    void init(const gsTensorBSplineBasis<d,T> & tbasis, index_t samples = 5)
    {
        m_deg = tbasis.maxDegree();
        m_breaks.clear();
        m_breaks.resize(1);
        m_breaks[0].resize(d);
        for (short_t j = 0; j < d; ++j)
            m_breaks[0][j] = tbasis.knots(j).breaks();

        point_t upperIndex;
        for (short_t j = 0; j < d; ++j)
            upperIndex[j] = static_cast<Z>(tbasis.knots(j).numElements());

        m_tree = Tree_t(TData_t(point_t::Zero(), upperIndex, 0));
        _classifyTree(samples);
    }

    /// Mirrors the HTensorBasis tree structure into m_tree.
    /// Populates m_breaks with per-level break vectors from the HTB's knot vectors.
    /// Sets leaf geometry (lowerCorner, upperCorner, level) to level-L element indices.
    /// Does NOT classify leaves (no sign assignment).
    ///
    /// Leaf coordinates are stored as level-L element indices (indexing into m_breaks[L][j]).
    /// This makes the trimmed domain self-contained and independent of the basis.
    /// Internal node split positions are stored in level-0 element indices for consistency
    /// with the paired traversal logic, though these are not used after construction.
    ///
    /// \param htbasis The hierarchical tensor basis to mirror.
    void _mirrorHTBTree(const gsHTensorBasis<d,T> & htbasis)
    {
        typedef gsHTree<d,index_t>                         HTree_t;
        typedef gsKdTree<d,index_t,gsHTreeData<d,index_t>> HNode_t;
        typedef gsVector<index_t,d>                        HPoint_t;

        const HTree_t & htree  = htbasis.tree();
        const unsigned  idxLv  = htree.getIndexLevel();
        const unsigned  maxLvl = htbasis.maxLevel();

        // Extract break vectors for all levels from the HTB knot vectors
        m_breaks.clear();
        m_breaks.resize(maxLvl + 1);
        for (unsigned L = 0; L <= maxLvl; ++L)
        {
            const gsTensorBSplineBasis<d,T> * tbL =
                dynamic_cast<const gsTensorBSplineBasis<d,T>*>(htbasis.getBases()[L]);
            GISMO_ASSERT(tbL != nullptr,
                         "gsHTensorBasis level-" << L << " basis is not a gsTensorBSplineBasis");
            m_breaks[L].resize(d);
            for (short_t j = 0; j < d; ++j)
                m_breaks[L][j] = tbL->knots(j).breaks();
        }

        // Root of cut-tree covers same level-0 element range as the HTB.
        // (Temporarily stored; will be overwritten when leaves are populated.)
        const point_t upp0 = htree.upperCornerIndex().template cast<Z>();
        m_tree = Tree_t(TData_t(point_t::Zero(), upp0, 0));

        // Paired traversal: mirror HTB tree structure into cut-tree
        typedef std::pair<const HNode_t*, Tree_t*> NodePair;
        std::vector<NodePair> stack;
        stack.reserve(128);
        stack.emplace_back(&htree, &m_tree);

        while (!stack.empty())
        {
            auto [htNode, cutNode] = stack.back();
            stack.pop_back();

            if (htNode->isLeaf())
            {
                const int L = htNode->nodeData().level();
                // Stride: how many global-index units per one level-L unit
                const index_t h = static_cast<index_t>(1) << (idxLv - L);

                const HPoint_t & gk1 = htNode->nodeData().lowerCorner();
                const HPoint_t & gk2 = htNode->nodeData().upperCorner();

                // Convert global-index-level coords to level-L local indices.
                // Leaf coordinates are stored as level-L element indices.
                point_t lcL, ucL;
                for (short_t j = 0; j < d; ++j)
                {
                    const index_t k1L = static_cast<index_t>(gk1[j]) / h;
                    const index_t k2L = static_cast<index_t>(gk2[j]) / h;
                    // Clamp to valid range within level-L break vector
                    const index_t nElemL = static_cast<index_t>(m_breaks[L][j].size()) - 1;
                    lcL[j] = static_cast<Z>(k1L);
                    ucL[j] = static_cast<Z>(std::min(k2L, nElemL));
                }

                cutNode->nodeData().lowerCorner() = lcL;
                cutNode->nodeData().upperCorner() = ucL;
                cutNode->nodeData().level() = static_cast<index_t>(L);
            }
            else
            {
                // Convert the HTB split position (global-index-level) to level-0 coords
                const int splitAxis  = htNode->axis;
                const index_t gPos   = static_cast<index_t>(htNode->pos);

                // We need to figure out the level of the children to compute the
                // level-0 split position. Use the left child's level as reference.
                // For internal nodes, we simply split the parent's level-0 range
                // at the mid-point corresponding to the HTB split.
                const point_t & parentLc = cutNode->nodeData().lowerCorner();
                const point_t & parentUc = cutNode->nodeData().upperCorner();

                // Determine the level this split operates at.
                // The HTB internal node splits at a position expressed in global-index coords.
                // The left child spans [parentGk1, gPos], right spans [gPos, parentGk2].
                // We compute the level-0 element index at gPos using the coarsest applicable level.
                // Use level of the cut-node (inherited from its construction) as a guide.
                const index_t parentL = cutNode->nodeData().level();
                const index_t h = static_cast<index_t>(1) << (idxLv - static_cast<unsigned>(parentL));
                const index_t kL  = gPos / h;              // level-parentL index
                const index_t c   = kL >> parentL;          // level-0 element index

                Z cutPos = static_cast<Z>(c);
                // Clamp to a valid split (must be strictly between parent corners)
                if (cutPos <= parentLc[splitAxis]) cutPos = parentLc[splitAxis] + static_cast<Z>(1);
                if (cutPos >= parentUc[splitAxis]) cutPos = parentUc[splitAxis] - static_cast<Z>(1);

                cutNode->split(static_cast<int>(splitAxis), cutPos);

                stack.emplace_back(htNode->left,  cutNode->left);
                stack.emplace_back(htNode->right, cutNode->right);
            }
        }
    }

    /// Initializes from a hierarchical tensor basis.
    ///
    /// First mirrors the HTensorBasis tree topology via _mirrorHTBTree(), then
    /// classifies all leaves via _classifyTreeAdaptive() with no refinement
    /// (since HTB structure is already fully specified).
    ///
    /// Leaf coordinates are stored as level-L element indices (indexing into m_breaks[L][j]),
    /// making the domain self-contained and independent of the original basis.
    ///
    /// \param htbasis The hierarchical tensor basis to initialize from.
    /// \param samples Number of quadrature points per direction for sign classification (default: 5).
    void init(const gsHTensorBasis<d,T> & htbasis,
              index_t samples = 5)
    {
        m_deg = htbasis.maxDegree();
        _mirrorHTBTree(htbasis);
        _classifyTreeAdaptive(samples, [](short_t sign, const gsVector<T,d>& u1, const gsVector<T,d>& u2) { return true; });
    }
};

} // namespace gismo
