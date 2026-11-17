/** @file htree_connect_boxes_probe.cpp

    @brief Probe of the private gsHTree::connect_Boxes through the public box queries.

    gsHTree::getBoxes() and gsHTree::getBoxesInLevelIndex() list the leaves of a
    hierarchical tree and then merge same-level leaves into larger boxes with the
    private function connect_Boxes(). That function contains a comparison of a
    box with itself, so a pair of mergeable boxes is merged only when the box
    that comes first in the list is the upper one along the merge axis.

    This driver asks whether the public queries ever return two boxes that
    should have been merged. It builds trees and hierarchical bases (aimed cases
    and exhaustive families), counts the unmerged pairs in the returned boxes
    and prints one line per case. It reports and never asserts: the return value
    is EXIT_SUCCESS whatever the outcome. The last line is
    "probe-verdict: reproduced" or "probe-verdict: not-reproduced".

    A verbatim copy of connect_Boxes and a corrected variant are run on the
    leaf lists of the aimed cases, in tree order and in reversed order, to show
    what the defect does independently of whether the tree can produce a
    list that exposes it.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): Hugo Verhelst
*/

#include <gismo.h>

#include <algorithm>
#include <functional>
#include <sstream>

using namespace gismo;

namespace
{

/// A box as stored by gsHTree::getBoxes_vec(): [lower(0..d-1), upper(0..d-1), level].
typedef std::vector<index_t> Box;

/// One call of gsHTree::insertBox(): corners in the index grid of level \a lvl.
struct Insertion
{
    index_t lvl;
    std::vector<index_t> lo;
    std::vector<index_t> up;
};

/// Insertion whose corners are given in level-0 element units; they are scaled to level \a lvl.
Insertion boxFromLevel0(index_t lvl, std::vector<index_t> lo, std::vector<index_t> up)
{
    for (size_t i = 0; i != lo.size(); ++i)
    {
        lo[i] <<= lvl;
        up[i] <<= lvl;
    }
    Insertion ins = {lvl, lo, up};
    return ins;
}

/// Insertion whose corners are already indices of level \a lvl.
Insertion boxAtLevel(index_t lvl, std::vector<index_t> lo, std::vector<index_t> up)
{
    Insertion ins = {lvl, lo, up};
    return ins;
}

/// Index level of the bare trees; insertions up to this level are allowed.
const unsigned treeIndexLevel = 4;

/// Set when any real query returned an unmerged pair.
bool reproduced = false;

/// Copy of gsHTree<d,Z>::connect_Boxes, src/gsDomain/gsHTree.hpp lines 915-985, with Z = index_t.
template<short_t d> void
connectBoxesAsInTree(std::vector<std::vector<index_t> > & boxes)
{
    bool change = true;
    while(change)
    {
        change =  false;
        size_t s = boxes.size();
        for(size_t i = 0; i < s; i++)
        {
            for(size_t j = i+1; j < s; j++)
            {
                if(boxes[i][2*d]==boxes[j][2*d]) // if( the levels are the same )
                {
                    unsigned nmCoordLo = 0;
                    unsigned nmCoordUp = 0;
                    unsigned nmCountUp = 0;
                    unsigned nmCountLo = 0; //...the "nm" is for non-matching

                    // Compare the lower and upper corners of the boxes
                    // coordinate-wise, and check if there are differences.
                    // If there are differences, count and store the coordinate
                    for(short_t k=0; k < d; k++)
                    {
                        if( boxes[i][k] != boxes[j][k] )
                        {
                            nmCountLo++;
                            nmCoordLo = k;
                        }

                        if( boxes[i][d+k] != boxes[j][d+k] )
                        {
                            nmCountUp++;
                            nmCoordUp = k;
                        }
                    }

                    // The boxes can only be merged if
                    // the lower and upper corners are the same,
                    // except in one coordinate direction.
                    if( nmCountLo == 1
                        && nmCountUp == 1
                        && nmCoordLo == nmCoordUp )
                    {

                        if( boxes[i][nmCoordLo] == boxes[j][d+nmCoordUp] )
                        {
                            // box i is "on top" of box j.
                            // It inherits the lower corner from box j:
                            boxes[i][nmCoordLo] = boxes[j][nmCoordLo];
                            boxes.erase( boxes.begin()+j );
                            s--;
                            j--;
                            change = true;
                        }

                        if( boxes[i][d+nmCoordUp] == boxes[i][nmCoordLo] )
                        {
                            // box i is "below" of box j.
                            // It inherits the upper corner from box j:
                            boxes[i][d+nmCoordUp] = boxes[j][d+nmCoordUp];
                            boxes.erase( boxes.begin()+j );
                            s--;
                            j--;
                            change = true;
                        }
                    }
                } // if boxes same level
            } // for j
        } // for i
    }
}

/// Variant of connectBoxesAsInTree() with the second test replaced by the intended
/// comparison against box j, guarded by an else branch so that it never runs on a
/// pair whose box j was already erased by the first test.
template<short_t d> void
connectBoxesFixed(std::vector<std::vector<index_t> > & boxes)
{
    bool change = true;
    while(change)
    {
        change =  false;
        size_t s = boxes.size();
        for(size_t i = 0; i < s; i++)
        {
            for(size_t j = i+1; j < s; j++)
            {
                if(boxes[i][2*d]==boxes[j][2*d])
                {
                    unsigned nmCoordLo = 0;
                    unsigned nmCoordUp = 0;
                    unsigned nmCountUp = 0;
                    unsigned nmCountLo = 0;

                    for(short_t k=0; k < d; k++)
                    {
                        if( boxes[i][k] != boxes[j][k] )
                        {
                            nmCountLo++;
                            nmCoordLo = k;
                        }

                        if( boxes[i][d+k] != boxes[j][d+k] )
                        {
                            nmCountUp++;
                            nmCoordUp = k;
                        }
                    }

                    if( nmCountLo == 1
                        && nmCountUp == 1
                        && nmCoordLo == nmCoordUp )
                    {
                        if( boxes[i][nmCoordLo] == boxes[j][d+nmCoordUp] )
                        {
                            boxes[i][nmCoordLo] = boxes[j][nmCoordLo];
                            boxes.erase( boxes.begin()+j );
                            s--;
                            j--;
                            change = true;
                        }
                        else if( boxes[i][d+nmCoordUp] == boxes[j][nmCoordLo] )
                        {
                            boxes[i][d+nmCoordUp] = boxes[j][d+nmCoordUp];
                            boxes.erase( boxes.begin()+j );
                            s--;
                            j--;
                            change = true;
                        }
                    }
                }
            }
        }
    }
}

/// Number of unordered pairs of boxes that a pairwise merge must have joined: same level,
/// equal corners in all axes but one, and touching along that axis.
template<short_t d> index_t
leftoverPairs(const std::vector<Box> & boxes)
{
    index_t count = 0;
    for (size_t a = 0; a != boxes.size(); ++a)
    {
        for (size_t b = a + 1; b != boxes.size(); ++b)
        {
            const Box & A = boxes[a];
            const Box & B = boxes[b];
            if (A[2*d] != B[2*d])
                continue;

            short_t differingAxes = 0;
            short_t axis = 0;
            for (short_t k = 0; k != d; ++k)
            {
                if (A[k] != B[k] || A[d+k] != B[d+k])
                {
                    ++differingAxes;
                    axis = k;
                }
            }
            if (differingAxes == 1
                && (A[d+axis] == B[axis] || B[d+axis] == A[axis]))
                ++count;
        }
    }
    return count;
}

/// The two public queries that run connect_Boxes.
const char * queryName(index_t q)
{
    return q == 0 ? "getBoxes" : "getBoxesInLevelIndex";
}

/// Run public query \a q on \a tree and return the boxes in the Box layout.
template<short_t d> std::vector<Box>
runQuery(const gsHTree<d,index_t> & tree, index_t q)
{
    gsMatrix<index_t> b1, b2;
    gsVector<index_t> level;
    if (q == 0)
        tree.getBoxes(b1, b2, level);
    else
        tree.getBoxesInLevelIndex(b1, b2, level);

    std::vector<Box> result(b1.rows());
    for (index_t r = 0; r != b1.rows(); ++r)
    {
        for (short_t k = 0; k != d; ++k)
            result[r].push_back(b1(r, k));
        for (short_t k = 0; k != d; ++k)
            result[r].push_back(b2(r, k));
        result[r].push_back(level[r]);
    }
    return result;
}

/// One line per box: "  lvl=<l> lo=(..) up=(..)".
template<short_t d> std::string
formatBoxes(const std::vector<Box> & boxes)
{
    std::ostringstream os;
    for (size_t r = 0; r != boxes.size(); ++r)
    {
        os << "  lvl=" << boxes[r][2*d] << " lo=(";
        for (short_t k = 0; k != d; ++k)
            os << (k ? "," : "") << boxes[r][k];
        os << ") up=(";
        for (short_t k = 0; k != d; ++k)
            os << (k ? "," : "") << boxes[r][d+k];
        os << ")\n";
    }
    return os.str();
}

/// The insertions as text, in the order they are applied.
std::string
formatInsertions(const std::vector<Insertion> & ins)
{
    std::ostringstream os;
    for (size_t i = 0; i != ins.size(); ++i)
    {
        os << (i ? "; " : "") << "l" << ins[i].lvl << " [";
        for (size_t k = 0; k != ins[i].lo.size(); ++k)
            os << (k ? "," : "") << ins[i].lo[k];
        os << "]-[";
        for (size_t k = 0; k != ins[i].up.size(); ++k)
            os << (k ? "," : "") << ins[i].up[k];
        os << "]";
    }
    return os.str();
}

/// Build a bare tree on a grid of \a grid level-0 elements by applying \a ins in order.
/// With \a compress the sibling leaves of equal level are merged afterwards, as
/// gsHTensorBasis::update_structure() does.
template<short_t d> void
buildTree(gsHTree<d,index_t> & tree, const std::vector<index_t> & grid,
          const std::vector<Insertion> & ins, bool compress)
{
    gsVector<index_t,d> upp;
    for (short_t k = 0; k != d; ++k)
        upp[k] = grid[k];
    tree.init(upp, treeIndexLevel);

    for (size_t i = 0; i != ins.size(); ++i)
    {
        gsVector<index_t,d> lo, up;
        for (short_t k = 0; k != d; ++k)
        {
            lo[k] = ins[i].lo[k];
            up[k] = ins[i].up[k];
        }
        tree.insertBox(lo, up, static_cast<int>(ins[i].lvl));
    }
    if (compress)
        tree.makeCompressed();
}

/// The leaves of \a tree in leaf-iterator order, with corners in the index grid of
/// getMaxInsLevel(): the list that getBoxes_vec() hands to connect_Boxes.
template<short_t d> std::vector<Box>
leafList(const gsHTree<d,index_t> & tree)
{
    const index_t maxIns = tree.getMaxInsLevel();
    std::vector<Box> list;
    for (typename gsHTree<d,index_t>::const_literator it = tree.beginLeafIterator();
         it.good(); it.next())
    {
        const gsVector<index_t,d> lo = it.lowerCorner();
        const gsVector<index_t,d> up = it.upperCorner();
        const index_t shift = maxIns - it.level();
        Box b;
        for (short_t k = 0; k != d; ++k)
            b.push_back(lo[k] << shift);
        for (short_t k = 0; k != d; ++k)
            b.push_back(up[k] << shift);
        b.push_back(it.level());
        list.push_back(b);
    }
    return list;
}

/// Hand-derived scenario: the insertions in a base order, the orders to run them in
/// (index permutations of the base) and the box count the query must return.
struct AimedCase
{
    std::string name;
    std::vector<index_t> grid;
    std::vector<Insertion> base;
    std::vector<std::pair<std::string, std::vector<size_t> > > orders;
    size_t expected;
};

/// Orders "asc" and "desc" for \a n insertions.
std::vector<std::pair<std::string, std::vector<size_t> > >
ascDesc(size_t n)
{
    std::vector<size_t> asc(n);
    for (size_t i = 0; i != n; ++i)
        asc[i] = i;
    std::vector<size_t> desc(asc.rbegin(), asc.rend());
    std::vector<std::pair<std::string, std::vector<size_t> > > orders;
    orders.push_back(std::make_pair(std::string("asc"), asc));
    orders.push_back(std::make_pair(std::string("desc"), desc));
    return orders;
}

/// Orders "asc", "desc" and the interleaved "1,3,0,2" for four insertions.
std::vector<std::pair<std::string, std::vector<size_t> > >
ascDescInterleaved()
{
    std::vector<std::pair<std::string, std::vector<size_t> > > orders = ascDesc(4);
    const size_t interleaved[] = {1, 3, 0, 2};
    orders.push_back(std::make_pair(std::string("1,3,0,2"),
                                    std::vector<size_t>(interleaved, interleaved + 4)));
    return orders;
}

/// Print the result of one query on one path of an aimed case.
template<short_t d> void
reportAimed(const AimedCase & c, const std::string & order, const char * path,
            const gsHTree<d,index_t> & tree)
{
    for (index_t q = 0; q != 2; ++q)
    {
        const std::vector<Box> boxes = runQuery<d>(tree, q);
        const index_t p = leftoverPairs<d>(boxes);
        if (p > 0)
            reproduced = true;
        gsInfo << "case=" << c.name << " order=" << order << " path=" << path
               << " query=" << queryName(q) << " leftoverPairs=" << p
               << " expected=" << c.expected << " got=" << boxes.size() << "\n"
               << formatBoxes<d>(boxes);
    }
}

/// Run one aimed case in every order on the three paths (compressed bare tree,
/// uncompressed bare tree, hierarchical basis), then the replica comparison on the
/// uncompressed tree.
template<short_t d> void
runAimed(const AimedCase & c)
{
    for (size_t o = 0; o != c.orders.size(); ++o)
    {
        const std::string & label = c.orders[o].first;
        std::vector<Insertion> ins;
        for (size_t i = 0; i != c.orders[o].second.size(); ++i)
            ins.push_back(c.base[c.orders[o].second[i]]);

        // Bare tree with compression of equal-level siblings.
        gsHTree<d,index_t> treeC;
        buildTree<d>(treeC, c.grid, ins, true);
        reportAimed<d>(c, label, "tree_c", treeC);

        // Bare tree without compression keeps the most non-sibling leaves.
        gsHTree<d,index_t> treeU;
        buildTree<d>(treeU, c.grid, ins, false);
        reportAimed<d>(c, label, "tree_u", treeU);

        // Hierarchical basis: level-0 grid of c.grid[0] elements per direction
        // (all aimed grids are cubes), refined through the public refineElements().
        gsKnotVector<> kv(0, 1, c.grid[0] - 1, 3);
        std::vector<gsKnotVector<> > kvs(d, kv);
        gsTensorBSplineBasis<d,real_t> tens(kvs);
        gsTHBSplineBasis<d,real_t> thb(tens);
        std::vector<index_t> flat;
        for (size_t i = 0; i != ins.size(); ++i)
        {
            flat.push_back(ins[i].lvl);
            flat.insert(flat.end(), ins[i].lo.begin(), ins[i].lo.end());
            flat.insert(flat.end(), ins[i].up.begin(), ins[i].up.end());
        }
        thb.refineElements(flat);
        reportAimed<d>(c, label, "thb", thb.tree());

        // Replica: the copy must reproduce the real getBoxes count on the tree's own
        // leaf order; on the reversed order the buggy and the corrected copy are compared.
        const std::vector<Box> list = leafList<d>(treeU);
        std::vector<Box> forward = list;
        connectBoxesAsInTree<d>(forward);
        const size_t realCount = runQuery<d>(treeU, 0).size();

        std::vector<Box> reversedBuggy(list.rbegin(), list.rend());
        std::vector<Box> reversedFixed = reversedBuggy;
        connectBoxesAsInTree<d>(reversedBuggy);
        connectBoxesFixed<d>(reversedFixed);
        gsInfo << "replica case=" << c.name << " order=" << label
               << " path=tree_u replicaAgrees=" << (forward.size() == realCount ? 1 : 0)
               << " reversed buggy=" << reversedBuggy.size()
               << " fixed=" << reversedFixed.size() << "\n";
    }
}

/// Receives each insertion list of a family.
typedef std::function<void(const std::vector<Insertion> &)> Visitor;
/// Calls a Visitor once for every member of a family.
typedef std::function<void(const Visitor &)> Generator;

/// Counters for one (compression flag, query) of a family.
struct FamilyStats
{
    FamilyStats() : cases(0), pairs(0), badCases(0) {}
    index_t cases;
    index_t pairs;
    index_t badCases;
    std::vector<std::string> offenders;
};

/// Run every member of a family on the bare tree, with and without compression, and
/// print one summary line per (compression flag, query).
template<short_t d> void
runFamily(const std::string & name, const std::vector<index_t> & grid,
          const Generator & generate, bool verbose)
{
    FamilyStats stats[2][2];
    gsHTree<d,index_t> tree;

    generate([&](const std::vector<Insertion> & ins)
    {
        for (index_t compressed = 0; compressed != 2; ++compressed)
        {
            buildTree<d>(tree, grid, ins, compressed != 0);
            for (index_t q = 0; q != 2; ++q)
            {
                const std::vector<Box> boxes = runQuery<d>(tree, q);
                const index_t p = leftoverPairs<d>(boxes);
                FamilyStats & st = stats[compressed][q];
                ++st.cases;
                st.pairs += p;
                if (p > 0)
                {
                    ++st.badCases;
                    reproduced = true;
                }
                if (verbose || (p > 0 && st.offenders.size() < 5))
                {
                    std::ostringstream os;
                    os << "  case ins=[" << formatInsertions(ins) << "] leftoverPairs=" << p
                       << " got=" << boxes.size() << "\n" << formatBoxes<d>(boxes);
                    if (verbose)
                        gsInfo << "family=" << name << " compressed=" << compressed
                               << " query=" << queryName(q) << "\n" << os.str();
                    else
                        st.offenders.push_back(os.str());
                }
            }
        }
    });

    for (index_t compressed = 0; compressed != 2; ++compressed)
    {
        for (index_t q = 0; q != 2; ++q)
        {
            const FamilyStats & st = stats[compressed][q];
            gsInfo << "family=" << name << " compressed=" << compressed
                   << " query=" << queryName(q) << " cases=" << st.cases
                   << " leftoverPairs=" << st.pairs << " expected=0 got=" << st.badCases << "\n";
            for (size_t i = 0; i != st.offenders.size(); ++i)
                gsInfo << st.offenders[i];
        }
    }
}

/// All axis-aligned rectangles of a level-0 grid with \a n elements per direction,
/// as level-1 insertions (corners in level-0 units, scaled to level 1).
std::vector<Insertion>
allRectangles(short_t dim, index_t n)
{
    std::vector<Insertion> rects;
    std::vector<index_t> lo(dim), up(dim);
    std::function<void(short_t)> rec = [&](short_t axis)
    {
        if (axis == dim)
        {
            rects.push_back(boxFromLevel0(1, lo, up));
            return;
        }
        for (index_t a = 0; a < n; ++a)
            for (index_t b = a + 1; b <= n; ++b)
            {
                lo[axis] = a;
                up[axis] = b;
                rec(axis + 1);
            }
    };
    rec(0);
    return rects;
}

/// All unit cells of the level-1 grid of a level-0 grid with \a n elements per direction.
std::vector<Insertion>
allLevel1Cells(short_t dim, index_t n)
{
    std::vector<Insertion> cells;
    std::vector<index_t> lo(dim), up(dim);
    std::function<void(short_t)> rec = [&](short_t axis)
    {
        if (axis == dim)
        {
            cells.push_back(boxAtLevel(1, lo, up));
            return;
        }
        for (index_t a = 0; a < 2 * n; ++a)
        {
            lo[axis] = a;
            up[axis] = a + 1;
            rec(axis + 1);
        }
    };
    rec(0);
    return cells;
}

/// Visit every ordered pair (including equal ones) of \a items.
Generator
orderedPairs(const std::vector<Insertion> & items)
{
    return [items](const Visitor & visit)
    {
        for (size_t a = 0; a != items.size(); ++a)
            for (size_t b = 0; b != items.size(); ++b)
            {
                std::vector<Insertion> ins;
                ins.push_back(items[a]);
                ins.push_back(items[b]);
                visit(ins);
            }
    };
}

/// Visit every set of 1, 2 or 3 distinct \a items in all insertion orders
/// (\a allOrders) or in ascending and reversed order only.
Generator
smallSets(const std::vector<Insertion> & items, bool allOrders)
{
    return [items, allOrders](const Visitor & visit)
    {
        const size_t n = items.size();
        if (allOrders)
        {
            for (size_t a = 0; a != n; ++a)
            {
                std::vector<Insertion> ins(1, items[a]);
                visit(ins);
                for (size_t b = 0; b != n; ++b)
                {
                    if (b == a)
                        continue;
                    ins.resize(1);
                    ins.push_back(items[b]);
                    visit(ins);
                    for (size_t c = 0; c != n; ++c)
                    {
                        if (c == a || c == b)
                            continue;
                        ins.resize(2);
                        ins.push_back(items[c]);
                        visit(ins);
                    }
                }
            }
            return;
        }

        for (size_t a = 0; a != n; ++a)
        {
            visit(std::vector<Insertion>(1, items[a]));
            for (size_t b = a + 1; b != n; ++b)
            {
                std::vector<Insertion> two;
                two.push_back(items[a]);
                two.push_back(items[b]);
                visit(two);
                visit(std::vector<Insertion>(two.rbegin(), two.rend()));
                for (size_t c = b + 1; c != n; ++c)
                {
                    std::vector<Insertion> three = two;
                    three.push_back(items[c]);
                    visit(three);
                    visit(std::vector<Insertion>(three.rbegin(), three.rend()));
                }
            }
        }
    };
}

/// Visit every rectangle (to level 1) followed by every unit cell of the level-2 grid.
Generator
rectangleThenLevel2Cell(const std::vector<Insertion> & rects, index_t level0Elements)
{
    return [rects, level0Elements](const Visitor & visit)
    {
        const index_t n2 = 4 * level0Elements;
        for (size_t a = 0; a != rects.size(); ++a)
            for (index_t i = 0; i != n2; ++i)
                for (index_t j = 0; j != n2; ++j)
                {
                    std::vector<Insertion> ins;
                    ins.push_back(rects[a]);
                    ins.push_back(boxAtLevel(2, std::vector<index_t>{i, j},
                                                std::vector<index_t>{i + 1, j + 1}));
                    visit(ins);
                }
    };
}

} // namespace

int main(int argc, char *argv[])
{
    // Options: the aimed cases always run; the exhaustive families can be skipped,
    // and every family member can be printed instead of the first few offenders.
    bool noFamilies = false;
    bool verbose = false;

    gsCmdLine cmd("Probe of gsHTree::connect_Boxes through getBoxes and getBoxesInLevelIndex.");
    cmd.addSwitch("noFamilies", "Run only the aimed cases, not the exhaustive families", noFamilies);
    cmd.addSwitch("verbose", "Print every member of every family", verbose);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // --------------- aimed cases, 2D: level-0 grid of 4x4 elements ---------------
    // Each case lists level-1 boxes in level-0 element units; the expected count is
    // derived by hand from the tree structure (level-1 pieces merge to one box,
    // the level-0 background splits along the first cut of the tree).
    const std::vector<index_t> grid4(2, 4);
    std::vector<AimedCase> cases2;
    {
        AimedCase c;
        c.grid = grid4;

        // Two level-1 halves of a strip, separated by the root cut at x = 2.
        c.name = "strip_x";
        c.base.clear();
        c.base.push_back(boxFromLevel0(1, std::vector<index_t>{0, 0}, std::vector<index_t>{2, 1}));
        c.base.push_back(boxFromLevel0(1, std::vector<index_t>{2, 0}, std::vector<index_t>{4, 1}));
        c.orders = ascDesc(c.base.size());
        c.expected = 2;
        cases2.push_back(c);

        // Two stacked level-1 cells in one column: level 1 merges, level 0 stays two boxes.
        c.name = "column_y_pair";
        c.base.clear();
        c.base.push_back(boxFromLevel0(1, std::vector<index_t>{0, 0}, std::vector<index_t>{1, 1}));
        c.base.push_back(boxFromLevel0(1, std::vector<index_t>{0, 1}, std::vector<index_t>{1, 2}));
        c.orders = ascDesc(c.base.size());
        c.expected = 3;
        cases2.push_back(c);

        // A full column of four level-1 cells, inserted bottom-up, top-down and interleaved.
        c.name = "column_y_full";
        c.base.clear();
        for (index_t j = 0; j != 4; ++j)
            c.base.push_back(boxFromLevel0(1, std::vector<index_t>{0, j}, std::vector<index_t>{1, j + 1}));
        c.orders = ascDescInterleaved();
        c.expected = 2;
        cases2.push_back(c);

        // A full row of four level-1 cells: a cascade of merges along x between non-siblings.
        c.name = "strip_x_cells";
        c.base.clear();
        for (index_t i = 0; i != 4; ++i)
            c.base.push_back(boxFromLevel0(1, std::vector<index_t>{i, 0}, std::vector<index_t>{i + 1, 1}));
        c.orders = ascDescInterleaved();
        c.expected = 2;
        cases2.push_back(c);
    }

    // --------------- aimed cases, 3D ---------------
    std::vector<AimedCase> cases3;
    {
        AimedCase c;
        c.grid = std::vector<index_t>(3, 4);

        // Two stacked level-1 cells in one z-column.
        c.name = "column_z_pair";
        c.base.clear();
        c.base.push_back(boxFromLevel0(1, std::vector<index_t>{0, 0, 0}, std::vector<index_t>{1, 1, 1}));
        c.base.push_back(boxFromLevel0(1, std::vector<index_t>{0, 0, 1}, std::vector<index_t>{1, 1, 2}));
        c.orders = ascDesc(c.base.size());
        c.expected = 4;
        cases3.push_back(c);

        // A full z-column of four level-1 cells.
        c.name = "column_z_full";
        c.base.clear();
        for (index_t k = 0; k != 4; ++k)
            c.base.push_back(boxFromLevel0(1, std::vector<index_t>{0, 0, k}, std::vector<index_t>{1, 1, k + 1}));
        c.orders = ascDesc(c.base.size());
        c.expected = 3;
        cases3.push_back(c);

        // The lower slab of a 2x2x2 grid, one level-1 cell pattern per (i,j).
        c.name = "slab_z_cells";
        c.grid = std::vector<index_t>(3, 2);
        c.base.clear();
        for (index_t i = 0; i != 2; ++i)
            for (index_t j = 0; j != 2; ++j)
                c.base.push_back(boxFromLevel0(1, std::vector<index_t>{i, j, 0}, std::vector<index_t>{i + 1, j + 1, 1}));
        c.orders = ascDesc(c.base.size());
        c.expected = 2;
        cases3.push_back(c);
    }

    for (size_t i = 0; i != cases2.size(); ++i)
        runAimed<2>(cases2[i]);
    for (size_t i = 0; i != cases3.size(); ++i)
        runAimed<3>(cases3[i]);

    // --------------- exhaustive families (bare tree only) ---------------
    // The aimed cases cannot cover every interaction of insertion order and
    // compression, so every pair, or every small set, of boxes on a small grid is
    // inserted and each returned list is checked for unmerged pairs.
    if (!noFamilies)
    {
        const std::vector<Insertion> rects2 = allRectangles(2, 4);
        const std::vector<Insertion> rects3 = allRectangles(3, 3);

        runFamily<2>("F2_rect_pairs", std::vector<index_t>(2, 4), orderedPairs(rects2), verbose);
        runFamily<2>("F2_cells3", std::vector<index_t>(2, 3), smallSets(allLevel1Cells(2, 3), true), verbose);
        runFamily<2>("F2_rect_then_l2", std::vector<index_t>(2, 4), rectangleThenLevel2Cell(rects2, 4), verbose);
        runFamily<3>("F3_rect_pairs", std::vector<index_t>(3, 3), orderedPairs(rects3), verbose);
        runFamily<3>("F3_cells3", std::vector<index_t>(3, 2), smallSets(allLevel1Cells(3, 2), false), verbose);
    }

    // --------------- verdict ---------------
    gsInfo << "probe-verdict: " << (reproduced ? "reproduced" : "not-reproduced") << "\n";
    return EXIT_SUCCESS;
}
