/** @file hierarchical_domain_partition.cpp

    @brief Collects the disjoint partition of a hierarchical domain into
    atomic cells (gsHTree leaves) and enumerates knot-span elements inside
    each cell. Validates the collection against the standard assembler loop.

    Each leaf omega_r of the gsHTree satisfies
        omega_r subset Omega^level \ Omega^{level+1},
    so the leaves form the disjoint partition
        Omega = union_l ( Omega^l \ Omega^{l+1} ).

    Usage:
        ./hierarchical_domain_partition [--plot]

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#include <gismo.h>

using namespace gismo;

// ============================================================
// Data structures
// ============================================================

/// One knot-span element [lower, upper] in parametric space.
struct ElementBox
{
    gsVector<real_t> lower;
    gsVector<real_t> upper;
};

/**
 * @brief One leaf of the gsHTree.
 *
 * Satisfies: omega subset Omega^level \ Omega^{level+1}.
 * Stores:
 *   - level         : hierarchy level of this leaf
 *   - lowerIndex    : lower index corner at level's own index system
 *   - upperIndex    : upper index corner at level's own index system
 *   - lowerParam    : parametric lower corner
 *   - upperParam    : parametric upper corner
 *   - elements      : all knot-span elements inside this leaf
 */
struct AtomicCell
{
    index_t            level;
    gsVector<index_t>  lowerIndex;
    gsVector<index_t>  upperIndex;
    gsVector<real_t>   lowerParam;
    gsVector<real_t>   upperParam;
    std::vector<ElementBox> elements;
};

// ============================================================
// Collect atomic cells
// ============================================================

/**
 * @brief Collect all gsHTree leaves as AtomicCells.
 *
 * For each leaf:
 *   1. Read leaf level and index corners (already in leaf-level index space).
 *   2. Convert index corners to parametric coords via the level's knot vector.
 *   3. Enumerate all knot-span elements by multi-index loop over alpha_d in
 *      [lowerIndex_d, upperIndex_d - 1].
 *
 * @note  Requires `!thb.manualLevels()`. For manually-leveled bases the leaf
 *        index corners need an additional dyadic→knot-index conversion (see
 *        gsHDomainIterator::updateLeaf) before calling kv(index).
 *
 * @tparam d  domain dimension
 * @tparam T  scalar type
 * @param  thb  THB- or HB-spline basis (standard uniform levels only)
 */
template<short_t d, class T>
std::vector<AtomicCell> collectAtomicCells(const gsHTensorBasis<d,T> & thb)
{
    GISMO_ASSERT( !thb.manualLevels(),
        "collectAtomicCells: manualLevels bases require a dyadic-to-knot-index "
        "conversion not implemented here." );

    std::vector<AtomicCell> cells;

    auto leaf = thb.tree().beginLeafIterator();

    for (; leaf.good(); leaf.next())
    {
        AtomicCell cell;
        cell.level      = leaf.level();
        cell.lowerIndex = leaf.lowerCorner(); // shifted to leaf-level index space
        cell.upperIndex = leaf.upperCorner();

        cell.lowerParam.resize(d);
        cell.upperParam.resize(d);

        // Parametric bounds from the level's knot vector
        for (short_t dd = 0; dd < d; ++dd)
        {
            const gsKnotVector<T> & kv =
                thb.tensorLevel(cell.level).component(dd).knots();

            cell.lowerParam[dd] = kv( cell.lowerIndex[dd] );
            cell.upperParam[dd] = kv( cell.upperIndex[dd] );
        }

        // Enumerate knot-span elements: alpha_d in [lowerIndex_d, upperIndex_d - 1]
        gsVector<index_t> alpha(d);
        alpha = cell.lowerIndex;

        bool more = true;
        while (more)
        {
            ElementBox elem;
            elem.lower.resize(d);
            elem.upper.resize(d);

            for (short_t dd = 0; dd < d; ++dd)
            {
                const gsKnotVector<T> & kv =
                    thb.tensorLevel(cell.level).component(dd).knots();

                elem.lower[dd] = kv( alpha[dd]     );
                elem.upper[dd] = kv( alpha[dd] + 1 );
            }

            cell.elements.push_back(elem);

            // Increment multi-index in lexicographic order
            more = false;
            for (short_t dd = d - 1; dd >= 0; --dd)
            {
                if (alpha[dd] < cell.upperIndex[dd] - 1)
                {
                    ++alpha[dd];
                    more = true;
                    break;
                }
                alpha[dd] = cell.lowerIndex[dd];
            }
        }

        cells.push_back(cell);
    }

    return cells;
}

// ============================================================
// main
// ============================================================

int main(int argc, char * argv[])
{
    bool plot = false;

    gsCmdLine cmd("Disjoint partition of a hierarchical domain into atomic cells.");
    cmd.addSwitch("plot", "Write ParaView files for the mesh", plot);
    try { cmd.getValues(argc, argv); } catch (int rv) { return rv; }

    // ------------------------------------------------------------------
    // Build a 2D THB-spline basis with L-shaped local refinement
    // (same pattern as hierarchical_quadrature_example.cpp)
    // ------------------------------------------------------------------
    gsKnotVector<> kv1(0, 1.0, 0, 3, 1);
    kv1.uniformRefine();
    kv1.uniformRefine();
    gsKnotVector<> kv2 = kv1;

    gsTensorBSplineBasis<2, real_t> tbsb(kv1, kv2);
    gsTHBSplineBasis<2, real_t> thb(tbsb);

    // L-shaped hierarchical refinement: 5 additional levels.
    // Bottom horizontal arm + left vertical arm at each level.
    {
        thb.uniformRefine(); // ensure level 1 exists
        for (index_t lvl = 1; lvl <= 5; ++lvl)
        {
            const index_t N = 4 * (1 << lvl);
            std::vector<index_t> boxes;
            // Horizontal arm: bottom 4 rows, full width
            boxes.insert(boxes.end(), {lvl, 0, 0, N, 4});
            // Vertical arm: left 4 columns, full height
            boxes.insert(boxes.end(), {lvl, 0, 0, 4, N});
            thb.refineElements(boxes);
        }
    }

    if (plot)
    {
        gsMesh<> mesh(thb);
        gsWriteParaview(mesh, "hierarchical_domain_mesh");
        gsInfo << "Mesh written to hierarchical_domain_mesh.vtp\n";
    }

    // ------------------------------------------------------------------
    // Step 1: collect atomic cells
    // ------------------------------------------------------------------
    const std::vector<AtomicCell> cells = collectAtomicCells<2, real_t>(thb);

    // ------------------------------------------------------------------
    // Step 2: print diagnostics
    // ------------------------------------------------------------------
    std::map<index_t, index_t> cellsPerLevel;
    std::map<index_t, index_t> elemsPerLevel;
    real_t totalArea = 0;
    index_t totalElems = 0;

    for (const AtomicCell & c : cells)
    {
        ++cellsPerLevel[c.level];
        elemsPerLevel[c.level] += (index_t)c.elements.size();
        totalElems              += (index_t)c.elements.size();

        for (const ElementBox & e : c.elements)
            totalArea += (e.upper[0] - e.lower[0]) * (e.upper[1] - e.lower[1]);
    }

    gsInfo << "Hierarchical domain partition diagnostics\n";
    gsInfo << "=========================================\n";
    gsInfo << "Total atomic cells : " << cells.size() << "\n";
    gsInfo << "Total elements     : " << totalElems   << "\n";

    for (const auto & p : cellsPerLevel)
        gsInfo << "  level " << p.first
               << " : atomic cells = "  << p.second
               << ", elements = "       << elemsPerLevel[p.first] << "\n";

    gsInfo << "Total parametric area : " << totalArea << "\n\n";

    // ------------------------------------------------------------------
    // Step 3: validate against domain()->allElements()
    // ------------------------------------------------------------------
    index_t refElems = 0;
    real_t  refArea  = 0;

    for (auto & elem : thb.domain()->allElements())
    {
        ++refElems;
        const gsVector<real_t> lo = elem.lowerCorner();
        const gsVector<real_t> hi = elem.upperCorner();
        refArea += (hi[0] - lo[0]) * (hi[1] - lo[1]);
    }

    gsInfo << "Reference (allElements)\n";
    gsInfo << "-----------------------\n";
    gsInfo << "  element count : " << refElems << "\n";
    gsInfo << "  total area    : " << refArea  << "\n\n";

    // ------------------------------------------------------------------
    // Assertions
    // ------------------------------------------------------------------
    GISMO_ENSURE( totalElems == refElems,
        "Element count mismatch: atomic cells give " << totalElems
        << " but allElements gives " << refElems );

    const real_t areaTol = 1e-10;
    GISMO_ENSURE( math::abs(totalArea - refArea) < areaTol,
        "Area mismatch: atomic cells give " << totalArea
        << " but allElements gives " << refArea );

    gsInfo << "Validation PASSED: element count and area match.\n";

    return 0;
}
