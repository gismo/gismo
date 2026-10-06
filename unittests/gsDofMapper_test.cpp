/** @file gsDofMapper_test.cpp

    @brief Tests for gismo::gsDofMapper (gsCore/gsDofMapper.h)

    The fixture tests (F1-F8) pin the complete observable state of a set of
    mappers as literal digests, so that a change of the internal storage of
    gsDofMapper cannot silently move an expectation along with the
    implementation.  The tests after them check individual properties, some
    against an oracle computed through the public interface.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
**/

#include "gismo_unittest.h"
#include <gsAssembler/gsDofMapperCreator.h>

#include <limits>
#include <set>

using namespace gismo;

namespace {

// Pins one observable string.  The tag is prefixed to both sides so that a
// failing assertion names the pinned area without any extra plumbing.
#define PIN(tag, expected, actual) \
    CHECK_EQUAL(std::string(tag "|") + (expected), std::string(tag "|") + (actual))

index_t nPatches(const gsDofMapper & m)
{ return static_cast<index_t>(m.numPatches()); }

std::string join(const gsVector<index_t> & v)
{
    std::ostringstream os;
    os << "[";
    for (index_t i = 0; i != v.size(); ++i)
        os << (i ? "," : "") << v[i];
    os << "]";
    return os.str();
}

std::string join(const std::vector<index_t> & v)
{
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i != v.size(); ++i)
        os << (i ? "," : "") << v[i];
    os << "]";
    return os.str();
}

std::string join(const std::vector<std::pair<index_t,index_t> > & v)
{
    std::ostringstream os;
    os << "[";
    for (size_t i = 0; i != v.size(); ++i)
        os << (i ? "," : "") << "(" << v[i].first << "," << v[i].second << ")";
    os << "]";
    return os.str();
}

// --- scalar/aggregate counters -------------------------------------------

// componentsSize, numComponents, numPatches, mapSize, size, freeSize,
// boundarySize, boundarySizeWithDuplicates, coupledSize, taggedSize,
// allFree, isPermutation, isFinalized.
std::string dumpCounts(const gsDofMapper & m)
{
    std::ostringstream os;
    os << "comps="   << m.componentsSize()
       << " ncomp="  << m.numComponents()
       << " npatch=" << m.numPatches()
       << " map="    << m.mapSize()
       << " size="   << m.size()
       << " free="   << m.freeSize()
       << " elim="   << m.boundarySize()
       << " cpld="   << m.coupledSize()
       << " tagged=" << m.taggedSize()
       << " allfree="<< (m.allFree() ? 1 : 0)
       << " perm="   << (m.isPermutation() ? 1 : 0)
       << " final="  << (m.isFinalized() ? 1 : 0);
    return os.str();
}

// boundarySizeWithDuplicates(), pinned on its own: it counts the stored
// entries above the global free count, i.e. the eliminated dofs with their
// multiplicity over patches and components.
std::string dumpElimDup(const gsDofMapper & m)
{
    std::ostringstream os;
    os << "elimdup=" << m.boundarySizeWithDuplicates();
    return os.str();
}

// Per-component size(c), freeSize(c), totalSize(c).
std::string dumpPerComponent(const gsDofMapper & m)
{
    std::ostringstream os;
    for (index_t c = 0; c != m.numComponents(); ++c)
        os << (c ? " " : "") << "c" << c
           << ":size=" << m.size(c)
           << ",free=" << m.freeSize(c)
           << ",total=" << m.totalSize(c);
    return os.str();
}

// offset(k) for every patch and patchSize(k,c) for every (patch,component).
std::string dumpLayout(const gsDofMapper & m)
{
    std::ostringstream os;
    os << "off=[";
    for (index_t k = 0; k != nPatches(m); ++k)
        os << (k ? "," : "") << m.offset(k);
    os << "] ps=";
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        os << (c ? " " : "") << "c" << c << ":[";
        for (index_t k = 0; k != nPatches(m); ++k)
            os << (k ? "," : "") << m.patchSize(k, c);
        os << "]";
    }
    return os.str();
}

// firstFreeIndex(0) and lastIndex().  See
// first_free_index_bounds_the_free_block_of_each_component for c>=1.
std::string dumpFirstLast(const gsDofMapper & m)
{
    std::ostringstream os;
    os << "first0=" << m.firstFreeIndex(0) << " last=" << m.lastIndex();
    return os.str();
}

// --- index maps -----------------------------------------------------------

std::string dumpAsVector(const gsDofMapper & m)
{
    std::ostringstream os;
    for (index_t c = 0; c != m.numComponents(); ++c)
        os << (c ? " " : "") << "c" << c << join(m.asVector(c));
    return os.str();
}

// mapIndex(n) over the full flat range [0, mapSize()).
std::string dumpMapIndex(const gsDofMapper & m)
{
    std::ostringstream os;
    os << "[";
    for (size_t n = 0; n != m.mapSize(); ++n)
        os << (n ? "," : "") << m.mapIndex(static_cast<index_t>(n));
    os << "]";
    return os.str();
}

// index(i,k,c) for every local dof of every patch and component.
std::string dumpIndex(const gsDofMapper & m)
{
    std::ostringstream os;
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        os << (c ? ";" : "") << "c" << c;
        for (index_t k = 0; k != nPatches(m); ++k)
        {
            os << "|";
            const index_t n = static_cast<index_t>(m.patchSize(k, c));
            for (index_t i = 0; i != n; ++i)
                os << (i ? "," : "") << m.index(i, k, c);
        }
    }
    return os.str();
}

// freeIndex(i,k,c) for every local dof.  Unlike index(), freeIndex does not
// add m_shift; on a shifted mapper this is the discriminator between the two.
std::string dumpFreeIndex(const gsDofMapper & m)
{
    std::ostringstream os;
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        os << (c ? ";" : "") << "c" << c;
        for (index_t k = 0; k != nPatches(m); ++k)
        {
            os << "|";
            const index_t n = static_cast<index_t>(m.patchSize(k, c));
            for (index_t i = 0; i != n; ++i)
                os << (i ? "," : "") << m.freeIndex(i, k, c);
        }
    }
    return os.str();
}

// localToGlobal(locals,k,globals,c) on the full local range of every patch.
std::string dumpLocalToGlobal(const gsDofMapper & m, const index_t c)
{
    std::ostringstream os;
    for (index_t k = 0; k != nPatches(m); ++k)
    {
        const index_t n = static_cast<index_t>(m.patchSize(k, c));
        gsMatrix<index_t> locals(n, 1), globals;
        for (index_t i = 0; i != n; ++i)
            locals(i, 0) = i;
        m.localToGlobal(locals, k, globals, c);
        os << (k ? " " : "") << "k" << k << ":"
           << globals.rows() << "x" << globals.cols() << "[";
        for (index_t i = 0; i != globals.rows(); ++i)
            os << (i ? "," : "") << globals(i, 0);
        os << "]";
    }
    return os.str();
}

// localToGlobal2(locals,k,globals,numFree,c) on the full local range of every
// patch: the returned numFree plus the full two-column result.  Column 0 is
// the position in locals and column 1 the global index; the free entries are
// packed from the top and the eliminated ones from the bottom.
std::string dumpLocalToGlobal2(const gsDofMapper & m, const index_t c)
{
    std::ostringstream os;
    for (index_t k = 0; k != nPatches(m); ++k)
    {
        const index_t n = static_cast<index_t>(m.patchSize(k, c));
        gsMatrix<index_t> locals(n, 1), globals;
        for (index_t i = 0; i != n; ++i)
            locals(i, 0) = i;
        index_t numFree = -1;
        m.localToGlobal2(locals, k, globals, numFree, c);
        os << (k ? " " : "") << "k" << k << ":nfree=" << numFree << " "
           << globals.rows() << "x" << globals.cols() << "[";
        for (index_t i = 0; i != globals.rows(); ++i)
            os << (i ? "," : "") << "(" << globals(i, 0) << "," << globals(i, 1) << ")";
        os << "]";
    }
    return os.str();
}

// is_free(i,k,c) / is_boundary(i,k,c) as one character per local dof.
// 'F' free, 'B' eliminated, '?' if the two predicates disagree.
std::string dumpFreeBoundaryFlags(const gsDofMapper & m)
{
    std::ostringstream os;
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        os << (c ? ";" : "");
        for (index_t k = 0; k != nPatches(m); ++k)
        {
            if (k) os << "|";
            const index_t n = static_cast<index_t>(m.patchSize(k, c));
            for (index_t i = 0; i != n; ++i)
            {
                const bool f = m.is_free    (i, k, c);
                const bool b = m.is_boundary(i, k, c);
                os << (f == b ? '?' : (f ? 'F' : 'B'));
            }
        }
    }
    return os.str();
}

// bindex(i,k,c) for the local dofs that are eliminated.
std::string dumpBindex(const gsDofMapper & m)
{
    std::ostringstream os;
    os << "[";
    bool first = true;
    for (index_t c = 0; c != m.numComponents(); ++c)
        for (index_t k = 0; k != nPatches(m); ++k)
        {
            const index_t n = static_cast<index_t>(m.patchSize(k, c));
            for (index_t i = 0; i != n; ++i)
                if (m.is_boundary(i, k, c))
                {
                    os << (first ? "" : ",") << c << ":" << k << ":" << i
                       << "->" << m.bindex(i, k, c);
                    first = false;
                }
        }
    os << "]";
    return os.str();
}

// is_free_index / is_boundary_index over the shifted global range, plus
// global_to_bindex on the eliminated ones.
std::string dumpGlobalFlags(const gsDofMapper & m, const index_t shift)
{
    std::ostringstream os;
    for (index_t gl = shift; gl != shift + m.size(); ++gl)
    {
        const bool f = m.is_free_index    (gl);
        const bool b = m.is_boundary_index(gl);
        os << (f == b ? '?' : (f ? 'F' : 'B'));
    }
    os << " gbi=[";
    bool first = true;
    for (index_t gl = shift; gl != shift + m.size(); ++gl)
        if (m.is_boundary_index(gl))
        {
            os << (first ? "" : ",") << gl << "->" << m.global_to_bindex(gl);
            first = false;
        }
    os << "]";
    return os.str();
}

// findBoundary(k,c) and findFree(k,c) for every patch.
std::string dumpFindBoundaryFree(const gsDofMapper & m, const index_t c = 0)
{
    std::ostringstream os;
    os << "bnd=";
    for (index_t k = 0; k != nPatches(m); ++k)
        os << (k ? "," : "") << join(m.findBoundary(k, c));
    os << " free=";
    for (index_t k = 0; k != nPatches(m); ++k)
        os << (k ? "," : "") << join(m.findFree(k, c));
    return os.str();
}

// The digests below that walk the global index range take the mapper's
// shift and query the shifted indices shift+g, g in [0,size()), but print g.
// A shifted mapper and its unshifted original therefore produce the same
// string, which is the whole contract of the shift: it relabels the global
// indices and changes nothing else.

// indexOnPatch(gl,k,local) over every global index and every patch.
std::string dumpIndexOnPatch(const gsDofMapper & m, const index_t shift = 0)
{
    std::ostringstream os;
    for (index_t gl = 0; gl != m.size(); ++gl)
    {
        os << (gl ? " " : "") << gl << ":";
        for (index_t k = 0; k != nPatches(m); ++k)
        {
            index_t local = -12345;
            if (m.indexOnPatch(shift + gl, k, local)) os << local;
            else                              os << "-";
            if (k + 1 != nPatches(m)) os << "/";
        }
    }
    return os.str();
}

// --- coupled queries ------------------------------------------------------

// is_coupled / is_coupled_index / cindex / findCoupled / findFreeUncoupled
// for component \a c.  The gflags column (is_coupled_index over the whole
// global range) is component-independent and therefore repeats verbatim for
// every c; it is kept in so that each component's literal stands on its own.
std::string dumpCoupledQueries(const gsDofMapper & m, const index_t c = 0,
                               const index_t shift = 0)
{
    std::ostringstream os;
    os << "flags=";
    for (index_t k = 0; k != nPatches(m); ++k)
    {
        if (k) os << "|";
        const index_t n = static_cast<index_t>(m.patchSize(k, c));
        for (index_t i = 0; i != n; ++i)
            os << (m.is_coupled(i, k, c) ? 'C' : '.');
    }
    os << " cidx=[";
    bool first = true;
    for (index_t k = 0; k != nPatches(m); ++k)
    {
        const index_t n = static_cast<index_t>(m.patchSize(k, c));
        for (index_t i = 0; i != n; ++i)
            if (m.is_coupled(i, k, c))
            {
                os << (first ? "" : ",") << k << ":" << i
                   << "->" << m.cindex(i, k, c);
                first = false;
            }
    }
    os << "] gflags=";
    for (index_t gl = 0; gl != m.size(); ++gl)
        os << (m.is_coupled_index(shift + gl) ? 'C' : '.');
    os << " findCoupled=";
    for (index_t k = 0; k != nPatches(m); ++k)
        os << (k ? "," : "") << join(m.findCoupled(k, -1, c));
    os << " findCoupledPairs=";
    for (index_t k = 0; k != nPatches(m); ++k)
        for (index_t j = 0; j != nPatches(m); ++j)
            os << ((k || j) ? "," : "") << k << j << join(m.findCoupled(k, j, c));
    os << " findFreeUncoupled=";
    for (index_t k = 0; k != nPatches(m); ++k)
        os << (k ? "," : "") << join(m.findFreeUncoupled(k, c));
    return os.str();
}

// inverseOnPatch(k) as a flat string.
std::string dumpInverseOnPatch(const gsDofMapper & m, const index_t k)
{
    const std::map<index_t,index_t> inv = m.inverseOnPatch(k);
    std::ostringstream os;
    os << "[";
    for (std::map<index_t,index_t>::const_iterator it = inv.begin();
         it != inv.end(); ++it)
        os << (it == inv.begin() ? "" : ",") << "(" << it->first << "," << it->second << ")";
    os << "]";
    return os.str();
}

// --- component / pre-image queries ---------------------------------------

std::string dumpComponentOf(const gsDofMapper & m, const index_t shift = 0)
{
    std::ostringstream os;
    os << "[";
    for (index_t gl = 0; gl != m.size(); ++gl)
        os << (gl ? "," : "") << m.componentOf(shift + gl);
    os << "]";
    return os.str();
}

std::string dumpPreImages(const gsDofMapper & m, const index_t shift = 0)
{
    std::ostringstream os;
    std::vector<std::pair<index_t,index_t> > pre;
    for (index_t gl = 0; gl != m.size(); ++gl)
    {
        m.preImage(shift + gl, pre);
        os << (gl ? " " : "") << gl << join(pre);
    }
    return os.str();
}

std::string dumpAnyPreImage(const gsDofMapper & m, const index_t shift = 0)
{
    std::ostringstream os;
    os << "[";
    for (index_t gl = 0; gl != m.size(); ++gl)
    {
        const std::pair<index_t,index_t> p = m.anyPreImage(shift + gl);
        os << (gl ? "," : "") << "(" << p.first << "," << p.second << ")";
    }
    os << "]";
    return os.str();
}

// --- tagged queries ------------------------------------------------------

// getTagged() is printed as stored, i.e. unshifted; see getTagged().
std::string dumpTagged(const gsDofMapper & m, const index_t shift = 0)
{
    std::ostringstream os;
    os << "tagged=" << join(m.getTagged())
       << " n=" << m.taggedSize()
       << " gflags=";
    for (index_t gl = 0; gl != m.size(); ++gl)
        os << (m.is_tagged_index(shift + gl) ? 'T' : '.');
    os << " flags=";
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        os << (c ? ";" : "");
        for (index_t k = 0; k != nPatches(m); ++k)
        {
            if (k) os << "|";
            const index_t n = static_cast<index_t>(m.patchSize(k, c));
            for (index_t i = 0; i != n; ++i)
                os << (m.is_tagged(i, k, c) ? 'T' : '.');
        }
    }
    os << " tindex=";
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        os << (c ? ";" : "");
        for (index_t k = 0; k != nPatches(m); ++k)
        {
            if (k) os << "|";
            const index_t n = static_cast<index_t>(m.patchSize(k, c));
            for (index_t i = 0; i != n; ++i)
                os << (i ? "," : "") << m.tindex(i, k, c);
        }
    }
    return os.str();
}

// --- fixtures -------------------------------------------------------------

// F1: one patch with six local dofs, one component, nothing coupled and
// nothing eliminated.  The simplest possible finalized mapper.
gsDofMapper singlePatchPlain()
{
    gsVector<index_t> sz(1);
    sz[0] = 6;
    gsDofMapper m(sz, 1);
    m.finalize();
    return m;
}

// F2: two patches of six local dofs each, one component.  Three matchDof
// calls build two coupling groups (the second one holds three local dofs),
// two plain free dofs are eliminated, and two dofs are tagged after
// finalize().
gsDofMapper twoPatchCoupledElim()
{
    gsVector<index_t> sz(2);
    sz[0] = 6; sz[1] = 6;
    gsDofMapper m(sz, 1);
    m.matchDof(0, 4, 1, 0);   // coupling group A: (0,4) ~ (1,0)
    m.matchDof(0, 5, 1, 1);   // coupling group B: (0,5) ~ (1,1)
    m.matchDof(1, 2, 0, 5);   // extends group B by (1,2)
    m.eliminateDof(0, 0);     // plain free dof of patch 0
    m.eliminateDof(3, 1);     // plain free dof of patch 1
    m.finalize();
    m.markTagged(1, 0);
    m.markTagged(4, 1);
    return m;
}

// F3: two patches of four and five local dofs, three components.  Each
// component gets a different matching and elimination pattern, so the three
// components end up with different free/coupled/eliminated counts.  This is
// what exercises the per-component prefix arrays and the multi-component
// re-offset pass at the end of finalize().
gsDofMapper threeCompUniform()
{
    gsVector<index_t> sz(2);
    sz[0] = 4; sz[1] = 5;
    gsDofMapper m(sz, 3);

    // component 0: one coupling pair, one eliminated dof
    m.matchDof(0, 3, 1, 0, 0);
    m.eliminateDof(0, 0, 0);

    // component 1: two coupling pairs, two eliminated dofs
    m.matchDof(0, 2, 1, 1, 1);
    m.matchDof(0, 3, 1, 2, 1);
    m.eliminateDof(4, 1, 1);
    m.eliminateDof(0, 0, 1);

    // component 2: no coupling, one eliminated dof
    m.eliminateDof(3, 1, 2);

    m.finalize();
    return m;
}

// F4: F3 with a nonzero global shift and a nonzero boundary shift.
gsDofMapper threeCompShifted()
{
    gsDofMapper m = threeCompUniform();
    m.setShift(100);
    m.setBoundaryShift(7);
    return m;
}

// F5: the identity/aliased layout: three patches, seven (already global)
// dofs, two components.
gsDofMapper identityMapper()
{
    gsDofMapper m;
    m.setIdentity(3, 7, 2);
    m.finalize();
    return m;
}

// F6: F3, then all coupled dofs are tagged and the free block of component 0
// is permuted by an explicit non-trivial permutation of length freeSize(0).
gsDofMapper permutedMapper()
{
    gsDofMapper m = threeCompUniform();
    m.markCoupledAsTagged();
    gsVector<index_t> perm(7);      // freeSize(0) == 7 on F3
    perm << 3, 0, 6, 1, 5, 2, 4;
    m.permuteFreeDofs(perm, 0);
    return m;
}

// F7: two patches of five local dofs, one component, and a single
// collapseDofs call gluing the local dofs 0, 2 and 4 of patch 0 together.
gsDofMapper collapsedMapper()
{
    gsVector<index_t> sz(2);
    sz[0] = 5; sz[1] = 5;
    gsDofMapper m(sz, 1);
    gsMatrix<unsigned> b(3, 1);
    b(0, 0) = 0; b(1, 0) = 2; b(2, 0) = 4;
    m.colapseDofs(0, b);
    m.finalize();
    return m;
}

// F8: a real basis fixture.  Two unit squares side by side sharing one
// interface, biquadratic, uniformly refined twice, run through the
// gsDofMapperCreator conforming path end to end.
gsMultiBasis<real_t> twoPatchBasis()
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
    gsMultiBasis<real_t> mb(mp);
    mb.degreeElevate(1);
    mb.uniformRefine(2);
    return mb;
}

gsDofMapper creatorTwoPatch()
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    return createMapper(mb, /*nComp=*/1, /*conforming=*/true, /*finalize=*/true);
}

// F9: a genuinely ragged patch-concatenated mapper, built through the
// per-component ragged constructor.  Component 0: patches of size 5,5
// (total 10).  Component 1: patches of size 2,4 (total 6) -- UNEQUAL to
// component 0's, on every patch.  hasDistinctComponentSpaces is declared
// false: the point of this fixture is the ragged per-component *sizes*,
// independent of the declared-distinct flag.
//
// raggedSetup() is F9 before finalize().  Local dof 2 of patch 0 in
// component 1 is one past that patch's end, and the flat storage slot it
// would address is local dof 0 of patch 1.
gsDofMapper raggedSetup()
{
    std::vector<gsVector<index_t> > sz(2);
    sz[0].resize(2); sz[0][0] = 5; sz[0][1] = 5;
    sz[1].resize(2); sz[1][0] = 2; sz[1][1] = 4;
    return gsDofMapper(sz, false);
}

gsDofMapper raggedPatchMapper()
{
    gsDofMapper m = raggedSetup();
    m.finalize();
    return m;
}

// F10: a genuinely ragged global-identity/aliased mapper, built through the
// per-component setIdentity overload.  Component 0 total 7, component 1
// total 10 -- UNEQUAL, same motivation as F9 for the aliased layout.
gsDofMapper raggedIdentityMapper()
{
    std::vector<size_t> dofsPerComponent;
    dofsPerComponent.push_back(7);
    dofsPerComponent.push_back(10);
    gsDofMapper m;
    m.setIdentity(3, dofsPerComponent);
    m.finalize();
    return m;
}

// --- the target configuration: the 2D Raviart-Thomas component pair ------
//
// One vector variable whose two components live in different tensor-product
// spline spaces on the same patch and the same mesh, degrees and
// regularities transposed between them:
//
//     component 0 = S^{3,2}_{2,1}      component 1 = S^{2,3}_{1,2}
//
// A 1D space of degree p and regularity r over E elements has dimension
// E(p-r)+r+1, so with p-r==1 in every direction here the component
// dimensions are (E0+3)(E1+2) and (E0+2)(E1+3).  On an ISOTROPIC mesh those
// two numbers coincide (4x4: 42 and 42) and only an anisotropic mesh
// separates them (4x2: 28 and 30).  That is why every fixture below comes in
// both mesh flavours, and why the distinctness of the component spaces has
// to be a declared property: on the isotropic mesh -- the default case, a
// unit square uniformly refined -- an RT mapper is size-indistinguishable
// from an ordinary uniform one, so no predicate built on cardinalities can
// tell them apart.

gsTensorBSplineBasis<2,real_t> rtComponentBasis(short_t p0, short_t p1,
                                                index_t e0, index_t e1)
{
    // Maximal regularity (interior multiplicity p-r == 1) in both directions.
    gsKnotVector<real_t> kv0(0.0, 1.0, e0-1, p0+1, 1, p0);
    gsKnotVector<real_t> kv1(0.0, 1.0, e1-1, p1+1, 1, p1);
    return gsTensorBSplineBasis<2,real_t>(kv0, kv1);
}

// One patch, one size vector per component, sizes taken from the two real
// bases so the dimension arithmetic above is proven rather than asserted.
// hasDistinctComponentSpaces is declared true: the components genuinely come
// from two different basis objects.
std::vector<gsVector<index_t> > rtSizes(index_t e0, index_t e1)
{
    const gsTensorBSplineBasis<2,real_t> b0 = rtComponentBasis(3, 2, e0, e1);
    const gsTensorBSplineBasis<2,real_t> b1 = rtComponentBasis(2, 3, e0, e1);

    std::vector<gsVector<index_t> > sz(2);
    sz[0].resize(1); sz[0][0] = b0.size();
    sz[1].resize(1); sz[1][0] = b1.size();
    return sz;
}

gsDofMapper rtMapper(index_t e0, index_t e1, bool declareDistinct = true)
{
    gsDofMapper m(rtSizes(e0, e1), declareDistinct);
    m.finalize();
    return m;
}

// F11: the RT pair on an ISOTROPIC 4x4 mesh -- 42 and 42, equal.
gsDofMapper rtIsotropicMapper() { return rtMapper(4, 4); }

// F12: the RT pair on an ANISOTROPIC 4x2 mesh -- 28 and 30, unequal.
gsDofMapper rtAnisotropicMapper() { return rtMapper(4, 2); }

} // anonymous namespace


SUITE(gsDofMapper_test)
{

// =========================================================================
// Default-constructed state
// =========================================================================

// The default-constructed observation: zero components, mapSize()==0 and,
// despite there being no patch data at all, numPatches()==1, the patch count
// the default constructor sets.
TEST(default_constructed_state)
{
    const gsDofMapper m;

    CHECK_EQUAL(0, m.numComponents());
    CHECK_EQUAL(0u, (unsigned)m.componentsSize());
    CHECK_EQUAL(0u, (unsigned)m.mapSize());
    CHECK_EQUAL(1u, (unsigned)m.numPatches());
    CHECK(!m.isFinalized());

    // size(), boundarySize() and coupledSize() guard with GISMO_ENSURE, which,
    // unlike GISMO_ASSERT (compiled out under NDEBUG), throws
    // std::runtime_error in Release builds too.  Hence plain CHECK_THROW and
    // not CHECK_THROW_IN_DEBUG.  The three "Ensure" banners these print on
    // std::cerr are expected output of this test.
    CHECK_THROW(m.size()        , std::runtime_error);
    CHECK_THROW(m.boundarySize(), std::runtime_error);
    CHECK_THROW(m.coupledSize() , std::runtime_error);
}

// =========================================================================
// F1 -- single patch, no coupling, no elimination
// =========================================================================

TEST(single_patch_plain)
{
    const gsDofMapper m = singlePatchPlain();
    PIN("f1.elimdup", "elimdup=0", dumpElimDup(m));
    PIN("f1.counts", "comps=1 ncomp=1 npatch=1 map=6 size=6 free=6 elim=0 cpld=0 tagged=0 allfree=1 perm=1 final=1", dumpCounts(m));
    PIN("f1.percomp", "c0:size=6,free=6,total=6", dumpPerComponent(m));
    PIN("f1.layout", "off=[0] ps=c0:[6]", dumpLayout(m));
    PIN("f1.firstlast", "first0=0 last=6", dumpFirstLast(m));

    PIN("f1.asvector", "c0[0,1,2,3,4,5]", dumpAsVector(m));
    PIN("f1.mapindex", "[0,1,2,3,4,5]", dumpMapIndex(m));
    PIN("f1.index", "c0|0,1,2,3,4,5", dumpIndex(m));
    PIN("f1.flags", "FFFFFF", dumpFreeBoundaryFlags(m));
    PIN("f1.bindex", "[]", dumpBindex(m));
    PIN("f1.gflags", "FFFFFF gbi=[]", dumpGlobalFlags(m, 0));
    PIN("f1.findbf", "bnd=[] free=[0,1,2,3,4,5]", dumpFindBoundaryFree(m));
    PIN("f1.onpatch", "0:0 1:1 2:2 3:3 4:4 5:5", dumpIndexOnPatch(m));

    // See three_comp_coupled_queries_per_component for the same digest on a
    // mapper with more than one component.
    PIN("f1.coupled", "flags=...... cidx=[] gflags=...... findCoupled=[] findCoupledPairs=00[] findFreeUncoupled=[0,1,2,3,4,5]", dumpCoupledQueries(m));

    PIN("f1.componentof", "[0,0,0,0,0,0]", dumpComponentOf(m));
    PIN("f1.preimage", "0[(0,0)] 1[(0,1)] 2[(0,2)] 3[(0,3)] 4[(0,4)] 5[(0,5)]", dumpPreImages(m));
    PIN("f1.anypreimage", "[(0,0),(0,1),(0,2),(0,3),(0,4),(0,5)]", dumpAnyPreImage(m));
    PIN("f1.anypreimages", "[(0,0),(0,1),(0,2),(0,3),(0,4),(0,5)]", join(m.anyPreImages(0)));

    PIN("f1.inverseonpatch", "[(0,0),(1,1),(2,2),(3,3),(4,4),(5,5)]",
        dumpInverseOnPatch(m, 0));

    // A single-component permutation: every position of the result is the
    // image of a local dof, so no sentinel is visible here.  See
    // inverse_as_vector_marks_other_components for the sentinel.
    PIN("f1.inverseasvector", "[0,1,2,3,4,5]", join(m.inverseAsVector(0)));
}

// =========================================================================
// F2 -- two patches, coupling, elimination, tagging
// =========================================================================

TEST(two_patch_coupled_elim)
{
    const gsDofMapper m = twoPatchCoupledElim();
    PIN("f2.elimdup", "elimdup=2", dumpElimDup(m));
    PIN("f2.counts", "comps=1 ncomp=1 npatch=2 map=12 size=9 free=7 elim=2 cpld=2 tagged=2 allfree=0 perm=0 final=1", dumpCounts(m));
    PIN("f2.percomp", "c0:size=9,free=7,total=12", dumpPerComponent(m));
    PIN("f2.layout", "off=[0,6] ps=c0:[6,6]", dumpLayout(m));
    PIN("f2.firstlast", "first0=0 last=7", dumpFirstLast(m));

    PIN("f2.asvector", "c0[7,0,1,2,5,6,5,6,6,8,3,4]", dumpAsVector(m));
    PIN("f2.mapindex", "[7,0,1,2,5,6,5,6,6,8,3,4]", dumpMapIndex(m));
    PIN("f2.index", "c0|7,0,1,2,5,6|5,6,6,8,3,4", dumpIndex(m));
    PIN("f2.flags", "BFFFFF|FFFBFF", dumpFreeBoundaryFlags(m));
    PIN("f2.bindex", "[0:0:0->0,0:1:3->1]", dumpBindex(m));
    PIN("f2.gflags", "FFFFFFFBB gbi=[7->0,8->1]", dumpGlobalFlags(m, 0));
    PIN("f2.findbf", "bnd=[0],[3] free=[1,2,3,4,5],[0,1,2,4,5]", dumpFindBoundaryFree(m));
    PIN("f2.onpatch", "0:1/- 1:2/- 2:3/- 3:-/4 4:-/5 5:4/0 6:5/1 7:0/- 8:-/3", dumpIndexOnPatch(m));

    PIN("f2.coupled", "flags=....CC|CCC... cidx=[0:4->0,0:5->1,1:0->0,1:1->1,1:2->1] gflags=.....CC.. findCoupled=[4,5],[0,1,2] findCoupledPairs=00[],01[4,5],10[0,1,2],11[] findFreeUncoupled=[1,2,3],[4,5]", dumpCoupledQueries(m));

    PIN("f2.componentof", "[0,0,0,0,0,0,0,0,0]", dumpComponentOf(m));
    PIN("f2.preimage", "0[(0,1)] 1[(0,2)] 2[(0,3)] 3[(1,4)] 4[(1,5)] 5[(0,4),(1,0)] 6[(0,5),(1,1),(1,2)] 7[(0,0)] 8[(1,3)]", dumpPreImages(m));
    PIN("f2.anypreimage", "[(0,1),(0,2),(0,3),(1,4),(1,5),(0,4),(0,5),(0,0),(1,3)]", dumpAnyPreImage(m));
    // One entry per global index, size() == 9, not one per stored local dof
    // (12 here).
    PIN("f2.anypreimages", "[(0,1),(0,2),(0,3),(1,4),(1,5),(0,4),(0,5),(0,0),(1,3)]", join(m.anyPreImages(0)));

    PIN("f2.tagged", "tagged=[0,3] n=2 gflags=T..T..... flags=.T....|....T. tindex=2,0,1,1,2,2|2,2,2,2,1,2", dumpTagged(m));

    // localToGlobal and localToGlobal2 are the assembly entry points and they
    // are pinned in their own right, not only transitively through index().
    // F2 carries both free and eliminated dofs, so the free-from-the-top /
    // boundary-from-the-bottom packing of localToGlobal2 is exercised.
    PIN("f2.l2g", "k0:6x1[7,0,1,2,5,6] k1:6x1[5,6,6,8,3,4]", dumpLocalToGlobal(m, 0));
    PIN("f2.l2g2", "k0:nfree=5 6x2[(1,0),(2,1),(3,2),(4,5),(5,6),(0,7)] k1:nfree=5 6x2[(0,5),(1,6),(2,6),(4,3),(5,4),(3,8)]", dumpLocalToGlobal2(m, 0));
}

// Known defect, pinned so that it does not change by accident:
// gsDofMapper::findTagged builds the intersection into a local std::list and
// then returns an untouched default-constructed gsVector, so it is always
// empty even when the patch carries tagged dofs.  Fixing findTagged must
// replace this test with one of the tagged local dofs.
TEST(find_tagged_returns_empty_defect)
{
    const gsDofMapper m = twoPatchCoupledElim();
    CHECK(m.taggedSize() > 0);            // the fixture really has tagged dofs
    CHECK_EQUAL(0, m.findTagged(0).size());
    CHECK_EQUAL(0, m.findTagged(1).size());
}

// =========================================================================
// F3 -- three components with different per-component patterns
// =========================================================================

TEST(three_comp_uniform)
{
    const gsDofMapper m = threeCompUniform();
    // boundarySizeWithDuplicates counts the eliminated dofs with their
    // multiplicity over patches, so it is bounded below by boundarySize().
    PIN("f3.elimdup", "elimdup=4", dumpElimDup(m));
    PIN("f3.counts", "comps=3 ncomp=3 npatch=2 map=27 size=24 free=20 elim=4 cpld=3 tagged=0 allfree=0 perm=0 final=1", dumpCounts(m));
    PIN("f3.percomp", "c0:size=8,free=7,total=9 c1:size=7,free=5,total=9 c2:size=9,free=8,total=9", dumpPerComponent(m));
    PIN("f3.layout", "off=[0,4] ps=c0:[4,5] c1:[4,5] c2:[4,5]", dumpLayout(m));
    PIN("f3.firstlast", "first0=0 last=20", dumpFirstLast(m));

    PIN("f3.asvector", "c0[20,0,1,6,6,2,3,4,5] c1[21,7,10,11,8,10,11,9,22] c2[12,13,14,15,16,17,18,23,19]", dumpAsVector(m));
    PIN("f3.mapindex", "[20,0,1,6,6,2,3,4,5,21,7,10,11,8,10,11,9,22,12,13,14,15,16,17,18,23,19]", dumpMapIndex(m));
    PIN("f3.index", "c0|20,0,1,6|6,2,3,4,5;c1|21,7,10,11|8,10,11,9,22;c2|12,13,14,15|16,17,18,23,19", dumpIndex(m));
    PIN("f3.flags", "BFFF|FFFFF;BFFF|FFFFB;FFFF|FFFBF", dumpFreeBoundaryFlags(m));
    PIN("f3.bindex", "[0:0:0->0,1:0:0->1,1:1:4->2,2:1:3->3]", dumpBindex(m));
    PIN("f3.gflags", "FFFFFFFFFFFFFFFFFFFFBBBB gbi=[20->0,21->1,22->2,23->3]", dumpGlobalFlags(m, 0));
    PIN("f3.findbf", "bnd=[0],[] free=[1,2,3],[0,1,2,3,4]", dumpFindBoundaryFree(m));
    PIN("f3.onpatch", "0:1/- 1:2/- 2:-/1 3:-/2 4:-/3 5:-/4 6:3/0 7:1/- 8:-/0 9:-/3 10:2/1 11:3/2 12:0/- 13:1/- 14:2/- 15:3/- 16:-/0 17:-/1 18:-/2 19:-/4 20:0/- 21:0/- 22:-/4 23:-/3", dumpIndexOnPatch(m));

    PIN("f3.componentof", "[0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,2,2,2,0,1,1,2]", dumpComponentOf(m));
    PIN("f3.preimage", "0[(0,1)] 1[(0,2)] 2[(1,1)] 3[(1,2)] 4[(1,3)] 5[(1,4)] 6[(0,3),(1,0)] 7[(0,1)] 8[(1,0)] 9[(1,3)] 10[(0,2),(1,1)] 11[(0,3),(1,2)] 12[(0,0)] 13[(0,1)] 14[(0,2)] 15[(0,3)] 16[(1,0)] 17[(1,1)] 18[(1,2)] 19[(1,4)] 20[(0,0)] 21[(0,0)] 22[(1,4)] 23[(1,3)]", dumpPreImages(m));
    PIN("f3.anypreimage", "[(0,1),(0,2),(1,1),(1,2),(1,3),(1,4),(0,3),(0,1),(1,0),(1,3),(0,2),(0,3),(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,4),(0,0),(0,0),(1,4),(1,3)]", dumpAnyPreImage(m));
    // See any_pre_images_agrees_with_any_pre_image for anyPreImages(c) on
    // this fixture; it is checked against anyPreImage rather than pinned as
    // a literal.
    PIN("f3.tagged", "tagged=[] n=0 gflags=........................ flags=....|.....;....|.....;....|..... tindex=0,0,0,0|0,0,0,0,0;0,0,0,0|0,0,0,0,0;0,0,0,0|0,0,0,0,0", dumpTagged(m));
}

// =========================================================================
// F4 -- F3 with a global and a boundary shift
// =========================================================================

// The shift relabels the global indices and changes nothing else; that
// relation is checked query by query in shift_is_a_relabelling.  Pinned here
// are the values the shifts themselves enter.
//
// Indices below the shift belong to no patch; see out_of_range_global_indices.
TEST(three_comp_shifted)
{
    const gsDofMapper m = threeCompShifted();
    PIN("f4.firstlast", "first0=100 last=120", dumpFirstLast(m));   // firstFreeIndex(0) == m_shift
    PIN("f4.asvector", "c0[120,100,101,106,106,102,103,104,105] c1[121,107,110,111,108,110,111,109,122] c2[112,113,114,115,116,117,118,123,119]", dumpAsVector(m));
    // freeIndex(i,k,c) returns the index WITHOUT m_shift; read against
    // f4.index just below, this literal is the pin on that distinction.
    PIN("f4.freeindex", "c0|20,0,1,6|6,2,3,4,5;c1|21,7,10,11|8,10,11,9,22;c2|12,13,14,15|16,17,18,23,19", dumpFreeIndex(m));
    PIN("f4.index", "c0|120,100,101,106|106,102,103,104,105;c1|121,107,110,111|108,110,111,109,122;c2|112,113,114,115|116,117,118,123,119", dumpIndex(m));
    PIN("f4.bindex", "[0:0:0->7,1:0:0->8,1:1:4->9,2:1:3->10]", dumpBindex(m));
    PIN("f4.gflags", "FFFFFFFFFFFFFFFFFFFFBBBB gbi=[120->7,121->8,122->9,123->10]", dumpGlobalFlags(m, 100));
}

// =========================================================================
// F5 -- identity / aliased layout
// =========================================================================

// The global-identity/aliased layout that setIdentity() produces: every
// patch carries the component's complete identity range.  patchSize(p,c) is
// the component total on every patch, offset(p,c) is zero for every real
// patch, and preImage/anyPreImage report patch 0 as the canonical preimage.
// So f5.layout, f5.index, f5.flags, f5.findbf, f5.onpatch and f5.tagged (whose
// is_tagged/tindex columns loop over patchSize per patch) show the full
// 7-entry range on each of the three patches.
TEST(identity_mapper)
{
    const gsDofMapper m = identityMapper();
    // boundarySizeWithDuplicates counts the eliminated dofs with their
    // multiplicity over patches, so it is bounded below by boundarySize().
    PIN("f5.elimdup", "elimdup=0", dumpElimDup(m));
    PIN("f5.counts", "comps=2 ncomp=2 npatch=3 map=14 size=14 free=14 elim=0 cpld=0 tagged=0 allfree=1 perm=1 final=1", dumpCounts(m));
    PIN("f5.percomp", "c0:size=7,free=7,total=7 c1:size=7,free=7,total=7", dumpPerComponent(m));
    PIN("f5.layout", "off=[0,0,0] ps=c0:[7,7,7] c1:[7,7,7]", dumpLayout(m));
    PIN("f5.firstlast", "first0=0 last=14", dumpFirstLast(m));
    PIN("f5.asvector", "c0[0,1,2,3,4,5,6] c1[7,8,9,10,11,12,13]", dumpAsVector(m));
    PIN("f5.mapindex", "[0,1,2,3,4,5,6,7,8,9,10,11,12,13]", dumpMapIndex(m));
    PIN("f5.index", "c0|0,1,2,3,4,5,6|0,1,2,3,4,5,6|0,1,2,3,4,5,6;c1|7,8,9,10,11,12,13|7,8,9,10,11,12,13|7,8,9,10,11,12,13", dumpIndex(m));
    PIN("f5.flags", "FFFFFFF|FFFFFFF|FFFFFFF;FFFFFFF|FFFFFFF|FFFFFFF", dumpFreeBoundaryFlags(m));
    PIN("f5.bindex", "[]", dumpBindex(m));
    PIN("f5.gflags", "FFFFFFFFFFFFFF gbi=[]", dumpGlobalFlags(m, 0));
    PIN("f5.findbf", "bnd=[],[],[] free=[0,1,2,3,4,5,6],[0,1,2,3,4,5,6],[0,1,2,3,4,5,6]", dumpFindBoundaryFree(m));
    PIN("f5.onpatch", "0:0/0/0 1:1/1/1 2:2/2/2 3:3/3/3 4:4/4/4 5:5/5/5 6:6/6/6 7:0/0/0 8:1/1/1 9:2/2/2 10:3/3/3 11:4/4/4 12:5/5/5 13:6/6/6", dumpIndexOnPatch(m));
    PIN("f5.componentof", "[0,0,0,0,0,0,0,1,1,1,1,1,1,1]", dumpComponentOf(m));
    PIN("f5.preimage", "0[(0,0)] 1[(0,1)] 2[(0,2)] 3[(0,3)] 4[(0,4)] 5[(0,5)] 6[(0,6)] 7[(0,0)] 8[(0,1)] 9[(0,2)] 10[(0,3)] 11[(0,4)] 12[(0,5)] 13[(0,6)]", dumpPreImages(m));
    PIN("f5.anypreimage", "[(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6),(0,0),(0,1),(0,2),(0,3),(0,4),(0,5),(0,6)]", dumpAnyPreImage(m));
    PIN("f5.tagged", "tagged=[] n=0 gflags=.............. flags=.......|.......|.......;.......|.......|....... tindex=0,0,0,0,0,0,0|0,0,0,0,0,0,0|0,0,0,0,0,0,0;0,0,0,0,0,0,0|0,0,0,0,0,0,0|0,0,0,0,0,0,0", dumpTagged(m));
}

// =========================================================================
// F6 -- F3 with markCoupledAsTagged and a permuted component-0 free block
// =========================================================================

// The cpld=2 in f6.counts: permuting a component destroys the ability to track ITS coupled dofs, so
// that component's own coupled count drops to zero and every later
// cumulative prefix loses exactly that count: F3 has coupled counts 1, 2, 0
// per component (cumulative 1, 3, 3), so permuting component 0 leaves
// cumulative 0, 2, 2 and coupledSize() == 2 -- component 1's two coupled
// dofs are still tracked.
//
// f6.tagged is the composition of the two tag operations.
// markCoupledAsTagged tags F3's coupled dofs, 6, 10 and 11; permuting
// component 0 then moves the one inside that component's free block [0,7),
// 6 -> perm[6] == 4, and leaves 10 and 11 -- which belong to component 1 --
// where they are.  See
// mark_coupled_as_tagged_tags_the_coupled_dofs and
// permute_free_dofs_keeps_other_components_tags.
TEST(permuted_mapper)
{
    const gsDofMapper m = permutedMapper();
    // boundarySizeWithDuplicates counts the eliminated dofs with their
    // multiplicity over patches, so it is bounded below by boundarySize().
    PIN("f6.elimdup", "elimdup=4", dumpElimDup(m));
    PIN("f6.counts", "comps=3 ncomp=3 npatch=2 map=27 size=24 free=20 elim=4 cpld=2 tagged=3 allfree=0 perm=0 final=1", dumpCounts(m));
    PIN("f6.percomp", "c0:size=8,free=7,total=9 c1:size=7,free=5,total=9 c2:size=9,free=8,total=9", dumpPerComponent(m));
    PIN("f6.layout", "off=[0,4] ps=c0:[4,5] c1:[4,5] c2:[4,5]", dumpLayout(m));
    PIN("f6.firstlast", "first0=0 last=20", dumpFirstLast(m));
    PIN("f6.asvector", "c0[20,3,0,4,4,6,1,5,2] c1[21,7,10,11,8,10,11,9,22] c2[12,13,14,15,16,17,18,23,19]", dumpAsVector(m));
    PIN("f6.mapindex", "[20,3,0,4,4,6,1,5,2,21,7,10,11,8,10,11,9,22,12,13,14,15,16,17,18,23,19]", dumpMapIndex(m));
    PIN("f6.index", "c0|20,3,0,4|4,6,1,5,2;c1|21,7,10,11|8,10,11,9,22;c2|12,13,14,15|16,17,18,23,19", dumpIndex(m));
    PIN("f6.flags", "BFFF|FFFFF;BFFF|FFFFB;FFFF|FFFBF", dumpFreeBoundaryFlags(m));
    PIN("f6.bindex", "[0:0:0->0,1:0:0->1,1:1:4->2,2:1:3->3]", dumpBindex(m));
    PIN("f6.gflags", "FFFFFFFFFFFFFFFFFFFFBBBB gbi=[20->0,21->1,22->2,23->3]", dumpGlobalFlags(m, 0));
    PIN("f6.findbf", "bnd=[0],[] free=[1,2,3],[0,1,2,3,4]", dumpFindBoundaryFree(m));
    PIN("f6.onpatch", "0:2/- 1:-/2 2:-/4 3:1/- 4:3/0 5:-/3 6:-/1 7:1/- 8:-/0 9:-/3 10:2/1 11:3/2 12:0/- 13:1/- 14:2/- 15:3/- 16:-/0 17:-/1 18:-/2 19:-/4 20:0/- 21:0/- 22:-/4 23:-/3", dumpIndexOnPatch(m));
    PIN("f6.componentof", "[0,0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,2,2,2,0,1,1,2]", dumpComponentOf(m));
    PIN("f6.preimage", "0[(0,2)] 1[(1,2)] 2[(1,4)] 3[(0,1)] 4[(0,3),(1,0)] 5[(1,3)] 6[(1,1)] 7[(0,1)] 8[(1,0)] 9[(1,3)] 10[(0,2),(1,1)] 11[(0,3),(1,2)] 12[(0,0)] 13[(0,1)] 14[(0,2)] 15[(0,3)] 16[(1,0)] 17[(1,1)] 18[(1,2)] 19[(1,4)] 20[(0,0)] 21[(0,0)] 22[(1,4)] 23[(1,3)]", dumpPreImages(m));
    PIN("f6.anypreimage", "[(0,2),(1,2),(1,4),(0,1),(0,3),(1,3),(1,1),(0,1),(1,0),(1,3),(0,2),(0,3),(0,0),(0,1),(0,2),(0,3),(1,0),(1,1),(1,2),(1,4),(0,0),(0,0),(1,4),(1,3)]", dumpAnyPreImage(m));
    PIN("f6.tagged", "tagged=[4,10,11] n=3 gflags=....T.....TT............ flags=...T|T....;..TT|.TT..;....|..... tindex=3,0,0,0|0,1,0,1,0;3,1,1,2|1,1,2,1,3;3,3,3,3|3,3,3,3,3", dumpTagged(m));
}

// =========================================================================
// F7 -- collapseDofs
// =========================================================================

TEST(collapsed_mapper)
{
    const gsDofMapper m = collapsedMapper();
    PIN("f7.elimdup", "elimdup=0", dumpElimDup(m));
    PIN("f7.counts", "comps=1 ncomp=1 npatch=2 map=10 size=8 free=8 elim=0 cpld=1 tagged=0 allfree=1 perm=0 final=1", dumpCounts(m));
    PIN("f7.percomp", "c0:size=8,free=8,total=10", dumpPerComponent(m));
    PIN("f7.layout", "off=[0,5] ps=c0:[5,5]", dumpLayout(m));
    PIN("f7.firstlast", "first0=0 last=8", dumpFirstLast(m));
    PIN("f7.asvector", "c0[7,0,7,1,7,2,3,4,5,6]", dumpAsVector(m));
    PIN("f7.mapindex", "[7,0,7,1,7,2,3,4,5,6]", dumpMapIndex(m));
    PIN("f7.index", "c0|7,0,7,1,7|2,3,4,5,6", dumpIndex(m));
    PIN("f7.flags", "FFFFF|FFFFF", dumpFreeBoundaryFlags(m));
    PIN("f7.bindex", "[]", dumpBindex(m));
    PIN("f7.gflags", "FFFFFFFF gbi=[]", dumpGlobalFlags(m, 0));
    PIN("f7.findbf", "bnd=[],[] free=[0,1,2,3,4],[0,1,2,3,4]", dumpFindBoundaryFree(m));
    PIN("f7.onpatch", "0:1/- 1:3/- 2:-/0 3:-/1 4:-/2 5:-/3 6:-/4 7:0/-", dumpIndexOnPatch(m));
    PIN("f7.coupled", "flags=C.C.C|..... cidx=[0:0->0,0:2->0,0:4->0] gflags=.......C findCoupled=[0,2,4],[] findCoupledPairs=00[],01[],10[],11[] findFreeUncoupled=[1,3],[0,1,2,3,4]", dumpCoupledQueries(m));
    PIN("f7.componentof", "[0,0,0,0,0,0,0,0]", dumpComponentOf(m));
    PIN("f7.preimage", "0[(0,1)] 1[(0,3)] 2[(1,0)] 3[(1,1)] 4[(1,2)] 5[(1,3)] 6[(1,4)] 7[(0,0),(0,2),(0,4)]", dumpPreImages(m));
    PIN("f7.anypreimage", "[(0,1),(0,3),(1,0),(1,1),(1,2),(1,3),(1,4),(0,0)]", dumpAnyPreImage(m));
    // One entry per global index, size() == 8; see f2.anypreimages.
    PIN("f7.anypreimages", "[(0,1),(0,3),(1,0),(1,1),(1,2),(1,3),(1,4),(0,0)]", join(m.anyPreImages(0)));
    PIN("f7.tagged", "tagged=[] n=0 gflags=........ flags=.....|..... tindex=0,0,0,0,0|0,0,0,0,0", dumpTagged(m));
}

// =========================================================================
// F8 -- the gsDofMapperCreator path on a real two-patch basis
// =========================================================================

TEST(creator_two_patch)
{
    const gsDofMapper m = creatorTwoPatch();
    PIN("f8.elimdup", "elimdup=0", dumpElimDup(m));
    PIN("f8.counts", "comps=1 ncomp=1 npatch=2 map=50 size=45 free=45 elim=0 cpld=5 tagged=0 allfree=1 perm=0 final=1", dumpCounts(m));
    PIN("f8.percomp", "c0:size=45,free=45,total=50", dumpPerComponent(m));
    PIN("f8.layout", "off=[0,25] ps=c0:[25,25]", dumpLayout(m));
    PIN("f8.firstlast", "first0=0 last=45", dumpFirstLast(m));

    PIN("f8.asvector", "c0[0,1,2,3,40,4,5,6,7,41,8,9,10,11,42,12,13,14,15,43,16,17,18,19,44,40,20,21,22,23,41,24,25,26,27,42,28,29,30,31,43,32,33,34,35,44,36,37,38,39]", dumpAsVector(m));
    PIN("f8.mapindex", "[0,1,2,3,40,4,5,6,7,41,8,9,10,11,42,12,13,14,15,43,16,17,18,19,44,40,20,21,22,23,41,24,25,26,27,42,28,29,30,31,43,32,33,34,35,44,36,37,38,39]", dumpMapIndex(m));
    PIN("f8.index", "c0|0,1,2,3,40,4,5,6,7,41,8,9,10,11,42,12,13,14,15,43,16,17,18,19,44|40,20,21,22,23,41,24,25,26,27,42,28,29,30,31,43,32,33,34,35,44,36,37,38,39", dumpIndex(m));
    PIN("f8.flags", "FFFFFFFFFFFFFFFFFFFFFFFFF|FFFFFFFFFFFFFFFFFFFFFFFFF", dumpFreeBoundaryFlags(m));
    PIN("f8.bindex", "[]", dumpBindex(m));
    PIN("f8.gflags", "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF gbi=[]", dumpGlobalFlags(m, 0));
    PIN("f8.findbf", "bnd=[],[] free=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24],[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24]", dumpFindBoundaryFree(m));
    PIN("f8.onpatch", "0:0/- 1:1/- 2:2/- 3:3/- 4:5/- 5:6/- 6:7/- 7:8/- 8:10/- 9:11/- 10:12/- 11:13/- 12:15/- 13:16/- 14:17/- 15:18/- 16:20/- 17:21/- 18:22/- 19:23/- 20:-/1 21:-/2 22:-/3 23:-/4 24:-/6 25:-/7 26:-/8 27:-/9 28:-/11 29:-/12 30:-/13 31:-/14 32:-/16 33:-/17 34:-/18 35:-/19 36:-/21 37:-/22 38:-/23 39:-/24 40:4/0 41:9/5 42:14/10 43:19/15 44:24/20", dumpIndexOnPatch(m));

    PIN("f8.coupled", "flags=....C....C....C....C....C|C....C....C....C....C.... cidx=[0:4->0,0:9->1,0:14->2,0:19->3,0:24->4,1:0->0,1:5->1,1:10->2,1:15->3,1:20->4] gflags=........................................CCCCC findCoupled=[4,9,14,19,24],[0,5,10,15,20] findCoupledPairs=00[],01[4,9,14,19,24],10[0,5,10,15,20],11[] findFreeUncoupled=[0,1,2,3,5,6,7,8,10,11,12,13,15,16,17,18,20,21,22,23],[1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,21,22,23,24]", dumpCoupledQueries(m));

    PIN("f8.componentof", "[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]", dumpComponentOf(m));
    PIN("f8.preimage", "0[(0,0)] 1[(0,1)] 2[(0,2)] 3[(0,3)] 4[(0,5)] 5[(0,6)] 6[(0,7)] 7[(0,8)] 8[(0,10)] 9[(0,11)] 10[(0,12)] 11[(0,13)] 12[(0,15)] 13[(0,16)] 14[(0,17)] 15[(0,18)] 16[(0,20)] 17[(0,21)] 18[(0,22)] 19[(0,23)] 20[(1,1)] 21[(1,2)] 22[(1,3)] 23[(1,4)] 24[(1,6)] 25[(1,7)] 26[(1,8)] 27[(1,9)] 28[(1,11)] 29[(1,12)] 30[(1,13)] 31[(1,14)] 32[(1,16)] 33[(1,17)] 34[(1,18)] 35[(1,19)] 36[(1,21)] 37[(1,22)] 38[(1,23)] 39[(1,24)] 40[(0,4),(1,0)] 41[(0,9),(1,5)] 42[(0,14),(1,10)] 43[(0,19),(1,15)] 44[(0,24),(1,20)]", dumpPreImages(m));
    PIN("f8.anypreimage", "[(0,0),(0,1),(0,2),(0,3),(0,5),(0,6),(0,7),(0,8),(0,10),(0,11),(0,12),(0,13),(0,15),(0,16),(0,17),(0,18),(0,20),(0,21),(0,22),(0,23),(1,1),(1,2),(1,3),(1,4),(1,6),(1,7),(1,8),(1,9),(1,11),(1,12),(1,13),(1,14),(1,16),(1,17),(1,18),(1,19),(1,21),(1,22),(1,23),(1,24),(0,4),(0,9),(0,14),(0,19),(0,24)]", dumpAnyPreImage(m));
    // One entry per global index, size() == 45; see f2.anypreimages.
    PIN("f8.anypreimages", "[(0,0),(0,1),(0,2),(0,3),(0,5),(0,6),(0,7),(0,8),(0,10),(0,11),(0,12),(0,13),(0,15),(0,16),(0,17),(0,18),(0,20),(0,21),(0,22),(0,23),(1,1),(1,2),(1,3),(1,4),(1,6),(1,7),(1,8),(1,9),(1,11),(1,12),(1,13),(1,14),(1,16),(1,17),(1,18),(1,19),(1,21),(1,22),(1,23),(1,24),(0,4),(0,9),(0,14),(0,19),(0,24)]", join(m.anyPreImages(0)));
}

// addShift adds to the current shift instead of replacing it, so after
// setShift(100) + addShift(5) the global numbering starts at 105.
TEST(add_shift_accumulates)
{
    gsDofMapper m = threeCompUniform();
    m.setShift(100);
    m.addShift(5);
    PIN("shift.accumulates", "first0=105 last=125", dumpFirstLast(m));
}

// =========================================================================
// Independent numbering oracle
// =========================================================================

// For an uncoupled, uneliminated, patch-concatenated multi-component mapper
// the finalized numbering is component-major:
//     index(i,p,c) == sum_{c'<c} totalSize(c') + offset(p) + i.
// This is not a recorded literal: it is an independent check on the
// component-major numbering that any reimplementation must keep satisfying.
TEST(component_major_numbering_oracle)
{
    gsVector<index_t> sz(2);
    sz[0] = 3; sz[1] = 4;
    gsDofMapper m(sz, 2);
    m.finalize();

    CHECK_EQUAL(0, m.boundarySize());
    CHECK_EQUAL(0, m.coupledSize());
    CHECK_EQUAL(14, m.size());

    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        index_t base = 0;
        for (index_t cp = 0; cp != c; ++cp)
            base += static_cast<index_t>(m.totalSize(cp));

        for (index_t p = 0; p != nPatches(m); ++p)
        {
            const index_t n = static_cast<index_t>(m.patchSize(p, c));
            for (index_t i = 0; i != n; ++i)
                CHECK_EQUAL(base + static_cast<index_t>(m.offset(p)) + i,
                            m.index(i, p, c));
        }
    }
}

// =========================================================================
// F9/F10 -- genuinely ragged mappers (unequal per-component sizes)
// =========================================================================
//
// Every other multi-component fixture in this file (F3/F4/F6 =
// threeCompUniform() and its derivatives, F5 = identityMapper()) has EQUAL
// per-component sizes -- only the coupling/elimination *pattern* differs --
// so a query that bounds its work by component 0's size, patchSize(k), where
// it needs the queried component's, patchSize(k,c), is only visible on
// these and in ragged_component_is_a_relabelled_scalar_mapper.

TEST(ragged_identity_indexOnPatch)
{
    const gsDofMapper m = raggedIdentityMapper();
    CHECK_EQUAL(2u, (unsigned)m.componentsSize());
    CHECK(gsDofMapper::GlobalIdentity == m.layout());
    CHECK_EQUAL(7u,  (unsigned)m.totalSize(0));
    CHECK_EQUAL(10u, (unsigned)m.totalSize(1));

    // Component 1, local index 7 (global 14) lies beyond component 0's total
    // (7) dofs.  Aliased layout: it is found -- with the SAME local index --
    // on every patch.
    index_t local = -12345;
    for (index_t k = 0; k != nPatches(m); ++k)
    {
        CHECK(m.indexOnPatch(14, k, local));
        CHECK_EQUAL(7, local);
    }

    // Component 0's own range.
    for (index_t k = 0; k != nPatches(m); ++k)
    {
        CHECK(m.indexOnPatch(0, k, local));
        CHECK_EQUAL(0, local);
    }
}


// =========================================================================
// F11/F12 -- the Raviart-Thomas pair, and the declared-distinctness contract
// =========================================================================

// Pins the dimension table this whole feature turns on: on an isotropic mesh the
// two RT component spaces have exactly the same dimension, and only an
// anisotropic mesh separates them.  If this ever stops holding, the fixtures
// below stop testing what they claim to test.
TEST(rt_component_dimensions)
{
    CHECK_EQUAL(42, rtComponentBasis(3, 2, 4, 4).size());
    CHECK_EQUAL(42, rtComponentBasis(2, 3, 4, 4).size());

    CHECK_EQUAL(28, rtComponentBasis(3, 2, 4, 2).size());
    CHECK_EQUAL(30, rtComponentBasis(2, 3, 4, 2).size());
}

// The single most important test in this file.  An isotropic-mesh RT mapper
// has identical per-patch cardinalities in both components, so every
// size-based predicate reports it "uniform" -- and the uniform evaluator
// would then assemble it with component 0's basis replicated, silently and
// with wrong results, on precisely the configuration this feature exists
// for.  The rejection must therefore come from the declared flag alone.
TEST(rt_isotropic_rejected_on_the_declared_flag_alone)
{
    const gsDofMapper m = rtIsotropicMapper();

    // Size-blind: nothing observable about the storage distinguishes this
    // mapper from an ordinary uniform two-component one.
    CHECK_EQUAL(2u, (unsigned)m.componentsSize());
    CHECK_EQUAL(42u, (unsigned)m.patchSize(0,0));
    CHECK_EQUAL(42u, (unsigned)m.patchSize(0,1));
    CHECK_EQUAL(m.patchSize(0,0), m.patchSize(0,1));
    CHECK_EQUAL(m.totalSize(0), m.totalSize(1));

    CHECK(m.hasDistinctComponentSpaces());
    CHECK(!m.hasUniformComponents());

    // Control: the very same storage shape, declared NOT distinct, is
    // accepted.  The two mappers are observably identical apart from the
    // declaration, which is what proves the rejection is not size-based.
    const gsDofMapper u = rtMapper(4, 4, /*declareDistinct=*/false);
    CHECK_EQUAL(m.patchSize(0,0), u.patchSize(0,0));
    CHECK_EQUAL(m.patchSize(0,1), u.patchSize(0,1));
    CHECK(!u.hasDistinctComponentSpaces());
    CHECK(u.hasUniformComponents());
}

// The anisotropic mesh is rejected twice over: by the declared flag, and --
// independently -- by the cardinality conjunct, which is what a mapper built
// from unequal component sizes trips even when nothing was declared.
TEST(rt_anisotropic_rejected_on_flag_and_on_cardinality)
{
    const gsDofMapper m = rtAnisotropicMapper();
    CHECK_EQUAL(28u, (unsigned)m.patchSize(0,0));
    CHECK_EQUAL(30u, (unsigned)m.patchSize(0,1));
    CHECK(m.hasDistinctComponentSpaces());
    CHECK(!m.hasUniformComponents());

    const gsDofMapper u = rtMapper(4, 2, /*declareDistinct=*/false);
    CHECK(!u.hasDistinctComponentSpaces());
    CHECK(!u.hasUniformComponents()); // rejected by cardinality alone
}

// The full truth table of the compatibility predicate over every fixture in
// this file.  Everything the current mainline builds must stay accepted --
// including a default-constructed mapper, which is the normal pre-init state
// the legacy boundary sees and must not fire on.
TEST(has_uniform_components_truth_table)
{
    const gsDofMapper empty;
    CHECK(!empty.hasDistinctComponentSpaces());
    CHECK(empty.hasUniformComponents());

    CHECK(singlePatchPlain()    .hasUniformComponents()); // 1 component
    CHECK(twoPatchCoupledElim() .hasUniformComponents()); // 1 component
    CHECK(threeCompUniform()    .hasUniformComponents()); // 3 equal comps
    CHECK(identityMapper()      .hasUniformComponents()); // uniform alias
    CHECK(creatorTwoPatch()     .hasUniformComponents()); // creator-built

    // None of the mainline paths declares distinctness.
    CHECK(!singlePatchPlain()   .hasDistinctComponentSpaces());
    CHECK(!threeCompUniform()   .hasDistinctComponentSpaces());
    CHECK(!identityMapper()     .hasDistinctComponentSpaces());
    CHECK(!creatorTwoPatch()    .hasDistinctComponentSpaces());

    // Ragged cardinalities are rejected under either layout.
    CHECK(!raggedPatchMapper()   .hasUniformComponents());
    CHECK(!raggedIdentityMapper().hasUniformComponents());
}

// A mapper whose components have equal totals but a different per-patch
// split is uniform by every total-size measure and still unusable: the
// cardinality conjunct is per patch, not per component total.
TEST(has_uniform_components_is_per_patch)
{
    std::vector<gsVector<index_t> > sz(2);
    sz[0].resize(2); sz[0][0] = 4; sz[0][1] = 6;
    sz[1].resize(2); sz[1][0] = 6; sz[1][1] = 4;
    gsDofMapper m(sz, /*hasDistinctComponentSpaces=*/false);
    m.finalize();

    CHECK_EQUAL(m.totalSize(0), m.totalSize(1)); // 10 == 10
    CHECK(!m.hasUniformComponents());
}

// =========================================================================
// Metadata carried through swap
// =========================================================================

TEST(swap_carries_layout_metadata)
{
    gsDofMapper a = rtAnisotropicMapper();   // 1 patch, concatenated, distinct
    gsDofMapper b = raggedIdentityMapper();  // 3 patches, identity, not distinct

    a.swap(b);

    CHECK_EQUAL(3u, (unsigned)a.numPatches());
    CHECK(gsDofMapper::GlobalIdentity == a.layout());
    CHECK(!a.hasDistinctComponentSpaces());
    CHECK_EQUAL(7u,  (unsigned)a.totalSize(0));
    CHECK_EQUAL(10u, (unsigned)a.totalSize(1));

    CHECK_EQUAL(1u, (unsigned)b.numPatches());
    CHECK(gsDofMapper::PatchConcatenated == b.layout());
    CHECK(b.hasDistinctComponentSpaces());
    CHECK_EQUAL(28u, (unsigned)b.patchSize(0,0));
    CHECK_EQUAL(30u, (unsigned)b.patchSize(0,1));
}

// =========================================================================
// setIdentity() as a reset of an already-populated mapper
// =========================================================================

// setIdentity() rebuilds the mapping, the counts and the layout, so it must
// also drop the tag list: its entries are global indices of the numbering
// being thrown away.  Left behind, they survive into a numbering that does
// not contain them, and taggedSize()/getTagged()/is_tagged_index() report
// dofs that do not exist.
TEST(set_identity_clears_stale_tags)
{
    gsDofMapper m = twoPatchCoupledElim();
    m.markCoupledAsTagged();
    CHECK(m.taggedSize() > 0);
    const index_t staleTag = m.getTagged().back();

    m.setIdentity(1, 3, 1);
    m.finalize();

    CHECK_EQUAL(0, m.taggedSize());
    CHECK(m.getTagged().empty());
    CHECK_EQUAL(3, m.size());
    // The stale tag was an index of the discarded numbering; nothing in the
    // new one may answer to it.
    for (index_t gl = 0; gl != m.size(); ++gl)
        CHECK(!m.is_tagged_index(gl));
    CHECK(!m.is_tagged_index(staleTag));
}

// Re-initialising a populated mapper with a smaller shape must leave nothing
// of the old one behind -- the reason every initialiser assigns rather than
// resizes.
TEST(set_identity_reinitializes_completely)
{
    gsDofMapper m = threeCompUniform();
    CHECK_EQUAL(3, m.numComponents());

    std::vector<size_t> dofs(2);
    dofs[0] = 4; dofs[1] = 6;
    m.setIdentity(2, dofs);
    m.finalize();

    CHECK_EQUAL(2, m.numComponents());
    CHECK_EQUAL(2u, (unsigned)m.numPatches());
    CHECK(gsDofMapper::GlobalIdentity == m.layout());
    CHECK(!m.hasDistinctComponentSpaces());
    CHECK_EQUAL(4u, (unsigned)m.totalSize(0));
    CHECK_EQUAL(6u, (unsigned)m.totalSize(1));
    CHECK_EQUAL(10u, (unsigned)m.mapSize());
    CHECK_EQUAL(10, m.size());
    CHECK_EQUAL(0, m.boundarySize());
    CHECK(!m.hasUniformComponents()); // ragged identity totals
}

// =========================================================================
// Construction-time validation
// =========================================================================

TEST(ragged_constructor_rejects_invalid_metadata)
{
    // No components at all.
    CHECK_THROW(gsDofMapper(std::vector<gsVector<index_t> >(), false),
                std::runtime_error);

    // No patches.
    {
        std::vector<gsVector<index_t> > sz(1);
        CHECK_THROW(gsDofMapper(sz, false), std::runtime_error);
    }

    // Components disagreeing about the patch count: ragged means different
    // sizes per (component,patch), never different patch sets.
    {
        std::vector<gsVector<index_t> > sz(2);
        sz[0].resize(2); sz[0][0] = 3; sz[0][1] = 3;
        sz[1].resize(3); sz[1][0] = 3; sz[1][1] = 3; sz[1][2] = 3;
        CHECK_THROW(gsDofMapper(sz, false), std::runtime_error);
    }

    // Negative patch size.
    {
        std::vector<gsVector<index_t> > sz(1);
        sz[0].resize(2); sz[0][0] = 3; sz[0][1] = -1;
        CHECK_THROW(gsDofMapper(sz, false), std::runtime_error);
    }
}

TEST(set_identity_rejects_invalid_metadata)
{
    gsDofMapper m;
    CHECK_THROW(m.setIdentity(0, 5, 1), std::runtime_error);   // no patches
    CHECK_THROW(m.setIdentity(-1, 5, 1), std::runtime_error);  // no patches
    CHECK_THROW(m.setIdentity(1, 5, 0), std::runtime_error);   // no components
    CHECK_THROW(m.setIdentity(1, std::vector<size_t>()), std::runtime_error);
}

// Dof counts are stored, and accumulated across components by finalize(), in
// index_t -- a build-time configurable type that is plain int by default and
// may be far narrower than size_t.  A count that overflows it must be
// rejected against index_t's range, and -- this is the operative part --
// BEFORE anything is allocated or narrowed: validating after the fact would
// mean either a wrapped negative count or a multi-gigabyte allocation on the
// way to the diagnostic.  Every case below therefore must throw without the
// mapper ever sizing its storage.
TEST(construction_rejects_counts_beyond_index_range)
{
    const size_t imax = static_cast<size_t>(std::numeric_limits<index_t>::max());

    // Cumulative over patches (one component).
    {
        gsVector<index_t> sz(2);
        sz[0] = std::numeric_limits<index_t>::max();
        sz[1] = std::numeric_limits<index_t>::max();
        CHECK_THROW(gsDofMapper(sz, 1), std::runtime_error);
    }

    // Cumulative over components: each component fits on its own, their sum
    // does not.  This is the case a per-component-only check misses.
    {
        gsVector<index_t> sz(1);
        sz[0] = static_cast<index_t>(imax/2 + 1);
        CHECK_THROW(gsDofMapper(sz, 2), std::runtime_error);
    }

    // Same, through the ragged constructor.
    {
        std::vector<gsVector<index_t> > sz(2);
        sz[0].resize(1); sz[0][0] = static_cast<index_t>(imax/2 + 1);
        sz[1].resize(1); sz[1][0] = static_cast<index_t>(imax/2 + 1);
        CHECK_THROW(gsDofMapper(sz, false), std::runtime_error);
    }

    // setIdentity takes size_t counts, so a single component can exceed
    // index_t's range on its own.  Only meaningful where index_t is no wider
    // than size_t, which is every supported configuration.
    if (sizeof(index_t) <= sizeof(size_t))
    {
        gsDofMapper m;
        CHECK_THROW(m.setIdentity(1, std::vector<size_t>(1, imax + 1)),
                    std::runtime_error);

        std::vector<size_t> dofs(2);
        dofs[0] = imax; dofs[1] = 1;    // each fits, the sum does not
        CHECK_THROW(m.setIdentity(1, dofs), std::runtime_error);
    }

    // The component count is narrowed to index_t by numComponents() and
    // passes every dof-count check when the components are empty, so it is
    // bounded on its own.  Only the scalar overload can be exercised here:
    // it must reject the count before allocating its per-component vector,
    // whereas reaching the other two entry points would take a vector of
    // imax+1 elements.  (On a build with a narrow index_t -- int8_t, 128
    // components -- those are reachable with ordinary inputs.)
    if (sizeof(index_t) <= sizeof(size_t))
    {
        gsDofMapper m;
        CHECK_THROW(m.setIdentity(1, 0, imax + 1), std::runtime_error);
    }
}


// =========================================================================
// Per-component prefix queries
// =========================================================================
//
// After finalize() the per-component count vectors m_numFreeDofs,
// m_numElimDofs and m_numCpldDofs are CUMULATIVE prefix sums, so a
// component's own count is always a prefix difference.  Reading a prefix
// where a difference is meant, or the last component's prefix where the
// queried component's is meant, is correct for component 0 and for a
// single-component mapper and wrong everywhere else -- which is why every
// test below needs a fixture with three components and asserts on the
// components >= 1.
//
// F3 (threeCompUniform) has, per component:
//
//   component | local storage per patch          | free | cpld | elim
//   ----------+----------------------------------+------+------+-----
//       0     | [20,0,1,6]  [6,2,3,4,5]          |   7  |   1  |  1
//       1     | [21,7,10,11] [8,10,11,9,22]      |   5  |   2  |  2
//       2     | [12,13,14,15] [16,17,18,23,19]   |   8  |   0  |  1
//
// so the free blocks are [0,7), [7,12), [12,20), the coupled dofs sit at the
// top of each ({6}, {10,11}, {}) and the eliminated ones are 20, {21,22} and
// 23.  Every literal below is read off that table.

TEST(three_comp_coupled_queries_per_component)
{
    const gsDofMapper m = threeCompUniform();

    PIN("f3.coupled.c0",
        "flags=...C|C.... cidx=[0:3->0,1:0->0] gflags=......C...CC............"
        " findCoupled=[3],[0] findCoupledPairs=00[],01[3],10[0],11[]"
        " findFreeUncoupled=[1,2],[1,2,3,4]",
        dumpCoupledQueries(m, 0));

    PIN("f3.coupled.c1",
        "flags=..CC|.CC.. cidx=[0:2->1,0:3->2,1:1->1,1:2->2]"
        " gflags=......C...CC............"
        " findCoupled=[2,3],[1,2] findCoupledPairs=00[],01[2,3],10[1,2],11[]"
        " findFreeUncoupled=[1],[0,3]",
        dumpCoupledQueries(m, 1));

    PIN("f3.coupled.c2",
        "flags=....|..... cidx=[] gflags=......C...CC............"
        " findCoupled=[],[] findCoupledPairs=00[],01[],10[],11[]"
        " findFreeUncoupled=[0,1,2,3],[0,1,2,4]",
        dumpCoupledQueries(m, 2));
}

// is_coupled_index, on its own and over the whole global range: exactly the
// three dofs that really are shared between the two patches, and nothing
// else.  Taking the cumulative m_numCpldDofs[c+1] as the band width instead
// of the component's own count widens component 1's band by component 0's
// one coupled dof and component 2's by all three, so free-uncoupled dofs of
// the later components are reported coupled.
TEST(is_coupled_index_uses_the_component_own_coupled_count)
{
    const gsDofMapper m = threeCompUniform();

    // Derived independently of the mapper's own bands: a coupled dof is a
    // free dof with more than one pre-image.
    std::vector<std::pair<index_t,index_t> > pre;
    for (index_t gl = 0; gl != m.size(); ++gl)
    {
        m.preImage(gl, pre);
        const bool shared = m.is_free_index(gl) && pre.size() > 1;
        CHECK_EQUAL(shared, m.is_coupled_index(gl));
    }

    CHECK_EQUAL(3, m.coupledSize());
}

TEST(three_comp_find_boundary_free_per_component)
{
    const gsDofMapper m = threeCompUniform();
    PIN("f3.findbf.c0", "bnd=[0],[] free=[1,2,3],[0,1,2,3,4]",   dumpFindBoundaryFree(m, 0));
    PIN("f3.findbf.c1", "bnd=[0],[4] free=[1,2,3],[0,1,2,3]",    dumpFindBoundaryFree(m, 1));
    PIN("f3.findbf.c2", "bnd=[],[3] free=[0,1,2,3],[0,1,2,4]",   dumpFindBoundaryFree(m, 2));
}

// =========================================================================
// firstFreeIndex(c)
// =========================================================================

// firstFreeIndex(c) and firstFreeIndex(c+1) bound the free indices of
// component c.  Checked against the free indices index(i,k,c) actually
// hands out, which is an independent oracle: it never consults the count
// vectors.  The block must contain every one of them and nothing else,
// and the last bound is the end of the free range, lastIndex().
namespace {
void checkFirstFreeIndexAgainstFreeIndices(const gsDofMapper & m)
{
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        std::set<index_t> free;
        for (index_t k = 0; k != nPatches(m); ++k)
        {
            const index_t n = static_cast<index_t>(m.patchSize(k, c));
            for (index_t i = 0; i != n; ++i)
            {
                const index_t gl = m.index(i, k, c);
                if (m.is_free_index(gl))
                    free.insert(gl);
            }
        }
        const index_t first = m.firstFreeIndex(c), end = m.firstFreeIndex(c+1);
        CHECK_EQUAL(static_cast<index_t>(free.size()), end - first);
        if (!free.empty())
        {
            CHECK_EQUAL(first,   *free.begin());
            CHECK_EQUAL(end - 1, *free.rbegin());
        }
    }
    CHECK_EQUAL(m.lastIndex(), m.firstFreeIndex(m.numComponents()));
}
} // anonymous namespace

TEST(first_free_index_bounds_the_free_block_of_each_component)
{
    checkFirstFreeIndexAgainstFreeIndices(threeCompUniform());
    checkFirstFreeIndexAgainstFreeIndices(threeCompShifted());
    checkFirstFreeIndexAgainstFreeIndices(identityMapper());
    checkFirstFreeIndexAgainstFreeIndices(raggedPatchMapper());

    // The literals for F3.  m_numFreeDofs[c]+m_numElimDofs[c] -- the
    // pre-renumbering base -- overshoots them by the eliminated dofs of the
    // earlier components once finalize() has moved the eliminated blocks
    // above all the free ones.
    const gsDofMapper m = threeCompUniform();
    CHECK_EQUAL(0,  m.firstFreeIndex(0));
    CHECK_EQUAL(7,  m.firstFreeIndex(1));
    CHECK_EQUAL(12, m.firstFreeIndex(2));
}

// A component that owns no dof at all has an empty free block, at the
// position where its dofs would have been numbered.
TEST(first_free_index_of_an_empty_component)
{
    std::vector<gsVector<index_t> > sz(3);
    sz[0].resize(1); sz[0][0] = 2;
    sz[1].resize(1); sz[1][0] = 0;   // empty component
    sz[2].resize(1); sz[2][0] = 3;
    gsDofMapper m(sz, /*hasDistinctComponentSpaces=*/false);
    m.finalize();

    CHECK_EQUAL(5, m.size());
    CHECK_EQUAL(5, m.freeSize());
    CHECK_EQUAL(0, m.boundarySize());
    CHECK_EQUAL(0, m.size(1));

    CHECK_EQUAL(0, m.firstFreeIndex(0));
    CHECK_EQUAL(2, m.firstFreeIndex(1));   // its (empty) free block
    CHECK_EQUAL(2, m.firstFreeIndex(2));
    checkFirstFreeIndexAgainstFreeIndices(m);
}

// A component whose dofs are all eliminated has an empty free block too.
// For component 0 that block sits at the start of the free range, so
// firstFreeIndex() is still the shift, not the start of component 0's
// eliminated block (which lies above every free block).
TEST(first_free_index_of_an_eliminated_only_component)
{
    gsVector<index_t> sz(1);
    sz[0] = 3;
    gsDofMapper m(sz, 2);
    for (index_t i = 0; i != 3; ++i)
        m.eliminateDof(i, 0, 0);
    m.finalize();

    CHECK_EQUAL(3, m.freeSize());
    CHECK_EQUAL(3, m.boundarySize());
    CHECK_EQUAL(0, m.freeSize(0));
    checkFirstFreeIndexAgainstFreeIndices(m);

    CHECK_EQUAL(0, m.firstFreeIndex(0));   // empty, at the start
    CHECK_EQUAL(0, m.firstFreeIndex(1));   // component 1's free block
    CHECK_EQUAL(3, m.firstFreeIndex(2));

    m.setShift(10);
    CHECK_EQUAL(10, m.firstFreeIndex());
    CHECK_EQUAL(13, m.lastIndex());
}

// Before finalize() the per-component free counts are not yet accumulated,
// so only component 0, whose block starts at the shift, has a defined start.
// Asking for a later component throws instead of returning the raw count.
TEST(first_free_index_before_finalize)
{
    CHECK_EQUAL(0, gsDofMapper().firstFreeIndex());

    gsVector<index_t> sz(2);
    sz << 3, 4;
    gsDofMapper m(sz, 2);
    m.setShift(10);
    CHECK_EQUAL(10, m.firstFreeIndex(0));
    CHECK_THROW(m.firstFreeIndex(1), std::runtime_error);
    CHECK_THROW(m.firstFreeIndex(2), std::runtime_error);

    m.finalize();
    CHECK_EQUAL(10, m.firstFreeIndex(0));
    CHECK_EQUAL(17, m.firstFreeIndex(1));
    CHECK_EQUAL(24, m.firstFreeIndex(2));
}

// The deprecated firstIndex(c) is firstFreeIndex(c).  In particular
// firstIndex() is the shift even when component 0 has no free dof:
// callers subtract it from a free index to get a position in the free
// range.
TEST(deprecated_first_index_is_first_free_index)
{
    gsVector<index_t> sz(1);
    sz[0] = 3;
    gsDofMapper m(sz, 2);
    for (index_t i = 0; i != 3; ++i)
        m.eliminateDof(i, 0, 0);
    m.finalize();
    m.setShift(10);

    CHECK_EQUAL(10, m.firstIndex());
    for (index_t c = 0; c <= m.numComponents(); ++c)
        CHECK_EQUAL(m.firstFreeIndex(c), m.firstIndex(c));
    const gsDofMapper u = threeCompUniform();
    for (index_t c = 0; c <= u.numComponents(); ++c)
        CHECK_EQUAL(u.firstFreeIndex(c), u.firstIndex(c));
}

// =========================================================================
// inverseAsVector(c)
// =========================================================================

// The inverse of asVector(c) over the whole global index space: every
// position that is not the image of a local dof of component c carries the
// -1 sentinel instead of whatever the allocation happened to hold.
TEST(inverse_as_vector_marks_other_components)
{
    gsVector<index_t> sz(2);
    sz[0] = 3; sz[1] = 4;
    gsDofMapper m(sz, 2);
    m.finalize();
    CHECK(m.isPermutation());

    PIN("inv.c0", "[0,1,2,3,4,5,6,-1,-1,-1,-1,-1,-1,-1]", join(m.inverseAsVector(0)));
    PIN("inv.c1", "[-1,-1,-1,-1,-1,-1,-1,0,1,2,3,4,5,6]", join(m.inverseAsVector(1)));

    // It really is the inverse of asVector on the component's own block.
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        const gsVector<index_t> fwd = m.asVector(c);
        const gsVector<index_t> inv = m.inverseAsVector(c);
        CHECK_EQUAL(m.size(), inv.size());
        for (index_t j = 0; j != fwd.size(); ++j)
            CHECK_EQUAL(j, inv[fwd[j]]);
    }
}

// =========================================================================
// inverseOnPatch(k)
// =========================================================================

// The inverse on a patch covers exactly the dofs that live on that patch,
// in every component.  Iterating a whole component vector from the patch
// offset instead attributes the following patches' dofs to the requested
// one and, for every patch but the first, reads past the end of the vector.
TEST(inverse_on_patch_is_bounded_to_the_patch)
{
    const gsDofMapper m = threeCompUniform();

    PIN("f3.inverseonpatch.k0",
        "[(0,1),(1,2),(6,3),(7,1),(10,2),(11,3),(12,0),(13,1),(14,2),(15,3),(20,0),(21,0)]",
        dumpInverseOnPatch(m, 0));
    PIN("f3.inverseonpatch.k1",
        "[(2,1),(3,2),(4,3),(5,4),(6,0),(8,0),(9,3),(10,1),(11,2),(16,0),(17,1),(18,2),(19,4),(22,4),(23,3)]",
        dumpInverseOnPatch(m, 1));

    // Independent oracle: the result is exactly { index(i,k,c) -> i }.
    for (index_t k = 0; k != nPatches(m); ++k)
    {
        const std::map<index_t,index_t> inv = m.inverseOnPatch(k);
        size_t expected = 0;
        for (index_t c = 0; c != m.numComponents(); ++c)
        {
            const index_t n = static_cast<index_t>(m.patchSize(k, c));
            expected += n;
            for (index_t i = 0; i != n; ++i)
            {
                const std::map<index_t,index_t>::const_iterator it =
                    inv.find(m.index(i, k, c));
                CHECK(it != inv.end());
                if (it != inv.end()) CHECK_EQUAL(i, it->second);
            }
        }
        // Coupled dofs appear once per local dof they carry on this patch,
        // so the map may be smaller than the total, never larger.
        CHECK(inv.size() <= expected);
    }
}

// Under the aliased layout every patch carries the complete inverse, so all
// patches return identical contents.
TEST(inverse_on_patch_is_aliased_under_global_identity)
{
    const gsDofMapper m = raggedIdentityMapper();
    const std::string k0 = dumpInverseOnPatch(m, 0);
    for (index_t k = 1; k != nPatches(m); ++k)
        PIN("f10.inverseonpatch", k0, dumpInverseOnPatch(m, k));
    CHECK_EQUAL(static_cast<size_t>(m.size()), m.inverseOnPatch(0).size());
}

// =========================================================================
// anyPreImages(c)
// =========================================================================

// anyPreImages indexes its result by the global dof value, which for every
// component but the first is larger than that component's own storage size:
// sizing the result by that storage size is an out-of-bounds write, not
// merely a wrong answer.  The result has exactly one entry per global index,
// and an index of another component carries (-1,-1) -- a -1 in the second
// slot as well, since 0 is a valid patch-local index.
namespace {
void checkAnyPreImagesAgainstAnyPreImage(const gsDofMapper & m, const index_t shift)
{
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        const std::vector<std::pair<index_t,index_t> > all = m.anyPreImages(c);
        CHECK_EQUAL(static_cast<size_t>(m.size()), all.size());

        for (index_t g = 0; g != m.size(); ++g)
        {
            if (m.componentOf(shift + g) == c)
                CHECK(m.anyPreImage(shift + g) == all[g]);
            else
                CHECK(std::make_pair(index_t(-1), index_t(-1)) == all[g]);
        }
    }
}
} // anonymous namespace

TEST(any_pre_images_agrees_with_any_pre_image)
{
    checkAnyPreImagesAgainstAnyPreImage(threeCompUniform(), 0);
    // Positions are unshifted, gl minus the shift, like inverseAsVector's.
    checkAnyPreImagesAgainstAnyPreImage(threeCompShifted(), 100);
}

// =========================================================================
// permuteFreeDofs(permutation, c)
// =========================================================================

// The permutation is component-local: it has one entry per free dof of the
// component it is applied to, and it may not move an index out of that
// component's free block.  Taking its length from the cumulative
// m_numFreeDofs[c+1] both demands a longer permutation than the component
// has free dofs and, indexed by an unrebased global index, would permute
// component c's dofs into its predecessors' bands.
TEST(permute_free_dofs_is_component_local)
{
    gsDofMapper m = threeCompUniform();
    CHECK_EQUAL(5, m.freeSize(1));

    gsVector<index_t> perm(5);
    perm << 4, 3, 2, 1, 0;
    m.permuteFreeDofs(perm, 1);

    // Component 1's free block is [7,12): every free index stays inside it,
    // and the two eliminated dofs (21 and 22) are untouched.
    PIN("perm.c1", "[21,11,8,7,10,8,7,9,22]", join(m.asVector(1)));

    // The other components are not touched at all.
    PIN("perm.c0", "[20,0,1,6,6,2,3,4,5]",          join(m.asVector(0)));
    PIN("perm.c2", "[12,13,14,15,16,17,18,23,19]",  join(m.asVector(2)));

    // Permuting a component destroys the tracking of ITS coupled dofs only;
    // component 0's single coupled dof is still counted.
    CHECK_EQUAL(1, m.coupledSize());
}

// markCoupledAsTagged tags the top of each component's own free block, which
// is where finalize() puts that component's coupled dofs.  Taking the band
// start from m_numFreeDofs[c+1]+m_numElimDofs[c] instead walks off the
// coupled dofs entirely -- the eliminated blocks are stacked above every
// component's free block, not interleaved with them -- and taking the band
// WIDTH from the cumulative m_numCpldDofs[c+1] then overruns, for the last
// component past size() itself.
TEST(mark_coupled_as_tagged_tags_the_coupled_dofs)
{
    gsDofMapper m = threeCompUniform();
    CHECK_EQUAL(0, m.taggedSize());       // nothing tagged beforehand
    m.markCoupledAsTagged();

    PIN("tagcpld.f3", "[6,10,11]", join(m.getTagged()));
    CHECK_EQUAL(m.coupledSize(), m.taggedSize());

    // The tagged set is now exactly the coupled set, component by component.
    for (index_t gl = 0; gl != m.size(); ++gl)
        CHECK_EQUAL(m.is_coupled_index(gl), m.is_tagged_index(gl));

    // Single component, one collapsed group: the one coupled dof, a valid
    // index below size() == 8.
    gsDofMapper c = collapsedMapper();
    c.markCoupledAsTagged();
    PIN("tagcpld.f7", "[7]", join(c.getTagged()));
    CHECK(c.getTagged().back() < c.size());
}

// The coupled dofs are added to whatever is already tagged, not substituted
// for it.
TEST(mark_coupled_as_tagged_keeps_existing_tags)
{
    gsDofMapper m = twoPatchCoupledElim();   // markTagged'd dofs 0 and 3
    PIN("tagcpld.before", "[0,3]", join(m.getTagged()));
    m.markCoupledAsTagged();
    PIN("tagcpld.after", "[0,3,5,6]", join(m.getTagged()));
}

// A component-local permutation moves only the tags that lie in that
// component's free block.  Rebuilding the tag list from the permuted
// component's storage -- which is the only storage the permutation walks --
// visits no other component, so every tag of every other component is
// silently deleted.  That breaks the documented "markCoupledAsTagged(), then
// permute, then use the tagged queries" workflow, which is the whole reason
// permuting is allowed to destroy the coupled bookkeeping in the first place.
TEST(permute_free_dofs_keeps_other_components_tags)
{
    gsDofMapper m = threeCompUniform();

    m.markTagged(1, 0, 0);   // component 0, global 0
    m.markTagged(1, 0, 1);   // component 1, global 7  -- the permuted one
    m.markTagged(1, 0, 2);   // component 2, global 13
    PIN("permtag.before", "[0,7,13]", join(m.getTagged()));

    gsVector<index_t> perm(5);
    perm << 4, 3, 2, 1, 0;
    m.permuteFreeDofs(perm, 1);

    // Component 1's free block is [7,12), so its tag 7 becomes 7+perm[0]==11.
    // The other two are outside the block and must be untouched.
    PIN("permtag.after", "[0,11,13]", join(m.getTagged()));
    CHECK_EQUAL(3, m.taggedSize());
    CHECK(m.is_tagged_index(0));
    CHECK(m.is_tagged_index(11));
    CHECK(m.is_tagged_index(13));
    CHECK(!m.is_tagged_index(7));
}

// An eliminated dof's tag is not in any free block, so it survives a
// permutation of its own component unchanged.
TEST(permute_free_dofs_keeps_eliminated_tags)
{
    gsDofMapper m = threeCompUniform();
    m.markTagged(0, 0, 0);   // component 0's eliminated dof, global 20

    gsVector<index_t> perm(7);
    perm << 3, 0, 6, 1, 5, 2, 4;
    m.permuteFreeDofs(perm, 0);

    PIN("permtag.elim", "[20]", join(m.getTagged()));
}

// =========================================================================
// Global shift
// =========================================================================
//
// Everything a mapper stores -- the dof numbering and the tags -- lives in
// the unshifted index space.  The shift is added to the global indices a
// query returns and subtracted from the global indices a query takes, and
// nothing else sees it.  A shifted mapper must therefore answer every query
// exactly as its unshifted original does, once the global indices on either
// side of the call are translated.  Multi-space assemblers shift every space
// but the first (gsExprAssembler::resetDimensions), so the translation is
// live, not theoretical.

namespace {
// Every query of \a s against \a m, translated by \a shift.  \a s must be
// \a m with setShift(shift) applied and nothing else.
void checkShiftIsARelabelling(const gsDofMapper & m, const gsDofMapper & s,
                              const index_t shift)
{
    CHECK_EQUAL(m.size(), s.size());
    std::vector<std::pair<index_t,index_t> > pm, ps;
    for (index_t g = 0; g != m.size(); ++g)
    {
        const index_t gl = shift + g;
        CHECK_EQUAL(m.is_free_index(g),     s.is_free_index(gl));
        CHECK_EQUAL(m.is_boundary_index(g), s.is_boundary_index(gl));
        CHECK_EQUAL(m.is_coupled_index(g),  s.is_coupled_index(gl));
        CHECK_EQUAL(m.is_tagged_index(g),   s.is_tagged_index(gl));
        CHECK_EQUAL(m.componentOf(g),       s.componentOf(gl));
        CHECK(m.anyPreImage(g) == s.anyPreImage(gl));
        m.preImage(g, pm);
        s.preImage(gl, ps);
        CHECK(pm == ps);
        for (index_t k = 0; k != nPatches(m); ++k)
        {
            index_t lm = -1, ls = -2;
            CHECK_EQUAL(m.indexOnPatch(g, k, lm), s.indexOnPatch(gl, k, ls));
            if (m.indexOnPatch(g, k)) CHECK_EQUAL(lm, ls);
        }
    }

    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        CHECK(m.anyPreImages(c) == s.anyPreImages(c));
        CHECK_EQUAL(m.firstFreeIndex(c) + shift, s.firstFreeIndex(c));
        for (index_t k = 0; k != nPatches(m); ++k)
        {
            CHECK(m.findBoundary(k, c)      == s.findBoundary(k, c));
            CHECK(m.findFree(k, c)          == s.findFree(k, c));
            CHECK(m.findFreeUncoupled(k, c) == s.findFreeUncoupled(k, c));
            for (index_t j = -1; j != nPatches(m); ++j)
                CHECK(m.findCoupled(k, j, c) == s.findCoupled(k, j, c));

            const index_t n = static_cast<index_t>(m.patchSize(k, c));
            for (index_t i = 0; i != n; ++i)
            {
                CHECK_EQUAL(m.index(i, k, c) + shift, s.index(i, k, c));
                CHECK_EQUAL(m.tindex(i, k, c),        s.tindex(i, k, c));
                CHECK_EQUAL(m.is_tagged(i, k, c),     s.is_tagged(i, k, c));
                CHECK_EQUAL(m.is_coupled(i, k, c),    s.is_coupled(i, k, c));
            }
        }
    }

    // inverseOnPatch is keyed by global index, so its keys carry the shift.
    for (index_t k = 0; k != nPatches(m); ++k)
    {
        const std::map<index_t,index_t> im = m.inverseOnPatch(k);
        const std::map<index_t,index_t> is = s.inverseOnPatch(k);
        CHECK_EQUAL(im.size(), is.size());
        for (std::map<index_t,index_t>::const_iterator it = im.begin();
             it != im.end(); ++it)
        {
            const std::map<index_t,index_t>::const_iterator jt = is.find(it->first + shift);
            CHECK(jt != is.end());
            if (jt != is.end()) CHECK_EQUAL(it->second, jt->second);
        }
    }
}
} // anonymous namespace

TEST(shift_is_a_relabelling)
{
    // Coupled, eliminated, three components.
    checkShiftIsARelabelling(threeCompUniform(), threeCompShifted(), 100);

    // Tags and a single component.
    gsDofMapper f2 = twoPatchCoupledElim();
    gsDofMapper f2s = f2;
    f2s.setShift(40);
    checkShiftIsARelabelling(f2, f2s, 40);

    // The aliased layout.
    gsDofMapper f5 = identityMapper();
    gsDofMapper f5s = f5;
    f5s.setShift(7);
    checkShiftIsARelabelling(f5, f5s, 7);
}

// markTagged, is_tagged_index, tindex and markCoupledAsTagged all agree on a
// shifted mapper: every tag is stored unshifted, so a markTagged'd dof is
// found by tindex and both ways of tagging produce lists in one index space.
TEST(tags_are_stored_unshifted)
{
    gsDofMapper m = threeCompShifted();
    m.markTagged(1, 0, 1);            // component 1, freeIndex 7, index 107
    PIN("shift.tagged", "[7]", join(m.getTagged()));
    CHECK(m.is_tagged(1, 0, 1));
    CHECK(m.is_tagged_index(107));
    CHECK(!m.is_tagged_index(7));     // an unshifted index is not a global one
    CHECK_EQUAL(0, m.tindex(1, 0, 1));

    // The tag stays with its dof when the shift changes afterwards.
    m.setShift(250);
    CHECK(m.is_tagged(1, 0, 1));
    CHECK(m.is_tagged_index(257));
    CHECK(!m.is_tagged_index(107));

    // markCoupledAsTagged on a shifted mapper tags exactly the coupled dofs,
    // in the same index space as markTagged.
    gsDofMapper c = threeCompShifted();
    c.markCoupledAsTagged();
    for (index_t gl = 100; gl != 100 + c.size(); ++gl)
        CHECK_EQUAL(c.is_coupled_index(gl), c.is_tagged_index(gl));
    PIN("shift.tagcpld", "[6,10,11]", join(c.getTagged()));
}

// permuteFreeDofs classifies the stored, unshifted values of the permuted
// component.  Classifying them with is_free_index, which expects a shifted
// index, counts the stored values of the first `shift` eliminated dofs as
// free and tries to permute them.  On a shifted mapper the whole F6 workflow
// must yield F6 plus the shift.
TEST(permute_free_dofs_on_a_shifted_mapper)
{
    gsDofMapper m = threeCompShifted();
    m.markCoupledAsTagged();
    gsVector<index_t> perm(7);
    perm << 3, 0, 6, 1, 5, 2, 4;
    m.permuteFreeDofs(perm, 0);

    const gsDofMapper f6 = permutedMapper();
    for (index_t c = 0; c != m.numComponents(); ++c)
    {
        const gsVector<index_t> expected = (f6.asVector(c).array() + 100).matrix();
        CHECK(expected == m.asVector(c));
    }
    CHECK(f6.getTagged() == m.getTagged());
    CHECK_EQUAL(f6.coupledSize(), m.coupledSize());
    PIN("shift.perm.tagged", dumpTagged(f6), dumpTagged(m, 100));
}

// A global index outside [shift, shift+size()) belongs to no dof of the
// mapper: the predicates answer false, and the queries that must return a
// component or a preimage throw.
TEST(out_of_range_global_indices)
{
    const gsDofMapper m = threeCompShifted();   // global range [100,124)
    const index_t outside[] = {-1, 0, 23, 99, 124, 1000};
    for (size_t t = 0; t != sizeof(outside)/sizeof(outside[0]); ++t)
    {
        const index_t gl = outside[t];
        CHECK(!m.is_free_index(gl));
        CHECK(!m.is_boundary_index(gl));
        CHECK(!m.is_coupled_index(gl));
        CHECK(!m.is_tagged_index(gl));
        for (index_t k = 0; k != nPatches(m); ++k)
            CHECK(!m.indexOnPatch(gl, k));
        std::vector<std::pair<index_t,index_t> > pre;
        // GISMO_ENSURE, so these throw in Release builds too.
        CHECK_THROW(m.componentOf(gl), std::runtime_error);
        CHECK_THROW(m.preImage(gl, pre), std::runtime_error);
        CHECK_THROW(m.anyPreImage(gl), std::runtime_error);
    }

    // The two ends of the range are inside it.
    CHECK(m.is_free_index(100));
    CHECK(m.is_boundary_index(123));
    CHECK_EQUAL(0, m.componentOf(100));
    CHECK_EQUAL(2, m.componentOf(123));
}

// The range checks must hold for every index_t value and every shift,
// including a negative one: computing gl - shift before checking the range
// is signed overflow for gl near either end of index_t, which is undefined
// behaviour rather than a wrong answer, and so is caught by the
// undefined-behaviour sanitizer rather than by the assertions below.
namespace {
void checkExtremesAreOutOfRange(const gsDofMapper & m, const index_t shift)
{
    const index_t lo = std::numeric_limits<index_t>::min();
    const index_t hi = std::numeric_limits<index_t>::max();
    const index_t probes[] = {lo, lo + 1, hi - 1, hi};
    for (size_t t = 0; t != sizeof(probes)/sizeof(probes[0]); ++t)
    {
        const index_t gl = probes[t];
        // Every probe is out of range unless the shift puts it inside.
        // Exact because every shift used here keeps shift+size()
        // representable, as index() itself requires.
        const bool inside = gl >= shift && gl < shift + m.size();
        if (inside) continue;
        CHECK(!m.is_free_index(gl));
        CHECK(!m.is_boundary_index(gl));
        CHECK(!m.is_coupled_index(gl));
        CHECK(!m.is_tagged_index(gl));
        for (index_t k = 0; k != nPatches(m); ++k)
            CHECK(!m.indexOnPatch(gl, k));
        CHECK_THROW(m.componentOf(gl), std::runtime_error);
    }
}
} // anonymous namespace

TEST(extreme_global_indices_and_shifts)
{
    const index_t lo = std::numeric_limits<index_t>::min();
    const index_t hi = std::numeric_limits<index_t>::max();
    const gsDofMapper f3 = threeCompUniform();
    const index_t n = f3.size();

    // Small positive and negative shifts, probed at both ends of index_t.
    const index_t shifts[] = {1, -1, -50, 0};
    for (size_t t = 0; t != sizeof(shifts)/sizeof(shifts[0]); ++t)
    {
        gsDofMapper s = f3;
        s.setShift(shifts[t]);
        checkExtremesAreOutOfRange(s, shifts[t]);
        checkShiftIsARelabelling(f3, s, shifts[t]);
    }

    // The largest shift for which every global index is representable: the
    // last dof sits at index_t max - 1 and max itself is one past the end.
    {
        gsDofMapper s = f3;
        s.setShift(hi - n);
        checkExtremesAreOutOfRange(s, hi - n);
        checkShiftIsARelabelling(f3, s, hi - n);
        CHECK(s.is_boundary_index(hi - 1));
        CHECK_EQUAL(2, s.componentOf(hi - 1));
        CHECK(!s.is_boundary_index(hi));
    }

    // The smallest shift: the first dof is index_t min itself.
    {
        gsDofMapper s = f3;
        s.setShift(lo);
        checkExtremesAreOutOfRange(s, lo);
        checkShiftIsARelabelling(f3, s, lo);
        CHECK(s.is_free_index(lo));
        CHECK_EQUAL(0, s.componentOf(lo));
    }
}

// =========================================================================
// Release-safe argument validation
// =========================================================================
//
// Every public entry point that would otherwise read or write outside the
// mapper's storage validates its arguments with GISMO_ENSURE, which throws
// std::runtime_error in Release builds as well -- hence plain CHECK_THROW
// throughout.  The exception is the assembly hot path -- the per-dof
// accessors index(), bindex(), cindex(), tindex(), freeIndex() and the
// per-element localToGlobal()/localToGlobal2() -- which check the same
// bounds in debug builds only (see hot_path_accessors_... and
// local_to_global_...).
//
// Batch and broadcast calls are validated in full before the first change,
// so a call that throws must leave the mapper exactly as it was.  That is
// checked by comparing digests of everything observable before and after;
// for a mapper still in setup, a finalized copy stands in for it.

namespace {

std::string finalDigest(const gsDofMapper & m)
{
    return dumpCounts(m) + "|" + dumpPerComponent(m) + "|" + dumpAsVector(m)
        + "|" + join(m.getTagged()) + "|" + dumpFirstLast(m);
}

std::string setupDigest(gsDofMapper m)
{
    m.finalize();
    return finalDigest(m);
}

gsMatrix<index_t> column(index_t a, index_t b)
{
    gsMatrix<index_t> m(2, 1);
    m << a, b;
    return m;
}

gsMatrix<unsigned> ucolumn(unsigned a, unsigned b)
{
    gsMatrix<unsigned> m(2, 1);
    m << a, b;
    return m;
}

} // anonymous namespace

// A local index must lie in its own patch's range, not merely somewhere in
// the component's storage: a bound taken from the component total lets an
// oversized index of one patch silently address the next patch's dofs.
TEST(oversized_local_index_is_rejected_per_patch)
{
    gsDofMapper m = raggedSetup();
    const std::string before = setupDigest(m);

    // Local 2 of (patch 0, component 1) would be local 0 of patch 1.
    CHECK_THROW(m.eliminateDof(2, 0, 1), std::runtime_error);
    CHECK_THROW(m.markCoupled(2, 0, 1), std::runtime_error);
    CHECK_THROW(m.matchDof(0, 2, 1, 3, 1), std::runtime_error);
    CHECK_THROW(m.matchDof(1, 3, 0, 2, 1), std::runtime_error);
    CHECK_THROW(m.markBoundary(0, column(0, 2), 1), std::runtime_error);
    CHECK_THROW(m.colapseDofs(0, ucolumn(0, 2), 1), std::runtime_error);
    // Local 5 of (patch 0, component 0) would be local 0 of patch 1.
    CHECK_THROW(m.eliminateDof(5, 0, 0), std::runtime_error);
    // One past the last patch, i.e. past the end of the storage.
    CHECK_THROW(m.eliminateDof(4, 1, 1), std::runtime_error);
    CHECK_THROW(m.eliminateDof(5, 1, 0), std::runtime_error);
    // Negative.
    CHECK_THROW(m.eliminateDof(-1, 0, 0), std::runtime_error);
    CHECK_THROW(m.matchDof(0, 0, 1, -1, 0), std::runtime_error);

    CHECK_EQUAL(before, setupDigest(m));

    // The last valid local index of every (patch,component) is accepted.
    m.eliminateDof(4, 0, 0);
    m.eliminateDof(4, 1, 0);
    m.eliminateDof(1, 0, 1);
    m.eliminateDof(3, 1, 1);
    m.finalize();
    CHECK_EQUAL(4, m.boundarySize());
}

// The Raviart-Thomas pair in both mesh flavours.  On the anisotropic mesh
// the local bound differs between the components (28 and 30 on one patch);
// on the isotropic mesh it coincides (42 and 42).
TEST(rt_local_bounds_are_per_component)
{
    {
        gsDofMapper m(rtSizes(4, 2), true);
        CHECK_EQUAL(28u, (unsigned)m.patchSize(0, 0));
        CHECK_EQUAL(30u, (unsigned)m.patchSize(0, 1));
        const std::string before = setupDigest(m);
        CHECK_THROW(m.eliminateDof(28, 0, 0), std::runtime_error);
        CHECK_THROW(m.eliminateDof(28, 0, -1), std::runtime_error);
        CHECK_THROW(m.eliminateDof(30, 0, 1), std::runtime_error);
        CHECK_EQUAL(before, setupDigest(m));
        m.eliminateDof(29, 0, 1);
        m.eliminateDof(27, 0, -1);
        m.finalize();
        CHECK_EQUAL(3, m.boundarySize());
    }
    {
        gsDofMapper m(rtSizes(4, 4), true);
        CHECK_EQUAL(42u, (unsigned)m.patchSize(0, 0));
        CHECK_EQUAL(42u, (unsigned)m.patchSize(0, 1));
        const std::string before = setupDigest(m);
        CHECK_THROW(m.eliminateDof(42, 0, -1), std::runtime_error);
        CHECK_THROW(m.eliminateDof(42, 0, 1), std::runtime_error);
        CHECK_EQUAL(before, setupDigest(m));
        m.eliminateDof(41, 0, -1);
        m.finalize();
        CHECK_EQUAL(2, m.boundarySize());
    }
}

// Under the aliased layout the local index is component-global: the bound
// is that component's own total, on every patch alike.
TEST(identity_local_bound_is_the_component_total)
{
    std::vector<size_t> dofs(2);
    dofs[0] = 7; dofs[1] = 10;
    gsDofMapper m;
    m.setIdentity(3, dofs);
    const std::string before = setupDigest(m);

    CHECK_THROW(m.eliminateDof(7, 0, 0), std::runtime_error);
    CHECK_THROW(m.eliminateDof(7, 2, 0), std::runtime_error);
    CHECK_THROW(m.eliminateDof(10, 1, 1), std::runtime_error);
    CHECK_THROW(m.eliminateDof(9, 0, -1), std::runtime_error);  // 9 >= 7
    CHECK_THROW(m.eliminateDof(0, 3, 1), std::runtime_error);   // no patch 3
    CHECK_EQUAL(before, setupDigest(m));

    m.eliminateDof(6, 0, 0);
    m.eliminateDof(9, 2, 1);    // accepted on any patch
    m.eliminateDof(9, 0, 1);    // the same dof again: nothing new
    m.finalize();
    CHECK_EQUAL(2, m.boundarySize());
}

TEST(invalid_component_and_patch_identifiers_in_setup)
{
    gsDofMapper m = raggedSetup();       // 2 components, 2 patches
    const std::string before = setupDigest(m);
    const index_t badComp[]  = {-2, 2, 100};
    const index_t badPatch[] = {-1, 2, 100};

    for (size_t t = 0; t != 3; ++t)
    {
        const index_t c = badComp[t], k = badPatch[t];
        CHECK_THROW(m.eliminateDof(0, 0, c), std::runtime_error);
        CHECK_THROW(m.eliminateDof(0, k, 0), std::runtime_error);
        CHECK_THROW(m.markCoupled(0, 0, c), std::runtime_error);
        CHECK_THROW(m.markCoupled(0, k, 0), std::runtime_error);
        CHECK_THROW(m.matchDof(0, 0, 1, 0, c), std::runtime_error);
        CHECK_THROW(m.matchDof(k, 0, 1, 0, 0), std::runtime_error);
        CHECK_THROW(m.matchDof(0, 0, k, 0, 0), std::runtime_error);
        CHECK_THROW(m.matchDofs(0, column(0,1), 1, column(0,1), c), std::runtime_error);
        CHECK_THROW(m.matchDofs(k, column(0,1), 1, column(0,1), 0), std::runtime_error);
        CHECK_THROW(m.matchDofs(0, column(0,1), k, column(0,1), 0), std::runtime_error);
        CHECK_THROW(m.markBoundary(0, column(0,1), c), std::runtime_error);
        CHECK_THROW(m.markBoundary(k, column(0,1), 0), std::runtime_error);
        CHECK_THROW(m.colapseDofs(0, ucolumn(0,1), c), std::runtime_error);
        CHECK_THROW(m.colapseDofs(k, ucolumn(0,1), 0), std::runtime_error);
    }
    CHECK_EQUAL(before, setupDigest(m));
}

TEST(invalid_component_and_patch_identifiers_in_queries)
{
    gsDofMapper m = raggedPatchMapper();     // finalized, 2 components, 2 patches
    const index_t badComp[]  = {-1, 2, 100};
    const index_t badPatch[] = {-1, 2, 100};
    gsMatrix<index_t> loc(1, 1), glob;
    loc << 0;
    index_t nf = 0;
    GISMO_UNUSED(nf);   // only read by the debug-only checks

    for (size_t t = 0; t != 3; ++t)
    {
        const index_t c = badComp[t], k = badPatch[t];
        CHECK_THROW(m.patchSize(0, c), std::runtime_error);
        CHECK_THROW(m.patchSize(k, 0), std::runtime_error);
        CHECK_THROW(m.totalSize(c), std::runtime_error);
        CHECK_THROW(m.offset(0, c), std::runtime_error);
        CHECK_THROW(m.offset(k, 0), std::runtime_error);
        CHECK_THROW(m.size(c), std::runtime_error);
        CHECK_THROW(m.freeSize(c), std::runtime_error);
        CHECK_THROW(m.asVector(c), std::runtime_error);
        CHECK_THROW(m.inverseAsVector(c), std::runtime_error);
        CHECK_THROW(m.anyPreImages(c), std::runtime_error);
        CHECK_THROW(m.inverseOnPatch(k), std::runtime_error);
        CHECK_THROW(m.indexOnPatch(0, k), std::runtime_error);
        CHECK_THROW(m.findBoundary(k, 0), std::runtime_error);
        CHECK_THROW(m.findBoundary(0, c), std::runtime_error);
        CHECK_THROW(m.findFree(k, 0), std::runtime_error);
        CHECK_THROW(m.findFree(0, c), std::runtime_error);
        CHECK_THROW(m.findCoupled(k, -1, 0), std::runtime_error);
        CHECK_THROW(m.findCoupled(0, -1, c), std::runtime_error);
        CHECK_THROW(m.findFreeUncoupled(k, 0), std::runtime_error);
        CHECK_THROW(m.findFreeUncoupled(0, c), std::runtime_error);
        CHECK_THROW(m.findTagged(k, 0), std::runtime_error);
        CHECK_THROW(m.findTagged(0, c), std::runtime_error);
        CHECK_THROW(m.markTagged(0, 0, c), std::runtime_error);
        CHECK_THROW(m.markTagged(0, k, 0), std::runtime_error);
        // Assembly hot path: debug-checked only, like index().
        CHECK_THROW_IN_DEBUG(m.localToGlobal(loc, k, glob, 0), std::logic_error);
        CHECK_THROW_IN_DEBUG(m.localToGlobal(loc, 0, glob, c), std::logic_error);
        CHECK_THROW_IN_DEBUG(m.localToGlobal2(loc, k, glob, nf, 0), std::logic_error);
        CHECK_THROW_IN_DEBUG(m.localToGlobal2(loc, 0, glob, nf, c), std::logic_error);
    }
    // The second patch of findCoupled() is -1 (any) or a real patch.
    CHECK_THROW(m.findCoupled(0, 2, 0), std::runtime_error);
    CHECK_THROW(m.findCoupled(0, -2, 0), std::runtime_error);
    CHECK_EQUAL(0, m.findCoupled(0, 1, 0).size());
    // markTagged checks the local index against its own patch too.
    CHECK_THROW(m.markTagged(2, 0, 1), std::runtime_error);
    // mapIndex() takes a flat index into [0, mapSize()).
    CHECK_THROW(m.mapIndex(-1), std::runtime_error);
    CHECK_THROW(m.mapIndex(static_cast<index_t>(m.mapSize())), std::runtime_error);
    CHECK_EQUAL(15, m.mapIndex(static_cast<index_t>(m.mapSize()) - 1));
    // firstFreeIndex(c) also accepts c == numComponents(): the end of the
    // last free block.
    CHECK_EQUAL(16, m.firstFreeIndex(2));
    CHECK_THROW(m.firstFreeIndex(3), std::runtime_error);
    CHECK_THROW(m.firstFreeIndex(-1), std::runtime_error);
    CHECK_EQUAL(0, m.taggedSize());
}

// A default-constructed mapper has no components, so every
// component-indexed query throws; firstFreeIndex() alone stays valid on it.
TEST(default_constructed_component_queries)
{
    const gsDofMapper m;
    CHECK_EQUAL(0, m.firstFreeIndex());
    CHECK_THROW(m.freeSize(0), std::runtime_error);
    CHECK_THROW(m.totalSize(0), std::runtime_error);
    CHECK_THROW(m.patchSize(0, 0), std::runtime_error);
    CHECK_THROW(m.offset(0), std::runtime_error);
    CHECK_THROW(m.asVector(0), std::runtime_error);
    CHECK_THROW(m.mapIndex(0), std::runtime_error);
}

// A broadcast (component -1) is validated for every component before any of
// them is changed.  The discriminating case is a local index that is valid
// in an EARLIER component and invalid in a LATER one: checking component by
// component while mutating would eliminate or match in the earlier one and
// only then throw.
TEST(broadcast_is_validated_for_every_component_before_mutating)
{
    {
        gsDofMapper m = raggedSetup();   // local 3 of patch 0: comp 0 yes, comp 1 no
        const std::string before = setupDigest(m);
        CHECK_THROW(m.eliminateDof(3, 0, -1), std::runtime_error);
        CHECK_THROW(m.markCoupled(3, 0, -1), std::runtime_error);
        CHECK_THROW(m.matchDof(0, 3, 1, 0, -1), std::runtime_error);
        CHECK_THROW(m.matchDof(1, 0, 0, 3, -1), std::runtime_error);
        CHECK_THROW(m.matchDofs(0, column(0,3), 1, column(0,1), -1), std::runtime_error);
        CHECK_THROW(m.markBoundary(0, column(0,3), -1), std::runtime_error);
        CHECK_THROW(m.colapseDofs(0, ucolumn(0,3), -1), std::runtime_error);
        CHECK_EQUAL(before, setupDigest(m));
    }
    {
        // The RT pair with the mesh transposed puts the larger component
        // first: 30 and 28 local dofs.
        gsDofMapper m(rtSizes(2, 4), true);
        CHECK_EQUAL(30u, (unsigned)m.patchSize(0, 0));
        CHECK_EQUAL(28u, (unsigned)m.patchSize(0, 1));
        const std::string before = setupDigest(m);
        CHECK_THROW(m.eliminateDof(29, 0, -1), std::runtime_error);
        CHECK_THROW(m.markBoundary(0, column(0,28), -1), std::runtime_error);
        CHECK_EQUAL(before, setupDigest(m));
    }
}

// A batch is validated entry by entry before the first entry is applied: a
// bad entry at the END is the one that exposes a validate-while-mutating
// implementation.
TEST(batch_is_validated_before_mutating)
{
    gsDofMapper m = raggedSetup();
    const std::string before = setupDigest(m);

    gsMatrix<index_t> b1(3, 1), b2(3, 1), bad(3, 1);
    b1  << 0, 1, 2;
    b2  << 0, 1, 2;
    bad << 0, 1, 4;           // 4 is past (patch 0, component 1)
    CHECK_THROW(m.matchDofs(1, b2, 0, bad, 1), std::runtime_error);
    CHECK_THROW(m.matchDofs(0, bad, 1, b2, 1), std::runtime_error);
    CHECK_THROW(m.markBoundary(0, bad, 1), std::runtime_error);
    {
        gsMatrix<unsigned> ub(3, 1);
        ub << 0, 1, 4;
        CHECK_THROW(m.colapseDofs(0, ub, 1), std::runtime_error);
    }
    // Unsigned entries beyond index_t's range are rejected, not narrowed
    // onto a valid-looking index.
    CHECK_THROW(m.colapseDofs(0, ucolumn(0, std::numeric_limits<unsigned>::max()), 0),
                std::runtime_error);
    // Mismatched lengths, and more than one column (only the first column
    // is ever read).
    CHECK_THROW(m.matchDofs(0, b1, 1, column(0,1), 0), std::runtime_error);
    {
        gsMatrix<index_t> wide(1, 2);
        wide << 0, 1;
        CHECK_THROW(m.matchDofs(0, wide, 1, wide, 0), std::runtime_error);
        gsMatrix<unsigned> uwide(1, 2);
        uwide << 0, 1;
        CHECK_THROW(m.colapseDofs(0, uwide, 0), std::runtime_error);
    }
    // Rows but no column: there is no first column to read.  Only a matrix
    // without rows is empty.
    CHECK_THROW(m.colapseDofs(0, gsMatrix<unsigned>(2, 0), 0), std::runtime_error);
    CHECK_THROW(m.markBoundary(0, gsMatrix<index_t>(2, 0), 0), std::runtime_error);
    // markBoundary() reads rows() entries: a row vector is rejected rather
    // than silently cut to its first entry.
    {
        gsMatrix<index_t> row(1, 2);
        row << 0, 1;
        CHECK_THROW(m.markBoundary(0, row, 0), std::runtime_error);
    }
    CHECK_EQUAL(before, setupDigest(m));

    // Fewer than two dofs collapse to nothing, and an empty boundary
    // eliminates nothing.
    m.markBoundary(0, gsMatrix<index_t>(0, 1), 0);
    m.markBoundary(0, gsMatrix<index_t>(), 0);
    m.colapseDofs(0, gsMatrix<unsigned>(), 0);
    m.colapseDofs(0, gsMatrix<unsigned>(1, 1).setZero(), 0);
    m.matchDofs(0, gsMatrix<index_t>(), 1, gsMatrix<index_t>(), 0);
    CHECK_EQUAL(before, setupDigest(m));

    // The valid batches go through.
    m.matchDofs(0, b1, 1, b2, 0);
    m.markBoundary(1, b2, 1);
    m.finalize();
    CHECK_EQUAL(3, m.coupledSize());
    CHECK_EQUAL(3, m.boundarySize());
}

// The setup mutators rewrite the setup-time encoding that finalize()
// replaces by the final numbering, so they are rejected afterwards -- and
// finalize() itself runs exactly once.  The post-finalize mutators
// markTagged(), markCoupledAsTagged() and permuteFreeDofs() are the
// opposite: they need the final numbering and are rejected before it.
TEST(setup_mutators_require_an_unfinalized_mapper)
{
    gsDofMapper m = twoPatchCoupledElim();
    const std::string before = finalDigest(m);

    CHECK_THROW(m.matchDof(0, 0, 1, 5), std::runtime_error);
    CHECK_THROW(m.matchDofs(0, column(0,1), 1, column(4,5), 0), std::runtime_error);
    CHECK_THROW(m.markCoupled(2, 0), std::runtime_error);
    CHECK_THROW(m.eliminateDof(2, 0), std::runtime_error);
    CHECK_THROW(m.markBoundary(0, column(2,3), 0), std::runtime_error);
    CHECK_THROW(m.colapseDofs(0, ucolumn(2,3), 0), std::runtime_error);
    CHECK_THROW(m.finalize(), std::runtime_error);
    CHECK_EQUAL(before, finalDigest(m));

    // Legitimate after finalize().
    m.markTagged(2, 0);
    m.markCoupledAsTagged();
    gsVector<index_t> id(m.freeSize(0));
    for (index_t i = 0; i != id.size(); ++i) id[i] = i;
    m.permuteFreeDofs(id, 0);

    // The queries below read the final numbering, which an unfinalized
    // mapper does not have: its storage holds the setup encoding instead.
    gsDofMapper s = raggedSetup();
    CHECK_THROW(s.markTagged(0, 0, 0), std::runtime_error);
    CHECK_THROW(s.markCoupledAsTagged(), std::runtime_error);
    // A permutation of the right length (component 0 holds 10 dofs, all of
    // them free before finalize()), so that nothing but the missing
    // finalize() can reject it.
    gsVector<index_t> id10(10);
    for (index_t i = 0; i != id10.size(); ++i) id10[i] = i;
    CHECK_THROW(s.permuteFreeDofs(id10, 0), std::runtime_error);
    // These two also reach size(), which throws on its own.
    CHECK_THROW(s.anyPreImages(0), std::runtime_error);
    CHECK_THROW(s.inverseAsVector(0), std::runtime_error);
    // These would otherwise answer from the setup encoding without a throw.
    CHECK_THROW(s.findBoundary(0, 0), std::runtime_error);
    CHECK_THROW(s.findFree(0, 0), std::runtime_error);
    CHECK_THROW(s.findCoupled(0, -1, 0), std::runtime_error);
    CHECK_THROW(s.findFreeUncoupled(0, 0), std::runtime_error);
    CHECK_THROW(s.findTagged(0, 0), std::runtime_error);
    CHECK_THROW(s.inverseOnPatch(0), std::runtime_error);
    CHECK_THROW(s.boundarySizeWithDuplicates(), std::runtime_error);
}

// permuteFreeDofs() takes a permutation of [0, freeSize(c)) and nothing
// else, and rejects anything else before rewriting a single dof.  A value
// out of range would move a dof out of its component's free block, or out
// of the mapper; a repeated value would map two dofs onto one index.
TEST(permute_free_dofs_rejects_non_permutations)
{
    gsDofMapper m = threeCompUniform();       // freeSize(0) == 7
    m.markCoupledAsTagged();
    const std::string before = finalDigest(m);

    gsVector<index_t> p(7);
    p << 3, 0, 6, 1, 5, 2, 7;                 // 7 is out of range
    CHECK_THROW(m.permuteFreeDofs(p, 0), std::runtime_error);
    p << 3, 0, 6, 1, 5, 2, -1;
    CHECK_THROW(m.permuteFreeDofs(p, 0), std::runtime_error);
    p << 3, 0, 6, 1, 5, 2, 3;                 // 3 twice, 4 missing
    CHECK_THROW(m.permuteFreeDofs(p, 0), std::runtime_error);
    CHECK_THROW(m.permuteFreeDofs(gsVector<index_t>(6), 0), std::runtime_error);
    CHECK_THROW(m.permuteFreeDofs(gsVector<index_t>(8), 0), std::runtime_error);
    p << 3, 0, 6, 1, 5, 2, 4;
    CHECK_THROW(m.permuteFreeDofs(p, 3), std::runtime_error);
    CHECK_THROW(m.permuteFreeDofs(p, -1), std::runtime_error);
    CHECK_EQUAL(before, finalDigest(m));

    m.permuteFreeDofs(p, 0);
    CHECK_EQUAL(finalDigest(permutedMapper()), finalDigest(m));
}

// Every index a mapper hands out -- and one past the last, which
// lastIndex() computes -- must be representable,
// so a shift that would push them past index_t's maximum is rejected when it
// is set rather than overflowing later in index().
TEST(shift_must_keep_every_index_representable)
{
    const index_t lo = std::numeric_limits<index_t>::min();
    const index_t hi = std::numeric_limits<index_t>::max();
    gsDofMapper m = threeCompUniform();
    const index_t n = m.size(), nb = m.boundarySize();

    m.setShift(hi - n);                       // the largest valid shift
    CHECK_EQUAL(hi, m.firstFreeIndex(0) + m.size());
    CHECK_THROW(m.setShift(hi - n + 1), std::runtime_error);
    CHECK_THROW(m.setShift(hi), std::runtime_error);
    CHECK_THROW(m.addShift(1), std::runtime_error);
    CHECK_EQUAL(hi - n, m.firstFreeIndex(0));     // unchanged
    m.addShift(-5);
    CHECK_EQUAL(hi - n - 5, m.firstFreeIndex(0));

    m.setShift(lo);                           // any negative shift is fine
    CHECK_THROW(m.addShift(-1), std::runtime_error);
    CHECK_THROW(m.addShift(lo), std::runtime_error);
    CHECK_EQUAL(lo, m.firstFreeIndex(0));
    m.addShift(1);
    CHECK_EQUAL(lo + 1, m.firstFreeIndex(0));

    m.setBoundaryShift(hi - nb);
    CHECK_EQUAL(hi - 1, m.bindex(3, 1, 2));   // the last boundary index
    CHECK_THROW(m.setBoundaryShift(hi - nb + 1), std::runtime_error);

    // Before finalize() the final size is not known yet; mapSize() bounds
    // it, and a shift accepted against that bound stays valid.
    gsDofMapper s = raggedSetup();
    const index_t ms = static_cast<index_t>(s.mapSize());
    CHECK_THROW(s.setShift(hi - ms + 1), std::runtime_error);
    s.setShift(hi - ms);
    s.eliminateDof(0, 0, 0);
    s.finalize();
    CHECK_EQUAL(hi - 1, s.index(0, 0, 0));   // the eliminated dof, numbered last
}

// The per-dof accessors stay unchecked in Release builds (they run once per
// local dof per element in every assembly), but debug builds apply the same
// patch-specific bound as the checked entry points.
TEST(hot_path_accessors_check_local_bounds_in_debug)
{
    const gsDofMapper m = raggedPatchMapper();
    CHECK_EQUAL(15, m.index(3, 1, 1));        // the last valid dof
    CHECK_THROW_IN_DEBUG(m.index(2, 0, 1), std::logic_error);   // would be (patch 1, local 0)
    CHECK_THROW_IN_DEBUG(m.index(5, 0, 0), std::logic_error);
    CHECK_THROW_IN_DEBUG(m.index(4, 1, 1), std::logic_error);
    CHECK_THROW_IN_DEBUG(m.index(-1, 0, 0), std::logic_error);
    CHECK_THROW_IN_DEBUG(m.index(0, 2, 0), std::logic_error);
    CHECK_THROW_IN_DEBUG(m.index(0, 0, 2), std::logic_error);
    CHECK_THROW_IN_DEBUG(m.freeIndex(2, 0, 1), std::logic_error);
    CHECK_THROW_IN_DEBUG(m.bindex(2, 0, 1), std::logic_error);
    CHECK_THROW_IN_DEBUG(m.cindex(2, 0, 1), std::logic_error);
    CHECK_THROW_IN_DEBUG(m.tindex(2, 0, 1), std::logic_error);
    CHECK_THROW_IN_DEBUG(m.is_free(2, 0, 1), std::logic_error);
    gsMatrix<index_t> loc(2, 1), glob;
    loc << 1, 2;
    CHECK_THROW_IN_DEBUG(m.localToGlobal(loc, 0, glob, 1), std::logic_error);
}

// localToGlobal(2) reads locals(i,0) for i < rows(), so they need exactly one
// column unless there are no rows at all.  Checked in debug builds only (the
// assembly hot path); an empty input is valid in every build.
TEST(local_to_global_requires_a_column)
{
    const gsDofMapper m = raggedPatchMapper();
    gsMatrix<index_t> glob;
    index_t nf = -1;
    CHECK_THROW_IN_DEBUG(m.localToGlobal(gsMatrix<index_t>(2, 0), 0, glob, 0), std::logic_error);
    CHECK_THROW_IN_DEBUG(m.localToGlobal2(gsMatrix<index_t>(2, 0), 0, glob, nf, 0), std::logic_error);
    {
        gsMatrix<index_t> wide(2, 2);
        wide.setZero();
        CHECK_THROW_IN_DEBUG(m.localToGlobal(wide, 0, glob, 0), std::logic_error);
        CHECK_THROW_IN_DEBUG(m.localToGlobal2(wide, 0, glob, nf, 0), std::logic_error);
    }
    m.localToGlobal(gsMatrix<index_t>(0, 1), 0, glob, 0);
    CHECK_EQUAL(0, glob.rows());
    m.localToGlobal2(gsMatrix<index_t>(0, 0), 0, glob, nf, 0);
    CHECK_EQUAL(0, glob.rows());
    CHECK_EQUAL(0, nf);
}

// localToGlobal2 resizes globals to two columns before it reads locals, so
// passing one matrix as both would destroy the input it is about to read.
// Rejected before anything is written; checked in debug builds only (the
// assembly hot path).
TEST(local_to_global2_rejects_aliased_arguments)
{
#ifndef NDEBUG
    const gsDofMapper m = raggedPatchMapper();
    gsMatrix<index_t> values(3, 1);
    values << 0, 1, 2;
    index_t nf = -1;
    CHECK_THROW(m.localToGlobal2(values, 0, values, nf, 0), std::logic_error);
    CHECK_EQUAL(3, values.rows());
    CHECK_EQUAL(1, values.cols());
    CHECK_EQUAL(0, values(0,0));
    CHECK_EQUAL(1, values(1,0));
    CHECK_EQUAL(2, values(2,0));
    CHECK_EQUAL(-1, nf);
#endif
}

namespace
{
// n dofs on one patch, each in a coupling group of its own: the largest
// number of coupling ids a component of n dofs can use.
void checkEveryDofCoupled(const index_t n)
{
    gsVector<index_t> sz(1);
    sz[0] = n;
    gsDofMapper m(sz, 1);
    for (index_t i = 0; i != n; ++i)
        m.markCoupled(i, 0);
    m.finalize();
    CHECK_EQUAL(n, m.size());
    CHECK_EQUAL(n, m.freeSize());
    CHECK_EQUAL(n, m.coupledSize());
    CHECK_EQUAL(0, m.boundarySize());
    for (index_t i = 0; i != n; ++i)
    {
        CHECK_EQUAL(i, m.index(i, 0));
        CHECK(m.is_coupled_index(i));
    }
}

// n dofs on one patch, every one eliminated: the largest number of
// elimination ids a mapper of n dofs can use.
void checkEveryDofEliminated(const index_t n)
{
    gsVector<index_t> sz(1);
    sz[0] = n;
    gsDofMapper m(sz, 1);
    gsMatrix<index_t> all(n, 1);
    for (index_t i = 0; i != n; ++i)
        all(i, 0) = i;
    m.markBoundary(0, all);
    m.finalize();
    CHECK_EQUAL(n, m.size());
    CHECK_EQUAL(0, m.freeSize());
    CHECK_EQUAL(n, m.boundarySize());
    for (index_t i = 0; i != n; ++i)
        CHECK_EQUAL(i, m.bindex(i, 0));
}
} // anonymous namespace

// Construction accepts up to max(index_t) dofs, so the ids handed out during
// setup must be representable up to that count as well.  At the limit itself
// only a narrow index_t (int8_t, int16_t) is affordable; wider builds check
// the same property on a small mapper.
TEST(setup_ids_are_representable_up_to_the_dof_count_limit)
{
    checkEveryDofCoupled(5);
    checkEveryDofEliminated(5);

    const index_t imax = std::numeric_limits<index_t>::max();
    if (static_cast<size_t>(imax) <= (static_cast<size_t>(1) << 16))
    {
        checkEveryDofCoupled(imax);
        checkEveryDofEliminated(imax);
    }
}

// =========================================================================
// Ragged storage against per-component scalar mappers
// =========================================================================

namespace
{

// Two patches.  Component 0 has 5 and 5 local dofs, component 1 has 3 and 4,
// so the components' offset rows differ on every patch.  Each component is
// glued across the patches and has eliminated dofs; in component 0 a whole
// coupled group is eliminated, in component 1 a group spans three dofs.
void raggedComp0Calls(gsDofMapper & m, index_t c)
{
    m.matchDof(0, 4, 1, 0, c);
    m.matchDof(0, 3, 1, 1, c);
    m.eliminateDof(0, 0, c);
    m.eliminateDof(3, 0, c);
}

void raggedComp1Calls(gsDofMapper & m, index_t c)
{
    m.matchDof(0, 2, 1, 0, c);
    m.matchDof(0, 2, 1, 3, c);
    m.markBoundary(1, column(1, 2), c);
    m.eliminateDof(0, 0, c);
}

gsVector<index_t> twoSizes(index_t a, index_t b)
{
    gsVector<index_t> v(2);
    v << a, b;
    return v;
}

// The ragged mapper r and, for each of its components c, the scalar mapper
// s[c] of that component's patch sizes with the same setup calls.  After
// finalize() the ragged numbering is component-major: the free blocks of all
// components first, then their eliminated blocks, and within each block a
// component is numbered as its scalar mapper numbers it.  relabel() is that
// correspondence, computed from the scalar mappers' counts alone.
struct RaggedOracle
{
    gsDofMapper r;
    gsDofMapper s[2];
    index_t freeOff[2], elimOff[2], cpldOff[2], nFree;

    RaggedOracle()
    {
        std::vector<gsVector<index_t> > sz(2);
        sz[0] = twoSizes(5, 5);
        sz[1] = twoSizes(3, 4);
        r = gsDofMapper(sz, false);
        s[0] = gsDofMapper(sz[0], 1);
        s[1] = gsDofMapper(sz[1], 1);
        raggedComp0Calls(r, 0);
        raggedComp0Calls(s[0], 0);
        raggedComp1Calls(r, 1);
        raggedComp1Calls(s[1], 0);
        r.finalize();
        s[0].finalize();
        s[1].finalize();

        freeOff[0] = elimOff[0] = cpldOff[0] = 0;
        freeOff[1] = s[0].freeSize();
        elimOff[1] = s[0].boundarySize();
        cpldOff[1] = s[0].coupledSize();
        nFree = s[0].freeSize() + s[1].freeSize();
    }

    /// The global index in r of global index \a g of s[c].
    index_t relabel(index_t c, index_t g) const
    {
        return g < s[c].freeSize() ? freeOff[c] + g
                                   : nFree + elimOff[c] + (g - s[c].freeSize());
    }
};

void checkSameVector(const gsVector<index_t> & expected, const gsVector<index_t> & actual)
{
    CHECK_EQUAL(expected.size(), actual.size());
    if (expected.size() != actual.size()) return;
    for (index_t i = 0; i != expected.size(); ++i)
        CHECK_EQUAL(expected[i], actual[i]);
}

} // anonymous namespace

// A component of a ragged mapper is numbered, and answers every
// per-component query, as the scalar mapper of its own patch sizes with the
// same setup calls does, up to RaggedOracle::relabel().  The two components'
// offset rows differ on every patch, so a query that reads another
// component's offsets or sizes answers differently from the oracle.
TEST(ragged_component_is_a_relabelled_scalar_mapper)
{
    const RaggedOracle o;
    const gsDofMapper & r = o.r;

    CHECK_EQUAL(o.nFree, r.freeSize());
    CHECK_EQUAL(o.s[0].boundarySize() + o.s[1].boundarySize(), r.boundarySize());
    CHECK_EQUAL(o.s[0].coupledSize()  + o.s[1].coupledSize(),  r.coupledSize());
    CHECK_EQUAL(o.s[0].boundarySizeWithDuplicates() + o.s[1].boundarySizeWithDuplicates(),
                r.boundarySizeWithDuplicates());

    std::vector<std::pair<index_t,index_t> > sp, rp;
    for (index_t c = 0; c != 2; ++c)
    {
        const gsDofMapper & s = o.s[c];
        CHECK_EQUAL(s.size(),      r.size(c));
        CHECK_EQUAL(s.freeSize(),  r.freeSize(c));
        CHECK_EQUAL(s.totalSize(), r.totalSize(c));

        const gsVector<index_t> sv = s.asVector(0);
        gsVector<index_t> expected(sv.size());
        for (index_t j = 0; j != sv.size(); ++j)
            expected[j] = o.relabel(c, sv[j]);
        checkSameVector(expected, r.asVector(c));

        for (index_t k = 0; k != 2; ++k)
        {
            CHECK_EQUAL(s.offset(k),    r.offset(k, c));
            CHECK_EQUAL(s.patchSize(k), r.patchSize(k, c));

            const index_t n = static_cast<index_t>(s.patchSize(k));
            gsMatrix<index_t> locals(n, 1), sg, rg;
            for (index_t i = 0; i != n; ++i)
            {
                locals(i, 0) = i;
                CHECK_EQUAL(o.relabel(c, s.index(i, k)), r.index(i, k, c));
                CHECK_EQUAL(s.is_free(i, k),     r.is_free(i, k, c));
                CHECK_EQUAL(s.is_boundary(i, k), r.is_boundary(i, k, c));
                CHECK_EQUAL(s.is_coupled(i, k),  r.is_coupled(i, k, c));
                if (s.is_boundary(i, k))
                    CHECK_EQUAL(o.elimOff[c] + s.bindex(i, k), r.bindex(i, k, c));
                if (s.is_coupled(i, k))
                    CHECK_EQUAL(o.cpldOff[c] + s.cindex(i, k), r.cindex(i, k, c));
            }

            s.localToGlobal(locals, k, sg);
            r.localToGlobal(locals, k, rg, c);
            CHECK_EQUAL(n, rg.rows());
            for (index_t i = 0; i != n && i < rg.rows(); ++i)
                CHECK_EQUAL(o.relabel(c, sg(i, 0)), rg(i, 0));

            index_t snf = -1, rnf = -2;
            s.localToGlobal2(locals, k, sg, snf);
            r.localToGlobal2(locals, k, rg, rnf, c);
            CHECK_EQUAL(snf, rnf);
            CHECK_EQUAL(n, rg.rows());
            for (index_t i = 0; i != n && i < rg.rows(); ++i)
            {
                CHECK_EQUAL(sg(i, 0), rg(i, 0));
                CHECK_EQUAL(o.relabel(c, sg(i, 1)), rg(i, 1));
            }

            checkSameVector(s.findBoundary(k),      r.findBoundary(k, c));
            checkSameVector(s.findFree(k),          r.findFree(k, c));
            checkSameVector(s.findFreeUncoupled(k), r.findFreeUncoupled(k, c));
            for (index_t j = -1; j != 2; ++j)
                checkSameVector(s.findCoupled(k, j), r.findCoupled(k, j, c));
        }

        // Queries taking a global index, at every global index of the
        // component.
        const std::vector<std::pair<index_t,index_t> > sa = s.anyPreImages(0);
        const std::vector<std::pair<index_t,index_t> > ra = r.anyPreImages(c);
        CHECK_EQUAL(static_cast<size_t>(r.size()), ra.size());
        for (index_t g = 0; g != s.size(); ++g)
        {
            const index_t gl = o.relabel(c, g);
            CHECK_EQUAL(c, r.componentOf(gl));
            s.preImage(g, sp);
            r.preImage(gl, rp);
            CHECK(sp == rp);
            CHECK(s.anyPreImage(g) == r.anyPreImage(gl));
            CHECK(sa[g] == ra[gl]);
            for (index_t k = 0; k != 2; ++k)
            {
                index_t sl = -1, rl = -2;
                CHECK_EQUAL(s.indexOnPatch(g, k, sl), r.indexOnPatch(gl, k, rl));
                if (s.indexOnPatch(g, k)) CHECK_EQUAL(sl, rl);
            }
        }
        // Every other component's index is left unclaimed.
        index_t claimed = 0;
        for (size_t g = 0; g != ra.size(); ++g)
            if (-1 != ra[g].first) ++claimed;
        CHECK_EQUAL(s.size(), claimed);
    }

    // inverseOnPatch covers every component, keyed by the ragged global index.
    for (index_t k = 0; k != 2; ++k)
    {
        const std::map<index_t,index_t> ri = r.inverseOnPatch(k);
        size_t n = 0;
        for (index_t c = 0; c != 2; ++c)
        {
            const std::map<index_t,index_t> si = o.s[c].inverseOnPatch(k);
            n += si.size();
            for (std::map<index_t,index_t>::const_iterator it = si.begin();
                 it != si.end(); ++it)
            {
                const std::map<index_t,index_t>::const_iterator jt =
                    ri.find(o.relabel(c, it->first));
                CHECK(jt != ri.end());
                if (jt != ri.end()) CHECK_EQUAL(it->second, jt->second);
            }
        }
        CHECK_EQUAL(n, ri.size());
    }

    // ... and a shift relabels the ragged mapper like any other.
    gsDofMapper shifted = r;
    shifted.setShift(100);
    checkShiftIsARelabelling(r, shifted, 100);
}

// inverseAsVector() needs a permutation, so it is checked on F9, which has
// neither coupling nor elimination: component 1's six dofs (patches of 2
// and 4) occupy global indices 10..15, after component 0's ten.
TEST(ragged_inverse_as_vector)
{
    const gsDofMapper m = raggedPatchMapper();
    const gsVector<index_t> v0 = m.inverseAsVector(0), v1 = m.inverseAsVector(1);
    CHECK_EQUAL(16, v0.size());
    CHECK_EQUAL(16, v1.size());
    for (index_t g = 0; g != 16; ++g)
    {
        CHECK_EQUAL(g < 10 ? g : -1,      v0[g]);
        CHECK_EQUAL(g < 10 ? -1 : g - 10, v1[g]);
    }
}

}
