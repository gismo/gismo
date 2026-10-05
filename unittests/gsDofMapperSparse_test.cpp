/** @file gsDofMapperSparse_test.cpp

    @brief Dense versus sparse storage of gsDofMapper: a sparse mapper built by
    the same calls as a dense one must answer every query identically.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):

**/

#include "gismo_unittest.h"
#include <gsAssembler/gsDofMapperCreator.h>
#include <gsMSplines/gsMappedBasis.h>

#include <algorithm>
#include <fstream>
#include <functional>
#include <random>
#include <set>

#if defined(__SANITIZE_ADDRESS__)
#  define SPARSE_TEST_ASAN 1
#elif defined(__has_feature)
#  if __has_feature(address_sanitizer)
#    define SPARSE_TEST_ASAN 1
#  endif
#endif

// RLIMIT_AS is the only portable way to make an O(N) allocation fail instead
// of succeeding lazily; it is unusable under AddressSanitizer, which reserves
// terabytes of shadow memory.
#if defined(__linux__) && !defined(SPARSE_TEST_ASAN)
#  define SPARSE_TEST_ADDRESS_GUARD 1
#  include <sys/mman.h>
#  include <sys/resource.h>
#  include <unistd.h>
#endif

using namespace gismo;

namespace {

typedef gsDofMapper::storage storage;

// Dump helpers copied from unittests/gsDofMapper_test.cpp, where they are file-local.
// They query size() and index() and therefore work on finalized mappers only.

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

// firstIndex(0) and lastIndex().  See first_index_is_the_smallest_index_of_
// the_component for firstIndex(c) with c>=1.
std::string dumpFirstLast(const gsDofMapper & m)
{
    std::ostringstream os;
    os << "first0=" << m.firstIndex(0) << " last=" << m.lastIndex();
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

bool sameVec(const gsVector<index_t> & a, const gsVector<index_t> & b)
{
    return a.size() == b.size() && (a.size() == 0 || (a.array() == b.array()).all());
}

// Compares every query of a finalized mapper \a a against the finalized
// mapper \a b that was built by the same calls, whatever their storage modes.
void compareQueries(const gsDofMapper & a, const gsDofMapper & b)
{
    CHECK(a.isFinalized() && b.isFinalized());

    CHECK_EQUAL(a.numComponents(), b.numComponents());
    CHECK_EQUAL(a.componentsSize(), b.componentsSize());
    CHECK_EQUAL(a.numPatches(), b.numPatches());
    CHECK_EQUAL(a.size(), b.size());
    CHECK_EQUAL(a.freeSize(), b.freeSize());
    CHECK_EQUAL(a.boundarySize(), b.boundarySize());
    CHECK_EQUAL(a.coupledSize(), b.coupledSize());
    CHECK_EQUAL(a.taggedSize(), b.taggedSize());
    CHECK_EQUAL(a.mapSize(), b.mapSize());
    CHECK_EQUAL(a.boundarySizeWithDuplicates(), b.boundarySizeWithDuplicates());
    CHECK_EQUAL(a.hasUniformComponents(), b.hasUniformComponents());
    CHECK_EQUAL(a.hasDistinctComponentSpaces(), b.hasDistinctComponentSpaces());
    CHECK_EQUAL(a.isPermutation(), b.isPermutation());
    CHECK_EQUAL(a.isFinalized(), b.isFinalized());
    CHECK_EQUAL(a.allFree(), b.allFree());
    CHECK_EQUAL(a.layout(), b.layout());
    CHECK(a.getTagged() == b.getTagged());
    CHECK_EQUAL(a.lastIndex(), b.lastIndex());

    const index_t nComp = a.numComponents();
    const index_t nPat  = static_cast<index_t>(a.numPatches());
    // lastIndex() = shift + freeSize()
    const index_t shift = a.lastIndex() - a.freeSize();

    for (index_t c = 0; c <= nComp; ++c)
        CHECK_EQUAL(a.firstIndex(c), b.firstIndex(c));

    for (index_t n = 0; n != static_cast<index_t>(a.mapSize()); ++n)
        CHECK_EQUAL(a.mapIndex(n), b.mapIndex(n));

    bool hasRemote = false, hasRemoteB = false;
    for (index_t c = 0; c != nComp; ++c)
        for (index_t k = 0; k != nPat; ++k)
            for (index_t i = 0; i != static_cast<index_t>(a.patchSize(k,c)); ++i)
            {
                hasRemote  = hasRemote  || a.is_remote(i,k,c);
                hasRemoteB = hasRemoteB || b.is_remote(i,k,c);
            }
    CHECK_EQUAL(hasRemote, hasRemoteB);

    for (index_t c = 0; c != nComp; ++c)
    {
        CHECK_EQUAL(a.totalSize(c), b.totalSize(c));
        CHECK_EQUAL(a.size(c), b.size(c));
        CHECK_EQUAL(a.freeSize(c), b.freeSize(c));
        CHECK(sameVec(a.asVector(c), b.asVector(c)));
        CHECK(a.anyPreImages(c) == b.anyPreImages(c));
        if (a.isPermutation())
            CHECK(sameVec(a.inverseAsVector(c), b.inverseAsVector(c)));

        for (index_t k = 0; k != nPat; ++k)
        {
            CHECK_EQUAL(a.patchSize(k,c), b.patchSize(k,c));
            CHECK_EQUAL(a.offset(k,c), b.offset(k,c));
            for (index_t i = 0; i != static_cast<index_t>(a.patchSize(k,c)); ++i)
            {
                CHECK_EQUAL(a.index(i,k,c), b.index(i,k,c));
                CHECK_EQUAL(a.freeIndex(i,k,c), b.freeIndex(i,k,c));
                CHECK_EQUAL(a.bindex(i,k,c), b.bindex(i,k,c));
                CHECK_EQUAL(a.cindex(i,k,c), b.cindex(i,k,c));
                CHECK_EQUAL(a.tindex(i,k,c), b.tindex(i,k,c));
                CHECK_EQUAL(a.is_free(i,k,c), b.is_free(i,k,c));
                CHECK_EQUAL(a.is_boundary(i,k,c), b.is_boundary(i,k,c));
                CHECK_EQUAL(a.is_coupled(i,k,c), b.is_coupled(i,k,c));
                CHECK_EQUAL(a.is_tagged(i,k,c), b.is_tagged(i,k,c));
            }

            CHECK(sameVec(a.findBoundary(k,c), b.findBoundary(k,c)));
            CHECK(sameVec(a.findFree(k,c), b.findFree(k,c)));
            CHECK(sameVec(a.findFreeUncoupled(k,c), b.findFreeUncoupled(k,c)));
            CHECK(sameVec(a.findTagged(k,c), b.findTagged(k,c)));
            for (index_t j = -1; j != nPat; ++j)
                CHECK(sameVec(a.findCoupled(k,j,c), b.findCoupled(k,j,c)));
        }
    }

    for (index_t k = 0; k != nPat; ++k)
        CHECK(a.inverseOnPatch(k) == b.inverseOnPatch(k));

    // one step beyond both ends of the shifted range
    for (index_t gl = shift - 2; gl != shift + a.size() + 2; ++gl)
    {
        const bool valid = gl >= shift && gl < shift + a.size();
        CHECK_EQUAL(a.is_free_index(gl), b.is_free_index(gl));
        CHECK_EQUAL(a.is_boundary_index(gl), b.is_boundary_index(gl));
        CHECK_EQUAL(a.is_coupled_index(gl), b.is_coupled_index(gl));
        CHECK_EQUAL(a.is_tagged_index(gl), b.is_tagged_index(gl));
        CHECK_EQUAL(a.is_remote_index(gl), b.is_remote_index(gl));
        if (valid)
        {
            std::vector<std::pair<index_t,index_t> > pa, pb;
            a.preImage(gl, pa);
            b.preImage(gl, pb);
            CHECK(pa == pb);
            CHECK(a.anyPreImage(gl) == b.anyPreImage(gl));
            CHECK_EQUAL(a.componentOf(gl), b.componentOf(gl));
        }
        for (index_t k = 0; k != nPat; ++k)
        {
            index_t la = -7, lb = -7;
            CHECK_EQUAL(a.indexOnPatch(gl, k, la), b.indexOnPatch(gl, k, lb));
            CHECK_EQUAL(la, lb);
            CHECK_EQUAL(a.indexOnPatch(gl, k), b.indexOnPatch(gl, k));
        }
    }

    PIN("counts", dumpCounts(a), dumpCounts(b));
    PIN("elimdup", dumpElimDup(a), dumpElimDup(b));
    PIN("perComponent", dumpPerComponent(a), dumpPerComponent(b));
    PIN("layout", dumpLayout(a), dumpLayout(b));
    PIN("firstLast", dumpFirstLast(a), dumpFirstLast(b));
    PIN("asVector", dumpAsVector(a), dumpAsVector(b));
    PIN("mapIndex", dumpMapIndex(a), dumpMapIndex(b));
    PIN("index", dumpIndex(a), dumpIndex(b));
    PIN("freeIndex", dumpFreeIndex(a), dumpFreeIndex(b));
    PIN("freeBoundaryFlags", dumpFreeBoundaryFlags(a), dumpFreeBoundaryFlags(b));
    PIN("bindex", dumpBindex(a), dumpBindex(b));
    PIN("globalFlags", dumpGlobalFlags(a, shift), dumpGlobalFlags(b, shift));
    PIN("indexOnPatch", dumpIndexOnPatch(a, shift), dumpIndexOnPatch(b, shift));
    PIN("componentOf", dumpComponentOf(a, shift), dumpComponentOf(b, shift));
    PIN("preImages", dumpPreImages(a, shift), dumpPreImages(b, shift));
    PIN("anyPreImage", dumpAnyPreImage(a, shift), dumpAnyPreImage(b, shift));
    PIN("tagged", dumpTagged(a, shift), dumpTagged(b, shift));
    for (index_t c = 0; c != nComp; ++c)
    {
        PIN("localToGlobal", dumpLocalToGlobal(a, c), dumpLocalToGlobal(b, c));
        PIN("localToGlobal2", dumpLocalToGlobal2(a, c), dumpLocalToGlobal2(b, c));
        PIN("findBoundaryFree", dumpFindBoundaryFree(a, c), dumpFindBoundaryFree(b, c));
        PIN("coupledQueries", dumpCoupledQueries(a, c, shift), dumpCoupledQueries(b, c, shift));
    }
    for (index_t k = 0; k != nPat; ++k)
        PIN("inverseOnPatch", dumpInverseOnPatch(a, k), dumpInverseOnPatch(b, k));

    std::ostringstream oa, ob;
    a.print(oa);
    b.print(ob);
    CHECK_EQUAL(oa.str(), ob.str());
}

// As compareQueries(), for a finalized dense mapper \a d against the sparse
// mapper \a s built by the same calls, including the storage modes.
void compareFinalized(const gsDofMapper & d, const gsDofMapper & s)
{
    CHECK(storage::dense  == d.storageMode());
    CHECK(storage::sparse == s.storageMode());
    compareQueries(d, s);
}

// Compares the queries that are available before finalize(), which must not
// be asked for size() or index().
void compareSetup(const gsDofMapper & a, const gsDofMapper & b)
{
    CHECK(!a.isFinalized() && !b.isFinalized());
    CHECK_EQUAL(a.numComponents(), b.numComponents());
    CHECK_EQUAL(a.numPatches(), b.numPatches());
    CHECK_EQUAL(a.mapSize(), b.mapSize());
    CHECK_EQUAL(a.layout(), b.layout());
    CHECK_EQUAL(a.hasUniformComponents(), b.hasUniformComponents());
    CHECK_EQUAL(a.hasDistinctComponentSpaces(), b.hasDistinctComponentSpaces());

    for (index_t c = 0; c != a.numComponents(); ++c)
    {
        CHECK_EQUAL(a.totalSize(c), b.totalSize(c));
        CHECK(sameVec(a.asVector(c), b.asVector(c)));
        for (index_t k = 0; k != static_cast<index_t>(a.numPatches()); ++k)
        {
            CHECK_EQUAL(a.patchSize(k,c), b.patchSize(k,c));
            CHECK_EQUAL(a.offset(k,c), b.offset(k,c));
        }
    }
    for (index_t n = 0; n != static_cast<index_t>(a.mapSize()); ++n)
        CHECK_EQUAL(a.mapIndex(n), b.mapIndex(n));
}

gsMultiBasis<real_t> twoPatchBasis()
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
    gsMultiBasis<real_t> mb(mp);
    mb.degreeElevate(1);
    mb.uniformRefine(2);
    return mb;
}

gsBoundaryConditions<real_t> allKindsBC()
{
    static gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west , condition_type::dirichlet, &g, 0, false, -1);
    bc.addCondition(1, boundary::north, condition_type::dirichlet, &g, 0, false,  0);
    bc.addCondition(1, boundary::east , condition_type::clamped  , &g, 0, false,  1);
    bc.addCondition(1, boundary::south, condition_type::collapsed, &g, 0, false, -1);
    bc.addCoupled(0, boundary::south, 0, boundary::north, 2, 0, -1);
    bc.addCornerValue(boundary::southeast, 0.0, 0, 0, -1);
    return bc;
}

gsDofMapper fixtureA(storage st)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    return createMapper(mb, allKindsBC(), 2, 0, /*conforming=*/true, /*finalize=*/true, st);
}

gsMatrix<index_t> col(std::initializer_list<index_t> l)
{
    gsMatrix<index_t> m(static_cast<index_t>(l.size()), 1);
    index_t i = 0;
    for (index_t v : l) m(i++, 0) = v;
    return m;
}

// Unfinalized mapper exercising every setup path, in the storage \a st.
gsDofMapper manualSetup(storage st, bool ragged)
{
    gsDofMapper m;
    if (ragged)
    {
        std::vector<gsVector<index_t> > sz(2);
        sz[0].resize(3); sz[0] << 20, 15, 12;
        sz[1].resize(3); sz[1] << 20, 16, 13;
        m = gsDofMapper(sz, true, st);
    }
    else
    {
        gsVector<index_t> sz(3);
        sz << 20, 15, 12;
        m = gsDofMapper(sz, 2, st);
    }

    m.matchDofs(0, col({15,16,17,18,19}), 1, col({0,1,2,3,4}), -1);
    m.matchDof(1, 14, 2, 0, 0);
    m.matchDof(1, 13, 2, 1, 1);
    m.markCoupled(3, 0, 1);
    m.markCoupled(3, 0, 0);

    gsMatrix<unsigned> b(4,1);
    b << 0, 1, 2, 3;
    m.colapseDofs(2, b, -1);
    gsMatrix<unsigned> b2(3,1);
    b2 << 5, 6, 7;
    m.colapseDofs(0, b2, 1);

    m.eliminateDof(10, 0, -1);
    m.eliminateDof(11, 1, 1);
    m.eliminateDof(19, 0, 0);          // a coupled dof
    m.markBoundary(1, col({2,3,4}), -1);
    m.markBoundary(2, col({8,9}), 0);
    m.matchDof(0, 10, 2, 11, 0);       // eliminated with free
    m.matchDof(0, 10, 1, 11, 1);       // eliminated with eliminated
    m.matchDof(0, 7, 1, 7, 0);         // free with free
    m.matchDof(0, 3, 0, 7, 0);         // coupled with coupled
    m.matchDof(2, 0, 2, 5, 0);         // collapsed with free
    return m;
}

// Four patches around one interior vertex, whose coupling group has four local dofs.
gsMultiBasis<real_t> fourPatchBasis()
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 1.0);
    gsMultiBasis<real_t> mb(mp);
    mb.degreeElevate(1);
    mb.uniformRefine(1);
    return mb;
}

// One scalar component; \a eliminate selects homogeneous Dirichlet conditions
// on every boundary side, otherwise the interfaces are glued and nothing is eliminated.
gsDofMapper fourPatchMapper(storage st, bool eliminate, bool finalize)
{
    gsMultiBasis<real_t> mb = fourPatchBasis();
    static gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    if (eliminate)
    {
        for (gsBoxTopology::const_biterator bit = mb.topology().bBegin();
             bit != mb.topology().bEnd(); ++bit)
            bc.addCondition(*bit, condition_type::dirichlet, &g, 0, false, -1);
        return createMapper(mb, bc, 1, 0, /*conforming=*/true, finalize, st);
    }
    return createMapper(mb, bc, dirichlet::none, iFace::glue, 1, 0, finalize, st);
}

gsBoundaryConditions<real_t> allKindsBC3()
{
    static gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    bc.addCondition(0, boundary::west , condition_type::dirichlet, &g, 0, false,  2);
    bc.addCondition(1, boundary::north, condition_type::dirichlet, &g, 0, false,  0);
    bc.addCondition(1, boundary::east , condition_type::clamped  , &g, 0, false,  1);
    bc.addCondition(1, boundary::south, condition_type::collapsed, &g, 0, false, -1);
    bc.addCoupled(0, boundary::south, 0, boundary::north, 2, 0, 0);
    bc.addCornerValue(boundary::southeast, 0.0, 0, 0, 2);
    return bc;
}

// Three components on two patches, every kind of condition, finalized.
gsDofMapper threeCompMapper(storage st)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    return createMapper(mb, allKindsBC3(), 3, 0, /*conforming=*/true, /*finalize=*/true, st);
}

// Marks three dofs of three different components as tagged.
void tagThree(gsDofMapper & m)
{
    m.markTagged(1, 0, 0);
    m.markTagged(3, 1, 1);
    m.markTagged(5, 0, 2);
}

// Two patches of four and five dofs, three components with different matchings
// and eliminations per component.
gsDofMapper manualThreeComp(storage st, bool finalize)
{
    gsVector<index_t> sz(2);
    sz << 4, 5;
    gsDofMapper m(sz, 3, st);
    m.matchDof(0, 3, 1, 0, 0);
    m.eliminateDof(0, 0, 0);
    m.matchDof(0, 2, 1, 1, 1);
    m.matchDof(0, 3, 1, 2, 1);
    m.eliminateDof(4, 1, 1);
    m.eliminateDof(0, 0, 1);
    m.eliminateDof(3, 1, 2);
    if (finalize) m.finalize();
    return m;
}

// One patch of three dofs, two components; every dof of component 0 is
// eliminated, so component 0 owns no free dof and the multi-component
// relabelling must still shift the later components.
gsDofMapper eliminatedOnlyComponent(storage st, bool finalize)
{
    gsVector<index_t> sz(1);
    sz[0] = 3;
    gsDofMapper m(sz, 2, st);
    for (index_t i = 0; i != 3; ++i)
        m.eliminateDof(i, 0, 0);
    m.matchDof(0, 0, 0, 2, 1);
    m.eliminateDof(1, 0, 1);
    if (finalize) m.finalize();
    return m;
}

// Ragged mapper whose component 1 has no dof at all.
gsDofMapper raggedEmptyComponent(storage st, bool finalize)
{
    std::vector<gsVector<index_t> > sz(3);
    sz[0].resize(2); sz[0] << 5, 4;
    sz[1].resize(2); sz[1] << 0, 0;
    sz[2].resize(2); sz[2] << 3, 6;
    gsDofMapper m(sz, /*hasDistinctComponentSpaces=*/false, st);
    m.matchDof(0, 4, 1, 0, 0);
    m.eliminateDof(0, 0, 0);
    m.matchDof(0, 2, 1, 5, 2);
    m.eliminateDof(2, 1, 2);
    m.markCoupled(1, 1, 2);
    m.eliminateDof(3, 1, 0);
    m.markCoupled(0, 0, 2);
    if (finalize) m.finalize();
    return m;
}

// Two components of different sizes on the same two patches.
gsDofMapper raggedCreateMapper(storage st)
{
    gsMultiBasis<real_t> mb = twoPatchBasis();
    gsMultiBasis<real_t> mb2 = twoPatchBasis();
    mb2.degreeElevate(1);
    std::vector<gsMultiBasis<real_t> > vmb;
    vmb.push_back(mb);
    vmb.push_back(mb2);
    return createMapper(vmb, allKindsBC(), 0, true, true, st);
}

// Identity layout with unequal totals per component. No dof is matched across
// patches: every patch aliases the same positions.
gsDofMapper raggedIdentity(storage st)
{
    gsDofMapper m;
    m.setIdentity(3, std::vector<size_t>{7, 10}, st);
    m.markBoundary(1, col({0, 4}), 0);
    m.eliminateDof(2, 0, 1);
    m.markCoupled(5, 0, 1);
    m.finalize();
    return m;
}

// Space size and storage mode of a mapper after it went through initSystem().
index_t numDofsAfterInitSystem(const gsMultiBasis<real_t> & mb, index_t nComp,
                               const gsDofMapper & m, storage & mode)
{
    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(mb);
    auto u = A.getSpace(mb, nComp);
    u.setupMapper(m);
    A.initSystem();
    mode = u.mapper().storageMode();
    return A.numDofs();
}


// --- sparse localize ------------------------------------------------------

typedef std::function<gsDofMapper(storage)> MapperBuilder;

struct NamedFixture
{
    const char *  name;
    MapperBuilder build;
    bool          hasCoupled;   // the coupled-only subset must be non-empty
};

// One scalar component, two patches of n0 and n1 dofs, n0 >= 8 and n1 >= 12.
// Marked positions: 0 and N-1 eliminated, 5 coupled alone, n0-1 and n0
// (patch 0 last, patch 1 first) one coupled pair.  N = n0+n1 is the number of
// positions; the ids below follow from the marks in position order.
gsDofMapper bigMarked(index_t n0, index_t n1, storage st, bool finalize = true)
{
    gsVector<index_t> sz(2);
    sz << n0, n1;
    gsDofMapper m(sz, 1, st);
    m.eliminateDof(0, 0, 0);
    m.matchDof(0, n0-1, 1, 0, 0);
    m.eliminateDof(n1-1, 1, 0);
    m.markCoupled(5, 0, 0);
    if (finalize) m.finalize();
    return m;
}

// index() of bigMarked() after finalize(), shift 0: regular ids in position
// order, then the coupled ids by first appearance, then the eliminated ones.
index_t expectedIndex(index_t i, index_t k, index_t n0, index_t n1)
{
    const index_t N = n0 + n1;
    if (k == 0)
    {
        if (i == 0)      return N - 3;
        if (i <= 4)      return i - 1;
        if (i == 5)      return N - 5;
        if (i <= n0 - 2) return i - 2;
        return N - 4;
    }
    if (i == 0)      return N - 4;
    if (i <= n1 - 2) return n0 + i - 4;
    return N - 2;
}

// After localize({1, 2, n0+6, N-5, N-4}).
index_t expectedLocal1(index_t i, index_t k, index_t n0, index_t n1)
{
    if (k == 0)
    {
        if (i == 0)      return 5;
        if (i == 2)      return 0;
        if (i == 3)      return 1;
        if (i == 5)      return 3;
        if (i == n0 - 1) return 4;
        return gsDofMapper::remoteDof();
    }
    if (i == 0)      return 4;
    if (i == 10)     return 2;
    if (i == n1 - 1) return 6;
    return gsDofMapper::remoteDof();
}

// After a further localize({1, 4}).
index_t expectedLocal2(index_t i, index_t k, index_t n0, index_t n1)
{
    if (k == 0)
    {
        if (i == 0)      return 2;
        if (i == 3)      return 0;
        if (i == n0 - 1) return 1;
        return gsDofMapper::remoteDof();
    }
    if (i == 0)      return 1;
    if (i == n1 - 1) return 3;
    return gsDofMapper::remoteDof();
}

index_t expectedAt(int stage, index_t i, index_t k, index_t n0, index_t n1)
{
    switch (stage)
    {
    case 0:  return expectedIndex(i, k, n0, n1);
    case 1:  return expectedLocal1(i, k, n0, n1);
    default: return expectedLocal2(i, k, n0, n1);
    }
}

std::vector<index_t> bigLocalSet(index_t n0, index_t n1)
{
    const index_t N = n0 + n1;
    return std::vector<index_t>{1, 2, n0 + 6, N - 5, N - 4};
}

// Checks index() and is_remote() of a bigMarked() mapper at the given
// (position, patch) pairs against the closed form of the given stage
// (0 finalized, 1 localized, 2 localized again).
void checkBigAt(const gsDofMapper & m, int stage, index_t n0, index_t n1,
                const std::vector<std::pair<index_t,index_t> > & pts)
{
    for (size_t q = 0; q != pts.size(); ++q)
    {
        const index_t i = pts[q].first, k = pts[q].second;
        const index_t e = expectedAt(stage, i, k, n0, n1);
        CHECK_EQUAL(e, m.index(i, k, 0));
        CHECK_EQUAL(e == gsDofMapper::remoteDof(), m.is_remote(i, k, 0));
    }
}

void checkBigCounts(const gsDofMapper & m, int stage, index_t n0, index_t n1)
{
    const index_t N = n0 + n1;
    switch (stage)
    {
    case 0:
        CHECK_EQUAL(N - 3, m.freeSize());
        CHECK_EQUAL(2, m.coupledSize());
        CHECK_EQUAL(2, m.boundarySize());
        CHECK_EQUAL(N - 1, m.size());
        break;
    case 1:
        CHECK_EQUAL(5, m.freeSize());
        CHECK_EQUAL(2, m.coupledSize());
        CHECK_EQUAL(7, m.size());
        CHECK_EQUAL(N - 6, m.boundarySizeWithDuplicates());
        break;
    default:
        CHECK_EQUAL(2, m.freeSize());
        CHECK_EQUAL(1, m.coupledSize());
        CHECK_EQUAL(4, m.size());
    }
}

std::vector<std::pair<index_t,index_t> > bigSamplePoints(index_t n0, index_t n1)
{
    std::vector<std::pair<index_t,index_t> > p;
    const index_t i0[] = {0, 1, 2, 3, 4, 5, 6, n0 - 2, n0 - 1};
    const index_t i1[] = {0, 1, 2, 10, n1 - 2, n1 - 1};
    for (size_t q = 0; q != sizeof(i0)/sizeof(i0[0]); ++q) p.push_back(std::make_pair(i0[q], index_t(0)));
    for (size_t q = 0; q != sizeof(i1)/sizeof(i1[0]); ++q) p.push_back(std::make_pair(i1[q], index_t(1)));
    return p;
}

std::vector<std::pair<index_t,index_t> > bigAllPoints(index_t n0, index_t n1)
{
    std::vector<std::pair<index_t,index_t> > p;
    for (index_t i = 0; i != n0; ++i) p.push_back(std::make_pair(i, index_t(0)));
    for (index_t i = 0; i != n1; ++i) p.push_back(std::make_pair(i, index_t(1)));
    return p;
}

struct BigBytes
{
    size_t setup, finalized, localized;
};

// Builds a bigMarked() mapper in storage \a st and takes it through
// finalize(), localize(bigLocalSet) and localize({1, 4}), checking the closed
// form at the sample positions after each step.  Only calls that are
// O(#marked + #runs) in sparse storage are made, so N may be 10^9.
BigBytes runBigMarked(index_t n0, index_t n1, storage st)
{
    const std::vector<std::pair<index_t,index_t> > pts = bigSamplePoints(n0, n1);
    BigBytes b;

    gsDofMapper m = bigMarked(n0, n1, st, false);
    CHECK(st == m.storageMode());
    b.setup = m.nBytes();

    m.finalize();
    b.finalized = m.nBytes();
    checkBigCounts(m, 0, n0, n1);
    checkBigAt(m, 0, n0, n1, pts);

    m.localize(bigLocalSet(n0, n1));
    b.localized = m.nBytes();
    checkBigCounts(m, 1, n0, n1);
    checkBigAt(m, 1, n0, n1, pts);

    m.localize(std::vector<index_t>{1, 4});
    checkBigCounts(m, 2, n0, n1);
    checkBigAt(m, 2, n0, n1, pts);
    return b;
}

#ifdef SPARSE_TEST_ADDRESS_GUARD
// Lowers the soft RLIMIT_AS to the current virtual size plus 64 MiB for the
// lifetime of the object.  Every legitimate allocation of the sparse path is a
// few KiB, whereas a full table of 10^9 positions needs at least 119 MiB (one
// bit each) and 3.7 GiB as index_t, so an O(N) allocation throws std::bad_alloc
// before a page of it is written.
class AddressSpaceGuard
{
public:
    AddressSpaceGuard() : m_restore(false)
    {
        unsigned long long pages = 0;
        std::ifstream statm("/proc/self/statm");
        statm >> pages;
        const rlim_t vsize = static_cast<rlim_t>(pages) * static_cast<rlim_t>(sysconf(_SC_PAGESIZE));
        if (0 != getrlimit(RLIMIT_AS, &m_old))
            return;
        rlim_t cap = vsize + (rlim_t(64) << 20);
        if (RLIM_INFINITY != m_old.rlim_max && cap > m_old.rlim_max)
            cap = m_old.rlim_max;
        if (RLIM_INFINITY == m_old.rlim_cur || cap < m_old.rlim_cur)
        {
            rlimit lim = m_old;
            lim.rlim_cur = cap;
            m_restore = (0 == setrlimit(RLIMIT_AS, &lim));
        }
    }

    ~AddressSpaceGuard()
    {
        if (m_restore)
            setrlimit(RLIMIT_AS, &m_old);
    }

    // True if a 128 MiB mapping is refused.  mmap, because a compiler may
    // elide an unused new or malloc and make the check vacuous.
    bool selfTest() const
    {
        void * p = mmap(nullptr, size_t(128) << 20, PROT_READ | PROT_WRITE,
                        MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (MAP_FAILED == p)
            return true;
        munmap(p, size_t(128) << 20);
        return false;
    }

private:
    rlimit m_old;
    bool   m_restore;
};
#endif

// --- subsets of free dofs --------------------------------------------------

// Sorted, unique, shifted ids of the free dofs of \a m, as index() returns
// them.  The free range starts at lastIndex() - freeSize(), which is not
// firstIndex() when component 0 owns no free dof.
//  0 empty, 1 all, 2 middle half, 3 every third, 4 coupled dofs only,
//  5 the last one, 6 pseudo-random half.
std::vector<index_t> pick(const gsDofMapper & m, int mode)
{
    const index_t n = m.freeSize();
    const index_t f = m.lastIndex() - n;
    std::vector<index_t> r;
    switch (mode)
    {
    case 0:
        break;
    case 1:
        for (index_t i = 0; i != n; ++i) r.push_back(f + i);
        break;
    case 2:
        for (index_t i = n / 4; i < n / 2; ++i) r.push_back(f + i);
        break;
    case 3:
        for (index_t i = 1; i < n; i += 3) r.push_back(f + i);
        break;
    case 4:
        for (index_t c = 0; c != m.numComponents(); ++c)
            for (index_t k = 0; k != static_cast<index_t>(m.numPatches()); ++k)
            {
                const gsVector<index_t> cp = m.findCoupled(k, -1, c);
                for (index_t q = 0; q != cp.size(); ++q)
                {
                    const index_t g = m.index(cp[q], k, c);
                    if (g >= f && g < f + n) r.push_back(g);
                }
            }
        std::sort(r.begin(), r.end());
        r.erase(std::unique(r.begin(), r.end()), r.end());
        break;
    case 5:
        if (n > 0) r.push_back(f + n - 1);
        break;
    default:
        {
            std::mt19937 rng(1234 + mode);
            for (index_t i = 0; i != n; ++i)
                if (rng() & 1u) r.push_back(f + i);
        }
    }
    return r;
}

std::vector<NamedFixture> localizeFixtures()
{
    std::vector<NamedFixture> f;
    NamedFixture x;
    x.name = "fixtureA"; x.hasCoupled = true;
    x.build = [](storage st) { return fixtureA(st); };
    f.push_back(x);
    x.name = "fourPatch_eliminated"; x.hasCoupled = false;
    x.build = [](storage st) { return fourPatchMapper(st, true, true); };
    f.push_back(x);
    x.name = "fourPatch_glue"; x.hasCoupled = true;
    x.build = [](storage st) { return fourPatchMapper(st, false, true); };
    f.push_back(x);
    x.name = "threeComp"; x.hasCoupled = true;
    x.build = [](storage st) { return threeCompMapper(st); };
    f.push_back(x);
    x.name = "threeComp_tagged"; x.hasCoupled = false;
    x.build = [](storage st) { gsDofMapper m = threeCompMapper(st); tagThree(m); return m; };
    f.push_back(x);
    x.name = "manualThreeComp"; x.hasCoupled = true;
    x.build = [](storage st) { return manualThreeComp(st, true); };
    f.push_back(x);
    x.name = "eliminatedOnlyComponent"; x.hasCoupled = false;
    x.build = [](storage st) { return eliminatedOnlyComponent(st, true); };
    f.push_back(x);
    x.name = "raggedEmptyComponent"; x.hasCoupled = false;
    x.build = [](storage st) { return raggedEmptyComponent(st, true); };
    f.push_back(x);
    x.name = "raggedCreateMapper"; x.hasCoupled = false;
    x.build = [](storage st) { return raggedCreateMapper(st); };
    f.push_back(x);
    x.name = "raggedIdentity"; x.hasCoupled = false;
    x.build = [](storage st) { return raggedIdentity(st); };
    f.push_back(x);
    x.name = "bigMarked_16_13"; x.hasCoupled = true;
    x.build = [](storage st) { return bigMarked(16, 13, st); };
    f.push_back(x);
    return f;
}

gsDofMapper buildShifted(const NamedFixture & fix, storage st, index_t shift)
{
    gsDofMapper m = fix.build(st);
    if (shift != 0) m.setShift(shift);
    return m;
}

// --- gsFeSpace::setMapperStorage -------------------------------------------

// Sink that records the pattern and the matrix passed by computePattern_into
// and assemble_into.  The assembler calls a sink from all OpenMP threads.
struct PatternSink
{
    std::vector<gsEigen::Triplet<real_t,index_t> > mat;
    std::set<std::pair<index_t,index_t> > pattern;

    void addMatrix(const gsVector<index_t> & rows, const gsVector<index_t> & cols,
                   const gsMatrix<real_t> & block)
    {
#       pragma omp critical (PatternSink)
        for (index_t j = 0; j != cols.size(); ++j)
            for (index_t i = 0; i != rows.size(); ++i)
                if (rows[i] >= 0 && cols[j] >= 0)
                    mat.push_back(gsEigen::Triplet<real_t,index_t>(rows[i], cols[j], block(i,j)));
    }

    void addRhs(const gsVector<index_t> &, const gsMatrix<real_t> &) {}

    void addPattern(const gsVector<index_t> & rows, const gsVector<index_t> & cols)
    {
#       pragma omp critical (PatternSink)
        for (index_t j = 0; j != cols.size(); ++j)
            for (index_t i = 0; i != rows.size(); ++i)
                if (rows[i] >= 0 && cols[j] >= 0)
                    pattern.insert(std::make_pair(rows[i], cols[j]));
    }

    // rows x cols matrix; duplicate triplets are summed
    gsSparseMatrix<real_t> matrix(index_t rows, index_t cols) const
    {
        gsSparseMatrix<real_t> K(rows, cols);
        K.setFromTriplets(mat.begin(), mat.end());
        return K;
    }
};

// Number of nonzeros of \a K whose (row, col) is not in \a pattern.
index_t missingFromPattern(const gsSparseMatrix<real_t> & K,
                           const std::set<std::pair<index_t,index_t> > & pattern)
{
    index_t missing = 0;
    for (index_t o = 0; o != K.outerSize(); ++o)
        for (gsSparseMatrix<real_t>::InnerIterator it(K, o); it; ++it)
            if (0 != it.value() && !pattern.count(std::make_pair(it.row(), it.col())))
                ++missing;
    return missing;
}

struct PoissonPass
{
    gsSparseMatrix<real_t> K;
    gsMatrix<real_t>       F;
    gsDofMapper            mapper;
    std::set<std::pair<index_t,index_t> > pattern;
};

// Assembles a Poisson problem on four patches with the mapper storage \a st
// requested through the space, checking the storage mode after every step.
PoissonPass poissonPass(storage st)
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 1.0);
    gsMultiBasis<real_t> mb(mp);
    mb.degreeElevate(1);
    mb.uniformRefine(1);
    gsFunctionExpr<real_t> f("1", 2), g("0", 2);
    gsBoundaryConditions<real_t> bc;
    for (gsBoxTopology::const_biterator bit = mb.topology().bBegin();
         bit != mb.topology().bEnd(); ++bit)
        bc.addCondition(*bit, condition_type::dirichlet, &g, 0, false, -1);
    bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(mb);
    gsExprAssembler<real_t>::geometryMap G = A.getMap(mp);
    auto u = A.getSpace(mb);
    CHECK(storage::dense == u.mapperStorage());
    u.setMapperStorage(st);
    CHECK(st == u.mapperStorage());
    u.setup(bc, dirichlet::interpolation, 0);
    CHECK(st == u.mapper().storageMode());
    A.initSystem();
    CHECK(st == u.mapper().storageMode());
    auto ff = A.getCoeff(f, G);
    A.assemble(igrad(u, G) * igrad(u, G).tr() * meas(G), u * ff * meas(G));
    CHECK(st == u.mapper().storageMode());

    PatternSink S;
    A.computePattern_into(S, igrad(u, G) * igrad(u, G).tr());

    PoissonPass r;
    r.pattern = S.pattern;
    r.K = A.matrix();
    r.F = A.rhs();
    r.mapper = u.mapper();
    return r;
}

// --- gsDofMapper::index_into -----------------------------------------------

// What the batch lookups were asked, so that a test can assert that its
// inputs reached the interesting branches.
struct BatchStats
{
    BatchStats() : queried(0), remote(0), nonAscending(0), empty(0), colPositive(0) { }
    size_t queried, remote, nonAscending, empty, colPositive;
    std::set<std::pair<index_t,index_t> > ends;  // (component*numPatches+patch, position)
};

// Entry i is index(act(i,col), k, c).
std::vector<index_t> referenceIndices(const gsDofMapper & m, const gsMatrix<index_t> & act,
                                      index_t col, index_t k, index_t c)
{
    std::vector<index_t> ref(act.rows());
    for (index_t i = 0; i != act.rows(); ++i)
        ref[i] = m.index(act(i,col), k, c);
    return ref;
}

// Calls index_into on column \a col of \a act and checks every entry against
// the per-index reference; guard entries around the output detect writes
// outside [0, act.rows()).
std::vector<index_t> checkBatch(const gsDofMapper & m, const gsMatrix<index_t> & act,
                                index_t col, index_t k, index_t c, BatchStats * st = nullptr)
{
    const index_t n = act.rows();
    const index_t guard = -987654;
    std::vector<index_t> buf(n + 2, guard);
    m.index_into(act, col, k, c, buf.data() + 1);
    CHECK_EQUAL(guard, buf[0]);
    CHECK_EQUAL(guard, buf[n + 1]);

    const std::vector<index_t> ref = referenceIndices(m, act, col, k, c);
    const std::vector<index_t> out(buf.begin() + 1, buf.end() - 1);
    CHECK_EQUAL(ref.size(), out.size());
    index_t bad = -1;
    for (index_t i = 0; bad < 0 && i != n; ++i)
        if (ref[i] != out[i]) bad = i;
    CHECK_EQUAL(index_t(-1), bad);
    if (bad >= 0)
        CHECK_EQUAL(ref[bad], out[bad]);

    if (st)
    {
        const index_t key = c * static_cast<index_t>(m.numPatches()) + k;
        if (0 == n) ++st->empty;
        if (col > 0) ++st->colPositive;
        bool desc = false;
        for (index_t i = 0; i != n; ++i)
        {
            ++st->queried;
            if (m.is_remote(act(i,col), k, c)) ++st->remote;
            if (act(i,col) == 0 || act(i,col) == m.patchSize(k, c) - 1)
                st->ends.insert(std::make_pair(key, act(i,col)));
            if (i + 1 < n && act(i+1,col) < act(i,col)) desc = true;
        }
        if (desc) ++st->nonAscending;
    }
    return out;
}

// Dense and sparse mapper give the reference and each other.
void checkBatchPair(const gsDofMapper & d, const gsDofMapper & s, const gsMatrix<index_t> & act,
                    index_t col, index_t k, index_t c, BatchStats & st)
{
    const std::vector<index_t> od = checkBatch(d, act, col, k, c, &st);
    const std::vector<index_t> os = checkBatch(s, act, col, k, c);
    CHECK(od == os);
}

// Single-column inputs for a patch of n >= 7 dofs: empty, single entries,
// whole patch ascending and descending, wrapping, permuted, a duplicate,
// strided.
std::vector<gsMatrix<index_t> > handMadeColumns(index_t n)
{
    std::vector<gsMatrix<index_t> > v;
    v.push_back(gsMatrix<index_t>(0, 1));
    v.push_back(col({0}));
    v.push_back(col({n-1}));
    gsMatrix<index_t> up(n, 1), down(n, 1);
    for (index_t i = 0; i != n; ++i) { up(i,0) = i; down(i,0) = n - 1 - i; }
    v.push_back(up);
    v.push_back(down);
    v.push_back(col({n-2, n-1, 0, 1}));
    v.push_back(col({5, 3, 4, 2}));
    v.push_back(col({2, 3, 3, 4}));
    v.push_back(col({0, 3, 6}));
    return v;
}

// Columns tied to the marked and run structure of bigMarked(16, 13).
std::vector<gsMatrix<index_t> > bigExtraColumns(index_t k, index_t n)
{
    std::vector<gsMatrix<index_t> > v;
    if (0 == k)
    {
        v.push_back(col({0,1,2,3,4,5,6}));
        v.push_back(col({n-2, n-1}));
    }
    else
        v.push_back(col({0,1,9,10,11,n-2,n-1}));
    return v;
}

// Feeds every input of the batch tests to the dense and sparse mapper: the
// hand-made columns, a packed matrix whose columns are used at col > 0, and,
// when \a bases is given (one multi-basis per component), the 7x7 grid of
// active_into columns of every patch.
void runBatchQueries(const gsDofMapper & d, const gsDofMapper & s,
                     const std::vector<gsMultiBasis<real_t> > & bases, bool big,
                     BatchStats & st)
{
    for (index_t c = 0; c != d.numComponents(); ++c)
        for (index_t k = 0; k != static_cast<index_t>(d.numPatches()); ++k)
        {
            const index_t n = d.patchSize(k, c);
            std::vector<gsMatrix<index_t> > cols = handMadeColumns(n);
            if (big)
            {
                const std::vector<gsMatrix<index_t> > e = bigExtraColumns(k, n);
                cols.insert(cols.end(), e.begin(), e.end());
            }
            for (size_t q = 0; q != cols.size(); ++q)
                checkBatchPair(d, s, cols[q], 0, k, c, st);

            gsMatrix<index_t> packed(4, 5);
            packed.col(0) << 0, 1, 2, 3;
            packed.col(1) << n-4, n-3, n-2, n-1;
            packed.col(2) << n-2, n-1, 0, 1;
            packed.col(3) << 5, 3, 4, 2;
            packed.col(4) << 2, 3, 3, 4;
            for (index_t j = 0; j != packed.cols(); ++j)
                checkBatchPair(d, s, packed, j, k, c, st);

            if (!bases.empty())
            {
                const gsBasis<real_t> & b = bases[c].basis(k);
                const gsMatrix<real_t> sup = b.support();
                gsMatrix<real_t> pts(2, 49);
                for (index_t a = 0; a != 7; ++a)
                    for (index_t e = 0; e != 7; ++e)
                    {
                        pts(0, 7*a+e) = sup(0,0) + (sup(0,1) - sup(0,0)) * e / 6;
                        pts(1, 7*a+e) = sup(1,0) + (sup(1,1) - sup(1,0)) * a / 6;
                    }
                gsMatrix<index_t> act;
                b.active_into(pts, act);
                CHECK_EQUAL(index_t(49), act.cols());
                for (index_t j = 0; j != act.cols(); ++j)
                    checkBatchPair(d, s, act, j, k, c, st);
            }
        }
}

// Every patch of every component was queried at position 0 and at its last
// position.
void checkEndsQueried(const gsDofMapper & m, const BatchStats & st)
{
    for (index_t c = 0; c != m.numComponents(); ++c)
        for (index_t k = 0; k != static_cast<index_t>(m.numPatches()); ++k)
        {
            const index_t key = c * static_cast<index_t>(m.numPatches()) + k;
            CHECK(st.ends.count(std::make_pair(key, index_t(0))) > 0);
            CHECK(st.ends.count(std::make_pair(key, m.patchSize(k, c) - 1)) > 0);
        }
}

struct BasisFixture
{
    const char * name;
    MapperBuilder build;
    std::vector<gsMultiBasis<real_t> > bases;   // one per component
};

std::vector<BasisFixture> batchFixtures()
{
    std::vector<BasisFixture> f;
    BasisFixture x;
    x.name = "fixtureA";
    x.build = [](storage st) { return fixtureA(st); };
    x.bases.assign(2, twoPatchBasis());
    f.push_back(x);
    x.name = "fourPatch_eliminated";
    x.build = [](storage st) { return fourPatchMapper(st, true, true); };
    x.bases.assign(1, fourPatchBasis());
    f.push_back(x);
    x.name = "fourPatch_glue";
    x.build = [](storage st) { return fourPatchMapper(st, false, true); };
    f.push_back(x);
    x.name = "threeComp";
    x.build = [](storage st) { return threeCompMapper(st); };
    x.bases.assign(3, twoPatchBasis());
    f.push_back(x);
    x.name = "raggedCreateMapper";
    x.build = [](storage st) { return raggedCreateMapper(st); };
    x.bases.assign(1, twoPatchBasis());
    x.bases.push_back(twoPatchBasis());
    x.bases[1].degreeElevate(1);
    f.push_back(x);
    return f;
}

// --- assembler with sinks ----------------------------------------------------

struct SinkResult
{
    gsSparseMatrix<real_t> K;
    std::set<std::pair<index_t,index_t> > pattern;
    index_t nTest, nTrial;
};

// Distinct test space on the same basis (hence the same actives) whose mapper
// differs from the trial mapper: Dirichlet on all sides for u, on the west
// sides only for v.
SinkResult distinctTestSpacePass(storage st)
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 1.0);
    gsMultiBasis<real_t> mb = fourPatchBasis();
    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bcAll, bcWest;
    for (gsBoxTopology::const_biterator bit = mb.topology().bBegin();
         bit != mb.topology().bEnd(); ++bit)
    {
        bcAll.addCondition(*bit, condition_type::dirichlet, &g, 0, false, -1);
        if (boundary::west == bit->side())
            bcWest.addCondition(*bit, condition_type::dirichlet, &g, 0, false, -1);
    }
    bcAll.setGeoMap(mp);
    bcWest.setGeoMap(mp);

    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(mb);
    gsExprAssembler<real_t>::geometryMap G = A.getMap(mp);
    auto u = A.getSpace(mb);
    auto v = A.getTestSpace(u, mb);
    u.setMapperStorage(st);
    v.setMapperStorage(st);
    u.setup(bcAll , dirichlet::interpolation, 0);
    v.setup(bcWest, dirichlet::interpolation, 0);
    A.initSystem();
    CHECK(st == u.mapper().storageMode());
    CHECK(st == v.mapper().storageMode());

    SinkResult r;
    r.nTest  = A.numTestDofs();
    r.nTrial = A.numDofs();
    CHECK(r.nTest != r.nTrial);

    PatternSink S;
    A.computePattern_into(S, igrad(v, G) * igrad(u, G).tr() * meas(G));
    r.pattern = S.pattern;
    PatternSink T;
    A.assemble_into(T, igrad(v, G) * igrad(u, G).tr() * meas(G));
    r.K = T.matrix(r.nTest, r.nTrial);
    return r;
}

// Vector-valued trial space (two components) with Dirichlet conditions on all sides.
SinkResult vectorValuedPass(storage st)
{
    gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 1, 1.0);
    gsMultiBasis<real_t> mb = twoPatchBasis();
    gsFunctionExpr<real_t> g("0", 2);
    gsBoundaryConditions<real_t> bc;
    for (gsBoxTopology::const_biterator bit = mb.topology().bBegin();
         bit != mb.topology().bEnd(); ++bit)
        bc.addCondition(*bit, condition_type::dirichlet, &g, 0, false, -1);
    bc.setGeoMap(mp);

    gsExprAssembler<real_t> A(1, 1);
    A.setIntegrationElements(mb);
    gsExprAssembler<real_t>::geometryMap G = A.getMap(mp);
    auto u = A.getSpace(mb, 2);
    u.setMapperStorage(st);
    u.setup(bc, dirichlet::interpolation, 0);
    A.initSystem();
    CHECK(st == u.mapper().storageMode());

    SinkResult r;
    r.nTest = r.nTrial = A.numDofs();
    PatternSink S;
    A.computePattern_into(S, u * u.tr() * meas(G));
    r.pattern = S.pattern;
    PatternSink T;
    A.assemble_into(T, u * u.tr() * meas(G));
    r.K = T.matrix(r.nTest, r.nTrial);
    return r;
}

// Dense and sparse sink results agree; every assembled entry is in the pattern.
void compareSinkResults(const SinkResult & d, const SinkResult & s)
{
    CHECK(!d.pattern.empty());
    CHECK(d.K.nonZeros() > 0);
    CHECK_EQUAL(d.nTest , s.nTest);
    CHECK_EQUAL(d.nTrial, s.nTrial);
    CHECK_EQUAL(d.nTest , d.K.rows());
    CHECK_EQUAL(d.nTrial, d.K.cols());
    CHECK(d.pattern == s.pattern);
    CHECK_EQUAL(0, missingFromPattern(d.K, d.pattern));
    CHECK_EQUAL(0, missingFromPattern(s.K, s.pattern));
    const real_t tol = 1e2 * std::numeric_limits<real_t>::epsilon();
    CHECK((d.K - s.K).norm() <= tol * d.K.norm());
}

} // namespace

SUITE(gsDofMapperSparse_test)
{
    // Fixture A: two conforming patches, two components, every kind of condition.
    TEST(fixtureA_dense_equals_sparse)
    {
        const gsDofMapper d = fixtureA(storage::dense);
        const gsDofMapper s = fixtureA(storage::sparse);
        CHECK(d.coupledSize() > 0);
        CHECK(d.boundarySize() > 0);
        CHECK(d.boundarySizeWithDuplicates() > d.boundarySize());
        compareFinalized(d, s);
        gsInfo << "fixture A nBytes: dense " << d.nBytes() << ", sparse " << s.nBytes() << "\n";
    }

    TEST(manual_setup_dense_equals_sparse)
    {
        for (int ragged = 0; ragged != 2; ++ragged)
        {
            gsDofMapper d = manualSetup(storage::dense , ragged != 0);
            gsDofMapper s = manualSetup(storage::sparse, ragged != 0);
            CHECK(storage::sparse == s.storageMode());

            // setup-time values, before finalize
            for (index_t c = 0; c != d.numComponents(); ++c)
                CHECK(sameVec(d.asVector(c), s.asVector(c)));
            CHECK_EQUAL(d.mapSize(), s.mapSize());

            d.finalize();
            s.finalize();
            CHECK(d.coupledSize() > 0);
            CHECK(d.boundarySize() > 0);

            d.markTagged(2, 0, 0);
            s.markTagged(2, 0, 0);
            d.markTagged(5, 1, 1);
            s.markTagged(5, 1, 1);
            d.markTagged(17, 0, 1);
            s.markTagged(17, 0, 1);
            CHECK(d.taggedSize() > 0);
            compareFinalized(d, s);

            d.markCoupledAsTagged();
            s.markCoupledAsTagged();
            d.setShift(7);
            s.setShift(7);
            d.setBoundaryShift(3);
            s.setBoundaryShift(3);
            compareFinalized(d, s);
        }
    }

    // Every createMapper overload forwards the storage mode.
    TEST(createMapper_overloads_forward_storage)
    {
        gsMultiBasis<real_t> mb = twoPatchBasis();
        const gsBoundaryConditions<real_t> bc = allKindsBC();
        const gsBoxTopology & topo = mb.topology();
        const storage sp = storage::sparse;

        // default argument is dense
        CHECK(storage::dense == createMapper(mb, bc, 2, 0, true, true).storageMode());
        CHECK(storage::dense == createMapper(mb, 1, true, true).storageMode());

        // 1: primary
        {
            gsDofMapper s = createMapper(mb, topo, bc, 2, 0, true, true, sp);
            CHECK(sp == s.storageMode());
            compareFinalized(createMapper(mb, topo, bc, 2, 0, true, true), s);
        }
        // 2: bases only
        {
            gsDofMapper s = createMapper(mb, 2, true, true, sp);
            CHECK(sp == s.storageMode());
            compareFinalized(createMapper(mb, 2, true, true), s);
        }
        // 3: bases + topology
        {
            gsDofMapper s = createMapper(mb, topo, 2, true, true, sp);
            CHECK(sp == s.storageMode());
            compareFinalized(createMapper(mb, topo, 2, true, true), s);
        }
        // 4: bases + bc
        {
            gsDofMapper s = createMapper(mb, bc, 2, 0, true, true, sp);
            CHECK(sp == s.storageMode());
            compareFinalized(createMapper(mb, bc, 2, 0, true, true), s);
        }
        // 5: strategies
        {
            gsDofMapper s = createMapper(mb, bc, dirichlet::elimination, iFace::glue, 2, 0, true, sp);
            CHECK(sp == s.storageMode());
            compareFinalized(createMapper(mb, bc, dirichlet::elimination, iFace::glue, 2, 0, true), s);

            gsDofMapper s2 = createMapper(mb, bc, dirichlet::none, iFace::glue, 2, 0, true, sp);
            CHECK(sp == s2.storageMode());
        }
        // 6: per-component function sets, shared pointer (delegates to the primary)
        {
            std::vector<const gsFunctionSet<real_t>*> shared(2, &mb);
            gsDofMapper s = createMapper(shared, topo, bc, 0, true, true, sp);
            CHECK(sp == s.storageMode());
            compareFinalized(createMapper(shared, topo, bc, 0, true, true), s);
        }
        // 6: distinct function sets (ragged path)
        {
            gsMultiBasis<real_t> mb2 = twoPatchBasis();
            std::vector<const gsFunctionSet<real_t>*> distinct;
            distinct.push_back(&mb);
            distinct.push_back(&mb2);
            gsDofMapper s = createMapper(distinct, topo, bc, 0, true, true, sp);
            CHECK(sp == s.storageMode());
            CHECK(s.hasDistinctComponentSpaces());
            compareFinalized(createMapper(distinct, topo, bc, 0, true, true), s);
        }
        // 7: vector of gsMultiBasis
        {
            std::vector<gsMultiBasis<real_t> > vmb(2, mb);
            gsDofMapper s = createMapper(vmb, bc, 0, true, true, sp);
            CHECK(sp == s.storageMode());
            compareFinalized(createMapper(vmb, bc, 0, true, true), s);
        }
        // mapped basis: the global-identity branch
        {
            auto geom = gsNurbsCreator<real_t>::BSplineSquare(2);
            gsMultiBasis<real_t> mbm(geom->basis());
            mbm.basis(0).uniformRefine(2);
            const index_t sz = mbm.basis(0).size();
            gsSparseMatrix<real_t> ident(sz, sz);
            ident.setIdentity();
            gsMappedBasis<2, real_t> mapB(mbm, ident);

            gsFunctionExpr<real_t> g("0", 2);
            gsBoundaryConditions<real_t> bcm;
            bcm.addCondition(0, boundary::west, condition_type::dirichlet, &g);

            gsDofMapper d = createMapper(mapB, bcm, 1, 0, false, true);
            gsDofMapper s = createMapper(mapB, bcm, 1, 0, false, true, sp);
            CHECK(storage::dense == d.storageMode());
            CHECK(sp == s.storageMode());
            CHECK(sameVec(d.asVector(0), s.asVector(0)));
            compareFinalized(d, s);
        }
    }

    TEST(setIdentity_sparse)
    {
        gsDofMapper d, s;
        d.setIdentity(1, std::vector<size_t>{4,6}, storage::dense);
        s.setIdentity(1, std::vector<size_t>{4,6}, storage::sparse);
        CHECK(storage::sparse == s.storageMode());
        CHECK(!d.hasUniformComponents());
        CHECK(!s.hasUniformComponents());
        CHECK_EQUAL(4u, (unsigned)s.totalSize(0));
        CHECK_EQUAL(6u, (unsigned)s.totalSize(1));
        d.finalize();
        s.finalize();
        for (index_t c = 0; c != 2; ++c)
            for (index_t i = 0; i != static_cast<index_t>(d.totalSize(c)); ++i)
                CHECK_EQUAL(d.index(i,0,c), s.index(i,0,c));
        compareFinalized(d, s);

        gsDofMapper sc;
        sc.setIdentity(1, 5, 2, storage::sparse);
        CHECK(storage::sparse == sc.storageMode());
        CHECK(sc.hasUniformComponents());

        // identity with several patches and eliminated dofs (global-identity layout)
        gsDofMapper d3, s3;
        d3.setIdentity(3, 9, 2, storage::dense);
        s3.setIdentity(3, 9, 2, storage::sparse);
        d3.markBoundary(1, col({0,4}), -1);
        s3.markBoundary(1, col({0,4}), -1);
        d3.matchDof(0, 2, 2, 6, 0);
        s3.matchDof(0, 2, 2, 6, 0);
        d3.finalize();
        s3.finalize();
        compareFinalized(d3, s3);

        // a full reset with the default argument returns to dense storage
        s3.setIdentity(2, 3);
        CHECK(storage::dense == s3.storageMode());
        s3.finalize();
        CHECK_EQUAL(3, s3.size());
        sc.setIdentity(2, 3, 1, storage::sparse);
        CHECK(storage::sparse == sc.storageMode());
    }

    TEST(copy_swap_preserve_storage)
    {
        gsDofMapper s = fixtureA(storage::sparse);
        gsDofMapper copy(s);
        CHECK(storage::sparse == copy.storageMode());
        compareFinalized(fixtureA(storage::dense), copy);

        gsDofMapper d = fixtureA(storage::dense);
        d.swap(copy);
        CHECK(storage::sparse == d.storageMode());
        CHECK(storage::dense == copy.storageMode());
        compareFinalized(copy, d);
    }

    // A sparse mapper installed in an assembler space survives initSystem().
    TEST(survives_initSystem)
    {
        gsMultiBasis<real_t> mb = twoPatchBasis();
        const gsDofMapper dense  = fixtureA(storage::dense);
        const gsDofMapper sparse = fixtureA(storage::sparse);

        gsExprAssembler<real_t> Ad(1, 1);
        Ad.setIntegrationElements(mb);
        auto ud = Ad.getSpace(mb, 2);
        ud.setupMapper(dense);
        Ad.initSystem();

        gsExprAssembler<real_t> As(1, 1);
        As.setIntegrationElements(mb);
        auto us = As.getSpace(mb, 2);
        us.setupMapper(sparse);
        As.initSystem();

        CHECK(storage::dense  == ud.mapper().storageMode());
        CHECK(storage::sparse == us.mapper().storageMode());
        CHECK_EQUAL(dense.freeSize(), Ad.numDofs());
        CHECK_EQUAL(Ad.numDofs(), As.numDofs());
    }

    TEST(localize_and_permute)
    {
        gsDofMapper s = fixtureA(storage::sparse);
        {
            gsDofMapper sl = s;
            gsDofMapper dl = fixtureA(storage::dense);
            const std::vector<index_t> l2g{0,1,2};
            sl.localize(l2g);
            dl.localize(l2g);
            CHECK(storage::sparse == sl.storageMode());
            for (index_t c = 0; c != dl.numComponents(); ++c)
                for (index_t k = 0; k != static_cast<index_t>(dl.numPatches()); ++k)
                    for (index_t i = 0; i != static_cast<index_t>(dl.patchSize(k, c)); ++i)
                        CHECK_EQUAL(dl.index(i, k, c), sl.index(i, k, c));
        }
        CHECK(storage::sparse == s.storageMode());

        // an invalid permutation throws and leaves the storage mode alone
        gsVector<index_t> bad = gsVector<index_t>::Zero(s.freeSize(0));
        CHECK_THROW(s.permuteFreeDofs(bad, 0), std::runtime_error);
        CHECK(storage::sparse == s.storageMode());

        gsDofMapper d = fixtureA(storage::dense);
        const index_t nFree = d.freeSize(0);
        gsVector<index_t> perm(nFree);
        for (index_t i = 0; i != nFree; ++i)
            perm[i] = nFree - 1 - i;
        d.permuteFreeDofs(perm, 0);
        s.permuteFreeDofs(perm, 0);
        CHECK(storage::dense == s.storageMode());
        for (index_t c = 0; c != 2; ++c)
            CHECK(sameVec(d.asVector(c), s.asVector(c)));
        CHECK_EQUAL(d.coupledSize(), s.coupledSize());
        CHECK_EQUAL(d.nBytes() > 0, s.nBytes() > 0);
    }


    // Four patches around an interior vertex, one scalar component.
    TEST(ncomp1_four_patch_interior_vertex)
    {
        for (int eliminate = 0; eliminate != 2; ++eliminate)
        {
            gsDofMapper d = fourPatchMapper(storage::dense , eliminate != 0, false);
            gsDofMapper s = fourPatchMapper(storage::sparse, eliminate != 0, false);
            CHECK(storage::sparse == s.storageMode());
            compareSetup(d, s);

            d.finalize();
            s.finalize();
            CHECK(d.coupledSize() > 0);
            if (eliminate)
                CHECK(d.boundarySize() > 0);
            else
                CHECK_EQUAL(0, d.boundarySize());

            // the interior vertex is one global dof with a preimage on each of the four patches
            size_t widest = 0;
            for (index_t gl = d.firstIndex(0); gl != d.firstIndex(0) + d.size(); ++gl)
            {
                std::vector<std::pair<index_t,index_t> > pre;
                d.preImage(gl, pre);
                widest = std::max(widest, pre.size());
            }
            CHECK_EQUAL(4u, (unsigned)widest);

            compareFinalized(d, s);
        }
    }

    TEST(ncomp3_all_condition_kinds_shift_and_tags)
    {
        gsDofMapper d = threeCompMapper(storage::dense);
        gsDofMapper s = threeCompMapper(storage::sparse);
        CHECK(d.coupledSize() > 0);
        CHECK(d.boundarySize() > 0);
        compareFinalized(d, s);

        tagThree(d);
        tagThree(s);
        CHECK_EQUAL(3, d.taggedSize());
        compareFinalized(d, s);

        d.setShift(1000);
        s.setShift(1000);
        d.addShift(5);
        s.addShift(5);
        d.setBoundaryShift(11);
        s.setBoundaryShift(11);
        compareFinalized(d, s);
    }

    TEST(ncomp3_manual_and_eliminated_only_component)
    {
        {
            gsDofMapper d = manualThreeComp(storage::dense , false);
            gsDofMapper s = manualThreeComp(storage::sparse, false);
            compareSetup(d, s);
            d.finalize();
            s.finalize();
            CHECK(d.coupledSize() > 0);
            CHECK(d.boundarySize() > 0);
            compareFinalized(d, s);
        }
        {
            gsDofMapper d = eliminatedOnlyComponent(storage::dense , false);
            gsDofMapper s = eliminatedOnlyComponent(storage::sparse, false);
            compareSetup(d, s);
            d.finalize();
            s.finalize();
            CHECK_EQUAL(0, d.freeSize(0));
            CHECK_EQUAL(0, s.freeSize(0));
            compareFinalized(d, s);
        }
    }

    TEST(ragged_with_empty_component)
    {
        gsDofMapper d = raggedEmptyComponent(storage::dense , false);
        gsDofMapper s = raggedEmptyComponent(storage::sparse, false);
        compareSetup(d, s);

        // a broadcast must validate the empty component, and a throwing call changes nothing
        CHECK_THROW(d.markBoundary(0, col({0}), -1), std::runtime_error);
        CHECK_THROW(s.markBoundary(0, col({0}), -1), std::runtime_error);
        compareSetup(d, s);

        d.finalize();
        s.finalize();
        CHECK_EQUAL(0u, (unsigned)s.totalSize(1));
        CHECK_EQUAL(0, s.size(1));
        compareFinalized(d, s);
    }

    TEST(ragged_createMapper_distinct_bases)
    {
        const gsDofMapper d = raggedCreateMapper(storage::dense);
        const gsDofMapper s = raggedCreateMapper(storage::sparse);
        CHECK(d.hasDistinctComponentSpaces());
        CHECK(s.hasDistinctComponentSpaces());
        CHECK(d.totalSize(0) != d.totalSize(1));
        compareFinalized(d, s);
    }

    TEST(identity_layout_ragged_and_uniform)
    {
        {
            const gsDofMapper d = raggedIdentity(storage::dense);
            const gsDofMapper s = raggedIdentity(storage::sparse);
            CHECK(gsDofMapper::GlobalIdentity == d.layout());
            CHECK(gsDofMapper::GlobalIdentity == s.layout());
            CHECK(d.totalSize(0) != d.totalSize(1));
            CHECK(d.boundarySize() > 0);
            compareFinalized(d, s);
        }
        {
            gsDofMapper d, s;
            d.setIdentity(2, 6, 3, storage::dense);
            s.setIdentity(2, 6, 3, storage::sparse);
            d.markBoundary(0, col({1, 2}), -1);
            s.markBoundary(0, col({1, 2}), -1);
            d.finalize();
            s.finalize();
            CHECK(gsDofMapper::GlobalIdentity == d.layout());
            CHECK(gsDofMapper::GlobalIdentity == s.layout());
            compareFinalized(d, s);

            d.markTagged(0, 0, 0);
            s.markTagged(0, 0, 0);
            CHECK_EQUAL(1, d.taggedSize());
            compareFinalized(d, s);
        }
    }

    // permuteFreeDofs converts a sparse mapper to dense storage; the result
    // must be the permuted dense mapper.
    TEST(permute_free_dofs_densifies_and_matches)
    {
        {
            gsDofMapper d = manualThreeComp(storage::dense , true);
            gsDofMapper s = manualThreeComp(storage::sparse, true);
            d.markTagged(1, 0, 0);
            s.markTagged(1, 0, 0);
            d.markTagged(2, 1, 1);
            s.markTagged(2, 1, 1);

            const index_t nFree = d.freeSize(1);
            CHECK(nFree > 1);
            gsVector<index_t> perm(nFree);
            for (index_t i = 0; i != nFree; ++i)
                perm[i] = nFree - 1 - i;
            d.permuteFreeDofs(perm, 1);
            s.permuteFreeDofs(perm, 1);
            CHECK(storage::dense == s.storageMode());
            compareQueries(d, s);
        }
        {
            gsDofMapper d = threeCompMapper(storage::dense);
            gsDofMapper s = threeCompMapper(storage::sparse);
            d.setShift(40);
            s.setShift(40);

            const index_t nFree = d.freeSize(0);
            CHECK(nFree > 1);
            gsVector<index_t> perm(nFree);
            for (index_t i = 0; i != nFree; ++i)
                perm[i] = (i + 3) % nFree;
            d.permuteFreeDofs(perm, 0);
            s.permuteFreeDofs(perm, 0);
            CHECK(storage::dense == s.storageMode());
            compareQueries(d, s);
        }
    }

    TEST(copy_assign_move_swap)
    {
        const gsDofMapper dA = fixtureA(storage::dense);
        const gsDofMapper dN = fourPatchMapper(storage::dense, true, true);

        // copy assignment
        {
            const gsDofMapper s = threeCompMapper(storage::sparse);
            gsDofMapper a;
            a = s;
            CHECK(storage::sparse == a.storageMode());
            CHECK(storage::sparse == s.storageMode());
            compareFinalized(threeCompMapper(storage::dense), a);
            compareQueries(s, a);
        }
        // move construction
        {
            gsDofMapper tmp = threeCompMapper(storage::sparse);
            gsDofMapper m2 = give(tmp);
            CHECK(storage::sparse == m2.storageMode());
            compareFinalized(threeCompMapper(storage::dense), m2);
        }
        // swap of two sparse mappers with different numbers of components
        {
            gsDofMapper a = fixtureA(storage::sparse);
            gsDofMapper b = fourPatchMapper(storage::sparse, true, true);
            CHECK_EQUAL(2, a.numComponents());
            CHECK_EQUAL(1, b.numComponents());
            a.swap(b);
            CHECK(storage::sparse == a.storageMode());
            CHECK(storage::sparse == b.storageMode());
            CHECK_EQUAL(1, a.numComponents());
            CHECK_EQUAL(2, b.numComponents());
            compareFinalized(dN, a);
            compareFinalized(dA, b);
        }
        // swap of an unfinalized sparse mapper with a finalized dense one
        {
            gsDofMapper u = fourPatchMapper(storage::sparse, true, false);
            gsDofMapper f = fixtureA(storage::dense);
            u.swap(f);
            CHECK(storage::dense == u.storageMode());
            CHECK(storage::sparse == f.storageMode());
            CHECK(u.isFinalized());
            CHECK(!f.isFinalized());
            compareQueries(dA, u);

            gsDofMapper twin = fourPatchMapper(storage::dense, true, false);
            compareSetup(twin, f);
            twin.finalize();
            f.finalize();
            compareFinalized(twin, f);
        }
    }

    // Calls that are guarded by GISMO_ENSURE throw in every build type and in
    // both storage modes.
    TEST(error_parity)
    {
        for (int sp = 0; sp != 2; ++sp)
        {
            const storage st = sp ? storage::sparse : storage::dense;
            gsDofMapper m = threeCompMapper(st);

            CHECK_THROW(m.mapIndex(static_cast<index_t>(m.mapSize())), std::runtime_error);
            CHECK_THROW(m.componentOf(m.firstIndex(0) + m.size()), std::runtime_error);
            CHECK_THROW(m.findBoundary(static_cast<index_t>(m.numPatches()), 0), std::runtime_error);
            CHECK_THROW(m.totalSize(m.numComponents()), std::runtime_error);
            CHECK_THROW(m.finalize(), std::runtime_error);

            gsVector<index_t> notPerm = gsVector<index_t>::Zero(m.freeSize(0));
            CHECK_THROW(m.permuteFreeDofs(notPerm, 0), std::runtime_error);
            CHECK(st == m.storageMode());

            gsDofMapper u = manualThreeComp(st, false);
            CHECK_THROW(u.markBoundary(0, col({4}), 0), std::runtime_error);
            CHECK_THROW(u.markBoundary(1, col({5}), 2), std::runtime_error);
        }
    }

    // anyPreImages() skips remote positions of a localized mapper, and
    // inverseAsVector() throws unless the mapper is a permutation.
    TEST(localized_remote_anyPreImages_inverseAsVector)
    {
        for (int sp = 0; sp != 2; ++sp)
        {
            const storage st = sp ? storage::sparse : storage::dense;
            const index_t n0 = 16, n1 = 13;
            gsDofMapper m = bigMarked(n0, n1, st);
            m.localize(bigLocalSet(n0, n1));
            CHECK(st == m.storageMode());

            bool hasRemote = false;
            for (index_t k = 0; k != static_cast<index_t>(m.numPatches()); ++k)
                for (index_t i = 0; i != static_cast<index_t>(m.patchSize(k,0)); ++i)
                    hasRemote = hasRemote || m.is_remote(i,k,0);
            CHECK(hasRemote);
            CHECK_EQUAL(7, m.size());
            CHECK(!m.isPermutation());

            const std::vector<std::pair<index_t,index_t> > r = m.anyPreImages(0);
            CHECK_EQUAL(static_cast<size_t>(m.size()), r.size());
            for (index_t v = 0; v != static_cast<index_t>(r.size()); ++v)
            {
                const index_t k = r[v].first, i = r[v].second;
                const bool okPatch = 0 <= k && k < static_cast<index_t>(m.numPatches());
                CHECK(okPatch);
                if (!okPatch) continue;
                const bool okDof = 0 <= i && i < static_cast<index_t>(m.patchSize(k,0));
                CHECK(okDof);
                if (!okDof) continue;
                CHECK_EQUAL(v, m.index(i,k,0));
                CHECK(!m.is_remote(i,k,0));
            }
            CHECK(std::make_pair(index_t(0), index_t(2))  == r[0]);
            CHECK(std::make_pair(index_t(0), index_t(3))  == r[1]);
            CHECK(std::make_pair(index_t(1), index_t(10)) == r[2]);
            CHECK(std::make_pair(index_t(0), index_t(5))  == r[3]);
            CHECK(std::make_pair(index_t(0), index_t(0))  == r[5]);
            CHECK(std::make_pair(index_t(1), index_t(12)) == r[6]);
            // Two positions share id 4; either is a valid representative.
            CHECK(std::make_pair(index_t(0), index_t(15)) == r[4] ||
                  std::make_pair(index_t(1), index_t(0))  == r[4]);

            CHECK_THROW(m.inverseAsVector(0), std::runtime_error);
            CHECK_THROW(bigMarked(n0, n1, st).inverseAsVector(0), std::runtime_error);

            gsVector<index_t> sz(2);
            sz << n0, n1;
            gsDofMapper p(sz, 1, st);
            p.finalize();
            CHECK(p.isPermutation());
            const gsVector<index_t> inv = p.inverseAsVector(0);
            CHECK_EQUAL(static_cast<index_t>(p.size()), inv.size());
            for (index_t k = 0; k != static_cast<index_t>(p.numPatches()); ++k)
                for (index_t i = 0; i != static_cast<index_t>(p.patchSize(k,0)); ++i)
                {
                    const index_t g = p.index(i,k,0);
                    const bool inRange = 0 <= g && g < inv.size();
                    CHECK(inRange);
                    if (inRange)
                        CHECK_EQUAL(static_cast<index_t>(p.offset(k,0)) + i, inv[g]);
                }
            CHECK((inv.array() != -1).all());
        }
    }

    TEST(initSystem_survival_ncomp1_and_3)
    {
        {
            const gsMultiBasis<real_t> mb = fourPatchBasis();
            storage md = storage::sparse, ms = storage::dense;
            const index_t nd = numDofsAfterInitSystem(mb, 1, fourPatchMapper(storage::dense , true, true), md);
            const index_t ns = numDofsAfterInitSystem(mb, 1, fourPatchMapper(storage::sparse, true, true), ms);
            CHECK(storage::dense  == md);
            CHECK(storage::sparse == ms);
            CHECK(nd > 0);
            CHECK_EQUAL(nd, ns);
        }
        {
            const gsMultiBasis<real_t> mb = twoPatchBasis();
            storage md = storage::sparse, ms = storage::dense;
            const index_t nd = numDofsAfterInitSystem(mb, 3, threeCompMapper(storage::dense ), md);
            const index_t ns = numDofsAfterInitSystem(mb, 3, threeCompMapper(storage::sparse), ms);
            CHECK(storage::dense  == md);
            CHECK(storage::sparse == ms);
            CHECK(nd > 0);
            CHECK_EQUAL(nd, ns);
        }
    }

    // After localize() every query of a sparse mapper equals the dense one,
    // for many local subsets of the free dofs.
    TEST(localize_dense_equals_sparse)
    {
        const std::vector<NamedFixture> fixtures = localizeFixtures();
        const index_t shifts[] = {0, 5};
        for (size_t q = 0; q != fixtures.size(); ++q)
            for (size_t sh = 0; sh != 2; ++sh)
                for (int mode = 0; mode != 7; ++mode)
                {
                    const NamedFixture & fix = fixtures[q];
                    gsDofMapper d = buildShifted(fix, storage::dense , shifts[sh]);
                    gsDofMapper s = buildShifted(fix, storage::sparse, shifts[sh]);
                    const std::vector<index_t> l = pick(d, mode);
                    if (4 == mode && fix.hasCoupled)
                        CHECK(!l.empty());

                    d.localize(l);
                    s.localize(l);
                    CHECK_EQUAL(static_cast<index_t>(l.size()), d.freeSize());
                    CHECK_EQUAL(static_cast<index_t>(l.size()), s.freeSize());
                    compareFinalized(d, s);

                    if (3 == mode && 5 == shifts[sh])
                    {
                        gsDofMapper s2 = s;
                        compareFinalized(d, s2);
                        gsDofMapper e;
                        e.swap(s2);
                        compareFinalized(d, e);

                        const index_t nFree = d.freeSize(0);
                        if (nFree > 0)
                        {
                            gsVector<index_t> perm(nFree);
                            for (index_t i = 0; i != nFree; ++i)
                                perm[i] = nFree - 1 - i;
                            d.permuteFreeDofs(perm, 0);
                            s.permuteFreeDofs(perm, 0);
                            CHECK(storage::dense == s.storageMode());
                            compareQueries(d, s);
                        }
                    }
                }
    }

    // Localizing a localized mapper (and once more) keeps sparse equal to dense.
    TEST(relocalize_dense_equals_sparse)
    {
        const std::vector<NamedFixture> fixtures = localizeFixtures();
        const index_t shifts[] = {0, 5};
        const int firstModes[]  = {1, 3, 6};
        const int secondModes[] = {0, 2, 3, 4};
        for (size_t q = 0; q != fixtures.size(); ++q)
            for (size_t sh = 0; sh != 2; ++sh)
                for (size_t a = 0; a != 3; ++a)
                    for (size_t b = 0; b != 4; ++b)
                    {
                        gsDofMapper d = buildShifted(fixtures[q], storage::dense , shifts[sh]);
                        gsDofMapper s = buildShifted(fixtures[q], storage::sparse, shifts[sh]);

                        const std::vector<index_t> l1 = pick(d, firstModes[a]);
                        d.localize(l1);
                        s.localize(l1);

                        const std::vector<index_t> l2 = pick(d, secondModes[b]);
                        d.localize(l2);
                        s.localize(l2);
                        CHECK_EQUAL(static_cast<index_t>(l2.size()), s.freeSize());
                        compareFinalized(d, s);

                        const std::vector<index_t> l3 = pick(d, 6);
                        d.localize(l3);
                        s.localize(l3);
                        CHECK_EQUAL(static_cast<index_t>(l3.size()), s.freeSize());
                        compareFinalized(d, s);
                    }
    }

    // A space that is asked for sparse mappers builds, keeps and assembles
    // with one; the result equals the dense one.
    TEST(feSpace_setMapperStorage_end_to_end)
    {
        const PoissonPass pd = poissonPass(storage::dense);
        const PoissonPass ps = poissonPass(storage::sparse);

        CHECK(storage::dense  == pd.mapper.storageMode());
        CHECK(storage::sparse == ps.mapper.storageMode());
        CHECK(pd.K.nonZeros() > 0);
        CHECK_EQUAL(pd.K.rows(), ps.K.rows());
        CHECK_EQUAL(pd.K.cols(), ps.K.cols());
        CHECK_EQUAL(pd.K.nonZeros(), ps.K.nonZeros());
        CHECK_EQUAL(pd.F.rows(), ps.F.rows());
        CHECK_EQUAL(pd.F.cols(), ps.F.cols());

        bool samePattern = pd.K.outerSize() == ps.K.outerSize();
        for (index_t o = 0; samePattern && o != pd.K.outerSize(); ++o)
        {
            gsSparseMatrix<real_t>::InnerIterator id(pd.K, o), is(ps.K, o);
            for (; id && is; ++id, ++is)
                samePattern = samePattern && id.row() == is.row() && id.col() == is.col();
            samePattern = samePattern && !id && !is;
        }
        CHECK(samePattern);

        CHECK(!pd.pattern.empty());
        CHECK(pd.pattern == ps.pattern);
        CHECK_EQUAL(0, missingFromPattern(pd.K, pd.pattern));
        CHECK_EQUAL(0, missingFromPattern(ps.K, ps.pattern));

        const real_t tol = 1e2 * std::numeric_limits<real_t>::epsilon();
        CHECK((pd.K - ps.K).norm() <= tol * pd.K.norm());
        CHECK((pd.F - ps.F).norm() <= tol * pd.F.norm());

        gsMultiPatch<real_t> mp = gsNurbsCreator<real_t>::BSplineSquareGrid(2, 2, 1.0);
        gsMultiBasis<real_t> mb(mp);
        mb.degreeElevate(1);
        mb.uniformRefine(1);

        // initSystem() builds the mapper itself when setup() was never called
        gsExprAssembler<real_t> B(1, 1);
        B.setIntegrationElements(mb);
        auto v = B.getSpace(mb);
        v.setMapperStorage(storage::sparse);
        B.initSystem();
        CHECK(storage::sparse == v.mapper().storageMode());

        // setupMapper() installs the mapper in its own mode and leaves the request alone
        gsExprAssembler<real_t> C(1, 1);
        C.setIntegrationElements(mb);
        auto w = C.getSpace(mb);
        w.setMapperStorage(storage::sparse);
        w.setupMapper(pd.mapper);
        CHECK(storage::dense  == w.mapper().storageMode());
        CHECK(storage::sparse == w.mapperStorage());
    }

    // index_into equals the per-index index() on fixtures with a basis, for
    // dense and sparse storage, finalized and localized.
    TEST(index_into_equals_index_basis_fixtures)
    {
        const std::vector<BasisFixture> fixtures = batchFixtures();
        const index_t shifts[] = {0, 5};
        const int modes[] = {-1, 0, 2, 3, 6};
        BatchStats total;
        size_t emptyTableCases = 0;
        for (size_t q = 0; q != fixtures.size(); ++q)
            for (size_t sh = 0; sh != 2; ++sh)
                for (size_t mi = 0; mi != 5; ++mi)
                {
                    const NamedFixture nf = {fixtures[q].name, fixtures[q].build, false};
                    gsDofMapper d = buildShifted(nf, storage::dense , shifts[sh]);
                    gsDofMapper s = buildShifted(nf, storage::sparse, shifts[sh]);
                    if (modes[mi] >= 0)
                    {
                        const std::vector<index_t> l = pick(d, modes[mi]);
                        d.localize(l);
                        s.localize(l);
                    }
                    BatchStats st;
                    runBatchQueries(d, s, fixtures[q].bases, false, st);
                    checkEndsQueried(d, st);
                    CHECK(st.queried > 0);
                    CHECK(st.nonAscending > 0);
                    CHECK(st.empty > 0);
                    CHECK(st.colPositive > 0);
                    if (modes[mi] >= 0)
                        CHECK(st.remote > 0);
                    else
                        CHECK_EQUAL(size_t(0), st.remote);
                    if (0 == modes[mi])
                        ++emptyTableCases;
                    total.remote += st.remote;
                    total.nonAscending += st.nonAscending;
                }
        CHECK(emptyTableCases > 0);
        CHECK(total.remote > 0);
        CHECK(total.nonAscending > 0);
    }

    // index_into on bigMarked(16, 13): finalized, localized by pick modes,
    // and localized by bigLocalSet and then {1, 4}.
    TEST(index_into_equals_index_big_marked)
    {
        const index_t shifts[] = {0, 5};
        const std::vector<gsMultiBasis<real_t> > noBases;
        size_t emptyTableCases = 0, remote = 0, nonAscending = 0;
        for (size_t sh = 0; sh != 2; ++sh)
        {
            for (int stage = 0; stage != 7; ++stage)
            {
                gsDofMapper d = bigMarked(16, 13, storage::dense);
                gsDofMapper s = bigMarked(16, 13, storage::sparse);
                // pick() lists shifted ids; bigLocalSet() and {1, 4} are unshifted
                if (stage < 5)
                {
                    d.setShift(shifts[sh]);
                    s.setShift(shifts[sh]);
                }
                switch (stage)
                {
                case 0: break;
                case 1: case 2: case 3: case 4:
                    {
                        const int modes[] = {0, 2, 3, 6};
                        const std::vector<index_t> l = pick(d, modes[stage-1]);
                        d.localize(l);
                        s.localize(l);
                        break;
                    }
                case 5:
                    d.localize(bigLocalSet(16, 13));
                    s.localize(bigLocalSet(16, 13));
                    break;
                default:
                    d.localize(bigLocalSet(16, 13));
                    s.localize(bigLocalSet(16, 13));
                    d.localize(std::vector<index_t>{1, 4});
                    s.localize(std::vector<index_t>{1, 4});
                }
                if (stage >= 5)
                {
                    d.setShift(shifts[sh]);
                    s.setShift(shifts[sh]);
                }
                BatchStats st;
                runBatchQueries(d, s, noBases, true, st);
                checkEndsQueried(d, st);
                CHECK(st.nonAscending > 0);
                CHECK(st.empty > 0);
                CHECK(st.colPositive > 0);
                if (stage >= 1)
                    CHECK(st.remote > 0);
                if (1 == stage)
                    ++emptyTableCases;
                remote += st.remote;
                nonAscending += st.nonAscending;
            }
        }
        CHECK(emptyTableCases > 0);
        CHECK(remote > 0);
        CHECK(nonAscending > 0);
    }

    // Pattern and matrix assembled into a sink do not depend on the mapper storage.
    TEST(assembler_pattern_dense_equals_sparse)
    {
        const SinkResult ad = distinctTestSpacePass(storage::dense);
        const SinkResult as = distinctTestSpacePass(storage::sparse);
        compareSinkResults(ad, as);

        const SinkResult bd = vectorValuedPass(storage::dense);
        const SinkResult bs = vectorValuedPass(storage::sparse);
        compareSinkResults(bd, bs);
    }

    // nBytes() of a sparse mapper is bounded independently of the number of
    // positions; on Linux builds without AddressSanitizer the same bounds are
    // checked at N = 10^9 in huge_sparse_mapper_bounded_memory, where the two
    // sparse sizes are compared.
    TEST(nBytes_bounds_setup_finalize_localize)
    {
        const index_t n0 = 600000, n1 = 400000;
        const size_t denseTable = static_cast<size_t>(n0 + n1) * sizeof(index_t);
        const size_t sparseCap = 16384;

        const BigBytes d = runBigMarked(n0, n1, storage::dense);
        const BigBytes s = runBigMarked(n0, n1, storage::sparse);
        const size_t db[] = {d.setup, d.finalized, d.localized};
        const size_t sb[] = {s.setup, s.finalized, s.localized};
        for (int stage = 0; stage != 3; ++stage)
        {
            CHECK(db[stage] >= denseTable);
            CHECK(sb[stage] <= sparseCap);
            CHECK(100 * sb[stage] < db[stage]);
        }
    }

    // A sparse mapper over 10^9 positions: the numbering is the dense one
    // (closed form, validated against dense at small size) and nothing on the
    // path allocates O(N) memory.
    TEST(huge_sparse_mapper_bounded_memory)
    {
#ifdef SPARSE_TEST_ADDRESS_GUARD
        {
            const index_t n0 = 600000000, n1 = 400000000;
            CHECK(static_cast<long long>(n0) + n1 - 2 < static_cast<long long>(gsDofMapper::remoteDof()));

            AddressSpaceGuard guard;
            const bool guarded = guard.selfTest();
            CHECK(guarded);
            if (guarded)
            {
                const BigBytes huge = runBigMarked(n0, n1, storage::sparse);
                const BigBytes small = runBigMarked(600000, 400000, storage::sparse);
                CHECK_EQUAL(small.setup, huge.setup);
                CHECK_EQUAL(small.finalized, huge.finalized);
                CHECK_EQUAL(small.localized, huge.localized);
                CHECK(huge.setup <= 16384);
                CHECK(huge.finalized <= 16384);
                CHECK(huge.localized <= 16384);
            }
        }
#endif

        // The closed form used above is the dense numbering.
        const index_t n0 = 16, n1 = 13;
        const std::vector<std::pair<index_t,index_t> > all = bigAllPoints(n0, n1);
        const std::vector<index_t> l = bigLocalSet(n0, n1);
        const std::vector<index_t> l2{1, 4};
        gsDofMapper d = bigMarked(n0, n1, storage::dense);
        gsDofMapper s = bigMarked(n0, n1, storage::sparse);
        for (int stage = 0; stage != 3; ++stage)
        {
            if (1 == stage) { d.localize(l);  s.localize(l);  }
            if (2 == stage) { d.localize(l2); s.localize(l2); }
            checkBigCounts(d, stage, n0, n1);
            checkBigCounts(s, stage, n0, n1);
            checkBigAt(d, stage, n0, n1, all);
            checkBigAt(s, stage, n0, n1, all);
            compareFinalized(d, s);
        }
    }
}
