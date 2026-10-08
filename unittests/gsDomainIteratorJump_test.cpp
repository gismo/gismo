/** @file gsDomainIteratorJump_test.cpp

    @brief Tests the O(D) jump of gsTensorDomainIterator::next(k), its use through
    gsCompositeDomain and gsIndexSubDomain, the serial alltoall/alltoallv and the
    MPI registration of (unsigned) long long.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include "gismo_unittest.h"

#include <gsDomain/gsTensorDomain.h>
#include <gsDomain/gsCompositeDomain.h>
#include <gsDomain/gsIndexSubDomain.h>

#include <algorithm>
#include <vector>

SUITE(gsDomainIteratorJump_test)
{
    // The reference for an iterator jump J = S; J += k is a copy of S advanced
    // by min(k, N-s) single steps (never past the end state), where s is the
    // position of S and N the number of elements. All comparisons below are
    // exact: both sides read the same knot values.

    /// Counters of the examined cases; used to prove that no check was vacuous.
    struct JumpStats
    {
        long total, inRange, pastEnd, continuation, bad;
        JumpStats() : total(0), inRange(0), pastEnd(0), continuation(0), bad(0) { }
    };

    /// Records a mismatch; the first one of a counter is printed with its context.
    static void noteBad(JumpStats & st, const char * fixture, size_t s, index_t k, const char * what)
    {
        if (0 == st.bad)
            gsWarn << "first mismatch: fixture " << fixture << ", start s=" << s
                   << ", jump k=" << k << ", category " << what << "\n";
        ++st.bad;
    }

    static bool sameVec(const gsVector<real_t> & a, const gsVector<real_t> & b)
    {
        return a.size() == b.size() && (a.array() == b.array()).all();
    }

    static gsKnotVector<real_t> kvA()
    { return gsKnotVector<real_t>(std::vector<real_t>{0, 0, 0, 0.1, 0.35, 0.35, 0.7, 1, 1, 1}); }

    static gsKnotVector<real_t> kvB()
    { return gsKnotVector<real_t>(std::vector<real_t>{-1, -1, -0.2, 0.5, 0.5, 2, 2}); }

    static gsKnotVector<real_t> kvC()
    { return gsKnotVector<real_t>(std::vector<real_t>{0, 0, 1, 1}); }

    /// Iterator at position s, reached by s single steps from beginAll().
    static gsDomain<real_t>::iterator advancedBySteps(const gsDomain<real_t> & dom, size_t s)
    {
        gsDomain<real_t>::iterator it = dom.beginAll();
        for (size_t j = 0; j < s; ++j)
            ++it;
        return it;
    }

    // ------------------------------------------------------------------
    // Tensor domains
    // ------------------------------------------------------------------

    /// Compares J = S+k with the reference for all starts s and jumps k in [0,N+2].
    template<int D>
    static void checkTensorJump(const gsDomain<real_t> & dom, const char * name,
                                size_t expectedN, JumpStats & st)
    {
        CHECK_EQUAL(expectedN, dom.numElements());
        gsDomain<real_t>::iterator it0 = dom.beginAll();
        CHECK((nullptr != dynamic_cast<gsTensorDomainIterator<real_t,D>*>(it0.get())));

        const size_t N = dom.numElements();
        gsDomain<real_t>::iterator S = dom.beginAll();
        for (size_t s = 0; s < N; ++s, ++S)
        {
            for (index_t k = 0; k <= static_cast<index_t>(N) + 2; ++k)
            {
                gsDomain<real_t>::iterator J = S;
                J += k;
                gsDomain<real_t>::iterator O = S;
                const size_t steps = std::min(static_cast<size_t>(k), N - s);
                for (size_t j = 0; j < steps; ++j)
                    ++O;

                ++st.total;
                const bool in = (s + k < N);
                if (J.id() != s + k || O.id() != s + steps)
                    noteBad(st, name, s, k, "id");
                if ((J.id() < N) != in)
                    noteBad(st, name, s, k, "in-range flag");
                if (!sameVec(J.lowerCorner(), O.lowerCorner()))
                    noteBad(st, name, s, k, "lowerCorner");
                if (J.isBoundaryElement() != O.isBoundaryElement())
                    noteBad(st, name, s, k, "isBoundaryElement");

                if (!in)
                {
                    ++st.pastEnd;
                    continue;
                }
                ++st.inRange;
                if (!sameVec(J.upperCorner(), O.upperCorner()))
                    noteBad(st, name, s, k, "upperCorner");
                if (dom.elementIndex(0, J.centerPoint()) != static_cast<index_t>(s + k))
                    noteBad(st, name, s, k, "elementIndex oracle");

                // The state left by the jump must continue like the reference.
                while (J.id() < N && O.id() < N)
                {
                    ++st.continuation;
                    if (!sameVec(J.lowerCorner(), O.lowerCorner()) ||
                        !sameVec(J.upperCorner(), O.upperCorner()))
                        noteBad(st, name, s, k, "continuation corners");
                    ++J;
                    ++O;
                }
                if (J.id() != N || O.id() != N ||
                    !sameVec(J.lowerCorner(), O.lowerCorner()))
                    noteBad(st, name, s, k, "continuation end state");
            }
        }
    }

    /// The past-the-end state reached by a jump is absorbing and equals the one
    /// reached by single steps. Only lowerCorner() and isBoundaryElement() are
    /// valid at that state.
    template<int D>
    static void checkTensorAbsorbing(const gsDomain<real_t> & dom, const char * name,
                                     size_t expectedN, JumpStats & st)
    {
        CHECK_EQUAL(expectedN, dom.numElements());
        const size_t N = dom.numElements();

        gsDomain<real_t>::iterator ref = advancedBySteps(dom, N);
        const gsVector<real_t> refLower = ref.lowerCorner();
        const bool refBdr = ref.isBoundaryElement();

        std::vector<size_t> starts;
        starts.push_back(0);
        starts.push_back(1);
        starts.push_back(N / 2);
        starts.push_back(N - 1);

        for (size_t q = 0; q < starts.size(); ++q)
        {
            const size_t s = starts[q];
            gsDomain<real_t>::iterator E = advancedBySteps(dom, s);
            E += static_cast<index_t>(N - s);
            ++st.total;
            if (E.id() != N || !sameVec(E.lowerCorner(), refLower) || E.isBoundaryElement() != refBdr)
                noteBad(st, name, s, static_cast<index_t>(N - s), "end state differs from stepping");

            for (index_t k = 1; k <= static_cast<index_t>(N) + 2; ++k)
            {
                gsDomain<real_t>::iterator J = E;
                J += k;
                ++st.pastEnd;
                if (J.id() != N + k)
                    noteBad(st, name, s, k, "id at end");
                if (!sameVec(J.lowerCorner(), refLower))
                    noteBad(st, name, s, k, "end lowerCorner changed");
                if (J.isBoundaryElement() != refBdr)
                    noteBad(st, name, s, k, "end isBoundaryElement changed");
            }
        }
    }

    /// mode 0: jump versus stepping; mode 1: past-the-end state (only if \a absorb).
    template<int D>
    static void visitTensor(int mode, const gsDomain<real_t> & dom, const char * name,
                            size_t N, bool absorb, JumpStats & st)
    {
        if (0 == mode)
            checkTensorJump<D>(dom, name, N, st);
        else if (absorb)
            checkTensorAbsorbing<D>(dom, name, N, st);
    }

    template<int D>
    static void visitDirect(int mode, const char * name, size_t N, bool absorb,
                            const std::vector<gsKnotVector<real_t> > & kvs, JumpStats & st)
    {
        std::vector<gsDomain<real_t>::Ptr> v;
        for (size_t i = 0; i < kvs.size(); ++i)
            v.push_back(memory::make_shared(new gsKnotVector<real_t>(kvs[i])));
        gsTensorDomain<real_t,D> dom(v);
        visitTensor<D>(mode, dom, name, N, absorb, st);
    }

    static void visitAllTensorFixtures(int mode, JumpStats & st)
    {
        std::vector<gsKnotVector<real_t> > kvs;

        kvs.assign(1, kvA());
        visitDirect<1>(mode, "direct D=1 {A}", 4, false, kvs, st);
        kvs.assign(1, kvC());
        visitDirect<1>(mode, "direct D=1 {C}", 1, false, kvs, st);

        {
            gsTensorBSplineBasis<2,real_t> b(kvA(), kvB());
            gsDomain<real_t>::Ptr dp = b.domain();
            visitTensor<2>(mode, *dp, "basis D=2 (A,B)", 12, true, st);
        }
        kvs.clear(); kvs.push_back(kvC()); kvs.push_back(kvA());
        visitDirect<2>(mode, "direct D=2 {C,A}", 4, false, kvs, st);
        kvs.clear(); kvs.push_back(kvB()); kvs.push_back(kvC());
        visitDirect<2>(mode, "direct D=2 {B,C}", 3, true, kvs, st);

        {
            gsTensorBSplineBasis<3,real_t> b(kvA(), kvB(), kvC());
            gsDomain<real_t>::Ptr dp = b.domain();
            visitTensor<3>(mode, *dp, "basis D=3 (A,B,C)", 12, true, st);
        }
        {
            gsTensorBSplineBasis<3,real_t> b(kvB(), kvC(), kvA());
            gsDomain<real_t>::Ptr dp = b.domain();
            visitTensor<3>(mode, *dp, "basis D=3 (B,C,A)", 12, true, st);
        }
        kvs.clear(); kvs.push_back(kvA()); kvs.push_back(kvA()); kvs.push_back(kvB());
        visitDirect<3>(mode, "direct D=3 {A,A,B}", 48, true, kvs, st);
        kvs.clear(); kvs.push_back(kvC()); kvs.push_back(kvA()); kvs.push_back(kvB());
        visitDirect<3>(mode, "direct D=3 {C,A,B}", 12, false, kvs, st);
    }

    TEST(TensorJumpMatchesStepping)
    {
        JumpStats st;
        visitAllTensorFixtures(0, st);
        CHECK(st.total > 0);
        CHECK(st.inRange > 0);
        CHECK(st.pastEnd > 0);
        CHECK(st.continuation > 0);
        CHECK_EQUAL(0, st.bad);
    }

    TEST(TensorJumpFromEndIsAbsorbing)
    {
        JumpStats st;
        visitAllTensorFixtures(1, st);
        CHECK(st.total > 0);
        CHECK(st.pastEnd > 0);
        CHECK_EQUAL(0, st.bad);
    }

    TEST(TensorJumpCostIsNotLinear)
    {
        // 1000^3 = 1e9 elements: below INT_MAX, but a loop of single steps would take seconds.
        std::vector<gsDomain<real_t>::Ptr> v;
        for (int i = 0; i < 3; ++i)
            v.push_back(memory::make_shared(new gsKnotVector<real_t>(0.0, 1.0, 999u, 3)));
        gsTensorDomain<real_t,3> dom(v);
        const size_t N = dom.numElements();
        CHECK_EQUAL(static_cast<size_t>(1000000000), N);

        gsDomain<real_t>::iterator it = dom.beginAll();
        gsStopwatch sw;
        it += static_cast<index_t>(N - 1);
        const double t = sw.stop();
        gsInfo << "jump by " << N - 1 << " elements took " << t << " s\n";

        CHECK(t < 0.1);
        CHECK_EQUAL(N - 1, it.id());
        CHECK_EQUAL(static_cast<index_t>(N - 1), dom.elementIndex(0, it.centerPoint()));
    }

    // ------------------------------------------------------------------
    // Multipatch domains
    // ------------------------------------------------------------------

    static void buildMb2(gsMultiBasis<real_t> & mb)
    {
        mb.addBasis(new gsTensorBSplineBasis<2,real_t>(kvA(), kvB()));
        mb.addBasis(new gsTensorBSplineBasis<2,real_t>(kvC(), kvC()));
        mb.addBasis(new gsTensorBSplineBasis<2,real_t>(kvB(), kvA()));
        mb.addBasis(new gsTensorBSplineBasis<2,real_t>(kvA(), kvC()));
    }

    static void buildMb3(gsMultiBasis<real_t> & mb)
    {
        mb.addBasis(new gsTensorBSplineBasis<3,real_t>(kvA(), kvB(), kvC()));
        mb.addBasis(new gsTensorBSplineBasis<3,real_t>(kvC(), kvB(), kvA()));
        mb.addBasis(new gsTensorBSplineBasis<3,real_t>(kvB(), kvB(), kvC()));
    }

    /// Cumulative element counts per patch, computed from the bases alone.
    static std::vector<size_t> patchOffsets(const gsMultiBasis<real_t> & mb)
    {
        std::vector<size_t> off(1, 0);
        for (size_t p = 0; p < mb.nBases(); ++p)
            off.push_back(off.back() + mb.basis(p).numElements());
        return off;
    }

    /// Patch p with off[p] <= e < off[p+1].
    static size_t patchOf(const std::vector<size_t> & off, size_t e)
    {
        return static_cast<size_t>(std::upper_bound(off.begin(), off.end(), e) - off.begin()) - 1;
    }

    static void checkCompositeJump(const gsDomain<real_t> & dom, const std::vector<size_t> & off,
                                   const char * name, size_t expectedN, JumpStats & st)
    {
        CHECK_EQUAL(expectedN, off.back());
        CHECK_EQUAL(expectedN, dom.numElements());
        gsDomain<real_t>::iterator it0 = dom.beginAll();
        CHECK(nullptr != dynamic_cast<gsCompositeDomainIterator<real_t>*>(it0.get()));

        const size_t N = dom.numElements();
        gsDomain<real_t>::iterator S = dom.beginAll();
        for (size_t s = 0; s < N; ++s, ++S)
        {
            for (index_t k = 0; k <= static_cast<index_t>(N) + 2; ++k)
            {
                gsDomain<real_t>::iterator J = S;
                J += k;
                gsDomain<real_t>::iterator O = S;
                const size_t steps = std::min(static_cast<size_t>(k), N - s);
                for (size_t j = 0; j < steps; ++j)
                    ++O;

                ++st.total;
                if (J.id() != s + k || O.id() != s + steps)
                    noteBad(st, name, s, k, "id");
                if ((J.id() < N) != (s + k < N))
                    noteBad(st, name, s, k, "in-range flag");
                if (s + k >= N)
                {
                    // The jump and the single steps leave different internal
                    // states at the end; only the id is meaningful there.
                    ++st.pastEnd;
                    continue;
                }
                ++st.inRange;

                const size_t e = s + k;
                const size_t p = patchOf(off, e);
                if (!sameVec(J.lowerCorner(), O.lowerCorner()) ||
                    !sameVec(J.upperCorner(), O.upperCorner()))
                    noteBad(st, name, s, k, "corners");
                if (J.patchIndex() != static_cast<index_t>(p) ||
                    J.localId() != static_cast<index_t>(e - off[p]))
                    noteBad(st, name, s, k, "patchIndex/localId");
                if (J.subdomainIndex() != O.subdomainIndex())
                    noteBad(st, name, s, k, "subdomainIndex");
                if (dom.elementIndex(J.patchIndex(), J.centerPoint()) != static_cast<index_t>(e))
                    noteBad(st, name, s, k, "elementIndex oracle");

                while (J.id() < N && O.id() < N)
                {
                    ++st.continuation;
                    if (!sameVec(J.lowerCorner(), O.lowerCorner()) ||
                        !sameVec(J.upperCorner(), O.upperCorner()))
                        noteBad(st, name, s, k, "continuation corners");
                    ++J;
                    ++O;
                }
                if (J.id() != N || O.id() != N)
                    noteBad(st, name, s, k, "continuation end id");
            }
        }
    }

    TEST(CompositeJumpMatchesStepping)
    {
        JumpStats st;
        {
            gsMultiBasis<real_t> mb;
            buildMb2(mb);
            const std::vector<size_t> off = patchOffsets(mb);
            CHECK(off == (std::vector<size_t>{0, 12, 13, 25, 29}));
            gsDomain<real_t>::Ptr dom = mb.domain();
            checkCompositeJump(*dom, off, "multipatch 2-D", 29, st);
        }
        {
            gsMultiBasis<real_t> mb;
            buildMb3(mb);
            const std::vector<size_t> off = patchOffsets(mb);
            CHECK(off == (std::vector<size_t>{0, 12, 24, 33}));
            gsDomain<real_t>::Ptr dom = mb.domain();
            checkCompositeJump(*dom, off, "multipatch 3-D", 33, st);
        }
        CHECK(st.total > 0);
        CHECK(st.inRange > 0);
        CHECK(st.pastEnd > 0);
        CHECK(st.continuation > 0);
        CHECK_EQUAL(0, st.bad);
    }

    // ------------------------------------------------------------------
    // gsIndexSubDomain
    // ------------------------------------------------------------------

    /// Data of one parent element, obtained by single steps through the parent.
    struct ElemInfo
    {
        gsVector<real_t> lower, upper;
        index_t patch, local;
    };

    static std::vector<ElemInfo> parentTable(const gsDomain<real_t> & dom)
    {
        std::vector<ElemInfo> tab;
        const size_t N = dom.numElements();
        gsDomain<real_t>::iterator it = dom.beginAll();
        for (size_t e = 0; e < N; ++e, ++it)
        {
            ElemInfo ei;
            ei.lower = it.lowerCorner();
            ei.upper = it.upperCorner();
            ei.patch = it.patchIndex();
            ei.local = static_cast<index_t>(it.localId());
            tab.push_back(ei);
        }
        return tab;
    }

    /// Strictly increasing list with a gap: e % 3 != 1, plus the last element.
    static std::vector<index_t> listL1(size_t N)
    {
        std::vector<index_t> l;
        for (size_t e = 0; e < N; ++e)
            if (e % 3 != 1 || e + 1 == N)
                l.push_back(static_cast<index_t>(e));
        return l;
    }

    /// Strictly increasing list with a gap, starting at 2: e % 4 >= 2.
    static std::vector<index_t> listL2(size_t N)
    {
        std::vector<index_t> l;
        for (size_t e = 0; e < N; ++e)
            if (e % 4 >= 2)
                l.push_back(static_cast<index_t>(e));
        return l;
    }

    static bool sameAsParent(const gsDomain<real_t>::iterator & it, const ElemInfo & ei)
    {
        return sameVec(it.lowerCorner(), ei.lower) && sameVec(it.upperCorner(), ei.upper) &&
            it.patchIndex() == ei.patch && static_cast<index_t>(it.localId()) == ei.local;
    }

    static void checkListProperties(const std::vector<index_t> & idx, const std::vector<ElemInfo> & tab)
    {
        const size_t N = tab.size();
        CHECK(idx.size() < N);
        CHECK(!idx.empty());
        CHECK(std::adjacent_find(idx.begin(), idx.end(), std::greater_equal<index_t>()) == idx.end());
        CHECK(tab[idx.front()].patch != tab[idx.back()].patch);
    }

    static void checkSubDomainWalk(const gsDomain<real_t>::Ptr & dom, const std::vector<ElemInfo> & tab,
                                   const std::vector<index_t> & idx, const char * name, JumpStats & st)
    {
        checkListProperties(idx, tab);
        gsIndexSubDomain<real_t> sub(dom, idx);
        const size_t M = idx.size();
        CHECK_EQUAL(M, sub.numElements());

        gsDomain<real_t>::iterator it0 = sub.beginAll();
        CHECK(nullptr != dynamic_cast<gsIndexSubDomainIterator<real_t>*>(it0.get()));

        const gsDomain<real_t>::iterator last = sub.endAll();
        size_t j = 0;
        for (gsDomain<real_t>::iterator it = sub.beginAll(); it != last && j <= M; ++it, ++j)
        {
            ++st.total;
            ++st.inRange;
            if (it.id() != j)
                noteBad(st, name, j, 1, "id");
            else if (!sameAsParent(it, tab[idx[j]]))
                noteBad(st, name, j, 1, "walk element");
        }
        if (j != M)
            noteBad(st, name, j, 1, "number of visited elements");
    }

    static void checkSubDomainJump(const gsDomain<real_t>::Ptr & dom, const std::vector<ElemInfo> & tab,
                                   const std::vector<index_t> & idx, const char * name, JumpStats & st)
    {
        gsIndexSubDomain<real_t> sub(dom, idx);
        const size_t M = idx.size();
        gsDomain<real_t>::iterator S = sub.beginAll();
        for (size_t s = 0; s < M; ++s, ++S)
        {
            for (index_t k = 0; k <= static_cast<index_t>(M) + 2; ++k)
            {
                gsDomain<real_t>::iterator J = S;
                J += k;
                ++st.total;
                if (J.id() != s + k)
                    noteBad(st, name, s, k, "id");
                if (s + k >= M)
                {
                    ++st.pastEnd;
                    continue;
                }
                ++st.inRange;
                if (!sameAsParent(J, tab[idx[s + k]]))
                    noteBad(st, name, s, k, "jumped element");
            }
        }
    }

    TEST(IndexSubDomainVisitsListedElements)
    {
        JumpStats st;
        for (int d = 2; d <= 3; ++d)
        {
            gsMultiBasis<real_t> mb;
            if (2 == d) buildMb2(mb); else buildMb3(mb);
            gsDomain<real_t>::Ptr dom = mb.domain();
            const std::vector<ElemInfo> tab = parentTable(*dom);
            const size_t N = tab.size();
            checkSubDomainWalk(dom, tab, listL1(N), 2 == d ? "2-D L1" : "3-D L1", st);
            checkSubDomainWalk(dom, tab, listL2(N), 2 == d ? "2-D L2" : "3-D L2", st);
        }
        CHECK(st.total > 0);
        CHECK_EQUAL(0, st.bad);
    }

    TEST(IndexSubDomainJumpMatchesStepping)
    {
        JumpStats st;
        for (int d = 2; d <= 3; ++d)
        {
            gsMultiBasis<real_t> mb;
            if (2 == d) buildMb2(mb); else buildMb3(mb);
            gsDomain<real_t>::Ptr dom = mb.domain();
            const std::vector<ElemInfo> tab = parentTable(*dom);
            const size_t N = tab.size();
            checkSubDomainJump(dom, tab, listL1(N), 2 == d ? "2-D L1" : "3-D L1", st);
            checkSubDomainJump(dom, tab, listL2(N), 2 == d ? "2-D L2" : "3-D L2", st);
        }
        CHECK(st.total > 0);
        CHECK(st.inRange > 0);
        CHECK(st.pastEnd > 0);
        CHECK_EQUAL(0, st.bad);
    }

    // ------------------------------------------------------------------
    // Serial alltoall / alltoallv
    // ------------------------------------------------------------------

    TEST(SerialAlltoallvHonoursDisplacements)
    {
        index_t send[8] = {10, 11, 12, 13, 14, 15, 16, 17};
        index_t recv[12];
        std::fill(recv, recv + 12, -7);
        int sc = 3, sd = 2, rc = 3, rd = 5;
        CHECK_EQUAL(0, gsSerialComm::alltoallv(send, &sc, &sd, recv, &rc, &rd));
        for (int i = 0; i < 12; ++i)
        {
            if (i >= 5 && i < 8)
                CHECK_EQUAL(12 + (i - 5), recv[i]);
            else
                CHECK_EQUAL(-7, recv[i]);
        }

        // Through an instance, with a 64-bit element type
        long long lsend[8] = {1LL << 40, 11, 12, 13, 14, 15, 16, -17};
        long long lrecv[12];
        std::fill(lrecv, lrecv + 12, -7LL);
        gsSerialComm ser;
        CHECK_EQUAL(0, ser.alltoallv(lsend, &sc, &sd, lrecv, &rc, &rd));
        for (int i = 0; i < 12; ++i)
        {
            if (i >= 5 && i < 8)
                CHECK_EQUAL(lsend[2 + (i - 5)], lrecv[i]);
            else
                CHECK_EQUAL(-7LL, lrecv[i]);
        }

        // Empty block: nothing is written
        index_t zrecv[12];
        std::fill(zrecv, zrecv + 12, -7);
        int zero = 0;
        CHECK_EQUAL(0, gsSerialComm::alltoallv(send, &zero, &sd, zrecv, &zero, &rd));
        for (int i = 0; i < 12; ++i)
            CHECK_EQUAL(-7, zrecv[i]);

        // A size-1 communicator behaves like the serial one
        gsMpi::init();
        gsMpiComm comm(gsMpi::localComm());
        CHECK_EQUAL(1, comm.size());
        index_t crecv[12];
        std::fill(crecv, crecv + 12, -7);
        CHECK_EQUAL(0, comm.alltoallv(send, &sc, &sd, crecv, &rc, &rd));
        CHECK(std::equal(recv, recv + 12, crecv));
    }

    TEST(SerialAlltoallCopiesSelfBlock)
    {
        index_t send[6] = {1, 2, 3, 4, 5, 6};
        index_t recv[8];
        std::fill(recv, recv + 8, -7);
        CHECK_EQUAL(0, gsSerialComm::alltoall(send, recv, 4, 4));
        for (int i = 0; i < 8; ++i)
            CHECK_EQUAL(i < 4 ? i + 1 : -7, recv[i]);

        gsMpi::init();
        gsMpiComm comm(gsMpi::localComm());
        CHECK_EQUAL(1, comm.size());
        index_t crecv[8];
        std::fill(crecv, crecv + 8, -7);
        CHECK_EQUAL(0, comm.alltoall(send, crecv, 4, 4));
        CHECK(std::equal(recv, recv + 8, crecv));
    }

#ifdef GISMO_WITH_MPI
    // ------------------------------------------------------------------
    // 64-bit integer registration with MPI
    // ------------------------------------------------------------------

    TEST(LongLongNativeMpiRegistration)
    {
        gsMpi::init();
        CHECK(MPITraits<long long>::getType() == MPI_LONG_LONG);
        CHECK(MPITraits<unsigned long long>::getType() == MPI_UNSIGNED_LONG_LONG);

        CHECK((Generic_MPI_Op<long long, std::plus<long long> >::get() == MPI_SUM));
        CHECK((Generic_MPI_Op<long long, std::multiplies<long long> >::get() == MPI_PROD));
        CHECK((Generic_MPI_Op<long long, _mpi_Min<long long> >::get() == MPI_MIN));
        CHECK((Generic_MPI_Op<long long, _mpi_Max<long long> >::get() == MPI_MAX));
        CHECK((Generic_MPI_Op<unsigned long long, std::plus<unsigned long long> >::get() == MPI_SUM));
        CHECK((Generic_MPI_Op<unsigned long long, std::multiplies<unsigned long long> >::get() == MPI_PROD));
        CHECK((Generic_MPI_Op<unsigned long long, _mpi_Min<unsigned long long> >::get() == MPI_MIN));
        CHECK((Generic_MPI_Op<unsigned long long, _mpi_Max<unsigned long long> >::get() == MPI_MAX));
    }

    template<typename T>
    static void checkSizeOneReductions(const T (&ref)[3])
    {
        gsMpiComm comm(gsMpi::localComm());
        CHECK_EQUAL(1, comm.size());
        T buf[3];

        std::copy(ref, ref + 3, buf);
        CHECK_EQUAL(0, comm.sum(buf, 3));
        CHECK(std::equal(ref, ref + 3, buf));
        CHECK_EQUAL(0, comm.min(buf, 3));
        CHECK(std::equal(ref, ref + 3, buf));
        CHECK_EQUAL(0, comm.max(buf, 3));
        CHECK(std::equal(ref, ref + 3, buf));
    }

    TEST(LongLongReductionsOnSizeOneComm)
    {
        gsMpi::init();
        const long long ll[3] = {1, 1LL << 40, -5};
        checkSizeOneReductions(ll);
        const unsigned long long ull[3] = {0, 1ULL << 63, 7};
        checkSizeOneReductions(ull);
        const int64_t i64[3] = {-1, 1LL << 50, 3};
        checkSizeOneReductions(i64);
    }
#endif

} // SUITE(gsDomainIteratorJump_test)
