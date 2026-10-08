/** @file gsDistributedPartition_test.cpp

    @brief Tests the distributed-partition building blocks in gsDistributedPartition.h

    - OrderedKeyDouble, OrderedKeyFloat: gsOrderedKey preserves the order and
      the equality of IEEE values (-0.0 == +0.0).
    - OrderedKeyLongDoubleUnavailable: only float and double have a key.
    - RadixSelectCount, RadixSelectWeighted: gsRadixSelect driven through
      gsAsMachine and runLockstep matches a std::sort oracle.
    - LockstepRejectsLengthMismatch, LockstepRejectsKindMismatch,
      LockstepRejectsEarlyDone: runLockstep enforces the collective protocol.
    - DistributedRcb*: the distributed RCB machine, run in lockstep, gives
      every rank exactly the element set that rcbSplit labels with that rank.
    - RunCollectivesSizeOneCommGivesAllIds: the communicator driver on a
      size-1 communicator.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
**/

#include "gismo_unittest.h"
#include <gsDomain/gsDistributedPartition.h>
#include <gsAssembler/gsDofMapperCreator.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <random>
#include <set>
#include <vector>

SUITE(gsDistributedPartition_test)
{
    // ------------------------------------------------------------------
    // Order-preserving keys
    // ------------------------------------------------------------------

    // True iff the key of a and b mirrors a < b and a == b.
    template<class T>
    static bool keyPairOk(T a, T b)
    {
        const uint64_t ea = gsOrderedKey<T>::encode(a), eb = gsOrderedKey<T>::encode(b);
        return ((a < b) == (ea < eb)) && ((a == b) == (ea == eb));
    }

    static double doubleFromBits(std::uint64_t b)
    {
        double x;
        std::memcpy(&x, &b, sizeof(x));
        return x;
    }

    static std::uint64_t doubleBits(double x)
    {
        std::uint64_t b;
        std::memcpy(&b, &x, sizeof(b));
        return b;
    }

    static std::uint32_t floatBits(float x)
    {
        std::uint32_t b;
        std::memcpy(&b, &x, sizeof(b));
        return b;
    }

    static float floatFromBits(std::uint32_t b)
    {
        float x;
        std::memcpy(&x, &b, sizeof(x));
        return x;
    }

    TEST(OrderedKeyDouble)
    {
        typedef std::numeric_limits<double> L;
        std::vector<double> sp;
        sp.push_back(-0.0);
        sp.push_back(0.0);
        sp.push_back(L::denorm_min());
        sp.push_back(-L::denorm_min());
        sp.push_back(doubleFromBits(0x0000123456789ABCULL));
        sp.push_back(-doubleFromBits(0x0000123456789ABCULL));
        sp.push_back(L::min());
        sp.push_back(-L::min());
        sp.push_back(1.0);
        sp.push_back(-1.0);
        sp.push_back(L::max());
        sp.push_back(L::lowest());
        sp.push_back(L::infinity());
        sp.push_back(-L::infinity());

        index_t cases = 0, bad = 0;
        for (size_t i = 0; i < sp.size(); ++i)
            for (size_t j = 0; j < sp.size(); ++j)
            {
                ++cases;
                if (!keyPairOk(sp[i], sp[j]))
                {
                    if (bad++ == 0)
                        gsWarn << "OrderedKeyDouble: special pair (" << i << ", " << j << ") disagrees\n";
                }
            }
        CHECK(gsOrderedKey<double>::encode(-0.0) == gsOrderedKey<double>::encode(0.0));

        std::mt19937_64 rng(20261007);
        for (index_t n = 0; n < 100000; ++n)
        {
            double a, b;
            do { a = doubleFromBits(rng()); } while (a != a);
            do
            {
                if (n % 5 == 0)      b = a;
                else if (n % 2 == 0) b = doubleFromBits(rng());
                else                 b = doubleFromBits(doubleBits(a) ^ (rng() & 0xFFULL));
            } while (b != b);
            ++cases;
            if (!keyPairOk(a, b))
            {
                if (bad++ == 0)
                    gsWarn << "OrderedKeyDouble: random pair " << a << ", " << b << " disagrees\n";
            }
        }
        gsInfo << "gsDistributedPartition OrderedKeyDouble: " << cases << " pair cases, " << bad << " mismatches\n";
        CHECK(cases > 0);
        CHECK_EQUAL(0, bad);
    }

    TEST(OrderedKeyFloat)
    {
        typedef std::numeric_limits<float> L;
        std::vector<float> sp;
        sp.push_back(-0.0f);
        sp.push_back(0.0f);
        sp.push_back(L::denorm_min());
        sp.push_back(-L::denorm_min());
        sp.push_back(floatFromBits(0x00123456U));
        sp.push_back(-floatFromBits(0x00123456U));
        sp.push_back(L::min());
        sp.push_back(-L::min());
        sp.push_back(1.0f);
        sp.push_back(-1.0f);
        sp.push_back(L::max());
        sp.push_back(L::lowest());
        sp.push_back(L::infinity());
        sp.push_back(-L::infinity());

        index_t cases = 0, bad = 0, badHi = 0;
        for (size_t i = 0; i < sp.size(); ++i)
        {
            if ((gsOrderedKey<float>::encode(sp[i]) >> 32) != 0) ++badHi;
            for (size_t j = 0; j < sp.size(); ++j)
            {
                ++cases;
                if (!keyPairOk(sp[i], sp[j]))
                {
                    if (bad++ == 0)
                        gsWarn << "OrderedKeyFloat: special pair (" << i << ", " << j << ") disagrees\n";
                }
            }
        }
        CHECK(gsOrderedKey<float>::encode(-0.0f) == gsOrderedKey<float>::encode(0.0f));

        std::mt19937_64 rng(20261008);
        for (index_t n = 0; n < 100000; ++n)
        {
            float a, b;
            do { a = floatFromBits((std::uint32_t)rng()); } while (a != a);
            do
            {
                if (n % 5 == 0)      b = a;
                else if (n % 2 == 0) b = floatFromBits((std::uint32_t)rng());
                else                 b = floatFromBits(floatBits(a) ^ (std::uint32_t)(rng() & 0xFFULL));
            } while (b != b);
            ++cases;
            if (!keyPairOk(a, b))
            {
                if (bad++ == 0)
                    gsWarn << "OrderedKeyFloat: random pair " << a << ", " << b << " disagrees\n";
            }
            if ((gsOrderedKey<float>::encode(a) >> 32) != 0) ++badHi;
        }
        gsInfo << "gsDistributedPartition OrderedKeyFloat: " << cases << " pair cases, " << bad
               << " mismatches, " << badHi << " keys with nonzero upper 32 bits\n";
        CHECK(cases > 0);
        CHECK_EQUAL(0, bad);
        CHECK_EQUAL(0, badHi);
    }

    TEST(OrderedKeyLongDoubleUnavailable)
    {
        CHECK(!gsOrderedKey<long double>::available);
        CHECK(gsOrderedKey<double>::available && gsOrderedKey<float>::available);
        gsInfo << "gsDistributedPartition OrderedKeyLongDoubleUnavailable: 3 availability cases compared\n";
    }

    // ------------------------------------------------------------------
    // Radix select
    // ------------------------------------------------------------------

    typedef gsRadixSelect<double> gsSel;
    typedef gsAsMachine<gsSel>    gsSelMachine;

    // One global element: select key, candidate group and weight.
    struct gsSelElem
    {
        gsPartitionKey key;
        index_t        group;
        int64_t        w;
    };

    struct gsSelCand
    {
        gsPartitionKey key;
        int64_t        w;
    };

    static bool selCandLess(const gsSelCand & a, const gsSelCand & b) { return a.key < b.key; }

    // Reference result per target: sort the candidates by key and scan.
    static void selOracle(const std::vector<gsSelElem> & el, gsSel::Mode mode,
                          const std::vector<index_t> & tg, const std::vector<int64_t> & tp,
                          std::vector<char> & found, std::vector<gsPartitionKey> & key)
    {
        found.assign(tg.size(), 0);
        key.assign(tg.size(), gsPartitionKey());
        for (size_t t = 0; t < tg.size(); ++t)
        {
            std::vector<gsSelCand> cand;
            for (size_t e = 0; e < el.size(); ++e)
                if (el[e].group == tg[t])
                {
                    gsSelCand c; c.key = el[e].key; c.w = el[e].w;
                    cand.push_back(c);
                }
            std::sort(cand.begin(), cand.end(), selCandLess);
            if (mode == gsSel::Count)
            {
                if (tp[t] >= 0 && tp[t] < (int64_t)cand.size())
                {
                    found[t] = 1;
                    key[t]   = cand[(size_t)tp[t]].key;
                }
            }
            else
            {
                int64_t prefix = 0;
                for (size_t c = 0; c < cand.size(); ++c)
                {
                    prefix += cand[c].w;
                    if (prefix >= tp[t])
                    {
                        found[t] = 1;
                        key[t]   = cand[c].key;
                        break;
                    }
                }
            }
        }
    }

    // Distinct ids spanning the whole index_t range: 0, 255/256, 2^16, 2^24,
    // their neighbours and max(), then random ids.
    static std::vector<index_t> selSparseIds(std::mt19937_64 & rng, index_t N)
    {
        const index_t big = std::numeric_limits<index_t>::max();
        const index_t fixed[] = { 0, big, 65536, (index_t)1 << 24, 255, big - 1, 65535,
                                  ((index_t)1 << 24) - 1, 1, 256, 65537, ((index_t)1 << 24) + 1, big / 2 };
        std::set<index_t> used;
        std::vector<index_t> ids;
        for (size_t i = 0; i < sizeof(fixed) / sizeof(fixed[0]) && (index_t)ids.size() < N; ++i)
        {
            used.insert(fixed[i]);
            ids.push_back(fixed[i]);
        }
        while ((index_t)ids.size() < N)
        {
            const index_t id = (index_t)(rng() % ((std::uint64_t)big + 1));
            if (used.insert(id).second)
                ids.push_back(id);
        }
        return ids;
    }

    // xkind 0: coordinates from a duplicate-rich pool plus random values, ids 0..N-1;
    //       1: every coordinate equal, ids 0..N-1;
    //       2: every coordinate equal, sparse ids over the index_t range;
    //       3: two coordinate values, sparse ids over the index_t range.
    // wkind 0: weights in [0,5], mostly zero; 1: all zero; 2: 0/1 with one weight of 1000.
    static std::vector<gsSelElem> selMakeElems(std::mt19937_64 & rng, index_t N, int xkind,
                                               bool multi, int wkind)
    {
        const double pool[] = { -0.0, 0.0, -1.5, 1.5, 1e-300, -1e300, 3.0,
                                std::numeric_limits<double>::denorm_min() };
        const index_t groupPool[] = { 0, 1, 2, 3, 4, -1, -1 };

        std::vector<index_t> ids;
        if (xkind >= 2)
            ids = selSparseIds(rng, N);
        else
        {
            ids.resize((size_t)N);
            for (index_t i = 0; i < N; ++i) ids[(size_t)i] = i;
        }
        std::shuffle(ids.begin(), ids.end(), rng);

        std::vector<gsSelElem> el((size_t)N);
        for (size_t i = 0; i < el.size(); ++i)
        {
            double x;
            if (xkind == 0)
                x = (rng() & 1) ? pool[rng() % 8]
                                : (double)((std::int64_t)(rng() % 20001) - 10000) / 100.0;
            else if (xkind == 3)
                x = (rng() & 1) ? -1.5 : 1.5;
            else
                x = 1.5;
            el[i].key.hi = gsOrderedKey<double>::encode(x);
            el[i].key.lo = (gsPartitionKey::lo_type)ids[i];
            el[i].group  = multi ? groupPool[rng() % 7] : 0;
            if (wkind == 0)      el[i].w = (rng() % 3 == 0) ? (int64_t)(rng() % 6) : 0;
            else if (wkind == 1) el[i].w = 0;
            else                 el[i].w = (int64_t)(rng() % 2);
        }
        if (wkind == 2 && N > 0)
            el[(size_t)(rng() % (std::uint64_t)N)].w = 1000;
        return el;
    }

    static void selMakeTargets(std::mt19937_64 & rng, const std::vector<gsSelElem> & el,
                               gsSel::Mode mode, bool multi, index_t P,
                               std::vector<index_t> & tg, std::vector<int64_t> & tp)
    {
        const int G = 7;
        std::vector<int64_t> cnt(G, 0), wsum(G, 0);
        for (size_t e = 0; e < el.size(); ++e)
            if (el[e].group >= 0 && el[e].group < G)
            {
                ++cnt[(size_t)el[e].group];
                wsum[(size_t)el[e].group] += el[e].w;
            }
        tg.clear(); tp.clear();
        if (!multi)
        {
            const int64_t n = (int64_t)el.size(), W = wsum[0];
            if (mode == gsSel::Count)
            {
                tp.push_back(0);
                tp.push_back(n - 1);
                tp.push_back(n > 0 ? (int64_t)(rng() % (std::uint64_t)n) : 0);
                tp.push_back(n);
                tp.push_back(-1);
            }
            else
            {
                tp.push_back(-3);
                tp.push_back(0);
                tp.push_back(1);
                tp.push_back(W >= 1 ? 1 + (int64_t)(rng() % (std::uint64_t)W) : 1);
                tp.push_back(W);
                tp.push_back(W + 1);
                for (index_t j = 1; j < P; ++j)
                    tp.push_back((W * j + P - 1) / P);
            }
            tg.assign(tp.size(), 0);
        }
        else
        {
            for (int pass = 0; pass < 2; ++pass)
                for (int g = 0; g < G; ++g)
                {
                    tg.push_back(g);
                    if (mode == gsSel::Count)
                        tp.push_back(pass == 0 ? (int64_t)(rng() % (std::uint64_t)(cnt[(size_t)g] + 2)) - 1 : 0);
                    else
                        tp.push_back(pass == 0 ? (int64_t)(rng() % (std::uint64_t)(wsum[(size_t)g] + 4)) - 2
                                               : wsum[(size_t)g]);
                }
        }
    }

    // Runs one select on P virtual ranks and compares every rank with the oracle.
    // The elements are spread over the ranks in slices N*r/P of the given order,
    // or all sit on the last rank (skewed).
    static void selectCase(const char * what, const std::vector<gsSelElem> & el, index_t P, bool skewed,
                           bool nullGroup, gsSel::Mode mode, const std::vector<index_t> & tg,
                           const std::vector<int64_t> & tp, bool vecOverload,
                           index_t & compared, index_t & bad, index_t & maxIssued)
    {
        const std::int64_t N = (std::int64_t)el.size();
        std::vector<std::vector<gsPartitionKey> > keys((size_t)P);
        std::vector<std::vector<index_t> >        grp((size_t)P);
        std::vector<std::vector<int64_t> >        wts((size_t)P);
        for (index_t r = 0; r < P; ++r)
        {
            const std::int64_t lo = skewed ? (r == P - 1 ? 0 : N) : N * r / P;
            const std::int64_t hi = skewed ? N : N * (r + 1) / P;
            for (std::int64_t e = lo; e < hi; ++e)
            {
                keys[(size_t)r].push_back(el[(size_t)e].key);
                grp[(size_t)r].push_back(el[(size_t)e].group);
                wts[(size_t)r].push_back(el[(size_t)e].w);
            }
        }

        std::vector<gsSel> sels((size_t)P);
        for (index_t r = 0; r < P; ++r)
            sels[(size_t)r].setup(mode, (index_t)keys[(size_t)r].size(), keys[(size_t)r].data(),
                                  nullGroup ? NULL : grp[(size_t)r].data(),
                                  mode == gsSel::Weighted ? wts[(size_t)r].data() : NULL, tg, tp);

        std::vector<gsSelMachine> ms;
        ms.reserve((size_t)P);
        for (index_t r = 0; r < P; ++r)
            ms.push_back(gsSelMachine(sels[(size_t)r]));
        index_t steps;
        if (vecOverload)
            steps = runLockstep(ms);
        else
        {
            std::vector<gsSelMachine *> ptrs((size_t)P);
            for (index_t r = 0; r < P; ++r) ptrs[(size_t)r] = &ms[(size_t)r];
            steps = runLockstep(ptrs);
        }

        std::vector<char> ofound;
        std::vector<gsPartitionKey> okey;
        selOracle(el, mode, tg, tp, ofound, okey);

        bool warned = false;
        for (index_t r = 0; r < P; ++r)
        {
            const gsSel & s = sels[(size_t)r];
            bool ok = s.finished() && s.numTargets() == (index_t)tg.size() && s.requestsIssued() == steps
                   && s.requestsIssued() <= (index_t)gsSel::maxPasses;
            maxIssued = std::max(maxIssued, s.requestsIssued());
            for (size_t t = 0; t < tg.size(); ++t)
            {
                ++compared;
                const bool match = (s.found((index_t)t) == (ofound[t] != 0))
                                && (!ofound[t] || s.key((index_t)t) == okey[t]);
                if (!match)
                {
                    ok = false;
                    if (!warned)
                        gsWarn << "radix select " << what << ": P=" << P << " N=" << N << " rank " << r
                               << " target " << t << " (group " << tg[t] << ", param " << tp[t]
                               << ") disagrees with the sorted oracle\n";
                    warned = true;
                }
            }
            if (!ok) ++bad;
        }
    }

    static void selectSweep(gsSel::Mode mode, index_t & runs, index_t & compared, index_t & bad,
                            index_t & maxIssued)
    {
        std::mt19937_64 rng(mode == gsSel::Count ? 4242 : 2424);
        const index_t Ps[] = { 1, 2, 3, 7, 16 };
        const index_t Ns[] = { 0, 1, 5, 40, 2000 };
        const int nW = (mode == gsSel::Weighted) ? 3 : 1;
        for (size_t ip = 0; ip < sizeof(Ps) / sizeof(Ps[0]); ++ip)
            for (size_t in = 0; in < sizeof(Ns) / sizeof(Ns[0]); ++in)
                for (int xkind = 0; xkind < 4; ++xkind)
                    for (int multi = 0; multi < 2; ++multi)
                        for (int wkind = 0; wkind < nW; ++wkind)
                            for (int skew = 0; skew < 2; ++skew)
                            {
                                const index_t P = Ps[ip], N = Ns[in];
                                const std::vector<gsSelElem> el = selMakeElems(rng, N, xkind, multi != 0, wkind);
                                std::vector<index_t> tg;
                                std::vector<int64_t> tp;
                                selMakeTargets(rng, el, mode, multi != 0, P, tg, tp);
                                const bool nullGroup = !multi && (rng() & 1);
                                const bool vecOverload = (rng() & 1) != 0;
                                char what[96];
                                std::snprintf(what, sizeof(what), "xkind %d multi %d wkind %d skew %d",
                                              xkind, multi, wkind, skew);
                                selectCase(what, el, P, skew != 0, nullGroup, mode, tg, tp, vecOverload,
                                           compared, bad, maxIssued);
                                ++runs;
                            }
    }

    TEST(RadixSelectCount)
    {
        index_t runs = 0, compared = 0, bad = 0, maxIssued = 0;
        selectSweep(gsSel::Count, runs, compared, bad, maxIssued);
        gsInfo << "gsDistributedPartition RadixSelectCount: " << runs << " runs, " << compared
               << " (rank, target) cases compared, " << bad << " mismatching ranks, max requestsIssued "
               << maxIssued << " (maxPasses " << (unsigned)gsSel::maxPasses << ")\n";
        CHECK(compared > 0);
        CHECK_EQUAL(0, bad);
        CHECK(maxIssued <= (index_t)gsSel::maxPasses);
    }

    TEST(RadixSelectWeighted)
    {
        index_t runs = 0, compared = 0, bad = 0, maxIssued = 0;
        selectSweep(gsSel::Weighted, runs, compared, bad, maxIssued);
        gsInfo << "gsDistributedPartition RadixSelectWeighted: " << runs << " runs, " << compared
               << " (rank, target) cases compared, " << bad << " mismatching ranks, max requestsIssued "
               << maxIssued << " (maxPasses " << (unsigned)gsSel::maxPasses << ")\n";
        CHECK(compared > 0);
        CHECK_EQUAL(0, bad);
        CHECK(maxIssued <= (index_t)gsSel::maxPasses);
    }

    // ------------------------------------------------------------------
    // Protocol enforcement of runLockstep
    // ------------------------------------------------------------------

    // Issues one collective of the given kind and length, then DONE. The first
    // request carries both 4-entry buffers, so that a driver without the check
    // under test misbehaves cleanly instead of reading out of bounds.
    struct gsToyMachine
    {
        typedef double Scalar;
        enum Mode { Sum, Min, DoneAtOnce };

        gsToyMachine(Mode mode, int len, int fill)
        : i64(4, fill), t(4, (double)fill), m_mode(mode), m_len(len), m_started(false)
        { }

        gsCollectiveRequest<double> step()
        {
            gsCollectiveRequest<double> r;
            if (!m_started)
            {
                m_started = true;
                if (m_mode == Sum)      r = gsCollectiveRequest<double>::sum(i64.data(), m_len);
                else if (m_mode == Min) r = gsCollectiveRequest<double>::min(t.data(), m_len);
                r.i64 = i64.data();
                r.t   = t.data();
                r.len = m_len;
            }
            return r;
        }

        std::vector<int64_t> i64;
        std::vector<double>  t;

    private:
        Mode m_mode;
        int  m_len;
        bool m_started;
    };

    // Equal lengths and kinds: one step, and the sum 1 + 2 + ... + P in every entry.
    static index_t toyControl(int P, gsToyMachine::Mode mode)
    {
        std::vector<gsToyMachine> ms;
        ms.reserve((size_t)P);
        for (int r = 0; r < P; ++r)
            ms.push_back(gsToyMachine(mode, 4, r + 1));
        const index_t steps = runLockstep(ms);
        index_t bad = (steps == 1) ? 0 : 1;
        for (int r = 0; r < P; ++r)
            for (int i = 0; i < 4; ++i)
            {
                if (mode == gsToyMachine::Sum && ms[(size_t)r].i64[(size_t)i] != P * (P + 1) / 2) ++bad;
                if (mode == gsToyMachine::Min && ms[(size_t)r].t[(size_t)i] != 1.0) ++bad;
            }
        return bad;
    }

    TEST(LockstepRejectsLengthMismatch)
    {
        index_t cases = 0;
        for (int P = 2; P <= 3; ++P)
        {
            std::vector<gsToyMachine> ms;
            ms.reserve((size_t)P);
            for (int r = 0; r < P; ++r)
                ms.push_back(gsToyMachine(gsToyMachine::Sum, r == 0 ? 3 : 4, r + 1));
            CHECK_THROW(runLockstep(ms), std::exception);

            std::vector<gsToyMachine> ml;
            ml.reserve((size_t)P);
            for (int r = 0; r < P; ++r)
                ml.push_back(gsToyMachine(gsToyMachine::Min, r == P - 1 ? 2 : 4, r + 1));
            CHECK_THROW(runLockstep(ml), std::exception);

            CHECK_EQUAL(0, toyControl(P, gsToyMachine::Sum));
            CHECK_EQUAL(0, toyControl(P, gsToyMachine::Min));
            cases += 4;
        }
        gsInfo << "gsDistributedPartition LockstepRejectsLengthMismatch: " << cases << " cases (throws and controls)\n";
        CHECK(cases > 0);
    }

    TEST(LockstepRejectsKindMismatch)
    {
        index_t cases = 0;
        for (int P = 2; P <= 3; ++P)
        {
            std::vector<gsToyMachine> ms;
            ms.reserve((size_t)P);
            for (int r = 0; r < P; ++r)
                ms.push_back(gsToyMachine(r == 0 ? gsToyMachine::Min : gsToyMachine::Sum, 3, r + 1));
            CHECK_THROW(runLockstep(ms), std::exception);

            std::vector<gsToyMachine> mm;
            mm.reserve((size_t)P);
            for (int r = 0; r < P; ++r)
                mm.push_back(gsToyMachine(r == P - 1 ? gsToyMachine::Min : gsToyMachine::Sum, 3, r + 1));
            CHECK_THROW(runLockstep(mm), std::exception);

            CHECK_EQUAL(0, toyControl(P, gsToyMachine::Sum));
            cases += 3;
        }
        gsInfo << "gsDistributedPartition LockstepRejectsKindMismatch: " << cases << " cases (throws and controls)\n";
        CHECK(cases > 0);
    }

    TEST(LockstepRejectsEarlyDone)
    {
        index_t cases = 0;
        for (int P = 2; P <= 3; ++P)
        {
            std::vector<gsToyMachine> ms;
            ms.reserve((size_t)P);
            for (int r = 0; r < P; ++r)
                ms.push_back(gsToyMachine(r == P - 1 ? gsToyMachine::DoneAtOnce : gsToyMachine::Sum, 4, r + 1));
            CHECK_THROW(runLockstep(ms), std::exception);

            std::vector<gsToyMachine> mf;
            mf.reserve((size_t)P);
            for (int r = 0; r < P; ++r)
                mf.push_back(gsToyMachine(r == 0 ? gsToyMachine::DoneAtOnce : gsToyMachine::Sum, 4, r + 1));
            CHECK_THROW(runLockstep(mf), std::exception);

            CHECK_EQUAL(0, toyControl(P, gsToyMachine::Sum));
            cases += 3;
        }
        gsInfo << "gsDistributedPartition LockstepRejectsEarlyDone: " << cases << " cases (throws and controls)\n";
        CHECK(cases > 0);
    }

    // ------------------------------------------------------------------
    // Distributed RCB == rcbSplit
    // ------------------------------------------------------------------

    // Re-exports the protected static RCB building blocks; never instantiated.
    struct gsRcbProbe : public gsGeometricPartitioner<real_t>
    {
        using gsGeometricPartitioner<real_t>::rcbBisect;
        using gsGeometricPartitioner<real_t>::rcbSplit;
    };

    // Fixed-seed 64-bit LCG (deterministic across platforms and stdlibs).
    static std::uint64_t rcbLcgNext(std::uint64_t & s)
    {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        return s >> 33;
    }

    // Centroid set + weights; all coordinates are exact binary fractions.
    struct gsRcbData
    {
        std::string name;
        bool weighted;
        gsMatrix<real_t> C;
        std::vector<index_t> w;
        bool spike;
        index_t N() const { return static_cast<index_t>(w.size()); }
    };

    // 48 x 48 grid: coordinate ties, a tied root axis choice, ~25% zero weights.
    static gsRcbData rcbGridData(bool weighted)
    {
        gsRcbData d;
        d.name = weighted ? "grid weighted" : "grid unweighted";
        d.weighted = weighted;
        d.spike = false;
        const index_t n = 48;
        d.C.resize(2, n * n);
        d.w.assign(n * n, 0);
        for (index_t j = 0; j < n; ++j)
            for (index_t i = 0; i < n; ++i)
            {
                const index_t id = i + n * j;
                d.C(0, id) = static_cast<real_t>(i);
                d.C(1, id) = static_cast<real_t>(j);
                d.w[id]    = (3 * i + 5 * j) % 4;
            }
        return d;
    }

    // Spike element index in the scatter data.
    static const index_t rcbSpikeId = 1000;

    // Clustered scatter on a dyadic lattice with exact duplicate points and a
    // unique minimum-x element rcbSpikeId. With `spike` all the weight sits
    // on that element.
    static gsRcbData rcbScatterData(bool weighted, bool spike)
    {
        gsRcbData d;
        d.name = spike ? "spike weighted"
               : (weighted ? "scatter weighted" : "scatter unweighted");
        d.weighted = weighted;
        d.spike = spike;
        const index_t N = 3000;
        d.C.resize(2, N);
        d.w.assign(N, 0);
        std::uint64_t seed = 12345;
        for (index_t id = 0; id < N; ++id)
        {
            const std::uint64_t a = rcbLcgNext(seed) % 256;
            const std::uint64_t b = rcbLcgNext(seed) % 256;
            d.C(0, id) = 0.5 + 1.5 * static_cast<real_t>(a * a) / 65536.0;
            d.C(1, id) = static_cast<real_t>(b) / 256.0;
            d.w[id]    = static_cast<index_t>(rcbLcgNext(seed) % 5);
            if (id % 7 == 6)
            {
                d.C(0, id) = d.C(0, id - 1);
                d.C(1, id) = d.C(1, id - 1);
            }
        }
        d.C(0, rcbSpikeId) = 0.0;
        d.C(1, rcbSpikeId) = 0.5;
        if (spike)
        {
            d.w.assign(N, 0);
            d.w[rcbSpikeId] = 1000;
        }
        return d;
    }

    // Weighted data whose weights are all zero: every node has target 0.
    static gsRcbData rcbZeroWeightData(const gsRcbData & base)
    {
        gsRcbData d = base;
        d.name = base.name + " zero weights";
        d.weighted = true;
        d.w.assign(d.N(), 0);
        return d;
    }

    // Ten points on a line with total weight 1 on the last one: the root target
    // 1*kL/k = 0 for P = 3.
    static gsRcbData rcbLineFixture()
    {
        gsRcbData d;
        d.name = "line target 0";
        d.weighted = true;
        d.spike = false;
        const index_t N = 10;
        d.C.resize(1, N);
        d.w.assign(N, 0);
        for (index_t e = 0; e < N; ++e) d.C(0, e) = static_cast<real_t>(e);
        d.w[N - 1] = 1;
        return d;
    }

    // Scatter coordinates with sparse 0/1 weights (about 1/16 of the elements weigh 1).
    static gsRcbData rcbSparseWeightData()
    {
        gsRcbData d = rcbScatterData(true, false);
        d.name = "scatter sparse 0/1 weights";
        std::uint64_t s = 777;
        for (index_t e = 0; e < d.N(); ++e)
            d.w[e] = (rcbLcgNext(s) % 16 == 0) ? 1 : 0;
        return d;
    }

    // N elements on `points` distinct points (point e % points), interleaved ids.
    static gsRcbData rcbTiedPointsData(index_t N, index_t points, bool weighted)
    {
        gsRcbData d;
        d.name = (points == 1 ? "one point" : "five points") + std::string(weighted ? " weighted" : " unweighted");
        d.weighted = weighted;
        d.spike = false;
        d.C.resize(2, N);
        d.w.assign(N, 0);
        for (index_t e = 0; e < N; ++e)
        {
            const index_t p = e % points;
            d.C(0, e) = static_cast<real_t>(p) + (points == 1 ? 2.0 : 0.0);
            d.C(1, e) = 0.5 * static_cast<real_t>((2 * p) % 5) + (points == 1 ? 3.0 : 0.0);
            d.w[e]    = e % 3;
        }
        return d;
    }

    // 2-D, N = 64: x = -0.0 / +0.0 alternating for e < 40, then 1 + e; y = 0.25 (e % 4).
    static gsRcbData rcbSignedZeroData(bool weighted)
    {
        gsRcbData d;
        d.name = weighted ? "signed zeros weighted" : "signed zeros unweighted";
        d.weighted = weighted;
        d.spike = false;
        const index_t N = 64;
        d.C.resize(2, N);
        d.w.assign(N, 0);
        for (index_t e = 0; e < N; ++e)
        {
            d.C(0, e) = (e < 40) ? ((e % 2) ? -0.0 : +0.0) : 1.0 + static_cast<real_t>(e);
            d.C(1, e) = 0.25 * static_cast<real_t>(e % 4);
            d.w[e]    = e % 3;
        }
        return d;
    }

    static gsRcbData rcbSingleElementData(bool weighted)
    {
        gsRcbData d;
        d.name = weighted ? "one element weighted" : "one element unweighted";
        d.weighted = weighted;
        d.spike = false;
        d.C.resize(2, 1);
        d.C(0, 0) = 3.0;
        d.C(1, 0) = 4.0;
        d.w.assign(1, 1);
        return d;
    }

    // N = 5 with a long y extent and a short x extent far from the origin.
    static gsRcbData rcbFiveElementData(bool weighted)
    {
        gsRcbData d;
        d.name = weighted ? "five elements weighted" : "five elements unweighted";
        d.weighted = weighted;
        d.spike = false;
        const real_t x[5] = { 100.0, 100.25, 100.5, 100.75, 101.0 };
        const real_t y[5] = { 0.0, 5.0, 1.0, 4.0, 2.0 };
        const index_t w[5] = { 1, 2, 0, 3, 1 };
        d.C.resize(2, 5);
        d.w.assign(5, 0);
        for (index_t e = 0; e < 5; ++e)
        {
            d.C(0, e) = x[e];
            d.C(1, e) = y[e];
            d.w[e]    = w[e];
        }
        return d;
    }

    static const index_t rcbParts[] = { 1, 2, 3, 5, 7, 16, 64 };
    static const size_t  rcbNumParts = sizeof(rcbParts) / sizeof(rcbParts[0]);

    // Distributes d over P virtual ranks (slices N*r/P), runs the distributed RCB
    // in lockstep and compares the own element set of every rank with the
    // elements rcbSplit labels with that rank. `bad` counts mismatching ranks.
    static void rcbRunCase(const gsRcbData & d, index_t P, index_t & cases, index_t & bad,
                           std::vector<index_t> * refOut = NULL)
    {
        const index_t N = d.N();
        const short_t gd = static_cast<short_t>(d.C.rows());

        std::vector<index_t> ref(N, -1), order(N);
        for (index_t i = 0; i < N; ++i) order[i] = i;
        gsRcbProbe::rcbSplit(d.C, d.w, d.weighted, ref, order, 0, N, P, 0);
        ++cases;

        std::vector<gsMatrix<real_t> >  Cs((size_t)P);
        std::vector<std::vector<index_t> > Ws((size_t)P);
        std::vector<index_t> lo((size_t)P), hi((size_t)P);
        for (index_t r = 0; r < P; ++r)
        {
            lo[(size_t)r] = static_cast<index_t>(static_cast<std::int64_t>(N) * r / P);
            hi[(size_t)r] = static_cast<index_t>(static_cast<std::int64_t>(N) * (r + 1) / P);
            if (hi[(size_t)r] > lo[(size_t)r])
            {
                Cs[(size_t)r] = d.C.middleCols(lo[(size_t)r], hi[(size_t)r] - lo[(size_t)r]);
                Ws[(size_t)r].assign(d.w.begin() + lo[(size_t)r], d.w.begin() + hi[(size_t)r]);
            }
        }

        std::vector<gsDistributedRcb<real_t> > ms;
        ms.reserve((size_t)P);
        for (index_t r = 0; r < P; ++r)
        {
            const index_t * wp = (d.weighted && !Ws[(size_t)r].empty()) ? Ws[(size_t)r].data() : NULL;
            ms.push_back(gsDistributedRcb<real_t>(Cs[(size_t)r], gd, lo[(size_t)r], NULL, wp,
                                                  d.weighted, P, r, N));
        }
        std::vector<gsDistributedRcb<real_t> *> ptrs((size_t)P);
        for (index_t r = 0; r < P; ++r) ptrs[(size_t)r] = &ms[(size_t)r];
        runLockstep(ptrs);

        std::vector<std::vector<index_t> > expect((size_t)P);
        bool refOk = true;
        for (index_t e = 0; e < N; ++e)
        {
            if (ref[e] < 0 || ref[e] >= P) { refOk = false; continue; }
            expect[(size_t)ref[e]].push_back(e);
        }
        if (!refOk)
            gsWarn << "rcb " << d.name << " P=" << P << ": rcbSplit left labels outside [0,P)\n";

        index_t levelBound = 0;
        while (((index_t)1 << levelBound) < P) ++levelBound;

        bool warned = false;
        std::int64_t total = 0;
        for (index_t r = 0; r < P; ++r)
        {
            const gsDistributedRcb<real_t> & m = ms[(size_t)r];
            const std::vector<index_t> & own = m.ownElements();
            total += (std::int64_t)own.size();

            bool ok = refOk && m.finished() && own == expect[(size_t)r]
                   && std::adjacent_find(own.begin(), own.end(), std::greater_equal<index_t>()) == own.end();
            const std::vector<index_t> & lab = m.localLabels();
            if (lab.size() != (size_t)(hi[(size_t)r] - lo[(size_t)r]))
                ok = false;
            else
                for (size_t i = 0; i < lab.size(); ++i)
                    if (lab[i] != ref[lo[(size_t)r] + (index_t)i]) ok = false;
            if (r == 0 && m.levels() > levelBound)
                ok = false;
            if (!ok)
            {
                if (!warned)
                    gsWarn << "rcb " << d.name << ": P=" << P << " rank " << r << " owns " << own.size()
                           << " elements, rcbSplit labels " << expect[(size_t)r].size()
                           << " with that rank (or levels/labels disagree)\n";
                warned = true;
                ++bad;
            }
        }
        if (total != N)
        {
            gsWarn << "rcb " << d.name << ": P=" << P << " own counts sum to " << total << ", not " << N << "\n";
            ++bad;
        }
        if (refOut) *refOut = ref;
    }

    static void rcbSweep(const std::vector<gsRcbData> & sets, index_t & cases, index_t & bad)
    {
        for (size_t c = 0; c < sets.size(); ++c)
            for (size_t ip = 0; ip < rcbNumParts; ++ip)
                rcbRunCase(sets[c], rcbParts[ip], cases, bad);
    }

    static void rcbReport(const char * test, index_t cases, index_t bad)
    {
        gsInfo << "gsDistributedPartition " << test << ": " << cases
               << " (data set, P) cases, " << bad << " mismatches\n";
    }

    TEST(DistributedRcbMatchesRcbSplit)
    {
        std::vector<gsRcbData> sets;
        sets.push_back(rcbGridData(false));
        sets.push_back(rcbGridData(true));
        sets.push_back(rcbScatterData(false, false));
        sets.push_back(rcbScatterData(true, false));
        sets.push_back(rcbScatterData(true, true));
        index_t cases = 0, bad = 0;
        rcbSweep(sets, cases, bad);
        rcbReport("DistributedRcbMatchesRcbSplit", cases, bad);
        CHECK(cases > 0);
        CHECK_EQUAL(0, bad);
    }

    TEST(DistributedRcbZeroTargetMatchesRcbSplit)
    {
        std::vector<gsRcbData> sets;
        sets.push_back(rcbZeroWeightData(rcbGridData(false)));
        sets.push_back(rcbZeroWeightData(rcbScatterData(false, false)));
        sets.push_back(rcbSparseWeightData());
        index_t cases = 0, bad = 0;
        rcbSweep(sets, cases, bad);

        // N = 10 on a line, P = 3: root target 0, labels {0,0,0,0,1,1,1,2,2,2}.
        const gsRcbData line = rcbLineFixture();
        std::vector<index_t> ref;
        rcbRunCase(line, 3, cases, bad, &ref);
        const index_t expected[10] = { 0, 0, 0, 0, 1, 1, 1, 2, 2, 2 };
        CHECK_EQUAL(10, (index_t)ref.size());
        for (size_t e = 0; e < ref.size() && e < 10; ++e)
            CHECK_EQUAL(expected[e], ref[e]);

        rcbReport("DistributedRcbZeroTargetMatchesRcbSplit", cases, bad);
        CHECK(cases > 0);
        CHECK_EQUAL(0, bad);
    }

    TEST(DistributedRcbTiesAndSignedZeros)
    {
        std::vector<gsRcbData> sets;
        sets.push_back(rcbTiedPointsData(200, 5, false));
        sets.push_back(rcbTiedPointsData(200, 5, true));
        sets.push_back(rcbTiedPointsData(100, 1, false));
        sets.push_back(rcbTiedPointsData(100, 1, true));
        index_t cases = 0, bad = 0;
        rcbSweep(sets, cases, bad);

        // A cut that falls inside the run of +-0 coordinates is decided by the id
        // tiebreak, which is only exercised if the run carries two labels.
        index_t cutInsideZeroRun = 0;
        for (int weighted = 0; weighted < 2; ++weighted)
        {
            const gsRcbData z = rcbSignedZeroData(weighted != 0);
            for (size_t ip = 0; ip < rcbNumParts; ++ip)
            {
                std::vector<index_t> ref;
                rcbRunCase(z, rcbParts[ip], cases, bad, &ref);
                std::set<index_t> labels(ref.begin(), ref.begin() + 40);
                if (labels.size() >= 2) ++cutInsideZeroRun;
            }
        }
        CHECK(cutInsideZeroRun > 0);

        rcbReport("DistributedRcbTiesAndSignedZeros", cases, bad);
        gsInfo << "gsDistributedPartition DistributedRcbTiesAndSignedZeros: " << cutInsideZeroRun
               << " (data set, P) cases with a cut inside the +-0 run\n";
        CHECK(cases > 0);
        CHECK_EQUAL(0, bad);
    }

    TEST(DistributedRcbTinyN)
    {
        std::vector<gsRcbData> sets;
        sets.push_back(rcbSingleElementData(false));
        sets.push_back(rcbSingleElementData(true));
        sets.push_back(rcbFiveElementData(false));
        sets.push_back(rcbFiveElementData(true));
        index_t cases = 0, bad = 0;
        rcbSweep(sets, cases, bad);

        // The y extent (5) exceeds the x extent (1): the root axis is y.
        const gsRcbData five = rcbFiveElementData(false);
        const real_t xExt = five.C.row(0).maxCoeff() - five.C.row(0).minCoeff();
        const real_t yExt = five.C.row(1).maxCoeff() - five.C.row(1).minCoeff();
        CHECK(yExt > xExt);

        rcbReport("DistributedRcbTinyN", cases, bad);
        CHECK(cases > 0);
        CHECK_EQUAL(0, bad);
    }

    // ------------------------------------------------------------------
    // runCollectives on a size-1 communicator
    // ------------------------------------------------------------------

    TEST(RunCollectivesSizeOneCommGivesAllIds)
    {
        gsMpi::init();
        const gsMpiComm comm(gsMpi::localComm());
        CHECK_EQUAL(1, comm.size());

        std::vector<gsRcbData> sets;
        sets.push_back(rcbGridData(true));
        sets.push_back(rcbScatterData(false, false));

        index_t cases = 0, bad = 0;
        for (size_t c = 0; c < sets.size(); ++c)
        {
            const gsRcbData & d = sets[c];
            const index_t N = d.N();
            const short_t gd = static_cast<short_t>(d.C.rows());
            const index_t * wp = d.weighted ? d.w.data() : NULL;

            gsDistributedRcb<real_t> viaComm(d.C, gd, 0, NULL, wp, d.weighted, 1, 0, N);
            runCollectives(viaComm, comm);

            gsDistributedRcb<real_t> viaLockstep(d.C, gd, 0, NULL, wp, d.weighted, 1, 0, N);
            std::vector<gsDistributedRcb<real_t> *> ptrs(1, &viaLockstep);
            runLockstep(ptrs);

            std::vector<index_t> all(N);
            for (index_t i = 0; i < N; ++i) all[i] = i;

            ++cases;
            if (!(viaComm.finished() && viaComm.ownElements() == all))
            {
                gsWarn << "runCollectives over a size-1 comm: " << d.name << " does not own all ids\n";
                ++bad;
            }
            if (!(viaComm.ownElements() == viaLockstep.ownElements()))
            {
                gsWarn << "runCollectives over a size-1 comm: " << d.name << " differs from runLockstep\n";
                ++bad;
            }
        }
        gsInfo << "gsDistributedPartition RunCollectivesSizeOneCommGivesAllIds: " << cases
               << " data sets, " << bad << " mismatches\n";
        CHECK(cases > 0);
        CHECK_EQUAL(0, bad);
    }

    // ------------------------------------------------------------------
    // Distributed curve machine == curveLabels
    // ------------------------------------------------------------------

    // Re-exports the protected static curve cut; never instantiated.
    struct gsDistCurveProbe : public gsGeometricPartitioner<real_t>
    {
        using gsGeometricPartitioner<real_t>::curveLabels;
    };

    typedef gsGeometricPartitioner<real_t> gsCurveGeo;

    // Non-uniform multi-patch geometry with its basis and DOF mapper. The
    // control-net box and the centroid box differ.
    struct gsCurveFixture
    {
        gsMultiPatch<real_t>         mp;
        gsMultiBasis<real_t>         mb;
        gsBoundaryConditions<real_t> bc;
        gsDofMapper                  mapper;

        gsCurveFixture(const gsMultiPatch<real_t> & geo, int refinements)
        : mp(geo), mb(mp)
        {
            for (int i = 0; i < refinements; ++i) mb.uniformRefine();
            mapper = createMapper(mb, bc, /*nComp=*/1, /*unk=*/0, /*conforming=*/true, /*finalize=*/true);
        }

        index_t numElements() const
        { return static_cast<index_t>(mb.domain()->numElements()); }
    };

    // Quadratic square whose centre control point is pulled to (2.5, 1.7),
    // plus a second patch: a bilinear square shifted to x in [3.1, 4.1] when
    // coincident is false, a copy of the first patch otherwise.
    static gsMultiPatch<real_t> curveSquares(bool coincident)
    {
        gsMultiPatch<real_t> mp;
        gsNurbsCreator<real_t>::TensorBSpline2Ptr a = gsNurbsCreator<real_t>::BSplineSquareDeg(2);
        a->coefs()(4, 0) = 2.5;
        a->coefs()(4, 1) = 1.7;
        mp.addPatch(*a);
        if (coincident)
            mp.addPatch(*a);
        else
            mp.addPatch(*gsNurbsCreator<real_t>::BSplineSquare(1.0, 3.1, 0.2));
        return mp;
    }

    // Quadratic cube whose centre control point (row 13) is pulled to
    // (2.2, 1.6, 2.9), plus an unpulled quadratic cube shifted by (3.3, 0.4, 0.7).
    static gsMultiPatch<real_t> curveCubes()
    {
        gsMultiPatch<real_t> mp;
        gsNurbsCreator<real_t>::TensorBSpline3Ptr a = gsNurbsCreator<real_t>::BSplineCube(2);
        a->coefs()(13, 0) = 2.2;
        a->coefs()(13, 1) = 1.6;
        a->coefs()(13, 2) = 2.9;
        mp.addPatch(*a);
        gsNurbsCreator<real_t>::TensorBSpline3Ptr b = gsNurbsCreator<real_t>::BSplineCube(2);
        b->coefs().col(0).array() += 3.3;
        b->coefs().col(1).array() += 0.4;
        b->coefs().col(2).array() += 0.7;
        mp.addPatch(*b);
        return mp;
    }

    // Keys of the columns of C, computed exactly as labelsByCurve does: one
    // gsSpaceFillingCurve on `box` (geoDim x 2, lower and upper corner).
    static std::vector<std::uint64_t> curveKeys(const gsMatrix<real_t> & C, const gsMatrix<real_t> & box,
                                                gsSpaceFillingCurve::Curve curve)
    {
        const gsSpaceFillingCurve sfc(box, curve);
        std::vector<std::uint64_t> keys(static_cast<size_t>(C.cols()));
        gsVector<real_t> pt(C.rows());
        for (index_t e = 0; e < C.cols(); ++e)
        {
            for (index_t d = 0; d < C.rows(); ++d) pt[d] = C(d, e);
            keys[(size_t)e] = sfc.encode(pt);
        }
        return keys;
    }

    // geoDim x 2 (lower, upper) row-wise extent of the centroids.
    static gsMatrix<real_t> centroidBox(const gsMatrix<real_t> & C)
    {
        gsMatrix<real_t> box(C.rows(), 2);
        box.col(0) = C.rowwise().minCoeff();
        box.col(1) = C.rowwise().maxCoeff();
        return box;
    }

    // Centroids (geoDim x N) and element weights from a serial partitioner.
    static gsMatrix<real_t> curveCentroids(gsCurveFixture & fx, gsCurveGeo::Strategy strat, index_t nparts,
                                           std::vector<index_t> * weights = NULL,
                                           std::vector<index_t> * labels = NULL)
    {
        gsCurveGeo::Options opts;
        opts.strategy      = strat;
        opts.keepCentroids = true;
        gsCurveGeo geo(fx.mp, fx.mb, fx.mapper, nparts, opts);
        geo.partition();
        if (weights) *weights = geo.elementWeights();
        if (labels)  *labels  = geo.labels();
        return geo.centroids();
    }

    static gsSpaceFillingCurve::Curve curveOf(gsCurveGeo::Strategy s)
    {
        return s == gsCurveGeo::hilbert ? gsSpaceFillingCurve::Hilbert : gsSpaceFillingCurve::Morton;
    }

    // Net-box keys of a fixture for one curve.
    static std::vector<std::uint64_t> fixtureKeys(gsCurveFixture & fx, gsSpaceFillingCurve::Curve curve)
    {
        gsMatrix<real_t> box;
        fx.mp.boundingBox(box);
        const gsCurveGeo::Strategy s = curve == gsSpaceFillingCurve::Hilbert ? gsCurveGeo::hilbert
                                                                              : gsCurveGeo::morton;
        return curveKeys(curveCentroids(fx, s, 2), box, curve);
    }

    // Result of one oracle-vs-machine comparison.
    struct gsCurveRun
    {
        std::vector<index_t>               labels; // [N] oracle label of every element
        std::vector<std::vector<index_t> > own;    // [P] own ids of every virtual rank
    };

    // Oracle labels of (keys, w) for P parts. O(N log N).
    static std::vector<index_t> curveOracle(const std::vector<std::uint64_t> & keys,
                                            const std::vector<index_t> & w, index_t P)
    {
        std::vector<index_t> labels;
        gsDistCurveProbe::curveLabels(keys, w, P, labels);
        return labels;
    }

    // Runs one gsDistributedCurve per virtual rank in lockstep on the slices
    // [N r / P, N (r+1) / P) of (keys, w); returns the own ids of every rank.
    static std::vector<std::vector<index_t> > runDistCurve(const std::vector<std::uint64_t> & keys,
                                                            const std::vector<index_t> & w, index_t P)
    {
        const index_t N = static_cast<index_t>(keys.size());
        std::vector<gsDistributedCurve<real_t> > ms;
        ms.reserve((size_t)P);
        for (index_t r = 0; r < P; ++r)
        {
            const index_t lo = static_cast<index_t>(static_cast<std::int64_t>(N) * r / P);
            const index_t hi = static_cast<index_t>(static_cast<std::int64_t>(N) * (r + 1) / P);
            const index_t n  = hi - lo;
            ms.push_back(gsDistributedCurve<real_t>(n, lo, n ? &keys[(size_t)lo] : NULL,
                                                    n ? &w[(size_t)lo] : NULL, P, r, N));
        }
        std::vector<gsDistributedCurve<real_t> *> ptrs((size_t)P);
        for (index_t r = 0; r < P; ++r) ptrs[(size_t)r] = &ms[(size_t)r];
        runLockstep(ptrs);

        std::vector<std::vector<index_t> > own((size_t)P);
        for (index_t r = 0; r < P; ++r)
            if (ms[(size_t)r].finished())
                own[(size_t)r] = ms[(size_t)r].ownElements();
        return own;
    }

    // Number of ranks whose own ids differ from {e : labels[e] == r}, plus one
    // if the own ids are not a disjoint cover of [0, N). Warns with `tag`.
    static index_t checkOwnVsLabels(const char * tag, const std::vector<std::vector<index_t> > & own,
                                    const std::vector<index_t> & labels, index_t P)
    {
        const index_t N = static_cast<index_t>(labels.size());
        index_t bad = 0;
        std::vector<std::vector<index_t> > expect((size_t)P);
        for (index_t e = 0; e < N; ++e)
            expect[(size_t)labels[(size_t)e]].push_back(e);

        std::vector<char> seen((size_t)N, 0);
        bool cover = (own.size() == (size_t)P);
        for (index_t r = 0; r < P && own.size() == (size_t)P; ++r)
        {
            const std::vector<index_t> & o = own[(size_t)r];
            if (o != expect[(size_t)r])
            {
                gsWarn << "distCurve " << tag << " P=" << P << " rank " << r << " owns " << o.size()
                       << " ids, the oracle labels " << expect[(size_t)r].size() << " with it\n";
                ++bad;
            }
            for (size_t i = 0; i < o.size(); ++i)
            {
                if (o[i] < 0 || o[i] >= N || seen[(size_t)o[i]]) cover = false;
                else seen[(size_t)o[i]] = 1;
            }
        }
        if (std::find(seen.begin(), seen.end(), 0) != seen.end()) cover = false;
        if (!cover)
        {
            gsWarn << "distCurve " << tag << " P=" << P << ": the own ids are not a disjoint cover of [0," << N << ")\n";
            ++bad;
        }
        return bad;
    }

    // Oracle plus machine for (keys, w, P); `bad` accumulates checkOwnVsLabels.
    static gsCurveRun curveCase(const char * tag, const std::vector<std::uint64_t> & keys,
                                const std::vector<index_t> & w, index_t P, index_t & cases, index_t & bad)
    {
        gsCurveRun run;
        run.labels = curveOracle(keys, w, P);
        run.own    = runDistCurve(keys, w, P);
        ++cases;
        bad += checkOwnVsLabels(tag, run.own, run.labels, P);
        return run;
    }

    static index_t maxOf(const std::vector<index_t> & v)
    { return *std::max_element(v.begin(), v.end()); }

    // Elements in (key, id) order.
    static std::vector<index_t> curveOrder(const std::vector<std::uint64_t> & keys)
    {
        std::vector<index_t> order(keys.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = static_cast<index_t>(i);
        std::sort(order.begin(), order.end(), [&keys](index_t a, index_t b)
                  { return keys[(size_t)a] < keys[(size_t)b] || (keys[(size_t)a] == keys[(size_t)b] && a < b); });
        return order;
    }

    // Number of parts that receive no element.
    static index_t emptyParts(const std::vector<index_t> & labels, index_t P)
    {
        std::vector<char> used((size_t)P, 0);
        for (size_t e = 0; e < labels.size(); ++e) used[(size_t)labels[e]] = 1;
        return static_cast<index_t>(std::count(used.begin(), used.end(), 0));
    }

    // A key set to sweep: a name, a dimension and the keys.
    struct gsCurveKeySet
    {
        std::string                name;
        std::vector<std::uint64_t> keys;
    };

    static std::vector<gsCurveKeySet> curveKeySets()
    {
        std::vector<gsCurveKeySet> sets;
        gsCurveFixture f2(curveSquares(false), 3), f3(curveCubes(), 2);
        for (int dim = 2; dim <= 3; ++dim)
            for (int c = 0; c < 2; ++c)
            {
                gsCurveKeySet s;
                const gsSpaceFillingCurve::Curve curve = c ? gsSpaceFillingCurve::Morton : gsSpaceFillingCurve::Hilbert;
                s.name = std::string(dim == 2 ? "2D " : "3D ") + (c ? "morton" : "hilbert");
                s.keys = fixtureKeys(dim == 2 ? f2 : f3, curve);
                sets.push_back(s);
            }
        return sets;
    }

    // Net-box and centroid-box checks of one fixture for both curves and P in {2,3,7}.
    // netBad: curveLabels(net keys) != geo.labels(), or machine != geo.labels();
    // differ: combinations where the centroid-box labels differ from geo.labels().
    static void curveNetBoxSweep(const char * name, gsCurveFixture & fx, index_t & cases, index_t & netBad,
                                 index_t & differ, real_t & boxGap)
    {
        gsMatrix<real_t> netBox;
        fx.mp.boundingBox(netBox);
        const gsCurveGeo::Strategy strats[] = { gsCurveGeo::hilbert, gsCurveGeo::morton };
        const index_t ps[] = { 2, 3, 7 };
        for (size_t is = 0; is < 2; ++is)
            for (size_t ip = 0; ip < 3; ++ip)
            {
                std::vector<index_t> w, geoLabels;
                const gsMatrix<real_t> C = curveCentroids(fx, strats[is], ps[ip], &w, &geoLabels);
                const gsMatrix<real_t> cenBox = centroidBox(C);
                boxGap = (std::max)(boxGap, (netBox - cenBox).cwiseAbs().maxCoeff()
                                            / (netBox.col(1) - netBox.col(0)).maxCoeff());

                const std::vector<std::uint64_t> netKeys = curveKeys(C, netBox, curveOf(strats[is]));
                const std::vector<std::uint64_t> cenKeys = curveKeys(C, cenBox, curveOf(strats[is]));
                ++cases;
                if (curveOracle(netKeys, w, ps[ip]) != geoLabels)
                {
                    gsWarn << "distCurve " << name << " P=" << ps[ip] << ": curveLabels(net-box keys) != geo.labels()\n";
                    ++netBad;
                }
                if (curveOracle(cenKeys, w, ps[ip]) != geoLabels) ++differ;
                netBad += checkOwnVsLabels(name, runDistCurve(netKeys, w, ps[ip]), geoLabels, ps[ip]);
            }
    }

    TEST(CurveKeysUseControlNetBox)
    {
        index_t cases = 0;
        for (int dim = 2; dim <= 3; ++dim)
        {
            gsCurveFixture fx(dim == 2 ? curveSquares(false) : curveCubes(), dim == 2 ? 3 : 2);
            index_t c = 0, netBad = 0, differ = 0;
            real_t boxGap = 0;
            curveNetBoxSweep(dim == 2 ? "2D net box" : "3D net box", fx, c, netBad, differ, boxGap);
            cases += c;
            gsInfo << "gsDistributedPartition CurveKeysUseControlNetBox " << dim << "D: N=" << fx.numElements()
                   << ", " << c << " (curve, P) cases, " << netBad << " mismatches, " << differ
                   << " centroid-box label sets differ, relative box gap " << boxGap << "\n";
            CHECK(fx.numElements() >= 64);
            CHECK_EQUAL(0, netBad);
            CHECK(differ >= 1);
            CHECK(boxGap > 0.05);
        }
        CHECK_EQUAL(12, cases);
    }

    TEST(DistCurveUnitWeights)
    {
        const std::vector<gsCurveKeySet> sets = curveKeySets();
        index_t cases = 0, bad = 0, vacuous = 0;
        for (size_t s = 0; s < sets.size(); ++s)
        {
            const std::vector<index_t> w(sets[s].keys.size(), 1);
            for (size_t ip = 0; ip < rcbNumParts; ++ip)
            {
                const index_t P = rcbParts[ip];
                const gsCurveRun run = curveCase(sets[s].name.c_str(), sets[s].keys, w, P, cases, bad);
                if (P > 1 && maxOf(run.labels) == 0) ++vacuous;
                if (P == 1 && run.own[0].size() != sets[s].keys.size()) ++bad;
            }
        }
        gsInfo << "gsDistributedPartition DistCurveUnitWeights: " << cases << " cases, " << bad << " mismatches\n";
        CHECK_EQUAL(4 * (index_t)rcbNumParts, cases);
        CHECK_EQUAL(0, bad);
        CHECK_EQUAL(0, vacuous);
    }

    TEST(DistCurveSkewedWeights)
    {
        const std::vector<gsCurveKeySet> sets = curveKeySets();
        index_t cases = 0, bad = 0, vacuous = 0;
        std::uint64_t seed = 12345;
        for (size_t s = 0; s < sets.size(); ++s)
        {
            std::vector<index_t> w(sets[s].keys.size());
            for (size_t e = 0; e < w.size(); ++e) w[e] = 1 + static_cast<index_t>(rcbLcgNext(seed) % 9);
            if (*std::max_element(w.begin(), w.end()) <= *std::min_element(w.begin(), w.end())) ++vacuous;
            for (size_t ip = 0; ip < rcbNumParts; ++ip)
            {
                const index_t P = rcbParts[ip];
                const gsCurveRun run = curveCase(sets[s].name.c_str(), sets[s].keys, w, P, cases, bad);
                if (P > 1 && maxOf(run.labels) == 0) ++vacuous;
            }
        }
        gsInfo << "gsDistributedPartition DistCurveSkewedWeights: " << cases << " cases, " << bad << " mismatches\n";
        CHECK_EQUAL(4 * (index_t)rcbNumParts, cases);
        CHECK_EQUAL(0, bad);
        CHECK_EQUAL(0, vacuous);
    }

    TEST(DistCurveAllZeroWeights)
    {
        const std::vector<gsCurveKeySet> sets = curveKeySets();
        // All-zero weights fall back to unit weights: the count split, with
        // part sizes floor(N/P) or ceil(N/P).
        index_t cases = 0, bad = 0, notCount = 0;
        for (size_t s = 0; s < sets.size(); ++s)
        {
            const std::vector<index_t> w(sets[s].keys.size(), 0), ones(sets[s].keys.size(), 1);
            const index_t N = static_cast<index_t>(w.size());
            for (size_t ip = 0; ip < rcbNumParts; ++ip)
            {
                const index_t P = rcbParts[ip];
                const gsCurveRun run = curveCase(sets[s].name.c_str(), sets[s].keys, w, P, cases, bad);
                if (run.labels != curveOracle(sets[s].keys, ones, P)) ++notCount;
                for (index_t r = 0; r < P; ++r)
                {
                    const index_t m = (index_t)run.own[(size_t)r].size();
                    if (m != N / P && m != (N + P - 1) / P) ++notCount;
                }
            }
        }
        gsInfo << "gsDistributedPartition DistCurveAllZeroWeights: " << cases << " cases, " << bad
               << " mismatches, " << notCount << " not count splits\n";
        CHECK_EQUAL(4 * (index_t)rcbNumParts, cases);
        CHECK_EQUAL(0, bad);
        CHECK_EQUAL(0, notCount);
    }

    TEST(DistCurveOneHeavyElement)
    {
        const std::vector<gsCurveKeySet> sets = curveKeySets();
        index_t cases = 0, bad = 0, noGap = 0;
        for (size_t s = 0; s < sets.size(); ++s)
        {
            const index_t N = static_cast<index_t>(sets[s].keys.size());
            std::vector<index_t> w(sets[s].keys.size(), 1);
            w[(size_t)curveOrder(sets[s].keys)[(size_t)(N / 3)]] = 10 * N;
            for (size_t ip = 0; ip < rcbNumParts; ++ip)
            {
                const index_t P = rcbParts[ip];
                const gsCurveRun run = curveCase(sets[s].name.c_str(), sets[s].keys, w, P, cases, bad);
                if (P >= 3 && emptyParts(run.labels, P) == 0) ++noGap;
            }
        }
        gsInfo << "gsDistributedPartition DistCurveOneHeavyElement: " << cases << " cases, " << bad << " mismatches\n";
        CHECK_EQUAL(4 * (index_t)rcbNumParts, cases);
        CHECK_EQUAL(0, bad);
        CHECK_EQUAL(0, noGap);
    }

    TEST(DistCurveZeroWeightRuns)
    {
        const std::vector<gsCurveKeySet> sets = curveKeySets();
        index_t cases = 0, bad = 0, edgeCuts = 0;
        for (size_t s = 0; s < sets.size(); ++s)
        {
            const index_t N = static_cast<index_t>(sets[s].keys.size());
            const std::vector<index_t> order = curveOrder(sets[s].keys);

            // Zero runs [start, start + len) in curve positions: one at position 0, one at the end.
            std::vector<std::pair<index_t, index_t> > runs;
            runs.push_back(std::make_pair((index_t)0, (index_t)3));
            for (index_t p = N / 8; p + 6 < N - 4; p += N / 8) runs.push_back(std::make_pair(p, (index_t)(1 + (p / 7) % 5)));
            runs.push_back(std::make_pair(N - 4, (index_t)4));
            std::vector<index_t> w(sets[s].keys.size(), 1);
            for (size_t k = 0; k < runs.size(); ++k)
                for (index_t i = 0; i < runs[k].second; ++i)
                    w[(size_t)order[(size_t)(runs[k].first + i)]] = 0;

            for (size_t ip = 0; ip < rcbNumParts; ++ip)
            {
                const index_t P = rcbParts[ip];
                const gsCurveRun run = curveCase(sets[s].name.c_str(), sets[s].keys, w, P, cases, bad);
                if (P == 1) continue;
                for (size_t k = 0; k < runs.size(); ++k)
                {
                    if (runs[k].first == 0) continue;
                    const index_t head = run.labels[(size_t)order[(size_t)runs[k].first]];
                    const index_t pred = run.labels[(size_t)order[(size_t)(runs[k].first - 1)]];
                    for (index_t i = 1; i < runs[k].second; ++i)
                        if (run.labels[(size_t)order[(size_t)(runs[k].first + i)]] != head) ++bad;
                    if (pred < head) ++edgeCuts;
                }
            }
        }
        gsInfo << "gsDistributedPartition DistCurveZeroWeightRuns: " << cases << " cases, " << bad
               << " mismatches, " << edgeCuts << " runs entered by a cut\n";
        CHECK_EQUAL(4 * (index_t)rcbNumParts, cases);
        CHECK_EQUAL(0, bad);
        CHECK(edgeCuts >= 1);
    }

    TEST(DistCurveDuplicateKeys)
    {
        std::vector<gsCurveKeySet> sets;
        {
            gsCurveFixture twin(curveSquares(true), 3);
            for (int c = 0; c < 2; ++c)
            {
                gsCurveKeySet s;
                s.name = c ? "twin patches morton" : "twin patches hilbert";
                s.keys = fixtureKeys(twin, c ? gsSpaceFillingCurve::Morton : gsSpaceFillingCurve::Hilbert);
                const size_t half = s.keys.size() / 2;
                for (size_t e = 0; e < half; ++e)
                    CHECK(s.keys[e] == s.keys[e + half]);
                sets.push_back(s);
            }
        }
        {
            gsCurveFixture f2(curveSquares(false), 3);
            gsCurveKeySet s;
            s.name = "key % 7";
            s.keys = fixtureKeys(f2, gsSpaceFillingCurve::Hilbert);
            for (size_t e = 0; e < s.keys.size(); ++e) s.keys[e] %= 7;
            sets.push_back(s);
        }

        index_t cases = 0, bad = 0, splitTies = 0;
        std::uint64_t seed = 777;
        for (size_t s = 0; s < sets.size(); ++s)
        {
            const std::vector<std::uint64_t> & keys = sets[s].keys;
            std::vector<index_t> skewed(keys.size());
            for (size_t e = 0; e < skewed.size(); ++e) skewed[e] = 1 + static_cast<index_t>(rcbLcgNext(seed) % 9);
            const std::vector<index_t> unit(keys.size(), 1);
            for (int kind = 0; kind < 2; ++kind)
                for (size_t ip = 0; ip < rcbNumParts; ++ip)
                {
                    const index_t P = rcbParts[ip];
                    const gsCurveRun run = curveCase(sets[s].name.c_str(), keys, kind ? skewed : unit, P, cases, bad);
                    if (P == 1) continue;
                    std::map<std::uint64_t, std::pair<index_t, index_t> > range; // key -> (min, max) label
                    for (size_t e = 0; e < keys.size(); ++e)
                    {
                        std::map<std::uint64_t, std::pair<index_t, index_t> >::iterator it = range.find(keys[e]);
                        if (it == range.end()) range[keys[e]] = std::make_pair(run.labels[e], run.labels[e]);
                        else
                        {
                            it->second.first  = (std::min)(it->second.first,  run.labels[e]);
                            it->second.second = (std::max)(it->second.second, run.labels[e]);
                        }
                    }
                    for (std::map<std::uint64_t, std::pair<index_t, index_t> >::iterator it = range.begin();
                         it != range.end(); ++it)
                        if (it->second.first != it->second.second) { ++splitTies; break; }
                }
        }
        gsInfo << "gsDistributedPartition DistCurveDuplicateKeys: " << cases << " cases, " << bad
               << " mismatches, " << splitTies << " cases with a tied key split over two parts\n";
        CHECK_EQUAL(3 * 2 * (index_t)rcbNumParts, cases);
        CHECK_EQUAL(0, bad);
        CHECK(splitTies >= 1);
    }

    TEST(DistCurveFewerElementsThanParts)
    {
        gsCurveFixture fx(curveSquares(false), 0);
        const index_t N = fx.numElements();
        CHECK(N >= 2 && N < 5);
        const index_t ps[] = { 5, 7, 16, 64 };
        index_t cases = 0, bad = 0, emptySlices = 0, emptyPartsSeen = 0;
        for (int c = 0; c < 2; ++c)
        {
            const std::vector<std::uint64_t> keys =
                fixtureKeys(fx, c ? gsSpaceFillingCurve::Morton : gsSpaceFillingCurve::Hilbert);
            for (int kind = 0; kind < 2; ++kind)
            {
                const std::vector<index_t> w(keys.size(), kind ? 0 : 1);
                for (size_t ip = 0; ip < 4; ++ip)
                {
                    const index_t P = ps[ip];
                    const gsCurveRun run = curveCase(c ? "tiny morton" : "tiny hilbert", keys, w, P, cases, bad);
                    for (index_t r = 0; r < P; ++r)
                        if (static_cast<std::int64_t>(N) * r / P == static_cast<std::int64_t>(N) * (r + 1) / P)
                            ++emptySlices;
                    emptyPartsSeen += emptyParts(run.labels, P);
                }
            }
        }
        gsInfo << "gsDistributedPartition DistCurveFewerElementsThanParts: N=" << N << ", " << cases << " cases, "
               << bad << " mismatches, " << emptySlices << " empty slices, " << emptyPartsSeen << " empty parts\n";
        CHECK_EQUAL(16, cases);
        CHECK_EQUAL(0, bad);
        CHECK(emptySlices >= 1);
        CHECK(emptyPartsSeen >= 1);
    }
}
