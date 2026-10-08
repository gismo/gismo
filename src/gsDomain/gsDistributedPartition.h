/** @file gsDistributedPartition.h

    @brief Building blocks of the distributed partitioners: a collective
    request protocol, an order-preserving key trait and a multi-target radix
    select.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

#include <cstring>
#include <limits>

/**
   @defgroup gsDistributedPartitionProtocol Collective request protocol

   A distributed algorithm is written as a resumable state machine that emits
   one collective request at a time (kind, length, buffer). A thin driver
   forwards each request to the communicator; a test harness runs P machines
   in lockstep in one process. The machine therefore never touches a
   communicator, and the sequence of collectives of all ranks is identical by
   construction as long as every decision of the machine depends on replicated
   inputs and reduced results only.

   Machine concept (top level):
   @code
   struct M
   {
       typedef <T> Scalar;
       gsCollectiveRequest<Scalar> step();
   };
   @endcode
   - The first step() starts the machine and receives no result. Every later
     call means: the request returned by the previous call has completed and
     its result is in that request's buffers.
   - step() returns a request of kind DONE when the machine has finished; any
     call after DONE returns DONE again.
   - All replicated decisions inside step() use only inputs that are
     identical on all ranks and reduced results.
   - A request points into storage owned by the machine (the data() of a
     member vector, never a scalar member); the storage stays valid and
     unresized until the next step().

   Sub-machine concept (an embeddable building block, e.g. gsRadixSelect):
   @code
   typedef <T> Scalar;
   bool step(gsCollectiveRequest<Scalar> & out); // true : out holds the next request, forward it
                                                 // false: finished, out untouched
   @endcode
   A sub-machine never produces DONE, and calls after the first false return
   false. A parent embeds it as a member, forwards its requests, and after
   false continues in the same step() call:
   @code
   gsCollectiveRequest<T> step()
   {
       gsCollectiveRequest<T> r;
       for (;;) switch (m_state)
       {
       case Box:      fillBox(); m_state = BoxDone; return gsCollectiveRequest<T>::min(m_box.data(), (int)m_box.size());
       case BoxDone:  m_sel.setup(...); m_state = Select; // fall through
       case Select:   if (m_sel.step(r)) return r; m_state = AfterSelect; break;
       case AfterSelect: useResults(m_sel); m_state = Finished; break;
       case Finished: return gsCollectiveRequest<T>::done();
       }
   }
   @endcode

   Rules for machine authors: the request kinds are int64 SUM, MIN/MAX of T
   (after canonicalising -0 to +0) and the final all-to-all exchange; there is
   no floating-point SUM (it is rank-order dependent). A rank without elements
   contributes numeric_limits<T>::max() / lowest() to MIN / MAX boxes, never 0.
   Active sets, numbers of passes and early exits are decided from reduced data
   only, and data-dependent failures are agreed through a reduction before the
   next collective, never thrown on a subset of ranks.
*/

namespace gismo
{

/// Kinds of collective a distributed partition machine can request.
struct gsCollective
{
    enum Kind
    {
        SUM_I64 = 0,     ///< in-place element-wise sum of len int64 values (i64[len])
        MIN_T,           ///< in-place element-wise min of len T values (t[len])
        MAX_T,           ///< in-place element-wise max of len T values (t[len])
        ALLTOALL_INT,    ///< one int per destination: sendCounts[len] -> recvCounts[len], len == nranks
        ALLTOALLV_INDEX, ///< index_t ids: sendIds/sendCounts/sendDispls -> recvIds/recvCounts/recvDispls, len == nranks
        DONE             ///< the machine has finished; no buffer
    };
};

/** @brief One collective request.

    SUM, MIN and MAX are in place, as gsMpiComm::sum/min/max(T* inout, int len);
    len == 0 is allowed (a no-op that must still be issued by all ranks).
    MIN / MAX buffers must already be +-0-canonicalised by the machine.
    ALLTOALL_INT sends one int to every destination (sendCounts[d]) and
    receives one from every source (recvCounts[s]). ALLTOALLV_INDEX sends
    sendCounts[d] ids starting at sendIds + sendDispls[d] to rank d and
    receives recvCounts[s] ids from rank s into recvIds + recvDispls[s]; the
    machine sizes recvIds (recvCapacity entries) from the recvCounts an earlier
    ALLTOALL_INT returned.

    All pointers point into storage owned by the issuing machine
    (std::vector::data() of a member, never a scalar member); they stay valid
    and unresized until the machine's next step(), which is when the result has
    been written in place.
*/
template<class T>
struct gsCollectiveRequest
{
    gsCollective::Kind kind;
    int        len;          ///< SUM/MIN/MAX: number of values; ALLTOALL*: number of ranks
    int64_t  * i64;          ///< SUM_I64 in/out [len]
    T        * t;            ///< MIN_T / MAX_T in/out [len]
    int      * sendCounts;   ///< ALLTOALL_INT in [len]; ALLTOALLV_INDEX in [len] (ids to destination d)
    int      * recvCounts;   ///< ALLTOALL_INT out [len]; ALLTOALLV_INDEX in [len] (ids from source s)
    int      * sendDispls;   ///< ALLTOALLV_INDEX in [len], offsets into sendIds
    int      * recvDispls;   ///< ALLTOALLV_INDEX in [len], offsets into recvIds
    index_t  * sendIds;      ///< ALLTOALLV_INDEX in
    index_t  * recvIds;      ///< ALLTOALLV_INDEX out
    int        recvCapacity; ///< ALLTOALLV_INDEX: size of the recvIds buffer (checked by runLockstep)

    /// Request of kind DONE: len = 0, all pointers null, recvCapacity = 0.
    gsCollectiveRequest()
    : kind(gsCollective::DONE), len(0), i64(NULL), t(NULL), sendCounts(NULL),
      recvCounts(NULL), sendDispls(NULL), recvDispls(NULL), sendIds(NULL),
      recvIds(NULL), recvCapacity(0)
    { }

    static gsCollectiveRequest done() { return gsCollectiveRequest(); }

    static gsCollectiveRequest sum(int64_t * inout, int len)
    {
        gsCollectiveRequest r;
        r.kind = gsCollective::SUM_I64; r.len = len; r.i64 = inout;
        return r;
    }

    static gsCollectiveRequest min(T * inout, int len)
    {
        gsCollectiveRequest r;
        r.kind = gsCollective::MIN_T; r.len = len; r.t = inout;
        return r;
    }

    static gsCollectiveRequest max(T * inout, int len)
    {
        gsCollectiveRequest r;
        r.kind = gsCollective::MAX_T; r.len = len; r.t = inout;
        return r;
    }

    static gsCollectiveRequest alltoallCounts(int * send, int * recv, int nranks)
    {
        gsCollectiveRequest r;
        r.kind = gsCollective::ALLTOALL_INT; r.len = nranks;
        r.sendCounts = send; r.recvCounts = recv;
        return r;
    }

    static gsCollectiveRequest alltoallvIds(index_t * sendIds, int * sendCounts, int * sendDispls,
                                            index_t * recvIds, int * recvCounts, int * recvDispls,
                                            int nranks, int recvCapacity)
    {
        gsCollectiveRequest r;
        r.kind = gsCollective::ALLTOALLV_INDEX; r.len = nranks;
        r.sendIds = sendIds; r.sendCounts = sendCounts; r.sendDispls = sendDispls;
        r.recvIds = recvIds; r.recvCounts = recvCounts; r.recvDispls = recvDispls;
        r.recvCapacity = recvCapacity;
        return r;
    }
};

/** @brief Runs a sub-machine as a top-level machine.

    The adapter holds a pointer to the sub-machine (it does not own it).
    step() forwards sub.step(r) and returns a DONE request once the
    sub-machine has finished.
*/
template<class Sub>
class gsAsMachine
{
public:
    typedef typename Sub::Scalar Scalar;

    explicit gsAsMachine(Sub & sub) : m_sub(&sub) { }

    gsCollectiveRequest<Scalar> step()
    {
        gsCollectiveRequest<Scalar> r;
        if (m_sub->step(r))
            return r;
        return gsCollectiveRequest<Scalar>::done();
    }

private:
    Sub * m_sub;
};

/** @brief Forwards each request of a machine to a communicator until DONE.

    Comm must provide rank(), size(), sum/min/max(X* inout, int len),
    alltoall(int*, int*, int, int) and
    alltoallv(index_t*, int*, int*, index_t*, int*, int*), as gsMpiComm and
    gsSerialComm do. The machine supplies the agreement between the ranks;
    this driver adds nothing.
*/
template<class Machine, class Comm>
void runCollectives(Machine & m, const Comm & comm)
{
    typedef typename Machine::Scalar T;
    for (;;)
    {
        gsCollectiveRequest<T> r = m.step();
        switch (r.kind)
        {
        case gsCollective::SUM_I64:
            comm.sum(r.i64, r.len);
            break;
        case gsCollective::MIN_T:
            comm.min(r.t, r.len);
            break;
        case gsCollective::MAX_T:
            comm.max(r.t, r.len);
            break;
        case gsCollective::ALLTOALL_INT:
            GISMO_ENSURE(r.len == (int)comm.size(),
                         "ALLTOALL_INT request length differs from the communicator size.");
            comm.alltoall(r.sendCounts, r.recvCounts, 1, 1);
            break;
        case gsCollective::ALLTOALLV_INDEX:
            GISMO_ENSURE(r.len == (int)comm.size(),
                         "ALLTOALLV_INDEX request length differs from the communicator size.");
            comm.alltoallv(r.sendIds, r.sendCounts, r.sendDispls,
                           r.recvIds, r.recvCounts, r.recvDispls);
            break;
        case gsCollective::DONE:
        default:
            return;
        }
    }
}

/** @brief In-process P-rank driver (testing and verification driver).

    ms[r] plays rank r of P = ms.size(). Each step calls step() on every
    machine in rank order 0..P-1 and then enforces the collective protocol:
    every rank issues the same kind (a mix of DONE and non-DONE is a mismatch)
    and the same length; reductions are carried out element-wise in rank order
    (SUM in int64; MIN / MAX starting from rank 0's value) and written to every
    machine's buffer; the all-to-alls require a consistent count matrix
    (sendCounts_s[d] == recvCounts_d[s]), non-negative counts and
    displacements, and receive ranges inside recvCapacity. Violations throw
    std::runtime_error (GISMO_ENSURE).

    @return the number of collective steps executed (DONE not counted).

    Complexity per step: O(P len) for reductions, O(P^2) + O(total ids) for the
    all-to-alls.
*/
template<class Machine>
index_t runLockstep(const std::vector<Machine*> & ms)
{
    typedef typename Machine::Scalar T;
    typedef gsCollectiveRequest<T> Req;
    const size_t P = ms.size();
    GISMO_ENSURE(P > 0, "runLockstep needs at least one machine.");

    index_t steps = 0;
    std::vector<Req> rq(P);
    for (;;)
    {
        for (size_t r = 0; r != P; ++r)
            rq[r] = ms[r]->step();

        for (size_t r = 1; r != P; ++r)
            GISMO_ENSURE(rq[r].kind == rq[0].kind,
                         "Collective kind mismatch between ranks.");
        if (rq[0].kind == gsCollective::DONE)
            return steps;

        switch (rq[0].kind)
        {
        case gsCollective::SUM_I64:
        {
            for (size_t r = 1; r != P; ++r)
                GISMO_ENSURE(rq[r].len == rq[0].len, "SUM_I64 length mismatch between ranks.");
            const int n = rq[0].len;
            if (n > 0)
            {
                std::vector<int64_t> acc(rq[0].i64, rq[0].i64 + n);
                for (size_t r = 1; r != P; ++r)
                    for (int i = 0; i != n; ++i)
                        acc[i] += rq[r].i64[i];
                for (size_t r = 0; r != P; ++r)
                    std::copy(acc.begin(), acc.end(), rq[r].i64);
            }
            break;
        }
        case gsCollective::MIN_T:
        case gsCollective::MAX_T:
        {
            const bool isMin = (rq[0].kind == gsCollective::MIN_T);
            for (size_t r = 1; r != P; ++r)
                GISMO_ENSURE(rq[r].len == rq[0].len, "MIN_T/MAX_T length mismatch between ranks.");
            const int n = rq[0].len;
            if (n > 0)
            {
                std::vector<T> acc(rq[0].t, rq[0].t + n);
                for (size_t r = 1; r != P; ++r)
                    for (int i = 0; i != n; ++i)
                    {
                        const T v = rq[r].t[i];
                        if (isMin ? (v < acc[i]) : (acc[i] < v))
                            acc[i] = v;
                    }
                for (size_t r = 0; r != P; ++r)
                    std::copy(acc.begin(), acc.end(), rq[r].t);
            }
            break;
        }
        case gsCollective::ALLTOALL_INT:
        {
            for (size_t r = 0; r != P; ++r)
                GISMO_ENSURE(rq[r].len == (int)P, "ALLTOALL_INT length differs from the number of ranks.");
            std::vector<int> mat(P * P); // mat[s*P+d] = sendCounts_s[d]
            for (size_t s = 0; s != P; ++s)
                for (size_t d = 0; d != P; ++d)
                    mat[s * P + d] = rq[s].sendCounts[d];
            for (size_t d = 0; d != P; ++d)
                for (size_t s = 0; s != P; ++s)
                    rq[d].recvCounts[s] = mat[s * P + d];
            break;
        }
        case gsCollective::ALLTOALLV_INDEX:
        {
            for (size_t r = 0; r != P; ++r)
                GISMO_ENSURE(rq[r].len == (int)P, "ALLTOALLV_INDEX length differs from the number of ranks.");
            for (size_t s = 0; s != P; ++s)
                for (size_t d = 0; d != P; ++d)
                {
                    const int cs = rq[s].sendCounts[d];
                    const int cr = rq[d].recvCounts[s];
                    GISMO_ENSURE(cs == cr, "ALLTOALLV_INDEX count matrix is inconsistent.");
                    GISMO_ENSURE(cs >= 0 && rq[s].sendDispls[d] >= 0 && rq[d].recvDispls[s] >= 0,
                                 "ALLTOALLV_INDEX negative count or displacement.");
                    GISMO_ENSURE((int64_t)rq[d].recvDispls[s] + cr <= (int64_t)rq[d].recvCapacity,
                                 "ALLTOALLV_INDEX receive range exceeds the receive buffer.");
                }
            for (size_t s = 0; s != P; ++s)
                for (size_t d = 0; d != P; ++d)
                {
                    const int c = rq[s].sendCounts[d];
                    if (c > 0)
                        std::copy(rq[s].sendIds + rq[s].sendDispls[d],
                                  rq[s].sendIds + rq[s].sendDispls[d] + c,
                                  rq[d].recvIds + rq[d].recvDispls[s]);
                }
            break;
        }
        case gsCollective::DONE:
        default:
            return steps;
        }
        ++steps;
    }
}

/// Forwards to the pointer version.
template<class Machine>
typename std::enable_if<!std::is_pointer<Machine>::value, index_t>::type
runLockstep(std::vector<Machine> & ms)
{
    std::vector<Machine*> p(ms.size());
    for (size_t r = 0; r != ms.size(); ++r)
        p[r] = &ms[r];
    return runLockstep(p);
}

/** @brief Order-preserving map T -> uint64_t.

    a < b <=> encode(a) < encode(b) and a == b <=> encode(a) == encode(b), on
    finite and infinite values. NaN is outside the domain (callers must ensure
    finiteness first).

    The primary template is available for every coefficient type (long double,
    multiprecision, rational, posit, ...) and fails at run time with
    GISMO_ENSURE, so that the distributed partition compiles for them while
    only the serial one is usable.

    Specialisations for IEEE float and double: the value is canonicalised
    (-0.0 -> +0.0 by `if (x == 0) x = 0;`), its IEEE bits b of width w
    (64 / 32) are read with std::memcpy and mapped to
    `(b >> (w-1)) ? ~b : (b | (1 << (w-1)))`: negative values flip all bits,
    non-negative values set the sign bit. A float key is zero-extended to
    uint64_t, which preserves the order. encode is branch-light and does not
    throw.
*/
template<class T> struct gsOrderedKey
{
    static const bool available = false;
    static uint64_t encode(const T &)
    {
        GISMO_ENSURE(false, "the distributed partition needs IEEE float or double coordinates; use the serial constructor");
        return 0;
    }
};

template<class T> const bool gsOrderedKey<T>::available;

template<> struct gsOrderedKey<double>
{
    static const bool available = true;
    static uint64_t encode(double x);
};

template<> struct gsOrderedKey<float>
{
    static const bool available = true;
    static uint64_t encode(float x);
};

inline uint64_t gsOrderedKey<double>::encode(double x)
{
    if (x == 0) x = 0;
    uint64_t b;
    std::memcpy(&b, &x, sizeof(b));
    return (b >> 63) ? ~b : (b | (uint64_t(1) << 63));
}

inline uint64_t gsOrderedKey<float>::encode(float x)
{
    if (x == 0) x = 0;
    uint32_t b;
    std::memcpy(&b, &x, sizeof(b));
    b = (b >> 31) ? ~b : (b | (uint32_t(1) << 31));
    return uint64_t(b);
}

/** @brief Full selection key: lexicographic (hi, lo).

    hi = gsOrderedKey of a coordinate or a space-filling-curve key,
    lo = element id as the unsigned type of index_t width. Ids are
    non-negative, so unsigned order equals id order, and keys of distinct
    elements are distinct. The composite key has 64 + 8 sizeof(index_t) bits;
    the radix select consumes it in 8-bit digits, most significant first:
    digit d (0-based) is byte 7-d of hi for d < 8, else byte
    sizeof(index_t)-1-(d-8) of lo.
*/
struct gsPartitionKey
{
    typedef std::make_unsigned<index_t>::type lo_type;
    uint64_t hi;
    lo_type  lo;
};

inline bool operator< (const gsPartitionKey & a, const gsPartitionKey & b)
{ return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo); }
inline bool operator==(const gsPartitionKey & a, const gsPartitionKey & b)
{ return a.hi == b.hi && a.lo == b.lo; }
inline bool operator!=(const gsPartitionKey & a, const gsPartitionKey & b) { return !(a == b); }
inline bool operator<=(const gsPartitionKey & a, const gsPartitionKey & b) { return !(b < a); }
inline bool operator> (const gsPartitionKey & a, const gsPartitionKey & b) { return b < a; }
inline bool operator>=(const gsPartitionKey & a, const gsPartitionKey & b) { return !(a < b); }

/** @brief Multi-target radix select, a sub-machine of the collective protocol.

    The candidates of target t are all elements, on all ranks, with
    elemGroup == targetGroup[t]; their keys are globally distinct (unique
    ids). Two selection modes:
    - Count: the candidate key of 0-based rank r in key order; not found if
      r < 0 or r >= the number of candidates.
    - Weighted: the smallest candidate key f with W(<= f) >= tau, where W is
      the inclusive prefix weight over the target's candidates in key order;
      not found if no such f exists (tau larger than the total candidate
      weight, including all-zero weights with tau >= 1, or no candidates).
      tau <= 0 is satisfied by the smallest candidate key (the target is
      treated as count rank 0).

    Typical uses: recursive coordinate bisection - one target per node, group =
    node index, Count, r = size of the left part, left = key < key(t); curve
    partitioning - one group, P-1 targets tau_j, Weighted, label =
    #{j : key > key(j)}.

    Algorithm: most-significant-digit radix over the composite key, all targets
    advanced together. Per target the state is the number of fixed digits d
    (equal for all active targets in a pass), the prefix (key with the first d
    digits fixed, the rest zero) and the remaining rank / threshold. The
    histogram rows are the distinct (group, prefix) pairs of the active
    targets, sorted by (group, masked key); targets with the same pair share a
    row. A row holds `bins` counts (Count) or `bins` counts followed by `bins`
    weight sums (Weighted). In a pass each element with a group masks its key
    to depth d, looks the pair up in the row table (binary search) and adds to
    the bin of its digit d. One SUM_I64 request carries
    [rows x rowLen][2 slots per pending target, in target order]. After the
    reduction each active target picks, from its row only, the smallest digit b
    with cumulative count > r (Count) or cumulative weight >= tau (Weighted)
    and subtracts the mass below b. A target whose chosen bucket has reduced
    count 1 before the last digit is pending: in the next request the unique
    rank owning that element writes its hi / lo into the two slots (all other
    ranks write 0), so the sum is exact and the key is known without further
    digit passes. At depth 0 the row totals decide `found`. Every decision
    (found, digit, exit, number of passes) is taken from reduced data and
    replicated inputs only, so all ranks run the same requests. Duplicate keys
    (non-unique ids) are detected, on reduced data so that all ranks fail
    together, only when the duplicated key is the one a target selects;
    elsewhere they go unnoticed, so callers must supply unique ids.

    Bounds: at most maxPasses SUM_I64 requests per run. Per request O(nLocal
    log G + G bins + nTargets bins) local work, G <= nTargets rows; request
    size <= nTargets rowLen + 2 nTargets int64 (rowLen = 256 Count / 512
    Weighted; e.g. curves with P = 256: 255 * 512 * 8 B ~ 1 MiB); memory
    O(G bins + nTargets) beyond the caller's arrays; no allocation sized by
    the global number of elements.

    @tparam T only fixes the request type, so the select can be embedded in a
    gsCollectiveRequest<T> machine.
*/
template<class T>
class gsRadixSelect
{
public:
    typedef T Scalar;
    enum Mode { Count = 0, Weighted = 1 };
    static const unsigned digitBits = 8;
    static const unsigned bins      = 1u << digitBits;                       // 256
    static const unsigned maxPasses = (64 + 8 * sizeof(index_t)) / digitBits; // 12 for 32-bit index_t, 16 for 64-bit

    /// Finished, zero targets: step() returns false at once.
    gsRadixSelect() : m_state(Idle), m_mode(Count), m_nLocal(0), m_keys(NULL),
                      m_group(NULL), m_weights(NULL), m_d(0), m_issued(0)
    { }

    /** @brief Sets up a run.

        Replicated inputs (identical on all ranks): mode, targetGroup,
        targetParam. Local inputs (caller-owned, must stay alive and unchanged
        until step() returns false): keys[nLocal], elemGroup[nLocal] (nullptr:
        every element is in group 0; a value no target names, e.g. -1, means
        "candidate of no target"), weights[nLocal] (Weighted only, int64 >= 0;
        nullptr in Count mode).

        targetGroup[t] is the candidate group of target t; targetParam[t] is
        the 0-based rank r (Count) or the threshold tau (Weighted).
    */
    void setup(Mode mode, index_t nLocal, const gsPartitionKey * keys, const index_t * elemGroup,
               const int64_t * weights, std::vector<index_t> targetGroup, std::vector<int64_t> targetParam)
    {
        GISMO_ASSERT(targetGroup.size() == targetParam.size(), "Target arrays differ in size.");
        GISMO_ASSERT(mode == Count || weights != NULL || nLocal == 0, "Weighted select needs weights.");
        m_mode = mode; m_nLocal = nLocal; m_keys = keys; m_group = elemGroup; m_weights = weights;
        m_tGroup.swap(targetGroup);
        m_rem.swap(targetParam);
        const size_t nT = m_tGroup.size();
        m_prefix.assign(nT, gsPartitionKey());
        for (size_t t = 0; t != nT; ++t) { m_prefix[t].hi = 0; m_prefix[t].lo = 0; }
        m_status.assign(nT, Active);
        m_found.assign(nT, 0);
        m_cnt.assign(nT, 1);
        m_result = m_prefix;
        m_row.assign(nT, 0);
        m_pendU.assign(nT, 0);
        for (size_t t = 0; t != nT; ++t)
        {
            if (mode == Count)
            {
                if (m_rem[t] < 0) m_status[t] = Done;
            }
            else if (m_rem[t] <= 0)
                m_rem[t] = 0;
            else
                m_cnt[t] = 0;
        }
        m_d = 0;
        m_issued = 0;
        m_state = nT ? Begin : Idle;
    }

    /// Sub-machine step, see \ref gsDistributedPartitionProtocol.
    bool step(gsCollectiveRequest<T> & out)
    {
        for (;;) switch (m_state)
        {
        case Idle:
            return false;
        case Begin:
            if (!buildPass())
            {
                m_state = Idle;
                return false;
            }
            ++m_issued;
            m_state = Reduced;
            out = gsCollectiveRequest<T>::sum(m_buf.data(), (int)m_buf.size());
            return true;
        case Reduced:
            consume();
            ++m_d;
            m_state = Begin;
            break;
        }
    }

    bool finished() const { return m_state == Idle; }
    index_t numTargets() const { return (index_t)m_tGroup.size(); }

    /// Whether target t has a result; valid once finished.
    bool found(index_t t) const { return m_found[t] != 0; }

    /// Key of target t; valid once finished and found(t).
    const gsPartitionKey & key(index_t t) const { return m_result[t]; }

    /// Number of SUM_I64 requests of the last run (diagnostics and tests).
    index_t requestsIssued() const { return m_issued; }

private:
    enum State  { Idle, Begin, Reduced };
    enum Status { Active, Pending, Done };

    struct RowKey { index_t g; gsPartitionKey k; };
    static bool rowLess(const RowKey & a, const RowKey & b)
    { return a.g < b.g || (a.g == b.g && a.k < b.k); }

    typedef gsPartitionKey::lo_type lo_type;

    static unsigned digitOf(const gsPartitionKey & k, unsigned d)
    {
        if (d < 8) return unsigned((k.hi >> (8 * (7 - d))) & 0xFFu);
        return unsigned((k.lo >> (8 * (sizeof(index_t) - 1 - (d - 8)))) & 0xFFu);
    }

    static void setDigit(gsPartitionKey & k, unsigned d, unsigned v)
    {
        if (d < 8) k.hi |= uint64_t(v) << (8 * (7 - d));
        else       k.lo |= lo_type(lo_type(v) << (8 * (sizeof(index_t) - 1 - (d - 8))));
    }

    /// Key with only the first d digits kept.
    static gsPartitionKey maskKey(const gsPartitionKey & k, unsigned d)
    {
        gsPartitionKey r;
        if (d >= 8)      r.hi = k.hi;
        else if (d == 0) r.hi = 0;
        else             r.hi = k.hi & (~uint64_t(0) << (64 - 8 * d));
        const unsigned kd = d > 8 ? d - 8 : 0;
        if (kd == 0)                      r.lo = 0;
        else if (kd >= sizeof(index_t))   r.lo = k.lo;
        else r.lo = lo_type(k.lo & lo_type(~lo_type(0) << (8 * (sizeof(index_t) - kd))));
        return r;
    }

    /// Builds the row tables and the request buffer of the pass at depth m_d from the local data;
    /// false if no target is active or pending. O(nLocal log G + G bins).
    bool buildPass()
    {
        const size_t nT = m_tGroup.size();
        m_rows.clear(); m_pendTab.clear(); m_pendList.clear();
        for (size_t t = 0; t != nT; ++t)
        {
            RowKey q; q.g = m_tGroup[t]; q.k = m_prefix[t];
            if (m_status[t] == Active)       m_rows.push_back(q);
            else if (m_status[t] == Pending) { m_pendTab.push_back(q); m_pendList.push_back((index_t)t); }
        }
        if (m_rows.empty() && m_pendList.empty())
            return false;

        sortUnique(m_rows);
        sortUnique(m_pendTab);
        for (size_t t = 0; t != nT; ++t)
        {
            RowKey q; q.g = m_tGroup[t]; q.k = m_prefix[t];
            if (m_status[t] == Active)
                m_row[t] = (index_t)(std::lower_bound(m_rows.begin(), m_rows.end(), q, rowLess) - m_rows.begin());
            else if (m_status[t] == Pending)
                m_pendU[t] = (index_t)(std::lower_bound(m_pendTab.begin(), m_pendTab.end(), q, rowLess) - m_pendTab.begin());
        }

        const size_t rowLen = (m_mode == Weighted) ? 2 * bins : bins;
        const size_t base   = m_rows.size() * rowLen;
        m_buf.assign(base + 2 * m_pendList.size(), 0);
        std::vector<int64_t> pendVal(2 * m_pendTab.size(), 0);

        const bool haveRows = !m_rows.empty(), havePend = !m_pendTab.empty();
        for (index_t e = 0; e != m_nLocal; ++e)
        {
            RowKey q;
            q.g = m_group ? m_group[e] : 0;
            q.k = maskKey(m_keys[e], m_d);
            if (haveRows)
            {
                typename std::vector<RowKey>::const_iterator it =
                    std::lower_bound(m_rows.begin(), m_rows.end(), q, rowLess);
                if (it != m_rows.end() && !rowLess(q, *it))
                {
                    const size_t off = size_t(it - m_rows.begin()) * rowLen;
                    const unsigned dg = digitOf(m_keys[e], m_d);
                    ++m_buf[off + dg];
                    if (m_mode == Weighted)
                        m_buf[off + bins + dg] += m_weights[e];
                }
            }
            if (havePend)
            {
                typename std::vector<RowKey>::const_iterator it =
                    std::lower_bound(m_pendTab.begin(), m_pendTab.end(), q, rowLess);
                if (it != m_pendTab.end() && !rowLess(q, *it))
                {
                    const size_t u = size_t(it - m_pendTab.begin());
                    int64_t h;
                    const uint64_t hi = m_keys[e].hi;
                    std::memcpy(&h, &hi, sizeof(h));
                    pendVal[2 * u]     = h;
                    pendVal[2 * u + 1] = (int64_t)m_keys[e].lo;
                }
            }
        }
        for (size_t k = 0; k != m_pendList.size(); ++k)
        {
            const size_t u = (size_t)m_pendU[m_pendList[k]];
            m_buf[base + 2 * k]     = pendVal[2 * u];
            m_buf[base + 2 * k + 1] = pendVal[2 * u + 1];
        }
        return true;
    }

    static void sortUnique(std::vector<RowKey> & v)
    {
        std::sort(v.begin(), v.end(), rowLess);
        v.erase(std::unique(v.begin(), v.end(), eqRow), v.end());
    }
    static bool eqRow(const RowKey & a, const RowKey & b) { return !rowLess(a, b) && !rowLess(b, a); }

    /// Advances all targets using the reduced buffer of the pass at depth m_d.
    void consume()
    {
        const size_t nT = m_tGroup.size();
        const size_t rowLen = (m_mode == Weighted) ? 2 * bins : bins;
        const size_t base   = m_rows.size() * rowLen;

        for (size_t k = 0; k != m_pendList.size(); ++k)
        {
            const size_t t = (size_t)m_pendList[k];
            const int64_t h = m_buf[base + 2 * k];
            uint64_t hi;
            std::memcpy(&hi, &h, sizeof(hi));
            m_result[t].hi = hi;
            m_result[t].lo = (lo_type)m_buf[base + 2 * k + 1];
            m_found[t]  = 1;
            m_status[t] = Done;
        }

        for (size_t t = 0; t != nT; ++t)
        {
            if (m_status[t] != Active)
                continue;
            const int64_t * c = m_buf.data() + size_t(m_row[t]) * rowLen;
            const int64_t * w = c + bins;
            int64_t cum = 0;
            unsigned b = 0;
            if (m_cnt[t])
                for (; b != bins; ++b) { if (cum + c[b] > m_rem[t]) break; cum += c[b]; }
            else
                for (; b != bins; ++b) { if (cum + w[b] >= m_rem[t]) break; cum += w[b]; }

            if (b == bins)
            {
                GISMO_ENSURE(m_d == 0, "radix select: inconsistent reduced histogram.");
                m_status[t] = Done;
                continue;
            }
            m_rem[t] -= cum;
            setDigit(m_prefix[t], m_d, b);
            const int64_t cnt = c[b];
            if (m_d + 1 == maxPasses)
            {
                GISMO_ENSURE(cnt == 1, "radix select: duplicate keys (ids not unique).");
                m_result[t] = m_prefix[t];
                m_found[t]  = 1;
                m_status[t] = Done;
            }
            else if (cnt == 1)
                m_status[t] = Pending;
        }
    }

private:
    State    m_state;
    Mode     m_mode;
    index_t  m_nLocal;
    const gsPartitionKey * m_keys;      // [nLocal], caller-owned
    const index_t        * m_group;     // [nLocal] or null, caller-owned
    const int64_t        * m_weights;   // [nLocal] or null, caller-owned
    unsigned m_d;                       // number of digits fixed in the current pass
    index_t  m_issued;                  // SUM_I64 requests issued in the current run

    std::vector<index_t>       m_tGroup;  // [nT] candidate group of each target
    std::vector<int64_t>       m_rem;     // [nT] remaining rank / threshold
    std::vector<char>          m_cnt;     // [nT] 1: count criterion, 0: weight criterion
    std::vector<char>          m_status;  // [nT] Status
    std::vector<char>          m_found;   // [nT]
    std::vector<gsPartitionKey> m_prefix; // [nT] key with the first m_d digits fixed
    std::vector<gsPartitionKey> m_result; // [nT]
    std::vector<index_t>       m_row;     // [nT] row of an active target in m_rows
    std::vector<index_t>       m_pendU;   // [nT] entry of a pending target in m_pendTab

    std::vector<RowKey>        m_rows;     // sorted unique (group, prefix) of the active targets
    std::vector<RowKey>        m_pendTab;  // sorted unique (group, prefix) of the pending targets
    std::vector<index_t>       m_pendList; // pending targets in target order
    std::vector<int64_t>       m_buf;      // request buffer: [rows x rowLen][2 slots per pending target]
};

template<class T> const unsigned gsRadixSelect<T>::digitBits;
template<class T> const unsigned gsRadixSelect<T>::bins;
template<class T> const unsigned gsRadixSelect<T>::maxPasses;

/** @brief Routing sub-machine: sends the id of local element i to rank dest[i].

    Afterwards ownIds() holds the received ids, sorted ascending, and the
    agreed coverage check (sum over the ranks of the own counts == numGlobal)
    has passed.

    Requests, in this order: ALLTOALL_INT (len nranks, the per-destination
    counts), ALLTOALLV_INDEX (len nranks, the ids), SUM_I64 (len 1, the
    coverage check). The local id of element i is idLo + i, or ids[i] when the
    ids overload of setup() is used.

    Complexity: O(nLocal + nranks) for the packing, O(own log own) for the
    sort of the received ids. Memory: O(nLocal + nranks + own); nothing is
    sized by the global number of elements.
*/
template<class T>
class gsOwnIdsExchange
{
public:
    typedef T Scalar;

    /// Finished: step() returns false.
    gsOwnIdsExchange() : m_state(Idle), m_idLo(0), m_ids(NULL), m_nLocal(0),
                         m_dest(NULL), m_nranks(0), m_numGlobal(0), m_expected(-1)
    { }

    /** @brief Sets up a run with the contiguous ids idLo .. idLo+nLocal-1.

        \a dest[nLocal] (caller-owned, alive and unchanged until step() returns
        false) holds the destination rank of each local element, in [0,
        nranks). \a expectedOwn >= 0: replicated expected own count (the caller
        has checked it against INT_MAX on all ranks); the receive buffer is
        sized from it and its equality with the received counts is enforced.
        expectedOwn < 0: the buffer is sized from the received counts.
    */
    void setup(index_t idLo, index_t nLocal, const index_t * dest, index_t nranks,
               int64_t numGlobal, int64_t expectedOwn = -1)
    { init(idLo, NULL, nLocal, dest, nranks, numGlobal, expectedOwn); }

    /// As above, with explicit ids[nLocal] (caller-owned, alive and unchanged until step() returns false).
    void setup(const index_t * ids, index_t nLocal, const index_t * dest, index_t nranks,
               int64_t numGlobal, int64_t expectedOwn = -1)
    { init(0, ids, nLocal, dest, nranks, numGlobal, expectedOwn); }

    /// Sub-machine step, see \ref gsDistributedPartitionProtocol.
    bool step(gsCollectiveRequest<T> & out)
    {
        for (;;) switch (m_state)
        {
        case Idle:
            return false;
        case Begin:
        {
            const size_t P = (size_t)m_nranks;
            m_sendCounts.assign(P, 0);
            m_recvCounts.assign(P, 0);
            m_sendDispls.assign(P, 0);
            m_recvDispls.assign(P, 0);
            for (index_t i = 0; i != m_nLocal; ++i)
            {
                GISMO_ASSERT(m_dest[i] >= 0 && m_dest[i] < m_nranks, "Destination rank out of range.");
                ++m_sendCounts[m_dest[i]];
            }
            m_state = Counts;
            out = gsCollectiveRequest<T>::alltoallCounts(m_sendCounts.data(), m_recvCounts.data(), (int)P);
            return true;
        }
        case Counts:
        {
            const size_t P = (size_t)m_nranks;
            int64_t tot = 0;
            for (size_t s = 0; s != P; ++s)
                tot += m_recvCounts[s];
            GISMO_ENSURE(m_expected < 0 || tot == m_expected,
                         "Received id count differs from the expected part size.");
            GISMO_ENSURE(tot <= (int64_t)std::numeric_limits<int>::max(),
                         "Own part exceeds INT_MAX elements.");
            int off = 0;
            for (size_t s = 0; s != P; ++s) { m_recvDispls[s] = off; off += m_recvCounts[s]; }
            off = 0;
            for (size_t d = 0; d != P; ++d) { m_sendDispls[d] = off; off += m_sendCounts[d]; }

            m_own.assign((size_t)tot, 0);
            m_sendIds.assign((size_t)m_nLocal, 0);
            std::vector<int> pos(m_sendDispls);
            for (index_t i = 0; i != m_nLocal; ++i)
                m_sendIds[pos[m_dest[i]]++] = m_ids ? m_ids[i] : m_idLo + i;
            m_state = Ids;
            out = gsCollectiveRequest<T>::alltoallvIds(m_sendIds.data(), m_sendCounts.data(), m_sendDispls.data(),
                                                       m_own.data(), m_recvCounts.data(), m_recvDispls.data(),
                                                       (int)P, (int)m_own.size());
            return true;
        }
        case Ids:
            std::sort(m_own.begin(), m_own.end());
            std::vector<index_t>().swap(m_sendIds);
            m_cov.assign(1, (int64_t)m_own.size());
            m_state = Cover;
            out = gsCollectiveRequest<T>::sum(m_cov.data(), 1);
            return true;
        case Cover:
            GISMO_ENSURE(m_cov[0] == m_numGlobal,
                         "The own parts of the ranks do not cover numGlobal elements.");
            m_state = Idle;
            break;
        }
    }

    /// Received ids, sorted ascending; valid once step() has returned false.
    std::vector<index_t> & ownIds() { return m_own; }
    const std::vector<index_t> & ownIds() const { return m_own; }

private:
    enum State { Idle, Begin, Counts, Ids, Cover };

    void init(index_t idLo, const index_t * ids, index_t nLocal, const index_t * dest,
              index_t nranks, int64_t numGlobal, int64_t expectedOwn)
    {
        GISMO_ASSERT(nranks >= 1 && nLocal >= 0, "Invalid exchange size.");
        m_idLo = idLo; m_ids = ids; m_nLocal = nLocal; m_dest = dest;
        m_nranks = nranks; m_numGlobal = numGlobal; m_expected = expectedOwn;
        m_own.clear();
        m_state = Begin;
    }

    State    m_state;
    index_t  m_idLo;
    const index_t * m_ids;     // [nLocal] or null (contiguous ids), caller-owned
    index_t  m_nLocal;
    const index_t * m_dest;    // [nLocal], caller-owned
    index_t  m_nranks;
    int64_t  m_numGlobal;
    int64_t  m_expected;       // replicated expected own count, or -1

    std::vector<int>     m_sendCounts; // [nranks] ids to every destination
    std::vector<int>     m_recvCounts; // [nranks] ids from every source
    std::vector<int>     m_sendDispls; // [nranks] exclusive prefix sums of m_sendCounts
    std::vector<int>     m_recvDispls; // [nranks] exclusive prefix sums of m_recvCounts
    std::vector<index_t> m_sendIds;    // [nLocal] ids grouped by destination, local order within a group
    std::vector<index_t> m_own;        // [own count] received ids, sorted after the exchange
    std::vector<int64_t> m_cov;        // [1] own count, then the global sum
};

/** @brief Distributed recursive coordinate bisection: the resumable machine
    (collective request protocol, \ref gsDistributedPartitionProtocol) whose
    result equals gsGeometricPartitioner<T>::rcbSplit(C, w, weighted, labels,
    order, 0, N, nparts, 0) on the concatenation of all ranks' slices.

    The number of parts equals the number of ranks: part p is owned by rank p.
    After DONE, ownElements() holds the ids e of the rank's part, i.e. those
    with labels[e] == rank, sorted ascending; localLabels() holds the label of
    every local element. The sets are bit-identical to the replicated
    reference for the same centroids and weights.

    <b>Algorithm.</b> The bisection tree is replicated: a node table of
    {label0, k, n} (k parts left to label, n global elements; at most
    2 nparts - 1 entries). Each local element knows its node. Per level, the
    active nodes (k >= 2, n >= 1) are processed together, one collective per
    kind for all of them:
    -# Box. MIN and MAX of the +-0-canonicalised per-rank coordinate extrema,
       A x geoDim values each, layout [slot*geoDim + d] (ranks without an
       element of a node contribute max() / lowest()). The axis is the
       longest one, ties to the lowest axis, exactly as rcbBisect.
    -# Weighted only. SUM_I64 (A) of the node weight sums gives target =
       totalW kL / k in int64. target == 0 selects the exact count split. The
       others run the 16 coordinate-bisection iterations of rcbBisect, each
       one SUM_I64 (B) of the left weights of the B cut nodes (layout: cut
       nodes in slot order), with every rank updating lo / hi / midc / cut
       from the reduced sums by the same expressions; then SUM_I64 (B) of the
       left counts. A cut with 0 or n elements on the left (n >= 2) falls
       back to the exact split.
    -# Exact split: the left part consists of the nL = (n kL + k/2)/k smallest
       elements in the total order (coordinate, id). The keys
       {gsOrderedKey::encode(coordinate), id} carry that order; one batched
       gsRadixSelect (Count mode, group = node, rank nL) returns the key of
       rank nL, and left = key < selected key.
    -# Children {label0, kL, nL} and {label0 + kL, k - kL, n - nL} are
       appended in slot order and the local elements move to them.
    Every decision (active nodes, axes, modes, number of select passes) is
    taken from replicated data and reduced results only. Finally the leaf
    labels are routed with gsOwnIdsExchange (ALLTOALL_INT, ALLTOALLV_INDEX,
    SUM_I64).

    <b>Requests.</b> First a SUM_I64 of len 3 {#local columns with a
    non-finite coordinate, #local ids outside [0, numGlobal), nLocal}; the
    reduced values are checked on every rank, so non-finite input, ids out of
    range and a total slice size different from numGlobal make all ranks throw
    (GISMO_ENSURE) from the same step. The ids must be distinct across all
    ranks, i.e. together a permutation of 0 .. numGlobal-1. Duplicate ids are
    not detected. Then per level l: MIN_T and MAX_T
    (A_l geoDim), unweighted: the select's s_l <= gsRadixSelect::maxPasses
    SUM_I64; weighted: SUM_I64 (A_l), if B_l > 0 then 17 SUM_I64 (B_l), and
    the select. At the end ALLTOALL_INT, ALLTOALLV_INDEX and a SUM_I64 of len
    1. The total is 4 + L (2 + s_l) unweighted and 4 + L (3 + 17 [B_l > 0] +
    s_l) weighted, L <= ceil(log2 nparts) levels.

    <b>Complexity.</b> Local work per level O(n_loc (2 + 16 weighted)) plus the
    select's O(n_loc log G) per pass (G <= A target rows); at most 4 +
    ceil(log2 nparts) (3 + 17 + maxPasses) collectives, each of at most
    max(A geoDim, (nparts/2) 256 + 2 nparts/2) values with A <= nparts/2.
    Memory: O(n_loc) for the keys, groups, node indices and send ids, O(nparts)
    for the node table and the per-level slot tables, O(nparts 256) for the
    select rows and the own part (received ids); nothing is sized by the
    global number of elements.

    <b>Buffers.</b> m_node[i]: node of local element i. m_boxMin / m_boxMax:
    [slot*geoDim + d]. m_wsum, m_slot*: [slot]. m_cLo, m_cHi, m_cMid, m_cWL,
    m_cTarget, m_cCut, m_cSlot: [cut node j]. m_keys, m_group: [local element].

    Only floating-point types with an IEEE order-preserving key (float,
    double) are supported; the constructor ensures gsOrderedKey<T>::available.

    The machine holds pointers into its request buffers once running: do not
    copy or move it after the first step() (before it, copying and moving are
    fine). The centroid matrix, ids and weights passed to the constructor must
    stay alive and unchanged until step() returns DONE.
*/
template<class T>
class gsDistributedRcb
{
public:
    typedef T Scalar;

    /** @brief Constructs the machine for one rank.

        \param centroids geoDim x nLocal, column i = centroid of local element
               i (column-major); an empty slice may pass a geoDim x 0 or a
               0 x 0 matrix, nLocal = centroids.cols()
        \param geoDim number of coordinates (replicated)
        \param idLo    local element i has the global id idLo + i if \a ids is null
        \param ids     else ids[i] (nLocal entries, idLo ignored); the ids of all ranks together
               must be distinct and lie in [0, numGlobal)
        \param weights nLocal entries (as gsGeometricPartitioner's weights), read
               only if \a weighted; may be null if !weighted
        \param weighted, nparts, numGlobal replicated; \a nparts is the number
               of ranks
        \param rank in [0, nparts)
    */
    gsDistributedRcb(const gsMatrix<T> & centroids, short_t geoDim,
                     index_t idLo, const index_t * ids,
                     const index_t * weights, bool weighted,
                     index_t nparts, index_t rank, int64_t numGlobal)
    : m_state(Start), m_cdata(NULL), m_gd(geoDim), m_nLocal((index_t)centroids.cols()),
      m_idLo(idLo), m_ids(ids), m_w(weights), m_weighted(weighted),
      m_nparts(nparts), m_rank(rank), m_numGlobal(numGlobal), m_levels(0), m_iter(0)
    {
        GISMO_ENSURE(gsOrderedKey<T>::available,
                     "The distributed partition needs IEEE float or double coordinates; use the serial constructor.");
        GISMO_ENSURE(nparts >= 1 && rank >= 0 && rank < nparts, "Invalid number of parts or rank.");
        GISMO_ENSURE(numGlobal >= 0, "Negative global element count.");
        GISMO_ENSURE(geoDim >= 1, "Invalid geometric dimension.");
        GISMO_ENSURE(centroids.cols() == 0 || centroids.rows() == (index_t)geoDim,
                     "Centroid matrix does not have geoDim rows.");
        GISMO_ENSURE(!weighted || m_nLocal == 0 || weights != NULL, "Weighted partition needs weights.");
        m_cdata = m_nLocal > 0 ? centroids.data() : NULL;
        Node root; root.label0 = 0; root.k = nparts; root.n = numGlobal; root.split = false;
        m_nodes.push_back(root);
        m_node.assign((size_t)m_nLocal, 0);
    }

    /// Machine concept, see \ref gsDistributedPartitionProtocol.
    gsCollectiveRequest<T> step()
    {
        gsCollectiveRequest<T> r;
        for (;;) switch (m_state)
        {
        case Start:
            fillValidity();
            m_state = Validity;
            return gsCollectiveRequest<T>::sum(m_val.data(), 3);

        case Validity:
            GISMO_ENSURE(m_val[0] == 0, "Non-finite centroid coordinate in a slice.");
            GISMO_ENSURE(m_val[1] == 0, "Element id out of range [0, numGlobal).");
            GISMO_ENSURE(m_val[2] == m_numGlobal, "The slices do not cover numGlobal elements.");
            m_state = LevelBegin;
            break;

        case LevelBegin:
            collectActive();
            if (m_active.empty())
            {
                startExchange();
                m_state = Exchange;
                break;
            }
            fillBox();
            m_state = BoxMin;
            return gsCollectiveRequest<T>::min(m_boxMin.data(), (int)m_boxMin.size());

        case BoxMin:
            m_state = BoxMax;
            return gsCollectiveRequest<T>::max(m_boxMax.data(), (int)m_boxMax.size());

        case BoxMax:
            chooseAxes();
            if (m_weighted)
            {
                fillWsum();
                m_state = Wsum;
                return gsCollectiveRequest<T>::sum(m_wsum.data(), (int)m_wsum.size());
            }
            m_mode.assign(m_active.size(), ModeExact);
            m_cSlot.clear();
            m_state = Prepare;
            break;

        case Wsum:
            planWeighted();
            m_iter = 0;
            m_state = m_cSlot.empty() ? Prepare : CutIter;
            break;

        case CutIter:
            fillCutWL();
            m_state = CutIterDone;
            return gsCollectiveRequest<T>::sum(m_cWL.data(), (int)m_cWL.size());

        case CutIterDone:
            for (size_t j = 0; j != m_cSlot.size(); ++j)
                if (m_cWL[j] < m_cTarget[j]) m_cLo[j] = m_cMid[j]; else m_cHi[j] = m_cMid[j];
            if (++m_iter != 16)
            {
                m_state = CutIter;
                break;
            }
            fillCutCount();
            m_state = CutCount;
            return gsCollectiveRequest<T>::sum(m_cWL.data(), (int)m_cWL.size());

        case CutCount:
            applyCutCounts();
            m_state = Prepare;
            break;

        case Prepare:
            if (prepareSelect())
            {
                m_state = Select;
                break;
            }
            m_state = AfterSelect;
            break;

        case Select:
            if (m_sel.step(r)) return r;
            m_state = AfterSelect;
            break;

        case AfterSelect:
            splitNodes();
            ++m_levels;
            m_state = LevelBegin;
            break;

        case Exchange:
            if (m_ex.step(r)) return r;
            m_state = Finished;
            break;

        case Finished:
            return gsCollectiveRequest<T>::done();
        }
    }

    /// True once step() has returned DONE.
    bool finished() const { return m_state == Finished; }

    /// Ids of the rank's part, sorted ascending (unique when the input ids are distinct across ranks); valid once finished.
    const std::vector<index_t> & ownElements() const { return m_ex.ownIds(); }

    /// Moves the own ids out (valid once finished); ownElements() is empty afterwards.
    std::vector<index_t> releaseOwnElements()
    {
        std::vector<index_t> r;
        r.swap(m_ex.ownIds());
        return r;
    }

    /// Label (== destination rank) of every local element, nLocal entries; valid once finished.
    const std::vector<index_t> & localLabels() const { return m_labels; }

    /// Number of bisection levels run (<= ceil(log2 nparts)).
    index_t levels() const { return m_levels; }

private:
    enum State { Start, Validity, LevelBegin, BoxMin, BoxMax, Wsum, CutIter, CutIterDone,
                 CutCount, Prepare, Select, AfterSelect, Exchange, Finished };
    enum Mode  { ModeExact, ModeCut, ModeAllLeft, ModeSelect };

    /// Node of the bisection tree: k parts labelled label0 .. label0+k-1 share n global elements;
    /// split is set once the node has been bisected.
    struct Node { index_t label0; index_t k; int64_t n; bool split; };

    T coord(short_t d, index_t i) const { return m_cdata[(size_t)i * (size_t)m_gd + (size_t)d]; }

    static bool isFinite(const T & x) { return (x - x) == T(0); }

    index_t idOf(index_t i) const { return m_ids ? m_ids[i] : m_idLo + i; }

    static T canon(T v) { if (v == T(0)) v = T(0); return v; }

    /// {#columns with non-finite coordinates, #ids out of range, nLocal}
    void fillValidity()
    {
        int64_t nonfinite = 0, bad = 0;
        for (index_t i = 0; i != m_nLocal; ++i)
        {
            bool nf = false;
            for (short_t d = 0; d != m_gd; ++d)
                if (!isFinite(coord(d, i))) nf = true;
            nonfinite += nf;
            const int64_t id = idOf(i);
            if (id < 0 || id >= m_numGlobal) ++bad;
        }
        m_val.assign(3, 0);
        m_val[0] = nonfinite; m_val[1] = bad; m_val[2] = (int64_t)m_nLocal;
    }

    /// Active nodes (k >= 2, n >= 1) in node order and the node -> slot table. O(nodes).
    void collectActive()
    {
        m_active.clear();
        m_slotOf.assign(m_nodes.size(), -1);
        for (size_t v = 0; v != m_nodes.size(); ++v)
            if (!m_nodes[v].split && m_nodes[v].k >= 2 && m_nodes[v].n >= 1)
            {
                m_slotOf[v] = (int)m_active.size();
                m_active.push_back((index_t)v);
            }
    }

    /// Local per-slot coordinate extrema. O(nLocal geoDim).
    void fillBox()
    {
        const size_t A = m_active.size(), gd = (size_t)m_gd;
        m_boxMin.assign(A * gd, std::numeric_limits<T>::max());
        m_boxMax.assign(A * gd, std::numeric_limits<T>::lowest());
        for (index_t i = 0; i != m_nLocal; ++i)
        {
            const int s = m_slotOf[m_node[i]];
            if (s < 0) continue;
            for (short_t d = 0; d != m_gd; ++d)
            {
                const T v = canon(coord(d, i));
                T & lo = m_boxMin[(size_t)s * gd + d];
                T & hi = m_boxMax[(size_t)s * gd + d];
                if (v < lo) lo = v;
                if (v > hi) hi = v;
            }
        }
    }

    /// Longest axis per slot from the reduced box (ties to the lowest axis).
    void chooseAxes()
    {
        const size_t A = m_active.size(), gd = (size_t)m_gd;
        m_axis.assign(A, 0); m_axLo.assign(A, T(0)); m_axHi.assign(A, T(0));
        for (size_t s = 0; s != A; ++s)
        {
            T bestE = -1;
            for (short_t d = 0; d != m_gd; ++d)
            {
                const T lo = m_boxMin[s * gd + d], hi = m_boxMax[s * gd + d];
                if (hi - lo > bestE)
                {
                    bestE = hi - lo; m_axis[s] = d; m_axLo[s] = lo; m_axHi[s] = hi;
                }
            }
        }
    }

    /// Local weight sum per slot (int64). O(nLocal).
    void fillWsum()
    {
        m_wsum.assign(m_active.size(), 0);
        for (index_t i = 0; i != m_nLocal; ++i)
        {
            const int s = m_slotOf[m_node[i]];
            if (s >= 0) m_wsum[s] += (int64_t)m_w[i];
        }
    }

    /// Modes from the reduced node weights: target == 0 -> exact split, else coordinate bisection.
    void planWeighted()
    {
        const size_t A = m_active.size();
        m_mode.assign(A, ModeExact);
        m_cSlot.clear(); m_cLo.clear(); m_cHi.clear(); m_cTarget.clear();
        for (size_t s = 0; s != A; ++s)
        {
            const Node & nd = m_nodes[m_active[s]];
            const index_t kL = (nd.k + 1) / 2;
            const int64_t target = (m_wsum[s] * (int64_t)kL) / (int64_t)nd.k;
            if (target == 0) continue;
            m_mode[s] = ModeCut;
            m_cSlot.push_back((int)s);
            m_cLo.push_back(m_axLo[s]);
            m_cHi.push_back(m_axHi[s]);
            m_cTarget.push_back(target);
        }
        m_cutIdx.assign(A, -1);
        for (size_t j = 0; j != m_cSlot.size(); ++j) m_cutIdx[m_cSlot[j]] = (int)j;
        m_cMid.assign(m_cSlot.size(), T(0));
        m_cCut.assign(m_cSlot.size(), T(0));
        m_cWL.assign(m_cSlot.size(), 0);
    }

    /// Local left weights of the midpoints of one bisection iteration. O(nLocal).
    void fillCutWL()
    {
        const size_t B = m_cSlot.size();
        for (size_t j = 0; j != B; ++j) m_cMid[j] = T(0.5) * (m_cLo[j] + m_cHi[j]);
        m_cWL.assign(B, 0);
        for (index_t i = 0; i != m_nLocal; ++i)
        {
            const int s = m_slotOf[m_node[i]];
            if (s < 0) continue;
            const int j = m_cutIdx[s];
            if (j >= 0 && coord(m_axis[s], i) < m_cMid[j])
                m_cWL[j] += (int64_t)m_w[i];
        }
    }

    /// Final cut values and local left counts. O(nLocal).
    void fillCutCount()
    {
        const size_t B = m_cSlot.size();
        for (size_t j = 0; j != B; ++j) m_cCut[j] = T(0.5) * (m_cLo[j] + m_cHi[j]);
        m_cWL.assign(B, 0);
        for (index_t i = 0; i != m_nLocal; ++i)
        {
            const int s = m_slotOf[m_node[i]];
            if (s < 0) continue;
            const int j = m_cutIdx[s];
            if (j >= 0 && coord(m_axis[s], i) < m_cCut[j])
                ++m_cWL[j];
        }
    }

    /// Reduced left counts: a degenerate cut (all or none left, n >= 2) becomes an exact split.
    void applyCutCounts()
    {
        m_nL.assign(m_active.size(), 0);
        for (size_t j = 0; j != m_cSlot.size(); ++j)
        {
            const size_t s = (size_t)m_cSlot[j];
            const int64_t n = m_nodes[m_active[s]].n, c = m_cWL[j];
            if ((c == 0 || c == n) && n >= 2)
                m_mode[s] = ModeExact;
            else
                m_nL[s] = c;
        }
    }

    /// Exact splits: nL = (n kL + k/2)/k clamped to [0, n]; nL < n becomes a select target.
    /// Returns whether a select has to run (replicated). O(nLocal + nodes).
    bool prepareSelect()
    {
        const size_t A = m_active.size();
        if (m_nL.size() != A) m_nL.assign(A, 0);
        m_tgtOf.assign(A, -1);
        std::vector<index_t> groups;
        std::vector<int64_t> params;
        for (size_t s = 0; s != A; ++s)
        {
            if (m_mode[s] != ModeExact) continue;
            const Node & nd = m_nodes[m_active[s]];
            const index_t kL = (nd.k + 1) / 2;
            int64_t nL = (nd.n * (int64_t)kL + (int64_t)(nd.k / 2)) / (int64_t)nd.k;
            if (nL < 0) nL = 0;
            if (nL > nd.n) nL = nd.n;
            m_nL[s] = nL;
            if (nL == nd.n)
                m_mode[s] = ModeAllLeft;
            else
            {
                m_mode[s] = ModeSelect;
                m_tgtOf[s] = (int)groups.size();
                groups.push_back(m_active[s]);
                params.push_back(nL);
            }
        }
        if (groups.empty()) return false;

        m_keys.assign((size_t)m_nLocal, gsPartitionKey());
        m_group.assign((size_t)m_nLocal, -1);
        for (index_t i = 0; i != m_nLocal; ++i)
        {
            const int s = m_slotOf[m_node[i]];
            if (s < 0 || m_mode[s] != ModeSelect) continue;
            m_keys[i].hi = gsOrderedKey<T>::encode(coord(m_axis[s], i));
            m_keys[i].lo = (gsPartitionKey::lo_type)idOf(i);
            m_group[i] = m_node[i];
        }
        m_sel.setup(gsRadixSelect<T>::Count, m_nLocal, m_keys.data(), m_group.data(), NULL,
                    groups, params);
        return true;
    }

    /// Appends the children of every active node and moves the local elements. O(nLocal + A).
    void splitNodes()
    {
        const size_t A = m_active.size();
        for (size_t s = 0; s != A; ++s)
            if (m_mode[s] == ModeSelect)
                GISMO_ENSURE(m_sel.found(m_tgtOf[s]), "RCB select found no element of the requested rank.");

        const index_t base = (index_t)m_nodes.size();
        for (size_t s = 0; s != A; ++s)
        {
            const Node nd = m_nodes[m_active[s]];
            m_nodes[m_active[s]].split = true;
            const index_t kL = (nd.k + 1) / 2;
            Node l; l.label0 = nd.label0;      l.k = kL;        l.n = m_nL[s];        l.split = false;
            Node r; r.label0 = nd.label0 + kL; r.k = nd.k - kL; r.n = nd.n - m_nL[s]; r.split = false;
            m_nodes.push_back(l);
            m_nodes.push_back(r);
        }
        for (index_t i = 0; i != m_nLocal; ++i)
        {
            const int s = m_slotOf[m_node[i]];
            if (s < 0) continue;
            bool left = true;
            switch (m_mode[s])
            {
            case ModeCut:
                left = coord(m_axis[s], i) < m_cCut[m_cutIdx[s]];
                break;
            case ModeSelect:
                left = m_keys[i] < m_sel.key(m_tgtOf[s]);
                break;
            default:
                break;
            }
            m_node[i] = base + 2 * (index_t)s + (left ? 0 : 1);
        }
    }

    /// Leaf labels, replicated part sizes and the start of the routing exchange. O(nLocal + nodes).
    void startExchange()
    {
        m_labels.assign((size_t)m_nLocal, 0);
        for (index_t i = 0; i != m_nLocal; ++i)
        {
            GISMO_ASSERT(m_nodes[m_node[i]].k == 1, "Element left in an inner node.");
            m_labels[i] = m_nodes[m_node[i]].label0;
        }
        std::vector<int64_t> partSize((size_t)m_nparts, 0);
        for (size_t v = 0; v != m_nodes.size(); ++v)
            if (m_nodes[v].k == 1)
                partSize[m_nodes[v].label0] += m_nodes[v].n;
        for (size_t p = 0; p != partSize.size(); ++p)
            GISMO_ENSURE(partSize[p] <= (int64_t)std::numeric_limits<int>::max(),
                         "A part exceeds INT_MAX elements.");
        if (m_ids)
            m_ex.setup(m_ids, m_nLocal, m_labels.data(), m_nparts, m_numGlobal, partSize[m_rank]);
        else
            m_ex.setup(m_idLo, m_nLocal, m_labels.data(), m_nparts, m_numGlobal, partSize[m_rank]);
    }

    State   m_state;
    const T * m_cdata;          // [geoDim x nLocal], column-major, caller-owned
    short_t m_gd;
    index_t m_nLocal;
    index_t m_idLo;
    const index_t * m_ids;      // [nLocal] or null, caller-owned
    const index_t * m_w;        // [nLocal] or null, caller-owned
    bool    m_weighted;
    index_t m_nparts, m_rank;
    int64_t m_numGlobal;
    index_t m_levels;
    int     m_iter;

    std::vector<int64_t> m_val;       // [3] validity sums
    std::vector<Node>    m_nodes;     // replicated node table
    std::vector<index_t> m_node;      // [nLocal] node of each local element
    std::vector<index_t> m_labels;    // [nLocal] leaf label of each local element

    std::vector<index_t> m_active;    // [A] active node indices
    std::vector<int>     m_slotOf;    // [nodes] slot of a node or -1
    std::vector<T>       m_boxMin, m_boxMax; // [A*geoDim], index slot*geoDim + d
    std::vector<short_t> m_axis;      // [A]
    std::vector<T>       m_axLo, m_axHi; // [A]
    std::vector<int64_t> m_wsum;      // [A]
    std::vector<char>    m_mode;      // [A] Mode
    std::vector<int64_t> m_nL;        // [A] left count
    std::vector<int>     m_tgtOf;     // [A] select target of a slot or -1

    std::vector<int>     m_cSlot;     // [B] slot of cut node j
    std::vector<int>     m_cutIdx;    // [A] cut index of a slot or -1
    std::vector<T>       m_cLo, m_cHi, m_cMid, m_cCut; // [B]
    std::vector<int64_t> m_cTarget;   // [B]
    std::vector<int64_t> m_cWL;       // [B] left weights, then left counts

    std::vector<gsPartitionKey> m_keys; // [nLocal]
    std::vector<index_t>        m_group; // [nLocal]
    gsRadixSelect<T>  m_sel;
    gsOwnIdsExchange<T> m_ex;
};

/** @brief Distributed space-filling-curve labelling: rank r ends with the
    sorted ids e with gsGeometricPartitioner<T>::curveLabels(keys, w, P)[e] ==
    r, without any rank holding O(N) data (collective request protocol, \ref
    gsDistributedPartitionProtocol).

    <b>Inputs.</b> This rank's contiguous id slice [idLo, idLo + nLocal), the
    curve key of every local element (computed by the caller exactly as the
    serial labelling does, so the machine needs no geometry), the element
    weights (index_t, >= 0), the number of ranks P (== nparts, part p is
    owned by rank p), this rank and the global element count N. Element i of
    the slice has the full key {hi = keys[i], lo = idLo + i}, the (key, id)
    total order of curveLabels.

    <b>Requests</b> (identical on all ranks; every branch depends on P, N or
    reduced data):
    -# SUM_I64 (len 3) {sum of weights, #negative weights, nLocal}. The reduced
       values are checked on every rank: a negative weight or slices that do
       not cover N make all ranks throw (GISMO_ENSURE) from the same step.
       W denotes the reduced weight sum.
    -# Skipped if P == 1 or N == 0. Otherwise a gsRadixSelect (one group)
       with the P-1 targets tau_j = ceil(W j / P), j = 1..P-1, computed
       as q j + (r j + P - 1) / P with q = W / P, r = W % P (exact int64, W j
       is never formed), at most gsRadixSelect::maxPasses SUM_I64. If W > 0
       it is Weighted with thresholds tau_j. If W == 0 every element weighs
       1 instead (W := N, the count split of rcbBisect) and it is Count with
       0-based ranks tau_j - 1: with unit weights the smallest key whose
       inclusive prefix reaches tau is the tau-th smallest. The selected keys
       f_1 <= ... <= f_{P-1} are the cuts.
    -# SUM_I64 (len P): the per-label element counts, i.e. the replicated
       part sizes. They are checked against INT_MAX and against N on every
       rank, and the own size is handed to the exchange.
    -# gsOwnIdsExchange of the labels: ALLTOALL_INT, ALLTOALLV_INDEX, SUM_I64
       (len 1, coverage).
    The total is at most 1 + maxPasses + 1 + 3 collectives (5 if the
    select is skipped).

    <b>Correctness.</b> Let A(e) = W(< e) be the sum of the weights of the
    elements with a full key smaller than that of e (exclusive prefix, in
    (key, id) order). It is the value of `acc` when the oracle loop of
    curveLabels reaches e. The loop only increases j, and its test
    P A >= W (j+1) is monotone in j (W >= 0), so it leaves
    @f[ \mathrm{label}(e) = \#\{ j \in [1,P-1] : P\,A(e) \ge W j \}. @f]
    - W == 0 (all weights zero): both this machine and the oracle use unit
      weights, A(e) = position of e in (key, id) order and W = N >= 1, so
      the argument below applies unchanged.
    - W > 0: A(e) is an integer, so P A(e) >= W j <=> A(e) >= ceil(W j / P)
      = tau_j. The identity ceil(W j / P) = q j + ceil(r j / P) holds for W =
      q P + r, and r j < P^2 keeps the arithmetic in int64.
    - Lemma. For 1 <= tau <= W and f(tau) = min{ f : W(<= f) >= tau } (the
      smallest key whose inclusive prefix weight reaches tau, which exists
      since tau <= W): A(e) >= tau <=> key(e) > f(tau).
      (<=) If key(e) > f, then W(< e) >= W(<= f) >= tau.
      (=>) If key(e) <= f, then W(< e) <= W(< f). Either f is the smallest
      key, so W(< f) = 0 < tau; or f has a predecessor f', and
      W(< f) = W(<= f') < tau by minimality of f.
    - Hence label(e) = #{ j : f_j < key(e) } = lower_bound(f, key(e)). It is
      lower_bound and not upper_bound: the element that IS f_j has
      W(< e) < tau_j and must not be counted for j.
    - Zero weights do not move A, so every element after f_j counts j
      whatever its weight. A heavy element h makes several consecutive f_j
      equal to key(h); the label jumps by several at the element after h and
      the parts in between are empty, as in the oracle. Duplicate curve keys
      are ordered by id in both.

    <b>Complexity.</b> Local work O(nLocal) for the totals and the keys,
    O(nLocal log P) per select pass (at most maxPasses passes), O(nLocal log
    P) for the labels, O(nLocal + m log m) for the exchange (m = own ids).
    Collectives: at most 1 + maxPasses + 1 + 3; the largest request is the
    select's, at most (P-1) 512 + 2 (P-1) int64 values. Memory: O(nLocal +
    m) + O(P 256); nothing is sized by the global number of elements.

    The machine holds pointers into its request buffers once running: do not
    copy or move it after the first step() (before it, copying and moving
    are fine). T only fixes the request type; no T-valued collective is
    issued, so the class works for every coefficient type.
*/
template<class T>
class gsDistributedCurve
{
public:
    typedef T Scalar;

    /** @brief Constructs the machine for one rank.

        \param nLocal  size of this rank's slice; local element i has the
               global id idLo + i
        \param idLo    first id of the slice
        \param keys    curve key of local element i, nLocal entries (may be
               null if nLocal == 0); copied
        \param weights weight of local element i (>= 0), nLocal entries (may
               be null if nLocal == 0); copied
        \param nparts  number of parts == number of ranks P (replicated)
        \param rank    this rank, in [0, nparts)
        \param numGlobal N (replicated)
    */
    gsDistributedCurve(index_t nLocal, index_t idLo, const uint64_t * keys, const index_t * weights,
                       index_t nparts, index_t rank, int64_t numGlobal)
    : m_state(Start), m_nLocal(nLocal), m_idLo(idLo), m_nparts(nparts), m_rank(rank),
      m_numGlobal(numGlobal), m_issued(0), m_totalW(0)
    {
        GISMO_ENSURE(nparts >= 1 && rank >= 0 && rank < nparts, "Invalid number of parts or rank.");
        GISMO_ENSURE(numGlobal >= 0, "Negative global element count.");
        GISMO_ENSURE(nLocal >= 0 && idLo >= 0, "Invalid slice.");
        GISMO_ENSURE(nLocal == 0 || (keys != NULL && weights != NULL), "Missing keys or weights.");
        m_keys.resize((size_t)nLocal);
        m_w64.resize((size_t)nLocal);
        for (index_t i = 0; i != nLocal; ++i)
        {
            m_keys[i].hi = keys[i];
            m_keys[i].lo = static_cast<gsPartitionKey::lo_type>(idLo + i);
            m_w64[i]     = static_cast<int64_t>(weights[i]);
        }
    }

    /// Machine concept, see \ref gsDistributedPartitionProtocol.
    gsCollectiveRequest<T> step()
    {
        gsCollectiveRequest<T> r;
        for (;;) switch (m_state)
        {
        case Start:
            m_tot.assign(3, 0);
            for (index_t i = 0; i != m_nLocal; ++i)
            {
                m_tot[0] += m_w64[i];
                m_tot[1] += (m_w64[i] < 0);
            }
            m_tot[2] = m_nLocal;
            m_state = Totals;
            return issue(gsCollectiveRequest<T>::sum(m_tot.data(), 3));

        case Totals:
            GISMO_ENSURE(m_tot[1] == 0, "Negative element weight.");
            GISMO_ENSURE(m_tot[2] == m_numGlobal, "The slices do not cover numGlobal elements.");
            m_totalW = m_tot[0];
            if (m_nparts > 1 && m_numGlobal > 0)
            {
                const bool unit = (m_totalW == 0);
                const int64_t W = unit ? m_numGlobal : m_totalW;
                const int64_t P = m_nparts, q = W / P, rm = W % P;
                std::vector<int64_t> tau((size_t)(P - 1));
                for (int64_t j = 1; j != P; ++j)
                    tau[j - 1] = q * j + (rm * j + P - 1) / P - (unit ? 1 : 0);
                m_sel.setup(unit ? gsRadixSelect<T>::Count : gsRadixSelect<T>::Weighted,
                            m_nLocal, m_keys.data(), NULL, unit ? NULL : m_w64.data(),
                            std::vector<index_t>((size_t)(P - 1), 0), give(tau));
                m_state = Select;
            }
            else
                m_state = Labels;
            break;

        case Select:
            if (m_sel.step(r))
                return issue(r);
            m_cuts.clear();
            for (index_t j = 1; j < m_nparts; ++j)
            {
                GISMO_ENSURE(m_sel.found(j - 1), "Curve cut not found.");
                m_cuts.push_back(m_sel.key(j - 1));
            }
            m_state = Labels;
            break;

        case Labels:
            computeLabels();
            m_state = Parts;
            return issue(gsCollectiveRequest<T>::sum(m_part.data(), (int)m_nparts));

        case Parts:
        {
            int64_t tot = 0;
            for (size_t p = 0; p != m_part.size(); ++p)
            {
                GISMO_ENSURE(m_part[p] <= (int64_t)std::numeric_limits<int>::max(),
                             "A part exceeds INT_MAX elements.");
                tot += m_part[p];
            }
            GISMO_ENSURE(tot == m_numGlobal, "The part sizes do not sum to numGlobal.");
            m_ex.setup(m_idLo, m_nLocal, m_localLabels.data(), m_nparts, m_numGlobal, m_part[m_rank]);
            m_state = Exchange;
            break;
        }

        case Exchange:
            if (m_ex.step(r))
                return issue(r);
            m_own.swap(m_ex.ownIds());
            m_state = Finished;
            break;

        case Finished:
            return gsCollectiveRequest<T>::done();
        }
    }

    /// True once step() has returned DONE.
    bool finished() const { return m_state == Finished; }

    /// Own element ids, sorted ascending and unique; valid once finished.
    const std::vector<index_t> & ownElements() const { return m_own; }

    /// Moves the own ids out; valid once finished, leaves ownElements() empty.
    std::vector<index_t> releaseOwnElements()
    {
        std::vector<index_t> out;
        out.swap(m_own);
        return out;
    }

    /// Label (== destination rank) of local element i, nLocal entries; valid once finished.
    const std::vector<index_t> & localLabels() const { return m_localLabels; }

    /// The cuts f_1 <= ... <= f_{P-1}; empty if P == 1 or N == 0. Valid once finished.
    const std::vector<gsPartitionKey> & cuts() const { return m_cuts; }

    /// Reduced total weight W; valid after the first reduction.
    int64_t totalWeight() const { return m_totalW; }

    /// Number of collective requests issued so far (DONE not counted).
    index_t requestsIssued() const { return m_issued; }

private:
    enum State { Start, Totals, Select, Labels, Parts, Exchange, Finished };

    gsCollectiveRequest<T> issue(const gsCollectiveRequest<T> & r)
    {
        ++m_issued;
        return r;
    }

    /// Labels of the local elements and the local per-label counts (m_part). O(nLocal log P + P).
    void computeLabels()
    {
        m_localLabels.assign((size_t)m_nLocal, 0);
        m_part.assign((size_t)m_nparts, 0);
        for (index_t i = 0; i != m_nLocal; ++i)
        {
            index_t l = 0;
            if (m_nparts > 1)
                l = (index_t)(std::lower_bound(m_cuts.begin(), m_cuts.end(), m_keys[i]) - m_cuts.begin());
            m_localLabels[i] = l;
            ++m_part[l];
        }
    }

    State   m_state;
    index_t m_nLocal;
    index_t m_idLo;
    index_t m_nparts, m_rank;
    int64_t m_numGlobal;
    index_t m_issued;
    int64_t m_totalW;

    std::vector<gsPartitionKey> m_keys;   // [nLocal] full key {curve key, id} of local element i
    std::vector<int64_t>        m_w64;    // [nLocal] weights
    std::vector<int64_t>        m_tot;    // [3] {sum of weights, #negative weights, nLocal}, then reduced
    std::vector<gsPartitionKey> m_cuts;   // [P-1] f_1 <= ... <= f_{P-1}, or empty
    std::vector<index_t>        m_localLabels; // [nLocal] label of local element i
    std::vector<int64_t>        m_part;   // [P] local per-label counts, then the reduced part sizes
    std::vector<index_t>        m_own;    // [own count] own ids, sorted ascending

    gsRadixSelect<T>    m_sel;
    gsOwnIdsExchange<T> m_ex;
};

} // namespace gismo
