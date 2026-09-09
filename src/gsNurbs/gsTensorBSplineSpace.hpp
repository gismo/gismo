/** @file gsTensorBSplineSpace.hpp

    @brief Implementation of tensor B-spline space descriptions.
*/

#pragma once

#include <gsNurbs/gsTensorBSplineSpace.h>
#include <gsNurbs/gsTensorBSpline.h>

#include <algorithm>
#include <vector>

namespace gismo
{

namespace internal
{

template <class T>
std::vector<T> uniqueKnotValues(const gsKnotVector<T>& kv)
{
    std::vector<T> values;
    values.reserve(static_cast<size_t>(kv.uSize()));
    for (index_t i = 0; i < static_cast<index_t>(kv.uSize()); ++i)
        values.push_back(kv.uValue(i));
    return values;
}

template <class T>
void appendUniqueKnotValues(std::vector<T>& values, const gsKnotVector<T>& kv)
{
    for (index_t i = 0; i < static_cast<index_t>(kv.uSize()); ++i)
    {
        const T xi = kv.uValue(i);
        if (std::find(values.begin(), values.end(), xi) == values.end())
            values.push_back(xi);
    }
}

template <class T>
std::vector<T> mergedUniqueKnotValues(const gsKnotVector<T>& a,
                                      const gsKnotVector<T>& b)
{
    std::vector<T> values = uniqueKnotValues(a);
    appendUniqueKnotValues(values, b);
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
    return values;
}

template <class T>
int multiplicityOrZero(const gsKnotVector<T>& kv, T xi)
{
    return kv.has(xi) ? kv.multiplicity(xi) : 0;
}

template <class T>
int productContinuity(const gsKnotVector<T>& kv, T xi, short_t productDegree)
{
    // If xi is not a breakpoint of this factor, the factor is polynomial
    // across xi and does not limit the continuity of the product there.
    return kv.has(xi)
         ? static_cast<int>(kv.degree()) - kv.multiplicity(xi)
         : static_cast<int>(productDegree);
}

template <class T>
gsKnotVector<T> knotVectorFromMultiplicities(const std::vector<T>& values,
                                             const std::vector<int>& mult,
                                             short_t degree)
{
    GISMO_ASSERT(values.size() == mult.size(),
                 "knotVectorFromMultiplicities: size mismatch.");
    std::vector<T> knots;
    size_t total = 0;
    for (size_t i = 0; i < mult.size(); ++i)
    {
        GISMO_ASSERT(mult[i] > 0,
                     "knotVectorFromMultiplicities: non-positive multiplicity.");
        total += static_cast<size_t>(mult[i]);
    }
    knots.reserve(total);
    for (size_t i = 0; i < values.size(); ++i)
        for (int k = 0; k < mult[i]; ++k)
            knots.push_back(values[i]);
    return gsKnotVector<T>(knots, degree);
}

template <class T>
gsKnotVector<T> degreeElevatedKnotVector(gsKnotVector<T> kv, short_t targetDegree)
{
    GISMO_ENSURE(targetDegree >= kv.degree(),
                 "Cannot lower degree while building a common tensor spline space.");
    const short_t delta = static_cast<short_t>(targetDegree - kv.degree());
    if (delta > 0)
        kv.degreeElevate(delta);
    return kv;
}

template <class T>
bool sameDomain(const gsKnotVector<T>& a, const gsKnotVector<T>& b)
{
    return a.first() == b.first() && a.last() == b.last();
}

template <class T>
bool sameBreaks(const gsKnotVector<T>& a, const gsKnotVector<T>& b)
{
    if (!sameDomain(a, b) || a.uSize() != b.uSize())
        return false;
    for (index_t i = 0; i < static_cast<index_t>(a.uSize()); ++i)
        if (a.uValue(i) != b.uValue(i))
            return false;
    return true;
}

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> elevatedSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec,
    short_t dir,
    short_t delta)
{
    std::vector< gsKnotVector<T> > kvs;
    kvs.reserve(static_cast<size_t>(d));
    for (short_t k = 0; k < d; ++k)
    {
        gsKnotVector<T> kv = spec.knots(k);
        if (k == dir && delta > 0)
            kv.degreeElevate(delta);
        kvs.push_back(kv);
    }
    return gsTensorBSplineSpaceSpec<d,T>(kvs);
}

} // namespace internal

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T>::gsTensorBSplineSpaceSpec()
{
}

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T>::gsTensorBSplineSpaceSpec(
    const std::vector<KnotVectorType>& knots)
    : m_knots(knots)
{
    GISMO_ENSURE(m_knots.size() == static_cast<size_t>(d),
                 "gsTensorBSplineSpaceSpec: expected " << d << " knot vectors, got "
                 << m_knots.size() << ".");
}

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T>
gsTensorBSplineSpaceSpec<d,T>::fromBasis(const Basis& basis)
{
    GISMO_ENSURE(!basis.isPeriodic(),
                 "gsTensorBSplineSpaceSpec currently supports non-periodic bases only.");
    std::vector<KnotVectorType> kvs;
    kvs.reserve(static_cast<size_t>(d));
    for (short_t k = 0; k < d; ++k)
        kvs.push_back(basis.knots(k));
    return gsTensorBSplineSpaceSpec(kvs);
}

template <short_t d, class T>
typename gsTensorBSplineSpaceSpec<d,T>::Basis
gsTensorBSplineSpaceSpec<d,T>::toBasis() const
{
    GISMO_ENSURE(this->isValid(), "gsTensorBSplineSpaceSpec is empty.");
    return Basis(m_knots);
}

template <short_t d, class T>
const typename gsTensorBSplineSpaceSpec<d,T>::KnotVectorType&
gsTensorBSplineSpaceSpec<d,T>::knots(short_t dir) const
{
    GISMO_ASSERT(dir >= 0 && dir < d, "Invalid tensor direction.");
    return m_knots[static_cast<size_t>(dir)];
}

template <short_t d, class T>
short_t gsTensorBSplineSpaceSpec<d,T>::degree(short_t dir) const
{
    return this->knots(dir).degree();
}

template <short_t d, class T>
index_t gsTensorBSplineSpaceSpec<d,T>::size() const
{
    return this->toBasis().size();
}

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> bezierSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec)
{
    std::vector< gsKnotVector<T> > kvs;
    kvs.reserve(static_cast<size_t>(d));
    for (short_t dir = 0; dir < d; ++dir)
    {
        const gsKnotVector<T>& src = spec.knots(dir);
        const std::vector<T> values = internal::uniqueKnotValues(src);
        std::vector<int> mult(values.size());
        const int bezMult = static_cast<int>(src.degree()) + 1;
        for (size_t i = 0; i < values.size(); ++i)
        {
            const bool boundary = (i == 0 || i + 1 == values.size());
            mult[i] = boundary ? src.multiplicity(values[i]) : bezMult;
        }
        kvs.push_back(internal::knotVectorFromMultiplicities(values, mult, src.degree()));
    }
    return gsTensorBSplineSpaceSpec<d,T>(kvs);
}

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> productSpace(
    const gsTensorBSplineSpaceSpec<d,T>& a,
    const gsTensorBSplineSpaceSpec<d,T>& b,
    gsTensorBSplineSpacePolicy policy)
{
    std::vector< gsKnotVector<T> > kvs;
    kvs.reserve(static_cast<size_t>(d));
    for (short_t dir = 0; dir < d; ++dir)
    {
        const gsKnotVector<T>& akv = a.knots(dir);
        const gsKnotVector<T>& bkv = b.knots(dir);
        GISMO_ENSURE(internal::sameDomain(akv, bkv),
                     "productSpace requires matching tensor-space domains.");

        const short_t p = akv.degree();
        const short_t q = bkv.degree();
        const short_t prodDegree = static_cast<short_t>(p + q);
        const std::vector<T> values = internal::mergedUniqueKnotValues(akv, bkv);
        std::vector<int> mult(values.size());

        for (size_t i = 0; i < values.size(); ++i)
        {
            const bool boundary = (i == 0 || i + 1 == values.size());
            if (boundary || policy == gsTensorBSplineSpacePolicy::Bezier)
            {
                mult[i] = static_cast<int>(prodDegree) + 1;
            }
            else
            {
                const int ca = internal::productContinuity(akv, values[i], prodDegree);
                const int cb = internal::productContinuity(bkv, values[i], prodDegree);
                mult[i] = static_cast<int>(prodDegree) - std::min(ca, cb);
                mult[i] = std::max(1, std::min(mult[i], static_cast<int>(prodDegree) + 1));
            }
        }
        kvs.push_back(internal::knotVectorFromMultiplicities(values, mult, prodDegree));
    }
    return gsTensorBSplineSpaceSpec<d,T>(kvs);
}

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> powerSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec,
    index_t power,
    gsTensorBSplineSpacePolicy policy)
{
    GISMO_ENSURE(power >= 1, "powerSpace expects power >= 1.");
    if (power == 1)
        return spec;

    gsTensorBSplineSpaceSpec<d,T> result = spec;
    for (index_t k = 1; k < power; ++k)
        result = productSpace(result, spec, policy);
    return result;
}

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> derivativeSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec,
    const gsVector<index_t,d>& orders)
{
    std::vector< gsKnotVector<T> > kvs;
    kvs.reserve(static_cast<size_t>(d));
    for (short_t dir = 0; dir < d; ++dir)
    {
        const gsKnotVector<T>& src = spec.knots(dir);
        const index_t order = orders[dir];
        const index_t sz = static_cast<index_t>(src.size());
        GISMO_ENSURE(order <= src.degree(),
                     "derivativeSpace order exceeds polynomial degree.");
        GISMO_ENSURE(sz > 2 * order,
                     "derivativeSpace cannot drop enough boundary knots.");

        std::vector<T> knots;
        knots.reserve(static_cast<size_t>(sz - 2 * order));
        for (index_t i = order; i < sz - order; ++i)
            knots.push_back(src[i]);
        kvs.push_back(gsKnotVector<T>(knots,
            static_cast<short_t>(src.degree() - order)));
    }
    return gsTensorBSplineSpaceSpec<d,T>(kvs);
}

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> divergenceSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec)
{
    std::vector< gsTensorBSplineSpaceSpec<d,T> > terms;
    terms.reserve(static_cast<size_t>(d));
    for (short_t dir = 0; dir < d; ++dir)
    {
        gsVector<index_t,d> orders;
        orders.setZero();
        orders[dir] = 1;
        terms.push_back(internal::elevatedSpace(
            derivativeSpace(spec, orders), dir, 1));
    }
    return commonSpace(terms);
}

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> laplacianSpace(
    const gsTensorBSplineSpaceSpec<d,T>& spec,
    gsTensorBSplineSpacePolicy policy)
{
    std::vector< gsTensorBSplineSpaceSpec<d,T> > terms;
    terms.reserve(static_cast<size_t>(d));
    for (short_t dir = 0; dir < d; ++dir)
    {
        gsVector<index_t,d> orders;
        orders.setZero();
        orders[dir] = 2;
        terms.push_back(internal::elevatedSpace(
            derivativeSpace(spec, orders), dir, 2));
    }
    gsTensorBSplineSpaceSpec<d,T> result = commonSpace(terms);
    return policy == gsTensorBSplineSpacePolicy::Bezier
         ? bezierSpace(result)
         : result;
}

template <short_t d, class T>
gsTensorBSplineSpaceSpec<d,T> commonSpace(
    const std::vector< gsTensorBSplineSpaceSpec<d,T> >& specs)
{
    GISMO_ENSURE(!specs.empty(), "commonSpace expects at least one input space.");

    std::vector< gsKnotVector<T> > kvs;
    kvs.reserve(static_cast<size_t>(d));
    for (short_t dir = 0; dir < d; ++dir)
    {
        short_t targetDegree = specs.front().degree(dir);
        for (size_t s = 1; s < specs.size(); ++s)
            targetDegree = std::max(targetDegree, specs[s].degree(dir));

        std::vector< gsKnotVector<T> > elevated;
        elevated.reserve(specs.size());
        std::vector<T> values;
        for (size_t s = 0; s < specs.size(); ++s)
        {
            gsKnotVector<T> kv =
                internal::degreeElevatedKnotVector(specs[s].knots(dir), targetDegree);
            if (s == 0)
                values = internal::uniqueKnotValues(kv);
            else
            {
                GISMO_ENSURE(internal::sameDomain(elevated.front(), kv),
                             "commonSpace requires matching tensor-space domains.");
                internal::appendUniqueKnotValues(values, kv);
            }
            elevated.push_back(kv);
        }

        std::sort(values.begin(), values.end());
        values.erase(std::unique(values.begin(), values.end()), values.end());

        std::vector<int> mult(values.size(), 0);
        for (size_t i = 0; i < values.size(); ++i)
            for (size_t s = 0; s < elevated.size(); ++s)
                mult[i] = std::max(mult[i],
                    internal::multiplicityOrZero(elevated[s], values[i]));

        kvs.push_back(internal::knotVectorFromMultiplicities(values, mult, targetDegree));
    }
    return gsTensorBSplineSpaceSpec<d,T>(kvs);
}

template <short_t d, class T>
gsSparseMatrix<T,RowMajor> transferMatrix(
    const gsTensorBSplineSpaceSpec<d,T>& source,
    const gsTensorBSplineSpaceSpec<d,T>& target)
{
    typedef gsTensorBSplineBasis<d,T> Basis;

    for (short_t dir = 0; dir < d; ++dir)
    {
        const short_t delta = static_cast<short_t>(target.degree(dir) - source.degree(dir));
        GISMO_ENSURE(delta >= 0,
                     "transferMatrix target degree is lower than source degree.");

        const gsKnotVector<T> elevatedSource =
            internal::degreeElevatedKnotVector(source.knots(dir), target.degree(dir));
        const gsKnotVector<T>& kvT = target.knots(dir);
        GISMO_ENSURE(internal::sameDomain(elevatedSource, kvT),
                     "transferMatrix target domain differs from source domain.");

        for (index_t i = 0; i < static_cast<index_t>(elevatedSource.uSize()); ++i)
        {
            const T xi = elevatedSource.uValue(i);
            GISMO_ENSURE(kvT.has(xi),
                         "transferMatrix target is missing a source breakpoint.");
            GISMO_ENSURE(kvT.multiplicity(xi) >= elevatedSource.multiplicity(xi),
                         "transferMatrix target has lower source multiplicity.");
        }
    }

    Basis srcBasis = source.toBasis();
    const Basis tgtBasis = target.toBasis();
    const index_t nSource = srcBasis.size();

    gsMatrix<T> coefs = gsMatrix<T>::Identity(nSource, nSource);
    gsTensorBSpline<d,T> lifted(srcBasis, coefs);

    for (short_t dir = 0; dir < d; ++dir)
    {
        const short_t delta = static_cast<short_t>(target.degree(dir) - lifted.degree(dir));
        GISMO_ENSURE(delta >= 0,
                     "transferMatrix target degree is lower than source degree.");
        if (delta > 0)
            lifted.degreeElevate(delta, dir);
    }

    for (short_t dir = 0; dir < d; ++dir)
    {
        const gsKnotVector<T>& kvT = target.knots(dir);
        const T first = kvT.first();
        const T last  = kvT.last();
        for (index_t i = 0; i < static_cast<index_t>(kvT.uSize()); ++i)
        {
            const T xi = kvT.uValue(i);
            if (xi <= first || xi >= last)
                continue;
            const int have = lifted.knots(dir).has(xi)
                ? lifted.knots(dir).multiplicity(xi) : 0;
            const int need = kvT.multiplicity(xi) - have;
            GISMO_ENSURE(need >= 0,
                         "transferMatrix target is not a superspace of source.");
            if (need > 0)
                lifted.insertKnot(xi, dir, need);
        }
    }

    for (short_t dir = 0; dir < d; ++dir)
    {
        GISMO_ENSURE(lifted.degree(dir) == target.degree(dir),
                     "transferMatrix degree mismatch after transfer construction.");
        const gsKnotVector<T>& have = lifted.knots(dir);
        const gsKnotVector<T>& want = target.knots(dir);
        GISMO_ENSURE(have.uSize() == want.uSize(),
                     "transferMatrix knot mismatch after transfer construction.");
        for (index_t i = 0; i < static_cast<index_t>(want.uSize()); ++i)
        {
            const T xi = want.uValue(i);
            GISMO_ENSURE(have.uValue(i) == xi &&
                         have.multiplicity(xi) == want.multiplicity(xi),
                         "transferMatrix multiplicity mismatch after transfer construction.");
        }
    }

    GISMO_ENSURE(lifted.coefs().rows() == tgtBasis.size(),
                 "transferMatrix output size mismatch.");
    return lifted.coefs().sparseView();
}

} // namespace gismo
