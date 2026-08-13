/** @file gsSpaceFillingCurve.h

    @brief Space-filling curves (Morton / Hilbert) for geometric partitioning

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>
#include <gsCore/gsMath.h>

#include <cstdint>
#include <vector>

namespace gismo
{

namespace internal
{

/// \brief Hilbert transpose -> axes (Skilling, "Programming the Hilbert
/// curve", AIP Conf. Proc. 707 (2004) 381). Transforms \a X in place.
/// \a X has \a n entries of \a b bits each. Not meant for n<=1.
inline void gsHilbertTransposeToAxes(uint64_t * X, int b, int n)
{
    uint64_t N = 2ULL << (b-1), P, Q, t;
    int i;
    // Gray decode by H ^ (H/2)
    t = X[n-1] >> 1;
    for (i = n-1; i > 0; i--) X[i] ^= X[i-1];
    X[0] ^= t;
    // Undo excess work
    for (Q = 2; Q != N; Q <<= 1)
    {
        P = Q - 1;
        for (i = n-1; i >= 0; i--)
            if (X[i] & Q) X[0] ^= P;                              // invert
            else { t = (X[0]^X[i]) & P;  X[0] ^= t;  X[i] ^= t; }  // exchange
    }
}

/// \brief Hilbert axes -> transpose (Skilling, "Programming the Hilbert
/// curve", AIP Conf. Proc. 707 (2004) 381). Transforms \a X in place.
/// \a X has \a n entries of \a b bits each. Not meant for n<=1.
inline void gsHilbertAxesToTranspose(uint64_t * X, int b, int n)
{
    uint64_t M = 1ULL << (b-1), P, Q, t;
    int i;
    // Inverse undo
    for (Q = M; Q > 1; Q >>= 1)
    {
        P = Q - 1;
        for (i = 0; i < n; i++)
            if (X[i] & Q) X[0] ^= P;                              // invert
            else { t = (X[0]^X[i]) & P;  X[0] ^= t;  X[i] ^= t; }  // exchange
    }
    // Gray encode
    for (i = 1; i < n; i++) X[i] ^= X[i-1];
    t = 0;
    for (Q = M; Q > 1; Q >>= 1)
        if (X[n-1] & Q) t ^= Q-1;
    for (i = 0; i < n; i++) X[i] ^= t;
}

} // namespace internal

/** \brief Morton (Z-order) and Hilbert space-filling curves for geometric
 *  partitioning.
 *
 * The static part is a purely integer bit-interleaving (Morton) resp.
 * Skilling axes-to-transpose (Hilbert) encoder/decoder between \a dim
 * lattice coordinates of \a bits bits each and a single uint64_t key, for
 * \a dim in {0,1,2,3}.
 *
 * The non-static (geometric) part maps a physical point inside a bounding
 * box to such a key: it quantises each coordinate to 2^bits levels and
 * drops degenerate axes (an axis whose extent is negligible relative to
 * the largest extent of the box is not encoded at all, so e.g. a surface
 * embedded in 3D encodes as a 2-D curve).
 *
 * \ingroup Utils
 */
class gsSpaceFillingCurve
{
public:

    /// Space-filling curve type
    enum Curve { Morton = 0, Hilbert = 1 };

    // ---------------- static, purely integer part ----------------------------

    /// Bits per axis for a curve of dimension \a dim: 0->0, 1->63, 2->31, 3->21.
    /// The geoDim-1 resolution (63 bits) is in practice limited by the
    /// real_t mantissa (about 53 bits for double); that is accepted, not a
    /// defect, since quantize() only ever produces values representable by
    /// the mantissa in the first place.
    static unsigned bitsPerAxis(short_t dim)
    {
        GISMO_ENSURE(dim >= 0 && dim <= 3, "gsSpaceFillingCurve supports dim in {0,1,2,3}, got "<<dim);
        static const unsigned budget[4] = {0u, 63u, 31u, 21u};
        return budget[dim];
    }

    /// Pack \a dim lattice coordinates of \a bits bits into one key
    /// (bit-interleave). The most-significant bit of X[0] is the
    /// most-significant bit of the key; within each bit level, axis 0
    /// comes first.
    static uint64_t interleave(const uint64_t * X, short_t dim, unsigned bits)
    {
        GISMO_ASSERT(dim >= 0 && dim <= 3, "dim must be in {0,1,2,3}");
        GISMO_ASSERT((unsigned)dim*bits <= 64, "dim*bits exceeds 64 bits");
        if (bits < 64)
            for (short_t k = 0; k != dim; ++k)
                GISMO_ASSERT(X[k] < (1ULL<<bits), "coordinate "<<k<<" exceeds "<<bits<<" bits");

        uint64_t key = 0;
        for (int j = (int)bits - 1; j >= 0; --j)
            for (short_t i = 0; i != dim; ++i)
                key = (key << 1) | ((X[i] >> j) & 1ULL);
        return key;
    }

    /// Inverse of interleave(); writes \a dim coordinates into \a X.
    static void deinterleave(uint64_t key, short_t dim, unsigned bits, uint64_t * X)
    {
        GISMO_ASSERT(dim >= 0 && dim <= 3, "dim must be in {0,1,2,3}");
        GISMO_ASSERT((unsigned)dim*bits <= 64, "dim*bits exceeds 64 bits");

        for (short_t i = 0; i != dim; ++i) X[i] = 0;
        for (int j = (int)bits - 1; j >= 0; --j)
            for (short_t i = 0; i != dim; ++i)
                X[i] |= ((key >> ((unsigned)dim*(unsigned)j + (unsigned)(dim-1-i))) & 1ULL) << j;
    }

    /// Morton (Z-order) encode: identical to interleave().
    static uint64_t mortonEncode(const uint64_t * coords, short_t dim, unsigned bits)
    {
        return interleave(coords, dim, bits);
    }

    /// Morton (Z-order) decode: identical to deinterleave().
    static void mortonDecode(uint64_t key, short_t dim, unsigned bits, uint64_t * coords)
    {
        deinterleave(key, dim, bits, coords);
    }

    /// Hilbert encode: Skilling axes-to-transpose on a local copy of
    /// \a coords, then interleave(). dim==0 and dim==1 are handled as
    /// special cases (Skilling's loops are not the identity at n==1).
    static uint64_t hilbertEncode(const uint64_t * coords, short_t dim, unsigned bits)
    {
        GISMO_ASSERT(dim >= 0 && dim <= 3, "dim must be in {0,1,2,3}");
        GISMO_ASSERT((unsigned)dim*bits <= 64, "dim*bits exceeds 64 bits");
        if (bits < 64)
            for (short_t k = 0; k != dim; ++k)
                GISMO_ASSERT(coords[k] < (1ULL<<bits), "coordinate "<<k<<" exceeds "<<bits<<" bits");

        if (dim == 0) return 0;
        if (dim == 1) return coords[0];

        uint64_t X[3];
        for (short_t i = 0; i != dim; ++i) X[i] = coords[i];
        internal::gsHilbertAxesToTranspose(X, (int)bits, (int)dim);
        return interleave(X, dim, bits);
    }

    /// Hilbert decode: deinterleave(), then Skilling transpose-to-axes.
    /// dim==0 and dim==1 are handled as special cases.
    static void hilbertDecode(uint64_t key, short_t dim, unsigned bits, uint64_t * coords)
    {
        GISMO_ASSERT(dim >= 0 && dim <= 3, "dim must be in {0,1,2,3}");
        GISMO_ASSERT((unsigned)dim*bits <= 64, "dim*bits exceeds 64 bits");

        if (dim == 0) return;
        if (dim == 1) { coords[0] = key; return; }

        deinterleave(key, dim, bits, coords);
        internal::gsHilbertTransposeToAxes(coords, (int)bits, (int)dim);
    }

    /// Dispatch on \a c.
    static uint64_t encode(Curve c, const uint64_t * coords, short_t dim, unsigned bits)
    {
        return (c == Hilbert) ? hilbertEncode(coords, dim, bits) : mortonEncode(coords, dim, bits);
    }

    /// Dispatch on \a c.
    static void decode(Curve c, uint64_t key, short_t dim, unsigned bits, uint64_t * coords)
    {
        if (c == Hilbert) hilbertDecode(key, dim, bits, coords);
        else              mortonDecode (key, dim, bits, coords);
    }

    // ---------------- geometric part -----------------------------------------

    /// Empty curve: geoDim()==0, curveDim()==0, encode() always returns 0.
    gsSpaceFillingCurve()
    : m_geoDim(0), m_curveDim(0), m_bits(0), m_curve(Hilbert)
    { }

    /** \brief Curve over the box \a box
     *
     * \param box    geoDim x 2 matrix; column 0 is the lower corner, column 1 the
     *               upper corner (the layout of gsMultiPatch::boundingBox)
     * \param curve  Morton or Hilbert
     * \param relTol relative tolerance of the degenerate-axis test
     */
    explicit gsSpaceFillingCurve(const gsMatrix<real_t> & box,
                                 Curve curve = Hilbert,
                                 real_t relTol = 1e-12)
    : m_curve(curve)
    {
        GISMO_ENSURE(box.cols() == 2, "box must have exactly two columns (lower, upper corner), got "<<box.cols());
        GISMO_ENSURE(box.rows() >= 1 && box.rows() <= 3, "box must have 1, 2 or 3 rows (geoDim), got "<<box.rows());
        GISMO_ENSURE(relTol >= 0, "relTol must be non-negative, got "<<relTol);

        m_geoDim = (short_t)box.rows();

        m_lo.resize(m_geoDim);
        m_extent.resize(m_geoDim);
        real_t maxExtent = 0;
        for (short_t i = 0; i != m_geoDim; ++i)
        {
            const real_t lo = box(i,0);
            const real_t hi = box(i,1);
            GISMO_ENSURE(hi >= lo, "box has hi["<<i<<"]="<<hi<<" < lo["<<i<<"]="<<lo);
            m_lo[i]     = lo;
            m_extent[i] = hi - lo;
            if (m_extent[i] > maxExtent) maxExtent = m_extent[i];
        }

        m_activeAxes.clear();
        for (short_t i = 0; i != m_geoDim; ++i)
            if (m_extent[i] > relTol * maxExtent)
                m_activeAxes.push_back(i);

        m_curveDim = (short_t)m_activeAxes.size();
        m_bits = bitsPerAxis(m_curveDim);
    }

    short_t  geoDim()   const { return m_geoDim; }   ///< rows of the box (0 for the empty curve)
    short_t  curveDim() const { return m_curveDim; }   ///< number of non-degenerate axes, 0..geoDim()
    unsigned bits()     const { return m_bits; }   ///< == bitsPerAxis(curveDim())
    Curve    curve()    const { return m_curve; }

    /// Indices of the non-degenerate axes, ascending; size() == curveDim()
    const std::vector<short_t> & activeAxes() const { return m_activeAxes; }

    /// Lattice coordinate of \a x along geometric axis \a axis; 0 if that axis is degenerate
    uint64_t quantize(real_t x, short_t axis) const
    {
        // Determine whether axis is active (linear scan; curveDim() <= 3)
        bool active = false;
        for (std::size_t k = 0; k != m_activeAxes.size(); ++k)
            if (m_activeAxes[k] == axis) { active = true; break; }
        if (!active) return 0;
        if (m_bits == 0) return 0;

        real_t t = (x - m_lo[axis]) / m_extent[axis]; // m_extent[axis] > 0 by construction
        if (t < 0) t = 0;
        if (t > 1) t = 1;

        const uint64_t L = 1ULL << m_bits;
        uint64_t q = (uint64_t)math::floor(t * (real_t)L);
        if (q > L-1) q = L-1;
        return q;
    }

    /// Key of column \a col of \a pts (pts.rows() >= geoDim())
    uint64_t encode(const gsMatrix<real_t> & pts, index_t col) const
    {
        GISMO_ASSERT(pts.rows() >= m_geoDim && col < pts.cols(), "point out of range for this curve");
        if (m_curveDim == 0) return 0;
        uint64_t c[3];
        for (short_t k = 0; k != m_curveDim; ++k)
            c[k] = quantize(pts(m_activeAxes[k], col), m_activeAxes[k]);
        return encode(m_curve, c, m_curveDim, m_bits);
    }

    /// Key of the point \a pt (pt.size() >= geoDim())
    uint64_t encode(const gsVector<real_t> & pt) const
    {
        GISMO_ASSERT(pt.size() >= m_geoDim, "point out of range for this curve");
        if (m_curveDim == 0) return 0;
        uint64_t c[3];
        for (short_t k = 0; k != m_curveDim; ++k)
            c[k] = quantize(pt[m_activeAxes[k]], m_activeAxes[k]);
        return encode(m_curve, c, m_curveDim, m_bits);
    }

    /// Centre of the lattice cell of \a key; \a pt is resized to geoDim().
    ///
    /// Guarantee: for any point p inside the box and any active axis i,
    /// |decode(encode(p))[i] - p[i]| <= 0.5 * extent[i] / 2^bits.
    void decode(uint64_t key, gsVector<real_t> & pt) const
    {
        pt.setZero(m_geoDim);
        for (short_t i = 0; i != m_geoDim; ++i)
            pt[i] = m_lo[i];
        if (m_curveDim == 0) return;

        uint64_t c[3];
        decode(m_curve, key, m_curveDim, m_bits, c);
        for (short_t k = 0; k != m_curveDim; ++k)
        {
            const short_t i = m_activeAxes[k];
            pt[i] = m_lo[i] + ( (real_t)c[k] + 0.5 ) * m_extent[i] / (real_t)(1ULL << m_bits);
        }
    }

private:

    short_t  m_geoDim;
    short_t  m_curveDim;
    unsigned m_bits;
    Curve    m_curve;

    std::vector<real_t>   m_lo;
    std::vector<real_t>   m_extent;
    std::vector<short_t>  m_activeAxes;
};

} // namespace gismo
