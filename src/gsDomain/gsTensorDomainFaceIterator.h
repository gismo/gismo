/** @file gsTensorDomainFaceIterator.h

    @brief Iterator over the interior faces of a per-direction break grid

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): H.M. Verhelst, A. Mantzaflaris
*/

#pragma once

#include <gsDomain/gsDomainFaceIterator.h>
#include <vector>

namespace gismo
{

/// Face predicate selecting every interior face of the grid.
/// Grids without a trimming concept report every element as interior.
struct AllFaces
{
    bool    operator()(size_t /*l*/, size_t /*r*/) const { return true;  }
    short_t sign      (size_t /*flat*/)            const { return -1;    }
};

/**
 * @brief Iterates the interior faces of a tensor-structured grid given by a
 * per-direction break vector, filtered by a \a FaceOp predicate.
 *
 * A \a FaceOp is a small copyable functor providing both
 * \code
 * bool    operator()(size_t leftFlatIndex, size_t rightFlatIndex) const; // keep this face?
 * short_t sign(size_t flatIndex) const;                                  // -1 / 0 / +1
 * \endcode
 * with flat indices being level-0 element indices, lexicographic with
 * direction 0 fastest. Signs are obtained through the predicate, giving
 * exactly one source of truth, checked at compile time.
 *
 * This class knows nothing about gsTensorDomain: it only consumes a break
 * grid, which lets gsTensorDomain include this header without a cycle and
 * lets other domains (e.g. gsTrimmedDomain) reuse it verbatim with their own
 * level-0 break vector and a sign-aware FaceOp.
 *
 * \ingroup Domain
 */
template<class T, int D, class FaceOp /*= AllFaces*/>
class gsTensorDomainFaceIterator : public gsDomainFaceIterator<T>
{
    typedef typename gsDomainIterator<T>::uPtr domainIter;
public:
    /// \param breaks breaks[j] = element boundaries in direction j, size N_j+1
    /// \param op     face predicate, see the FaceOp concept above
    explicit gsTensorDomainFaceIterator(const std::vector< std::vector<T> > & breaks,
                                        FaceOp op = FaceOp())
    :
    m_breaks(breaks),
    m_nelem(D),
    m_stride(D),
    m_op(give(op)),
    m_dir(0),
    m_good(false)
    {
        GISMO_ASSERT((index_t)breaks.size()==D,
                     "gsTensorDomainFaceIterator: expected "<<D<<" break vectors, got "<<breaks.size());
        size_t stride = 1;
        for (short_t j = 0; j < D; ++j)
        {
            GISMO_ASSERT(m_breaks[j].size()>=2,
                         "gsTensorDomainFaceIterator: direction "<<j<<" needs at least 2 breaks.");
            m_nelem[j]  = m_breaks[j].size()-1;
            m_stride[j] = stride;
            stride *= m_nelem[j];
        }
        reset();
    }

    gsTensorDomainFaceIterator(const gsTensorDomainFaceIterator & other) = default;
    domainIter clone() const override
    { return domainIter(new gsTensorDomainFaceIterator(*this)); }

    /// Number of faces the same (breaks, op) pair yields
    static size_t numFaces(const std::vector< std::vector<T> > & breaks,
                           FaceOp op = FaceOp())
    {
        gsTensorDomainFaceIterator it(breaks, op);
        size_t n = 0;
        for (; it.good(); it.next()) ++n;
        return n;
    }

    /// True while the iterator stands on a valid face
    bool good() const { return m_good; }

    void next() override
    {
        while (m_good)
        {
            m_good = advanceRaw();
            if (m_good && m_op(leftFlat(), rightFlat())) break;
        }
        if (m_good) updateSide();
    }

    void next(index_t increment) override
    {
        for (index_t i = 0; i!=increment && m_good; ++i)
            next();
    }

    void reset() override
    {
        m_dir = 0;
        while (m_dir < D && m_nelem[m_dir] < 2) ++m_dir;
        m_cur.setZero();
        m_good = (m_dir < D);
        if (m_good && !m_op(leftFlat(), rightFlat()))
            next();
        else if (m_good)
            updateSide();
    }

    gsVector<T> lowerCorner() const override
    {
        GISMO_ASSERT(m_good, "gsTensorDomainFaceIterator: lowerCorner() called past the end.");
        gsVector<T> result(D);
        for (short_t j = 0; j < D; ++j)
            result[j] = (j==m_dir) ? m_breaks[j][m_cur[j]+1] : m_breaks[j][m_cur[j]];
        return result;
    }

    gsVector<T> upperCorner() const override
    {
        GISMO_ASSERT(m_good, "gsTensorDomainFaceIterator: upperCorner() called past the end.");
        gsVector<T> result(D);
        // Identical to lowerCorner() in every component: the face has zero
        // extent, degenerating to the shared boundary at m_cur[m_dir]+1.
        for (short_t j = 0; j < D; ++j)
            result[j] = m_breaks[j][m_cur[j]+1];
        return result;
    }

    const T getPerpendicularCellSize() const override
    {
        GISMO_ASSERT(m_good, "gsTensorDomainFaceIterator: getPerpendicularCellSize() called past the end.");
        return m_breaks[m_dir][m_cur[m_dir]+1] - m_breaks[m_dir][m_cur[m_dir]];
    }

    const T getPerpendicularCellSizeRight() const override
    {
        GISMO_ASSERT(m_good, "gsTensorDomainFaceIterator: getPerpendicularCellSizeRight() called past the end.");
        return m_breaks[m_dir][m_cur[m_dir]+2] - m_breaks[m_dir][m_cur[m_dir]+1];
    }

    size_t leftElementId() const override
    {
        GISMO_ASSERT(m_good, "gsTensorDomainFaceIterator: leftElementId() called past the end.");
        return leftFlat();
    }

    size_t rightElementId() const override
    {
        GISMO_ASSERT(m_good, "gsTensorDomainFaceIterator: rightElementId() called past the end.");
        return rightFlat();
    }

    short_t leftSign() const override
    {
        GISMO_ASSERT(m_good, "gsTensorDomainFaceIterator: leftSign() called past the end.");
        return m_op.sign(leftFlat());
    }

    short_t rightSign() const override
    {
        GISMO_ASSERT(m_good, "gsTensorDomainFaceIterator: rightSign() called past the end.");
        return m_op.sign(rightFlat());
    }

private:

    size_t leftFlat()  const
    { size_t f=0; for (short_t j=0;j<D;++j) f += (size_t)m_cur[j]*m_stride[j]; return f; }

    size_t rightFlat() const { return leftFlat() + m_stride[m_dir]; }

    /// Raw cursor advance -- the odometer over the current direction's face
    /// lattice. Component j ranges over [0, N_j-1], except the normal
    /// direction, which ranges over [0, N_dir-2] (the left element of a face
    /// is never the last one in its direction).
    bool advanceRaw()
    {
        for (short_t j = 0; j < D; ++j)
        {
            const index_t lim = (j==m_dir) ? (index_t)m_nelem[j]-1 : (index_t)m_nelem[j];
            if (++m_cur[j] < lim) return true;
            m_cur[j] = 0;
        }
        // multi-index wrapped: move on to the next direction that has faces
        do { ++m_dir; } while (m_dir < D && m_nelem[m_dir] < 2);
        if (m_dir >= D) return false;
        m_cur.setZero();
        return true;
    }

    /// Refreshes the stored side, since the normal direction changes as
    /// iteration crosses from one direction's face lattice to the next.
    void updateSide()
    {
        const index_t p = this->m_pside.patch;
        this->m_pside = patchSide(p, boxSide(m_dir, true));
    }

    std::vector< std::vector<T> > m_breaks;  // size D, m_breaks[j].size() == N_j+1
    std::vector<size_t> m_nelem;             // N_j = m_breaks[j].size()-1
    std::vector<size_t> m_stride;            // stride_0 = 1, stride_j = stride_{j-1}*N_{j-1}
    FaceOp   m_op;
    short_t  m_dir;                          // current normal direction
    gsVector<index_t,D> m_cur;               // multi-index of the LEFT element
    bool     m_good;
};

} // namespace gismo
