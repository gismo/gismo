/** @file gsFiberMatrix.h

    @brief A specialized sparse matrix class which stores separately
    each fiber.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

namespace gismo
{

/**
 * \brief A specialized sparse matrix class which stores each row
 *  as a separate sparse vector.
 *
 *  This allows efficient row resizing and insertion
 *  operations, particularly for knot insertion algorithms.
 *
 *  \subsection lazy Lazy column allocation
 *
 *  When constructed via \a resizeLazy, fibers are allocated on first
 *  insertion (allocate-on-first-touch) instead of eagerly up-front.
 *  Untouched fibers cost one null pointer (8 B). This is
 *  memory-lean for subdomain assembly where only a fraction of the
 *  global DOF columns receive entries.
 *
 *  Thread-safety: lazy allocation during pattern construction is
 *  race-free because \c _pattern::push already serializes per column
 *  via \c omp_set_lock before any insertion. During value assembly
 *  (\c _eval::push), \c coeffRef only touches pattern-existing entries
 *  (fiber already allocated). \c assembleJacobian 's pattern-free
 *  insertions run under \c omp critical — also safe.
 *
 *  Eager mode (the default, via \a resize) is bit-identical to the
 *  pre-lazy behaviour.
 */
template <class T, int Major = ColMajor> // RowMajor==0, ColMajor==1
class gsFiberMatrix
{
    static constexpr bool IsRowMajor = (Major==RowMajor);
public:
    typedef gsSparseVector<T> Fiber;
    typedef typename Fiber::iterator iterator;

    struct RowBlockXpr;

    gsFiberMatrix()
    : m_innerSize(0), m_lazy(false), m_reserveHint(0)
    { }

    gsFiberMatrix(index_t rows, index_t cols)
    : m_innerSize(0), m_lazy(false), m_reserveHint(0)
    {
        resize(rows, cols);
    }

    gsFiberMatrix(const gsFiberMatrix& other)
    : m_fibers(other.m_fibers.size()),
      m_innerSize(other.m_innerSize),
      m_lazy(other.m_lazy),
      m_reserveHint(other.m_reserveHint)
    {
        for (size_t i = 0; i < m_fibers.size(); ++i)
            m_fibers[i] = other.m_fibers[i] ? new Fiber( *other.m_fibers[i] ) : nullptr;
    }

    gsFiberMatrix(const RowBlockXpr& rowxpr)
    : m_fibers(rowxpr.num),
      m_innerSize(rowxpr.mat.m_innerSize),
      m_lazy(rowxpr.mat.m_lazy),
      m_reserveHint(rowxpr.mat.m_reserveHint)
    {
        for (index_t i = 0; i < rowxpr.num; ++i)
            m_fibers[i] = rowxpr.mat.m_fibers[rowxpr.start + i]
                            ? new Fiber( *rowxpr.mat.m_fibers[rowxpr.start + i] )
                            : nullptr;
    }

    ~gsFiberMatrix()
    {
        clear();
    }

    iterator begin(index_t j) const
    {
        GISMO_ASSERT(m_fibers[j], "Fiber "<<j<<" is null (lazy not yet allocated)");
        return iterator(*m_fibers[j]);
    }

#if EIGEN_HAS_RVALUE_REFERENCES
    gsFiberMatrix(gsFiberMatrix&& other)
    : m_fibers(give(other.m_fibers)),
      m_innerSize(other.m_innerSize),
      m_lazy(other.m_lazy),
      m_reserveHint(other.m_reserveHint)
    { other.m_innerSize = 0; other.m_lazy = false; other.m_reserveHint = 0; }

    /// Assignment operator
    gsFiberMatrix& operator= ( const gsFiberMatrix& other )
    {
        clear();
        m_fibers.resize(other.m_fibers.size());
        m_innerSize   = other.m_innerSize;
        m_lazy        = other.m_lazy;
        m_reserveHint = other.m_reserveHint;
        for (size_t i = 0; i < m_fibers.size(); ++i)
            m_fibers[i] = other.m_fibers[i] ? new Fiber( *other.m_fibers[i] ) : nullptr;
        return *this;
    }

    /// Move assignment operator
    gsFiberMatrix& operator= ( gsFiberMatrix&& other )
    {
        clear();
        m_fibers      = give(other.m_fibers);
        m_innerSize   = other.m_innerSize;
        m_lazy        = other.m_lazy;
        m_reserveHint = other.m_reserveHint;
        other.m_innerSize = 0;
        other.m_lazy      = false;
        other.m_reserveHint = 0;
        return *this;
    }
#else
    gsFiberMatrix& operator= (gsFiberMatrix other)
    {
        this->swap( other );
        return *this;
    }
#endif

    gsFiberMatrix& operator= (const RowBlockXpr& rowxpr)
    {
        gsFiberMatrix temp(rowxpr);
        this->swap( temp );
        return *this;
    }

    inline index_t fibers() const { return static_cast<index_t>(m_fibers.size()); }

    /** \returns the size of the storage major dimension,
     * i.e., the number of columns for a columns major matrix, and the number of rows otherwise */
    inline index_t outerSize() const
    { return static_cast<index_t>(m_fibers.size()); }

    /** \returns the size of the inner dimension according to the storage order,
      * i.e., the number of rows for a columns major matrix, and the number of cols otherwise */
    inline index_t innerSize() const
    { return m_innerSize; }

    void setZero()
    {
        for( auto & f : m_fibers)
            if (f) f->setZero();
    }

    /** \returns the number of rows of the matrix */
    inline index_t rows() const { return IsRowMajor ? outerSize() : innerSize(); }

    /** \returns the number of columns of the matrix */
    inline index_t cols() const { return IsRowMajor ? innerSize() : outerSize(); }

    Fiber& fiber(index_t i)
    {
        GISMO_ASSERT( i>=0 && i<outerSize(), "Invalid fiber: "<<i);
        GISMO_ENSURE(m_fibers[i], "Fiber "<<i<<" is null (lazy not yet allocated)");
        return *m_fibers[i];
    }
    const Fiber& fiber(index_t i) const
    {
        GISMO_ASSERT( i>=0 && i<outerSize(), "Invalid fiber: "<<i);
        GISMO_ENSURE(m_fibers[i], "Fiber "<<i<<" is null (lazy not yet allocated)");
        return *m_fibers[i];
    }

    Fiber& row(index_t i)
    {
        GISMO_ASSERT( i>=0 && i<rows(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows());
        GISMO_ENSURE(IsRowMajor, "Cannot access row in col-major fiber matrix");
        return fiberRef(i);
    }

    const Fiber& row(index_t i) const
    {
        GISMO_ASSERT( i>=0 && i<rows(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows());
        GISMO_ENSURE(IsRowMajor, "Cannot access row in col-major fiber matrix");
        GISMO_ENSURE(m_fibers[i], "Row "<<i<<" is null (lazy not yet allocated)");
        return *m_fibers[i];
    }

    Fiber& col(index_t i)
    {
        GISMO_ASSERT(i>=0 && i<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<cols()"<<"="<<cols());
        GISMO_ENSURE(!IsRowMajor, "Cannot access col in col-major fiber matrix");
        return fiberRef(i);
    }

    const Fiber& col(index_t i) const
    {
        GISMO_ASSERT(i>=0 && i<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<cols()"<<"="<<cols());
        GISMO_ENSURE(!IsRowMajor, "Cannot access col in row-major fiber matrix");
        GISMO_ENSURE(m_fibers[i], "Col "<<i<<" is null (lazy not yet allocated)");
        return *m_fibers[i];
    }

    T coeff(index_t i, index_t j) const
    {
        GISMO_ASSERT( i>=0 && i<rows() && j>=0 && j<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows()<<"  &&  "<<j<<">=0 && "<<j<<"<cols()"<<"="<<cols() );
        if (!IsRowMajor) std::swap(i,j);
        if (!m_fibers[i]) return T(0);
        return m_fibers[i]->coeff(j);
    }

    T & coeffRef(index_t i, index_t j)
    {
        GISMO_ASSERT( i>=0 && i<rows() && j>=0 && j<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows()<<"  &&  "<<j<<">=0 && "<<j<<"<cols()"<<"="<<cols() );
        if (!IsRowMajor) std::swap(i,j);
        return fiberRef(i).coeffRef(j);
    }

    T & insert(index_t i, index_t j)
    {
        GISMO_ASSERT( i>=0 && i<rows() && j>=0 && j<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows()<<"  &&  "<<j<<">=0 && "<<j<<"<cols()"<<"="<<cols() );
        if (!IsRowMajor) std::swap(i,j);
        return fiberRef(i).insert(j);
    }

    void insertExplicitZero(index_t i, index_t j)
    {
        GISMO_ASSERT( i>=0 && i<rows() && j>=0 && j<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows()<<"  &&  "<<j<<">=0 && "<<j<<"<cols()"<<"="<<cols() );
        if (!IsRowMajor) std::swap(i,j);
        fiberRef(i).data().atWithInsertion(j);
    }

    bool isExplicitZero(index_t i, index_t j) const
    {
        GISMO_ASSERT( i>=0 && i<rows() && j>=0 && j<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows()<<"  &&  "<<j<<">=0 && "<<j<<"<cols()"<<"="<<cols() );
        if (!IsRowMajor) std::swap(i,j);
        if (!m_fibers[i]) return true; // null fiber → no explicit zero slot present
        auto & vdata = m_fibers[i]->data();
        const index_t jj = vdata.searchLowerIndex(j);
        return ((jj==vdata.size()) || (vdata.index(jj)!=j));
    }

    void clear()
    {
        for (int i = 0; i < fibers(); ++i)
            delete m_fibers[i]; // safe on nullptr
        m_fibers.clear();
        m_innerSize   = 0;
        m_lazy        = false;
        m_reserveHint = 0;
    }

    //void prune()

    void swap(gsFiberMatrix& other)
    {
        m_fibers.swap( other.m_fibers );
        std::swap(m_innerSize,   other.m_innerSize);
        std::swap(m_lazy,        other.m_lazy);
        std::swap(m_reserveHint, other.m_reserveHint);
    }

    void setIdentity(index_t n)
    {
        GISMO_ASSERT( n >= 0, "n must be positive." );

        resize(n, n);

        for (index_t i = 0; i < n; ++i)
            m_fibers[i]->insert(i) = (T)(1.0);
    }

    void assignZero()
    {
        for (auto & fb : m_fibers)
            if (fb) std::fill(fb->valuePtr(), fb->valuePtr() + fb->nonZeros(), (T)0.);
    }

    /// \brief Eager resize: allocates every fiber up-front (default; bit-identical to pre-lazy behaviour).
    void resize(index_t rows, index_t cols)
    {
        GISMO_ASSERT( rows >= 0 && cols >= 0, "Invalid row/col in resize.");
        if (!IsRowMajor) std::swap(rows,cols);

        clear();
        m_fibers.resize(rows);
        for (index_t i = 0; i < rows; ++i)
            m_fibers[i] = new Fiber(cols);
        m_innerSize = cols;
        m_lazy      = false;
    }

    /// \brief Lazy resize: fibers are allocated on first touch (null until then).
    ///
    /// Untouched fibers cost one null pointer (8 B). Memory-lean for
    /// subdomain assembly where only a fraction of global DOF columns
    /// receive entries. Thread-safe under the pattern-then-value
    /// discipline of gsExprAssembler (see class docstring).
    void resizeLazy(index_t rows, index_t cols)
    {
        GISMO_ASSERT( rows >= 0 && cols >= 0, "Invalid row/col in resizeLazy.");
        if (!IsRowMajor) std::swap(rows,cols);

        clear();
        m_fibers.resize(rows, nullptr);
        m_innerSize = cols;
        m_lazy      = true;
    }

    void reservePerColumn(index_t nz)
    {
        if (m_lazy)
        {
            m_reserveHint = nz;
            for (index_t i = 0; i < fibers(); ++i)
                if (m_fibers[i]) m_fibers[i]->reserve(nz);
        }
        else
        {
            for (index_t i = 0; i < fibers(); ++i)
                m_fibers[i]->reserve(nz);
        }
    }

    template<typename Cont>
    void reserve(const Cont &nz)
    {
        GISMO_ASSERT(m_fibers.size()==(size_t)nz.size(), "Wrong size in nonzero vector.");
        for (index_t i = 0; i < fibers(); ++i)
            if (m_fibers[i]) m_fibers[i]->reserve(nz[i]);
    }

    void conservativeResize(index_t newRows, index_t newCols)
    {
        if (!IsRowMajor) std::swap(newRows,newCols);

        const index_t oldRows = fibers();
        const index_t m = std::min(oldRows, newRows);
        for (index_t i = 0; i < m; ++i)
            if (m_fibers[i]) m_fibers[i]->conservativeResize(newCols);

        if (newCols > m_innerSize) m_innerSize = newCols;
        resizeFibers(newRows);
    }

    void resizeFibers(index_t newSz)
    {
        const index_t oldSz = fibers();

        // delete any fibers which will be removed, if any
        for (index_t i = newSz; i < oldSz; ++i)
            delete m_fibers[i]; // safe on nullptr

        m_fibers.resize(newSz);

        // allocate newly added fibers, if any
        for (index_t i = oldSz; i < newSz; ++i)
            m_fibers[i] = m_lazy ? nullptr : new Fiber(m_innerSize);
    }

    void duplicateRow(index_t k) //..
    {
        GISMO_ASSERT(IsRowMajor &&  0 <= k && k < fibers(), "k out of bounds.");

        //todo, something like: m_fibers.insert(m_fibers.begin()+k, new Fiber( fiber(k) );

        // add one new fiber
        m_fibers.resize(m_fibers.size()+1);

        // shift rows [k+1,...) down to [k+2,...)
        for (index_t i = fibers() - 1; i > k + 1; --i)
            m_fibers[i] = m_fibers[i-1];

        // allocate new row
        m_fibers[k+1] = m_fibers[k] ? new Fiber( *m_fibers[k] ) : nullptr;
    }

    // row expressions //..

    RowBlockXpr topRows(index_t num)       { return RowBlockXpr(*this, 0, num); }
    const RowBlockXpr topRows(index_t num) const { return RowBlockXpr(*this, 0, num); }

    RowBlockXpr bottomRows(index_t num)       { return RowBlockXpr(*this, fibers() - num, num); }
    const RowBlockXpr bottomRows(index_t num) const { return RowBlockXpr(*this, fibers() - num, num); }

    RowBlockXpr middleRows(index_t start, index_t num)        { return RowBlockXpr(*this, start, num); }
    const RowBlockXpr middleRows(index_t start, index_t num) const  { return RowBlockXpr(*this, start, num); }

    index_t nonZeros() const
    {
        index_t nnz = 0;
        for (index_t i = 0; i < fibers(); ++i)
            if (m_fibers[i]) nnz += m_fibers[i]->nonZeros();
        return nnz;
    }

    gsVector<index_t> nonZerosPerFiber() const
    {
        gsVector<index_t> result(fibers());
        for (size_t i = 0; i != m_fibers.size(); ++i)
            result[i] = m_fibers[i] ? m_fibers[i]->nonZeros() : 0;
        return result;
    }

    gsSparseMatrix<T> toSparseMatrix() const
    {
        gsSparseMatrix<T> rvo;
        toSparseMatrix_into(rvo);
        return rvo;
    }

    template <class Derived>
    void toSparseMatrix_into(gsEigen::SparseMatrixBase<Derived>& m) const
    {
        m.derived().resize( rows(), cols() );
        m.derived().reserve( nonZerosPerFiber() );
        for (index_t i = 0; i < fibers(); ++i)
        {
            if (!m_fibers[i]) continue;
            for (typename Fiber::InnerIterator it(*m_fibers[i]); it; ++it)
                m.derived().insert(IsRowMajor?i:it.index(), IsRowMajor?it.index():i) = it.value();
        }
        m.derived().makeCompressed();
    }

    struct RowBlockXpr
    {
        RowBlockXpr(const gsFiberMatrix& _mat, index_t _start, index_t _num)
        : mat(const_cast<gsFiberMatrix&>(_mat)), start(_start), num(_num)
        {
            // HACK: We cast away the constness of the matrix, otherwise we would need two versions of
            // this expression class.
            // It's still safe because the row block methods in gsFiberMatrix above return the proper constness.
            GISMO_ASSERT( 0 <= num && 0 <= start   , "Invalid block.");
            GISMO_ASSERT( start < mat.rows()       , "Invalid block.");
            GISMO_ASSERT( start + num <= mat.rows(), "Invalid block.");
        }

        gsFiberMatrix & mat;
        index_t start, num;

        RowBlockXpr& operator= (const RowBlockXpr& other)
        {
            GISMO_ASSERT(num == other.num, "Wrong size in assignment.");
            for (index_t i = 0; i < num; ++i)
                mat.row(start + i) = other.mat.row(other.start + i);
            return *this;
        }

        RowBlockXpr& operator= (const gsFiberMatrix & other)
        {
            GISMO_ASSERT(num == other.rows(), "Wrong size in assignment.");
            for (index_t i = 0; i < num; ++i)
                mat.row(start + i) = other.row(i);
            return *this;
        }
    };

private:
    /// Allocate-on-first-touch: returns the fiber at \a outer, allocating
    /// it (with the reservation hint) if it is null. In eager mode this
    /// is a null-check + return — no allocation occurs.
    Fiber& fiberRef(index_t outer)
    {
        GISMO_ASSERT(outer >= 0 && outer < outerSize(), "Invalid fiber index");
        if (!m_fibers[outer])
        {
            m_fibers[outer] = new Fiber(m_innerSize);
            if (m_reserveHint > 0)
                m_fibers[outer]->reserve(m_reserveHint);
        }
        return *m_fibers[outer];
    }

    std::vector< Fiber* > m_fibers;
    index_t m_innerSize    = 0;  ///< inner dimension (fiber length); explicit so it survives null fibers
    bool    m_lazy         = false; ///< true when fibers may be null (allocated on first touch)
    index_t m_reserveHint  = 0;  ///< reservation hint applied at first-touch allocation (lazy mode)

}; // class gsFiberMatrix

} // namespace gismo