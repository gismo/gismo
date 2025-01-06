/** @file gsSparseRows.h

    @brief A specialized sparse matrix class which stores separately

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, A. Mantzaflaris
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
 */
template <class T, bool IsRowMajor = true>
class gsSparseRows // todo: rename as gsFiberMatrix
{
public:
    typedef gsEigen::SparseVector<T> Fiber;

    struct RowBlockXpr;

    gsSparseRows()
    { }

    gsSparseRows(index_t rows, index_t cols)
    : m_fibers(rows)
    {
        for (index_t i = 0; i < rows; ++i)
            m_fibers[i] = new Fiber(cols);
    }

    gsSparseRows(const gsSparseRows& other)
    : m_fibers(other.rows())
    {
        for (int i = 0; i < rows(); ++i)
            m_fibers[i] = new Fiber( *other.m_fibers[i] );
    }

    gsSparseRows(const RowBlockXpr& rowxpr)
    : m_fibers(rowxpr.num)
    {
        for (index_t i = 0; i < rowxpr.num; ++i)
            m_fibers[i] = new Fiber( *rowxpr.mat.m_fibers[rowxpr.start + i] );
    }

    ~gsSparseRows()
    {
        clear();
    }

    gsSparseRows& operator= (const gsSparseRows other)
    {
        this->swap( other );
        return *this;
    }

    gsSparseRows& operator= (const RowBlockXpr& rowxpr)
    {
        gsSparseRows temp(rowxpr);
        this->swap( temp );
        return *this;
    }

    index_t fibers() const { return m_fibers.size(); }

    index_t rows() const { return m_fibers.size(); }
    index_t cols() const { return (m_fibers.size() > 0) ? m_fibers[0]->size() : 0; }

    Fiber& fiber(index_t i)              { return *m_fibers[i]; }
    const Fiber& fiber(index_t i) const { return *m_fibers[i]; }

    Fiber& row(index_t i)             { return *m_fibers[i]; }
    const Fiber& row(index_t i) const { return *m_fibers[i]; }

    T & coef(index_t i, index_t j)
    {
        if (!IsRowMajor) std::swap(i,j);
        return m_fibers[i]->coeff(j);
    }

    T & coeffRef(index_t i, index_t j)
    {
        if (!IsRowMajor) std::swap(i,j);
        return m_fibers[i]->coeffRef(j);
    }

    void clear()
    {
        for (int i = 0; i < fibers(); ++i)
            delete m_fibers[i];
        m_fibers.clear();
    }

    void swap(gsSparseRows& other)
    {
        m_fibers.swap( other.m_fibers );
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
        for (index_t i = 0; i < fibers(); ++i)
            std::fill(m_fibers[i]->valuePtr(), m_fibers[i]->valuePtr() + nonZeros(), (T)0.);
    }

    void resize(index_t rows, index_t cols)
    {
        GISMO_ASSERT( rows >= 0 && cols >= 0, "Invalid row/col in resize.");
        if (!IsRowMajor) std::swap(rows,cols);

        clear();
        m_fibers.resize(rows);
        for (index_t i = 0; i < rows; ++i)
            m_fibers[i] = new Fiber(cols);
    }

    void reservePerColumn(index_t nz)
    {
        for (index_t i = 0; i < fibers(); ++i)
            m_fibers[i]->reserve(nz);
    }

    void conservativeResize(index_t newRows, index_t newCols)
    {
        GISMO_ASSERT(rows() > 0 && cols() != newCols, "Cannot resize columns -- not implemented");
        if (!IsRowMajor) std::swap(newRows,newCols);

        const index_t oldRows = fibers();
        resizeFibers(newRows);

        // allocate newly added rows, if any
        for (index_t i = oldRows; i < newRows; ++i)
            m_fibers[i] = new Fiber(newCols);
    }

    void duplicateRow(index_t k) //..
    {
        GISMO_ASSERT( 0 <= k && k < fibers(), "k out of bounds.");

        //todo, something like: m_fibers.insert(m_fibers.begin()+k, new Fiber( fiber(k) );
        
        // add one new fiber
        resizeFibers( fibers() + 1 );

        // shift rows [k+1,...) down to [k+2,...)
        for (index_t i = fibers() - 1; i > k + 1; --i)
            m_fibers[i] = m_fibers[i-1];

        // allocate new row
        m_fibers[k+1] = new Fiber( row(k) );
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
            nnz += m_fibers[i]->nonZeros();
        return nnz;
    }

    template <class Derived>
    void toSparseMatrix(gsEigen::SparseMatrixBase<Derived>& m) const
    {
        m.derived().resize( rows(), cols() );
        m.derived().reserve( nonZeros() );
        for (index_t i = 0; i < fibers(); ++i)
        {
            for (typename Fiber::InnerIterator it(*m_fibers[i]); it; ++it)
                m.derived().insert(i, it.index()) = it.value();
        }
        m.derived().makeCompressed();
    }

    template <class Derived>
    void toSparseMatrixTransposed(gsEigen::SparseMatrixBase<Derived>& m) const
    {
        m.derived().resize(cols(), rows() );
        m.derived().reserve( nonZeros() );
        for (index_t i = 0; i < fibers(); ++i)
        {
            for (typename Fiber::InnerIterator it(*m_fibers[i]); it; ++it)
                m.derived().coeffRef(it.index(), i) = it.value();
        }
        m.derived().makeCompressed();
    }

    struct RowBlockXpr
    {
        RowBlockXpr(const gsSparseRows& _mat, index_t _start, index_t _num)
        : mat(const_cast<gsSparseRows&>(_mat)), start(_start), num(_num)
        {
            // HACK: We cast away the constness of the matrix, otherwise we would need two versions of
            // this expression class.
            // It's still safe because the row block methods in gsSparseRows above return the proper constness.
            GISMO_ASSERT( 0 <= num && 0 <= start   , "Invalid block.");
            GISMO_ASSERT( start < mat.rows()       , "Invalid block.");
            GISMO_ASSERT( start + num <= mat.rows(), "Invalid block.");
        }

        gsSparseRows & mat;
        index_t start, num;

        RowBlockXpr& operator= (const RowBlockXpr& other)
        {
            GISMO_ASSERT(num == other.num, "Wrong size in assignment.");
            for (index_t i = 0; i < num; ++i)
                mat.row(start + i) = other.mat.row(other.start + i);
            return *this;
        }

        RowBlockXpr& operator= (const gsSparseRows& other)
        {
            GISMO_ASSERT(num == other.rows(), "Wrong size in assignment.");
            for (index_t i = 0; i < num; ++i)
                mat.row(start + i) = other.row(i);
            return *this;
        }
    };

private:
    std::vector< Fiber* > m_fibers;
  
    /// Change the number of fibers without allocating newly added rows
    void resizeFibers(index_t newRows)
    {
        // delete fibers which will be removed from the array
        // (does nothing if newRows >= fibers())
        for (index_t i = newRows; i < fibers(); ++i)
            delete m_fibers[i];

        m_fibers.resize(newRows);
    }

};

} // namespace gismo
