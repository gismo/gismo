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
template <class T>
class gsSparseRows // rename as gsFiberMatrix
{
public:
    typedef gsEigen::SparseVector<T> Row;

    struct RowBlockXpr;

    gsSparseRows()
    { }

    gsSparseRows(index_t rows, index_t cols)
    : m_rows(rows)
    {
        for (index_t i = 0; i < rows; ++i)
            m_rows[i] = new Row(cols);
    }

    gsSparseRows(const gsSparseRows& other)
    : m_rows(other.rows())
    {
        for (int i = 0; i < rows(); ++i)
            m_rows[i] = new Row( *other.m_rows[i] );
    }

    gsSparseRows(const RowBlockXpr& rowxpr)
    : m_rows(rowxpr.num)
    {
        for (index_t i = 0; i < rowxpr.num; ++i)
            m_rows[i] = new Row( *rowxpr.mat.m_rows[rowxpr.start + i] );
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

    index_t rows() const { return m_rows.size(); }
    index_t cols() const { return (m_rows.size() > 0) ? m_rows[0]->size() : 0; }

    Row& row(index_t i)             { return *m_rows[i]; }
    const Row& row(index_t i) const { return *m_rows[i]; }

    T & coef(index_t i, index_t j) { return m_rows[i]->coeff(j); }

    T & coeffRef(index_t i, index_t j) { return m_rows[i]->coeffRef(j); }

    void clear()
    {
        for (int i = 0; i < rows(); ++i)
            delete m_rows[i];
        m_rows.clear();
    }

    void swap(gsSparseRows& other)
    {
        m_rows.swap( other.m_rows );
    }

    void setIdentity(index_t n)
    {
        GISMO_ASSERT( n >= 0, "n must be positive." );

        resize(n, n);

        for (index_t i = 0; i < n; ++i)
            m_rows[i]->insert(i) = (T)(1.0);
    }

    void assignZero()
    {
        for (index_t i = 0; i < rows(); ++i)
            std::fill(m_rows[i]->valuePtr(), m_rows[i]->valuePtr() + nonZeros(), (T)0.);
    }

    void resize(index_t rows, index_t cols)
    {
        GISMO_ASSERT( rows >= 0 && cols >= 0, "Invalid row/col in resize.");

        clear();
        m_rows.resize(rows);
        for (index_t i = 0; i < rows; ++i)
            m_rows[i] = new Row(cols);
    }

    void reservePerColumn(index_t nz)
    {
        for (index_t i = 0; i < rows(); ++i)
            m_rows[i]->reserve(nz);
    }

    void conservativeResize(index_t newRows, index_t newCols)
    {
        GISMO_ASSERT(rows() > 0 && cols() != newCols, "Cannot resize columns -- not implemented");

        const index_t oldRows = rows();
        resizeRows(newRows);

        // allocate newly added rows, if any
        for (index_t i = oldRows; i < newRows; ++i)
            m_rows[i] = new Row(newCols);
    }

    void duplicateRow(index_t k)
    {
        GISMO_ASSERT( 0 <= k && k < rows(), "k out of bounds.");

        // add one new row
        resizeRows( rows() + 1 );

        // shift rows [k+1,...) down to [k+2,...)
        for (index_t i = rows() - 1; i > k + 1; --i)
            m_rows[i] = m_rows[i-1];

        // allocate new row
        m_rows[k+1] = new Row( row(k) );
    }

    // row expressions

    RowBlockXpr topRows(index_t num)       { return RowBlockXpr(*this, 0, num); }
    const RowBlockXpr topRows(index_t num) const { return RowBlockXpr(*this, 0, num); }

    RowBlockXpr bottomRows(index_t num)       { return RowBlockXpr(*this, rows() - num, num); }
    const RowBlockXpr bottomRows(index_t num) const { return RowBlockXpr(*this, rows() - num, num); }

    RowBlockXpr middleRows(index_t start, index_t num)        { return RowBlockXpr(*this, start, num); }
    const RowBlockXpr middleRows(index_t start, index_t num) const  { return RowBlockXpr(*this, start, num); }

    index_t nonZeros() const
    {
        index_t nnz = 0;
        for (index_t i = 0; i < rows(); ++i)
            nnz += m_rows[i]->nonZeros();
        return nnz;
    }

    template <class Derived>
    void toSparseMatrix(gsEigen::SparseMatrixBase<Derived>& m) const
    {
        m.derived().resize( rows(), cols() );
        m.derived().reserve( nonZeros() );
        for (index_t i = 0; i < rows(); ++i)
        {
            for (typename Row::InnerIterator it(*m_rows[i]); it; ++it)
                m.derived().insert(i, it.index()) = it.value();
        }
        m.derived().makeCompressed();
    }

    template <class Derived>
    void toSparseMatrixTransposed(gsEigen::SparseMatrixBase<Derived>& m) const
    {
        m.derived().resize(cols(), rows() );
        m.derived().reserve( nonZeros() );
        for (index_t i = 0; i < rows(); ++i)
        {
            for (typename Row::InnerIterator it(*m_rows[i]); it; ++it)
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
    std::vector< Row* > m_rows;

    /// Change the number of rows without allocating newly added rows
    void resizeRows(index_t newRows)
    {
        // delete rows which will be removed from the array
        // (does nothing if newRows >= rows())
        for (index_t i = newRows; i < rows(); ++i)
            delete m_rows[i];

        m_rows.resize(newRows);
    }

};

} // namespace gismo
