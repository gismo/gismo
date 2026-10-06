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
    { }

    gsFiberMatrix(index_t rows, index_t cols)
    : m_fibers(IsRowMajor?rows:cols)
    {
        for (size_t i = 0; i < m_fibers.size(); ++i)
            m_fibers[i] = new Fiber(cols);
    }

    gsFiberMatrix(const gsFiberMatrix& other)
    : m_fibers(other.outerSize())
    {
        for (size_t i = 0; i < m_fibers.size(); ++i)
            m_fibers[i] = new Fiber( *other.m_fibers[i] );
    }

    gsFiberMatrix(const RowBlockXpr& rowxpr)
    : m_fibers(rowxpr.num)
    {
        for (index_t i = 0; i < rowxpr.num; ++i)
            m_fibers[i] = new Fiber( *rowxpr.mat.m_fibers[rowxpr.start + i] );
    }

    ~gsFiberMatrix()
    {
        clear();
    }

    iterator begin(index_t j) const { return iterator(*m_fibers[j]); }

#if EIGEN_HAS_RVALUE_REFERENCES
    gsFiberMatrix(gsFiberMatrix&& other) : m_fibers(give(other.m_fibers)) {}

    /// Assignment operator
    gsFiberMatrix& operator= ( const gsFiberMatrix& other )
    {
        clear();
        m_fibers.resize(other.outerSize());
        for (size_t i = 0; i < m_fibers.size(); ++i)
            m_fibers[i] = new Fiber( *other.m_fibers[i] );
        return *this;
    }

    /// Move assignment operator
    gsFiberMatrix& operator= ( gsFiberMatrix&& other )
    {
        clear();
        m_fibers = give(other.m_fibers);
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

    inline index_t fibers() const { return m_fibers.size(); }

    /** \returns the size of the storage major dimension,
     * i.e., the number of columns for a columns major matrix, and the number of rows otherwise */
    inline index_t outerSize() const
    { return m_fibers.size(); }

    /** \returns the size of the inner dimension according to the storage order,
      * i.e., the number of rows for a columns major matrix, and the number of cols otherwise */
    inline index_t innerSize() const
    //    { return ( m_fibers.empty() ? 0 : m_fibers.front()->size() ); }
    { return ( m_fibers.size()>0 ? m_fibers.front()->size() : 0 ); }

    void setZero()
    {
        for( auto & f : m_fibers)
            f->setZero();
    }

    /** \returns the number of rows of the matrix */
    inline index_t rows() const { return IsRowMajor ? outerSize() : innerSize(); }

    /** \returns the number of columns of the matrix */
    inline index_t cols() const { return IsRowMajor ? innerSize() : outerSize(); }

    Fiber& fiber(index_t i)              { return *m_fibers[i]; }
    const Fiber& fiber(index_t i) const { return *m_fibers[i]; }

    Fiber& row(index_t i)
    {
        GISMO_ASSERT( i>=0 && i<rows(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows());
        GISMO_ENSURE(IsRowMajor, "Cannot access row in col-major fiber matrix");
        return *m_fibers[i];
    }

    const Fiber& row(index_t i) const
    {
        GISMO_ASSERT( i>=0 && i<rows(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows());
        GISMO_ENSURE(IsRowMajor, "Cannot access row in col-major fiber matrix");
        return *m_fibers[i];
    }

    Fiber& col(index_t i)
    {
        GISMO_ASSERT(i>=0 && i<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<cols()"<<"="<<cols());
        GISMO_ENSURE(!IsRowMajor, "Cannot access col in col-major fiber matrix");
        return *m_fibers[i];
    }

    const Fiber& col(index_t i) const
    {
        GISMO_ASSERT(i>=0 && i<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<cols()"<<"="<<cols());
        GISMO_ENSURE(!IsRowMajor, "Cannot access col in row-major fiber matrix");
        return *m_fibers[i];
    }

    T & coef(index_t i, index_t j)
    {
        GISMO_ASSERT( i>=0 && i<rows() && j>=0 && j<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows()<<"  &&  "<<j<<">=0 && "<<i<<"<cols()"<<"="<<cols() );
        if (!IsRowMajor) std::swap(i,j);
        return m_fibers[i]->coeff(j);
    }

    T & coeffRef(index_t i, index_t j)
    {
        GISMO_ASSERT( i>=0 && i<rows() && j>=0 && j<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows()<<"  &&  "<<j<<">=0 && "<<i<<"<cols()"<<"="<<cols() );
        if (!IsRowMajor) std::swap(i,j);
        return m_fibers[i]->coeffRef(j);
    }

    void insertExplicitZero(index_t i, index_t j)
    {
        GISMO_ASSERT( i>=0 && i<rows() && j>=0 && j<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows()<<"  &&  "<<j<<">=0 && "<<i<<"<cols()"<<"="<<cols() );
        if (!IsRowMajor) std::swap(i,j);
        m_fibers[i]->data().atWithInsertion(j);
    }

    bool isExplicitZero(index_t i, index_t j) const
    {
        GISMO_ASSERT( i>=0 && i<rows() && j>=0 && j<cols(), "Invalid element: "<<i<<">=0 && "<<i<<"<rows()"<<"="<<rows()<<"  &&  "<<j<<">=0 && "<<i<<"<cols()"<<"="<<cols() );
        if (!IsRowMajor) std::swap(i,j);
        auto & vdata = m_fibers[i]->data();
        const index_t jj = vdata.searchLowerIndex(j);
        return ((jj==vdata.size()) || (vdata.index(jj)!=j));
    }

    void clear()
    {
        for (int i = 0; i < fibers(); ++i)
            delete m_fibers[i];
        m_fibers.clear();
    }

    //void prune()
    
    void swap(gsFiberMatrix& other)
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
        for (auto & fb : m_fibers)
            std::fill(fb->valuePtr(), fb->valuePtr() + fb->nonZeros(), (T)0.);
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

    template<typename Cont>
    void reserve(const Cont &nz)
    {
        GISMO_ASSERT(m_fibers.size()==(size_t)nz.size(), "Wrong size in nonzero vector.");
        for (index_t i = 0; i < fibers(); ++i)
            m_fibers[i]->reserve(nz[i]);
    }

    void conservativeResize(index_t newRows, index_t newCols)
    {
        if (!IsRowMajor) std::swap(newRows,newCols);

        const index_t oldRows = fibers();

        // delete any fibers which will be removed, if any
        for (index_t i = newRows; i < oldRows; ++i)
            delete m_fibers[i];

        m_fibers.resize(newRows);

        // allocate newly added fibers, if any
        for (index_t i = oldRows; i < newRows; ++i)
            m_fibers[i] = new Fiber(newCols);

        const index_t m = std::min(oldRows, newRows);
        for (index_t i = 0; i < m; ++i)
            m_fibers[i]->conservativeResize(newCols);
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

    gsVector<index_t> nonZerosPerFiber() const
    {
        gsVector<index_t> result(fibers());
        for (size_t i = 0; i != m_fibers.size(); ++i)
            result[i] = m_fibers[i]->nonZeros();
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
