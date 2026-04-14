/** @file gsKroneckerMatrix.h

    @brief

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s):
*/

#pragma once

#include <gsSolver/gsIterativeSolver.h>
#include <gsMatrix/gsFiberMatrix.h>
#include <gsTensor/gsTensorTools.h>


namespace gismo
{

template<class T, int d>
class gsKroneckerMatrixIter
{
    typedef gsFiberMatrix<T,ColMajor> FiberMatrix;
    typedef std::vector<FiberMatrix> cwiseMat;
    typedef typename gsFiberMatrix<T>::iterator cIter;
public:

    gsKroneckerMatrixIter(const cwiseMat & _cwise)
    : m_cwise(_cwise)
    {
        const index_t n = m_cwise.size();
        GISMO_ASSERT(d==-1 || d==n, "Invalid input");

        // assuming square matrix blocks:
        m_rank = m_cwise.front().cols() / m_cwise.front().rows();
        //gsDebugVar(m_rank);

        gsVector<index_t,d> csz(n);
        m_rit.resize( m_rank * n);
        for(index_t i = 0; i!= n; ++i)
        {
            csz  [i] = m_cwise[i].cols() / m_rank;
            m_rit[i] = m_cwise[i].begin(0);
            for(size_t k = n; k!=m_rit.size(); k+=n )
                m_rit[k+i] = m_cwise[i].begin((k/n)*csz[i]);
        }

        m_col = gsGridIterator<index_t,CUBE,d>(csz);
    }

    operator bool() const {return m_col;}

    // const point & operator*() const {return m_cur;}
    // const point * operator->() const {return &m_cur;}

    inline gsKroneckerMatrixIter & operator++()
    {
        const index_t n = m_col->size();// as member ?
        const index_t rkn = n * m_rank;  // as member ?
        index_t t;
        for (t = 0; t != n; ++t)
        {
            // assuming common structure for all skeleton matrices
            for( index_t k = n; k!=rkn; k+=n) ++m_rit[k+t];
            if ( ++m_rit[t] )
                break;
            else
                if ( t + 1 == n )
                { ++t; break;}
                else
                for(size_t k = 0; k!=rkn; k+=n)
                    m_rit[k+t] = m_cwise[t].begin((k/n)*m_cwise[t].rows() + m_col->at(t));
        }

        // RETURN (t == n);
        if (t == n) // row end
        {
            ++m_col;
            for(index_t i = 0; i!=n; ++i)
                for(size_t k = 0; k!=m_rit.size(); k+=n)
                    m_rit[k+i] = m_cwise[i].begin( (k/n)*m_cwise[i].rows() + m_col->at(i) );
        }

        return *this;
    }

    index_t rowIndex(const index_t i) const {return m_rit[i].index();}
    index_t colIndex(const index_t i) const {return m_col->at(i);}
    //index_t cwiseValue(const index_t i) const {return m_rit[i].value();} //! rank 1

    T value(index_t & nr, index_t & nc) const
    {
        const index_t n = m_col->size();
        const index_t l = m_col->size()-1;

        T tmp = m_rit[l].value();
        nr    = m_rit[l].index();
        nc    = m_col->at(l);

        // row/col indices and first component value
        for(index_t i = l-1; i>=0; --i)
        {
            nr   = nr * m_cwise[i].rows()    + m_rit[i].index();
            nc   = nc * (m_col.upper()[i]+1) + m_col->at(i);
            tmp *= m_rit[i].value();
        }

        // add up upto rank
        for( size_t k = n; k!=m_rit.size(); k+=n )
        {
            T tmpr = m_rit[k].value();
            for(index_t i = 1; i!=n; ++i)
                tmpr *= m_rit[k+i].value();
            tmp += tmpr;
        }

        return tmp;
    }

    T value() const
    {
        const index_t & n = m_col->size();

        T tmp = m_rit[0].value();
        for(index_t i = 1; i<n; ++i)
            tmp *= m_rit[i].value();

        for( size_t k = n; k!=m_rit.size(); k+=n )
        {
            T tmpr = m_rit[k].value();
            for(index_t i = 1; i!=n; ++i)
                tmpr *= m_rit[k+i].value();
            tmp += tmpr;
        }
        return tmp;
    }

private:
    const cwiseMat & m_cwise;
    index_t m_rank;

    std::vector<cIter> m_rit;       ///< Row
    gsGridIterator<index_t,CUBE,d> m_col; ///< Column
};

template<class T, short_t d = -1> class gsKroneckerMatrix;

/**
   A Matrix represented by a Canonical representation (ie. a sum of
   Kronecker products)

   In theory, this matrix is a (i1,..id, j1,..,jd) tensor,
   ie. 2d-tensor of canonical rank R

   Taking a matricization with respect to i's and j's we obtain the
   canonical representation of rank R of the 2-tensor M_{ij}

   M_{ij} = \sum_{r=1..R}  S_{i1,j1}^r * ... *  S_{id,jd}^r

   Since the skeleton matrices correspond to indices (i1,j1) etc,
   these have also some rank. However, the (matrix) rank of these
   matrices is full (they are non-singular). Moreover (in the
   applications which we are interested in) the rank of these matrices
   have elements which are integrals, which are obtained by
   quadrature:

   S_{i1,j1} = \sum_{nodes x_n1 on [0,1]} w1_{n1} * w2_{n1} * B_{i1}(x_n1)  * B_{j1}^T(x_n1)

   This makes our final matrix being a huge guy with rank

   M_{ij} = \sum_{r=1..R} \sum_{nodes x_{n1,..,n_d}}

   w1_{n1} * .. * wd_{nd} # B_{i1}(x_n1) *..* B_{id}(x_nd) #  \omega^r(x1)*..* \omega^r(xd) # B_{j1}^T(x_n1) * ... * B_{j1}^T(x_n)



 */
template<class T, short_t d>
class gsKroneckerMatrix : public gsKroneckerMatrix<T,-1>
{
public:
    typedef gsFiberMatrix<T, ColMajor> FiberMatrix;
    typedef std::vector<FiberMatrix> cwiseMat;
    typedef gsKroneckerMatrixIter<T, d> iterator;

    typedef memory::unique_ptr<gsKroneckerMatrix> uPtr;

public:

    /// Initialize a low rank matrix by coordinate matrices
    gsKroneckerMatrix(cwiseMat cwise)
            : m_sz(cwise.size())
    {
        m_cwise.swap(cwise);
        for (size_t i = 0; i != m_cwise.size(); ++i)
            m_sz(i) = m_cwise[i].rows();
        m_rank = m_cwise.front().cols() / m_sz[0];// assuming square blocks

        m_dbd.resize(1);
        m_dbd(0) = 0;
    }

    gsKroneckerMatrix()
            : m_sz(0), m_rank(0)
    {}

    short_t dim() const { return d; }

    /// Order of the LR matrix (number of matrices in the Kronecker product)
    index_t order() const
    { return m_sz.size(); }

    /// Kronecker rank (number of Kronecker products in the sum)
    index_t krRank() const
    { return m_rank; }

    /// Rows of the global matrix
    index_t rows() const
    { return m_sz.prod(); }

    /// Columns of the global matrix
    index_t cols() const
    { return m_sz.prod(); }

    /// Size (rows or columns) of Kronecker factor matrices
    /// NOTE: the class assumes that the matrices are square
    const gsVector<index_t, d> &dimensions() const
    { return m_sz; }

private:// temporary memory for matrix-vector product
    mutable gsMatrix<T, -1, -1, ColMajor> q0, q1;
public:

    /// Computes the matrix-vector product M*in and writes the result in \out
    void apply(const gsMatrix <T> &in, gsMatrix <T> &out) const
    {
        static const index_t n = 1; // in.rows();
        GISMO_ASSERT(in.cols() == 1, "single rhs for now, got " << in.cols());

        const index_t sz = m_sz.prod();

        out.setZero(sz, n);
        for (index_t r = 0; r != m_rank; ++r)
        {
            q0 = in;

            for (unsigned i = 0; i != m_cwise.size(); ++i)
            {
                const index_t sz_i = m_cwise[i].rows();//cols
                const index_t r_i = sz / sz_i;
                const index_t rsz_i = r * sz_i;

                q0.resize(sz_i, n * r_i);
                q1.resize(r_i, n * sz_i);

                /*// Multiple rhs
                  for ( index_t k = 0; k!=n; ++k)
                  q1.middleCols(k*sz_i, sz_i).noalias() =
                  (  m_cwise[i].block(0, rsz_i, sz_i, sz_i)
                  * q0.middleCols(k*r_i, r_i) ).transpose();
                  */
                //TODO Add block to fiberMatrix
                //q1.noalias() = (m_cwise[i].block(0, rsz_i, sz_i, sz_i) * q0).transpose();

                q1.noalias() = (m_cwise[i].toSparseMatrix() // HIGHLY INEFFICIENT
                                .block(0, rsz_i, sz_i, sz_i) * q0).transpose();
                // complexity: (2p+1) sz_i^d

                q1.swap(q0);
            }
            q0.resize(sz, n);
            out += q0;
        }

        ///*
        // Hardcode all-Dirichlet boundaries as diagonal ones
        // destroys symmetry:
        for (index_t i = 0; i != m_dbd.rows(); ++i)
        {
            const index_t k = m_dbd.at(i);
            out.row(k) = in.row(k);
            // to regain symmetry:
            // add: for all non-dirichet entries j, substruct the matrix entry: out(j) -= Mat(m_dbd.at(i),j);
        }

    }

    //----------------
    gsMatrix<unsigned> m_dbd;
    //T penalty;
    //----------------


    /// Retuns the \a rowId, \a colId coefficient of the global matrix
    inline const T coeff(index_t rowId, index_t colId) const
    {
        gsVector<index_t, d> u(m_sz.size()), v(m_sz.size());
        for (size_t i = 0; i != m_cwise.size(); ++i)
        {
            rowId -= u[i] = rowId % m_sz[i];
            rowId /= m_sz[i];
            colId -= v[i] = colId % m_sz[i];
            colId /= m_sz[i];
        }
        return coeff(u, v);

        /*
        T result = 0;
        for (index_t r = 0; r != m_rank; ++r)
        {
            index_t row = rowId, col = colId;
            T tmp = 1;
            for (unsigned i = 0; i != m_cwise.size(); ++i)
            {
                const index_t cr = row % m_sz[i];
                row -= cr; row /= m_sz[i];
                const index_t cc = col % m_sz[i];
                col -= cc; col /= m_sz[i];

                tmp *= m_cwise[i](cr, r*m_sz[i] + cc );
            }
            result += tmp;
        }
        return result;
        */
    }

    /// Same as coeff(i,j) but indices given in tensor form
    inline const T coeff(const gsVector<index_t, d> &rowId,
                         const gsVector<index_t, d> &colId) const
    {
        T result = 0;
        for (index_t r = 0; r != m_rank; ++r)
        {
            T tmp = 1;
            for (unsigned i = 0; i != m_cwise.size(); ++i)
                tmp *= m_cwise[i](rowId[i], r * m_sz[i] + colId[i]);
            result += tmp;
        }
        return result;
    }

    /// Computes the \a rowId column of the global matrix and stores
    /// it in \a result
    void col_into(index_t rowId, gsSparseVector <T> &result) const
    {
        typedef typename gsSparseMatrix<T>::iterator cIter;

        const index_t n = (d == -1 ? m_sz.size() : d);// how to do static??
        gsVector<index_t, d> u(n), rstr(n);
        const index_t rkn = n * m_rank;
        std::vector<cIter> r(rkn);
        for (size_t i = 0; i != n; ++i)
        {
            rowId -= u[i] = rowId % m_sz[i];
            rowId /= m_sz[i];
            for (index_t k = 0; k != rkn; k += n)
                r[k + i] = m_cwise[i].begin((k / n) * m_sz[i] + u[i]);
        }
        const index_t nc = m_sz.prod();
        result.resize(nc);
        index_t nzPerCol = 1;
        for (index_t i = 0; i != n; ++i)
            nzPerCol *= bandWidth(m_cwise[i], m_sz[i]);
        result.reserve(nzPerCol);

        tensorStrides(m_sz, rstr);
        index_t t;
        do // for each row
        {
            index_t rr = r[0].index();// row index in the result
            T tmp = r[0].value();
            for (index_t i = 1; i != n; ++i)
            {
                tmp *= r[i].value();
                rr += rstr[i] * r[i].index();
            }
            for (index_t k = n; k != rkn; k += n)
            {
                T tmpr = r[k].value();
                for (index_t i = 1; i != n; ++i)
                    tmpr *= r[k + i].value();
                tmp += tmpr;
            }

            result.insert(rr) = tmp;

            // start next position
            for (t = 0; t != n; ++t)
            {
                // assuming common structure for all skeleton matrices
                for (index_t k = n; k != rkn; k += n) ++r[k + t];
                if (++r[t])
                    break;
                else
                {
                    if (t + 1 == n)
                    {
                        ++t;
                        break;
                    }
                    else
                        for (index_t k = 0; k != rkn; k += n)
                            r[k + t] = m_cwise[t].begin((k / n) * m_sz[t] + u[t]);
                }
            }
            // end next position

        } while (t != n);
        // result.makeCompressed();
    }

    /// Computes the diagonal of the global matrix
    void diagonal_into(gsVector <T> &result) const
    {
        result.setZero(m_sz.prod()); // the matrix diagonal

        T tmp;
        gsGridIterator<index_t, CUBE, d> it(m_sz);
        for (index_t c = 0; it; ++it, ++c)
        {
            for (index_t r = 0; r != m_rank; ++r)
            {
                tmp = 1;
                for (unsigned i = 0; i != m_cwise.size(); ++i)
                    tmp *= m_cwise[i].coeff(it->at(i), r * m_sz[i] + it->at(i));
                result[c] += tmp;
            }
        }

        ///*
        // Hardcode all-Dirichlet boundaries
        for (index_t i = 0; i != m_dbd.rows(); ++i)
            result[m_dbd.at(i)] = 1;
        //*/

    }

    FiberMatrix toFiberMatrix() const
    {
#ifdef DO_TIMING
        timer.restart();
#endif

        std::vector<index_t> sparsity;
        _reserveKroneckerProduct(sparsity);

#ifdef DO_TIMING
        pattern += timer.stop();
        timer.restart();
#endif

        FiberMatrix result(cols(), cols());
        result.reserve(sparsity);

#ifdef DO_TIMING
        reserve += timer.stop();
        timer.restart();
#endif

        _gsKroneckerProductSum(result);

#ifdef DO_TIMING
        kronecker += timer.stop();
#endif
        return result;
    }

    gsSparseMatrix<T> toSparseMatrix() const
    {
        return toFiberMatrix().toSparseMatrix();
    }

    struct wipeRowCol
    {
        wipeRowCol(index_t _s) : s(_s)
        {}

        // false --> remove, true --> keep
        bool operator()(const index_t &row, const index_t &col, real_t &value) const
        {
            return (row != 0 && col % s != 0) && (row + 1 != s && (col + 1) % s != 0);
        }

        const index_t s;
    };

    void setExtremeZeros()
    {
        for (unsigned i = 0; i != m_cwise.size(); ++i)
            m_cwise[i].prune(wipeRowCol(m_sz[i]));
    }

    void eliminateBoundary()
    {
        for (index_t i = 0; i != m_cwise.size(); i++)
        {

        }
    }


    const cwiseMat &skeleton() const
    { return m_cwise; }


private:
    cwiseMat m_cwise;
    //const cwiseMat & m_cwise;
    gsVector<index_t, d> m_sz;
    index_t m_rank;

// #define DO_INNER_TIMING
    void _gsKroneckerProductSum(FiberMatrix & result) const
    {
    /* // to test: is this faster ? --> currently NOT: make static dim!
    index_t n = 1;
    for( size_t i = 0; i!= m.size(); ++i)
    n *= m[i].rows();// row sizes (used for columns as well)

    //gsLRSparseMatrix<T,D> lrm(m);
    result.resize(n, n);
    reserveKroneckerProduct_V0(m, result);

    index_t r, c;
    for(typename gsLRSparseMatrix<T>::iterator it(m); it; ++it)
    {
        const T val = it.value(r,c);
        result.insert(r, c) = val;
    }
    //*/

    ///*
        // idea 1: make specific versions for 2D, 3D using for loops --> loop unrolling
        // idea 2: cache nested multiplications within the computation

        typedef typename gsFiberMatrix<T>::iterator cIter;
        const index_t n = (d==-1 ? m_cwise.size() : d);// how to do static??
        GISMO_ASSERT(n>1, "Need at least 2 matrices");

        gsVector<index_t,d> rsz(n), rstr(n);

        for( index_t i = 0; i!= n; ++i)
            rsz [i] = m_cwise[i].rows();// row sizes (used for columns as well)

        const index_t rk = m_cwise.front().cols() / rsz[0];// assuming square blocks
        const index_t rkn = n * rk;

        tensorStrides(rsz,rstr);

        gsVector<index_t,d> c(n);
        std::vector<cIter> r(rkn);
        //gsVector<cIter,D> r(rkn);
        c.setZero();
        index_t t;

        //gsGridIterator<cIter,CUBE,D,true> it(r);
        index_t cc = -1;

        do // for each column
        {
            // const index_t cc = c.dot(rstr); // rows==cols
            ++cc;

            for( index_t i = 0; i!=n; ++i)
                for( index_t k = 0; k!=rkn; k+=n)
                    r[k+i] = m_cwise[i].begin((k/n)*rsz[i] + c[i]);

            auto& fiber = result.fiber(cc);
            do // for each row
            {
                index_t rr = r[0].index();// row index in the result
                T tmp      = r[0].value();
                for( index_t i = 1; i!=n; ++i)
                {
                    tmp *= r[i].value();
                    rr  += rstr[i] * r[i].index();
                }
#ifdef DO_INNER_TIMING
                innerTimer.restart();
#endif
                for( index_t k = n; k!=rkn; k+=n)
                {
                    T tmpr = r[k].value();
                    for( index_t i = 1; i!=n; ++i)
                        tmpr *= r[k+i].value();
                    tmp += tmpr;
                }
#ifdef DO_INNER_TIMING
                product += innerTimer.stop();
#endif

                fiber.insertBack(rr) =  tmp;

                // GISMO_ASSERT(tmp - tmp ==tmp - tmp, "Inf");
                // GISMO_ASSERT(tmp == tmp, "Nan");

                // start next position
                for (t = 0; t != n; ++t)
                {
                    // assuming common structure for all skeleton matrices
                    for( index_t k = n; k!=rkn; k+=n) ++r[k+t];
                    if ( ++r[t] )
                        break;
                    else
                        if (t == n - 1)
                        { ++t; break;}
                        else
                            for( index_t k = 0; k!=rkn; k+=n)
                            r[k+t] = m_cwise[t].begin((k/n)*rsz[t] + c[t]);
                }
                // end next position

            } while (t!=n);

        } while (nextLexicographic(c, rsz));//assuming square blocks
//*/
    }
#undef DO_INNER_TIMING

    void _reserveKroneckerProduct(std::vector<index_t> & result) const
    {

        typedef typename gsFiberMatrix<T>::iterator cIter;
        const index_t n = (d==-1 ? m_cwise.size() : d);// how to do static??
        GISMO_ASSERT(n>1, "Need at least 2 matrices");

        gsVector<index_t,d> rsz(n);

        for( index_t i = 0; i!= n; ++i)
            rsz [i] = m_cwise[i].rows();// row sizes (used for columns as well)

        result.resize( rsz.prod(), 0);


        gsVector<index_t,d> c(n);
        std::vector<cIter> r(n);
        c.setZero();
        index_t t;

        //gsGridIterator<cIter,CUBE,D,true> it(r);
        index_t cc = -1;

        do // for each column
        {
            // const index_t cc = c.dot(rstr); // rows==cols
            ++cc;

            for( index_t i = 0; i!=n; ++i)
                r[i] = m_cwise[i].begin(c[i]);

            do // for each row
            {
                ++result[cc];
                // GISMO_ASSERT(tmp - tmp ==tmp - tmp, "Inf");
                // GISMO_ASSERT(tmp == tmp, "Nan");

                // start next position
                for (t = 0; t != n; ++t)
                {
                    if ( ++r[t] )
                        break;
                    else
                        if (t == n - 1)
                        { ++t; break;}
                        else
                            r[t] = m_cwise[t].begin(c[t]);
                }
                // end next position

            } while (t!=n);

        } while (nextLexicographic(c, rsz));//assuming square blocks
    }

};

// template< class T, int d = -1>
//     class gsLRMatrix : public gsLinearOperator


template<class T>
class gsKroneckerMatrix<T,-1> : public gsLinearOperator<T>
{
protected:
    index_t m_rank; ///< Hadamard rank of the matrix
public:
    typedef gsFiberMatrix<T> FiberMatrix;
    typedef std::vector<FiberMatrix> cwiseMat;
    typedef std::map<index_t, cwiseMat> blockColumnRep; ///< Contains block associated to a column
    typedef std::vector<blockColumnRep> blockRowColumnRep;
    typedef memory::unique_ptr<gsKroneckerMatrix > uPtr;

    virtual ~gsKroneckerMatrix() = default;

    virtual void apply(const gsMatrix <T> &in, gsMatrix <T> &out) const = 0;

    virtual void diagonal_into(gsVector <T> &result) const = 0;
            
    virtual short_t dim() const = 0;

    index_t rank() const { return m_rank;}

    virtual index_t cols() const = 0;
    virtual index_t rows() const = 0;
            
    gsSparseMatrix<T> toSparseMatrix() const
    { return toFiberMatrix().toSparseMatrix(); }

    virtual FiberMatrix toFiberMatrix() const
    {
        // static dispatch vs virtual function? what is faster here?
        switch( this->dim() )
        {
        case 2: return static_cast<const gsKroneckerMatrix<T,2>&>(*this).toFiberMatrix();
        case 3: return static_cast<const gsKroneckerMatrix<T,3>&>(*this).toFiberMatrix();
        case 4: return static_cast<const gsKroneckerMatrix<T,4>&>(*this).toFiberMatrix();
        default:
            GISMO_ERROR("Dimension error: "<< this->dim());
        }
    }

    virtual FiberMatrix toFiberMatrixBlockApproach() const
    {
        // static dispatch vs virtual function? what is faster here?
        switch( this->dim() )
        {
        case 2: return static_cast<const gsKroneckerMatrix<T,2>&>(*this).toFiberMatrixBlockApproach();
        case 3: return static_cast<const gsKroneckerMatrix<T,3>&>(*this).toFiberMatrixBlockApproach();
        case 4: return static_cast<const gsKroneckerMatrix<T,4>&>(*this).toFiberMatrixBlockApproach();
        default:
            GISMO_ERROR("Dimension error: "<< this->dim());
        }
    }


    virtual size_t gradingParameter() const
    {
        // static dispatch vs virtual function? what is faster here?
        switch( this->dim() )
        {
        case 2: return static_cast<const gsKroneckerMatrix<T,2>&>(*this).gradingParameter();
        case 3: return static_cast<const gsKroneckerMatrix<T,3>&>(*this).gradingParameter();
        case 4: return static_cast<const gsKroneckerMatrix<T,4>&>(*this).gradingParameter();
        default:
            GISMO_ERROR("Dimension error: "<< this->dim());
        }
    }

    virtual const cwiseMat & skeleton() const
    {
        // static dispatch vs virtual function? what is faster here?
        switch( this->dim() )
        {
        case 2: return static_cast<const gsKroneckerMatrix<T,2>&>(*this).skeleton();
        case 3: return static_cast<const gsKroneckerMatrix<T,3>&>(*this).skeleton();
        case 4: return static_cast<const gsKroneckerMatrix<T,4>&>(*this).skeleton();
        default:
            GISMO_ERROR("Dimension error: "<< this->dim());
        }
    }

    virtual const FiberMatrix & truncation() const
    {
        // static dispatch vs virtual function? what is faster here?
        switch( this->dim() )
        {
        case 2: return static_cast<const gsKroneckerMatrix<T,2>&>(*this).truncation();
        case 3: return static_cast<const gsKroneckerMatrix<T,3>&>(*this).truncation();
        case 4: return static_cast<const gsKroneckerMatrix<T,4>&>(*this).truncation();
        default:
            GISMO_ERROR("Dimension error: "<< this->dim());
        }
    }

};


/**
   Operator expressing a diagonal matrix

   Mostly used for Jacobi (diagonal) preconditionners
 */
template< class T>
class gsDiagonalOp : public gsLinearOperator<T>
{
public:

    typedef memory::shared_ptr<gsDiagonalOp> Ptr;
    typedef memory::unique_ptr<gsDiagonalOp> uPtr;

    /// Constructor by a vector (the diagonal entries
    gsDiagonalOp(const gsVector<T> & diag)
    : m_diag(diag)
    { }

    /// Make function returning a smart pointer
    static uPtr make(const gsVector<T> & diag)
    { return uPtr( new gsDiagonalOp(diag) ); }


    void apply(const gsMatrix<T> & in, gsMatrix<T> & out) const
    {
        out = m_diag.array() * in.array() ;
    }

    const gsVector<T> & vector() {return m_diag; }

    index_t rows() const {return m_diag.rows();}

    index_t cols() const {return m_diag.rows();}

private:
    const gsVector<T> & m_diag;
};



} // namespace gismo
