/** @file gsMatrixBlockView.h

    @brief Wraps a matrix and attaches a block structure to it.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

// Assuming that the Eigen library is included already

namespace gismo {

/** \brief Represents a block-view of the given matrix
 *
 * The blocks are references to the matrix segments. Each block
 * can be seen a a standalone matrix.
 *
 * \tparam T coefficient type
 * \ingroup Matrix
 */

// isConst .. to do
template <typename MatrixType, bool isConst = false>
class gsMatrixBlockView
{
public:
    typedef Eigen::Block<MatrixType>   block_t    ;
    typedef Eigen::Block<MatrixType> * block_ptr_t;

    typedef Eigen::Matrix<index_t,Eigen::Dynamic, 1, Eigen::ColMajor> Vector_t;

public:

    gsMatrixBlockView() : m_rowSize(0), m_colSize(0)
    { }

    /** \brief Creates a block-view of the given matrix
     *
     * The blocks are references to the matrix segments. Each block
     * can be seen a a standalone matrix.
     *
     * \param[in] matrix   the matrix to be "partitioned" into blocks
     * \param[in] rowSizes the sizes of the row blocks, must sum up to matrix.rows()
     * \param[in] colSizes the sizes of the column blocks, must sum up to matrix.cols()
     */
    gsMatrixBlockView(MatrixType & matrix,
                      const Vector_t & rowSizes,
                      const Vector_t & colSizes)
    : m_rowSize(rowSizes.size()),
      m_colSize(colSizes.size())
    {
        GISMO_ASSERT( rowSizes.sum() == matrix.rows() &&
                      colSizes.sum() == matrix.cols() ,
                      "Invalid block structure.");

        // NB: SparseMatrix blocks are not assignable, therefore we
        // use pointers
        block_ptr_t tmp;

        m_blocks.reserve(m_rowSize*m_colSize);
        index_t row, col(0);

        for ( index_t j = 0; j<m_colSize; ++j )
        {
            row = 0;
            for ( index_t i = 0; i<m_rowSize; ++i )
            {
                tmp = new block_t( matrix.block(
                                   row        , col        ,
                                   rowSizes[i], colSizes[j])
                                 );

                m_blocks.push_back(tmp);

                row += rowSizes[i];
            }
            col += colSizes[j];
        }
    }

    /** \brief Creates a block-view of the given matrix, using only one column piece
     *
     * This constructor is suitable for vectors,
     * The blocks are references to the matrix segments. Each block
     * can be seen a a standalone matrix.
     *
     * \param[in] matrix the matrix to be "partitioned" into blocks
     * \param[in] rowSizes the sizes of the row blocks, must sum up to matrix.rows()
     */
    gsMatrixBlockView(MatrixType & matrix,
                      const Vector_t & rowSizes)
    : m_rowSize(rowSizes.size()),
      m_colSize(1)
    {
        GISMO_ASSERT( rowSizes.sum() == matrix.rows(),
                      "Invalid block structure.");

        // NB: SparseMatrix blocks are not assignable, therefore we
        // use pointers
        block_ptr_t tmp;

        m_blocks.reserve(m_rowSize*m_colSize);
        const index_t cols = matrix.cols();
        index_t row(0);

        for ( index_t i = 0; i<m_rowSize; ++i )
        {
            tmp = new block_t( matrix.block(row, 0, rowSizes[i], cols ) );

            m_blocks.push_back(tmp);

            row += rowSizes[i];
        }
    }

    /// Combatibility constuctor giving the \a matrix as a block
    gsMatrixBlockView(const MatrixType & matrix)
    : m_rowSize(1),
      m_colSize(1)
    {
        m_blocks.push_back( new block_t(matrix.topRows(matrix.rows())) );
    }

    /// Copy constructor
    gsMatrixBlockView (const gsMatrixBlockView & other)
    {
        freeAll(m_blocks);
        m_blocks.clear();
        m_rowSize = other.m_rowSize;
        m_colSize = other.m_colSize;

        block_ptr_t tmp;

        for ( typename std::vector<block_ptr_t>::const_iterator it
                  = other.m_blocks.begin(); it != other.m_blocks.end(); ++it )
        {
                tmp = new block_t(**it);
                m_blocks.push_back(tmp);
        }
    }

    ~gsMatrixBlockView()
     {
         freeAll(m_blocks);
     }

public:

    /// Returns the number of blocks
    size_t numBlocks()    const { return m_blocks.size(); }

    /// Returns the number of row-blocks
    size_t numRowBlocks() const { return m_rowSize;       }

    /// Returns the number of col-blocks
    size_t numColBlocks() const { return m_colSize;       }

    /// Returns the block indexed \a i (row) and \a j (column)
    block_t & operator () ( index_t i, index_t j = 0) const
    {
        GISMO_ASSERT( i < m_rowSize && j < m_colSize ,
                      "Invalid block requested.");

        return *m_blocks[j*m_rowSize+i];
    }


    gsMatrixBlockView & operator = (gsMatrixBlockView  other)
    {
        if(this != &other)
        {
            m_blocks.swap(other.m_blocks);
            m_rowSize = other.m_rowSize;
            m_colSize = other.m_colSize;
        }

        return *this;
    }

    /// Overwrites the contents of block (i,j) with matrix \a other
	template<typename OtherDerived>
    void assign(index_t i, index_t j, const Eigen::EigenBase<OtherDerived>& other)
    {
        GISMO_ASSERT( i < m_rowSize && j < m_colSize ,
                      "Assign to invalid block requested.");

        // Works for dense Eigen matrices
        *m_blocks[j*m_colSize+i] = other;
    }

    /// Prints the block structure
    friend std::ostream & operator << (std::ostream & os, const gsMatrixBlockView & mv)
    {
        os<< "Matrix block-view size: "<< mv.m_rowSize <<"x"<< mv.m_colSize
          <<" blocks. Structure:\n";

        for ( index_t i = 0; i<mv.m_rowSize; ++i )
        {
            for ( index_t j = 0; j<mv.m_colSize; ++j )
            {
                block_t & bl = mv(i,j);
                os<< bl.rows() <<"x"<<bl.cols() <<"  ";
            }
            os<<"\n";
        }

        return os;
    }

private:

    std::vector<block_ptr_t> m_blocks;

    index_t m_rowSize;

    index_t m_colSize;

};





} // namespace gismo
