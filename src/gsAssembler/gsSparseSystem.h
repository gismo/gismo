/** @file gsSparseSystem.h

    @brief Class representing a sparse linear system (with dense
    right-hand side(s)), indexed by sets of degrees of freedom

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, C. Hofer, K. Rafetseder
*/

#pragma once

#include <gsCore/gsStdVectorRef.h>

namespace gismo
{

/**
    @brief A sparse linear system indexed by sets of degrees of
    freedom

    \ingroup Assembler
*/
template<typename T, bool symm>
class gsSparseSystem
{
protected:
    typedef gsStdVectorRef<gsDofMapper> DofMapperRef;

    typedef std::vector<gsDofMapper>    DofMappers;

    typedef typename gsSparseMatrix<T>::BlockView matBlockView;

    typedef typename gsMatrix<T>::BlockView rhsBlockView;

protected:

    /// @brief the system matrix
    gsSparseMatrix<T> m_matrix;

    /// @brief the right hand side of the system
    gsMatrix<T> m_rhs;


    // -- Structure

    /// @brief the mappers for the different blocks. The assignment of the
    /// mappers to the column (and row) blocks is done via the members
    /// \a m_row and \a m_col. This also allows for different mappers for columns
    /// and rows. There are no assumptions on the ordering of \a m_mappers.
    ///  The relation of column mappers and the used multibasis is done via \a m_cvar.
    DofMappers m_mappers;

    /// @brief map between row blocks and index of \a m_mappers, i.e.
    ///  row block i is described by m_mappers[m_row[i]].
    gsVector<index_t> m_row;

    /// @brief map between column blocks and index of \a m_mappers i.e.
    ///  col block j is described by m_mappers[m_col[j]].
    gsVector<index_t> m_col;

    /// @brief strides for the row blocks (shifting of mapped indices).
    /// The mapper do not have this information anymore.
    gsVector<index_t> m_rstr;

    /// @brief strides for the column blocks (shifting of mapped indices).
    /// The mapper do not have this information anymore.
    gsVector<index_t> m_cstr;

    /// @brief map between column blocks and used multibasis
    /// (which unknown/component correlate to which multibasis), i.e. column block i
    /// uses multibasis m_cvar[i]. So this allows e.g. a single multibasis for several components.
    gsVector<index_t> m_cvar;

    gsVector<index_t> m_dims;

public:

    gsSparseSystem()
    { }

    /**
     * @brief gsSparseSystem Constuctor for the sparse System only considering one block (scalar equation),
     * described by the one DofMapper \a mapper.
     * @param[in] mapper the one DofMapper discribing the discretization.
     */
    gsSparseSystem(gsDofMapper & mapper)
        : m_mappers(1),
          m_row    (1),
          m_col    (1),
          m_rstr   (1),
          m_cstr   (1),
          m_cvar   (1),
          m_dims   (1)
    {
        m_row [0] =  m_col [0] =
                m_rstr[0] =  m_cstr[0] =
                m_cvar[0] = 0;

        m_dims[0] = 1;

        m_mappers.front().swap(mapper);

        m_matrix.resize( m_mappers.front().freeSize() ,
                         m_mappers.front().freeSize() );
    }


    /**
     * @brief gsSparseSystem Constructor of a sparse System specified by the number of unknows for each
     * block, given by \a dims. It is assumed that the used Basis (m_bases in gsAssembler) has a one to
     * one relation with the given \a mappers.
     * E.g.: dims = {1,3,1,2}, then we consider a matrix with 1+3+1+2 = 7 x 7 blockstructure.  Then
     * a) mappers.size() == 7: In the first case, row and column mappers are identical, so there is no
     *    difference bettween column and row blocks.
     *    - Column block 0     uses mapper 0 and corresponds to m_bases[0]
     *    - Column block 1,2,3 uses mapper 1 and corresponds to m_bases[1]
     *    - Column block 4     uses mapper 2 and corresponds to m_bases[2]
     *    - Column block 5,6   uses mapper 3 and corresponds to m_bases[3]     (in C++ indexing)
     *    Identical maps for Row blocks
     * b) mappers.size() == 14: In the second case, row and column mappers are not identical.
     *    - Row block 0        uses mapper 0
     *    - Row block 1,2,3    uses mapper 1
     *    - Row block 4        uses mapper 2
     *    - Row block 5,6      uses mapper 3
     *    - Column block 0     uses mapper 4 and corresponds to m_bases[0]
     *    - Column block 1,2,3 uses mapper 5 and corresponds to m_bases[1]
     *    - Column block 4     uses mapper 6 and corresponds to m_bases[2]
     *    - Column block 5,6   uses mapper 7 and corresponds to m_bases[3]     (in C++ indexing)
     * Other cases are not handled!
     * @param mappers a vector of gsDofMapper with size == dims.sum() or 2*dims.sum() according to the
     *                definition above
     * @param dims Defines how many unknown are determined by a certain dofMapper.
     */
    gsSparseSystem(DofMappers & mappers,
                   const gsVector<index_t> & dims)
        : m_row(dims.sum()),
          m_col(dims.sum()),
          m_rstr(dims.sum()),
          m_cstr(dims.sum()),
          m_dims(dims.cast<index_t>())
    {
        const index_t d = dims.size();
        const index_t s = dims.sum();
        const index_t ms = mappers.size();

        GISMO_ASSERT(ms==s ||ms==2*s, "Connot deduce block structure");

        m_mappers.swap(mappers);

        //Calculate the map for blocks to mappers
        index_t k=0;
        for(index_t i=0;i<d;++i)
            for(index_t j=0; j<dims[i];++j)
            {
                m_row[k]=i;
                ++k;
            }
        //At first glance, row blocks and column blocks have the same mapper
        m_col = m_row;

        if (ms == 1 )
        {
            m_row.setZero();
            m_col.setZero();
            m_cvar.setZero(1);
        }
        else if ( ms == 2*s )
            m_col.array() += s; //mappers s+1... 2s are then the column mappers

        //assumes that the mappers are ordered as the bases in gsAssembler and are starting from m_bases[0]!
        m_cvar = m_row.cast<index_t>();

        //Fill the strides
        m_rstr[0] = m_cstr[0] = 0;
        for (index_t r = 1; r < d; ++r) // for all row-blocks
            m_rstr[r] = m_rstr[r-1] + m_mappers[m_row[r-1]].freeSize();
        for (index_t c = 1; c < d; ++c) // for all col-blocks
            m_cstr[c] = m_cstr[c-1] + m_mappers[m_col[c-1]].freeSize();

        m_matrix.resize( m_rstr.at(d-1) + m_mappers[m_row[d-1]].freeSize() ,
                m_cstr.at(d-1) + m_mappers[m_col[d-1]].freeSize() );

    }

    /**
     * @brief gsSparseSystem Constructor for the sparse system, with given number of row and column blocks.
     *        It is assumed that the mappers have a one to one correspondence with the blocks, i.e.
     *        a) mappers.size() == rows + cols ==> the first \a row entries of \a mappers correspond to
     *           the row blocks and the last \a cols entries of \a mappers correspond to the column blocks
     *           Moreover, it is assumes that the basis for the solution space (column blocks) are indexed
     *           starting by 0. i.e. mappers[cols],...,mappers[rows+cols-1] correspond in a one to one relation
     *           to m_bases[0],...,m_bases[cols-1] in gsAssembler!!!!
     *        b) mappers.size() == rows == cols ==> row and column mappers are identical, and
     *           \a mappers[i] is the mapper for i-th column (and row) block
     *        c) mappers.size() == 1 ==> only one block
     * @param mappers the set of mappers, with restrictions above
     * @param rows the number of row blocks
     * @param cols the number of column blocks
     */
    gsSparseSystem(DofMappers & mappers,
                   const index_t rows,
                   const index_t cols)
        : m_row (gsVector<index_t>::LinSpaced(rows,0,rows-1)),
          m_col (gsVector<index_t>::LinSpaced(cols,0,cols-1)),
          m_rstr(rows),
          m_cstr(cols),
          m_dims(cols)
    {
        GISMO_ASSERT( rows > 0 && cols > 0, "Block dimensions must be positive");

        m_mappers.swap(mappers);

        if ( static_cast<index_t>(m_mappers.size()) == rows + cols )
        {
            m_col.array() += rows;
            m_cvar = m_row.cast<index_t>();//assume 1-1 bases
        }
        else if ( static_cast<index_t>(m_mappers.size()) == rows )
        {
            GISMO_ENSURE( rows == cols, "Dof Mapper vector does not match block dimensions");
            m_cvar = m_row.cast<index_t>();//assume 1-1 bases
        }
        else if ( m_mappers.size() == 1 )
        {
            m_row.setZero();
            m_col.setZero();
            m_cvar.setZero(1);
        }
        else
        {
            GISMO_ERROR("Cannot deduce block structure.");
        }

        m_dims.setOnes();

        m_rstr[0] = m_cstr[0] = 0;
        for (index_t r = 1; r < rows; ++r) // for all row-blocks
            m_rstr[r] = m_rstr[r-1] + m_mappers[m_row[r-1]].freeSize();
        for (index_t c = 1; c < cols; ++c) // for all col-blocks
            m_cstr[c] = m_cstr[c-1] + m_mappers[m_col[c-1]].freeSize();

        m_matrix.resize( m_rstr.at(rows-1) + m_mappers[m_row[rows-1]].freeSize() ,
                m_cstr.at(cols-1) + m_mappers[m_col[cols-1]].freeSize() );
    }

    /**
     * @brief gsSparseSystem The constructor for the sparse system given the set of mappers, may be different
     * for each column and row block, the assignment of mappers to row blocks and column block, and the relation
     * which Multibases correspond to which unknown. This constructor allows for the most possible freedome to
     * design your block system.
     *
     * What is not possible here, cannot be done!
     *
     * a fancy example: rowInd = {1,1,0,3}, colInd={0,2,2}, colvar={1,0,0}.
     * this leads to a 4x3 block structure matrix
     *
     *    - Row block 0,1      uses mapper 1
     *    - Row block 2        uses mapper 0
     *    - Row block 3        uses mapper 3
     *    - Column block 0     uses mapper 0 and corresponds to m_bases[1]
     *    - Column block 1,2   uses mapper 2 and corresponds to m_bases[0]   (in C++ indexing)
     *
     * @param[in] mappers the set of mappers, need not be a one to one relation with the blocks.
     * @param[in] rowInd assignment of row blocks to mappers, e.g. row block i uses mapper rowInd[i] of the
     *            set \a mappers.
     * @param[in] colInd assignment of column blocks to mappers, e.g. column block j uses mapper colInd[j] of the
     *            set \a mappers. (need not be the same mappers as for the row blocks)
     * @param[in] colvar assignment of unknowns to the used bases (i.e. column blocks to MultiBasis)
     */
    gsSparseSystem(DofMappers & mappers,
                   const gsVector<index_t> & rowInd,
                   const gsVector<index_t> & colInd,
                   const gsVector<index_t> & colvar)
        : m_row (rowInd),
          m_col (colInd),
          m_rstr((index_t)rowInd.size()),
          m_cstr((index_t)colInd.size()),
          m_dims(colInd.size())
        // ,m_cvar(colvar) //<< Bug
    {
        m_dims.setOnes();
        m_cvar = colvar;
        const index_t rows = m_row.size();
        const index_t cols = m_col.size();
        GISMO_ASSERT( rows > 0 && cols > 0,
                      "Block dimensions must be positive");

        m_mappers.swap(mappers);

        m_rstr[0] = m_cstr[0] = 0;
        for (index_t r = 1; r < rows; ++r) // for all row-blocks
            m_rstr[r] = m_rstr[r-1] + m_mappers[m_row[r-1]].freeSize();
        for (index_t c = 1; c < cols; ++c) // for all col-blocks
            m_cstr[c] = m_cstr[c-1] + m_mappers[m_col[c-1]].freeSize();

        m_matrix.resize( m_rstr.at(rows-1) + m_mappers[m_row[rows-1]].freeSize() ,
                m_cstr.at(cols-1) + m_mappers[m_col[cols-1]].freeSize() );
    }

    /**
     * @brief swap swaps the content of the Sparse System with the other given one
     * @param[in] other the other sparse system
     */
    void swap(gsSparseSystem & other)
    {
        m_matrix .swap(other.m_matrix );
        m_rhs    .swap(other.m_rhs    );
        m_mappers.swap(other.m_mappers);
        m_row    .swap(other.m_row    );
        m_col    .swap(other.m_col    );
        m_rstr   .swap(other.m_rstr   );
        m_cstr   .swap(other.m_cstr   );
        m_cvar   .swap(other.m_cvar   );
        m_dims   .swap(other.m_dims   );
    }

    /**
     * @brief reserve reserves the memory for the sparse matrix and the rhs.
     * @param[in] nz Non-zeros per column for the sparse matrix
     * @param [in] numRhs number of columns
     */
    void reserve(const index_t nz, const index_t numRhs)
    {
        GISMO_ASSERT( 0 != m_mappers.size(), "Sparse system was not initialized");
        if ( 0 != m_matrix.cols() )
        {
            m_matrix.reservePerColumn(nz);
            if ( 0 != numRhs )
                m_rhs.setZero(m_matrix.cols(), numRhs);
        }
    }

    /**
     * @brief Reserves the memory for the sparse matrix and the rhs,
     * based on the polynomial degree of the first basis-piece in
     * /em mb, as well as the input options bdA, bdB and bdO
     *
     * At each column approximately bdA * deg + dbB non-zero entries
     * are expected. An extra amount of memory of bdO percent is
     * allocated, in order to speedup the process.
     * @param mb
     * @param opt
     * @param [in] numRhs number of columns
     */
    void reserve(const gsMultiBasis<T> & mb, const gsOptionList & opt,
                 const index_t numRhs)
    {
        reserve(numColNz(mb,opt), numRhs);
    }

    /// @brief Provides an estimation of the number of non-zero matrix
    /// entries per column. This value can be used for sparse matrix
    /// memory allocation
    index_t numColNz(const gsMultiBasis<T> & mb, const gsOptionList & opt) const
    {
        // Pick up values from options
        const T bdA       = opt.getReal("bdA");
        const index_t bdB = opt.getInt("bdB");
        const T bdO       = opt.getReal("bdO");
        const gsBasis<T> & b = mb[0];
        T nz = 1;
        for (short_t i = 0; i != b.dim(); ++i)
            nz *= bdA * b.degree(i) + bdB;
        return cast<T,short_t>(nz*(1.0+bdO));
    }

    /// @brief set everything to zero
    void setZero()
    {
        m_matrix.setZero();
        m_rhs   .setZero();
    }

    /// @brief the number of matrix columns
    index_t cols() const { return m_matrix.cols(); }

    /// @brief the number of matrix rows
    index_t rows() const { return m_matrix.rows(); }

    /// @brief the rows of the right-hand side vector
    index_t rhsRows() const { return m_rhs.rows(); }

    /// @brief the number of right-hand side vectors
    index_t rhsCols() const { return m_rhs.cols(); }

public: /* Accessors */

    /// @brief Access the system Matrix
    const gsSparseMatrix<T> & matrix() const
    { return m_matrix; }

    /// @brief Access the system Matrix
    gsSparseMatrix<T> & matrix()
    { return m_matrix; }

    /// @brief Access the right hand side
    const gsMatrix<T> & rhs() const
    { return m_rhs; }

    /// @brief Access the right hand side
    gsMatrix<T> & rhs()
    { return m_rhs; }

    /// @brief returns a block view of the matrix, easy way to extract single blocks
    matBlockView blockView()
    {
        gsVector<index_t> rowSizes(m_row.size()), colSizes(m_col.size());

        for (index_t r = 0; r != rowSizes.size(); ++r) // for all row-blocks
            rowSizes[r] = m_mappers[m_row[r]].freeSize();

        for (index_t c = 0; c != colSizes.size(); ++c) // for all col-blocks
            colSizes[c] = m_mappers[m_col[c]].freeSize();

        return m_matrix.blockView(rowSizes,colSizes);
    }

    /**
     * @brief returns a block view of the matrix, where you can choose which blocks (of gsSparseSystem) should be combined to one block
     * @param[in] numRowBlocksNew number of row blocks
     * @param[in] numColBlocksNew number of column blocks
     * @param[in] rowBlocksNew a vector defining the row blocks
     * @param[in] colBlocksNew a vector defining the column blocks
     *
     * example:
     * assume in gsSparseSystem we have 4 row blocks
     *
     * rowBlocksNew = {0,1,1,2}
     * yields a block system with only 3 row blocks, where the blocks 1 and 2 of gsSparseSystem are combined to one block
     */


    matBlockView blockView(size_t numRowBlocksNew, size_t numColBlocksNew, const gsVector<index_t>& rowBlocksNew, const gsVector<index_t>& colBlocksNew)
    {
        gsVector<index_t> rowSizes(numRowBlocksNew), colSizes(numColBlocksNew);
        rowSizes.setZero();
        colSizes.setZero();

        for (index_t r = 0; r != m_row.size(); ++r) // for all row-blocks
            rowSizes[rowBlocksNew[r]] += m_mappers[m_row[r]].freeSize();

        for (index_t c = 0; c != m_col.size(); ++c) // for all col-blocks
            colSizes[colBlocksNew[c]] += m_mappers[m_col[c]].freeSize();

        return m_matrix.blockView(rowSizes, colSizes);
    }

    /**
     * @brief returns a block view of the rhs, where you can choose which blocks (of gsSparseSystem) should be combined to one block
     * @param[in] numRowBlocksNew number of row blocks
     * @param[in] rowBlocksNew a vector defining the row blocks
     */

    rhsBlockView blockViewRhs(size_t numRowBlocksNew, const gsVector<index_t>& rowBlocksNew)
    {
        gsVector<index_t> rowSizes(numRowBlocksNew), colSizes(1);
        rowSizes.setZero();
        colSizes[0]=1;

        for (index_t r = 0; r != m_row.size(); ++r) // for all row-blocks
            rowSizes[rowBlocksNew[r]] += m_mappers[m_row[r]].freeSize();

        return m_rhs.blockView(rowSizes, colSizes);
    }

    /// @brief returns the number of column blocks
    index_t numColBlocks() const {return m_col.size();}

    /// @brief returns the number of row blocks
    index_t numRowBlocks() const {return m_row.size();}

    /**
     * @brief rowMapper returns the dofMapper for row block \a r
     * @param[in] r the index of the row block
     * @return the corresponding DofMapper
     */
    const gsDofMapper & rowMapper(const index_t r) const
    { return m_mappers[m_row[r]]; }

    /**
     * @brief rowMapper returns the dofMapper for row block \a r
     * @param[in] r the index of the row block
     * @return the corresponding DofMapper
     */
    gsDofMapper & rowMapper(const index_t r)
    { return m_mappers[m_row[r]]; }

    /**
     * @brief colMapper returns the dofMapper for column block \a c
     * @param[in] c the index of the column block
     * @return the corresponding DofMapper
     */
    const gsDofMapper & colMapper(const index_t c) const
    {
        return m_mappers[m_col[c]];
    }

    /**
     * @brief colMapper returns the dofMapper for column block \a c
     * @param[in] c the index of the column block
     * @return the corresponding DofMapper
     */
    gsDofMapper & colMapper(const index_t c)
    { return m_mappers[m_col[c]]; }

    /**
     * @brief colBasis returns the index of the Basis used for column block \a c
     * @param[in] c the index of the column block
     * @return the index of the used basis
     */
    index_t colBasis(const index_t c) const // better name ?
    { return m_cvar[c]; }

    /// @brief return the number of components for the given component
    index_t unkSize(const index_t unk) const
    {return m_dims[unk];}

    /// @brief returns the number of unknowns
    index_t numUnknowns() const {return m_dims.size(); }

    /// @brief returns all dof Mappers.
    /// \note the result is not a one to one relation with the blocks.
    const DofMappers & dofMappers() const
    { return m_mappers; }

    /// @brief checks if the system was initialized
    bool initialized() const
    {
        return 0 != m_col.size();
    }

    /// @brief returns true if only half of the matrix is stored, due to its symmetry
    bool symmetry() const
    {
        return symm;
    }

    /*
    gsRefVector<gsDofMapper> colMappers()
    {
        return gsRefVector<gsDofMapper>(m_mappers, m_col);
    }

    gsRefVector<gsDofMapper> rowMappers()
    {
        return gsRefVector<gsDofMapper>(m_mappers, m_row);
    }
    */


public: /* mapping patch-local to global indices */

    /**
     * @brief mapRowIndices Maps a set of basis indices by the corresponding dofMapper.
     * \note that the result is not the position in the sparse system, since the shifts
     * are not included.
     * @param[in] actives the set of basis indices
     * @param[in] patchIndex the patch under consideration
     * @param[out] result the mapped indices
     * @param[in] r the considered row block
     */
    void mapRowIndices(const gsMatrix<index_t> & actives,
                       const index_t patchIndex,
                       gsMatrix<index_t> & result,
                       const index_t r = 0) const
    {
        m_mappers[m_row.at(r)].localToGlobal(actives, patchIndex, result);
    }

    /**
     * @brief mapColIndices Maps a set of basis indices by the corresponding dofMapper.
     * \note that the result is not the position in the sparse system, since the shifts
     * due to the block structure are not included.
     * @param[in] actives the set of basis indices
     * @param[in] patchIndex the patch under consideration
     * @param[out] result the mapped indices
     * @param[in] c the considered column block
     */
    void mapColIndices(const gsMatrix<index_t> & actives,
                       const index_t patchIndex,
                       gsMatrix<index_t> & result,
                       const index_t c = 0) const
    {
        m_mappers[m_col.at(c)].localToGlobal(actives, patchIndex, result);
    }

    /**
     * @brief mapToGlobalRowIndex Maps a single basis index to the final position in the system,
     * i.e. including the right shift due to the block structure
     * @param[in] active the basis index
     * @param[in] patchIndex the patch index
     * @param[out] result the mapped index
     * @param[in] r the considered row block
     */
    void mapToGlobalRowIndex(const index_t active,
                             const index_t patchIndex,
                             index_t & result,
                             const index_t r = 0) const
    {
        result = m_mappers[m_row.at(r)].index(active, patchIndex)+m_rstr[r];
    }

    /**
     * @brief mapToGlobalColIndex Maps a single basis index to the final position in the system,
     * i.e. including the right shift due to the block structure
     * @param[in] active the basis index
     * @param[in] patchIndex the patch index
     * @param[out] result the mapped index
     * @param[in] c the considered column block
     */
    void mapToGlobalColIndex(const index_t active,
                             const index_t patchIndex,
                             index_t & result,
                             const index_t c = 0) const
    {
        result = m_mappers[m_col.at(c)].index(active, patchIndex)+m_cstr[c];
    }

public: /* Add local contributions to system matrix */

    /**
     * @brief pushToMatrix pushes the local matrix for an element to the global system,
     * \note
     * 1. the same index set is assumed for row and column block
     * 2. eliminated dofs are incorporated in the right way
     * 3. no assembling is done for the rhs
     * @param[in] localMat the local matrix
     * @param[in] actives the mapped index of basis functions, without shifts!
     * @param[in] eliminatedDofs the values for the dofs, which are removed from the system
     * @param[in] r the row block
     * @param[in] c the column block
     */
    void pushToMatrix(const gsMatrix<T>  & localMat,
                      const gsMatrix<index_t> & actives,
                      const gsMatrix<T> & eliminatedDofs,
                      const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive = actives.rows();
        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];
        GISMO_ASSERT( &rowMap == &m_mappers[m_col.at(c)], "Error");

        for (index_t i = 0; i != numActive; ++i)
        {
            const index_t ii = m_rstr.at(r) + actives.at(i); // N_i

            if ( rowMap.is_free_index(actives.at(i)) )
            {
                for (index_t j = 0; j != numActive; ++j)
                {
                    const index_t jj = m_cstr.at(c) + actives.at(j); // N_j

                    if ( rowMap.is_free_index( actives.at(j)) )
                    {
                        // If matrix is symmetric, we store only lower
                        // triangular part
                        if ( (!symm) || jj <= ii )
                            m_matrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                    else if(0!=eliminatedDofs.size())
                    {
                        m_rhs.row(ii).noalias() -= localMat(i, j) *
                            eliminatedDofs.row( rowMap.global_to_bindex(actives.at(j)) );
                    }
                }
            }
        }
    }

    /**
     * @brief pushToMatrix pushes the local matrix for an element to the global system,
     * \note
     * 1. different index sets are used for row and column block
     * 2. eliminated dofs are incorporated in the right way
     * 3. no assembling is done for the rhs
     * @param[in] localMat the local matrix
     * @param[in] actives_i the mapped index of row - basis functions, without shifts!
     * @param[in] actives_j the mapped index of column - basis functions, without shifts!
     * @param[in] eliminatedDofs_j the values for the dofs (corresponding to the columns), which are
     *            removed from the system
     * remark: if no column dofs are eliminated (only row dofs are eliminated), you can just give an empty matrix, because it will not be used
     * @param[in] r the row block
     * @param[in] c the column block
     */
    void pushToMatrix(const gsMatrix<T>  & localMat,
                      const gsMatrix<index_t> & actives_i,
                      const gsMatrix<index_t> & actives_j,
                      const gsMatrix<T> & eliminatedDofs_j,
                      const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive_i = actives_i.rows();
        const index_t numActive_j = actives_j.rows();

        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];
        const gsDofMapper & colMap = m_mappers[m_col.at(c)];

        for (index_t i = 0; i != numActive_i; ++i)
        {
            const index_t ii = m_rstr.at(r) + actives_i.at(i); // N_i

            if ( rowMap.is_free_index(actives_i.at(i)) )
            {
                for (index_t j = 0; j != numActive_j; ++j)
                {
                    const index_t jj = m_cstr.at(c) + actives_j.at(j); // N_j

                    if ( colMap.is_free_index(actives_j.at(j)) )
                    {
                        // If matrix is symmetric, we store only lower
                        // triangular part
                        if ( (!symm) || jj <= ii )
                            m_matrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                    else
                    {
                        m_rhs.row(ii).noalias() -= localMat(i, j) *
                                eliminatedDofs_j.row( colMap.global_to_bindex(actives_j.at(j)) );
                    }
                }
            }
        }
    }


    /**
     * @brief pushToMatrixAllFree pushes the local matrix for an element to the global system,
     * \note
     * 1. the same index set is assumed for row and column block
     * 2. no checks are done if an index is eliminated or not
     * 3. no assembling is done for the rhs
     * @param[in] localMat the local matrix
     * @param[in] actives the mapped index of basis functions, without shifts!
     * @param[in] r the row block
     * @param[in] c the column block
     */
    void pushToMatrixAllFree(const gsMatrix<T>  & localMat,
                             const gsMatrix<index_t> & actives,
                             const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive = actives.rows();
        GISMO_ASSERT( &m_mappers[m_row.at(r)] == &m_mappers[m_col.at(c)], "Error");

        for (index_t i = 0; i != numActive; ++i)
        {
            const index_t ii = m_rstr.at(r) + actives.at(i); // N_i

            for (index_t j = 0; j != numActive; ++j)
            {
                const index_t jj = m_cstr.at(c) + actives.at(j); // N_j

                // If matrix is symmetric, we store only lower
                // triangular part
                if ( (!symm) || jj <= ii )
                    m_matrix.coeffRef(ii, jj) += localMat(i, j);


            }
        }
    }


    /**
     * @brief pushToMatrixAllFree pushes the local matrix for an element to the global system,
     * \note
     * 1. different index sets are used for row and column block
     * 2. no checks are done if an index is eliminated or not
     * 3. no assembling is done for the rhs
     * @param[in] localMat the local matrix
     * @param[in] actives_i the mapped index of row - basis functions, without shifts!
     * @param[in] actives_j the mapped index of column - basis functions, without shifts!
     * @param[in] r the row block
     * @param[in] c the column block
     */
    void pushToMatrixAllFree(const gsMatrix<T>  & localMat,
                             const gsMatrix<index_t> & actives_i,
                             const gsMatrix<index_t> & actives_j,
                             const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive_i = actives_i.rows();
        const index_t numActive_j = actives_j.rows();

        for (index_t i = 0; i != numActive_i; ++i)
        {
            const index_t ii = m_rstr.at(r) + actives_i.at(i); // N_i

            for (index_t j = 0; j != numActive_j; ++j)
            {
                const index_t jj = m_cstr.at(c) + actives_j.at(j); // N_j

                // If matrix is symmetric, we store only lower
                // triangular part
                if ( (!symm) || jj <= ii )
                    m_matrix.coeffRef(ii, jj) += localMat(i, j);
            }
        }
    }

    /**
     * @brief pushToMatrix pushes one local matrix consisting of several blocks corresponding to blocks of the global system
     * \note
     * 1. Usefull for bilinear forms depending on vector valued functions
     * 2. different index sets are used for row and column blocks
     * 3. eliminated dofs are incorporated in the right way
     * 4. assume identical row and column mappers for the global system, therefore only one vector of mapped index sets is given
     * @param[in] localMat local system matrix
     * @param[in] actives_vec a vector of mapped index sets (for ALL blocks of the global system), accessed via \a actives_vec[\a r_vec(i)]
     * @param[in] eliminatedDofs a vector of values for the dofs (corresponding to the columns) that are eliminated from the system
     *            (for ALL blocks of the global system), accessed via \a eliminatedDofs[\a r_vec(i)]
     * @param[in] r_vec a vector of row block indices to which the local matrix is pushed
     * @param[in] c_vec a vector of column block indices to which the local matrix is pushed
     */
    void pushToMatrix(const gsMatrix<T> & localMat,
                      const std::vector<gsMatrix<index_t> >& actives_vec,
                      const std::vector<gsMatrix<T> > & eliminatedDofs,
                      const gsVector<index_t> & r_vec,
                      const gsVector<index_t> & c_vec)
    {
        int rstrLocal = 0;
        int cstrLocal = 0;

        for (index_t r_ind = 0; r_ind != r_vec.size(); ++r_ind) // for row-blocks
        {
            size_t r = r_vec(r_ind);
            const gsDofMapper & rowMap    = m_mappers[m_row.at(r)];
            const index_t numActive_i = actives_vec[r].rows();

            for (index_t c_ind = 0; c_ind != c_vec.size(); ++c_ind) // for col-blocks
            {
                size_t c = c_vec(c_ind);
                const gsDofMapper & colMap    = m_mappers[m_col.at(c)];
                const index_t numActive_j = actives_vec[c].rows();
                const gsMatrix<T> & eliminatedDofs_j = eliminatedDofs[c];

                for (index_t i = 0; i != numActive_i; ++i) // N_i
                {
                    const int ii =  m_rstr.at(r) + actives_vec[r].at(i); // row index global matrix
                    const int iiLocal = rstrLocal + i;                   // row index local matrix

                    if ( rowMap.is_free_index(actives_vec[r].at(i)) )
                    {

                        for (index_t j = 0; j != numActive_j; ++j) // N_j
                        {
                            const int jj =  m_cstr.at(c) + actives_vec[c].at(j); // column index global matrix
                            const int jjLocal = cstrLocal + j;                   // column index local matrix

                            if ( colMap.is_free_index(actives_vec[c].at(j)) )
                            {
                                // If matrix is symmetric, we store only lower
                                // triangular part
                                if ( (!symm) || jj <= ii )
                                    m_matrix.coeffRef(ii, jj) += localMat(iiLocal, jjLocal);
                            }
                            else // Fixed DoF
                            {
                                m_rhs.row(ii).noalias() -= localMat(iiLocal, jjLocal) * eliminatedDofs_j.row( colMap.global_to_bindex(actives_vec[c].at(j)));
                            }
                        }
                    }
                }
                cstrLocal += numActive_j;
            }
            cstrLocal = 0;
            rstrLocal += numActive_i;
        }

    }



public: /* Add local contributions to system right-hand side */

    /**
     * @brief pushToRhs pushes the local rhs for an element to the global system
     * \note checks are done if an index is eliminated or not
     * @param[in] localRhs the local right hand side matrix/vector
     * @param[in] actives the corresponding mapped index of basis functions without shifts
     * @param[in] r the row block associated to
     */
    void pushToRhs(const gsMatrix<T> & localRhs,
                   const gsMatrix<index_t> & actives,
                   const size_t r = 0)
    {
        const gsDofMapper & mapper = m_mappers[m_row.at(r)];
        const index_t    numActive = actives.rows();

        for (index_t i = 0; i != numActive; ++i)
        {
            const index_t ii =  m_rstr.at(r) + actives.at(i);
            if ( mapper.is_free_index(actives.at(i)) )
            {
                m_rhs.row(ii) += localRhs.row(i);
            }
        }
    }

    /**
     * @brief pushToRhsAllFree pushes the local rhs for an element to the global system
     * \note no checks are done if an index is eliminated or not
     * @param[in] localRhs the local right hand side matrix/vector
     * @param[in] actives the corresponding mapped index of basis functions without shifts
     * @param[in] r the row block associated to
     */
    void pushToRhsAllFree(const gsMatrix<T> & localRhs,
                          const gsMatrix<index_t> & actives,
                          const size_t r = 0)
    {
        const index_t    numActive = actives.rows();

        for (index_t i = 0; i != numActive; ++i)
        {
            const index_t ii =  m_rstr.at(r) + actives.at(i);
            m_rhs.row(ii) += localRhs.row(i);
        }
    }


    /**
     * @brief pushToRhs pushes one local rhs consisting of several blocks corresponding to blocks of the global system
     * \note Usefull for rhs depending on a vector valued function
     * @param[in] localRhs local system matrix
     * @param[in] actives_vec a vector of mapped index sets (for ALL row blocks of the global system), accessed via \a actives_vec[\a r_vec(i)]
     * @param[in] r_vec a vector of row block indices to which the local matrix is pushed
     */

    void pushToRhs(const gsMatrix<T> & localRhs,
                   const std::vector<gsMatrix<index_t> >& actives_vec,
                   const gsVector<index_t> & r_vec)
    {
        int rstrLocal = 0;

        for (index_t r_ind = 0; r_ind != r_vec.size(); ++r_ind) // for row-blocks
        {
            index_t r = r_vec(r_ind);
            const gsDofMapper & rowMap    = m_mappers[m_row.at(r)];
            const index_t numActive_i = actives_vec[r].rows();


            for (index_t i = 0; i != numActive_i; ++i) // N_i
            {
                const int ii =  m_rstr.at(r) + actives_vec[r].at(i); // row index global matrix
                const int iiLocal = rstrLocal + i;                   // row index local matrix

                if ( rowMap.is_free_index(actives_vec[r].at(i)) )
                {
                    m_rhs.row(ii) += localRhs.row(iiLocal);
                }
            }
            rstrLocal += numActive_i;
        }

    }

public: /* Add local contributions to system matrix and right-hand side */

    /**
     * @brief push pushes the local system matrix and rhs for an element to the global system,
     * \note
     * 1. the same index set is assumed for row and column block
     * 2. eliminated dofs are incorporated in the right way
     * @param[in] localMat the local system matrix
     * @param[in] localRhs the local rhs matrix/vector
     * @param[in] actives the mapped index of basis functions, without shifts!
     * @param[in] eliminatedDofs the values for the dofs, which are removed from the system
     * @param[in] r the row block
     * @param[in] c the column block
     */
    void push(const gsMatrix<T> & localMat,
              const gsMatrix<T> & localRhs,
              const gsMatrix<index_t> & actives,
              const gsMatrix<T> & eliminatedDofs,
              const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive = actives.rows();
        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];

        GISMO_ASSERT( &rowMap == &m_mappers[m_col.at(c)], "Error");
        GISMO_ASSERT( m_matrix.cols() == m_rhs.rows(), "gsSparseSystem is not allocated");
        //Assert eliminatedDofs.rows() == rowMap.boundarySize()

        for (index_t i = 0; i != numActive; ++i)
        {
            const int ii =  m_rstr.at(r) + actives(i);
            if ( rowMap.is_free_index(actives.at(i)) )
            {
                m_rhs.row(ii) += localRhs.row(i);

                for (index_t j = 0; j < numActive; ++j)
                {
                    const int jj =  m_cstr.at(c) + actives(j);
                    if ( rowMap.is_free_index(actives.at(j)) )
                    {
                        // If matrix is symmetric, we store only lower
                        // triangular part
                        if ( (!symm) || jj <= ii )
                            m_matrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                    else // if ( mapper.is_boundary_index(jj) ) // Fixed DoF?
                    {
                        m_rhs.row(ii).noalias() -= localMat(i, j) *
                                eliminatedDofs.row( rowMap.global_to_bindex(actives.at(j)) );
                    }
                }
            }
        }
    }

    /**
     * @brief push pushes the local system matrix and rhs for an element to the global system,
     * \note
     * 1. different index sets can be used for row and column blocks
     * 2. eliminated dofs are incorporated in the right way
     * @param[in] localMat the local system matrix
     * @param[in] localRhs the local rhs matrix/vector
     * @param[in] actives_i the mapped index of row basis functions, without shifts!
     * @param[in] actives_j the mapped index of column basis functions, without shifts!
     * @param[in] eliminatedDofs_j the values for the dofs (corresponding to the columns), which are
     *                         removed from the system
     * @param[in] r the row block
     * @param[in] c the column block
     */
    void push(const gsMatrix<T> & localMat,
              const gsMatrix<T> & localRhs,
              const gsMatrix<index_t> & actives_i,
              const gsMatrix<index_t> & actives_j,
              const gsMatrix<T> & eliminatedDofs_j,
              const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive_i = actives_i.rows();
        const index_t numActive_j = actives_j.rows();
        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];
        const gsDofMapper & colMap = m_mappers[m_col.at(c)];

        GISMO_ASSERT( m_matrix.cols() == m_rhs.rows(), "gsSparseSystem is not allocated");
        //Assert eliminatedDofs.rows() == rowMap.boundarySize()

        for (index_t i = 0; i != numActive_i; ++i)
        {
            const int ii =  m_rstr.at(r) + actives_i.at(i);
            if ( rowMap.is_free_index(actives_i.at(i)) )
            {
                m_rhs.row(ii) += localRhs.row(i);

                for (index_t j = 0; j < numActive_j; ++j)
                {
                    const int jj =  m_cstr.at(c) + actives_j.at(j);
                    if ( colMap.is_free_index(actives_j.at(j)) )
                    {
                        // If matrix is symmetric, we store only lower
                        // triangular part
                        if ( (!symm) || jj <= ii )
                            m_matrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                    else // if ( mapper.is_boundary_index(jj) ) // Fixed DoF?
                    {
                        m_rhs.row(ii).noalias() -= localMat(i, j) *
                                eliminatedDofs_j.row( colMap.global_to_bindex(actives_j.at(j)) );
                    }
                }
            }
        }
    }


    /**
     * @brief pushAllFree pushes the local system matrix and rhs for an element to the global system,
     * \note
     * 1. the same index set is assumed for row and column block
     * 2. no checks are done if an index is eliminated or not
     *
     * Use this functions if you already know that your block has no eliminated dofs
     * @param[in] localMat the local system matrix
     * @param[in] localRhs the local rhs matrix/vector
     * @param[in] actives the mapped index of basis functions, without shifts!
     * @param[in] r the row block
     * @param[in] c the column block
     */
    void pushAllFree(const gsMatrix<T>  & localMat,
                     const gsMatrix<T>  & localRhs,
                     const gsMatrix<index_t> & actives,
                     const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive = actives.rows();
        //GISMO_ASSERT( &m_mappers[m_row.at(r)] == &m_mappers[m_col.at(c)], "Error");

        for (index_t j=0; j!=numActive; ++j)
        {
            const unsigned jj = m_cstr.at(c) + actives(j);
            m_rhs.row(jj) += localRhs.row(j);
            for (index_t i=0; i!=numActive; ++i)
            {
                const unsigned ii = m_rstr.at(r) + actives(i);
                // If matrix is symmetric, we store only lower
                // triangular part
                if ( (!symm) || jj <= ii )
                    m_matrix( ii, jj ) += localMat(i,j);
            }
        }
    }


    /**
     * @brief pushToMatrix pushes one local matrix and rhs consisting of several blocks corresponding to blocks of the global system
     * \note
     * 1. Usefull for bilinear forms depending on vector valued functions
     * 2. different index sets are used for row and column blocks
     * 3. eliminated dofs are incorporated in the right way
     * 4. assume identical row and column mappers for the global system, therefore only one vector of mapped index sets is given
     * @param[in] localMat local system matrix
     * @param[in] localRhs local rhs vector
     * @param[in] actives_vec a vector of mapped index sets (for ALL blocks of the global system), accessed via \a actives_vec[\a r_vec(i)]
     * @param[in] eliminatedDofs a vector of values for the dofs (corresponding to the columns) that are eliminated from the system
     *            (for ALL blocks of the global system), accessed via \a eliminatedDofs[\a r_vec(i)]
     * @param[in] r_vec a vector of row block indices to which the local matrix is pushed
     * @param[in] c_vec a vector of column block indices to which the local matrix is pushed
     */

    void push(const gsMatrix<T> & localMat,
              const gsMatrix<T> & localRhs,
              const std::vector<gsMatrix<index_t> >& actives_vec,
              const std::vector<gsMatrix<T> > & eliminatedDofs,
              const gsVector<index_t> & r_vec,
              const gsVector<index_t> & c_vec)
    {
        int rstrLocal = 0;
        int cstrLocal = 0;

        for (index_t r_ind = 0; r_ind != r_vec.size(); ++r_ind) // for row-blocks
        {
            index_t r = r_vec(r_ind);
            const gsDofMapper & rowMap    = m_mappers[m_row.at(r)];
            const index_t numActive_i = actives_vec[r].rows();

            for (index_t c_ind = 0; c_ind != c_vec.size(); ++c_ind) // for col-blocks
            {
                index_t c = c_vec(c_ind);
                const gsDofMapper & colMap    = m_mappers[m_col.at(c)];
                const index_t numActive_j = actives_vec[c].rows();
                const gsMatrix<T> & eliminatedDofs_j = eliminatedDofs[c];

                for (index_t i = 0; i != numActive_i; ++i) // N_i
                {
                    const int ii =  m_rstr.at(r) + actives_vec[r].at(i); // row index global matrix
                    const int iiLocal = rstrLocal + i;                   // row index local matrix

                    if ( rowMap.is_free_index(actives_vec[r].at(i)) )
                    {
                        // rhs should not be pushed for each col-block (but only once)
                        if(c_ind == 0)
                            m_rhs.row(ii) += localRhs.row(iiLocal);

                        for (index_t j = 0; j != numActive_j; ++j) // N_j
                        {
                            const int jj =  m_cstr.at(c) + actives_vec[c].at(j); // column index global matrix
                            const int jjLocal = cstrLocal + j;                   // column index local matrix

                            if ( colMap.is_free_index(actives_vec[c].at(j)) )
                            {
                                // If matrix is symmetric, we store only lower
                                // triangular part
                                if ( (!symm) || jj <= ii )
                                    m_matrix.coeffRef(ii, jj) += localMat(iiLocal, jjLocal);
                            }
                            else // Fixed DoF
                            {
                                m_rhs.row(ii).noalias() -= localMat(iiLocal, jjLocal) * eliminatedDofs_j.row( colMap.global_to_bindex(actives_vec[c].at(j)));
                            }
                        }
                    }
                }
                cstrLocal += numActive_j;
            }
            cstrLocal = 0;
            rstrLocal += numActive_i;
        }

    }


    // AM: incomplete, not working yet
    /**
     * @brief push pushes one local system matrix and one local rhs for an element to
     *        the global system for several blocks
     * \note
     * 1. the same index set is assumed for row and column blocks
     * 2. eliminated dofs are incorporated in the right way
     * 3. several blocks can be pused with one function call, but the same local
     *    matrix is assumed
     * 4. if the size of the vector is smaller, than the number of blocks, then only the
     *    the first \a actives.size() blocks are filled.
     * @param[in] localMat local system matrices
     * @param[in] localRhs local rhs vector/matrices
     * @param[in] actives a vector of index sets (one per block)
     * @param[in] eliminatedDofs a vector of values for the dofs (corresponding to the columns), which are
     *           removed from the system, (one per block)
     */
    void push(const gsMatrix<T> & localMat,
              const gsMatrix<T> & localRhs,
              const std::vector<gsMatrix<index_t> >& actives,
              const std::vector<gsMatrix<T> > & eliminatedDofs)
    {
        GISMO_ASSERT( m_matrix.cols() == m_rhs.rows(), "gsSparseSystem is not allocated");

        for (size_t r = 0; r != actives.size(); ++r) // for all row-blocks
        {
            const gsDofMapper & rowMap    = m_mappers[m_row.at(r)];
            const index_t numRowActive = actives[r].rows();

            for (size_t c = 0; c != actives.size(); ++c) // for all col-blocks
            {
                const gsMatrix<T> & fixedDofs = eliminatedDofs[m_col.at(c)];

                for (index_t i = 0; i != numRowActive; ++i)
                {
                    const int ii =  m_rstr.at(r) + actives[r].at(i);
                    if ( rowMap.is_free_index(actives[r].at(i)) )
                    {
                        m_rhs.row(ii) += localRhs.row(i + r * numRowActive); //  + c *
                        const index_t numColActive = actives[c].rows();

                        for (index_t j = 0; j < numColActive; ++j)
                        {
                            const int jj =  m_cstr.at(c) + actives[c].at(j);
                            if ( rowMap.is_free_index(actives[c].at(j)) )
                            {
                                // If matrix is symmetric, we store only lower
                                // triangular part
                                if ( (!symm) || jj <= ii )
                                    m_matrix.coeffRef(ii, jj) += localMat(i + r * numRowActive,
                                                                          j + c * numRowActive); //  + c * ..
                            }
                            else // if ( mapper.is_boundary_index(jj) ) // Fixed DoF?
                            {
                                m_rhs.at(ii) -= localMat(i + r * numRowActive, j + c * numRowActive) *  //  + c *..
                                    fixedDofs.coeff( rowMap.global_to_bindex(actives[c].at(j)), 0 );
                            }
                        }
                    }
                }
            }
        }
    }

    /**
     * @brief push pushes several local system matrices and local rhs for an element to
     *        the global system for several blocks
     * \note
     * 1. the same index set is assumed for row and column blocks
     * 2. eliminated dofs are incorporated in the right way
     * 3. several blocks can be pused with one function call
     *
     * <hr>\b Parameters
     * \n[in]\b localMat a vector local system matrices
     * \n[in]\b localRhs a vector local rhs vector/matrices
     * \n[in]\b actives a vector of index sets (one per block)
     * \n[in]\b fixedDofs
     * \n[in]\b r the row blocks where the matrices should be pushed
     * \n[in]\b c the colum blocks where the matrices shouldbe pushed
     */
    void push(const std::vector<gsMatrix<T> > &,
              const std::vector<gsMatrix<T> > &,
              const std::vector<gsMatrix<index_t> > &,
              const std::vector<gsMatrix<T> > &,
              const gsVector<index_t> &, const gsVector<index_t> &)
    {
        GISMO_NO_IMPLEMENTATION
    }


    GISMO_DEPRECATED
    /**
     * @brief pushToMatrix pushes the local matrix for an element to the global system,
     * \note
     *  1. the same index set is assumed for row and column block
     *  2. dofs wich are eliminated are NOT moved to the rhs, they are just ignored in
     *        the assembly (only works for homogeneous BC)
     *  3. no assembling is done for the rhs
     * @param[in] localMat the local matrix
     * @param[in] actives the mapped index of basis functions, without shifts!
     * @param[in] r the row block
     * @param[in] c the column block
     */
    void pushToMatrix(const gsMatrix<T>  & localMat,
                      const gsMatrix<index_t> & actives,
                      const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive = actives.rows();
        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];
        GISMO_ASSERT( &rowMap == &m_mappers[m_col.at(c)], "Error");

        for (index_t i = 0; i != numActive; ++i)
        {
            const int ii = m_rstr.at(r) + actives.at(i); // N_i

            if ( rowMap.is_free_index(actives.at(i)) )
            {
                for (index_t j = 0; j != numActive; ++j)
                {
                    const int jj = m_cstr.at(c) + actives.at(j); // N_j

                    if ( rowMap.is_free_index( actives.at(j)) )
                        // If matrix is symmetric, we store only lower
                        // triangular part
                        if ( (!symm) || jj <= ii )
                            m_matrix.coeffRef(ii, jj) += localMat(i, j);
                }
            }
        }
    }

    GISMO_DEPRECATED
    /**
     * @brief pushToMatrix pushes the local matrix for an element to the global system,
     * \note
     * 1. the different index sets are used for row and column block
     * 2. dofs wich are eliminated are NOT moved to the rhs, they are just ignored in
     *    the assembly (only works for homogeneous BC)
     * 3. no assembling is done for the rhs
     * @param[in] localMat the local matrix
     * @param[in] actives_i the mapped index of row - basis functions, without shifts!
     * @param[in] actives_j the mapped index of column - basis functions, without shifts!
     * @param[in] r the row block
     * @param[in] c the column block
     */
    void pushToMatrix(const gsMatrix<T>  & localMat,
                      const gsMatrix<index_t> & actives_i,
                      const gsMatrix<index_t> & actives_j,
                      const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive_i = actives_i.rows();
        const index_t numActive_j = actives_j.rows();

        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];
        const gsDofMapper & colMap = m_mappers[m_col.at(c)];

        for (index_t i = 0; i != numActive_i; ++i)
        {
            const int ii = m_rstr.at(r) + actives_i.at(i); // N_i

            if ( rowMap.is_free_index(actives_i.at(i)) )
            {
                for (index_t j = 0; j != numActive_j; ++j)
                {
                    const int jj = m_cstr.at(c) + actives_j.at(j); // N_j

                    if ( colMap.is_free_index(actives_j.at(j)) )
                        // If matrix is symmetric, we store only lower
                        // triangular part
                        if ( (!symm) || jj <= ii )
                            m_matrix.coeffRef(ii, jj) += localMat(i, j);
                }
            }
        }
    }

    GISMO_DEPRECATED
    /**
     * @brief push pushes the local system matrix and rhs for an element to the global system,
     * \note
     * 1. the same index set is assumed for row and column block
     * 2. dofs wich are eliminated are NOT moved to the rhs, they are just ignored in
     *    the assembly (only works for homogeneous BC)
     * @param[in] localMat the local system matrix
     * @param[in] localRhs the local rhs matrix/vector
     * @param[in] actives the mapped index of basis functions, without shifts!
     * @param[in] r the row block
     * @param[in] c the column block
     */
    void push(const gsMatrix<T>  & localMat,
              const gsMatrix<T>  & localRhs,
              const gsMatrix<index_t> & actives,
              const size_t r = 0, const size_t c = 0)
    {
        GISMO_ASSERT( m_matrix.cols() == m_rhs.rows(), "gsSparseSystem is not allocated");

        const index_t numActive = actives.rows();
        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];
        GISMO_ASSERT( &rowMap == &m_mappers[m_col.at(c)], "Error");

        for (index_t i = 0; i != numActive; ++i)
        {
            const int ii =  m_rstr.at(r) + actives(i);
            if ( rowMap.is_free_index(actives(i)) )
            {
                m_rhs.row(ii) += localRhs.row(i);

                for (index_t j = 0; j != numActive; ++j)
                {
                    const int jj =  m_cstr.at(c) + actives(j);
                    if ( rowMap.is_free_index(actives(j)) )
                    {
                        // If matrix is symmetric, we store only lower
                        // triangular part
                        if ( (!symm) || jj <= ii )
                            m_matrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                }
            }
        }
    }


public:
    void pushSparse(const gsSparseMatrix<T> & localMat,
                    const gsMatrix<T> & localRhs,
                    const gsMatrix<index_t> & actives_i,
                    const gsMatrix<index_t> & actives_j,
                    const gsMatrix<T> & eliminatedDofs_j,
                    const size_t r = 0, const size_t c = 0)
    {
        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];
        const gsDofMapper & colMap = m_mappers[m_col.at(c)];

        GISMO_ASSERT( m_matrix.cols() == m_rhs.rows(), "gsSparseSystem is not allocated");
        //Assert eliminatedDofs.rows() == rowMap.boundarySize()

        for (index_t i = 0; i<localMat.outerSize(); ++i)
            for (typename gsSparseMatrix<T>::iterator it(localMat,i); it; ++it)
            {
                const int ii =  m_rstr.at(r) + actives_i.at(it.row());
                if ( rowMap.is_free_index(actives_i.at(it.row())) )
                {
                    const int jj =  m_cstr.at(c) + actives_j.at(it.col());
                    if ( colMap.is_free_index(actives_j.at(it.col())) )
                    {
                        // If matrix is symmetric, we store only lower
                        // triangular part
                        if ( (!symm) || jj <= ii )
                            m_matrix.coeffRef(ii, jj) += it.value();
                    }
                    else // if ( mapper.is_boundary_index(jj) ) // Fixed DoF?
                    {
                        m_rhs.row(ii).noalias() -= it.value() *
                                eliminatedDofs_j.row( colMap.global_to_bindex(jj) );
                    }
                }
            }

        for(index_t i=0; i<actives_i.rows();++i)
        {
            const int ii =  m_rstr.at(r) + actives_i.at(i);
            if ( rowMap.is_free_index(actives_i.at(i)) )
                m_rhs.row(ii) += localRhs.row(i);
        }
    }

};  // class gsSparseSystem


template<class T>
std::ostream &operator<<(std::ostream &os, gsSparseSystem<T> & ss)
{
    os << ss.blockView();
    os << "Right-hand side: "
       << ss.rhs().rows() <<"x"<< ss.rhs().cols() <<"\n";
    return os;
}


} // namespace gismo
