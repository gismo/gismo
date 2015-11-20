/** @file gsSparseSystem.h

    @brief Class representing a sparse linear system (with dense
    right-hand side(s)), indexed by sets of degrees of freedom

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

//#include <gsCore/gsRefVector.h>

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

    gsSparseMatrix<T> m_matrix;

    gsMatrix<T> m_rhs;


    // -- Structure
    DofMappers m_mappers;

    gsVector<size_t> m_row;

    gsVector<size_t> m_col;    

    gsVector<index_t> m_rstr;

    gsVector<index_t> m_cstr;

    gsVector<index_t> m_cvar;
    //

public:

    gsSparseSystem() 
    { }

    gsSparseSystem(gsDofMapper & mapper)
    : m_mappers(1), 
      m_row    (1),
      m_col    (1),
      m_rstr   (1),
      m_cstr   (1),
      m_cvar   (1)
    {
        m_row [0] =  m_col [0] = 
        m_rstr[0] =  m_cstr[0] = 
        m_cvar[0] = 0;

        m_mappers.front().swap(mapper);
        
        m_matrix.resize( m_mappers.front().freeSize() , 
                         m_mappers.front().freeSize() );        
    }

    gsSparseSystem(DofMappers & mappers,
                   const gsVector<unsigned> & dims)
    : m_row (gsVector<size_t>::LinSpaced(dims.size(),0,dims.size()-1)),
      m_col (gsVector<size_t>::LinSpaced(dims.size(),0,dims.size()-1)),
      m_rstr(dims.sum()),
      m_cstr(dims.sum())
    {
        m_mappers.swap(mappers);

        const index_t d = dims.size();
        const index_t s = dims.sum();

        if ( static_cast<index_t>(m_mappers.size()) == 2*d )
        {
            m_col.array() += d;
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
            GISMO_ASSERT(m_mappers.size() == d, "Cannot deduce block structure.");
        }

        m_rstr[0] = m_cstr[0] = 0;
        for (index_t r = 1; r < d; ++r) // for all row-blocks
            m_rstr[r] = m_rstr[r-1] + m_mappers[m_row[r-1]].freeSize();
        for (index_t c = 1; c < d; ++c) // for all col-blocks
            m_cstr[c] = m_cstr[c-1] + m_mappers[m_col[c-1]].freeSize();

        m_matrix.resize( m_rstr.at(d-1) + m_mappers[m_row[d-1]].freeSize() , 
                         m_cstr.at(d-1) + m_mappers[m_col[d-1]].freeSize() );

    }

    gsSparseSystem(DofMappers & mappers, 
                   const index_t rows, 
                   const index_t cols)
    : m_row (gsVector<size_t>::LinSpaced(rows,0,rows-1)),
      m_col (gsVector<size_t>::LinSpaced(cols,0,cols-1)),
      m_rstr(rows),
      m_cstr(cols)
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

        m_rstr[0] = m_cstr[0] = 0;
        for (index_t r = 1; r < rows; ++r) // for all row-blocks
            m_rstr[r] = m_rstr[r-1] + m_mappers[m_row[r-1]].freeSize();
        for (index_t c = 1; c < cols; ++c) // for all col-blocks
            m_cstr[c] = m_cstr[c-1] + m_mappers[m_col[c-1]].freeSize();

        m_matrix.resize( m_rstr.at(rows-1) + m_mappers[m_row[rows-1]].freeSize() , 
                         m_cstr.at(cols-1) + m_mappers[m_col[cols-1]].freeSize() );
    }

    gsSparseSystem(DofMappers & mappers, 
                   const gsVector<index_t> & rowInd, 
                   const gsVector<index_t> & colInd, 
                   const gsVector<index_t> & colvar)
    : m_row (rowInd),
      m_col (colInd),
      m_rstr(rowInd.size()),
      m_cstr(colInd.size()),
      m_cvar(colvar)
    {
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
    }
    
    /// Non-zeros per column for the sparse matrix, number of columns
    /// for the right-hand side
    void reserve(const index_t nz, const index_t numRhs)
    {
        GISMO_ASSERT( 0 != m_mappers.size(), "Sparse system was not initialized");
        m_matrix.reservePerColumn(nz); 
        m_rhs.setZero(m_matrix.cols(), numRhs);
    }

    void setZero()
    {
        m_matrix.setZero();
        m_rhs   .setZero();
    }

    index_t cols()
    {
        return m_rhs.rows();
    }

public: /* Accessors */

    const gsSparseMatrix<T> & matrix() const 
    { return m_matrix; }

    gsSparseMatrix<T> & matrix() 
    { return m_matrix; }

    const gsMatrix<T> & rhs() const 
    { return m_rhs; }

    gsMatrix<T> & rhs() 
    { return m_rhs; }

    matBlockView blockView()
    {
        gsVector<index_t> rowSizes(m_row.size()), colSizes(m_row.size());

        for (index_t r = 0; r != rowSizes.size(); ++r) // for all row-blocks
            rowSizes[r] = m_mappers[r].freeSize();

        for (size_t c = 0; c != colSizes.size(); ++c) // for all col-blocks
            colSizes[c] = m_mappers[c].freeSize();

        return m_matrix.blockView(rowSizes,colSizes); 
    }

    const gsDofMapper & rowMapper(const index_t r) const
    { return m_mappers[m_row[r]]; }

    gsDofMapper & rowMapper(const index_t r)
    { return m_mappers[m_row[r]]; }

    const gsDofMapper & colMapper(const index_t c) const
    {
        return m_mappers[m_col[c]]; 
    }
    
    bool initialized() const
    {
        return 0 != m_col.size();
    }

    bool symmetry() const
    {
        return symm;
    }

    // Returns the mapper index for column block \a c
    gsDofMapper & colMapper(const index_t c)
    { return m_mappers[m_col[c]]; }

    // Returns the basis index for column block \a c
    index_t colBasis(const index_t c) // better name ?
    { return m_cvar[c]; }

    const DofMappers & dofMappers() const
    { return m_mappers; }

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

    void mapRowIndices(const gsMatrix<unsigned> & actives,
                       const index_t patchIndex, 
                       gsMatrix<unsigned> & result,
                       const size_t r = 0) const 
    {
        m_mappers[m_row.at(r)].localToGlobal(actives, patchIndex, result);
    }
    
    void mapColIndices(gsMatrix<unsigned> & actives,
                       const index_t patchIndex, 
                       gsMatrix<unsigned> & result,
                       const size_t c = 0) const 
    {
        m_mappers[m_col.at(c)].localToGlobal(actives, patchIndex, result);
    }
 
public: /* Add local contributions to system matrix */

    // Same index sets for row and column
    void pushToMatrix(const gsMatrix<T>  & localMat, 
                      const gsMatrix<unsigned> & actives,
                      const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive = actives.rows();
        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];
        GISMO_ASSERT( &rowMap == &m_mappers[m_col.at(c)], "Error");

        for (index_t i = 0; i != numActive; ++i)
        {
            const int ii = m_rstr.at(r) + actives.at(i); // N_i
            
            if ( rowMap.is_free_index(ii) )
            {
                for (index_t j = 0; j != numActive; ++j)
                {
                    const int jj = m_cstr.at(c) + actives.at(j); // N_j

                    if ( rowMap.is_free_index(jj) )
                        // If matrix is symmetric, we store only lower
                        // triangular part
                        if ( (!symm) || jj <= ii ) 
                            m_matrix.coeffRef(ii, jj) += localMat(i, j);
                }
            }
        }
    }

    // Different index sets for row and column
    void pushToMatrix(const gsMatrix<T>  & localMat, 
                     const gsMatrix<unsigned> & actives_i,
                     const gsMatrix<unsigned> & actives_j,
                     const size_t r = 0, const size_t c = 0)
    {
        GISMO_NO_IMPLEMENTATION
    }

public: /* Add local contributions to system right-hand side */

    void pushToRhs(const gsMatrix<T> & localRhs, 
                  const gsMatrix<unsigned> & actives, 
                  const size_t c = 0)
    {
        const gsDofMapper & mapper = m_mappers[m_col.at(c)];
        const index_t    numActive = actives.rows();

        for (index_t i = 0; i != numActive; ++i)
        {
            const int ii =  m_cstr.at(c) + actives.at(i);
            if ( mapper.is_free_index(ii) )
            {
                m_rhs.row(ii) += localRhs.row(i);
            }
        }

    }

public: /* Add local contributions to system matrix and right-hand side */

    // Same index set for row and column block
    void push(const gsMatrix<T>  & localMat, 
              const gsMatrix<T>  & localRhs,
              const gsMatrix<unsigned> & actives,
              const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive = actives.rows();
        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];
        GISMO_ASSERT( &rowMap == &m_mappers[m_col.at(c)], "Error");

        for (index_t i = 0; i != numActive; ++i)
        {
            const int ii =  m_rstr.at(r) + actives(i);
            if ( rowMap.is_free_index(ii) )
            {
                m_rhs.row(ii) += localRhs.row(i);
                
                for (index_t j = 0; j != numActive; ++j)
                {
                    const int jj =  m_cstr.at(c) + actives(j);
                    if ( rowMap.is_free_index(jj) )
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

    // Same index set for row and column block, no checks for
    // free/boundary dofs
    void pushAllFree(const gsMatrix<T>  & localMat, 
                     const gsMatrix<T>  & localRhs,
                     const gsMatrix<unsigned> & actives,
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
    
    // Same index set for row and column block, with fixed DoF values
    void push(const gsMatrix<T> & localMat, 
              const gsMatrix<T> & localRhs,
              const gsMatrix<unsigned> & actives,
              const gsMatrix<T> & eliminatedDofs,
              const size_t r = 0, const size_t c = 0)
    {
        const index_t numActive = actives.rows();
        const gsDofMapper & rowMap = m_mappers[m_row.at(r)];
        
        GISMO_ASSERT( &rowMap == &m_mappers[m_col.at(c)], "Error");
        //Assert eliminatedDofs.rows() == rowMap.boundarySize()

        for (index_t i = 0; i != numActive; ++i)
        {
            const int ii =  m_rstr.at(r) + actives(i);
            if ( rowMap.is_free_index(ii) )
            {
                m_rhs.row(ii) += localRhs.row(i);
                
                for (index_t j = 0; j < numActive; ++j)
                {
                    const int jj =  m_cstr.at(c) + actives(j);
                    if ( rowMap.is_free_index(jj) )
                    {
                        // If matrix is symmetric, we store only lower
                        // triangular part
                        if ( (!symm) || jj <= ii ) 
                            m_matrix.coeffRef(ii, jj) += localMat(i, j);
                    }
                    else // if ( mapper.is_boundary_index(jj) ) // Fixed DoF?
                    {
                        m_rhs.row(ii).noalias() -= localMat(i, j) * 
                            eliminatedDofs.row( rowMap.global_to_bindex(jj) );
                    }
                }
            }
        }
    }

    // Different index sets for row and column blocks
    void push(const gsMatrix<T> & localMat, 
              const gsMatrix<T> & localRhs,
              const gsMatrix<unsigned> & actives_i,
              const gsMatrix<unsigned> & actives_j,
              const gsMatrix<T> & eliminatedDofs_j,
              const size_t r = 0, const size_t c = 0)
    {
        GISMO_NO_IMPLEMENTATION
    }

    // Local matrix and rhs with the same structure as the global system, 
    // same index sets -- one per block
    // AM: incomplete, not working yet
    void push(const gsMatrix<T> & localMat, 
              const gsMatrix<T> & localRhs,
              const std::vector<gsMatrix<unsigned> >& actives, 
              const std::vector<gsMatrix<T> > & eliminatedDofs)
    {
        for (size_t r = 0; r != actives.size(); ++r) // for all row-blocks
        {
            const gsDofMapper & rowMap    = m_mappers[m_row.at(r)];
            const gsMatrix<T> & fixedDofs = eliminatedDofs[m_row.at(r)];
            const index_t numRowActive = actives[r].rows();

            for (size_t c = 0; c != actives.size(); ++c) // for all col-blocks
            {
                for (index_t i = 0; i != numRowActive; ++i)
                {
                    const int ii =  m_rstr.at(r) + actives[r].at(i);
                    if ( rowMap.is_free_index(ii) )
                    {
                        m_rhs.row(ii) += localRhs.row(i); //  + c * 
                        const index_t numColActive = actives[c].rows();
  
                        for (index_t j = 0; j < numColActive; ++j)
                        {
                            const int jj =  m_cstr.at(c) + actives[c].at(j);
                            if ( rowMap.is_free_index(jj) )
                            {
                                // If matrix is symmetric, we store only lower
                                // triangular part
                                if ( (!symm) || jj <= ii ) 
                                    m_matrix.coeffRef(ii, jj) += localMat(i, j); //  + c * ..
                            }
                            else // if ( mapper.is_boundary_index(jj) ) // Fixed DoF?
                            {
                                m_rhs.row(ii).noalias() -= localMat(i, j) *  //  + c *..
                                    fixedDofs.row( rowMap.global_to_bindex(jj) );
                            }
                        }
                    }
                }
            }   
        }
    }

    // Seveal local matrices/rhs, same row/col indexing
    void push(const std::vector<gsMatrix<T> > & localMat, 
              const std::vector<gsMatrix<T> > & localRhs,
              const std::vector<gsMatrix<unsigned> > & actives,
              const std::vector<gsMatrix<T> > & fixedDofs,
              const gsVector<size_t> & r, const gsVector<size_t> & c)
    {
        GISMO_NO_IMPLEMENTATION
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
