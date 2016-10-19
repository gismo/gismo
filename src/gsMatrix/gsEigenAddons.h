
/** @file gsEigenAddons.h

    @brief Extends functionality of the Eigen library

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

namespace Eigen { 

namespace internal {

// *** General case adjugate implementation

template<typename MatrixType, typename ResultType, int Size = MatrixType::RowsAtCompileTime>
struct compute_adjugate
{
    static inline void run(const MatrixType& matrix, ResultType& result)
    {
        eigen_assert( matrix.rows() == matrix.cols() && "Matrix is not square.");

        switch( matrix.rows() )
        {
        case 1:
            compute_adjugate<MatrixType, ResultType, 1>::run(matrix,result);
            break;
        case 2:
            compute_adjugate<MatrixType, ResultType, 2>::run(matrix,result);
            break;
        case 3:
            compute_adjugate<MatrixType, ResultType, 3>::run(matrix,result);
            break;
        case 4:
            compute_adjugate<MatrixType, ResultType, 4>::run(matrix,result);
            break;
        default:
            eigen_assert( false && "Not implemented.");
            break;
        };
    }

};

// ***  Size 1 implementation

template<typename MatrixType, typename ResultType>
struct compute_adjugate<MatrixType, ResultType, 1>
{
    static inline void run(const MatrixType& matrix, ResultType& result)
    {
        result.coeffRef(0,0) = typename ResultType::Scalar(1);
    }
};

// ***  Size 2 implementation

template<typename MatrixType, typename ResultType>
struct compute_adjugate<MatrixType, ResultType, 2>
{
    static inline void run(const MatrixType& matrix, ResultType& result)
    {
        result.coeffRef(0,0) =  matrix.coeff(1,1);
        result.coeffRef(1,0) = -matrix.coeff(1,0);
        result.coeffRef(0,1) = -matrix.coeff(0,1);
        result.coeffRef(1,1) =  matrix.coeff(0,0);
    }
};

// ***  Size 3 implementation

template<typename MatrixType, typename ResultType>
struct compute_adjugate<MatrixType, ResultType, 3>
{
    static inline void run(const MatrixType& matrix, ResultType& result)
    {
        result.coeffRef(0,0) =  cofactor_3x3<MatrixType,0,0>(matrix);
        result.coeffRef(0,1) =  cofactor_3x3<MatrixType,1,0>(matrix);
        result.coeffRef(0,2) =  cofactor_3x3<MatrixType,2,0>(matrix);
        result.coeffRef(1,0) =  cofactor_3x3<MatrixType,0,1>(matrix);
        result.coeffRef(1,1) =  cofactor_3x3<MatrixType,1,1>(matrix);
        result.coeffRef(1,2) =  cofactor_3x3<MatrixType,2,1>(matrix);
        result.coeffRef(2,0) =  cofactor_3x3<MatrixType,0,2>(matrix);
        result.coeffRef(2,1) =  cofactor_3x3<MatrixType,1,2>(matrix);
        result.coeffRef(2,2) =  cofactor_3x3<MatrixType,2,2>(matrix);
    }
};

// ***  Size 4 implementation

template<typename MatrixType, typename ResultType>
struct compute_adjugate<MatrixType, ResultType, 4>
{
    static inline void run(const MatrixType& matrix, ResultType& result)
    {
        result.coeffRef(0,0) =  cofactor_4x4<MatrixType,0,0>(matrix);
        result.coeffRef(1,0) = -cofactor_4x4<MatrixType,0,1>(matrix);
        result.coeffRef(2,0) =  cofactor_4x4<MatrixType,0,2>(matrix);
        result.coeffRef(3,0) = -cofactor_4x4<MatrixType,0,3>(matrix);
        result.coeffRef(0,2) =  cofactor_4x4<MatrixType,2,0>(matrix);
        result.coeffRef(1,2) = -cofactor_4x4<MatrixType,2,1>(matrix);
        result.coeffRef(2,2) =  cofactor_4x4<MatrixType,2,2>(matrix);
        result.coeffRef(3,2) = -cofactor_4x4<MatrixType,2,3>(matrix);
        result.coeffRef(0,1) = -cofactor_4x4<MatrixType,1,0>(matrix);
        result.coeffRef(1,1) =  cofactor_4x4<MatrixType,1,1>(matrix);
        result.coeffRef(2,1) = -cofactor_4x4<MatrixType,1,2>(matrix);
        result.coeffRef(3,1) =  cofactor_4x4<MatrixType,1,3>(matrix);
        result.coeffRef(0,3) = -cofactor_4x4<MatrixType,3,0>(matrix);
        result.coeffRef(1,3) =  cofactor_4x4<MatrixType,3,1>(matrix);
        result.coeffRef(2,3) = -cofactor_4x4<MatrixType,3,2>(matrix);
        result.coeffRef(3,3) =  cofactor_4x4<MatrixType,3,3>(matrix);
    }
};

// *** Adjugate expression

template<typename MatrixType>
struct adjugate_impl : public ReturnByValue<adjugate_impl<MatrixType> >
{
    typedef typename MatrixType::Index Index;
    typedef typename internal::eval<MatrixType>::type MatrixTypeNested;
    typedef typename remove_all<MatrixTypeNested>::type MatrixTypeNestedCleaned;
    MatrixTypeNested m_matrix;

    adjugate_impl(const MatrixType& matrix)
    : m_matrix(matrix)
    {}

    inline Index rows() const { return m_matrix.rows(); }
    inline Index cols() const { return m_matrix.cols(); }

    template<typename Dest> inline void evalTo(Dest& dst) const
    {
        static const int Size = EIGEN_PLAIN_ENUM_MIN(MatrixType::ColsAtCompileTime,Dest::ColsAtCompileTime);
        EIGEN_ONLY_USED_FOR_DEBUG(Size);
        eigen_assert(( (Size<=1) || (Size>4) || (extract_data(m_matrix)!=extract_data(dst)))
                     && "Aliasing problem detected in adjugate(), you need to do adjugate().eval() here.");

        compute_adjugate<MatrixTypeNestedCleaned, Dest>::run(m_matrix, dst);
    }
};

template<typename MatrixType>
struct traits<adjugate_impl<MatrixType> >
{
    typedef typename MatrixType::PlainObject ReturnType;
};

} // namespace internal


// ***  Implementation of MatrixBase::adjugate()

template<typename Derived>
inline const internal::adjugate_impl<Derived>
MatrixBase<Derived>::adjugate() const
{
    eigen_assert(rows() == cols());
    return internal::adjugate_impl<Derived>(derived());
}

template<typename Derived>
inline void MatrixBase<Derived>::adjugateInPlace()
{
    derived() = adjugate().eval();
}


/**
  * \class BlockDiag
  *
  * \brief Expression for block diagonal replication of a matrix or vector
  *
  * \param MatrixType the type of the object we are replicating
  *
  */
namespace internal {
template<typename MatrixType,int NumBlocks>
struct traits<BlockDiag<MatrixType,NumBlocks> >
 : traits<MatrixType>
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename traits<MatrixType>::StorageKind StorageKind;
  typedef typename traits<MatrixType>::XprKind XprKind;
  enum {
    Factor = (NumBlocks==Dynamic) ? Dynamic : NumBlocks*NumBlocks
  };
  typedef typename nested<MatrixType,Factor>::type MatrixTypeNested;
  typedef typename remove_reference<MatrixTypeNested>::type _MatrixTypeNested;
  enum {
    RowsAtCompileTime = NumBlocks==Dynamic || int(MatrixType::RowsAtCompileTime)==Dynamic
                      ? Dynamic
                      : NumBlocks * MatrixType::RowsAtCompileTime,
    ColsAtCompileTime = NumBlocks==Dynamic || int(MatrixType::ColsAtCompileTime)==Dynamic
                      ? Dynamic
                      : NumBlocks * MatrixType::ColsAtCompileTime,
   //FIXME we don't propagate the max sizes !!!
    MaxRowsAtCompileTime = RowsAtCompileTime,
    MaxColsAtCompileTime = ColsAtCompileTime,
    IsRowMajor = MaxRowsAtCompileTime==1 && MaxColsAtCompileTime!=1 ? 1
               : MaxColsAtCompileTime==1 && MaxRowsAtCompileTime!=1 ? 0
               : (MatrixType::Flags & RowMajorBit) ? 1 : 0,
    Flags = (_MatrixTypeNested::Flags & HereditaryBits & ~RowMajorBit) | (IsRowMajor ? RowMajorBit : 0),
    CoeffReadCost = _MatrixTypeNested::CoeffReadCost
  };
};
}

template<typename MatrixType,int NumBlocks> class BlockDiag
  : public internal::dense_xpr_base< BlockDiag<MatrixType,NumBlocks> >::type
{
    typedef typename internal::traits<BlockDiag>::MatrixTypeNested MatrixTypeNested;
    typedef typename internal::traits<BlockDiag>::_MatrixTypeNested _MatrixTypeNested;
  public:

    typedef typename internal::dense_xpr_base<BlockDiag>::type Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(BlockDiag)

    template<typename OriginalMatrixType>
    inline explicit BlockDiag(const OriginalMatrixType& a_matrix)
      : m_matrix(a_matrix), m_numBlocks(NumBlocks)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<typename internal::remove_const<MatrixType>::type,OriginalMatrixType>::value),
                          THE_MATRIX_OR_EXPRESSION_THAT_YOU_PASSED_DOES_NOT_HAVE_THE_EXPECTED_TYPE)
      eigen_assert(NumBlocks!=Dynamic);
    }

    template<typename OriginalMatrixType>
    inline BlockDiag(const OriginalMatrixType& a_matrix, Index numBlocks)
      : m_matrix(a_matrix), m_numBlocks(numBlocks)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<typename internal::remove_const<MatrixType>::type,OriginalMatrixType>::value),
                          THE_MATRIX_OR_EXPRESSION_THAT_YOU_PASSED_DOES_NOT_HAVE_THE_EXPECTED_TYPE)
    }

    inline Index rows() const { return m_matrix.rows() * m_numBlocks.value(); }
    inline Index cols() const { return m_matrix.cols() * m_numBlocks.value(); }

    inline Scalar coeff(Index rowId, Index colId) const
    {
      if ( rowId / m_matrix.rows() !=  colId / m_matrix.cols() )
        return Scalar(0);
      // try to avoid using modulo; this is a pure optimization strategy
      const Index actual_row  = internal::traits<MatrixType>::RowsAtCompileTime==1 ? 0
                            : NumBlocks==1 ? rowId
                            : rowId%m_matrix.rows();
      const Index actual_col  = internal::traits<MatrixType>::ColsAtCompileTime==1 ? 0
                            : NumBlocks==1 ? colId
                            : colId%m_matrix.cols();

      return m_matrix.coeff(actual_row, actual_col);
    }
    template<int LoadMode>
    inline PacketScalar packet(Index rowId, Index colId) const
    {
     
      if ( rowId / m_matrix.rows() !=  colId / m_matrix.cols() )
          GISMO_ERROR("not implemented");
          
      const Index actual_row  = internal::traits<MatrixType>::RowsAtCompileTime==1 ? 0
                            : NumBlocks==1 ? rowId
                            : rowId%m_matrix.rows();
      const Index actual_col  = internal::traits<MatrixType>::ColsAtCompileTime==1 ? 0
                            : NumBlocks==1 ? colId
                            : colId%m_matrix.cols();

      return m_matrix.template packet<LoadMode>(actual_row, actual_col);
    }

    const _MatrixTypeNested& nestedExpression() const
    { 
      return m_matrix; 
    }

  protected:
    MatrixTypeNested m_matrix;
    const internal::variable_if_dynamic<Index, NumBlocks> m_numBlocks;
};

/**
  * \return an expression of the replication of \c *this
  *
  * Example: \include MatrixBase_blockDiag_int_int.cpp
  * Output: \verbinclude MatrixBase_blockDiag_int_int.out
  *
  * \sa VectorwiseOp::blockDiag(), DenseBase::blockDiag<int,int>(), class BlockDiag
  */
template<typename Derived>
const typename MatrixBase<Derived>::BlockDiagReturnType
MatrixBase<Derived>::blockDiag(Index numBlocks) const
{
  return BlockDiag<Derived,Dynamic>(derived(),numBlocks);
}




/**
  * \class BlockTranspose
  *
  * \brief Expression of block-wise transposition of a tiled matrix
  *
  * \param MatrixType the type of the object we are replicating
  *
  */
namespace internal {
template<typename MatrixType,int NumBlocks>
struct traits<BlockTranspose<MatrixType,NumBlocks> >
 : traits<MatrixType>
{
  typedef typename MatrixType::Scalar Scalar;
  typedef typename traits<MatrixType>::StorageKind StorageKind;
  typedef typename traits<MatrixType>::XprKind XprKind;
  enum {
    Factor = (NumBlocks==Dynamic) ? Dynamic : NumBlocks*NumBlocks
  };
  typedef typename nested<MatrixType,Factor>::type MatrixTypeNested;
  typedef typename remove_reference<MatrixTypeNested>::type _MatrixTypeNested;
  enum {
    RowsAtCompileTime = NumBlocks==Dynamic || int(MatrixType::RowsAtCompileTime)==Dynamic
                      ? Dynamic
                      : MatrixType::ColsAtCompileTime / NumBlocks,
    ColsAtCompileTime = NumBlocks==Dynamic || int(MatrixType::ColsAtCompileTime)==Dynamic
                      ? Dynamic
                      : NumBlocks * MatrixType::RowsAtCompileTime,
   //FIXME we don't propagate the max sizes !!!
    MaxRowsAtCompileTime = RowsAtCompileTime,
    MaxColsAtCompileTime = ColsAtCompileTime,
    IsRowMajor = MaxRowsAtCompileTime==1 && MaxColsAtCompileTime!=1 ? 1
               : MaxColsAtCompileTime==1 && MaxRowsAtCompileTime!=1 ? 0
               : (MatrixType::Flags & RowMajorBit) ? 1 : 0,
    Flags = (_MatrixTypeNested::Flags & HereditaryBits & ~RowMajorBit) | (IsRowMajor ? RowMajorBit : 0),
    CoeffReadCost = _MatrixTypeNested::CoeffReadCost
  };
};
}

template<typename MatrixType,int NumBlocks> class BlockTranspose
  : public internal::dense_xpr_base< BlockTranspose<MatrixType,NumBlocks> >::type
{
    typedef typename internal::traits<BlockTranspose>::MatrixTypeNested MatrixTypeNested;
    typedef typename internal::traits<BlockTranspose>::_MatrixTypeNested _MatrixTypeNested;
  public:

    typedef typename internal::dense_xpr_base<BlockTranspose>::type Base;
    EIGEN_DENSE_PUBLIC_INTERFACE(BlockTranspose)

    template<typename OriginalMatrixType>
    inline explicit BlockTranspose(const OriginalMatrixType& a_matrix)
      : m_matrix(a_matrix), m_numBlocks(NumBlocks)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<typename internal::remove_const<MatrixType>::type,OriginalMatrixType>::value),
                          THE_MATRIX_OR_EXPRESSION_THAT_YOU_PASSED_DOES_NOT_HAVE_THE_EXPECTED_TYPE)
      eigen_assert(NumBlocks!=Dynamic);
    }

    template<typename OriginalMatrixType>
    inline BlockTranspose(const OriginalMatrixType& a_matrix, Index numBlocks)
      : m_matrix(a_matrix), m_numBlocks(numBlocks)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<typename internal::remove_const<MatrixType>::type,OriginalMatrixType>::value),
                          THE_MATRIX_OR_EXPRESSION_THAT_YOU_PASSED_DOES_NOT_HAVE_THE_EXPECTED_TYPE)
    }

    inline Index rows() const { return m_matrix.cols() / m_numBlocks.value(); }
    // return m_rows;

    inline Index cols() const { return m_matrix.rows() * m_numBlocks.value(); }

    inline Index numBlocks() const { return m_numBlocks.value(); }
    // return m_matrix.cols() / rows();
    
    inline Scalar coeff(Index rowId, Index colId) const
    {
        const Index r = colId % m_matrix.rows();                  // actual row
        const Index c = (colId/m_matrix.rows()) * rows() + rowId; // actual col
        return m_matrix.coeff(r, c);
    }
    template<int LoadMode>
    inline PacketScalar packet(Index rowId, Index colId) const
    {
        const Index r = colId % m_matrix.rows();                  // actual row
        const Index c = (colId/m_matrix.rows()) * rows() + rowId; // actual col
        return m_matrix.template packet<LoadMode>(r, c);
    }

    const _MatrixTypeNested& nestedExpression() const
    { 
      return m_matrix; 
    }

  protected:
    MatrixTypeNested m_matrix;
    const internal::variable_if_dynamic<Index, NumBlocks> m_numBlocks;
};

/**
  * \return an expression of block-wise transposed tiled matrix
  *
  */
template<typename Derived>
const typename MatrixBase<Derived>::BlockTransposeReturnType
MatrixBase<Derived>::blockTranspose(Index numBlocks) const
{
  return BlockTranspose<Derived,Dynamic>(derived(),numBlocks);
}


} // namespace Eigen
