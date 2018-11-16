/** @file BlockDiag.h

    @brief BlockDiag extension for Eigen matrix objects

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


namespace Eigen { 

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
  typedef typename ref_selector<MatrixType>::type MatrixTypeNested;
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

    // FIXME enable DirectAccess with negative strides?
    Flags = IsRowMajor ? RowMajorBit : 0
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
    typedef typename internal::remove_all<MatrixType>::type NestedExpression;
    
    template<typename OriginalMatrixType>
    EIGEN_DEVICE_FUNC
    inline explicit BlockDiag(const OriginalMatrixType& a_matrix)
      : m_matrix(a_matrix), m_numBlocks(NumBlocks)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<typename internal::remove_const<MatrixType>::type,OriginalMatrixType>::value),
                          THE_MATRIX_OR_EXPRESSION_THAT_YOU_PASSED_DOES_NOT_HAVE_THE_EXPECTED_TYPE)
      eigen_assert(NumBlocks!=Dynamic);
    }

    template<typename OriginalMatrixType>
    EIGEN_DEVICE_FUNC
    inline BlockDiag(const OriginalMatrixType& a_matrix, Index numBlocks)
      : m_matrix(a_matrix), m_numBlocks(numBlocks)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<typename internal::remove_const<MatrixType>::type,OriginalMatrixType>::value),
                          THE_MATRIX_OR_EXPRESSION_THAT_YOU_PASSED_DOES_NOT_HAVE_THE_EXPECTED_TYPE)
    }

    EIGEN_DEVICE_FUNC
    inline Index rows() const { return m_matrix.rows() * m_numBlocks.value(); }

    EIGEN_DEVICE_FUNC
    inline Index cols() const { return m_matrix.cols() * m_numBlocks.value(); }
/*
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
*/
    EIGEN_DEVICE_FUNC
    const _MatrixTypeNested& nestedExpression() const
    { 
      return m_matrix; 
    }

  protected:
    MatrixTypeNested m_matrix;
    const internal::variable_if_dynamic<Index, NumBlocks> m_numBlocks;
};

namespace internal {
template<typename ArgType, int NumBlocks> 
struct unary_evaluator<BlockDiag<ArgType, NumBlocks> >
  : evaluator_base<BlockDiag<ArgType, NumBlocks> >
{
  typedef BlockDiag<ArgType, NumBlocks> XprType;
  typedef typename XprType::CoeffReturnType CoeffReturnType;
  enum {
    Factor = (NumBlocks==Dynamic) ? Dynamic : NumBlocks*NumBlocks
  };
  typedef typename internal::nested_eval<ArgType,Factor>::type ArgTypeNested;
  typedef typename internal::remove_all<ArgTypeNested>::type ArgTypeNestedCleaned;
  
  enum {
    CoeffReadCost = evaluator<ArgTypeNestedCleaned>::CoeffReadCost,
    LinearAccessMask = XprType::IsVectorAtCompileTime ? LinearAccessBit : 0,
    Flags = (evaluator<ArgTypeNestedCleaned>::Flags & (HereditaryBits|LinearAccessMask) & ~RowMajorBit) | (traits<XprType>::Flags & RowMajorBit),
    
    Alignment = evaluator<ArgTypeNestedCleaned>::Alignment
  };

  EIGEN_DEVICE_FUNC explicit unary_evaluator(const XprType& xpr)
    : m_arg(xpr.nestedExpression()),
      m_argImpl(m_arg),
      m_rows(xpr.nestedExpression().rows()),
      m_cols(xpr.nestedExpression().cols())
    {}
 
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  CoeffReturnType coeff(Index rowId, Index colId) const
  {
      if ( rowId / m_rows.value() !=  colId / m_cols.value() )
        return CoeffReturnType(0);
      // try to avoid using modulo; this is a pure optimization strategy
      const Index actual_row  = internal::traits<XprType>::RowsAtCompileTime==1 ? 0
                            : NumBlocks==1 ? rowId
                            : rowId%m_rows.value();
      const Index actual_col  = internal::traits<XprType>::ColsAtCompileTime==1 ? 0
                            : NumBlocks==1 ? colId
                            : colId%m_cols.value();

      return m_argImpl.coeff(actual_row, actual_col);
  }
  
  template<int LoadMode, typename PacketType>
  EIGEN_STRONG_INLINE
  PacketType packet(Index rowId, Index colId) const
  {
      assert( rowId / m_rows.value() !=  colId / m_cols.value() &&
              "Not implemented");
          
      const Index actual_row  = internal::traits<XprType>::RowsAtCompileTime==1 ? 0
                            : NumBlocks==1 ? rowId
                            : rowId%m_rows.value();
      const Index actual_col  = internal::traits<XprType>::ColsAtCompileTime==1 ? 0
                            : NumBlocks==1 ? colId
                            : colId%m_cols.value();

      return m_argImpl.template packet<LoadMode>(actual_row, actual_col);
  }
   
protected:
  ArgTypeNested m_arg;
  evaluator<ArgTypeNestedCleaned> m_argImpl;
  const variable_if_dynamic<Index, ArgType::RowsAtCompileTime> m_rows;
  const variable_if_dynamic<Index, ArgType::ColsAtCompileTime> m_cols;
};

} // namespace internal

/** \memberof Eigen::MatrixBase
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

} // namespace eigen
