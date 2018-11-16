/** @file BlockTranspose.h

    @brief BlockTranspose extension for Eigen matrix objects

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/


namespace Eigen { 

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
  typedef typename ref_selector<MatrixType>::type MatrixTypeNested;
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

    // FIXME enable DirectAccess with negative strides?
    Flags = IsRowMajor ? RowMajorBit : 0
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
    typedef typename internal::remove_all<MatrixType>::type NestedExpression;
    
    template<typename OriginalMatrixType>
    EIGEN_DEVICE_FUNC
    inline explicit BlockTranspose(const OriginalMatrixType& a_matrix)
      : m_matrix(a_matrix), m_numBlocks(NumBlocks)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<typename internal::remove_const<MatrixType>::type,OriginalMatrixType>::value),
                          THE_MATRIX_OR_EXPRESSION_THAT_YOU_PASSED_DOES_NOT_HAVE_THE_EXPECTED_TYPE)
      eigen_assert(NumBlocks!=Dynamic);
    }

    template<typename OriginalMatrixType>
    EIGEN_DEVICE_FUNC
    inline BlockTranspose(const OriginalMatrixType& a_matrix, Index numBlocks)
      : m_matrix(a_matrix), m_numBlocks(numBlocks)
    {
      EIGEN_STATIC_ASSERT((internal::is_same<typename internal::remove_const<MatrixType>::type,OriginalMatrixType>::value),
                          THE_MATRIX_OR_EXPRESSION_THAT_YOU_PASSED_DOES_NOT_HAVE_THE_EXPECTED_TYPE)
    }

    EIGEN_DEVICE_FUNC
    inline Index rows() const { return m_matrix.cols() / m_numBlocks.value(); }
    // return m_rows;

    EIGEN_DEVICE_FUNC
    inline Index cols() const { return m_matrix.rows() * m_numBlocks.value(); }

    EIGEN_DEVICE_FUNC
    inline Index numBlocks() const { return m_numBlocks.value(); }
    // return m_matrix.cols() / rows();

    /*
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
struct unary_evaluator<BlockTranspose<ArgType, NumBlocks> >
  : evaluator_base<BlockTranspose<ArgType, NumBlocks> >
{
  typedef BlockTranspose<ArgType, NumBlocks> XprType;
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
      m_cols(xpr.nestedExpression().cols()),
      m_str(m_cols.value() / xpr.numBlocks())
  {}
 
  EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE
  CoeffReturnType coeff(Index rowId, Index colId) const
  {
      const Index r = colId % m_rows.value();                 // actual row
      const Index c = (colId/m_rows.value()) * m_str + rowId; // actual col
      return m_argImpl.coeff(r, c);
  }
  
  template<int LoadMode, typename PacketType>
  EIGEN_STRONG_INLINE
  PacketType packet(Index rowId, Index colId) const
  {
      const Index r = colId % m_rows.value();                 // actual row
      const Index c = (colId/m_rows.value()) * m_str + rowId; // actual col
      return m_argImpl.template packet<LoadMode>(r, c);
  }
   
protected:
  ArgTypeNested m_arg;
  evaluator<ArgTypeNestedCleaned> m_argImpl;
  const variable_if_dynamic<Index, ArgType::RowsAtCompileTime> m_rows;
  const variable_if_dynamic<Index, ArgType::ColsAtCompileTime> m_cols;
  const Index m_str ;
};

} // namespace internal

/** \memberof Eigen::MatrixBase
  * \return an expression of block-wise transposed tiled matrix
  */
template<typename Derived>
const typename MatrixBase<Derived>::BlockTransposeReturnType
MatrixBase<Derived>::blockTranspose(Index numBlocks) const
{
  return BlockTranspose<Derived,Dynamic>(derived(),numBlocks);
}


}
