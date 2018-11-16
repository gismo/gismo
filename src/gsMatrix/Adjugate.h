/** @file Adjugate.h

    @brief Adjugate extension for Eigen matrix objects

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
    static inline void run(const MatrixType& , ResultType& result)
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

/// @memberof Eigen::MatrixBase
template<typename Derived>
inline const internal::adjugate_impl<Derived>
MatrixBase<Derived>::adjugate() const
{
    eigen_assert(rows() == cols());
    return internal::adjugate_impl<Derived>(derived());
}

/// @memberof Eigen::MatrixBase
template<typename Derived>
inline void MatrixBase<Derived>::adjugateInPlace()
{
    derived() = adjugate().eval();
}

}
