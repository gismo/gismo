/** @file gsMatrixAddons.h

    @brief Provides extra member declarations to the Eigen MatrixBase class

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

inline const internal::adjugate_impl<Derived> adjugate() const;

inline void adjugateInPlace();
    
typedef BlockDiag<Derived,Dynamic> BlockDiagReturnType;
inline const BlockDiagReturnType blockDiag(Index rowFactor) const;

typedef BlockTranspose<Derived,Dynamic> BlockTransposeReturnType;
inline const BlockTransposeReturnType blockTranspose(Index rowFactor) const;

template<typename IndicesType>
const RowSelection<Derived,IndicesType> selectRows(const IndicesType & ind) const;


/**
  * \brief Simple (inplace) Gauss elimination without any pivoting
  */
void gaussElim()
{
    Derived & M = derived();
    index_t piv = 0;
    const index_t nr = M.rows();
    const index_t nc = M.cols();
    for (index_t r=0; r!=nr; ++r)
    {
        piv = 0;
        while (piv!= nc && 0 == M(r, piv)) ++piv;               
        if (piv == nc ) continue;
        
        const index_t br = nr-r-1;
        const index_t bc = nc-piv-1;
        M.block(r+1, piv+1, br, bc).noalias() -=
            M.col(piv).tail(br) * M.row(r).tail(bc) / M(r, piv);
        M.col(piv).tail(br).setZero();
    }
}


/**
  * \brief Inversion for small matrices using Cramer's Rule
  */
inline Matrix<Scalar, Dynamic, Dynamic> cramerInverse() const
{
    const Derived & M = derived();
    eigen_assert(M.rows() == M.cols() && "Matrix is not square.");

    Matrix<Scalar, Dynamic, Dynamic> rvo1(M.rows(), M.rows());
    switch (M.rows())
    {
        case 1:
            rvo1 = M.template topLeftCorner<1, 1>().inverse();
            break;
        case 2:
            rvo1 = M.template topLeftCorner<2, 2>().inverse();
            break;
        case 3:
            rvo1 = M.template topLeftCorner<3, 3>().inverse();
            break;
        case 4:
            rvo1 = M.template topLeftCorner<4, 4>().inverse();
            break;
        default:
            eigen_assert(false && "Not implemented.");
    };
    return rvo1;
}

/**
  * \brief Inplace inversion for small matrices using Cramer's Rule
  */
inline void cramerInverseInPlace()
{
//    derived() = cramerInverse().eval();
    derived().swap(cramerInverse().eval());
}