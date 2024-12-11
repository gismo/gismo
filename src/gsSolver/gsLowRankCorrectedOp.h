/** @file gsLowRankCorrectedOp.h

    @brief Provides the inverse of \f$Ainv^{-1} + U Q^{-1} V^T\f$ with \f$U Q^{-1} V^T\f$ being a low rank matrix using the Sherman Morrisson Woodburry formula

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): S. Takacs
*/

#pragma once

#include <gsCore/gsExport.h>
#include <gsSolver/gsLinearOperator.h>

namespace gismo
{

/// @brief Class for representing the inverse of a low rank corrected operator
///
/// Represents the inverse of \f$Ainv^{-1} + U Q^{-1} V^T\f$ with \f$U Q^{-1} V^T\f$ being a low rank matrix
///
/// \ingroup Solver
template<typename T = real_t>
class GISMO_EXPORT gsLowRankCorrectedOp : public gsLinearOperator<T>
{
public:

    /// Shared pointer for gsLowRankCorrectedOp
    typedef memory::shared_ptr<gsLowRankCorrectedOp> Ptr;

    /// Unique pointer for gsLowRankCorrectedOp
    typedef memory::unique_ptr<gsLowRankCorrectedOp> uPtr;
    
    /// Base class
    typedef typename gsLinearOperator<T>::Ptr BasePtr;

    /// Constructor
    gsLowRankCorrectedOp(const BasePtr& Ainv, const gsSparseMatrix<T>& Q, const gsSparseMatrix<T>& U, const gsSparseMatrix<T>& V);
    
    /// Make command returning a smart pointer
    static uPtr make(const BasePtr& Ainv, const gsSparseMatrix<T>& Q, const gsSparseMatrix<T>& U, const gsSparseMatrix<T>& V)
    {
        return memory::make_unique( new gsLowRankCorrectedOp( Ainv, Q, U, V ) );
    }

    virtual void apply(const gsMatrix<T> & input, gsMatrix<T> & result) const;

    virtual index_t rows() const    { return m_Ainv->rows(); }
    virtual index_t cols() const    { return m_Ainv->cols(); }

private:
    // disable copy constructor and assignment operator since for now we have no good
    // way of cloning the linear operators
    gsLowRankCorrectedOp(const gsLowRankCorrectedOp& other);
    gsLowRankCorrectedOp& operator=(const gsLowRankCorrectedOp& other);

private:

    gsSparseMatrix<T> m_U;
    gsSparseMatrix<T> m_V;
    BasePtr m_Ainv;
    BasePtr m_Winv;
    
};

} // end of namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsLowRankCorrectedOp.hpp)
#endif
