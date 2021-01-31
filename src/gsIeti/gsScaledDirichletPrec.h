/** @file gsScaledDirichletPrec.h

    @brief This class represents the sclaed Dirichlet preconditioner.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): S. Takacs
*/

#pragma once

#include <gsSolver/gsLinearOperator.h>
#include <gsSolver/gsMatrixOp.h>
#include <gsSolver/gsBlockOp.h>
#include <gsSolver/gsProductOp.h>
#include <gsSolver/gsSumOp.h>
#include <gsSolver/gsAdditiveOp.h>

#define DEBUGVAR(a) gsInfo << "  " << #a << ": " << a << std::endl
#define DEBUGMATRIX(a) gsInfo << "  " << #a << ": " << a.rows() << " x " << a.cols() << std::endl

namespace gismo
{

/** @brief   This class represents the sclaed Dirichlet preconditioner.

    \ingroup Solver
*/
template< typename T >
class gsScaledDirichletPrec
{
    typedef typename gsLinearOperator<T>::Ptr OpPtr;
    typedef gsSparseMatrix<T,RowMajor> Transfer;
    typedef memory::shared_ptr<Transfer> TransferPtr;
    typedef gsMatrix<T> Matrix;
public:

//TODO: needed?
/*    void setupSparseLUSolvers()
    {
        const size_t sz = localMatrixOps.size();
        localSolverOps.clear();
        localSolverOps.reserve(sz);
        for (size_t i=0; i<sz; ++i)
        {
            gsMatrixOp< gsSparseMatrix<T> >* matop
              = dynamic_cast< gsMatrixOp< gsSparseMatrix<T> >* >(localMatrixOps[i].get());
            GISMO_ENSURE( matop, "gsScaledDirichletPrec::setupSparseLUSolvers requires the "
              "local systems in localMatrixOps to by of type gsMatrixOp< gsSparseMatrix<T> >." );
            localSolverOps.push_back(makeSparseLUSolver(gsSparseMatrix<T>(matop->matrix())));
        }
    }
*/    
    
    index_t numberOfLagrangeMultipliers() const
    {
        GISMO_ASSERT( jumpMatrices.size()>0, "gsScaledDirichletPrec: Number of Lagrange multipliers "
            "can only be determined if there are jump matrices.");
        return jumpMatrices[0]->rows();
    }
    
    void setupMultiplicityScaling()
    {
        GISMO_ASSERT( jumpMatrices.size() == localMatrixOps.size(),
            "TODO" );

        const size_t pnr = jumpMatrices.size();
        localScaling.clear();
        localScaling.reserve(pnr);

        for (size_t k=0; k<pnr; ++k)
        {
            const size_t sz = localMatrixOps[k]->rows();
            gsMatrix<T> sc(sz,1);
            for (index_t i=0; i<sz; ++i)
              sc(i,0) = 1;

            Transfer & jm = *(jumpMatrices[k]);

            for (index_t i=0; i<jm.outerSize(); ++i){
                for (typename Transfer::InnerIterator it(jm, i); it; ++it) {
                    const index_t c = it.col();
                    sc(c,0) += 1;
                }
            }
            localScaling.push_back(give(sc));
        }
    }

    OpPtr secaledDirichletPreconditioner() const
    {
        GISMO_ASSERT( jumpMatrices.size() == localScaling.size(),
            "gsScaledDirichletPrec::secaledDirichletPreconditioner needs the localScaling matrices given. "
            "Forgot to call setupMultiplicityScaling()?" );
        GISMO_ASSERT( jumpMatrices.size() == localSchurOps.size(),
            "gsScaledDirichletPrec::secaledDirichletPreconditioner needs the localSchur operators given." );

        const size_t pnr = jumpMatrices.size();

        std::vector<OpPtr> scalingOps;
        scalingOps.reserve(pnr);

        for (size_t k=0; k<pnr; ++k)
        {
            const size_t sz = localMatrixOps[k]->rows();
            gsSparseMatrix<T> scaling(sz,sz);
            for (size_t i=0; i<sz; ++i)
                scaling(i,i) = 1./localScaling[k](i,0);
            scalingOps.push_back(makeMatrixOp(scaling.moveToPtr()));
        }

        typename gsSumOp<T>::Ptr result = gsSumOp<T>::make();

        for (size_t i=0; i<pnr; ++i)
        {
            typename gsProductOp<T>::Ptr local = gsProductOp<T>::make();
            local->addOperator(makeMatrixOp(jumpMatrices[i]->transpose()));
            local->addOperator(scalingOps[i]);
            local->addOperator(localSchurOps[i]);
            local->addOperator(scalingOps[i]);
            local->addOperator(makeMatrixOp(jumpMatrices[i]));
            result->addOperator(local);
        }

        return result;
    }

public:
    std::vector< OpPtr >        localMatrixOps;
    std::vector< OpPtr >        localSchurOps;
    std::vector< Matrix >       localScaling;
    std::vector< TransferPtr >  jumpMatrices;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsScaledDirichletPrec.hpp)
#endif
