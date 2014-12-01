/** @file gsDivConSolution.h

    @brief Recontructs the solution vector field from the divergence preserving transformation.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Sogn
*/

#pragma once

#include <gsCore/gsFunction.h>
#include <gsCore/gsGeometry.h>

//----TODO----
//document


namespace gismo
{


template<class T, class Basis_t>
class gsDivConSolution : public gsFunction<T>
{

public:

    /// @name Constructors
    /// @{

    /// Default empty constructor
    gsDivConSolution() { }

    /// Constructor which copies the given coefficient matrix \a SolutionCoefs.
    gsDivConSolution(const gsMatrix<T> & solutionCoefs, const gsGeometry<T> & geo, std::vector< gsBasis<T> *> const & basis) :
        m_solutionCoeff(solutionCoefs), m_geometry(geo), m_basis(basis)
    {
        componentShifts.resize(targetDim());
        componentShifts[0] = 0;
        for (index_t k = 1; k < targetDim(); ++k)
        {
            componentShifts[k] = componentShifts[k-1] + m_basis[k-1]->size();
        }
    }


    /// @}

public:

    /*/** @name Evaluation functions

        These functions allow to evaluate the geometry as well as its derivatives
        at one or multiple points of the parameter space.
        All evaluation functions of gsFunction, from which gsGeometry derives,
        are also supported.

        \note
        These functions generally do not have to be overridden in
        derived classes since the basis type will provide the proper implementation.

        @{
    */



    /// Evaluates the solution into result

    void eval_into(const gsMatrix<T>& u, gsMatrix<T>& result) const
    {
        //JS2: 2D Verified, 3D not Verified(Tested)

        const index_t TarDim = targetDim();
        GISMO_ASSERT(m_geometry.geoDim() == TarDim, "Geometric dimention and target dimention not matching!");

        const index_t numPts = u.cols();
        gsMatrix<T> B, jacobi;
        gsMatrix<unsigned> ind;

        result.setZero( TarDim, numPts );

        std::auto_ptr< gsGeometryEvaluator<T> > geoEval (m_geometry.evaluator(NEED_VALUE | NEED_MEASURE | NEED_GRAD_TRANSFORM));
        geoEval->evaluateAt(u);
        jacobi = geoEval->jacobians();

        for (index_t comp = 0; comp< TarDim; ++comp)
        {
            m_basis[comp]->eval_into(u, B);
            m_basis[comp]->active_into(u, ind);
            for ( index_t j=0; j< numPts ; j++ ) // for all points (columns of u)
            {
                T det = (jacobi.block(0, j*TarDim, TarDim, TarDim)).determinant();
                for ( index_t i = 0; i < ind.rows(); ++i ) // for all non-zero basis functions
                {
                    result.col(j) += jacobi.block(0, comp+j*TarDim, TarDim, 1)
                            * m_solutionCoeff(ind(i,j) + componentShifts[comp], 0) * B(i,j)/det;//geoEval.measure(j);
                }
            }
        }
    }

    /// @}
    /*************************************************************************/

    /// @name Accessors
    /// @{

    /// \brief Returns the basis.
    std::vector<gsBasis<T> * > basis() {return m_basis;}

    /// Dimension \em n of the physical space
    int geoDim() const {return m_geometry.geoDim();}

    /// Dimension \em n of the coefficients (control points)
    int coefDim() const { return m_geometry.coefDim(); }

    /// Dimension of the absent physical space (overriding gsFunction::targetDim())
    int targetDim() const { return m_basis.size();}

    /// Dimension \em d of the parameter domain.
    virtual int parDim() const { return m_basis[0]->dim(); }


    /// Dimension \em d of the parameter domain (overriding gsFunction::domainDim()).
    int domainDim() const { return m_basis[0]->dim(); }
    /// @}



protected:

    // The coefficient of the solution field (at one pacht), size of basis X targetDim.
    gsMatrix<T> m_solutionCoeff;

    // The geometry, need this for the jacobian.
    const gsGeometry<T>  & m_geometry;//(&)

    // The basis
    std::vector< gsBasis<T> * >m_basis;

    // shift index for components in m_solutionCoeff;
    std::vector<index_t> componentShifts;


}; // class gsDivConSolution


}
