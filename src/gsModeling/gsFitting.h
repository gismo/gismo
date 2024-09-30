/** @file gsFitting.h

    @brief Provides declaration of data fitting algorithms by least
    squares approximation.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): M. Kapl, G. Kiss, A. Mantzaflaris, D. Mokris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <vector>
#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMappedSpline.h>

namespace gismo
{

/**
  @brief 
   Class for performing a least squares fit of a parametrized point cloud with a gsGeometry.
    
   \ingroup Modeling
**/
template<class T>
class gsFitting
{
public:
    /// default constructor
    gsFitting()
    {
        m_basis = NULL;
        m_result= NULL ;
    }

    /// constructor
    gsFitting(gsMatrix<T> const & param_values, 
              gsMatrix<T> const & points, 
              gsBasis<T>  & basis);

        /// constructor
    gsFitting(gsMatrix<T> const & param_values, 
              gsMatrix<T> const & points,
              gsVector<index_t>  offset,
              gsMappedBasis<2,T>  & mbasis) ;

    /// Destructor
    virtual ~gsFitting();

public:

    /// Computes the least squares fit for a gsBasis
    void compute(T lambda = 0);

    void parameterCorrection(T accuracy = 1e-8,
                             index_t maxIter = 10,
                             T tolOrth = 1e-6);

    /// Computes the euclidean error for each point
    void computeErrors();

    /// Computes the maximum norm error for each point
    void computeMaxNormErrors();

    /// Computes the approximation error of the fitted curve to the original point cloud
    void computeApproxError(T & error, int type = 0) const;

    ///return the errors for each point
    void get_Error(std::vector<T>& errors, int type = 0) const;

    /// Returns the minimum point-wise error from the pount cloud (or zero if not fitted)
    T minPointError() const { return m_min_error; }

    /// Returns the maximum point-wise error from the pount cloud (or zero if not fitted)
    T maxPointError() const { return m_max_error; }

    /// Return the errors for each point
    const std::vector<T> & pointWiseErrors() const
    {
        return m_pointErrors;
    }

    /// Computes the number of points below the error threshold (or zero if not fitted)
    size_t numPointsBelow(T threshold) const
    { 
        const size_t result=
            std::count_if(m_pointErrors.begin(), m_pointErrors.end(), 
                          GS_BIND2ND(std::less<T>(), threshold));
        return result; 
    }

    /// Computes the least squares fit for a gsBasis
    void iterativeCompute( T const & tolerance, unsigned const & num_iters = 10);

    /// Adds to the matrix A_mat terms for minimization of second derivative, weighted
    /// with parameter lambda.
    void applySmoothing(T lambda, gsSparseMatrix<T> & A_mat);
    gsSparseMatrix<T> smoothingMatrix(T lambda) const;
    /// Assembles system for the least square fit.
    void assembleSystem(gsSparseMatrix<T>& A_mat, gsMatrix<T>& B);


public:

    /// gives back the computed approximation
    gsGeometry<T> * result() const { return m_result; }

    /// gives back the computed approximation for multipatch geometry
    const gsMappedSpline<2,T> & mresult() const { return m_mresult; }

    /// Returns the basis of the approximation
    const gsBasis<T> & getBasis() const {return *static_cast<const gsBasis<T>*>(m_basis);}

    void setBasis(gsBasis<T> & basis) {m_basis=&basis;}

    /// returns the parameter values
    gsMatrix<T> & getreturnParamValues() {return m_param_values;}
    gsMatrix<T> & returnParamValues() {return m_param_values;}

    /// returns the points
    gsMatrix<T> returnPoints() const {return m_points;}

    /// Sets constraints that the coefficients of the resulting
    /// geometry have to conform to. More precisely, denoting the
    /// coefficient vector by \a x, it enforces
    ///  \a lhs * \a x = \a rhs.
    void setConstraints(const gsSparseMatrix<T>& lhs, const gsMatrix<T>& rhs)
    {
	m_constraintsLHS = lhs;
	m_constraintsRHS = rhs;
    }

    /// Sets constraints on that the coefficients of the resulting geometry have to conform to.
    /// \param indices indices (in the coefficient vector) of the prescribed coefficients.
    /// \param coefs prescribed coefficients.
    void setConstraints(const std::vector<index_t>& indices,
			const std::vector<gsMatrix<T> >& coefs);

    /// Sets constraints in such a way that the previous values at \a
    /// fixedSides of the geometry remain intact.
    void setConstraints(const std::vector<boxSide>& fixedSides);

    /// Set constraints in such a way that the resulting geometry on
    /// each of \a fixedSides will coincide with the corresponding
    /// curve in \a fixedCurves.
    void setConstraints(const std::vector<boxSide>& fixedSides,
            const std::vector<gsBSpline<T> >& fixedCurves);
    void setConstraints(const std::vector<boxSide>& fixedSides,
            const std::vector<gsGeometry<T> * >& fixedCurves);


private:
    /// Extends the system of equations by taking constraints into account.
    void extendSystem(gsSparseMatrix<T>& A_mat, gsMatrix<T>& m_B);

protected:

    //gsOptionList

    /// the parameter values of the point cloud
    gsMatrix<T> m_param_values;

    /// the points of the point cloud
    gsMatrix<T> m_points;

    // Patch offsets
    gsVector<index_t> m_offset;

    /// Pointer keeping the basis
    gsFunctionSet<T> * m_basis;

    /// Pointer keeping the resulting geometry
    gsGeometry<T> * m_result;

    /// Pointer keeping the resulting multipatch geometry
    gsMappedSpline<2,T>  m_mresult;

    // All point-wise errors
    std::vector<T> m_pointErrors;

    mutable T m_last_lambda;

    /// Maximum point-wise error
    T m_max_error;

    /// Minimum point-wise error
    T m_min_error;

    /// Left hand-side of the constraints that the coefficients of the
    /// resulting geometry have to conform to.
    /// This corresponds to matrix D in Prautzch, Boehm, Paluszny:
    /// Bezier and B-spline techniques, Section 4.7.
    gsSparseMatrix<T> m_constraintsLHS;

    /// Right hand-side of the constraints that the coefficients of the
    /// resulting geometry have to conform to.
    /// This corresponds to vector q in Prautzch, Boehm, Paluszny:
    /// Bezier and B-spline techniques, Section 4.7.
    gsMatrix<T>       m_constraintsRHS;

private:
    //void applySmoothing(T lambda, gsMatrix<T> & A_mat);

}; // class gsFitting


#ifdef GISMO_WITH_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsKnotVector
   */
  void pybind11_init_gsFitting(pybind11::module &m);

#endif // GISMO_WITH_PYBIND11


}// namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsFitting.hpp)
#endif
