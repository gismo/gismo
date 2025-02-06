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
#include <gsMSplines/gsMappedBasis.h>
#include <gsMSplines/gsMappedSpline.h>

namespace gismo
{

/**
  @brief
   Class for performing a fit of a parametrized point cloud with a gsGeometry.

   \ingroup Modeling
**/
template<class T>
class gsFitting
{
public:
    /// default constructor
    gsFitting()
    {
        m_basis = nullptr;
        m_result= nullptr;
    }

    /** @brief gsFitting: Main constructor of the fitting class
    * @param param_values a matrix containing the parameter values that parametrize the \a points
    * @param points matrix containing the points to be fitted
    * @param basis basis to use for fitting
    */
    gsFitting(gsMatrix<T> const & param_values,
              gsMatrix<T> const & points,
              gsBasis<T>  & basis);

    /// @brief gsFitting: Main constructor of the fitting class for multi-patch
    gsFitting(gsMatrix<T> const & param_values,
              gsMatrix<T> const & points,
              gsVector<index_t>  offset,
              gsMappedBasis<2,T>  & mbasis) ;

    /// Destructor
    virtual ~gsFitting();

public:

    /** @brief compute: Computes the coefficients of the spline geometry via penalized least squares
    * @param lambda smoothing weight
    */
    void compute(T lambda = 0);

    /** @brief updateGeometry: Updates the fitted geometry with new coefficients and parameters
    * @param coefficients the new coefficients
    * @param parameters the new parameters
    */
    void updateGeometry(gsMatrix<T> coefficients, gsMatrix<T> parameters);


    /** @brief initializeGeometry: Initializes the fitted geometry with given coefficients and parameters
    * @param coefficients the input coefficients
    * @param parameters the input parameters
    */
    void initializeGeometry(const gsMatrix<T> & coefficients, const gsMatrix<T> & parameters);

    /// choose the method for computing the coefficients: TDM, PDM, HDM with different blending weights
    enum tdm_method
    {
        tdm_boundary_pdm, // TDM
        pdm, // PDM
        hybrid_pdm_tdm_boundary_pdm, // HDM with constant weight
        hybrid_error_pdm_tdm_boundary_pdm, // HDM with weight based on point-wise error
        hybrid_curvature_pdm_tdm_boundary_pdm, // HDM with weight based on curvature
    };

    /**
     * @brief compute_hdm: computes the coefficients of the spline geometry via Hybrid Distance Minimization (HDM)
     * @param lambda smoothing weight
     * @param mu weight for the PDM
     * @param sigma weight for the TDM
     * @param interpIdx vector containing the number of interior points and the indices of the boundary points
     * @param method method for computing the blending weights
     */
    void compute_tdm(T lambda, T mu, T sigma, const std::vector<index_t> & interpIdx,
                     tdm_method method = hybrid_curvature_pdm_tdm_boundary_pdm);


    /// check if the given parameter \a parameter is a corner of the domain \a parametric_domain
    bool is_corner(gsMatrix<T> & parametric_domain, gsVector<T> & parameter);



    /// check if the given parameter \a parameter is within the cell \a element
    bool is_point_within_cell(const gsMatrix<T>& parameter, const gsMatrix<T>& element);
    /// check if the given parameter  \a x, \a y is within the cell \a element; same as \ref is_point_within_cell, but different input format
    bool is_point_within_cell(const T x, const T y, const gsMatrix<T>& element);

    /// check if the given parameter \a parameter is inside the support \a support
    /// difference with \ref is_point_inside_cell in the inclusion of the left and right interval extremes.
    bool is_point_inside_support(const gsMatrix<T>& parameter, const gsMatrix<T>& support);
    /// check if the given parameter \a x, \a y is inside the support \a support; same as \ref is_point_inside_support, but different input format
    bool is_point_inside_support(const T x, const T y, const gsMatrix<T>& support);


    /** @brief parameterCorrection: globally apply \a maxIter steps of parameter correction to the least squares fitted geometry
     * @param accuracy accuracy of the closest point computation
     * @param maxIter maximum number of parameter correction steps
     * @param tolOrth orthogonality tolerance
     */
    void parameterCorrection(T accuracy = 1e-8,
                             index_t maxIter = 10,
                             T tolOrth = 1e-6);

    /** @brief parameterProjectionSepBoundary: project the points onto the fitted geometry, separating interior and boundary points
     * @param accuracy accuracy of the closest point computation, for the foot-point projection
     * @param interpIdx vector containing the number of interior points and the indices of the boundary points
     */
    void parameterProjectionSepBoundary(T accuracy,const std::vector<index_t>& interpIdx);

    /** @brief parameterCorrectionSepBoundary_pdm: apply \a maxIter steps of parameter correction for PDM method, separating interior and boundary points
     * @param accuracy accuracy of the closest point computation
     * @param maxIter maximum number of parameter correction steps
     * @param sepIndex vector containing the number of interior points and the indices of the boundary points
     */
    void parameterCorrectionSepBoundary_pdm(T accuracy, index_t maxIter, const std::vector<index_t>& sepIndex);

    /** @brief parameterCorrectionSepBoundary_tdm: apply \a maxIter steps of parameter correction for HDM method, separating interior and boundary points
    * @param accuracy accuracy of the closest point computation
    * @param maxIter maximum number of parameter correction steps
    * @param mu weight for PDM
    * @param sigma weight for TDM
    * @param sepIndex vector containing the number of interior points and the indices of the boundary points
    * @param method method for computing the blending weights
    */
    void parameterCorrectionSepBoundary_tdm(T accuracy, index_t maxIter, T mu, T sigma, const std::vector<index_t>& sepIndex, tdm_method method = hybrid_curvature_pdm_tdm_boundary_pdm);

    /// Computes the point-wise errors in euclidean norm as well as the max and min errors,
    /// and updates the member variables \a m_pointErrors, \a m_max_error, \a m_min_error;
    /// different from \ref computeMaxNormErrors(), where the error is computed in inifinity/maximum norm
    void computeErrors();

    /// Compute point-wise error in euclidean norm between the fitted geometry at the parameters \a parameters and the input point cloud \a points
    /// similar to \ref computeErrors(), but different input and output format; it does not update the member variables
    gsMatrix<T> pointWiseErrors(const gsMatrix<> & parameters,const gsMatrix<> & points);

    /// Computes min, max and mse errors in euclidean norms between the fitted geometry at the parameters \a param_values and the input point cloud \a points
    /// it does not update the member variables
    std::vector<T> computeErrors(const gsMatrix<> & param_values,const gsMatrix<> & points);

    /// Computes the point-wise errors in infinity/maximum norm as well as the max and min errors,
    /// and updates the member variables \a m_pointErrors, \a m_max_error, \a m_min_error;
    /// different from \ref computeErrors(), where the error is computed in euclidean norm
    void computeMaxNormErrors();

    /// Computes the approximation error \a error of the fitted geometry to the original point cloud
    /// \a type = 0: sum of squares, \a type = 1: sum of absolute values (l_1 norm)
    void computeApproxError(T & error, int type = 0) const;

    /// Compute the point-wise errors for each point
    /// \a type = 0: point-wise infinity/maximum norm
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

    /// Adds to the matrix A_mat terms for minimization of second derivative, weighted with parameter lambda.
    void applySmoothing(T lambda, gsSparseMatrix<T> & A_mat);
    gsSparseMatrix<T> smoothingMatrix(T lambda) const;
    /// Assembles system for the least square fit.
    void assembleSystem(gsSparseMatrix<T>& A_mat, gsMatrix<T>& B);

    /** @brief compute_normals: Computes the normals of the fitted geometry at the input parameter values \a params_int
    * @param num_int index of the input parameter values
    * @param params_int input parameter values
    * @param N_int matrix containing the normals
    */
    void compute_normals(const index_t & num_int, const gsMatrix<T> & params_int, gsSparseMatrix<T> & N_int);

    /// vector of length \a num_int containing all the point-wise errors; store also the max err value in \a max_err_int
    gsMatrix<T> fill_pointWiseErrors(const index_t & num_int, T & max_err_int);

    /// compute the principal curvatures (c1, c2) at the given parameters \a params
    gsMatrix<T> principal_curvatures(const gsMatrix<T> & params);

    /// vector of length \a num_int containing rho = 1/max(c1, c2), where c1, c2 are the principal curvature values computed at every parametric point \a params_int
    gsMatrix<T> inverse_principal_curvatures(const index_t & num_int, const gsMatrix<T> & params_int);

    /** @brief blending_weights: computes the blending weights \a mu and \a sigma for the balance mu * PDM + sigma * TDM in the HDM method
    * @param N_int matrix containing the nomals at the parameters \a params_int
    * @param num_int indeces of the interior parameters
    * @param mu weight for PDM
    * @param sigma weight for TDM
    * @param params_int input parameter values
    * @param method method for computing the blending weights: constant, based on point-wise error, based on curvature
    * @param NNT output matrix containing the normals and the blending weights
     */
    void blending_weights(const gsSparseMatrix<T> & N_int, const index_t & num_int, const T & mu, const T & sigma,
                          const gsMatrix<T> & params_int, tdm_method method, gsSparseMatrix<T> & NNT);

    /** @brief assembleSystem: assembles the linear system for the Hybrid Distance Minimization method
    * @param points_int interior points
    * @param params_int interior parameters
    * @param points_bdy boundary points
    * @param params_bdy boundary parameters
    * @param num_basis dimension of the basis
    * @param NNT matrix containing the normals and the blending weights
    * @param A_tilde output system matrix
    * @param rhs output right-hand side vector
    */
    void assembleSystem(const gsMatrix<T> & points_int, const gsMatrix<T> & params_int,
                        const gsMatrix<T> & points_bdy, const gsMatrix<T> & params_bdy,
                        const index_t & num_basis, const gsSparseMatrix<T> & NNT,
                        gsSparseMatrix<T> & A_tilde, gsMatrix<T> & rhs);

public:

    /// gives back the computed approximation
    gsGeometry<T> * result() const { return m_result; }

    /// gives back the computed approximation for multipatch geometry
    const gsMappedSpline<2,T> & mresult() const { return m_mresult; }

    /// Returns the basis of the approximation
    const gsBasis<T> & getBasis() const {return *static_cast<const gsBasis<T>*>(m_basis);}

    void setBasis(gsBasis<T> & basis) {m_basis=&basis;}

    /// returns the parameter values
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

    /// Initialize the parametric domain of the point cloud
    void initParametricDomain()
    {
        m_uMin = m_param_values.row(0).minCoeff();
        m_uMax = m_param_values.row(0).maxCoeff();
        m_vMin = m_param_values.row(1).minCoeff();
        m_vMax = m_param_values.row(1).maxCoeff();

        gsInfo << "Parametric domain: ["
               << m_uMin << ", " << m_uMax << "] x ["
               << m_vMin << ", " << m_vMax << "]" << std::endl;
    }

    /// Returns the smoothing weight used in the last fitting
    T lambda() const {return m_last_lambda;}


private:
    /// Extends the system of equations by taking constraints into account.
    void extendSystem(gsSparseMatrix<T>& A_mat, gsMatrix<T>& m_B);

protected:

    /// Assembles 3xblock collocation matrix.
    void assembleBlockB(const gsMatrix<T>& points,
                        const gsMatrix<T>& params,
                        index_t num_basis,
                        gsSparseMatrix<T>& result) const
    {
        GISMO_UNUSED(points);
        GISMO_UNUSED(num_basis);
        //index_t num_pts = points.rows();
        gsSparseMatrix<T> sparseColloc = m_result->basis().collocationMatrix(params);
        threeOnDiag(sparseColloc, result);
    }

    /// Assembles the right hand side vectors for PDM/TDM.
    void assembleBlockX(const gsMatrix<T>& points,
                        gsMatrix<T>& result) const
    {
        result.resize(points.rows() * 3, 1);
        result << points.col(0), points.col(1), points.col(2);
    }

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

    // Interior c1 and c2 curvature, as column of the matrix
    gsMatrix<T> m_pointCurvature;

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

    T m_uMin, m_uMax, m_vMin, m_vMax;

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
