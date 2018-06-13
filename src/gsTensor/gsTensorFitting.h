/** @file gsTensorFitting.h

    @brief Provides declaration of an iterative fitting method.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Mayr
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>

namespace gismo
{
	
	/**
	 * @brief
	 * Class for performing an iterative fitting method, 
	 * including a fitting step, a knot removing step and a projection step.
	 **/ 	
	
	template<class T>
	class gsTensorFitting
	{
		public:
		
		/// default constructor
		gsTensorFitting(){}
		
		/** constructor
		 * @param[in] data_points	3x(m_1*m_2)-matrix
		 * @param[in] m_1	number of data points in first direction
		 * @param[in] m_2	number of data points in second direction
		 * @param[in] degree	maximal degree
		 * @param[in] error_bound	bound for knot removal
		 */ 
		gsTensorFitting(const gsMatrix<T>& data_points,
						const index_t& m_1, 
						const index_t& m_2,
						const index_t& degree, 
						const T& error_bound); 
		
		/// Destructor
		virtual ~gsTensorFitting()
		{
			
		}
		
		/// set the tolerance (error bound) for removing knots
		void setTolerance(T& error_bound)
		{
			E_ = error_bound; 
		}
		
		/// return the current tolerance (error bound) for removing knots
		T tolerance() const { return E_; }
		
		/// set the maximal degree
		void setMaxDegree(index_t& max_degree)
		{
			p_ = max_degree; 
		}
		
		/// return the maximal degree
		index_t maxDegree() const { return p_; }
		
		/// returns the computed TensorBSpline
		gsTensorBSpline<2,T>& result() {return Tensor_; }
		
		/// returns the meanPointError; throws if .compute() has not been called yet
		T meanPointError() const {return mean_error_; }
		
		/// returns the maxPointError; throws if .compute() has not been called yet
		T maxPointError() const {return max_error_; }
		
		///computes the BSplineSurface
		void compute(); 
		
		private:
		
		//methods
		
		//----------------------------------------------------------------------
		//constructing knot vector and parameter values
		/** computes distance of two points in 3d
		* @param[in] a	first point
		* @param[in] b	second point
		* @return d	euclidean distance between a and b
		* @warning	assumes a.size = b.size = 3
		*/
		T dist_3d(const gsVector<T>& a, const gsVector<T>& b); 
		
		/** computation of parameter values for approximation
		* cf nurbs book, p. 365, (9.6)
		* @param[in] dir	indicating if direction 0 or 1
		* @param[out] u_k	size = n_p and filled with 0s required, containing n_p parameter values
		*/
		void param_vals_centripetal(const int dir, std::vector<T>& u_k);
		
		/** constructs knot vector out of the values contained in another vector 
		* @param[in] vec	contains knot-values u_0, .., u_m
		* @param[in] m	size of vec
		* @return knot vector with 0, 0, u_1, .., u_m-1, 1, 1
		*/ 
		gsKnotVector<T> get_knotvec_from_vec(const std::vector<T>& vec, const int m); 
		
		//----------------------------------------------------------------------
		//remove knots from surface
		/** computes for a given basis function the range of parameter indices (not the range of knots)
		* as required on p.429 of nurbs book
		* @param[in] idx_tmp_begin	indicates at which basis functions to start, idx=0,...,n_t-1 (independent of direction)
		* @param[in] idx_tmp_end	indicates at which basis function to end, st i < idx_tmp_end; last basis function is n_t-1
		* @param[in] kv_u	knot-vector
		* @param[in] n_t	number of basis functions in kv_u
		* @param[in] uk	vector containing parameter values
		* @param[out] idx_param_vals	vector of size 2*n_t, idx_param_vals[i]=a and idx_param_vals[i+1]=b implies that parameters uk(j) j=a,..,b-1 lie in support of basis function i
		*/ 
		void get_param_range(const size_t idx_tmp_begin, const size_t idx_tmp_end, const gsKnotVector<T>& kv_u, const size_t n_t, std::vector<T>& uk, std::vector<index_t>& idx_param_vals); 
		
		/** gets the repeated index of a knot vector (ie if u_i = .. = u_{i+k}, then returns i+k) out of the unique_idx
		* @param[in] kv	knot vector
		* @param[in] unique_idx	ie index without repetitions
		* @param[out] repeated_idx	in case of multiplicities, the last index 
		*/ 
		void unique_to_repeated_idx(const gsKnotVector<T>& kv, const index_t unique_idx, index_t& repeated_idx); 
		
		/** computes the evaluation of basis functions at points given by p_vals
		* @param[in] p_vals	evaluation points
		* @param[in] m_1	number of evaluation points
		* @param[in] kv_u	knot vector
		* @param[in] n_t	number of basis functions in kv_u
		* @param[in] start_idx	indicates at which basisfunction the evaluation should start (starts at 0)
		* @param[in] end_idx	indicates which basisfunction is the last one to evaluate (so if n basisfunction exists, it has to hold end_idx<=(n-1)!) 
		* @param[out] basis_evals	row j of matrix corresponds to basis function j evaluated at all parameter values p_vals, j=0,..,n-1
		*/ 
		void compute_basis_evals(std::vector<T>& p_vals, const index_t m_1, const gsKnotVector<T>& kv_u, const size_t n_t, const size_t start_idx, const size_t end_idx, gsMatrix<T>& basis_evals); 
		
		/** creates a vector containing the removal bounds from function above for each distinct interior point in kv
		* @param[in] n_t	number of control points in this direction
		* @param[in] n_o	number of control points in other direction
		* @param[in] p	degree of spline curve in this direction
		* @param[in] P	3x(n_t*n_o)-matrix containing control points
		* @param[in] kv	knot vector containing all knots in this direction
		* @param[in] dir	indicates which direction is 'this', dir==0 and dir==1 possible
		* @param[in] start_idx	>=(p+1) as repeated index (including multiplicities), indicates for which interior point the error bounds should be computed
		* @param[in] end_idx	<=(n_t-1) as repeated index (including multiplicities), indicates for which interior point the error bounds should be computed
		* @param[out] Br_mat	matrix of size n_ox(kv.uSize()-2), where each unique interior knot corresponds to a column
		*/
		void matrix_removal_bounds(const size_t n_t, const size_t n_o, const index_t p, const gsMatrix<T>& P, const gsKnotVector<T>& kv, const int dir, const size_t start_idx, const size_t end_idx, gsMatrix<T>& Br_mat); 
		
		/**computation of knot removal error bound
		* cf. nurbs book p. 428, Algorithm 9.8
		* @param[in] p	degree of spline curve
		* @param[in] kv	knot vector, size =n+p+1
		* @param[in] Points	3xn-matrix containing control points
		* @param[in] r	index of interior point u which might be removed; r is index with multiplicites, starting with 0
		* @param[in] s	multiplicity of point u which might be removed
		* @param[out] Br	knot removal error bound
		*/ 
		void get_removal_bounds(const index_t p, const gsKnotVector<T>& kv, const gsMatrix<T>& Points, const index_t r, const unsigned s, T& Br); 
		
		/** computes bound for each interior knot by computing t_ij = basis_evals.transpose*Br_mat and setting bound_vec = sum_i t_ij
		* @param[in] Br_mat	matrix of size m_qxinterior_knots, containing all B^j_r
		* @param[in] m_2	number of basis functions in other direction
		* @param[in] interior_knots	number of unique interior knots, ie number of unique knots -2
		* @param[out] bound_vec	vector of size interior_knots, containing bounds for all unique interior knots
		* @param[in] basis_evals_v	matrix with m_q rows, containing the evaluations of basis-functions of other direction
		*/
		void get_bound_surface(const gsMatrix<T>& Br_mat, const size_t m_2, const size_t interior_knots, std::vector<T>& bound_vec, const gsMatrix<T>& basis_evals_v); 
		
		/** computes new control points for a knot-vector where point at position r is removed
		* algorithm is similar to that at nurbs book p185, but without checking if knot can be removed 
		* because this algorithm will only be called in remove_knots_from_surface where there is a check before
		* @param[in] p	degree of curve
		* @param[in] r	index of knot to be removed, where r is the index with multiplicities
		* @param[in] s	multiplicity of knot to be removed
		* @param[in] kv	knot vector defining splines
		* @param[in,out] Points	containing old and new control points
		* @param[out] remove_idx	column remove_idx of Points is labelled as removable (but won't be removed, so size of Points does not change)
		* @warning no knot gets actually removed from kv
		*/
		void remove_knot(const index_t p, const index_t r, const unsigned s, const gsKnotVector<T>& kv, gsMatrix<T>& Points, index_t& remove_idx); 
		
		/** removes knots from surface such that error does not change
		* cf nurbs book, p. 429, Algorithm A9.9 -> modified for surfaces
		* @param[in] kv_u	knot vector from which knots are removed
		* @param[in] p	degree of kv_u
		* @param[in] n_p	number of basis functions in kv_u
		* @param[in] kv_v	knot vector in other direction
		* @param[in] n_q	number of basis functions in kv_v
		* @param[in] Points	matrix of control points (each row is a point) needed to describe current geometry
		* @param[in] p_vals	size = m_p, parameter values in direction of kv_u
		* @param[in] m_p	number of data points in direction of kv_u
		* @param[in] p_vals_v	size = m_q, parameter values in direction of kv_v
		* @param[in] m_q	number of data points in direction of kv_v
		* @param[in,out] error_vals	m_qxm_p-matrix containing error-values for each tensor-parameter-value (uk,vk) where uk in p_vals and vk in p_vals_v
		* @param[in] dir	indicates in which direction knots are removed: dir=0 and dir=1 possible
		* @param[out] n_new	new number of basis-functions and control points, replacing n_p
		* @param[out] kv_new	new knot vector for spline-basis, replacing kv_u
		* @param[out] P_new	new matrix of control points, initially the same size as Points, replacing Points
		* @warning if dir=0: TensorBSpline is (kv_u,kv_v,Points), if dir=1: TensorBSpline is (kv_v, kv_u, Points.transpose())
		*/ 
		void remove_knots_from_surface(const gsKnotVector<T>& kv_u, const index_t p, const size_t n_p, const gsKnotVector<T>& kv_v, const size_t n_q, const gsMatrix<T>& Points, const std::vector<T>& p_vals, const index_t m_1, const std::vector<T>& p_vals_v, const index_t m_2, 
						gsMatrix<T>& error_vals, const int dir, size_t& n_new, gsKnotVector<T>& kv_new, gsMatrix<T>& P_new); 
		
		
		//----------------------------------------------------------------------
		//functions for creating a basis matrix and solving the linear least squares system

		/** checks if knotspan and evaluation points satisfy the condition of Lemma XIV.2 of 	
		* de Boor, practical guide to splines -> gives invertibility
		* @param[in] kv	knot vector
		* @param[in] p	degree of kv
		* @param[in] n	number of basis functions in kv
		* @param[in] u_k	evaluation points
		* @param[out] is_bad	is true if the condition is not satisfied and the matrix thus singular
		*/ 
		void check_knot_span(const gsKnotVector<T>& kv, const index_t p, const size_t n, const std::vector<T>& u_k, bool &is_bad); 
		
		/** constructs basis for a given knot vector and then evaluates this basis at points u_k
		* also uses this basis to construct R which is part of the righthandside needed later in least squares
		* @param[in] kv	knot vector 
		* @param[in] u_k	size = m, containing parameter values
		* @param[in] m	number of points & length of u_k
		* @param[in] n	number of basis functions
		* @param[in] p	degree of knot vector kv
		* @param[out] B	matrix containing b-spline-values after evaluation; if N is matrix from (9.66), then B=N
		* @param[in] n_p	number of data points in other direction
		* @param[in] Q	3x(m*n_p)-matrix of data points, each col is a point
		* @param[out] R	((m-2)*n_px3)-matrix containing righthandside needed in least squares, like (9.67), p.412, in nurbs book
		* @param[in] stp	either 0 or 1 for first resp second direction
		* @param[out] is_bad	indicate if system in fitting-step will be solvable, = true if not solvable
		* @warning	only stp=0 or stp=1 allowed
		* @warning	if is_bad=true, then no values for R or B will be computed
		*/ 
		void make_basis_matrix(const gsKnotVector<T>& kv, const std::vector<T>& u_k, const index_t m, const size_t n, const index_t p, gsMatrix<T>& B, const size_t n_p, const gsMatrix<T>& Q, gsMatrix<T>& R, const int stp, bool& is_bad); 
		
		/** does the fitting step by calling make_basis_matrix and solving the least squares system
		* @param[in] n_p	number of curves in 'other' direction
		* @param[in] n_q	number of desired control points along 'this' direction
		* @param[in] m_q	number of given data points along 'this' direction
		* @param[in] kv_v	knot vector used to construct basis
		* @param[in] q	degree of kv_v
		* @param[in] v_k	vector containing parameter values, size(v_k) = m_q
		* @param[in] P_p	3x(n_p*m_q)-matrix containing approximation points - each col is a point
		* @param[out] P_points	3x(n_p*n_q)-matrix containing control points after fitting
		* @param[in] stp	indicates if first (stp=0) or second (stp=1) fitting step is made
		* @param[out] is_bad	= true if system is not invertible -> in this case, no points are computed
		* @warning	no check if stp==1 or stp==2 only
		*/
		void fitting_step(const size_t n_p, const size_t n_q, const index_t m_q, const gsKnotVector<T>& kv_v, const index_t q, const std::vector<T>& v_k, const gsMatrix<T>& P_p, gsMatrix<T>& P_points, const int stp, bool& is_bad); 
		
		//----------------------------------------------------------------------
		//measuring the distance between approximation and desired points
		/** checks if any of convergence criteria for point projection is satisfied
		* cf. nurbs book, p.233, (6.8)
		* @param[in] R	distance between evaluated surface and desired point 
		* @param[in] k	contains f(u,v) and g(u,v), see p. 232, (6.5)
		* @param[in] Su	first derivatives in u direction
		* @param[in] Sv	first derivatives in v direction
		* @param[in] eps1	epsilon for point coincidence
		* @param[in] eps2	epsilon for zero cosine check
		* @return bool-variable which is true if either point coincidence or zero cosine is satisfied
		*/
		bool check_conv_crit(const gsMatrix<T>& R, const std::vector<T>& k, const std::vector<T>& Su, const std::vector<T>& Sv, const T& eps1, const T& eps2); 
		
		/** checks if point uk is in range [0,1]
		* if not, then uk is set to 0 resp. 1
		* @param[in,out] uk
		*/
		void check_new_vals(T& uk); 
		
		/** checks if parameters do not change significantly, cf nurbs book, p. 234, 4.
		* @param[in] diff_u	i.e. difference between u_i+1 and u_i
		* @param[in] diff_v	i.e. difference between v_i+1 and v_i
		* @param[in] Su	(3x1)-matrix containing first derivatives in u-dir, evaluated
		* @param[in] Sv	(3x1)-matrix containing first derivatives in v-dir, evaluated
		* @param[in] eps1	epsilon for check
		* @return bool-variable which is true if condition holds
		*/
		bool check_change(const T& diff_u, const T& diff_v, const std::vector<T>& Su, const std::vector<T>& Sv, const T& eps1);
		
		/** projects a given point Q onto the surface defined by Tensor, using (u0,v0) as a starting point
		* @param[in] u0	starting point in one direction
		* @param[in] v0	starting point in other direction
		* @param[in] Tensor	B-Spline surface, needed for evaluation of points and first and second derivatives
		* @param[in] Q	point (as one row of a matrix) which gets projected onto a surface
		* @param[out] uk	parameter value in one direction which yields Q
		* @param[out] vk	parameter value in other direction which yields Q
		*/
		void project_to_surface(const T& u0, const T& v0, const gsTensorBSpline<2,T>& Tensor, const gsMatrix<T>& Q, T& uk, T& vk); 
		
		/** calls function project_to_surface for every given point to compute the minimal distance stored in diff_k
		* and also computes new parameter-values
		* @param[in,out] u_k	parameter values in first direction
		* @param[in] m_p	size of u_k
		* @param[in,out] v_k	parameter values in second direction
		* @param[in] m_q	size of v_k
		* @param[in] U	knot vector in first direction
		* @param[in] V	knot_vector in second direction
		* @param[in] P_points	matrix containing control points (each col is a point) to describe current geometry
		* @param[in] Q	3x(m_p*m_q)-matrix of data points whose distance to current geometry is measured
		* @param[out] diff_k	(m_qxm_p)-matrix containing the distance in l2-norm of each point to the geometry
		* @warning no check, if P_points has matching size
		*/ 
		void projection_step(std::vector<T>& u_k, const index_t m_p, std::vector<T>& v_k, const index_t m_q, const gsKnotVector<T>& U, const gsKnotVector<T>& V, const gsMatrix<T>& P_points, const gsMatrix<T>& Q, gsMatrix<T>& diff_k); 
		
		/** searches maximum of all entries in matrix M
		* @param[in] M	matrix
		* @param[out] _row	row-index of maximum
		* @param[out] _col	column-index of maximum
		* @param[out] _value	maximum
		*/ 
		void get_max_in_matrix(const gsMatrix<T>& M, index_t& _row, index_t& _col, T& _value); 
		
		//----------------------------------------------------------------------
		//members 
		
		/// data points
		gsMatrix<T> Q_; 
		
		/// number of data points in 1st direction
		index_t m_p_; 
		
		/// number of data points in 2nd direction
		index_t m_q_; 
		
		///desired degree
		index_t p_; 
		
		///error bound
		T E_; 
		
		///the created Tensor B-Spline
		gsTensorBSpline<2,T> Tensor_; 
				
		/// maximal error 
		T max_error_; 
		
		///mean of errors for every data point in Q_
		T mean_error_;
		
	}; //class gsTensorFitting
	
}//namespace gismo


#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsTensorFitting.hpp)
#endif
