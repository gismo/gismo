/** @file gsTensorFitting.h

    @brief Provides declaration of an iterative fitting method.

    This file is part of the G+Smo library.
    
    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Mayr
*/

#include <gsNurbs/gsTensorBSpline.h>

namespace gismo
{
	
	template<class T>
	gsTensorFitting<T>::gsTensorFitting(const gsMatrix<T>& data_points,
						const index_t& m_1, 
						const index_t& m_2,
						const index_t& degree, 
						const T& error_bound)
	{
		Q_ = data_points; 
		m_p_ = m_1; 
		m_q_ = m_2; 
		p_ = degree; 
		E_ = error_bound; 
	}
		
	template<class T>
	void gsTensorFitting<T>::compute()
	{
		//check, if data fits
		GISMO_ASSERT(Q_.rows()==3 && Q_.cols()==(m_p_*m_q_), "Q does not have correct size"); //a point in Q should correspond to a row
	
		//compute centripetal parameter values u_k and v_k and knot vectors from the parameter values
		std::vector<T> u_k(m_p_, 0.), v_k(m_q_, 0.); 
		gsKnotVector<T> U, V; 
		#pragma omp parallel sections shared(u_k, v_k, U, V)
		{
			#pragma omp section
			{
				param_vals_centripetal(0, u_k);
				U = get_knotvec_from_vec(u_k, m_p_); 
			}
			#pragma omp section
			{
				param_vals_centripetal(1, v_k);
				V = get_knotvec_from_vec(v_k, m_q_); 
			}
		}
	
		//as a start, set the degrees to 1, knot vector as the parameter values and the control points as the given points Q
		// -> this corresponds to interpolation
	
		//set the error e_kl (which corresponds to the points of Q) to zero
		gsMatrix<T> error_vals(m_q_,m_p_); 
		error_vals.setZero(); 

		//set the number of control points (which is now =m_p resp =m_q)
		size_t n_p = m_p_; 
		size_t n_q = m_q_; 
	
		//set the initial control points, P should now have a point in each row
		gsMatrix<T> P = Q_; 

		size_t n_new, n_new2; 
		gsKnotVector<T> kv_new, kv_new2; 
		gsMatrix<T> P_new, P_new2; 
		bool restart = false; 
		bool is_bad = false; 
		std::string message(""); 
		gsMatrix<T> tmp_Points(3,n_q*n_p); 
	
		index_t kv_deg = 1; 
	
		//dir indicates in which direction the fitting-removing takes place
		//for-loop for degrees in both directions
		for (index_t deg = 1; deg<=p_; ++deg)
		{
			//knot removing in p-direction
			remove_knots_from_surface(U, kv_deg, n_p, V, n_q, P, u_k, m_p_, v_k, m_q_, error_vals, 0, n_new, kv_new, P_new); 
					
			//knot removing in q-direction
			remove_knots_from_surface(V, kv_deg, n_q, kv_new, n_new, P_new, v_k, m_q_, u_k, m_p_, error_vals, 1, n_new2, kv_new2, P_new2);  
				
			if ((deg < p_) && (!restart) )
			{
				//check if enough knots have been removed
				if ( (n_new + 1 + (kv_new.uSize() - 2)) > (size_t)m_p_)
				{
					//then number of knots after degree-elevating would be too large -> restart
					restart = true; 
				
					message = message + "\nin direction = 0, deg = " + std::to_string(deg) + ": number of control points does exceed number of given points -> restart \n"; 
				}
				else if ( (n_new2 + 1 + (kv_new2.uSize() - 2)) > (size_t)m_q_)
				{
					//then number of knots after degree-elevating would be too large -> restart
					restart = true; 
				
					message = message + "\nin in direction = 1, deg = " + std::to_string(deg) + ": number of control points does exceed number of given points -> restart \n";
				}
				else
				{
					//fitting-step in both directions is possible!
					//-> update values
					n_p = n_new; 
					U = kv_new; 
				
					//increase all multiplicities in U by 1 -> then degree should be deg+1
					U.degreeElevate(1); 
				
					//update the number of control points resp. basis functions
					n_p += (U.uSize()-2) + 1; 
				
					//fitting to the points in Q
					is_bad = false; 
					fitting_step(m_q_, n_p, m_p_, U, kv_deg+1, u_k, Q_, tmp_Points, 0, is_bad);
				
					//check if basis-matrix is singular (=bad) -> in this case, tmp_Points can't be used
					if (is_bad)
					{
						message = message + "\nin direction = 0, deg = " + std::to_string(deg) + ": basis matrix is singular -> restart \n";
					
						restart = true; 
					}
					else
					{
						//P_points is larger than needed
						P = tmp_Points.block(0,0,3,n_p*m_q_); 
						
						//update values
						n_q = n_new2; 
						V = kv_new2; 
					
						//increase all multiplicites in kv_new2 by 1
						V.degreeElevate(1); 
					
						//update the number of control points resp. basis functions
						n_q += 1 + (V.uSize() - 2); 
					
						//fitting to the points which were computed in first fitting-step
						is_bad = false; 
						fitting_step(n_p, n_q, m_q_, V, kv_deg+1, v_k, P, tmp_Points, 1, is_bad); 
					
						//check if basis-matrix in fitting_step is singular
						if (is_bad)
						{
							message = message + "\nin direction = 1, deg = " + std::to_string(deg) + ": basis matrix is singular -> restart \n";
						
							restart = true; 
						}
						else
						{
							//P_points2 is larger than needed
							P = tmp_Points.block(0,0,3,n_p*n_q); 
						
							++kv_deg; 
						
							//project all points in Q to current surface to get distances and update parameter values
							projection_step(u_k, m_p_, v_k, m_q_,  U, V, P, Q_, error_vals);
													
							//parameter-values and error_vals are updated in projection_step
						}
					}
				}
		
				if (restart)
				{
					//If restart is true, elevate the degree by 1 of the previous surface
					// in all directions 
					message = message + "\n restart \n"; 
				
					//check if degree-elevating would exceed the number of given points in any direction
					if ( ((n_new + 1 + (kv_new.uSize()-2))>(size_t)m_p_) || ( (n_new2 + 1 + (kv_new2.uSize()-2))>(size_t)m_q_) )
					{
						//just return current values
						U = kv_new; 
						V = kv_new2; 
						n_p = n_new; 
						n_q = n_new2; 
						P = P_new2; 
					}
					else
					{	
						gsTensorBSplineBasis<2,T> basis(kv_new, kv_new2); 
						gsTensorBSpline<2,T> Tensor(basis, P_new2.transpose()); 
					
						//elevate the degree
						Tensor.degreeElevate(1,-1); 
					
						U = Tensor.knots(0); 
						V = Tensor.knots(1); 
					
						n_p = U.size() - U.degree() - 1; 
						n_q = V.size() - V.degree() - 1; 
					
						P = Tensor.coefs().transpose(); 
					
						++kv_deg; 
					
						//project all points in Q to current surface to get distances and update parameter values
						projection_step(u_k, m_p_, v_k, m_q_, U, V, P, Q_, error_vals); 
					
						//error-values and parameter values are already updated in projection_step										
														
						restart = false; 
					}
				}
			} 
		}
		
		gsKnotVector<T> kv_u, kv_v; 
		gsMatrix<T> Points; 
	
		if ( (n_p == n_new) && (n_q == n_new2) ) //in case of restart: goes in here
		{
			//no knots have been removed -> return these values
			kv_u = kv_new; 
			kv_v = kv_new2; 
			Points = P_new2; 
		}
		else
		{
			restart = false; 
		
			//Update values
			n_p = n_new; 
			U = kv_new; 
		
			//fitting to the points in Q
			is_bad = false; 
			fitting_step(m_q_, n_p, m_p_, U, kv_deg, u_k, Q_, tmp_Points, 0, is_bad); 
				
			//check if basis-matrix is singular (=bad) -> in this case, P_points can't be used
			if (is_bad)
			{
				message = message + "\nin direction = 0, after loop: basis matrix is singular -> restart \n";
			
				restart = true; 
			}
			else
			{
				//P_points is larger than needed
				P = tmp_Points.block(0,0,3,n_new*m_q_); 
			
				//update values
				n_q = n_new2; 
				V = kv_new2; 
			
				//fitting to the points which were computed in first fitting-step
				is_bad = false; 
				fitting_step(n_p, n_q, m_q_, V, kv_deg, v_k, P, tmp_Points, 1, is_bad); 
				
				//check if basis-matrix in fitting_step is singular
				if (is_bad)
				{
					message = message + "\nin direction = 1, after loop: basis matrix is singular -> restart \n";
				
					restart = true; 
				}
				else
				{
					//P_points2 is larger than needed
					P = tmp_Points.block(0,0,3,n_p*n_q); 
						
					//project all points in Q to current surface to get distances and update parameter values
					projection_step(u_k, m_p_, v_k, m_q_, U, V, P, Q_, error_vals); 
								
					//last knot-removing in p-direction
					remove_knots_from_surface(U, kv_deg, n_p, V, n_q, P, u_k, m_p_, v_k, m_q_, error_vals,0, n_new, kv_new, P_new); 
						
					//last knot-removing in q-direction
					remove_knots_from_surface(V, kv_deg, n_q, kv_new, n_new, P_new, v_k, m_q_, u_k, m_p_, error_vals, 1, n_new2, kv_new2, P_new2); 
						
					//return these values
					kv_u = kv_new;
					kv_v = kv_new2; 
					Points = P_new2; 
					n_p = n_new; 
					n_q = n_new2; 
				}
			}
		
			if (restart)
			{
				message = message + "\n restart at end!\n"; 
				//return the values from last knot-removing step in loop
				kv_u = kv_new; 
				kv_v = kv_new2; 
				Points = P_new2; 
				n_p = n_new; 
				n_q = n_new2; 
			}
		}
		
		projection_step(u_k, m_p_, v_k, m_q_, kv_u, kv_v, Points, Q_, error_vals); 
		T val_;
		index_t i_, j_; 
		get_max_in_matrix(error_vals,i_,j_,val_);

		T mean_diff = 0.; 
		for (index_t i=0; i<m_q_; ++i)
		{
			for (index_t j=0; j<m_p_; ++j)
			{
				mean_diff += error_vals(i,j); 
			}
		}
	
		mean_diff = mean_diff/(m_p_*m_q_); 
		
		mean_error_ = mean_diff; 
		max_error_ = val_; 
	
		gsTensorBSpline<2,T> temp_Tensor(kv_u,kv_v,Points.transpose()); 
		Tensor_ = std::move(temp_Tensor); 

		gsInfo << message << "\n"; 
	
		return; 
	}
	
	
	template<class T>
	T gsTensorFitting<T>::dist_3d(const gsVector<T>& a, const gsVector<T>& b)
	{
		GISMO_ASSERT( (( a.size() == 3) && (a.size() == b.size())), "parameters do not have correct size"); 
		T d = 0.; 
	
		for (int i=0; i<3; i++)
		{
			d += (a(i) - b(i))*(a(i) - b(i));
		}
		d = std::sqrt(d); 
	
		return std::move(d); 
	}

	template<class T>
	void gsTensorFitting<T>::param_vals_centripetal(const int dir, std::vector<T>& u_k)
	{
		//compute centripetal method as on p. 365
		//compute parameter values using chord length; cf nurbs book, p. 364
	
		std::vector<T> dist_tmp(m_p_-1); 
		for (index_t l=0; l<m_q_; ++l)
		{
			T d_p = 0.; 
			for (index_t k=1; k<m_p_; ++k)
			{
				if (dir==0)
				{
					dist_tmp[k-1] = std::sqrt(dist_3d(Q_.col(l*m_p_ + k),Q_.col(l*m_p_ + k-1)));
				}
				else
				{
					dist_tmp[k-1] = std::sqrt(dist_3d(Q_.col(k*m_q_ + l),Q_.col((k-1)*m_q_ + l)));
				}
			
				d_p += dist_tmp[k-1];
			}
		
			T u_l = 0; 
			for (index_t k=1; k<m_p_; ++k)
			{
				u_l += dist_tmp[k-1]/d_p;
				//construct u_k like suggested in nurbs book, p. 376
				u_k[k] += u_l; 
			}			
		}
	
		for (index_t k=0; k<(m_p_-1); ++k)
		{
			u_k[k] /= m_q_;
		}
		u_k[m_p_-1] = 1.; 
	
		return; 
	}

	template<class T>
	gsKnotVector<T> gsTensorFitting<T>::get_knotvec_from_vec(const std::vector<T>& vec, const int m)
	{
		std::vector<T> knot_vec(m+2,0.);
		std::copy(vec.begin() + 1, vec.begin() + m-1, knot_vec.begin() + 2); 
		knot_vec[m] = 1.; 
		knot_vec[m+1] = 1.; 
	
		gsKnotVector<T> U(knot_vec,1); 

		return std::move(U); 
	}

	template<class T>
	void gsTensorFitting<T>::check_knot_span(const gsKnotVector<T>& kv, const index_t p, const size_t n, const std::vector<T>& u_k, bool &is_bad)
	{	
		bool is_good; 
		bool all_is_good = true; 
	
		size_t start_idx = 1; 
		size_t k = start_idx; 
		
		//go through all knot-spans and find a parameter value uk in it
		T left, right; 
	
		for (size_t i=0; i<n; ++i) //basis functions 0 and n are not evaluated	
		{
			left = kv.at(i); 
			right = kv.at(i+p+1); 
			
			is_good = false;
		
			k = start_idx; 
			while ( ( k<(u_k.size()-1)) && (!is_good) )
			{
				if ( ( left + 1e-2< u_k[k] ) && (u_k[k] < right - 1e-2))
				{
					is_good = true; 
				}
				
				++k; 
			}
			start_idx = k; 
		
			all_is_good = all_is_good && is_good; 
		}
	
		is_bad = !all_is_good; 
	
		return; 
	}

	//with check_knot_span from de Boor
	template<class T>
	void gsTensorFitting<T>::make_basis_matrix(const gsKnotVector<T>& kv, const std::vector<T>& u_k, const index_t m, const size_t n, const index_t p, gsMatrix<T>& B, const size_t n_p, const gsMatrix<T>& Q, gsMatrix<T>& R, const int stp, bool& is_bad)
	{
		GISMO_ASSERT( ( stp == 0) || (stp == 1) , "input parameter taking a forbidden value!"); 
		//check if size of u_k is correct
		GISMO_ASSERT(u_k.size() == (size_t)m, "input parameter has wrong size!"); 
		//check if size of Q is consistent
		GISMO_ASSERT((Q.rows()==3) && ((size_t)Q.cols()==(m*n_p)), "input parameter has wrong size!"); 
		//check if R has correct size
		GISMO_ASSERT((R.cols()==3) && ((size_t)R.rows()==(m-2)*n_p), "input parameter has wrong size!"); 
	
		//check knot span
		is_bad = false; 
		check_knot_span(kv, p, n, u_k, is_bad); 
	
		if ( !is_bad)
		{
			//construct matrix consisting of evaluations of basis at points u_k
			gsBSplineBasis<T> bsb(kv);
	
			//have to write u_k as matrix, otherwise no correct evaluation
			gsMatrix<T> uk_mat(1,m-2); 
			for (index_t i=1; i<(m-1); i++)
			{
				uk_mat(0,i-1) = u_k[i]; 
			}
		
			//evaluate bsplinebasis at parameter values, except u[0] and u[end]=u[m-1]
			//B_u = (m-2)x(n-2) matrix where col i-1 corresopnds to basis function i
			gsMatrix<T> tmp;  
			for (size_t i=1; i<(n-1); i++)
			{
				//"evaluate a single basis function i at points uk_mat" (out of gismo-docu)
				tmp = bsb.evalSingle(i,uk_mat); 
				B.col(i-1) = tmp.transpose();
			}
		
			gsMatrix<T> b0 = bsb.evalSingle(0,uk_mat); 
			gsMatrix<T> bn = bsb.evalSingle(n-1,uk_mat);
	
			size_t idx=0; 
			size_t idx0=0; 
			size_t idxm=0; 
	
			//construct R used for solving least squares
			//cf to (9.63), p. 411 in nurbs book
			gsMatrix<T> R_x(m-2,1), R_y(m-2,1), R_z(m-2,1); 

			for (size_t j=0; j<n_p; j++)
			{
				for (size_t i=0; i<(size_t)(m-2); i++)
				{
					switch (stp)
					{
						case 1:
						{
							idx = (i+1)*n_p + j; 
							idx0 = 0*n_p + j; 
							idxm = (m-1)*n_p + j; 
							break;
						}
						case 0:
						{
							idx = j*m + i+1;
							idx0 = j*m + 0; 
							idxm = j*m + m-1;  
							break;
						}
					}
					R_x(i,0) = Q(0,idx) - b0(0,i)*Q(0,idx0) - bn(0,i)*Q(0,idxm); 
					R_y(i,0) = Q(1,idx) - b0(0,i)*Q(1,idx0) - bn(0,i)*Q(1,idxm); 
					R_z(i,0) = Q(2,idx) - b0(0,i)*Q(2,idx0) - bn(0,i)*Q(2,idxm); 
				}
			
				R.block(j*(m-2),0,m-2,1) = R_x; 
				R.block(j*(m-2),1,m-2,1) = R_y; 
				R.block(j*(m-2),2,m-2,1) = R_z; 
			}
		}
		
		return; 
	}

	template<class T>
	void gsTensorFitting<T>::fitting_step(const size_t n_p, const size_t n_q, const index_t m_q, const gsKnotVector<T>& kv_v, const index_t q, const std::vector<T>& v_k, const gsMatrix<T>& P_p, gsMatrix<T>& P_points, const int stp, bool& is_bad)
	{
		GISMO_ASSERT( ( stp == 0) || (stp == 1) ,"input parameter is taking forbidden value!"); 
	
		//check, if P_points (matrix for output-points) has correct size
		GISMO_ASSERT((P_points.rows()==3) && ((size_t)P_points.cols() >= n_q*n_p) , "input parameter has wrong size!"); 
		//check if input-matrix P_p has correct size
		GISMO_ASSERT(((size_t)P_p.cols()==m_q*n_p) && (P_p.rows()==3), "input parameter has wrong size!"); 
		//check if parameter-values have correct size
		GISMO_ASSERT(v_k.size()==(size_t)m_q, "input parameter has wrong size!"); 
	
		P_points.setZero(); 
		
		if (n_q>2)
		{
			gsMatrix<> B_v(m_q-2,n_q-2); //basis matrix
			gsMatrix<> R_q((m_q-2)*n_p,3); //modified R_q
		
			make_basis_matrix(kv_v, v_k, m_q, n_q, q, B_v, n_p, P_p, R_q, stp, is_bad);
		
			if ( !is_bad)
			{
				gsMatrix<> N = B_v.transpose()*B_v; //size=(n_q-2)x(n_q-2)
			
				Eigen::LDLT<Eigen::Matrix<T,-1,-1> > B_solve(N); 
				Eigen::ComputationInfo e_info = B_solve.info();
		
				if (e_info != 0)
				{
					is_bad = true; 
				}
				else
				{
					//only complete fitting-step if basis-matrix was solvable (=is_bad) and 
					// if solving was successful (e_info =0)
		
					gsVector<T> sol_x, sol_y, sol_z; //solutions for linear system of equations
					gsMatrix<T> R((n_q-2)*n_p, 3); //righthandside for LDLT-factorization
				
					gsMatrix<T> B_tr = B_v.transpose(); 
				
					switch (stp)
					{
						case 0:
						{
							size_t idx; 
							
							for (size_t i=0; i<n_p; i++)
							{
								//first point stays the same
								P_points.col(i*n_q) = P_p.col(m_q*i);
						
								//next points are computed by solving least square system
								R = B_tr*(R_q.block(i*(m_q-2),0,m_q-2,3) );
			 
								sol_x = B_solve.solve(R.col(0)); 
								sol_y = B_solve.solve(R.col(1)); 
								sol_z = B_solve.solve(R.col(2)); 
						
								for (size_t k=1; k<(n_q-1); k++)
								{
									idx = i*n_q + k; 
									P_points(0,idx) = sol_x(k-1); 
									P_points(1,idx) = sol_y(k-1); 
									P_points(2,idx) = sol_z(k-1); 
								}
			
								//last point stays the same again
								P_points.col(n_q*(i+1)-1) = P_p.col(m_q*(i+1) - 1);
							}
						
							break; 
						}
						case 1:
						{
							size_t idx; 
							
							for (size_t i=0; i<n_p; i++)
							{
								//first point stays the same
								P_points.col(i) = P_p.col(i);
						
								//next points are computed by solving least square system
								R = B_tr*(R_q.block(i*(m_q-2),0,m_q-2,3) );
			 
								sol_x = B_solve.solve(R.col(0)); 
								sol_y = B_solve.solve(R.col(1)); 
								sol_z = B_solve.solve(R.col(2)); 
						
								for (size_t k=1; k<(n_q-1); k++)
								{
									idx = n_p*k + i; 
								
									P_points(0,idx) = sol_x(k-1); 
									P_points(1,idx) = sol_y(k-1); 
									P_points(2,idx) = sol_z(k-1); 
								}
			
								//last point stays the same again
								P_points.col(n_p*(n_q-1) + i) = P_p.col((m_q-1)*n_p + i);
							}
						
							break; 
						}
					}
				}
			}
		}
		else //n_q=2
		{
			size_t idx1, idx2, idx_a, idx_b; 

			for (size_t j=0; j<n_p; j++)
			{
				if (stp==2)
				{
					idx1 = 0*n_p + j; 
					idx2 = (m_q-1)*n_p + j; 
				
					idx_a = 0*2 + j; 
					idx_b = 1*2 + j; 				
				}
				else 
				{
					idx1 = j*m_q + 0; 
					idx2 = j*m_q + m_q-1; 
				
					idx_a = j*2 + 0; 
					idx_b = j*2 + 1; 
				}
			
				//only two points in this direction, ie only first and last point 
				//-> therefore just take first and last point without changing them
				P_points.col(idx_a) = P_p.col(idx1); 
				P_points.col(idx_b) = P_p.col(idx2); 
			}
		}
		
		return; 
	}

	template<class T>
	void gsTensorFitting<T>::project_to_surface(const T& u0, const T& v0, const gsTensorBSpline<2,T>& Tensor, const gsMatrix<T>& Q, T& uk, T& vk)
	{
		GISMO_ASSERT( (Q.rows() == 3) && (Q.cols() == 1) , "input parameter has wrong size!"); 
	
		//u0 and v0 are parameter values from u_k and v_k -> used here as starting point
		// for the projection onto the surface which is described by Tensor
		
		//start
		T utmp = u0; 
		T vtmp = v0; 
	
		gsMatrix<T> eval_pts(2,1);
		gsMatrix<T> result_d2; 
		gsMatrix<T> result_d; 	
		gsMatrix<T> R(3,1); 
	
		std::vector<T> Su(3), Sv(3);
		std::vector<T> Suu(3), Svv(3), Suv(3);  
		std::vector<T> kappa(2); 
		std::vector<T> delta(2); 

		bool check_zero = false;  //prevent from working with a zero-matrix
		bool check_conv = false; 
		bool check_ch = false; 
		bool check_if_worse = false; 
	
		int counter = 0;  

		//Iterate until any condition is satisfied
		while (check_conv==false && check_ch==false && check_zero==false && check_if_worse==false && counter<1000)
		{
			++counter; 
		
			eval_pts(0,0) = utmp; 
			eval_pts(1,0) = vtmp; 
	
			Tensor.deriv2_into(eval_pts, result_d2); 
			//result_d2 now contains dxxf1(uk,vk),dyyf1(uk,vk),dxyf1(uk,vk),dxxf2(uk,vk),..)^T
		
			Tensor.deriv_into(eval_pts, result_d); 
			//result_d: (duf1(uk,vk),dvf1(uk,vk),duf2(uk,vk),dvf2(uk,vk),duf3(uk,vk),dvf3(uk,vk))^T
		
			Tensor.eval_into(eval_pts, R); 
			R -= Q; //Q is supposed to be also a 1x3-matrix
	
			//result_d2 and result_d are not in the correct order for this projection, 
			//therefore make new matrices Suu,Svv,etc as in nurbs book, p.233
			Su[0] = result_d(0,0); 
			Su[1] = result_d(2,0); 
			Su[2] = result_d(4,0); 
		
			Sv[0] = result_d(1,0); 
			Sv[1] = result_d(3,0); 
			Sv[2] = result_d(5,0); 
			
			//kappa_i as on p. 233
			//kappa(0) = <-R,Su>, kappa(1) = <-R, Sv>
			kappa[0] = -R(0,0)*Su[0] - R(1,0)*Su[1] - R(2,0)*Su[2]; 
			kappa[1] = -R(0,0)*Sv[0] - R(1,0)*Sv[1] - R(2,0)*Sv[2]; 
		
			//check if convergence criteria are already satisfied
			//cf nurbs book, p. 233, (6.8)
			T eps1 = 0.00000001; 
			T eps2 = 0.00001; 
		
			//check for point-coincidence and zero-cosine
			check_conv = check_conv_crit(R, kappa, Su, Sv, eps1, eps2); 
		
			if (check_conv) //if point-coindidence or zero-cosine is satisfied
			{
				uk = utmp; 
				vk = vtmp; 
			}
			else //check_conv=false -> compute new points
			{
				//result_d2 and result_d are not in the correct order for this projection, 
				//therefore make new matrices Suu,Svv,etc as in nurbs book, p.233
				Suu[0] = result_d2(0,0); 
				Suu[1] = result_d2(3,0); 
				Suu[2] = result_d2(6,0); 
			
				Svv[0] = result_d2(1,0); 
				Svv[1] = result_d2(4,0); 
				Svv[2] = result_d2(7,0); 
			
				Suv[0] = result_d2(2,0); 
				Suv[1] = result_d2(5,0); 
				Suv[2] = result_d2(8,0); 
			
				//J_i as on p.233
				//J is symmetric
				T j11 = Su[0]*Su[0] + Su[1]*Su[1] + Su[2]*Su[2] + R(0,0)*Suu[0] + R(1,0)*Suu[1] + R(2,0)*Suu[2];
				T j12 = Su[0]*Sv[0] + Su[1]*Sv[1] + Su[2]*Sv[2] + R(0,0)*Suv[0] + R(1,0)*Suv[1] + R(2,0)*Suv[2]; 
				T j22 = Sv[0]*Sv[0] + Sv[1]*Sv[1] + Sv[2]*Sv[2] + R(0,0)*Svv[0] + R(1,0)*Svv[1] + R(2,0)*Svv[2];  
			
				T det_J = (j11*j22 - j12*j12); 

				if ( std::abs(det_J) < 1e-10 )
				{
					check_zero = true; 
					uk = utmp; 
					vk = vtmp; 
				}
				else
				{
					//solve system J*delta = kappa
					delta[0] = (j22*kappa[0] - j12*kappa[1])/det_J; 
					delta[1] = (-j12*kappa[0] + j11*kappa[1])/det_J; 
					uk = utmp + delta[0]; 
					check_new_vals(uk); 
				
					vk = vtmp + delta[1]; 
					check_new_vals(vk); 
				
					//If check_ch=true, then point is found
					//condition on change of parameters is satisfied -> exiting iteration
					check_ch = check_change(uk-utmp, vk-vtmp, Su, Sv, eps1); 
						
					//new check
					gsMatrix<T> Rtmp; 
					eval_pts(0,0) = uk; 
					eval_pts(1,0) = vk; 
					Tensor.eval_into(eval_pts, Rtmp); 
					Rtmp -= Q;
					T Rtmp_norm = std::sqrt(Rtmp(0,0)*Rtmp(0,0) + Rtmp(1,0)*Rtmp(1,0) + Rtmp(2,0)*Rtmp(2,0)); 
					T R_norm = std::sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0) + R(2,0)*R(2,0)); 
					
					if (Rtmp_norm > R_norm)
					{
						check_if_worse = true; 
					
						//return old values if the new values give a worse result
						uk = utmp; 
						vk = vtmp; 
					}
					else
					{
						if (!check_ch)
						{
							utmp = uk; 
							vtmp = vk;
						}
					}	
				}
			}
		
		} 
	
		return; 
	}

	
	template<class T>
	bool gsTensorFitting<T>::check_conv_crit(const gsMatrix<T>& R, const std::vector<T>& k, const std::vector<T>& Su, const std::vector<T>& Sv, const T& eps1, const T& eps2)
	{
		GISMO_ASSERT( (R.rows()==3) && (R.cols()==1) , "input parameter has wrong size!"); 
		GISMO_ASSERT( (Su.size() == 3) && (Sv.size() == 3) , "input parameter has wrong size!"); 
	
		//R = S(ui,vi) - P
		//Su = S_u(ui,vi)
		//Sv = S_v(ui,vi)
		//k = { <Su,R>, <Sv,R> }
	
		bool check_yn; 
	
		//point coincidence
		T R_tmp = std::sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0) + R(2,0)*R(2,0)); 
		bool check1 = (R_tmp<eps1); 
	
		//zero cosine
		T Su_tmp = std::sqrt(Su[0]*Su[0] + Su[1]*Su[1] + Su[2]*Su[2]);
		T Sv_tmp = std::sqrt(Sv[0]*Sv[0] + Sv[1]*Sv[1] + Sv[2]*Sv[2]); 
		
		bool check2u = ( std::abs(k[0]) <= (eps2*Su_tmp*R_tmp) ); 
		bool check2v = ( std::abs(k[1]) <= (eps2*Sv_tmp*R_tmp) ); 

		//iteration halts if any condition is satisfied
		check_yn = check1 || (check2u && check2v) ; 
	
		return check_yn; 
	}

	template<class T>
	void gsTensorFitting<T>::check_new_vals(T& uk)
	{
		if (uk < 0)
		{
			uk = 0; 
		}
		else if (uk > 1)
		{
			uk = 1; 
		}

		return; 
	}

	template<class T>
	bool gsTensorFitting<T>::check_change(const T& diff_u, const T& diff_v, const std::vector<T>& Su, const std::vector<T>& Sv, const T& eps1)
	{
		GISMO_ASSERT( (Su.size() == 3) && (Sv.size() == 3) , "input parameter has wrong size!"); 
	
		std::vector<T> Su_tmp(3), Sv_tmp(3); 
		for (int i=0; i<3; ++i)
		{
			Su_tmp[i] = Su[i]*diff_u; 
			Sv_tmp[i] = Sv[i]*diff_v; 
		}
		
		T abs_S = std::sqrt( (Su_tmp[0] + Sv_tmp[0])*(Su_tmp[0] + Sv_tmp[0]) 
							+ (Su_tmp[1] + Sv_tmp[1])*(Su_tmp[1] + Sv_tmp[1])
							+ (Su_tmp[2] + Sv_tmp[2])*(Su_tmp[2] + Sv_tmp[2])
							); 

		return (abs_S<eps1); 
	}

	template<class T>
	void gsTensorFitting<T>::projection_step(std::vector<T>& u_k, const index_t m_p, std::vector<T>& v_k, const index_t m_q, const gsKnotVector<T>& U, const gsKnotVector<T>& V, const gsMatrix<T>& P_points, const gsMatrix<T>& Q, gsMatrix<T>& diff_k)
	{	
		//check if sizes are consistent
		GISMO_ASSERT((diff_k.rows()==m_q) && (diff_k.cols()==m_p), "input parameter has wrong size!"); 
		GISMO_ASSERT((Q.cols()==m_q*m_p) && (Q.rows()==3), "input parameter has wrong size!"); 
	
		//set uk_new and vk_new to 0
		std::vector<T> uk_new(m_p,0.); 
		std::vector<T> vk_new(m_q,0.); 
	
		//construct values for projection
		gsTensorBSpline<2,T> Tensor(U,V,P_points.transpose()); 
	
		std::vector<T> uk_proj(m_p*m_q), vk_proj(m_p*m_q); 

		#pragma omp parallel for shared(uk_proj, vk_proj, Tensor, u_k, v_k) collapse(2) schedule(static)
		for (index_t l=0; l<m_q; l++)
		{
			for (index_t k=0; k<m_p; k++)
			{
				project_to_surface(u_k[k], v_k[l], Tensor, Q.col(l*m_p + k), uk_proj[k*m_q+l], vk_proj[l*m_p+k]);
			}
		}
	
		//recompute parameter values
		u_k[0] = 0.; 
		for (index_t i=1; i<m_p-1; ++i)
		{
			u_k[i] = 0.; 
			for (index_t j=0; j<m_q; ++j)
			{
				u_k[i] += uk_proj[i*m_q+j]; 
			}
			u_k[i] /= m_q; 
		}
		u_k[m_p-1] = 1.; 
	
		v_k[0] = 0.; 
		for (index_t j=1; j<m_q-1; ++j)
		{
			v_k[j] = 0.; 
			for (index_t i=0; i<m_p; ++i)
			{
				v_k[j] += vk_proj[j*m_p+i]; 
			}
			v_k[j] /= m_p; 
		}
		v_k[m_q-1] = 1.; 

		gsMatrix<T> R(3,1); 
		gsMatrix<T> eval_pts(2,1);
		
		#pragma omp parallel for firstprivate(eval_pts, R) shared(Tensor, u_k, v_k, Q) collapse(2) schedule(static)
		for (index_t k=0; k<m_p; ++k)
		{
			for (index_t l=0; l<m_q; ++l)
			{
				eval_pts(0,0) = u_k[k];
				eval_pts(1,0) = v_k[l];
			
				//evaluate points and compute distance to original point
				Tensor.eval_into(eval_pts,R);  
				R -= Q.col(l*m_p + k);
				diff_k(l,k) = std::sqrt(R(0,0)*R(0,0) + R(1,0)*R(1,0) + R(2,0)*R(2,0)); 	
			}
		}
			
		return; 
	}

	//here r is index with multiplicities
	template<class T>
	void gsTensorFitting<T>::get_removal_bounds(const index_t p, const gsKnotVector<T>& kv, const gsMatrix<T>& Points, const index_t r, const unsigned s, T& Br)
	{
		//only remove interior points -> check if r is interior point
		GISMO_ASSERT((r>p), "input parameter has wrong size!"); 

		index_t i = r - p; 
		index_t j = r - s; 
		index_t p_1 = p + 1; 
		T u = kv.at(r); 

		T alpha_i, alpha_j; 
		gsMatrix<T> P_tmp(3,p-s+3); 
	
		P_tmp.col(0) = Points.col(i-1); 
		P_tmp.col(p-s+2) = Points.col(j+1); 

		index_t idx_i = 1; 
		index_t idx_j = p-s+1; 
	
		while (j-i>0)
		{
			alpha_i = (u - kv.at(i))/(kv.at(i+p_1) - kv.at(i)); 
			alpha_j = (u - kv.at(j))/(kv.at(j+p_1) - kv.at(j)); 

			P_tmp.col(idx_i) = (Points.col(i) - (1-alpha_i)*P_tmp.col(idx_i-1))/alpha_i; 
			P_tmp.col(idx_j) = (Points.col(j) - alpha_j*P_tmp.col(idx_j+1))/(1-alpha_j);
		
			++i; 
			--j; 
			++idx_i; 
			--idx_j; 
		}
		
		if (j-i<0)
		{
			//equation (9.82)
			Br = dist_3d(P_tmp.col(idx_i-1),P_tmp.col(idx_j+1)); 
		}
		else //j-i==0
		{
			//equation (9.80)
			alpha_i = (u - kv.at(i))/(kv.at(i+p_1) - kv.at(i)); 
			Br = dist_3d(Points.col(i), alpha_i*P_tmp.col(idx_i+1) + (1-alpha_i)*P_tmp.col(idx_i-1)); 		
		}
		return; 
	}

	template<class T>
	void gsTensorFitting<T>::matrix_removal_bounds(const size_t n_t, const size_t n_o, const index_t p, const gsMatrix<T>& P, const gsKnotVector<T>& kv, const int dir, const size_t start_idx, const size_t end_idx, gsMatrix<T>& Br_mat)
	{	
		//check if P has correct size
		GISMO_ASSERT(((size_t)P.cols()==(n_t*n_o)) && (P.rows()==3), "input parameter has wrong size!"); 
		//check, if knot-vector and n_t and p fit together
		GISMO_ASSERT(kv.size() == (size_t)(n_t + p + 1), "input parameter has wrong size!"); 
		//check, if start_idx and end_idx are in idx-range of knot vector
		GISMO_ASSERT((start_idx>=(size_t)(p+1)) && (end_idx<=n_t), "input parameter has wrong size!"); 
		GISMO_ASSERT( (dir == 0) || (dir == 1) ,"input parameter is taking a forbidden value"); 
	
		T value = kv.at(start_idx); 
		typename gsKnotVector<T>::uiterator uit_start = kv.uFind(value); 
	
		value = kv.at(end_idx+1); 
		typename gsKnotVector<T>::uiterator uit_end = kv.uFind(value);
	
		if ( end_idx == n_t)
		{
			uit_end = kv.uend() - 1;  
		}
		else
		{
			++uit_end; 
		}
	
		gsMatrix<T> P_tmp(3,n_t); 
		index_t gen_idx_start; 
		if (start_idx==(size_t)(p+1))
		{
			gen_idx_start = p; 
		}
		else //start_idx>p+1
		{
			index_t uidx = (uit_start - kv.ubegin()) - 1; //unique index before unique start-idx
			unique_to_repeated_idx(kv, uidx, gen_idx_start);
		}

		//stepping through all curves 

		switch (dir)
		{
			case 0: 
			{
				for (size_t l=0; l<n_o; l++)
				{
					P_tmp = P.block(0,l*n_t,3,n_t); 
		
					index_t gen_idx = gen_idx_start; 
						
					//iterating through unique interior knots between unique_start and unique_end
					for (typename gsKnotVector<T>::uiterator uit = uit_start; uit != uit_end; ++uit)
					{
						index_t k = uit - kv.ubegin(); 
						unsigned mult = kv.u_multiplicityIndex(k); 
						gen_idx += mult; 
						get_removal_bounds(p, kv, P_tmp, gen_idx, mult, Br_mat(l,k-1));
					}			
				}
				break; 
			}
			case 1: 
			{
				for (size_t l=0; l<n_o; l++)
				{
			
					for (size_t k=0; k<n_t; k++)
					{
						P_tmp.col(k) = P.col(k*n_o + l); 
					}
				
					index_t gen_idx = gen_idx_start; 
						
					//iterating through unique interior knots between unique_start and unique_end
					for (typename gsKnotVector<T>::uiterator uit = uit_start; uit != uit_end; ++uit)
					{
						index_t k = uit - kv.ubegin(); 
						unsigned mult = kv.u_multiplicityIndex(k); 
						gen_idx += mult; 
						get_removal_bounds(p, kv, P_tmp, gen_idx, mult, Br_mat(l,k-1));
					}		
				}
				break; 
			}
		}
	
		return; 
	}

	//here idx_tmp is index of basis function, so only repeated_idx from knot vector can be used for that
	template<class T>
	void gsTensorFitting<T>::get_param_range(const size_t idx_tmp_begin, const size_t idx_tmp_end, const gsKnotVector<T>& kv_u, const size_t n_t, std::vector<T>& uk, std::vector<index_t>& idx_param_vals)
	{
		GISMO_ASSERT((idx_tmp_begin>=0) && (idx_tmp_end<=n_t), "input parameter has wrong size!"); 
		GISMO_ASSERT(idx_param_vals.size() == 2*n_t, "input parameter has wrong size!"); 
	
		gsBSplineBasis<T> basis(kv_u); 
		gsMatrix<T> M_tmp; 
		typename std::vector<T>::iterator it;
	
		for (size_t i=idx_tmp_begin; i < idx_tmp_end; ++i)
		{
			M_tmp = basis.support(i);
	
			it = std::lower_bound(uk.begin(), uk.end(), M_tmp(0)); 
			idx_param_vals[2*i] = it - uk.begin(); 
		
			it = std::upper_bound(uk.begin(), uk.end(), M_tmp(1)); 
			idx_param_vals[2*i+1] = it - uk.begin(); 
		}
	
		return; 
	}

	template<class T>
	void gsTensorFitting<T>::get_bound_surface(const gsMatrix<T>& Br_mat, const size_t m_2, const size_t interior_knots, std::vector<T>& bound_vec, const gsMatrix<T>& basis_evals_v)
	{	
		GISMO_ASSERT( ((size_t)Br_mat.rows() == m_2) && ((size_t)Br_mat.cols() == interior_knots), "input parameter has wrong size!"); 
		GISMO_ASSERT((size_t)bound_vec.size()==interior_knots, "input parameter has wrong size!"); 
		GISMO_ASSERT((size_t)basis_evals_v.rows()==m_2, "input parameter has wrong size!"); 
	
		//m_2 = number of evaluation points of other direction
		//interior_knots = number of interior knots in this direction
	
		gsMatrix<T> tmp = basis_evals_v.transpose()*Br_mat; 
				
		//write sum of tmp into bound_vec
		for (size_t i=0; i<interior_knots; ++i)
		{
			bound_vec[i] = tmp.col(i).sum(); 
		}
	
		return; 
	}

	template<class T>
	void gsTensorFitting<T>::remove_knot(const index_t p, const index_t r, const unsigned s, const gsKnotVector<T>& kv, gsMatrix<T>& Points, index_t& remove_idx)
	{
		//check if only interior knots will be removed
		GISMO_ASSERT((r>p), "input parameter has wrong size!");  
	
		GISMO_ASSERT( (Points.rows() == 3) , "input parameter has wrong size!"); 
		
		index_t i = r - p; 
		index_t j = r - s; 
		index_t p_1 = p + 1; 
		T u = kv.at(r); 
		index_t off = i - 1; 

		T alpha_i, alpha_j; 
		gsMatrix<T> P_tmp(3,p-s+3); 
		
		//initial temporarily points
		P_tmp.col(0) = Points.col(off); 
		P_tmp.col(j+1-off) = Points.col(j+1); 
	
		index_t idx_i = 1; 
		index_t idx_j = p-s+1; 
	
		while (j-i>0)
		{
			alpha_i = (u - kv.at(i))/(kv.at(i+p_1) - kv.at(i)); 
			alpha_j = (u - kv.at(j))/(kv.at(j+p_1) - kv.at(j)); 

			P_tmp.col(idx_i) = (Points.col(i) - (1-alpha_i)*P_tmp.col(idx_i-1))/alpha_i; 
			P_tmp.col(idx_j) = (Points.col(j) - alpha_j*P_tmp.col(idx_j+1))/(1-alpha_j);
		
			++i; 
			--j; 
			++idx_i; 
			--idx_j; 
		}
	
		//no check if knot is removable
	
		i = r - p; 
		j = r - s; 		
		while (j-i>0)
		{
			Points.col(i) = P_tmp.col(i-off); 
			Points.col(j) = P_tmp.col(j-off); 

			++i; 
			--j; 
		}
	 
		//remove row i
		remove_idx = i; 
		
		return; 	
	}
	
	template<class T>
	void gsTensorFitting<T>::unique_to_repeated_idx(const gsKnotVector<T>& kv, const index_t unique_idx, index_t& repeated_idx)
	{
		T value = kv.uValue(unique_idx); 
	
		typename gsKnotVector<T>::iterator it = kv.iFind(value); 
		repeated_idx = it - kv.begin(); 

		return; 
	}
	
	template<class T>
	void gsTensorFitting<T>::compute_basis_evals(std::vector<T>& p_vals, const index_t m_1, const gsKnotVector<T>& kv_u, const size_t n_t, const size_t start_idx, const size_t end_idx, gsMatrix<T>& basis_evals)
	{
		gsBSplineBasis<T> basis(kv_u); 
		GISMO_ASSERT(((size_t)basis_evals.rows()==n_t) && (basis_evals.cols()==m_1), "input parameter has wrong size!"); 
	
		//make a matrix out of parameter-values, because a vector is not accepted
		//it suffices to compute everything for the first row, because the values don't change
		gsMatrix<T> p_vals_mat(1,m_1);
		p_vals_mat.row(0) = gsAsVector<T>(p_vals); 

		gsMatrix<T> tmp(1,m_1); 
		for (size_t j=start_idx; j<=end_idx; j++)
		{
			tmp = basis.evalSingle(j,p_vals_mat);
			basis_evals.row(j) = std::move(tmp); //row j of matrix corresponds to basis function j evaluated at all parameter values
		}

		return; 
	}

	template<class T>
	void gsTensorFitting<T>::remove_knots_from_surface(const gsKnotVector<T>& kv_u, const index_t p, const size_t n_p, const gsKnotVector<T>& kv_v, const size_t n_q, const gsMatrix<T>& Points, const std::vector<T>& p_vals, const index_t m_1, const std::vector<T>& p_vals_v, const index_t m_2, 
						gsMatrix<T>& error_vals, const int dir, size_t& n_new, gsKnotVector<T>& kv_new, gsMatrix<T>& P_new)
	{
			
		GISMO_ASSERT( (dir == 0) || (dir == 1) ,"input parameter is taking a forbidden value!"); 
		GISMO_ASSERT( (error_vals.rows() == m_2) && (error_vals.cols() == m_1) , "input parameter has wrong size!"); 
	
		if ( kv_u.uSize() <= 2)	
		{
			kv_new = kv_u; 
			n_new = n_p; 
			P_new = Points; 
		}
		else
		{
	
		//knot-vectors for tensor grid
		gsKnotVector<> kv1, kv2; 
	
		if (dir == 0)
		{
			kv1 = kv_u; 
			kv2 = kv_v; 
		}
		else
		{
			kv1 = kv_v; 
			kv2 = kv_u; 
		}
		
		//p = degree of kv_u, 
		//n_p = number of control points of kv_u, 
		//n_q = number of control points of kv_v
	
		size_t n_t = n_p; 
		size_t n_o = n_q; 
	
		//m_1 = number of data points in direction of kv_u
		//m_2 = number of data points in direction of kv_v
	
		//check if control points have correct size
		GISMO_ASSERT(((size_t)Points.cols()==(n_t*n_o)) && (Points.rows()==3), "input parameter has wrong size!"); 
	
		std::vector<T> uk(p_vals); 
		std::vector<T> vk(p_vals_v); 
	
		//compute for each basis function (in both directions) the range of parameter indices
		//only need to compute these ranges for the first curve, as it does not change
		//compute evaluations of basis functions in both directions
	
		std::vector<index_t> idx_param_vals(2*n_t), idx_param_vals_v(2*n_o);
		gsMatrix<T> basis_evals(n_t,m_1), basis_evals_v(n_o,m_2); 
		#pragma omp parallel sections shared(idx_param_vals, idx_param_vals_v, basis_evals, basis_evals_v, uk, vk)
		{
			#pragma omp section
			{
				get_param_range(0, n_t, kv_u, n_t, uk, idx_param_vals); 
				compute_basis_evals(uk, m_1, kv_u, n_t, 0, n_t-1, basis_evals);
			}
			#pragma omp section
			{
				get_param_range(0, n_o, kv_v, n_o, vk, idx_param_vals_v);
				compute_basis_evals(vk, m_2, kv_v, n_o, 0, n_o-1, basis_evals_v);
			}
		}
	
		size_t num_uinterior = kv_u.uSize() - 2; 

		//compute error bounds for all curves in direction 'dir'
		gsMatrix<> Br_mat(n_o,num_uinterior);  
		matrix_removal_bounds(n_t, n_o, p, Points, kv_u, dir, p+1, n_t, Br_mat); 
	
		//compute error bounds for interior points as sum over all curves 
		std::vector<T> Br_vec(num_uinterior);  
		get_bound_surface(Br_mat, n_o, num_uinterior, Br_vec, basis_evals_v);

		bool keep_iterating = true; 
		bool knot_is_removable = true; 
		typename std::vector<T>::iterator it; 
		gsMatrix<T> new_error(m_2,m_1); 
		gsMatrix<T> temp_error = error_vals;
		T Br_value;
		gsMatrix<T> P_tmp(3,n_t); 
		index_t repeated_idx; 
		std::vector<index_t> remove_indices(n_o); 
	
		//set new values to initial values 
		kv_new = kv_u; 
		P_new = Points; 
		n_new = n_t; 
	
		const T big_number = 1000000000; 

		while (keep_iterating)
		{
			knot_is_removable = true; 

			//find knot with smallest B_r and set unique index of knot and multiplicity of this knot
			it = std::min_element(Br_vec.begin(), Br_vec.end()); 
			index_t unique_idx = (it - Br_vec.begin() ) + 1; 
			Br_value = *it; 
			unsigned mult = kv_new.u_multiplicityIndex(unique_idx); 
		
			unique_to_repeated_idx(kv_new,unique_idx,repeated_idx); 
		
			//check if minimum is really minimum
			if (std::abs(Br_value-big_number)<1e-14)
			{
				keep_iterating = false; //need to step out immediately
			}
			else
			{
				//compute error bounds as in (9.87) and (9.89), p434 nurbs book
				index_t tmp_int = (p+mult)%2; 
				index_t k, index_case; 
				T alpha; 
				switch (tmp_int)
				{
					case 0:
					{
						k = (p+mult)/2; 
				
						alpha = 0.; //alpha basically not needed here
				
						index_case = repeated_idx - k;
				
						break; 
					}
					case 1:
					{
						k = (p+mult+1)/2; 
					
						index_case = repeated_idx - k + 1;
					
						alpha = (kv_new.at(repeated_idx)-kv_new.at(index_case))/(kv_new.at(index_case+p+1)-kv_new.at(index_case));
					
						break; 
					}
					default:
					{
						GISMO_ASSERT( ( tmp_int == 0) || (tmp_int == 1) , "input parameter is taking a forbidden value!"); 
					}
				}//end switch
								
				//computing error bounds---------------------------------------
				//go through parameter values, that lie in support of basis function index_case
				//need here repeated index, because this refers to the basis functions
			
				//set new_error to 0
				new_error.setZero(); 
						
				switch (dir)
				{
					case 0: 
					{
						index_t dist = idx_param_vals[2*index_case+1] - idx_param_vals[2*index_case]; 
						std::vector<bool> is_removable_vec(dist,true);
						index_t start_idx = idx_param_vals[2*index_case];  
													
						for (index_t i=idx_param_vals[2*index_case]; i<idx_param_vals[2*index_case+1]; ++i)
						{		 
							//go through all basis functions in other direction
							for (size_t j=0; j<n_o; j++)
							{ 
								//and now go through all parameter values in other direction which lie in support of these basis functions
								for (index_t l=idx_param_vals_v[2*j]; l<idx_param_vals_v[2*j+1]; ++l)
								{
									//Need here unique_idx, because Br_mat only exists for unique interior indices
									new_error(l,i) += basis_evals_v(j,l)*Br_mat(j,unique_idx-1); 
								}
							}
				
							//here repeated_idx is used, because basis_evals refers to the index of basis functions
							temp_error.col(i) = error_vals.col(i) + (1-alpha)*basis_evals(index_case,i)*new_error.col(i); 
		
							//checks if knot is removable, i.e. if temp_error(i) <=E_bound
							//If for all i the knot is removable, then true should be the result
							for (index_t j=0; j<m_2; ++j)
							{
								is_removable_vec[i-start_idx] = is_removable_vec[i-start_idx] && (temp_error(j,i)<=E_); 
							}
						}
					
						for (index_t i=idx_param_vals[2*index_case]; i<idx_param_vals[2*index_case+1]; ++i)
						{
							knot_is_removable = knot_is_removable && is_removable_vec[i-start_idx]; 
						}
					
						break; 
					}
					case 1: 
					{
						index_t dist = idx_param_vals[2*index_case+1] - idx_param_vals[2*index_case]; 
						std::vector<bool> is_removable_vec(dist,true);
						index_t start_idx = idx_param_vals[2*index_case];  
					
						for (index_t i=idx_param_vals[2*index_case]; i<idx_param_vals[2*index_case+1]; i++)
						{	 
							//go through all basis functions in other direction
							for (size_t j=0; j<n_o; j++)
							{ 
								//and now go through all parameter values in other direction which lie in support of these basis functions
								for (index_t l=idx_param_vals_v[2*j]; l<idx_param_vals_v[2*j+1]; l++)
								{
									//Need here unique_idx, because Br_mat only exists for unique interior indices
									new_error(i,l) += basis_evals_v(j,l)*Br_mat(j,unique_idx-1); 
								}
							}
				
							//here repeated_idx is used, because basis_evals refers to the index of basis functions
							temp_error.row(i) = error_vals.row(i) + (1-alpha)*basis_evals(index_case,i)*new_error.row(i); 
						
							//checks if knot is removable, i.e. if temp_error(i) <=E_bound
							//If for all i the knot is removable, then true should be the result
							for (index_t j=0; j<m_2; j++)
							{
								is_removable_vec[i-start_idx] = is_removable_vec[i-start_idx] && (temp_error(i,j)<=E_); 
							}
						}
					
						for (index_t i=idx_param_vals[2*index_case]; i<idx_param_vals[2*index_case+1]; ++i)
						{
							knot_is_removable = knot_is_removable && is_removable_vec[i-start_idx]; 
						}
						
						break; 
					}
				}
			
				//Now check if knot can be removed
				if (knot_is_removable)
				{
					switch (dir)
					{
						case 0: 
						{	
							//update error-values
							index_t dist_param_vals = idx_param_vals[2*index_case+1] - idx_param_vals[2*index_case]; 
							error_vals.block(0,idx_param_vals[2*index_case],m_2,dist_param_vals) = temp_error.block(0,idx_param_vals[2*index_case],m_2,dist_param_vals);  
						
							//compute new control points for after knot-removal
							for (size_t i=0; i<n_o; ++i)
							{
								P_tmp = P_new.block(0,i*n_new,3,n_new);
						
								//remove knot -> warning, a row in P_tmp gets removed during function
								remove_knot(p,repeated_idx,mult,kv_new, P_tmp, remove_indices[i]);
					
								//now remove also row in P_new and store back the entries of P_tmp into P_new
								P_new.block(0,i*n_new,3,n_new) = std::move(P_tmp); 
							}
						
							//need to update control points
							gsMatrix<T> P(3, (n_new-1)*n_o);
							for (size_t i=0; i<n_o; ++i)
							{
								size_t ni = n_new - remove_indices[i] - 1; 
								P.block( 0, i*(n_new-1), 3, remove_indices[i]) = P_new.block(0, i*n_new, 3, remove_indices[i]); 
								P.block( 0, i*(n_new-1) + remove_indices[i], 3, ni) = P_new.block(0, i*n_new + remove_indices[i]+1, 3, ni); 
							}
							P_new = std::move(P); 
						
							break; 
						}
						case 1: 
						{
							//update error-values
							index_t dist_param_vals = idx_param_vals[2*index_case+1] - idx_param_vals[2*index_case]; 
							error_vals.block(idx_param_vals[2*index_case],0,dist_param_vals, m_2) = temp_error.block(idx_param_vals[2*index_case], 0, dist_param_vals, m_2);  
						
							//to get correct size for P_tmp
							P_tmp = P_new.block(0,0,3,n_new); 
						
							//compute new control points for after knot-removal							
							for (size_t i=0; i<n_o; ++i)
							{
								for (size_t k=0; k<n_new; k++)
								{
									P_tmp.col(k) = P_new.col(k*n_o + i); 
								}
							
								//remove knot -> warning, a row in P_tmp gets removed during function
								remove_knot(p,repeated_idx,mult,kv_new, P_tmp, remove_indices[i]);
					
								//now remove also row in P_new and store back the entries of P_tmp into P_new
								for (size_t k=0; k<n_new; k++)
								{
									P_new.col(k*n_o + i) = P_tmp.col(k); 
								}
							}
						
							//need to update control points
							gsMatrix<> P(3, (n_new-1)*n_o);
							for (size_t i=0; i<n_o; ++i)
							{								
								for (size_t j=0; j<(size_t)remove_indices[i]; ++j)
								{
									P.col(j*n_o + i) = P_new.col(j*n_o + i); 
								}	
								for (size_t j=(size_t)remove_indices[i]; j<(n_new-1); ++j)
								{
									P.col(j*n_o + i) = P_new.col((j+1)*n_o + i); 
								}
							}
							P_new = std::move(P); 
						
							break; 
						}
					}
				
					//the corresponding knot still has to be removed
					T knot = kv_new.at(repeated_idx); 
					kv_new.remove(knot); 
					--n_new; 

					//if no more internal knots are left, then break
					if (kv_new.uSize()==2) //only the two outer knots are left
					{
						keep_iterating = false;
					}
					else
					{
						//erase values because now there is one basis function less
						idx_param_vals.erase( idx_param_vals.begin() + 2*(repeated_idx-p), idx_param_vals.begin() + 2*(repeated_idx-p+1)); 
					
						//according to (9.84), update only parameters (with new numbers) repeated_idx-p-1, ..., repeated_idx-mult
						get_param_range(repeated_idx-p-1, repeated_idx-mult+1, kv_new, n_new, uk, idx_param_vals); 
					
						//also need to recompute the basis evaluations, because basis-functions have changed
						gsMatrix<T> basis_evals_tmp(n_new,m_1); 
						index_t nrp = n_new - (repeated_idx-p); 
						basis_evals_tmp.block(0,0,repeated_idx-p,m_1) = basis_evals.block(0,0,repeated_idx-p,m_1); 
						basis_evals_tmp.block(repeated_idx-p,0, nrp, m_1) = basis_evals.block(repeated_idx-p+1,0,nrp,m_1); 
					
						basis_evals = std::move(basis_evals_tmp); 
					
						compute_basis_evals(uk, m_1, kv_new, n_new, repeated_idx-p-1, repeated_idx-mult, basis_evals);
					
						//compute new error bounds for the relevant knots
						//Br won't necessarily get smaller as a knot with multiplicity 2 is still there with mult=1 after removing
						if (mult==1) //only one occurence of knot BEFORE removal -> knot is now gone
						{
							//then Br is now smaller
							//need unique_idx here, because Br_mat only exists for unique_interior knots
							Br_mat.removeCol(unique_idx-1); //should be the right knot
						
							//n_new is already updated! 
							matrix_removal_bounds(n_new,n_o,p,P_new,kv_new, dir, std::max(repeated_idx-p,p+1), std::min((size_t)(repeated_idx+p-mult+1),n_new-1), Br_mat);
							Br_vec.erase(Br_vec.begin() + unique_idx-1); 
							get_bound_surface(Br_mat, n_o, kv_new.uSize()-2, Br_vec, basis_evals_v);
						}
						else //more occurencies of knot -> knot is still there
						{
							matrix_removal_bounds(n_new,n_o,p,P_new,kv_new, dir, std::max(repeated_idx-p,p+1), std::min((size_t)(repeated_idx+p-mult+1),n_new-1), Br_mat);
							get_bound_surface(Br_mat, n_o, kv_new.uSize()-2, Br_vec, basis_evals_v);
						}
					}
				}
				else //knot can't be removed
				{
					Br_vec[unique_idx-1] = big_number; 
				}
			
			}//end else
		
		}//end of while-loop
	
		}
		
		return; 
	}

	template<class T>
	void gsTensorFitting<T>::get_max_in_matrix(const gsMatrix<T>& M, index_t& _row, index_t& _col, T& _value)
	{
		T max_val = -10000000000; 
	
		#pragma parallel for collapse(2) reduction(max:max_val) schedule(static)
		for (index_t i=0; i<M.rows(); i++)
		{
			for (index_t j=0; j<M.cols(); j++)
			{
				if (M(i,j)>=max_val)
				{
					max_val = M(i,j); 
					_row = i; 
					_col = j; 
				}
			}
		}
		
		_value = max_val; 
		
		return; 
	}
	
} //namespace gismo
