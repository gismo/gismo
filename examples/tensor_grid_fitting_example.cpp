/** @file tensor_grid_fitting_example.cpp

    @brief Demonstrates the use of the gsTensorFitting-class

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): J. Mayr
*/


#include <gismo.h>

using namespace gismo;

int main(int argc, char *argv[])
{
	
	//Options with default parameters
	double E_bound = 0.3; //0.3; 
	int p = 4; 
	
	std::string fn = "fitting/grid_surface.xml";
	
	
	//! [Read data]
    // Tensor grid surface fitting
    // Expected input is a file with matrices with:
    // id 1:  Q_t     -- tensor grid of data points, size 3xN, where N = n*m
    // 					a column at position i*n+j in Q_t corresponds to a 
	//					data point Q(i,j) in an nxm-matrix Q containing a tensor grid of data points
    // id 0:  [n,m] -- number of rows and columns in matrix 
    gsFileData<> fd_in(fn);
    gsMatrix<> Q; 
	gsMatrix<> n_coords(2,1); 
    fd_in.getId<gsMatrix<> >(0, n_coords );
    fd_in.getId<gsMatrix<> >(1, Q);
    //! [Read data]
	
	gsTensorFitting<real_t> tensor_fit(Q, n_coords(0,0), n_coords(1,0), p, E_bound); 
	
	tensor_fit.compute(); 
	
	gsTensorBSpline<2,real_t> surf = tensor_fit.result(); 
	
	gsInfo << "result: " << surf << "\n"; 
	
	gsWriteParaview(surf,"tensor_grid_fitting_result"); 
	
	
	return 0; 
}	
