/** @file gsBemLaplace.h

    @brief Provides declaration of Basis abstract interface.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Falini, A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsForwardDeclarations.h>
#include <gsSolver/gsBemUtils.h>

namespace gismo
{
  /** 
      Bem laplace solver
  */


template<class T>
class gsBemLaplace
{

public:
    
    /// Default empty constructor
    gsBemLaplace() { };
    
    /// Constructor
    gsBemLaplace( gsCurveLoop<T> * loop, gsFunction<T> * boundary, gsBasis<T> * basis = NULL )
    {
        m_boundary_fun.push_back(boundary); // m_boundary_fun.push_back( boundary);
        m_pdomain = new gsPlanarDomain<T>(loop);
        m_basis.push_back(basis) ;
	    
	    // if ( basis)
	    //  m_basis.push_back(basis);
	    // else
	    // m_basis = & loop->singleCurve.basis();
	    //    
    }
    
    /// Constructor
    gsBemLaplace( gsPlanarDomain<T> * pl, std::vector<gsFunction<T> *> boundaries,
                  std::vector< gsBasis<T> *> basis);

    ~gsBemLaplace() 
     { 
         
     }
    
public:
    
    gsBemLaplace * clone() const
    { return new gsBemLaplace(*this); }
    
    /// Solves the laplace problem and returns a pointer to the solution
    /// \param parametric_bc If set to true, the boundary conditions are
    /// considered to be expressed in the parameter domain
    gsBemSolution<T> * solve(bool parametric_bc = false);
    
    //gsField<T> * outerBoundary();
    
    // gsBasis<T>    & basis()  { return *m_basis; }
    
// Data members
private:
    
    /// The green function
    gsGreenFunction2d<T> m_green_fun;
    
    std::vector<gsFunction<T> *>  m_boundary_fun;
    gsPlanarDomain<T> *m_pdomain;
    std::vector< gsBasis<T> *> m_basis;

public:
    // Needed since m_green_fun is 16B aligned
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    
}; // class gsClass
    
    
//////////////////////////////////////////////////
//////////////////////////////////////////////////
   
    
} // namespace gismo


#ifndef GISMO_HEADERS_ONLY
#include GISMO_HPP_HEADER(gsBemLaplace.hpp)
#endif
