/** @file gsAssemblerOptions.h

    @brief Provides assembler and solver options.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris
*/

#pragma once

namespace gismo
{
struct dirichlet
{	
	enum values
	{
        homogeneous   = 1, ///< Assume homogeneous Dirichlet conditions

        interpolation = 2, ///< Compute Dirichlet DoFs by using interpolation on the boundary
        
        l2Projection  = 3, ///< Compute Dirichlet DoFs by using L2 projection on the boundary
        
        user          = 10 ///< User will provide values of the Dirichlet dofs
    };

	enum strategy
	{
        elimination  = 11, ///< Compute Dirichlet DoFs by using interpolation on the boundary

        nitsche      = 12, ///< Enforce the boundary condition weakly by a penalty term
	
        penalize     = 13, ///< Penalize the diagonal at the position of Dirichlet DoFs,

        //nullifyrow    = 14,
        
        /// Compute Dirichlet DoFs in the normal direction (for a vector valued function),
        /// The tangential component are handled with the Nitsche method.
        eliminatNormal = 14,

        none         = 0 ///<< Do absolutely nothing for Dirichlet boundary conditions.
	};
};

struct iFace
{	
	enum strategy
	{
        /// Glue patches together by merging DoFs across an
        /// interface into one. This only works for conforming
	    /// interfaces.
	    conforming = 1,
        glue       = 1,
	    
	    /// Use discontinuous Galerkin-like coupling between
	    /// adjacent patches.
	    dg = 2,

	    /// Do absolutely nothing for coupling the interfaces.
	    none = 0
	};
};

struct transform
{	
	enum type
	{
	    Hgrad = 1, // covariant, inverse_composition
	    Hdiv  = 2, // Piola 
	    Hcurl = 3
	};
};

// for mixed formulations
struct space
{	
	enum type
	{
	    taylorHood    = 1,
	    raviartThomas = 2,

	    none          = 0
	};
};



struct gsAssemblerOptions
{
public:
    // Default constructor
    gsAssemblerOptions()
    : dirValues    (dirichlet::l2Projection ), 
      dirStrategy  (dirichlet::nitsche      ), 
      intStrategy  (iFace    ::glue         ),
      transformType(transform::Hgrad        ),
      spaceType    (space    ::taylorHood   )
    { }

public:

    dirichlet::values  dirValues;

    dirichlet::strategy  dirStrategy;

    iFace::strategy      intStrategy;

    transform::type      transformType;
    
    space::type          spaceType;

};

} // namespace gismo
