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
	enum strategy
	{
        /// Apply homogeneous dirichlet conditions
        homogeneous  = 1,

        /// Compute Dirichlet dofs by using interpolation on the boundary
        interpolation  = 2,
        elimination  = 2,

        /// Compute Dirichlet dofs by using least-squares fitting on the boundary
        leastSquares = 3,

        /// Compute Dirichlet dofs by using least-squares fitting on the boundary
        l2Projection = 4,

        /// Fixed values provided by the user.
        fixed        = 5,

        /// Enforce the boundary condition weakly by a penalty
        /// term. Not compatible/ignores eliminate==true
        nitsche      = 11,
	
        /// Penalize the diagonal at the position of Dirichlet dofs,
        /// Not compatible/ignores eliminate==true 
        penalize     = 12,
        

        /// Compute Dirichlet dofs in the normal direction (for a vector valued function),
        /// The tangential component are handled with the nitche method.
        eliminatNormal = 17,

        /// Do absolutely nothing for Dirichlet boundary conditions.
        none         = 0
	};

	/* ///If true, the dirichet dofs are part of the system Dofs,
	   /// else they are eliminated a priori from the system
	enum eliminate
	{
	    no  = 0,
	    yes = 1
	};
    */
};

struct iFace
{	
	enum strategy
	{
        /// Glue patches together by merging dofs across an
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

// To do: add more options 
struct gsAssemblerOptions
{
public:
    gsAssemblerOptions()
    : dirStrategy(dirichlet::nitsche), intStrategy(iFace::glue)
    { }

    gsAssemblerOptions(dirichlet::strategy _dirStrategy,
                       iFace::strategy     _intStrategy)
    : dirStrategy(_dirStrategy), intStrategy(_intStrategy)
    { }

public:

    dirichlet::strategy  dirStrategy;

    iFace::strategy      intStrategy;
};

} // namespace gismo
