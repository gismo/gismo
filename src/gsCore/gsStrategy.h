/** @file gsStrategy.h

    @brief Discretization strategy enumerations (Dirichlet enforcement,
    interface coupling, transforms, mixed spaces), shared by bases,
    assemblers and boundary-condition handling.

    Moved here from gsAssembler/gsAssemblerOptions.h: these enums are used
    inline by core classes (e.g. gsMultiBasis) and are conceptually core
    discretization strategies, not assembler options.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): A. Mantzaflaris, C. Hofer
*/

#pragma once

namespace gismo
{

struct dirichlet
{
    enum strategy
    {
        elimination  = 11, ///< Enforce Dirichlet BCs by eliminating them from the system

        penalize     = 13, ///< Penalize the diagonal at the position of Dirichlet DoFs,

        nitsche      = 12, ///< Enforce the boundary condition weakly by a penalty term

        /// Compute Dirichlet DoFs in the normal direction (for a vector valued function),
        /// The tangential component are handled with the Nitsche method.
        eliminatNormal = 14,

        none         = 0 ///<< Do absolutely nothing for Dirichlet boundary conditions.
    };

    enum values
    {
        homogeneous   = 100, ///< Assume homogeneous Dirichlet conditions

        interpolation = 101, ///< Compute Dirichlet DoFs by using interpolation on the boundary

        l2Projection  = 102, ///< Compute Dirichlet DoFs by using L2 projection on the boundary

        user          = 103 ///< User will provide values of the Dirichlet dofs
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

        /// Use enhanced smoothness splines between interfaces of adjacent patches.
        smooth = 3,

        /// Do absolutely nothing for coupling the interfaces.
        none = 0
    };

};

/*
    enum iFaceTopology
    {
        nested   = 1,

        clamped  = 2,
    }
*/

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
struct discreteSpace
{
    enum type
    {
        taylorHood    = 1,
        //instead of raviartThomas, there should be nested_Space and BubbleElement
        //raviartThomas should go away here.
        raviartThomas = 2,

        none          = 0
    };
};

} // namespace gismo
