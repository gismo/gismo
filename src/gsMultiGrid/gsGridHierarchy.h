/** @file gsGridHierarchy.h

    @brief Coarsening algorithms for knot vectors and bases.

    This file is part of the G+Smo library.

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.

    Author(s): C. Hofreither, S. Takacs
*/

#pragma once

#include <gsCore/gsDofMapper.h>
#include <gsNurbs/gsKnotVector.h>
#include <gsIO/gsOptionList.h>

namespace gismo
{

/** @brief
 *  Grid Hierarchy
 *
 *  This class allows to construct a grid hierarchy and stors a grid hierarchy (vector of
 *  bases, local transfer matrices and transfer matrices). The transfer matrices can then
 *  be used for the setup of a \a gsMultiGridOp object.
 *
 *  @ingroup Solver
**/
template< typename T >
class gsGridHierarchy
{

public:

    /// @brief This function sets up a multigrid hierarchy by uniform refinement
    ///
    /// @param mBasis                    The gsMultiBasis to be refined (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param assemblerOptions          A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
    /// @param levels                    The number of levels
    /// @param numberOfKnotsToBeInserted The number of knots to be inserted, defaulted to 1
    /// @param multiplicityOfKnotsToBeInserted The multiplicity of the knots to be inserted, defaulted to 1
    /// @param unk                       Since the gsBoundaryCondition object can obtain data for systems
    ///                                  of PDEs, we have to provide information concerning which unknown
    ///                                  we are refering to.
    static gsGridHierarchy buildByRefinement(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& assemblerOptions,
        index_t levels,
        index_t numberOfKnotsToBeInserted = 1,
        index_t multiplicityOfKnotsToBeInserted = 1,
        index_t unk = 0
        );

    /// @brief This function sets up a multigrid hierarchy by uniform refinement
    ///
    /// @param mBasis                    The gsMultiBasis to be refined (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param options                   A gsOptionList defining the necessary infomation
    /// Takes boundary conditions for unknown 0.
    static gsGridHierarchy buildByRefinement(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& options
        )
    {
        return gsGridHierarchy::buildByRefinement(
            give(mBasis),
            boundaryConditions,
            options,
            options.askInt( "Levels", 3 ),
            options.askInt( "NumberOfKnotsToBeInserted", 1 ),
            options.askInt( "MultiplicityOfKnotsToBeInserted", 1 ),
            0
        );
    }


    /// @brief This function sets up a grid hierarchy by coarsening
    ///
    /// @param mBasis                    The gsMultiBasis to be coarsened (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param assemblerOptions          A gsOptionList defining a "DirichletStrategy" and a "InterfaceStrategy"
    /// @param levels                    The maximum number of levels
    /// @param degreesOfFreedom          Number of dofs in the coarsest grid in the grid hierarchy
    /// @param unk                       Since the gsBoundaryCondition object can obtain data for systems
    ///                                  of PDEs, we have to provide information concerning which unknown
    ///                                  we are refering to.
    ///
    /// The algorithm terminates if either the number of levels is reached or the number of degrees of freedom
    /// is below the given threshold.
    ///
    static gsGridHierarchy buildByCoarsening(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& assemblerOptions,
        index_t levels,
        index_t degreesOfFreedom = 0,
        index_t unk = 0
        );

    /// @brief This function sets up a grid hierarchy by coarsening
    ///
    /// @param mBasis                    The gsMultiBasis to be coarsened (initial basis)
    /// @param boundaryConditions        The boundary conditions
    /// @param options                   A gsOptionList defining the necessary infomation
    ///
    /// Takes boundary conditions for unkown 0.
    static gsGridHierarchy buildByCoarsening(
        gsMultiBasis<T> mBasis,
        const gsBoundaryConditions<T>& boundaryConditions,
        const gsOptionList& options
    )
    {
        return gsGridHierarchy::buildByCoarsening(
            give(mBasis),
            boundaryConditions,
            options,
            options.askInt( "Levels", 3 ),
            options.askInt( "DegreesOfFreedom", 0 ),
            0
        );
    }

    /// Get the default options
    static gsOptionList defaultOptions()
    {
        gsOptionList opt;
        opt.addInt( "DirichletStrategy", "Method for enforcement of Dirichlet BCs [11..14]", 11 );
        opt.addInt( "InterfaceStrategy", "Method of treatment of patch interfaces [0..3]", 1  );
        opt.addInt( "Levels", "Number of levels to be constructed in the grid hierarchy", 3 );
        opt.addInt( "DegreesOfFreedom",   "Number of dofs in the coarsest grid in the grid hierarchy (only buildByCoarsening)", 0 );
        opt.addInt( "NumberOfKnotsToBeInserted", "The number of knots to be inserted (only buildByRefinement)", 1 );
        opt.addInt( "MultiplicityOfKnotsToBeInserted",   "The multiplicity of the knots to be inserted (only buildByRefinement)", 1 );
        return opt;
    }

    /// Reset the object (to save memory)
    void clear()
    {
        m_mBases.clear();
        m_transferMatrices.clear();
    }

    /// Get the vector of \a gsMultiBasis objects (by reference)
    const std::vector< gsMultiBasis<T> >& getMultiBases() const
    { return m_mBases; }

    /// Get the vector of \a gsMultiBasis objects
    gsGridHierarchy& moveMultiBasesTo( std::vector< gsMultiBasis<T> >& o )
    { o = give(m_mBases); return *this; }

    /// Get the vector of transfer matrices (by reference)
    const std::vector< gsSparseMatrix<T, RowMajor> >& getTransferMatrices() const
    { return m_transferMatrices; }

    /// Get the vector of transfer matrices
    gsGridHierarchy& moveTransferMatricesTo( std::vector< gsSparseMatrix<T, RowMajor> >& o )
    { o = give(m_transferMatrices); return *this; }

private:
    std::vector< gsMultiBasis<T> > m_mBases;
    std::vector< gsSparseMatrix<T, RowMajor> > m_transferMatrices;
};

} // namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsGridHierarchy.hpp)
#endif
